import os
import re
import csv
import math
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =========================================================
# CONFIG
# =========================================================

GUIDE_LEN = 20
MAX_MISMATCH = 3

TARGET_MULTI_DIR = os.path.join(os.getcwd(), "data", "targets_multistrain")
print("Looking for FASTA files in:", TARGET_MULTI_DIR)
PLASMID_DIR = "data/plasmids"
GENOME_DIR = "data/genomes"
RESULTS_DIR = "results_panstrain"

GENE_FILES = {
    "blaKPC": "blaKPC_multi.fasta",
    "blaNDM1": "blaNDM1_multi.fasta",
    "mcr1": "mcr1_multi.fasta",
    "mecA": "mecA_multi.fasta",
}

GENE_META = {
    "blaKPC":  {"label": "blaKPC",  "resistance": "Carbapenem"},
    "blaNDM1": {"label": "blaNDM1", "resistance": "Carbapenem"},
    "mcr1":    {"label": "mcr1",    "resistance": "Colistin"},
    "mecA":    {"label": "mecA",    "resistance": "Methicillin / MRSA"},
}


# =========================================================
# BASIC HELPERS
# =========================================================

def reverse_complement(seq):
    table = str.maketrans("ATCGNatcgn", "TAGCNtagcn")
    return seq.translate(table)[::-1]


def parse_fasta(filepath):
    records = []
    header = None
    seq_parts = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts).upper()))

    return records


def load_all_fasta_from_dir(directory):
    seqs = []
    if not os.path.isdir(directory):
        return seqs

    for fname in os.listdir(directory):
        if fname.lower().endswith((".fa", ".fasta", ".fna")):
            path = os.path.join(directory, fname)
            for header, seq in parse_fasta(path):
                seqs.append({
                    "source_file": fname,
                    "header": header,
                    "sequence": seq
                })
    return seqs


def gc_content(seq):
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return round(100.0 * gc / len(seq), 2)


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


# =========================================================
# PAM / GUIDE DISCOVERY
# =========================================================

def find_spcas9_guides(sequence, guide_len=20):
    """
    Find NGG PAMs on both strands from a reference sequence.
    """
    results = []
    seq = sequence.upper()
    rc = reverse_complement(seq)
    n = len(seq)

    # forward strand
    for i in range(n - guide_len - 2):
        spacer = seq[i:i + guide_len]
        pam = seq[i + guide_len:i + guide_len + 3]
        if re.fullmatch(r"[ATCG]GG", pam) and "N" not in spacer:
            results.append({
                "position": i,
                "strand": "+",
                "spacer": spacer,
                "pam": pam
            })

    # reverse strand
    for i in range(len(rc) - guide_len - 2):
        spacer = rc[i:i + guide_len]
        pam = rc[i + guide_len:i + guide_len + 3]
        if re.fullmatch(r"[ATCG]GG", pam) and "N" not in spacer:
            genomic_pos = n - (i + guide_len + 3)
            results.append({
                "position": genomic_pos,
                "strand": "-",
                "spacer": spacer,
                "pam": pam
            })

    return results


# =========================================================
# ON-TARGET SCORING
# =========================================================

def score_on_target(spacer, pam):
    score = 50.0
    gc = gc_content(spacer)

    if 40 <= gc <= 70:
        score += 18
    elif 30 <= gc < 40 or 70 < gc <= 80:
        score += 6
    else:
        score -= 12

    if "TTTT" in spacer:
        score -= 20
    elif "TTT" in spacer:
        score -= 8

    if "GGGG" in spacer:
        score -= 8

    if spacer[-1] == "G":
        score += 4

    if spacer[0] != "T":
        score += 2

    pam_bonus = {
        "TGG": 10,
        "AGG": 8,
        "CGG": 8,
        "GGG": 5
    }
    score += pam_bonus.get(pam, 4)

    for nuc in "ATCG":
        if nuc * 5 in spacer:
            score -= 10

    return max(0, min(100, round(score, 1)))


# =========================================================
# OFF-TARGET SCREENING
# =========================================================

def count_near_matches_in_sequence(guide, sequence, max_mismatch=3):
    counts = {0: 0, 1: 0, 2: 0, 3: 0}
    seq = sequence.upper()
    rc = reverse_complement(seq)
    L = len(guide)

    for i in range(len(seq) - L - 2):
        candidate = seq[i:i + L]
        pam = seq[i + L:i + L + 3]
        if not re.fullmatch(r"[ATCG]GG", pam):
            continue
        d = hamming_distance(guide, candidate)
        if d <= max_mismatch:
            counts[d] += 1

    for i in range(len(rc) - L - 2):
        candidate = rc[i:i + L]
        pam = rc[i + L:i + L + 3]
        if not re.fullmatch(r"[ATCG]GG", pam):
            continue
        d = hamming_distance(guide, candidate)
        if d <= max_mismatch:
            counts[d] += 1

    return counts


def aggregate_offtarget_counts(guide, background_records, max_mismatch=3):
    total = {0: 0, 1: 0, 2: 0, 3: 0}
    for rec in background_records:
        counts = count_near_matches_in_sequence(guide, rec["sequence"], max_mismatch)
        for k in total:
            total[k] += counts[k]
    return total


def offtarget_penalty(hit_counts):
    penalty = (
        hit_counts[0] * 50 +
        hit_counts[1] * 20 +
        hit_counts[2] * 8 +
        hit_counts[3] * 3
    )
    return round(penalty, 1)


def specificity_score_from_penalty(penalty):
    score = 100 * math.exp(-penalty / 50.0)
    return round(max(0, min(100, score)), 1)


# =========================================================
# CONSERVATION / PAN-STRAIN SCORING
# =========================================================

def scan_guide_in_sequence(guide, sequence, max_mismatch=1):
    """
    Returns best hit for guide in one strain sequence, requiring NGG PAM.
    Checks both strands. Returns dict or None.
    """
    seq = sequence.upper()
    rc = reverse_complement(seq)
    L = len(guide)

    best = None

    # forward
    for i in range(len(seq) - L - 2):
        candidate = seq[i:i + L]
        pam = seq[i + L:i + L + 3]
        if not re.fullmatch(r"[ATCG]GG", pam):
            continue
        d = hamming_distance(guide, candidate)
        if d <= max_mismatch:
            record = {
                "mismatches": d,
                "pam": pam,
                "strand": "+",
                "position": i
            }
            if best is None or d < best["mismatches"]:
                best = record

    # reverse
    for i in range(len(rc) - L - 2):
        candidate = rc[i:i + L]
        pam = rc[i + L:i + L + 3]
        if not re.fullmatch(r"[ATCG]GG", pam):
            continue
        d = hamming_distance(guide, candidate)
        if d <= max_mismatch:
            genomic_pos = len(seq) - (i + L + 3)
            record = {
                "mismatches": d,
                "pam": pam,
                "strand": "-",
                "position": genomic_pos
            }
            if best is None or d < best["mismatches"]:
                best = record

    return best


def conservation_profile(guide, strain_records, max_mismatch=1):
    """
    Evaluate guide across all strains of a gene.
    """
    total = len(strain_records)
    perfect = 0
    one_mm_or_better = 0
    pam_supported = 0

    heatmap_values = []  # 1=perfect, 0.5=1mm, 0=no hit

    for header, seq in strain_records:
        hit = scan_guide_in_sequence(guide, seq, max_mismatch=max_mismatch)
        if hit is None:
            heatmap_values.append(0.0)
            continue

        pam_supported += 1

        if hit["mismatches"] == 0:
            perfect += 1
            one_mm_or_better += 1
            heatmap_values.append(1.0)
        elif hit["mismatches"] == 1:
            one_mm_or_better += 1
            heatmap_values.append(0.5)

    perfect_frac = perfect / total if total else 0
    one_mm_frac = one_mm_or_better / total if total else 0
    pam_frac = pam_supported / total if total else 0

    # weighted conservation score
    conservation_score = (
        70 * perfect_frac +
        20 * max(0, one_mm_frac - perfect_frac) +
        10 * pam_frac
    )
    conservation_score = round(min(100, conservation_score), 1)

    return {
        "total_strains": total,
        "perfect_match_strains": perfect,
        "one_mismatch_or_better_strains": one_mm_or_better,
        "pam_supported_strains": pam_supported,
        "perfect_match_fraction": round(perfect_frac, 4),
        "one_mismatch_or_better_fraction": round(one_mm_frac, 4),
        "pam_supported_fraction": round(pam_frac, 4),
        "conservation_score": conservation_score,
        "heatmap_values": heatmap_values
    }


# =========================================================
# BACKGROUND LOADING
# =========================================================

def load_background_sequences():
    background = []
    for folder in [PLASMID_DIR, GENOME_DIR]:
        seqs = load_all_fasta_from_dir(folder)
        for rec in seqs:
            rec["source_dir"] = folder
            background.append(rec)
    return background


# =========================================================
# ANALYSIS
# =========================================================

def classify_final_score(score):
    if score >= 80:
        return "Excellent"
    elif score >= 65:
        return "Good"
    elif score >= 50:
        return "Moderate"
    else:
        return "Poor"


def analyze_multistrain_gene(gene_name, fasta_path, background_records, top_n=20):
    strain_records = parse_fasta(fasta_path)
    if not strain_records:
        return None

    ref_header, ref_sequence = strain_records[0]
    guides = find_spcas9_guides(ref_sequence, guide_len=GUIDE_LEN)

    print(f"\nAnalyzing {gene_name}")
    print(f"Reference header: {ref_header[:70]}")
    print(f"Number of strains: {len(strain_records)}")
    print(f"Reference guide candidates: {len(guides)}")

    results = []

    for g in guides:
        spacer = g["spacer"]
        pam = g["pam"]

        on_target = score_on_target(spacer, pam)

        off_counts = aggregate_offtarget_counts(
            spacer,
            background_records,
            max_mismatch=MAX_MISMATCH
        )

        penalty = offtarget_penalty(off_counts)
        specificity = specificity_score_from_penalty(penalty)

        cons = conservation_profile(spacer, strain_records, max_mismatch=1)

        # integrated paper-style score
        final_score = round(
            0.45 * on_target +
            0.30 * specificity +
            0.25 * cons["conservation_score"],
            1
        )

        results.append({
            "gene": gene_name,
            "position": g["position"],
            "strand": g["strand"],
            "spacer": spacer,
            "pam": pam,
            "gc_content": gc_content(spacer),
            "on_target_score": on_target,
            "offtarget_hits_0mm": off_counts[0],
            "offtarget_hits_1mm": off_counts[1],
            "offtarget_hits_2mm": off_counts[2],
            "offtarget_hits_3mm": off_counts[3],
            "offtarget_penalty": penalty,
            "specificity_score": specificity,
            "total_strains": cons["total_strains"],
            "perfect_match_strains": cons["perfect_match_strains"],
            "one_mismatch_or_better_strains": cons["one_mismatch_or_better_strains"],
            "pam_supported_strains": cons["pam_supported_strains"],
            "perfect_match_fraction": cons["perfect_match_fraction"],
            "one_mismatch_or_better_fraction": cons["one_mismatch_or_better_fraction"],
            "pam_supported_fraction": cons["pam_supported_fraction"],
            "conservation_score": cons["conservation_score"],
            "final_score": final_score,
            "classification": classify_final_score(final_score),
            "heatmap_values": cons["heatmap_values"]
        })

    results.sort(key=lambda x: x["final_score"], reverse=True)

    stats = {
        "total_guides": len(results),
        "excellent": sum(1 for r in results if r["classification"] == "Excellent"),
        "good": sum(1 for r in results if r["classification"] == "Good"),
        "moderate": sum(1 for r in results if r["classification"] == "Moderate"),
        "poor": sum(1 for r in results if r["classification"] == "Poor"),
        "best_score": max((r["final_score"] for r in results), default=0),
        "mean_score": round(sum(r["final_score"] for r in results) / len(results), 1) if results else 0,
        "best_conservation": max((r["conservation_score"] for r in results), default=0),
        "mean_conservation": round(sum(r["conservation_score"] for r in results) / len(results), 1) if results else 0
    }

    print(f"Best final score: {stats['best_score']}")
    print(f"Mean final score: {stats['mean_score']}")
    print(f"Best conservation: {stats['best_conservation']}")

    return {
        "gene": gene_name,
        "strain_records": strain_records,
        "results": results,
        "top_results": results[:top_n],
        "stats": stats
    }


# =========================================================
# CSV EXPORT
# =========================================================

def save_csv_outputs(all_data):
    os.makedirs(RESULTS_DIR, exist_ok=True)

    fieldnames = [
        "gene", "position", "strand", "spacer", "pam", "gc_content",
        "on_target_score",
        "offtarget_hits_0mm", "offtarget_hits_1mm", "offtarget_hits_2mm", "offtarget_hits_3mm",
        "offtarget_penalty", "specificity_score",
        "total_strains", "perfect_match_strains", "one_mismatch_or_better_strains",
        "pam_supported_strains", "perfect_match_fraction",
        "one_mismatch_or_better_fraction", "pam_supported_fraction",
        "conservation_score", "final_score", "classification"
    ]

    # full results
    full_path = os.path.join(RESULTS_DIR, "all_panstrain_guide_candidates.csv")
    with open(full_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for gene in all_data:
            for row in all_data[gene]["results"]:
                clean_row = {k: v for k, v in row.items() if k in fieldnames}
                writer.writerow(clean_row)

    # top per gene
    top_path = os.path.join(RESULTS_DIR, "top20_per_gene_panstrain.csv")
    with open(top_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for gene in all_data:
            for row in all_data[gene]["top_results"]:
                clean_row = {k: v for k, v in row.items() if k in fieldnames}
                writer.writerow(clean_row)

    # top global
    global_rows = []
    for gene in all_data:
        global_rows.extend(all_data[gene]["results"])
    global_rows.sort(key=lambda x: x["final_score"], reverse=True)

    global_path = os.path.join(RESULTS_DIR, "top30_global_panstrain_guides.csv")
    with open(global_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in global_rows[:30]:
            clean_row = {k: v for k, v in row.items() if k in fieldnames}
            writer.writerow(clean_row)

    # summary
    summary_path = os.path.join(RESULTS_DIR, "summary_statistics_panstrain.csv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "gene", "resistance", "total_guides", "excellent", "good", "moderate", "poor",
            "best_score", "mean_score", "best_conservation", "mean_conservation"
        ])
        for gene in all_data:
            s = all_data[gene]["stats"]
            writer.writerow([
                gene,
                GENE_META.get(gene, {}).get("resistance", ""),
                s["total_guides"],
                s["excellent"],
                s["good"],
                s["moderate"],
                s["poor"],
                s["best_score"],
                s["mean_score"],
                s["best_conservation"],
                s["mean_conservation"]
            ])

    print(f"\nSaved: {full_path}")
    print(f"Saved: {top_path}")
    print(f"Saved: {global_path}")
    print(f"Saved: {summary_path}")


# =========================================================
# FIGURES
# =========================================================

def make_score_distribution(all_data):
    plt.figure(figsize=(10, 6))
    for gene in all_data:
        scores = [r["final_score"] for r in all_data[gene]["results"]]
        plt.hist(scores, bins=25, alpha=0.45, label=gene)
    plt.xlabel("Pan-Guide Final Score")
    plt.ylabel("Guide Count")
    plt.title("Pan-Strain Guide Score Distribution Across MDR Genes")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "Figure1_PanStrainScoreDistribution.png"), dpi=180)
    plt.close()


def make_top_global_plot(all_data):
    global_rows = []
    for gene in all_data:
        global_rows.extend(all_data[gene]["results"])
    global_rows.sort(key=lambda x: x["final_score"], reverse=True)
    top = global_rows[:20]

    if not top:
        return

    labels = [f"{r['gene']}:{r['position']}" for r in top][::-1]
    vals = [r["final_score"] for r in top][::-1]

    plt.figure(figsize=(12, 8))
    plt.barh(labels, vals)
    plt.xlabel("Final Pan-Guide Score")
    plt.ylabel("Guide")
    plt.title("Top 20 Pan-Strain Guides Across MDR Genes")
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "Figure2_TopPanStrainGuides_Global.png"), dpi=180)
    plt.close()


def make_gene_comparison_plot(all_data):
    genes = list(all_data.keys())
    mean_scores = [all_data[g]["stats"]["mean_score"] for g in genes]
    mean_cons = [all_data[g]["stats"]["mean_conservation"] for g in genes]

    x = range(len(genes))

    plt.figure(figsize=(10, 6))
    plt.bar([i - 0.2 for i in x], mean_scores, width=0.4, label="Mean Final Score")
    plt.bar([i + 0.2 for i in x], mean_cons, width=0.4, label="Mean Conservation")
    plt.xticks(list(x), genes)
    plt.ylabel("Score")
    plt.title("Mean Final Score vs Mean Conservation by Gene")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "Figure3_PanStrainGeneComparison.png"), dpi=180)
    plt.close()


def make_specificity_vs_conservation_plot(all_data):
    plt.figure(figsize=(10, 6))
    for gene in all_data:
        rows = all_data[gene]["top_results"][:20]
        xvals = [r["conservation_score"] for r in rows]
        yvals = [r["specificity_score"] for r in rows]
        plt.scatter(xvals, yvals, label=gene, alpha=0.7)

    plt.xlabel("Conservation Score")
    plt.ylabel("Specificity Score")
    plt.title("Specificity vs Conservation of Top Pan-Strain Guides")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "Figure4_Specificity_vs_Conservation.png"), dpi=180)
    plt.close()


def make_coverage_plot(all_data):
    genes = list(all_data.keys())
    frac90 = []
    frac100 = []

    for gene in genes:
        rows = all_data[gene]["results"]
        n = len(rows) if rows else 1
        n90 = sum(1 for r in rows if r["perfect_match_fraction"] >= 0.90)
        n100 = sum(1 for r in rows if r["perfect_match_fraction"] >= 1.00)
        frac90.append(round(100 * n90 / n, 1))
        frac100.append(round(100 * n100 / n, 1))

    x = range(len(genes))

    plt.figure(figsize=(10, 6))
    plt.bar([i - 0.2 for i in x], frac90, width=0.4, label="Guides covering ≥90% strains")
    plt.bar([i + 0.2 for i in x], frac100, width=0.4, label="Guides covering 100% strains")
    plt.xticks(list(x), genes)
    plt.ylabel("Percent of Guides")
    plt.title("Pan-Strain Guide Coverage by Gene")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "Figure5_GuideCoverageByGene.png"), dpi=180)
    plt.close()


def make_heatmap_for_gene(gene, gene_data, top_n=20):
    rows = gene_data["top_results"][:top_n]
    strain_records = gene_data["strain_records"]

    if not rows or not strain_records:
        return

    matrix = [r["heatmap_values"] for r in rows]
    labels = [f"{r['position']}:{r['spacer'][:8]}..." for r in rows]
    strain_labels = [h[:18] for h, _ in strain_records]

    plt.figure(figsize=(max(8, len(strain_labels) * 0.4), max(6, len(labels) * 0.35)))
    plt.imshow(matrix, aspect="auto", interpolation="nearest")
    plt.colorbar(label="Match quality (1=perfect, 0.5=1 mismatch, 0=no hit)")
    plt.yticks(range(len(labels)), labels)
    plt.xticks(range(len(strain_labels)), strain_labels, rotation=90)
    plt.xlabel("Strains")
    plt.ylabel("Top Guides")
    plt.title(f"{gene} Top Guide Conservation Heatmap")
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, f"Heatmap_{gene}_TopGuides.png"), dpi=180)
    plt.close()


def make_figures(all_data):
    os.makedirs(RESULTS_DIR, exist_ok=True)

    make_score_distribution(all_data)
    make_top_global_plot(all_data)
    make_gene_comparison_plot(all_data)
    make_specificity_vs_conservation_plot(all_data)
    make_coverage_plot(all_data)

    for gene in all_data:
        make_heatmap_for_gene(gene, all_data[gene], top_n=20)

    print(f"Saved pan-strain figures in: {RESULTS_DIR}/")


# =========================================================
# REPORT
# =========================================================

def print_report(all_data):
    print("\n" + "=" * 100)
    print("PAN-STRAIN CRISPR GUIDE ANALYSIS REPORT")
    print("=" * 100)
    print(f"{'Gene':<10} {'Resistance':<22} {'Guides':<10} {'Best':<8} {'Mean':<8} {'BestCons':<10} {'MeanCons'}")
    print("-" * 100)

    for gene in all_data:
        s = all_data[gene]["stats"]
        resistance = GENE_META.get(gene, {}).get("resistance", "")
        print(
            f"{gene:<10} {resistance:<22} {s['total_guides']:<10} {s['best_score']:<8} "
            f"{s['mean_score']:<8} {s['best_conservation']:<10} {s['mean_conservation']}"
        )
    print("=" * 100)


# =========================================================
# MAIN
# =========================================================

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    background_records = load_background_sequences()
    print(f"Loaded {len(background_records)} background FASTA records for off-target screening.")

    all_data = {}

    for gene, fname in GENE_FILES.items():
        fpath = os.path.join(TARGET_MULTI_DIR, fname)
        if not os.path.exists(fpath):
            print(f"Missing multi-strain FASTA: {fpath}")
            continue

        result = analyze_multistrain_gene(gene, fpath, background_records, top_n=20)
        if result:
            all_data[gene] = result

    if not all_data:
        print("No genes analyzed. Add multi-strain FASTA files first.")
        return

    save_csv_outputs(all_data)
    make_figures(all_data)
    print_report(all_data)

    print(f"\nDone. Open '{RESULTS_DIR}/' for pan-strain outputs.")


if __name__ == "__main__":
    main()
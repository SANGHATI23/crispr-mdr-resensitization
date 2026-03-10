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

TARGET_DIR = "data/targets"
PLASMID_DIR = "data/plasmids"
GENOME_DIR = "data/genomes"
RESULTS_DIR = "results"

GENE_META = {
    "blaKPC":  {"label": "blaKPC",  "resistance": "Carbapenem"},
    "blaNDM1": {"label": "blaNDM1", "resistance": "Carbapenem"},
    "mcr1":    {"label": "mcr1",    "resistance": "Colistin"},
    "mecA":    {"label": "mecA",    "resistance": "Methicillin / MRSA"},
}


# =========================================================
# BASIC SEQUENCE HELPERS
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
# PAM SCANNING
# =========================================================

def find_spcas9_guides(sequence, guide_len=20):
    """
    Find NGG PAMs on both strands.
    Returns list of dicts.
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
    """
    Count candidate off-target sites adjacent to SpCas9-like PAMs
    on both strands. Returns dict mismatch_count -> hits.
    """
    counts = {0: 0, 1: 0, 2: 0, 3: 0}
    seq = sequence.upper()
    rc = reverse_complement(seq)
    L = len(guide)

    # forward strand scanning
    for i in range(len(seq) - L - 2):
        candidate = seq[i:i + L]
        pam = seq[i + L:i + L + 3]
        if not re.fullmatch(r"[ATCG]GG", pam):
            continue
        d = hamming_distance(guide, candidate)
        if d <= max_mismatch:
            counts[d] += 1

    # reverse strand scanning
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
        seq = rec["sequence"]
        counts = count_near_matches_in_sequence(guide, seq, max_mismatch=max_mismatch)
        for k in total:
            total[k] += counts[k]

    return total


def offtarget_penalty(hit_counts):
    """
    Heavier penalty for exact/near-exact matches.
    """
    penalty = (
        hit_counts[0] * 50 +
        hit_counts[1] * 20 +
        hit_counts[2] * 8 +
        hit_counts[3] * 3
    )
    return round(penalty, 1)


def specificity_score_from_penalty(penalty):
    """
    Convert penalty into 0-100 specificity score.
    Higher specificity is better.
    """
    score = 100 * math.exp(-penalty / 50.0)
    return round(max(0, min(100, score)), 1)


# =========================================================
# BACKGROUND LOADING
# =========================================================

def load_background_sequences():
    background = []

    for folder in [TARGET_DIR, PLASMID_DIR, GENOME_DIR]:
        seqs = load_all_fasta_from_dir(folder)
        for rec in seqs:
            rec["source_dir"] = folder
            background.append(rec)

    return background


# =========================================================
# ANALYSIS
# =========================================================

def classify_final_score(score):
    if score >= 70:
        return "Excellent"
    elif score >= 55:
        return "Good"
    elif score >= 40:
        return "Moderate"
    else:
        return "Poor"


def analyze_gene(gene_name, fasta_path, background_records, top_n=10):
    records = parse_fasta(fasta_path)
    if not records:
        return None

    header, sequence = records[0]
    guides = find_spcas9_guides(sequence, guide_len=GUIDE_LEN)

    results = []

    for g in guides:
        spacer = g["spacer"]
        pam = g["pam"]

        on_target = score_on_target(spacer, pam)

        # off-target search across all loaded backgrounds
        hit_counts = aggregate_offtarget_counts(
            spacer,
            background_records,
            max_mismatch=MAX_MISMATCH
        )

        # subtract one intended perfect match if present
        if hit_counts[0] > 0:
            hit_counts[0] -= 1

        penalty = offtarget_penalty(hit_counts)
        specificity = specificity_score_from_penalty(penalty)

        final_score = round((0.6 * on_target) + (0.4 * specificity), 1)

        results.append({
            "gene": gene_name,
            "position": g["position"],
            "strand": g["strand"],
            "spacer": spacer,
            "pam": pam,
            "gc_content": gc_content(spacer),
            "on_target_score": on_target,
            "offtarget_hits_0mm": hit_counts[0],
            "offtarget_hits_1mm": hit_counts[1],
            "offtarget_hits_2mm": hit_counts[2],
            "offtarget_hits_3mm": hit_counts[3],
            "offtarget_penalty": penalty,
            "specificity_score": specificity,
            "final_score": final_score,
            "classification": classify_final_score(final_score)
        })

    results.sort(key=lambda x: x["final_score"], reverse=True)

    stats = {
        "total_guides": len(results),
        "excellent": sum(1 for r in results if r["classification"] == "Excellent"),
        "good": sum(1 for r in results if r["classification"] == "Good"),
        "moderate": sum(1 for r in results if r["classification"] == "Moderate"),
        "poor": sum(1 for r in results if r["classification"] == "Poor"),
        "best_score": max([r["final_score"] for r in results], default=0),
        "mean_score": round(sum(r["final_score"] for r in results) / len(results), 1) if results else 0
    }

    print(f"\nAnalyzing {gene_name}")
    print(f"Sequence header: {header[:60]}")
    print(f"Guide candidates: {len(results)}")
    print(f"Best score: {stats['best_score']}")
    print(f"Mean score: {stats['mean_score']}")

    return {
        "gene": gene_name,
        "results": results,
        "top_results": results[:top_n],
        "stats": stats
    }


# =========================================================
# CSV EXPORT
# =========================================================

def save_csv_outputs(all_data):
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # full results
    full_path = os.path.join(RESULTS_DIR, "all_guide_candidates.csv")
    fieldnames = [
        "gene", "position", "strand", "spacer", "pam", "gc_content",
        "on_target_score",
        "offtarget_hits_0mm", "offtarget_hits_1mm", "offtarget_hits_2mm", "offtarget_hits_3mm",
        "offtarget_penalty", "specificity_score", "final_score", "classification"
    ]

    with open(full_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for gene in all_data:
            for row in all_data[gene]["results"]:
                writer.writerow(row)

    # top 10 per gene
    top10_path = os.path.join(RESULTS_DIR, "top10_per_gene.csv")
    with open(top10_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for gene in all_data:
            for row in all_data[gene]["top_results"]:
                writer.writerow(row)

    # top 20 global
    global_rows = []
    for gene in all_data:
        global_rows.extend(all_data[gene]["results"])
    global_rows.sort(key=lambda x: x["final_score"], reverse=True)

    top20_global_path = os.path.join(RESULTS_DIR, "top20_global_guides.csv")
    with open(top20_global_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in global_rows[:20]:
            writer.writerow(row)

    # summary statistics
    summary_path = os.path.join(RESULTS_DIR, "summary_statistics.csv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "gene", "resistance", "total_guides",
            "excellent", "good", "moderate", "poor",
            "best_score", "mean_score"
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
                s["mean_score"]
            ])

    print(f"\nSaved: {full_path}")
    print(f"Saved: {top10_path}")
    print(f"Saved: {top20_global_path}")
    print(f"Saved: {summary_path}")


# =========================================================
# FIGURES
# =========================================================

def make_figures(all_data):
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Figure 1: score distributions
    plt.figure(figsize=(10, 6))
    for gene in all_data:
        scores = [r["final_score"] for r in all_data[gene]["results"]]
        plt.hist(scores, bins=25, alpha=0.45, label=gene)
    plt.xlabel("Final Score")
    plt.ylabel("Guide Count")
    plt.title("Final Guide Score Distribution Across MDR Genes")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "Figure1_ScoreDistribution.png"), dpi=180)
    plt.close()

    # Figure 2: top global guides
    global_rows = []
    for gene in all_data:
        global_rows.extend(all_data[gene]["results"])
    global_rows.sort(key=lambda x: x["final_score"], reverse=True)
    top = global_rows[:20]

    if top:
        labels = [f"{r['gene']}:{r['position']}" for r in top][::-1]
        vals = [r["final_score"] for r in top][::-1]
        plt.figure(figsize=(12, 8))
        plt.barh(labels, vals)
        plt.xlabel("Final Score")
        plt.ylabel("Guide")
        plt.title("Top 20 Guides Across All MDR Genes")
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, "Figure2_TopGuides_Global.png"), dpi=180)
        plt.close()

    # Figure 3: gene comparison
    genes = list(all_data.keys())
    mean_scores = [all_data[g]["stats"]["mean_score"] for g in genes]
    best_scores = [all_data[g]["stats"]["best_score"] for g in genes]

    x = range(len(genes))
    plt.figure(figsize=(10, 6))
    plt.bar([i - 0.2 for i in x], mean_scores, width=0.4, label="Mean Score")
    plt.bar([i + 0.2 for i in x], best_scores, width=0.4, label="Best Score")
    plt.xticks(list(x), genes)
    plt.ylabel("Score")
    plt.title("Mean vs Best Guide Scores by Gene")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "Figure3_GeneComparison.png"), dpi=180)
    plt.close()

    # Figure 4: on-target vs specificity
    plt.figure(figsize=(10, 6))
    for gene in all_data:
        rows = all_data[gene]["top_results"][:20]
        xvals = [r["specificity_score"] for r in rows]
        yvals = [r["on_target_score"] for r in rows]
        plt.scatter(xvals, yvals, label=gene, alpha=0.7)
    plt.xlabel("Specificity Score")
    plt.ylabel("On-Target Score")
    plt.title("On-Target Efficiency vs Specificity (Top Guides)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "Figure4_OnTarget_vs_OffTarget.png"), dpi=180)
    plt.close()

    print(f"Saved: {os.path.join(RESULTS_DIR, 'Figure1_ScoreDistribution.png')}")
    print(f"Saved: {os.path.join(RESULTS_DIR, 'Figure2_TopGuides_Global.png')}")
    print(f"Saved: {os.path.join(RESULTS_DIR, 'Figure3_GeneComparison.png')}")
    print(f"Saved: {os.path.join(RESULTS_DIR, 'Figure4_OnTarget_vs_OffTarget.png')}")


# =========================================================
# REPORT
# =========================================================

def print_report(all_data):
    print("\n" + "=" * 80)
    print("CRISPR MDR GUIDE ANALYSIS REPORT")
    print("=" * 80)
    print(f"{'Gene':<10} {'Resistance':<22} {'Guides':<10} {'Best':<8} {'Mean':<8} {'Excellent'}")
    print("-" * 80)

    for gene in all_data:
        s = all_data[gene]["stats"]
        resistance = GENE_META.get(gene, {}).get("resistance", "")
        print(f"{gene:<10} {resistance:<22} {s['total_guides']:<10} {s['best_score']:<8} {s['mean_score']:<8} {s['excellent']}")

    print("=" * 80)


# =========================================================
# MAIN
# =========================================================

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    target_files = {
        "blaKPC": os.path.join(TARGET_DIR, "blaKPC.fasta"),
        "blaNDM1": os.path.join(TARGET_DIR, "blaNDM1.fasta"),
        "mcr1": os.path.join(TARGET_DIR, "mcr1.fasta"),
        "mecA": os.path.join(TARGET_DIR, "mecA.fasta"),
    }

    background_records = load_background_sequences()
    print(f"Loaded {len(background_records)} background FASTA records for off-target screening.")

    all_data = {}

    for gene, path in target_files.items():
        if not os.path.exists(path):
            print(f"Missing target FASTA: {path}")
            continue

        result = analyze_gene(gene, path, background_records, top_n=10)
        if result:
            all_data[gene] = result

    if not all_data:
        print("No genes analyzed.")
        return

    save_csv_outputs(all_data)
    make_figures(all_data)
    print_report(all_data)

    print(f"\nDone. Open '{RESULTS_DIR}/' for the upgraded outputs.")


if __name__ == "__main__":
    main()
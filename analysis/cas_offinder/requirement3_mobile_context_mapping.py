#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

HIT_TABLE = PROJECT_ROOT / "results_cas_offinder/expanded_panel/functional_annotation/final_requirement2_database_refined/requirement2_database_refined_hit_annotations.csv"
GFF_DIR = PROJECT_ROOT / "data_cas_offinder/expanded_panel/annotations"
OUT_DIR = PROJECT_ROOT / "results_cas_offinder/expanded_panel/mobile_context_mapping"

OUT_DIR.mkdir(parents=True, exist_ok=True)

WINDOW_BP = 10000


MOBILE_TERMS = [
    "transposase", "transposon", "insertion sequence", "insertion element",
    "integrase", "recombinase", "resolvase", "invertase",
    "mobilization", "mobilisation", "relaxase", "conjugal",
    "conjugative", "tn3", "tn7", "tn4401", "is family"
]

PLASMID_TERMS = [
    "plasmid", "replication initiation protein", "replicase",
    "partition protein", "relaxase", "conjugal transfer",
    "conjugative transfer", "mobilization protein"
]

INTEGRON_TERMS = [
    "integron", "inti1", "inti2", "inti3", "gene cassette", "attc"
]

SCCMEC_TERMS = [
    "sccmec", "staphylococcal cassette chromosome",
    "cassette chromosome", "meca", "mecr1", "meci",
    "ccra", "ccrb", "ccrc", "orfX".lower()
]

AMR_TERMS = [
    "resistance", "antibiotic", "beta-lactamase", "carbapenemase",
    "multidrug", "efflux", "colistin", "polymyxin",
    "mcr", "bla", "ndm", "kpc", "meca", "aac", "aph", "aad",
    "erm", "tet", "sul", "qnr"
]


def parse_gff_attributes(attr):
    d = {}
    if pd.isna(attr):
        return d

    for item in str(attr).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip().replace("%20", " ")

    return d


def extract_accession(path):
    m = re.search(r"(GC[AF]_\d+\.\d+)", path.name)
    if m:
        return m.group(1)
    return path.stem


def contains_any(text, terms):
    text = str(text).lower()
    return any(t in text for t in terms)


def infer_seqid(text):
    if pd.isna(text):
        return ""
    return str(text).replace(">", "").split()[0]


def pick_col(df, candidates, required=True):
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]

    if required:
        raise ValueError(
            f"Missing required column. Tried {candidates}. Available: {list(df.columns)}"
        )
    return None


def load_gff_folder():
    files = list(GFF_DIR.glob("*.gff")) + list(GFF_DIR.glob("*.gff3"))

    if not files:
        raise FileNotFoundError(f"No GFF/GFF3 files found in {GFF_DIR}")

    rows = []

    for gff_file in files:
        accession = extract_accession(gff_file)

        with open(gff_file, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue

                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9:
                    continue

                seqid, source, feature_type, start, end, score, strand, phase, attrs = parts

                try:
                    start = int(start)
                    end = int(end)
                except ValueError:
                    continue

                parsed = parse_gff_attributes(attrs)

                gene = (
                    parsed.get("gene")
                    or parsed.get("Name")
                    or parsed.get("locus_tag")
                    or parsed.get("ID")
                    or ""
                )

                product = (
                    parsed.get("product")
                    or parsed.get("Note")
                    or parsed.get("function")
                    or ""
                )

                rows.append({
                    "accession": accession,
                    "seqid": seqid,
                    "feature_type": feature_type,
                    "feature_start": start,
                    "feature_end": end,
                    "feature_strand": strand,
                    "feature_gene": gene,
                    "feature_product": product,
                    "feature_attributes": attrs,
                    "gff_file": str(gff_file.relative_to(PROJECT_ROOT))
                })

    gff = pd.DataFrame(rows)

    if gff.empty:
        raise RuntimeError("GFF files loaded, but no feature rows parsed.")

    gff["feature_text"] = (
        gff["feature_gene"].fillna("").astype(str)
        + " "
        + gff["feature_product"].fillna("").astype(str)
        + " "
        + gff["feature_attributes"].fillna("").astype(str)
    ).str.lower()

    return gff


def classify_window(window):
    if window.empty:
        return {
            "nearby_feature_count": 0,
            "nearby_cds_count": 0,
            "nearby_mobile_feature_count": 0,
            "nearby_plasmid_feature_count": 0,
            "nearby_integron_feature_count": 0,
            "nearby_sccmec_feature_count": 0,
            "nearby_amr_feature_count": 0,
            "mobile_context_flag": 0,
            "plasmid_context_flag": 0,
            "integron_context_flag": 0,
            "sccmec_context_flag": 0,
            "resistance_island_context_flag": 0,
            "operon_proxy_flag": 0,
            "context_class": "No_nearby_GFF_features",
            "nearby_gene_products": ""
        }

    text = window["feature_text"].fillna("").astype(str)

    mobile = text.apply(lambda x: contains_any(x, MOBILE_TERMS))
    plasmid = text.apply(lambda x: contains_any(x, PLASMID_TERMS))
    integron = text.apply(lambda x: contains_any(x, INTEGRON_TERMS))
    sccmec = text.apply(lambda x: contains_any(x, SCCMEC_TERMS))
    amr = text.apply(lambda x: contains_any(x, AMR_TERMS))

    cds = window["feature_type"].astype(str).str.lower().eq("cds")

    mobile_count = int(mobile.sum())
    plasmid_count = int(plasmid.sum())
    integron_count = int(integron.sum())
    sccmec_count = int(sccmec.sum())
    amr_count = int(amr.sum())

    resistance_island = int(amr_count > 0 and (mobile_count > 0 or integron_count > 0))

    operon_proxy = 0
    cds_window = window[cds].copy()

    if len(cds_window) >= 2:
        same_strand_max = cds_window["feature_strand"].value_counts().max()
        span = int(cds_window["feature_end"].max() - cds_window["feature_start"].min())

        if same_strand_max >= 2 and span <= 5000:
            operon_proxy = 1

    labels = []

    if mobile_count > 0:
        labels.append("Mobile_element_context")
    if plasmid_count > 0:
        labels.append("Plasmid_associated_context")
    if integron_count > 0:
        labels.append("Integron_like_context")
    if sccmec_count > 0:
        labels.append("SCCmec_like_context")
    if resistance_island:
        labels.append("Resistance_island_like_context")
    if operon_proxy:
        labels.append("Operon_proxy_context")

    if not labels:
        labels.append("No_dedicated_mobile_context_signal")

    nearby_products = (
        window[["feature_gene", "feature_product", "feature_type"]]
        .fillna("")
        .astype(str)
        .drop_duplicates()
        .head(12)
        .apply(lambda r: f"{r['feature_gene']}|{r['feature_product']}|{r['feature_type']}", axis=1)
        .tolist()
    )

    return {
        "nearby_feature_count": int(len(window)),
        "nearby_cds_count": int(cds.sum()),
        "nearby_mobile_feature_count": mobile_count,
        "nearby_plasmid_feature_count": plasmid_count,
        "nearby_integron_feature_count": integron_count,
        "nearby_sccmec_feature_count": sccmec_count,
        "nearby_amr_feature_count": amr_count,
        "mobile_context_flag": int(mobile_count > 0),
        "plasmid_context_flag": int(plasmid_count > 0),
        "integron_context_flag": int(integron_count > 0),
        "sccmec_context_flag": int(sccmec_count > 0),
        "resistance_island_context_flag": resistance_island,
        "operon_proxy_flag": operon_proxy,
        "context_class": ";".join(labels),
        "nearby_gene_products": " || ".join(nearby_products)
    }


def main():
    print("Requirement 3: mobile-context mapping started")
    print(f"Project root: {PROJECT_ROOT}")

    if not HIT_TABLE.exists():
        raise FileNotFoundError(f"Hit table not found: {HIT_TABLE}")

    hits = pd.read_csv(HIT_TABLE)
    print(f"Loaded hit table: {hits.shape}")

    gff = load_gff_folder()
    print(f"Loaded GFF features: {gff.shape}")

    accession_col = pick_col(hits, ["accession", "genome_accession", "assembly_accession"], required=False)
    position_col = pick_col(hits, ["position_1based", "position", "genomic_position"], required=True)
    target_gene_col = pick_col(hits, ["target_gene", "gene", "amr_gene"], required=False)
    guide_col = pick_col(hits, ["guide_id", "selected_guide_id", "guide"], required=False)
    contig_col = pick_col(hits, ["seqid", "contig", "contig_id", "contig_description"], required=False)

    if accession_col is None:
        for possible in ["genome_file", "config_file", "source_file"]:
            if possible in hits.columns:
                hits["accession_inferred"] = hits[possible].astype(str).str.extract(r"(GC[AF]_\d+\.\d+)")
                accession_col = "accession_inferred"
                break

    if accession_col is None:
        raise ValueError("Could not find or infer accession column.")

    if contig_col is None:
        raise ValueError("Could not find contig/seqid/contig_description column.")

    hits["seqid_for_mapping"] = hits[contig_col].apply(infer_seqid)

    if target_gene_col is None:
        hits["target_gene_inferred"] = "unknown_gene"
        target_gene_col = "target_gene_inferred"

    if guide_col is None:
        hits["guide_id_inferred"] = "unknown_guide"
        guide_col = "guide_id_inferred"

    hits[position_col] = pd.to_numeric(hits[position_col], errors="coerce")
    hits = hits.dropna(subset=[position_col]).copy()
    hits[position_col] = hits[position_col].astype(int)

    gff_groups = {
        key: group
        for key, group in gff.groupby(["accession", "seqid"], dropna=False)
    }

    mapped_rows = []

    for _, hit in hits.iterrows():
        accession = str(hit[accession_col])
        seqid = str(hit["seqid_for_mapping"])
        pos = int(hit[position_col])

        window_start = max(1, pos - WINDOW_BP)
        window_end = pos + WINDOW_BP

        window = gff_groups.get((accession, seqid), pd.DataFrame())

        if window.empty:
            window = gff[gff["seqid"].astype(str).eq(seqid)]

        if not window.empty:
            window = window[
                (window["feature_end"] >= window_start)
                & (window["feature_start"] <= window_end)
            ]

        context = classify_window(window)

        row = hit.to_dict()
        row["requirement3_window_bp"] = WINDOW_BP
        row["requirement3_window_start"] = window_start
        row["requirement3_window_end"] = window_end
        row["requirement3_accession_used"] = accession
        row["requirement3_seqid_used"] = seqid

        row.update(context)
        mapped_rows.append(row)

    mapped = pd.DataFrame(mapped_rows)

    hit_out = OUT_DIR / "requirement3_hit_level_mobile_context_mapping.csv"
    mapped.to_csv(hit_out, index=False)

    guide_summary = (
        mapped.groupby([target_gene_col, guide_col], dropna=False)
        .agg(
            total_hits=(position_col, "count"),
            mobile_context_hits=("mobile_context_flag", "sum"),
            plasmid_context_hits=("plasmid_context_flag", "sum"),
            integron_context_hits=("integron_context_flag", "sum"),
            sccmec_context_hits=("sccmec_context_flag", "sum"),
            resistance_island_like_hits=("resistance_island_context_flag", "sum"),
            operon_proxy_hits=("operon_proxy_flag", "sum"),
        )
        .reset_index()
    )

    guide_summary["requirement3_mobile_context_burden"] = (
        guide_summary["mobile_context_hits"] * 2
        + guide_summary["plasmid_context_hits"] * 3
        + guide_summary["integron_context_hits"] * 4
        + guide_summary["sccmec_context_hits"] * 4
        + guide_summary["resistance_island_like_hits"] * 5
        + guide_summary["operon_proxy_hits"] * 1
    )

    def decision(x):
        if x == 0:
            return "Low_mobile_context_signal"
        elif x <= 10:
            return "Moderate_mobile_context_signal"
        elif x <= 50:
            return "Elevated_mobile_context_signal"
        else:
            return "High_mobile_context_signal_review_or_deprioritize"

    guide_summary["requirement3_mobile_context_decision"] = guide_summary[
        "requirement3_mobile_context_burden"
    ].apply(decision)

    guide_out = OUT_DIR / "requirement3_guide_level_mobile_context_summary.csv"
    guide_summary.to_csv(guide_out, index=False)

    context_counts = (
        mapped["context_class"]
        .value_counts()
        .reset_index()
    )
    context_counts.columns = ["context_class", "hit_count"]

    context_out = OUT_DIR / "requirement3_context_class_counts.csv"
    context_counts.to_csv(context_out, index=False)

    summary_out = OUT_DIR / "requirement3_mobile_context_mapping_summary.txt"

    with open(summary_out, "w") as f:
        f.write("Requirement 3 Mobile Context Mapping Summary\n")
        f.write("===========================================\n\n")
        f.write(f"Input hit table: {HIT_TABLE.relative_to(PROJECT_ROOT)}\n")
        f.write(f"GFF folder: {GFF_DIR.relative_to(PROJECT_ROOT)}\n")
        f.write(f"Window size: +/- {WINDOW_BP} bp\n")
        f.write(f"Total mapped hits: {len(mapped)}\n")
        f.write(f"Total guide rows summarized: {len(guide_summary)}\n\n")

        f.write("Hit-level context signals:\n")
        for col in [
            "mobile_context_flag",
            "plasmid_context_flag",
            "integron_context_flag",
            "sccmec_context_flag",
            "resistance_island_context_flag",
            "operon_proxy_flag",
        ]:
            f.write(f"- {col}: {int(mapped[col].sum())}\n")

        f.write("\nGuide-level mobile-context decision counts:\n")
        for k, v in guide_summary["requirement3_mobile_context_decision"].value_counts().items():
            f.write(f"- {k}: {v}\n")

        f.write("\nOutput files:\n")
        f.write(f"- {hit_out.relative_to(PROJECT_ROOT)}\n")
        f.write(f"- {guide_out.relative_to(PROJECT_ROOT)}\n")
        f.write(f"- {context_out.relative_to(PROJECT_ROOT)}\n")

    print(f"Wrote: {hit_out.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {guide_out.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {context_out.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {summary_out.relative_to(PROJECT_ROOT)}")

    print("\nTop guide-level mobile-context burden:")
    print(
        guide_summary.sort_values(
            ["requirement3_mobile_context_burden", "total_hits"],
            ascending=[False, False]
        ).head(20).to_string(index=False)
    )

    print("\nRequirement 3 Step 1 completed.")


if __name__ == "__main__":
    main()

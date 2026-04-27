import subprocess
import pandas as pd
from pathlib import Path

TARGETS = Path("results_external/crisot/crisot_master_targets_23nt.txt")
GENOME = "mdr_reference_combined.fa"
OUTDIR = Path("results_external/crisot/master_batch")
OUTDIR.mkdir(parents=True, exist_ok=True)

summary_rows = []

targets = [x.strip().upper() for x in TARGETS.read_text().splitlines() if x.strip()]

for i, target in enumerate(targets, start=1):
    out_file = OUTDIR / f"crisot_{i:03d}.csv"

    cmd = [
        "python",
        "external_tools/CRISOT/CRISOT.py",
        "opti",
        "--tar", target,
        "--genome", GENOME,
        "--out", str(out_file),
    ]

    print(f"[{i}/{len(targets)}] Scoring {target}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        summary_rows.append({
            "Target sequence": target,
            "error": result.stderr,
        })
        continue

    if out_file.exists():
        try:
            df = pd.read_csv(out_file)
            if len(df) > 0:
                row = df.iloc[0].to_dict()
                row["Target sequence"] = target
                summary_rows.append(row)
            else:
                summary_rows.append({"Target sequence": target, "error": "empty output"})
        except Exception as e:
            summary_rows.append({"Target sequence": target, "error": str(e)})
    else:
        summary_rows.append({"Target sequence": target, "error": "output file not created"})

summary = pd.DataFrame(summary_rows)
summary.to_csv("results_external/crisot/crisot_master_rescore_output.csv", index=False)

print("Saved: results_external/crisot/crisot_master_rescore_output.csv")
print(summary.head())

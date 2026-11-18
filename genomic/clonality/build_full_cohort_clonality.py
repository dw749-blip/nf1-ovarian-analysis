"""
NF1 LOH/deletion clonality across the FULL TCGA-OV cohort (all NF1-altered
samples in NF1_expression_functional_groups_redo.tsv), split by event type,
for the clonality summary panel.

Reads (read-only):
  ../NF1_expression_functional_groups_redo.tsv
  chr17_segments.tsv   (hg19 ABSOLUTE chr17 segments, already downloaded)

Outputs (this folder):
  all_cohort_sample_list.tsv
  all_cohort_segments_raw.tsv
  all_cohort_clonality_results.tsv
"""
import pandas as pd

NF1_HG19_START = 29421945
NF1_HG19_END = 29704695
PARENT = "../"

func = pd.read_csv(PARENT + "NF1_expression_functional_groups_redo.tsv", sep="\t",
                    header=0, index_col=0)
func = func.rename_axis("sample_id").reset_index()

cohort = func[func["sample_id"].str.startswith("TCGA")].copy()

# Combine the two mutation+LOH subtypes into one category for the figure,
# keep a separate mutation_subtype column for reference.
cohort["mutation_subtype"] = cohort["group"].where(
    cohort["group"].isin(["Missense + LOH", "Truncating + LOH"]))
cohort["event_type"] = cohort["group"].replace({
    "Missense + LOH": "Mutation + LOH",
    "Truncating + LOH": "Mutation + LOH",
})

cohort.to_csv("all_cohort_sample_list.tsv", sep="\t", index=False)
print(f"Wrote all_cohort_sample_list.tsv: {len(cohort)} samples")
print(cohort["event_type"].value_counts())

seg = pd.read_csv("chr17_segments.tsv", sep="\t", header=0)
seg = seg[seg["Sample"].isin(cohort["sample_id"])]
overlap = seg[(seg["Start"] <= NF1_HG19_END) & (seg["End"] >= NF1_HG19_START)].copy()
overlap = overlap.sort_values(["Sample", "Start"])
overlap.to_csv("all_cohort_segments_raw.tsv", sep="\t", index=False)
print(f"\nWrote all_cohort_segments_raw.tsv: {len(overlap)} segments overlapping NF1 "
      f"for {overlap['Sample'].nunique()} of {len(cohort)} samples")

missing = sorted(set(cohort["sample_id"]) - set(overlap["Sample"]))
if missing:
    print(f"NOTE: {len(missing)} samples have no ABSOLUTE call overlapping NF1: {missing}")


def classify(row):
    hscn1, hscn2 = row["Modal_HSCN_1"], row["Modal_HSCN_2"]
    is_loh = (hscn1 == 0) or (hscn2 == 0)
    if not is_loh:
        return pd.Series({"LOH_call": False, "clonality": "no LOH", "ccf_lost_allele": None})
    flags, ccfs = [], []
    if hscn1 == 0:
        flags.append(row["Subclonal_HSCN_a1"]); ccfs.append(row["Cancer_cell_frac_a1"])
    if hscn2 == 0:
        flags.append(row["Subclonal_HSCN_a2"]); ccfs.append(row["Cancer_cell_frac_a2"])
    clonality = "subclonal" if any(f == 1 for f in flags) else "clonal"
    return pd.Series({"LOH_call": True, "clonality": clonality, "ccf_lost_allele": min(ccfs)})


calls = overlap.join(overlap.apply(classify, axis=1))
calls["overlap_bp"] = calls.apply(
    lambda r: min(r["End"], NF1_HG19_END) - max(r["Start"], NF1_HG19_START), axis=1)
calls["n_segments_this_sample"] = calls.groupby("Sample")["Sample"].transform("count")
best = calls.sort_values("overlap_bp", ascending=False).drop_duplicates("Sample", keep="first")

result = cohort.merge(best, left_on="sample_id", right_on="Sample", how="left")
result = result.drop(columns=["Sample"]).rename(
    columns={"NF1_functional_status": "functional_status", "group": "genotype_group"})

cols = ["sample_id", "functional_status", "genotype_group", "event_type",
        "mutation_subtype", "expression", "Chromosome", "Start", "End",
        "n_segments_this_sample", "Modal_HSCN_1", "Modal_HSCN_2", "Modal_Total_CN",
        "Homozygous_deletion", "LOH_call", "clonality", "ccf_lost_allele", "solution"]
result = result[cols]
result.to_csv("all_cohort_clonality_results.tsv", sep="\t", index=False)

print("\n=== Summary: full cohort ===")
print(f"Total: {len(result)}   No ABSOLUTE call at NF1: {result['LOH_call'].isna().sum()}")
print(result["clonality"].value_counts(dropna=False))
print("\nBy event_type:")
print(result.groupby(["event_type", "clonality"], dropna=False).size().unstack(fill_value=0))

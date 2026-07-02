import pandas as pd
import numpy as np
import os


def calculate_pi_in_windows(hap_data: pd.DataFrame, window_size: int = 1000, step: int = 500):

    nonhap_cols = ["contig", "gene_id", "site_pos", "site_type"]
    hap_cols = [c for c in hap_data.columns if c not in nonhap_cols]

    if len(hap_cols) == 0:
        raise ValueError("No haplotype columns found.")

    # Ensure haplotype columns are numeric
    hap_data[hap_cols] = hap_data[hap_cols].apply(pd.to_numeric, errors="coerce")

    pi_results = []
    sliding_window_results = []

    for contig, contig_data in hap_data.groupby("contig", sort=False):

        # Calculate per-site pi
        non_na_count = contig_data[hap_cols].notna().sum(axis=1)
        n_correction = non_na_count.where(non_na_count > 1, other=pd.NA) / (
            non_na_count.where(non_na_count > 1, other=pd.NA) - 1
        )

        p = contig_data[hap_cols].mean(axis=1, skipna=True)
        pi_raw = 2 * p * (1 - p)
        pi_at_S = pi_raw * n_correction

        pi_init = contig_data[nonhap_cols].copy()
        pi_init["pi_at_S"] = pi_at_S

        min_pos = int(contig_data["site_pos"].min())
        max_pos = int(contig_data["site_pos"].max())
        half = window_size // 2

        # SNP-centered windows
        for pos in contig_data["site_pos"]:

            left = max(int(pos) - half, min_pos)
            right = min(int(pos) + half, max_pos)

            mask = (pi_init["site_pos"] >= left) & (pi_init["site_pos"] <= right)

            window_pi_sum = pi_init.loc[mask, "pi_at_S"].sum(min_count=1)
            window_len = right - left + 1

            gene = pi_init.loc[pi_init["site_pos"] == pos, "gene_id"].iloc[0]

            pi_val = window_pi_sum / window_len

            pi_results.append(
                {
                    "contig": contig,
                    "midpt": pos,
                    "gene_id": gene,
                    "left_coord": left,
                    "right_coord": right,
                    "pi": pi_val,
                }
            )

        # Sliding windows
        for start in range(min_pos, max_pos + 1, step):

            left = start
            right = min(start + window_size - 1, max_pos)

            mask = (pi_init["site_pos"] >= left) & (pi_init["site_pos"] <= right)

            window_pi_sum = pi_init.loc[mask, "pi_at_S"].sum()
            window_len = right - left + 1

            pi_val = window_pi_sum / window_len

            sliding_window_results.append(
                {
                    "contig": contig,
                    "left_coord": left,
                    "right_coord": right,
                    "pi": pi_val,
                }
            )

    return pd.DataFrame(pi_results), pd.DataFrame(sliding_window_results)


def main():

    input_file = "/u/project/ngarud/Garud_lab/HMP_haplos_proc/haplotypes/Ruminococcus_bromii_62047/FP929051_haplotypes.csv"

    if not os.path.exists(input_file):
        raise FileNotFoundError("C_diff.txt not found in current directory")

    print("Reading haplotype data...")

    hap_data = pd.read_csv(input_file)

    if hap_data.empty:
        raise ValueError("Input file is empty.")

    window_size = 1000
    step = 500

    print("Calculating pi...")

    pi_windows, sliding_windows = calculate_pi_in_windows(
        hap_data, window_size=window_size, step=step
    )

    pi_windows.to_csv("R_bromii_pi_windows.tsv", sep="\t", index=False)
    sliding_windows.to_csv("R_bromii_sliding_windows.tsv", sep="\t", index=False)

    print("Finished!")


if __name__ == "__main__":
    main()

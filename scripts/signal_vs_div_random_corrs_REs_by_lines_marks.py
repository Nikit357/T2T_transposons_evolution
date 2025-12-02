import os

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm  # For progress bars

import warnings
warnings.filterwarnings('ignore')
def pearsonr_ci(x, y, alpha=0.05):
    """
    Calculates Pearson r, p-value, and 95% CI using Fisher z-transformation.
    This works on all versions of SciPy.
    """
    # Calculate r and p-value
    r, p = stats.pearsonr(x, y)
    n = len(x)

    # Fisher z-transformation
    z = np.arctanh(r)

    # Standard error of z
    se_z = 1.0 / np.sqrt(n - 3)

    # Get the z-score for the desired confidence level
    z_crit = stats.norm.ppf(1 - alpha / 2)  # This is 1.96 for 95% CI

    # Calculate CI for z
    z_low = z - z_crit * se_z
    z_high = z + z_crit * se_z

    # Back-transform z CI to r CI
    r_low = np.tanh(z_low)
    r_high = np.tanh(z_high)

    return r, p, r_low, r_high


# --- Helper Class for Colored Logging ---
class Logger:
    """A simple logger for pretty, colored console output."""

    BLUE = "\033[94m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    ENDC = "\033[0m"

    @staticmethod
    def info(message):
        print(f"{Logger.BLUE}INFO: {message}{Logger.ENDC}")

    @staticmethod
    def success(message):
        print(f"{Logger.GREEN}SUCCESS: {message}{Logger.ENDC}")

    @staticmethod
    def warning(message):
        print(f"{Logger.YELLOW}WARNING: {message}{Logger.ENDC}")

    @staticmethod
    def step(message):
        print(f"  -> {message}")


cell_lines = [
    "Caco-2",
    "SJCRH30",
    "BE2C",
    "SJSA1",
    "HAP-1",
    "HL-60",
    "MG63",
    "VCaP",
    "C4-2B",
    "RWPE2",
    "22Rv1",
    "RWPE1",
]
modifications = [
    "CTCF",
    "H3K27ac",
    "H3K4me1",
    "H3K9me3",
    "H3K36me3",
    "H3K27me3",
    "H3K4me3",
]

family_names = pd.read_csv("repeat_masker_family_name_counts.csv")["Unnamed: 0"]
individual_names = pd.read_csv("repeat_masker_individual_name_counts.csv")["Unnamed: 0"]

# Name	Family	Class
# --- Main Script ---
Logger.info("Starting batch plot generation...")
families_stats_list = []

individuals_stats_list = []
# Use tqdm for the outer loop progress bar
for cell_line in tqdm(cell_lines, desc="Processing Cell Lines"):

    # Use a nested tqdm for the inner loop
    for modification in tqdm(
        modifications, desc=f"  Mods for {cell_line}", leave=False
    ):

        filepath = f"./epigenomic_files/{cell_line}.{modification}.chm13v2.0.mapped_on_repeat_masker.bedGraph"
        # --- Data Loading ---
        Logger.step(f"Loading data from: {filepath}")
        try:
            sample_enrichment = pd.read_csv(
                filepath,
                sep="\t",
                header=None,
            )
            Logger.success("File loaded successfully.")
        except FileNotFoundError:
            Logger.warning(f"File not found, skipping: {filepath}")
            continue  # Skip to the next modification
        except Exception as e:
            Logger.warning(
                f"An unexpected error occurred while reading {filepath}: {e}"
            )
            continue

        sample_enrichment.columns = [
            "Chromosome",
            "Start",
            "End",
            "Divergence",
            "Name",
            "Family",
            "Class",
            "Signal",
        ]
        sample_enrichment["Signal"] = pd.to_numeric(
            sample_enrichment["Signal"], errors="coerce"
        )
        sample_enrichment["Divergence"] = pd.to_numeric(
            sample_enrichment["Divergence"], errors="coerce"
        )
        for family_name in tqdm(family_names, desc="Processing Family Names"):

            # Filter the data for the current family
            query_local = f'Family == "{family_name}"'
            sample_enrichment_local = sample_enrichment.query(query_local)

            # Calculate correlation
            clean_data = sample_enrichment_local.dropna(subset=["Divergence", "Signal"])
            clean_data['Signal'] = clean_data['Signal'].sample(frac=1, random_state=1).set_axis(clean_data.index)
            if len(clean_data) > 3:
                corr_coef, p_val, ci_low, ci_high = pearsonr_ci(
                    clean_data["Divergence"], clean_data["Signal"]
                )
            else:
                corr_coef, p_val, ci_low, ci_high = None, None, None, None

            families_stats_list.append(
                {
                    "cell_line": cell_line,
                    "modification": modification,
                    "family_name": family_name,
                    "r": corr_coef,
                    "p_value": p_val,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                }
            )

        for individual_name in tqdm(
            individual_names, desc="Processing Individual Names"
        ):

            # Filter the data for the current name
            query_local = f'Name == "{individual_name}"'
            sample_enrichment_local = sample_enrichment.query(query_local)

            # Calculate correlation
            clean_data = sample_enrichment_local.dropna(subset=["Divergence", "Signal"])
            clean_data['Signal'] = clean_data['Signal'].sample(frac=1, random_state=1).set_axis(clean_data.index)
            if len(clean_data) > 3:
                corr_coef, p_val, ci_low, ci_high = pearsonr_ci(
                    clean_data["Divergence"], clean_data["Signal"]
                )
            else:
                corr_coef, p_val, ci_low, ci_high = None, None, None, None

            individuals_stats_list.append(
                {
                    "cell_line": cell_line,
                    "modification": modification,
                    "individual_name": individual_name,
                    "r": corr_coef,
                    "p_value": p_val,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                }
            )


correlation_df_families = pd.DataFrame(families_stats_list)
correlation_df_individuals = pd.DataFrame(individuals_stats_list)

correlation_df_families.to_csv(
    "Signal_vs_divergence_correlation_by_lines_mods_families_detailed_random.csv", index=True
)


correlation_df_individuals.to_csv(
    "Signal_vs_divergence_correlation_by_lines_mods_individuals_detailed_random.csv",
    index=True,
)
Logger.info("All processing complete!")

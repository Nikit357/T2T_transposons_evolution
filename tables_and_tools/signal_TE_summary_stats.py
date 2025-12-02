import os

import numpy as np
import pandas as pd
from tqdm import tqdm  # For progress bars

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
Logger.info("Starting batch statistics calculation...")
means_df_list_classes = []
means_df_list_families = []
means_df_list_individuals = []

sds_df_list_classes = []
sds_df_list_families = []
sds_df_list_individuals = []

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


        # --- Means ---
        class_means = sample_enrichment.groupby("Class")["Signal"].mean().to_dict()
        class_means.update({"cell_line": cell_line, "modification": modification})
        means_df_list_classes.append(class_means)

        family_means = sample_enrichment.groupby("Family")["Signal"].mean().to_dict()
        family_means.update({"cell_line": cell_line, "modification": modification})
        means_df_list_families.append(family_means)

        individual_means = sample_enrichment.groupby("Name")["Signal"].mean().to_dict()
        individual_means.update({"cell_line": cell_line, "modification": modification})
        means_df_list_individuals.append(individual_means)

        # --- Standard Deviations ---
        class_sds = sample_enrichment.groupby("Class")["Signal"].std().to_dict()
        class_sds.update({"cell_line": cell_line, "modification": modification})
        sds_df_list_classes.append(class_sds)

        family_sds = sample_enrichment.groupby("Family")["Signal"].std().to_dict()
        family_sds.update({"cell_line": cell_line, "modification": modification})
        sds_df_list_families.append(family_sds)

        individual_sds = sample_enrichment.groupby("Name")["Signal"].std().to_dict()
        individual_sds.update({"cell_line": cell_line, "modification": modification})
        sds_df_list_individuals.append(individual_sds)

# --- CORRECTED DataFrame construction ---
# Use pd.DataFrame() to build from a list of dicts (rows)
# Use index=False so you don't save the 0, 1, 2... row numbers
Logger.info("Assembling and saving DataFrames...")

pd.DataFrame(means_df_list_classes).to_csv(
    "Average_signal_by_lines_mods_TE_classes.csv", index=False
)

pd.DataFrame(means_df_list_families).to_csv(
    "Average_signal_by_lines_mods_TE_families.csv", index=False
)

pd.DataFrame(means_df_list_individuals).to_csv(
    "Average_signal_by_lines_mods_TE_individuals.csv", index=False
)

pd.DataFrame(sds_df_list_classes).to_csv(
    "SD_signal_by_lines_mods_TE_classes.csv", index=False
)

pd.DataFrame(sds_df_list_families).to_csv(
    "SD_signal_by_lines_mods_TE_families.csv", index=False
)

pd.DataFrame(sds_df_list_individuals).to_csv(
    "SD_signal_by_lines_mods_TE_individuals.csv", index=False
)

Logger.success("All processing complete!")

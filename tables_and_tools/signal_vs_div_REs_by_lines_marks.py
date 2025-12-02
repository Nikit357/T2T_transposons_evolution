import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from tqdm import tqdm  # For progress bars
plt.rcParams['svg.fonttype'] = 'none'

# --- Helper Class for Colored Logging ---
class Logger:
    """A simple logger for pretty, colored console output."""
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    ENDC = '\033[0m'

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
    "C4-2B",
    "RWPE1",
    "VCaP",
    "SJCRH30",
    "SJSA1",
    "HAP-1",
    "BE2C",
    "HL-60",
    "MG63",
    "Caco-2",
#    "C4-2B",
    "RWPE2",
    "22Rv1",
#    "RWPE1",
]
modifications = [
    "CTCF",
    "H3K27ac",
    "H3K4me1",
    "H3K9me3",
    "H3K36me3",
    "H3K27me3",
    "H3K4me3",
#    "H3K27ac",
]
class_names = ["LINE", "LTR", "SINE", "DNA", "Retroposon", "RC"]

# --- Main Script ---
Logger.info("Starting batch plot generation...")
correlation_df_list = []
# Use tqdm for the outer loop progress bar
for cell_line in tqdm(cell_lines, desc="Processing Cell Lines"):
    
    # Use a nested tqdm for the inner loop
    for modification in tqdm(modifications, desc=f"  Mods for {cell_line}", leave=False):
        
        filepath = f"./epigenomic_files/{cell_line}.{modification}.chm13v2.0.mapped_on_repeat_masker.bedGraph"
        save_path_png = f"plots/Epigenetic_signal_vs_divergence_by_class/Divergence_vs_Signal_by_classes_{cell_line}_{modification}_tidy.png"
        save_path_svg = f"plots/Epigenetic_signal_vs_divergence_by_class/Divergence_vs_Signal_by_classes_{cell_line}_{modification}_tidy.svg"        
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
            continue # Skip to the next modification
        except Exception as e:
            Logger.warning(f"An unexpected error occurred while reading {filepath}: {e}")
            continue

        sample_enrichment.columns = [
            "Chromosome", "Start", "End", "Divergence", "Name",
            "Family", "Class", "Signal",
        ]
        sample_enrichment["Signal"] = pd.to_numeric(sample_enrichment["Signal"], errors="coerce")
        sample_enrichment["Divergence"] = pd.to_numeric(sample_enrichment["Divergence"], errors="coerce")

        # --- Main Plotting Logic ---
        Logger.step("Generating jointplot grid...")
        correlations = {'cell_line': cell_line, 'modification': modification}
        figure = plt.figure(figsize=(16, 10))
        class_names_grid = np.array(class_names).reshape((2, 3))
        
        # Create the main 2x3 grid for the entire figure
        outer_grid = GridSpec(2, 3, figure=figure, wspace=0.1, hspace=0.1)
        
        # --- Store references to the first axis in each row/column ---
        shared_x_axes = [None] * 3  # One for each column
        shared_y_axes = [None] * 2  # One for each row
        
        for r in range(2):
            for c in range(3):
                class_name = class_names_grid[r, c]
                tqdm.write(f"Plotting: {class_name} (Row {r}, Col {c})")
                
                # For each cell in the outer grid, create a nested GridSpec
                inner_grid = GridSpecFromSubplotSpec(
                    2, 2, subplot_spec=outer_grid[r, c],
                    width_ratios=[5, 1], height_ratios=[1, 5],
                    wspace=0.05, hspace=0.05,
                )
        
                # --- Create axes with sharing logic ---
                ax_main = figure.add_subplot(
                    inner_grid[1, 0], 
                    sharex=shared_x_axes[c], 
                    sharey=shared_y_axes[r]
                )
                ax_hist_x = figure.add_subplot(inner_grid[0, 0], sharex=ax_main)
                ax_hist_y = figure.add_subplot(inner_grid[1, 1], sharey=ax_main)
                
                # --- Store the first axis of each row/col as reference ---
                if r == 0:
                    shared_x_axes[c] = ax_main
                if c == 0:
                    shared_y_axes[r] = ax_main
        
                # Filter the data for the current class
                query_local = f'Class == "{class_name}"'
                sample_enrichment_local = sample_enrichment.query(query_local)
        
                # --- Plot the data ---
                sns.histplot(
                    data=sample_enrichment_local, x="Divergence", y="Signal",
                    bins=48, pmax=0.9, cmap="viridis", ax=ax_main,
                )
                sns.histplot(data=sample_enrichment_local, x="Divergence", ax=ax_hist_x, bins=50)
                sns.histplot(data=sample_enrichment_local, y="Signal", ax=ax_hist_y, bins=50)
        
                # --- START: MODIFIED TIDY-UP BLOCK ---
                # Tidy up: Remove all ticks, labels, and borders from marginal plots
                
                # Top histogram (ax_hist_x)
                ax_hist_x.set_ylabel("")
                ax_hist_x.set_xlabel("")
                ax_hist_x.tick_params(axis="both", bottom=False, left=False, 
                                      labelbottom=False, labelleft=False)
                for spine in ax_hist_x.spines.values():
                    spine.set_visible(False)
                    
                # Right histogram (ax_hist_y)
                ax_hist_y.set_xlabel("")
                ax_hist_y.set_ylabel("")
                ax_hist_y.tick_params(axis="both", bottom=False, left=False, 
                                      labelbottom=False, labelleft=False)
                for spine in ax_hist_y.spines.values():
                    spine.set_visible(False)
                # --- END: MODIFIED TIDY-UP BLOCK ---
                
                # --- Tidy up: Main plot ticks/labels (for shared axes) ---
                if c > 0: # Not in the first column
                    ax_main.tick_params(axis="y", labelleft=False)
                    ax_main.set_ylabel("")
                else:
                    ax_main.set_ylabel("Signal", fontsize=12) # Only show Y label on first col
        
                if r < 1: # Not in the bottom row
                    ax_main.tick_params(axis="x", labelbottom=False)
                    ax_main.set_xlabel("")
                else:
                    ax_main.set_xlabel("Divergence", fontsize=12) # Only show X label on bottom row
        
                # Calculate correlation
                clean_data = sample_enrichment_local.dropna(subset=["Divergence", "Signal"])
                if len(clean_data) > 1:
                    corr_coef = np.corrcoef(clean_data["Divergence"], clean_data["Signal"])[0, 1]
                    corr_str = f"r={corr_coef:.3f}"
                else:
                    corr_str = "r=N/A"
        
                # Set the title for this "jointplot" cell
                ax_hist_x.set_title(f"{class_name}, {corr_str}", fontsize=16, pad=10)
                correlations.update({class_name: corr_coef})

        correlation_df_list.append(pd.Series(correlations))

        # Final figure adjustments
        figure.suptitle(f"Divergence vs. Signal for {cell_line} - {modification}", fontsize=20, y=0.98)
        figure.tight_layout(rect=[0, 0.03, 1, 0.95])

        # Save the final figure
        Logger.step(f"Saving plot to: {save_path_png}")
        # Ensure the directory exists
        os.makedirs(os.path.dirname(save_path_png), exist_ok=True)
        plt.savefig(save_path_png)
        plt.savefig(save_path_svg)
        plt.close(figure) # IMPORTANT: Close the figure to free up memory
        Logger.success("Plot saved.")
#correlation_df = pd.concat(correlation_df_list)
#correlation_df.to_csv('Signal_vs_divergence_correlation_by_lines_mods_classes.csv')
Logger.info("All processing complete!")

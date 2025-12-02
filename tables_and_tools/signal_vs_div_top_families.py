import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from tqdm import tqdm  # For progress bars

correlation_df_families = pd.read_csv(
    "Signal_vs_divergence_correlation_by_lines_mods_families.csv", index_col=0
).T
family_names = list(correlation_df_families.columns[2:])
correlation_df_families[family_names] = correlation_df_families[family_names].apply(
    pd.to_numeric, errors="coerce"
)
corrs_for_test = correlation_df_families.melt(
    id_vars=["cell_line", "modification"], value_vars=family_names
)
corrs_for_test.columns = ["cell_line", "modification", "TE_family", "correlation"]
corrs_for_test["transformed_corr"] = np.arctanh(corrs_for_test.correlation)

corrs_for_test = corrs_for_test.sort_values('correlation')

top_corrs = corrs_for_test[-10:].reset_index(drop=True)


plt.rcParams["svg.fonttype"] = "none"


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
    def s_success(message):
        return f"{Logger.GREEN}SUCCESS: {message}{Logger.ENDC}"

    @staticmethod
    def s_warning(message):
        return f"{Logger.YELLOW}WARNING: {message}{Logger.ENDC}"


save_path_svg = f"plots/Divergence_vs_Signal_by_families_top10_correlations.svg"

figure = plt.figure(figsize=(30, 10))

N_ROWS = 2
N_COLS = 5

# Create the main 2x5 grid for the entire figure
# MODIFIED: Increased wspace and hspace for independent labels
outer_grid = GridSpec(N_ROWS, N_COLS, figure=figure, wspace=0.4, hspace=0.5)

# --- REMOVED: Shared axis lists are no longer needed ---
# shared_x_axes = [None] * N_COLS
# shared_y_axes = [None] * N_ROWS

# --- IMPROVEMENT: Add a cache for data loading ---
data_cache = {}

# --- IMPROVEMENT: Use a more Pythonic loop ---
for i, (idx, row) in enumerate(tqdm(top_corrs.iterrows(), total=len(top_corrs), desc="Plotting Top 10")):
    
    # Calculate grid position
    r = i // N_COLS
    c = i % N_COLS

    # Get data from the row
    family_name = row['TE_family']
    cell_line = row['cell_line']
    modification = row['modification']

    filepath = f"./epigenomic_files/{cell_line}.{modification}.chm13v2.0.mapped_on_repeat_masker.bedGraph"

    # --- Check cache before loading ---
    if filepath in data_cache:
        tqdm.write(f"  -> Found in cache: {filepath}")
        sample_enrichment = data_cache[filepath]
    else:
        # --- Data Loading ---
        tqdm.write(f"  -> Loading data from: {filepath}")
        try:
            sample_enrichment = pd.read_csv(
                filepath,
                sep="\t",
                header=None,
            )
            data_cache[filepath] = sample_enrichment
            tqdm.write(Logger.s_success("File loaded and cached."))
            
        except FileNotFoundError:
            tqdm.write(Logger.s_warning(f"File not found, skipping: {filepath}"))
            ax = figure.add_subplot(outer_grid[r, c])
            ax.axis('off')
            continue
        except Exception as e:
            tqdm.write(Logger.s_warning(
                f"An unexpected error occurred while reading {filepath}: {e}"
            ))
            ax = figure.add_subplot(outer_grid[r, c])
            ax.axis('off')
            continue

        # --- Process data (only needs to be done once, on load) ---
        sample_enrichment.columns = [
            "Chromosome", "Start", "End", "Divergence",
            "Name", "Family", "Class", "Signal",
        ]
        sample_enrichment["Signal"] = pd.to_numeric(
            sample_enrichment["Signal"], errors="coerce"
        )
        sample_enrichment["Divergence"] = pd.to_numeric(
            sample_enrichment["Divergence"], errors="coerce"
        )
        data_cache[filepath] = sample_enrichment

    # --- Query must be done every time ---
    query_local = f'Family == "{family_name}"'
    sample_enrichment_local = sample_enrichment.query(query_local)

    # For each cell in the outer grid, create a nested GridSpec
    inner_grid = GridSpecFromSubplotSpec(
        2, 2, subplot_spec=outer_grid[r, c],
        width_ratios=[5, 1], height_ratios=[1, 5],
        wspace=0.05, hspace=0.05,
    )

    # --- Create axes with sharing logic ---
    # MODIFIED: ax_main no longer shares with other main plots
    ax_main = figure.add_subplot(inner_grid[1, 0])
    
    # --- NEW: Set background color for main plot ---
    ax_main.set_facecolor('black')

    # Marginal plots *still* share with their own main plot
    ax_hist_x = figure.add_subplot(inner_grid[0, 0], sharex=ax_main)
    ax_hist_y = figure.add_subplot(inner_grid[1, 1], sharey=ax_main)

    # --- REMOVED: Storing shared axes references ---

    # --- Plot the data ---
    sns.histplot(
        data=sample_enrichment_local, x="Divergence", y="Signal",
        bins=48, pmax=0.9, cmap="viridis", ax=ax_main,
    )
    sns.histplot(
        data=sample_enrichment_local, x="Divergence", ax=ax_hist_x, bins=50
    )
    sns.histplot(
        data=sample_enrichment_local, y="Signal", ax=ax_hist_y, bins=50
    )

    # --- Tidy up: Marginal plots (unchanged) ---
    ax_hist_x.set_ylabel("")
    ax_hist_x.set_xlabel("")
    ax_hist_x.tick_params(axis="both", bottom=False, left=False,
                          labelbottom=False, labelleft=False)
    for spine in ax_hist_x.spines.values():
        spine.set_visible(False)

    ax_hist_y.set_xlabel("")
    ax_hist_y.set_ylabel("")
    ax_hist_y.tick_params(axis="both", bottom=False, left=False,
                          labelbottom=False, labelleft=False)
    for spine in ax_hist_y.spines.values():
        spine.set_visible(False)

    # --- Tidy up: Main plot ticks/labels ---
    # MODIFIED: All plots get labels and white text
    ax_main.set_ylabel(
        f"{cell_line}\n{modification}", fontsize=10, color='white'
    )
    ax_main.set_xlabel("Divergence", fontsize=12, color='white')
    ax_main.tick_params(axis='x', colors='white')
    ax_main.tick_params(axis='y', colors='white')
    for spine in ax_main.spines.values():
        spine.set_color('white')

    # Calculate correlation
    clean_data = sample_enrichment_local.dropna(
        subset=["Divergence", "Signal"]
    )
    if len(clean_data) > 1:
        corr_coef = np.corrcoef(
            clean_data["Divergence"], clean_data["Signal"]
        )[0, 1]
        corr_str = f"r={corr_coef:.3f}"
    else:
        corr_coef = np.nan # Set to nan
        corr_str = "r=N/A"

    # Set the title
    ax_hist_x.set_title(f"{family_name}\n{corr_str}", fontsize=14, pad=10)


# --- Ensure plots directory exists ---
os.makedirs(os.path.dirname(save_path_svg), exist_ok=True)
figure.savefig(save_path_svg, format="svg", facecolor='white', transparent=False)
plt.show()

Logger.success(f"Plot saved to {save_path_svg}")

import os
import re

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy import stats
from supervenn import supervenn
from tqdm import tqdm


def calculate_weighted_averages(df):
    """
    Calculates weighted averages for numeric columns using hit length as weight.
    Length is calculated as abs(rep_end - rep_start).
    """

    # Columns expected to contain numeric data or "."
    numeric_cols = [
        #"sw_score",
        "perc_div",
        #"perc_del",
        #"perc_ins",
        "ts_tv_ratio",
        #"transitions",
        #"transversions",
        #"gap_init_rate",
        #"avg_gap_size",
        "kimura_divCpGMod",
        #"cpg_sites",
        "kimura_unadjusted",
    ]

    def process_weighted_avg(row):
        # 1. Calculate Weights (Length of each hit)
        try:
            starts = [
                float(x) if x != "." else np.nan
                for x in str(row["rep_start"]).split(",")
            ]
            ends = [
                float(x) if x != "." else np.nan for x in str(row["rep_end"]).split(",")
            ]

            # Weight = |end - start|
            weights = np.abs(np.array(ends) - np.array(starts))
        except (ValueError, TypeError):
            # If coordinates are missing or malformed, return row as is
            return row

        # 2. Iterate through each numeric column to calculate the average
        for col in numeric_cols:
            if col not in row or row[col] is None or str(row[col]) == "None":
                row[f"{col}_weighted_average"] = np.nan
                continue

            try:
                # Convert comma string to numeric array, treating "." as NaN
                vals = np.array(
                    [float(x) if x != "." else np.nan for x in str(row[col]).split(",")]
                )

                # Mask out NaNs from both values and weights
                mask = ~np.isnan(vals) & ~np.isnan(weights)

                valid_vals = vals[mask]
                valid_weights = weights[mask]

                if len(valid_vals) > 0 and np.sum(valid_weights) > 0:
                    weighted_avg = np.average(valid_vals, weights=valid_weights)
                    row[f"{col}_weighted_average"] = weighted_avg
                else:
                    row[f"{col}_weighted_average"] = np.nan
            except (ValueError, ZeroDivisionError):
                row[f"{col}_weighted_average"] = np.nan

        return row

    # Apply calculations row-wise
    return df.apply(process_weighted_avg, axis=1)


repeat_masker_with_div_filtered = pd.read_csv("repeat_masker_with_div_filtered.csv")
repeat_masker_with_div_filtered_average = calculate_weighted_averages(
    repeat_masker_with_div_filtered
)
repeat_masker_with_div_filtered_average.to_csv(
    "repeat_masker_with_div_filtered_average.csv"
)

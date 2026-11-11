#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 16:28:04 2025

@author: B. Mary (ICA-CSIC)
"""

import pandas as pd
from datetime import datetime

data = {
    "Plot": [1, 2, 5, 6, 9, 12, 13, 3, 4, 7, 8, 10, 11],
    "Treatment": [
        "Control", "Control", "Control", "Control", "Control", "Control", "Control",
        "Fascines", "Fascines",
        "Mulching", "Mulching", "Mulching", "Mulching"
    ],
    "Initial_LAI": [0.109, 0.115, 0.079, 0.095, 0.103, 0.094, 0.082,
                    0.087, 0.078, 0.138, 0.180, 0.097, 0.074],
    "Final_LAI": [0.168, 0.177, 0.205, 0.190, 0.191, 0.193, 0.235,
                  0.193, 0.185, 0.261, 0.266, 0.238, 0.220],
    "Recovery_rate_LAI_per_year_MB": [0.044, 0.055, 0.067, 0.066, 0.061, 0.061, 0.069,
                                   0.066, 0.068, 0.079, 0.074, 0.073, 0.073],
    # "Category_recovery_rate": ["Low", "Low", "Medium", "Medium", "Low", "Low", "Medium",
    #                            "Medium", "Medium", "High", "High", "High", "High"],
    "Recovery_percent_MB": [53.9, 54.2, 160.6, 99.8, 86.1, 105.8, 186.3,
                         120.9, 137.7, 88.7, 47.5, 145.8, 198.0],
    # "Category_recovery_percent": ["Low", "Low", "High", "Medium", "Low", "Medium", "High",
    #                               "Medium", "Medium", "Medium", "Low", "High", "High"],
    "Variability_STD": [0.069, 0.078, 0.093, 0.088, 0.083, 0.083, 0.100,
                        0.090, 0.096, 0.112, 0.112, 0.094, 0.101]
}

df = pd.DataFrame(data)
print(df)

#%%

# Time span in years (Aug 2020 → Dec 2023)
# Define dates
start_date = datetime(2020, 8, 1)    # August 1, 2020
end_date = datetime(2023, 12, 1)     # December 1, 2023

# Calculate difference in days
delta_days = (end_date - start_date).days

# Convert to years (approximate, using 365.25 days per year to account for leap years)
time_years = delta_days / 365.25

print("Time span in years:", time_years)

# Recalculate Recovery Rate (LAI/year)
df["Recovery_rate_LAI_per_year"] = (df["Final_LAI"] - df["Initial_LAI"]) / time_years

# Recalculate Recovery Percent (%)
df["Recovery_percent"] = (df["Final_LAI"] - df["Initial_LAI"]) / df["Initial_LAI"] * 100

# Categorize Recovery Rate (example thresholds)
def categorize_rate(rate):
    if rate < 0.06:
        return "Low"
    elif rate < 0.075:
        return "Medium"
    else:
        return "High"

df["Category_recovery_rate"] = df["Recovery_rate_LAI_per_year"].apply(categorize_rate)

# Categorize Recovery Percent (example thresholds)
def categorize_percent(percent):
    if percent < 80:
        return "Low"
    elif percent < 150:
        return "Medium"
    else:
        return "High"

df["Category_recovery_percent"] = df["Recovery_percent"].apply(categorize_percent)

# Show updated DataFrame
print(df)

#%%

df["Recovery_rate_LAI_per_year"] - df["Recovery_rate_LAI_per_year_MB"] 
df["Recovery_percent"] - df["Recovery_percent_MB"] 

#%%

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(8,6))

# Horizontal bar plot with standard deviation
sns.barplot(
    x="Recovery_rate_LAI_per_year",
    y="Treatment",
    data=df,
    ci="sd",          # show standard deviation as error bars
    palette="Set2",
    orient="h"
)

# Overlay individual plot points
sns.stripplot(
    x="Recovery_rate_LAI_per_year",
    y="Treatment",
    data=df,
    color="black",
    size=6,
    jitter=True,
    orient="h"
)

plt.xlabel("Recovery Rate (LAI/year)")
plt.ylabel("Treatment")
plt.title("Mean Recovery Rate per Treatment (with Variability)")
plt.tight_layout()
plt.show()


#%%

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# Assume df already has 'Recovery_rate_LAI_per_year' and 'Treatment'

# 1️⃣ Calculate Tukey HSD
tukey = pairwise_tukeyhsd(endog=df["Recovery_rate_LAI_per_year"], 
                          groups=df["Treatment"], alpha=0.05)
print(tukey)  # Optional: see results

# 2️⃣ Map significance letters manually based on Tukey results
# Example: Control (a), Fascines (ab), Mulching (b)
sig_letters = {
    "Control": "a",
    "Fascines": "ab",
    "Mulching": "b"
}

# 3️⃣ Plot mean recovery rate per treatment with SD and significance letters
plt.figure(figsize=(8,6))
sns.barplot(
    x="Recovery_rate_LAI_per_year",
    y="Treatment",
    data=df,
    ci="sd",
    palette="Set2",
    orient="h"
)

# Overlay individual plot points
sns.stripplot(
    x="Recovery_rate_LAI_per_year",
    y="Treatment",
    data=df,
    color="black",
    size=6,
    jitter=True,
    orient="h"
)

# Add significance letters
for i, treatment in enumerate(df["Treatment"].unique()):
    mean_val = df[df["Treatment"]==treatment]["Recovery_rate_LAI_per_year"].mean()
    plt.text(mean_val + 0.002, i, sig_letters[treatment], fontsize=14, fontweight='bold')

plt.xlabel("Recovery Rate (LAI/year)")
plt.ylabel("Treatment")
plt.title("Mean Recovery Rate per Treatment with Significance")
plt.tight_layout()
plt.show()


#%%

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Assuming 'df' is your DataFrame
# Set Seaborn style
sns.set(style="whitegrid")

# 1️⃣ LAI Change by Plot (Initial vs Final)
plt.figure(figsize=(10,6))
df_melted = df.melt(id_vars=["Plot","Treatment"], value_vars=["Initial_LAI","Final_LAI"],
                    var_name="LAI_Type", value_name="LAI")
sns.barplot(x="Plot", y="LAI", hue="LAI_Type", data=df_melted, palette="Set2")
plt.title("Initial vs Final LAI per Plot")
plt.ylabel("LAI")
plt.xlabel("Plot")
plt.legend(title="")
plt.show()

# 2️⃣ Recovery Rate by Treatment (Boxplot)
plt.figure(figsize=(8,5))
sns.boxplot(x="Treatment", y="Recovery_rate_LAI_per_year", data=df, palette="Set3")
sns.stripplot(x="Treatment", y="Recovery_rate_LAI_per_year", data=df, color="black", size=6, jitter=True)
plt.title("Recovery Rate (LAI/year) by Treatment")
plt.ylabel("Recovery rate (LAI/year)")
plt.show()

# 3️⃣ Recovery Percentage vs Initial LAI (Scatter)
plt.figure(figsize=(8,6))
sns.scatterplot(x="Initial_LAI", y="Recovery_percent", hue="Treatment",
                size="Variability_STD", sizes=(50,200), data=df, palette="Set1", alpha=0.8)
plt.title("Recovery % vs Initial LAI (Bubble size = Variability)")
plt.xlabel("Initial LAI")
plt.ylabel("Recovery (%)")
plt.show()



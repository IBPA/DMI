import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl
from sklearn.decomposition import PCA

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
})

os.makedirs("results/PCA", exist_ok=True)

concentration = pd.read_csv("compound_conc.csv")
product = pd.read_csv("products.csv")

concentration = concentration.rename(columns={"Unnamed: 0": "SampleID"})
product = product.rename(columns={"DMD ID": "SampleID"})

merged = concentration.merge(
    product[["SampleID", "Region Code", "Brand Code", "Time Point Code"]],
    on="SampleID"
)

feature_cols = [c for c in concentration.columns if c != "SampleID"]
X = merged[feature_cols].fillna(0)
X_centered = X - X.mean(axis=0)

pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_centered)

merged["PC1"] = X_pca[:, 0]
merged["PC2"] = X_pca[:, 1]

region_code_to_name = {
    1: "Northeast",
    2: "Mid-South",
    3: "West",
    4: "South Central",
    5: "Plains",
    6: "California-1",
    7: "Southeast",
    8: "Great Lakes",
    9: "California-2",
    10: "California-3",
}

region_to_color = {
    1: "#8DD3C7",
    2: "#FFFFB3",
    3: "#BEBADA",
    4: "#FB8072",
    5: "#80B1D3",
    6: "#FDB462",
    7: "#B3DE69",
    8: "#FCCDE5",
    9: "#D9D9D9",
    10: "#BC80BD",
}

brand_code_to_label = {1: "Brand A", 2: "Brand B"}
brand_code_to_marker = {1: "o", 2: "^"}

time_code_to_label = {1: "T1", 2: "T2", 3: "T3"}
time_values = sorted(merged["Time Point Code"].dropna().unique())
time_to_alpha_nominal = {1: 0.5, 2: 0.75, 3: 1.0}
time_to_alpha = {t: time_to_alpha_nominal.get(t, 0.75) for t in time_values}

fig = plt.figure(figsize=(7, 5))
ax = plt.gca()

ax.grid(False)

fig.patch.set_facecolor("white")
ax.set_facecolor("white")

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

for spine in ["left", "bottom"]:
    ax.spines[spine].set_linewidth(0.7)
    ax.spines[spine].set_color("#a0a0a0")

ax.tick_params(axis="both", length=3, width=0.7, color="#a0a0a0")

for _, row in merged.iterrows():
    region_code = row["Region Code"]
    brand_code = row["Brand Code"]
    time_code = row["Time Point Code"]
    alpha = time_to_alpha.get(time_code, 0.75)
    lw = 0.4 + 0.2 * alpha
    ax.scatter(
        row["PC1"],
        row["PC2"],
        s=32,
        edgecolors="black",
        linewidths=lw,
        color=region_to_color[region_code],
        marker=brand_code_to_marker.get(brand_code, "o"),
        alpha=alpha,
    )

regions = sorted(merged["Region Code"].dropna().unique())
region_handles = [
    mlines.Line2D(
        [], [],
        marker="o",
        linestyle="",
        markerfacecolor=region_to_color[r],
        markeredgecolor="black",
        markeredgewidth=0.8,
        markersize=5,
        label=region_code_to_name.get(r, str(r)),
    )
    for r in regions
]

brands = sorted(merged["Brand Code"].dropna().unique())
brand_handles = [
    mlines.Line2D(
        [], [],
        marker=brand_code_to_marker.get(b, "o"),
        linestyle="",
        markerfacecolor="white",
        markeredgecolor="black",
        markeredgewidth=0.8,
        markersize=5,
        label=brand_code_to_label.get(b, f"Brand {b}"),
    )
    for b in brands
]

time_handles = [
    mlines.Line2D(
        [], [],
        marker="o",
        linestyle="",
        markerfacecolor="black",
        markeredgecolor="black",
        markeredgewidth=0.4 + 0.2 * time_to_alpha[t],
        markersize=5,
        alpha=time_to_alpha[t],
        label=time_code_to_label.get(t, f"T{t}"),
    )
    for t in time_values
]

leg_brand = ax.legend(
    handles=brand_handles,
    title="Brand",
    bbox_to_anchor=(1.02, 1.0),
    loc="upper left",
    borderaxespad=0.0,
    frameon=False,
)
leg_brand._legend_box.align = "left"
ax.add_artist(leg_brand)

leg_time = ax.legend(
    handles=time_handles,
    title="Time point",
    bbox_to_anchor=(1.02, 0.85),
    loc="upper left",
    borderaxespad=0.0,
    frameon=False,
)
leg_time._legend_box.align = "left"
ax.add_artist(leg_time)

leg_region = ax.legend(
    handles=region_handles,
    title="Region",
    bbox_to_anchor=(1.02, 0.65),
    loc="upper left",
    borderaxespad=0.0,
    frameon=False,
)
leg_region._legend_box.align = "left"

ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}% variance)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}% variance)")
ax.set_title("PCA of 60 Milk Samples")

ax.ticklabel_format(style="plain", axis="both")
ax.xaxis.get_offset_text().set_visible(False)
ax.yaxis.get_offset_text().set_visible(False)

plt.tight_layout(rect=[0, 0, 0.78, 1])
plt.savefig("pca_60_milk_samples.svg", format="svg")
plt.close()
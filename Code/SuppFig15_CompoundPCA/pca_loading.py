import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
molecules = pd.read_csv("molecules.csv", encoding="latin1")

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

loadings = pca.components_.T
pc1_load = loadings[:, 0]
pc2_load = loadings[:, 1]
feature_names = np.array(feature_cols)

mol_map = (
    molecules
    .set_index("DMD ID")["Molecule Name"]
    .to_dict()
)

mapped_names = np.array([mol_map.get(f, f) for f in feature_names])

idx_pc1 = np.argsort(np.abs(pc1_load))[::-1][:10]
idx_pc2 = np.argsort(np.abs(pc2_load))[::-1][:10]

pc1_vals = pc1_load[idx_pc1]
pc2_vals = pc2_load[idx_pc2]
names_pc1 = mapped_names[idx_pc1]
names_pc2 = mapped_names[idx_pc2]

order_pc1 = np.argsort(pc1_vals)
order_pc2 = np.argsort(pc2_vals)

pc1_vals_sorted = pc1_vals[order_pc1]
pc2_vals_sorted = pc2_vals[order_pc2]
names_pc1_sorted = names_pc1[order_pc1]
names_pc2_sorted = names_pc2[order_pc2]

fig, axes = plt.subplots(1, 2, figsize=(9, 6))

def signed_colors(vals):
    return ["#c43c3c" if v > 0 else "#1aa6b7" for v in vals]

ax1 = axes[0]
y_pos1 = np.arange(len(names_pc1_sorted))
ax1.barh(y_pos1, pc1_vals_sorted, color=signed_colors(pc1_vals_sorted))
ax1.set_yticks(y_pos1)
ax1.set_yticklabels(names_pc1_sorted)
ax1.invert_yaxis()
ax1.axvline(0, color="#555555", linewidth=0.8)
ax1.set_xlabel("Loading (PC 1)")
ax1.set_ylabel("Compound")
ax1.set_title("Loadings of PC 1")

ax1.text(-0.1, 1.02, "A", transform=ax1.transAxes,
         fontsize=11, fontweight="bold", ha="right", va="bottom")

ax2 = axes[1]
y_pos2 = np.arange(len(names_pc2_sorted))
ax2.barh(y_pos2, pc2_vals_sorted, color=signed_colors(pc2_vals_sorted))
ax2.set_yticks(y_pos2)
ax2.set_yticklabels(names_pc2_sorted)
ax2.invert_yaxis()
ax2.axvline(0, color="#555555", linewidth=0.8)
ax2.set_xlabel("Loading (PC 2)")
ax2.set_title("Loadings of PC 2")

ax2.text(-0.1, 1.02, "B", transform=ax2.transAxes,
         fontsize=11, fontweight="bold", ha="right", va="bottom")

for ax in axes:
    ax.set_facecolor("white")
    for spine in ["top", "right"]:
        spine_obj = ax.spines[spine]
        spine_obj.set_visible(False)
    for spine in ["left", "bottom"]:
        spine_obj = ax.spines[spine]
        spine_obj.set_linewidth(0.7)
        spine_obj.set_color("#a0a0a0")
    ax.tick_params(axis="both", length=3, width=0.7, color="#a0a0a0")

plt.tight_layout()
plt.savefig("pca_top10_loadings_PC1_PC2.svg", format="svg")
plt.close()
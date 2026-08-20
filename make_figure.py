#!/usr/bin/env python3
"""
Generate a comparison figure for the essential-gene review manuscript.

Figure: four method paradigms scored across six qualitative dimensions
(with values 0-5). Produces a compact radar + horizontal strength chart.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Categories and dimensions
cats = ["Expression-robust\n(scEssentials)",
        "Functional\n(DepMap/CRISPR)",
        "Graph/topology\n(EPGAT, FluxGAT)",
        "Generative / scCRISPR\n(scLAMBDA, PerturbNet)"]
dims = ["Data\navailability",
        "Interpret-\nability",
        "Genomic\ncoverage",
        "Cross-\nspecies",
        "Personal/\ncontext\ndetail",
        "Single-cell\nresolution"]

# Scores 0-5
data = np.array([
    [4, 5, 2, 3, 2, 5],   # Expression-robust
    [1, 5, 5, 2, 3, 2],   # Functional
    [3, 3, 4, 4, 3, 3],   # Graph/topology
    [2, 3, 3, 3, 4, 5],   # Generative/scCRISPR
])

N = len(dims)
angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
angles += angles[:1]

# --- Panel A: radar ---
fig = plt.figure(figsize=(13, 5.6))
ax = fig.add_subplot(121, polar=True)
colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
for i, d in enumerate(data):
    vals = d.tolist() + d[:1].tolist()
    ax.plot(angles, vals, color=colors[i], linewidth=2, label=cats[i].replace("\n"," ")[:28])
    ax.fill(angles, vals, color=colors[i], alpha=0.15)
ax.set_xticks(angles[:-1])
ax.set_xticklabels(dims, fontsize=6.5)
ax.set_ylim(0, 5)
ax.set_yticks([1,2,3,4,5])
ax.set_yticklabels(["1","2","3","4","5"], fontsize=6)
ax.set_ylim(0,5)
ax.set_title("(A) Method paradigm comparison", fontsize=10, pad=18)
ax.grid(True)

# --- Panel B: grouped bar ---
ax2 = fig.add_subplot(122)
x = np.arange(len(dims))
width = 0.2
for i, d in enumerate(data):
    ax2.bar(x + i*width, d, width, label=cats[i].replace("\n"," ")[:28], color=colors[i], alpha=0.85)
ax2.set_xticks(x + width*1.5)
ax2.set_xticklabels([d.replace("\n","\n") for d in dims], fontsize=7)
ax2.set_ylabel("Qualitative score (0-5)", fontsize=8)
ax2.set_ylim(0,5.3)
ax2.set_title("(B) Dimension scores", fontsize=10)
ax2.legend(fontsize=6, loc="upper right", framealpha=0.9)

plt.tight_layout()
out = "/data/liuwenting/tools/essentialgene/manuscript_figure_comparison.png"
plt.savefig(out, dpi=200, bbox_inches="tight")
print("saved", out)

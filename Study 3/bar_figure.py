
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

mpl.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})

# network colours matched to Panel B brain figure 
NET = {
    "Default Mode": "#D34F4F",   # red
    "Control":      "#E8943A",   # orange
    "Limbic":       "#EFC11F",   # yellow
    "Sensorimotor": "#5B9BD5",   # blue
}

# (region, network, t, p_fdr, MNI [x,y,z])
data = [
    ("Superior lateral occipital cortex (R)", "Default Mode", -4.7, 0.020, [51, -59, 44]),
    ("Superior lateral occipital cortex (L)", "Default Mode", -3.6, 0.025, [-39, -80, 31]),
    ("Angular gyrus (R)",                     "Default Mode", -3.5, 0.025, [54, -50, 28]),
    ("Temporal pole (R)",                     "Default Mode", -3.5, 0.027, [47, 13, -30]),
    ("Anterior fusiform cortex (L)",          "Limbic",       -3.4, 0.034, [-29, -6, -39]),
    ("Posterior cingulate gyrus (L)",         "Control",      -3.3, 0.025, [-5, -29, 28]),
    ("Parahippocampal gyrus (R)",             "Limbic",       -3.3, 0.034, [25, -11, -32]),
    ("Temporal pole (R)",                     "Limbic",       -3.2, 0.034, [30, 9, -38]),
    ("Posterior inferior temporal gyrus (L)", "Limbic",       -3.2, 0.034, [-45, -20, -30]),
    ("Posterior cingulate gyrus (R)",         "Control",      -3.1, 0.049, [5, -24, 31]),
    ("Posterior inferior temporal gyrus (R)", "Limbic",       -3.0, 0.049, [47, -12, -35]),
    ("Planum temporale (R)",                  "Sensorimotor", -3.0, 0.049, [64, -23, 8]),
    ("Superior lateral occipital cortex (L)", "Control",       3.3, 0.027, [-35, -62, 48]),
]

data = sorted(data, key=lambda d: d[2]) # most negative at top
regions = [d[0] for d in data]
y = list(range(len(data)))

fig, ax = plt.subplots(figsize=(8.6, 6.2), dpi=300)

COORD_X = 4.6 # x position of the MNI column
DARK = "#222629"

for yi, (reg, net, t, p, mni) in zip(y, data):
    ax.barh(yi, t, color=NET[net], edgecolor="white", linewidth=0.6, height=0.74, zorder=3)
    # t, p inside the bar near the baseline
    if t < 0:
        ax.text(-0.12, yi, f"{t:.1f}, p={p:.3f}", va="center", ha="right",
                fontsize=7.6, color=DARK, zorder=5)
    else:
        ax.text(0.12, yi, f"{t:.1f}, p={p:.3f}", va="center", ha="left",
                fontsize=7.6, color=DARK, zorder=5)
    # MNI coordinates as an aligned right-hand column
    ax.text(COORD_X, yi, f"[{mni[0]}, {mni[1]}, {mni[2]}]", va="center", ha="left",
            fontsize=8, color="#5d6469", zorder=5)

ax.text(COORD_X, -0.9, "MNI (x, y, z)", va="center", ha="left",
        fontsize=8.5, color="#5d6469", style="italic")

ax.axvline(0, color=DARK, lw=1.1, zorder=4)
ax.set_yticks(y)
ax.set_yticklabels(regions, fontsize=9)
ax.invert_yaxis()
ax.set_xlim(-5.8, 7.0)
ax.set_xticks([-4, -2, 0, 2])
ax.set_xlabel("t value", fontsize=11)
ax.set_ylabel("Region", fontsize=11)
ax.tick_params(axis="y", length=0)
for s in ["top", "right", "left"]:
    ax.spines[s].set_visible(False)
ax.grid(axis="x", color="#ebedee", lw=0.8, zorder=0)

# panel label
ax.text(-0.22, 1.02, "B", transform=ax.transAxes, fontsize=15, fontweight="bold",
        va="bottom", ha="left")

# horizontal legend at bottom
order = ["Sensorimotor", "Limbic", "Control", "Default Mode"]
handles = [Patch(facecolor=NET[n], edgecolor="white", label=n) for n in order]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.10),
          ncol=4, frameon=False, fontsize=9.5, handlelength=1.2, columnspacing=1.8)

plt.tight_layout()
plt.savefig("/home/claude/bar_fig.png", dpi=300, bbox_inches="tight", facecolor="white")


#!/usr/bin/env python3
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle, Rectangle
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Figures"
OUT.mkdir(parents=True, exist_ok=True)

blue, orange, green = "#2F6B9A", "#D9822B", "#2E8B57"
dark, pale = "#273444", "#F3F6F8"

fig, ax = plt.subplots(figsize=(110 / 25.4, 20 / 25.4), dpi=300)
ax.set_xlim(0, 110)
ax.set_ylim(0, 20)
ax.axis("off")
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Lineage branch
ax.add_patch(Circle((8.5, 10), 6.4, facecolor=blue, edgecolor="white", lw=0.8))
ax.text(8.5, 10, "hMSC", color="white", ha="center", va="center", fontsize=10, weight="bold")
for y, label, color in [(14.1, "Adipocyte", orange), (5.9, "Osteoblast", green)]:
    ax.add_patch(FancyArrowPatch((15.0, 10), (20.5, y), arrowstyle="-|>", mutation_scale=7,
                                 lw=1.0, color="#777777"))
    ax.add_patch(FancyBboxPatch((20.5, y - 2.6), 21, 5.2,
                                boxstyle="round,pad=0.25,rounding_size=1.2",
                                facecolor=color, edgecolor="none"))
    ax.text(31, y, label, color="white", ha="center", va="center", fontsize=10, weight="bold")

# Result-centered compartment remodeling summary.
ax.text(55.3, 17.2, "777 shifts", ha="center", va="center",
        fontsize=10, weight="bold", color="#8B2F2F")
for i in range(10):
    ax.add_patch(Rectangle((46 + i * 1.85, 8.0), 1.55, 6.0,
                           facecolor=blue if i in (0, 1, 5, 8) else orange,
                           edgecolor="white", lw=0.25))
ax.text(55.3, 4.6, "12/60 · 2.21×", ha="center", va="center",
        fontsize=10, weight="bold", color=dark)

ax.add_patch(FancyArrowPatch((65.5, 10.7), (72.5, 10.7), arrowstyle="-|>",
                             mutation_scale=8, lw=1.1, color="#777777"))

# Multi-omic prioritization: emphasize biological results rather than missingness.
ax.add_patch(FancyBboxPatch((72.8, 1.7), 35.8, 16.6,
                            boxstyle="round,pad=0.35,rounding_size=1.5",
                            facecolor=pale, edgecolor="#B8C2CC", lw=0.8))
ax.text(90.7, 14.7, "91 candidates", ha="center", va="center",
        fontsize=10, weight="bold", color="#8B2F2F")
ax.text(90.7, 10.2, "SFRP4 · SPP1", ha="center", va="center",
        fontsize=10, weight="bold", color=dark)
ax.text(90.7, 5.7, "TNFRSF11B", ha="center", va="center",
        fontsize=10, weight="bold", color=dark)

fig.subplots_adjust(0, 0, 1, 1)
png = OUT / "graphical_abstract.png"
tiff = OUT / "graphical_abstract.tiff"
fig.savefig(png, dpi=300, facecolor="white")
plt.close(fig)

with Image.open(png) as im:
    im = im.convert("RGB").resize((1300, 236), Image.Resampling.LANCZOS)
    im.save(png, dpi=(300, 300))
    rgb = Image.new("RGB", im.size, "white")
    if im.mode == "RGBA":
        rgb.paste(im, mask=im.getchannel("A"))
    else:
        rgb.paste(im.convert("RGB"))
    rgb.save(tiff, dpi=(300, 300), compression="tiff_lzw")

print(png)

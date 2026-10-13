#!/usr/bin/env python3
"""Convert the 600-dpi PNGs of Figures 1-5 (made with FIG_DPI=600) to LZW TIFF in figures/tiff/ (MDPI: >= 600 dpi)."""
import os, glob
from PIL import Image
D = os.path.dirname(os.path.abspath(__file__)); os.makedirs(f"{D}/tiff", exist_ok=True)
for n, pat in ((1, "Fig1_workflow"), (2, "Fig2_screen"), (3, "Fig3_trees"), (4, "Fig4_islands"), (5, "Fig5_boletales_genera")):
    src = f"{D}/{pat}.png"
    if not os.path.exists(src): print("missing", src); continue
    im = Image.open(src); dpi = im.info.get("dpi", (600, 600))
    im.convert("RGB").save(f"{D}/tiff/Figure{n}.tif", compression="tiff_lzw", dpi=(600, 600))
    print(f"Figure{n}.tif {im.size} dpi_in={dpi}")

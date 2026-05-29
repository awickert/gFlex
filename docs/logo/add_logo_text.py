#!/usr/bin/env python
"""Composite 'gFlex' text onto the rendered logo PNG.

Run after blender_logo.py:
    python docs/logo/add_logo_text.py

Reads  docs/_static/logo.png  (Blender render, RGBA)
Writes docs/_static/logo.png  — final logo with PRIMARY_FONT
"""

import os
from pathlib import Path  # used for FONT_PATH
from PIL import Image, ImageDraw, ImageFont

REPO    = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTDIR  = os.path.join(REPO, "docs", "_static")
LOGO_IN = os.path.join(OUTDIR, "logo.png")   # Blender render (no text)

FONT_PATH = str(Path.home() / ".local/share/fonts/AlbertSans/AlbertSans-Bold.ttf")


def add_text(img, font_path, font_size, pad_x, pad_y):
    font   = ImageFont.truetype(font_path, font_size)
    out    = img.copy()
    draw   = ImageDraw.Draw(out)
    stroke = max(2, font_size // 18)
    draw.text((pad_x, pad_y), "gFlex", font=font,
              fill=(255, 255, 255, 255),
              stroke_width=stroke,
              stroke_fill=(0, 0, 0, 255))
    return out


base  = Image.open(LOGO_IN).convert("RGBA")
W, H  = base.size
fsize = int(H * 0.12)
px    = int(W * 0.03)
py    = int(H * 0.03)

add_text(base, FONT_PATH, fsize, px, py).save(LOGO_IN)
print(f"✓  logo.png written  ({W}×{H}, {fsize}px)")

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).parent
BEFORE = ROOT / "before"
AFTER = ROOT / "after"
OUTPUT = ROOT / "comparison"
OUTPUT.mkdir(parents=True, exist_ok=True)
CASES = ("urea", "acetanilide", "primary-amine", "secondary-amine")

try:
    FONT = ImageFont.truetype("DejaVuSans.ttf", 16)
except OSError:
    FONT = ImageFont.load_default()

for case in CASES:
    for theme in ("light", "dark"):
        before = Image.open(BEFORE / f"{case}-before-{theme}.png").convert("RGB")
        after = Image.open(AFTER / f"{case}-after-{theme}.png").convert("RGB")
        if before.size != after.size:
            raise ValueError(f"before/after size mismatch: {case} {theme}")
        header = 34
        gap = 2
        comparison = Image.new(
            "RGB",
            (before.width + after.width + gap, before.height + header),
            "#d7dce2",
        )
        draw = ImageDraw.Draw(comparison)
        draw.text((12, 8), f"BEFORE — {case} — {theme}", fill="#1f2937", font=FONT)
        draw.text(
            (before.width + gap + 12, 8),
            f"AFTER — {case} — {theme}",
            fill="#1f2937",
            font=FONT,
        )
        comparison.paste(before, (0, header))
        comparison.paste(after, (before.width + gap, header))
        comparison.save(OUTPUT / f"{case}-{theme}-before-after.png")

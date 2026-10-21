from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).parent
MENUS = ("energy-palette", "mo-presets", "ligand-presets")


def _font():
    try:
        return ImageFont.truetype("DejaVuSans.ttf", 14)
    except OSError:
        return ImageFont.load_default()


def _compare(theme: str, scale: str, menu: str) -> None:
    suffix = f"-{theme}.png"
    before_dir = ROOT / "before" / ("scale-1.5" if scale == "1.5" else "")
    after_dir = ROOT / "after" / ("scale-1.5" if scale == "1.5" else "")
    before = Image.open(before_dir / f"{menu}{suffix}").convert("RGBA")
    after = Image.open(after_dir / f"{menu}{suffix}").convert("RGBA")

    pad = 12
    header = 30
    gap = 16
    font = _font()
    before_heading = f"ANTES · {theme} · {scale}x"
    after_heading = f"DESPUÉS · {theme} · {scale}x"
    measure = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    heading_width = max(
        measure.textbbox((0, 0), before_heading, font=font)[2],
        measure.textbbox((0, 0), after_heading, font=font)[2],
    )
    panel_width = max(before.width, after.width, heading_width) + 2 * pad
    panel_height = max(before.height, after.height) + 2 * pad
    canvas = Image.new(
        "RGBA",
        (2 * panel_width + gap, panel_height + header + pad),
        (238, 242, 247, 255),
    )
    draw = ImageDraw.Draw(canvas)
    draw.text((pad, 7), before_heading, fill=(24, 36, 56, 255), font=font)
    draw.text(
        (panel_width + gap + pad, 7),
        after_heading,
        fill=(24, 36, 56, 255),
        font=font,
    )
    canvas.alpha_composite(before, (pad, header + pad))
    canvas.alpha_composite(after, (panel_width + gap + pad, header + pad))

    out_dir = ROOT / "comparison"
    out_dir.mkdir(parents=True, exist_ok=True)
    scale_suffix = "" if scale == "1.0" else "-scale-1.5"
    canvas.convert("RGB").save(
        out_dir / f"{menu}-{theme}{scale_suffix}-before-after.png"
    )


def main() -> None:
    for theme in ("light", "dark"):
        for scale in ("1.0", "1.5"):
            for menu in MENUS:
                _compare(theme, scale, menu)


if __name__ == "__main__":
    main()

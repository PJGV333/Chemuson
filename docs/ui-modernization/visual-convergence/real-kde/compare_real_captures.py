"""Coarse region/ink comparison for real Qt Wayland grabs vs golden captures."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
GOLD = HERE.parent / "captures"

def read(path):
    return np.array(Image.open(path).convert("RGB"))

def measure(img, box, dark=False):
    x0,y0,x1,y1=box
    a=img[y0:y1,x0:x1]
    # Select saturated/dark-enough ink, excluding neutral backgrounds.
    if dark:
        mask=(a.max(axis=2)<170) & ((a.max(axis=2)-a.min(axis=2))>10)
    else:
        mask=(a.max(axis=2)<150) & ((a.max(axis=2)-a.min(axis=2))>10)
    ys,xs=np.where(mask)
    if not len(xs): return {"ink":0,"bbox":None}
    return {"ink":int(len(xs)),"bbox":[int(xs.min()),int(ys.min()),int(xs.max()),int(ys.max())],
            "region_wh":[x1-x0,y1-y0]}

def main():
 out={}
 for theme in ("light","dark"):
  real=read(HERE/f"real-qt-wayland-{theme}.png")
  gold=read(GOLD/f"prod-full-{theme}.png")
  # normalize screenshot layout to the known 1440x900 reference via resize;
  # evaluate shell regions only (the canvas intentionally differs).
  R=np.array(Image.fromarray(real).resize((1440,900),Image.Resampling.LANCZOS))
  regions={"appbar_logo":(50,0,230,54),"rail_first6":(0,54,58,420),
           "rail_last_icons":(0,580,58,866),"appbar_controls":(1240,0,1440,54)}
  out[theme]={"real_size":list(Image.open(HERE/f"real-qt-wayland-{theme}.png").size),"regions":{}}
  for name,box in regions.items():
   out[theme]["regions"][name]={"real":measure(R,box,theme=="dark"),"gold":measure(gold,box,theme=="dark")}
   out[theme]["regions"][name]["mean_abs_rgb_diff"]=float(np.abs(R[box[1]:box[3],box[0]:box[2]].astype(int)-gold[box[1]:box[3],box[0]:box[2]].astype(int)).mean())
 print(json.dumps(out,indent=2))
 (HERE/"real-capture-comparison.json").write_text(json.dumps(out,indent=2)+"\n")
if __name__=='__main__': main()

# supersplat-snap

## 🎥 Demo

Click the thumbnail below to watch the demo:

[![Demo Video](https://img.youtube.com/vi/W5FI3u79UrY/maxresdefault.jpg)](https://www.youtube.com/shorts/W5FI3u79UrY)

A **Splat Align** tool for [SuperSplat](https://github.com/playcanvas/supersplat) — the open-source 3D Gaussian Splat editor by PlayCanvas.

Adds a point-correspondence alignment tool that lets you place matching control points across two or more splat files and automatically compute the best-fit transform to align them.

> Built as a patch on top of SuperSplat. You clone SuperSplat separately and drop these files in — no fork needed.

---

## What it does

- Place **matching control points** on two or more splats (A1↔B1, A2↔B2, etc.)
- Computes an optimal **rotation + translation + scale** transform using the [Umeyama similarity method](https://en.wikipedia.org/wiki/Kabsch_algorithm) — least-squares best fit across all point pairs simultaneously
- **All secondary splats align to Splat A** (the first one you place points on)
- Points are **editable** after placement — drag X/Y/Z values to fine-tune
- After aligning, points move with the splat so you can **re-align iteratively**
- Fully **undoable** with Ctrl+Z
- Works with `.splat` and `.ply` files
- Supports **3+ splats** (B, C, D… all align to A)

## Tips for good alignment

> ⚠️ **Do not place all points on a flat plane** (e.g. all on the ground). Three coplanar points cannot define a unique 3D rotation. Place points at **different heights**:
> - Top of a structure or pole
> - A point on the ground
> - Top of a bench or wall
> - A tree base vs a higher branch

The more spread out in 3D space, the more accurate the alignment. 4–6 well-distributed points is ideal.

---

## Setup
You will need Node.JS installed also, download [here](https://nodejs.org/en/download)

### 1. Clone and set up SuperSplat

```bash
git clone https://github.com/playcanvas/supersplat.git
cd supersplat
npm install
```

### 2. Download and apply this patch

```bash
# Clone this repo somewhere
git clone https://github.com/VFX-Soup/supersplat-snap.git

# Copy the new/modified files into your supersplat folder
cp supersplat-align/src/tools/align-tool.ts  supersplat/src/tools/
cp supersplat-align/src/camera.ts             supersplat/src/
cp supersplat-align/src/main.ts               supersplat/src/
cp supersplat-align/src/ui/bottom-toolbar.ts  supersplat/src/ui/
```

### 3. Run

```bash
cd supersplat
npm run develop
```

Open [http://localhost:3000](http://localhost:3000) and hard-refresh (Ctrl+Shift+R).

---

## Usage

1. Load two or more splat files via the SuperSplat **Open** menu
2. Click the **⊕** (Align) icon in the bottom toolbar — it appears after the ruler/measure tool
3. In the **Scene Manager** (left panel), select **Splat A** (your base/reference)
4. Click on the viewport to place points **A1, A2, A3…** on recognisable features
5. Select **Splat B** in Scene Manager
6. Click matching features to place **B1, B2, B3…** in the same order
7. Repeat for any additional splats (C, D…)
8. Click **Align** — B (and C, D…) will transform to match A
9. Check the residual error shown in the panel. If it's high, try adding more points or adjusting point positions using the X/Y/Z drag inputs
10. Export as normal using SuperSplat's built-in export

### Panel controls

| Control | Action |
|---------|--------|
| Click viewport | Place a point on the selected splat |
| Drag X/Y/Z input | Scrub value (Shift = fine, Alt = coarse) |
| Undo Last | Remove the most recently placed point |
| Clear | Remove all points |
| Align | Compute and apply the transform |

---

## Files changed from SuperSplat

| File | Type | What changed |
|------|------|-------------|
| `src/tools/align-tool.ts` | **New** | The entire align tool (~400 lines) |
| `src/camera.ts` | Modified | Added `intersectSplat()` method — picks depth on a specific splat only (essential for overlapping splats) |
| `src/main.ts` | Modified | +3 lines: import + register AlignTool |
| `src/ui/bottom-toolbar.ts` | Modified | +5 lines: button + click handler + active state + tooltip |

---

## Compatibility

Tested against SuperSplat **v2.27.3**. If you're on a different version and the modified files cause issues, you can apply the changes to `main.ts` and `bottom-toolbar.ts` manually — the changes are minimal and described in the table above.

---

## Licence

`align-tool.ts` is original work, released under the **MIT Licence**.

The other modified files are based on SuperSplat source code © PlayCanvas, also MIT licensed.

```
MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software to deal in the Software without restriction, including the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
```

---

## Credits

- [SuperSplat](https://github.com/playcanvas/supersplat) by [PlayCanvas](https://playcanvas.com) — the editor this patches into
- Alignment uses the [Umeyama (1991)](https://web.stanford.edu/class/cs273/refs/umeyama.pdf) least-squares similarity transform method

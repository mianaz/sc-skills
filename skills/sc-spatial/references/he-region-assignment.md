# H&E segmentation → Visium spot-to-region assignment

Segment tissue regions from H&E images, then assign each Visium spot to its enclosing region. Useful for TLS detection, tumor/stroma annotation, or necrosis delineation. Per-sample threshold tuning is required.

(source: sc-paper-distill/papers/2026-pan-cancer-TLS.md M2)

```python
from skimage import color, morphology, measure, segmentation
from skimage.transform import rescale
import json, numpy as np
from skimage import io

img = io.imread("tissue_hires_image.png")
lum = color.rgb2gray(img)

# Threshold (tune per sample)
thre_low, thre_high = 0.2, 0.4
mask = np.logical_and(lum > thre_low, lum < thre_high)

# Morphological cleanup
mask = morphology.remove_small_holes(
    morphology.remove_small_objects(mask, min_size=500), area_threshold=600)
mask = morphology.erosion(mask, morphology.disk(5))

# Label connected components + rescale to spot coordinates
labels = measure.label(mask, connectivity=2)
scale_factor = json.load(open("scalefactors_json.json"))["tissue_hires_scalef"]
labels_scaled = rescale(labels, 1/scale_factor, anti_aliasing=False, order=0)
labels_expanded = segmentation.expand_labels(
    measure.label(labels_scaled), distance=50)

# Assign spots: spot_df has columns X, Y (pixel coords from tissue_positions)
spot_df["region"] = [labels_expanded[y, x] for y, x in zip(spot_df.Y, spot_df.X)]
```

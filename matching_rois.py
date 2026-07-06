import matplotlib.pyplot as plt

from helpers import *
from preprocessing import *
from track2p_preprocessing import *
from suite2p_preprocessing import *
from figures import *
from cellreg_preprocessing import *
from plot_metrics import *

#do everything with non deconvolved data (f-0.7*fneu)
# animal = 'EC_EBc_02'
# path = r'E:\dark_drift\data'
# track2p_folder = 'track2p-grat-functional'
# suite2p_obj = suite2pPreprocessing(animal,'grat', path, track2p_folder, ntheta = 8, fps = 20, tracked_cells = False, roi_detection = 'functional',deconvolved = True, zscore_threshold=1, std_threshold=None)
#

animals = ['EC_EB_08','EC_EBc_02']
animal_dobs = {'EC_EBc_02': '20250429',
                'EC_EB_01': '20250331',
               'EC_EB_08': '20250401',
               'EC_EB_03': '20250331',}

good_tracked_cells = {'EC_EBc_02': np.array([2, 4, 6, 8, 11, 12, 15,16, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 37, 41, 42, 45, 46, 51, 52, 53, 54, 55, 56, 57, 58, 60, 63, 64, 66, 67, 68, 71, 74, 75, 78, 79, 83, 85, 86, 88, 92, 93, 97, 99, 101, 102, 105, 106, 107, 108, 110, 115, 116, 117, 118, 119, 120, 122, 124,125,126, 128, 129,  131, 132, 135, 137, 138, 139, 140, 141, 144, 145, 147, 154, 159, 166, 168, 172, 176, 178, 183, 187, 188, 193, 194, 195, 200, 201, 204, 206, 207, 213, 215, 216, 217, 219, 227, 230, 231, 232, 233, 235, 238, 241, 246, 247, 249, 260, 265, 266, 270, 271, 273, 279, 280, 283, 286, 293, 295, 296, 299, 305, 306, 307, 319, 320, 321, 322, 324, 327, 328, 329]),
                      'EC_EB_08': np.array([0, 6, 7, 9, 11, 16, 25, 29, 30, 39, 40, 64, 65, 70, 73, 78,82, 83, 84, 92, 98, 102, 104, 110, 113, 119, 122,123, 126, 150, 153, 157, 161, 164])}

path = r'E:\dark_drift\data'
data_object = batchProcessing(animals, animal_dobs, path,
                              stim = 'grat',
                              tracked_cells = True,
                              checked_tracked_cells = good_tracked_cells,
                              roi_detection = 'functional',
                              ntheta = 8,
                              fps = 20,
                              zscore_threshold=0.8,
                              deconvolved = True)

animal = 'EC_EBc_02'
track_obj = data_object.dat[animal].track2p_obj

# build_roi_masks(ops, stat, iscell)
# overlay_masks_on_mean(self.meanimg,
#                       self.m[0][:,:,:20],
#                       labels= self.m[1][:20])

# FOR EACH ROI TRACKED, PLOT WHERE IT IS

good_tracked_cells = {'EC_EBc_02': np.array([2, 4, 6, 8, 11, 12, 15,16, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 37, 41, 42, 45, 46, 51, 52, 53, 54, 55, 56, 57, 58, 60, 63, 64, 66, 67, 68, 71, 74, 75, 78, 79, 83, 85, 86, 88, 92, 93, 97, 99, 101, 102, 105, 106, 107, 108, 110, 115, 116, 117, 118, 119, 120, 122, 124,125,126, 128, 129,  131, 132, 135, 137, 138, 139, 140, 141, 144, 145, 147, 154, 159, 166, 168, 172, 176, 178, 183, 187, 188, 193, 194, 195, 200, 201, 204, 206, 207, 213, 215, 216, 217, 219, 227, 230, 231, 232, 233, 235, 238, 241, 246, 247, 249, 260, 265, 266, 270, 271, 273, 279, 280, 283, 286, 293, 295, 296, 299, 305, 306, 307, 319, 320, 321, 322, 324, 327, 328, 329]),
                      'EC_EB_08': np.array([0, 6, 7, 9, 11, 16, 25, 29, 30, 39, 40, 64, 65, 70, 73, 78,82, 83, 84, 92, 98, 102, 104, 110, 113, 119, 122,123, 126, 150, 153, 157, 161, 164])}

fig, ax = plt.subplots (2,2, figsize = (15,10))
ax = ax.ravel()
for i_day in range(len(track_obj.meanImg)):
    ax[i_day].imshow(track_obj.meanImg[i_day], cmap = 'gray')
    ax[i_day].axis ('off')

    masks, labels_idx = build_roi_masks_tracked(track_obj.all_ops[i_day], track_obj.all_stat_t2p[i_day])

    start, end = 60,120 # to do
    # cells = np.array([24])
    # overlay_masks_on_mean(track_obj.meanImg[i_day],
    #                       masks[:,:,start:end],
    #                       labels= labels_idx[start:end],
    #                       alpha = 0.15,
    #                       ax = ax[i_day])
    cells = good_tracked_cells [animal]
    overlay_masks_on_mean(track_obj.meanImg[i_day],
                          masks[:,:,:],
                          labels= None,
                          alpha = 0.15,
                          ax = ax[i_day])


    # overlay_masks_on_mean(track_obj.meanImg[i_day],
    #                       masks[:,:,good_tracked_cells['EC_EBc_02']],
    #                       labels= None,
    #                       alpha = 0.2,
    #                       draw_outlines=True,
    #                       ax = ax[i_day])

plt.suptitle(animal)
plt.tight_layout()
plt.show()

indices = good_tracked_cells['EC_EBc_02']
n_cells = indices.max() + 1
binary_array = np.zeros(n_cells, dtype=int)
binary_array[indices] = 1
plt.plot (np.cumsum(binary_array))


def build_roi_masks_tracked(ops, stat):
    """
    Build boolean masks for ROIs listed in `keep_idx`.
    Returns (mask_stack [Ly,Lx,N], labels [N])
    """
    Ly, Lx = ops['Ly'], ops['Lx']
    #print(stat)
    keep_idx = np.arange(len(stat))

    mask_stack = np.zeros((Ly, Lx, len(keep_idx)), dtype=bool)
    for i, k in enumerate(keep_idx):
        y = np.asarray(stat[k]['ypix']).astype(int)
        x = np.asarray(stat[k]['xpix']).astype(int)
        m = np.zeros((Ly, Lx), dtype=bool)
        valid = (y >= 0) & (y < Ly) & (x >= 0) & (x < Lx)
        m[y[valid], x[valid]] = True
        mask_stack[..., i] = m
    return mask_stack, np.asarray(keep_idx)

import numpy as np
import matplotlib.pyplot as plt

def overlay_masks_on_mean(mean_img,
                          mask_stack,
                          labels=None,
                          *,
                          ax=None,
                          alpha=0.25,
                          draw_outlines=False,
                          draw_fill=True,
                          label_every=1,
                          one_indexed=False,
                          dot_centroid=False,
                          dot_size=20):
    """
    Overlay ROI masks on a mean image.
    mask_stack: (Ly, Lx, N) boolean
    If ax is provided, draws into it; otherwise creates a new fig/ax.
    """
    Ly, Lx = mean_img.shape
    cmap = "gray"

    # If user didn’t provide an axis, make a fresh one
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7), dpi=120)
    else:
        fig = ax.figure

    ax.imshow(mean_img, cmap=cmap, interpolation="nearest")
    ax.set_axis_off()

    # colors for N ROIs
    N = mask_stack.shape[-1]
    hues = np.linspace(0, 1, max(1, N), endpoint=False)
    colors = plt.cm.hsv(hues)[:, :3]  # RGB only

    # --------- FILLED OVERLAY (single draw) ---------
    if draw_fill and N > 0:
        lab = np.full((Ly, Lx), -1, dtype=int)
        for i in range(N):
            m = mask_stack[..., i]
            lab[m] = i
        overlay = np.zeros((Ly, Lx, 4), dtype=float)
        covered = lab >= 0
        overlay[covered, :3] = colors[lab[covered]]
        overlay[covered, 3] = alpha
        ax.imshow(overlay, interpolation="nearest")

    # --------- OUTLINES (single draw) ---------
    if draw_outlines:
        if not draw_fill:  # build lab if not already made
            lab = np.full((Ly, Lx), -1, dtype=int)
            for i in range(N):
                lab[mask_stack[..., i]] = i

        edge = np.zeros((Ly, Lx), dtype=bool)
        for dy, dx in ((-1,0),(1,0),(0,-1),(0,1)):
            nbr = np.roll(lab, shift=(dy, dx), axis=(0,1))
            edge |= (lab != nbr)
        edge &= (lab >= 0)

        if edge.any():
            edge_rgba = np.zeros((Ly, Lx, 4), dtype=float)
            edge_rgba[..., :3] = colors[lab.clip(min=0)]
            edge_rgba[..., 3] = edge.astype(float)
            ax.imshow(edge_rgba, interpolation="nearest")

    # --------- LABELS ---------
    if labels is not None and label_every > 0:
        yy, xx = np.indices((Ly, Lx))
        for i in range(0, N, label_every):
            m = mask_stack[..., i]
            if not m.any():
                continue
            y = yy[m].mean() + 8
            x = xx[m].mean() + 8
            text = str(int(labels[i]) + (1 if one_indexed else 0))
            ax.text(x, y, text,
                    ha="center", va="center", fontsize=7, color="w",
                    bbox=dict(facecolor="k", alpha=0.55, boxstyle="round,pad=0.15"),
                    clip_on=True)
            if dot_centroid:
                ax.scatter([x], [y], s=dot_size, c=[colors[i]],
                           edgecolors="k", linewidths=0.3)

    return fig, ax


def create_roi_mask(stat, ops):
    masks = []
    Ly, Lx = ops['Ly'], ops['Lx']
    for roi in stat:
        mask = np.zeros((Ly, Lx))
        ypix = roi['ypix']
        xpix = roi['xpix']
        weights = roi['lam']
        mask[ypix, xpix] = weights
        masks.append(mask)
    return np.array(masks)


def build_roi_masks(ops, stat, iscell):
    """
    Build boolean masks for ROIs listed in `keep_idx`.
    Returns (mask_stack [Ly,Lx,N], labels [N])
    """
    Ly, Lx = ops['Ly'], ops['Lx']
    keep_idx = np.where(iscell[:, 0] == 1)[0]

    mask_stack = np.zeros((Ly, Lx, len(keep_idx)), dtype=bool)
    for i, k in enumerate(keep_idx):
        y = np.asarray(stat[k]['ypix']).astype(int)
        x = np.asarray(stat[k]['xpix']).astype(int)
        m = np.zeros((Ly, Lx), dtype=bool)
        valid = (y >= 0) & (y < Ly) & (x >= 0) & (x < Lx)
        m[y[valid], x[valid]] = True
        mask_stack[..., i] = m
    return mask_stack, np.asarray(keep_idx)
import numpy as np

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

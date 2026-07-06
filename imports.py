import numpy as np
import matplotlib
# matplotlib.use('TkAgg')
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.cm as cm
from scipy.spatial.distance import cdist

import os
os.environ["QT_API"] = "pyqt5"      # helps qtpy pick PyQt5 if qtpy is present
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt


import os
import os.path
import ast
import random
from scipy.io import savemat
from itertools import chain
from scipy.stats import mannwhitneyu
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D  # required for 3D projection
from matplotlib.cm import get_cmap
import scipy
from scipy.io import loadmat
from scipy.io import savemat
from scipy.stats import sem
from scipy.stats import ttest_ind
from scipy.linalg import svd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import zscore
from scipy.stats import pearsonr
from datetime import datetime, timedelta

from joblib import load, dump
import imageio
import cv2
import random
import math
from types import SimpleNamespace
import colorednoise as cn
from PIL import Image
from tqdm import tqdm

from scipy.ndimage import median_filter
from scipy import signal
from scipy.io import loadmat
from scipy.ndimage import median_filter
from scipy import signal
from scipy.io import loadmat
from skimage import exposure
from pathlib import Path
from scipy.stats import ks_2samp
from scipy.stats import linregress
from matplotlib.patches import Ellipse
import matplotlib.colors as mcolors
from scipy.stats import ks_2samp, ttest_ind
from scipy.stats import kruskal


from sklearn.preprocessing import minmax_scale
from sklearn.preprocessing import StandardScaler
import skimage
import skimage.transform as transform
from skimage.measure import block_reduce
import skimage.transform as transform
from skvideo.io import vreader

from sklearn.decomposition import PCA
from itertools import combinations
from sklearn.preprocessing import normalize

import re
from datetime import datetime
import umap.umap_ as umap
from scipy.signal import find_peaks, peak_widths
import tifffile as tiff
import imageio.v2 as imageio
import glob
from scipy.ndimage import uniform_filter1d  # temporal smoothing
from sbxreader import sbx_memmap


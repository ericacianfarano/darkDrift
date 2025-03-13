import numpy as np
import matplotlib
# matplotlib.use('TkAgg')
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

import os
import os.path
import ast
import random

import scipy
from scipy.io import loadmat
from scipy.io import savemat
from scipy.stats import sem
from scipy.stats import ttest_ind
from scipy.linalg import svd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import zscore
from scipy.stats import pearsonr

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


from sklearn.preprocessing import minmax_scale
from sklearn.preprocessing import StandardScaler
import skimage
import skimage.transform as transform
from skimage.measure import block_reduce
import skimage.transform as transform
from skvideo.io import vreader


import re


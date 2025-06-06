from helpers import *
from preprocessing import *
from track2p_preprocessing import *
from suite2p_preprocessing import *
from figures import *

#do everything with non deconvolved data (f-0.7*fneu)

animal = 'EC_GCaMP6s_09'
path = r'I:\dark_drift_trial\data'
track2p_folder = 'track2p - 10 weeks'
suite2p_obj = suite2pPreprocessing(animal, path, track2p_folder, ntheta = 8, fps = 20, deconvolved = True)

# figures for proposal
fov_across_days(suite2p_obj, brightness = 0.5, contrast = 2.2)
hist_osi_angle (suite2p_obj)        # use deconvolved
plot_response(suite2p_obj, cell_i = 5) # use fluorescence traces
rasters_across_days(suite2p_obj)
polar_plots_across_days_rois(suite2p_obj, cells_to_plot = [10,21,75,103])
polar_plots_across_days(suite2p_obj, cell = 5)
plot_rois_across_days (suite2p_obj, suite2p_obj.track2p_obj, 65, n_cells_to_plot = 2, brightness = 0.5, contrast = 1) # just keep reruning this if it doesn't plot


plot_raw_responses (suite2p_obj, animal)
plot_response(suite2p_obj, cell_i = 5)

# to upload
# >> git status
# >> git add . (stage all updated files)
# >> git add filename.py (stage a specific file)
# >> git commit -m "Describe your change here"
# git push
# to stage and commit all in one line: git commit -am "Quick update"


plot_corr(suite2p_obj,n = 1000, across_days = True)



# step 1) process calcium imaging sbx files with suite2p
    # in anaconda prompt / terminal: type 'conda activate suite2p' to activate suite2p enviroinment
    # in anaconda prompt / terminal: type 'suite2p' to open the gui

# step 2) once an entire dataset is collected, process suite2p files with track2p to find same cells across days
    # in pycharm terminal: conda activate track2p
    # in pycharm terminal: set paths & settings > type 'python -m run_track2p' to run track2p algorithm
    # output will be saved in track2p folder
    # run track2pPreprocessing class (in track2p_preprocessing.py file) to process and load matched/tracked cells across days

    # track2p_obj = track2pPreprocessing(path, animal, 'track2p',single_plane = True, deconvolved = True)

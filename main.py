from helpers import *
from preprocessing import *
from track2p_preprocessing import *
from suite2p_preprocessing import *
from figures import *

#do everything with non deconvolved data (f-0.7*fneu)

# animal = 'EC_dark_01'
# path = r'I:\dark_drift\data'
# track2p_folder = 'track2p'
# suite2p_obj = suite2pPreprocessing(animal, path, track2p_folder, ntheta = 8, fps = 20, deconvolved = True)

animals = ['EC_dark_01', 'EC_dark_03', 'EC_dark_05', 'EC_dark_08'] #['EC_dark_01', 'EC_dark_03'] #'EC_dark_05', 'EC_dark_08']

animal_dobs = {'EC_dark_01': '20250331',
               'EC_dark_03': '20250331',
               'EC_dark_05': '20250401',
               'EC_dark_08': '20250401',
               }

path = r'I:\dark_drift\data'
data_object = batchProcessing(animals, animal_dobs, path,
                              stim = 'grat',
                              tracked_cells = True,
                              ntheta = 8,
                              fps = 20,
                              deconvolved = True)

animal = 'EC_dark_01'
if (data_object.tracked_cells) and (data_object.stim == 'grat'):
    for animal in data_object.list_animals:
        fov_across_days(data_object, brightness = 0.5, contrast = 2.2)
        plot_response(data_object, animal)
        rasters_across_days(data_object, animal)
        polar_plots_across_days(data_object, animal)
        polar_plots_across_days_rois(data_object, animal)
        #plot_corr(data_object, animal, n=1000, across_days=False) # includes within-day shuffle
        # add activity threshold

elif (not data_object.tracked_cells) and (data_object.stim == 'grat'):
    hist_osi_angle(data_object)  # use deconvolved
    polar_plots(data_object, animal) #tuning curves of all rois

elif (data_object.stim == 'spon')
    var_explained, slope_eigenvals = spontaneous_analysis_slope(data_object)
    del slope_eigenvals['dark']['P70_s']
    plot_slope_eigenvals(slope_eigenvals)


plot_rois_across_days (data_object,animal, 65, n_cells_to_plot = 2, brightness = 0.5, contrast = 1) # just keep reruning this if it doesn't plot
plot_raw_responses (data_object, animal)
plot_response(data_object,animal, cell_i = 5)
plot_corr(data_object, animal, n = 3, across_days = False)

# stuff to analyze: # cells that are sileneced/inhibited by visual stims, are cells that are spontaneously active just as active with gratings?
# to do: verify tracking of ROIs, output all of the tuning curves and polar plots, automatially save correlation plots

# to upload
# >> git status
# >> git add . (stage all updated files)
# >> git add filename.py (stage a specific file)
# >> git commit -m "Describe your change here"
# git push
# to stage and commit all in one line: git commit -am "Quick update"





# step 1) process calcium imaging sbx files with suite2p
    # in anaconda prompt / terminal: type 'conda activate suite2p' to activate suite2p enviroinment
    # in anaconda prompt / terminal: type 'suite2p' to open the gui

# step 2) once an entire dataset is collected, process suite2p files with track2p to find same cells across days
    # in pycharm terminal: conda activate track2p
    # in pycharm terminal: set paths & settings > type 'python -m run_track2p' to run track2p algorithm
    # output will be saved in track2p folder
    # run track2pPreprocessing class (in track2p_preprocessing.py file) to process and load matched/tracked cells across days

    # track2p_obj = track2pPreprocessing(path, animal, 'track2p',single_plane = True, deconvolved = True)

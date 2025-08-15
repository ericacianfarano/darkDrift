from helpers import *
from preprocessing import *
from track2p_preprocessing import *
from suite2p_preprocessing import *
from figures import *
from cellreg_preprocessing import *
from plot_metrics import *

#do everything with non deconvolved data (f-0.7*fneu)

# animal = 'EC_dark_01'
# path = r'I:\dark_drift\data'
# track2p_folder = 'track2p'
# suite2p_obj = suite2pPreprocessing(animal, path, track2p_folder, ntheta = 8, fps = 20, deconvolved = True)

animals = ['EC_dark_01', 'EC_dark_03', 'EC_dark_05', 'EC_dark_08', 'EC_ctrl_10', 'EC_ctrl_11', 'EC_ctrl_12'] #['EC_dark_01', 'EC_dark_03'] #'EC_dark_05', 'EC_dark_08']

animal_dobs = {'EC_dark_01': '20250331',
               'EC_dark_03': '20250331',
               'EC_dark_05': '20250401',
               'EC_dark_08': '20250401',
               'EC_dark_14': '20250430',
               'EC_ctrl_10': '20250429',
               'EC_ctrl_11': '20250429',
               'EC_ctrl_12': '20250430',
               }

# lets add a condition where we only load X animals when tracking is good and analysis is tracking based
# vs just loading data from a particular day for each animal, if tracking doesnt matter
# or loading all data we possibly have for example spontaneous analysis

path = r'E:\dark_drift\data'
data_object = batchProcessing(animals, animal_dobs, path,
                              stim = 'grat',
                              tracked_cells = False,
                              roi_detection = 'functional',
                              ntheta = 8,
                              fps = 20,
                              zscore_threshold=0.8,
                              deconvolved = True)


animal = 'EC_dark_01'
if (data_object.tracked_cells) and (data_object.stim == 'grat'):
    for animal in data_object.list_animals:
        #fov_across_days(data_object, brightness = 0.5, contrast = 2.2)
        data_object.dat[animal].corr_matrix = plot_corr(data_object, animal, thresholded_cells=0, n=1,
                                                        across_days=False)  # includes within-day shuffle
        plot_response(data_object, animal)
        rasters_across_days(data_object, animal)
        polar_plots_across_days(data_object, animal)
        polar_plots_across_days_rois(data_object, animal)

        plot_corr_change(corr)

        plot_response(data_object, animal)
        # look at cells that are responsive at least once + cell that are responsive all days
        # add activity threshold


elif (not data_object.tracked_cells) and (data_object.stim == 'grat'):
    for animal in data_object.list_animals:
        #hist_osi_angle(data_object.dat[animal])  # use deconvolved
        #polar_plots(data_object, animal) #tuning curves of all rois
        responsiveness(data_object, 'P70')
        responsive_cells_across_days(data_object)
        trial_by_trial_reliability(data_object, 'P70')
        plot_avg_response(data_object)
        plot_avg_response_time(data_object)


elif (data_object.stim == 'spon')
    var_explained, slope_eigenvals = spontaneous_analysis_slope(data_object)
    del slope_eigenvals['dark']['P70_s']
    del slope_eigenvals['ctrl']['P70_s']
    plot_slope_eigenvals(slope_eigenvals)


def plot_corr_change(corr_matrix):
    '''
    :param corr_matrix: output of plot_corr function
    :return:
    '''
    upper_diag = np.array([np.diag(corr_matrix[..., i], k=1) for i in range(corr_matrix.shape[-1])])

    fig, ax = plt.subplots()
    for vec in upper_diag:
        ax.plot(np.arange(upper_diag.shape[-1]), vec, c = 'grey', alpha = 0.6)
    ax.plot(np.arange(upper_diag.shape[-1]), upper_diag.mean(axis= 0), c='blue')
    ax.set_xticks(np.arange(upper_diag.shape[-1]))
    ax.set_xticklabels(['P70×P77', 'P77×P84', 'P84×91'])
    ax.set_xlabel ('Week comparison')
    ax.set_ylabel('Pearson Correlation')
    plt.show()


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

import matplotlib.pyplot as plt
from helpers import *
from preprocessing import *
from track2p_preprocessing import *
from suite2p_preprocessing import *
from figures import *
from cellreg_preprocessing import *
from plot_metrics import *
from animal_info import *
from population_helpers import *

#do everything with non deconvolved data (f-0.7*fneu)
# load_stim = 'grat'
# tracked_cells = True
# tracking_across_days = True # if False, track within 1 day
# group = 'LBc'

load_stim = 'grat'
tracking_across_days = 1 # if False, track within 1 day
tracked_cells = False
group = None

if tracking_across_days and tracked_cells: # track cells across weeks
    animals = ['EC_EB_03', 'EC_EB_14', 'EC_EB_23', 'EC_EB_28', 'EC_EB_29', 'EC_EB_27', 'EC_EB_31', 'EC_EB_32','EC_EB_33',
                         'EC_EBc_02', 'EC_EBc_03', 'EC_EBc_04', 'EC_EBc_05', 'EC_EBc_06', 'EC_EBc_07', 'EC_EBc_09',
                         'EC_LB_01', 'EC_LB_02', 'EC_LB_03', 'EC_LB_05', 'EC_LB_06', 'EC_LB_09', 'EC_LB_10', 'EC_LB_11',
                         'EC_LBc_02','EC_LBc_03', 'EC_LBc_06', 'EC_LBc_09', 'EC_LBc_11', 'EC_LBc_14', 'EC_LBc_16', 'EC_LBc_19', 'EC_LBc_20'] # still need to do LB group

elif not tracked_cells:     # don't track cells
    animals = ['EC_EB_01', 'EC_EB_03', 'EC_EB_05', 'EC_EB_08', 'EC_EB_10', 'EC_EB_13', 'EC_EB_14', 'EC_EB_21','EC_EB_22','EC_EB_23','EC_EB_24','EC_EB_27','EC_EB_28','EC_EB_29','EC_EB_31', 'EC_EB_32', 'EC_EB_33',
               'EC_EBc_01', 'EC_EBc_02', 'EC_EBc_03', 'EC_EBc_04','EC_EBc_05','EC_EBc_06','EC_EBc_07', 'EC_EBc_09','EC_EBc_11',
               'EC_LB_01','EC_LB_02','EC_LB_03','EC_LB_05','EC_LB_06','EC_LB_07','EC_LB_09','EC_LB_10', 'EC_LB_11',
               'EC_LBc_02','EC_LBc_03', 'EC_LBc_06', 'EC_LBc_09', 'EC_LBc_11', 'EC_LBc_13', 'EC_LBc_14', 'EC_LBc_16', 'EC_LBc_19', 'EC_LBc_20']

elif not tracking_across_days and tracked_cells: # track cells within a day
    animals = ['EC_EB_01', 'EC_EB_03', 'EC_EB_05', 'EC_EB_08', 'EC_EB_10', 'EC_EB_13', 'EC_EB_14', 'EC_EB_21','EC_EB_22','EC_EB_23','EC_EB_24','EC_EB_27','EC_EB_28','EC_EB_29','EC_EB_31','EC_EB_32','EC_EB_33',
               'EC_EBc_01', 'EC_EBc_02', 'EC_EBc_03', 'EC_EBc_04', 'EC_EBc_05', 'EC_EBc_06', 'EC_EBc_07', 'EC_EBc_09', 'EC_EBc_11',
               'EC_LB_01', 'EC_LB_02', 'EC_LB_03', 'EC_LB_05', 'EC_LB_06', 'EC_LB_07', 'EC_LB_09', 'EC_LB_10',
               'EC_LBc_02', 'EC_LBc_03', 'EC_LBc_06', 'EC_LBc_09', 'EC_LB_11', 'EC_LBc_13','EC_LBc_14','EC_LBc_16', 'EC_LBc_14', 'EC_LBc_19','EC_LBc_20']

if group is not None: # only load animals within a specific group. otherwise, load all animals
    animals = [a for a in animals if (group + '_') in a]


# animals = ['EC_EB_27', 'EC_EB_32', 'EC_EBc_02','EC_EBc_04', 'EC_LB_05','EC_LB_09', 'EC_LBc_14','EC_LB_16']
path = r'E:\dark_drift\data'
data_object = batchProcessing(animals, animal_dobs, path,
                              stim = load_stim,
                              tracked_cells = tracked_cells,
                              tracking_across_days = tracking_across_days, # if true,track across days. if false, track within a day
                              roi_detection = 'functional',
                              ntheta = 8,
                              fps = 20,
                              zscore_threshold=0.8,
                              deconvolved = True,
                              plot_fov = False)

trac = data_object.dat['EC_EBc_09'].dat_subject['20251210']['spon1']['zscored_traces'][np.array([4,7,9,11,17,23,26,34]), 800:2000]

n_traces, T = trac.shape

# --- sampling rate (CHANGE if needed)
fs = 20  # Hzz

# --- time axis in seconds
time = np.arange(T) / fs

# --- normalize each trace
trac_norm = trac / np.max(np.abs(trac), axis=1, keepdims=True)

# --- vertical spacing
spacing = 2.0

# --- colors
colors = plt.cm.plasma(np.linspace(0, 0.93, n_traces))

plt.figure(figsize=(5, 5))

for i in range(n_traces):
    offset = i * spacing
    plt.plot(time, trac_norm[i] + offset, color=colors[i], linewidth=1.5)

# aesthetics
plt.xlabel('Time (s)')
plt.yticks([])
plt.title('Neural Traces (normalized + offset)')

plt.tight_layout()
plt.show()


corr_path = r'E:\dark_drift\figures\correlation distribution'
eb_aligned, eb_shuffled = np.load(os.path.join(corr_path, 'corr dist all animals EB aligned.npy')), np.load(os.path.join(corr_path, 'corr dist all animals EB shuffled.npy'))
ebc_aligned, ebc_shuffled = np.load(os.path.join(corr_path, 'corr dist all animals EBc aligned.npy')), np.load(os.path.join(corr_path, 'corr dist all animals EBc shuffled.npy'))
lb_aligned, lb_shuffled = np.load(os.path.join(corr_path, 'corr dist all animals LB aligned.npy')), np.load(os.path.join(corr_path, 'corr dist all animals LB shuffled.npy'))
lbc_aligned, lbc_shuffled = np.load(os.path.join(corr_path, 'corr dist all animals LBc aligned.npy')), np.load(os.path.join(corr_path, 'corr dist all animals LBc shuffled.npy'))

plt.figure()
x = np.arange(eb_aligned[0, :].mean(axis=-1).shape[0])
plt.scatter(x, eb_aligned[0, :].mean(axis = -1), c = 'blue', label = 'EB')
plt.plot(x, eb_aligned[0, :].mean(axis = -1), c = 'blue')
plt.fill_between(x, eb_aligned[0, :].mean(axis = -1) - sem (eb_aligned[0, :], axis = -1), eb_aligned[0, :].mean(axis = -1) + sem (eb_aligned[0, :], axis = -1), color='blue', alpha=0.2)
plt.scatter(x, ebc_aligned[0, :].mean(axis = -1), c = 'black', label = 'EBc')
plt.plot(x, ebc_aligned[0, :].mean(axis = -1), c = 'black')
plt.fill_between(x, ebc_aligned[0, :].mean(axis = -1) - sem (ebc_aligned[0, :], axis = -1), ebc_aligned[0, :].mean(axis = -1) + sem (ebc_aligned[0, :], axis = -1), color='black', alpha=0.2)
plt.scatter(x, lb_aligned[0, :].mean(axis = -1), c = 'red', label = 'LB')
plt.plot(x, lb_aligned[0, :].mean(axis = -1), c = 'red')
plt.fill_between(x, lb_aligned[0, :].mean(axis = -1) - sem (lb_aligned[0, :], axis = -1), lb_aligned[0, :].mean(axis = -1) + sem (lb_aligned[0, :], axis = -1), color='red', alpha=0.2)
plt.scatter(x, lbc_aligned[0, :].mean(axis = -1), c = 'green', label = 'LBc')
plt.plot(x, lbc_aligned[0, :].mean(axis = -1), c = 'green')
plt.fill_between(x, lbc_aligned[0, :].mean(axis = -1) - sem (lbc_aligned[0, :], axis = -1), lbc_aligned[0, :].mean(axis = -1) + sem (lbc_aligned[0, :], axis = -1), color='green', alpha=0.2)
xtick_labels = ['W0-W0', 'W0-W1', 'W0-W2', 'W0-W3']
plt.xticks(x, xtick_labels)
plt.ylabel('Correlation to W0')
plt.legend()
plt.show()

plt.figure()

x = np.arange(np.diagonal(eb_aligned.mean(axis = -1), offset = 1).shape[0])

plt.scatter(x, np.diagonal(eb_aligned.mean(axis = -1), offset = 1), c = 'blue', label = 'EB')
plt.plot(x, np.diagonal(eb_aligned.mean(axis = -1), offset = 1), c = 'blue')
plt.fill_between(x, np.diagonal(eb_aligned.mean(axis = -1), offset = 1) - sem (np.diagonal(eb_aligned, offset = 1).T, axis = -1), np.diagonal(eb_aligned.mean(axis = -1), offset = 1) + sem (np.diagonal(eb_aligned, offset = 1).T, axis = -1), color='blue', alpha=0.2)

plt.scatter(x, np.diagonal(ebc_aligned.mean(axis = -1), offset = 1), c = 'black', label = 'EBc')
plt.plot(x, np.diagonal(ebc_aligned.mean(axis = -1), offset = 1), c = 'black')
plt.fill_between(x, np.diagonal(ebc_aligned.mean(axis = -1), offset = 1) - sem (np.diagonal(ebc_aligned, offset = 1).T, axis = -1), np.diagonal(ebc_aligned.mean(axis = -1), offset = 1) + sem (np.diagonal(ebc_aligned, offset = 1).T, axis = -1), color='black', alpha=0.2)

plt.scatter(x, np.diagonal(lb_aligned.mean(axis = -1), offset = 1), c = 'red', label = 'LB')
plt.plot(x, np.diagonal(lb_aligned.mean(axis = -1), offset = 1), c = 'red')
plt.fill_between(x, np.diagonal(lb_aligned.mean(axis = -1), offset = 1) - sem (np.diagonal(lb_aligned, offset = 1).T, axis = -1), np.diagonal(lb_aligned.mean(axis = -1), offset = 1) + sem (np.diagonal(lb_aligned, offset = 1).T, axis = -1), color='red', alpha=0.2)

plt.scatter(x, np.diagonal(lbc_aligned.mean(axis = -1), offset = 1), c = 'green', label = 'LBc')
plt.plot(x, np.diagonal(lbc_aligned.mean(axis = -1), offset = 1), c = 'green')
plt.fill_between(x, np.diagonal(lbc_aligned.mean(axis = -1), offset = 1) - sem (np.diagonal(lbc_aligned, offset = 1).T, axis = -1), np.diagonal(lbc_aligned.mean(axis = -1), offset = 1) + sem (np.diagonal(lbc_aligned, offset = 1).T, axis = -1), color='green', alpha=0.2)

xtick_labels = ['W0-W1', 'W1-W2', 'W2-W3']
plt.xticks(x, xtick_labels)
plt.ylabel('Correlation between weeks')
plt.legend()
plt.show()


groups = np.unique([a.split('_')[1] for a in data_object.list_animals])
plot_group_corr_matrix(data_object, 'EB', animals)
plot_group_corr_matrix(data_object, 'EBc', animals)

# animal = 'EC_EB_08'
# print("You're gay!")
if (data_object.tracked_cells) and (data_object.stim == 'grat'):
    for animal in data_object.list_animals:
        # fov_across_days(data_object, brightness = 0.5, contrast = 2.2)


        data_object.dat[animal].corr_matrix, data_object.dat[animal].null_corr_matrix  = plot_corr(data_object, animal, thresholded_cells=0,
                                                        across_days=False)  # includes within-day shuffle
        plot_response(data_object, animal)
        rasters_across_days(data_object, animal)
        polar_plots_across_days(data_object, animal)
        polar_plots_across_days_rois(data_object, animal)

        plot_corr_change(corr)

        plot_response(data_object, animal)
        # look at cells that are responsive at least once + cell that are responsive all days
        # add activity threshold


# ALL cells, gratings
elif (not data_object.tracked_cells) and (data_object.stim == 'grat'):

    for animal in data_object.dat.keys():
        if data_object.dat[animal].dat_subject:
            plot_ori_projection(data_object, animal, grat=True, n_components_plot=3, plot_time=False)
            plot_ori_projection(data_object, animal, grat=True, n_components_plot=3, plot_time=True)

    for metric in ['OSI', 'DSI', 'preferred_orientation', 'preferred_direction']: # use deconvolved traces
        hist_osi_angle_all_animals(data_object, metric)
        hist_osi_angle_all_animals_outlines(data_object, metric)
    # for recording_day in ['P70', 'P77', 'P84', 'P91','P105','P112','P119','P126']:
        # responsiveness(data_object, recording_day)
        # trial_by_trial_reliability(data_object)
    responsive_cells_across_days(data_object)
    plot_avg_response(data_object, 'EB')
    plot_avg_response(data_object, 'LB')
    plot_avg_response_time(data_object)

    reliability_x_osi(data_object)
    trial_by_trial_reliability_days(data_object)
    for animal in data_object.list_animals:
        # polar_plots(data_object, animal) #tuning curves of all rois
        #trial_by_trial_reliability(data_object, 'P70')


elif (data_object.stim == 'spon'):
    var_explained, slope_eigenvals = spontaneous_analysis_slope(data_object, rec_type='spon')
    # del slope_eigenvals['EB']['P70_s']
    # del slope_eigenvals['EBc']['P70_s']
    # del slope_eigenvals['LB']['P105_s']
    #del slope_eigenvals['LB']['P70_s']
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

import matplotlib.pyplot as plt
from imports import *
from helpers import *
#from preprocessing import *
from track2p_preprocessing import *
from cellreg_preprocessing import *
from behavioural_analysis import *
import pandas as pd
import os
import numpy as np


def load_dlc_pose_estimates(path_to_behav, likelihood_threshold=0.9):
    """
    Load DLC pose estimates from CSV using header names.

    Parameters
    ----------
    path_to_behav : str
        Path to behavioral folder containing DLC CSV file.
    likelihood_threshold : float
        DLC points with likelihood below this threshold are set to NaN.

    Returns
    -------
    pupil_coords : np.ndarray
        Shape (n_frames, 8), columns:
        [top_pupil_x, top_pupil_y,
         right_pupil_x, right_pupil_y,
         bottom_pupil_x, bottom_pupil_y,
         left_pupil_x, left_pupil_y]

    face_coords : np.ndarray
        Shape (n_frames, 26), columns:
        [top_eye_x, top_eye_y,
         right_eye_x, right_eye_y,
         bottom_eye_x, bottom_eye_y,
         left_eye_x, left_eye_y,
         nose_bridge_x, nose_bridge_y,
         top_nose_x, top_nose_y,
         tip_nose_x, tip_nose_y,
         bottom_nose_x, bottom_nose_y,
         mouth_opening_x, mouth_opening_y,
         chin_x, chin_y,
         whisker_top_x, whisker_top_y,
         whisker_left_x, whisker_left_y,
         whisker_right_x, whisker_right_y]
    """

    dlc_file_name = [file for file in os.listdir(path_to_behav) if file.endswith('.csv')][0]
    dlc_file_path = os.path.join(path_to_behav, dlc_file_name)

    df = pd.read_csv(dlc_file_path, header=[0, 1, 2], index_col=0)

    scorer = df.columns.get_level_values(0)[0]

    pupil_parts = [
        "top_pupil",
        "right_pupil",
        "bottom_pupil",
        "left_pupil",
    ]

    face_parts = [
        "top_eye",
        "right_eye",
        "bottom_eye",
        "left_eye",
        "nose_bridge",
        "top_nose",
        "tip_nose",
        "bottom_nose",
        "mouth_opening",
        "chin",
        "top_whisker",
        "left_whisker",
        "right_whisker",
    ]

    def extract_bodyparts(df, scorer, bodyparts, likelihood_threshold):
        coords_list = []

        for bp in bodyparts:
            try:
                x = df[(scorer, bp, "x")].copy()
                y = df[(scorer, bp, "y")].copy()
                likelihood = df[(scorer, bp, "likelihood")]

                x[likelihood < likelihood_threshold] = np.nan
                y[likelihood < likelihood_threshold] = np.nan

                coords_list.append(x.to_numpy())
                coords_list.append(y.to_numpy())

            except KeyError:
                raise KeyError(f"Bodypart '{bp}' not found in DLC CSV.")

        return np.column_stack(coords_list)

    pupil_coords = extract_bodyparts(df, scorer, pupil_parts, likelihood_threshold)
    face_coords = extract_bodyparts(df, scorer, face_parts, likelihood_threshold)

    return pupil_coords, face_coords

class track2pPreprocessing:

    def __init__(self, main_path, subject, tracked_cells, roi_detection, track2p_folder_name, single_plane = True, day = None):

        self.main_path = main_path
        self.subject = subject
        self.tracked_cells = tracked_cells
        self.roi_detection = roi_detection
        self.track2p_folder_name = track2p_folder_name
        self.single_plane = single_plane
        self.meanImg = []
        self.day = day
        #self.deconvolved = deconvolved      # whether or not we want to use deconvolved spikes (deconvolved = True) or raw fluorescence (deconvolved = False)
        self.load_track2p_output()
        if self.tracked_cells:
            self.track_cells()


    def load_track2p_output (self):
        '''
        load track2p output

        :param main_path: path to the data folder, for example : r'E:\DriftScape\Data'
        :param subject: animal name, for example 'EC_GECO_09'
        :param track2p_folder_name: folder that stores track2p output folder, for example 'track2p-3days'
        :param single_plane (bool): whether or not we have a single plane
        :generates:
            t2p_match_mat: matrix containing indices of matches neurons present across all sessions for a given plane
                # Neurons not successfully tracked on all days will contain None values in the matrix
            track_ops_dict
            track_ops
        '''

        # directory containing the 'track2p' output folder
        if self.day is None: # if we don't specify a day, we are loading the t2p output for the entire period
            t2p_save_path = os.path.join(self.main_path, self.subject)
        else: # if we do specify a day, navigate into that day's folder and take the t2p output on that day
            t2p_save_path = os.path.join(self.main_path, self.subject, self.day)

        # the plane to process
        if self.single_plane:
            self.plane = 'plane0'
            # for multiple planes, just loop through them (will need to reconfigure some code)

        # load plane0_match_mat.npy (the match matrix) > matrix with indices of matches neurons present across all sessions for a given plane
        # Neurons not successfully tracked on all days will contain None values in the matrix
        self.t2p_match_mat = np.load(os.path.join(t2p_save_path, self.track2p_folder_name, f'{self.plane}_match_mat.npy'), allow_pickle=True)

        # load track_ops.npy (contains settings/suite2p paths)
        self.track_ops_dict = np.load(os.path.join(t2p_save_path, self.track2p_folder_name, 'track_ops.npy'), allow_pickle=True).item()
        self.track_ops = SimpleNamespace(**self.track_ops_dict) # create dummy object from the track_ops dictionary

        # To get neurons tracked across all days > only take the rows of the matrices containing no None values in t2p_match_mat matrix
        self.t2p_match_mat_allday = self.t2p_match_mat[~np.any(self.t2p_match_mat == None, axis=1), :]
        print(f'There are {self.t2p_match_mat_allday.shape[0]} cells tracked across {self.t2p_match_mat_allday.shape[1]} days')

        # track_ops.npy ('settings file') contains all the paths to suite2p folders used for track2p algorithm
        # print paths to datasets used
        print('Datasets used for t2p:')


        for ds_path in self.track_ops.all_ds_path:
            print('\t', ds_path)


    def track_cells(self):
        '''
        Find cells that are tracked & present in all recordings ('matched cells')

        Load the activity of matched cells by:
            - looping through the datasets, loading suite2p files (ops.npy, stat.npy, iscell.npy, and F.npy)
            - filtering stat.npy and fluo.npy by the track2p iscell threshold (classical suite2p)
            - filtering stat.npy and fluo.npy by the appropriate indices from the matrix of neurons matched on all days

        Output: data structure where the indices of cells are matched within the stat and fluo objects
        '''

        iscell_thr = self.track_ops.iscell_thr  # use the same threshold as when running the algo (to be consistent with indexing)

        #self.['meanImg'] = []
        self.all_stat_t2p = []
        self.all_f_t2p_spikes = [] # deconvolved traces, spikes
        self.all_f_t2p_traces = [] # fluorescence
        self.all_ops = []  # ops dont change

        for (i, ds_path) in enumerate(self.track_ops.all_ds_path):

            # look in the G drive instead
            # ds_path = ds_path.replace('I:', 'E:')

            # i changed animal names. unless i rerun track2p, it doesnt update in the file
            old_animal_name = ds_path.split('\\')[3]
            ds_path = ds_path.replace(old_animal_name,self.subject)

            ds_path = Path(ds_path)

            # if not ds_path.exists():
            #     if ds_path.drive == 'G:':
            #         ds_path = Path(str(ds_path).replace('G:', 'E:', 1))
            #     elif ds_path.drive == 'E:':
            #         ds_path = Path(str(ds_path).replace('E:', 'G:', 1))

            ops = np.load(os.path.join(ds_path, f'suite2p {self.roi_detection}', self.plane, 'ops.npy'), allow_pickle=True).item()    #options and intermediate outputs (dictionary)
            stat = np.load(os.path.join(ds_path, f'suite2p {self.roi_detection}', self.plane, 'stat.npy'), allow_pickle=True)         # list of statistics computed for each cell (ROIs by 1)
            f = np.load(os.path.join(ds_path, f'suite2p {self.roi_detection}', self.plane, 'F.npy'), allow_pickle=True)               # array of fluorescence traces (ROIs by timepoints)
            fneu = np.load(os.path.join(ds_path, f'suite2p {self.roi_detection}', self.plane, 'Fneu.npy'), allow_pickle=True)         # array of neuropil fluorescence traces (ROIs by timepoints)
            spikes = np.load(os.path.join(ds_path, f'suite2p {self.roi_detection}', self.plane, 'spks.npy'), allow_pickle=True)       # array of deconvolved traces (ROIs by timepoints)
            iscell = np.load(os.path.join(ds_path, f'suite2p {self.roi_detection}', self.plane, 'iscell.npy'), allow_pickle=True)     # specifies whether an ROI is a cell, first column is 0/1, and second column is probability that the ROI is a cell based on the default classifier

            self.meanImg.append(ops['meanImg'])

            # if self.deconvolved: # use spikes
            #load spikes (deconvolved)
            if self.track_ops.iscell_thr is None:
                stat_iscell = stat[iscell[:, 0] == 1]
                f_iscell_spikes = spikes[iscell[:, 0] == 1, :]
            else:
                stat_iscell = stat[iscell[:, 1] > iscell_thr]
                f_iscell_spikes = spikes[iscell[:, 1] > iscell_thr, :]

            # load raw fluorescence traces
            # when using fluorescence (f) > need to subtract neuropil (F - 0.7*Fneu) > to get real fluorescence signal
            f -= 0.7 * fneu
            if self.track_ops.iscell_thr == None:
                stat_iscell = stat[iscell[:, 0] == 1]
                f_iscell_traces = f[iscell[:, 0] == 1, :]

            else:
                stat_iscell = stat[iscell[:, 1] > iscell_thr]
                f_iscell_traces = f[iscell[:, 1] > iscell_thr, :]

            stat_t2p = stat_iscell[self.t2p_match_mat_allday[:, i].astype(int)]
            f_t2p_traces = f_iscell_traces[self.t2p_match_mat_allday[:, i].astype(int), :]
            f_t2p_spikes = f_iscell_spikes[self.t2p_match_mat_allday[:, i].astype(int), :]

            self.all_stat_t2p.append(stat_t2p)
            self.all_f_t2p_spikes.append(f_t2p_spikes)
            self.all_f_t2p_traces.append(f_t2p_traces)
            self.all_ops.append(ops)

class suite2pPreprocessing:

    def __init__ (self, animal, stim, path, track2p_folder, ntheta = 8, fps = 10, track_cells_within_day = True, roi_detection = 'functional', deconvolved = False, zscore_threshold = None, std_threshold = None):

        self.animal = animal
        self.save_path = str(Path(path).parents[0] / 'figures')
        #self.path = os.path.join(path, animal, track2p_folder)
        self.track2p_folder = track2p_folder
        self.path = path
        self.deconvolved = deconvolved              # True = load deconvolved traces (spikes) ; False = load raw fluorescence traces
        self.track_cells_within_day = track_cells_within_day
        self.roi_detection = roi_detection
        self.fps = fps
        self.ntheta = ntheta
        self.stim = stim
        self.zscore_threshold = zscore_threshold
        self.std_threshold = std_threshold
        self.tracked_cells = self.track_cells_within_day

        if self.track_cells_within_day:

            self.track2p_obj = track2pPreprocessing(self.path, self.animal, True, self.roi_detection,
                                                    self.track2p_folder)  # load track2p folder
            self.paths = self.track2p_obj.track_ops.all_ds_path  # list of paths to data

            # take all the days in the file, but skip the first two weeks, which have no tracking data >  list of days (strings)
            self.days = sorted([x for x in os.listdir(os.path.join(path, animal)) if x.isdigit()])
            self.dat_subject = {str(day): {} for day in self.days}  # dictionary of day (key) and empty dictionary (value) pairs
            self.load_dat()

    def load_dat(self):

        for i_day, day in enumerate(self.days):

            day = str(day)

            if self.track_cells_within_day:

                for i_recording, recording_path in enumerate(self.paths):

                    recording = recording_path.split(os.sep)[5]
                    self.dat_subject[day][recording] = {}

                    suite2p_path = os.path.join(recording_path, f'suite2p {self.roi_detection}','plane0')
                    responses, iscell, ops, stat = suite2p_files(suite2p_path, response_type='fluorescence')
                    spikes, _, _, _ = suite2p_files(suite2p_path, response_type='deconvolved')
                    self.meanimg = load_ops_meanimg(ops, 'meanImg')

                    self.dat_subject[day][recording]['traces'] = self.track2p_obj.all_f_t2p_traces[i_recording]
                    self.dat_subject[day][recording]['deconvolved'] = self.track2p_obj.all_f_t2p_spikes[i_recording]

                    self.dat_subject[day][recording]['zscored_traces'] = zscore(self.dat_subject[day][recording]['traces'], axis=1)
                    self.dat_subject[day][recording]['zscored_deconvolved'] = zscore(self.dat_subject[day][recording]['deconvolved'], axis=1)

                    self.dat_subject[day][recording]['log'] = {}
                    logpath = os.path.join(os.path.dirname(recording_path), 'logfiles')
                    logfile = [file for file in os.listdir(logpath) if file.endswith('_log.txt')][0]
                    path_to_log_file = os.path.join(logpath, logfile)
                    self.dat_subject[day][recording]['log']['stim_dict'], self.dat_subject[day][recording]['log']['list_stim_types'], self.dat_subject[day][recording]['log']['list_stim_names'] = load_stims (self, path_to_log_file)

                    # store location of each cells, in array of shape (n_cells, xPos, yPos, z_pos)
                    # since cells are in same pos for all recordings (because of track2p), just look at the first day
                    if i_recording == 0:
                        z_coord = int(recording_path.split(os.sep)[-2].split('_')[1]) # assuming the z position stored in the filename
                        self.dat_subject[day]['cell_xyz_pos'] = np.zeros((self.track2p_obj.all_stat_t2p[0].shape[0], 3))
                        for cell in np.arange(self.track2p_obj.all_stat_t2p[0].shape[0]):
                            self.dat_subject[day]['cell_xyz_pos'][cell] = np.array(
                                [np.round(self.track2p_obj.all_stat_t2p[0][cell]['xpix'].mean(), 2),
                                 np.round(self.track2p_obj.all_stat_t2p[0][cell]['ypix'].mean(), 2),
                                 z_coord])

                    # load ttl / neural information
                    ttl_path = os.path.join(recording_path, f'{recording_path.split(os.sep)[-2]}.mat')
                    event_id = np.squeeze(loadmat(ttl_path)['info']['event_id'][0][0])
                    self.dat_subject[day][recording]['ttl_data'] = np.squeeze(loadmat(ttl_path)['info']['frame'][0][0])#[event_id == 1]

                    # p = Path (recording_path)
                    behav_path = os.path.join(Path (recording_path).parent, 'behav')
                    self.dat_subject[day][recording]['pupil_pos'], self.dat_subject[day][recording]['face_pos'] = load_dlc_pose_estimates(behav_path, 0.7)


                    # print(self.dat_subject[day][recording]['pupil_pos'].shape, self.dat_subject[day][recording]['zscored_traces'].shape)
                    x_center, y_center, width, height, angle = fit_ellipse(self.dat_subject[day][recording]['pupil_pos'])

                    # plt.figure()
                    # t = 0
                    # plt.plot()
                    #
                    # fig, ax = plt.subplots()
                    # # plot_ellipse(ax, x_center, y_center, width, height, angle, points=None, color='red')
                    #
                    # [plot_ellipse(ax,
                    #              float(x_center[t]), float(y_center[t]),
                    #              float(width[t]), float(height[t]),
                    #              float(angle[t]),
                    #              color='red') for t in np.arange(self.dat_subject[day][recording]['pupil_pos'].shape[0])]
                    #
                    # plt.show()

                    # fig, ax = plt.subplots()
                    # for t in np.arange(300):
                    # # for t in np.arange(self.dat_subject[day][recording]['pupil_pos'].shape[0]):
                    #     # ellipse
                    #     ell = Ellipse(
                    #         (float(x_center[t]), float(y_center[t])),
                    #         float(width[t]), float(height[t]),
                    #         angle=float(np.degrees(angle[t])),
                    #         edgecolor="red", facecolor="none", linewidth=2
                    #     )
                    #     ax.add_patch(ell)
                    #
                    #     # plot the 4 points for that frame (so you SEE something)
                    #     data = self.dat_subject[day][recording]['pupil_pos']
                    #     pts = data[t].reshape(4, 2)  # [[topx,topy],[rightx,righty],[bottomx,bottomy],[leftx,lefty]]
                    #     ax.scatter(pts[:, 0], pts[:, 1], s=20)
                    #
                    #     # ax.set_xlim(300,1000)
                    #     # ax.set_ylim(300,1000)
                    #
                    #     ax.set_aspect("equal", adjustable="box")
                    #     plt.gca().invert_yaxis()  # optional if these are image coordinates (0,0 top-left)
                    # plt.show()



                    # #load ttl data, correct ttl data, loads dictionary that holds neural data for ttls values
                    get_behav_response_ttls(self, day, recording)

                    # min_n_frames = np.min([self.dat_subject[day][recording]['traces'].shape[1], self.dat_subject[day][recording]['pupil_pos'].shape[0], self.dat_subject[day][recording]['face_pos'].shape[0]])
                    # print(min_n_frames)

                    # self.dat_subject[day][recording]['eye_matrix'] = self.dat_subject[day][recording]['pupil_pos']
                    # self.dat_subject[day][recording]['face_matrix'] = self.dat_subject[day][recording]['face_pos']

                    del self.dat_subject[day][recording]['responses_ttls_whole']
                    del self.dat_subject[day][recording]['zscored_responses_ttls_whole']

                    # rename the keys in the dictionary to be the names of the stims
                    self.dat_subject[day][recording]['responses_ttls'] = {
                        self.dat_subject[day][recording]['dict_stim_names'][k]: v for k, v in self.dat_subject[day][recording]['responses_ttls'].items()}
                    self.dat_subject[day][recording]['zscored_responses_ttls'] = {
                        self.dat_subject[day][recording]['dict_stim_names'][k]: v for k, v in self.dat_subject[day][recording]['zscored_responses_ttls'].items()}

                    # self.dat_subject[day][recording]['dict_stim_names']
                    #build_parameter_matrix(self, day, response_window='whole', zscore=True)
                    # store_metrics(self, day, recording, response_window = 'whole', zscore = True)

                    # responses_ttls houses the response to each image, with the stim names as values

                    # the neural and behavioral data is originally in a giant matrix, and we separate it out according to which stimulus is shown at that time
                    # we then just want to keep the data (without dictionary keys) but retain the sequential information
                    # so use list comprehension to lump everything together into one giant matrix of shape
                    # neural dat: (n_stims, n_cells, n_timepoints)
                    # eye / face matrix (n_stims, n_timepoints, n_keypoints)
                    # matrix just lumps everything together > we are grouping the neural/behav data accor
                    self.dat_subject[day][recording]['matrix'] = np.stack(
                        [arr.squeeze(0) for k, arr in self.dat_subject[day][recording]['responses_ttls'].items() if
                         'Wait' not in k], axis=0)
                    self.dat_subject[day][recording]['zscored_matrix'] = np.stack(
                        [arr.squeeze(0) for k, arr in self.dat_subject[day][recording]['zscored_responses_ttls'].items() if
                         'Wait' not in k], axis=0)
                    self.dat_subject[day][recording]['eye_matrix'] = np.stack(
                        [arr.squeeze(0) for k, arr in self.dat_subject[day][recording]['pupil_keypoints'].items() if
                         'Wait' not in k], axis=0)
                    self.dat_subject[day][recording]['face_matrix'] = np.stack(
                        [arr.squeeze(0) for k, arr in self.dat_subject[day][recording]['face_keypoints'].items() if
                         'Wait' not in k], axis=0)

        # combining data from all the chunks into one array
        matrices = [
            self.dat_subject[day][recording]['matrix']
            for day in self.days
            for recording in self.dat_subject[day] if 'chunk' in recording]

        self.matrix = np.concatenate(matrices, axis=0)

        eye_behav_matrices = [
            self.dat_subject[day][recording]['eye_matrix']
            for day in self.days
            for recording in self.dat_subject[day] if 'chunk' in recording]

        self.eye_behav_matrices = np.concatenate(eye_behav_matrices, axis=0)

        face_behav_matrices = [
            self.dat_subject[day][recording]['face_matrix']
            for day in self.days
            for recording in self.dat_subject[day] if 'chunk' in recording]

        self.face_behav_matrices = np.concatenate(face_behav_matrices, axis=0)

        self.stims = [stim
            for day in self.days
            for recording in self.dat_subject[day] if 'chunk' in recording
            for stim in self.dat_subject[day][recording]['log']['list_stim_names']
            if 'Wait' not in stim]




path = r'I:\sensorium\data'
animal = 'EC_EB_40' # EC_EBc_06 was the pilot control
day = '20260413'

data_object = suite2pPreprocessing(animal,
                              stim = 'chunk',
                              path = path,
                              track2p_folder= os.path.join(path, animal, day, 'track2p-chunk-functional'),
                              track_cells_within_day = True,
                              roi_detection = 'functional',
                              ntheta = 8,
                              fps = 20,
                              zscore_threshold=0.8,
                              deconvolved = True)

np.save(fr"{path}\{animal}\{day}\dat_matrices\pupil_keypoints.npy", data_object.eye_behav_matrices)
np.save(fr"{path}\{animal}\{day}\dat_matrices\face_keypoints.npy", data_object.face_behav_matrices)
np.save(fr"{path}\{animal}\{day}\dat_matrices\responses.npy", data_object.matrix)
np.save(fr"{path}\{animal}\{day}\dat_matrices\stims.npy", data_object.stims)
np.save(fr"{path}\{animal}\{day}\dat_matrices\cell_locations.npy", data_object.dat_subject[day]['cell_xyz_pos'])

# pupil of shape n_trials, n_frames, n_keypoints:
        # [top_pupil_x, top_pupil_y,
        #  right_pupil_x, right_pupil_y,
        #  bottom_pupil_x, bottom_pupil_y,
        #  left_pupil_x, left_pupil_y]
pupil = np.load(fr"{path}\{animal}\{day}\dat_matrices\pupil_keypoints.npy", allow_pickle=True)

# face of shape n_trials, n_frames, n_keypoints:
        # [top_eye_x, top_eye_y,
        #  right_eye_x, right_eye_y,
        #  bottom_eye_x, bottom_eye_y,
        #  left_eye_x, left_eye_y,
        #  nose_bridge_x, nose_bridge_y,
        #  top_nose_x, top_nose_y,
        #  tip_nose_x, tip_nose_y,
        #  bottom_nose_x, bottom_nose_y,
        #  mouth_opening_x, mouth_opening_y,
        #  chin_x, chin_y,
        #  whisker_top_x, whisker_top_y,
        #  whisker_left_x, whisker_left_y,
        #  whisker_right_x, whisker_right_y]
face = np.load(fr"{path}\{animal}\{day}\dat_matrices\face_keypoints.npy", allow_pickle=True)

# responses of shape n_trials, n_cells, n_frames
responses = np.load(fr"{path}\{animal}\{day}\dat_matrices\responses.npy", allow_pickle=True)

# stims of shape n_trials
stims = np.load(fr"{path}\{animal}\{day}\dat_matrices\stims.npy", allow_pickle=True)

#cell_locations of shape n_cells, 3 [x,y,z]
cell_locations = np.load(fr"{path}\{animal}\{day}\dat_matrices\cell_locations.npy", allow_pickle=True)

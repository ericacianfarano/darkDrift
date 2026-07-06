# demo from https://github.com/juremaj/track2p/blob/main/demo_t2p_ouputs.ipynb
from imports import *

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
            ds_path = ds_path.replace('I:', 'E:')

            # i changed animal names. unless i rerun track2p, it doesnt update in the file
            old_animal_name = ds_path.split('\\')[3]
            ds_path = ds_path.replace(old_animal_name,self.subject)

            ds_path = Path(ds_path)

            if not ds_path.exists():
                if ds_path.drive == 'G:':
                    ds_path = Path(str(ds_path).replace('G:', 'E:', 1))
                elif ds_path.drive == 'E:':
                    ds_path = Path(str(ds_path).replace('E:', 'G:', 1))

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

############################################################################################################################
# VISUALIZE ROIs & ACTIVITY OF MATCHED CELLS
# extract ROI information from all_stat
    # index by the day to get stat_t2p from all_stat2p (this is the sorted stat object for that day).
    # get the ROI info by indexing stat_t2p by the index of the cell match (because of resorting we use the same index across days).

def plot_rois_across_days (data_object, wind_value, n_cells_to_plot = 8):
    '''
    Randomly sample n_cells_to_plot ROIs
    Plot ROI of cells across all recordings days

    :param wind_value: how much of each surround we wish to see
        - larger wind_value provides more zoomed out view of ROI
        - smaller wind_value provides more zoomed in view of ROI
    :param n_cells_to_plot: number of cells to randomly choose for plotting
    '''

    # randomly choose n_cells_to_plot cell ROIs to plot across days
    n_tracked_cells = data_object.t2p_match_mat_allday.shape[0]
    cells_to_plot = np.random.randint(0, n_tracked_cells, n_cells_to_plot)

    fig, ax = plt.subplots(n_cells_to_plot, len(data_object.track_ops.all_ds_path), figsize = (5,8))
    for i, cell_idx in enumerate(cells_to_plot):
        for i_day, path in enumerate(data_object.track_ops.all_ds_path):
            mean_img = data_object.all_ops[i_day]['meanImg']
            stat_t2p = data_object.all_stat_t2p[i_day]
            median_coord = stat_t2p[cell_idx]['med']

            # plot a short window around the ROI centroid
            ax[i, i_day].imshow(mean_img[int(median_coord[0])-wind_value:int(median_coord[0])+wind_value, int(median_coord[1])-wind_value:int(median_coord[1])+wind_value], cmap='gray')
            ax[i, i_day].scatter(wind_value, wind_value)
            ax[i, i_day].axis('off')

            if i ==0:
                ax[i, i_day].set_title(path.split('\\')[4])

    plt.suptitle(path.split('\\')[3])
    fig.tight_layout()
    plt.show()

def plot_cell_activity_days (data_object):
    '''
    Plot the activity/trace/fluorescence of cell 'cell_idx' across all days

    :param cell_idx:
    '''
    # randomly choose ROI to plot
    n_tracked_cells = data_object.t2p_match_mat_allday.shape[0]
    cell_idx = np.random.randint(0, n_tracked_cells)

    fig, ax = plt.subplots( len(data_object.track_ops.all_ds_path),1, figsize = (10,6))
    for i_day, path in enumerate(data_object.track_ops.all_ds_path):
        ax[i_day].plot(data_object.all_f_t2p[i_day][cell_idx, :])
        ax[i_day].set_xlabel('Frame')
        ax[i_day].set_ylabel('F')
        ax[i_day].set_title(path.split('\\')[4])
    plt.suptitle(f'ROI # {cell_idx}')
    fig.tight_layout()
    plt.show()

def rasters_across_days (data_object):
    '''
    Visualize the raster plots
        - data is already sorted such that the rows represent the same cell across days
        - so we just need to loop through all_f_t2p and plot each element
    '''
    fig, ax = plt.subplots(len(data_object.track_ops.all_ds_path),1, figsize = (10,6))
    for i_day, path in enumerate(data_object.track_ops.all_ds_path):
        f_plot = zscore(data_object.all_f_t2p[i_day], axis=1)
        ax[i_day].imshow(f_plot, aspect='auto', cmap='Greys', vmin=0, vmax=1.96)
        ax[i_day].set_xlabel('Frame')
        ax[i_day].set_ylabel('ROI')
        ax[i_day].set_title(path.split('\\')[4])
    plt.suptitle(f'Raster plots - all ROIs tracked across days')
    fig.tight_layout()
    plt.show()

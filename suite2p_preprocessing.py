from imports import *
from helpers import *
#from preprocessing import *
from track2p_preprocessing import *

class suite2pPreprocessing:

    def __init__ (self, animal, path, track2p_folder, ntheta = 8, grating_repeats = 2, grating_blocks = 4, fps = 10, deconvolved = False):

        self.animal = animal
        self.save_path = str(Path(path).parents[0] / 'figures')
        self.path = os.path.join(path, animal, track2p_folder)
        self.deconvolved = deconvolved              # True = load deconvolved traces (spikes) ; False = load raw fluorescence traces
        self.fps = fps
        self.ntheta = ntheta
        self.grating_repeats = grating_repeats
        self.grating_blocks = grating_blocks
        self.track2p_obj = track2pPreprocessing(path, animal, track2p_folder, deconvolved=self.deconvolved)  # load track2p folder

        self.paths = self.track2p_obj.track_ops.all_ds_path               # list of paths to data

        # look in the G drive, instead of the E drive
        self.paths = [p.replace('E:', 'H:') for p in self.paths]

        self.days = [path_str.split('\\')[4] for path_str in self.paths]        # list of days (strings)
        self.dat_subject = {day: {} for day in self.days}    # dictionary of day (key) and empty dictionary (value) pairs
        self.load_dat(self.track2p_obj)
        #self.analysis()

    def load_dat(self, track2p_object):

        for i_day, day in enumerate(self.days):

            # load suite2p roi data processed by track2p
            self.dat_subject[day]['tracked_fluorescence'] = track2p_object.all_f_t2p[i_day]
            self.dat_subject[day]['zscored_tracked_fluorescence'] = zscore(self.dat_subject[day]['tracked_fluorescence'], axis=1)

            # load log file information
            self.dat_subject[day]['log'] = {}
            logpath = os.path.join(os.path.dirname(self.paths[i_day]), 'logfiles')

            logpath = logpath.replace('E:', 'G:')

            logfile = [file for file in os.listdir(logpath) if file.endswith('_log.txt')][0]
            path_to_log_file = os.path.join(logpath, logfile)

            # not sure what to do with this yet
            self.dat_subject[day]['log']['stim_dict'], self.dat_subject[day]['log']['list_stim_types'], self.dat_subject[day]['log']['list_stim_names'] = load_stims (self, path_to_log_file)

            # load behaviour video (deal with this later)
            # behav_path = os.path.join(os.path.dirname(self.paths[i_day]), 'behav_files', f'{self.paths[i_day].split('\\')[-2]}_eye.mj2')
            # want to read through it, bin, and downsample it, and rewrite it to use in facemap (https://github.com/MouseLand/facemap)

            # load ttl / neural information
            path_char = '\\'
            ttl_path = os.path.join(self.paths[i_day], f'{self.paths[i_day].split(path_char)[-2]}.mat')
            event_id = np.squeeze(loadmat(ttl_path)['info']['event_id'][0][0])
            self.dat_subject[day]['ttl_data'] = np.squeeze(loadmat(ttl_path)['info']['frame'][0][0])#[event_id == 1]

            # WANT TO Z SCORE ACCORDING TO BASELINE HERE

            parsed = np.array([parse_grating_array(
                self.dat_subject[day]['log']['stim_dict'][grating_block]) for grating_block in
                self.dat_subject[day]['log']['stim_dict'].keys() if 'Grating' in grating_block])

            self.dat_subject[day]['thetas'] = parsed[:, :, 0]
            self.dat_subject[day]['n_repeats'] = self.dat_subject[day]['thetas'].shape[0]
            self.dat_subject[day]['n_thetas'] = self.dat_subject[day]['thetas'].shape[1]

            #load ttl data, correct ttl data, loads dictionary that holds neural data for ttls values
            get_neuronal_responses_ttls(self, day)

            #build_parameter_matrix(self, day, response_window='whole', zscore=True)
            #store_metrics(self, day, response_window='moving', zscore=True)

            store_metrics(self, day, response_window = 'whole', zscore = True)

            # calculating orientation selectivity > average responses to the same orientation within the same 'block'
            #responses_gratings(self, day, n_grating_repeats=self.grating_repeats)

            #self.ttl_d, self.ttls, self.stimulus_responses = stim_responses(self.fps, self.ttl_data, self.stim_dict, self.dat_subject[day]['tracked_fluorescence'].T)

    def analysis (self):
        get_complex_num(self)

# animal = 'EC_GCaMP6s_06'
# path = r'I:\dark_drift_trial\data'
# track2p_folder = 'track2p'
# suite2p_obj = suite2pPreprocessing(animal, path, track2p_folder, ntheta = 8,  grating_repeats = 6, fps = 20, deconvolved = True)

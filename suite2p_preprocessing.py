from imports import *
from helpers import *
#from preprocessing import *
from track2p_preprocessing import *

class batchProcessing:
    def __init__(self, list_animals, animal_dobs, path, ntheta = 8, fps = 20,  stim = 'grat', tracked_cells = False, deconvolved = True):
        self.list_animals = list_animals
        self.animal_dobs = animal_dobs
        self.path = path
        self.ntheta = ntheta
        self.fps = fps
        self.stim = stim # either 'grat' or 'spon'
        self.tracked_cells = tracked_cells
        self.deconvolved = deconvolved
        self.dat = {animal: suite2pPreprocessing(animal, self.path, os.path.join(self.path, animal, f'track2p-{self.stim}', 'track2p'), ntheta = self.ntheta, fps = self.fps, tracked_cells=self.tracked_cells, deconvolved = self.deconvolved) for animal in self.list_animals}


class suite2pPreprocessing:

    def __init__ (self, animal, path, track2p_folder, ntheta = 8, fps = 10, tracked_cells = False, deconvolved = False):

        self.animal = animal
        self.save_path = str(Path(path).parents[0] / 'figures')
        #self.path = os.path.join(path, animal, track2p_folder)
        self.deconvolved = deconvolved              # True = load deconvolved traces (spikes) ; False = load raw fluorescence traces
        self.tracked_cells = tracked_cells
        self.fps = fps
        self.ntheta = ntheta
        self.track2p_obj = track2pPreprocessing(path, animal, self.tracked_cells, track2p_folder)  # load track2p folder
        self.paths = self.track2p_obj.track_ops.all_ds_path               # list of paths to data

        # look in the G drive, instead of the E drive
        #self.paths = [p.replace('E:', 'H:') for p in self.paths]

        self.days = np.unique([path_str.split('\\')[4] for path_str in self.paths])        # list of days (strings)
        self.dat_subject = {day: {} for day in self.days}    # dictionary of day (key) and empty dictionary (value) pairs
        self.load_dat(self.track2p_obj)
        #self.analysis()

    def load_dat(self, track2p_object):

        for i_day, day in enumerate(self.days):
            for recording_path in [p for p in self.paths if day in p]: #self.paths[i_day].split(os.sep)[5]:

                recording = recording_path.split(os.sep)[5]
                self.dat_subject[day][recording] = {}

                # load suite2p information
                # load both deconvolved and fluorescence traces
                suite2p_path = os.path.join(recording_path, 'suite2p','plane0')
                responses, iscell, ops, stat = suite2p_files(suite2p_path, response_type='fluorescence')
                spikes, _, _, _ = suite2p_files(suite2p_path, response_type='deconvolved')

                # store suite2p information (all cells)
                self.dat_subject[day][recording]['meanImg'] = ops['meanImg']

                if self.tracked_cells: # load suite2p roi data processed by track2p
                    self.dat_subject[day][recording]['traces'] = track2p_object.all_f_t2p_traces[i_day]
                    self.dat_subject[day][recording]['deconvolved'] = track2p_object.all_f_t2p_spikes[i_day]

                else: # load suite2p data, for all cells
                    self.dat_subject[day][recording]['traces'] = responses[iscell == 1]
                    self.dat_subject[day][recording]['deconvolved'] = spikes[iscell == 1]

                self.dat_subject[day][recording]['zscored_traces'] = zscore(self.dat_subject[day][recording]['traces'], axis=1)
                self.dat_subject[day][recording]['zscored_deconvolved'] = zscore(self.dat_subject[day][recording]['deconvolved'], axis=1)

                # load behaviour video (deal with this later)
                # behav_path = os.path.join(os.path.dirname(self.paths[i_day]), 'behav_files', f'{self.paths[i_day].split('\\')[-2]}_eye.mj2')
                # want to read through it, bin, and downsample it, and rewrite it to use in facemap (https://github.com/MouseLand/facemap)

                # only load log/ttl file if recording is 'grat'
                if 'spon' not in recording:

                    # load log file information
                    self.dat_subject[day][recording]['log'] = {}
                    logpath = os.path.join(os.path.dirname(recording_path), 'logfiles')

                    logpath = logpath.replace('E:', 'G:')

                    logfile = [file for file in os.listdir(logpath) if file.endswith('_log.txt')][0]
                    path_to_log_file = os.path.join(logpath, logfile)

                    # not sure what to do with this yet
                    self.dat_subject[day][recording]['log']['stim_dict'], self.dat_subject[day][recording]['log']['list_stim_types'], self.dat_subject[day][recording]['log']['list_stim_names'] = load_stims (self, path_to_log_file)

                    # load ttl / neural information
                    ttl_path = os.path.join(recording_path, f'{recording_path.split(os.sep)[-2]}.mat')
                    event_id = np.squeeze(loadmat(ttl_path)['info']['event_id'][0][0])
                    self.dat_subject[day][recording]['ttl_data'] = np.squeeze(loadmat(ttl_path)['info']['frame'][0][0])#[event_id == 1]

                    # WANT TO Z SCORE ACCORDING TO BASELINE HERE

                    parsed = np.array([parse_grating_array(
                        self.dat_subject[day][recording]['log']['stim_dict'][grating_block]) for grating_block in
                        self.dat_subject[day][recording]['log']['stim_dict'].keys() if 'Grating' in grating_block])

                    self.dat_subject[day][recording]['thetas'] = parsed[:, :, 0]
                    self.dat_subject[day][recording]['n_repeats'] = self.dat_subject[day][recording]['thetas'].shape[0]
                    self.dat_subject[day][recording]['n_thetas'] = self.dat_subject[day][recording]['thetas'].shape[1]

                    #load ttl data, correct ttl data, loads dictionary that holds neural data for ttls values
                    get_neuronal_responses_ttls(self, day, recording)

                    #build_parameter_matrix(self, day, response_window='whole', zscore=True)
                    #store_metrics(self, day, response_window='moving', zscore=True)

                    store_metrics(self, day, recording, response_window = 'whole', zscore = True)

                    # calculating orientation selectivity > average responses to the same orientation within the same 'block'
                    #responses_gratings(self, day, n_grating_repeats=self.grating_repeats)

                    #self.ttl_d, self.ttls, self.stimulus_responses = stim_responses(self.fps, self.ttl_data, self.stim_dict, self.dat_subject[day]['tracked_fluorescence'].T)

    def analysis (self):
        get_complex_num(self)

# animal = 'EC_GCaMP6s_06'
# path = r'I:\dark_drift_trial\data'
# track2p_folder = 'track2p'
# suite2p_obj = suite2pPreprocessing(animal, path, track2p_folder, ntheta = 8,  grating_repeats = 6, fps = 20, deconvolved = True)

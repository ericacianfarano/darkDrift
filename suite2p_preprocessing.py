import matplotlib.pyplot as plt

from imports import *
from helpers import *
#from preprocessing import *
from track2p_preprocessing import *
from cellreg_preprocessing import *

class batchProcessing:
    def __init__(self, list_animals, animal_dobs, path, ntheta = 8, fps = 20,  stim = 'grat', tracked_cells = False, tracking_across_days = True, checked_tracked_cells = None, roi_detection = 'functional', deconvolved = True, plot_fov = False, zscore_threshold = None, std_threshold = None):
        self.list_animals = list_animals
        self.animal_dobs = animal_dobs
        self.path = path
        self.ntheta = ntheta
        self.fps = fps
        self.stim = stim # either 'grat' or 'spon'
        self.tracked_cells = tracked_cells
        self.tracking_across_days = tracking_across_days
        self.checked_tracked_cells = checked_tracked_cells
        self.roi_detection = roi_detection
        self.deconvolved = deconvolved
        self.plot_fov = plot_fov
        self.zscore_threshold = zscore_threshold
        self.std_threshold = std_threshold
        # if stim is None and (not tracking_across_days): # want to load track2p for a specifiec day
        #     self.dat = {animal: suite2pPreprocessing(animal, self.stim, self.path, os.path.join(self.path, animal, f'track2p-{self.roi_detection}'), ntheta = self.ntheta, fps = self.fps, tracked_cells=self.tracked_cells, tracking_across_days = self.tracking_across_days, checked_tracked_cells = self.checked_tracked_cells, roi_detection = self.roi_detection, deconvolved = self.deconvolved, zscore_threshold = self.zscore_threshold, std_threshold = self.std_threshold) for animal in self.list_animals}
        # else: # want to load track2p for all the days
        #     self.dat = {animal: suite2pPreprocessing(animal, self.stim, self.path, os.path.join(self.path, animal, f'track2p-{self.stim}-{self.roi_detection}'), ntheta = self.ntheta, fps = self.fps, tracked_cells=self.tracked_cells, tracking_across_days = self.tracking_across_days, checked_tracked_cells = self.checked_tracked_cells, roi_detection = self.roi_detection, deconvolved = self.deconvolved, zscore_threshold = self.zscore_threshold, std_threshold = self.std_threshold) for animal in self.list_animals}
        # #self.analysis()
        if stim is None and (not tracking_across_days): # want to load track2p for a specifiec day
            self.dat = {animal: suite2pPreprocessing(animal, self.stim, self.path, f'track2p-{self.roi_detection}', ntheta = self.ntheta, fps = self.fps, tracked_cells=self.tracked_cells, tracking_across_days = self.tracking_across_days, checked_tracked_cells = self.checked_tracked_cells, roi_detection = self.roi_detection, deconvolved = self.deconvolved, plot_fov = self.plot_fov, zscore_threshold = self.zscore_threshold, std_threshold = self.std_threshold) for animal in self.list_animals}
        else: # want to load track2p for all the days
            self.dat = {animal: suite2pPreprocessing(animal, self.stim, self.path, f'track2p-{self.stim}-{self.roi_detection}', ntheta = self.ntheta, fps = self.fps, tracked_cells=self.tracked_cells, tracking_across_days = self.tracking_across_days, checked_tracked_cells = self.checked_tracked_cells, roi_detection = self.roi_detection, deconvolved = self.deconvolved, plot_fov = self.plot_fov, zscore_threshold = self.zscore_threshold, std_threshold = self.std_threshold) for animal in self.list_animals}
        #self.analysis()


    # def analysis(self):
    #
    #     if self.stim == 'grat':
    #         for animal in self.list_animals:
    #
    #             # take cells that pass treshold at least once > get matrix shape (n_days, n_days, n_cells)
    #             self.dat_subject['corr_matrix'] = corr_vector(self, animal, thresholded_cells=1, null_distribution=False, across_days=False)


class suite2pPreprocessing:

    def __init__ (self, animal, stim, path, track2p_folder, ntheta = 8, fps = 10, tracked_cells = False, tracking_across_days = True, checked_tracked_cells = None, roi_detection = 'functional', deconvolved = False, plot_fov = False, zscore_threshold = None, std_threshold = None):

        self.animal = animal
        self.save_path = str(Path(path).parents[0] / 'figures')
        #self.path = os.path.join(path, animal, track2p_folder)
        self.track2p_folder = track2p_folder
        self.path = path
        self.deconvolved = deconvolved              # True = load deconvolved traces (spikes) ; False = load raw fluorescence traces
        self.plot_fov = plot_fov
        self.tracked_cells = tracked_cells
        self.tracking_across_days = tracking_across_days
        if checked_tracked_cells is not None:
            self.checked_tracked_cells = checked_tracked_cells[animal]
        else:
            self.checked_tracked_cells = checked_tracked_cells
        self.roi_detection = roi_detection
        self.fps = fps
        self.ntheta = ntheta
        self.stim = stim
        self.zscore_threshold = zscore_threshold
        self.std_threshold = std_threshold

        if (not self.tracked_cells) or (self.tracking_across_days):

            animal_path = Path(self.path) / self.animal

            if (not animal_path.exists()) or (not animal_path.is_dir()):
                if Path(self.path).drive == 'I:':
                    self.path = self.path.replace('I:', 'E:')
                elif Path(self.path).drive == 'E:':
                    self.path = self.path.replace('E:', 'G:')
                elif Path(self.path).drive == 'G:':
                    self.path = self.path.replace('G:', 'E:')

            if self.tracked_cells:
                self.track2p_obj = track2pPreprocessing(self.path, self.animal, self.tracked_cells, self.roi_detection, self.track2p_folder)  # load track2p folder
                self.paths = self.track2p_obj.track_ops.all_ds_path               # list of paths to data
            else:
                self.paths = get_paths(self.stim, os.path.join(self.path, self.animal))
                print('Datasets used (using all cells for processing): \n')
                print("\n".join("\t" + item for item in self.paths))

            # look in the E drive, instead of the I drive
            self.paths = [p.replace('I:', 'E:') for p in self.paths]

            if self.paths == []:
                self.paths = [p.replace('E:', 'G:') for p in self.paths]
            if self.paths == []:
                self.paths = [p.replace('G:', 'E:') for p in self.paths]


            self.days = np.unique([path_str.split('\\')[4] for path_str in self.paths])        # list of days (strings)
            # print(self.days)
            self.dat_subject = {str(day): {} for day in self.days}    # dictionary of day (key) and empty dictionary (value) pairs

            if self.tracked_cells:
                self.load_dat(self.track2p_obj)
            else:
                self.load_dat()

        elif (not self.tracking_across_days): # track cells within days

            # take all the days in the file, but skip the first two weeks, which have no tracking data >  list of days (strings)
            self.days = sorted([x for x in os.listdir(os.path.join(path, animal)) if x.isdigit()])[2:]
            self.dat_subject = {str(day): {} for day in self.days}  # dictionary of day (key) and empty dictionary (value) pairs
            self.load_dat()

    def load_dat(self, track2p_object = None):

        for i_day, day in enumerate(self.days):

            day = str(day)

            if self.tracked_cells and (not self.tracking_across_days):  # track cells within days
                track2p_folder_day = os.path.join(self.path, self.animal, day, f'track2p-{self.roi_detection}')
                track2p_obj = track2pPreprocessing(self.path, self.animal, self.tracked_cells, self.roi_detection, track2p_folder_day, day = day)  # load track2p folder
                self.paths = track2p_obj.track_ops.all_ds_path  # list of paths to data
                self.paths = [p.replace('I:', 'E:') for p in self.paths]

                for i_recording, recording_path in enumerate(self.paths):
                    recording = recording_path.split(os.sep)[5]
                    self.dat_subject[day][recording] = {}

                    suite2p_path = os.path.join(recording_path, f'suite2p {self.roi_detection}','plane0')
                    responses, iscell, ops, stat = suite2p_files(suite2p_path, response_type='fluorescence')
                    spikes, _, _, _ = suite2p_files(suite2p_path, response_type='deconvolved')
                    self.meanimg = load_ops_meanimg(ops, 'meanImg')

                    self.dat_subject[day][recording]['traces'] = track2p_obj.all_f_t2p_traces[i_recording]
                    self.dat_subject[day][recording]['deconvolved'] = track2p_obj.all_f_t2p_spikes[i_recording]

                    self.dat_subject[day][recording]['zscored_traces'] = zscore(self.dat_subject[day][recording]['traces'], axis=1)
                    self.dat_subject[day][recording]['zscored_deconvolved'] = zscore(self.dat_subject[day][recording]['deconvolved'], axis=1)

                    if 'spon' not in recording: # only load log/ttl file if recording is 'grat' and not 'spon'

                        self.dat_subject[day][recording]['log'] = {}
                        logpath = os.path.join(os.path.dirname(recording_path), 'logfiles')
                        logfile = [file for file in os.listdir(logpath) if file.endswith('_log.txt')][0]
                        path_to_log_file = os.path.join(logpath, logfile)
                        self.dat_subject[day][recording]['log']['stim_dict'], self.dat_subject[day][recording]['log']['list_stim_types'], self.dat_subject[day][recording]['log']['list_stim_names'] = load_stims (self, path_to_log_file)
                        # load ttl / neural information
                        ttl_path = os.path.join(recording_path, f'{recording_path.split(os.sep)[-2]}.mat')
                        event_id = np.squeeze(loadmat(ttl_path)['info']['event_id'][0][0])
                        self.dat_subject[day][recording]['ttl_data'] = np.squeeze(loadmat(ttl_path)['info']['frame'][0][0])#[event_id == 1]

                        parsed = np.array([parse_grating_array(
                            self.dat_subject[day][recording]['log']['stim_dict'][grating_block]) for grating_block in
                            self.dat_subject[day][recording]['log']['stim_dict'].keys() if 'Grating' in grating_block])

                        self.dat_subject[day][recording]['thetas'] = parsed[:, :, 0]
                        self.dat_subject[day][recording]['n_repeats'] = self.dat_subject[day][recording]['thetas'].shape[0]
                        self.dat_subject[day][recording]['n_thetas'] = self.dat_subject[day][recording]['thetas'].shape[1]

                        #load ttl data, correct ttl data, loads dictionary that holds neural data for ttls values
                        get_neuronal_responses_ttls(self, day, recording)

                        #build_parameter_matrix(self, day, response_window='whole', zscore=True)

                        store_metrics(self, day, recording, response_window = 'whole', zscore = True)


            else:
                for recording_path in [p for p in self.paths if day in p]: #self.paths[i_day].split(os.sep)[5]:

                    old_animal_name = recording_path.split('\\')[3]
                    recording_path = recording_path.replace(old_animal_name, self.animal)

                    if not Path(recording_path).exists():
                        if Path(recording_path).drive == 'G:':
                            recording_path = str(Path(str(recording_path).replace('G:', 'E:', 1)))
                        elif recording_path.drive == 'E:':
                            recording_path = str(Path(str(recording_path).replace('E:', 'G:', 1)))


                    recording = recording_path.split(os.sep)[5]
                    self.dat_subject[day][recording] = {}

                    # load suite2p information
                    # load both deconvolved and fluorescence traces
                    suite2p_path = os.path.join(recording_path, f'suite2p {self.roi_detection}','plane0')

                    responses, iscell, ops, stat = suite2p_files(suite2p_path, response_type='fluorescence')
                    spikes, _, _, _ = suite2p_files(suite2p_path, response_type='deconvolved')

                    self.meanimg = load_ops_meanimg(ops, 'meanImg')
                    imE = load_ops_meanimg(ops, 'meanImgE')
                    im = load_ops_meanimg(ops, 'meanImg')

                    if self.plot_fov:
                        # out_path = fr"E:\dark_drift\figures\cellpose_images\{self.animal}_{day}.png"
                        # plt.imsave(out_path, im, cmap="gray")
                        plt.figure(figsize = (13,9))
                        plt.imshow(im,cmap = 'gray')
                        plt.title(self.animal +' '+ day)
                        plt.axis('off')
                        # plt.imsave(fr'E:\dark_drift\figures\cellpose_images\{self.animal},{day}')
                        # plt.close('all')
                        plt.tight_layout()
                        plt.show()
                        #self.m = build_roi_masks(ops, stat, iscell)

                        # overlay_masks_on_mean(self.meanimg,
                        #                       self.m[0][:,:,:20],
                        #                       labels= self.m[1][:20])

                        #shape n_rois, x_dim, y_dim
                        #self.dat_subject[day][recording]['roi_mask'] = create_roi_mask(stat, ops)

                    # store suite2p information (all cells)
                    self.dat_subject[day][recording]['meanImg'] = ops['meanImg']

                    if self.tracked_cells: # load suite2p roi data processed by track2p
                        if self.checked_tracked_cells is None:
                            self.dat_subject[day][recording]['traces'] = track2p_object.all_f_t2p_traces[i_day]
                            self.dat_subject[day][recording]['deconvolved'] = track2p_object.all_f_t2p_spikes[i_day]
                        else:
                            self.dat_subject[day][recording]['traces'] = track2p_object.all_f_t2p_traces[i_day][self.checked_tracked_cells]
                            self.dat_subject[day][recording]['deconvolved'] = track2p_object.all_f_t2p_spikes[i_day][self.checked_tracked_cells]

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

                        #logpath = logpath.replace('E:', 'G:')

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

                        store_metrics(self, day, recording, response_window = 'whole', zscore = True)

                        # calculating orientation selectivity > average responses to the same orientation within the same 'block'
                        #responses_gratings(self, day, n_grating_repeats=self.grating_repeats)

                        #self.ttl_d, self.ttls, self.stimulus_responses = stim_responses(self.fps, self.ttl_data, self.stim_dict, self.dat_subject[day]['tracked_fluorescence'].T)

    def analysis (self):
        get_complex_num(self)

        # if self.stim == 'grat':
        #     # take cells that pass treshold at least once > get matrix shape (n_days, n_days, n_cells)
        #     self.dat_subject ['corr_matrix'] = corr_vector(self, self.animal, thresholded_cells=1, null_distribution=False, across_days=False)


# animal = 'EC_GCaMP6s_06'
# path = r'I:\dark_drift_trial\data'
# track2p_folder = 'track2p'
# suite2p_obj = suite2pPreprocessing(animal, path, track2p_folder, ntheta = 8,  grating_repeats = 6, fps = 20, deconvolved = True)

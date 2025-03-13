from helpers import *

'''
Notes on the data files collected from each experiment

file organization: animal > day > recording (r1/r2/r3) > subrecording (e.g., r1_ZZZ_000) > behav/experiments/log folders

    - rX_YYY_ZZZ_eye.mat > behavioural video data (face + eye/pupil)
    - rX_YYY_ZZZ.mat > metadata; contains scan rate and other info, including ttls (stored in 'frames' variable)
    - rX_YYY_ZZZ_realtime.mat > contains data from ROIs collected during scanning (don't need)
    - rX_YYY_ZZZ.sbx > scanbox 2p data (don't need as we use this to run suite2p)
    - F.npy > array of fluorescence traces (ROIs by timepoints)
    - Fneu.npy > array of neuropil fluorescence traces (ROIs by timepoints)
    - spks.npy > array of deconvolved traces (ROIs by timepoints)
    - stat.npy > list of statistics computed for each cell (ROIs by 1)
    - ops.npy > options and intermediate outputs (dictionary)
    - iscell.npy > specifies whether an ROI is a cell, first column is 0/1, and second column is probability that the ROI is a cell based on the default classifier
    - yearmonthday-hourminutesecond_log.txt > log file specifying which stimuli are presented (includes specific info on each stimulus presented, including time of presentation)
'''

class DataAnalysis:

    def __init__(self, nat_im_path, list_animals, list_days = False, fps = 10, save_sf = True):

        self.list_animals = list_animals
        self.list_days = list_days
        self.nat_im_path = nat_im_path
        self.fps = fps              # sampling rate of the data acquisition
        self.save_sf = save_sf      # whether or not we want to process/save the spatial footprints

        self.stim_dict = {}
        self.responses_cells = {}
        self.ttl_data = {}
        self.ttls = {}
        self.stimulus_responses = {}
        self.list_stim_types = {}
        self.list_stim_names = {}
        self.preprocessing()

        self.grating_responses = {}
        self.orientations = {}
        #self.analysis()

    def preprocessing (self):
        for animal in self.list_animals:
            self.stim_dict[animal], self.responses_cells[animal],self.ttl_data[animal] = {}, {}, {}
            self.list_stim_types[animal], self.list_stim_names[animal] = {}, {}
            self.ttls[animal], self.stimulus_responses[animal] = {}, {}

            # if we do not pass self.list_days as input, then process all days in the folder
            if not self.list_days:
                #self.list_days = self.list_days
                #else:
                self.list_days = [file for file in os.listdir(os.path.join(self.nat_im_path, animal)) if ((not file.endswith('.mat')) and (not file.endswith('.png')) and (not file.endswith('.tif')))]

            for day in tqdm(self.list_days, desc = f'Processing: {animal}'):
                self.stim_dict[animal][day], self.responses_cells[animal][day],self.ttl_data[animal][day] = {}, {}, {}
                self.list_stim_types[animal][day], self.list_stim_names[animal][day] = {}, {}
                self.ttls[animal][day], self.stimulus_responses[animal][day] = {}, {}

                recordings = [file for file in os.listdir(os.path.join(self.nat_im_path, animal, day)) if ((not file.endswith('.mat')) and (not file.endswith('.png')) and (not file.endswith('.tif')))]
                for recording in recordings:
                    sub_file = [file for file in os.listdir(os.path.join(self.nat_im_path, animal, day, recording)) if ((not file.endswith('.mat')) and (not file.endswith('.png')) and (not file.endswith('.tiff')))][0]

                    # load log file information
                    log_file = [file for file in os.listdir(os.path.join(self.nat_im_path, animal, day, recording, sub_file, 'log')) if file.endswith('_log.txt')][0]
                    log_path = os.path.join(self.nat_im_path, animal, day, recording, sub_file, 'log', log_file)
                    self.stim_dict[animal][day][recording], self.list_stim_types[animal][day][recording], self.list_stim_names[animal][day][recording] = load_stims (log_path)

                    # load behaviour video
                    behav_path = os.path.join(self.nat_im_path, animal, day, recording, sub_file, 'behav_files', f'{sub_file}_eye.mj2')
                    # want to read through it, bin, and downsample it, and rewrite it to use in facemap (https://github.com/MouseLand/facemap)

                    # load suite2p information
                    suite2p_path = os.path.join(self.nat_im_path, animal, day, recording, sub_file, 'experiments', 'suite2p', 'plane0')
                    responses, iscell, ops, stat = suite2p_files (suite2p_path, response_type ='fluorescence')
                    #if you load in F, you need to subtract the neuropil (F - 0.7*Fneu) > this is the real fluorescence signal
                    #ttls > load mat

                    # #only keep responses of cells that we manually chose in the suite2p GUI
                    # self.responses_cells[animal][day][recording] = np.array([cell for i_cell, cell in enumerate(responses) if iscell[i_cell]])

                    # # create & save spatial footprint arrays in MATLAB file so that we can perform cell registration across multiple sessions in large-scale calcium imaging data
                    # if self.save_sf:
                    #     spatial_footprints = masked_spatial_footprints (ops, stat)
                    #     savemat(os.path.join(os.path.split(self.nat_im_path)[0], 'Spatial Footprints', 'ROIs', f'spatial footprint {animal}_{day}_{recording}.mat'), {f'spatial_footprints': spatial_footprints})
                    #
                    #     #generate spatial footprint images
                    #     generate_sf(animal, day, recording, stat, ops, self.nat_im_path)
                    #
                    # # load ttl / neural information
                    # ttl_path = os.path.join(self.nat_im_path, animal, day, recording, sub_file, 'experiments', f'{sub_file}.mat')
                    # event_id = np.squeeze(loadmat(ttl_path)['info']['event_id'][0][0])
                    # self.ttl_data[animal][day][recording] = np.squeeze(loadmat(ttl_path)['info']['frame'][0][0])#[event_id == 1]
                    #
                    # #print(day, recording)
                    # self.ttl_data[animal][day][recording], self.ttls[animal][day][recording], self.stimulus_responses[animal][day][recording] = stim_responses(self.fps, self.ttl_data[animal][day][recording], self.stim_dict[animal][day][recording], self.responses_cells[animal][day][recording])

    # def analysis(self):
    #
    #     for animal in self.list_animals:
    #         self.grating_responses[animal] = {}
    #         self.orientations[animal] = {}
    #         for day in list(self.ttls[animal].keys()):
    #             self.grating_responses[animal][day] = {}
    #             self.orientations[animal][day] = {}
    #             for recording in list(self.ttls[animal][day].keys()):
    #                 stim_dict_animal = self.stim_dict[animal][day][recording]
    #                 ttl_data_animal = self.ttl_data[animal][day][recording]
    #                 neural_responses_animal = self.responses_cells[animal][day][recording]
    #
    #                 self.grating_responses[animal][day][recording], self.orientations[animal][day][recording] = responses_gratings(stim_dict_animal, ttl_data_animal, neural_responses_animal, self.fps)
    #

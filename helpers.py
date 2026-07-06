import matplotlib.pyplot as plt

from imports import *

def load_stims (self, path_to_log):
    '''
    Reads all lines in log file
    :param path_to_log:
    :return:
    - stim_dict: ex, {'Wait 0': ['Wait600.00'], 'GratingStim 0': ['grating_143_SF0.02_TF1.00', 'grating_323_SF0.02_TF1.00',...]...}
    - list_stim_types: ex, ['Wait', 'GratingStim', 'GratingStim', 'GratingStim', 'GratingStim', 'GratingStim', 'GratingStim', ...]
    - list_stim_names ex, ['Wait600.00', 'grating_143_SF0.02_TF1.00', 'grating_323_SF0.02_TF1.00', 'grating_359_SF0.02_TF1.00',...]
    '''

    # read through all the line in the log file, and add each line to 'list_stims'
    list_stims = []
    with open(path_to_log) as f:
        lines = f.readlines()  # all the lines in the log file

        for i, l in enumerate(lines):  # enumerate through all lines in log file
            if i > 0:  # ignore first line, which just has info on experiment
                l_new = ast.literal_eval(l)
                list_stims.append(l_new)  # store in list

    # store the name and type of each stimulus in sepearate lists
    list_stim_types = [stim['type'] for stim in list_stims]
    list_stim_names = [stim['name'] for stim in list_stims]

    # create a dictionary where we sort the stimulus block together with the names of the stimuli (wait stims together, and grating stims together, in blocks of 10)
    stim_dict = {}
    wait_count, grating_count, image_count = 0, 0, 0
    for i_stim, stim in enumerate(list_stim_types):

        if 'Wait' in stim:
            stim_dict[f'{stim} {wait_count}'] = [list_stim_names[i_stim]]
            wait_count += 1

        elif ('Grating' in stim):
            if (grating_count%self.ntheta) == 0:
                list_gratings = []
                for i_grating in np.arange(i_stim, i_stim+self.ntheta):
                    list_gratings.append(list_stim_names[i_grating])
                stim_dict[f'{stim} {grating_count//self.ntheta}'] = list_gratings
            grating_count += 1

        elif ('Image' in stim):
            stim_dict[f'{stim} {image_count}'] = [list_stim_names[i_stim]]
            image_count += 1

    return stim_dict, list_stim_types, list_stim_names

def check_ttls (obj, ttl_arr, stim_names_list):
    '''
    Check if ttls are aligned

    :param ttl_arr:
    :return: corrected ttl_array, if necessary
    '''
    # if NOT FFF stims
    #if not (('black.png' in stim_names_list) and ('white.png' in stim_names_list)):
    #print(ttl_arr)
    difference = np.diff(ttl_arr)

    n_ttls_needed = len(stim_names_list) * 2

    # since each recording is flanked by a long wait time, the first and last 'difference' period should be very long
    # check that the wait times are each over 400 frames long and that all ttl periods are less than 200 frames long (20fps * 5s = 100frames) (i.e., all stim periods are less than 10 seconds)
    if (difference[0] > 400) and (difference[-1] > 400) and (np.all(np.diff(ttl_arr)[1:-1] < (obj.fps * 10))):
        if len(ttl_arr) + 1 == n_ttls_needed: # missing ONE ttl between end of wait and first stim
            print('Missing one TTL - appending')
            ttl_arr = np.insert(ttl_arr, 2, ttl_arr[1])
        else:
            print('TTL array is correct')

    # have 2 less ttls than we need (missing start of wait and end of wait)
    elif (len(ttl_arr) + 2 == n_ttls_needed):
        print('Missing two TTLs - appending')
        #add duplicate 'start'
        ttl_arr = np.insert(ttl_arr, 0, ttl_arr[0])

        # there is a ttl missing, so we'd want to add a 0 at the start
        ttl_arr = np.insert(ttl_arr, 0, 0)

    # have 1 less ttls than we need (missing start of first stim / end of wait)
    elif ((len(ttl_arr) + 1 == n_ttls_needed) and (ttl_arr[1]!=ttl_arr[2])):
        print('Missing one TTL - appending')
        #add duplicate 'start'
        ttl_arr = np.insert(ttl_arr, 1, ttl_arr[1])

    # have 1 MORE ttl than we need (extra start?)
    elif ((len(ttl_arr) - 1 == n_ttls_needed) and (ttl_arr[1]!=ttl_arr[2])):
        print('One extra TTL - removing')
        #remove duplicate 'start'
        ttl_arr = ttl_arr[1:]

    # else:
    #     # have 1 less ttls than we need (missing start of first stim / end of wait)
    #     n_ttls_needed = len(stim_names_list) * 2
    #     if ((len(ttl_arr) + 1 == n_ttls_needed) and (ttl_arr[1]!=ttl_arr[2])):
    #         #add duplicate 'start'
    #         ttl_arr = np.insert(ttl_arr, 1, ttl_arr[1])
    #     # if the first ttl period is longer than 4s
    #     if np.diff(ttl_arr)[0] > obj.fps * 4:
    #         ttl_arr[0] = ttl_arr[1] - obj.fps * 4

    # if we're cutting out the start and the end of the recording, the ttls need to start at '0'
    # if we dont cut out the start and end, dont shift TTLs
    return ttl_arr# - ttl_arr[0]

def zip_ttl_data(obj, object_datsubject_day):
    '''
    Put ttl data into a tuple of (stim_type, stim_name, [ttl_start, ttl_stop], ...)

    :param object_datsubject_day: example: suite2p_obj.dat_subject['20231106']
    :return:
    '''

    stim_names = object_datsubject_day['log']['list_stim_names']                            # list of stimuli names, ex: [ 'Wait600.00', 'ILSVRC2012_val_00000385.JPEG', 'grating_240_SF0.02_TF1.00'...]
    stim_types = object_datsubject_day['log']['list_stim_types']                            # list of stimuli types, ex: [ 'Wait', 'ImageStim', 'GratingStim'...]

    object_datsubject_day['ttl_data'] = check_ttls(obj, object_datsubject_day['ttl_data'], stim_names)       # correct TTL data
    ttl = object_datsubject_day['ttl_data']                                                 # 1d array of TTLs, after being corrected for data collection errors, ex:      [ 102 6086 6086 6094...]
    zipped_ttl_dat = []
    for i, (stim_name, stim_type) in enumerate(zip(stim_names, stim_types)):
        zipped_ttl_dat.append((stim_type, stim_name, [ttl[i * 2], ttl[(2 * i) + 1]]))

    return zipped_ttl_dat
def get_ttl_dicts (obj, obj_dat_subject_day):
    '''

    :param obj_dat_subject_day: example: suite2p_obj.dat_subject['20231106']
    :return:
        - dict_stim_names (dict): dictionary with stimulus type (key) and list of stimuli names (value) for each 'block'
        - dict_stim_ttls (dict): dictionary with stimulus type (key) and array of ttls (value) for each 'block' (shape n_stim_presentations x 2)
    '''

    ttl_tuple = zip_ttl_data(obj, obj_dat_subject_day)

    obj_dat_subject_day['dict_stim_names'] = {}
    obj_dat_subject_day['dict_stim_ttls'] = {}

    wait_count, grating_count, image_count = 0, 0, 0
    for i_zip in range(len(ttl_tuple)):

        if 'Wait' == ttl_tuple[i_zip][0]:
            type = ttl_tuple[i_zip][0]  # stimulus type, eg 'WaitStim'
            obj_dat_subject_day['dict_stim_names'][f'{type}_{wait_count}'] = ttl_tuple[i_zip][1] #ttl_tuple[i_zip][1] # specifications, eg Grating_90_sf0.02, etc
            obj_dat_subject_day['dict_stim_ttls'][f'{type}_{wait_count}'] = [ttl_tuple[i_zip][2]]  #ttl_tuple[i_zip][2] #list of ttls
            wait_count += 1

        elif 'ImageStim' == ttl_tuple[i_zip][0]:
            type = ttl_tuple[i_zip][0]  # stimulus type, eg 'WaitStim'
            obj_dat_subject_day['dict_stim_names'][f'{type}_{image_count}'] = ttl_tuple[i_zip][1] #ttl_tuple[i_zip][1] # specifications, eg Grating_90_sf0.02, etc
            obj_dat_subject_day['dict_stim_ttls'][f'{type}_{image_count}'] = [ttl_tuple[i_zip][2]]  #ttl_tuple[i_zip][2] #list of ttls
            image_count += 1

        # grating stim, and the remainder is 1 (to account for the wait)
        elif (i_zip % obj_dat_subject_day['n_thetas']) == 1:
            z_names = [ttl_tuple[i][1] for i in range(i_zip, i_zip + obj_dat_subject_day['n_thetas'])]
            z_ttls = [ttl_tuple[i][2] for i in range(i_zip, i_zip + obj_dat_subject_day['n_thetas'])]

            type = ttl_tuple[i_zip][0]                  # stimulus type, eg 'GratingStim'
            obj_dat_subject_day['dict_stim_names'][f'{type}_{grating_count}'] = z_names # specifications, eg Grating_90_sf0.02, etc
            obj_dat_subject_day['dict_stim_ttls'][f'{type}_{grating_count}'] = z_ttls #list of ttls
            grating_count += 1

    return obj_dat_subject_day['dict_stim_ttls']

# def zscore_baseline (obj, ):
#     self.dat[animal][day][sub_file]['baseline'] = self.dat[animal][day][sub_file]['responses'][:, :self.fps * 9]
#     baseline_mean = self.dat[animal][day][sub_file]['baseline'].mean(axis=-1, keepdims=True)
#     baseline_std = self.dat[animal][day][sub_file]['baseline'].std(axis=-1, keepdims=True)
#     self.dat[animal][day][sub_file]['baseline_zscore_responses_ttls'] =
#     (self.dat[animal][day][sub_file]['responses_ttls'] - baseline_mean) / baseline_std

def get_neuronal_responses_ttls (obj, day, recording):
    '''
    Indexes the tracked neuronal responses according to the ttl values for each stimulus

    :param obj: ex: suite2p_obj
    :param day (str): yearmonthday, ex: '20231106'
    :return: void - modifies obj in place
    '''

    # load corrected ttl data into a dictionary with stim type (key) and array of ttls (value) (shape: n_stim_presentations x 2)
    dict_stim_ttls = get_ttl_dicts(obj, obj.dat_subject[day][recording])

    fluorescence = obj.dat_subject[day][recording]['deconvolved']
    zscore_fluorescence = obj.dat_subject[day][recording]['zscored_deconvolved']
    fps = obj.fps

    # create a dictionary that holds the neural activity for each pair of TTLs (in chronological order)
    obj.dat_subject[day][recording]['responses_ttls'] = {}
    obj.dat_subject[day][recording]['responses_ttls_whole'] = {}
    obj.dat_subject[day][recording]['zscored_responses_ttls'] = {}
    obj.dat_subject[day][recording]['zscored_responses_ttls_whole'] = {}

    for stim in dict_stim_ttls.keys():
        obj.dat_subject[day][recording]['responses_ttls'][stim] = []
        obj.dat_subject[day][recording]['responses_ttls_whole'][stim] = []
        obj.dat_subject[day][recording]['zscored_responses_ttls'][stim] = []
        obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim] = []

        for [start_ttl, end_ttl] in dict_stim_ttls[stim]:
            if 'Wait' in stim:  # take entire wait period
                obj.dat_subject[day][recording]['responses_ttls'][stim].append(fluorescence[:, start_ttl:end_ttl])
                obj.dat_subject[day][recording]['responses_ttls_whole'][stim].append(fluorescence[:, start_ttl:end_ttl])
                obj.dat_subject[day][recording]['zscored_responses_ttls'][stim].append(zscore_fluorescence[:, start_ttl:end_ttl])
                obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim].append(zscore_fluorescence[:, start_ttl:end_ttl])
            elif 'Grating' in stim:  # 1s static + 3s moving + 1s off > only use 3 sec moving
                obj.dat_subject[day][recording]['responses_ttls'][stim].append(fluorescence[:, start_ttl + (fps * 1): start_ttl + (fps * 4)])
                obj.dat_subject[day][recording]['zscored_responses_ttls'][stim].append(zscore_fluorescence[:, start_ttl + (fps * 1):  start_ttl + (fps * 4)])

                # 1 s before stim + 1s static + 3s moving
                obj.dat_subject[day][recording]['responses_ttls_whole'][stim].append(fluorescence[:, start_ttl - (fps * 1): start_ttl + (fps * 4)])
                obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim].append(zscore_fluorescence[:, start_ttl - (fps * 1): start_ttl + (fps * 4)])

            elif 'Image' in stim:  # 300-500 ms of grey screen (jitter) + 0.5s on STIM + 1.5s off + 300-500 ms of grey screen (jitter)
                # for some reason ttl starts, then 16 frames elapses, then second ttl comes
                # should be 6 frames before / after
                # want to include 300ms for baseline, then 0.8s after stim start (includes stim presentaiton for 0.5s and 0.3 grey)
                # print('time diff in ttls', end_ttl - start_ttl, start_ttl, end_ttl)
                # print('actual diff in ttls', int(start_ttl + (fps * 0.8)) - int(start_ttl - (fps * 0.3)),int(start_ttl - (fps * 0.3)), int(start_ttl + (fps * 0.8)))
                obj.dat_subject[day][recording]['responses_ttls'][stim].append(fluorescence[:, int(start_ttl - (fps * 0.3)) : int(start_ttl + (fps * 0.8))])
                obj.dat_subject[day][recording]['zscored_responses_ttls'][stim].append(fluorescence[:, int(start_ttl - (fps * 0.3)) : int(start_ttl + (fps * 0.8))])

        # broadcasting list into array
        if 'Grating' in stim:  # in order to have an array of shape stimulus x cells x time, everything needs to be the same shape
            #max_size = min([l.shape[-1] for l in obj.dat_subject[day]['responses_ttls'][stim]])
            #obj.dat_subject[day]['responses_ttls'][stim] = np.array([n[:, :max_size] for n in obj.dat_subject[day]['responses_ttls'][stim]])
            obj.dat_subject[day][recording]['responses_ttls'][stim] = np.array([n for n in obj.dat_subject[day][recording]['responses_ttls'][stim]])
            obj.dat_subject[day][recording]['responses_ttls_whole'][stim] = np.array([n for n in obj.dat_subject[day][recording]['responses_ttls_whole'][stim]])
            obj.dat_subject[day][recording]['zscored_responses_ttls'][stim] = np.array([n for n in obj.dat_subject[day][recording]['zscored_responses_ttls'][stim]])
            obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim] = np.array([n for n in obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim]])

        elif 'Image' in stim:
            obj.dat_subject[day][recording]['responses_ttls'][stim] = np.array(obj.dat_subject[day][recording]['responses_ttls'][stim])
            # obj.dat_subject[day][recording]['responses_ttls_whole'][stim] = np.array(obj.dat_subject[day][recording]['responses_ttls_whole'][stim])
            obj.dat_subject[day][recording]['zscored_responses_ttls'][stim] = np.array(obj.dat_subject[day][recording]['zscored_responses_ttls'][stim])
            # obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim] = np.array(obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim])

        elif 'Wait' in stim:  # outputs a singleton dimension, so we want to use np.squeeze
            obj.dat_subject[day][recording]['responses_ttls'][stim] = np.squeeze(np.array(obj.dat_subject[day][recording]['responses_ttls'][stim]))
            obj.dat_subject[day][recording]['zscored_responses_ttls'][stim] = np.squeeze(np.array(obj.dat_subject[day][recording]['responses_ttls'][stim]))
            obj.dat_subject[day][recording]['responses_ttls_whole'][stim] = np.squeeze(np.array(obj.dat_subject[day][recording]['responses_ttls_whole'][stim]))
            obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim] = np.squeeze(np.array(obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim]))

        obj.dat_subject[day][recording]['responses_ttls'][stim] = np.array(obj.dat_subject[day][recording]['responses_ttls'][stim])
        obj.dat_subject[day][recording]['responses_ttls_whole'][stim] = np.array(obj.dat_subject[day][recording]['responses_ttls_whole'][stim])
        obj.dat_subject[day][recording]['zscored_responses_ttls'][stim] = np.array(obj.dat_subject[day][recording]['zscored_responses_ttls'][stim])
        obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim] = np.array(obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim])



def get_behav_response_ttls (obj, day, recording):
    '''
    Indexes the tracked neuronal responses according to the ttl values for each stimulus

    :param obj: ex: suite2p_obj
    :param day (str): yearmonthday, ex: '20231106'
    :return: void - modifies obj in place
    '''

    # load corrected ttl data into a dictionary with stim type (key) and array of ttls (value) (shape: n_stim_presentations x 2)
    dict_stim_ttls = get_ttl_dicts(obj, obj.dat_subject[day][recording])

    # print(dict_stim_ttls)

    pupil = obj.dat_subject[day][recording]['pupil_pos']
    face = obj.dat_subject[day][recording]['face_pos']
    fluorescence = obj.dat_subject[day][recording]['deconvolved']
    zscore_fluorescence = obj.dat_subject[day][recording]['zscored_deconvolved']
    fps = obj.fps

    # create a dictionary that holds the neural activity for each pair of TTLs (in chronological order)
    obj.dat_subject[day][recording]['pupil_keypoints']= {}
    obj.dat_subject[day][recording]['face_keypoints'] = {}
    obj.dat_subject[day][recording]['responses_ttls'] = {}
    obj.dat_subject[day][recording]['responses_ttls_whole'] = {}
    obj.dat_subject[day][recording]['zscored_responses_ttls'] = {}
    obj.dat_subject[day][recording]['zscored_responses_ttls_whole'] = {}

    for stim in dict_stim_ttls.keys():
        obj.dat_subject[day][recording]['pupil_keypoints'][stim]= []
        obj.dat_subject[day][recording]['face_keypoints'][stim] = []
        obj.dat_subject[day][recording]['responses_ttls'][stim] = []
        obj.dat_subject[day][recording]['responses_ttls_whole'][stim] = []
        obj.dat_subject[day][recording]['zscored_responses_ttls'][stim] = []
        obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim] = []

        for [start_ttl, end_ttl] in dict_stim_ttls[stim]:
            if 'Wait' in stim:  # take entire wait period
                obj.dat_subject[day][recording]['pupil_keypoints'][stim].append(pupil[start_ttl:end_ttl])
                obj.dat_subject[day][recording]['face_keypoints'][stim].append(face[start_ttl:end_ttl])
                obj.dat_subject[day][recording]['responses_ttls'][stim].append(fluorescence[:, start_ttl:end_ttl])
                obj.dat_subject[day][recording]['responses_ttls_whole'][stim].append(fluorescence[:, start_ttl:end_ttl])
                obj.dat_subject[day][recording]['zscored_responses_ttls'][stim].append(zscore_fluorescence[:, start_ttl:end_ttl])
                obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim].append(zscore_fluorescence[:, start_ttl:end_ttl])
            elif 'Grating' in stim:  # 1s static + 3s moving + 1s off > only use 3 sec moving
                obj.dat_subject[day][recording]['pupil_keypoints'][stim].append(pupil[start_ttl + (fps * 1): start_ttl + (fps * 4)])
                obj.dat_subject[day][recording]['face_keypoints'][stim].append(face[start_ttl + (fps * 1): start_ttl + (fps * 4)])

                obj.dat_subject[day][recording]['responses_ttls'][stim].append(fluorescence[:, start_ttl + (fps * 1): start_ttl + (fps * 4)])
                obj.dat_subject[day][recording]['zscored_responses_ttls'][stim].append(zscore_fluorescence[:, start_ttl + (fps * 1):  start_ttl + (fps * 4)])

                # 1 s before stim + 1s static + 3s moving
                obj.dat_subject[day][recording]['responses_ttls_whole'][stim].append(fluorescence[:, start_ttl - (fps * 1): start_ttl + (fps * 4)])
                obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim].append(zscore_fluorescence[:, start_ttl - (fps * 1): start_ttl + (fps * 4)])

            elif 'Image' in stim:  # 300-500 ms of grey screen (jitter) + 0.5s on STIM + 1.5s off + 300-500 ms of grey screen (jitter)
                # for some reason ttl starts, then 16 frames elapses, then second ttl comes
                # should be 6 frames before / after
                # want to include 300ms for baseline, then 0.8s after stim start (includes stim presentaiton for 0.5s and 0.3 grey)
                # print('time diff in ttls', end_ttl - start_ttl, start_ttl, end_ttl)
                # print('actual diff in ttls', int(start_ttl + (fps * 0.8)) - int(start_ttl - (fps * 0.3)),int(start_ttl - (fps * 0.3)), int(start_ttl + (fps * 0.8)))
                obj.dat_subject[day][recording]['pupil_keypoints'][stim].append(
                    pupil[int(start_ttl - (fps * 0.3)): int(start_ttl + (fps * 0.8))])
                obj.dat_subject[day][recording]['face_keypoints'][stim].append(
                    face[int(start_ttl - (fps * 0.3)): int(start_ttl + (fps * 0.8))])
                obj.dat_subject[day][recording]['responses_ttls'][stim].append(fluorescence[:, int(start_ttl - (fps * 0.3)) : int(start_ttl + (fps * 0.8))])
                obj.dat_subject[day][recording]['zscored_responses_ttls'][stim].append(fluorescence[:, int(start_ttl - (fps * 0.3)) : int(start_ttl + (fps * 0.8))])

        # broadcasting list into array
        if 'Grating' in stim:  # in order to have an array of shape stimulus x cells x time, everything needs to be the same shape
            #max_size = min([l.shape[-1] for l in obj.dat_subject[day]['responses_ttls'][stim]])
            #obj.dat_subject[day]['responses_ttls'][stim] = np.array([n[:, :max_size] for n in obj.dat_subject[day]['responses_ttls'][stim]])
            obj.dat_subject[day][recording]['pupil_keypoints'][stim] = np.array(
                [n for n in obj.dat_subject[day][recording]['pupil_keypoints'][stim]])
            obj.dat_subject[day][recording]['face_keypoints'][stim] = np.array(
                [n for n in obj.dat_subject[day][recording]['face_keypoints'][stim]])
            obj.dat_subject[day][recording]['responses_ttls'][stim] = np.array([n for n in obj.dat_subject[day][recording]['responses_ttls'][stim]])
            obj.dat_subject[day][recording]['responses_ttls_whole'][stim] = np.array([n for n in obj.dat_subject[day][recording]['responses_ttls_whole'][stim]])
            obj.dat_subject[day][recording]['zscored_responses_ttls'][stim] = np.array([n for n in obj.dat_subject[day][recording]['zscored_responses_ttls'][stim]])
            obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim] = np.array([n for n in obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim]])

        elif 'Image' in stim:
            obj.dat_subject[day][recording]['pupil_keypoints'][stim] = np.array(
                obj.dat_subject[day][recording]['pupil_keypoints'][stim])
            obj.dat_subject[day][recording]['face_keypoints'][stim] = np.array(
                obj.dat_subject[day][recording]['face_keypoints'][stim])
            obj.dat_subject[day][recording]['responses_ttls'][stim] = np.array(obj.dat_subject[day][recording]['responses_ttls'][stim])
            # obj.dat_subject[day][recording]['responses_ttls_whole'][stim] = np.array(obj.dat_subject[day][recording]['responses_ttls_whole'][stim])
            obj.dat_subject[day][recording]['zscored_responses_ttls'][stim] = np.array(obj.dat_subject[day][recording]['zscored_responses_ttls'][stim])
            # obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim] = np.array(obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim])

        elif 'Wait' in stim:  # outputs a singleton dimension, so we want to use np.squeeze
            obj.dat_subject[day][recording]['responses_ttls'][stim] = np.squeeze(np.array(obj.dat_subject[day][recording]['responses_ttls'][stim]))
            obj.dat_subject[day][recording]['pupil_keypoints'][stim] = np.squeeze(
                np.array(obj.dat_subject[day][recording]['pupil_keypoints'][stim]))
            obj.dat_subject[day][recording]['face_keypoints'][stim] = np.squeeze(
                np.array(obj.dat_subject[day][recording]['face_keypoints'][stim]))
            obj.dat_subject[day][recording]['zscored_responses_ttls'][stim] = np.squeeze(np.array(obj.dat_subject[day][recording]['responses_ttls'][stim]))
            obj.dat_subject[day][recording]['responses_ttls_whole'][stim] = np.squeeze(np.array(obj.dat_subject[day][recording]['responses_ttls_whole'][stim]))
            obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim] = np.squeeze(np.array(obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim]))

        obj.dat_subject[day][recording]['pupil_keypoints'][stim] = np.array(
            obj.dat_subject[day][recording]['pupil_keypoints'][stim])
        obj.dat_subject[day][recording]['face_keypoints'][stim] = np.array(
            obj.dat_subject[day][recording]['face_keypoints'][stim])
        obj.dat_subject[day][recording]['responses_ttls'][stim] = np.array(obj.dat_subject[day][recording]['responses_ttls'][stim])
        obj.dat_subject[day][recording]['responses_ttls_whole'][stim] = np.array(obj.dat_subject[day][recording]['responses_ttls_whole'][stim])
        obj.dat_subject[day][recording]['zscored_responses_ttls'][stim] = np.array(obj.dat_subject[day][recording]['zscored_responses_ttls'][stim])
        obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim] = np.array(obj.dat_subject[day][recording]['zscored_responses_ttls_whole'][stim])


def parse_grating_string(s):
    '''
    translate grating string to extract orientation, TF, and SF
    :param s:
    :return:
    '''
    match = re.search(r'grating_(\d+)_SF(\d+\.\d+)_TF(\d+\.\d+)', s)
    if match:
        # return {
        #     'deg': int(match.group(1)),
        #     'SF': float(match.group(2)),
        #     'TF': float(match.group(3))
        # }
        deg = int(match.group(1))
        sf = float(f"{float(match.group(2)):.3f}")  # Format to 2 decimal places
        tf = float(f"{float(match.group(3)):.3f}")  # Format to 2 decimal places
        sf = 0.005 if sf == 0.01 else sf        # for some reasont the 0.005 SF string format is wrong
        return [deg, sf, tf]
    return None

def parse_grating_array (arr):
    np.set_printoptions(suppress=True, precision=4)
    return np.array([parse_grating_string(s) for s in arr])

def build_parameter_matrix(object, day, recording, response_window = 'whole', zscore = True):
    '''
    :param object:  data_object
    :param animal:
    :param day:
    :param subfile:
    :param response_window: either 'whole' (1s wait + 1s static + 3s moving + 1s wait) or 'moving' (3 s moving)
    :param zscore: True or False: whether or not we want z scored responses in our calculations
    :return:
    '''

    # each of shape (n_repeats, x nORI) > parameters for stimuli shown for each repeat
    thetas_stim = object.dat_subject[day][recording]['thetas']

    # array of orientations shown (shape n_orientations), in increasing order
    thetas = np.unique(thetas_stim)

    # dictionary - 1 key per repeat. each repeat contains array of shape (nSF x nOri, timepoints)
    if response_window == 'whole':      # 1s wait + 1s static + 3s moving + 1s wait
        if zscore:
            responses_ttls = object.dat_subject[day][recording]['zscored_responses_ttls_whole']
        else:
            responses_ttls = object.dat_subject[day][recording]['responses_ttls_whole']
    elif response_window == 'moving':   # 3s moving
        if zscore:
            responses_ttls = object.dat_subject[day][recording]['zscored_responses_ttls']
        else:
            responses_ttls = object.dat_subject[day][recording]['responses_ttls']

    # array shape (n_repeats, nORI, n_cells, timepoints)
    gratings_responses = np.array([responses_ttls[key] for key in responses_ttls.keys() if 'Grating' in key])

    n_repeats = gratings_responses.shape[0]
    n_orientation = len(thetas)
    cells = gratings_responses.shape[2]
    timepoints = gratings_responses.shape[-1]

    # the reordered array> contains response matrix for each repeat, with sorted orientation and SF
    responses_ordered = np.zeros((n_repeats, n_orientation, cells, timepoints))

    for i_repeat in range(n_repeats):
        for i_theta, theta in enumerate(thetas):
            # find where this sf/orientation combination exists in thetas_repeat
            idx = np.where(thetas_stim[i_repeat] == theta)[0][0]

            # store response at this index
            responses_ordered[i_repeat, i_theta] = gratings_responses[i_repeat, idx]

    # store in object
    object.dat_subject[day][recording][f'param_matrix_{response_window}{zscore*"_zscore"}'] = responses_ordered

    return responses_ordered, thetas


def build_tuning_curves (object, day, recording, zscore = False):
    '''
    :param object:
    :param animal:
    :param day:
    :param zscore: fine to z-score when building TCs, but can't z-score if trying to calculate the complex phase
    :return:
        - tuning_curve: array of shape ((n_repeats, n_orientation/n_SF, cells))
        - theta: array of shape n_theta, theta of gratings in increasing order
    '''

    # response: shape ((n_repeats, n_orientation, cells, timepoints)) > ordered responses
    # 1s wait + 1s static + 3s moving + 1s static
    response, thetas = build_parameter_matrix(object, day, recording, response_window='whole', zscore=zscore)

    # baseline period is the average over the first second of stimulus (WAIT period) > shape ((n_repeats, n_orientation/n_SF, cells))
    baseline = response[:, :, :, :object.fps].mean(axis=-1)

    # moving period is average of the first 2 seconds moving > shape ((n_repeats, n_orientation/n_SF, cells))
    tuning_curve_moving = response[:, :, :, 2 * (object.fps):4 * (object.fps)].mean(axis=-1)

    # baseline subtracted tuning curve, ensuring that all responses are non-negative> shape (on_repeats, n_orientation/n_SF, cells))
    tuning_curve = np.maximum((tuning_curve_moving - baseline), 0)

    return response, tuning_curve, thetas




def deltaf_over_f(obj, param_matrix):
    '''
    calculate delta f / f (Orientations x Repeats)

    :param param_matrix: shape (n_repeats, n_ori, n_cells, n_timepoints)
    :return: dF/F matrix shape (n_repeats, n_ori, n_cells, n_timepoints)
    '''

    # take the mean over the first second of the resposne, across all dimensions except cells
    # baseline_F is of shape n_cells
    baseline_F = param_matrix[..., :obj.fps].mean(axis=(0, 1, -1), keepdims=True)

    return (param_matrix - baseline_F) / baseline_F


def zscore_baseline(obj, param_matrix):
    '''
    z-score the parameter matrix according to the baseline period for each trial (ori x repeats)

    :param param_matrix: shape (n_repeats, n_ori, n_cells, n_timepoints)
    :return: z scored matrix shape (n_repeats, n_ori, n_cells, n_timepoints)
    '''

    # take first second of response (baseline) > mean/std across time/orientation/repeat axis
    # baseline_mean and baseline_std are of shape (n_cells)
    baseline_mean = param_matrix[..., :obj.fps].mean(axis=(0, 1, -1), keepdims=True)
    baseline_std = param_matrix[..., :obj.fps].std(axis=(0, 1, -1), keepdims=True)

    # z score the response, and then average across repeats (shape = ((n_repeats, n_ori, n_cells, n_timepoints)
    return ((param_matrix - baseline_mean) / baseline_std)

def get_cell_by_cell_threshold (obj, param_matrix, zscore_threshold = None, std_threshold = None): #TO DO !!!!!
    '''
    :param param_matrix: shape = ((n_repeats, n_ori, n_cells, n_timepoints))
    :param zscore_threshold:
    :param std_threshold:
    :return: cell_theshold, shape ((n_cells))
    '''
    if zscore_threshold:  # compute z score of entire trace, then pick some global threshold (for all cells) of how many std above baseline

        z_scored_response = zscore_baseline(obj, param_matrix) # array of responses, shape ((n_repeats, n_ori, ..., n_cells, timepoints))
        cell_threshold = np.ones((z_scored_response.shape[-2])) * zscore_threshold

    elif std_threshold: # compute individual (cell-by-cell) threshold
        # either compute DC offset, or deltaF/

        # DELTA F / F, shape = ((n_repeats, n_ori, n_cells, n_timepoints
        dFoF_response = deltaf_over_f(obj, param_matrix) # array of responses, shape ((n_repeats, n_ori, ..., n_cells, timepoints))

        # first second of the response (baseline period) > std over orientations, repeats and timepoints > multiply each cell by threshold
        cell_threshold = dFoF_response[..., :obj.fps].std(axis=(0, 1, -1)) * std_threshold
    #
    # else:
    #     return None
    return cell_threshold


def zscore_thresholding (object, day, recording, zscore_threshold = None, std_threshold = None):
    '''
    Determines which cells were responsive

    :param object:
    :param zscore_threshold: 2
    :param std_threshold: 1.5

    '''

    # since we're z scoring according to the baseline period, want to grab raw responses so we don't z score twice
    # shape (n_repeats, n_ori, n_cells, timepoints)
    parameter_matrix, _ = build_parameter_matrix(object, day, recording, response_window='whole', zscore=False)

    # ZSCORED RESPONSES, shape = ((n_repeats, n_ori, n_cells, n_timepoints
    z_scored_response = zscore_baseline(object, parameter_matrix)

    # DELTA F / F, shape = ((n_repeats, n_ori, n_cells, n_timepoints
    dFoF_response = deltaf_over_f (object, parameter_matrix)

    cell_threshold = get_cell_by_cell_threshold(object, parameter_matrix, zscore_threshold=zscore_threshold, std_threshold=std_threshold)

    # shape n_cells > take the z scored responses during the 'moving grating' period > average (n_ori, n_cells, n_timepoints) over orientations/repeats to get shape (n_cells)
    #z_scored_response_mean = z_scored_response[..., object.fps * 2: object.fps * 5].mean(axis=(0,-1))
    # shape (n_ori, n_cells)
    z_scored_response_mean = z_scored_response[..., object.fps * 2: object.fps * 5].mean(axis=(0,-1))
    z_scored_response_median = np.median(z_scored_response[..., object.fps * 2: object.fps * 5], axis=(0,-1))

    # see if the mean response during moving exceeds the cell std threshold
    #cell_exceeds_threshold = z_scored_response_mean > cell_threshold
    # shape (n_ori, n_cells) > (n_cells)
    # cell_exceeds_threshold = (z_scored_response_mean > cell_threshold[None,:]).any(axis=0) & (
    #             z_scored_response_median > cell_threshold[None,:] * 0.4).any(axis=0)

    # shape (n_ori, n_cells) > (n_cells)
    cell_exceeds_threshold = (z_scored_response_mean > cell_threshold[None,:]).any(axis=0)

    return z_scored_response, dFoF_response, cell_exceeds_threshold


def responsive_cells(obj):

    groups = np.unique([a.split('_')[1] for a in obj.list_animals])

    # {'EB': {animal 1: [], animal 2:[]}],'EBc': {animal 1: [] ...}, 'LB':{animal 1: [] ...}, 'LBc': {animal 1: [] ...}}
    percent_responsive_groups = {g: {} for g in groups}

    for i_animal, animal in enumerate(obj.dat):

        group = animal.split('_')[1]
        days_recordings = [(day, subfile) for day in obj.dat[animal].dat_subject for subfile in
                           obj.dat[animal].dat_subject[day] if 'grat' in subfile]

        # if there are no gratings for the animal, skip it
        if not days_recordings:
            continue

        p_resp = {}

        for i, (day, subfile) in enumerate(days_recordings):

            timepoint = calculate_animal_age(obj.animal_dobs[animal], day)

            # shape n_Cells, n_timepoints
            percent_cells = 100 * obj.dat[animal].dat_subject[day][subfile]['thresholded_cells'].sum() / len(obj.dat[animal].dat_subject[day][subfile]['thresholded_cells'])
            p_resp[timepoint] = percent_cells

        percent_responsive_groups [group] [animal] = p_resp

    return percent_responsive_groups


def trial_reliability(response_matrix):
    """
    Compute average trial-by-trial reliability (correlation) for each orientation and cell.

    :param response_matrix: array of shape (n_repeats, n_orientations, n_cells, n_timepoints)
    :return: reliability matrix of shape (n_orientations, n_cells)
    """
    n_repeats, n_orientations, n_cells, n_timepoints = response_matrix.shape
    reliability = np.full((n_orientations, n_cells), np.nan) # array full of nans of shape (n_ori, n_cells)

    for o in range(n_orientations):
        for c in range(n_cells):
            trials = response_matrix[:, o, c, :]  # shape (n_repeats, n_timepoints)

            # correlating mean of first half of trials with mean of second half of trials
            n_repeats = trials.shape[0]
            half_repeats = n_repeats//2
            trials_half1, trials_half2 = trials[:half_repeats].mean(axis = 0),  trials[half_repeats:].mean(axis = 0) # each of shape n_timepoints
            trials = np.vstack((trials_half1, trials_half2))
            corr_matrix = np.corrcoef(trials)  # shape (n_repeats, n_repeats) > correlate each pair of trials

            # if we are correlating each trial with every other trial
            #corr_matrix = np.corrcoef(trials)  # shape (n_repeats, n_repeats) > correlate each pair of trials
            iu = np.triu_indices(trials.shape[0], k=1) # indices of upper triangle of (n_repeats x n_repeats) matrix, exclude self-corr (diagonals)
            reliability[o, c] = np.nanmean(corr_matrix[iu]) # take mean of upper triangle

    return reliability

def trial_reliability_3(response_matrix):
    """
    Split repeats down the middle and compute one tuning-reliability value per cell.

    Parameters
    ----------
    response_matrix : array
        Shape (n_repeats, n_orientations, n_cells, n_timepoints)

    Returns
    -------
    reliability : array
        Shape (n_cells,)
    """
    resp = response_matrix.mean(axis=-1)   # (n_repeats, n_orientations, n_cells)
    n_repeats, n_orientations, n_cells = resp.shape
    reliability = np.full(n_cells, np.nan)

    half = n_repeats // 2

    for c in range(n_cells):
        trials = resp[:, :, c]             # (n_repeats, n_orientations)

        half1 = trials[:half].mean(axis=0) # (n_orientations,)
        half2 = trials[half:].mean(axis=0) # (n_orientations,)

        if np.std(half1) > 0 and np.std(half2) > 0:
            reliability[c] = np.corrcoef(half1, half2)[0, 1]

    return reliability

def trial_reliability_2(response_matrix):
    """
    response_matrix: (n_repeats, n_orientations, n_cells, n_timepoints)
    returns: (n_cells,)
    """
    resp = response_matrix.mean(axis=-1)  # (n_repeats, n_orientations, n_cells)
    n_repeats, n_orientations, n_cells = resp.shape

    reliability = np.full(n_cells, np.nan)

    for c in range(n_cells):
        trials = resp[:, :, c]   # (n_repeats (10), n_orientations (8))
        corr = np.corrcoef(trials) # shape n_repeats, n_repeats > similarity matrix showing how close each repeat (all orientations) is to every other one
        iu = np.triu_indices_from(corr, k=1)
        reliability[c] = np.nanmean(corr[iu]) # average upper triangular values from the similarity matrix

    return reliability # shape n_cells

def complex_phase_from_tuning (tuning_curves, thetas):
    '''
    :param tuning_curves: shape ((n_repeats, n_orientation/n_SF, cells))
    :param thetas: shape ((n_orientation))
    :param sfs: shape ((n_SF))
    :return:
    '''

    # shape n_orientations
    theta_rad = np.deg2rad(thetas)

    # mean response to the stimulus across repeats > shape: (n_cells, n_orientations/n_SF)
    # # HOW DO I GET OSI FOR AVERAGE ACROSS ORIENTATINO)
    responses_theta = tuning_curves.mean(axis=0).T

    # shape (n_cells);  NB: np.exp outputs shape n_orientations
    # in python, j is the imaginary component
    complex_ori = (responses_theta * np.exp(1j * 2 * theta_rad)).sum(axis=-1) / (responses_theta.sum(axis=-1))
    # note that e^(i2theta) doubles the angles because orientations are periodic over 180deg, not 360deg.
    # this ensures that responses to 0deg and 180deg are treated as equivalent (represent same orientation)
    # consequently, the resultant vector encodes the doubled angle of the preferred orientation.
    # note that to calculate direction selectivity you wouldnt have the 2 factor

    # shape (n_cells) > magnitude (absolute value) of the complex number > represents the amplitude of the response (how strongly neuron responds) = OSI (should be between 0 & 1)
    osi = np.abs(complex_ori)
    osi = np.nan_to_num(osi, nan=0) #if the OSI is nan, put it to 0

    # the argument (angle / preferred angle) of the complex number encodes the preferred orientation (the orientation to which the neuron responds most strongly).
    pref_orientation = np.rad2deg(np.angle(complex_ori))/2  # preferred angle
    # divide by 2 because of the i*theta*2 in the original equation
    pref_orientation = np.nan_to_num(pref_orientation, nan=0)  # if the pref_orientation is nan, put it to 0

    complex_direction = (responses_theta * np.exp(1j * theta_rad)).sum(axis=-1) / (responses_theta.sum(axis=-1))
    dsi = np.abs(complex_direction)  # ranges between 0 and 1, where 0 = no direction selectivity, 1 = perfect direction selectivity.
    pref_direction = np.rad2deg(np.angle(complex_direction))  # neuron’s preferred direction.
    dsi = np.nan_to_num(dsi, nan=0)  # if the pref_orientation is nan, put it to 0
    pref_direction = np.nan_to_num(pref_direction, nan=0)  # if the pref_orientation is nan, put it to 0

    return (complex_ori, osi, pref_orientation), (complex_direction, dsi, pref_direction)

def store_metrics(object, day, recording, response_window = 'whole', zscore = True):
    '''
    Calculates and stores several important metrics to be plotted by plot_parameter_matrix
    '''

    # param_matrix: shape ((n_repeats, n_orientation, n_sf, cells, timepoints))
    param_matrix, object.dat_subject[day][recording]['orientations'] = build_parameter_matrix(object, day, recording, response_window=response_window, zscore=zscore)

    _, object.dat_subject[day][recording]['tuning_curves'], _ = build_tuning_curves(object, day, recording, zscore=False)

    orientation, direction = complex_phase_from_tuning(object.dat_subject[day][recording]['tuning_curves'], object.dat_subject[day][recording]['orientations'])

    _, object.dat_subject[day][recording]['OSI'], object.dat_subject[day][recording]['preferred_orientation'] = orientation
    _, object.dat_subject[day][recording]['DSI'], object.dat_subject[day][recording]['preferred_direction'] = direction

    object.dat_subject[day][recording]['zscored'], object.dat_subject[day][recording]['dFoF'], object.dat_subject[day][recording]['thresholded_cells'] =  zscore_thresholding(object, day, recording, zscore_threshold=object.zscore_threshold, std_threshold=object.std_threshold)


def corr_vector (obj, animal, thresholded_cells = 0, null_distribution = False, n = 1000, across_days = False):
    '''
    :param obj:
    :param animal (str):
    :param thresholded_cells: either 0 (take all cells for the analysis), 1 (cells that pass treshold at least once), 2 (cells that pass threshold on all days)
    :param null_distribution (bool):
    :param n (int):
    :param across_days: True if we want to get 1 correlation value across all days(average the matrix). False if we want to compare separate days
    :return: correlation matrix, shape (n_days, n_days, n_cells)

    NB: deconvolved must be true
    '''

    dat_object = obj.dat[animal]

    #only looking in days with gratings
    days_recordings = [(day, subfile) for day in dat_object.dat_subject for subfile in dat_object.dat_subject[day] if 'grat' in subfile]
    days = [d[0] for d in days_recordings]
    #print(obj.checked_tracked_cells, type(obj.checked_tracked_cells))
    if thresholded_cells == 0: # take all cells, regardless of whether or not they pass threshold at any point
        # if obj.checked_tracked_cells is None:
        #     n_cells = dat_object.track2p_obj.track_ops_dict['n_tracked']
        #     cell_indices = range(n_cells)
        # else:
        # n_cells = len(obj.dat[animal].checked_tracked_cells) # dat_object.track2p_obj.track_ops_dict['n_tracked']
        # cell_indices = range(n_cells)
        n_cells = dat_object.track2p_obj.track_ops_dict['n_tracked']
        cell_indices = range(n_cells)


    elif thresholded_cells == 1: # only take cells that pass threshold on a minimum of 1 day
        cell_indices = np.where(np.array([obj.dat[animal].dat_subject[day][subfile]['thresholded_cells'] for (day, subfile) in days_recordings]).any(axis = 0))[0]
        n_cells = len(cell_indices)
    elif thresholded_cells == 2:
        cell_indices = np.where(np.array([obj.dat[animal].dat_subject[day][subfile]['thresholded_cells'] for (day, subfile) in days_recordings]).all(axis=0))[0]
        n_cells = len(cell_indices)
    else: # default is to just take all cells
        n_cells = dat_object.track2p_obj.track_ops_dict['n_tracked']
        cell_indices = range(n_cells)

     # if using deconvolved traces, z scoring doesnt make sense (cant have neg values). only use z score if deconvovled = false
    #param_matrix_key = 'param_matrix_whole' if dat_object.deconvolved else 'param_matrix_whole_zscore'
    param_matrix_key = 'zscored'

    if across_days:     # 1 correlation value per cell
        corr_values = np.zeros(n_cells)
    else:               # 1 correlation value per day per cell
        corr_values = np.zeros((len(days), len(days), n_cells))

    for i_cell, cell_num in enumerate(tqdm(cell_indices, desc = 'Calculating correlation distribution for each cell')):
        correlation_matrix = np.zeros((len(days), len(days)))
        for i, day_i in enumerate(days):
            for j, day_j in enumerate(days):

                if i != j: # if we are comparing two different days
                    if null_distribution:
                        # param_matrix_whole_zscore is : shape (n_repeats, n_orientations, n_cells, n_timepoints)
                        # response vector i > (shape n_orientations) > average response across time
                        vector_i = dat_object.dat_subject[day_i]['grat'][param_matrix_key].mean(axis=(0,-1))[:, cell_num]

                        null_dist = np.zeros(n)
                        # for the null distribution: shuffle cells across days (or just pick random cells)
                        for i_n in range(n):
                            # output of interleave_responses is : shape (n_orientations, n_cells, n_timepoints)
                            # response vector j > (shape n_orientations) > average response across time
                            vector_j = dat_object.dat_subject[day_i]['grat'][param_matrix_key].mean(axis=(0,-1))[:, np.random.choice(cell_indices)]

                            #roll_by = np.random.randint(obj.dat_subject[day_i]['mean_ordered_grat_responses'].shape[-1])
                            #vector_j = obj.dat_subject[day_i]['mean_ordered_grat_responses'][:, cell].mean(axis=-1)

                            corr_coef, p_value = pearsonr(vector_i, vector_j)
                            null_dist[i_n] = corr_coef

                        correlation_matrix[i, j] = null_dist.mean()

                    else:
                        # response vector (shape n_orientations) > average across time (average response)
                        vector_i = dat_object.dat_subject[day_i]['grat'][param_matrix_key].mean(axis=(0,-1))[:,cell_num]
                        vector_j = dat_object.dat_subject[day_j]['grat'][param_matrix_key].mean(axis=(0,-1))[:,cell_num]

                        corr_coef, p_value = pearsonr(vector_i, vector_j)
                        correlation_matrix[i, j] = corr_coef
                elif i == j: #comparing the same day, want to do a split-half reliability metrix
                    if null_distribution:

                        param_data = dat_object.dat_subject[day_i]['grat'][param_matrix_key]  # (n_repeats, n_orientations, n_cells, n_timepoints)

                        # Split data into first half and second half along the first axis (repeats/trials)
                        half = param_data.shape[0] // 2
                        vector_first_half = param_data[:half, :, cell_num, :].mean(axis=(0, -1))  # (n_orientations,)

                        null_dist = np.zeros(n)
                        # for the null distribution: shuffle cells across days (or just pick random cells)
                        for i_n in range(n):
                            vector_second_half = param_data[half:, :, np.random.choice(cell_indices), :].mean(axis=(0, -1))  # (n_orientations,)

                            # Compute Pearson correlation between halves
                            corr_coef, _ = pearsonr(vector_first_half, vector_second_half)
                            null_dist[i_n] = corr_coef

                        correlation_matrix[i, j] = null_dist.mean()

                    else:
                        param_data = dat_object.dat_subject[day_i]['grat'][param_matrix_key]  # (n_repeats, n_orientations, n_cells, n_timepoints)

                        # Split data into first half and second half along the first axis (repeats/trials)
                        half = param_data.shape[0] // 2
                        vector_first_half = param_data[:half, :, cell_num, :].mean(axis=(0, -1))  # (n_orientations,)
                        vector_second_half = param_data[half:, :, cell_num, :].mean(axis=(0, -1))  # (n_orientations,)

                        # Compute Pearson correlation between halves
                        corr_coef, _ = pearsonr(vector_first_half, vector_second_half)
                        correlation_matrix[i, j] = corr_coef

        if across_days:
            # corr_values is of shape n_cells
            # correlation_matrix is of shape (n_days, n_days)
            # want to exclude within-day comparisons
            corr_values[i_cell] = correlation_matrix[np.triu_indices(correlation_matrix.shape[0], k=1)].mean()

        else:
            # corr_values is of shape (n_days, n_cells)
            # correlation_matrix is of shape (n_days, n_days) > only take 1st row because it compares day 1 with all other days > changed to taking whole array (also want to compare day2 with itself)
            corr_values[:, :, i_cell] = np.array(correlation_matrix)

    return corr_values


#######################################


def plot_tuning_curves_days (object, n_cells_to_plot = 5):
    n_cells = object.track2p_obj.track_ops.n_tracked
    days = object.days
    cells = random.sample(range(0,n_cells), n_cells_to_plot)

    # variables and dependencies for colour mapping
    plasma = plt.get_cmap('plasma')
    cNorm  = colors.Normalize(vmin=0, vmax=object.grating_blocks)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)
    scalarMap.set_array([])

    fig, ax = plt.subplots (n_cells_to_plot, len(days)+1, figsize = (4*len(days),1.7*n_cells_to_plot), width_ratios= [3] * len(days) + [1])

    for i_cell, cell in enumerate(cells):
        max_y = np.array([object.dat_subject[day]['mean_ordered_grat_responses'][:, :, cell].mean(axis=-1).max() for day in days]).max()

        for i_day, day in enumerate(days):
            response_arr = object.dat_subject[day]['mean_ordered_grat_responses'].mean(axis = -1)
            orientations = object.dat_subject[day]['mean_ordered_grat']

            # mean response across frames for each block / orientation / cell
            [ax[i_cell, i_day].plot(orientations[block], response_arr[block, :, cell], label=f'Cell {cell + 1}', c = scalarMap.to_rgba(block)) for block in range(response_arr.shape[0])]
            ax[i_cell,i_day].set_xlabel('Orientation (deg)') if cell == cells[-1] else None
            ax[i_cell,i_day].set_ylabel(f'Response (ROI {cell})') if i_day == 0 else None
            ax[i_cell,i_day].set_ylim([0, max_y + 3])
            ax[i_cell,i_day].set_title(day) if cell == cells[0] else None
            ax[i_cell, i_day].set_xticks([]) if cell != cells[-1] else None
        ax[i_cell, -1].axis('off')

    plt.tight_layout()
    cax = fig.add_axes([0.92, 0.13, 0.02, 0.8])
    cb = plt.colorbar(scalarMap, cax = cax, ticks=range(4),orientation='vertical')
    cb.set_label('Grating Block', labelpad=10)# , y=1.05, rotation=0)
    plt.show()


def get_complex_num(obj):
    '''
    calculate tuning curves > sum [r_theta * e ^ (i2theta)] / sum [r_theta]
    theta = direction of the kth condition (array of all the directions)
    r_theta = neuronal responses to the directions

    average the responses
    '''

    for day in obj.dat_subject.keys():

        n_offsets = 2

        orientations = np.array([np.array([obj.dat_subject[day]['mean_ordered_grat'][repeat] for repeat in np.arange(start,obj.dat_subject[day]['mean_ordered_grat'].shape[0], n_offsets)]).mean(axis = 0) for start in range(n_offsets)])
        theta = np.deg2rad(np.array([val for pair in zip(orientations[0], orientations[1]) for val in pair])) # > shape: n_orientations

        obj.dat_subject[day]['theta'] = theta

        r = np.array([np.array([obj.dat_subject[day]['mean_ordered_grat_responses'][repeat] for repeat in np.arange(start,obj.dat_subject[day]['mean_ordered_grat_responses'].shape[0], n_offsets)]).mean(axis = 0) for start in range(n_offsets)])
        responses_theta = np.array([val for pair in zip(r[0], r[1]) for val in pair]).mean(axis=-1).T # mean response to the stimulus (average over frames) > shape: n_cells x n_orientations
        # np.exp outputs shape n_orientations

        obj.dat_subject[day]['tuning_curves'] = responses_theta

        #complex_nums = np.zeros(responses_theta.shape[0], dtype = np.complex64)

        # in python, j is the imaginary component
        complex_nums = (responses_theta * np.exp(1j * 2 * theta)).sum(axis=-1) / (responses_theta.sum(axis=-1))  # > shape n_cells

        obj.dat_subject[day]['complex'] = complex_nums

        # magnitude (absolute value) of the complex number > represents the amplitude of the response (how strongly neuron responds) = OSI (should be between 0 & 1)
        obj.dat_subject[day]['OSI'] = np.abs(complex_nums)

        # the argument (angle / preferred angle) of the complex number encodes the preferred orientation (the orientation to which the neuron responds most strongly).
        obj.dat_subject[day]['pref_orientation'] = np.rad2deg(np.angle(complex_nums))  # preferred angle
        # divide by 2 because of the i*theta*2 in the original equation

def add_headers(
    fig,
    *,
    row_headers=None,
    col_headers=None,
    row_pad=1,
    col_pad=5,
    rotate_row_headers=True,
    **text_kwargs
):
    # Based on https://stackoverflow.com/a/25814386

    axes = fig.get_axes()

    for ax in axes:
        sbs = ax.get_subplotspec()

        # Putting headers on cols
        if (col_headers is not None) and sbs.is_first_row():
            ax.annotate(
                col_headers[sbs.colspan.start],
                xy=(0.5, 1),
                xytext=(0, col_pad),
                xycoords="axes fraction",
                textcoords="offset points",
                ha="center",
                va="baseline",
                **text_kwargs,
            )

        # Putting headers on rows
        if (row_headers is not None) and sbs.is_first_col():
            ax.annotate(
                row_headers[sbs.rowspan.start],
                xy=(0, 0.5),
                xytext=(-ax.yaxis.labelpad - row_pad, 0),
                xycoords=ax.yaxis.label,
                textcoords="offset points",
                ha="right",
                va="center",
                rotation=rotate_row_headers * 90,
                **text_kwargs,
            )

def hist_osi_angle (obj):
    '''
    plots histogram distribution of OSIs/preferred orientation for each day
    :param obj:
    :return:
    '''

    fig, ax = plt.subplots(len(obj.dat_subject.keys()),2, figsize = (6.5,8))
    print(obj.dat_subject.keys())
    for i, day in enumerate(obj.dat_subject.keys()):

        ax[i,0].hist(obj.dat_subject[day]['grat']['preferred_orientation'], color =  'r', alpha = 0.5, label = f'Cell count, {day}', bins = np.linspace(-180,180,25))
        ax[i, 1].hist(obj.dat_subject[day]['grat']['OSI'], color =  'b', alpha = 0.5, label = f'Cell count, {day}', bins = np.linspace(0,1,25))

        # ax[i,0].set_ylim([0, 6])
        # ax[i, 0].set_xlim([-100, 100])
        # ax[i, 1].set_ylim([0, 18])
        ax[i, 0].set_xlabel('Orientation (deg)')
        ax[i, 1].set_xlabel('OSI')

    font_kwargs = dict(fontsize="large")
    add_headers(fig, col_headers=["Preferred Orientation", "Orientation Selectivity Index"], row_headers=['Cell count \n (' + day + ')' for day in list(obj.dat_subject.keys())], **font_kwargs)

    plt.tight_layout()
    plt.show()

#hist_osi_angle (suite2p_obj)

def plot_change_osi_angle(obj):
    n_cells = obj.track2p_obj.track_ops.n_tracked
    ncols = 13
    nrows = int(np.ceil(n_cells / ncols))
    fig, ax = plt.subplots(nrows, ncols, subplot_kw={'projection': 'polar'}, figsize=(18, 9))
    ax = ax.ravel()

    max_mag = np.array([np.array([obj.dat_subject[day]['OSI'][cell] for day in obj.dat_subject.keys()]) for cell in
                        range(n_cells)]).max()

    for cell in range(nrows * ncols):

        # polar plots need to be plotted in radians
        if cell < n_cells:
            r = np.array([obj.dat_subject[day]['OSI'][cell] for day in obj.dat_subject.keys()])
            theta = np.deg2rad(np.array([obj.dat_subject[day]['pref_orientation'][cell] for day in obj.dat_subject.keys()]))

            # subtract the starting angle so all cells's starting preferred angle is at 0 degrees
            #theta- theta[0]
            ax[cell].scatter(theta , r, s=20, c=np.arange(theta.shape[0]), cmap='plasma')
            ax[cell].plot(theta, r, c='black', alpha=0.5)
            ax[cell].set_thetamin(-90)
            ax[cell].set_thetamax(90)

            ax[cell].set_thetagrids([0, 90, -90], y=0.2, labels=['0', '\u03c0' + '/2', '3' + '\u03c0' + '/2'], fontsize=8)  # labels = ['0', '','\u03c0','']
            #
            # ax[cell].set_thetagrids([0, 90, 180, 270], y=0.2,
            #                         labels=['0', '\u03c0' + '/2', '\u03c0', '3' + '\u03c0' + '/2'],
            #                         fontsize=8)  # labels = ['0', '','\u03c0','']
            #ax[cell].set_rmax(max_mag)
            ax[cell].set_rlabel_position(45)
            ax[cell].tick_params(axis='y', labelsize=8)
            ax[cell].grid(True)

        if cell >= n_cells:
            ax[cell].axis('off')

    #plt.tight_layout(pad=0.9)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle('Change in OSI/pref angle across days')
    plt.show()

def tuning_curves_change (object):
    n_cells = object.track2p_obj.track_ops.n_tracked
    days = object.days
    ncols = 7
    nrows = int(np.ceil(n_cells / ncols))

    # variables and dependencies for colour mapping
    plasma = plt.get_cmap('plasma')
    cNorm  = colors.Normalize(vmin=0, vmax=object.grating_blocks)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)
    scalarMap.set_array([])

    fig, ax = plt.subplots (nrows, ncols, sharex = True, sharey = False, figsize = (19,10))
    ax = ax.ravel()

    for i_cell in range(ncols*nrows):

        if i_cell < n_cells:
            max_y = np.array([object.dat_subject[day]['mean_ordered_grat_responses'][:, :, i_cell].mean(axis=-1).max() for day in days]).max()

            responses_days = np.array([object.dat_subject[day]['tuning_curves'][i_cell] for day in days])
            orientations = np.linspace(0, 330, 12)

            # plot mean tuning curve across each day
            [ax[i_cell].plot(orientations, responses_days[i_day], c=scalarMap.to_rgba(i_day)) for i_day in range (len(days))]
            ax[i_cell].set_xticks(np.arange(0,360,90))
            ax[i_cell].tick_params(axis='x', labelsize=8)
            ax[i_cell].tick_params(axis='y', labelsize=8)

    #         ax[i_cell,i_day].set_ylim([0, max_y + 3])

        if i_cell >= n_cells:
            ax[i_cell].axis('off')
    fig.suptitle('Tuning curves across days',fontsize = 15, y = 1)
    fig.supxlabel('Orientation (deg)',fontsize = 11)
    fig.supylabel('Response', x = 0, ha = 'left',fontsize = 11)
    plt.tight_layout(pad = 0.4)
    plt.show()

def plot_curve_polar(object, roi):

    fig = plt.figure(figsize = (9,3))

    days = object.days

    # plotting line plot
    # variables and dependencies for colour mapping
    plasma = plt.get_cmap('plasma')
    cNorm  = colors.Normalize(vmin=0, vmax=len(days))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)
    scalarMap.set_array([])

    responses_days = np.array([object.dat_subject[day]['tuning_curves'][roi] for day in days])
    orientations = np.linspace(0, 330, 12)

    # plot mean tuning curve across each day
    ax1 = fig.add_subplot(121)
    p1 = [ax1.plot(orientations, responses_days[i_day], linewidth = 2, alpha = 0.8, c=scalarMap.to_rgba(i_day)) for i_day in range (len(days))]
    ax1.set_xticks(np.arange(0,360,90))
    ax1.tick_params(axis='x', labelsize=8)
    ax1.tick_params(axis='y', labelsize=8)
    ax1.set_xlabel('Orientation (deg)', fontsize = 9)
    ax1.set_ylabel('Response', fontsize = 9)
    ax1.set_title('Tuning Curves', fontsize = 11)

    # plotting polar plot
    r = np.array([object.dat_subject[day]['OSI'][roi] for day in object.dat_subject.keys()])
    theta = np.deg2rad(np.array([object.dat_subject[day]['pref_orientation'][roi] for day in object.dat_subject.keys()]))

    ax2 = fig.add_subplot(122, polar=True)
    #p2 = ax2.scatter(theta, r, s=25, c=np.arange(theta.shape[0]), cmap='plasma')
    p2 = [ax2.scatter(theta[i_day], r[i_day], s=25, color=scalarMap.to_rgba(i_day)) for i_day in range (len(days))]
    ax2.plot(theta, r, c='black', alpha=0.5)
    ax2.set_thetamin(-90)
    ax2.set_thetamax(90)
    ax2.set_thetagrids([0, 90, -90], y=0.05, labels=['0', '\u03c0' + '/2', '3' + '\u03c0' + '/2'], fontsize=8)  # labels = ['0', '','\u03c0','']
    ax2.set_rlabel_position(45)
    ax2.set_rticks(np.round(np.linspace(0, r.max(), 4),2))
    ax2.tick_params(axis='y', labelsize=8)
    ax2.set_title('Preferred Orientation \n & Magnitude ', fontsize =  10)
    ax2.grid(True)

    # ax3 = fig.add_subplot(122)
    # p3 = ax3.colorbar(scalarMap, cax = ax3, ticks=range(4),orientation='vertical')
    # #p3.set_label('Day', labelpad=10)# , y=1.05, rotation=0)
    cax = fig.add_axes([0.92, 0.13, 0.02, 0.8])
    cb = plt.colorbar(scalarMap, cax = cax, ticks=np.arange(4),orientation='vertical')
    cb.set_label('Day', labelpad=8)# , y=1.05, rotation=0)

    plt.tight_layout(pad = 0.2)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle(f'ROI # {roi}', fontsize = 15)
    plt.show()
    plt.savefig(fr'C:\Users\erica\OneDrive\Desktop\roi{roi}.png')

def diff_bw_angles(theta1, theta2):
    diff = theta2 - theta1
    if diff < -180:
        diff += 360
    elif diff > 180:
        diff -= 360
    return np.abs(diff)

def diff_bw_angles_arr(theta1_arr, theta2_arr):
    # both arrays of same shape
    diff_arr = np.zeros_like(theta1_arr)
    for i in range(theta1_arr.shape[0]):
        diff_arr[i] = diff_bw_angles(theta1_arr[i], theta2_arr[i])
    return diff_arr


def difference_across_days(obj):
    '''
    calculate the average change in preferred direction, and the average change in orientation selectivity for the whole population

    :param obj:
    :return:
    '''
    n_cells = obj.track2p_obj.track_ops.n_tracked
    days = list(obj.dat_subject.keys())

    fig, ax = plt.subplots(1, 2, figsize=(6, 3))
    osi_diff = [np.array([np.abs(obj.dat_subject[days[i]]['OSI'] - obj.dat_subject[days[i + i_apart]]['OSI']) for i in
                          range(len(days) - i_apart)]).mean(axis=-1) for i_apart in range(1, 4)]
    pd_diff = [np.array([np.abs(diff_bw_angles_arr(obj.dat_subject[days[i]]['pref_orientation'],
                                                   obj.dat_subject[days[i + i_apart]]['pref_orientation'])) for i in
                         range(len(days) - i_apart)]).mean(axis=-1) for i_apart in range(1, 4)]

    for i, day_apart in enumerate(osi_diff):
        ax[0].scatter([i + 1] * day_apart.shape[0], day_apart, c='grey', alpha=0.7, s=30)
        # scatter[0].bar([i + 1], day_apart.mean(), c='deeppink', alpha=0.7, s=50)
        ax[0].bar([i + 1], day_apart.mean(), color='deeppink', alpha=0.55)
        ax[0].set_ylabel(r'$\Delta$' + ' OSI')
        ax[0].set_xlabel(r'$\Delta$' + ' time (days)')

    for i, day_apart in enumerate(pd_diff):
        ax[1].scatter([i + 1] * day_apart.shape[0], day_apart, c='grey', alpha=0.7, s=30)
        # ax[1].scatter([i + 1], day_apart.mean(), c='darkviolet', alpha=0.7, s=50)
        ax[1].bar([i + 1], day_apart.mean(), color='darkviolet', alpha=0.55)
        ax[1].set_ylabel(r'$\Delta$' + ' PD')
        ax[1].set_xlabel(r'$\Delta$' + ' time (days)')
    plt.tight_layout()
    plt.show()



def plot_cell_diff (obj, pd_or_osi):
    '''
    plot the OSI/PD of each ROI on day X against the OSI/PD of each ROI on day X+1 or X+2 or X+3
    :param obj:
    :param pd_or_osi (str): 'pref_orientation' or 'OSI
    :return:
    '''
    fig, ax = plt.subplots(4, 4, figsize=(10, 10))

    if pd_or_osi == 'pref_orientation':
        pr_dir = np.array([obj.dat_subject[day]['pref_orientation'] for day in obj.dat_subject.keys()])

        for i_day, iday in enumerate(obj.dat_subject.keys()):
            for k_day, kday in enumerate(obj.dat_subject.keys()):

                if i_day != k_day:
                    osi_i, osi_k = obj.dat_subject[iday]['pref_orientation'], obj.dat_subject[kday]['pref_orientation']

                    buffer = 0
                    ax[i_day, k_day].scatter(osi_i, osi_k, s = 25, alpha = 0.6, c = 'deeppink')
                    # the r value is the correlation coefficient
                    corr = np.round(pearsonr(osi_i, osi_k)[0], 2)
                    ax[i_day, k_day].text(-70,70, f'r = {corr}')
                    ax[i_day, k_day].set_title((i_day, k_day))
                    ax[i_day, k_day].set_xlabel(f'Angle (deg) (day {i_day})')
                    ax[i_day, k_day].set_title(r'$\Delta$' + ' Preferred Orientation', fontsize = 10)
                    ax[i_day, k_day].set_ylabel(f'Angle (deg) (day {k_day})')
                    ax[i_day,k_day].plot ([0-buffer, pr_dir.max()+buffer],[0-buffer, pr_dir.max()+buffer],'--', alpha = 0.6, zorder = 0, c = 'grey')
                    ax[i_day, k_day].set_box_aspect(1)
                else:
                    ax[i_day,k_day].axis('off')

    elif pd_or_osi == 'OSI':
        osis = np.array([obj.dat_subject[day]['OSI'] for day in obj.dat_subject.keys()])

        for i_day, iday in enumerate(obj.dat_subject.keys()):
            for k_day, kday in enumerate(obj.dat_subject.keys()):

                if i_day != k_day:
                    osi_i, osi_k = obj.dat_subject[iday]['OSI'], obj.dat_subject[kday]['OSI']

                    ax[i_day, k_day].scatter(osi_i, osi_k, s=25, alpha=0.6, c='darkviolet')
                    # the r value is the correlation coefficient
                    corr = np.round(pearsonr(osi_i, osi_k)[0], 2)
                    ax[i_day, k_day].text(0.05, 0.6, f'r = {corr}')
                    ax[i_day, k_day].set_title((i_day, k_day))
                    ax[i_day, k_day].set_xlabel(f'OSI (day {i_day})')
                    ax[i_day, k_day].set_title(r'$\Delta$' + ' Orientation Selectivity Index', fontsize=10)
                    ax[i_day, k_day].set_ylabel(f'OSI (day {k_day})')
                    ax[i_day, k_day].set_xlim([0 - osis.max() / 12, osis.max() + osis.max() / 12])
                    ax[i_day, k_day].set_ylim([0 - osis.max() / 12, osis.max() + osis.max() / 12])
                    ax[i_day, k_day].plot([0 - osis.max() / 12, osis.max() + osis.max() / 12],
                                          [0 - osis.max() / 12, osis.max() + osis.max() / 12], '--', alpha=0.6,
                                          zorder=0, c='grey')
                    ax[i_day, k_day].set_box_aspect(1)
                else:
                    ax[i_day, k_day].axis('off')

    plt.tight_layout()
    plt.show()


################################



def repair_opsfile(suite2p_path_list):
    '''
    Modifies each 'ops.py' file in 'suite2p_path_list'
    Saves the fixed file as 'ops.npy' and saves the old version as 'ops old.npy'

    Running cellReg in suite2p (anatomical_detection = 3) outputs a buggy  ops.npy file, with a Zero-Division Error
    After getting suite2p outputs, run this function to fix these bugs, then we can select/refine the ROI selection

    :param suite2p_path_list: list of paths to the 'plane0' file where each ops.npy file is located
    :return: Void
    '''

    for path in suite2p_path_list:
        ops = np.load(os.path.join(path, 'ops.npy'), allow_pickle=True).item()
        np.save(os.path.join(path, 'ops old.npy'), ops, allow_pickle=True)  #save original ops file as 'ops old. npy'
        ops['diameter'][0], ops['diameter'][1] = 1, 1                       # change the buggy values in the file (zero division error)
        np.save(os.path.join(path, 'ops.npy'), ops, allow_pickle=True)      #save new ops file as 'ops.npy'

def masked_spatial_footprints(ops, stat):
    '''
    From the suite2p output, creates an array of shape (n_ROIs x height x width)
    such that each ROI represents a spatial mask with the spatial profile of that ROI/cell

    To be able to utilize the CellReg (ZivLab) toolbox that allows tracking of the same cells across days

    :param ops: ops.npy file
    :param stat: stat.npy file
    :return: array of shape (n_ROIs x height x width) with the spatial profile of each ROI/cell
    '''
    n_rois = stat.shape[0]

    # create set of ROIs found, as a spatial mask of shape (n_rois x height x width)
    height, width = ops['Ly'], ops['Lx']  # Image size in pixels
    spatial_footprint_array = np.zeros((n_rois, height, width))
    for n in range(n_rois):
        non_overlapping_indices = np.argwhere(~stat[n]['overlap'])[:,
                                  0]  # indices of the pixels in the ROI that aren't overlapping with other cells
        ypix, xpix = stat[n]['ypix'][~stat[n]['overlap']], stat[n]['xpix'][
            ~stat[n]['overlap']]  # y and x pixels of that cell/ROI
        spatial_footprint_array[n, ypix, xpix] = stat[n]['lam'][non_overlapping_indices]

    return spatial_footprint_array

#def generate_sf (animal_name, date, session, stat_file, ops_file, im_path):
def generate_sf (animal_name, date, stat_file, ops_file):
    im = np.zeros((ops_file['Ly'], ops_file['Lx']))
    ncells = stat_file.shape[0]
    for n in range(ncells):
        ypix = stat_file[n]['ypix'][~stat_file[n]['overlap']] # y-pixels of cell that are non-overlapping with other cells
        xpix = stat_file[n]['xpix'][~stat_file[n]['overlap']] # x-pixels of cell that are non-overlapping with other cells
        im[ypix,xpix] = n+1
    im[im == 0] = np.nan                    # make zero-values into nans (for the purpose of the visualization)

    fig, ax = plt.subplots(figsize = (7,7))
    ax.imshow(im, cmap = 'plasma')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.title(f'Spatial Footprint ({animal_name}: {date})')
    #save_path = os.path.split(im_path)[0]
    #plt.savefig (os.path.join(save_path,'Spatial Footprints', 'Images', f'Spatial Footprint ({animal_name}: {date}, #{session})'))
    plt.show()
    #plt.close()

def suite2p_files (suite2p_path, response_type = 'deconvolved'):
    '''
    :param response_type: 'deconvolved' or 'fluorescence' - specifies if we want deconvolved (binarized) traces or raw fluorescnece values

    :param suite2p_path: path to the suite2p plane0 folder
    :return: neuronal responses file, is cell file, ops file, stat file

    - F.npy > array of fluorescence traces (ROIs by timepoints)
    - Fneu.npy > array of neuropil fluorescence traces (ROIs by timepoints)
    - spks.npy > array of deconvolved traces (ROIs by timepoints)
    - stat.npy > list of statistics computed for each cell (ROIs by 1)
    - ops.npy > options and intermediate outputs (dictionary)
    -
    '''

    if response_type == 'deconvolved': # spks.npy > array of deconvolved traces (ROIs by timepoints) (generally cleaner)
        responses_file = np.load(os.path.join(suite2p_path, 'spks.npy'), allow_pickle=True)  # array of fluorescence traces (cells x timepoints)
    elif response_type == 'fluorescence': # F.npy > array of fluorescence traces (ROIs by timepoints)
        F = np.load(os.path.join(suite2p_path, 'F.npy'),allow_pickle=True)  # array of fluorescence traces (cells x timepoints)
        Fneu = np.load(os.path.join(suite2p_path, 'Fneu.npy'), allow_pickle=True)  # array of fluorescence traces (cells x timepoints)

        responses_file = F - (0.7 * Fneu)  # array of fluorescence traces (cells x timepoints)

    #iscell.npy > specifies whether an ROI is a cell, first column is 0/1, and second column is probability that the ROI is a cell based on the default classifier
    iscell_file = np.load(os.path.join(suite2p_path, 'iscell.npy'), allow_pickle=True)[:,0]  # specifies whether an ROI is a cell, either 0/1

    # ops.npy > options and intermediate outputs (dictionary)
    ops_file = np.load(os.path.join(suite2p_path, 'ops.npy'), allow_pickle=True).item()

    #stat.npy > list of statistics computed for each cell (ROIs by 1)
    stat_file = np.load(os.path.join(suite2p_path, 'stat.npy'), allow_pickle=True)

    return responses_file, iscell_file, ops_file, stat_file


# def get_sbx_suite2p_path (path_to_data, animal_name, date_yymmdd):
#     '''
#     :param path_to_data: ex: '/Volumes/Erica1/TrenholmLab/RepresentationalDrift/Data'
#     :param animal_name: ex: 'EC_GECO_11'
#     :param date_yymmdd: ex: '20231206'
#     :return: returns the path to the scanbox file. to view a scanbox file, in the terminal, type f'sbxv
#

def save_scanbox_img (path_to_img, file_name):
    '''
    Saves a .mat imaging file (file_name) located in path_to_image to .png
    '''

    mat = scipy.io.loadmat(os.path.join(path_to_img, file_name))['img']
    plt.imsave(os.path.join(path_to_img, os.path.splitext(file_name)[0] + '.png'), mat, cmap='gray')


# save_scanbox_img (r'I:\dark_drift\data\EC_dark_03\20250526\spon0\spon0_000_000\location', 'wdf_26.37deg.mat')

def save_scanbox_images (path_to_folder):
    '''
    Converts all .mat files taken during an imaging session (located in path_to_folder) to .png
    '''

    to_convert = [f for f in os.listdir(path_to_folder)
                  if f.lower().endswith('.mat') and os.path.isfile(os.path.join(path_to_folder, f))]

    for file in to_convert:
        save_scanbox_img(path_to_folder, file)


#save_scanbox_images (r'I:\dark_drift\data\EC_dark_03\20250526\spon0\spon0_000_000\location')

############## spontaneous analysis helpers

def compute_coactivation_matrix(data_trials_all, valid_epochs, response_threshold):
    '''
    :param data_trials_all:
    :param valid_epochs:
    :param response_threshold:
    :return: symmetric coactivation_matrix, where high values mean two neurons often fire together across epochs.
    '''
    n_cells = data_trials_all.shape[1]
    coactivation_matrix = np.zeros((n_cells, n_cells), dtype=int)

    for epoch in data_trials_all[valid_epochs]:  # shape: (n_cells, n_timepoints)
        active_cells = np.where(np.any(epoch > response_threshold, axis=1))[0]

        # Count self-activation
        for i in active_cells:
            coactivation_matrix[i, i] += 1

        # Count co-activation
        for i, j in combinations(active_cells, 2):
            coactivation_matrix[i, j] += 1
            coactivation_matrix[j, i] += 1

    return coactivation_matrix

def perform_pca (data_trials_normalized, n_components = 10):
    '''
    :param data_trials_normalized: array shape (n_trials, n_features), each trial (each row) is normalized to unit length (l2 normalization)
    '''
    # performing PCA & projection
    pca = PCA(n_components=n_components, svd_solver='full') # full forces the lapack solver
    principal_components = pca.fit_transform(data_trials_normalized)
    explained_variance_ratio = pca.explained_variance_ratio_ # eigenvalues
    return pca, principal_components, explained_variance_ratio

def filter_active_epochs (data_trials_array, response_threshold = 3, cell_threshold = 4):
    '''
    Remove epochs with no activity or only 1-2 cells active
    :param data_trials_array: numpy array of shape (n_cells, n_chunks)
    :param threshold: z-score threshold that cells need to be to be considered active
    :return:
    '''
    # check whether each cell goes above response threshold in each epoch
    # then count number of cells (per epoch) that exceed response threshold (shape n_epochs)
    cells_above_threshold = np.any(data_trials_array > response_threshold, axis=-1).sum(axis = -1)

    # valid epochs are only those with > 2 cells co-active (shape n_epochs)
    valid_epochs = cells_above_threshold >= cell_threshold

    return valid_epochs

def active_epochs_idx (epochs_array, response_threshold = 3, min_active_cells = 2):
    '''
    Return indices of epochs that have some mean level of activity within the epoch
    :param epochs_array: numpy array of shape (n_cells, n_epochs, n_timepoints)
    :param threshold: z-score threshold that cells need to be to be considered active
    :return:
    '''
    # consider the mean activity across the epoch
    # boolean checks that each cell's mean activity in an epoch goes over response_threshold > shape (n_cells, n_epochs)
    # then we check, across each epoch, that there is at min_active_cells that is over the threshold
    epoch_mask = np.sum(epochs_array.mean(axis=-1) > response_threshold, axis=0) >= min_active_cells

    # return the boolean (T/F) mask with T indices that pass this test
    return epoch_mask



def decay_eigenspectra(eigenvalues):
    '''
    Fit power-law in log-log space

    Power-law: y = Ax**-alpha, where x is the rank, alpha is the slope (rate of decay), A is constant (intercept in log-log space)
    - to solve for alpha > take log of both sides
    - > log(y) = log(Ax**-alpha)
    - > log(y) = log(A) + log(x**-alpha)
    - > log(y) = log(A) -alpha*log(x)
    - > log(y) ~= -alpha*log(x)
    - > log(y)/log(x) ~= -alpha

    :param eigenvalues: variance explained output from PCA
    :return: -slope
    '''
    x = np.arange(1, len(eigenvalues) + 1)
    log_x = np.log10(x)
    log_y = np.log10(eigenvalues + 1e-10)  # add epsilon to avoid log(0)
    slope, intercept, r_value, _, _ = linregress(log_x, log_y)
    return -slope

def binarize_array(data, threshold=4):
    return (data >= threshold).astype(int)

def sort_key(label):
    num = int(re.findall(r'\d+', label)[0])
    return (num, 1 if '_s' in label else 0)

def plot_slope_eigenvals(power_law_slopes):
    ncols = 2
    fig, ax = plt.subplots(
        nrows=int(np.ceil(len(power_law_slopes)/ncols)),
        ncols=ncols, figsize=(12, 7), sharey=True
    )
    ax = ax.ravel()
    colours = {'EB': 'blue', 'EBc': 'black', 'LB': 'red', 'LBc':'green'}

    max_days = max(
        len(days)
        for group_dict in power_law_slopes.values()
        for days in group_dict.values())

    if len(power_law_slopes) == 1:
        ax = [ax]

    for i_group, group in enumerate(power_law_slopes):
        group_days = sorted(
            {d for days in power_law_slopes[group].values() for d in days.keys()},
            key=sort_key
        )

        ax[i_group].plot(group_days, [np.nan]*len(group_days), alpha=0)

        for animal in power_law_slopes[group]:
            sorted_items = sorted(power_law_slopes[group][animal].items(),
                                  key=lambda x: sort_key(x[0]))
            x = [k for k, _ in sorted_items]
            y = [v for _, v in sorted_items]
            ax[i_group].plot(x, y, marker='o', color='gray', alpha=0.4)

        # group mean ± SEM
        values_by_day = {}
        for animal, days in power_law_slopes[group].items():
            for day, val in days.items():
                values_by_day.setdefault(day, []).append(val)

        mean_vals = {day: np.mean(vals) for day, vals in values_by_day.items()}
        sem_vals  = {day: np.std(vals, ddof=1)/np.sqrt(len(vals)) for day, vals in values_by_day.items()}

        x_mean = group_days
        y_mean = [mean_vals[d] for d in group_days if d in mean_vals]
        y_sem  = [sem_vals[d]  for d in group_days if d in sem_vals]

        ax[i_group].errorbar(x_mean[:len(y_mean)], y_mean, yerr=y_sem, fmt='d-', capsize=4, alpha = 0.6,
                             color=colours.get(group, 'gray'), markersize=8, label=group)

        ax[i_group].set_xlabel('Recording session')
        ax[i_group].set_ylabel('Power-law slope')
        ax[i_group].set_title(f'{group} group')
        ax[i_group].grid(axis='y', linestyle='--', alpha=0.5)
        ax[i_group].tick_params(axis='x', rotation=45)
        ax[i_group].set_xlim([-0.5, max_days])

    plt.tight_layout()

    # folder_path = os.path.join(obj.dat[animal].save_path, f'slope eigenvals')
    # if not os.path.exists(folder_path):
    #     os.makedirs(folder_path)
    # plt.savefig(os.path.join(folder_path, f'slope eigenvals.svg'))
    # plt.savefig(os.path.join(folder_path, f'slope eigenvals.png'))


    plt.show()

def calculate_animal_age(dob, imaging_date):
    '''
    Returns the age of the animal born on DOB, on imaging_date (in post-natal days, and weeks)
    :param dob: format 'YYYY/MM/DD'
    :param imaging_date: 'YYYY/MM/DD'
    :return:
    '''
    # Convert input strings to date objects
    dob = datetime.strptime(dob, "%Y%m%d")
    imaging_date = datetime.strptime(imaging_date, "%Y%m%d")

    # Calculate the difference between the dates in days
    age_in_days = (imaging_date - dob).days

    # Calculate weeks and remaining days
    weeks = age_in_days // 7
    days = age_in_days % 7

    return f"P{age_in_days}"

def calculate_day_animal_age(dob, target_postnatal_day):
    '''
    Returns day that the animal born on DOB will be of age target_postnatal_day
    :param dob: format 'YYYYMMDD'
    :param target_postnatal_day: INT
    :return:
    '''
    # Convert input strings to date objects
    dob = datetime.strptime(dob, "%Y%m%d")
    formatted_dob = dob.strftime("%Y%m%d")

    dob_plus_target = dob + timedelta(days=target_postnatal_day)
    formatted_date = dob_plus_target.strftime("%Y%m%d")

    print( f"Animal born on {formatted_dob} will be P{target_postnatal_day} on {formatted_date})")


def calculate_day_animal_dob(target_postnatal_day, target_recording_day):
    '''
    Returns DOB of an animal that we want to be of age target_postnatal_day by target_recording_day
    :param target_recording_day (str): format 'YYYYMMDD'
    :param target_postnatal_day (int):
    :return: None (prints string)
    '''
    # Convert input strings to date objects
    target_recording_day = datetime.strptime(target_recording_day, "%Y%m%d")
    target_day_minus_age = target_recording_day - timedelta(days=target_postnatal_day)
    formatted_dob = target_day_minus_age.strftime("%Y%m%d")

    print( f"Animal born on {formatted_dob} will be P{target_postnatal_day} on {target_recording_day})")

#calculate_day_animal_dob(91,'20250916') # want animal to be P91 on '2025/09/16'

def variance_explained(var_explained_dict, log=False, timepoints = False):
    """
    Plot variance explained of different groups in 'var_explained_dict'.
    If log=True, plot log10(PC index) vs log10(variance explained).
    """

    plt.figure(figsize=(8, 5))
    if timepoints:
        colors = {'Control_0': 'cornflowerblue', 'Control_1': 'darkblue', 'RD1_0': 'salmon', 'RD1_1': 'firebrick'}
        x = np.arange(1, var_explained_dict['Control_0'].shape[1] + 1)
    else:
        colors = {'Control': 'black', 'RD1': 'red', 'GNAT': 'green'}
        x = np.arange(1, var_explained_dict['Control'].shape[1] + 1)
    x_plot = np.log10(x) if log else x

    for group, data in var_explained_dict.items():
        data = np.array(data)
        if log:
            data = np.log10(data + 1e-10)
        mean = data.mean(axis=0)
        error = sem(data, axis=0)

        c = colors.get(group, 'gray')
        plt.plot(x_plot, data.T, alpha=0.3, c=c)
        plt.plot(x_plot, mean, alpha=0.8, c=c, linewidth=2.5, label=group)
        plt.fill_between(x_plot, mean - error, mean + error, color=c, alpha=0.2)

    plt.legend(fontsize=14)
    plt.xlabel('log10(PC index)' if log else 'PC index')
    plt.ylabel('log10(variance explained)' if log else 'Variance explained')
    plt.title('Log-Log Variance Decay' if log else 'Variance explained by PCs')
    #plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()


def compute_event_features_per_cell(data, z_threshold=3, fps=10):
    """
    Compute spontaneous event rates, amplitudes, and durations from z-scored traces.

    Parameters:
    - data: array of shape (n_cells, n_timepoints)
    - z_threshold: threshold to detect peaks
    - fps: imaging frame rate (frames per second)

    Returns:
    - rates: array of shape (n_cells,) with event rates in Hz
    - all_amplitudes: list of all peak amplitudes across cells
    - all_durations: list of all event durations (in seconds) across cells
    """
    n_cells, n_timepoints = data.shape
    duration = n_timepoints / fps  # in seconds
    rates = np.zeros(n_cells)
    all_amplitudes = []
    all_durations = []

    for cell_idx in range(n_cells):
        trace = data[cell_idx, :]
        peaks, properties = find_peaks(trace, height=z_threshold)
        rates[cell_idx] = len(peaks) / duration

        # Collect amplitudes (z-scored peak heights)
        all_amplitudes.extend(properties["peak_heights"])

        # Collect durations (half-max width in seconds)
        widths = peak_widths(trace, peaks, rel_height=0.5)[0] / fps
        all_durations.extend(widths)

    return rates, all_amplitudes, all_durations


def plot_spon_event_rate(event_rates_groups, timepoints=False):
    '''
    Plotting spontaneous activity rates
    :param event_rates_groups: dictionary with group data
    :param timepoints:
    :return:
    '''
    plt.figure(figsize=(5, 4))

    if timepoints:
        colors = {'Control_0': 'cornflowerblue', 'Control_1': 'darkblue', 'RD1_0': 'salmon', 'RD1_1': 'firebrick'}
    else:
        colors = {'Control': 'black', 'RD1': 'red', 'GNAT': 'green'}

    if timepoints:
        event_rates_groups = {k: v for k, v in event_rates_groups.items() if 'RD1' in k}

    all_counts = np.concatenate([event_rates_groups[g] for g in event_rates_groups])
    bins = np.linspace(all_counts.min(), all_counts.max(), 30)  # 45 bins

    for g in event_rates_groups.keys():
        plt.hist(event_rates_groups[g], bins = bins,histtype ='step', linewidth = 2.5, density = True, label = g, color = colors[g], alpha = 1)
    plt.legend()
    plt.title('Spontaneous activity rates')
    plt.xlabel('Spontaneous Event rates (Hz)')
    plt.ylabel('Probability density')
    plt.show()

    # Optional: KS test between first two groups
    group1, group2 = list(event_rates_groups.keys())[:2]
    stat, pval = ks_2samp(event_rates_groups[group1], event_rates_groups[group2])
    print(f"KS test between {group1} and {group2}: D = {stat:.3f}, p = {pval:.3e}")

    stat, pval = mannwhitneyu(event_rates_groups[group1], event_rates_groups[group2],alternative='two-sided')
    print(f"Mann-Whitney U test: U = {stat:.3f}, p = {pval:.3e}")

def plot_spon_event_properties(amplitudes, durations, timepoints=False):
    '''
    Plotting spontaneous activity rates
    :param event_rates_groups: dictionary with group data
    :param timepoints:
    :return:
    '''
    fig, ax = plt.subplots (1,2,figsize=(10, 5))

    if timepoints:
        colors = {'Control_0': 'cornflowerblue', 'Control_1': 'darkblue', 'RD1_0': 'salmon', 'RD1_1': 'firebrick'}
    else:
        colors = {'Control': 'black', 'RD1': 'red', 'GNAT': 'green'}

    if timepoints: # only plot the rd1 data
        amplitudes = {k: v for k, v in amplitudes.items() if 'RD1' in k}
        durations = {k: v for k, v in durations.items() if 'RD1' in k}

    # amplitudes:
    all_counts = np.concatenate([amplitudes[g] for g in amplitudes])
    bins = np.linspace(all_counts.min(), 15, 30)  #all_counts.max()
    for g in amplitudes.keys():
        ax[0].hist(amplitudes[g], bins = bins,histtype ='step', linewidth = 2.5, density = True, label = g, color = colors[g], alpha = 0.7)
    ax[0].legend()
    ax[0].set_title('Spontaneous event amplitudes')
    ax[0].set_xlabel('Spontaneous Event amplitudes (z-score)')
    ax[0].set_ylabel('Probability density')
    #plt.show()
    # stats between first two groups
    group1, group2 = list(amplitudes.keys())[:2]
    stat, pval = ks_2samp(amplitudes[group1], amplitudes[group2])
    print(f"AMPLITUDE: KS test between {group1} and {group2}: D = {stat:.3f}, p = {pval:.3e}")
    stat, pval = mannwhitneyu(amplitudes[group1], amplitudes[group2],alternative='two-sided')
    print(f"AMPLITUDE: Mann-Whitney U test: U = {stat:.3f}, p = {pval:.3e}")

    # event durations:
    all_counts = np.concatenate([durations[g] for g in durations])
    bins = np.linspace(all_counts.min(), 7, 30)  # all_counts.max()
    for g in durations.keys():
        ax[1].hist(durations[g], bins = bins,histtype ='step', linewidth = 2.5, density = True, label = g, color = colors[g], alpha = 0.7)
    ax[1].legend()
    ax[1].set_title('Spontaneous event duration')
    ax[1].set_xlabel('Spontaneous Event duration')
    ax[1].set_ylabel('Probability density')

    # stats between first two groups
    group1, group2 = list(durations.keys())[:2]
    stat, pval = ks_2samp(durations[group1], durations[group2])
    print(f"DURATION: KS test between {group1} and {group2}: D = {stat:.3f}, p = {pval:.3e}")
    stat, pval = mannwhitneyu(durations[group1], durations[group2],alternative='two-sided')
    print(f"DURATION: Mann-Whitney U test: U = {stat:.3f}, p = {pval:.3e}")

    plt.show()

def plot_fourrier(dom_frequs_dict, power_spec, timepoints=False):
    '''
    Plotting spontaneous activity rates
    :param dom_frequs_dict: dictionary with group data
    :param timepoints:
    :return:
    '''
    fig, ax = plt.subplots (1,2,figsize=(10, 5))

    if timepoints:
        colors = {'Control_0': 'cornflowerblue', 'Control_1': 'darkblue', 'RD1_0': 'salmon', 'RD1_1': 'firebrick'}
    else:
        colors = {'Control': 'black', 'RD1': 'red', 'GNAT': 'green'}

    if timepoints: # only plot the rd1 data
        dom_frequs_dict = {k: v for k, v in dom_frequs_dict.items() if 'RD1' in k}
        power_spec = {k: v for k, v in power_spec.items() if 'RD1' in k}

    # dom freq:
    all_counts = np.concatenate([dom_frequs_dict[g] for g in dom_frequs_dict])
    #bins = np.linspace(all_counts.min(), 15, 30)  #all_counts.max()
    bins = np.linspace(all_counts.min(), all_counts.max(), 30)  # all_counts.max()
    for g in dom_frequs_dict.keys():
        ax[0].hist(dom_frequs_dict[g], bins = bins,histtype ='step', linewidth = 2.5, density = True, label = g, color = colors[g], alpha = 0.7)
    ax[0].legend()
    ax[0].set_title('Fourier')
    ax[0].set_xlabel('Dominant frequency per cell')
    ax[0].set_ylabel('Probability density')
    #plt.show()
    # stats between first two groups
    group1, group2 = list(dom_frequs_dict.keys())[:2]
    stat, pval = ks_2samp(dom_frequs_dict[group1], dom_frequs_dict[group2])
    print(f"dom_frequs_dict: KS test between {group1} and {group2}: D = {stat:.3f}, p = {pval:.3e}")
    stat, pval = mannwhitneyu(dom_frequs_dict[group1], dom_frequs_dict[group2],alternative='two-sided')
    print(f"dom_frequs_dict: Mann-Whitney U test: U = {stat:.3f}, p = {pval:.3e}")

    # power spectrum
    all_counts = np.concatenate([power_spec[g] for g in power_spec])
    #bins = np.linspace(all_counts.min(), 15, 30)  #all_counts.max()
    bins = np.linspace(all_counts.min(), data_object.fps/2, 30)  # all_counts.max()
    for g in power_spec.keys():
        plt.plot(fft_freqs, power.mean(axis=0))

        ax[1].hist(power_spec[g], bins = bins,histtype ='step', linewidth = 2.5, density = True, label = g, color = colors[g], alpha = 0.7)
    ax[1].legend()
    ax[1].set_title('Fourier')
    ax[1].set_xlabel('Power spectrum (av. across cells)')
    ax[1].set_ylabel('Probability density')
    #plt.show()
    # stats between first two groups
    group1, group2 = list(power_spec.keys())[:2]
    stat, pval = ks_2samp(power_spec[group1], power_spec[group2])
    print(f"dom_frequs_dict: KS test between {group1} and {group2}: D = {stat:.3f}, p = {pval:.3e}")
    stat, pval = mannwhitneyu(power_spec[group1], power_spec[group2],alternative='two-sided')
    print(f"dom_frequs_dict: Mann-Whitney U test: U = {stat:.3f}, p = {pval:.3e}")
    plt.show()

def coactive_cells_per_frame(coactive_counts_groups, timepoints=False):
    '''
    Plotting spontaneous activity rates
    :param event_rates_groups: dictionary with group data
    :param timepoints:
    :return:
    '''
    plt.figure(figsize=(5, 4))

    if timepoints:
        colors = {'Control_0': 'cornflowerblue', 'Control_1': 'darkblue', 'RD1_0': 'salmon', 'RD1_1': 'firebrick'}
    else:
        colors = {'Control': 'black', 'RD1': 'red', 'GNAT': 'green'}

    if timepoints:
        coactive_counts_groups = {k: v for k, v in coactive_counts_groups.items() if 'RD1' in k}

        all_counts = np.concatenate([coactive_counts_groups[g] for g in coactive_counts_groups])
        #all_counts = np.concatenate([coactive_counts_groups[g] for g in [g for g in coactive_counts.keys() if 'Control' not in g]])
        bins = np.linspace(all_counts.min(), all_counts.max(), 50)  # 45 bins

        for g in coactive_counts_groups.keys(): #[g for g in coactive_counts.keys() if 'Control' not in g]:
            plt.hist(coactive_counts_groups[g], bins=bins,histtype ='step', linewidth = 2.5, density=True, label=g, color=colors[g], alpha=1)
    else:
        all_counts = np.concatenate([coactive_counts_groups[g] for g in coactive_counts_groups])
        bins = np.linspace(all_counts.min(), all_counts.max(), 110)  # 45 bins

        for g in coactive_counts_groups.keys():
            plt.hist(coactive_counts_groups[g], bins = bins, histtype ='step', linewidth = 2.5,density = True, label = g, color = colors[g], alpha =1)
    plt.legend()
    plt.title('Co-active cells per frame')
    plt.xlim([0,30])
    plt.xlabel('Number of co-active cells')
    plt.ylabel('Probability density')
    plt.show()

    # Optional: KS test between first two groups
    group1, group2 = list(coactive_counts_groups.keys())[:2]
    stat, pval = ks_2samp(coactive_counts_groups[group1], coactive_counts_groups[group2])
    print(f"KS test between {group1} and {group2}: D = {stat:.3f}, p = {pval:.3e}")

    stat, pval = mannwhitneyu(coactive_counts_groups[group1], coactive_counts_groups[group2],alternative='two-sided')
    print(f"Mann-Whitney U test: U = {stat:.3f}, p = {pval:.3e}")


def get_earliest_file_time(folder, use_mtime=True):
    """
    Return the earliest timestamp (creation or modification) among files inside the folder.
    If the folder is empty or contains no files, return a large value to sort it last.
    """
    timestamps = []
    for root, _, files in os.walk(folder):
        for f in files:
            fpath = os.path.join(root, f)
            try:
                time = os.path.getmtime(fpath) if use_mtime else os.path.getctime(fpath)
                timestamps.append(time)
            except Exception:
                pass  # skip problematic files
    return min(timestamps) if timestamps else float('inf')


def get_paths(stim, root_dir_animal):
    '''
    :param stim: 'spon' or 'grat'
    :param root_dir_animal: example, f'I:\dark_drift\data\{animal}'
    :return:

    Get list of directory paths, up until experiments folder
    '''
    # root_dir = f'I:\dark_drift\data\{animal}'

    files = []
    for dirpath, dirnames, filenames in os.walk(root_dir_animal):
        if os.path.basename(dirpath) == 'experiments':
            if stim in dirpath.split('\\')[5]:
                files.append(dirpath)
    files.sort(key=lambda f: get_earliest_file_time(f, use_mtime=True))

    return files


def spontaneous_analysis_slope(obj, rec_type = 'spon', coactive_epochs=True, chunk_size=3):
    '''
    spontaneous analysis: slope of eigenspectrum across time
    :param obj:
    :param coactive_epochs:
    :param chunk_size:
    :return:
    '''

    groups = np.unique([a.split('_')[1] for a in obj.list_animals])

    # {'EB': {animal 1: [], animal 2:[]}],'EBc': {animal 1: [] ...}, 'LB':{animal 1: [] ...}, 'LBc': {animal 1: [] ...}}
    var_explained = {g: {} for g in groups}
    eigenvals_slope = {g: {} for g in groups}

    for i_animal, animal in enumerate(obj.dat):

        group = animal.split('_')[1]
        days_recordings = [(day, subfile) for day in obj.dat[animal].dat_subject for subfile in
                           obj.dat[animal].dat_subject[day] if rec_type in subfile]

        vexpl = {}
        eslope = {}

        for i, (day, subfile) in enumerate(days_recordings):

            timepoint = calculate_animal_age(obj.animal_dobs[animal], day) #+ '_s' * int(subfile[-1])
            if rec_type=='spon' and  (timepoint == 'P70' or timepoint == 'P105'): # if animals are seeing for the firs time (P70 or P105), add the _s for spon1
                timepoint += '_s' * int(subfile[-1])

            # shape n_Cells, n_timepoints
            spon_arr = obj.dat[animal].dat_subject[day][subfile]['zscored_traces']

            if coactive_epochs:
                n_full_chunks = spon_arr.shape[1] // chunk_size

                # shape > (n_trials, n_features, n_timepoints_per_epoch) (n_epochs, n_cells, n_timepoints)
                # data_trials = np.array([spon_arr[:,i:i + chunk_size] for i in range(0, n_full_chunks * chunk_size, chunk_size)]).reshape (n_full_chunks, -1)
                data_trials_all = np.array(
                    [spon_arr[:, i:i + chunk_size] for i in range(0, n_full_chunks * chunk_size, chunk_size)])

                # only take epochs that have at least 'cell_threshold' co-active cells that each go above 'response_threshold' > (n_epochs, n_cells, n_timepoints)
                valid_epochs = filter_active_epochs(data_trials_all, response_threshold=3, cell_threshold=5)
                # plot_raster(data_trials_all, valid_epochs, n_trials_to_plot=60, n_cells_to_plot=400, threshold=4)

                # then average over each epoch (time) to get average response> (n_epochs, n_cells)
                data_trials = data_trials_all[valid_epochs].mean(axis=-1)

            else:
                data_trials = spon_arr.T
            # Normalize each trial's activity pattern (vector) to unit length (L2 norm)
            # each trial (each row) is normalized to unit length (its L2 norm is 1)
            norms = np.linalg.norm(data_trials, axis=1, keepdims=True)
            data_trials_normalized = data_trials / norms

            # performing PCA & projection
            # explained variance ratio > eigenvalues
            pca, principal_components, explained_variance_ratio = perform_pca(data_trials_normalized, n_components=50)

            # vexpl.append(explained_variance_ratio)
            # eslope.append(decay_eigenspectra(explained_variance_ratio))
            #

            vexpl[timepoint] = explained_variance_ratio
            eslope[timepoint] = decay_eigenspectra(explained_variance_ratio)

        var_explained [group] [animal] = vexpl
        eigenvals_slope [group] [animal] = eslope

    return var_explained, eigenvals_slope

def load_ops_meanimg(ops, type):
    """Load a mean image from ops.npy with several fallbacks.
    type: meanImg, meanImgE"""

    if isinstance(ops, np.ndarray) and ops.shape == ():  # saved as a dict via np.save
        ops = ops.item()
    elif isinstance(ops, np.ndarray) and ops.size == 1:
        ops = ops.ravel()[0].item()

    # Common keys Suite2p writes
    # for k in ("meanImg", "meanImgE", "refImg", "refImgE", "meanImg_chan2"):
    k = type
    if k in ops and ops[k] is not None and np.size(ops[k]) > 0:
        img = np.asarray(ops[k])
        if img.ndim == 3:
            # pick first channel/slice if needed
            img = img[..., 0]
        return img

def build_roi_masks(ops, stat, iscell):
    """
    Build boolean masks for ROIs listed in `keep_idx`.
    Returns (mask_stack [Ly,Lx,N], labels [N])
    """
    Ly, Lx = ops['Ly'], ops['Lx']
    keep_idx = np.where(iscell == 1)[0]

    mask_stack = np.zeros((Ly, Lx, len(keep_idx)), dtype=bool)
    for i, k in enumerate(keep_idx):
        y = np.asarray(stat[k]['ypix']).astype(int)
        x = np.asarray(stat[k]['xpix']).astype(int)
        m = np.zeros((Ly, Lx), dtype=bool)
        valid = (y >= 0) & (y < Ly) & (x >= 0) & (x < Lx)
        m[y[valid], x[valid]] = True
        mask_stack[..., i] = m
    return mask_stack, np.asarray(keep_idx)



def overlay_masks_on_mean(mean_img,
                          mask_stack,
                          labels=None,
                          *,
                          alpha=0.25,
                          draw_outlines=True,
                          draw_fill=True,
                          label_every=1,
                          one_indexed=False,
                          dot_centroid=False,
                          dot_size=20):
    """
    Overlay ROI masks on a mean image.
    - invert_gray: use 'gray_r' (white background) instead of 'gray'
    - draw_outlines: draw mask boundaries (tries scikit-image; falls back to pixel edges)
    - draw_fill: semi-transparent colored blobs
    - labels: array of ROI indices to display; None = no labels
    - label_every: label every Nth ROI to reduce clutter
    - one_indexed: show 1-based labels
    - dot_centroid: draw a dot at the centroid
    """
    Ly, Lx = mean_img.shape
    cmap = "gray"

    fig, ax = plt.subplots(figsize=(7, 7), dpi=120)
    ax.imshow(mean_img, cmap=cmap, interpolation="nearest")
    ax.set_axis_off()

    n = mask_stack.shape[-1]
    hues = np.linspace(0, 1, max(1, n), endpoint=False)
    colors = plt.cm.hsv(hues)

    from numpy import nonzero

    # optional: outlines via scikit-image
    def draw_outline(mask, color):
        try:
            from skimage import measure
            for contour in measure.find_contours(mask.astype(float), 0.5):
                ax.plot(contour[:, 1], contour[:, 0], '-', linewidth=0.6, color=color)
        except Exception:
            # fallback: crude pixel edge plotting
            yy, xx = nonzero(mask)
            ax.plot(xx, yy, ',', color=color, alpha=0.8)

    for i in range(n):
        m = mask_stack[..., i]
        col = colors[i]

        if draw_fill:
            # show as semi-transparent blob
            # (imshow is faster & nicer than scattering every pixel)
            im = np.zeros((*mean_img.shape, 4), dtype=float)
            im[..., :3] = col[:3]
            im[..., 3] = (m.astype(float) * alpha)
            ax.imshow(im, interpolation="nearest")

        if draw_outlines:
            draw_outline(m, col)

        if labels is not None and (i % max(1, label_every) == 0):
            # centroid for labeling
            ys, xs = nonzero(m)
            if len(xs):
                cx = xs.mean()
                cy = ys.mean()
                text = str(int(labels[i]) + (1 if one_indexed else 0))
                ax.text(cx, cy, text,
                        ha="center", va="center", fontsize=7, color="w",
                        bbox=dict(facecolor="k", alpha=0.55, boxstyle="round,pad=0.15"),
                        clip_on=True)
                if dot_centroid:
                    ax.scatter([cx], [cy], s=dot_size, c=[col], edgecolors="k", linewidths=0.3)

    fig.tight_layout()
    return fig, ax


def overlay_masks_on_mean(mean_img,
                          mask_stack,
                          labels=None,
                          *,
                          alpha=0.25,
                          draw_outlines=False,
                          draw_fill=True,
                          label_every=1,
                          one_indexed=False,
                          dot_centroid=False,
                          dot_size=20):
    """
    Faster overlay: composites all fills/outlines into 1–2 images, then draws once.
    mask_stack: (Ly, Lx, N) boolean
    """
    Ly, Lx = mean_img.shape
    cmap = "gray"

    # colors for N ROIsI
    N = mask_stack.shape[-1]
    hues = np.linspace(0, 1, max(1, N), endpoint=False)
    colors = plt.cm.hsv(hues)[:, :3]  # RGB only

    fig, ax = plt.subplots(figsize=(7, 7), dpi=120)
    ax.imshow(mean_img, cmap=cmap, interpolation="nearest")
    ax.set_axis_off()

    # --------- FILLED OVERLAY (single draw) ---------
    if draw_fill:
        # Label image: per-pixel ROI index, -1 if none
        lab = np.full((Ly, Lx), -1, dtype=int)
        # (Assuming Suite2p masks don't overlap; first-come keeps, else prefer last)
        # Use last-wins (usually fine and avoids branches)
        for i in range(N):
            m = mask_stack[..., i]
            lab[m] = i

        overlay = np.zeros((Ly, Lx, 4), dtype=float)
        covered = lab >= 0
        if covered.any():
            overlay[covered, :3] = colors[lab[covered]]
            overlay[covered, 3] = alpha
            ax.imshow(overlay, interpolation="nearest")

    # --------- OUTLINES (single draw) ---------
    if draw_outlines:
        # Compute edges from the label image quickly (no skimage)
        # Boundary if pixel has a different neighbor label (or neighbor is -1)
        # If we didn’t build lab above (draw_fill False), build a minimal one now
        if not draw_fill:
            lab = np.full((Ly, Lx), -1, dtype=int)
            for i in range(N):
                lab[mask_stack[..., i]] = i

        edge = np.zeros((Ly, Lx), dtype=bool)
        # compare with 4-neighbors using roll (fast)
        for dy, dx in ((-1,0),(1,0),(0,-1),(0,1)):
            nbr = np.roll(lab, shift=(dy, dx), axis=(0,1))
            edge |= (lab != nbr)

        edge &= (lab >= 0)
        if edge.any():
            edge_rgba = np.zeros((Ly, Lx, 4), dtype=float)
            edge_rgba[..., :3] = colors[lab.clip(min=0)]  # color by ROI
            edge_rgba[..., 3] = edge.astype(float)        # alpha 1 on edges, 0 elsewhere
            ax.imshow(edge_rgba, interpolation="nearest")

    # --------- LABELS (vectorized centroids; 1 draw per label) ---------
    if labels is not None and label_every > 0:
        # Compute centroids from masks via pixel means (fast & accurate)
        # Build index arrays once
        yy, xx = np.indices((Ly, Lx))
        shown = 0
        for i in range(0, N, label_every):
            m = mask_stack[..., i]
            if not m.any():
                continue
            y = yy[m].mean() + 8
            x = xx[m].mean() + 8
            text = str(int(labels[i]) + (1 if one_indexed else 0))
            ax.text(x, y, text,
                    ha="center", va="center", fontsize=7, color="w",
                    bbox=dict(facecolor="k", alpha=0.55, boxstyle="round,pad=0.15"),
                    clip_on=True)
            if dot_centroid:
                ax.scatter([x], [y], s=dot_size, c=[colors[i]], edgecolors="k", linewidths=0.3)
            shown += 1

    fig.tight_layout()
    return fig, ax

# for EC_EBc_02 manual registration:
# CELLS 1-10
# # cell 1: 1,3,2,1
# # cell 2

def scree_info(pca, thresh=0.80):
    """Return explained variance ratios, cumulative, and n_keep to reach thresh."""
    explained_var_ratio = pca.explained_variance_ratio_
    cum = np.cumsum(explained_var_ratio)
    n_keep = int(np.searchsorted(cum, thresh) + 1) # first index where cumulative sum of expl. var > threshold.
    return explained_var_ratio, cum, n_keep

def plot_scree(ax, explained_var_ratio, cum, n_keep, title=None, label = None, color = None):
    #fig, ax = plt.subplots()
    xs = np.arange(1, len(explained_var_ratio) + 1)
    #ax.bar(xs, explained_var_ratio, alpha=0.6)
    ax.plot(xs, explained_var_ratio, marker='o', color=color, label = label, alpha = 0.7) # can also plot cum_sum
    # ax.axhline(0.80, ls='--', lw=1, label='80%')
    # ax.axvline(n_keep, ls='--', lw=1, color='k', label=f'n_keep={n_keep}')
    ax.set_xlabel('PC')
    ax.set_ylabel('Explained Variance')
    ax.set_title(title)
    ax.set_ylim(0, 0.3)
    ax.set_xlim(0, 20)
    ax.legend(frameon=False)

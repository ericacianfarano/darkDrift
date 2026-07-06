from track2p.t2p import run_t2p
from track2p.ops.default import DefaultTrackOps
import os

# C:\Users\erica\.conda\envs\track2p\Lib\site-packages\track2p

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


if __name__ == '__main__':

    across_days = False

    if across_days:

        animals = ['EC_EB_40']
        drive = 'I'

        for animal in animals:
            for stim in ['chunk']:
                for roi_detection in ['functional']:# 'anatomical']:

                    print(f'{animal}: {stim} stimulus, {roi_detection} suite2p detection')

                    root_dir = f'{drive}:\sensorium\data\{animal}'
                    files = []

                    for dirpath, dirnames, filenames in os.walk(root_dir):
                        if os.path.basename(dirpath) == 'experiments':
                            if stim in dirpath.split('\\')[5]:
                                files.append(dirpath)

                    # sort folder based on when it was created > want chronological ordering
                    files.sort(key=lambda f: get_earliest_file_time(f, use_mtime=True))

                    # print(files)
                    # Get default parameters
                    track_ops = DefaultTrackOps()
                    track_ops.s2p_folder_name = f"suite2p {roi_detection}"
                    track_ops.t2p_save_folder_name = f"track2p-{stim}-{roi_detection}"

                    track_ops.all_ds_path = files

                    #track_ops.save_path = r'E:\DriftScape\Data\EC_GECO_09\test' # path where to save the outputs of algorithm (a 'track2p' folder will be created where figures for visualisation and matrices of matches would be saved)
                    track_ops.save_path = fr'{drive}:\dark_drift\data\{animal}'#\track2p-{stim}-{roi_detection}'  # path where to save the outputs of algorithm (a 'track2p' folder will be created where figures for visualisation and matrices of matches would be saved)

                    # if the registration is poor, try cropping out black edges
                    if (animal == 'EC_EBc_02') or (animal == 'EC_EBc_03') or (animal == 'EC_EB_01') or (animal == 'EC_EBc_09') or (animal == 'EC_LBc_12')\
                            or (animal == 'EC_EB_24'):
                        track_ops.border_crop = 30
                    elif animal == 'EC_EB_21':
                        track_ops.border_crop = 30
                        track_ops.reg_chan= 1
                        # track_ops.transform_type = 'rigid'
                        # track_ops.iou_dist_thr = 2
                    elif animal == 'EC_EB_22':
                        track_ops.border_crop = 10
                        track_ops.reg_chan= 1
                        # track_ops.transform_type = 'rigid'
                        track_ops.iou_dist_thr = 10
                    # elif animal == 'EC_EB_28':
                        # track_ops.border_crop = 30
                        # track_ops.reg_chan= 1
                        # track_ops.transform_type = 'rigid'
                        # track_ops.iou_dist_thr = 10
                        # running rigid helps align last w5-w6 weeks. but affine helps align w3-w5 (all else unchanged)
                    elif ((animal == 'EC_EBc_04') or (animal == 'EC_EBc_05') or (animal == 'EC_EBc_06') or (animal == 'EC_EBc_07')
                          or (animal == 'EC_LBc_02') or (animal == 'EC_LBc_03') or (animal == 'EC_LBc_06')
                          or (animal == 'EC_LB_01') or (animal == 'EC_LB_02') or (animal == 'EC_LB_03') or (animal == 'EC_LB_05') or (animal == 'EC_LB_06') or (animal == 'EC_LB_09') or (animal == 'EC_LB_10') or (animal == 'EC_LB_11')
                    ):
                        track_ops.border_crop = 0
                        track_ops.iou_dist_thr = 16
                    elif (animal == 'EC_LBc_13'):
                        track_ops.border_crop = 20 # need to crop spon 30 cm
                    else: # registration is good without cropping too much
                        track_ops.border_crop = 0

                    # track_ops.transform_type = 'rigid'
                    track_ops.iou_dist_thr = 16  # 4   # decreasing this makes matches more conservative (less matches but they need to overlap more strongly). default is 16

                    # something to implement is despite the crop, align rois on the whole image.

                    # if the registration is poor, use a different regisntration method
                    if (animal == 'EC_EB_01'):
                        track_ops.transform_type = 'rigid'  ############# used for EB_01
                        track_ops.iou_dist_thr = 4   # decreasing this makes matches more conservative (less matches but they need to overlap more strongly). default is 16
                    if (animal == 'EC_EB_13'): # tried normal transform: 0/10/20/30/50/70/100 crop, did not work
                        track_ops.border_crop = 30
                        # track_ops.transform_type = 'rigid'  ############# used for EB_01
                        # track_ops.iou_dist_thr = 4   # decreasing this makes matches more conservative (less matches but they need to overlap more strongly). default is 16
                    track_ops.thr_remove_zeros = True ####################################
                    # track_ops.sat_perc = 99.9  ##################
                    # track_ops.thr_method = 'min' ###################
                    # track_ops.iou_dist_thr = 5
                    # if tracking is poor
                    print(track_ops.iou_dist_thr, track_ops.matching_method, track_ops.win_size)
                    #track_ops.matching_method = 'dist'  # instead of 'iou'

                    #track_ops.win_size = 32 # default is 48

                    #track_ops.reg_chan = 0 # channel to use for registration (0=functional, 1=anatomical) (use 0 if only recording gcamp!)
                    if roi_detection == 'anatomical':
                        track_ops.iscell_thr = 0 # set this to 0 (basically take all cells)
                        track_ops.reg_chan = 0
                    elif roi_detection == 'functional':
                        track_ops.iscell_thr = 0.02 # small positive to get rid of garbage
                        track_ops.reg_chan = 0

                        if animal == 'EC_EB_10' or animal == 'EC_EB_31': # only want right half of FOV > only want cells we manually classified as cells
                            track_ops.iscell_thr = None
                        #track_ops.matching_method = 'cent'
                        #track_ops.transform_type = 'nonrigid' #nonrigid or affine
                        #track_ops.iou_dist_thr = 30

                    #print(track_ops)
                    #print(track_ops.save_path)
                    # print all the settings / parameters used for running the algorithm
                    for attr, value in track_ops.__dict__.items():
                        print(attr, '=', value)

                    # Run the algorithm
                    run_t2p(track_ops)

    ############################################################## track2p on all sessions within a day
    else: #not across days, aligning within a day
        animals = ['EC_LB_12','EC_LB_13', 'EC_LB_16'] #'EC_LBc_09', 'EC_LBc_11','EC_LBc_13','EC_LBc_14','EC_LBc_16',
        # animals = ['EC_EB_01', 'EC_EB_03', 'EC_EB_05', 'EC_EB_08', 'EC_EB_10', 'EC_EB_13', 'EC_EB_14',
        #            'EC_EBc_01', 'EC_EBc_02', 'EC_EBc_03', 'EC_EBc_04', 'EC_EBc_05', 'EC_EBc_06', 'EC_EBc_07',
        #            'EC_EBc_09', 'EC_EBc_11',
        #            'EC_LB_01', 'EC_LB_02', 'EC_LB_03', 'EC_LB_05', 'EC_LB_06', 'EC_LB_07', 'EC_LB_09', 'EC_LB_10',
        #            'EC_LBc_02', 'EC_LBc_03', 'EC_LBc_06']
        drive = 'I'
        for animal in animals:
            # root_dir = f'{drive}:\dark_drift\data\{animal}'
            root_dir = f'{drive}:\sensorium\data\{animal}'
            for day in [x for x in os.listdir(root_dir) if x.isdigit()]:
                files = []

                for stim in ['chunk']: #['grat', 'spon']:
                    for roi_detection in ['functional']:# 'anatomical']:

                        for dirpath, dirnames, filenames in os.walk(root_dir):
                            if os.path.basename(dirpath) == 'experiments' and (day in dirpath) and (stim in dirpath):
                                if stim in dirpath.split('\\')[5]:
                                    files.append(dirpath)


                    # sort folder based on when it was created > want chronological ordering
                    files.sort(key=lambda f: get_earliest_file_time(f, use_mtime=True))

                if len (files) > 1: # if we have +1 recording

                    # Get default parameters
                    track_ops = DefaultTrackOps()
                    track_ops.s2p_folder_name = f"suite2p {roi_detection}"
                    track_ops.t2p_save_folder_name = f"track2p-{roi_detection}"

                    track_ops.all_ds_path = files

                    #track_ops.save_path = r'E:\DriftScape\Data\EC_GECO_09\test' # path where to save the outputs of algorithm (a 'track2p' folder will be created where figures for visualisation and matrices of matches would be saved)
                    # track_ops.save_path = fr'E:\dark_drift\data\{animal}\{day}'#\track2p-{stim}-{roi_detection}'  # path where to save the outputs of algorithm (a 'track2p' folder will be created where figures for visualisation and matrices of matches would be saved)
                    track_ops.save_path = fr'{drive}:\sensorium\data\{animal}\{day}'  # \track2p-{stim}-{roi_detection}'  # path where to save the outputs of algorithm (a 'track2p' folder will be created where figures for visualisation and matrices of matches would be saved)


                    # if the registration is poor, try cropping out black edges
                    if (animal == 'EC_EBc_02') or (animal == 'EC_EBc_03') or (animal == 'EC_EB_01') or (animal == 'EC_EBc_09'):
                        track_ops.border_crop = 30
                    elif ((animal == 'EC_EBc_04') or (animal == 'EC_EBc_05') or (animal == 'EC_EBc_06') or (animal == 'EC_EBc_07')
                          or (animal == 'EC_LBc_02') or (animal == 'EC_LBc_03') or (animal == 'EC_LBc_06')
                          or (animal == 'EC_LB_01') or (animal == 'EC_LB_02') or (animal == 'EC_LB_03') or (animal == 'EC_LB_05') or (animal == 'EC_LB_06') or (animal == 'EC_LB_09') or (animal == 'EC_LB_10') or (animal == 'EC_LB_11')
                    ):
                        track_ops.border_crop = 0
                        track_ops.iou_dist_thr = 16
                    else: # registration is good without cropping too much
                        track_ops.border_crop = 0

                    # track_ops.transform_type = 'rigid'
                    track_ops.iou_dist_thr = 16  # 4   # decreasing this makes matches more conservative (less matches but they need to overlap more strongly). default is 16

                    # something to implement is despite the crop, align rois on the whole image.

                    # if the registration is poor, use a different regisntration method
                    if (animal == 'EC_EB_01'):
                        track_ops.transform_type = 'rigid'  ############# used for EB_01
                        track_ops.iou_dist_thr = 4   # decreasing this makes matches more conservative (less matches but they need to overlap more strongly). default is 16
                    if (animal == 'EC_EB_13'): # tried normal transform: 0/10/20/30/50/70/100 crop, did not work
                        track_ops.border_crop = 30
                        # track_ops.transform_type = 'rigid'  ############# used for EB_01
                        # track_ops.iou_dist_thr = 4   # decreasing this makes matches more conservative (less matches but they need to overlap more strongly). default is 16
                    track_ops.thr_remove_zeros = True ####################################
                    # track_ops.sat_perc = 99.9  ##################
                    # track_ops.thr_method = 'min' ###################
                    # track_ops.iou_dist_thr = 5
                    # if tracking is poor
                    print(track_ops.iou_dist_thr, track_ops.matching_method, track_ops.win_size)
                    #track_ops.matching_method = 'dist'  # instead of 'iou'

                    #track_ops.win_size = 32 # default is 48

                    #track_ops.reg_chan = 0 # channel to use for registration (0=functional, 1=anatomical) (use 0 if only recording gcamp!)
                    if roi_detection == 'anatomical':
                        track_ops.iscell_thr = 0 # set this to 0 (basically take all cells)
                        track_ops.reg_chan = 0
                    elif roi_detection == 'functional':
                        track_ops.iscell_thr = 0.05 # small positive to get rid of garbage
                        track_ops.reg_chan = 0

                        if animal == 'EC_EB_10': # only want right half of FOV > only want cells we classified as cells
                            track_ops.iscell_thr = None
                        #track_ops.matching_method = 'cent'
                        #track_ops.transform_type = 'nonrigid' #nonrigid or affine
                        #track_ops.iou_dist_thr = 30

                    #print(track_ops)
                    #print(track_ops.save_path)
                    # print all the settings / parameters used for running the algorithm
                    for attr, value in track_ops.__dict__.items():
                        print(attr, '=', value)

                    # Run the algorithm
                    run_t2p(track_ops)
        #############################################################################################################
'''    # FOR TESTING T2P
    animal = 'EC_ctrl_11'
    stim = 'grat'
    roi_detection = 'functional'

    print(f'{animal}: {stim} stimulus, {roi_detection} suite2p detection')

    # stim = 'grat' # 'spon'
    # roi_detection = 'functional'  # or functional
    root_dir = f'E:\dark_drift\data\{animal}'
    files = []

    for dirpath, dirnames, filenames in os.walk(root_dir):
        if os.path.basename(dirpath) == 'experiments':
            if stim in dirpath:
                files.append(dirpath)

    # folder creation time
    #files.sort(key=os.path.getctime)
    files.sort(key=lambda f: get_earliest_file_time(f, use_mtime=True))

    #files = files[-2:]

    # Get default parameters
    track_ops = DefaultTrackOps()
    #print(track_ops)
    track_ops.s2p_folder_name = f"suite2p {roi_detection}"
    track_ops.t2p_save_folder_name = f"track2p-{stim}-{roi_detection}"

    track_ops.all_ds_path = files

    #track_ops.save_path = r'E:\DriftScape\Data\EC_GECO_09\test' # path where to save the outputs of algorithm (a 'track2p' folder will be created where figures for visualisation and matrices of matches would be saved)
    track_ops.save_path = fr'E:\dark_drift\data\{animal}'#\track2p-{stim}-{roi_detection}'  # path where to save the outputs of algorithm (a 'track2p' folder will be created where figures for visualisation and matrices of matches would be saved)

    #track_ops.reg_chan = 0 # channel to use for registration (0=functional, 1=anatomical) (use 0 if only recording gcamp!)
    if roi_detection == 'anatomical':
        track_ops.iscell_thr = 0 # set this to 0 (basically take all cells)
        track_ops.reg_chan = 0
    elif roi_detection == 'functional':
        track_ops.iscell_thr = 0.1 # small positive to get rid of garbage
        track_ops.reg_chan = 0.
        #track_ops.matching_method = 'cent'
        track_ops.transform_type = 'nonrigid' #nonrigid or affine
        # track_ops.iou_dist_thr = 1
    #
    # #print(track_ops)
    # #print(track_ops.save_path)
    # print all the settings / parameters used for running the algorithm
    for attr, value in track_ops.__dict__.items():
        print(attr, '=', value)

    # Run the algorithm
    run_t2p(track_ops)'''

# to run from terminal: python -m run_track2p

# changes to files
    # 1) in 'loaders.py' file (in track2p > io), in line 22:
            # print(f'Loading ROIs for plane{plane_idx} in dataset {track_ops_path.split("/")[-2]}')
        # changed to
            # sess_id = track_ops_path.split("/")[-2]
            # print(f'Loading ROIs for plane{plane_idx} in dataset {sess_id}')
        # the issue was trying to get the "rx_xxx_000' string from the path,
        # but we couldn't separate the path according to '/', we can only separate
        # it with '\\' (but can't use this expresion with f-statement so we store in variable)

    # 2) in 'elastix.py' file, on line 40-42
        # original code:
            # all_roi_array_reg[:,:,i] = roi_array_reg
        # changed to
            # dim0, dim1 = roi_array_reg.shape
            # all_roi_array_reg[:dim0,:dim1,i] = roi_array_reg

        # the issue is we were trying to store roi_array_reg (shape 800 x 796)
        # in all_roi_array_reg array (shape 800 x 800), so we were getting a
        # valueError / broadcasting mismatch. so by storing roi_array_reg shape
        # and indexing to those dimensions, we fix the error

    # 3) in 'output.py' file, on line 184
        # original code:
            # this_ax = axs[i] if track_ops.nplanes==1 else axs[j][i]
        # changed to
            # added if/else where [ (track_ops.nplanes == 1) and (len(all_ds_thr_met) == 1)]
            # axs[i] changed to ax

    #4 ) in 'loaders.py' file> line 32 in get_all_roi_array_from_stat()
        # original code:
            # n_ypix = track_ops.all_ds_avg_ch1[0][0].shape[0]
        # changed to;
            # n_ypix = track_ops.all_ds_avg_ch1[0][0].shape[1]

    #5) Also changed all instances of 'suite2p' to a new variable stored in the ops file called  track_ops.s2p_folder_name
    # where we can specify whether it should run on functional or anatomical suite2p files
    # changed in loaders, in t2p, suite2p loaders, etc


# note that itk wont work with python 3.9, only 3.8


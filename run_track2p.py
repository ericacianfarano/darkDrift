from track2p.t2p import run_t2p
from track2p.ops.default import DefaultTrackOps
import os

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

    animal = 'EC_dark_08'
    stim = 'spon' # 'spon'
    root_dir = f'I:\dark_drift\data\{animal}'
    files = []

    for dirpath, dirnames, filenames in os.walk(root_dir):
        if os.path.basename(dirpath) == 'experiments':
            if stim in dirpath:
                files.append(dirpath)

    # folder creation time
    #files.sort(key=os.path.getctime)
    files.sort(key=lambda f: get_earliest_file_time(f, use_mtime=True))

    # animal = 'EC_GCaMP6s_09'
    #
    # # subfiles
    # if animal == 'EC_GCaMP6s_06':
    #     last_digits = [0,0,0,0,1,0]
    # elif animal == 'EC_GCaMP6s_09':
    #     last_digits = [1,1,0,0,2,0,0]

    # Get default parameters
    track_ops = DefaultTrackOps()

    # for now parameters are defined manually:
    # track_ops.all_ds_path = [                           # list of paths to datasets containing a `suite2p` folder
    #             fr'I:\dark_drift_trial\data\{animal}\20250213\big100_ori\big100_ori_000_00{last_digits[0]}\experiments',
    #             #fr'I:\dark_drift_trial\data\{animal}\20250221\big100_ori\big100_ori_000_00{last_digits[1]}\experiments',
    #             #fr'I:\dark_drift_trial\data\{animal}\20250303\big100_ori\big100_ori_000_00{last_digits[2]}\experiments',
    #             #fr'I:\dark_drift_trial\data\{animal}\20250311\big100_ori\big100_ori_000_00{last_digits[3]}\experiments',
    #             fr'I:\dark_drift_trial\data\{animal}\20250319\big100_ori\big100_ori_000_00{last_digits[4]}\experiments',
    #             fr'I:\dark_drift_trial\data\{animal}\20250425\big100_ori\big100_ori_000_00{last_digits[5]}\experiments',
    #             fr'I:\dark_drift_trial\data\{animal}\20250523\big100_ori\big100_ori_000_00{last_digits[6]}\experiments'
    #         ]

    track_ops.all_ds_path = files
    # # # for now parameters are defined manually:
    # track_ops.all_ds_path = [                           # list of paths to datasets containing a `suite2p` folder
    #             r'E:\DriftScape\Data\EC_GECO_09\20231106\r2\r2_205_000\experiments',
    #             r'E:\DriftScape\Data\EC_GECO_09\20231107\r1\r1_228_000\experiments',
    #             r'E:\DriftScape\Data\EC_GECO_09\20231108\r1\r1_192_000\experiments',
    #             r'E:\DriftScape\Data\EC_GECO_09\20231109\r1\r1_210_000\experiments'
    #         ]

    #track_ops.save_path = r'E:\DriftScape\Data\EC_GECO_09\test' # path where to save the outputs of algorithm (a 'track2p' folder will be created where figures for visualisation and matrices of matches would be saved)
    track_ops.save_path = fr'I:\dark_drift\data\{animal}\track2p-{stim}'  # path where to save the outputs of algorithm (a 'track2p' folder will be created where figures for visualisation and matrices of matches would be saved)

    track_ops.reg_chan = 0 # channel to use for registration (0=functional, 1=anatomical) (use 0 if only recording gcamp!)
    track_ops.iscell_thr = 0.1 # set this to 0 (basically take all cells)

    # #print(track_ops)
    # #print(track_ops.save_path)
    # # print all the settings / parameters used for running the algorithm
    # for attr, value in track_ops.__dict__.items():
    #     print(attr, '=', value)

    # Run the algorithm
    run_t2p(track_ops)

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



# note that itk wont work with python 3.9, only 3.8


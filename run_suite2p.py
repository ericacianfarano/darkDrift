import os
import numpy as np
import glob
from suite2p import run_s2p
import imageio
from sbxreader import sbx_memmap
import tifffile as tiff, os

def convert_sbx_tiff (path_to_exp):
    '''
    function created because i accidentally saved two channels. sujite2p shits its pants when you make it process both

    so convert sbx file to tiff (only 1 channel), so suite2p only needs to grab one of the channels (cant do this with sbx files)

    if this occurs, you need to uncomment the specific OPS settings in the suite 2p, but also switch files.append(dirpath) with files.append(os.path.join(dirpath, 'tiffs')
    also, you need to change F = dat[i:min(i+chunk, N), 0, 0] to F = dat[i:min(i+chunk, N), 0, 1]

    :param sbx_path: path to scanbox file
    :param out_dir: where to save tiffs

    :return:
    '''
    sbx_path = os.path.join (path_to_exp, path_to_exp.split('\\')[6]+'.sbx')
    out_dir = os.path.join(path_to_exp, 'tiffs')

    print(sbx_path, out_dir)

    os.makedirs(out_dir, exist_ok=True)

    dat = sbx_memmap(sbx_path)          # (frames, planes, channels, Ly, Lx)
    N = dat.shape[0]
    chunk = 2000
    i = k = 0
    while i < N:
        F = dat[i:min(i+chunk, N), 0, 1]  # plane 0, channel index 1 = PMT2 -> use dat[i:min(i+chunk, N), 0, 1] if you accidentally record extra PMT
        tiff.imwrite(os.path.join(out_dir, f"chunk_{k:05d}.tif"), F, photometric='minisblack')
        i += chunk; k += 1

# convert_sbx_tiff ( fr"E:\dark_drift\data\EC_LBc_09\20260225\grat\grat_000_003\experiments")
# convert_sbx_tiff ( fr"E:\dark_drift\data\EC_LBc_09\20260225\spon0\spon0_000_000\experiments")
# convert_sbx_tiff ( fr"E:\dark_drift\data\EC_EB_08\20250610\grat\grat_000_000\experiments")
# convert_sbx_tiff ( fr"E:\dark_drift\data\EC_LB_10\20251016\spon0\spon0_000_003\experiments")

def chunks_to_mp4(tiff_dir, fps=20, vmin=None, vmax=None):
    """
    Take chunked TIFF stacks (chunk_00000.tif, chunk_00001.tif, ...)
    and write them as an MP4 movie.

    Parameters
    ----------
    tiff_dir : str > path to experiments
        Directory containing chunk_*.tif files.
    out_mp4 : str
        Output MP4 file path.
    fps : int or float
        Frames per second for the video.
    vmin, vmax : int or None
        Intensity range to map to 0–255. If None, they will be estimated
        from the data of the first chunk.
    """

    out_mp4 = os.path.join(tiff_dir, 'vid.mp4')

    # find chunks in order
    tiff_paths = sorted(glob.glob(os.path.join(tiff_dir,'tiffs', "chunk_*.tif")))
    if not tiff_paths:
        raise FileNotFoundError(f"No chunk_*.tif found in {tiff_dir}")

    print(f"Found {len(tiff_paths)} chunk files")

    # if needed, estimate vmin/vmax from first chunk
    if vmin is None or vmax is None:
        first_stack = tiff.imread(tiff_paths[0])  # (n_frames, y, x)
        vmin = np.percentile(first_stack, 1) if vmin is None else vmin
        vmax = np.percentile(first_stack, 99) if vmax is None else vmax
        print(f"Auto intensity range: vmin={vmin:.1f}, vmax={vmax:.1f}")

    def to_uint8(frame):
        frame = np.clip((frame - vmin) / (vmax - vmin) * 255, 0, 255)
        return frame.astype(np.uint8)

    # open video writer (let imageio pick the right plugin, usually ffmpeg)
    with imageio.get_writer(out_mp4, fps=fps, codec="libx264") as writer:
        for p in tiff_paths:
            print("Adding frames from:", os.path.basename(p))
            stack = tiff.imread(p)  # (n_frames, y, x)
            for frame in stack:
                writer.append_data(to_uint8(frame))

    print("Done, wrote:", out_mp4)

# chunks_to_mp4(fr"E:\dark_drift\data\EC_EB_14\20251124\grat\grat_000_000\experiments", fps=30, vmin=None, vmax=None, smooth_window = 3)
# chunks_to_mp4 ( fr"E:\dark_drift\data\EC_LBc_09\20260225\grat\grat_000_003\experiments")

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

    animals = ['EC_LB_12','EC_LB_13','EC_LB_16'] #'EC_LBc_03' #need to run ctrl 12
    drive = 'K'
    day =None# '20251107' #None
    stim = None # 'spon' or 'grat', if None, will run both

    for animal in animals:
        for roi_detection in ['functional']:

            print(f'Running suite2p for {animal} {day} {stim} {roi_detection}')

            # load ops file
            if roi_detection == 'functional':
                ops = np.load(fr'C:\Users\erica\suite2p_thy1gcamp6s.npy', allow_pickle=True).item()
            elif roi_detection == 'anatomical':
                ops = np.load(fr'C:\Users\erica\suite2p_anatomical_cellpose.npy', allow_pickle=True).item()
                # tried cellprob = 0.5 and flow = 1.0, but got few rois with clean edge
                # tried cell prob = -1 and flow = 2. good boundaries. not that many more rois
                # flow threshold doesnt seem tobe doing much
                ops['cellprob_threshold'] = 2 #tried 0.5
                # ops['flow_threshold'] = 1.5  # tighten boundaries (default 1.5 is looser) # tried 1.0, cell look good but didnt get many ROIS
                # ops['preclassify'] = 0.0
                # ops['use_builtin_classifier'] = False
                ops['pretrained_model'] = r"C:\Users\erica\.cellpose\models\cyto2torch_0"
                print(ops)

                #ops['cellprob_threshold'] = 0.5  # try 0.5; if still big, 0.6–0.7
                #ops['flow_threshold'] = 1.0  # tighten boundaries (default 1.5 is looser)
                # ops['spatial_hp_cp'] = 50  # suppress broad halos (try 30–80 range)
                # ops['pretrained_model'] = 'cyto2'  # often hugs soma edges better on 2P
                #ops['diameter'] = 5

                #ops['flow_threshold'] = 0.5

            # modify ops file

            ops['save_mat'] = 1
            ops['fs'] = 15
            ops['fast_disk'] = 'C:\\'
            ops['input_format'] = 'sbx'
            ops['save_folder'] = f'suite2p {roi_detection}'
            ops['delete_bin'] = True
            ops['keep_movie_raw'] = False
            ops['reg_tif'] = False
            ops['batch_size'] = 150

            # when running on TIFF files. if running on sbx, comment this out
            # ops ['input_format'] = 'tiff'
            # ops ['nchannels'] = 1
            # ops ['functional_chan'] = 1
            # ops ['align_by_chan'] = 1
            # ops ['save_chan2'] = False
            # ops ['nplanes'] = 1
            # ops ['delete_bin'] = True
            # ops ['keep_movie_raw'] = False
            # ops ['reg_tif'] = False
            # ops ['batch_size'] = 250
            #########################################

            # root_dir = f'{drive}:\dark_drift\data\{animal}'
            root_dir = f'{drive}:\sensorium\data\{animal}'
            # root_dir = fr'{drive}:\vision_restored\{animal}'

            print(root_dir)

            #print(ops)

            files = []

            if day: # if a day is specified, only process files in that day
                for dirpath, dirnames, filenames in os.walk(os.path.join(root_dir, day)):
                    if os.path.basename(dirpath) == 'experiments':
                        if stim:
                            if stim in dirpath and ('extra_rec' not in dirpath):
                                # files.append(os.path.join(dirpath, 'tiffs'))
                                files.append(dirpath)
                        else:
                            if ('extra_rec' not in dirpath):
                                # files.append(os.path.join(dirpath, 'tiffs'))
                                files.append(dirpath)

            else: # if a day isnt specified, process all files
                for dirpath, dirnames, filenames in os.walk(root_dir):
                    if os.path.basename(dirpath) == 'experiments':
                        if stim:
                            if stim in dirpath:
                                files.append(dirpath)
                        else:
                            files.append(dirpath)


            # sort folders according to creation time
            files.sort(key=lambda f: get_earliest_file_time(f, use_mtime=True))
            print('files', files)
            for file in files:
                db = {'data_path':[file]}
                opsEND = run_s2p(ops=ops, db=db)

                print("Outputs saved to:", opsEND['save_path'])  # e.g., E:\s2p_out\suite2p
                print("Ops file:", opsEND['ops_path'])


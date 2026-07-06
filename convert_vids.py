import tifffile as tiff, os
from imports import *

def convert_sbx_tiff (dir_path):
    '''
    function created because i accidentally saved two channels. sujite2p shits its pants when you make it process both

    so convert sbx file to tiff (only 1 channel), so suite2p only needs to grab one of the channels (cant do this with sbx files)

    :param sbx_path: path to scanbox file
    :param out_dir: where to save tiffs

    :return:
    '''

    path_to_exp = os.path.join(dir_path, 'experiments')

    sbx_path = os.path.join(path_to_exp, path_to_exp.split('\\')[-2] + '.sbx')
    # sbx_path = os.path.join (path_to_exp, path_to_exp.split('\\')[6]+'.sbx')
    out_dir = os.path.join(path_to_exp, 'tiffs')

    print(sbx_path, out_dir)

    os.makedirs(out_dir, exist_ok=True)

    dat = sbx_memmap(sbx_path)          # (frames, planes, channels, Ly, Lx)
    N = dat.shape[0]
    chunk = 2000
    i = k = 0
    while i < N:
        F = dat[i:min(i+chunk, N), 0, 0]  # plane 0, channel index 1 = PMT2 -> use dat[i:min(i+chunk, N), 0, 1] if you accidentally record extra PMT
        tiff.imwrite(os.path.join(out_dir, f"chunk_{k:05d}.tif"), F, photometric='minisblack')
        i += chunk; k += 1
#
# convert_sbx_tiff ( fr"E:\dark_drift\data\EC_EB_08\20250610\grat\grat_000_000\experiments")
# convert_sbx_tiff ( fr"E:\dark_drift\data\EC_LB_10\20251016\spon0\spon0_000_003\experiments")
def chunks_to_mp4(
    dir,
    fps=20,
    vmin=None,
    vmax=None,
    smooth_window=1,
    add_stim_bar=True,    # we will rename later if you want
    stim_duration_s=4.0,
    border_px=6,          # thickness of green border
):

    tiff_dir = os.path.join(dir, 'experiments')

    out_mp4 = os.path.join(tiff_dir, 'vid.mp4')

    # ---- build stim frame set from TTLs ----
    stim_frame_set = None
    if add_stim_bar:
        stim_frame_set = make_stim_window_frames(
            r_path=dir,
            fps=fps,
            duration_s=stim_duration_s,
        )
        print(f"Stim overlay: {len(stim_frame_set)} frames in stim windows")

    crop_start = max(list(stim_frame_set)[0] - 80, 0)

    # re-index ttls frames
    stim_frame_set = {f - crop_start for f in stim_frame_set if f >= crop_start}

    # ---- find TIFF chunks ----
    tiff_paths = sorted(glob.glob(os.path.join(tiff_dir, 'tiffs', "chunk_*.tif")))
    if not tiff_paths:
        raise FileNotFoundError(f"No chunk_*.tif found in {os.path.join(tiff_dir, 'tiffs')}")

    print(f"Found {len(tiff_paths)} chunk files")

    # ---- estimate vmin/vmax ----
    if vmin is None or vmax is None:
        first_stack = tiff.imread(tiff_paths[0])
        if vmin is None:
            vmin = np.percentile(first_stack, 1)
        if vmax is None:
            vmax = np.percentile(first_stack, 99)
        print(f"Auto intensity range: vmin={vmin:.1f}, vmax={vmax:.1f}")

    def to_uint8(frame):
        frame = np.clip((frame - vmin) / (vmax - vmin) * 255, 0, 255)
        return frame.astype(np.uint8)

    global_frame_idx = 0  # original frame index
    out_idx = 0           # index in cropped movie

    with imageio.get_writer(out_mp4, fps=fps, codec="libx264") as writer:
        for p in tiff_paths:
            print("Adding frames from:", os.path.basename(p))
            stack = tiff.imread(p).astype(np.float32)

            if smooth_window > 1:
                stack = uniform_filter1d(stack, size=smooth_window, axis=0, mode='nearest')

            for frame in stack:
                # only start writing once we hit crop_start
                if global_frame_idx >= crop_start:
                    g = to_uint8(frame)
                    frame_rgb = np.stack([g, g, g], axis=-1)

                    if stim_frame_set is not None and out_idx in stim_frame_set:
                        frame_rgb[:border_px, :, :]  = [0, 255, 0]
                        frame_rgb[-border_px:, :, :] = [0, 255, 0]
                        frame_rgb[:, :border_px, :]  = [0, 255, 0]
                        frame_rgb[:, -border_px:, :] = [0, 255, 0]

                    writer.append_data(frame_rgb)
                    out_idx += 1

                global_frame_idx += 1

    print("Total frames IN 2p movie:", global_frame_idx)
    print("Total frames written to cropped 2p movie:", out_idx)
    print("Done, wrote:", out_mp4)


def make_stim_window_frames(r_path, fps = 20, duration_s=4.0):
    """
    stim_onsets : array-like of global frame indices where TTL starts
    fps        : frames per second of imaging
    duration_s : how long after onset to mark (seconds)
    """

    recording_path = os.path.join(r_path, 'experiments')
    ttl_path = os.path.join(recording_path, f'{recording_path.split(os.sep)[-2]}.mat')
    ttl_info = loadmat(ttl_path)['info'][0][0]

    event_id = np.squeeze(ttl_info['event_id'])
    frame_idx = np.squeeze(ttl_info['frame'])

    # e.g. event_id == 1 means stimulus is present
    stim_frames = frame_idx[event_id == 1]

    stim_onsets = np.asarray(stim_frames, dtype=int)[1:] #exclude first ttl (which is just start of rec

    n_after = int(round(duration_s * fps))

    stim_frame_set = set()
    for onset in stim_onsets:
        # include onset itself and the next n_after frames
        for f in range(onset, onset + n_after):
            stim_frame_set.add(f)
    return stim_frame_set

import imageio.v2 as imageio
import cv2
import imageio.v2 as imageio
import numpy as np
import os
from scipy.io import loadmat  # already used elsewhere

def behav_with_ttl(
    dir,
    stim_duration_s=4.0,
    border_px=12,
    out_name=None,
):
    """
    Open the behaviour video in `dir/behav` and save a copy with
    a green border whenever TTL is ON (for stim_duration_s seconds
    after each onset).

    Parameters
    ----------
    dir : str
        Root recording directory, e.g. ...\\EC_EB_01\\20250609\\grat\\grat_000_000
        Must contain:
            - 'experiments'  (with TTL .mat)
            - 'behav'        (with <recording_root>_eye.MJ2)
    stim_duration_s : float
        Seconds after each TTL onset to display the border.
    border_px : int
        Border thickness in pixels.
    out_name : str or None
        Name of output MP4 inside behav folder. If None, will use
        '<recording_root>_eye_ttl.mp4'
    """

    behav_dir = os.path.join(dir, 'behav')

    # recording_root: e.g. 'grat_000_000'
    recording_root = os.path.basename(dir)

    # input behaviour video: <recording_root>_eye.MJ2
    in_video = os.path.join(behav_dir, f"{recording_root}_eye.MJ2")
    if not os.path.exists(in_video):
        raise FileNotFoundError(f"Behaviour video not found: {in_video}")

    if out_name is None:
        out_name = f"{recording_root}_eye_ttl.mp4"
    out_video = os.path.join(behav_dir, out_name)

    print(f"Behaviour video in:  {in_video}")
    print(f"TTL-annotated video: {out_video}")

    # ---- open with OpenCV ----
    cap = cv2.VideoCapture(in_video)
    if not cap.isOpened():
        raise IOError(f"Could not open behaviour video: {in_video}")

    fps = 20

    # ---- build stim frame set using same TTL logic ----
    stim_frame_set = make_stim_window_frames(
        r_path=dir,
        fps=fps,
        duration_s=stim_duration_s,
    )

    crop_start = max(list(stim_frame_set)[0] - 80, 0)

    # re-index ttls frames
    stim_frame_set = {f - crop_start for f in stim_frame_set if f >= crop_start}

    print(f"Stim overlay: {len(stim_frame_set)} frames in stim windows")

    writer = imageio.get_writer(out_video, fps=fps, codec="libx264")

    frame_idx = 0  # original behaviour frame index
    out_idx = 0    # cropped index

    try:
        while True:
            ret, frame_bgr = cap.read()
            if not ret:
                break

            if frame_idx >= crop_start:
                frame_rgb = cv2.cvtColor(frame_bgr, cv2.COLOR_BGR2RGB)

                frame_yuv = cv2.cvtColor(frame_rgb, cv2.COLOR_RGB2YUV)
                frame_yuv[:, :, 0] = clahe.apply(frame_yuv[:, :, 0])
                frame_rgb = cv2.cvtColor(frame_yuv, cv2.COLOR_YUV2RGB)

                if out_idx in stim_frame_set:
                    h, w, _ = frame_rgb.shape
                    b = min(border_px, h // 2, w // 2)
                    frame_rgb[:b, :, :]  = [0, 255, 0]
                    frame_rgb[-b:, :, :] = [0, 255, 0]
                    frame_rgb[:, :b, :]  = [0, 255, 0]
                    frame_rgb[:, -b:, :] = [0, 255, 0]

                writer.append_data(frame_rgb)
                out_idx += 1

            frame_idx += 1
    finally:
        cap.release()
        writer.close()

    print("Total frames IN behav movie:", frame_idx)
    print("Total frames written to cropped behav movie:", out_idx)
    print("Done, wrote behaviour video with TTL overlay.")

convert_sbx_tiff ( fr"H:\vision_restored\EC_GCaMP6s_25\20251218\cell+2.8\cell+2.8_000_010")

clahe = cv2.createCLAHE(clipLimit=7.0, tileGridSize=(8,8))
for recording_path in [ fr"H:\vision_restored\EC_GCaMP6s_25\20251218\cell+2.8\cell+2.8_000_010"]:
    # convert_sbx_tiff(recording_path)
    chunks_to_mp4(recording_path, fps=20, vmin=None, vmax=None, smooth_window = 5)

    # behaviour
    behav_with_ttl(
        recording_path
    )

# p = fr'H:\test\grat_000_003'
# convert_sbx_tiff(p)

#
# # e.g. event_id == 1 means stimulus is present
# stim_frames = frame_idx[event_id == 1]
#
# ttl_path = os.path.join(recording_path, f'{recording_path.split(os.sep)[-2]}.mat')
# event_id = np.squeeze(loadmat(ttl_path)['info']['event_id'][0][0])
# self.dat_subject[day][recording]['ttl_data'] = np.squeeze(loadmat(ttl_path)['info']['frame'][0][0])  # [event_id == 1]
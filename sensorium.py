from imports import *
import cv2

# save_path = r'H:\sensorium\dat1\images_png'
# path = r"K:\sensorium\train_normal_1_test_shuffled\train_normal_1_test_shuffled"
path = r'K:\sensorium\stims\EC_LB_16\train_hi_1_test_shuffled'
save_path = r'K:\sensorium\stims\EC_LB_16\chunk_stims_p105'

###### shuffle the images, partition them into 4 different chunks, and save each chunk to a folder
n_ims = len(os.listdir(path)) - 1 #-1 to accunt for the metadata file
indices = np.arange(n_ims)
# np.random.shuffle(indices)         # shuffles in-place
chunks = np.array_split(indices, 5) # split indices into 4 equal (ish) arrays > each array has indices of images

# img_dir = r'H:\sensorium\dat1\images'
# list of file names that end in png in the path folder
files_sorted = sorted(
    [f for f in os.listdir(path) if f.lower().endswith('.png')],
    key=lambda x: int(os.path.splitext(x)[0]))

for i_chunk, chunk in enumerate(chunks):

    print(f'Processing chunk {i_chunk} of {len(chunks)-1}')

    chunk_save_path = fr'{save_path}\chunk{i_chunk}'
    os.makedirs(chunk_save_path, exist_ok=True)

    im_names_chunk = np.array(files_sorted)[chunk] # only grab images in that chunk

    for image_name in im_names_chunk:
        im_path = os.path.join(path, image_name)

        # read image as grayscale
        im = cv2.imread(im_path, cv2.IMREAD_GRAYSCALE)

        if im is None:
            print(f"Could not read {im_path}")
            continue

        img_up = cv2.resize(
            im,
            (1920, 1080),
            interpolation=cv2.INTER_LINEAR
        )

        # save using same filename
        plt.imsave(
            os.path.join(chunk_save_path, image_name),
            img_up,
            cmap="gray"
        )

for folder_chunk in os.listdir(save_path):
    folder_path = os.path.join(save_path, folder_chunk)

    if os.path.isdir(folder_path):
        for i in range(6):
            os.makedirs(os.path.join(folder_path, str(i)), exist_ok=True)

for chunk in chunks:
    im_chunk = ims[chunk]

    for img in ims[chunk]: # array is shape (size_chunk, x_dim, y_dim)
        # upsampling each image to (1080 x 1920) pixels to match monitor dimensions
        img_up = cv2.resize(
            img,
            (1920, 1080),  # (width, height)
            interpolation=cv2.INTER_LINEAR
        )

        plt.imsave(
            os.path.join(save_path, f"{image}.png"),
            img_up,
            cmap="gray"
        )

[len(c) for c in chunks]

for image in os.listdir(path):

    # store each image (size is 144 x 256)
    img = np.squeeze(np.load(os.path.join(path, image)))

    # upsampling each image to (1080 x 1920) pixels to match monitor dimensions
    img_up = cv2.resize(
        img,
        (1920, 1080),  # (width, height)
        interpolation=cv2.INTER_LINEAR
    )

    plt.imsave(
        os.path.join(save_path,f"{image}.png"),
        img_up,
        cmap="gray"
    )



img = np.squeeze(np.load(os.path.join(path, '964.npy')))
# plt.figure(figsize = (16,9))
# plt.imshow(np.squeeze(img), cmap = 'gray')
# plt.show()

# upsampling each image to (1080 x 1920) pixels to match monitor dimensions
img_up = cv2.resize(
    img,
    (1920, 1080),                    # (width, height)
    interpolation=cv2.INTER_LINEAR
)

# plt.figure(figsize = (16,9))
# plt.imshow(np.squeeze(img_up), cmap = 'gray')
# plt.show()

ims = np.zeros((len(os.listdir(path)), 144, 256))
for i, image in enumerate(os.listdir(path)):
    ims[i] = np.squeeze(np.load(os.path.join(path, image)))

for chunk in chunks:
    im_chunk = ims[chunk]

tiers = np.load(r'H:\sensorium\dat1\tiers.npy')

validation_indices = np.where(tiers == "validation")[0]
train_indices = np.where(tiers == "train")[0]
test_indices = np.where(tiers == "test")[0]

print(len(validation_indices), len(train_indices), len(test_indices))

###########

import os
import numpy as np
import hashlib
from collections import defaultdict

img_dir = r'H:\sensorium\dat1\images'

hash_to_files = defaultdict(list)

files_sorted = sorted(os.listdir(img_dir), key=lambda x: int(x[:-4]))

for fname in files_sorted:

    i = int(fname[:-4])

    if i in train_indices:

        img = np.squeeze(np.load(os.path.join(img_dir, fname)))

        # hash raw bytes (fast + exact)
        h = hashlib.sha256(img.tobytes()).hexdigest()
        hash_to_files[h].append(fname)
num_images = sum(len(v) for v in hash_to_files.values())
num_unique = len(hash_to_files)

print(f"Total images: {num_images}")
print(f"Unique images: {num_unique}")
print(f"Repeated images: {num_images - num_unique}")


img_dir = r'H:\sensorium\dat1\images'
files_sorted = sorted(os.listdir(img_dir), key=lambda x: int(x[:-4]))
test_ims = np.zeros((len(test_indices), 144, 256))

count_idx = 0
for fname in files_sorted:

    i = int(fname[:-4]) #image number

    if i in test_indices:
        test_ims[count_idx] = np.squeeze(np.load(os.path.join(img_dir, fname)))
        count_idx +=1

for i, im in enumerate(test_ims):
    num_repeat = 0
    for comp_im in test_ims:
        if np.array_equal(im, comp_im):
            num_repeat +=1

    if num_repeat != 10:
        print(i, num_repeat)



#################################################
import numpy as np
import os

suite2p_path = r'I:\sensorium\data\EC_LB_12\20260618\chunk0\chunk0_175_000\experiments\suite2p functional\plane0'
rois_to_delete = [750]

files = [
    "F.npy",
    "Fneu.npy",
    "iscell.npy",
    "spks.npy",
    "stat.npy",
]

for file in files:
    arr = np.load(os.path.join(suite2p_path,file), allow_pickle=True)
    arr = np.delete(arr, rois_to_delete, axis=0)
    np.save(os.path.join(suite2p_path,file), arr, allow_pickle=True)

#
# stat = np.load(fr"{path_to_suite2p}\stat.npy", allow_pickle=True)
# F = np.load(fr"{path_to_suite2p}\F.npy", allow_pickle=True)
# iscell = np.load(fr"{path_to_suite2p}\iscell.npy", allow_pickle=True)
# spks = np.load(fr"{path_to_suite2p}\spks.npy", allow_pickle=True)
# fneu = np.load(fr"{path_to_suite2p}\Fneu.npy", allow_pickle=True)

# fneu0 = np.delete(fneu, rois_to_delete,0)
# F0 = np.delete(F, rois_to_delete,0)
# stat0 = np.delete(stat, rois_to_delete,0)
# iscell0 = np.delete(iscell, rois_to_delete,0)
# spks0 = np.delete(spks, rois_to_delete,0)
#
# np.save(fr"{path_to_suite2p}\fneu.npy", fneu0, allow_pickle=True)
# np.save(fr"{path_to_suite2p}\F.npy", F0 ,allow_pickle=True)
# np.save(fr"{path_to_suite2p}\iscell.npy", iscell0, allow_pickle=True)
# np.save(fr"{path_to_suite2p}\stat.npy", stat0, allow_pickle=True)
# np.save(fr"{path_to_suite2p}\spks.npy", spks0, allow_pickle=True)


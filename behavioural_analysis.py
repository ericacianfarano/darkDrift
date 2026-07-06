from helpers import *

def load_dlc_pose_estimates(path_to_behav, likelihood_threshold):
    '''
    :param path_to_behav: path to behavioural folder, ex: r'G:\vision_restored\data\EC_RD1_04\20240626\brs\brs_000_000\behav'
    :param likelihood_threshold (float): flt between 0 and 1. DLC likelihood values below likelihood_threshold will be naned out
    :return: poses (array): DLC pose estimates of shape n_timepoints x n_poses
        [top_pupilX, top_pupilY, right_pupilX, right_pupilY, bottom_pupilX, bottom_pupilY, left_pupilX, left_pupilY,]
    '''
    dlc_file_name = [file for file in os.listdir(path_to_behav) if file.endswith('.csv')][0]
    dlc_file_path = os.path.join(path_to_behav, dlc_file_name)

    # skip the first 2 rows, and the first column
    pose_estimation_csv = np.genfromtxt(dlc_file_path, delimiter=',', skip_header=1)[2:, 1:]

    # extract indices where liklihood values are less than the set threshold
    likelihoods = pose_estimation_csv[:, 2::3]  # columns of the likelihood values
    subthreshold_likelihood = np.where(np.all(likelihoods < likelihood_threshold, axis=1))[0]

    # delete likelihood values from poses > delete every third columns starting from the second column
    poses = np.delete(pose_estimation_csv, np.s_[2::3], axis=1)

    # nan out values that are less than the likelihood threshold
    poses[subthreshold_likelihood, :] *= np.nan

    return poses

def fit_ellipse(data):
    """
    Fit an ellipse to given points data.

    Parameters:
    - data: numpy array of shape (n_timepoints, 8) where columns are [topX, topY, rightX, rightY, bottomX, bottomY, leftX, leftY].

    Returns:
    - x_center, y_center, width, height, angle: Parameters of the fitted ellipse.
    """
    # Extract points from data array
    topX, topY = data[:, 0], data[:, 1]
    rightX, rightY = data[:, 2], data[:, 3]
    bottomX, bottomY = data[:, 4], data[:, 5]
    leftX, leftY = data[:, 6], data[:, 7]

    # Center of the ellipse
    x_center = (leftX + rightX) / 2
    y_center = (topY + bottomY) / 2

    # Horizontal and vertical diameters
    height = np.sqrt((rightX - leftX)**2 + (rightY - leftY)**2)
    width = np.sqrt((topX - bottomX)**2 + (topY - bottomY)**2)

    # Angle of rotation (assuming you have the angle or calculate it)
    angle = np.arctan2(topY - bottomY, topX - bottomX)

    return x_center, y_center, width, height, angle

def plot_ellipse(ax, x_center, y_center, width, height, angle, points=None, color = 'red'):
    """
    Plot an ellipse on a given axis.

    Parameters:
    - ax: Matplotlib axis object where the ellipse will be plotted.
    - x_center, y_center, width, height, angle: Parameters of the ellipse.
    - points: Optional, numpy array of shape (4, 2) containing the top, bottom, left, and right points coordinates.

    Returns:
    - None
    """
    # Create the ellipse
    ellipse = Ellipse((x_center, y_center), width, height, angle=np.degrees(angle), edgecolor=color, fc='None')

    # Add the ellipse to the plot
    ax.add_patch(ellipse)

    # Plot the points if provided
    if points is not None:
        for i in range(4):
            ax.scatter(points[i*2], points[(i*2)+1], color=color)

    ax.set_aspect('equal', adjustable='box')


def moving_average(data, window_size = 5):
    '''
    Smoothing to remove noise > moving average filter

    Move the window over the data one point at a time, calculating the average of the data points within the window at each step.
    :param data:
    :param window_size:
    :return:
    '''

    # np.convolve convolves the data array with the array of ones of length 'window_size', normalizing by window_size > calcualte avg of each window
    # 'valid' > so output array only has values where window fully overlaps with input data > to avoid boundary effects
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

def interpolate_nans(signal):
    """Interpolate NaN values in a signal."""
    nans = np.isnan(signal)
    not_nans = ~nans
    indices = np.arange(len(signal))
    interpolator = interp1d(indices[not_nans], signal[not_nans], bounds_error=False, fill_value="extrapolate")
    signal[nans] = interpolator(indices[nans])
    return signal

def analyze_frequency (signal, sampling_rate):
    '''
    Fast Fourier Transform (FFT) > frequency analysis to identify dominant frequencies (i.e., repetitive movements)
    :param signal:
    :param sampling_rate: n fps, acquisition rate, number of frames per second
    :return:
        - freq (array): > array of frequency values
        - np.abs(fourier) (array) > magnitude of each frequency component, which represents the amplitude of that frequency in the signal
    '''

    n = len(signal)                         # N samples in the signal array
    fourier = fft(signal)                   # Compute FFT of the signal > complex array that represents amplitude & phase of each frequency component of the signal
    freq = fftfreq(n, 1 / sampling_rate)    # array of frequencie values that corresponds to the FFT result/components in 'fourier''
    return freq, np.abs(fourier)            # frequencies & magnitude of the FFT result, respectively.

def interpolate_nans_1d(signal):
    """Interpolate NaN values in a signal."""
    nans = np.isnan(signal)
    not_nans = ~nans
    indices = np.arange(len(signal))
    interpolator = interp1d(indices[not_nans], signal[not_nans], bounds_error=False, fill_value="extrapolate")
    signal[nans] = interpolator(indices[nans])
    return signal



def plot_ellipse_on_frame(frame, x_center, y_center, width, height, angle, color=(0, 0, 255), thickness=2):
    """
    Draw an ellipse on an image frame using OpenCV.

    Parameters:
    - frame: The image frame where the ellipse will be drawn.
    - x_center, y_center: Center of the ellipse.
    - width, height: Width and height of the ellipse.
    - angle: Rotation angle of the ellipse in degrees.
    - color: Color of the ellipse (BGR format).
    - thickness: Thickness of the ellipse outline.

    Returns:
    - frame: The image frame with the ellipse drawn on it.
    """
    # Convert width and height from ellipse parameters to OpenCV format
    center = (int(x_center), int(y_center))
    axes = (int(width / 2), int(height / 2))
    angle = int(angle)  # Convert angle to integer

    # Draw the ellipse on the frame
    cv2.ellipse(frame, center, axes, angle, 0, 360, color, thickness)
    cv2.circle(frame, center, radius=1, color=(255, 255, 0), thickness=1)

    return frame

def ellipse_behav_video(path_to_folder, poses, start_idx, end_idx):
    animal, day, folder_name = path_to_folder.split('\\')[3], path_to_folder.split('\\')[4], path_to_folder.split('\\')[-2]
    input_file = os.path.join(path_to_folder, f'{folder_name}_eye.mj2')
    output_file = os.path.join(path_to_folder, f'ellipse_{animal}-{day}-{folder_name}_eye.avi')

    if not os.path.exists(input_file):
        print(f'Input path {input_file} does not exist')
    else:
        print(f'Creating {output_file} video')

        # Open the MJ2 video file
        cap = cv2.VideoCapture(input_file)

        # Get the width and height of frames
        width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
        height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
        fps = cap.get(cv2.CAP_PROP_FPS)

        fourcc = cv2.VideoWriter_fourcc(*'XVID')
        out = cv2.VideoWriter(output_file, fourcc, fps, (width, height))

        # since 'poses' is only the data after the first ttl and before the last ttl, we want to index 'poses'/ellipse data differently than we need to index the raw frames from the video
        frame_num = 0
        within_frame_num = 0

        x_center, y_center, width, height, angle = fit_ellipse(poses)

        x_center_proc = interpolate_nans_1d(x_center)
        y_center_proc = interpolate_nans_1d(y_center)
        width_proc = interpolate_nans_1d(width)
        height_proc = interpolate_nans_1d(height)
        angle_proc = interpolate_nans_1d(angle)

        while (cap.isOpened()):
            ret, frame = cap.read()
            if ret:

                if (frame_num >= start_idx) and (frame_num < end_idx):
                    # plot the ellipse on the frame
                    frame_with_ellipse = plot_ellipse_on_frame(frame, x_center_proc[within_frame_num], y_center_proc[within_frame_num], width_proc[within_frame_num], height_proc[within_frame_num], angle_proc[within_frame_num])
                    out.write(frame_with_ellipse)
                    within_frame_num += 1

                frame_num += 1

            else:
                break

        # Release everything if job is finished
        cap.release()
        out.release()

    cv2.destroyAllWindows()



def nystagmus (object, poses):
    fps = object.fps
    x_center, y_center,width, height, angle = fit_ellipse(poses)
    x_center_interp = interpolate_nans_1d(x_center)
    y_center_interp = interpolate_nans_1d(y_center)

    plt.figure()
    plt.scatter(x_center_interp, y_center_interp, c=np.arange(x_center.shape[0]), cmap='plasma')
    plt.plot(x_center_interp, y_center_interp, c='black', alpha=0.2)
    plt.show()

    plt.figure()
    plt.scatter(np.arange(x_center.shape[0]), x_center_interp, c='red', alpha=0.5)
    plt.title('x pos')
    plt.show()

    plt.figure(figsize = (7,3))
    plt.plot(np.arange(x_center.shape[0]), y_center_interp, c='blue', alpha=0.5)
    plt.title('y pos')
    plt.show()

    fig = plt.figure(figsize = (7,3))
    ax = plt.axes(projection='3d')
    ax.scatter(x_center_interp, y_center_interp, np.arange(x_center.shape[0]), c=np.arange(x_center.shape[0]), cmap='plasma')
    plt.show()

    # # Calculate velocities
    dt = 1/fps # seconds
    velocity_x = np.diff(x_center_interp) / dt  # velocity components in x
    velocity_y = np.diff(y_center_interp) / dt  # velocity components in y
    velocity = np.sqrt(velocity_x**2 + velocity_y**2) # velocity magnitude

    # do fourier on velocities
    frequency_x, fourier_x = analyze_frequency(velocity_x, fps)
    frequency_y, fourier_y = analyze_frequency(velocity_y, fps)
    frequency, fourier = analyze_frequency(velocity, fps)
    ###########
    plt.figure(figsize=(12, 9))
    plt.subplot(3, 1, 1)
    plt.plot(frequency_x, fourier_x)
    plt.title('Frequency Analysis of X Velocity')
    plt.xlabel('Frequency (Hz)')
    plt.xlim([0,10])
    plt.ylabel('Amplitude')
    plt.yscale('log')
    ###########
    plt.subplot(3, 1, 2)
    plt.plot(frequency_y, fourier_y)
    plt.title('Frequency Analysis of Y Velocity')
    plt.xlabel('Frequency (Hz)')
    plt.xlim([0, 10])
    plt.ylabel('Amplitude')
    plt.yscale('log')
    ##########
    plt.subplot(3, 1, 3)
    plt.plot(frequency, fourier)
    plt.title('Frequency Analysis of Velocity')
    plt.xlabel('Frequency (Hz)')
    plt.xlim([0, 10])
    plt.ylabel('Amplitude')
    plt.yscale('log')
    plt.tight_layout()
    plt.show()

    # do fourier on position
    frequency_x, fourier_x = analyze_frequency(x_center_interp, fps)
    frequency_y, fourier_y = analyze_frequency(y_center_interp, fps)
    ###########
    plt.figure(figsize=(12, 6))
    plt.subplot(2, 1, 1)
    plt.plot(frequency_x, fourier_x)
    plt.title('Frequency Analysis of X Position')
    plt.xlabel('Frequency (Hz)')
    plt.xlim([0,10])
    plt.ylabel('Amplitude')
    plt.yscale('log')
    ###########
    plt.subplot(2, 1, 2)
    plt.plot(frequency_y, fourier_y)
    plt.title('Frequency Analysis of Y Position')
    plt.xlabel('Frequency (Hz)')
    plt.xlim([0, 10])
    plt.ylabel('Amplitude')
    plt.yscale('log')
    plt.tight_layout()
    plt.show()

#nystagmus (data_object, data_object.dat['EC_RD1_02']['20240627']['brs_000_000']['poses'])





from imports import *
from helpers import *
from suite2p_preprocessing import *

animals = ['EC_dark_01', 'EC_dark_03']
path = r'I:\dark_drift\data'
data_object = batchProcessing(animals, path, ntheta = 8, fps = 20, deconvolved = True)

# SPONTANEOUS ANALYSIS ACROSS GROUPS

'''coactive_epochs = True
chunk_size = 3 # n frames

var_explained = {'Control': [], 'RD1':[]}
eigenvals_slope = {'Control': [], 'RD1':[]}
corr_values = {'Control': [], 'RD1':[]}
# var_explained = {'Control': [], 'RD1':[], 'GNAT':[]}
# eigenvals_slope = {'Control': [], 'RD1':[], 'GNAT':[]}

for i_animal, animal in enumerate([animal for animal in animals_days.keys() if 'GNAT' not in animal]):

    group = 'Control' if 'GCaMP6s' in animal else animal.split('_')[1]

    # list of tuples, with each entry being (day, subfile)
    if group != 'GNAT':
        days_recordings = [(day, subfile) for day in animals_days[animal] for subfile in data_object.dat[animal][day].keys() if ('SFxTF' in subfile) and ('big' in subfile)][:2]
    else:
        days_recordings = [(day, subfile) for day in animals_days[animal] for subfile in
                           data_object.dat[animal][day].keys() if ('big' in subfile)][:2]

    for i, (day, subfile) in enumerate(days_recordings):

        # shape n_Cells, n_timepoints
        spon_arr = np.squeeze(data_object.dat[animal][day][subfile]['zscored_responses_ttls']['Wait_1'])

        if coactive_epochs:
            n_full_chunks = spon_arr.shape[1] // chunk_size

            # shape > (n_trials, n_features, n_timepoints_per_epoch) (n_epochs, n_cells, n_timepoints)
            #data_trials = np.array([spon_arr[:,i:i + chunk_size] for i in range(0, n_full_chunks * chunk_size, chunk_size)]).reshape (n_full_chunks, -1)
            data_trials_all = np.array([spon_arr[:, i:i + chunk_size] for i in range(0, n_full_chunks * chunk_size, chunk_size)])

            # only take epochs that have at least 'cell_threshold' co-active cells that each go above 'response_threshold' > (n_epochs, n_cells, n_timepoints)
            valid_epochs = filter_active_epochs(data_trials_all, response_threshold=3, cell_threshold = 5)
            #plot_raster(data_trials_all, valid_epochs, n_trials_to_plot=60, n_cells_to_plot=400, threshold=4)

            # then average over each epoch (time) to get average response> (n_epochs, n_cells)
            data_trials = data_trials_all[valid_epochs].mean (axis = -1)

        else:
            data_trials = spon_arr.T
        # Normalize each trial's activity pattern (vector) to unit length (L2 norm)
        # each trial (each row) is normalized to unit length (its L2 norm is 1)
        norms = np.linalg.norm(data_trials, axis=1, keepdims=True)
        data_trials_normalized = data_trials / norms

        # performing PCA & projection
        # explained variance ratio > eigenvalues
        pca, principal_components, explained_variance_ratio = perform_pca (data_trials_normalized, n_components = 50)

        var_explained[group].append(explained_variance_ratio)
        eigenvals_slope[group].append(decay_eigenspectra(explained_variance_ratio))

        corr_matrix = np.corrcoef(data_trials.T)
        pca = PCA(n_components=2)
        reduced = pca.fit_transform(corr_matrix)
        plt.figure(figsize=(5, 5))
        colors = {'Control': 'black', 'RD1': 'red', 'GNAT': 'green'}
        plt.scatter(reduced[:, 0], reduced[:, 1], alpha = 0.6, s = 15,color = colors[group])
        plt.xlabel('PC 1')
        plt.ylabel('PC 2')
        plt.title(f'PCA of pearson correlation {group}')
        plt.show()

for group in var_explained:
    var_explained[group] = np.array(var_explained[group])
    eigenvals_slope[group] = np.array(eigenvals_slope[group])
    #corr_matrix[group] = np.array(corr_matrix[group])

variance_explained (var_explained, log = True)
plot_slope_eigenvals(eigenvals_slope)



event_rates = {'Control': [], 'RD1':[]}
event_amplitudes = {'Control': [], 'RD1':[]}
event_durations = {'Control': [], 'RD1':[]}
dominant_frequencies_per_cell = {'Control': [], 'RD1':[]}
power_spectrum = {'Control': [], 'RD1':[]}

coactive_counts = {'Control': [], 'RD1':[]}
z_thresh = 3
for i_animal, animal in enumerate([animal for animal in animals_days.keys() if 'GNAT' not in animal]):
    group = 'Control' if 'GCaMP6s' in animal else animal.split('_')[1]

    # list of tuples, with each entry being (day, subfile)
    if group != 'GNAT':
        days_recordings = [(day, subfile) for day in animals_days[animal] for subfile in data_object.dat[animal][day].keys() if ('SFxTF' in subfile) and ('big' in subfile)][:2]
    else:
        days_recordings = [(day, subfile) for day in animals_days[animal] for subfile in
                           data_object.dat[animal][day].keys() if ('big' in subfile)][:2]

    for i, (day, subfile) in enumerate(days_recordings):
        spon_arr = np.squeeze(data_object.dat[animal][day][subfile]['zscored_responses_ttls']['Wait_1']) # shape n_Cells, n_timepoints
        #event_rates[group].append(compute_event_rates_per_cell(spon_arr, z_threshold=z_thresh))

        rates, all_amplitudes, all_durations = compute_event_features_per_cell(spon_arr, z_threshold=z_thresh)
        event_rates[group].append(rates)
        event_amplitudes[group].append(all_amplitudes)
        event_durations[group].append(all_durations)

        active = spon_arr > z_thresh  # shape (n_cells, n_timepoints) > binarized
        coactive_counts[group].append(active.sum(axis=0))  # number of active cells at each timepoint

        # fourier
        n_cells, n_timepoints = spon_arr.shape
        sampling_rate = data_object.fps
        # Compute FFT
        fft_vals = np.fft.rfft(spon_arr, axis=1)  # Only positive frequencies
        fft_freqs = np.fft.rfftfreq(n_timepoints, d=1 / sampling_rate)
        # Compute power spectrum
        power = np.abs(fft_vals) ** 2  # shape: (n_cells, n_freqs)
        # Get dominant frequency per cell
        dominant_freqs = fft_freqs[np.argmax(power[:, 1:], axis=1) + 1]  # skip DC (index 0)
        dominant_frequencies_per_cell[group].append(dominant_freqs)
        power_spectrum[group].append(power.mean(axis = 0))

for group in event_rates:
    event_rates[group] = np.array([item for sublist in event_rates[group] for item in sublist])
    coactive_counts[group] = np.array([item for sublist in coactive_counts[group] for item in sublist])
    event_amplitudes[group] = np.array([item for sublist in event_amplitudes[group] for item in sublist])
    event_durations[group] = np.array([item for sublist in event_durations[group] for item in sublist])
    dominant_frequencies_per_cell[group] = np.array([item for sublist in dominant_frequencies_per_cell[group] for item in sublist])
    power_spectrum[group] = np.array([item for sublist in power_spectrum[group] for item in sublist])
    #eigenvals_slope[group] = np.array(eigenvals_slope[group])

plot_spon_event_rate(event_rates, timepoints=False)
coactive_cells_per_frame(coactive_counts, timepoints=False)
plot_spon_event_properties(event_amplitudes, event_durations, timepoints=False)
plot_fourrier(dominant_frequencies_per_cell,power_spectrum, timepoints=False)'''

### SPONTANEOUS ANALYSIS ACROSS TIME



#variance_explained (var_explained, log = True, timepoints = True)

##### stponaneous activity rates at different time points

event_rates = {'Control_0': [],'Control_1': [], 'RD1_0':[], 'RD1_1':[]}
coactive_counts = {'Control_0': [],'Control_1': [], 'RD1_0':[], 'RD1_1':[]}
event_amplitudes= {'Control_0': [],'Control_1': [], 'RD1_0':[], 'RD1_1':[]}
event_durations= {'Control_0': [],'Control_1': [], 'RD1_0':[], 'RD1_1':[]}
dominant_frequencies_per_cell= {'Control_0': [],'Control_1': [], 'RD1_0':[], 'RD1_1':[]}
z_thresh = 3
for i_animal, animal in enumerate(animals_days.keys()):#[animal for animal in animals_days.keys() if 'GNAT' not in animal]):
    g = 'Control' if 'GCaMP6s' in animal else animal.split('_')[1]
    days_recordings = [(day, subfile) for day in animals_days[animal] for subfile in data_object.dat[animal][day].keys() if ('SFxTF' in subfile) and ('big' in subfile)][:2]

    for i, (day, subfile) in enumerate(days_recordings):
        group = g + '_' +str(i)
        spon_arr = np.squeeze(data_object.dat[animal][day][subfile]['zscored_responses_ttls']['Wait_1']) # shape n_Cells, n_timepoints
        #event_rates[group].append(compute_event_rates_per_cell(spon_arr, z_threshold=z_thresh))

        rates, all_amplitudes, all_durations = compute_event_features_per_cell(spon_arr, z_threshold=z_thresh)
        event_rates[group].append(rates)
        event_amplitudes[group].append(all_amplitudes)
        event_durations[group].append(all_durations)

        active = spon_arr > z_thresh  # shape (n_cells, n_timepoints) > binarized
        coactive_counts[group].append(active.sum(axis=0))  # number of active cells at each timepoint

        # fourier
        n_cells, n_timepoints = spon_arr.shape
        sampling_rate = data_object.fps
        # Compute FFT
        fft_vals = np.fft.rfft(spon_arr, axis=1)  # Only positive frequencies
        fft_freqs = np.fft.rfftfreq(n_timepoints, d=1 / sampling_rate)
        # Compute power spectrum
        power = np.abs(fft_vals) ** 2  # shape: (n_cells, n_freqs)
        # Get dominant frequency per cell
        dominant_freqs = fft_freqs[np.argmax(power[:, 1:], axis=1) + 1]  # skip DC (index 0)
        dominant_frequencies_per_cell[group].append(dominant_freqs)

for group in event_rates:
    event_rates[group] = np.array([item for sublist in event_rates[group] for item in sublist])
    coactive_counts[group] = np.array([item for sublist in coactive_counts[group] for item in sublist])
    event_amplitudes[group] = np.array([item for sublist in event_amplitudes[group] for item in sublist])
    event_durations[group] = np.array([item for sublist in event_durations[group] for item in sublist])
    dominant_frequencies_per_cell[group] = np.array([item for sublist in dominant_frequencies_per_cell[group] for item in sublist])
    #eigenvals_slope[group] = np.array(eigenvals_slope[group])

plot_spon_event_rate(event_rates, timepoints=True)
coactive_cells_per_frame(coactive_counts, timepoints=True)
plot_spon_event_properties(event_amplitudes, event_durations, timepoints=True)
plot_dominant_frequencies(dominant_frequencies_per_cell, timepoints=True)


################################################################################################




def plot_pca_projection(animal, principal_components, n_repeats, n_orientations, trial_type = 'evoked'):
    """
    Plots the projection of trials onto the first two principal components,
    with color indicating stimulus orientation.

    Args:
        animal (str): Animal identifier.
        principal_components (np.ndarray): Principal components of shape (n_trials, n_components).
        n_repeats (int): Number of stimulus repeats.
        n_orientations (int): Number of stimulus orientations.
        trial_type (str): either 'evoked' or 'spontaneous'
    """

    if trial_type == 'evoked':
        # Assign colors according to orientation
        orientations = np.linspace(0, 315, n_orientations)  # 0, 45, 90, ..., 315
        orientation_ids = np.repeat(np.arange(n_orientations), n_repeats)
        colors = orientations[orientation_ids]

        fig = plt.figure(figsize=(5, 4))
        #ax = fig.add_subplot(111, projection = '3d')
        ax = fig.add_subplot(111)

        scatter = ax.scatter(
            principal_components[:, 0],
            principal_components[:, 1],
            #principal_components[:, 2],
            c=colors,
            cmap='plasma',
            alpha=0.8
        )

        # color bar showing orientation
        cbar = plt.colorbar(scatter, ax=ax, ticks=orientations)
        cbar.set_label('Orientation (degrees)')
        cbar.set_ticks(orientations)
        cbar.set_ticklabels([f'{int(o)}°' for o in orientations])


    elif trial_type == 'spontaneous':

        fig = plt.figure(figsize=(5, 4))
        #ax = fig.add_subplot(111, projection='3d')
        ax = fig.add_subplot(111)

        scatter = ax.scatter(
            principal_components[:, 0],
            principal_components[:, 1],
            #principal_components[:, 2],
            c=np.arange(principal_components.shape[0]),
            cmap='plasma',
            alpha=0.8
        )

    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')
    ax.set_title(f'{animal} - PCA Projection \n of single-trial grating-evoked activity')

    plt.show()
    return fig



def valid_epochs(array, coactive_epochs=True, chunk_size=3):
    '''
    spontaneous analysis: slope of eigenspectrum across time
    :param obj:
    :param coactive_epochs:
    :param chunk_size:
    :return:
    '''

    # array should be z-scored trackes

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

# then average over each epoch (time) to get average response> (n_epochs, n_cells)
data_trials = data_trials_all [valid_epochs].mean(axis=-1)

print(f'{animal}, {day}, {data_trials.shape[0] / data_trials_all.shape[0] * 100}% of trials have >4 coactive cells')

# Normalize each trial's activity pattern (vector) to unit length (L2 norm)
# each trial (each row) is normalized to unit length (its L2 norm is 1)
norms = np.linalg.norm(data_trials, axis=1, keepdims=True)
data_trials_normalized = data_trials / norms

pca, principal_components, explained_variance_ratio = perform_pca (data_trials_normalized.T, n_components = 2)

plot_pca_projection('EC_GCaMP6s_09', principal_components, 5, 8, trial_type = 'spontaneous')

# EMBEDDING THE POPULATION VECTORS
#calcualte each population vector
arr = data_object.dat['EC_LBc_03'].dat_subject['20251017']['grat']['param_matrix_whole_zscore'].mean(axis = -1) # shape n_repeats, n_orientations, n_cells
pop_vec = arr.reshape(-1, arr.shape[-1]) # shape n_repeats * n_orientations , n_cells

# normalize each vector to unit length
norms = np.linalg.norm(pop_vec, axis=1, keepdims=True)
pop_vec_norm = pop_vec / np.maximum(norms, 1e-12)

pca = PCA(n_components=3)
ori_labels = data_object.dat['EC_LBc_03'].dat_subject['20251017']['grat']['orientations'][np.tile(np.arange(arr.shape[1]), arr.shape[0])]

# # Example if you have 8 orientations spaced 45° apart
# n_orientations = 8
# ori_deg = np.arange(0, 360, 360 / n_orientations)
# # Collapse 0°–180°, 45°–225°, etc.
# ori_collapsed = (ori_deg % 180)
# # Repeat for each trial
# ori_labels = np.tile(ori_collapsed, arr.shape[0])

scores = pca.fit_transform(pop_vec_norm)
plt.figure(figsize=(5,5))
plt.scatter(scores[:,0], scores[:,1],c=ori_labels, cmap='plasma', s=40, alpha = 0.7)
plt.axis('equal')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('population vecs')
plt.colorbar(label='Orientation (deg)')
plt.show()

############# doing the analyssi8s across all aninmals

def make_ori_projection_group (obj, group, postnatal_day, orientation = True):

    all_animals_traces = []
    # all_animals_ori_labels = []

    # put all animals in group and postnatal day into one giant array
    for animal in [a for a in data_object.dat if group + '_' in a]:

        for i, day in enumerate([d for d in data_object.dat[animal].dat_subject.keys() if calculate_animal_age(obj.animal_dobs[animal], d) == postnatal_day]):

            # traces_epochs, ori_epochs = frames_epochs(experiment_object, epoch_len=20, grating=True)

            experiment_object = obj.dat[animal].dat_subject[day]['grat']

            #if each epoch is an actual trial > embedding population vectors
            traces_epochs = experiment_object['param_matrix_whole_zscore'][...,20:].mean(axis=-1)  # shape n_repeats, n_orientations, n_cells, remove first second (20 frames) where screen grey
            if orientation:
                ori_epochs = experiment_object['orientations'][np.tile(np.arange(traces_epochs.shape[1]), traces_epochs.shape[0])] #% 180
            else:
                ori_epochs = experiment_object['orientations'][
                    np.tile(np.arange(traces_epochs.shape[1]), traces_epochs.shape[0])]   % 180
            traces_epochs = traces_epochs.reshape(-1, traces_epochs.shape[-1])  # shape n_repeats * n_orientations , n_cells
            #remove each epoch’s global mean across cells
            traces_epochs = traces_epochs - traces_epochs.mean(axis=1, keepdims=True)

            norms = np.linalg.norm(traces_epochs, axis=1, keepdims=True)
            normalized_traces = traces_epochs / np.maximum(norms, 1e-12)

            all_animals_traces.append(normalized_traces)
            # all_animals_ori_labels.append(ori_epochs)

    # all_animals_traces
    # pseudo-stacking population vectors, as if all animals were recorded simulataneously. should work as the grating data is organized/sorted
    # each array in the dic should be shape (n_ori x n_trials) x (n_total_cells_across_animals)

    all_animals_traces = np.hstack(all_animals_traces)
    # all_animals_ori_labels = np.hstack(all_animals_ori_labels)

    pca = PCA(svd_solver='full') #keep all the PCs

    scores = pca.fit_transform(all_animals_traces)  # (n_epochs, 3)
    epoch_idx = np.arange(len(scores))  # 0..n_epochs-1
    norm_t = epoch_idx / epoch_idx.max()

    return ori_epochs, scores, pca

def plot_umap_animals(obj, group, orientation = True):
    n_days = 4
    ncols = int(np.ceil(np.sqrt(n_days)))
    nrows = int(np.ceil(n_days / ncols))

    if 'EB' in group:
        keys = ['P70', 'P77', 'P84', 'P91']
        animal_groups = ['EB', 'EBc']
    elif 'LB' in group:
        keys = ['P105', 'P112', 'P119', 'P126']
        animal_groups = ['LB', 'LBc']
    print(animal_groups)
    for g in animal_groups:
        print(g)
        fig_umap, ax_umap = plt.subplots(nrows=1, ncols=n_days, figsize=(4 * n_days, 3))
        for i, postnatal_day in enumerate(keys):

            ori_epochs, scores, pca_obj = make_ori_projection_group(obj, g, postnatal_day, orientation = orientation)

            explained_var_ratio, cum, n_keep = scree_info(pca_obj, thresh=0.80)
            feats = scores[:, :n_keep] # shape pop_vec, n_components kept
            plot_umap(ax_umap[i], feats, colors=ori_epochs, title=f'P{postnatal_day}', n_neighbors=15, min_dist=0.1, metric='cosine',random_state=0)

        unique_oris = np.arange(0, 360, 45)
        cbar = fig_umap.colorbar(
            plt.cm.ScalarMappable(
                cmap=plt.cm.Paired,
                norm=plt.Normalize(vmin=0, vmax=len(unique_oris) - 1)
            ),
            ax=ax_umap[i],
            shrink=0.8,
            pad=0.02
        )

        # --- Properly label the discrete ticks ---
        cbar.set_ticks(np.arange(len(unique_oris)))
        cbar.set_ticklabels(unique_oris.astype(int))
        cbar.set_label('Orientation (°)')

        fig_umap.suptitle(f'{g} UMAP')

        folder_path = os.path.join(os.path.dirname(data_object.path), 'figures', f'umap_gratings_all_animals')
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        plt.savefig(os.path.join(folder_path,f'umap {g} {orientation * "orientation"}.svg'))
        plt.savefig(os.path.join(folder_path,f'umap {g} {orientation * "orientation"}.png'))

        #return ori_epochs
        plt.show(block=False)  # shows both figures without blocking

plot_umap_animals(data_object,'EB', orientation = 1)
plot_umap_animals(data_object,'EB', orientation = 0)
plot_umap_animals(data_object,'LB', orientation = 1)
plot_umap_animals(data_object,'LB', orientation = 0)

for animal in data_object.dat.keys():
    if data_object.dat[animal].dat_subject:
        plot_ori_projection(data_object, animal, grat = True, n_components_plot=3)

#####################################################


for day in data_object.dat[animal].dat_subject.keys():
    plot_ori_projection(data_object.dat[animal].dat_subject[day]['grat'], grat = True)
    #plot_ori_projection(data_object.dat[animal].dat_subject[day]['spon1'], grat = False)



plt.figure(figsize=(8,8))
plt.scatter(scores[:,1], scores[:,2], c=ori_epochs, cmap='plasma', s=15 + 35 * norm_t, alpha=0.2 + 0.8 * norm_t)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('PC1'); plt.ylabel('PC2'); plt.title('Epoch-wise population patterns (PC1–PC2)')
plt.colorbar(label='Orientation (deg)')
plt.show()

'''
# ---- 3D scatter ----
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(
    scores[:, 1], scores[:, 2], scores[:, 3],
    c=ori_epochs, cmap='plasma',
    s=15 + 35 * norm_t, alpha=0.2 + 0.8 * norm_t
)
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
ax.set_title('Epoch-wise population patterns (3D PCA projection)')
# Colorbar
cb = fig.colorbar(sc, ax=ax, pad=0.1)
cb.set_label('Orientation (deg)')
plt.show()
'''

def compute_event_features_per_cell(data, z_threshold=3, fps=20):
    """
    Compute spontaneous event rates, amplitudes, and durations from z-scored traces.

    Parameters:
    - data: array of shape (n_cells, n_timepoints)
    - z_threshold: threshold to detect peaks
    - fps: imaging frame rate (frames per second)

    Returns:
    - rates: array of shape (n_cells,) with event rates in Hz > peaks per second
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
        rates[cell_idx] = len(peaks) / duration # num_peaks*fps / n_frames = #peaks/sec

        # amplitudes (z-scored peak heights)
        all_amplitudes.extend(properties["peak_heights"])

        # durations (half-max width in seconds)
        widths = peak_widths(trace, peaks, rel_height=0.5)[0] / fps
        all_durations.extend(widths)

    #shape (n_cells),
    return rates, all_amplitudes, all_durations

def event_features_groups (obj):

    groups = np.unique([a.split('_')[1] for a in obj.list_animals])

    # {'EB': {animal 1: [], animal 2:[]}],'EBc': {animal 1: [] ...}, 'LB':{animal 1: [] ...}, 'LBc': {animal 1: [] ...}}
    d_rates = {g: {} for g in groups}
    d_amplitudes = {g: {} for g in groups}
    d_durations = {g: {} for g in groups}

    for i_animal, animal in enumerate(obj.dat):

        group = animal.split('_')[1]
        days_recordings = [(day, subfile) for day in obj.dat[animal].dat_subject for subfile in
                           obj.dat[animal].dat_subject[day] if 'spon' in subfile]

        d_animal_rates = {}
        d_animal_amplitudes = {}
        d_animal_durations = {}

        for i, (day, subfile) in enumerate(days_recordings):

            timepoint = calculate_animal_age(obj.animal_dobs[animal], day) #+ '_s' * int(subfile[-1])
            #if timepoint == 'P70' or timepoint == 'P105': # if animals are seeing for the firs time (P70 or P105), add the _s for spon1
            timepoint += '_s' * int(subfile[-1])

            spon_arr = data_object.dat[animal].dat_subject[day][subfile]['zscored_traces']
            rates, amplitudes, durations = compute_event_features_per_cell(spon_arr, z_threshold=3, fps=20)

            d_animal_rates[timepoint] =  np.array(rates).mean()
            d_animal_amplitudes[timepoint]= np.array(amplitudes).mean()
            d_animal_durations[timepoint]= np.array(durations).mean()

        d_rates[group][animal] = d_animal_rates
        d_amplitudes[group][animal] = d_animal_amplitudes
        d_durations[group][animal] = d_animal_durations

    return (d_rates, d_amplitudes, d_durations)


def plot_event_rates(obj):

    triplex = event_features_groups(obj)

    labels = ['spontaneous event rates (Hz)', 'spontaneous event amplitudes', 'spontaneous event durations (Hz)']

    for i, metric in enumerate(triplex): # rates, amplitude, duration

        ncols = 2
        fig, ax = plt.subplots(
            nrows=int(np.ceil(len(metric)/ncols)),
            ncols=ncols, figsize=(12, 7), sharey=True
        )
        ax = ax.ravel()
        colours = {'EB': 'blue', 'EBc': 'black', 'LB': 'red', 'LBc':'green'}

        max_days = max(
            len(days)
            for group_dict in metric.values()
            for days in group_dict.values())

        if len(metric) == 1:
            ax = [ax]

        for i_group, group in enumerate(metric):
            group_days = sorted(
                {d for days in metric[group].values() for d in days.keys()},
                key=sort_key
            )

            ax[i_group].plot(group_days, [np.nan]*len(group_days), alpha=0)

            for animal in metric[group]:
                sorted_items = sorted(metric[group][animal].items(),
                                      key=lambda x: sort_key(x[0]))
                x = [k for k, _ in sorted_items]
                y = [v for _, v in sorted_items]
                ax[i_group].plot(x, y, marker='o', color='gray', alpha=0.4)

            # group mean ± SEM
            values_by_day = {}
            for animal, days in metric[group].items():
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
            ax[i_group].set_ylabel(labels[i])
            ax[i_group].set_title(f'{group} group')
            ax[i_group].grid(axis='y', linestyle='--', alpha=0.5)
            ax[i_group].tick_params(axis='x', rotation=45)
            ax[i_group].set_xlim([-0.5, max_days])

        plt.tight_layout()
        plt.show()

################# PROJECTING SPON INTO EVOKED SPACE ##################


def project_spon (experiment_object, grat = True):

    if grat:
        # traces_epochs, ori_epochs = frames_epochs(experiment_object, epoch_len=20, grating=True)

        #if each epoch is an actual trial > embedding population vectors
        traces_epochs = experiment_object['param_matrix_whole_zscore'][...,20:].mean(axis=-1)  # shape n_repeats, n_orientations, n_cells, remove first second (20 frames) where screen grey
        ori_epochs = experiment_object['orientations'][np.tile(np.arange(traces_epochs.shape[1]), traces_epochs.shape[0])] #% 180
        traces_epochs = traces_epochs.reshape(-1, traces_epochs.shape[-1])  # shape n_repeats * n_orientations , n_cells
        # print(traces_epochs.shape, ori_epochs.shape)

        # # collapsing across repeats > still collapses in the PCA but not in the UMAP
        # traces_epochs = experiment_object['param_matrix_whole_zscore'][...,20:].mean(axis=(0,-1))  # shape n_orientations, n_cells, remove first second (20 frames) where screen grey
        # #print(traces_epochs.shape)
        # ori_epochs = experiment_object['orientations']#[np.tile(traces_epochs.shape[0])]# % 180
        # #print(ori_epochs)
        # traces_epochs = traces_epochs.reshape(-1, traces_epochs.shape[-1])  # shape n_orientations , n_cells
        # #print(traces_epochs.shape, ori_epochs.shape)

    else:
        traces_epochs =  frames_epochs(experiment_object, epoch_len=20, grating=False)

    #remove each epoch’s global mean across cells
    traces_epochs = traces_epochs - traces_epochs.mean(axis=1, keepdims=True)

    norms = np.linalg.norm(traces_epochs, axis=1, keepdims=True)
    normalized_traces = traces_epochs / np.maximum(norms, 1e-12)

    pca = PCA(svd_solver='full') #keep all the PCs

    scores = pca.fit_transform(normalized_traces)  # (n_epochs, 3)
    #print(scores.shape)
    epoch_idx = np.arange(len(scores))  # 0..n_epochs-1
    norm_t = epoch_idx / epoch_idx.max()

    return traces_epochs, ori_epochs, scores, pca

animal = 'EC_EBc_07'
for day in data_object.dat[animal].dat_subject:

     # shape n_repeats, n_orientations, n_cells (remove first second (20 frames) where screen grey)
    evoked_arr =  data_object.dat[animal].dat_subject[day]['grat']['param_matrix_whole_zscore'][..., 20:].mean(axis=-1)

    ori_epochs = data_object.dat[animal].dat_subject[day]['grat']['orientations'][
                      np.tile(np.arange(evoked_arr.shape[1]), evoked_arr.shape[0])] % 180

    # shape n_repeats * n_orientations , n_cells
    evoked_arr = evoked_arr.reshape(-1, evoked_arr.shape[-1])  # shape n_repeats * n_orientations , n_cells

    # remove global mean per sample (across cells)
    evoked_arr -= evoked_arr.mean(axis=1, keepdims=True) # shape (n_repeats * n_orientations, n_cells) remove global mean across cells for each pop vector, (80, 1) broadcasts across cells

    # remove mean activity per neuron (across epochs) (for PCA centering) > shape n_repeats * n_orientations , n_cells
    mu_evoked = evoked_arr.mean(axis=0, keepdims=True) #subtract mean to center data, shape (1, n_cells)
    evoked_arr -= evoked_arr.mean(axis=0, keepdims=True) # shape (n_repeats * n_orientations, n_cells) remove mean for each neuron, (1, n_cells)

    norms_evoked = np.linalg.norm(evoked_arr, axis=1, keepdims=True)
    normalized_evoked = evoked_arr / np.maximum(norms_evoked, 1e-12)
    pca = PCA(n_components=None, svd_solver='full')  # keep all the PCs
    scores_evoked = pca.fit_transform(normalized_evoked)  # (n_epochs, n_pcs),

    for stim in [s for s in data_object.dat[animal].dat_subject[day] if 'spon' in s]:
        spon_arr = data_object.dat[animal].dat_subject[day][stim]['zscored_traces'].T # shape timepoints, n_cells

        # chunk spon_arr into epochs > spontaneous samples/epochs should match those of evoked samples (which are trial-averaged response windows)
        epoch_length = 20 # match evoked averaging window
        nT = (spon_arr.shape[0] // epoch_length) * epoch_length # num frames that fit perfectly into window length
        spon_arr_epochs = spon_arr[:nT].reshape(-1, epoch_length, spon_arr[:nT].shape[1]).mean(axis=1)  # shape (n_epochs, n_cells)

        # remove global mean per sample (across cells)
        spon_arr_epochs -= spon_arr_epochs.mean(axis=1, keepdims=True)  # (n_epochs, 1) > mean is shape (n_epochs, 1)

        # center spon using evoked mean (mean per neuron), shape (n_epochs, n_cells)
        spon_arr_epochs -= mu_evoked

        norms_spon = np.linalg.norm(spon_arr_epochs, axis=1, keepdims=True)
        normalized_spon = spon_arr_epochs / np.maximum(norms_spon, 1e-12)

        scores_spon = pca.transform(normalized_spon)

    xmin = min(scores_spon[:,0].min(), scores_evoked[:,0].min())
    xmax = max(scores_spon[:,0].max(), scores_evoked[:,0].max())
    ymin = min(scores_spon[:,1].min(), scores_evoked[:,1].min())
    ymax = max(scores_spon[:,1].max(), scores_evoked[:,1].max())
    plt.figure()
    plt.hist2d(scores_spon[:,0], scores_spon[:,1], bins=60, range=[[xmin, xmax], [ymin, ymax]])
    plt.colorbar(label="spont count")
    # plt.scatter(scores_evoked[:,0], scores_evoked[:,1], s=30, color="white", edgecolors="k", zorder=3)
    sc = plt.scatter(
         scores_evoked[:, 1],
         scores_evoked[:, 2],
         c=ori_epochs,
         cmap="plasma",
         s=35,
         zorder=3
     )
    cb = plt.colorbar(sc)
    cb.set_label("orientation (deg, mod 180)")
    plt.xlim(xmin, xmax); plt.ylim(ymin, ymax)
    plt.xlabel("PC1 (evoked)"); plt.ylabel("PC2 (evoked)")
    plt.title(day)
    plt.show()



#### projecting into spon space  ####
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

animal = 'EC_EBc_07'

for day in data_object.dat[animal].dat_subject:

    # ---------------------------
    # EVOKED: single-trial patterns (repeat x orientation), already z-scored
    # ---------------------------
    # grat_arr: (n_repeats, n_orientations, n_cells, n_time=100)
    grat_arr = data_object.dat[animal].dat_subject[day]['grat']['param_matrix_whole_zscore']

    # average evoked response over frames 20:100 (skip first 20 grey frames)
    # evoked_arr: (n_repeats, n_orientations, n_cells)
    evoked_arr = grat_arr[..., 20:].mean(axis=-1)

    ori_epochs =  data_object.dat[animal].dat_subject[day]['grat']['orientations'][
        np.tile(np.arange(evoked_arr.shape[1]), evoked_arr.shape[0])] % 180

    # flatten to (n_trials, n_cells), where n_trials = n_repeats * n_orientations
    evoked_arr = evoked_arr.reshape(-1, evoked_arr.shape[-1])

    # ---------------------------
    # SPONT: define spontaneous manifold (PCA on spontaneous activity)
    # ---------------------------
    for stim in [s for s in data_object.dat[animal].dat_subject[day] if 'spon' in s]:

        # spon_arr: (n_timepoints, n_cells), already z-scored
        spon_arr = data_object.dat[animal].dat_subject[day][stim]['zscored_traces'].T

        # make spontaneous "samples" as window-averages over same duration as evoked window (100-20 = 80 frames)
        epoch_length = 40
        nT = (spon_arr.shape[0] // epoch_length) * epoch_length
        spon_arr_epochs = spon_arr[:nT].reshape(-1, epoch_length, spon_arr.shape[1]).mean(axis=1)  # (n_epochs_spon, n_cells)

        # ---------------------------
        # Match the paper's top panel:
        # - Fit PCA on spontaneous activity
        # - Project (unit-length) single-trial evoked patterns into spont PC space
        #
        # Centering: use spontaneous per-neuron mean as the reference mean for BOTH spont and evoked.
        # (No per-sample "global mean across cells" subtraction here, because the paper text doesn't mention it.)
        # ---------------------------
        mu_spon = spon_arr_epochs.mean(axis=0, keepdims=True)          # (1, n_cells)

        # center spont (for PCA)
        spon_centered = spon_arr_epochs - mu_spon

        # spon_centered-=spon_arr_epochs.mean(axis=1, keepdims=True)

        # PCA on spont
        pca = PCA(n_components=None, svd_solver='full')
        scores_spon = pca.fit_transform(spon_centered)                 # (n_epochs_spon, K)

        # center evoked using the SAME spont mean
        evoked_centered = evoked_arr - mu_spon

        # normalize evoked patterns to unit length (paper explicitly says this)
        norms_evoked = np.linalg.norm(evoked_centered, axis=1, keepdims=True)
        normalized_evoked = evoked_centered / np.maximum(norms_evoked, 1e-12)

        # project evoked into spont PC space
        scores_evoked = pca.transform(normalized_evoked)               # (n_trials, K)

        # ---------------------------
        # Plot TOP PANEL: evoked points over spont PC1/PC2 space
        # ---------------------------
        xmin = min(scores_spon[:, 0].min(), scores_evoked[:, 0].min())
        xmax = max(scores_spon[:, 0].max(), scores_evoked[:, 0].max())
        ymin = min(scores_spon[:, 1].min(), scores_evoked[:, 1].min())
        ymax = max(scores_spon[:, 1].max(), scores_evoked[:, 1].max())
        # xmin = min(scores_evoked[:, 0].min())
        # xmax = max(scores_evoked[:, 0].max())
        # ymin = min(scores_evoked[:, 1].min())
        # ymax = max(scores_evoked[:, 1].max())


        plt.figure()
        plt.hist2d(scores_spon[:, 0], scores_spon[:, 1], bins=60, range=[[xmin, xmax], [ymin, ymax]])
        # plt.colorbar(label="spont count")
        # plt.scatter(scores_evoked[:, 0], scores_evoked[:, 1], s=25, color="white", edgecolors="k", linewidths=0.5, zorder=3)
        sc = plt.scatter(
            scores_evoked[:, 1],
            scores_evoked[:, 2],
            c=ori_epochs,
            cmap="plasma",
            s=35,
            zorder=3
        )
        cb = plt.colorbar(sc)
        cb.set_label("orientation (deg, mod 180)")
        # plt.xlim(xmin, xmax)
        # plt.ylim(ymin, ymax)
        plt.xlabel("PC1 (spont PCA space)")
        plt.ylabel("PC2 (spont PCA space)")
        plt.title(f"{animal} | {day} | {stim} : evoked in spont PC space")
        plt.show()


def chunk_data (animal, day, window_size = 20):

    # for day in data_object.dat[animal].dat_subject:

    # arr: (n_cells, n_timepoints)
    grat_arr = data_object.dat[animal].dat_subject[day]['grat']['zscored_traces']

    # get which orientaiton being shown in each frame. we can average later as if there is a 'mix' it will nan it out automatically
    ori_per_frame, segments = build_ori_per_frame(
        dict_stim_ttls= data_object.dat[animal].dat_subject[day]['grat']['dict_stim_ttls'],
        dict_stim_names= data_object.dat[animal].dat_subject[day]['grat']['dict_stim_names'],
        n_timepoints=grat_arr.shape[-1],
        collapse_180=True,  # set True if you want pure orientation (0≡180, etc.)
        fps = 20
    )

    nT = (grat_arr.shape[1] // window_size) * window_size # n timepoints to keep for analysis (how many fit cleanly into window)
    grat_arr_epochs = grat_arr[:, :nT].reshape(grat_arr.shape[0], -1, window_size) # (n_cells, n_epochs, n_timepoints)

    ori_per_frame_epochs = ori_per_frame[:nT].reshape(-1, window_size).mean(axis= -1) # shape n_epochs
    # # only take epochs that have at least 'cell_threshold' co-active cells that each go above 'response_threshold' > (n_epochs, n_cells, n_timepoints)
    # valid_epochs = filter_active_epochs(data_trials_all, response_threshold=3, cell_threshold=5)
    # # plot_raster(data_trials_all, valid_epochs, n_trials_to_plot=60, n_cells_to_plot=400, threshold=4)
    #
    # # then average over each epoch (time) to get average response> (n_epochs, n_cells)
    # data_trials = data_trials_all[valid_epochs].mean(axis=-1)

    whole_dat = grat_arr_epochs.copy()
    whole_dat_labels = ori_per_frame_epochs #np.zeros((whole_dat.shape[1]))+2 # shape (n_epochs), 2 if grat

    # TRY: filtering out silent epochs, averaging over epochs

    for stim in [s for s in data_object.dat[animal].dat_subject[day] if 'spon' in s]:

        # arr: (n_cells, n_timepoints)
        spon_arr = data_object.dat[animal].dat_subject[day][stim]['zscored_traces']

        nT = (spon_arr.shape[1] // window_size) * window_size  # n timepoints to keep for analysis (how many fit cleanly into window)
        spon_arr_epochs = spon_arr[:, :nT].reshape(spon_arr.shape[0], -1,window_size)  # (n_cells, n_epochs, n_timepoints)

        spon_labels = np.zeros((spon_arr_epochs.shape[1])) + -100*(-1+int(stim[-1]))  # shape (n_epochs), 2 if grat, 0 if spon0, 1 if spon1

        whole_dat = np.concatenate((whole_dat, spon_arr_epochs), axis=1)
        whole_dat_labels = np.concatenate((whole_dat_labels, spon_labels), axis=0)

    return whole_dat.mean(axis = -1), whole_dat_labels


def plot_umap(ax, feats, colors=None, title='UMAP', n_neighbors=15, min_dist=0.1, metric='cosine', random_state=0):
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=random_state)
    emb = reducer.fit_transform(feats)  # (n_epochs, 2)
    sc = ax.scatter(emb[:, 0], emb[:, 1], alpha = 0.7, c=colors if colors is not None else 'tab:blue',cmap = 'Paired', s=12)
    ax.set_xlabel('UMAP-1')
    ax.set_ylabel('UMAP-2')
    ax.set_title(title)
    return sc


def plot_chunks (animal, window_size):

    fig_umap, ax_umap = plt.subplots(nrows=1, ncols=4, figsize=(16,4))

    for i_day, day in enumerate(data_object.dat[animal].dat_subject):

        postnatal_day = calculate_animal_age(data_object.animal_dobs[animal], day)

        whole_dat, whole_dat_labels = chunk_data(animal,day, window_size)

        # mu_global = whole_dat.mean(axis=0, keepdims=True)
        mu_cell = whole_dat.mean(axis=1, keepdims=True)
        # whole_dat-=mu_global
        whole_dat -= mu_cell
        # print(mu_cell.shape)

        pca = PCA(n_components=None, svd_solver='full')
        scores = pca.fit_transform(whole_dat.T)  # (n_epochs, n_PCs)

        norms_scores = np.linalg.norm(scores, axis=1, keepdims=True)
        normalized_scores = scores / np.maximum(norms_scores, 1e-12)
        #
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection="3d")
        # sc = ax.scatter(
        #     normalized_scores[:, 0],  # PC1
        #     normalized_scores[:, 1],  # PC2
        #     normalized_scores[:, 2],  # PC3
        #     c=whole_dat_labels,
        #     cmap="viridis",
        #     s=35,
        #     alpha=0.5
        # )
        # cb = fig.colorbar(sc, ax=ax)
        # cb.set_label("stim")
        # ax.set_xlabel("PC1")
        # ax.set_ylabel("PC2")
        # ax.set_zlabel("PC3")
        # plt.show()

        explained_var_ratio, cum, n_keep = scree_info(pca, thresh=0.80)
        feats = normalized_scores[:, :n_keep]  # shape pop_vec, n_components kept

        s=plot_umap(ax_umap[i_day], feats, colors=whole_dat_labels, title=f'P{postnatal_day}', n_neighbors=15, min_dist=0.1,
                  metric='cosine', random_state=0)

    fig_umap.colorbar(s)
    fig_umap.suptitle(animal)

    folder_path = os.path.join(os.path.dirname(data_object.path), 'figures', f'umap_spon_grat')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'umap_spon_grat {animal}.svg'))
    plt.savefig(os.path.join(folder_path, f'umap_spon_grat {animal}.png'))

    plt.show()

for animal in data_object.dat.keys():
    plot_chunks (animal, 20)

############## trying the umap analysis for the LBc across days to see if we can get consistent grating mapping

def build_ori_per_frame(dict_stim_ttls, dict_stim_names, n_timepoints, collapse_180=False, fps = 20):
    """
    Returns:
      ori_per_frame: (n_timepoints,) float array, deg per frame, NaN where no grating
      segments: list of (start, end, deg, key) for sanity checks
    """
    ori_per_frame = np.full(n_timepoints, np.nan, dtype=float)
    segments = []

    # helper: parse degrees from a name like 'grating_270_SF0.02_TF1.00'
    def parse_deg(name):
        m = re.search(r'grating_(-?\d+(?:\.\d+)?)', name)
        if not m:
            return None
        return float(m.group(1))

    # iterate keys in chronological order (numeric suffix if present)
    def key_order(k):
        m = re.search(r'_(\d+)$', k)
        return int(m.group(1)) if m else 10_000

    for key in sorted(dict_stim_ttls.keys(), key=key_order):
        ttls = dict_stim_ttls[key]   # list of [start, end] for this repeat/type
        names = dict_stim_names.get(key, None)

        # skip non-grating blocks (e.g., 'Wait_0')
        if names is None or key.lower().startswith('wait'):
            continue

        if len(ttls) != len(names):
            raise ValueError(f"Length mismatch for {key}: {len(ttls)} intervals vs {len(names)} names")

        for (start, end), name in zip(ttls, names):
            deg = parse_deg(name)
            if deg is None:
                # not a grating (or unparsable) → leave as NaN
                continue
            if collapse_180: # merge direction pairs togehter. 0 & 180 degrees are the same
                deg = deg % 180.0

            # assume [start, end) (end exclusive); clip to array bounds just in case
            s = max(0, int(start)+(fps)) # exclude the first static second > ttl starts at static stim
            e = min(n_timepoints, int(start)+(fps*4)) # only include first 3 seconds of movement
            if e > s:
                ori_per_frame[s:e] = deg
                segments.append((s, e, deg, key))

    return ori_per_frame, segments

def chunk_data_grat(animal, day, window_size=20):
    # for day in data_object.dat[animal].dat_subject:

    # arr: (n_cells, n_timepoints)
    grat_arr = data_object.dat[animal].dat_subject[day]['grat']['zscored_traces']

    # get which orientaiton being shown in each frame. we can average later as if there is a 'mix' it will nan it out automatically
    ori_per_frame, segments = build_ori_per_frame(
        dict_stim_ttls=data_object.dat[animal].dat_subject[day]['grat']['dict_stim_ttls'],
        dict_stim_names=data_object.dat[animal].dat_subject[day]['grat']['dict_stim_names'],
        n_timepoints=grat_arr.shape[-1],
        collapse_180=True,  # set True if you want pure orientation (0≡180, etc.)
        fps=20
    )

    nT = (grat_arr.shape[
              1] // window_size) * window_size  # n timepoints to keep for analysis (how many fit cleanly into window)
    grat_arr_epochs = grat_arr[:, :nT].reshape(grat_arr.shape[0], -1, window_size)  # (n_cells, n_epochs, n_timepoints)

    ori_per_frame_epochs = ori_per_frame[:nT].reshape(-1, window_size).mean(axis=-1)  # shape n_epochs
    # # only take epochs that have at least 'cell_threshold' co-active cells that each go above 'response_threshold' > (n_epochs, n_cells, n_timepoints)
    # valid_epochs = filter_active_epochs(data_trials_all, response_threshold=3, cell_threshold=5)
    # # plot_raster(data_trials_all, valid_epochs, n_trials_to_plot=60, n_cells_to_plot=400, threshold=4)
    #
    # # then average over each epoch (time) to get average response> (n_epochs, n_cells)
    # data_trials = data_trials_all[valid_epochs].mean(axis=-1)

    # turn 'inactive' epochs into a different colour
    active_indices = active_epochs_idx(grat_arr_epochs, response_threshold=2.5, min_active_cells=2) # shape n_epochs

    ori_per_frame_epochs[~active_indices & ~np.isnan(ori_per_frame_epochs)]=-100 # if epoch is 'inactive' and it is also not already nan, make it -100 (mask is shape n_epochs)

    # ori_per_frame_epochs[np.isnan(ori_per_frame_epochs)] = -100 # put all leftover nans at -200

    return grat_arr_epochs.mean(axis=-1), ori_per_frame_epochs


def plot_umap(ax, feats, colors=None, title='UMAP', n_neighbors=15, min_dist=0.1, metric='cosine', random_state=0, colour_map = 'discrete'):
    '''
    :param colour_map: either 'discrete' or 'continuous'
    :return:
    '''
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=random_state)
    emb = reducer.fit_transform(feats)  # (n_epochs, 2)
    if colour_map == 'discrete':
        c = 'Paired'
    elif colour_map == 'continuous':
        c = 'plasma'
    sc = ax.scatter(emb[:, 0], emb[:, 1], alpha=0.7, c=colors if colors is not None else 'tab:blue', cmap=c,
                    s=12)
    # ax.plot(emb[:, 0], emb[:, 1], alpha=0.5, c='black')
    ax.set_xlabel('UMAP-1')
    ax.set_ylabel('UMAP-2')
    ax.set_title(title)
    return sc


def plot_chunks_grat(animal, window_size):
    fig_umap, ax_umap = plt.subplots(nrows=1, ncols=4, figsize=(16, 4))
    pca = PCA(n_components=None, svd_solver='full')

    for i_day, day in enumerate(data_object.dat[animal].dat_subject):
        postnatal_day = calculate_animal_age(data_object.animal_dobs[animal], day)

        grat_dat, grat_dat_labels = chunk_data_grat(animal, day, window_size)
        # print('total', len(grat_dat_labels))
        # print('not nan', np.sum(~np.isnan(grat_dat_labels)))
        # grat_dat_labels*=np.nan

        # grat_dat = grat_dat[:, ~np.isnan(grat_dat_labels)]
        # grat_dat_labels = grat_dat_labels[~np.isnan(grat_dat_labels)]

        # # mu_global = whole_dat.mean(axis=0, keepdims=True)
        # mu_cell = grat_dat.mean(axis=1, keepdims=True)
        # # whole_dat-=mu_global
        # grat_dat -= mu_cell
        # # print(mu_cell.shape)

        if i_day == 0:
            scores = pca.fit_transform(grat_dat.T)  # (n_epochs, n_PCs)
        else:
            scores = pca.transform(grat_dat.T)  # (n_epochs, n_PCs)

        norms_scores = np.linalg.norm(scores, axis=1, keepdims=True)
        normalized_scores = scores / np.maximum(norms_scores, 1e-12)
        #
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection="3d")
        # sc = ax.scatter(
        #     normalized_scores[:, 0],  # PC1
        #     normalized_scores[:, 1],  # PC2
        #     normalized_scores[:, 2],  # PC3
        #     c=whole_dat_labels,
        #     cmap="viridis",
        #     s=35,
        #     alpha=0.5
        # )
        # cb = fig.colorbar(sc, ax=ax)
        # cb.set_label("stim")
        # ax.set_xlabel("PC1")
        # ax.set_ylabel("PC2")
        # ax.set_zlabel("PC3")
        # plt.show()
        explained_var_ratio, cum, n_keep = scree_info(pca, thresh=0.8)
        # print(n_keep)
        feats = normalized_scores[:, :n_keep]  # shape pop_vec, n_components kept
        # print(feats.shape)
        # grat_dat_labels = feats[:, 2] # np.linalg.norm(feats, axis=1) > color by population magnitude
        s = plot_umap(ax_umap[i_day], feats, colors=grat_dat_labels, title=f'P{postnatal_day}', n_neighbors=20,
                      min_dist=0.1,
                      metric='euclidean', random_state=0, colour_map = 'Discrete')
    fig_umap.colorbar(s)
    fig_umap.suptitle(animal)

    folder_path = os.path.join(os.path.dirname(data_object.path), 'figures', f'umap_whole_grat')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'umap_whole_grat {animal}.svg'))
    plt.savefig(os.path.join(folder_path, f'umap_whole_grat {animal}.png'))

    plt.show()


for animal in [d for d in data_object.dat.keys() if 'EB_' in d]:
    plot_chunks_grat(animal, 10)

############################## chunk animal date across days togehter ##################### and put into PCA #@###############


def chunk_data (animal, window_size = 20):

    whole_dat = np.array([])
    whole_dat_labels = np.array([])

    for i_day, day in enumerate(data_object.dat[animal].dat_subject):

        ###################### if we want to take proper grating epoch #################
        # param_matrix shape (n_repeats, n_orientations, n_cells, n_time=100) > single-trial patterns
        # exclude static grating portion, only take moving grating portion
        # then take mean across time axis (try also taking the max)
        grat_arr = data_object.dat[animal].dat_subject[day]['grat']['param_matrix_whole_zscore'][..., data_object.fps:]

        ori_epochs = data_object.dat[animal].dat_subject[day]['grat']['orientations'][
                         np.tile(np.arange(grat_arr.shape[1]), grat_arr.shape[0])] % 180

        # flatten to (n_trials, n_cells), where n_trials = n_repeats * n_orientations, shape (n_cells, n_epochs or trials (n_repeats*n_orientations), n_timepoints)
        grat_arr_epochs = grat_arr.reshape(grat_arr.shape[2], grat_arr.shape[0]*grat_arr.shape[1], -1)

        ####################### if we want to chunk raw array ########################
        # arr: (n_cells, n_timepoints)
        # grat_arr = data_object.dat[animal].dat_subject[day]['grat']['zscored_traces']
        #
        # # get which orientaiton being shown in each frame. we can average later as if there is a 'mix' it will nan it out automatically
        # ori_per_frame, segments = build_ori_per_frame(
        #     dict_stim_ttls= data_object.dat[animal].dat_subject[day]['grat']['dict_stim_ttls'],
        #     dict_stim_names= data_object.dat[animal].dat_subject[day]['grat']['dict_stim_names'],
        #     n_timepoints=grat_arr.shape[-1],
        #     collapse_180=True,  # set True if you want pure orientation (0≡180, etc.)
        #     fps = 20
        # )
        #
        # nT = (grat_arr.shape[1] // window_size) * window_size # n timepoints to keep for analysis (how many fit cleanly into window)
        # grat_arr_epochs = grat_arr[:, :nT].reshape(grat_arr.shape[0], -1, window_size) # (n_cells, n_epochs, n_timepoints)
        #
        # ori_per_frame_epochs = ori_per_frame[:nT].reshape(-1, window_size).mean(axis= -1) # shape n_epochs
        # # only take epochs that have at least 'cell_threshold' co-active cells that each go above 'response_threshold' > (n_epochs, n_cells, n_timepoints)
        # valid_epochs = filter_active_epochs(data_trials_all, response_threshold=3, cell_threshold=5)
        # # plot_raster(data_trials_all, valid_epochs, n_trials_to_plot=60, n_cells_to_plot=400, threshold=4)
        #
        # # then average over each epoch (time) to get average response> (n_epochs, n_cells)
        # data_trials = data_trials_all[valid_epochs].mean(axis=-1)

        ##################################################################################################

        # turn 'inactive' epochs into a different colour
        # active_indices = active_epochs_idx(grat_arr_epochs, response_threshold=1.4,
        #                                    min_active_cells=2)  # shape n_epochs
        #
        # ori_epochs[~active_indices & ~np.isnan(ori_epochs)] = np.nan  # if epoch is 'inactive' and it is also not already nan, make it -100 (mask is shape n_epochs)
        # print(np.unique(ori_per_frame_epochs))

        if i_day == 0:
            whole_dat = grat_arr_epochs.copy()
            whole_dat_labels = ori_epochs.copy()

            day_labels = np.zeros_like(ori_epochs)
        else:
            print('a')
            # whole_dat = np.concatenate((whole_dat, grat_arr_epochs), axis=1)
            # whole_dat_labels = np.concatenate((whole_dat_labels, ori_epochs), axis=0)
            # day_labels = np.concatenate((day_labels, np.zeros_like(ori_epochs)+i_day), axis = 0)

    # day_labels[np.isnan(whole_dat_labels)] *= np.nan
    #whole dat > shape (n_cells, n_epochs)
    return whole_dat.mean(axis = -1), whole_dat_labels, day_labels


def plot_umap(ax, feats, colors=None, title='UMAP', n_neighbors=15, min_dist=0.1, metric='cosine', random_state=0, colour_map = 'discrete'):
    '''
    :param colour_map: either 'discrete' or 'continuous'
    :return:
    '''
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=random_state)
    emb = reducer.fit_transform(feats)  # (n_epochs, 2)
    if colour_map == 'discrete':
        c = 'Paired'
    elif colour_map == 'continuous':
        c = 'plasma'
    sc = ax.scatter(emb[:, 0], emb[:, 1], alpha=0.7, c=colors if colors is not None else 'tab:blue', cmap=c,
                    s=12)
    # ax.plot(emb[:, 0], emb[:, 1], alpha=0.5, c='black')
    ax.set_xlabel('UMAP-1')
    ax.set_ylabel('UMAP-2')
    ax.set_title(title)
    return sc, emb


def plot_chunks_grat(animal, window_size, umap_metric, color_map = 'grat'):
    '''

    :param animal:
    :param window_size:
    :param umap_metric: either 'cosine' or 'euclidean' (euclidean better for ring)
    :param color_map: either 'grat' or 'day'
    :return:
    '''
    fig_umap, ax_umap = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    pca = PCA(n_components=None, svd_solver='full')

    # grat_dat is shape (n_cells, n_epochs)
    grat_dat, grat_dat_labels, day_labels = chunk_data(animal, window_size)
    # grat_dat_labels*=np.nan

    # grat_dat = grat_dat[:, ~np.isnan(grat_dat_labels)]
    # grat_dat_labels = grat_dat_labels[~np.isnan(grat_dat_labels)]

    # # mu_global = whole_dat.mean(axis=0, keepdims=True)
    # mu_cell = grat_dat.mean(axis=1, keepdims=True)
    # # whole_dat-=mu_global
    # grat_dat -= mu_cell
    # # print(mu_cell.shape)

    # #remove each epoch’s global mean across cells
    grat_dat = grat_dat - grat_dat.mean(axis=1, keepdims=True)
    #
    norms = np.linalg.norm(grat_dat, axis=0, keepdims=True)
    normalized_traces = grat_dat / np.maximum(norms, 1e-12)

    # #remove each epoch’s global mean across cells
    # traces_epochs = traces_epochs - traces_epochs.mean(axis=0, keepdims=True)
    #
    # norms = np.linalg.norm(traces_epochs, axis=1, keepdims=True)
    # normalized_traces = traces_epochs / np.maximum(norms, 1e-12)

    scores = pca.fit_transform(normalized_traces.T)  # (n_epochs, n_PCs)

    print(scores.shape)


    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(
        scores[:, 0],  # PC1
        scores[:, 1],  # PC2
        scores[:, 2],  # PC3
        c=grat_dat_labels,
        cmap="plasma",
        s=35,
        alpha=0.5
    )
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("stim")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")
    plt.show()
    explained_var_ratio, cum, n_keep = scree_info(pca, thresh=0.8)
    feats = scores[:, :n_keep]  # shape pop_vec, n_components kept

    # norms_feats = np.linalg.norm(feats, axis=1, keepdims=True)
    # normalized_feats = feats / np.maximum(norms_feats, 1e-12)

    # grat_dat_labels = feats[:, 2] # np.linalg.norm(feats, axis=1) > color by population magnitude
    if color_map == 'grat':
        labels = grat_dat_labels
    elif color_map=='day':
        labels = day_labels
    s, embedding = plot_umap(ax_umap, feats, colors=labels, n_neighbors=20, min_dist=0.1, metric=umap_metric, random_state=0, colour_map = 'continuous')

    # unique_labels, centroids = centroid_metrics(normalized_feats, grat_dat_labels, day_labels)
    centroids=0
    unique_labels = 0
    # [plt.scatter(centroids[i_label, 0], centroids[i_label, 1], c = 'black', s = 100, alpha = 0.6) for i_label, label in enumerate(unique_labels)]
    fig_umap.colorbar(s)
    fig_umap.suptitle(animal)

    folder_path = os.path.join(os.path.dirname(data_object.path), 'figures', f'umap_whole_grat')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'umap_grat_day {animal}.svg'))
    plt.savefig(os.path.join(folder_path, f'umap_grat_day {animal}.png'))

    plt.show()

    return feats, grat_dat_labels, day_labels, unique_labels, centroids


for animal in [d for d in data_object.dat.keys() if 'LBc_06' in d]:
    feats, grat_dat_labels, day_labels, unique_labels, centroids = plot_chunks_grat(animal, 10, umap_metric = 'cosine', color_map='grat')

    # comparing 1 orientation cluster across days > where second number is the same
    days = np.unique([day for (day, o) in unique_labels])
    orientations_shown = np.unique([o for (day, o) in unique_labels])
    # distance_matrix = np.zeros((len(orientations_shown), len(days), len(days)))
    # for i, orientation in enumerate(orientations_shown):  # iterate through orientations
    #     ori_cent_across_days = np.array([cent
    #                                      for (day, ori), cent in zip(unique_labels, centroids)
    #                                      if ori == orientation])
    #
    #     distance_matrix[i] = cdist(ori_cent_across_days, ori_cent_across_days, metric="euclidean")
    #
    # plt.figure()
    # plt.title(animal)
    # plt.imshow(distance_matrix.mean(axis = 0), cmap= 'plasma', vmin = 0, vmax = 0.8)
    # plt.show()


# comparing distance of the orientation clusters within each day:

def centroid_metrics (features, grat_labels, day_labels):
    '''
    :param features: umap output, shape (n_epochs, n_features)
    :param labels: either
            - grat_dat_labels: shape (n_features)
            - day_labels: shape (n_features)
    :return:
    '''
    # feats shape (n_epochs, n_features)
    # grat_dat_labels shape (n_features)
    # day_labels shape (n_features)

    #gives tuple of (DAY, ORIENTATION CLUSTER)
    unique_labels = [[(int(d), int(g)) for d in np.unique(day_labels) if ~np.isnan(d)] for g in np.unique(grat_labels) if ~np.isnan(g)] # list of unique labels (either day or orientation)
    unique_labels = [pair for sublist in unique_labels for pair in sublist]

    # labels == l: gives boolean mask > (shape n_epochs)
    # we then want to index feat array on this mask to grab only epochs that have that label > (shape n_epochs, feats)
    # then we can average across epochs > (shape feats)
    # then we do this for each label in unique labels > we get 1 centroid per label
    # shapes: label (day/ori), features
    centroids = np.vstack([features[(grat_labels == g) & (day_labels == d)].mean(axis=0) for (d,g) in unique_labels])

    return unique_labels, centroids

###################################################################################


def make_ori_projection (animal, ori = True):

    for i_day, day in enumerate(data_object.dat[animal].dat_subject.keys()):
        experiment_object = data_object.dat[animal].dat_subject[day]['grat']

        # traces_epochs, ori_epochs = frames_epochs(experiment_object, epoch_len=20, grating=True)

        #if each epoch is an actual trial > embedding population vectors
        traces_epochs = experiment_object['param_matrix_whole_zscore'][...,data_object.fps:data_object.fps*4].mean(axis=-1)  # shape n_repeats, n_orientations, n_cells, remove first second (20 frames) where screen grey
        if ori:
            ori_epochs = experiment_object['orientations'][np.tile(np.arange(traces_epochs.shape[1]), traces_epochs.shape[0])] % 180 # shape n_epochs
        else:
            ori_epochs = experiment_object['orientations'][
                             np.tile(np.arange(traces_epochs.shape[1]), traces_epochs.shape[0])]  # shape n_epochs
        traces_epochs = traces_epochs.reshape(-1, traces_epochs.shape[-1])  # shape n_repeats * n_orientations  (n_epochs), n_cells

        #remove each epoch’s global mean across cells
        traces_epochs = traces_epochs - traces_epochs.mean(axis=1, keepdims=True)

        # norms = np.linalg.norm(traces_epochs, axis=1, keepdims=True)
        # normalized_traces = traces_epochs / np.maximum(norms, 1e-12)
        #
        # pca = PCA(svd_solver='full') #keep all the PCs
        #
        # scores = pca.fit_transform(normalized_traces)  # (n_epochs (1 trial, 1 orientation), n_components)
        # #print(scores.shape)
        # epoch_idx = np.arange(len(scores))  # 0..n_epochs-1
        # norm_t = epoch_idx / epoch_idx.max()

        if i_day == 0:
            whole_dat = traces_epochs.copy()
            whole_dat_labels = ori_epochs.copy()
            day_labels = np.zeros_like(ori_epochs)
        else:
            whole_dat = np.concatenate((whole_dat, traces_epochs), axis=0)
            whole_dat_labels = np.concatenate((whole_dat_labels, ori_epochs), axis=0)
            day_labels = np.concatenate((day_labels, np.zeros_like(ori_epochs)+i_day), axis = 0)

    return whole_dat, whole_dat_labels, day_labels

def plot_umap(ax, feats, colors=None, title='UMAP', n_neighbors=15, min_dist=0.1, metric='cosine', random_state=0, colour_map = 'discrete'):
    '''
    :param colour_map: either 'discrete' or 'continuous'
    :return:
    '''
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=random_state)
    emb = reducer.fit_transform(feats)  # (n_epochs, 2)
    if colour_map == 'discrete':
        c = 'Paired'
    elif colour_map == 'continuous':
        c = 'plasma'
    sc = ax.scatter(emb[:, 0], emb[:, 1], alpha=0.7, c=colors if colors is not None else 'tab:blue', cmap=c,
                    s=12)
    # ax.plot(emb[:, 0], emb[:, 1], alpha=0.5, c='black')
    ax.set_xlabel('UMAP-1')
    ax.set_ylabel('UMAP-2')
    ax.set_title(title)
    return sc, emb


def plot_ori_projection(animal, labels = 'grat', ori = True):

    # fig_scatter = plt.figure(figsize=(5 * ncols, 5 * nrows))
    # fig_scree, ax_scree = plt.subplots(figsize=(4,4))
    fig_umap, ax_umap = plt.subplots(nrows = 1,ncols = 1, figsize=(3, 3))

    for i, day in enumerate(data_object.dat[animal].dat_subject.keys()):

        # traces_epochs is shape (n_epochs, n_cells)
        # ori_epochs is shape (n_epochs)
        traces_epochs, ori_epochs, day_labels = make_ori_projection(animal, ori = ori)

        pca = PCA(svd_solver='full') #keep all the PCs
        scores = pca.fit_transform(traces_epochs)  # (n_epochs (1 trial, 1 orientation), n_components)
        # epoch_idx = np.arange(len(scores))  # 0..n_epochs-1
        # norm_t = epoch_idx / epoch_idx.max()

        explained_var_ratio, cum, n_keep = scree_info(pca, thresh=0.80)
        #ax_scree = fig_scree.add_subplot(nrows, ncols, i + 1)
        # plot_scree(ax_scree, explained_var_ratio, cum, n_keep,  title=f'{animal}', label = f'P{postnatal_day}', color = colors_days[i])

        feats = scores[:, :n_keep] # shape pop_vec, n_components kept

        norms = np.linalg.norm(feats, axis=1, keepdims=True)
        feats = feats / np.maximum(norms, 1e-12)


        if labels == 'grat':
            c = ori_epochs
            cmap = 'discrete'
        elif labels == 'day':
            c = day_labels
            cmap = 'continuous'
        plot_umap(ax_umap, feats, colors=c, title=f'UMAP {animal}', n_neighbors=15, min_dist=0.1, metric='cosine', random_state=0,colour_map=cmap)

        # if n_components_plot == 3:
        #     ax_sc = fig_scatter.add_subplot(nrows, ncols, i + 1, projection='3d')
        #     sc = ax_sc.scatter(
        #         scores[:, 0], scores[:, 1], scores[:, 2],
        #         c=(ori_epochs if grat else 'blue'),
        #         cmap=('plasma' if grat else None), s=30
        #     )
        #     ax_sc.set_zlabel('PC3')
        #     ax_sc.view_init(elev=25, azim=40)
        #     ax_sc.set_box_aspect((1, 1, 1))
        # else:
        #     ax_sc = fig_scatter.add_subplot(nrows, ncols, i + 1)
        #     sc = ax_sc.scatter(
        #         scores[:, 0], scores[:, 1],
        #         c=(ori_epochs if grat else 'blue'),
        #         cmap=('plasma' if grat else None), s=30
        #     )
        #
        # ax_sc.set_xlabel('PC1'); ax_sc.set_ylabel('PC2')
        # ax_sc.set_title(f'{animal} — {day} ({"Gratings" if grat else "Spontaneous"})')
        # if grat:
        #     fig_scatter.colorbar(sc, ax=ax_sc, pad=0.05, shrink=0.7, label='Orientation (deg)')
        #
        # if n_components_plot == 3:
        #     ax_sc.set_zlabel('PC3')
        #     ax_sc.view_init(elev=25, azim=40)
        #     ax_sc.set_box_aspect((1, 1, 1))
        #

    # cbar = fig_umap.colorbar(
    #     plt.cm.ScalarMappable(
    #         cmap=plt.cm.Paired,
    #         norm=plt.Normalize(vmin=0, vmax=len(unique_oris) - 1)
    #     ),
    #     ax=ax_umap[i],
    #     shrink=0.8,
    #     pad=0.02
    # )

    # # --- Properly label the discrete ticks ---
    # cbar.set_ticks(np.arange(len(unique_oris)))
    # cbar.set_ticklabels(unique_oris.astype(int))
    # cbar.set_label('Orientation (°)')
    #
    # fig_scatter.tight_layout()
    # fig_scree.tight_layout()
    # fig_umap.suptitle(f'{animal} UMAP')

    # folder_path = os.path.join(os.path.dirname(obj.path), 'figures', f'umap_gratings')
    # if not os.path.exists(folder_path):
    #     os.makedirs(folder_path)
    # plt.savefig(os.path.join(folder_path, f'umap {animal}.svg'))
    # plt.savefig(os.path.join(folder_path, f'umap {animal}.png'))

    # plt.show(block=False)  # shows both figures without blocking

for animal in data_object.dat.keys():
    plot_ori_projection(animal, labels = 'grat', ori = True)
    plot_ori_projection(animal, labels='grat', ori = False)
    plot_ori_projection(animal, labels='day')


#
# def plot_ori_projection(obj, animal, grat=True, n_components_plot = 3):
#     days = list(obj.dat[animal].dat_subject.keys())
#     n_days = len(days)
#     ncols = int(np.ceil(np.sqrt(n_days)))
#     nrows = int(np.ceil(n_days / ncols))
#
#     fig_scatter = plt.figure(figsize=(5 * ncols, 5 * nrows))
#     fig_scree, ax_scree = plt.subplots(figsize=(4,4))
#     fig_umap, ax_umap = plt.subplots(nrows = 1,ncols = n_days, figsize=(4*n_days, 3))
#
#     if n_days == 1:
#         ax_umap = [ax_umap]
#
#     for i, day in enumerate(days):
#
#         postnatal_day = calculate_animal_age(obj.animal_dobs[animal], day)
#         cmap = plt.cm.Paired
#         unique_oris = np.arange(0,360,45)
#
#         cmap_days = plt.cm.plasma
#         colors_days = [cmap_days(i / max(1, n_days - 1)) for i in range(n_days)]
#
#         traces_epochs, ori_epochs, scores, pca_obj = make_ori_projection(
#             obj.dat[animal].dat_subject[day]['grat'], grat=True
#         )
#
#         explained_var_ratio, cum, n_keep = scree_info(pca_obj, thresh=0.80)
#         #ax_scree = fig_scree.add_subplot(nrows, ncols, i + 1)
#         plot_scree(ax_scree, explained_var_ratio, cum, n_keep,  title=f'{animal}', label = f'P{postnatal_day}', color = colors_days[i])
#
#         feats = scores[:, :n_keep] # shape pop_vec, n_components kept
#         plot_umap(ax_umap[i], feats, colors=ori_epochs, title=f'P{postnatal_day}', n_neighbors=15, min_dist=0.1, metric='cosine',random_state=0)
#
#         if n_components_plot == 3:
#             ax_sc = fig_scatter.add_subplot(nrows, ncols, i + 1, projection='3d')
#             sc = ax_sc.scatter(
#                 scores[:, 0], scores[:, 1], scores[:, 2],
#                 c=(ori_epochs if grat else 'blue'),
#                 cmap=('plasma' if grat else None), s=30
#             )
#             ax_sc.set_zlabel('PC3')
#             ax_sc.view_init(elev=25, azim=40)
#             ax_sc.set_box_aspect((1, 1, 1))
#         else:
#             ax_sc = fig_scatter.add_subplot(nrows, ncols, i + 1)
#             sc = ax_sc.scatter(
#                 scores[:, 0], scores[:, 1],
#                 c=(ori_epochs if grat else 'blue'),
#                 cmap=('plasma' if grat else None), s=30
#             )
#
#         ax_sc.set_xlabel('PC1'); ax_sc.set_ylabel('PC2')
#         ax_sc.set_title(f'{animal} — {day} ({"Gratings" if grat else "Spontaneous"})')
#         if grat:
#             fig_scatter.colorbar(sc, ax=ax_sc, pad=0.05, shrink=0.7, label='Orientation (deg)')
#
#         if n_components_plot == 3:
#             ax_sc.set_zlabel('PC3')
#             ax_sc.view_init(elev=25, azim=40)
#             ax_sc.set_box_aspect((1, 1, 1))
#
#
#     cbar = fig_umap.colorbar(
#         plt.cm.ScalarMappable(
#             cmap=plt.cm.Paired,
#             norm=plt.Normalize(vmin=0, vmax=len(unique_oris) - 1)
#         ),
#         ax=ax_umap[i],
#         shrink=0.8,
#         pad=0.02
#     )
#
#     # --- Properly label the discrete ticks ---
#     cbar.set_ticks(np.arange(len(unique_oris)))
#     cbar.set_ticklabels(unique_oris.astype(int))
#     cbar.set_label('Orientation (°)')
#
#     fig_scatter.tight_layout()
#     fig_scree.tight_layout()
#     fig_umap.suptitle(f'{animal} UMAP')
#
#     folder_path = os.path.join(os.path.dirname(obj.path), 'figures', f'umap_gratings')
#     if not os.path.exists(folder_path):
#         os.makedirs(folder_path)
#     plt.savefig(os.path.join(folder_path, f'umap {animal}.svg'))
#     plt.savefig(os.path.join(folder_path, f'umap {animal}.png'))
#
#
#     plt.show(block=False)  # shows both figures without blocking
#
# for animal in data_object.dat.keys():
#     if data_object.dat[animal].dat_subject:
#         plot_ori_projection(data_object, animal, grat = True, n_components_plot=3)
#

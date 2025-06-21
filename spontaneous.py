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



import matplotlib.pyplot as plt

from helpers import *
from imports import *

def plot_avg_response (object):
    '''
    :param object:
    :param recording_day (str): PXX, postnatal day we want to analyze
    :return:
    '''

    response = {k: {} for k in np.unique([a.split ('_')[1] for a in object.list_animals])}

    for animal in object.list_animals:
        group = animal.split('_')[1]

        for day in object.dat[animal].dat_subject.keys():

            postnatal_day = calculate_animal_age(object.animal_dobs[animal], day)
            responsive_indices = object.dat[animal].dat_subject[day]['grat']['thresholded_cells']
            matrix = object.dat[animal].dat_subject[day]['grat']['zscored'][:, :, responsive_indices, :].mean(axis = (0,1,2))  # shape n_repeats, n_ori, n_cells, n_timepoints

            if postnatal_day not in response[group]:
                response[group][postnatal_day] = [matrix]
            else:
                response[group][postnatal_day].append(matrix)

    # PLOTTING
    # ensure consistent sorting
    postnatal_days = sorted(set(day for group_data in response.values() for day in group_data)) # get consistent sorting
    colours = {'ctrl': 'black', 'dark': 'red'}

    fig, ax = plt.subplots (len(postnatal_days),1, figsize = (6,9))
    if len(postnatal_days) == 1:
        ax = [ax]  # Ensure ax is iterable

    for i_day, postnatal_day in enumerate(postnatal_days):
        all_group_data = {}
        for group in response:
            if postnatal_day in response[group]:
                data = np.array(response[group][postnatal_day])  # shape: n_animals, n_timepoints
                all_group_data[group] = data

                # plot individual animals
                for animal in data:
                    ax[i_day].plot(np.arange(data.shape[1]), animal, color=colours.get(group, 'gray'), alpha=0.2)
                # plot group average
                ax[i_day].plot(np.arange(data.shape[1]), data.mean(axis=0), color=colours.get(group, 'gray'), label=group)

        # statistical comparison
        if set(all_group_data.keys()) == {'ctrl', 'dark'}:
            ctrl = all_group_data['ctrl']
            dark = all_group_data['dark']

            if ctrl.shape[1] == dark.shape[1]:  # ensure same n_timepoints
                pvals = np.array([ttest_ind(ctrl[:, t], dark[:, t], equal_var=False).pvalue for t in range(ctrl.shape[1])])
                sig_mask = pvals < 0.05

                # plot asterisks
                y_max = max(ctrl.max(), dark.max())
                for t in np.where(sig_mask)[0]:
                    ax[i_day].text(t, y_max + 0.1, '*', ha='center', va='bottom', fontsize=10, color='k')

        ax[i_day].set_ylim ([min(ctrl.min(), dark.min()), max(ctrl.max(), dark.max()) + 0.25 * max(ctrl.max(), dark.max())])
        ax[i_day].set_ylabel('Response (z-scored)')
        ax[i_day].set_xticks(np.arange(ctrl.shape[1])[::object.fps])
        ax[i_day].set_xticklabels( [int(a/object.fps) for a in (np.arange(ctrl.shape[1])[::object.fps])])
        ax[i_day].axvline(x=object.fps, color='k', linestyle=':', linewidth=1)
        ax[i_day].axvspan(object.fps * 2, object.fps * 5, color='gray', alpha=0.2)
        ax[i_day].set_title(f"Postnatal Day {postnatal_day}")
        ax[i_day].legend()

    ax[-1].set_xlabel("Time (s)")
    plt.tight_layout()

    plt.show()

def plot_avg_response_time (object):
    '''
    :param object:
    :param recording_day (str): PXX, postnatal day we want to analyze
    :return:
    '''

    response = {k: {} for k in np.unique([a.split ('_')[1] for a in object.list_animals])}

    for animal in object.list_animals:
        group = animal.split('_')[1]

        for day in object.dat[animal].dat_subject.keys():

            postnatal_day = calculate_animal_age(object.animal_dobs[animal], day)
            responsive_indices = object.dat[animal].dat_subject[day]['grat']['thresholded_cells']
            matrix = object.dat[animal].dat_subject[day]['grat']['zscored'][:, :, responsive_indices, :].mean(axis = (0,1,2))  # shape n_repeats, n_ori, n_cells, n_timepoints

            if postnatal_day not in response[group]:
                response[group][postnatal_day] = [matrix]
            else:
                response[group][postnatal_day].append(matrix)

    for group in response.keys():
        for day in response[group].keys():
            response[group][day] = np.array(response[group][day]).mean(axis = 0) # average across animals for each day

    # so now response is a dictionary with group / day, with an array of the average response across animals

    # PLOTTING
    # ensure consistent sorting
    postnatal_days = sorted(set(day for group_data in response.values() for day in group_data)) # get consistent sorting
    colors = plt.cm.plasma(np.linspace(0, 0.8, len(postnatal_days)))

    fig, ax = plt.subplots (len(response.keys()),1, figsize = (6,9), sharey = True, sharex = True)

    if len(response.keys()) == 1:
        ax = [ax]  # Ensure ax is iterable

    for i_group, group in enumerate(response.keys()):
        for i_day, postnatal_day in enumerate(response[group].keys()):
            ax[i_group].plot(np.arange(len(response[group][postnatal_day])), response[group][postnatal_day], color=colors[i_day], alpha=0.8, label = postnatal_day)


        #ax[i_group].set_ylim ([min(ctrl.min(), dark.min()), max(ctrl.max(), dark.max()) + 0.25 * max(ctrl.max(), dark.max())])
        ax[i_group].set_ylabel('Response (z-scored)')
        ax[i_group].set_xticks(np.arange(object.fps*5)[::object.fps])
        ax[i_group].set_xticklabels( [int(a/object.fps) for a in (np.arange(object.fps*5)[::object.fps])])
        ax[i_group].axvline(x=object.fps, color='k', linestyle=':', linewidth=1)
        ax[i_group].axvspan(object.fps * 2, object.fps * 5, color='gray', alpha=0.2)
        ax[i_group].set_title(f"{group}")
        ax[i_group].legend()

    ax[-1].set_xlabel("Time (s)")

    plt.suptitle('Average response over time')
    plt.tight_layout()
    plt.show()


def trial_by_trial_reliability (object, recording_day = None):
    '''
    :param object:
    :param recording_day (str): PXX, postnatal day we want to analyze
    :return:
    '''
    reliabilities = {k: [] for k in np.unique([a.split ('_')[1] for a in object.list_animals])}
    for animal in object.list_animals:
        for day in [d for d in object.dat[animal].dat_subject.keys() if calculate_animal_age(object.animal_dobs[animal], d) == recording_day]:

            group = animal.split('_')[1]
            responsive_indices = object.dat[animal].dat_subject[day]['grat']['thresholded_cells']
            matrix = object.dat[animal].dat_subject[day]['grat']['param_matrix_whole'][:,:, responsive_indices,object.fps*1:object.fps*5] # shape n_repeats, n_ori, n_cells, n_timepoints
            r = trial_reliability(matrix).max(axis = 0) # shape (n_cells)
            reliabilities[group].extend(r)

    fig, ax = plt.subplots (nrows = len(reliabilities), ncols = 1, figsize = (4,7), sharey = True, sharex =True)
    for i, group in enumerate(reliabilities):
        ax[i].hist (reliabilities[group], bins = np.linspace(np.concatenate(list(reliabilities.values())).min(),np.concatenate(list(reliabilities.values())).max(),30), density = True)
        ax[i].set_title(group)
        ax[i].set_xlabel('correlation')
    plt.suptitle('trial-by-trial correlation')
    plt.show()

def trial_by_trial_reliability_days (object):
    '''
    :param object:
    :param recording_day (str): PXX, postnatal day we want to analyze
    :return:
    '''

    d = {k: {} for k in np.unique([a.split ('_')[1] for a in object.list_animals])}

    for animal in object.list_animals:
        group = animal.split('_')[1]

        for day in object.dat[animal].dat_subject.keys():
            postnatal_day = calculate_animal_age(object.animal_dobs[animal], day)
            responsive_indices = object.dat[animal].dat_subject[day]['grat']['thresholded_cells']
            matrix = object.dat[animal].dat_subject[day]['grat']['zscored'][:,:, responsive_indices,object.fps*1:object.fps*5] # shape n_repeats, n_ori, n_cells, n_timepoints

            #trial_reliability(matrix) is of shape n_ori, n_cells
            r = trial_reliability(matrix).max(axis = 0) # shape (n_cells)

            if postnatal_day not in d[group]:
                d[group][postnatal_day] = [r]
            else:
                d[group][postnatal_day].append(r)

    # combine all data from all animals into a single array that we can plot
    for group in d.keys():
        for day in d[group].keys():
            d[group][day] = np.concatenate([arr.ravel() for arr in d[group][day]])#.reshape(-1) # average across animals for each day

    postnatal_days = sorted(set(day for group_data in d.values() for day in group_data))
    fig, ax = plt.subplots (nrows = len(d.keys()), ncols = len(postnatal_days), figsize = (10,7), sharex = True, sharey = True)

    for i_group, group in enumerate(d.keys()):
        for i_day, day in enumerate(d[group].keys()):

            ax[i_group, i_day].hist(d[group][day], bins=np.linspace(0,  1, 30), density=True)
            ax[i_group, i_day].set_title(f'{group}, {day}')
            ax[i_group, i_day].set_xlabel('correlation')

    plt.suptitle('trial-by-trial correlation')
    plt.tight_layout()
    plt.show()


def responsiveness (object, recording_day):
    '''
    For 'recording_day' plot the percentage of responsive cells for each group, as a bar plot

    :param object:
    :param recording_day (str): PXX, postnatal day we want to analyze
    :return:
    '''
    percent_responsive = {k: [] for k in np.unique([a.split ('_')[1] for a in object.list_animals])}
    for animal in object.list_animals:
        for day in [d for d in object.dat[animal].dat_subject.keys() if (calculate_animal_age(object.animal_dobs[animal], d) == recording_day) ]:# or (calculate_animal_age(object.animal_dobs[animal], d) == 'P77') or (calculate_animal_age(object.animal_dobs[animal], d) == 'P84')]:

            group = animal.split('_')[1]
            responsive_cells = object.dat[animal].dat_subject[day]['grat']['thresholded_cells']
            percent = 100*responsive_cells.sum()/len(responsive_cells)
            percent_responsive[group].append(percent)
    for group in percent_responsive:
        percent_responsive[group] = np.array(percent_responsive[group])


    fig, ax = plt.subplots(figsize=(4, 4))
    for i, group in enumerate(percent_responsive.keys()):
        ax.bar(i, percent_responsive[group].mean(), color = 'grey', alpha = 0.4, label=group)
        ax.scatter([i]*len(percent_responsive[group]), percent_responsive[group], color='grey', alpha=0.8, label=group)

        #error = sem(percent_responsive[group], nan_policy='omit')
        ax.bar([i], percent_responsive[group].mean(), yerr=sem(percent_responsive[group], nan_policy='omit'), capsize=5, color='grey', alpha=0.4)

    ax.set_xlabel('Group')
    ax.set_ylabel('% cells')
    ax.set_title(f'Responsive Cells ({recording_day})')
    ax.set_xticks(np.arange(len(percent_responsive)))
    ax.set_xticklabels(list(percent_responsive.keys()))
    plt.tight_layout()

    stat, p = ttest_ind(percent_responsive['ctrl'], percent_responsive['dark'], equal_var=False)
    print(f"t-test: t={stat:.3f}, p={p:.4f}")

    stat, p = mannwhitneyu(percent_responsive['ctrl'], percent_responsive['dark'], alternative='two-sided')
    print(f"Mann-Whitney U test: U={stat:.3f}, p={p:.4f}")

    if p < 0.05:
        x1, x2 = 0, 1  # positions of the two bars
        y, h, col = max(percent_responsive['ctrl'].mean(), percent_responsive['dark'].mean()) + 2, 2, 'k'
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        ax.text((x1 + x2) * 0.5, y + h + 0.5, '*', ha='center', va='bottom', color=col, fontsize=16)

    plt.show()


def responsive_cells_across_days (object):
    '''
    Plot the percentage of responsive cells, across all days (scatter/line plot)

    :param object:
    '''

    dict_responsive_cells = responsive_cells(object)


    fig, ax = plt.subplots(nrows=len(dict_responsive_cells), ncols=1, figsize=(7, 4 * len(dict_responsive_cells)), sharex=True, sharey=True)

    colours = {'dark': 'blue', 'ctrl': 'black', 'light': 'red'}

    if len(dict_responsive_cells) == 1:
        ax = [ax]

    for i_group, group in enumerate(dict_responsive_cells):
        # print(i_group, group)
        session_keys = list(dict_responsive_cells[group].keys())
        x_pos = np.arange(1, len(session_keys) + 1)
        session_vals = list(dict_responsive_cells[group].values())

        max_animals = max(len(vals) for vals in session_vals)
        for i_animal in range(max_animals):
            y_vals = []
            x_vals = []
            for i, vals in enumerate(session_vals):
                if i_animal < len(vals):
                    y_vals.append(vals[i_animal])
                    jitter = np.random.normal(0, 0.05)
                    x_vals.append(x_pos[i] + jitter)
            ax[i_group].plot(x_vals, y_vals, marker='o', color='gray', alpha=0.4)

        # Plot mean points without jitter
        means = [np.mean(vals) for vals in session_vals]
        ax[i_group].plot(x_pos, means, marker='d', markersize=10, color=colours.get(group, 'gray'), label=group)
        ax[i_group].plot(x_pos, means, color=colours.get(group, 'gray'))
        ax[i_group].set_xticks(x_pos)
        ax[i_group].set_xticklabels(session_keys)
        ax[i_group].set_xlabel('Recording session')
        ax[i_group].set_ylabel('% responsive cells')
        ax[i_group].set_title(f'{group} group')
        ax[i_group].grid(axis='y', linestyle='--', alpha=0.5)

    # Make sure x-tick labels are visible on all axes
    for axis in ax:
        axis.tick_params(labelbottom=True)  # force labels on

    session_keys = list(dict_responsive_cells[list(dict_responsive_cells.keys())[0]].keys())
    x_pos = np.arange(1, len(session_keys) + 1)

    for axis in ax:
        axis.set_xticks(x_pos)
        axis.set_xticklabels(session_keys)

    plt.tight_layout()
    plt.show()



def hist_osi_angle (obj, metric):
    '''
    plots histogram distribution of OSIs/preferred orientation for each day ACROSS ALL ANIMALS

    :param obj:
    :param metric: 'OSI' or 'DSI' or 'preferred_orientation' or 'preferred_direction'
    '''

    d = {k: {} for k in np.unique([a.split ('_')[1] for a in obj.list_animals])}

    for animal in obj.list_animals:
        group = animal.split('_')[1]

        for day in obj.dat[animal].dat_subject.keys():

            postnatal_day = calculate_animal_age(obj.animal_dobs[animal], day)
            responsive_indices = obj.dat[animal].dat_subject[day]['grat']['thresholded_cells']
            osi = obj.dat[animal].dat_subject[day]['grat'][metric][responsive_indices]  # shape n_cells

            if postnatal_day not in d[group]:
                d[group][postnatal_day] = [osi]
            else:
                d[group][postnatal_day].append(osi)

    # combine all data from all animals into a single array that we can plot
    for group in d.keys():
        for day in d[group].keys():
            d[group][day] = np.concatenate([arr.ravel() for arr in d[group][day]])#.reshape(-1) # average across animals for each day

    postnatal_days = sorted(set(day for group_data in d.values() for day in group_data))

    for post_d in postnatal_days:
        ctrl = np.asarray(d['ctrl'][post_d])
        dark = np.asarray(d['dark'][post_d])
        stat, p = ks_2samp(ctrl, dark, alternative='two-sided', mode='auto') # Two-sample KS test (two-sided by default)
        print(f"KS test (two-sided) ({post_d}): D={stat:.4f}, p={p:.4e}")

    colours = {'ctrl': 'black', 'dark': 'red'}
    fig, ax = plt.subplots (nrows = len(d.keys()), ncols = len(postnatal_days), figsize = (9,5), sharex=True, sharey=True)

    for i_group, group in enumerate(d.keys()):
        for i_day, day in enumerate(d[group].keys()):

            if metric == 'preferred_orientation':
                bins = np.linspace(-180,180,8)
            elif metric == 'OSI' or metric == 'DSI':
                bins = np.linspace(0, 1, 25)

            weights = np.ones_like(d[group][day]) * (100.0 / d[group][day].size)
            print(weights)
            ax[i_group, i_day].hist(
                d[group][day], bins=bins, weights=weights,
                density=False, color=colours.get(group, 'gray'), alpha=0.5
            )

            #ax[i_group, i_day].hist (d[group][day], density = True, color = colours.get(group, 'gray'), alpha = 0.5)# bins = bins)

            # ax[i_group, i_day].set_ylim([0, 6])
            # ax[i_group, i_day].set_xlim([-100, 100])
            # ax[i_group, i_day].set_ylim([0, 18])
            ax[i_group, i_day].set_label('cell count')
            ax[i_group, i_day].set_title(f'{day} ({group})')

    # font_kwargs = dict(fontsize="large")
    # add_headers(fig, col_headers=["Preferred Orientation", "Orientation Selectivity Index"], row_headers=['Cell count \n (' + day + ')' for day in list(obj.dat_subject.keys())], **font_kwargs)
    plt.suptitle(metric)
    plt.tight_layout()
    plt.show()
hist_osi_angle (data_object, 'OSI')
def reliability_x_osi (object):
    '''
    :param object:
    :param recording_day (str): PXX, postnatal day we want to analyze
    :return:
    '''

    reliability = {k: {} for k in np.unique([a.split ('_')[1] for a in object.list_animals])}
    osis = {k: {} for k in np.unique([a.split('_')[1] for a in object.list_animals])}

    for animal in object.list_animals:
        group = animal.split('_')[1]

        for day in object.dat[animal].dat_subject.keys():
            postnatal_day = calculate_animal_age(object.animal_dobs[animal], day)
            responsive_indices = object.dat[animal].dat_subject[day]['grat']['thresholded_cells']
            matrix = object.dat[animal].dat_subject[day]['grat']['zscored'][:,:, responsive_indices,object.fps*1:object.fps*5] # shape n_repeats, n_ori, n_cells, n_timepoints
            r = trial_reliability(matrix).max(axis = 0) # shape (n_cells)

            osi = object.dat[animal].dat_subject[day]['grat']['OSI'][responsive_indices]

            if postnatal_day not in reliability[group]:
                reliability[group][postnatal_day] = [r]
                osis[group][postnatal_day] = [osi]
            else:
                reliability[group][postnatal_day].append(r)
                osis[group][postnatal_day].append(osi)

    # combine all data from all animals into a single array that we can plot
    for group in reliability.keys():
        for day in reliability[group].keys():
            reliability[group][day] = np.concatenate([arr.ravel() for arr in reliability[group][day]])#.reshape(-1) # average across animals for each day
            osis[group][day] = np.concatenate([arr.ravel() for arr in osis[group][day]])  # .reshape(-1) # average across animals for each day

    postnatal_days = sorted(set(day for group_data in reliability.values() for day in group_data))
    colours = {'ctrl': 'black', 'dark': 'red'}
    fig, ax = plt.subplots (nrows = len(reliability.keys()), ncols = len(postnatal_days), figsize = (10,7), sharex = True, sharey = True)

    for i_group, group in enumerate(reliability.keys()):
        for i_day, day in enumerate(reliability[group].keys()):

            ax[i_group, i_day].scatter(osis[group][day], reliability[group][day], color = colours[group], alpha = 0.5)
            ax[i_group, i_day].set_title(f'{group}, {day}')
            ax[i_group, i_day].set_ylabel('reliability')
            ax[i_group, i_day].set_xlabel('OSI')

    plt.suptitle('trial-by-trial correlation x OSI')
    plt.tight_layout()
    plt.show()

#reliability_x_osi(data_object)

#trial_by_trial_reliability(data_object)
#trial_by_trial_reliability_days(data_object)
# hist_osi_angle (data_object, 'OSI')
# hist_osi_angle (data_object, 'DSI')
# hist_osi_angle (data_object, 'preferred_orientation')
# hist_osi_angle (data_object, 'preferred_direction')


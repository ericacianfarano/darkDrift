import matplotlib.pyplot as plt

from helpers import *
from imports import *

def plot_avg_response (object, group_type):
    '''
    :param object:
    :param recording_day (str): PXX, postnatal day we want to analyze
    :param group_type (str): EB or LB
    :return:
    '''

    response = {k: {} for k in np.unique([a.split ('_')[1] for a in object.list_animals if group_type in a ])}
    sub_list_animals = [animal for animal in object.list_animals if group_type in animal]

    for animal in sub_list_animals:
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
    colours = {'EBc': 'black', 'EB': 'blue', 'LBc': 'black', 'LB': 'red'}

    fig, ax = plt.subplots (len(postnatal_days),1, figsize = (6,len(postnatal_days)*2))
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
        #if set(all_group_data.keys()) == {'EBc', 'EB'}:
        #print([k for k in all_group_data.keys() if 'c' in k],[k for k in all_group_data.keys() if 'c' not in k])
        ctrl = all_group_data[[k for k in all_group_data.keys() if 'c' in k][0]]
        dark = all_group_data[[k for k in all_group_data.keys() if 'c' not in k][0]]

        if ctrl.shape[1] == dark.shape[1]:  # ensure same n_timepoints
            pvals = np.array([ttest_ind(ctrl[:, t], dark[:, t], equal_var=False).pvalue for t in range(ctrl.shape[1])])
            sig_mask = pvals < 0.05

            # plot asterisks
            y_max = max(ctrl.max(), dark.max())
            for t in np.where(sig_mask)[0]:
                ax[i_day].text(t, y_max + 0.1, '*', ha='center', va='bottom', fontsize=10, color='k')

        #ax[i_day].set_ylim ([min(ctrl.min(), dark.min()), max(ctrl.max(), dark.max()) + 0.25 * max(ctrl.max(), dark.max())])
        ax[i_day].set_ylabel('Response (z-scored)')
        ax[i_day].set_xticks(np.arange(ctrl.shape[1])[::object.fps])
        ax[i_day].set_xticklabels( [int(a/object.fps) for a in (np.arange(ctrl.shape[1])[::object.fps])])
        ax[i_day].axvline(x=object.fps, color='k', linestyle=':', linewidth=1)
        ax[i_day].axvspan(object.fps * 2, object.fps * 5, color='gray', alpha=0.2)
        ax[i_day].set_title(f"Postnatal Day {postnatal_day}")
        ax[i_day].legend()

    ax[-1].set_xlabel("Time (s)")
    plt.tight_layout()

    folder_path = os.path.join(os.path.dirname(object.path), 'figures', 'population_all_animals', 'avg_response')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'avg_response {group_type}.svg'))
    plt.savefig(os.path.join(folder_path, f'avg_response {group_type}.png'))

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
    colors = plt.cm.plasma(np.linspace(0, 0.8, 4))

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

    folder_path = os.path.join(os.path.dirname(object.path), 'figures', 'population_all_animals', 'avg_response')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'avg_response_all.svg'))
    plt.savefig(os.path.join(folder_path, f'avg_response_all.png'))


    plt.show()


def trial_by_trial_reliability (object, recording_day):
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
        ax[i].set_title(group + ' '+recording_day)
        ax[i].set_xlabel('correlation')
    plt.suptitle('trial-by-trial correlation')

    folder_path = os.path.join(os.path.dirname(object.path), 'figures', 'population_all_animals', 'reliability')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'reliability {recording_day}.svg'))
    plt.savefig(os.path.join(folder_path, f'reliability {recording_day}.png'))

    plt.show()

def trial_by_trial_reliability_days (object):
    '''
    :param object:
    :param recording_day (str): PXX, postnatal day we want to analyze
    :return:
    '''

    d = {k: {} for k in np.unique([a.split ('_')[1] for a in object.list_animals])}
    colours = {'EBc': 'black', 'EB': 'blue', 'LBc': 'green', 'LB': 'red'}

    for animal in object.list_animals:
        group = animal.split('_')[1]

        for day in object.dat[animal].dat_subject.keys():
            postnatal_day = calculate_animal_age(object.animal_dobs[animal], day)
            responsive_indices = object.dat[animal].dat_subject[day]['grat']['thresholded_cells']
            matrix = object.dat[animal].dat_subject[day]['grat']['zscored'][:,:, responsive_indices,object.fps*1:object.fps*5] # shape n_repeats, n_ori, n_cells, n_timepoints

            #trial_reliability(matrix) is of shape n_ori, n_cells > so take the 'best' reliability across orientations
            r = trial_reliability(matrix).max(axis = 0) # shape (n_cells)
            # r = trial_reliability_3(matrix)  # shape (n_cells)

            if postnatal_day not in d[group]:
                d[group][postnatal_day] = [r]
            else:
                d[group][postnatal_day].append(r)

    # combine all data from all animals into a single array that we can plot
    for group in d.keys():
        for day in d[group].keys():
            d[group][day] = np.concatenate([arr.ravel() for arr in d[group][day]])#.reshape(-1) # average across animals for each day

    postnatal_days = sorted(set(day for group_data in d.values() for day in group_data))
    fig, ax = plt.subplots (nrows = len(d.keys()), ncols = 4, figsize = (10,7), sharex = True, sharey = True)

    for i_group, group in enumerate(d.keys()):
        for i_day, day in enumerate(d[group].keys()):

            ax[i_group, i_day].hist(d[group][day], bins=np.linspace(0,  1, 15), color = colours[group], alpha = 0.5, density=True)
            ax[i_group, i_day].set_title(f'{group}, {day}')
            ax[i_group, i_day].set_xlabel('correlation')

    plt.suptitle('trial-by-trial correlation')
    plt.tight_layout()

    folder_path = os.path.join(os.path.dirname(object.path), 'figures', 'population_all_animals', 'reliability')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'reliability days.svg'))
    plt.savefig(os.path.join(folder_path, f'reliability days.png'))

    plt.show()


def responsiveness (object, recording_day):
    '''
    For 'recording_day' plot the percentage of responsive cells for each group, as a bar plot

    :param object:
    :param recording_day (str): PXX, postnatal day we want to analyze
    :return:
    '''
    if int(recording_day[1:]) <= 91:
        group_string = 'EB'
    elif int(recording_day[1:]) >= 105:
        group_string = 'LB'

    percent_responsive = {
        a.split('_')[1]: []
        for a in object.list_animals
        if group_string in a}

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

    # stat, p = ttest_ind(percent_responsive['EBc'], percent_responsive['EB'], equal_var=False)
    # print(f"t-test: t={stat:.3f}, p={p:.4f}")
    #
    # stat, p = mannwhitneyu(percent_responsive['EBc'], percent_responsive['EB'], alternative='two-sided')
    # print(f"Mann-Whitney U test: U={stat:.3f}, p={p:.4f}")
    #
    # if p < 0.05:
    #     x1, x2 = 0, 1  # positions of the two bars
    #     y, h, col = max(percent_responsive['EBc'].mean(), percent_responsive['EB'].mean()) + 2, 2, 'k'
    #     ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
    #     ax.text((x1 + x2) * 0.5, y + h + 0.5, '*', ha='center', va='bottom', color=col, fontsize=16)
    folder_path = os.path.join(os.path.dirname(object.path), 'figures', 'population_all_animals', 'responsiveness')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'responsiveness {recording_day}.svg'))
    plt.savefig(os.path.join(folder_path, f'responsiveness {recording_day}.png'))

    plt.show()

def responsive_cells_across_days(object):

    dict_responsive_cells = responsive_cells(object)

    ncols = 2
    fig, ax = plt.subplots(
        nrows=int(np.ceil(len(dict_responsive_cells)/ncols)),
        ncols=ncols, figsize=(8, 5), sharey=True
    )
    ax = ax.ravel()
    colours = {'EB': 'blue', 'EBc': 'black', 'LB': 'red', 'LBc':'green'}

    max_days = max(
        len(days)
        for group_dict in dict_responsive_cells.values()
        for days in group_dict.values())

    if len(dict_responsive_cells) == 1:
        ax = [ax]

    for i_group, group in enumerate(dict_responsive_cells):
        group_days = sorted(
            {d for days in dict_responsive_cells[group].values() for d in days.keys()},
            key=sort_key
        )

        ax[i_group].plot(group_days, [np.nan]*len(group_days), alpha=0)

        for animal in dict_responsive_cells[group]:
            sorted_items = sorted(dict_responsive_cells[group][animal].items(),
                                  key=lambda x: sort_key(x[0]))
            x = [k for k, _ in sorted_items]
            y = [v for _, v in sorted_items]
            ax[i_group].plot(x, y, marker='o', color='gray', alpha=0.4)

        # group mean ± SEM
        values_by_day = {}
        for animal, days in dict_responsive_cells[group].items():
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
        ax[i_group].set_ylabel('Responsive Cells (%)')
        ax[i_group].set_title(f'{group} group')
        ax[i_group].grid(axis='y', linestyle='--', alpha=0.5)
        ax[i_group].tick_params(axis='x', rotation=45)
        ax[i_group].set_xlim([-0.5, max_days])

    plt.tight_layout()

    folder_path = os.path.join(os.path.dirname(object.path), 'figures', 'population_all_animals', 'responsiveness')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'responsiveness all_animals.svg'))
    plt.savefig(os.path.join(folder_path, f'responsiveness all_animals.png'))

    plt.show()


def hist_osi_angle_all_animals (obj, metric):
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
    postnatal_days = {}
    for group in d.keys():
        for day in d[group].keys():
            d[group][day] = np.concatenate([arr.ravel() for arr in d[group][day]])#.reshape(-1) # average across animals for each day
        postnatal_days[group] = sorted(set(day for day in d[group]))

    # for post_d in postnatal_days:
    #     ctrl = np.asarray(d['EBc'][post_d])
    #     dark = np.asarray(d['EB'][post_d])
    #     stat, p = ks_2samp(ctrl, dark, alternative='two-sided', mode='auto') # Two-sample KS test (two-sided by default)
    #     print(f"KS test (two-sided) ({post_d}): D={stat:.4f}, p={p:.4e}")

    colours = {'EBc': 'black', 'EB': 'blue', 'LBc': 'green', 'LB': 'red'}
    fig, ax = plt.subplots (nrows = len(d.keys()), ncols = max([len(l) for l in postnatal_days.values()]), figsize = (10,8), sharex=True, sharey=True)

    for i_group, group in enumerate(d.keys()):
        for i_day, day in enumerate(d[group].keys()):

            if (metric == 'preferred_orientation') or (metric == 'preferred_direction'):
                bins = np.linspace(-180,180,8)
            elif metric == 'OSI' or metric == 'DSI':
                bins = np.linspace(0, 1, 20)

            weights = np.ones_like(d[group][day]) * (100.0 / d[group][day].size)
            ax[i_group, i_day].hist(
                d[group][day], bins=bins, weights=weights,
                density=False, color=colours.get(group, 'gray'), alpha=0.5
            )
            #ax[i_group, i_day].hist (d[group][day], density = True, color = colours.get(group, 'gray'), alpha = 0.5)# bins = bins)
            # ax[i_group, i_day].set_ylim([0, 6])
            # ax[i_group, i_day].set_xlim([-100, 100])
            # ax[i_group, i_day].set_ylim([0, 18])
            # ax[i_group, i_day].set_ylabel(f'{group} \n cell count')
            if i_day == 0:
                ax[i_group, i_day].set_ylabel(
                    rf'$\bf{{{group}}}$' + '\n % cells',
                    fontsize=14)
            ax[i_group, i_day].set_xlabel(metric)
            ax[i_group, i_day].set_title(f'{day}')

    # font_kwargs = dict(fontsize="large")
    # add_headers(fig, col_headers=["Preferred Orientation", "Orientation Selectivity Index"], row_headers=['Cell count \n (' + day + ')' for day in list(obj.dat_subject.keys())], **font_kwargs)
    plt.suptitle(metric)
    plt.tight_layout()
    folder_path = os.path.join(os.path.dirname(obj.path), 'figures', 'population_all_animals','histograms')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'hist_osi_across_days_{metric}.svg'))
    plt.savefig(os.path.join(folder_path, f'hist_osi_across_days_{metric}.png'))

    plt.show()

def hist_osi_angle_all_animals_outlines (obj, metric):
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
    postnatal_days = {}
    for group in d.keys():
        for day in d[group].keys():
            d[group][day] = np.concatenate([arr.ravel() for arr in d[group][day]])#.reshape(-1) # average across animals for each day
        postnatal_days[group] = sorted(set(day for day in d[group]))

    fig, ax = plt.subplots (nrows = 1, ncols =len(d.keys()) , figsize = (10,4), sharex=True, sharey=True)
    cmap = get_cmap('plasma')

    for i_group, group in enumerate(d.keys()):

        days = sorted(d[group].keys())
        colors = cmap(np.linspace(0, 0.9, len(days)))

        if (metric == 'preferred_orientation') or (metric == 'preferred_direction'):
            bins = np.linspace(-180, 180, 8)
        elif metric == 'OSI' or metric == 'DSI':
            bins = np.linspace(0, 1, 20)

        for c, day in zip(colors, days):

            # this transforms the y axis into 'percentage of cells'
            weights = np.ones_like(d[group][day]) * (100.0 / d[group][day].size)
            ax[i_group].hist(
                d[group][day], bins=bins, weights=weights, histtype = 'step', linewidth = 2, label = f'{day}, {len(d[group][day])} cells',
                density=False, color=c, alpha=0.5)

            ax[i_group].axvline(
                np.mean(d[group][day]), color=c, linestyle='--', linewidth=1.5, alpha = 0.7)

            ax[i_group].set_xlabel(metric)
            ax[i_group].set_ylabel('% cells')
            ax[i_group].set_title(f'{group}')

        ax[i_group].legend(frameon=False, fontsize=8)

    # font_kwargs = dict(fontsize="large")
    # add_headers(fig, col_headers=["Preferred Orientation", "Orientation Selectivity Index"], row_headers=['Cell count \n (' + day + ')' for day in list(obj.dat_subject.keys())], **font_kwargs)
    plt.suptitle(metric)
    plt.tight_layout()
    folder_path = os.path.join(os.path.dirname(obj.path), 'figures', 'population_all_animals','histograms')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'hist_osi_across_days_outlines_{metric}.svg'))
    plt.savefig(os.path.join(folder_path, f'hist_osi_across_days_outlines_{metric}.png'))
    plt.show()

    return d, postnatal_days

def plot_mean_hist_metrics (obj, metric):
    d, postnatal_days = hist_osi_angle_all_animals_outlines(obj, metric)
    # plt.close('all')

    fig, ax = plt.subplots (nrows = 1, ncols =2 , figsize = (10,4), sharex=True, sharey=True)
    colours = {'EBc': 'black', 'EB': 'blue', 'LBc': 'green', 'LB': 'red'}
    ticks = ['W0', 'W1', 'W2', 'W3']

    for i_group, group in enumerate(d.keys()):

        if 'EB' in group:
            ax_num = 0
        else:
            ax_num = 1

        days = sorted(d[group].keys())

        means = np.array([d[group][day].mean() for day in d[group].keys()])
        sems = np.array([sem(d[group][day]) for day in d[group].keys()])

        ax[ax_num].plot(np.arange(len(means)), means, c = colours[group], label = group)

        ax[ax_num].fill_between(
            np.arange(len(means)),
            means - sems,
            means + sems,
            color=colours[group],
            alpha=0.3
        )

    plt.suptitle(metric)
    ax[0].legend()
    ax[1].legend()
    ax[0].set_xticks(np.arange(len(ticks)), ticks)
    ax[1].set_xticks(np.arange(len(ticks)), ticks)
    plt.show()



def hist_osi_angle_all_animals_outlines_2 (obj, metric):
    '''
    plots histogram distribution of OSIs/preferred orientation for each day ACROSS ALL ANIMALS

    :param obj:
    :param metric: 'OSI' or 'DSI' or 'preferred_orientation' or 'preferred_direction'
    '''

    d = {k: {} for k in np.unique([a.split ('_')[1] for a in obj.list_animals])}

    for animal in obj.list_animals:
        group = animal.split('_')[1]

        for day in list(obj.dat[animal].dat_subject.keys())[::3]:

            postnatal_day = calculate_animal_age(obj.animal_dobs[animal], day)
            responsive_indices = obj.dat[animal].dat_subject[day]['grat']['thresholded_cells']
            osi = obj.dat[animal].dat_subject[day]['grat'][metric][responsive_indices]  # shape n_cells

            if postnatal_day not in d[group]:
                d[group][postnatal_day] = [osi]
            else:
                d[group][postnatal_day].append(osi)

    # combine all data from all animals into a single array that we can plot
    postnatal_days = {}
    for group in d.keys():
        for day in d[group].keys():
            d[group][day] = np.concatenate([arr.ravel() for arr in d[group][day]])#.reshape(-1) # average across animals for each day
        postnatal_days[group] = sorted(set(day for day in d[group]))

    fig, ax = plt.subplots (nrows = 1, ncols =len(d.keys()) , figsize = (13,4), sharex=True, sharey=True)
    cmap = get_cmap('plasma')

    for i_group, group in enumerate(d.keys()):

        days = sorted(d[group].keys())
        colors = cmap(np.linspace(0, 0.7, len(days)))

        if (metric == 'preferred_orientation') or (metric == 'preferred_direction'):
            bins = np.linspace(-180, 180, 8)
        elif metric == 'OSI' or metric == 'DSI':
            bins = np.linspace(0, 1, 20)

        for c, day in zip(colors, days):

            # this transforms the y axis into 'percentage of cells'
            weights = np.ones_like(d[group][day]) * (100.0 / d[group][day].size)
            ax[i_group].hist(
                d[group][day], bins=bins, weights=weights, histtype = 'step', linewidth = 2, label = f'{day}, {len(d[group][day])} cells',
                density=False, color=c, alpha=0.5)

            ax[i_group].axvline(
                np.mean(d[group][day]), color=c, linestyle='--', linewidth=1.5, alpha = 0.7)

            ax[i_group].set_xlabel(metric)
            ax[i_group].set_ylabel('% cells')
            ax[i_group].set_title(f'{group}')

        ax[i_group].legend(frameon=False, fontsize=8)

    # font_kwargs = dict(fontsize="large")
    # add_headers(fig, col_headers=["Preferred Orientation", "Orientation Selectivity Index"], row_headers=['Cell count \n (' + day + ')' for day in list(obj.dat_subject.keys())], **font_kwargs)
    plt.suptitle(metric)
    plt.tight_layout()
    # folder_path = os.path.join(os.path.dirname(obj.path), 'figures', 'population_all_animals','histograms')
    # if not os.path.exists(folder_path):
    #     os.makedirs(folder_path)
    # plt.savefig(os.path.join(folder_path, f'hist_osi_across_days_outlines_{metric}.svg'))
    # plt.savefig(os.path.join(folder_path, f'hist_osi_across_days_outlines_{metric}.png'))
    plt.show()


def hist_osi_angle_all_animals_outlines_3 (obj, metric):
    '''
    plots histogram distribution of OSIs/preferred orientation for each day ACROSS ALL ANIMALS

    :param obj:
    :param metric: 'OSI' or 'DSI' or 'preferred_orientation' or 'preferred_direction'
    '''

    d = {k: {} for k in np.unique([a.split ('_')[1] for a in obj.list_animals])}

    for animal in obj.list_animals:
        group = animal.split('_')[1]

        for day in list(obj.dat[animal].dat_subject.keys())[0:1]:

            postnatal_day = calculate_animal_age(obj.animal_dobs[animal], day)
            responsive_indices = obj.dat[animal].dat_subject[day]['grat']['thresholded_cells']
            osi = obj.dat[animal].dat_subject[day]['grat'][metric][responsive_indices]  # shape n_cells

            if postnatal_day not in d[group]:
                d[group][postnatal_day] = [osi]
            else:
                d[group][postnatal_day].append(osi)

    # combine all data from all animals into a single array that we can plot
    postnatal_days = {}
    for group in d.keys():
        for day in d[group].keys():
            d[group][day] = np.concatenate([arr.ravel() for arr in d[group][day]])#.reshape(-1) # average across animals for each day
        postnatal_days[group] = sorted(set(day for day in d[group]))

    print(f"\n--- STATISTICAL TESTS ---, {metric}")

    ks_stat, ks_p = ks_2samp(d['EB']['P70'], d['EBc']['P70'])
    t_stat, t_p = ttest_ind(d['EB']['P70'], d['EBc']['P70'], equal_var=False)
    kw_stat, kw_p = kruskal(d['EB']['P70'], d['EBc']['P70'])

    print('EB vs EBc')
    print(f"    KS test: stat={ks_stat:.3f}, p={ks_p:.3e}")
    print(f"    t-test:  stat={t_stat:.3f}, p={t_p:.3e}")
    print(f"    KW test: stat={kw_stat:.3f}, p={kw_p:.3e}")

    ks_stat, ks_p = ks_2samp(d['LB']['P105'], d['LBc']['P105'])
    t_stat, t_p = ttest_ind(d['LB']['P105'], d['LBc']['P105'], equal_var=False)
    kw_stat, kw_p = kruskal(d['LB']['P105'], d['LBc']['P105'])

    print('LB vs LBc')
    print(f"    KS test: stat={ks_stat:.3f}, p={ks_p:.3e}")
    print(f"    t-test:  stat={t_stat:.3f}, p={t_p:.3e}")
    print(f"    KW test: stat={kw_stat:.3f}, p={kw_p:.3e}")

    fig, ax = plt.subplots (nrows = 1, ncols =2 , figsize = (9,4), sharex=True, sharey=True)
    cmap = get_cmap('plasma')

    for i_group, group in enumerate(d.keys()):

        days = sorted(d[group].keys())
        colors = cmap(np.linspace(0, 0.7, 2))

        if (metric == 'preferred_orientation') or (metric == 'preferred_direction'):
            bins = np.linspace(-180, 180, 8)
        elif metric == 'OSI' or metric == 'DSI':
            bins = np.linspace(0, 1, 20)

        for c, day in zip(colors, days):


            if 'EB' in group:
                i_plot = 0
                lab = 'EB-EBc'

            elif 'LB' in group:
                lab = 'LB-LBc'
                i_plot = 1

            col = 'red'
            if 'c' in group:
                col = 'black'
            # this transforms the y axis into 'percentage of cells'
            weights = np.ones_like(d[group][day]) * (100.0 / d[group][day].size)
            ax[i_plot].hist(
                d[group][day], bins=bins, weights=weights, histtype = 'step', linewidth = 2, label = f'{group}, {len(d[group][day])} cells',
                density=False, color=col, alpha=0.5)

            ax[i_plot].axvline(
                np.mean(d[group][day]), color=col, linestyle='--', linewidth=1.5, alpha = 0.7)

            ax[i_plot].set_xlabel(metric)
            ax[i_plot].set_ylabel('% cells')
            ax[i_plot].set_title(f'{day}')

        ax[i_plot].legend(frameon=False, fontsize=8)

    # font_kwargs = dict(fontsize="large")
    # add_headers(fig, col_headers=["Preferred Orientation", "Orientation Selectivity Index"], row_headers=['Cell count \n (' + day + ')' for day in list(obj.dat_subject.keys())], **font_kwargs)
    plt.suptitle(metric)
    plt.tight_layout()
    # folder_path = os.path.join(os.path.dirname(obj.path), 'figures', 'population_all_animals','histograms')
    # if not os.path.exists(folder_path):
    #     os.makedirs(folder_path)
    # plt.savefig(os.path.join(folder_path, f'hist_osi_across_days_outlines_{metric}.svg'))
    # plt.savefig(os.path.join(folder_path, f'hist_osi_across_days_outlines_{metric}.png'))
    plt.show()



# hist_osi_angle (data_object, 'OSI')
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
    colours = {'EBc': 'black', 'EB': 'blue', 'LBc': 'black', 'LB': 'red'}
    fig, ax = plt.subplots (nrows = len(reliability.keys()), ncols = 4, figsize = (10,7), sharex = True, sharey = True)

    for i_group, group in enumerate(reliability.keys()):
        for i_day, day in enumerate(reliability[group].keys()):

            ax[i_group, i_day].scatter(osis[group][day], reliability[group][day], color = colours[group], alpha = 0.5)
            ax[i_group, i_day].set_title(f'{group}, {day}')
            ax[i_group, i_day].set_ylabel('reliability')
            ax[i_group, i_day].set_xlabel('OSI')

    plt.suptitle('trial-by-trial correlation x OSI')
    plt.tight_layout()
    folder_path = os.path.join(os.path.dirname(object.path), 'figures', 'population_all_animals','reliability')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'reliability x osi.svg'))
    plt.savefig(os.path.join(folder_path, f'reliability x osi.png'))
    plt.show()

#reliability_x_osi(data_object)

#trial_by_trial_reliability(data_object)
#trial_by_trial_reliability_days(data_object)
# hist_osi_angle (data_object, 'OSI')
# hist_osi_angle (data_object, 'DSI')
# hist_osi_angle (data_object, 'preferred_orientation')
# hist_osi_angle (data_object, 'preferred_direction')

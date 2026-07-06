from retmap import *

path = r'E:\dark_drift\data'
for animal in ['EC_EBc_05']:
    for day in ['20251202']:
        for subfile in [item for item in os.listdir(os.path.join(path, animal, day)) if
                        os.path.isdir(os.path.join(path, animal, day, item)) and ('retmap' in item)]:
            # config_path = fr'E:\retmap\{animal}\{day}\{subfile}\config.txt'
            path_to_wdf = fr'E:\dark_drift\data\{animal}\{day}\{subfile}\wdf'

            mode = "1"
            load_maps
            run_retmap(path_to_wdf, mode, min_area=250, animal_dob='20250902',
                       animal=animal, day=day, draw_mask=False)

from retmap import *

if __name__ == '__main__':

    # run in console
    from retmap import *

    animals_days = {
                    'EC_dark_01': ['20250609', '20250616'],
                    'EC_dark_03': ['20250609', '20250616'],
                    'EC_dark_05': ['20250610', '20250617'],
                    'EC_dark_08': ['20250610', '20250617']
                    }

    animals_birthdate = {'EC_dark_01': '20250331',
                         'EC_dark_03': '20250331',
                         'EC_dark_05': '20250401',
                         'EC_dark_08': '20250401',
                   }


    path = r'I:\dark_drift\data'
    for animal in animals_days.keys():
        for day in animals_days[animal]:
            for subfile in [item for item in os.listdir(os.path.join(path, animal, day)) if
                            os.path.isdir(os.path.join(path, animal, day, item)) and ('retmap' in item)]:
                print(animal, day, subfile)
                config_path = os.path.join (path, animal, day, 'retmap','wdf','config.txt')# fr'I:\retmap\{animal}\{day}\{subfile}\config.txt'
                mode = "2"
                run_retmap(config_path, mode, min_area = 250, animal = animal, day = day, draw_mask = False)

    plt.close('all')
    # total visual area across different groups
    total_visual_area_g, largest_neg_patch_area_g, mean_amplitude_g = {}, {}, {}
    for animal in animals_days:
        g = animal[3:-3]

        if (g not in total_visual_area_g):
            total_visual_area_g[g] = [total_visual_area[animal]]
            largest_neg_patch_area_g[g] = [largest_neg_patch_area[animal]]
            mean_amplitude_g[g] = [mean_amplitude[animal]]
        else:
            total_visual_area_g[g].append(total_visual_area[animal])
            largest_neg_patch_area_g[g].append(largest_neg_patch_area[animal])
            mean_amplitude_g[g].append(mean_amplitude[animal])

    fig, ax = plt.subplots (1,3, figsize = (12,4))
    group_names = ['GCaMP6s', 'GNAT', 'RD1', 'RD1opto']
    for i_subplot, metric, in enumerate([total_visual_area_g, largest_neg_patch_area_g, mean_amplitude_g]):

        for i, group in enumerate(group_names):
            arr = np.array(metric[group])

            if group == 'GCaMP6s':
                colour = 'black'
            if group == 'GNAT':
                colour = 'tomato'
            if group == 'RD1':
                colour = 'firebrick'
            if group == 'RD1opto':
                colour = 'blue'

            ax[i_subplot].scatter ([i]*len(arr) + np.random.uniform(-0.2, 0.2, len(arr)), arr, c = colour, alpha = 0.5, s = 30)
            ax[i_subplot].bar ([i], arr.mean(), color = colour, alpha = 0.6)

            if i_subplot == 0:
                ax[i_subplot].set_title ('total visual area')
            elif i_subplot == 1:
                ax[i_subplot].set_title ('largest negative patch')
            elif i_subplot == 2:
                ax[i_subplot].set_title ('mean pixel amplitude')

        # significance tests
        control, gnat, rd1, opto = metric['GCaMP6s'], metric['GNAT'], metric['RD1'], metric['RD1opto']
        stat, p = kruskal(control, gnat, rd1, opto)
        print(f'KW H-statistic: {stat:.3f}, p-value: {p:.3f}')

        if p < 0.05:  # follow up with testing pairwise comparisons (with correction for multiple comparisons)
            # Do pairwise comparisons manually:
            print('mannwhitney two-sided test, with bonferroni multiple comparison correction')
            pvals = [
                mannwhitneyu(control, gnat, alternative='two-sided').pvalue,
                mannwhitneyu(control, rd1, alternative='two-sided').pvalue,
                mannwhitneyu(control, opto, alternative='two-sided').pvalue,
                mannwhitneyu(rd1, opto, alternative='two-sided').pvalue
            ]

            # Apply Bonferroni correction manually (3 comparisons):
            _, pvals_corrected, _, _ = multipletests(pvals, alpha=0.05, method='bonferroni')

            print(pvals_corrected)

            # Define comparisons in same order as pvals
            comparisons = [
                ('GCaMP6s', 'GNAT'),
                ('GCaMP6s', 'RD1'),
                ('GCaMP6s', 'RD1opto'),
                ('RD1', 'RD1opto')
            ]

            # Coordinates for group positions on x-axis
            group_coords = {'GCaMP6s': 0, 'GNAT':1, 'RD1': 2, 'RD1opto': 3}

            max_val = max(np.max(control), np.max(gnat), np.max(rd1), np.max(opto)) * 1.05
            print(max_val)
            step = max_val * 0.05  # vertical spacing between significance bars
            h = max_val

            # Loop through each comparison and plot if significant
            for (group1, group2), p_val in zip(comparisons, pvals_corrected):
                if p_val < 0.05:
                    x1, x2 = group_coords[group1], group_coords[group2]
                    y = h
                    h += step  # update height for next line if needed

                    print(x1, x2, y)

                    # Decide number of stars
                    if p_val < 0.001:
                        stars = '***'
                    elif p_val < 0.01:
                        stars = '**'
                    else:
                        stars = '*'

                    # Draw the line and stars
                    ax[i_subplot].plot([x1, x1, x2, x2], [y, y + step / 2, y + step / 2, y], lw=1.5, c='k')
                    ax[i_subplot].text((x1 + x2) * 0.5, y + step * 0.6, stars,
                                       ha='center', va='bottom', fontsize=14)

        ax[i_subplot].set_xticks(range(len(np.unique([animal[3:-3] for animal in animals_days.keys()]))))
        ax[i_subplot].set_xticklabels(np.unique([animal[3:-3] for animal in animals_days.keys()]))
        ax[i_subplot].set_ylabel('# pixels')
    plt.show()


    # over time
    from retmap_modified import *
    screen = 'big100'
    animals_days = {
                    'EC_dark_01': ['20250609', '20250616'],
                    'EC_dark_03': ['20250609', '20250616'],
                    'EC_dark_05': ['20250610', '20250617'],
                    'EC_dark_08': ['20250610', '20250617']
                    }

    animals_birthdate = {'EC_dark_01': '20250331',
                         'EC_dark_03': '20250331',
                         'EC_dark_05': '20250401',
                         'EC_dark_08': '20250401',
                   }

    path = r'E:\dark_drift\data'
    for animal in ['EC_EBc_05']:
        for day in ['20251202']:
            for subfile in [item for item in os.listdir(os.path.join(path, animal, day)) if
                            os.path.isdir(os.path.join(path, animal, day, item)) and ('retmap' in item)]:
                config_path = fr'E:\retmap\{animal}\{day}\{subfile}\config.txt'
                actual_path = fr'E:\dark_drift\data\{animal}\{day}\{subfile}\wdf'

                mode = "2"
                run_retmap(config_path, actual_path, mode, min_area=250, animal_dob = '20250902', animal=animal, day=day, draw_mask=False)
    plt.close('all')

    # Define colors for each group
    group_colors = {
        'GCaMP6s': 'black',
        'GNAT': 'tomato',
        'RD1': 'firebrick'
    }

    fig, ax = plt.subplots (1,2, figsize = (8,4))
    group_names = ['GCaMP6s', 'GNAT', 'RD1']
    for i_subplot, metric, in enumerate([total_visual_area_days, largest_neg_patch_area_days]):

        for subject, data in metric.items():
            # Assume subject is in the format "EC_GROUP_xx" (e.g., "EC_GCaMP6s_06")
            parts = subject.split('_')
            group = parts[1]  # e.g., "GCaMP6s", "GNAT", "RD1", etc.

            # Extract age and value pairs; convert age from "P86" to 86.
            ages = []
            values = []
            for age_str, val in data.items():
                ages.append(int(age_str.lstrip('P')))
                values.append(val)

            # Sort the pairs by age
            ages = np.array(ages)
            values = np.array(values)
            order = np.argsort(ages)
            ages_sorted = ages[order]
            values_sorted = values[order]

            # Plot the trajectory for this animal
            ax[i_subplot].scatter(ages_sorted, values_sorted, linestyle='-',
                     color=group_colors.get(group, 'gray'), label=subject.split('_')[1], s = 20, alpha = 0.6)
            ax[i_subplot].plot(ages_sorted, values_sorted, linestyle='-',
                     color=group_colors.get(group, 'gray'), alpha = 0.7)

        ax[i_subplot].set_xlabel("Age (Postnatal days)")
        ax[i_subplot].set_ylabel("# pixels")
        if i_subplot==0:
            ax[i_subplot].set_title("Visually responsive area")
        elif i_subplot ==1:
            ax[i_subplot].set_title("Largest negative spot")
        #ax[i_subplot].legend(fontsize=8, loc='best', ncol=2)
        legend_handles = [mpatches.Patch(color=group_colors[group], label=group) for group in group_names]
        ax[i_subplot].legend(handles=legend_handles, fontsize=8, loc='best')
    plt.show()


    # for animal in [item for item in os.listdir(os.path.join(path))]:
    #     for day in [item for item in os.listdir(os.path.join(path, animal)) if
    #                 os.path.isdir(os.path.join(path, animal, item))]:
    #         for subfile in [item for item in os.listdir(os.path.join(path, animal, day)) if
    #                         os.path.isdir(os.path.join(path, animal, day, item))]:
    #             print(animal, day, subfile)
    #             config_path = fr'I:\retmap\{animal}\{day}\{subfile}\config.txt'
    #             mode = "2"
    #             run_retmap(config_path, mode)

    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="File path to config file")
    parser.add_argument("-m", "--mode", help="Choose analysis mode: 1 - Create complex fields\t 2 - Create and plot sign map")
    args = parser.parse_args()

    run_retmap(args.config, args.mode)


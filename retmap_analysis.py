from retmap import *
from helpers import *
from animal_info import *

animals = ['EC_EB_01', 'EC_EB_03', 'EC_EB_05', 'EC_EB_08', 'EC_EB_10', 'EC_EB_13', 'EC_EB_14', 'EC_EB_21', 'EC_EB_22', 'EC_EB_23', 'EC_EB_24', 'EC_EB_27', 'EC_EB_28', 'EC_EB_29', 'EC_EB_31', 'EC_EB_32',# 'EC_EB_33',
           'EC_EBc_02', 'EC_EBc_03', 'EC_EBc_04', 'EC_EBc_05', 'EC_EBc_06', 'EC_EBc_07', 'EC_EBc_09',
           'EC_LB_01', 'EC_LB_02', 'EC_LB_03', 'EC_LB_05', 'EC_LB_06','EC_LB_09', 'EC_LB_10', 'EC_LB_11',
           'EC_LBc_02', 'EC_LBc_03', 'EC_LBc_06', 'EC_LBc_09', 'EC_LBc_11', 'EC_LBc_13', 'EC_LBc_14', 'EC_LBc_16', 'EC_LBc_19', 'EC_LBc_20']

# animals = ['EC_EB_24', 'EC_EB_29',
#           'EC_EBc_04',
#            'EC_LB_11',
#            'EC_LBc_14', 'EC_LBc_16']

# # best example animals
# animals = ['EC_EB_31', 'EC_EBc_02', 'EC_LB_02', 'EC_LBc_06']

class batchRetmap:
    def __init__(self, list_animals, animal_dobs, path):
        self.list_animals = list_animals
        self.animal_dobs = animal_dobs
        self.path = path
        self.dat = {animal: retmapPreprocessing(animal, self.path) for animal in self.list_animals}

class retmapPreprocessing:

    def __init__ (self, animal, path):

        self.animal = animal
        self.save_path = str(Path(path).parents[0] / 'figures')

        animal_path = Path(path) / self.animal

        if not animal_path.exists():
            if animal_path.drive == 'I:':
                animal_path = Path(str(animal_path).replace('I:', 'E:', 1))
            elif animal_path.drive == 'E:':
                animal_path = Path(str(animal_path).replace('E:', 'G:', 1))
            elif animal_path.drive == 'G:':
                animal_path = Path(str(animal_path).replace('G:', 'E:', 1))

        # add date folders with 'retmap' subfile
        self.days = sorted(
            p.name for p in animal_path.iterdir()
            if p.is_dir()
            and len(p.name) == 8
            and p.name.isdigit()
            and (p / "retmap").is_dir()
        )

        # print(animal, self.days)

        self.animal_path = str(animal_path)

path = r'E:\dark_drift\data'
data_object = batchRetmap(animals, animal_dobs, path)

# make animal_days dictionary which houses each retmap imaging session for all the animals
animals_days = {}
for animal in data_object.dat.keys():
    animals_days[animal] = data_object.dat[animal].days


def batch_retmap(dat_object, animals_days_dict, birth_dates=None, draw_mask=False, draw_lines = False):
    for animal in animals_days_dict.keys():
        for day in animals_days_dict[animal][:1]:
            path_to_retmap = os.path.join(dat_object.dat[animal].animal_path, day, 'retmap', 'wdf')

            config_path = Path(os.path.join(path_to_retmap, 'config.txt'))

            mode = "2"
            if birth_dates:
                run_retmap(config_path, mode, min_area=200, animal_dob=birth_dates[animal], animal=animal,
                           day=day, draw_mask=draw_mask, draw_lines = draw_lines)

    plt.close('all')


batch_retmap(data_object, animals_days, birth_dates=animal_dobs, draw_mask=False)



# Define colors for each group
group_colors = {
    'EB': 'blue',
    'EBc': 'black',
    'LB': 'red',
    'LBc': 'green'
}

def plot_stats (group_to_plot):
    '''
    :param group_to_plot:  either 'EB' or 'LB'
    '''
    fig, ax = plt.subplots (1,2, figsize = (14,4))
    group_names = [o for o in list(group_colors.keys()) if group_to_plot in o]
    for i_subplot, metric, in enumerate([total_visual_area_days, largest_neg_patch_area_days]):

        for subject, data in metric.items():
            # Assume subject is in the format "EC_GROUP_xx" (e.g., "EC_GCaMP6s_06")
            parts = subject.split('_')
            group = parts[1]  # e.g., "GCaMP6s", "GNAT", "RD1", etc.

            if group in group_names:
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


plot_stats ('EB')
plot_stats ('LB')

# for animal in [item for item in os.listdir(os.path.join(path))]:
#     for day in [item for item in os.listdir(os.path.join(path, animal)) if
#                 os.path.isdir(os.path.join(path, animal, item))]:
#         for subfile in [item for item in os.listdir(os.path.join(path, animal, day)) if
#                         os.path.isdir(os.path.join(path, animal, day, item))]:
#             print(animal, day, subfile)
#             config_path = fr'I:\retmap\{animal}\{day}\{subfile}\config.txt'
#             mode = "2"
#             run_retmap(config_path, mode)

# parser = argparse.ArgumentParser()
# parser.add_argument("config", help="File path to config file")
# parser.add_argument("-m", "--mode", help="Choose analysis mode: 1 - Create complex fields\t 2 - Create and plot sign map")
# args = parser.parse_args()
#
# run_retmap(args.config, args.mode)


if __name__ == '__main__':

    from retmap import *

    screen = 'big100'
    project = 'ethogram'  # restored or ethogram
    over_time = False
    path = r'E:\retmap'


    if project == 'ethogram':
        animals_days = {
            'EC_phpeb_11': ['20241111'],
            'EC_phpeb_12': ['20250121'],
            'EC_phpeb_13': ['20250121'],
            'EC_GNAT_03': ['20240923'],
            'EC_GNAT_04': ['20241003'],
            'EC_GNAT_05': ['20240923'],
            'EC_GNAT_06': ['20240923']
        }

        # animals_days = {'EC_phpeb_11': ['20241113'],
        #                 # 'EC_GCaMP6s_12': ['20250123'],
        #                 'EC_phpeb_13': ['20250123'],
        #                 'EC_GNAT_03': ['20240924'],
        #                 'EC_GNAT_05': ['20240924'],
        #                 'EC_GNAT_06': ['20240924'],
        #                 'EC_GNAT_04': ['20241004'],
        #                 }

    elif project == 'restored' and over_time:  # looking across rd1s and gnats, over time
        animals_days = {
            'EC_GCaMP6s_06': ['20241121', '20250408'],  # maybe remove day 1
            'EC_GCaMP6s_09': ['20241009', '20241121', '20250408'],
            'EC_GNAT_06': ['20240923', '20250404'],
            'EC_GNAT_03': ['20240923', '20250408'],
            # 'EC_GNAT_05': ['20240923', '20250408'],
            'EC_RD1_05': ['20240903', '20240918', '20241121'],
            'EC_RD1_06': ['20241001', '20241204', '20241218', '20250107', '20250408'],
            'EC_RD1_08': ['20241001', '20241204', '20241218', '20250107', '20250408'],
            'EC_RD1_09': ['20241006', '20241219', '20250109'],
            'EC_RD1_10': ['20241006', '20241219', '20250109', '20250408'],
        }
        animals_birthdate = {'EC_GCaMP6s_06': '20240624',
                             'EC_GCaMP6s_09': '20240414',
                             'EC_GNAT_06': '20240621',
                             'EC_GNAT_03': '20240621',
                             # 'EC_GNAT_05': '20240621',
                             'EC_RD1_05': '20240606',
                             'EC_RD1_06': '20240704',
                             # 'EC_RD1_07': '20240704',
                             'EC_RD1_08': '20240704',
                             'EC_RD1_09': '20240704',
                             'EC_RD1_10': '20240704'}
    elif project == 'restored':
        animals_days = {'EC_GCaMP6s_05': ['20240917'],
                        'EC_GCaMP6s_06': ['20241121'],
                        'EC_GCaMP6s_08': ['20240918'],
                        'EC_GCaMP6s_09': ['20241121'],
                        'EC_RD1_05': ['20241121'],
                        'EC_RD1_06': ['20250107'],
                        'EC_RD1_07': ['20241204'],
                        'EC_RD1_08': ['20250107'],
                        'EC_RD1_09': ['20250109'],
                        'EC_RD1_10': ['20250109'],
                        'EC_RD1opto_04': ['20241111'],
                        'EC_RD1opto_02': ['20241111'],
                        'EC_RD1opto_05': ['20241118'],
                        'EC_RD1opto_03': ['20241119'],
                        'EC_RD1opto_08': ['20250305'],
                        'EC_RD1opto_10': ['20250305'],
                        'EC_GNAT_03': ['20240923'],
                        'EC_GNAT_04': ['20241003'],
                        'EC_GNAT_05': ['20240923'],
                        'EC_GNAT_06': ['20250404']
                        }
        animals_birthdate = {'EC_GCaMP6s_05': '20240624',
                             'EC_GCaMP6s_06': '20240624',
                             'EC_GCaMP6s_08': '20240414',
                             'EC_GCaMP6s_09': '20240414',
                             'EC_GNAT_03': '20240621',
                             'EC_GNAT_04': '20240621',
                             'EC_GNAT_05': '20240621',
                             'EC_GNAT_06': '20240621',
                             'EC_RD1_05': '20240606',
                             'EC_RD1_06': '20240704',
                             'EC_RD1_07': '20240704',
                             'EC_RD1_08': '20240704',
                             'EC_RD1_09': '20240704',
                             'EC_RD1_10': '20240704',
                             'EC_RD1opto_04': '20240606',
                             'EC_RD1opto_02': '20240517',
                             'EC_RD1opto_05': '20240619',
                             'EC_RD1opto_03': '20240606',
                             'EC_RD1opto_08': '20240725',
                             'EC_RD1opto_10': '20240731',
                             }

    batch_retmap(animals_days, screen, path, draw_mask=False, draw_lines = False)

    # plt.close('all')

    # quantifying total visual area, largest negative patch, mean pixel amplitude > across diff groups
    total_visual_area_g, largest_neg_patch_area_g, mean_amplitude_g, largest_pos_patch_area_g = {}, {}, {}, {}
    for animal in animals_days:
        g = animal.split('_')[1]

        if (g not in total_visual_area_g):
            total_visual_area_g[g] = [total_visual_area[animal]]
            largest_neg_patch_area_g[g] = [largest_neg_patch_area[animal]]
            largest_pos_patch_area_g[g] = [largest_pos_patch_area[animal]]
            mean_amplitude_g[g] = [mean_amplitude[animal]]
        else:
            total_visual_area_g[g].append(total_visual_area[animal])
            largest_neg_patch_area_g[g].append(largest_neg_patch_area[animal])
            largest_pos_patch_area_g[g].append(largest_pos_patch_area[animal])
            mean_amplitude_g[g].append(mean_amplitude[animal])

    fig, ax = plt.subplots (1,4, figsize = (14,4))
    group_names = ['phpeb', 'GNAT']
    for i_subplot, metric, in enumerate([total_visual_area_g, largest_neg_patch_area_g, mean_amplitude_g, largest_pos_patch_area_g]):

        for i, group in enumerate(group_names):
            arr = np.array(metric[group])

            if group == 'phpeb':
                colour = 'black'
            if group == 'GNAT':
                colour = 'green'

            ax[i_subplot].scatter ([i]*len(arr) + np.random.uniform(-0.2, 0.2, len(arr)), arr, c = colour, alpha = 0.5, s = 30)
            ax[i_subplot].bar ([i], arr.mean(), color = colour, alpha = 0.6)

            if i_subplot == 0:
                ax[i_subplot].set_title ('total visual area')
            elif i_subplot == 1:
                ax[i_subplot].set_title ('largest negative patch')
            elif i_subplot == 2:
                ax[i_subplot].set_title ('mean pixel amplitude')
            elif i_subplot == 3:
                ax[i_subplot].set_title ('largest positive patch')

        # significance tests
        control, gnat = metric['phpeb'], metric['GNAT']
        stat, p = kruskal(control, gnat)
        print(f'KW H-statistic: {stat:.3f}, p-value: {p:.3f}')

        mw_pval = mannwhitneyu(control, gnat, alternative='two-sided').pvalue
        print(f'Mann-Whitney p value {mw_pval}')

        # if p < 0.05:  # follow up with testing pairwise comparisons
        #     # Do pairwise comparisons manually:
        #     print('mannwhitney two-sided test')
        #     pvals = [mannwhitneyu(control, gnat, alternative='two-sided').pvalue]
        #
        #     print(pvals)
        #
        #     comparisons = [('GCaMP6s', 'GNAT')] # Define comparisons in same order as pvals
        #     group_coords = {'GCaMP6s': 0, 'GNAT':1}  # Coordinates for group positions on x-axis
        #
        #     max_val = max(np.max(control), np.max(gnat)) * 1.05
        #     step = max_val * 0.05  # vertical spacing between significance bars
        #     h = max_val
        #
        #     # Loop through each comparison and plot if significant
        #     for (group1, group2), p_val in zip(comparisons, pvals):
        #         if p_val < 0.05:
        #             x1, x2 = group_coords[group1], group_coords[group2]
        #             y = h
        #             h += step  # update height for next line if needed
        #
        #             print(x1, x2, y)
        #
        #             # plot the number of stars according to significance value
        #             if p_val < 0.001:
        #                 stars = '***'
        #             elif p_val < 0.01:
        #                 stars = '**'
        #             else:
        #                 stars = '*'
        #
        #             ax[i_subplot].plot([x1, x1, x2, x2], [y, y + step / 2, y + step / 2, y], lw=1.5, c='k')
        #             ax[i_subplot].text((x1 + x2) * 0.5, y + step * 0.6, stars,
        #                                ha='center', va='bottom', fontsize=14)

        ax[i_subplot].set_xticks(range(len(np.unique([animal[3:-3] for animal in animals_days.keys()]))))
        ax[i_subplot].set_xticklabels(['phpeb','gnat'])
        ax[i_subplot].set_ylabel('# pixels')
    plt.tight_layout()
    plt.show()

groups = list(dict.fromkeys([k.split('_')[1] for k in total_visual_area.keys()]))
plt.figure(figsize = (5,4))
for i, group in enumerate(groups):
    colors = {'EB': 'blue', 'EBc': 'black', 'LBc': 'black', 'LB': 'red'}
    group_dat = []
    for animal in [a for a in total_visual_area if group +'_' in a]:
        area = total_visual_area[animal]
        group_dat.append(area)

        plt.scatter(i, area, color=colors[group], alpha = 0.5)
    # print(group, np.array(group_dat))
    plt.bar(i, np.mean(np.array(group_dat)), color=colors[group])
plt.show()


for d in [total_visual_area, largest_neg_patch_area, largest_pos_patch_area]:

    from scipy.stats import sem

    # consistent order
    groups = ['EBc', 'EB', 'LBc', 'LB']

    colors = {
        'EB': 'royalblue',
        'EBc': 'black',
        'LBc': 'dimgray',
        'LB': 'firebrick'
    }

    fig, ax = plt.subplots(figsize=(5, 4))

    for i, group in enumerate(groups):

        # get values for this group
        group_dat = [
            d[a]
            for a in d
            if f'_{group}_' in a
        ]

        group_dat = np.array(group_dat)

        # jittered x positions
        x_jitter = np.random.normal(i, 0.06, size=len(group_dat))

        # scatter points
        ax.scatter(
            x_jitter,
            group_dat,
            color=colors[group],
            alpha=0.6,
            s=70,
            edgecolor='white',
            linewidth=0.8,
            zorder=3
        )

        # mean ± SEM
        mean = np.mean(group_dat)
        error = sem(group_dat)

        ax.errorbar(
            i,
            mean,
            yerr=error,
            fmt='_',
            color='black',
            capsize=4,
            markersize=25,
            linewidth=2.5,
            zorder=4
        )

    # formatting
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups)

    ax.set_ylabel('Total Visual Area')
    ax.set_xlim(-0.5, len(groups)-0.5)

    # cleaner style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(axis='both', labelsize=11)

    plt.tight_layout()
    plt.show()


#########
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem, ttest_ind

for d, ylabel in zip(
    [total_visual_area, largest_neg_patch_area, largest_pos_patch_area],
    ['Total Visual Area', 'Largest Negative Patch', 'Largest Positive Patch']
):

    groups = ['EBc', 'EB', 'LBc', 'LB']

    colors = {
        'EB': 'royalblue',
        'EBc': 'black',
        'LBc': 'black',
        'LB': 'firebrick'
    }

    fig, ax = plt.subplots(figsize=(4, 4))

    # store data for stats
    group_values = {}

    for i, group in enumerate(groups):

        # get values for this group
        group_dat = np.array([
            d[a]
            for a in d
            if f'_{group}_' in a
        ])

        group_values[group] = group_dat

        # jittered x positions
        x_jitter = np.random.normal(i, 0.06, size=len(group_dat))

        # scatter points
        ax.scatter(
            x_jitter,
            group_dat,
            color=colors[group],
            alpha=0.6,
            s=70,
            edgecolor='white',
            linewidth=0.8,
            zorder=3
        )

        # mean ± SEM
        mean = np.mean(group_dat)
        error = sem(group_dat)

        ax.errorbar(
            i,
            mean,
            yerr=error,
            fmt='_',
            color='black',
            capsize=4,
            markersize=25,
            linewidth=2.5,
            zorder=4
        )

    # -------------------------
    # STATISTICS
    # -------------------------

    # Welch's t-tests
    stat_eb, p_eb = ttest_ind(
        group_values['EB'],
        group_values['EBc'],
        equal_var=False
    )

    stat_lb, p_lb = ttest_ind(
        group_values['LB'],
        group_values['LBc'],
        equal_var=False
    )

    print(f'\n{ylabel}')
    print(f'EB vs EBc: p = {p_eb:.4f}')
    print(f'LB vs LBc: p = {p_lb:.4f}')

    # -------------------------
    # SIGNIFICANCE BARS
    # -------------------------

    ymax = max([
        np.max(v) for v in group_values.values()
    ])

    y_range = ymax - min([
        np.min(v) for v in group_values.values()
    ])

    def add_sig_bar(ax, x1, x2, y, p):

        # significance label
        if p < 0.001:
            sig = '***'
        elif p < 0.01:
            sig = '**'
        elif p < 0.05:
            sig = '*'
        else:
            sig = 'n.s.'

        ax.plot(
            [x1, x1, x2, x2],
            [y, y + 0.03*y_range,
             y + 0.03*y_range, y],
            c='black',
            lw=1.5
        )

        ax.text(
            (x1 + x2)/2,
            y + 0.04*y_range,
            sig,
            ha='center'
        )

    # EBc vs EB
    add_sig_bar(ax, 0, 1, ymax + 0.08*y_range, p_eb)

    # LBc vs LB
    add_sig_bar(ax, 2, 3, ymax + 0.18*y_range, p_lb)

    # formatting
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups)

    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.5, len(groups)-0.5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(axis='both', labelsize=11)

    plt.tight_layout()
    plt.show()
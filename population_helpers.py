from imports import *
from helpers import *
from suite2p_preprocessing import *


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
            s = max(0, int(start) + fps) # exclude the first second wehere there is a grey screen
            e = min(n_timepoints, int(end))
            if e > s:
                ori_per_frame[s:e] = deg
                segments.append((s, e, deg, key))

    return ori_per_frame, segments


def make_epochs(n_time, epoch_len= 10):
    '''
    :param n_time:
    :param epoch_len: n_frames to chunk each epoch
    :param step_s:
    :return:
    '''
    L = int(round(epoch_len))
    starts = np.arange(0, n_time - L + 1, L)
    return np.array([np.array((s, s+L)) for s in starts])

def frames_epochs (experiment_object, epoch_len = 10, grating = True):
    '''
    :param traces: array shape n_cells x n_timepoints
    :param epoch_len: n frames in each epoch
    :param grating: whether or not we want to analyze grating data
    :return:
    '''

    traces = experiment_object['zscored_traces']  # (n_cells, n_timepoints)
    n_time = traces.shape[1]

    epoch_idx = make_epochs(traces.shape[1], epoch_len=epoch_len)

    if grating:

        ori_per_frame, segments = build_ori_per_frame(
            dict_stim_ttls=experiment_object['dict_stim_ttls'],
            dict_stim_names=experiment_object['dict_stim_names'],
            n_timepoints=n_time,
            collapse_180=True,  # set True if you want pure orientation (0≡180, etc.)
            fps = 20
        )

        traces_epochs = []
        ori_epochs = []
        for s,e in epoch_idx:
            ori_slice = ori_per_frame[s:e]
            if (len(np.unique(ori_slice)) == 1) and ~np.isnan(ori_slice).any(): # if there is only 1 orientation represented in the chunk and there are NO nans in the chunk
                traces_epochs.append(traces[:, s:e])
                ori_epochs.append(ori_slice)

        traces_epochs = np.array(traces_epochs)
        ori_epochs = np.array(ori_epochs)

        #ori_epochs = (ori_epochs % 180.0).astype(float) # merge direction paris together. 0 & 180 degree are the same

        # output is:shape (n_epochs, n_cells) & shape (n_epochs)
        return traces_epochs.mean(axis = -1), (ori_epochs.mean(axis = -1) % 180.0).astype(float)

    else:
        traces_epochs = []
        for s, e in epoch_idx:
            traces_epochs.append(traces[:, s:e])
        traces_epochs = np.array(traces_epochs)

        # output is:shape (n_epochs, n_cells)
        return traces_epochs.mean(axis=-1)


def make_ori_projection (experiment_object, grat = True):

    if grat:
        # traces_epochs, ori_epochs = frames_epochs(experiment_object, epoch_len=20, grating=True)

        #if each epoch is an actual trial > embedding population vectors
        traces_epochs = experiment_object['param_matrix_whole_zscore'][...,20:].mean(axis=-1)  # shape n_repeats, n_orientations, n_cells, remove first second (20 frames) where screen grey
        ori_epochs = experiment_object['orientations'][np.tile(np.arange(traces_epochs.shape[1]), traces_epochs.shape[0])] #% 180
        repeat_epochs = np.repeat(np.arange(traces_epochs.shape[0]), traces_epochs.shape[1])  # COLOR ACCORDING TO REPEAT NUM
        traces_epochs = traces_epochs.reshape(-1, traces_epochs.shape[-1])  # shape n_repeats * n_orientations , n_cells

        # print(repeat_epochs)

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

    scores = pca.fit_transform(normalized_traces)  # (n_epochs (1 trial, 1 orientation), n_components)
    #print(scores.shape)
    epoch_idx = np.arange(len(scores))  # 0..n_epochs-1
    norm_t = epoch_idx / epoch_idx.max()

    return traces_epochs, ori_epochs,repeat_epochs, scores, pca

def plot_umap(ax, feats, colors=None, title='UMAP', n_neighbors=15, min_dist=0.1, metric='cosine', random_state=0):
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=random_state)
    emb = reducer.fit_transform(feats)  # (n_epochs, 2)
    sc = ax.scatter(emb[:, 0], emb[:, 1], c=colors if colors is not None else 'tab:blue',cmap = 'Paired', s=12)
    # colors are : 0, 45, 90, 135, 180, 225, 270, 315
    # i want colors to correspond to :
    #0 & 180: dodgerblue, mediumblue
    # 45 & 225: violet, darkviolet
    # 90 & 270: lightgreen, green
    # 135 & 315: salmon, orangered
    ax.set_xlabel('UMAP-1')
    ax.set_ylabel('UMAP-2')
    ax.set_title(title)
    return sc

# def plot_umap(ax, feats, colors=None, title='UMAP', n_neighbors=15, min_dist=0.1, metric='cosine', random_state=0):
#     reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=random_state)
#     emb = reducer.fit_transform(feats)  # (n_epochs, 2)
#
#     if colors is not None:
#         ori_color_map = {
#             0: 'dodgerblue',
#             180: 'mediumblue',
#             45: 'violet',
#             225: 'darkviolet',
#             90: 'lightgreen',
#             270: 'green',
#             135: 'salmon',
#             315: 'red'
#         }
#
#         # ori_color_map = {
#         #     0: 'black',
#         #     1: 'darkviolet',
#         #     2: 'violet',
#         #     3: 'mediumblue',
#         #     4: 'blue',
#         #     5: 'green',
#         #     6: 'lightgreen',
#         #     7: 'yellow',
#         #     8: 'orange',
#         #     9: 'red',
#         #     10: 'salmon'
#         # }
#
#         point_colors = [ori_color_map[o] for o in colors]
#     else:
#         point_colors = 'tab:blue'
#
#     sc = ax.scatter(emb[:, 0], emb[:, 1], c=point_colors, s=12)
#     ax.set_xlabel('UMAP-1')
#     ax.set_ylabel('UMAP-2')
#     ax.set_title(title)
#
#     return sc, ori_color_map

import numpy as np
import umap
import matplotlib.pyplot as plt

def plot_umap(ax, feats, colors=None, title='UMAP', n_neighbors=15, min_dist=0.1, metric='cosine', random_state=0, plot_time=False):

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state
    )
    emb = reducer.fit_transform(feats)  # (n_epochs, 2)

    ori_color_map = None

    if colors is not None:
        if plot_time:
            sc = ax.scatter(emb[:, 0], emb[:, 1], c=colors, cmap='plasma', s=12)
        else:
            ori_color_map = {
                0: 'dodgerblue',
                180: 'mediumblue',
                45: 'violet',
                225: 'darkviolet',
                90: 'lightgreen',
                270: 'green',
                135: 'salmon',
                315: 'red'
            }
            point_colors = [ori_color_map[int(o)] for o in colors]
            sc = ax.scatter(emb[:, 0], emb[:, 1], c=point_colors, s=12)
    else:
        sc = ax.scatter(emb[:, 0], emb[:, 1], c='tab:blue', s=12)

    ax.set_xlabel('UMAP-1')
    ax.set_ylabel('UMAP-2')
    ax.set_title(title)

    return sc, ori_color_map

def plot_ori_projection(obj, animal, grat=True, n_components_plot = 3, plot_time = False):
    print(animal)
    days = list(obj.dat[animal].dat_subject.keys())
    n_days = len(days)

    ncols = int(np.ceil(np.sqrt(n_days)))
    nrows = int(np.ceil(n_days / ncols))

    if not plot_time:
        fig_scatter = plt.figure(figsize=(5 * ncols, 5 * nrows))
        fig_scree, ax_scree = plt.subplots(figsize=(4,4))
    fig_umap, ax_umap = plt.subplots(nrows = 1,ncols = n_days, figsize=(4*n_days, 3))

    if n_days == 1:
        ax_umap = [ax_umap]

    for i, day in enumerate(days):
        postnatal_day = calculate_animal_age(obj.animal_dobs[animal], day)
        # cmap = plt.cm.Paired
        unique_oris = np.arange(0,360,45)

        cmap_days = plt.cm.plasma
        colors_days = [cmap_days(i / max(1, n_days - 1)) for i in range(n_days)]

        traces_epochs, ori_epochs, repeat_epochs, scores, pca_obj = make_ori_projection(
            obj.dat[animal].dat_subject[day]['grat'], grat=True
        )

        explained_var_ratio, cum, n_keep = scree_info(pca_obj, thresh=0.80)
        #ax_scree = fig_scree.add_subplot(nrows, ncols, i + 1)

        if not plot_time:
            plot_scree(ax_scree, explained_var_ratio, cum, n_keep,  title=f'{animal}', label = f'P{postnatal_day}', color = colors_days[i])

        feats = scores[:, :n_keep] # shape pop_vec, n_components kept

        first_last = first_last_distances(feats, ori_epochs, repeat_epochs)

        # print("First vs last trial distances:")
        # for ori, d in first_last.items():
        #     print(f"{ori}°: {d:.3f}")

        opp_dists = opposing_orientation_distances(feats, ori_epochs)

        # print("\nOpposing orientation distances:")
        # for pair, d in opp_dists.items():
        #     print(f"{pair}: {d:.3f}")

        spread_dist = within_cluster_spread(feats, ori_epochs, metric="euclidean")

        # print(ori_epochs), cmap previously ori_color_map
        if plot_time:
            colour_mapping = repeat_epochs
        else:
            colour_mapping = ori_epochs

        sc, ori_color_map = plot_umap(ax_umap[i], feats, colors=colour_mapping, title=f'P{postnatal_day}', n_neighbors=15, min_dist=0.1, metric='cosine',random_state=0, plot_time = plot_time)

        # sc, ori_color_map = plot_umap(
        #     ax_umap[i],
        #     feats,
        #     colors=cmap,
        #     title=f'P{postnatal_day}',
        #     n_neighbors=15,
        #     min_dist=0.1,
        #     metric='cosine',
        #     random_state=0,
        #     plot_time=True
        # )

        # print(ori_color_map)
        if not plot_time:
            if n_components_plot == 3:
                ax_sc = fig_scatter.add_subplot(nrows, ncols, i + 1, projection='3d')
                sc = ax_sc.scatter(
                    scores[:, 0], scores[:, 1], scores[:, 2],
                    c=(ori_epochs if grat else 'blue'),
                    cmap=('plasma' if grat else None), s=30
                )
                ax_sc.set_zlabel('PC3')
                ax_sc.view_init(elev=25, azim=40)
                ax_sc.set_box_aspect((1, 1, 1))
            else:
                ax_sc = fig_scatter.add_subplot(nrows, ncols, i + 1)
                sc = ax_sc.scatter(
                    scores[:, 0], scores[:, 1],
                    c=(ori_epochs if grat else 'blue'),
                    cmap=('plasma' if grat else None), s=30
                )

            ax_sc.set_xlabel('PC1'); ax_sc.set_ylabel('PC2')
            ax_sc.set_title(f'{animal} — {day} ({"Gratings" if grat else "Spontaneous"})')
            if grat:
                fig_scatter.colorbar(sc, ax=ax_sc, pad=0.05, shrink=0.7, label='Orientation (deg)')

            if n_components_plot == 3:
                ax_sc.set_zlabel('PC3')
                ax_sc.view_init(elev=25, azim=40)
                ax_sc.set_box_aspect((1, 1, 1))

    if plot_time:
        cbar = fig_umap.colorbar(sc, ax=ax_umap, shrink=0.8, pad=0.02)
        cbar.set_label('Repeat')
    else:
        cmap = mcolors.ListedColormap([ori_color_map[int(o)] for o in unique_oris])
        cbar = fig_umap.colorbar(
            plt.cm.ScalarMappable(
                cmap=cmap,
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

    if not plot_time:
        fig_scatter.tight_layout()
        fig_scree.tight_layout()
    fig_umap.suptitle(f'{animal} UMAP')

    folder_path = os.path.join(os.path.dirname(obj.path), 'figures', f'umap_gratings')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'umap_{animal}{"_time"*plot_time}.svg'))
    plt.savefig(os.path.join(folder_path, f'umap_{animal}{"_time"*plot_time}.png'))


    plt.show(block=False)  # shows both figures without blocking

    return first_last, opp_dists, spread_dist

# #
# dict_first_last_dist = {} #0, 45, 90, 135, 180, 225, 270,315
# dict_inter_cluster_distance = {} #0-180, 45-225, 90-270, 135-315
# dict_within_cluster_spread = {} #00, 45, 90, 135, 180, 225, 270,315
# for animal in data_object.dat.keys():
#     group = animal.split('_')[1]
#     first_last_dist, inter_cluster_distance, spread_dist = plot_ori_projection(data_object, animal)
#     # dict_first_last_dist[group]
#     # dict_inter_cluster_distance[group]
#     plt.close('all')
#
#     if group not in dict_first_last_dist:
#         dict_first_last_dist[group] = []
#         dict_first_last_dist[group].append(list(first_last_dist.values()))
#     else:
#         dict_first_last_dist[group].append(list(first_last_dist.values()))
#
#     if group not in dict_inter_cluster_distance:
#         dict_inter_cluster_distance[group] = []
#         dict_inter_cluster_distance[group].append(list(inter_cluster_distance.values()))
#     else:
#         dict_inter_cluster_distance[group].append(list(inter_cluster_distance.values()))
#
#     if group not in dict_within_cluster_spread:
#         dict_within_cluster_spread[group] = []
#         dict_within_cluster_spread[group].append(list(spread_dist.values()))
#     else:
#         dict_within_cluster_spread[group].append(list(spread_dist.values()))
#
#     plt.close('all')
#
# # dict_first_last_dist: each entry is a group, with shape (n_animals, 8 directions) > mean across directions
# # dict_inter_cluster_distance: each entry is a group, with shape (n_animals, 4 comparisons) > mean across comparisons
# for group in dict_first_last_dist:
#     dict_first_last_dist[group] = np.array(dict_first_last_dist[group]).mean(axis= -1)
#     dict_inter_cluster_distance[group] = np.array(dict_inter_cluster_distance[group]).mean(axis= -1)
#     dict_within_cluster_spread[group] = np.array(dict_within_cluster_spread[group]).mean(axis=-1)
#
# from scipy.stats import mannwhitneyu, ttest_ind
#
# print("\n--- STATISTICAL TESTS ---")
#
# # EB vs EBc
# eb = dict_inter_cluster_distance['EB']
# ebc = dict_inter_cluster_distance['EBc']
#
# u_stat, u_p = mannwhitneyu(eb, ebc, alternative='two-sided')
# t_stat, t_p = ttest_ind(eb, ebc, equal_var=False)
#
# print('EB vs EBc')
# print(f'    Mann-Whitney U: stat={u_stat:.3f}, p={u_p:.3e}')
# print(f'    Welch t-test:  stat={t_stat:.3f}, p={t_p:.3e}')
#
#
# # LB vs LBc
# lb = dict_inter_cluster_distance['LB']
# lbc = dict_inter_cluster_distance['LBc']
#
# u_stat, u_p = mannwhitneyu(lb, lbc, alternative='two-sided')
# t_stat, t_p = ttest_ind(lb, lbc, equal_var=False)
#
# print('LB vs LBc')
# print(f'    Mann-Whitney U: stat={u_stat:.3f}, p={u_p:.3e}')
# print(f'    Welch t-test:  stat={t_stat:.3f}, p={t_p:.3e}')

from scipy.spatial.distance import cdist
import numpy as np

def first_last_distances(feats, ori_labels, repeat_labels, metric="euclidean"):
    ori_labels = np.asarray(ori_labels)
    repeat_labels = np.asarray(repeat_labels)

    unique_oris = np.sort(np.unique(ori_labels))
    dists = {}

    for ori in unique_oris:
        mask = ori_labels == ori

        feats_ori = feats[mask]
        repeats_ori = repeat_labels[mask]

        # find first and last trial
        first_idx = np.argmin(repeats_ori)
        last_idx = np.argmax(repeats_ori)

        f = feats_ori[first_idx][None, :]
        l = feats_ori[last_idx][None, :]

        dist = cdist(f, l, metric=metric)[0, 0]
        dists[int(ori)] = dist

    return dists

from scipy.spatial.distance import cdist
import numpy as np

def within_cluster_spread(feats, ori_labels, metric="euclidean"):
    ori_labels = np.asarray(ori_labels)

    unique_oris = np.sort(np.unique(ori_labels))
    spread = {}

    for ori in unique_oris:
        mask = ori_labels == ori
        feats_ori = feats[mask]

        # centroid of this orientation cluster
        centroid = feats_ori.mean(axis=0, keepdims=True)

        # distance of each trial to centroid
        dists = cdist(feats_ori, centroid, metric=metric).ravel()

        # average spread around centroid
        spread[int(ori)] = np.mean(dists)

    return spread

from scipy.spatial.distance import cdist

def opposing_orientation_distances(feats, ori_labels, metric="euclidean"):
    ori_labels = np.asarray(ori_labels)
    unique_oris = np.sort(np.unique(ori_labels))

    # compute centroids
    centroids = {
        ori: feats[ori_labels == ori].mean(axis=0)
        for ori in unique_oris
    }

    pairs = [(0,180), (45,225), (90,270), (135,315)]
    dists = {}

    for a, b in pairs:
        if a in centroids and b in centroids:
            d = cdist(
                centroids[a][None, :],
                centroids[b][None, :]
            )[0, 0]
            dists[f"{a}-{b}"] = d

    return dists

def plot_group_summary(metric_dict, ylabel, title):
    groups = list(metric_dict.keys())
    x = np.arange(len(groups))

    plt.figure(figsize=(6,4))

    for i, group in enumerate(groups):
        vals = metric_dict[group]

        # jittered individual animals
        jitter = np.random.normal(0, 0.05, size=len(vals))
        plt.scatter(np.full_like(vals, i) + jitter, vals, alpha=0.6)

        # mean + SEM
        mean = np.mean(vals)
        sem = np.std(vals) / np.sqrt(len(vals))

        plt.errorbar(i, mean, yerr=sem, fmt='o', capsize=5, color='black')

    plt.xticks(x, groups)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.show()
#
# plot_group_summary(
#     dict_first_last_dist,
#     ylabel="Mean first-last distance",
#     title="Within-orientation drift (averaged across orientations)"
# )
#
# plot_group_summary(
#     dict_inter_cluster_distance,
#     ylabel="Mean opposing-orientation distance",
#     title="Orientation separation (averaged across pairs)"
# )

############## plot the gratings projection across dyas in the same space ############

def make_ori_projection_days (dat_obj, animal):

    # if each epoch is an actual trial > embedding population vectors

    # shape (n_days, n_repeats, n_ori, n_cells) > took mean across moving timepoint
    traces_epochs_days = np.array([dat_obj.dat[animal].dat_subject[day]['grat']['param_matrix_whole_zscore'][...,20:].mean(axis=-1) for day in dat_obj.dat[animal].dat_subject.keys()])

    # shape (n_days, n_repeats*n_ori) > grating/day labels
    ori_epochs_days = np.array([dat_obj.dat[animal].dat_subject[day]['grat']['orientations'][np.tile(np.arange(traces_epochs_days.shape[2]), traces_epochs_days.shape[1])] for day in dat_obj.dat[animal].dat_subject.keys()])
    repeat_epochs_days = np.array([np.repeat(np.arange(traces_epochs_days.shape[1]), traces_epochs_days.shape[2]) for day in dat_obj.dat[animal].dat_subject.keys()])
    day_label_days = np.array([np.repeat(day, traces_epochs_days.shape[1]*traces_epochs_days.shape[2]) for day in range(len(dat_obj.dat[animal].dat_subject.keys()))])

    # shape to be (n_days*n_repeats*n_orientations, n_cells)
    traces_epochs_days = traces_epochs_days.reshape(traces_epochs_days.shape[0] * traces_epochs_days.shape[1] * traces_epochs_days.shape[2], -1)

    # shape to be (n_days*n_repeats*n_orientations)
    ori_epochs_days = ori_epochs_days.reshape(ori_epochs_days.shape[0] * ori_epochs_days.shape[1])
    repeat_epochs_days = repeat_epochs_days.reshape(repeat_epochs_days.shape[0] * repeat_epochs_days.shape[1])
    day_label_days = day_label_days.reshape(day_label_days.shape[0] * day_label_days.shape[1])
    # might have to flatten other arrays!!!!!

    # traces_epochs = experiment_object['param_matrix_whole_zscore'][...,20:].mean(axis=-1)  # shape n_repeats, n_orientations, n_cells, remove first second (20 frames) where screen grey
    # ori_epochs = experiment_object['orientations'][np.tile(np.arange(traces_epochs.shape[1]), traces_epochs.shape[0])] #% 180
    # repeat_epochs = np.repeat(np.arange(traces_epochs.shape[0]), traces_epochs.shape[1])  # COLOR ACCORDING TO REPEAT NUM
    # traces_epochs = traces_epochs.reshape(-1, traces_epochs.shape[-1])  # shape n_repeats * n_orientations , n_cells

    #remove each epoch’s global mean across cells
    traces_epochs = traces_epochs_days - traces_epochs_days.mean(axis=1, keepdims=True)

    norms = np.linalg.norm(traces_epochs, axis=1, keepdims=True)
    normalized_traces = traces_epochs / np.maximum(norms, 1e-12)

    pca = PCA(svd_solver='full') #keep all the PCs

    scores = pca.fit_transform(normalized_traces)  # (n_epochs (1 trial, 1 orientation), n_components)
    #print(scores.shape)
    epoch_idx = np.arange(len(scores))  # 0..n_epochs-1
    norm_t = epoch_idx / epoch_idx.max()

    return traces_epochs, ori_epochs_days, repeat_epochs_days, day_label_days, scores, pca

def plot_umap_days(ax, feats,  colors_mapping=None, title='UMAP', n_neighbors=15, min_dist=0.1, metric='cosine', random_state=0, plot_time=False):

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state
    )
    emb = reducer.fit_transform(feats)  # (n_epochs, 2)

    ori_color_map = None

    if colors_mapping is not None:
        if plot_time:

            import matplotlib.colors as mcolors

            def truncate_colormap(cmap, minval=0.0, maxval=0.95, n=256):
                return mcolors.LinearSegmentedColormap.from_list(
                    f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',
                    cmap(np.linspace(minval, maxval, n))
                )

            cmap_trunc = truncate_colormap(plt.cm.plasma, 0.0, 0.95)
            cmap_trunc = truncate_colormap(plt.cm.viridis, 0.0, 0.95)

            sc = ax.scatter(emb[:, 0], emb[:, 1], c=colors_mapping, cmap=cmap_trunc, s=12)
        else:
            ori_color_map = {
                0: 'dodgerblue',
                180: 'mediumblue',
                45: 'violet',
                225: 'darkviolet',
                90: 'lightgreen',
                270: 'green',
                135: 'salmon',
                315: 'red'
            }
            point_colors = [ori_color_map[int(o)] for o in colors_mapping]

            sc = ax.scatter(emb[:, 0], emb[:, 1], c=point_colors, s=12)

    else:
        sc = ax.scatter(emb[:, 0], emb[:, 1], c='tab:blue', s=12)

    ax.set_xlabel('UMAP-1')
    ax.set_ylabel('UMAP-2')
    ax.set_title(title)

    return sc, ori_color_map


def plot_umap_days_2(
    ax, feats, colors_mapping=None, days_mapping=None,
    title='UMAP', n_neighbors=15, min_dist=0.1,
    metric='cosine', random_state=0, plot_time=False):

    '''
    colors are orientation, but size and alpha epends on the day
    '''

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state
    )
    emb = reducer.fit_transform(feats)

    ori_color_map = None

    # default point size / alpha
    point_sizes = 12
    point_alphas = np.ones(emb.shape[0])

    # map days to size + alpha
    if days_mapping is not None:
        days_mapping = np.asarray(days_mapping)
        unique_days = np.unique(days_mapping)

        sizes = np.linspace(10, 80, len(unique_days))
        alphas = np.linspace(1, 0.1, len(unique_days))  # early faint, late solid

        day_to_size = {d: s for d, s in zip(unique_days, sizes)}
        day_to_alpha = {d: a for d, a in zip(unique_days, alphas)}

        point_sizes = np.array([day_to_size[d] for d in days_mapping])
        point_alphas = np.array([day_to_alpha[d] for d in days_mapping])

    if colors_mapping is not None:

        if plot_time:
            import matplotlib.colors as mcolors

            def truncate_colormap(cmap, minval=0.0, maxval=0.95, n=256):
                return mcolors.LinearSegmentedColormap.from_list(
                    f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',
                    cmap(np.linspace(minval, maxval, n))
                )

            cmap_trunc = truncate_colormap(plt.cm.plasma, 0.0, 0.95)
            # cmap_trunc = truncate_colormap(plt.cm.viridis, 0.0, 0.95)

            sc = ax.scatter(
                emb[:, 0],
                emb[:, 1],
                c=colors_mapping,
                cmap=cmap_trunc,
                s=point_sizes,
                alpha=0.8
            )

        else:
            import matplotlib.colors as mcolors

            ori_color_map = {
                0: 'dodgerblue',
                180: 'mediumblue',
                45: 'violet',
                225: 'darkviolet',
                90: 'lightgreen',
                270: 'green',
                135: 'salmon',
                315: 'red'
            }

            point_colors = [ori_color_map[int(o)] for o in colors_mapping]

            point_colors_rgba = [
                (*mcolors.to_rgba(c)[:3], a)
                for c, a in zip(point_colors, point_alphas)
            ]

            sc = ax.scatter(
                emb[:, 0],
                emb[:, 1],
                c=point_colors_rgba,
                s=point_sizes
            )

    else:
        sc = ax.scatter(
            emb[:, 0],
            emb[:, 1],
            c='tab:blue',
            s=point_sizes,
            alpha=0.8
        )

    ax.set_xlabel('UMAP-1')
    ax.set_ylabel('UMAP-2')
    ax.set_title(title)

    return sc, ori_color_map

def plot_ori_projection_days(obj, animal, plot_time = False, multi_dimension = False):

    fig_umap, ax_umap = plt.subplots(nrows = 1,ncols = 1, figsize=(6, 5))
    ax_umap = ax_umap
    n_days = len(list(obj.dat[animal].dat_subject.keys()))

    # cmap = plt.cm.Paired
    unique_oris = np.arange(0,360,45)

    max_val = 0.8  # or 0.85 if you want even softer
    colors_days = [plt.cm.plasma(i / (n_days - 1) * max_val) for i in range(n_days)]

    traces_epochs, ori_epochs, repeat_epochs, day_label_days, scores, pca_obj = make_ori_projection_days(obj, animal)

    explained_var_ratio, cum, n_keep = scree_info(pca_obj, thresh=0.80)

    feats = scores[:, :n_keep] # shape pop_vec, n_components kept

    if plot_time:
        colour_mapping = repeat_epochs#day_label_days # repeat_epochs
    else:
        colour_mapping = ori_epochs

    # feats = feats [:80]
    # colour_mapping = colour_mapping[:80]

    if multi_dimension:
        sc, ori_color_map = plot_umap_days_2(
            ax_umap,
            feats,
            colors_mapping=ori_epochs,  # color = orientation
            days_mapping=day_label_days,  # size/alpha = day
            title='map',
            n_neighbors=15,
            min_dist=0.1,
            metric='cosine',
            random_state=0,
            plot_time=False
        )
    else:
        sc, ori_color_map = plot_umap_days(ax_umap, feats, colors_mapping=colour_mapping, n_neighbors=15, min_dist=0.1, metric='cosine',random_state=0, plot_time = plot_time)

    if plot_time:
        cbar = fig_umap.colorbar(sc, ax=ax_umap, shrink=0.8, pad=0.02)
        cbar.set_label('Day')
    else:
        cmap = mcolors.ListedColormap([ori_color_map[int(o)] for o in unique_oris])
        cbar = fig_umap.colorbar(
            plt.cm.ScalarMappable(
                cmap=cmap,
                norm=plt.Normalize(vmin=0, vmax=len(unique_oris) - 1)
            ),
            ax=ax_umap,
            shrink=0.8,
            pad=0.02
        )

        # --- Properly label the discrete ticks ---
        cbar.set_ticks(np.arange(len(unique_oris)))
        cbar.set_ticklabels(unique_oris.astype(int))
        cbar.set_label('Orientation (°)')

    fig_umap.suptitle(f'{animal} UMAP')

    folder_path = os.path.join(os.path.dirname(obj.path), 'figures', f'umap_gratings_days')
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'umap_{animal}{"_time"*plot_time}{"_multi"*multi_dimension}.svg'))
    plt.savefig(os.path.join(folder_path, f'umap_{animal}{"_time"*plot_time}{"_multi"*multi_dimension}.png'))

    plt.show(block=False)  # shows both figures without blocking



def plot_avg_response_repeats (object, group_type):
    '''
    :param object:
    :param recording_day (str): PXX, postnatal day we want to analyze
    :param group_type (str): EB or LB or EBc or LBc
    :return:
    '''

    group_type += '_'
    # lets do one group type at a time

    # so each entry in response is group, then postnatal day, then array of shape (n_animals, n_repeats, n_cells)
    response = {k: {} for k in np.unique([a.split ('_')[1] for a in object.list_animals if group_type in a ])}
    # response = []
    sub_list_animals = [animal for animal in object.list_animals if group_type in animal]

    for animal in sub_list_animals:
        group = animal.split('_')[1]

        for day in object.dat[animal].dat_subject.keys():

            postnatal_day = calculate_animal_age(object.animal_dobs[animal], day)
            responsive_indices = object.dat[animal].dat_subject[day]['grat']['thresholded_cells']

            # originally, shape n_repeats, n_ori, n_cells, n_timepoints, after mean shape is (n_repeats, n_timepoints)
            matrix = object.dat[animal].dat_subject[day]['grat']['zscored'][:, :, responsive_indices, :].mean(axis = (1,2))

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
                data = np.array(response[group][postnatal_day])  # shape: n_animals, n_repeats, n_timepoints
                all_group_data[group] = data

                mean_dat = data.mean(axis = 0) # average across all animals, shape (n_repeats, n_timepoints)

                n_repeats = mean_dat.shape[0]
                colors = plt.cm.plasma(np.linspace(0, 0.9, n_repeats))  # avoid bright yellow

                for i in range(n_repeats):
                    ax[i_day].plot(mean_dat[i], color=colors[i])

                # plot individual animals
                # for animal in data:
                #     ax[i_day].plot(np.arange(data.shape[1]), animal, color=colours.get(group, 'gray'), alpha=0.2)
                #
                # plot group average
                # ax[i_day].plot(np.arange(data.shape[1]), data.mean(axis=0), color=colours.get(group, 'gray'), label=group)

        #ax[i_day].set_ylim ([min(ctrl.min(), dark.min()), max(ctrl.max(), dark.max()) + 0.25 * max(ctrl.max(), dark.max())])
        ax[i_day].set_ylabel('Response (z-scored)')
        # ax[i_day].set_xticks(np.arange(ctrl.shape[1])[::object.fps])
        # ax[i_day].set_xticklabels( [int(a/object.fps) for a in (np.arange(ctrl.shape[1])[::object.fps])])
        # ax[i_day].axvline(x=object.fps, color='k', linestyle=':', linewidth=1)
        # ax[i_day].axvspan(object.fps * 2, object.fps * 5, color='gray', alpha=0.2)
        ax[i_day].set_title(f"Postnatal Day {postnatal_day}")
        ax[i_day].legend()

    ax[-1].set_xlabel("Time (s)")
    plt.tight_layout()

    # folder_path = os.path.join(os.path.dirname(object.path), 'figures', 'population_all_animals', 'avg_response')
    # if not os.path.exists(folder_path):
    #     os.makedirs(folder_path)
    # plt.savefig(os.path.join(folder_path, f'avg_response {group_type}.svg'))
    # plt.savefig(os.path.join(folder_path, f'avg_response {group_type}.png'))

    plt.show()


def responsive_cells_across_repeats(object, group):

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


def responsive_cells(obj):

    groups = np.unique([a.split('_')[1] for a in obj.list_animals])

    # {'EB': {animal 1: [], animal 2:[]}],'EBc': {animal 1: [] ...}, 'LB':{animal 1: [] ...}, 'LBc': {animal 1: [] ...}}
    percent_responsive_groups = {g: {} for g in groups}

    for i_animal, animal in enumerate(obj.dat):

        group = animal.split('_')[1]
        days_recordings = [(day, subfile) for day in obj.dat[animal].dat_subject for subfile in
                           obj.dat[animal].dat_subject[day] if 'grat' in subfile]

        # if there are no gratings for the animal, skip it
        if not days_recordings:
            continue

        p_resp = {}

        for i, (day, subfile) in enumerate(days_recordings):

            timepoint = calculate_animal_age(obj.animal_dobs[animal], day)

            # shape n_Cells, n_timepoints
            percent_cells = 100 * obj.dat[animal].dat_subject[day][subfile]['thresholded_cells'].sum() / len(obj.dat[animal].dat_subject[day][subfile]['thresholded_cells'])
            p_resp[timepoint] = percent_cells

        percent_responsive_groups [group] [animal] = p_resp

    return percent_responsive_groups


# def get_repeat_by_cell_threshold (obj, param_matrix, zscore_threshold = None):
#     '''
#     :param param_matrix: shape = ((n_repeats, n_ori, n_cells, n_timepoints))
#     :param zscore_threshold:
#     :param std_threshold:
#     :return: cell_theshold, shape ((n_cells))
#     '''
#
#     z_scored_response = zscore_baseline(obj, param_matrix) # array of responses, shape ((n_repeats, n_ori, ..., n_cells, timepoints))
#     cell_threshold = np.ones((z_scored_response.shape[-2])) * zscore_threshold
#
#     # else:
#     #     return None
#     return cell_threshold


def zscore_thresholding_repeats (object, animal, day, zscore_threshold = None, std_threshold = None):
    '''
    Determines which cells were responsive

    :param object:
    :param zscore_threshold: 2
    :param std_threshold: 1.5

    '''

    # shape (n_repeats, n_ori, n_cells, timepoints)
    z_scored_response = object.dat[animal].dat_subject[day]['grat']['param_matrix_whole_zscore']

    # shape n_cells
    cell_threshold = get_cell_by_cell_threshold(object, z_scored_response, zscore_threshold=zscore_threshold, std_threshold=std_threshold)

    # shape n_cells > take the z scored responses during the 'moving grating' period > average (n_ori, n_cells, n_timepoints) over orientations/repeats to get shape (n_cells)
    #z_scored_response_mean = z_scored_response[..., object.fps * 2: object.fps * 5].mean(axis=(0,-1))
    # shape (n_repeats, n_ori, n_cells)
    z_scored_response_mean = z_scored_response[..., object.fps * 2: object.fps * 5].mean(axis=(-1))
    z_scored_response_median = np.median(z_scored_response[..., object.fps * 2: object.fps * 5], axis=(-1))

    # shape (n_repeat, n_ori, n_cells) > (n_repeat, n_cells)
    cell_exceeds_threshold = np.any(z_scored_response_mean > cell_threshold, axis=1)

    return cell_exceeds_threshold



def responsive_cells_repeats(obj):

    groups = np.unique([a.split('_')[1] for a in obj.list_animals])

    # {'EB': {animal 1: [], animal 2:[]}],'EBc': {animal 1: [] ...}, 'LB':{animal 1: [] ...}, 'LBc': {animal 1: [] ...}}
    percent_responsive_groups = {g: {} for g in groups}

    for i_animal, animal in enumerate(obj.dat):

        group = animal.split('_')[1]
        days_recordings = [(day, subfile) for day in obj.dat[animal].dat_subject for subfile in
                           obj.dat[animal].dat_subject[day] if 'grat' in subfile]

        # if there are no gratings for the animal, skip it
        if not days_recordings:
            continue

        p_resp = {}

        for i, (day, subfile) in enumerate(days_recordings):

            timepoint = calculate_animal_age(obj.animal_dobs[animal], day)

            # shape (n_repeats, n_cells)
            responsive_cells_repeats = zscore_thresholding_repeats (obj, animal, day, zscore_threshold = 0.8, std_threshold = 0.8)

            # shape (n_repeats)
            percent_responsive_cells_repeats = np.array([100*(responsive_cells_repeats[i_repeat].sum()/responsive_cells_repeats[i_repeat].shape[0]) for i_repeat in range(responsive_cells_repeats.shape[0])])

            p_resp[timepoint] = percent_responsive_cells_repeats

        percent_responsive_groups [group] [animal] = p_resp

    return percent_responsive_groups

def group_responsive_repeats (obj):

    # dictionary: group > animal > postnatal day > array of shape n_repeats
    responsive_repeats = responsive_cells_repeats(obj)

    # collapse animal level and stack across animals within a group
    grouped = {}
    for group, animals in responsive_repeats.items():
        grouped[group] = {}

        for animal, days in animals.items():
            for day, arr in days.items():

                # initialize list for this group/day
                if day not in grouped[group]:
                    grouped[group][day] = []

                grouped[group][day].append(arr)

        # stack arrays for each day
        for day in grouped[group]:
            grouped[group][day] = np.stack(grouped[group][day], axis=0)

    groups = list(grouped.keys())
    n_groups = len(groups)

    fig, axes = plt.subplots(
        nrows=1,
        ncols=n_groups,
        figsize=(5 * n_groups, 4),
        sharey=True
    )

    axes = np.atleast_1d(axes)

    for ax, group in zip(axes, groups):

        days = sorted(grouped[group].keys())
        colors = plt.cm.plasma(np.linspace(0, 0.9, len(days)))

        for i, day in enumerate(days):
            data = grouped[group][day]   # shape (n_animals, n_repeats)

            mean_per_repeat = np.nanmean(data, axis=0)
            sem_per_repeat = np.nanstd(data, axis=0) / np.sqrt(data.shape[0])

            repeats = np.arange(data.shape[1])

            ax.plot(
                repeats,
                mean_per_repeat,
                color=colors[i],
                label=f'P{day}'
            )

            ax.fill_between(
                repeats,
                mean_per_repeat - sem_per_repeat,
                mean_per_repeat + sem_per_repeat,
                color=colors[i],
                alpha=0.2
            )

        ax.set_title(group)
        ax.set_xlabel('Repeat')
        ax.legend(title='Day', fontsize=8)

    axes[0].set_ylabel('Number of active cells')

    fig.suptitle('Active cells across repeats and days')
    plt.tight_layout()
    plt.show()



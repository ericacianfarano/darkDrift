import matplotlib.pyplot as plt
import numpy as np
from imports import *
from helpers import *

def adjust_image (img, brightness = 0.5, contrast = 2.2):
    img = (img - img.min()) / (img.max() - img.min())  # normalize to range [0, 1]
    adjusted_img = (img - 0.5) * contrast + 0.5 + brightness
    adjusted_img = np.clip(adjusted_img, 0, 1)  # clip to avoid overflow

    return adjusted_img

def fov_across_days (obj, brightness = 0.5, contrast = 2.2):

    for i_day in range(len(obj.dat_subject.keys())):
        plt.figure(figsize = (10,5))
        plt.imshow(adjust_image(obj.track2p_obj.meanImg[i_day][10:-10, 30:-30], brightness=brightness, contrast=contrast), cmap = 'gray')
        plt.axis('off')

        folder_path = os.path.join(obj.save_path, 'FOV_across_days', obj.animal)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        plt.savefig(os.path.join(folder_path, f'day {i_day}.svg'))
        plt.savefig(os.path.join(folder_path, f'day {i_day}.png'))
        plt.show()

def polar_plots_across_days(obj, cell):
    '''
    plot orientation tuning curves as circular polar plots

    for each ROI, plot tuning curve for each day separately
    :param obj:
    :return:
    '''

    n_cells = obj.track2p_obj.track_ops.n_tracked
    days = list(obj.dat_subject.keys())

    # variables and dependencies for colour mapping
    plasma = plt.get_cmap('plasma')
    cNorm  = colors.Normalize(vmin=0, vmax=n_cells+5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)
    scalarMap.set_array([])

    fig, ax = plt.subplots(1, len(days), subplot_kw={'projection': 'polar'}, figsize=(7, 3))

    rmax = np.array([obj.dat_subject[day]['tuning_curves'].mean(axis = 0)[:,cell] / obj.dat_subject[day]['tuning_curves'].mean(axis = 0)[:,cell].sum() for day in days]).max()
    #max = np.array([suite2p_obj.dat_subject[day]['tuning_curves'][cell] for day in days]).max()
    for i_day, day in enumerate(days):

        # polar plots need to be plotted in radians
        # subtract the starting angle so all cells's starting preferred angle is at 0 degrees
        #theta- theta[0]

        # r = vector of responses for each direction(response vector)
        r = obj.dat_subject[day]['tuning_curves'].mean(axis = 0)[:,cell]
        r /= r.sum() # normalizing responses so they're between 0 and 1
        theta = np.deg2rad(obj.dat_subject[day]['orientations'])

        #to join the last point and first point
        idx = np.arange(r.shape[0] + 1)
        idx[-1] = 0

        # plotting
        ax[i_day].plot(theta, r, linewidth = 3.5, color=scalarMap.to_rgba(cell), alpha = 0.6)
        ax[i_day].plot(theta[idx], r[idx], linewidth = 3.5,color=scalarMap.to_rgba(cell), alpha = 0.6)
        ax[i_day].set_thetagrids([0, 90, 180, 270], y=0.1,
                                labels=['0', '\u03c0' + '/2', '\u03c0', '3' + '\u03c0' + '/2'],
                                fontsize=10)  # labels = ['0', '','\u03c0','']
        ax[i_day].set_rmax(rmax)
        ax[i_day].set_rlabel_position(45) # r is normalized response
        ax[i_day].tick_params(axis='y', labelsize=8)
        ax[i_day].set_rticks([])#set_rticks(np.round(np.linspace(0, rmax, 2),1))
        ax[i_day].grid(True)
        ax[i_day].set_title(f'Day {i_day}', fontsize = 12)

    fig.tight_layout(rect=[0, 0.06, 1, 0.80])
    plt.suptitle(f'Tuning Curves  \n ROI #{cell}', fontsize = 13)

    folder_path = os.path.join(obj.save_path, 'roi_polar_plots', obj.animal)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'polar_roi_{cell}.svg'))
    plt.savefig(os.path.join(folder_path, f'polar_roi_{cell}.png'))
    plt.show()

def polar_plots_across_days_rois(obj, cells_to_plot = [10,21,75,103]):
    '''
    plot orientation tuning curves as circular polar plots
    plot a handful of ROIs, each ROI in a separate plot (all days on same polar plot)

    :param obj:
    :return:
    '''

    n_cells = obj.track2p_obj.track_ops.n_tracked
    days = list(obj.dat_subject.keys())

    # variables and dependencies for colour mapping
    plasma = plt.get_cmap('plasma')
    cNorm  = colors.Normalize(vmin=0, vmax=len(days)-0.5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)
    scalarMap.set_array([])

    fig, ax = plt.subplots(1, len(cells_to_plot), subplot_kw={'projection': 'polar'}, figsize=(10, 3))

    for i_cell, cell in enumerate(cells_to_plot):

        rmax = np.array([obj.dat_subject[day]['tuning_curves'].mean(axis = 0)[:, cell] / obj.dat_subject[day]['tuning_curves'].mean(axis = 0)[:, cell].sum() for day in days]).max()
        #max = np.array([suite2p_obj.dat_subject[day]['tuning_curves'][cell] for day in days]).max()
        for i_day, day in enumerate(days):

            # polar plots need to be plotted in radians
            # subtract the starting angle so all cells's starting preferred angle is at 0 degrees
            #theta- theta[0]

            # r = vector of responses for each direction(response vector)
            r = obj.dat_subject[day]['tuning_curves'].mean(axis = 0)[:, cell] # mean across repeats
            r /= r.sum() # normalizing responses so they're between 0 and 1
            theta = np.deg2rad(obj.dat_subject[day]['orientations'])

            #to join the last point and first point
            idx = np.arange(r.shape[0] + 1)
            idx[-1] = 0

            # plotting > rad labels=['0', '\u03c0' + '/2', '\u03c0', '3' + '\u03c0' + '/2'],
            ax[i_cell].plot(theta, r, linewidth = 3.5, color=scalarMap.to_rgba(i_day), alpha = 0.35)
            ax[i_cell].plot(theta[idx], r[idx], linewidth = 3.5,color=scalarMap.to_rgba(i_day), alpha = 0.35)

        ax[i_cell].set_thetagrids([0, 90, 180, 270], y=0.13,
                                labels=['0', '\u03c0' + '/2', '\u03c0', '3' + '\u03c0' + '/2'],
                                fontsize=12)  # labels = ['0', '','\u03c0','']
        ax[i_cell].set_rmax(rmax)
        ax[i_cell].set_rlabel_position(45) # r is normalized response
        ax[i_cell].tick_params(axis='y', labelsize=8)
        #ax[i_day].set_rticks(np.round(np.linspace(0, rmax, 2),1))
        ax[i_cell].set_rticks([])
        ax[i_cell].grid(True)
        #ax[i_day].set_title(f'ROI #{cell}', fontsize = 10)

    #plt.tight_layout()#pad=0.9)
    fig.tight_layout(rect=[0, 0.03, 1.1, 0.95])
    plt.suptitle(f'Tuning Curves ', fontsize = 12)
    cbar = fig.colorbar(scalarMap, ax=ax, orientation='vertical', pad=0.1)
    cbar.set_label('Days', fontsize=14)
    cbar.set_ticks(np.linspace(0, len(days)-1, len(days)))
    cbar.set_ticklabels([int(num) for num in np.linspace(0, len(days)-1, len(days)) + 1 ])
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylim(0, len(days)-1)

    folder_path = os.path.join(obj.save_path, 'roi_polar_plot_roi_across', obj.animal)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'rois_days.svg'))
    plt.savefig(os.path.join(folder_path, f'rois_days.png'))
    plt.show()

def plot_raw_responses (object, animal, svd = False):
    '''
    Plot the trace of 20 random cells + the start of each TTL (grey dotted line)
    :param object:
    :param animal:
    :return:
    '''
    n_cells_to_plot = 8
    first_minutes = 5
    for day in object.dat_subject.keys():
        data_array = object.dat_subject[day]['zscored_tracked_fluorescence']#[:, 12*60*object.fps:13.5*60*object.fps]#[:, :first_minutes*60*object.fps]
        if svd:
            u, s, vt = np.linalg.svd(data_array)
            data_array = u[:, :2].reshape(-1, 2) @ np.diag(s[:2]) @ vt[:2].reshape(2, -1)
        ttl = object.dat_subject[day]['ttl_data']
        cells = [int(num) for num in np.linspace(0, data_array.shape[0]-1, n_cells_to_plot)]
        fig, ax = plt.subplots(figsize=(8.5, 7))
        colors = plt.cm.turbo(np.linspace(0, 1, data_array.shape[0]))
        colors = plt.cm.plasma(np.linspace(0, 0.9, data_array.shape[0]))
        fig.tight_layout()
        global_counter = 0
        for i_cell, cell in enumerate(cells):
            plt.plot( gaussian_filter1d(global_counter+data_array[cell,:]/data_array[cell,:].max(), sigma = 0.1), color=colors[cell], alpha = 0.8, linewidth = 2.2)
            # plt.plot(global_counter + data_array[cell, :] / data_array[cell, :].max(),
            #          color=colors[cell], alpha=0.8)
            global_counter += 0.5
            #[ax.axvline(x, c='grey', ls='--', linewidth=0.8, alpha=0.4) for x in np.unique(ttl)]
            #[ax.axvline(x, c= 'grey', ls = '--', linewidth = 0.8, alpha = 0.4) for x in np.unique(ttl) if x < (first_minutes*60*object.fps)]
        ax.set_yticks([])
        ax.set_xticks([])
        ax.axis('off')
        ax.set_title(f'{animal} ({day})', fontsize=15)
        #ax.set_xlabel('Time', fontsize=14)
        #ax.set_ylabel('Cell #', fontsize=14)
        plt.subplots_adjust(top=0.95)
        plt.xlim([3 * 60 * object.fps, 5 * 60 * object.fps])
        plt.tight_layout()
        #plt.xlim([12*60*object.fps,14*60*object.fps])
        plt.show()
        # if not os.path.exists(os.path.join(object.save_path, 'responses')):
        #     os.mkdir(os.path.join(object.save_path, 'responses'))
        #
        # plt.savefig(os.path.join(object.save_path,'responses', f'responses-{animal}_{day}_{session}.png' ))
        #
        # if not object.show_plots:
        #     plt.close('all')

def plot_rois_across_days (suite2pobj, data_object, wind_value, n_cells_to_plot = 8, brightness = 0.5, contrast = 2.2):
    '''
    Randomly sample n_cells_to_plot ROIs
    Plot ROI of cells across all recordings days

    :param wind_value: how much of each surround we wish to see
        - larger wind_value provides more zoomed out view of ROI
        - smaller wind_value provides more zoomed in view of ROI
    :param n_cells_to_plot: number of cells to randomly choose for plotting
    '''

    # randomly choose n_cells_to_plot cell ROIs to plot across days
    n_tracked_cells = data_object.t2p_match_mat_allday.shape[0]
    cells_to_plot = np.random.randint(0, n_tracked_cells, n_cells_to_plot)

    fig, ax = plt.subplots(n_cells_to_plot, len(data_object.track_ops.all_ds_path), figsize = (12,8))
    for i, cell_idx in enumerate(cells_to_plot):
        for i_day, path in enumerate(data_object.track_ops.all_ds_path):
            mean_img = data_object.all_ops[i_day]['meanImg']
            stat_t2p = data_object.all_stat_t2p[i_day]
            median_coord = stat_t2p[cell_idx]['med']

            # plot a short window around the ROI centroid
            img = mean_img[int(median_coord[0])-wind_value:int(median_coord[0])+wind_value, int(median_coord[1])-wind_value:int(median_coord[1])+wind_value]
            ax[i, i_day].imshow(adjust_image(img, brightness=brightness, contrast=contrast), cmap='gray')
            #ax[i, i_day].scatter(wind_value, wind_value)
            ax[i, i_day].axis('off')

            if i ==0:
                ax[i, i_day].set_title(path.split('\\')[4])

    plt.suptitle(path.split('\\')[3])
    fig.tight_layout()
    folder_path = os.path.join(suite2pobj.save_path, 'rois_across_days', suite2pobj.animal)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'day {i_day}.svg'))
    plt.savefig(os.path.join(folder_path, f'day {i_day}.png'))
    plt.show()

def rasters_across_days(object):
    #ndays x n orientations x n_Cells x n_timepoints
    arr = np.array([object.dat_subject[day]['param_matrix_whole'].mean(axis = 0) for day in object.dat_subject.keys()])
    cells_to_plot = np.array([7,35,46,60])# 12])
    ncols = len(object.dat_subject.keys())
    fig, ax = plt.subplots (len(cells_to_plot), ncols, figsize = (10,6))
    for i_cell, cell in enumerate(cells_to_plot):
        for i_day, day in enumerate(object.dat_subject.keys()):
            orientations = object.dat_subject[day]['orientations']
            #ax[i, i_day].imshow(arr[:, i, :], vmin = arr.min(), vmax = arr.max())
            ax[i_cell, i_day].imshow(arr[i_day, :, cell, :], aspect='auto', cmap='Greys', vmin = arr[:, :, i_cell, :].min(), vmax = arr[:, :, cell, :].max() )#, vmin=0, vmax=1.96)
            ax[i_cell, i_day].set_xlabel('Time since stim onset (s)', fontsize = 8)
            ax[i_cell, i_day].set_ylabel(f'Ori (deg)', fontsize = 8)
            ax[i_cell, i_day].set_title(f'ROI {i_cell} (day {i_day})')
            if cell == cells_to_plot[-1]:
                ax[i_cell, i_day].set_xticks((np.arange(arr.shape[-1])[::int(object.fps)]))
                ax[i_cell, i_day].set_xticklabels([int(np.round(x)) for x in (np.arange(-object.fps, arr.shape[-1]-object.fps)[::int(object.fps)])/int(object.fps)], fontsize =13)
            else:
                ax[i_cell, i_day].set_xticks([])
            ax[i_cell, i_day].set_yticks([])
            ax[i_cell, i_day].set_yticks(np.arange(arr.shape[1])[::4])
            ax[i_cell, i_day].set_yticklabels ([int(x) for x in orientations[::4]])
            ax[i_cell, i_day].axvline (int(object.fps), c = 'red', alpha = 0.4)
            ax[i_cell, i_day].axvline(2*int(object.fps), c='red', alpha = 0.4)
    plt.tight_layout()

    folder_path = os.path.join(object.save_path, 'rasters_across_days', object.animal)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'day {i_day}.svg'))
    plt.savefig(os.path.join(folder_path, f'day {i_day}.png'))
    plt.show()

def plot_response(object, cell_i = 0):

    '''
    plot response to all 12 orientations rows, for 4 rois (columns)
    :param object:
    :return:
    '''

    plasma = plt.get_cmap('plasma')
    cNorm = colors.Normalize(vmin=0, vmax=len(object.dat_subject.keys()) + 0.5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)
    scalarMap.set_array([])

    n_thetas = int(np.array([object.dat_subject[day]['n_thetas'] for day in object.dat_subject.keys()]).mean())

    fig, ax = plt.subplots (nrows = n_thetas, ncols = len(object.dat_subject.keys()), figsize = (10,5), sharey = True)

    for i_day, day in enumerate (object.dat_subject.keys()):

        thetas = object.dat_subject[day]['orientations']
        responses = object.dat_subject[day]['param_matrix_whole_zscore'][:,:,cell_i, :]

        for i_theta, theta in enumerate(thetas):

            ax[i_theta, i_day].plot (responses[:,i_theta].mean(axis = 0), c=scalarMap.to_rgba(i_day), linewidth = 2, alpha=0.8, zorder = 1)
            [ax[i_theta, i_day].plot(responses[i_response, i_theta], c='grey', alpha=0.4, zorder = 0) for i_response in range(responses.shape[0])]

            ax[i_theta, i_day].set_yticks([])
            ax[i_theta, i_day].axvline(object.fps, c = 'black', linestyle = 'dotted', alpha = 0.4)
            ax[i_theta, i_day].axvline(object.fps*2, c='black', linestyle = 'dashed', alpha = 0.6)

            if i_theta ==0:
                ax[i_theta, i_day].set_title(f'Day {i_day}')
            if i_day == 0:
                ax[i_theta, i_day].set_ylabel(f'{int(theta)}\u00b0', rotation = 0,  labelpad=20)
                ax[i_theta, i_day].yaxis.set_label_coords(-0.1, 0.4)  # Adjust the vertical position (0.4 moves it down)
            if theta == thetas[-1]:
                ax[i_theta, i_day].set_xticks((np.arange(responses.shape[-1])[::int(object.fps)]))
                ax[i_theta, i_day].set_xticklabels([int(np.round(x)) for x in (np.arange(-object.fps, responses.shape[-1]-object.fps)[::int(object.fps)])/int(object.fps)], fontsize =10)
                ax[i_theta, i_day].set_xlabel('Time since stim onset (s)')
            else:
                ax[i_theta, i_day].set_xticks([])

    plt.suptitle(f'ROI {cell_i}')

    folder_path = os.path.join(object.save_path, 'response_matrix', object.animal)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'day {i_day}.svg'))
    plt.savefig(os.path.join(folder_path, f'day {i_day}.png'))
    plt.show()


def plot_corr (obj, n = 1000, cumulative = False, across_days = False):
    # 0 > within day (0 with 0)
    # 1 > across day (0 with 1)
    # 2 > across day (0 with 2)
    # 3 > across day (0 with 3)
    # 4 > across day (0 with 4)

    if across_days:

        # just plot histogram
        # these are both of shape n_cells
        vec, vec_null = corr_vector(obj, null_distribution=False, across_days=across_days), corr_vector(obj, null_distribution=True, n=n, across_days=across_days)
        min_corr = np.min((vec.min(), vec_null.min()))

        if not cumulative:
            bin_edges = np.linspace(min(vec.min(), vec_null.min()), max(vec.max(), vec_null.max()), 10)

            #fig, ax = plt.subplots(1,2, sharey = True,figsize = (8,5))
            plt.figure(figsize=(5, 5))

            for i in range(2):
                plt.xlabel('Pearson (r) correlation', fontsize = 14)
                plt.ylabel('Cell count', fontsize = 14)
                plt.xlim([min_corr - 0.05, 1])

            plt.title('Correlation distribution (across days)', fontsize = 15)
            plt.hist(vec, color='lightsalmon', alpha = 0.8, bins = bin_edges,label = 'Aligned')
            plt.hist(vec_null,color='gray', alpha = 0.8,bins =bin_edges, label = 'Shuffled')
            plt.tick_params(axis='both', labelsize=12)
            plt.legend(fontsize=14)
            plt.tight_layout()
            plt.show()

        else :
            # interpretation: at Pearson = ____, the CDF shows 50% of the data is below or equal to it.
            # sort pearson correlation values in ascending order
            vec_sorted = np.sort(vec)
            vec_null_sorted = np.sort(vec_null)

            # get cumulative proportions > generate evenly spaced points between 0 and 1 > assigns a fraction (% of data) to each sorted value
            cum_prop = np.linspace(0, 1, len(vec_sorted))
            cum_prop_null = np.linspace(0, 1, len(vec_null_sorted))

            plt.figure()
            plt.plot(vec_sorted, cum_prop, c='lightsalmon', linewidth = 3, alpha=0.8, label='Aligned')
            plt.plot(vec_null_sorted, cum_prop_null, c='gray',  linewidth = 3,alpha=0.8, label='Shuffled')
            plt.tick_params(axis='both', labelsize=12)
            plt.legend(fontsize=14)
            plt.xlabel("Pearson (r) correlation", fontsize  = 14)
            plt.ylabel("Proportion of cells", fontsize  = 14)
            plt.title("Cumulative correlation distribution", fontsize  = 16)
            plt.show()

    else: # across days is false
        vec_null = corr_vector(obj, n = n, null_distribution=True, across_days = across_days)
        vec = corr_vector(obj, null_distribution=False, across_days=across_days)

        fig, ax = plt.subplots(1, len(obj.days), sharey = True, sharex = True, figsize = (10,8))

        for i in range(0, len(obj.days)):
            print(f'day {i}')
            stat, p_val = ks_2samp(vec[i], vec_null[i])
            print(f"KS test: D = {stat:.3f}, p = {p_val:.3e}")

        for i in range(0, len(obj.days)):
            ax[i-1].hist(vec[i], color='lightsalmon', alpha = 0.8, bins = np.linspace(vec.min(), 1, 10), label = 'Aligned')

        for i in range(0, len(obj.days)):
            ax[i].set_xlabel('Pearson (r) correlation')
            ax[i].set_ylabel('Cell count')
            ax[i].set_title(f'day 1 - day{i+1}')
            ax[i].hist(vec_null[i], color='gray',alpha = 0.6, bins = np.linspace(vec.min(), 1, 10), label = 'Shuffled')
        plt.suptitle('Correlation distribution')
        plt.tight_layout()
        plt.legend()
        plt.show()


def hist_osi_angle (obj):
    '''
    plots histogram distribution of OSIs/preferred orientation for each day
    :param obj:
    :param metric: 'preferred_orientation', 'OSI', 'preferred_direction', 'DSI'
    :return:
    '''
    for metric in ['OSI', 'preferred_orientation', 'DSI', 'preferred_direction']:
        if 'preferred' in metric:
            fig, ax = plt.subplots(1, len(obj.dat_subject.keys()), subplot_kw={'projection': 'polar'}, figsize=(2.5*len(obj.dat_subject.keys()), 4))
            c = 'indigo'
            rmax = np.array([np.histogram(np.deg2rad(obj.dat_subject[day][metric]),
                                          bins=np.linspace(int(np.round(np.deg2rad(obj.dat_subject[day][metric]).min())),
                                                           int(np.round(np.deg2rad(obj.dat_subject[day][metric]).max())),
                                                           25))[0].max() for day in obj.dat_subject.keys()]).max()
        else:
            fig, ax = plt.subplots(1, len(obj.dat_subject.keys()), figsize = (2.5*len(obj.dat_subject.keys()),3.5), sharey = True)
            c = 'mediumorchid'

        for i, day in enumerate(obj.dat_subject.keys()):
            if 'preferred' in metric:
                data = np.deg2rad(obj.dat_subject[day][metric])
                min_val, max_val = int(np.round(data.min())), int(np.round(data.max()))
                bins = np.linspace(min_val, max_val, 25)

                counts, bin_edges = np.histogram(data, bins=bins)

                # Calculate bin centers to use as the angular values (theta)
                theta = (bin_edges[:-1] + bin_edges[1:]) / 2
                radii = counts  # Histogram counts represent the radius values
                width = np.diff(bin_edges)  # Width of each bin for plotting

                # Plot the histogram on the polar axis

                if 'orientation' in metric:
                    mask = (theta >= np.deg2rad(-90)) & (theta <= np.deg2rad(90))
                    ax[i].bar(theta[mask], radii[mask], width=width[mask], bottom=0.0, color=c, alpha=0.7, label=f'Cell count, {day}')

                else:
                    ax[i].bar(theta, radii, width=width, bottom=0.0, color=c, alpha=0.7, label=f'Cell count, {day}')

                # Optionally, adjust the limits and titles
                ax[i].set_ylim(0, max(radii) + 10)  # Adjust the radial axis limit

                ax[i].set_thetagrids([0, 90, 180, 270], y=0.1,
                                         labels=['0', '\u03c0' + '/2', '\u03c0', '3' + '\u03c0' + '/2'],
                                         fontsize=10)  # labels = ['0', '','\u03c0','']
                ax[i].set_rmax(rmax)
                ax[i].set_rlabel_position(45)  # r is normalized response
                ax[i].tick_params(axis='y', labelsize=8)
                ax[i].set_rticks(np.round(np.linspace(0, rmax, 2),1))  # set_rticks(np.round(np.linspace(0, rmax, 2),1))
                ax[i].grid(True)


            else:
                data = obj.dat_subject[day][metric]
                min_val, max_val = int(np.round(data.min())), int(np.round(data.max()))
                bins = np.linspace(min_val, max_val, 25)
                ax[i].hist(obj.dat_subject[day][metric], color =  c, alpha = 0.7, label = f'Cell count, {day}', bins = bins)
                ax[i].set_xlim([min_val-(max_val/10), max_val+(max_val/10)])
                ax[i].set_xticks ([int(x) for x in np.linspace(min_val, max_val, 2)])

            ax[i].set_title (f'Day {i + 1}')

        plt.suptitle(metric)
        plt.tight_layout()

        folder_path = os.path.join(obj.save_path, 'hist_metrics', obj.animal)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        plt.savefig(os.path.join(folder_path, f'{metric}, day {i}.svg'))
        plt.savefig(os.path.join(folder_path, f'{metric}, day {i}.png'))
        plt.show()


##################################### need to adapt the following functions `

# generate_sf (animal, '20231106', np.load(os.path.join(r'G:\DriftScape\Data\EC_GECO_09\20231106\r2\r2_205_000\experiments\suite2p\plane0', 'stat.npy'), allow_pickle=True), np.load(os.path.join(r'G:\DriftScape\Data\EC_GECO_09\20231106\r2\r2_205_000\experiments\suite2p\plane0', 'ops.npy'), allow_pickle=True).item())
# generate_sf (animal, '20231107', np.load(os.path.join(r'G:\DriftScape\Data\EC_GECO_09\20231107\r1\r1_228_000\experiments\suite2p\plane0', 'stat.npy'), allow_pickle=True), np.load(os.path.join(r'G:\DriftScape\Data\EC_GECO_09\20231107\r1\r1_228_000\experiments\suite2p\plane0', 'ops.npy'), allow_pickle=True).item())
# generate_sf (animal, '20231108', np.load(os.path.join(r'G:\DriftScape\Data\EC_GECO_09\20231108\r1\r1_192_000\experiments\suite2p\plane0', 'stat.npy'), allow_pickle=True), np.load(os.path.join(r'G:\DriftScape\Data\EC_GECO_09\20231108\r1\r1_192_000\experiments\suite2p\plane0', 'ops.npy'), allow_pickle=True).item())
# generate_sf (animal, '20231109', np.load(os.path.join(r'G:\DriftScape\Data\EC_GECO_09\20231109\r1\r1_210_000\experiments\suite2p\plane0', 'stat.npy'), allow_pickle=True), np.load(os.path.join(r'G:\DriftScape\Data\EC_GECO_09\20231109\r1\r1_210_000\experiments\suite2p\plane0', 'ops.npy'), allow_pickle=True).item())

def tuning_curves_change (object):
    n_cells = object.track2p_obj.track_ops.n_tracked
    days = object.days
    ncols = 7
    nrows = int(np.ceil(n_cells / ncols))

    # variables and dependencies for colour mapping
    plasma = plt.get_cmap('plasma')
    cNorm  = colors.Normalize(vmin=0, vmax=object.grating_blocks)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)
    scalarMap.set_array([])

    for i_cell in range(n_cells):

        max_y = np.array([object.dat_subject[day]['mean_ordered_grat_responses'][:, :, i_cell].mean(axis=-1).max() for day in days]).max()

        responses_days = np.array([object.dat_subject[day]['tuning_curves'][i_cell] for day in days])
        orientations = np.linspace(0, 330, 12)

        fig, ax = plt.subplots(figsize=(7, 2))

        # plot mean tuning curve across each day
        [ax.plot(orientations, responses_days[i_day], c=scalarMap.to_rgba(i_day)) for i_day in range (len(days))]
        ax.set_xticks(np.arange(0,360,90))
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)

        plt.title('Tuning curves across days',fontsize = 15, y = 1)
        plt.xlabel('Orientation (deg)',fontsize = 11)
        plt.ylabel('Response', x = 0, ha = 'left',fontsize = 11)

        # Add colorbar
        cb = fig.colorbar(scalarMap, ax=ax, orientation='vertical')
        cb.set_label('Day', labelpad=8)
        cb.set_ticks(np.arange(len(days)))  # Adjust ticks to match the number of days

        #plt.colorbar()
        plt.tight_layout(pad=0.4)
        plt.show()



def plot_change_osi_angle(obj):
    n_cells = obj.track2p_obj.track_ops.n_tracked
    ncols = 13
    nrows = int(np.ceil(n_cells / ncols))

    max_mag = np.array([np.array([obj.dat_subject[day]['OSI'][cell] for day in obj.dat_subject.keys()]) for cell in
                        range(n_cells)]).max()

    for cell in [1,20,24, 28, 29, 33, 46, 48, 53]:

        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(3, 3))

        # polar plots need to be plotted in radians
        r = np.array([obj.dat_subject[day]['OSI'][cell] for day in obj.dat_subject.keys()])
        theta = np.deg2rad(np.array([obj.dat_subject[day]['pref_orientation'][cell] for day in obj.dat_subject.keys()]))

        # subtract the starting angle so all cells's starting preferred angle is at 0 degrees
        #theta- theta[0]
        ax.scatter(theta , r, s=60, c=np.arange(theta.shape[0]), cmap='plasma')
        ax.plot(theta, r, c='black', alpha=0.4)
        ax.set_thetamin(-180)
        ax.set_thetamax(180)
        ax.set_thetagrids([0, 90, -90, 180], y=0.05, labels=['0', '\u03c0' + '/2', '3' + '\u03c0' + '/2','\u03c0'], fontsize=10)  # labels = ['0', '','\u03c0','']
        ax.set_rlabel_position(45)
        ax.tick_params(axis='y', labelsize=10)
        ax.grid(True)

        #plt.tight_layout(pad=0.9)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.suptitle('Change in OSI/pref angle across days')
        plt.show()

#plot_change_osi_angle (suite2p_obj)







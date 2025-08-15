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

    for animal in obj.dat:

        for i_day, day in enumerate(obj.dat[animal].dat_subject.keys()):
            plt.figure(figsize = (10,5))
            plt.imshow(adjust_image(obj.dat[animal].track2p_obj.meanImg[i_day][10:-10, 30:-30], brightness=brightness, contrast=contrast), cmap = 'gray')
            plt.title(day)
            plt.axis('off')

            folder_path = os.path.join(obj.dat[animal].save_path, 'FOV_across_days', animal)
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            plt.savefig(os.path.join(folder_path, f'day {day}.svg'))
            plt.savefig(os.path.join(folder_path, f'day {day}.png'))
            #plt.show()
            plt.close()

def polar_plots(obj,animal):
    '''
    plot orientation tuning curves as circular polar plots
    :param obj:
    :return:
    '''

    days_recordings = [(day, subfile) for day in obj.dat[animal].dat_subject for subfile in obj.dat[animal].dat_subject[day] if 'grat' in subfile]
    days = [d[0] for d in days_recordings]
    n_cells = obj.dat[animal].track2p_obj.track_ops.n_tracked

    # variables and dependencies for colour mapping
    plasma = plt.get_cmap('plasma')
    cNorm  = colors.Normalize(vmin=0, vmax=n_cells+5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)
    scalarMap.set_array([])

    for i_day, day in enumerate(days):
        for cell in np.arange(n_cells):

            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(4, 4))

            rmax = np.array([obj.dat[animal].dat_subject[day]['grat']['tuning_curves'].mean(axis = 0)[:,cell] / obj.dat[animal].dat_subject[day]['grat']['tuning_curves'].mean(axis = 0)[:,cell].sum() for day in days]).max()
            #max = np.array([suite2p_obj.dat_subject[day]['tuning_curves'][cell] for day in days]).max()

            # polar plots need to be plotted in radians
            # subtract the starting angle so all cells's starting preferred angle is at 0 degrees
            #theta- theta[0]

            # r = vector of responses for each direction(response vector)
            r = obj.dat[animal].dat_subject[day]['grat']['tuning_curves'].mean(axis = 0)[:,cell]
            r /= r.sum() # normalizing responses so they're between 0 and 1
            theta = np.deg2rad(obj.dat[animal].dat_subject[day]['grat']['orientations'])

            #to join the last point and first point
            idx = np.arange(r.shape[0] + 1)
            idx[-1] = 0

            # plotting
            ax.plot(theta, r, linewidth = 3.5, color=scalarMap.to_rgba(cell), alpha = 0.6)
            ax.plot(theta[idx], r[idx], linewidth = 3.5,color=scalarMap.to_rgba(cell), alpha = 0.6)
            ax.set_thetagrids([0, 90, 180, 270], y=0.1,
                                    labels=['0', '\u03c0' + '/2', '\u03c0', '3' + '\u03c0' + '/2'],
                                    fontsize=10)  # labels = ['0', '','\u03c0','']
            ax.set_rmax(rmax)
            ax.set_rlabel_position(45) # r is normalized response
            ax.tick_params(axis='y', labelsize=8)
            ax.set_rticks([])#set_rticks(np.round(np.linspace(0, rmax, 2),1))
            ax.grid(True)
            ax.set_title(f'{calculate_animal_age(obj.animal_dobs[animal], day)}', fontsize = 12)

            fig.tight_layout(rect=[0, 0.06, 1, 0.80])
            plt.suptitle(f'ROI #{cell}', fontsize = 13)
            folder_path = os.path.join(obj.dat[animal].save_path, 'polar_plots_untracked',calculate_animal_age(obj.animal_dobs[animal], day), animal)
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            plt.savefig(os.path.join(folder_path, f'untracked tuning curve {cell}.svg'))
            plt.savefig(os.path.join(folder_path, f'untracked tuning curve {cell}.png'))
            #plt.show()
            plt.close()

def polar_plots_across_days(obj,animal):
    '''
    plot orientation tuning curves as circular polar plots

    for each ROI, plot tuning curve for each day separately
    :param obj:
    :return:
    '''

    days_recordings = [(day, subfile) for day in obj.dat[animal].dat_subject for subfile in obj.dat[animal].dat_subject[day] if 'grat' in subfile]
    days = [d[0] for d in days_recordings]
    n_cells = obj.dat[animal].track2p_obj.track_ops.n_tracked

    # variables and dependencies for colour mapping
    plasma = plt.get_cmap('plasma')
    cNorm  = colors.Normalize(vmin=0, vmax=n_cells+5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)
    scalarMap.set_array([])

    for cell in np.arange(n_cells):

        fig, ax = plt.subplots(1, len(days), subplot_kw={'projection': 'polar'}, figsize=(7, 3))

        rmax = np.array([obj.dat[animal].dat_subject[day]['grat']['tuning_curves'].mean(axis = 0)[:,cell] / obj.dat[animal].dat_subject[day]['grat']['tuning_curves'].mean(axis = 0)[:,cell].sum() for day in days]).max()
        #max = np.array([suite2p_obj.dat_subject[day]['tuning_curves'][cell] for day in days]).max()
        for i_day, day in enumerate(days):

            # polar plots need to be plotted in radians
            # subtract the starting angle so all cells's starting preferred angle is at 0 degrees
            #theta- theta[0]

            # r = vector of responses for each direction(response vector)
            r = obj.dat[animal].dat_subject[day]['grat']['tuning_curves'].mean(axis = 0)[:,cell]
            r /= r.sum() # normalizing responses so they're between 0 and 1
            theta = np.deg2rad(obj.dat[animal].dat_subject[day]['grat']['orientations'])

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
            ax[i_day].set_title(f'{calculate_animal_age(obj.animal_dobs[animal], day)}', fontsize = 12)

        fig.tight_layout(rect=[0, 0.06, 1, 0.80])
        plt.suptitle(f'ROI #{cell}', fontsize = 13)

        folder_path = os.path.join(obj.dat[animal].save_path, 'roi_polar_plots', animal)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        plt.savefig(os.path.join(folder_path, f'polar_roi_{cell}.svg'))
        plt.savefig(os.path.join(folder_path, f'polar_roi_{cell}.png'))
        #plt.show()
        plt.close()

def polar_plots_across_days_rois(obj, animal):
    '''
    plot orientation tuning curves as circular polar plots
    plot a handful of ROIs, each ROI in a separate plot (all days on same polar plot)

    :param obj:
    :return:
    '''

    days_recordings = [(day, subfile) for day in obj.dat[animal].dat_subject for subfile in obj.dat[animal].dat_subject[day] if 'grat' in subfile]
    days = [d[0] for d in days_recordings]
    n_cells = obj.dat[animal].track2p_obj.track_ops.n_tracked

    # variables and dependencies for colour mapping
    plasma = plt.get_cmap('plasma')
    cNorm  = colors.Normalize(vmin=0, vmax=len(days)-0.5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)
    scalarMap.set_array([])

    for cell in range(n_cells):

        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(5, 5))

        rmax = np.array([obj.dat[animal].dat_subject[day]['grat']['tuning_curves'].mean(axis = 0)[:, cell] / obj.dat[animal].dat_subject[day]['grat']['tuning_curves'].mean(axis = 0)[:, cell].sum() for day in days]).max()
        #max = np.array([suite2p_obj.dat_subject[day]['tuning_curves'][cell] for day in days]).max()
        for i_day, day in enumerate(days):

            # polar plots need to be plotted in radians
            # subtract the starting angle so all cells's starting preferred angle is at 0 degrees
            #theta- theta[0]

            # r = vector of responses for each direction(response vector)
            r = obj.dat[animal].dat_subject[day]['grat']['tuning_curves'].mean(axis = 0)[:, cell] # mean across repeats
            r /= r.sum() # normalizing responses so they're between 0 and 1
            theta = np.deg2rad(obj.dat[animal].dat_subject[day]['grat']['orientations'])

            #to join the last point and first point
            idx = np.arange(r.shape[0] + 1)
            idx[-1] = 0

            # plotting > rad labels=['0', '\u03c0' + '/2', '\u03c0', '3' + '\u03c0' + '/2'],
            ax.plot(theta, r, linewidth = 3.5, color=scalarMap.to_rgba(i_day), alpha = 0.35)
            ax.plot(theta[idx], r[idx], linewidth = 3.5,color=scalarMap.to_rgba(i_day), alpha = 0.35)

            ax.set_thetagrids([0, 90, 180, 270], y=0.03,
                                    labels=['0', '\u03c0' + '/2', '\u03c0', '3' + '\u03c0' + '/2'],
                                    fontsize=12)  # labels = ['0', '','\u03c0','']
            ax.set_rmax(rmax)
            ax.set_rlabel_position(45) # r is normalized response
            ax.tick_params(axis='y', labelsize=8)
            #ax[i_day].set_rticks(np.round(np.linspace(0, rmax, 2),1))
            ax.set_rticks([])
            ax.grid(True)
            #ax.set_title(f'{calculate_animal_age(obj.animal_dobs[animal], day)}', fontsize = 10, pad = 15)

        #plt.tight_layout()#pad=0.9)
        #fig.tight_layout(rect=[0, 0.03, 1.1, 0.95])
        plt.suptitle(f'Cell #{cell}', fontsize = 12)
        cbar = fig.colorbar(scalarMap, ax=ax, orientation='vertical', pad=0.1)
        cbar.set_label('Days', fontsize=14)
        cbar.set_ticks(np.linspace(0, len(days)-1, len(days)))
        # cbar.set_ticklabels([int(num) for num in np.linspace(0, len(days)-1, len(days)) + 1 ])
        cbar.set_ticklabels([calculate_animal_age(obj.animal_dobs[animal], day) for day in days])
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.set_ylim(0, len(days)-1)

        folder_path = os.path.join(obj.dat[animal].save_path, 'roi_polar_plot_roi_across', animal)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        plt.savefig(os.path.join(folder_path, f'{animal}, roi {cell}.svg'))
        plt.savefig(os.path.join(folder_path, f'{animal}, roi {cell}.png'))
        #plt.show()
        plt.close()

def plot_raw_responses (object, animal, svd = False):
    '''
    Plot the trace of 20 random cells + the start of each TTL (grey dotted line)
    :param object:
    :param animal:
    :return:
    '''

    n_cells_to_plot = 8
    first_minutes = 5
    for day in object.dat[animal].dat_subject.keys():
        for recording in object.dat[animal].dat_subject[day]:
            data_array = object.dat[animal].dat_subject[day][recording]['zscored_traces']#[:, 12*60*object.fps:13.5*60*object.fps]#[:, :first_minutes*60*object.fps]
            if svd:
                u, s, vt = np.linalg.svd(data_array)
                data_array = u[:, :2].reshape(-1, 2) @ np.diag(s[:2]) @ vt[:2].reshape(2, -1)
            if 'grat' in recording:
                ttl = object.dat[animal].dat_subject[day][recording]['ttl_data']
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
            ax.set_title(f'{animal} ({day}, {recording})', fontsize=15)
            #ax.set_xlabel('Time', fontsize=14)
            #ax.set_ylabel('Cell #', fontsize=14)
            plt.subplots_adjust(top=0.95)
            plt.xlim([3 * 60 * object.fps, 5 * 60 * object.fps])
            plt.tight_layout()
            #plt.xlim([12*60*object.fps,14*60*object.fps])
            plt.show()
            if not os.path.exists(os.path.join(object.dat[animal].save_path, 'responses')):
                os.mkdir(os.path.join(object.dat[animal].save_path, 'responses'))

            plt.savefig(os.path.join(object.dat[animal].save_path,'responses', f'responses-{animal}_{day}_{recording}.png' ))


def plot_rois_across_days (obj, animal, wind_value, n_cells_to_plot = 8, brightness = 0.5, contrast = 2.2):
    '''
    Randomly sample n_cells_to_plot ROIs
    Plot ROI of cells across all recordings days

    :param wind_value: how much of each surround we wish to see
        - larger wind_value provides more zoomed out view of ROI
        - smaller wind_value provides more zoomed in view of ROI
    :param n_cells_to_plot: number of cells to randomly choose for plotting
    '''
    suite2pobj = obj.dat[animal]
    data_object = suite2pobj.track2p_obj

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

def rasters_across_days(object, animal):
    '''
    Plot each cell's response on each recording day
    '''

    days_recordings = [(day, subfile) for day in object.dat[animal].dat_subject for subfile in object.dat[animal].dat_subject[day] if 'grat' in subfile]
    days = [d[0] for d in days_recordings]
    folder_path = os.path.join(object.dat[animal].save_path, 'rasters_across_days', animal)

    #ndays x n orientations x n_Cells x n_timepoints
    arr = np.array([object.dat[animal].dat_subject[day]['grat']['param_matrix_whole'].mean(axis = 0) for day in days])

    for i_cell in range(object.dat[animal].track2p_obj.track_ops.n_tracked):
        fig, ax = plt.subplots(nrows=len(days), ncols=1, figsize=(5, 2 * len(days)))
        for i_day, day in enumerate(days):
            orientations = object.dat[animal].dat_subject[day]['grat']['orientations']
            #ax[i, i_day].imshow(arr[:, i, :], vmin = arr.min(), vmax = arr.max())
            #ax[i_day].imshow(arr[i_day, :, i_cell, :], aspect='auto', cmap='Greys', vmin = arr[:, :, i_cell, :].min(), vmax = arr[:, :, i_cell, :].max() )#, vmin=0, vmax=1.96)
            ax[i_day].imshow(arr[i_day, :, i_cell, :], aspect='auto', cmap='Greys', vmin=arr[i_day, :, i_cell, :].min(),
                             vmax=arr[i_day, :, i_cell, :].max())  # , vmin=0, vmax=1.96)

            ax[i_day].set_xlabel('Time since stim onset (s)', fontsize = 8)
            ax[i_day].set_ylabel(f'Ori (deg)', fontsize = 8)
            ax[i_day].set_title(f'{calculate_animal_age(object.animal_dobs[animal], day)}')
            ax[i_day].set_xticks((np.arange(arr.shape[-1])[::int(object.fps)]))
            ax[i_day].set_xticklabels([int(np.round(x)) for x in (np.arange(-object.fps, arr.shape[-1]-object.fps)[::int(object.fps)])/int(object.fps)], fontsize =10)
            ax[i_day].set_yticks([])
            ax[i_day].set_yticks(np.arange(arr.shape[1])[::4])
            ax[i_day].set_yticklabels ([int(x) for x in orientations[::4]])
            ax[i_day].axvline (int(object.fps), c = 'red', alpha = 0.4)
            ax[i_day].axvline(2*int(object.fps), c='red', alpha = 0.4)
        plt.suptitle(f'ROI {i_cell}')

        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        plt.savefig(os.path.join(folder_path, f'{animal}, day {i_day}, cell {i_cell}.svg'))
        plt.savefig(os.path.join(folder_path, f'{animal}, day {i_day}, cell {i_cell}.png'))
        plt.tight_layout()
        plt.close()



def plot_response(object, animal):

    '''
    plot fluorescence trace response to all orientations (rows), for each recordng (grating) day (columns)
    :param object:
    :return:
    '''

    days_recordings = [(day, subfile) for day in object.dat[animal].dat_subject for subfile in object.dat[animal].dat_subject[day] if 'grat' in subfile]
    days = [d[0] for d in days_recordings]
    folder_path = os.path.join(object.dat[animal].save_path, 'response_matrix', animal)

    plasma = plt.get_cmap('plasma')
    cNorm = colors.Normalize(vmin=0, vmax=len(days) + 0.5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)
    scalarMap.set_array([])

    n_thetas = object.ntheta

    for cell_i in range(object.dat[animal].track2p_obj.track_ops.n_tracked):

        fig, ax = plt.subplots (nrows = n_thetas+1, ncols = len(days), figsize = (6*len(days),5), sharey = True)

        # add axis for correlation figure
        right_ax = fig.add_axes([0.92, 0.3, 0.05, 0.4])  # [left, bottom, width, height] in figure coords
        upper_diag = np.array([np.diag(object.dat[animal].corr_matrix[..., i], k=1) for i in range(object.dat[animal].corr_matrix.shape[-1])])
        right_ax.plot(upper_diag[cell_i], marker='o')
        right_ax.set_title('Corr', fontsize=8)
        right_ax.tick_params(axis='both', labelsize=6)
        right_ax.set_ylim(-1, 1)  # optional: set fixed y-axis range for correlation

        for i_day, (day, subfile) in enumerate(days_recordings):

            thetas = object.dat[animal].dat_subject[day][subfile]['orientations']
            responses = object.dat[animal].dat_subject[day][subfile]['zscored'][:,:,cell_i, :]

            for i_theta, theta in enumerate(thetas):

                ax[i_theta, i_day].plot (responses[:,i_theta].mean(axis = 0), c=scalarMap.to_rgba(i_day), linewidth = 2, alpha=0.8, zorder = 1)
                [ax[i_theta, i_day].plot(responses[i_response, i_theta], c='grey', alpha=0.4, zorder = 0) for i_response in range(responses.shape[0])]

                ax[i_theta, i_day].set_yticks([])
                ax[i_theta, i_day].axvline(object.fps, c = 'black', linestyle = 'dotted', alpha = 0.4)
                ax[i_theta, i_day].axvline(object.fps*2, c='black', linestyle = 'dashed', alpha = 0.6)

                if i_theta ==0 : # set title if in top row

                    if object.dat[animal].dat_subject[day][subfile]['thresholded_cells'][cell_i]:
                        ax[i_theta, i_day].set_title(f'{calculate_animal_age(object.animal_dobs[animal], day)}', color='green', fontweight='bold')
                    else:
                        ax[i_theta, i_day].set_title(f'{calculate_animal_age(object.animal_dobs[animal], day)}', color='r')

                if i_day == 0:
                    ax[i_theta, i_day].set_ylabel(f'{int(theta)}\u00b0', rotation = 0,  labelpad=20)
                    ax[i_theta, i_day].yaxis.set_label_coords(-0.1, 0.4)  # Adjust the vertical position (0.4 moves it down)
                if theta == thetas[-1]:
                    ax[i_theta, i_day].set_xticks((np.arange(responses.shape[-1])[::int(object.fps)]))
                    ax[i_theta, i_day].set_xticklabels([int(np.round(x)) for x in (np.arange(-object.fps, responses.shape[-1]-object.fps)[::int(object.fps)])/int(object.fps)], fontsize =10)
                    ax[i_theta, i_day].set_xlabel('Time since stim onset (s)')
                else:
                    ax[i_theta, i_day].set_xticks([])

        if np.array([object.dat[animal].dat_subject[day][subfile]['thresholded_cells'][cell_i] for (day, subfile) in days_recordings]).any():
            plt.suptitle(f'Cell {cell_i}', fontsize=16, color='green', fontweight='bold')
        else:
            plt.suptitle(f'Cell {cell_i}', fontsize=16, color='r')

        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        plt.savefig(os.path.join(folder_path, f'{animal}, cell {cell_i}.svg'))
        plt.savefig(os.path.join(folder_path, f'{animal}, cell {cell_i}.png'))
        #plt.show()
        plt.close()


def plot_corr (obj, animal, thresholded_cells = 0, n = 1000, cumulative = False, across_days = False):
    '''
    :param obj:
    :param animal (str):
    :param thresholded_cells: either
        - 0 (take all cells for the analysis)
        - 1 (cells that pass treshold at least once)
        - 2 (cells that pass threshold on all days)
        - 3 cells that pass thresholds on both days we are considering
        - 4 cells that pass threshold on only 1 of the two days we are considering
    :param null_distribution (bool):
    :param n (int):
    :param across_days: True if we want to get 1 correlation value across all days(average the matrix). False if we want to compare separate days
    :return: correlation matrix, shape (n_days, n_days, n_cells)

    NB: deconvolved must be true
    '''
    # 0 > within day (0 with 0)
    # 1 > across day (0 with 1)
    # 2 > across day (0 with 2)
    # 3 > across day (0 with 3)
    # 4 > across day (0 with 4)

    dat_object = obj.dat[animal].dat_subject

    days_recordings = [(day, subfile) for day in dat_object for subfile in dat_object[day] if 'grat' in subfile]
    days = [d[0] for d in days_recordings]

    if across_days:

        # just plot histogram
        # these are both of shape n_cells
        vec, vec_null = corr_vector(obj, animal,thresholded_cells = thresholded_cells, null_distribution=False, across_days=across_days), corr_vector(obj, animal,thresholded_cells=thresholded_cells, null_distribution=True, n=n, across_days=across_days)
        min_corr = np.min((vec.min(), vec_null.min()))

        if not cumulative:
            bin_edges = np.linspace(min(vec.min(), vec_null.min()), max(vec.max(), vec_null.max()), 10)

            #fig, ax = plt.subplots(1,2, sharey = True,figsize = (8,5))
            plt.figure(figsize=(5, 5))

            for i in range(2):
                plt.xlabel('Pearson (r) correlation', fontsize = 14)
                plt.ylabel('Cell count', fontsize = 14)
                plt.xlim([min_corr - 0.05, 1])

            plt.title(f'Correlation distribution (across days) {animal}', fontsize = 15)
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
            plt.title(f"Cumulative correlation distribution {animal}", fontsize  = 16)
            plt.show()

    else: # across days is false
        # each of these arrays is shape (n_days x n_days x n_cells)
        # on-diagonal ([0,0], [1,1]...) > within day comparisons (split-half)
        # off-diagonal ([0,1], [0,2], [1,0]...) > across-week comparison
        vec_null = corr_vector(obj, animal, thresholded_cells = thresholded_cells, n = n, null_distribution=True, across_days = across_days)
        vec = corr_vector(obj, animal, thresholded_cells = thresholded_cells, null_distribution=False, across_days=across_days)

        #vec_thesholded = np.zeros_like(vec)

        # vec is of shape ((n_days, n_days, n_cells))

        fig, ax = plt.subplots(len(days), len(days), sharey = True, sharex = True, figsize = (10,8))

        for i in range(len(days)):
            for k in range(len(days)):


                if thresholded_cells == 0:
                    vec_thresholded, vec_null_thresholded = vec, vec_null
                elif thresholded_cells == 3: # take cells that are responsive on both days
                    thresholded_cells_idx = dat_object[days[i]]['grat']['thresholded_cells'] & dat_object[days[k]]['grat']['thresholded_cells']
                    vec_thresholded, vec_null_thresholded = vec[..., thresholded_cells_idx], vec_null[..., thresholded_cells_idx]
                elif thresholded_cells == 4:
                    # take cells that are responsive on at least one day
                    thresholded_cells_idx = dat_object[days[i]]['grat']['thresholded_cells'] | dat_object[days[k]]['grat']['thresholded_cells']
                    vec_thresholded, vec_null_thresholded = vec[..., thresholded_cells_idx], vec_null[..., thresholded_cells_idx]

                #vec_thesholded[i,k] = vec[i, k, thresholded_cells_idx]

                if obj.tracked_cells:
                    ax[i,k].hist(vec_thresholded[i, k], color='lightsalmon', alpha=0.8, bins=np.linspace(vec_thresholded.min(), 1, 10), label='Aligned')
                    ax[i,k].hist(vec_null_thresholded[i, k], color='gray', alpha=0.6, bins=np.linspace(vec_thresholded.min(), 1, 10), label='Shuffled')
                    ax[i,k].set_xlabel('Pearson (r) correlation')
                    ax[i,k].set_ylabel('Cell count')
                    ax[i,k].set_title(f'{calculate_animal_age(obj.animal_dobs[animal], days[i])} x {calculate_animal_age(obj.animal_dobs[animal], days[k])}')
                elif i==k: # if not tracking cells over time, only plot the diagonal (within-day corr)
                    ax[i,k].hist(vec_thresholded[i, k], color='lightsalmon', alpha=0.8, bins=np.linspace(vec_thresholded.min(), 1, 10), label='Aligned')
                    ax[i,k].hist(vec_null_thresholded[i, k], color='gray', alpha=0.6, bins=np.linspace(vec_thresholded.min(), 1, 10), label='Shuffled')
                    ax[i,k].set_xlabel('Pearson (r) correlation')
                    ax[i,k].set_ylabel('Cell count')
                    ax[i,k].set_title(f'{calculate_animal_age(obj.animal_dobs[animal], days[i])} x {calculate_animal_age(obj.animal_dobs[animal], days[k])}')
                else:
                    ax[i,k].set_axis_off

        plt.suptitle(f'Correlation distribution {animal}')
        plt.tight_layout()
        plt.legend()
        plt.show()

    if thresholded_cells==0:
        t = 'all_cells'
    elif thresholded_cells==1:
        t = 'thresholded'
    elif thresholded_cells==2:
        t = 'thresholded_all_days'
    else:
        t = ''

    folder_path = os.path.join(obj.dat[animal].save_path, f'correlation distribution{", tracked_cells"*obj.tracked_cells}', animal)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plt.savefig(os.path.join(folder_path, f'corr dist{"across days"*across_days}{", cumulative"*cumulative}{", tracked_cells"*obj.tracked_cells} {t}.svg'))
    plt.savefig(os.path.join(folder_path, f'corr dist{"across days"*across_days}{", cumulative"*cumulative}{", tracked_cells"*obj.tracked_cells} {t}.png'))

    return vec

def hist_osi_angle (obj):
    '''
    plots histogram distribution of OSIs/preferred orientation for each day
    :param obj:
    :param metric: 'preferred_orientation', 'OSI', 'preferred_direction', 'DSI'
    :return:
    '''
    for animal in obj.dat:

        subject_data = obj.dat[animal].dat_subject
        days_recordings = [(day, subfile) for day in subject_data for subfile in subject_data[day] if 'grat' in subfile]

        for metric in ['OSI', 'preferred_orientation', 'DSI', 'preferred_direction']:

            if 'preferred' in metric: # +1 is temporary, to change when we get more recordings
                fig, ax = plt.subplots(1, len(days_recordings), subplot_kw={'projection': 'polar'}, figsize=(4*len(days_recordings), 4))
                c = 'indigo'

            else:
                fig, ax = plt.subplots(1, len(days_recordings), figsize = (4*len(days_recordings),3.5), sharey = True)
                c = 'mediumorchid'

            for i, (day, subfile) in enumerate(days_recordings):

                rmax = np.array([np.histogram(np.deg2rad(subject_data[day][subfile][metric]),
                                              bins=np.linspace(int(np.round(np.deg2rad(subject_data[day][subfile][metric]).min())),
                                                               int(np.round(np.deg2rad(subject_data[day][subfile][metric]).max())),
                                                               25))[0].max() for day in [tup[0] for tup in days_recordings]]).max()

                if 'preferred' in metric:
                    data = np.deg2rad(subject_data[day][subfile][metric])
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
                    data = subject_data[day][subfile][metric]
                    min_val, max_val = int(np.round(data.min())), int(np.round(data.max()))
                    bins = np.linspace(min_val, max_val, 25)
                    ax[i].hist(subject_data[day][subfile][metric], color =  c, alpha = 0.7, label = f'Cell count, {day}', bins = bins)
                    ax[i].set_xlim([min_val-(max_val/10), max_val+(max_val/10)])
                    ax[i].set_xticks ([int(x) for x in np.linspace(min_val, max_val, 2)])

                ax[i].set_title (calculate_animal_age(obj.animal_dobs[animal], day))

            plt.suptitle(f'{animal}, {metric}')
            plt.tight_layout()

            folder_path = os.path.join(obj.dat[animal].save_path, 'hist_metrics', animal)
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            plt.savefig(os.path.join(folder_path, f'{metric} {animal}.svg'))
            plt.savefig(os.path.join(folder_path, f'{metric} {animal}.png'))
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







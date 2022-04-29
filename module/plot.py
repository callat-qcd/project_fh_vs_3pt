import matplotlib.pyplot as plt
import numpy as np
import gvar as gv
from matplotlib.ticker import FormatStrFormatter

plt.rcParams.update({"text.usetex": True})

grey = "#808080" # nstates = 1
red = "#FF6F6F" # nstates = 2
peach = "#FF9E6F" # nstates = 3
orange = "#FFBC6F" # nstates = 4
sunkist = "#FFDF6F"
yellow = "#FFEE6F"
lime = "#CBF169"
green = "#5CD25C" # nstates = 5
turquoise = "#4AAB89"
blue = "#508EAD" # nstates = 6
grape = "#635BB1"
violet = "#7C5AB8" # nstates = 7
fuschia = "#C3559F"

color_list = [red, peach, orange, green, blue, violet, grey] # for stability plots
psychedelic_list = [grey, fuschia, violet, grape, blue, turquoise, green, lime, yellow, sunkist, orange, peach, red] # for psychedelic moose

marker_list = ['8', 'h', 'D', 's', 'o', '*']

fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
fig_size  = (fig_width, fig_width / gr)
plt_axes  = [0.12,0.15,0.87,0.83]
gridspec_tmax = {'height_ratios': [3, 1], 'left': 0.12, 'right': 0.99, 'bottom': 0.15, 'top': 0.98}
gridspec_tmin = {'height_ratios': [4, 1, 1], 'left': 0.12, 'right': 0.99, 'bottom': 0.15, 'top': 0.98}
gridspec_tmin_div = {'height_ratios': [3, 1], 'left': 0.12, 'right': 0.99, 'bottom': 0.16, 'top': 0.98}
gridspec_prior_width = {'height_ratios': [3, 1], 'left': 0.12, 'right': 0.99, 'bottom': 0.15, 'top': 0.98}

plt.rcParams['figure.figsize'] = fig_size

fs_p = {"fontsize": 18} # font size of text, label, ticks
gA_fs_p = {"fontsize": 20} # font size of ga label
ls_p = {"labelsize": 18}
errorp = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1} # circle
errorb = {"markersize": 5, "linestyle": "none", "capsize": 3, "elinewidth": 1} # solid circle

linspace_num = 100 # bigger, more smoothy
plot_space = 0.05 # tau_c_plot, m_z_plot
omega_imp_a09 = 0.08730 # converse lattice to fm

gA_ylim = [1.01, 1.349]

# labels
ra_label = r"$R_{A_3}(t_{\rm sep}, \tau)$"
tau_label = r"$(\tau - t_{\rm sep}/2) (a_{09})$"
fm_tau_label = r"$(\tau - t_{\rm sep}/2) ({\rm fm})$"
fha_label = r"${\rm FH}_{A_3}(t_{\rm sep}, \tau_c)$"
tsep_label = r'$t_{\textrm{sep}} (a_{09})$'
tsep_fm_label = r'$t_{\textrm{sep}} ({\rm fm})$'
meff_label = r"$m_{\textrm{eff}}$"
zeff_label = r"$10^3 z_{\textrm{eff}}$"
ga_label = r"$\mathring{g}_A$"
e0_label = r"$E_{0}$"
q_label = r"$Q$"
w_label = r"$w$"

tmax_label = r"$t_{\rm sep}^{\rm max}$"
tmin_label = r"$t_{\rm sep}^{\rm min}$"

def legend_without_duplicate_labels(ax, loc, ncol):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), loc=loc, ncol=ncol, handletextpad=0, columnspacing=0, **fs_p)

def moose_plot(data_avg_dict_completed, fit_result, fitter, pt3_data_range, tau_c_2, divide_n):
    gA = {}
    gA_err = {}
    gA_fit = {}
    gA_fit_err = {}

    gA_tsep = []
    gA_tau = []
    gV_tsep = []
    gV_tau = []

    tau_c_1 = 1 # tau_c for tsep in range(2, divide_n)
    # tau_c_2 : tau_c for tsep in range(divide_n, 15)

    for i in pt3_data_range:
        gA['tsep_'+str(i)] = []
        gA_err['tsep_'+str(i)] = []

        for j in range(tau_c_1, i - tau_c_1 + 1):
            temp = ( data_avg_dict_completed['pt3_A3_tsep_'+str(i)][j] + data_avg_dict_completed['pt3_A3_tsep_'+str(i)][i-j] ) / (2 * data_avg_dict_completed['pt2_tsep_'+str(i)])

            gA['tsep_'+str(i)].append(temp.mean)
            gA_err['tsep_'+str(i)].append(temp.sdev)

        for j in np.linspace(tau_c_1, i - tau_c_1, linspace_num):
            gA_tsep.append(i)
            gA_tau.append(j)
            gV_tsep.append(i)
            gV_tau.append(j)

    pt2_fitter = fitter.pt2_fit_function(np.array([t for t in pt3_data_range]), fit_result.p)['pt2']

    pt3_gA_fitter = fitter.pt3_fit_function(np.array(gA_tsep), np.array(gV_tsep), np.array(gA_tau), np.array(gV_tau), fit_result.p)['pt3_A3']

    index = 0

    for i in pt3_data_range:
        gA_fit['tsep_'+str(i)] = []
        gA_fit_err['tsep_'+str(i)] = []

        for j in range(linspace_num):
            index = int((i-pt3_data_range[0])*linspace_num + j)

            temp = pt3_gA_fitter[index] / pt2_fitter[i - pt3_data_range[0]]
            gA_fit['tsep_'+str(i)].append(temp.mean)
            gA_fit_err['tsep_'+str(i)].append(temp.sdev)

    ########### here start a fig ###########
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    for i in pt3_data_range:
        idx = i - 2 # psychedelic_list start from tsep = 2
        color = psychedelic_list[idx]

        x_errorbar = np.arange(tau_c_1 - i/2, i/2 - tau_c_1 + 1) * omega_imp_a09
        x_fill = np.linspace(tau_c_1 - i/2, i/2 - tau_c_1, linspace_num) * omega_imp_a09

        if 2 < i < divide_n:
            ax.errorbar(x_errorbar, np.array(gA['tsep_' + str(i)]), yerr=np.array(gA_err['tsep_' + str(i)]), marker='o', color=color, **errorp)

            gA_fit_y1 = np.array(gA_fit['tsep_'+str(i)]) + np.array(gA_fit_err['tsep_'+str(i)])
            gA_fit_y2 = np.array(gA_fit['tsep_'+str(i)]) - np.array(gA_fit_err['tsep_'+str(i)])
            ax.fill_between(x_fill, gA_fit_y1, gA_fit_y2, color=color, alpha=0.3)

        elif i == 2:
            ax.errorbar(np.array([0]), np.array(gA['tsep_' + str(i)]), yerr=np.array(gA_err['tsep_' + str(i)]), marker='s', color=grey, **errorp)

        elif i >= divide_n:
            for idx in range( len(x_errorbar) ):
                if (tau_c_2 - i/2) * omega_imp_a09 <= x_errorbar[idx] <= (i/2 - tau_c_2) * omega_imp_a09:
                    ax.errorbar(np.array(x_errorbar[idx]), np.array(gA['tsep_' + str(i)][idx]), yerr=np.array(gA_err['tsep_' + str(i)][idx]), marker='o', color=color, **errorp)

                else:
                    ax.errorbar(np.array(x_errorbar[idx]), np.array(gA['tsep_' + str(i)][idx]), yerr=np.array(gA_err['tsep_' + str(i)][idx]), marker='s', color=grey, **errorp)

            gA_fit_y1 = np.array(gA_fit['tsep_'+str(i)]) + np.array(gA_fit_err['tsep_'+str(i)])
            gA_fit_y2 = np.array(gA_fit['tsep_'+str(i)]) - np.array(gA_fit_err['tsep_'+str(i)])

            x_fill_in = []
            x_fill_out_1 = []
            x_fill_out_2 = []

            fit_y1_in = []
            fit_y2_in = []
            fit_y1_out_1 = []
            fit_y2_out_1 = []
            fit_y1_out_2 = []
            fit_y2_out_2 = []

            for idx in range( len(x_fill) ):
                if (tau_c_2 - i/2) * omega_imp_a09 <= x_fill[idx] <= (i/2 - tau_c_2) * omega_imp_a09:
                    x_fill_in.append(x_fill[idx])
                    fit_y1_in.append(gA_fit_y1[idx])
                    fit_y2_in.append(gA_fit_y2[idx])

                elif (tau_c_2 - i/2) * omega_imp_a09 > x_fill[idx]:
                    x_fill_out_1.append(x_fill[idx])
                    fit_y1_out_1.append(gA_fit_y1[idx])
                    fit_y2_out_1.append(gA_fit_y2[idx])

                elif x_fill[idx] > (i/2 - tau_c_2) * omega_imp_a09:
                    x_fill_out_2.append(x_fill[idx])
                    fit_y1_out_2.append(gA_fit_y1[idx])
                    fit_y2_out_2.append(gA_fit_y2[idx])

            ax.fill_between(np.array(x_fill_in), np.array(fit_y1_in), np.array(fit_y2_in), color=color, alpha=0.3)
            ax.fill_between(np.array(x_fill_out_1), np.array(fit_y1_out_1), np.array(fit_y2_out_1), color=grey, alpha=0.3)
            ax.fill_between(np.array(x_fill_out_2), np.array(fit_y1_out_2), np.array(fit_y2_out_2), color=grey, alpha=0.3)

    # plot best fit result of ga
    best_x = np.linspace(-6.5*omega_imp_a09, 6.5*omega_imp_a09, 100)
    best_y1 = (fit_result.p['A3_00'].mean + fit_result.p['A3_00'].sdev) * np.ones(len(best_x))
    best_y2 = (fit_result.p['A3_00'].mean - fit_result.p['A3_00'].sdev) * np.ones(len(best_x))

    ax.fill_between(best_x, best_y1, best_y2, color=grey, alpha=0.3)

    x_lim = [-6.5, 6.5]
    ax.set_xlim([num*omega_imp_a09 for num in x_lim])
    ax.set_xlabel(fm_tau_label, **fs_p)
    ax.set_ylim(gA_ylim)
    ax.set_ylabel(ra_label, **fs_p)
    ax.tick_params(direction='in', **ls_p)

    ax_x = ax.twiny() # x-axis on the top
    ax_x.set_xlim(x_lim)
    ax_x.set_xlabel(tau_label, labelpad=-35, **fs_p)
    ax_x.tick_params(direction='in', **ls_p)
    ax_x.tick_params(axis="x", pad=-20)
    plt.savefig(f"./new_plots/psychedelic_moose.pdf", transparent=True)
    plt.show()

def moose_23_late_plot(data_avg_dict_completed, fit_result, fitter, pt3_data_range):
    gA = {}
    gA_err = {}
    gA_fit = {}
    gA_fit_err = {}

    gA_tsep = []
    gA_tau = []
    gV_tsep = []
    gV_tau = []

    tau_c_1 = 1 # tau_c for plot including grey part

    for i in pt3_data_range:
        gA['tsep_'+str(i)] = []
        gA_err['tsep_'+str(i)] = []

        for j in range(tau_c_1, i - tau_c_1 + 1):
            temp = ( data_avg_dict_completed['pt3_A3_tsep_'+str(i)][j] + data_avg_dict_completed['pt3_A3_tsep_'+str(i)][i-j] ) / (2 * data_avg_dict_completed['pt2_tsep_'+str(i)])

            gA['tsep_'+str(i)].append(temp.mean)
            gA_err['tsep_'+str(i)].append(temp.sdev)

        for j in np.linspace(tau_c_1, i - tau_c_1, linspace_num):
            gA_tsep.append(i)
            gA_tau.append(j)
            gV_tsep.append(i)
            gV_tau.append(j)

    pt2_fitter = fitter.pt2_fit_function(np.array([t for t in pt3_data_range]), fit_result.p)['pt2']

    pt3_gA_fitter = fitter.pt3_fit_function(np.array(gA_tsep), np.array(gV_tsep), np.array(gA_tau), np.array(gV_tau), fit_result.p)['pt3_A3']

    index = 0

    for i in pt3_data_range:
        gA_fit['tsep_'+str(i)] = []
        gA_fit_err['tsep_'+str(i)] = []

        for j in range(linspace_num):
            index = int((i - pt3_data_range[0]) / 2 * linspace_num + j) # /2 because 12 - 10 = 2 but 12 is the second one

            temp = pt3_gA_fitter[index] / pt2_fitter[int( (i - pt3_data_range[0]) / 2) ]
            gA_fit['tsep_'+str(i)].append(temp.mean)
            gA_fit_err['tsep_'+str(i)].append(temp.sdev)

    ########### here start a fig ###########
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    for i in pt3_data_range:
        idx = i - 2 # psychedelic_list start from tsep = 2

        x_errorbar = np.arange(tau_c_1 - i/2, i/2 - tau_c_1 + 1) * omega_imp_a09
        x_fill = np.linspace(tau_c_1 - i/2, i/2 - tau_c_1, linspace_num) * omega_imp_a09

        if i == 10:
            tau_c_2 = 3.5
            color = blue
            marker = 'p'

        elif i == 12:
            tau_c_2 = 3.5
            color = green
            marker = 'D'

        elif i == 14:
            tau_c_2 = 4.5 # 5-0.5 to make color change at half way
            color = red
            marker = 'o'

        for idx in range( len(x_errorbar) ):
            if (tau_c_2 - i/2) * omega_imp_a09 <= x_errorbar[idx] <= (i/2 - tau_c_2) * omega_imp_a09:
                ax.errorbar(np.array(x_errorbar[idx]), np.array(gA['tsep_' + str(i)][idx]), yerr=np.array(gA_err['tsep_' + str(i)][idx]), marker=marker, color=color, label=r'$t_{\rm sep} = %d$' %i , **errorp)

            else:
                ax.errorbar(np.array(x_errorbar[idx]), np.array(gA['tsep_' + str(i)][idx]), yerr=np.array(gA_err['tsep_' + str(i)][idx]), marker=marker, color=grey, **errorp)

        gA_fit_y1 = np.array(gA_fit['tsep_'+str(i)]) + np.array(gA_fit_err['tsep_'+str(i)])
        gA_fit_y2 = np.array(gA_fit['tsep_'+str(i)]) - np.array(gA_fit_err['tsep_'+str(i)])

        x_fill_in = []
        x_fill_out_1 = []
        x_fill_out_2 = []

        fit_y1_in = []
        fit_y2_in = []
        fit_y1_out_1 = []
        fit_y2_out_1 = []
        fit_y1_out_2 = []
        fit_y2_out_2 = []

        for idx in range( len(x_fill) ):
            if (tau_c_2 - i/2) * omega_imp_a09 <= x_fill[idx] <= (i/2 - tau_c_2) * omega_imp_a09:
                x_fill_in.append(x_fill[idx])
                fit_y1_in.append(gA_fit_y1[idx])
                fit_y2_in.append(gA_fit_y2[idx])

            elif (tau_c_2 - i/2) * omega_imp_a09 > x_fill[idx]:
                x_fill_out_1.append(x_fill[idx])
                fit_y1_out_1.append(gA_fit_y1[idx])
                fit_y2_out_1.append(gA_fit_y2[idx])

            elif x_fill[idx] > (i/2 - tau_c_2) * omega_imp_a09:
                x_fill_out_2.append(x_fill[idx])
                fit_y1_out_2.append(gA_fit_y1[idx])
                fit_y2_out_2.append(gA_fit_y2[idx])

        ax.fill_between(np.array(x_fill_in), np.array(fit_y1_in), np.array(fit_y2_in), color=color, alpha=0.3)
        ax.fill_between(np.array(x_fill_out_1), np.array(fit_y1_out_1), np.array(fit_y2_out_1), color=grey, alpha=0.3)
        ax.fill_between(np.array(x_fill_out_2), np.array(fit_y1_out_2), np.array(fit_y2_out_2), color=grey, alpha=0.3)

    # plot best fit result of ga
    best_x = np.linspace(-6.5*omega_imp_a09, 6.5*omega_imp_a09, 100)
    best_y1 = (fit_result.p['A3_00'].mean + fit_result.p['A3_00'].sdev) * np.ones(len(best_x))
    best_y2 = (fit_result.p['A3_00'].mean - fit_result.p['A3_00'].sdev) * np.ones(len(best_x))

    ax.fill_between(best_x, best_y1, best_y2, color=grey, alpha=0.3)

    x_lim = [-6.5, 6.5]
    ax.set_xlim([num*omega_imp_a09 for num in x_lim])
    ax.set_xlabel(fm_tau_label, **fs_p)
    ax.set_ylim(gA_ylim)
    ax.set_ylabel(ra_label, **fs_p)
    ax.tick_params(direction='in', **ls_p)

    legend_without_duplicate_labels(ax, 'lower center', 3)

    ax_x = ax.twiny() # x-axis on the top
    ax_x.set_xlim(x_lim)
    ax_x.set_xlabel(tau_label, labelpad=-35, **fs_p)
    ax_x.tick_params(direction='in', **ls_p)
    ax_x.tick_params(axis="x", pad=-20)
    plt.savefig(f"./new_plots/psychedelic_moose_23_late.pdf", transparent=True)
    plt.show()

def tau_c_plot(file_path, f_range, d_range, fit_result): # fit range and data range
    from module.prepare_data import Prepare_data
    from module.fit import Fit
    from module.prior_setting import prior_ho_width_1
    prior = prior_ho_width_1

    file_name = 'a09m310_e_gA_srcs0-15_full_tau.h5'
    pt2_nstates = 5
    pt3_nstates = pt2_nstates
    sum_nstates = 5
    include_2pt = True
    include_3pt = True
    include_sum = True

    pt2_data_range = [0, 96]
    pt3_data_range = [2, 15]

    prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

    data_avg_dict = prepare_data.read_data_with_average()

    plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    for sum_tau_cut in range(f_range[0], f_range[1]):
        fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut, include_2pt, include_3pt, include_sum)

        tmin = max(2, 2*sum_tau_cut) - 0.5
        tmax = 23

        x_fill = np.arange(tmin, tmax, plot_space)[:int(-1 / plot_space)]

        pt2_fitter = fitter.pt2_fit_function(np.arange(tmin, tmax, plot_space), fit_result.p)['pt2']

        sum_A3_fitter = fitter.summation_same_can(np.arange(tmin, tmax, plot_space), np.arange(tmin, tmax, plot_space), fit_result.p)['sum_A3']

        temp = []

        for i in range(len(np.arange(tmin, tmax, plot_space)) - int(1 / plot_space)):
            temp.append(sum_A3_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i])


        y1 = np.array([val.mean for val in temp]) + np.array([val.sdev for val in temp])
        y2 = np.array([val.mean for val in temp]) - np.array([val.sdev for val in temp])

        ax.fill_between(x_fill * omega_imp_a09, y1, y2, color=color_list[sum_tau_cut], alpha=0.3)

    for sum_tau_cut in range(d_range[0], d_range[1]):

        data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)

        sum_data_list = []

        tmin = max(2, 2*sum_tau_cut)

        for t in range(tmin, 14):
            sum_data_list.append(data_avg_dict_completed['sum_A3_fit_'+str(t)])

        temp_mean = np.array([val.mean for val in sum_data_list])
        temp_sdev = np.array([val.sdev for val in sum_data_list])

        if sum_tau_cut == 1:
            ax.errorbar(np.arange(tmin, 14) * omega_imp_a09, temp_mean, yerr=temp_sdev, label=r'$\tau_c=%d$' %sum_tau_cut, marker=marker_list[sum_tau_cut], color=color_list[sum_tau_cut], mfc=color_list[sum_tau_cut], **errorb)

        else:
            ax.errorbar(np.arange(tmin, 14) * omega_imp_a09, temp_mean, yerr=temp_sdev, label=r'$\tau_c=%d$' %sum_tau_cut, marker=marker_list[sum_tau_cut], color=color_list[sum_tau_cut], **errorp)

    # plot best fit result of ga
    best_x = np.linspace(1 * omega_imp_a09, 23 * omega_imp_a09, 100)
    best_y1 = (fit_result.p['A3_00'].mean + fit_result.p['A3_00'].sdev) * np.ones(len(best_x))
    best_y2 = (fit_result.p['A3_00'].mean - fit_result.p['A3_00'].sdev) * np.ones(len(best_x))

    ax.fill_between(best_x, best_y1, best_y2, color=grey, alpha=0.3)

    ax.set_ylim(gA_ylim)
    ax.set_xlim([1 * omega_imp_a09, 22 * omega_imp_a09])
    ax.set_ylabel(fha_label, **gA_fs_p)
    ax.set_xlabel(tsep_fm_label, **fs_p)
    ax.text(0.525, 1.295, tsep_label, **fs_p)
    ax.tick_params(direction='in', **ls_p)
    legend_without_duplicate_labels(ax, 'lower right', 3)

    ax1 = ax.twiny()
    ax1.set_xlim([1, 22])
    ax1.tick_params(axis="x", pad=-20)
    ax1.tick_params(direction='in', **ls_p)

    plot_range = str(f_range[0]) + '_to_' + str(f_range[1]-1)
    plt.savefig(f"./new_plots/tau_c_sum_plot_{plot_range}.pdf", transparent=True)
    plt.show()

def FH_R_plot(file_path, f_range, d_range, fit_result, combined_best_data_avg_dict_completed): # fit range and data range
    from module.prepare_data import Prepare_data
    from module.fit import Fit
    from module.prior_setting import prior_ho_width_1
    prior = prior_ho_width_1

    file_name = 'a09m310_e_gA_srcs0-15_full_tau.h5'
    pt2_nstates = 5
    pt3_nstates = pt2_nstates
    sum_nstates = 5
    include_2pt = True
    include_3pt = True
    include_sum = True

    pt2_data_range = [0, 96]
    pt3_data_range = [2, 15]

    prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

    data_avg_dict = prepare_data.read_data_with_average()



    plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    for sum_tau_cut in range(f_range[0], f_range[1]):
        fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut, include_2pt, include_3pt, include_sum)

        tmin = max(2, 2*sum_tau_cut) - 0.5
        tmax = 23

        x_fill = np.arange(tmin, tmax, plot_space)[:int(-1 / plot_space)]

        pt2_fitter = fitter.pt2_fit_function(np.arange(tmin, tmax, plot_space), fit_result.p)['pt2']

        sum_A3_fitter = fitter.summation_same_can(np.arange(tmin, tmax, plot_space), np.arange(tmin, tmax, plot_space), fit_result.p)['sum_A3']

        temp = []

        for i in range(len(np.arange(tmin, tmax, plot_space)) - int(1 / plot_space)):
            temp.append(sum_A3_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i])


        y1 = np.array([val.mean for val in temp]) + np.array([val.sdev for val in temp])
        y2 = np.array([val.mean for val in temp]) - np.array([val.sdev for val in temp])

        ax.fill_between(x_fill * omega_imp_a09, y1, y2, color=color_list[sum_tau_cut], alpha=0.3)

    for sum_tau_cut in range(d_range[0], d_range[1]):

        data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)

        sum_data_list = []

        tmin = max(2, 2*sum_tau_cut)

        for t in range(tmin, 14):
            sum_data_list.append(data_avg_dict_completed['sum_A3_fit_'+str(t)])

        temp_mean = np.array([val.mean for val in sum_data_list])
        temp_sdev = np.array([val.sdev for val in sum_data_list])

        if sum_tau_cut == 1:
            ax.errorbar(np.arange(tmin, 14) * omega_imp_a09, temp_mean, yerr=temp_sdev, label=r'$\tau_c=%d$' %sum_tau_cut, marker=marker_list[sum_tau_cut], color=color_list[sum_tau_cut], mfc=color_list[sum_tau_cut], **errorb)

        else:
            ax.errorbar(np.arange(tmin, 14) * omega_imp_a09, temp_mean, yerr=temp_sdev, label=r'$\tau_c=%d$' %sum_tau_cut, marker=marker_list[sum_tau_cut], color=color_list[sum_tau_cut], **errorp)

    ## R(tau=t/2) part
    gA = {}
    gA_err = {}

    for t in range(pt3_data_range[0], pt3_data_range[1]):
        if t % 2 == 0:
            gA['tsep_'+str(t)] = []
            gA_err['tsep_'+str(t)] = []

            tau = int(t / 2)
            temp = ( combined_best_data_avg_dict_completed['pt3_A3_tsep_'+str(t)][tau] ) / (combined_best_data_avg_dict_completed['pt2_tsep_'+str(t)])

            gA['tsep_'+str(t)].append(temp.mean)
            gA_err['tsep_'+str(t)].append(temp.sdev)

            color = psychedelic_list[t-2]

            ax.errorbar(np.array([t]) * omega_imp_a09, np.array(gA['tsep_'+str(t)]), yerr=np.array(gA_err['tsep_'+str(t)]), marker='o', color=color, label=r'$R_{A_3}(t_{\rm sep}, \tau=t_{\rm sep}/2)$', **errorp)

    # plot best fit result of ga
    best_x = np.linspace(1 * omega_imp_a09, 23 * omega_imp_a09, 100)
    best_y1 = (fit_result.p['A3_00'].mean + fit_result.p['A3_00'].sdev) * np.ones(len(best_x))
    best_y2 = (fit_result.p['A3_00'].mean - fit_result.p['A3_00'].sdev) * np.ones(len(best_x))

    ax.fill_between(best_x, best_y1, best_y2, color=grey, alpha=0.3)

    ax.set_ylim(gA_ylim)
    ax.set_xlim([1 * omega_imp_a09, 22 * omega_imp_a09])
    ax.set_ylabel(fha_label, **gA_fs_p)
    ax.set_xlabel(tsep_fm_label, **fs_p)
    ax.text(0.525, 1.295, tsep_label, **fs_p)
    ax.tick_params(direction='in', **ls_p)
    legend_without_duplicate_labels(ax, 'lower right', 1)

    ax1 = ax.twiny()
    ax1.set_xlim([1, 22])
    ax1.tick_params(axis="x", pad=-20)
    ax1.tick_params(direction='in', **ls_p)

    plt.savefig(f"./new_plots/FH_R_plot.pdf", transparent=True)
    plt.show()

def excited_states_plot(pt3_list, ratio_list):
    [pt3_gA_tra_fit, pt3_gA_sca_fit, pt3_gA_gs_fit, pt3_gA_all_fit, pt3_gA_tra_data, pt3_gA_sca_data, pt3_gA_gs_data, pt3_gA_all_data] = pt3_list

    [ratio_gA_tra_fit, ratio_gA_sca_fit, ratio_gA_gs_fit, ratio_gA_all_fit, ratio_gA_tra_data, ratio_gA_sca_data, ratio_gA_gs_data, ratio_gA_all_data] = ratio_list

    plot_xmax = 3

    pt3_tsep_start = 2
    pt3_tsep_end = 15 # sum_end = pt3_end - 1

    pt3_tsep_data = np.arange(pt3_tsep_start, pt3_tsep_end, 2) # errorbar do not need plot density

    pt3_tsep_fit = np.linspace(0, 35, linspace_num) # till 3fm

    plot_gap = (35 - 0) / (linspace_num - 1)

    sum_tsep_data = np.arange(pt3_tsep_start, pt3_tsep_end-1)

    sum_tsep_fit = np.linspace(0, 35, linspace_num)[:-int(1/plot_gap)]

    pt3_tsep_data = pt3_tsep_data * omega_imp_a09
    pt3_tsep_fit = pt3_tsep_fit * omega_imp_a09
    sum_tsep_data = sum_tsep_data * omega_imp_a09
    sum_tsep_fit = sum_tsep_fit * omega_imp_a09


    # transite ratio to sum-sub
    sum_gA_tra_fit = []
    sum_gA_sca_fit = []
    sum_gA_gs_fit = []
    sum_gA_all_fit = []

    sum_list = [sum_gA_tra_fit, sum_gA_sca_fit, sum_gA_gs_fit, sum_gA_all_fit]

    ratio_list = [ratio_gA_tra_fit, ratio_gA_sca_fit, ratio_gA_gs_fit, ratio_gA_all_fit]

    for i in range(len(sum_list)):
        for t in range(len(sum_tsep_fit)):
            temp = (ratio_list[i][t+int(1/plot_gap)] - ratio_list[i][t]) # here gA[t+1] - gA[t] != 1
            sum_list[i].append(temp)

    [sum_gA_tra_fit, sum_gA_sca_fit, sum_gA_gs_fit, sum_gA_all_fit] = [np.array(lis) for lis in sum_list]


    sum_gA_tra_data = []
    sum_gA_sca_data = []
    sum_gA_gs_data = []
    sum_gA_all_data = []

    sum_list = [sum_gA_tra_data, sum_gA_sca_data, sum_gA_gs_data, sum_gA_all_data]

    ratio_list = [ratio_gA_tra_data, ratio_gA_sca_data, ratio_gA_gs_data, ratio_gA_all_data]

    for i in range(len(sum_list)):
        for t in range(len(sum_tsep_data)):
            temp = (ratio_list[i][t+1] - ratio_list[i][t]) # here gA[t+1] - gA[t] != 1
            sum_list[i].append(temp)

    [sum_gA_tra_data, sum_gA_sca_data, sum_gA_gs_data, sum_gA_all_data] = [np.array(lis) for lis in sum_list]


    ###### 3pt/2pt plot ######
    fig = plt.figure(figsize=fig_size)
    ax  = plt.axes(plt_axes)

    pt2_es_fit = pt3_gA_all_fit - pt3_gA_gs_fit - pt3_gA_tra_fit - pt3_gA_sca_fit
    pt2_es_data = pt3_gA_all_data - pt3_gA_gs_data - pt3_gA_tra_data - pt3_gA_sca_data


    temp_mean = np.array([val.mean for val in (pt3_gA_all_fit - pt3_gA_gs_fit) / pt3_gA_gs_fit])
    temp_sdev = np.array([val.sdev for val in (pt3_gA_all_fit - pt3_gA_gs_fit) / pt3_gA_gs_fit])
    ax.fill_between(pt3_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, facecolor=grey, edgecolor='white', hatch='\ ', alpha=0.5, label=r'$\rm 3pt$')

    temp_mean = np.array([val.mean for val in (pt3_gA_sca_fit + pt2_es_fit) / pt3_gA_gs_fit])
    temp_sdev = np.array([val.sdev for val in (pt3_gA_sca_fit + pt2_es_fit) / pt3_gA_gs_fit])
    ax.fill_between(pt3_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, facecolor='None', edgecolor=grey, hatch='\\', alpha=0.8, label=r'$\rm 2pt+sc$')


    plot_fit_list = [pt3_gA_tra_fit, pt3_gA_sca_fit, pt2_es_fit]
    plot_data_list = [pt3_gA_tra_data, pt3_gA_sca_data, pt2_es_data]
    label_list = [r'$\rm{3pt\ tr}$', r'$\rm{3pt\ sc}$', r'$\rm{2pt\ es}$']
    color_list = [blue, red, grape]

    for i in range(len(plot_fit_list)):
        temp_mean = np.array([val.mean for val in plot_fit_list[i] / pt3_gA_gs_fit])
        temp_sdev = np.array([val.sdev for val in plot_fit_list[i] / pt3_gA_gs_fit])
        ax.fill_between(pt3_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, color=color_list[i], alpha=0.3, label=label_list[i])

        temp_mean = np.array([val.mean for val in plot_data_list[i] / pt3_gA_gs_data])
        temp_sdev = np.array([val.sdev for val in plot_data_list[i] / pt3_gA_gs_data])
        ax.errorbar(pt3_tsep_data, temp_mean, yerr=temp_sdev, marker='o', color=color_list[i], **errorp)

    ax.plot(np.linspace(0, 3.8, 10), np.ones(10) * 0., color='k', linestyle='-', linewidth=1)
    ax.plot(np.linspace(0, 3.8, 10), np.ones(10) * 0.01, color=red, linestyle='--', linewidth=1)
    ax.plot(np.linspace(0, 3.8, 10), np.ones(10) * -0.01, color=red, linestyle='--', linewidth=1)

    ax.set_xlabel(tsep_fm_label, **fs_p)
    ax.tick_params(direction='in', **ls_p)

    ax.text(20 * omega_imp_a09, 0.03, r'$\frac{R^{\rm es}_{A_3}(t_{\rm sep}, \tau=t_{\rm sep}/2)} {\mathring{g}_A}$', fontsize=27)

    ax.set_ylim([-0.151, 0.149])
    ax.set_xlim([0.01, plot_xmax + 0.8])

    ax_x = ax.twiny()
    ax_x.set_xlim([0.01 / omega_imp_a09, 3.8 / omega_imp_a09])
    ax_x.set_xlabel(tsep_label, labelpad=-35, **fs_p)
    ax_x.tick_params(direction='in', **ls_p)
    ax_x.tick_params(axis="x", pad=-20)

    legend_without_duplicate_labels(ax, 'lower right', 3)

    plt.savefig(f"./new_plots/ratio_excited_states.pdf", transparent=True)
    plt.show()


    ###### sum-sub plot ######
    fig = plt.figure(figsize=fig_size)
    ax  = plt.axes(plt_axes)

    pt2_es_fit = sum_gA_all_fit - sum_gA_gs_fit - sum_gA_tra_fit - sum_gA_sca_fit
    pt2_es_data = sum_gA_all_data - sum_gA_gs_data - sum_gA_tra_data - sum_gA_sca_data


    temp_mean = np.array([val.mean for val in (sum_gA_all_fit - sum_gA_gs_fit) / sum_gA_gs_fit])
    temp_sdev = np.array([val.sdev for val in (sum_gA_all_fit - sum_gA_gs_fit) / sum_gA_gs_fit])
    ax.fill_between(sum_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, facecolor=grey,edgecolor='white', hatch='\ ', alpha=0.5, label=r'$\rm fh$')

    temp_mean = np.array([val.mean for val in (sum_gA_sca_fit + pt2_es_fit) / sum_gA_gs_fit])
    temp_sdev = np.array([val.sdev for val in (sum_gA_sca_fit + pt2_es_fit) / sum_gA_gs_fit])
    ax.fill_between(sum_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, facecolor='None',edgecolor=grey, hatch='\\', alpha=0.8, label=r'$\rm 2pt+sc$')


    plot_fit_list = [sum_gA_tra_fit, sum_gA_sca_fit, pt2_es_fit]
    plot_data_list = [sum_gA_tra_data, sum_gA_sca_data, pt2_es_data]
    label_list = [r'$\rm{fh\ tr}$', r'$\rm{fh\ sc}$', r'$\rm{2pt\ es}$']
    color_list = [blue, red, grape]

    for i in range(len(plot_fit_list)):
        temp_mean = np.array([val.mean for val in plot_fit_list[i] / sum_gA_gs_fit])
        temp_sdev = np.array([val.sdev for val in plot_fit_list[i] / sum_gA_gs_fit])
        ax.fill_between(sum_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, color=color_list[i], alpha=0.3, label=label_list[i])

        temp_mean = np.array([val.mean for val in plot_data_list[i] / sum_gA_gs_data])
        temp_sdev = np.array([val.sdev for val in plot_data_list[i] / sum_gA_gs_data])
        ax.errorbar(sum_tsep_data, temp_mean, yerr=temp_sdev, marker='o', color=color_list[i], **errorp)

    ax.plot(np.linspace(0, 3.8, 10), np.ones(10) * 0., color='k', linestyle='-', linewidth=1)
    ax.plot(np.linspace(0, 3.8, 10), np.ones(10) * 0.01, color=red, linestyle='--', linewidth=1)
    ax.plot(np.linspace(0, 3.8, 10), np.ones(10) * -0.01, color=red, linestyle='--', linewidth=1)

    ax.set_xlabel(tsep_fm_label, **fs_p)
    ax.tick_params(direction='in', **ls_p)

    ax.text(20 * omega_imp_a09, 0.03, r'$\frac{{\rm FH}^{\rm es}_{A_3}(t_{\rm sep}, \tau_c=1)} {\mathring{g}_A}$', fontsize=27)

    ax.set_ylim([-0.151, 0.149])
    ax.set_xlim([0.01, plot_xmax + 0.8])

    ax_x = ax.twiny()
    ax_x.set_xlim([0.01 / omega_imp_a09, 3.8 / omega_imp_a09])
    ax_x.set_xlabel(tsep_label, labelpad=-35, **fs_p)
    ax_x.tick_params(direction='in', **ls_p)
    ax_x.tick_params(axis="x", pad=-20)

    legend_without_duplicate_labels(ax, 'lower right', 3)

    plt.savefig(f"./new_plots/fh_excited_states.pdf", transparent=True)
    plt.show()

def m_z_plot(data_avg_dict_completed, fit_result, fitter):
    pt2_data_range = [0, 96]

    m0_eff = []
    m0_eff_err = []
    m0_eff_fit = []
    m0_eff_fit_err = []

    x_errorbar = np.arange(pt2_data_range[0], pt2_data_range[1]-1) * omega_imp_a09
    x_fill = np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)[:int(-1/plot_space)] * omega_imp_a09

    for i in range(pt2_data_range[0], pt2_data_range[1]-1):
        temp = gv.log(data_avg_dict_completed['pt2_tsep_'+str(i)] / data_avg_dict_completed['pt2_tsep_'+str(i+1)])
        m0_eff.append(temp.mean)
        m0_eff_err.append(temp.sdev)

    fig = plt.figure(figsize=fig_size)
    ax  = plt.axes(plt_axes)
    ax.errorbar(x_errorbar, np.array(m0_eff), yerr=np.array(m0_eff_err), marker='o', color="k", **errorp)

    pt2_fitter = fitter.pt2_fit_function(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space), fit_result.p)['pt2']

    for i in range(len(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)) - int(1/plot_space)):
        temp = gv.log(pt2_fitter[i] / pt2_fitter[i+ int(1/plot_space)])
        m0_eff_fit.append(temp.mean)
        m0_eff_fit_err.append(temp.sdev)

    m0_eff_fit_y1 = []
    m0_eff_fit_y2 = []

    for i in range(len(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)) - int(1/plot_space)):
        m0_eff_fit_y1.append(m0_eff_fit[i] + m0_eff_fit_err[i])
        m0_eff_fit_y2.append(m0_eff_fit[i] - m0_eff_fit_err[i])
    ax.fill_between(x_fill, np.array(m0_eff_fit_y1), np.array(m0_eff_fit_y2), color=blue, alpha=0.3, label='fit')

    # grey band of prior
    ax.fill_between(x_fill, np.ones(len(x_fill)) * (0.5-0.05), np.ones(len(x_fill)) * (0.5+0.05), color=grey, alpha=0.3, label='prior' )

    x_lim = [2, 25.5]
    ax.set_xlim([num*omega_imp_a09 for num in x_lim])
    ax.set_xlabel(tsep_fm_label, **fs_p)

    ax.set_ylim([0.43, 0.62])
    ax.set_ylabel(meff_label, **fs_p)

    ax1 = ax.twiny()
    ax1.set_xlim(x_lim)
    ax1.set_xlabel(tsep_label, labelpad=-35, **fs_p)
    ax1.tick_params(axis="x", pad=-20)

    ax.tick_params(direction='in', **ls_p)
    ax1.tick_params(direction='in', **ls_p)
    ax.yaxis.set_major_formatter(FormatStrFormatter('$%.2f$'))

    plt.savefig(f"./new_plots/meff_23s.pdf", transparent=True)
    plt.show()


    #zeff
    zeff = []

    for i in range(pt2_data_range[0], pt2_data_range[1]-1):
        meff = gv.log(data_avg_dict_completed['pt2_tsep_'+str(i)] / data_avg_dict_completed['pt2_tsep_'+str(i+1)])
        zeff.append(np.sqrt(data_avg_dict_completed['pt2_tsep_'+str(i)]*np.exp(meff*i)))

    fig = plt.figure(figsize=fig_size)
    ax  = plt.axes(plt_axes)
    ax.errorbar(x_errorbar, [i.mean * 1000 for i in zeff], yerr=[i.sdev * 1000 for i in zeff], marker='o', color="k", **errorp)

    z0_eff_fit = []
    z0_eff_fit_err = []

    x = np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)[:int(-1/plot_space)]
    pt2_fitter = fitter.pt2_fit_function(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space), fit_result.p)['pt2']

    for idx, i in enumerate(range(len(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)) - int(1/plot_space))):
        meff = gv.log(pt2_fitter[i] / pt2_fitter[i+ int(1/plot_space)])
        zeff = np.sqrt(pt2_fitter[i] * np.exp(meff*x[idx]))
        z0_eff_fit.append(zeff.mean)
        z0_eff_fit_err.append(zeff.sdev)

    z0_eff_fit_y1 = []
    z0_eff_fit_y2 = []

    for i in range(len(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)) - int(1/plot_space)):
        z0_eff_fit_y1.append(z0_eff_fit[i] + z0_eff_fit_err[i])
        z0_eff_fit_y2.append(z0_eff_fit[i] - z0_eff_fit_err[i])

    ax.fill_between(x_fill, np.array(z0_eff_fit_y1) * 1000, np.array(z0_eff_fit_y2) * 1000, color=blue, alpha=0.3, label='fit')

    # grey band of prior
    ax.fill_between(x_fill, np.ones(len(x_fill)) * 1000 * (0.00034-0.00034), np.ones(len(x_fill)) * 1000 * (0.00034+0.00034), color=grey, alpha=0.3, label='prior' )

    x_lim = [2, 25.5]
    ax.set_xlim([num*omega_imp_a09 for num in x_lim])
    ax.set_xlabel(tsep_fm_label, **fs_p)

    ax.set_ylim([-0.2, 0.79])
    ax.set_ylabel(zeff_label, labelpad=-5, **fs_p)

    ax1 = ax.twiny()
    ax1.set_xlim(x_lim)
    ax1.set_xlabel(tsep_label, labelpad=-35, **fs_p)

    ax.tick_params(direction='in', **ls_p)

    ax1.tick_params(direction='in', **ls_p)
    ax1.tick_params(axis="x", pad=-20)

    plt.savefig(f"./new_plots/zeff_23s.pdf", transparent=True)
    plt.show()

def sum_gA_plot(data_avg_dict_completed, fit_result, fitter):
    pt3_data_range = [2, 15]
    pt2_nstates = 5
    pt3_nstates = pt2_nstates
    sum_nstates = 5
    gap = 1 # summation(t+gap) - summation(t)

    gA = []
    gA_err = []

    gA_fit = []
    gA_fit_err = []

    x_errorbar = np.arange(pt3_data_range[0], pt3_data_range[1]-gap) * omega_imp_a09
    x_fill = np.arange(pt3_data_range[0], pt3_data_range[1], plot_space)[:int(-gap / plot_space)] * omega_imp_a09

    for i in range(pt3_data_range[0], pt3_data_range[1]-gap):
        temp1 = data_avg_dict_completed['sum_A3_fit_'+str(i)]
        gA.append(temp1.mean)
        gA_err.append(temp1.sdev)


    fig = plt.figure(figsize=fig_size)
    ax  = plt.axes(plt_axes)

    x_odd = []
    x_even = []
    gA_odd = []
    gA_even = []
    gA_err_odd = []
    gA_err_even = []

    for i in range(len(x_errorbar)):
        t = round(x_errorbar[i] / omega_imp_a09)

        if t%2 == 1:
            x_odd.append(t)
            gA_odd.append(gA[i])
            gA_err_odd.append(gA_err[i])

        elif t%2 == 0:
            x_even.append(t)
            gA_even.append(gA[i])
            gA_err_even.append(gA_err[i])

    ax.errorbar(np.array(x_odd) * omega_imp_a09, np.array(gA_odd), yerr=np.array(gA_err_odd), marker='o', color='k', **errorp)
    ax.errorbar(np.array(x_even) * omega_imp_a09, np.array(gA_even), yerr=np.array(gA_err_even), marker='o', color='k', **errorp)#mfc=blue, **errorb)


    if sum_nstates == pt2_nstates:
        pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), fit_result.p)['pt2']
        sum_A3_fitter = fitter.summation_same_can(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), fit_result.p)['sum_A3']

    else:
        pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), fit_result.p, sum_nstates)['pt2']
        sum_A3_fitter = fitter.summation(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), fit_result.p)['sum_A3']

    for i in range(len(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space)) - int(gap / plot_space)):
        temp1 = (sum_A3_fitter[i+int(gap / plot_space)] / pt2_fitter[i+int(gap / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i]) / gap
        gA_fit.append(temp1.mean)
        gA_fit_err.append(temp1.sdev)

    gA_fit_y1 = np.array(gA_fit) + np.array(gA_fit_err)
    gA_fit_y2 = np.array(gA_fit) - np.array(gA_fit_err)

    ax.fill_between(x_fill, gA_fit_y1, gA_fit_y2, color=blue, alpha=0.3, label='fit')

    # grey band of prior
    ax.fill_between(np.linspace(0, 1.2, 10), np.ones(10) * (1.2-0.2), np.ones(10) * (1.2+0.2), color=grey, alpha=0.3, label='prior' )

    x_lim = [1.5, 13.5]
    ax1 = ax.twiny()
    ax1.set_xlim(x_lim)
    ax1.set_xlabel(tsep_label, labelpad=-35, **fs_p)
    ax1.tick_params(axis="x", pad=-20)
    ax.tick_params(axis='both', which='major', direction='in', **ls_p)
    ax1.tick_params(axis='both', which='major', direction='in', **ls_p)

    ax.set_xlim([num*omega_imp_a09 for num in x_lim])
    ax.set_xlabel(tsep_fm_label, **fs_p)

    ax.set_ylim(0.9, 1.49)
    ax.set_ylabel(ga_label, **gA_fs_p)

    plt.savefig(f"./new_plots/soaeff_23s.pdf", transparent=True)
    plt.show()

def tmin_tmax_combined_plot(x, value):
    best_gA = value['gA'][1]
    best_gA_err = value['gA_err'][1]
    best_Q = value['Q'][1]

    t_range = [2, 15]
    nstate = 5

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)
    ax1.set_ylabel(ga_label, **gA_fs_p)
    ax1.set_ylim(gA_ylim)

    ax1.tick_params(direction='in', **ls_p)

    ax1.vlines(7.5, gA_ylim[0], gA_ylim[1], colors = "k", linestyles = "dashed")

    ax1.errorbar(np.array(x), np.array(value['gA']), yerr=np.array(value['gA_err']), marker='o', color=color_list[nstate-2], **errorp)

    #best fit
    ax1.fill_between(np.arange(1.5, 15.5), (best_gA + best_gA_err)*np.ones([14]), (best_gA - best_gA_err)*np.ones([14]), color=color_list[nstate-2], alpha=0.2)

    ax1.errorbar(np.array([3]), np.array([best_gA]), yerr=np.array([best_gA_err]), marker='o', mfc=color_list[nstate-2], color=color_list[nstate-2], **errorb)

    ax1.errorbar(np.array([14]), np.array([best_gA]), yerr=np.array([best_gA_err]), marker='o', mfc=color_list[nstate-2], color=color_list[nstate-2], **errorb)


    ax2.set_ylabel(q_label, **fs_p)
    ax2.set_ylim([0, 1.1])

    ax2.vlines(7.5, 0, 1.1, colors = "k", linestyles = "dashed")

    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    ax2.scatter(np.array(x), np.array(value['Q']), marker='o', c='None', edgecolors=color_list[nstate-2])

    #best fit
    ax2.scatter(np.array([3]), np.array([best_Q]), marker='o', c=color_list[nstate-2])

    ax2.scatter(np.array([14]), np.array([best_Q]), marker='o', c=color_list[nstate-2])

    ax2.tick_params(direction='in', **ls_p)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(tmin_label+'\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ '+tmax_label, labelpad=1.5, **fs_p)
    plt.xlim([1.5, 14.5])

    plt.savefig('./new_plots/23s_gA-both_tmin_tmax_combined.pdf', transparent=True)

def tmax_scattered_plot(best_n, tmax_name, situation_list, save_name):

    value={}
    value['Q']=[]
    value['gA']=[]
    value['gA_err']=[]

    x=[]

    for situation in situation_list:
        tmax_dict = {}
        tmax_dict['2pt'] = situation.pt2_tmax
        tmax_dict['3pt'] = situation.pt3_A3_tsep_max
        tmax_dict['sum'] = situation.sum_A3_tsep_max

        x.append(tmax_dict[tmax_name])  # here is the varying parameter

    x.sort()
    for i in range(len(x)):
        for situation in situation_list:
            tmax_dict = {}
            tmax_dict['2pt'] = situation.pt2_tmax
            tmax_dict['3pt'] = situation.pt3_A3_tsep_max
            tmax_dict['sum'] = situation.sum_A3_tsep_max
            if tmax_dict[tmax_name] == x[i]:
                value['Q'].append(situation.Q_value)
                value['gA'].append(situation.A300)
                value['gA_err'].append(situation.A300_err)

    print(x)

    #####################################################################################

    gA_mean = 1.2534405806490314
    gA_sdev = 0.01930206389474567

    # gA - tmax
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)
    ax1.set_ylabel(ga_label, **gA_fs_p)
    ax1.set_ylim(gA_ylim)
    ax1.tick_params(direction='in', **ls_p)

    ax1.errorbar(np.array(x), np.array(value['gA']), yerr=np.array(value['gA_err']), marker='o', color=color_list[best_n-2], **errorp)

    ax1.fill_between(np.arange(4.5, 15.5, 1), (gA_mean + gA_sdev) * np.ones([11]), (gA_mean - gA_sdev) * np.ones([11]), color=color_list[best_n-2], alpha=0.2)


    ax2.set_ylabel(q_label, **fs_p)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange(4.5, 15.5, 1), 0.1 * np.ones([11]), 'r--')

    ax2.scatter(np.array(x), np.array(value['Q']), marker='o', c='None', edgecolors=color_list[best_n-2])

    ax2.tick_params(direction='in', **ls_p)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(tmax_label, **fs_p)
    plt.xlim([4.5, 14.5])

    plt.savefig('./new_plots/23s_gA-'+save_name+'_tmax.pdf', transparent=True)

def late_23_tau_inc_plot(A3, A3_err, Q):
    # stability plot of late tsep [10, 12, 14] # varying nstates and tau_cut
    n_ga = [1.25343903, 0.01930155]
    y1a = n_ga[0] - n_ga[1]
    y2a = n_ga[0] + n_ga[1]

    # nstate = 1, 2, 3, 4

    #####################gA
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)
    ax1.set_ylabel(ga_label, **gA_fs_p)
    ax1.set_ylim(gA_ylim)
    ax1.tick_params(axis='both', direction='in', which='major', **ls_p)

    for i in range(4):
        ax1.errorbar(np.array([i-0.2]), np.array(A3['opt+1'][i]), yerr=np.array(A3_err['opt+1'][i]), marker='^', color=green, label=r'$\tau_{\rm inc} = \tau_{\rm inc}^{\rm opt} +1$', **errorp)
        ax1.errorbar(np.array([i]), np.array(A3['opt'][i]), yerr=np.array(A3_err['opt'][i]), marker='o', color=blue, label=r'$\tau_{\rm inc} = \tau_{\rm inc}^{\rm opt}$', **errorp)
        ax1.errorbar(np.array([i+0.2]), np.array(A3['opt-1'][i]), yerr=np.array(A3_err['opt-1'][i]), marker='P', color=peach, label=r'$\tau_{\rm inc} = \tau_{\rm inc}^{\rm opt} -1$', **errorp)

    ax1.errorbar(np.array([1]), np.array(A3['opt'][1]), yerr=np.array(A3_err['opt'][1]), marker='o', mfc=blue, color=blue, **errorb)

    ax1.fill_between(x = [-1, 4], y1 = [y1a, y1a], y2 = [y2a, y2a], alpha = 0.3, color=grape)
    ax1.fill_between(x = [-1, 4], y1 = A3['opt'][1]+A3_err['opt'][1], y2 = A3['opt'][1]-A3_err['opt'][1], alpha = 0.3, color=blue)

    legend_without_duplicate_labels(ax1, 'lower center', 3)

    ax2.set_ylabel(q_label, **fs_p)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange(0 - 0.5, 4 + 0.5, 1), 0.1 * np.ones([5]), 'r--')

    for i in range(4):
        ax2.scatter(np.array([i]), np.array(Q['opt'][i]), marker='o', c='None', edgecolors=blue)
        ax2.scatter(np.array([i-0.2]), np.array(Q['opt+1'][i]), marker='^', c='None', edgecolors=green)
        ax2.scatter(np.array([i+0.2]), np.array(Q['opt-1'][i]), marker='P', c='None', edgecolors=peach)

    ax2.scatter(np.array([1]), np.array(Q['opt'][1]), marker='o', c=blue, edgecolors=blue)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xticks([0, 1, 2, 3], [r'$1$', r'$2$', r'$3$', r'$4$'])
    plt.xlabel(r'$\rm nstates$', **fs_p)
    plt.xlim([-0.5, 3.5])
    ax2.tick_params(axis='both', direction='in', which='major', **ls_p)
    fig.savefig("./new_plots/ga_23_late_tsep_101214_taucut.pdf", transparent=True)

def tmin_plot(n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, fit_name, xlabel):
    value={}
    value['Q']=[]
    value['logGBF']=[]
    value['gA']=[]
    value['gA_err']=[]
    x = []

    for n_ in range(n_range[1]): # n_ represents index of n in n_range
        value['Q'].append([])
        value['logGBF'].append([])
        value['gA'].append([])
        value['gA_err'].append([])
        x.append([])

    for n in range(n_range[0], n_range[1]): # n represents nstates value
        for situation in situation_list:
            nstate_dict = {}
            nstate_dict['2pt'] = situation.pt2_nstates
            nstate_dict['3pt'] = situation.pt3_nstates
            nstate_dict['sum'] = situation.sum_nstates

            if nstate_dict[nstate_name] == n:
                tmin_dict = {}
                tmin_dict['2pt'] = situation.pt2_tmin
                tmin_dict['3pt_gA'] = situation.pt3_A3_tsep_min
                tmin_dict['sum_gA'] = situation.sum_A3_tsep_min

                x[n].append(tmin_dict[tmin_name])  # here is the varying parameter

        x[n].sort() # fix n, search for all situations to complete x[n], then sort x[n]

        for i in range(len(x[n])):
            for situation in situation_list:
                nstate_dict = {}
                nstate_dict['2pt'] = situation.pt2_nstates
                nstate_dict['3pt'] = situation.pt3_nstates
                nstate_dict['sum'] = situation.sum_nstates

                tmin_dict = {}
                tmin_dict['2pt'] = situation.pt2_tmin
                tmin_dict['3pt_gA'] = situation.pt3_A3_tsep_min
                tmin_dict['sum_gA'] = situation.sum_A3_tsep_min

                if nstate_dict[nstate_name] == n and tmin_dict[tmin_name] == x[n][i]:
                    value['Q'][n].append(situation.Q_value)
                    value['logGBF'][n].append(situation.log_GBF)
                    value['gA'][n].append(situation.A300)
                    value['gA_err'][n].append(situation.A300_err)


    print(x)

    best_n_ = best_n - n_range[0] # n_ represents index, n represents value
    best_t_ = best_t - t_range[0] # t_ represents index, t represents value

    #####################################################################################
    # gA - tmin
    #fig=plt.figure()

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)

    # ax1
    ax1.set_ylabel(ga_label, **gA_fs_p)
    ax1.set_ylim(gA_ylim)


    for n in range(n_range[0], n_range[1]):
        ax1.errorbar(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['gA'][n]), yerr=np.array(value['gA_err'][n]), marker='o', color=color_list[n-2], label = r'$n_{\rm s} = %d$' %n, **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), (value['gA'][best_n][best_t_]+value['gA_err'][best_n][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['gA'][best_n][best_t - t_range[0]]-value['gA_err'][best_n][best_t - t_range[0]])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n - 2], alpha=0.2)

    # best fit
    ax1.errorbar(np.array([best_t + (best_n-4)*0.1]), np.array([value['gA'][best_n][best_t_]]), yerr=np.array([value['gA_err'][best_n][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)
    print("Best fit")
    print(np.array([value['gA'][best_n][best_t_]]), np.array([value['gA_err'][best_n][best_t_]]))


    legend_without_duplicate_labels(ax1, 'lower center', 4)

    # ax2
    ax2.set_ylabel(q_label, **fs_p)
    ax2.set_ylim([0.01, 1.1])
    ax2.set_yticks([0.1,0.75])

    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        ax2.scatter(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['Q'][n]), marker='o', c='None', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([value['Q'][best_n][best_t_]]), marker='o', c=color_list[best_n-2])


    # ax3
    ax3.set_ylabel(w_label, **fs_p)
    ax3.set_ylim([0.01, 1.1])

    #ax3.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.3 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    log_max = {}

    for t in range(t_range[0], t_range[1]):
        t_ = t - t_range[0]
        logGBF_list = []
        for n in range(n_range[0], n_range[1]):
            logGBF_list.append(value['logGBF'][n][t_])

        log_max['t='+str(t)] = max(logGBF_list)
        w_norm = np.sum( np.array( [np.exp(lgbf-max(logGBF_list)) for lgbf in logGBF_list]))
        print(w_norm, [np.exp(lgbf-max(logGBF_list)) for lgbf in logGBF_list])
        if t == best_t:
            w_norm_best = w_norm

        for n in range(n_range[0], n_range[1]):
            #w = np.exp(value['logGBF'][n][t_] - log_max['t='+str(t)])
            w = np.exp(value['logGBF'][n][t_] - log_max['t='+str(t)]) / w_norm
            ax3.scatter(np.array([t + (n-4)*0.1]), np.array([w]), marker='o', c='None', edgecolors=color_list[n-2])

    # best fit
    print(value['logGBF'][best_n][best_t_])
    ax3.scatter(np.array([best_t + (best_n-4)*0.1]),
        #np.array([np.exp( value['logGBF'][best_n][best_t_] - log_max['t='+str(best_t)] )]),
        np.array([np.exp( value['logGBF'][best_n][best_t_] - log_max['t='+str(best_t)]) / w_norm_best]),
        marker='o', c=color_list[best_n-2])
    ax3.set_yticks([0.1,0.75])

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **fs_p)
    plt.xlim([t_range[0] - 0.5, t_range[1] - 0.5])
    ax1.tick_params(axis='both', direction='in', which='major', **ls_p)
    ax2.tick_params(axis='both', direction='in', which='major', **ls_p)
    ax3.tick_params(axis='both', direction='in', which='major', **ls_p)

    plt.savefig('./new_plots/'+fit_name+'_gA-'+tmin_name+'_tmin.pdf', transparent=True)

def tmax_plot(t_range, best_n, best_t, tmax_name, situation_list, fit_name, xlabel):

    value={}
    value['Q']=[]
    value['gA']=[]
    value['gA_err']=[]

    x=[]

    for situation in situation_list:
        tmax_dict = {}
        tmax_dict['2pt'] = situation.pt2_tmax
        tmax_dict['3pt'] = situation.pt3_A3_tsep_max
        tmax_dict['sum'] = situation.sum_A3_tsep_max

        x.append(tmax_dict[tmax_name])  # here is the varying parameter

    x.sort()
    for i in range(len(x)):
        for situation in situation_list:
            tmax_dict = {}
            tmax_dict['2pt'] = situation.pt2_tmax
            tmax_dict['3pt'] = situation.pt3_A3_tsep_max
            tmax_dict['sum'] = situation.sum_A3_tsep_max
            if tmax_dict[tmax_name] == x[i]:
                value['Q'].append(situation.Q_value)
                value['gA'].append(situation.A300)
                value['gA_err'].append(situation.A300_err)

    print(x)

    best_t_ = best_t - t_range[0]

    #####################################################################################
    #####################################################################################

    # gA - tmax
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)
    ax1.set_ylabel(ga_label, **gA_fs_p)
    ax1.set_ylim(gA_ylim)

    ax1.tick_params(direction='in', **ls_p)

    print('DEBUG:')
    print(np.array(x))
    print(np.array(value['gA']))
    ax1.errorbar(np.array(x)-1, np.array(value['gA']), yerr=np.array(value['gA_err']), marker='o', color=color_list[best_n-2], **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1)-1, (value['gA'][best_t_]+value['gA_err'][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['gA'][best_t_]-value['gA_err'][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n-2], alpha=0.2)
    # best fit
    ax1.errorbar(np.array([best_t])-1, np.array([value['gA'][best_t_]]), yerr=np.array([value['gA_err'][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)


    ax2.set_ylabel(q_label, **fs_p)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1)-1, 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    ax2.scatter(np.array(x)-1, np.array(value['Q']), marker='o', c='None', edgecolors=color_list[best_n-2])

    # best fit
    ax2.scatter(np.array([best_t])-1, np.array([value['Q'][best_t_]]), marker='o', c=color_list[best_n-2])

    ax2.tick_params(direction='in', **ls_p)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **fs_p)
    plt.xlim([t_range[0] - 1.5, t_range[1] - 1.5])

    plt.savefig('./new_plots/'+fit_name+'_gA-'+tmax_name+'_tmax.pdf', transparent=True)

def best_fit_comparison_plot(A3, A3_err, Q):
    n_ga = [1.266, 0.011]
    y1a = n_ga[0] - n_ga[1]
    y2a = n_ga[0] + n_ga[1]

    #####################gA
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)
    ax1.set_ylabel(ga_label, **gA_fs_p)
    ax1.set_ylim(gA_ylim)
    ax1.tick_params(axis='both', which='major', direction='in', **ls_p)

    for i in range(3):
        ax1.errorbar(np.array([i]), np.array(A3[i]), yerr=np.array(A3_err[i]), marker='o', color="k", **errorp)

    ax1.fill_between(x = [-1, 3], y1 = [y1a, y1a], y2 = [y2a, y2a], alpha = 0.2, color="k")

    ax2.set_ylabel(q_label, **fs_p)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange(0 - 0.5, 3 + 0.5, 1), 0.1 * np.ones([4]), 'r--')

    for i in range(3):
        ax2.scatter(np.array([i]), np.array(Q[i]), marker='o', c='None', edgecolors="k")

    plt.subplots_adjust(wspace=0, hspace=0)

    plt.xticks([0, 1, 2], [r'$23\rm s$', r'$2\rm s$', r'$23$'], **fs_p)
    plt.xlim([-0.5, 2.5])
    ax2.tick_params(axis='both', which='major', direction='in', **ls_p)
    plt.savefig("./new_plots/ga_summary.pdf", transparent=True)

def tmin_div_plot(n_range, t_range, best_n, best_t, tmin_name, situation_list, fit_name, xlabel):
    value={}
    value['Q']=[]
    value['gA']=[]
    value['gA_err']=[]
    x = []

    for n_ in range(n_range[1]):
        value['Q'].append([])
        value['gA'].append([])
        value['gA_err'].append([])
        x.append([])


    for n in range(n_range[0], n_range[1]):
        n_ = n - n_range[0]
        for t in range(t_range[n_][0], t_range[n_][1]):
            for situation in situation_list[n_]:
                tmin_dict = {}
                tmin_dict['2pt'] = situation.pt2_tmin
                tmin_dict['3pt_gA'] = situation.pt3_A3_tsep_min
                tmin_dict['sum_gA'] = situation.sum_A3_tsep_min

                if tmin_dict[tmin_name] == t:
                    value['Q'][n].append(situation.Q_value)
                    value['gA'][n].append(situation.A300)
                    value['gA_err'][n].append(situation.A300_err)

                    x[n].append(tmin_dict[tmin_name])  # here is the varying parameter

    print(x)

    best_n_ = best_n - n_range[0]
    best_t_ = []
    for n_ in range(len(best_t)):
        best_t_.append(best_t[n_] - t_range[n_][0])

    plot_tmin = t_range[n_range[1] - n_range[0] - 1][0]
    plot_tmax = t_range[0][1]

    #####################################################################################
    # gA - tmin

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmin_div)

    # ax1
    ax1.set_ylabel(ga_label, **gA_fs_p)
    ax1.set_ylim(gA_ylim)

    for n in range(n_range[0], n_range[1]):

        col = color_list[n-2]
        alp = 1
        if col == grey:
            col = 'k'
            alp = 0.7

        n_ = n - n_range[0]
        ax1.errorbar(np.arange(t_range[n_][0], t_range[n_][1]) + (n-2) * 0.2, np.array(value['gA'][n]), yerr=np.array(value['gA_err'][n]), marker=marker_list[n], color=col, label=r'$n_{\rm s} = %d$' %n, alpha=alp, **errorp)


    # best fit
    ax1.fill_between(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['gA'][best_n][best_t_[best_n_]]+value['gA_err'][best_n][best_t_[best_n_]])*np.ones([plot_tmax - plot_tmin + 1]), (value['gA'][best_n][best_t_[best_n_]]-value['gA_err'][best_n][best_t_[best_n_]])*np.ones([plot_tmax - plot_tmin + 1]), color=color_list[best_n - 2], alpha=0.2)
    for n_ in range(n_range[1] - n_range[0]):
        n = n_ + n_range[0]

        col = color_list[n-2]
        alp = 1
        if col == grey:
            col = 'k'
            alp = 0.7

        if n != best_n:
            ax1.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['gA'][n][best_t_[n_]]+value['gA_err'][n][best_t_[n_]])*np.ones([plot_tmax - plot_tmin + 1]), color = col, alpha=alp, linestyle='--')
            ax1.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['gA'][n][best_t_[n_]]-value['gA_err'][n][best_t_[n_]])*np.ones([plot_tmax - plot_tmin + 1]), color = col, alpha=alp, linestyle='--')

        ax1.errorbar(np.array([best_t[n_] + (n - 2)*0.2]), np.array([value['gA'][n][best_t_[n_]]]), yerr=np.array([value['gA_err'][n][best_t_[n_]]]), marker=marker_list[n], mfc=col, color=col, alpha=alp, **errorb)


    legend_without_duplicate_labels(ax1, 'lower center', 3)

    print("Best fit")
    print(np.array([value['gA'][best_n][best_t_[best_n_]]]), np.array([value['gA_err'][best_n][best_t_[best_n_]]]))


    # ax2
    ax2.set_ylabel(q_label, **fs_p)
    ax2.set_ylim([0, 1.1])

    ax2.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), 0.1 * np.ones([plot_tmax - plot_tmin + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):

        col = color_list[n-2]
        alp = 1
        if col == grey:
            col = 'k'
            alp = 0.7

        n_ = n - n_range[0]
        ax2.scatter(np.arange(t_range[n_][0], t_range[n_][1]) + (n-2) * 0.2, np.array(value['Q'][n]), marker=marker_list[n], c='None', alpha=alp, edgecolors=col)

    # best fit
    for n_ in range(n_range[1] - n_range[0]):
        n = n_ + n_range[0]

        col = color_list[n-2]
        alp = 1
        if col == grey:
            col = 'k'
            alp = 0.7

        ax2.scatter(np.array([best_t[n_] + (n-2)*0.2]), np.array([value['Q'][n][best_t_[n_]]]), marker=marker_list[n], c=col, alpha=alp)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **fs_p)
    plt.xlim([plot_tmin - 0.5, plot_tmax - 0.5])
    ax1.tick_params(axis='both', direction='in', which='major', **ls_p)
    ax2.tick_params(axis='both', direction='in', which='major', **ls_p)

    plt.savefig('./new_plots/'+fit_name+'_gA-'+tmin_name+'_tmin.pdf', transparent=True)

def tau_cut_plot(A3, A3_err, Q, n, fit_name):
    #####################gA
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)

    ax1.set_ylabel(ga_label, **gA_fs_p)
    ax1.set_ylim(gA_ylim)

    ax1.errorbar(np.arange(0, 5), np.array(A3), yerr=np.array(A3_err), marker='o', color=color_list[n-2], **errorp)

    ax1.fill_between(np.arange(0 - 0.5, 5 + 0.5, 1), (A3[1] + A3_err[1])*np.ones([6]), (A3[1] - A3_err[1])*np.ones([6]), color=color_list[n-2], alpha=0.2)

    if fit_name == '23s':
        best_i = 1

    elif fit_name == '23':
        best_i = 0

    # best fit
    ax1.errorbar(np.array([best_i]), np.array([A3[best_i]]), yerr=np.array([A3_err[best_i]]), marker='o', mfc=color_list[n-2], color=color_list[n-2], **errorb)


    ax2.set_ylabel(q_label, **fs_p)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange(0 - 0.5, 5 + 0.5, 1), 0.1 * np.ones([6]), 'r--')

    ax2.scatter(np.arange(0, 5), np.array(Q), marker='o', c='None', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([best_i]), np.array([Q[best_i]]), marker='o', c=color_list[n-2])

    plt.subplots_adjust(wspace=0, hspace=0)


    if fit_name == '23s':
        plt.xticks([0, 1, 2, 3, 4], [r'$\tau_{\rm inc}^{\rm opt} + 1$', r'$\tau_{\rm inc}^{\rm opt}$', r'$\tau_{\rm inc}^{\rm opt} - 1$', r'$\tau_{\rm inc}^{\rm opt} - 2$', r'$\tau_{\rm inc}^{\rm opt} - 3$'])

    elif fit_name == '23':
        plt.xticks([0, 1, 2, 3, 4], [r'$\tau_{\rm inc}^{\rm opt}$', r'$\tau_{\rm inc}^{\rm opt} - 1$', r'$\tau_{\rm inc}^{\rm opt} - 2$', r'$\tau_{\rm inc}^{\rm opt} - 3$', r'$\tau_{\rm inc}^{\rm opt} - 4$'])

    plt.xlim([-0.5, 4.5])
    ax1.tick_params(axis='both', which='major', direction='in', **ls_p)
    ax2.tick_params(axis='both', which='major', direction='in', **ls_p)

    plt.savefig("./new_plots/"+fit_name+"_gA_tins.pdf", transparent=True)

def late_23_tmin_plot(E0_ylim, n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, fit_name, xlabel):
    value={}
    value['Q']=[]
    value['logGBF']=[]
    value['E0']=[]
    value['E0_err']=[]
    x = []

    for n_ in range(n_range[1]):
        value['Q'].append([])
        value['logGBF'].append([])
        value['E0'].append([])
        value['E0_err'].append([])
        x.append([])


    for n in range(n_range[0], n_range[1]):
        for t in range(t_range[0], t_range[1]):
            for situation in situation_list:
                nstate_dict = {}
                nstate_dict['2pt'] = situation.pt2_nstates
                nstate_dict['3pt'] = situation.pt3_nstates
                nstate_dict['sum'] = situation.sum_nstates

                tmin_dict = {}
                tmin_dict['2pt'] = situation.pt2_tmin
                tmin_dict['3pt_gA'] = situation.pt3_A3_tsep_min
                tmin_dict['3pt_gV'] = situation.pt3_V4_tsep_min
                tmin_dict['sum_gA'] = situation.sum_A3_tsep_min
                tmin_dict['sum_gV'] = situation.sum_V4_tsep_min

                if nstate_dict[nstate_name] == n and tmin_dict[tmin_name] == t:
                    value['Q'][n].append(situation.Q_value)
                    value['logGBF'][n].append(situation.log_GBF)
                    value['E0'][n].append(situation.E0)
                    value['E0_err'][n].append(situation.E0_err)

                    x[n].append(tmin_dict[tmin_name])  # here is the varying parameter

    print(x)

    best_n_ = best_n - n_range[0]
    best_t_ = best_t - t_range[0]

    #####################################################################################
    # E0 - tmin
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)

    # ax1
    ax1.set_ylabel(e0_label, **fs_p)
    ax1.set_ylim(E0_ylim)

    for n in range(n_range[0], n_range[1]):
        ax1.errorbar(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['E0'][n]), yerr=np.array(value['E0_err'][n]), marker=marker_list[n], color=color_list[n-2], label=r'$n_{\rm s} = %d$' %n, **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), (value['E0'][best_n][best_t_]+value['E0_err'][best_n][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['E0'][best_n][best_t - t_range[0]]-value['E0_err'][best_n][best_t - t_range[0]])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n - 2], alpha=0.2)

    # best fit
    ax1.errorbar(np.array([best_t + (best_n-4)*0.1]), np.array([value['E0'][best_n][best_t_]]), yerr=np.array([value['E0_err'][best_n][best_t_]]), marker=marker_list[best_n], mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)

    legend_without_duplicate_labels(ax1, 'lower center', 4)

    # ax2
    ax2.set_ylabel(q_label, **fs_p)
    ax2.set_ylim([0, 1.1])
    ax2.set_yticks([0.1,0.75])

    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        ax2.scatter(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['Q'][n]), marker=marker_list[n], c='None', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([value['Q'][best_n][best_t_]]), marker=marker_list[best_n], c=color_list[best_n-2])


    # ax3
    ax3.set_ylabel(w_label, **fs_p)
    ax3.set_ylim([0, 1.3])

    #ax3.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.3 * np.ones([t_range[1] - t_range[0] + 1]), 'b--')
    #ax3.axhline(0.3, )

    log_max = {}

    for t in range(t_range[0], t_range[1]):
        t_ = t - t_range[0]
        logGBF_list = []
        for n in range(n_range[0], n_range[1]):
            logGBF_list.append(value['logGBF'][n][t_])

        #log_max['t='+str(t)] = max(logGBF_list)
        w_norm = np.sum( np.array( [np.exp(lgbf) for lgbf in logGBF_list]))
        if t == best_t:
            w_norm_best = w_norm

        w_list = []
        for n in range(n_range[0], n_range[1]):
            #w = np.exp(value['logGBF'][n][t_] - log_max['t='+str(t)])
            w = np.exp(value['logGBF'][n][t_]) / w_norm

            ax3.scatter(np.array([t + (n-4)*0.1]), np.array([w]), marker=marker_list[n], c='None', edgecolors=color_list[n-2])

    # best fit
    ax3.scatter(np.array([best_t + (best_n-4)*0.1]),
        #np.array([np.exp( value['logGBF'][best_n][best_t_] - log_max['t='+str(best_t)] )]),
        np.array([np.exp( value['logGBF'][best_n][best_t_]) / w_norm_best]),
        marker=marker_list[best_n], c=color_list[best_n-2])
    ax3.set_yticks([0.1,0.75])

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **fs_p)
    plt.xlim([t_range[0] - 0.5, t_range[1] - 0.5])
    plt.xticks([t for t in range(3, 17, 2)])
    ax1.tick_params(axis='both', direction='in', which='major', **ls_p)
    ax2.tick_params(axis='both', direction='in', which='major', **ls_p)
    ax3.tick_params(axis='both', direction='in', which='major', **ls_p)

    plt.savefig('./new_plots/'+fit_name+'_E0-'+tmin_name+'_tmin.pdf', transparent=True)

def prior_width_plot(value):
    # gA - prior width
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_prior_width)
    ax1.set_ylabel(ga_label, **gA_fs_p)
    ax1.set_ylim(gA_ylim)
    ax1.tick_params(axis='both', which='major', direction='in', **ls_p)

    mask = [0, 1, 3, 4]
    ax1.errorbar([-0.05, 0.95, 2.95, 3.95], np.array(value['A3'])[mask], yerr=np.array(value['A3_err'])[mask],label=r'$\rm{g.s.}$', marker="P", color=red, **errorp)
    ax1.errorbar([0.05, 1.05, 3.05, 4.05], np.array(value['A3_'])[mask], yerr=np.array(value['A3_err_'])[mask], label=r'$\rm{g.s.\ and\ 1\ e.s.}$', marker="^", color=blue, **errorp)

    ax1.fill_between(np.arange(- 0.5, 5.5, 1), (value['A3'][2]+value['A3_err'][2])*np.ones([6]), (value['A3'][2]-value['A3_err'][2])*np.ones([6]), color=green, alpha=0.3)

    # highest logGBF
    ax1.errorbar(np.array([2]), np.array([value['A3'][2]]), yerr=np.array([value['A3_err'][2]]), color=green, marker="o", label=r"$\rm{best\ fit}$", **errorb)

    legend_without_duplicate_labels(ax1, 'lower center', 3)

    ax2.set_ylabel(q_label, **fs_p)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange( - 0.5, 5 + 0.5, 1), 0.1 * np.ones([6]), 'r--')

    ax2.scatter([-0.05, 0.95, 2.95, 3.95], np.array(value['Q'])[mask], marker='P', c='None', edgecolors=red)
    ax2.scatter([0.05, 1.05, 3.05, 4.05], np.array(value['Q_'])[mask], marker='^', c='None', edgecolors=blue)

    # highest logGBF
    ax2.scatter(np.array([2]), np.array([value['Q'][2]]), marker="o", c=green)

    ax2.tick_params(axis='both', direction='in', which='major', **ls_p)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(r'$ \rm prior\ width$', **fs_p)
    plt.xticks([0, 1, 2, 3, 4], [r'$0.1\sigma$', r'$0.5\sigma$', r'$1\sigma$', r'$2\sigma$', r'$10\sigma$'], **fs_p)
    plt.xlim([- 0.5, 4.5])

    plt.savefig('./new_plots/'+'23s_gA-prior_width.pdf', transparent=True)
    plt.show()

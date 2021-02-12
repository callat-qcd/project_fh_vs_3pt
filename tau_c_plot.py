# %%
import h5py as h5 
import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  
import os

plt.rcParams.update({"text.usetex": True})
import matplotlib as mpl
mpl.pyplot.ion()
#%matplotlib inline

# %%
from module.prepare_data import Prepare_data
from module.fit import Fit
from module.prior_setting import prior_ho_width_1 
prior = prior_ho_width_1

plt.rcParams.update({"text.usetex": True})
fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
fig_size  = (fig_width, fig_width / gr)
plt_axes  = [0.15,0.15,0.845,0.845]
fs_text   = 18 # font size of text
fs_leg    = 16 # legend font size
tick_size = 16 # tick size
plt.rcParams['figure.figsize'] = fig_size

errorp = {"markersize": 5, "linestyle": "none", "capsize": 3, "elinewidth": 1}

grey = "#808080" # nstates = 1
red = "#FF6F6F" # nstates = 2
peach = "#FF9E6F"
orange = "#FFBC6F" # nstates = 3
sunkist = "#FFDF6F"
yellow = "#FFEE6F" # nstates = 4
lime = "#CBF169"
green = "#5CD25C" # nstates = 5
turquoise = "#4AAB89"
blue = "#508EAD" # nstates = 6
grape = "#635BB1"
violet = "#7C5AB8" # nstates = 7
fuschia = "#C3559F"

color_list = [grey, red, orange, green, blue, grape]

# %%
def tau_c_sum_plot(f_range, d_range, fit_result):

    plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    for sum_tau_cut in range(f_range[0], f_range[1]):

        fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

        plot_space = 0.05 # for fill_between plot

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

        ax.fill_between(x_fill, y1, y2, color=color_list[sum_tau_cut], alpha=0.3)

    for sum_tau_cut in range(d_range[0], d_range[1]):

        data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)
        
        sum_data_list = []

        tmin = max(2, 2*sum_tau_cut)

        for t in range(tmin, 14):
            sum_data_list.append(data_avg_dict_completed['sum_A3_fit_'+str(t)])

        temp_mean = np.array([val.mean for val in sum_data_list])
        temp_sdev = np.array([val.sdev for val in sum_data_list])

        if sum_tau_cut == 1:
            ax.errorbar(np.arange(tmin, 14), temp_mean, yerr=temp_sdev, label=r'$\tau_c=%d$' %sum_tau_cut, marker='o', color=color_list[sum_tau_cut], mfc=color_list[sum_tau_cut], **errorp)

        else:
            ax.errorbar(np.arange(tmin, 14), temp_mean, yerr=temp_sdev, label=r'$\tau_c=%d$' %sum_tau_cut, marker='o', color=color_list[sum_tau_cut], mfc='none', **errorp)

    ax.set_ylim(1.0, 1.349)
    ax.set_xlim([1, 23])
    ax.set_xlabel(r'$t$', fontsize=fs_text)

    ax.tick_params(direction='in', labelsize=tick_size)

    ax.legend(loc='lower right', ncol=3, columnspacing=0, handletextpad=0.1, fontsize=fs_leg)
    plot_range = str(f_range[0]) + '_to_' + str(f_range[1]-1)
    plt.savefig(f"./new_plots/tau_c_sum_plot_{plot_range}.pdf", transparent=True)
    plt.show()

# %%
#############################################
#################### 2pt+3pt+sum best fit
#############################################
file_name = 'a09m310_e_gA_srcs0-15_full_tau.h5'
file_path = os.getcwd() + '/' + file_name # just put the data file inside the 'my_project' folder, if no database, data file should be put in the same path as this file.

pt2_data_range = [0, 96]
pt3_data_range = [2, 15]

prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

data_avg_dict = prepare_data.read_data_with_average()

# do the fit

################ set parameters ###################

pt2_t = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])
pt2_nstates = 5

pt3_A3_t = [3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14]
pt3_A3_tau = [1, 1, 2, 1, 2, 1, 2, 3, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 7]
pt3_V4_t = pt3_A3_t
pt3_V4_tau = pt3_A3_tau

pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]
pt3_nstates = pt2_nstates

sum_A3 = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
sum_V4 = sum_A3 
sum_nstates = 5
sum_tau_cut = 1

include_2pt = True
include_3pt = True
include_sum = True

save = False

data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)

best_p0 = {'E0': 0.49007432827585923, 'log(dE1)': -1.2182574657830274, 'z0': 0.0003202622719326246, 'z1': 0.00033086411929253656, 'z0_ps': 0.003, 'z1_ps': 1.4290420366321432e-21, 
                   'log(dE2)': -1.0203037715038503, 'z2': -0.0003420981842067054, 'z2_ps': 4.77757294807751e-19, 'log(dE3)': -0.6763116611503818, 'z3': 0.0006114436301814257, 'z3_ps': 6.18600096063414e-20, 
                   'log(dE4)': 0.1096154276143707, 'z4': 0.00030912500545967415, 'z4_ps': -2.1747716630250984e-21, 'log(dE5)': -1.25, 'z5': -1.5592681273181436e-22, 'z5_ps': 1.3797601616142584e-19, 
                   'log(E_fh)': -1.2512330692654803, 'z_fh_ss': 1.747909186937848e-05, 'z_fh_ps': 4.363143347633416e-21, 'A3_00': 1.2564901354721003, 'V4_00': 1.0231038334230558, 'A3_01': -0.2679609024371228, 
                   'V4_01': -0.004166047120256729, 'A3_11': 0.8761219755499188, 'V4_11': 0.9902000810575383, 'A3_02': -0.1521913122970933, 'V4_02': 0.0014475402638539437, 'A3_12': -0.06373623540219797, 
                   'V4_12': 0.01132917934978911, 'A3_22': 0.7453490951343201, 'V4_22': 1.0632482453026812, 'A3_03': 0.10049661764180362, 'V4_03': 0.08038556296566231, 'A3_13': 0.1361847606270459, 
                   'V4_13': 0.07099939052507677, 'A3_23': 0.25916215890865607, 'V4_23': -0.01133073437078887, 'A3_33': 0.865312735370937, 'V4_33': 0.7748558903412537, 'A3_04': 0.6204441395578857, 
                   'V4_04': 0.3522708502720259, 'A3_14': 0.40763003829174205, 'V4_14': 0.335117716953269, 'A3_24': -1.101075754867819, 'V4_24': -0.3733418857549809, 'A3_34': 0.7190174286711979, 
                   'V4_34': 0.6261439927183979, 'A3_44': -0.25761311272489074, 'V4_44': 0.9973793333624172, 'A3_05': -8.257694779526328e-21, 'V4_05': 4.986761175050056e-20, 'A3_15': -3.247769578767987e-21, 
                   'V4_15': 2.124892576801777e-20, 'A3_25': -5.58683860326533e-20, 'V4_25': 1.6748959838190552e-20, 'A3_35': 7.295422809984337e-21, 'V4_35': 2.214178969103922e-20, 'A3_45': 3.5577787778307475e-20, 
                   'V4_45': -2.5297406179805884e-20, 'A3_55': 1.5137185489740753e-20, 'V4_55': 1.0, 'fh_A3_0': 0.002666305722764564, 'fh_V4_0': 0.01628080958292751, 'd0_ss_A3': -1.5e-06, 'd0_ss_V4': 1.3e-06, 
                   'd0_ps_A3': -9e-06, 'd0_ps_V4': 7.5e-06, 'fh_A3_1': -0.0024127278385662476, 'fh_V4_1': 0.003973560962251952, 'd1_ss_A3': 4.059923399570421e-26, 'd1_ss_V4': -6.166269907949912e-26, 'd1_ps_A3': 2.3105448329306423e-26, 
                   'd1_ps_V4': 4.902765855922857e-26, 'fh_A3_2': -0.022933720252869036, 'fh_V4_2': -0.055108665759743665, 'd2_ss_A3': -5.052876947576615e-26, 'd2_ss_V4': -1.7883348339254526e-27, 'd2_ps_A3': 9.668146857259591e-27, 
                   'd2_ps_V4': -2.664874890850679e-25, 'fh_A3_3': -0.025212744548757944, 'fh_V4_3': -0.09513944044174277, 'd3_ss_A3': -4.139294809856933e-26, 'd3_ss_V4': 3.1366216019351247e-26, 'd3_ps_A3': 6.126928367111804e-26, 
                   'd3_ps_V4': 8.925147564672736e-27, 'fh_A3_4': -0.0007628417809514719, 'fh_V4_4': -0.0026882359647512925, 'd4_ss_A3': -1.8337248632272437e-26, 'd4_ss_V4': 1.0839295795776121e-30, 'd4_ps_A3': -8.40639490383513e-31, 
                   'd4_ps_V4': -1.8447470978983475e-30, 'd5_ss_A3': -4.996925396047227e-31, 'd5_ss_V4': 1.0374204770912712e-30, 'd5_ps_A3': -2.4636016481713015e-30, 'd5_ps_V4': -5.605193857299268e-45, 'fh_A3_5': -9.616184647013097e-06, 'fh_V4_5': 0.9999986504760016}

fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

fit_result = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0)[0]

print(fit_result)

# %%
f_range = [0, 2]
d_range = [0, 2]

tau_c_sum_plot(f_range, d_range, fit_result)

# %%

# %%
import h5py as h5 
import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  
import hashlib
import os
from scipy import interpolate

plt.rcParams.update({"text.usetex": True})
import matplotlib as mpl
mpl.pyplot.ion()
# %matplotlib inline # for jupyter

# %%
from module.prepare_data import Prepare_data
from module.fit import Fit
from module.prior_setting import prior_ho_width_1
prior = prior_ho_width_1

# %%
### standard colors

red = "#FF6F6F"
peach = "#FF9E6F"
orange = "#FFBC6F"
sunkist = "#FFDF6F"
yellow = "#FFEE6F"
lime = "#CBF169"
green = "#5CD25C"
turquoise = "#4AAB89"
blue = "#508EAD"
grape = "#635BB1"
violet = "#7C5AB8"
fuschia = "#C3559F"

# style

fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
figsize  = (fig_width, fig_width / gr)
gridspec_tmax = {'height_ratios': [3, 1], 'left': 0.12, 'right': 0.99, 'bottom': 0.15, 'top': 0.98}
q_label = r"$Q$"
oa00_label = r"$\mathring{g}_A$"
ov00_label = r"$O^V_{00}$"
c2pt_tmin = r"$t_{\rm sep}^{\rm min}:C_2$"
gA_ylim=[1.01, 1.349]
gV_ylim=[1.0051, 1.03]
textp = {"fontsize": 18}
labelp = {"labelsize": 18}
errorp = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1}
errorb = {"markersize": 5, "linestyle": "none", "capsize": 3, "elinewidth": 1}
aspect=[0.12,0.15,0.87,0.83]

plt.rcParams['figure.figsize'] = figsize

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), loc='lower center', ncol=3, fontsize=14)


# %%
### 2pt+3pt [10, 12, 14]  # tau cut + 1
# stability plot of late tsep [10, 12, 14] # varying nstates and tau_cut
file_name = 'a09m310_e_gA_srcs0-15.h5'
file_path = os.getcwd() + '/' + file_name # just put the data file inside the 'my_project' folder, if no database, data file should be put in the same path as this file.

pt2_data_range = [0, 96]
pt3_data_range = [2, 15]

prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

data_avg_dict = prepare_data.read_data_with_average()

# do the fit

#############################################################
################ set parameters #############################

pt3_A3_t = []
pt3_A3_tau = []
pt3_V4_t = []
pt3_V4_tau = []

pt3_A3_t = [10, 12, 14] # [10, 12, 12, 14, 14]
pt3_A3_tau = [5, 6, 7] # [5, 5, 6, 6, 7]
pt3_V4_t = pt3_A3_t
pt3_V4_tau = pt3_A3_tau

pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]

sum_A3 = np.array([])
sum_V4 = sum_A3 
sum_nstates = 0
sum_tau_cut = 1

use_p0 = False
include_2pt = True
include_3pt = True
include_sum = False

# for nstates = [1,5]
pt2_t_list = [
    [],
    [12, 13, 14, 15, 16, 17],
    [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
] 

##############################################################
##############################################################

A3_p1 = []
A3_err_p1 = []
V4_p1 = []
V4_err_p1 = []
Q_p1 = []

for nstates in range(1, 5):
    pt2_nstates = nstates
    pt3_nstates = pt2_nstates
    pt2_t = np.array(pt2_t_list[nstates])

    data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)


    if use_p0 == False:
        best_p0 = 0

    best_p0 = {'E0': 0.4926861952896794, 'log(dE1)': -1.273594898116278, 'z0': 0.0003260892031046884, 'z1': 0.000263341199546721, 'z0_ps': 0.002999999999999449, 'z1_ps': 8.742252813171513e-17,
               'log(dE2)': -1.2923954972884248, 'z2': -0.00041038719837353235, 'z2_ps': 1.3527493699158273e-19, 'log(dE3)': -1.2363538425685927, 'z3': -0.00021858370411089373, 'z3_ps': -2.0481201222876496e-19,
               'log(dE4)': -0.8211404557969483, 'z4': 0.00022859544792623758, 'z4_ps': -2.269364858601806e-19, 'log(E_fh)': -1.2500000000000864, 'z_fh_ss': -2.1292822629149773e-20, 'z_fh_ps': 1.7147835093259457e-19,
               'A3_00': 1.2661290957780242, 'V4_00': 1.0264985362386028, 'A3_01': -0.3516778615231407, 'V4_01': -0.04468534689899525, 'A3_11': 0.5133635952637498, 'V4_11': 1.0495939630022617, 'A3_02': -0.16492505182541714,
               'V4_02': -0.10206490798498477, 'A3_12': -0.42836059937108667, 'V4_12': -0.03879580934696547, 'A3_22': 0.490651712437167, 'V4_22': 0.9461073860857595, 'A3_03': 0.3915796384039534, 'V4_03': 0.40932562996043054,
               'A3_13': 0.1994825622768483, 'V4_13': 0.005796290661118874, 'A3_23': 0.3126615943412394, 'V4_23': -0.3775397962601796, 'A3_33': 0.09835455689160501, 'V4_33': 0.9985879807579602, 'A3_04': 1.0032935731651338,
               'V4_04': 0.8348906017782508, 'A3_14': 0.01354959867275492, 'V4_14': 0.21364847516473137, 'A3_24': -0.10285387912945705, 'V4_24': -0.243947645810037, 'A3_34': -0.07604449549348145, 'V4_34': -0.04355265360537515,
               'A3_44': 0.018277411532817922, 'V4_44': 1.0002794679730667, 'fh_A3_0': 0.0, 'fh_V4_0': 0.0, 'd0_ss_A3': -1.5e-06, 'd0_ss_V4': 1.3e-06, 'd0_ps_A3': -9e-06, 'd0_ps_V4': 7.5e-06, 'fh_A3_1': 0.0, 'fh_V4_1': 0.0,
               'd1_ss_A3': 0.0, 'd1_ss_V4': 0.0, 'd1_ps_A3': 0.0, 'd1_ps_V4': 0.0, 'fh_A3_2': 0.0, 'fh_V4_2': 0.0, 'd2_ss_A3': 0.0, 'd2_ss_V4': 0.0, 'd2_ps_A3': 0.0, 'd2_ps_V4': 0.0, 'fh_A3_3': 0.0, 'fh_V4_3': 0.0,
               'd3_ss_A3': 0.0, 'd3_ss_V4': 0.0, 'd3_ps_A3': 0.0, 'd3_ps_V4': 0.0, 'd4_ss_A3': 0.0, 'd4_ss_V4': 0.0, 'd4_ps_A3': 0.0, 'd4_ps_V4': 0.0, 'fh_A3_4': 0.0, 'fh_V4_4': 1.0}

    fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

    fit_result = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0)[0]

    best_p0 = {key: fit_result.p[key].mean for key in fit_result.p}   

    print(fit_result.format(100)) 
    
    A3_p1.append(fit_result.p['A3_00'].mean)
    A3_err_p1.append(fit_result.p['A3_00'].sdev)
    V4_p1.append(fit_result.p['V4_00'].mean)
    V4_err_p1.append(fit_result.p['V4_00'].sdev)
    Q_p1.append(fit_result.Q)

# %%
### 2pt+3pt [10, 12, 14] # tau cut
# stability plot of late tsep [10, 12, 14] # varying nstates and tau_cut
file_name = 'a09m310_e_gA_srcs0-15.h5'
file_path = os.getcwd() + '/' + file_name # just put the data file inside the 'my_project' folder, if no database, data file should be put in the same path as this file.

pt2_data_range = [0, 96]
pt3_data_range = [2, 15]

prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

data_avg_dict = prepare_data.read_data_with_average()

# do the fit

#############################################################
################ set parameters #############################

pt3_A3_t = []
pt3_A3_tau = []
pt3_V4_t = []
pt3_V4_tau = []

pt3_A3_t = [10, 12, 12, 14, 14] # [10, 12, 12, 14, 14]
pt3_A3_tau = [5, 5, 6, 6, 7] # [5, 5, 6, 6, 7]
pt3_V4_t = pt3_A3_t
pt3_V4_tau = pt3_A3_tau

pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]

sum_A3 = np.array([])
sum_V4 = sum_A3 
sum_nstates = 0
sum_tau_cut = 1

use_p0 = False
include_2pt = True
include_3pt = True
include_sum = False

# for nstates = [1,5]
pt2_t_list = [
    [],
    [12, 13, 14, 15, 16, 17],
    [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
] 

##############################################################
##############################################################

A3 = []
A3_err = []
V4 = []
V4_err = []
Q = []

for nstates in range(1, 5):
    pt2_nstates = nstates
    pt3_nstates = pt2_nstates
    pt2_t = np.array(pt2_t_list[nstates])

    data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)


    if use_p0 == False:
        best_p0 = 0

    best_p0 = {'E0': 0.4926861952896794, 'log(dE1)': -1.273594898116278, 'z0': 0.0003260892031046884, 'z1': 0.000263341199546721, 'z0_ps': 0.002999999999999449, 'z1_ps': 8.742252813171513e-17,
               'log(dE2)': -1.2923954972884248, 'z2': -0.00041038719837353235, 'z2_ps': 1.3527493699158273e-19, 'log(dE3)': -1.2363538425685927, 'z3': -0.00021858370411089373, 'z3_ps': -2.0481201222876496e-19,
               'log(dE4)': -0.8211404557969483, 'z4': 0.00022859544792623758, 'z4_ps': -2.269364858601806e-19, 'log(E_fh)': -1.2500000000000864, 'z_fh_ss': -2.1292822629149773e-20, 'z_fh_ps': 1.7147835093259457e-19,
               'A3_00': 1.2661290957780242, 'V4_00': 1.0264985362386028, 'A3_01': -0.3516778615231407, 'V4_01': -0.04468534689899525, 'A3_11': 0.5133635952637498, 'V4_11': 1.0495939630022617, 'A3_02': -0.16492505182541714,
               'V4_02': -0.10206490798498477, 'A3_12': -0.42836059937108667, 'V4_12': -0.03879580934696547, 'A3_22': 0.490651712437167, 'V4_22': 0.9461073860857595, 'A3_03': 0.3915796384039534, 'V4_03': 0.40932562996043054,
               'A3_13': 0.1994825622768483, 'V4_13': 0.005796290661118874, 'A3_23': 0.3126615943412394, 'V4_23': -0.3775397962601796, 'A3_33': 0.09835455689160501, 'V4_33': 0.9985879807579602, 'A3_04': 1.0032935731651338,
               'V4_04': 0.8348906017782508, 'A3_14': 0.01354959867275492, 'V4_14': 0.21364847516473137, 'A3_24': -0.10285387912945705, 'V4_24': -0.243947645810037, 'A3_34': -0.07604449549348145, 'V4_34': -0.04355265360537515,
               'A3_44': 0.018277411532817922, 'V4_44': 1.0002794679730667, 'fh_A3_0': 0.0, 'fh_V4_0': 0.0, 'd0_ss_A3': -1.5e-06, 'd0_ss_V4': 1.3e-06, 'd0_ps_A3': -9e-06, 'd0_ps_V4': 7.5e-06, 'fh_A3_1': 0.0, 'fh_V4_1': 0.0,
               'd1_ss_A3': 0.0, 'd1_ss_V4': 0.0, 'd1_ps_A3': 0.0, 'd1_ps_V4': 0.0, 'fh_A3_2': 0.0, 'fh_V4_2': 0.0, 'd2_ss_A3': 0.0, 'd2_ss_V4': 0.0, 'd2_ps_A3': 0.0, 'd2_ps_V4': 0.0, 'fh_A3_3': 0.0, 'fh_V4_3': 0.0,
               'd3_ss_A3': 0.0, 'd3_ss_V4': 0.0, 'd3_ps_A3': 0.0, 'd3_ps_V4': 0.0, 'd4_ss_A3': 0.0, 'd4_ss_V4': 0.0, 'd4_ps_A3': 0.0, 'd4_ps_V4': 0.0, 'fh_A3_4': 0.0, 'fh_V4_4': 1.0}

    fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

    fit_result = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0)[0]

    best_p0 = {key: fit_result.p[key].mean for key in fit_result.p}   

    print(fit_result.format(100)) 
    
    A3.append(fit_result.p['A3_00'].mean)
    A3_err.append(fit_result.p['A3_00'].sdev)
    V4.append(fit_result.p['V4_00'].mean)
    V4_err.append(fit_result.p['V4_00'].sdev)
    Q.append(fit_result.Q)
    
# %%
### 2pt+3pt [10, 12, 14] # tau cut - 1
# stability plot of late tsep [10, 12, 14] # varying nstates and tau_cut
file_name = 'a09m310_e_gA_srcs0-15.h5'
file_path = os.getcwd() + '/' + file_name # just put the data file inside the 'my_project' folder, if no database, data file should be put in the same path as this file.

pt2_data_range = [0, 96]
pt3_data_range = [2, 15]

prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

data_avg_dict = prepare_data.read_data_with_average()

# do the fit

#############################################################
################ set parameters #############################

pt3_A3_t = []
pt3_A3_tau = []
pt3_V4_t = []
pt3_V4_tau = []

pt3_A3_t = [10, 10, 12, 12, 12, 14, 14, 14] # [10, 12, 12, 14, 14]
pt3_A3_tau = [4, 5, 4, 5, 6, 5, 6, 7] # [5, 5, 6, 6, 7]
pt3_V4_t = pt3_A3_t
pt3_V4_tau = pt3_A3_tau

pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]

sum_A3 = np.array([])
sum_V4 = sum_A3 
sum_nstates = 0
sum_tau_cut = 1

use_p0 = False
include_2pt = True
include_3pt = True
include_sum = False

# for nstates = [1,5]
pt2_t_list = [
    [],
    [12, 13, 14, 15, 16, 17],
    [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
] 

##############################################################
##############################################################

A3_m1 = []
A3_err_m1 = []
V4_m1 = []
V4_err_m1 = []
Q_m1 = []

for nstates in range(1, 5):
    pt2_nstates = nstates
    pt3_nstates = pt2_nstates
    pt2_t = np.array(pt2_t_list[nstates])

    data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)


    if use_p0 == False:
        best_p0 = 0

    best_p0 = {'E0': 0.4926861952896794, 'log(dE1)': -1.273594898116278, 'z0': 0.0003260892031046884, 'z1': 0.000263341199546721, 'z0_ps': 0.002999999999999449, 'z1_ps': 8.742252813171513e-17,
               'log(dE2)': -1.2923954972884248, 'z2': -0.00041038719837353235, 'z2_ps': 1.3527493699158273e-19, 'log(dE3)': -1.2363538425685927, 'z3': -0.00021858370411089373, 'z3_ps': -2.0481201222876496e-19,
               'log(dE4)': -0.8211404557969483, 'z4': 0.00022859544792623758, 'z4_ps': -2.269364858601806e-19, 'log(E_fh)': -1.2500000000000864, 'z_fh_ss': -2.1292822629149773e-20, 'z_fh_ps': 1.7147835093259457e-19,
               'A3_00': 1.2661290957780242, 'V4_00': 1.0264985362386028, 'A3_01': -0.3516778615231407, 'V4_01': -0.04468534689899525, 'A3_11': 0.5133635952637498, 'V4_11': 1.0495939630022617, 'A3_02': -0.16492505182541714,
               'V4_02': -0.10206490798498477, 'A3_12': -0.42836059937108667, 'V4_12': -0.03879580934696547, 'A3_22': 0.490651712437167, 'V4_22': 0.9461073860857595, 'A3_03': 0.3915796384039534, 'V4_03': 0.40932562996043054,
               'A3_13': 0.1994825622768483, 'V4_13': 0.005796290661118874, 'A3_23': 0.3126615943412394, 'V4_23': -0.3775397962601796, 'A3_33': 0.09835455689160501, 'V4_33': 0.9985879807579602, 'A3_04': 1.0032935731651338,
               'V4_04': 0.8348906017782508, 'A3_14': 0.01354959867275492, 'V4_14': 0.21364847516473137, 'A3_24': -0.10285387912945705, 'V4_24': -0.243947645810037, 'A3_34': -0.07604449549348145, 'V4_34': -0.04355265360537515,
               'A3_44': 0.018277411532817922, 'V4_44': 1.0002794679730667, 'fh_A3_0': 0.0, 'fh_V4_0': 0.0, 'd0_ss_A3': -1.5e-06, 'd0_ss_V4': 1.3e-06, 'd0_ps_A3': -9e-06, 'd0_ps_V4': 7.5e-06, 'fh_A3_1': 0.0, 'fh_V4_1': 0.0,
               'd1_ss_A3': 0.0, 'd1_ss_V4': 0.0, 'd1_ps_A3': 0.0, 'd1_ps_V4': 0.0, 'fh_A3_2': 0.0, 'fh_V4_2': 0.0, 'd2_ss_A3': 0.0, 'd2_ss_V4': 0.0, 'd2_ps_A3': 0.0, 'd2_ps_V4': 0.0, 'fh_A3_3': 0.0, 'fh_V4_3': 0.0,
               'd3_ss_A3': 0.0, 'd3_ss_V4': 0.0, 'd3_ps_A3': 0.0, 'd3_ps_V4': 0.0, 'd4_ss_A3': 0.0, 'd4_ss_V4': 0.0, 'd4_ps_A3': 0.0, 'd4_ps_V4': 0.0, 'fh_A3_4': 0.0, 'fh_V4_4': 1.0}

    fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

    fit_result = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0)[0]

    best_p0 = {key: fit_result.p[key].mean for key in fit_result.p}   

    print(fit_result.format(100)) 
    
    A3_m1.append(fit_result.p['A3_00'].mean)
    A3_err_m1.append(fit_result.p['A3_00'].sdev)
    V4_m1.append(fit_result.p['V4_00'].mean)
    V4_err_m1.append(fit_result.p['V4_00'].sdev)
    Q_m1.append(fit_result.Q)
    
# %%
# stability plot of late tsep [10, 12, 14] # varying nstates and tau_cut
n_ga = [1.25343903, 0.01930155]
n_gv = [1.02243618, 0.00086969]
y1a = n_ga[0] - n_ga[1]
y2a = n_ga[0] + n_ga[1]
y1v = n_gv[0] - n_gv[1]
y2v = n_gv[0] + n_gv[1]    
    
# nstate = 1, 2, 3, 4    

#####################gA
fig=plt.figure()
plt.rcParams['figure.figsize'] = figsize
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)      
ax1.set_ylabel(oa00_label, fontsize=20)
ax1.set_ylim(gA_ylim)
ax1.tick_params(axis='both', which='major', **labelp)

for i in range(4):
    ax1.errorbar(np.array([i-0.2]), np.array(A3_m1[i]), yerr=np.array(A3_err_m1[i]), marker='^', color=green, **errorp, label=r'$\tau_c = \tau_c^{\rm opt} -1$')
    ax1.errorbar(np.array([i]), np.array(A3[i]), yerr=np.array(A3_err[i]), marker='o', color=blue, **errorp, label=r'$\tau_c = \tau_c^{\rm opt}$')
    ax1.errorbar(np.array([i+0.2]), np.array(A3_p1[i]), yerr=np.array(A3_err_p1[i]), marker='x', color=peach, **errorp, label=r'$\tau_c = \tau_c^{\rm opt} +1$')
    
ax1.errorbar(np.array([1]), np.array(A3[1]), yerr=np.array(A3_err[1]), marker='o', mfc=blue, color=blue, **errorb)
    
ax1.fill_between(x = [-1, 4], y1 = [y1a, y1a], y2 = [y2a, y2a], alpha = 0.3, color=yellow)
ax1.fill_between(x = [-1, 4], y1 = A3[1]+A3_err[1], y2 = A3[1]-A3_err[1], alpha = 0.3, color=blue)

legend_without_duplicate_labels(ax1)

print( A3[1])
print( A3_err[1])

ax2.set_ylabel(q_label, **textp)
ax2.set_ylim([0, 1.1])
ax2.plot(np.arange(0 - 0.5, 4 + 0.5, 1), 0.1 * np.ones([5]), 'r--')

for i in range(4):
    ax2.scatter(np.array([i]), np.array(Q[i]), marker='o', c='', edgecolors=blue)
    ax2.scatter(np.array([i-0.2]), np.array(Q_m1[i]), marker='o', c='', edgecolors=green)
    ax2.scatter(np.array([i+0.2]), np.array(Q_p1[i]), marker='o', c='', edgecolors=peach)
    
ax2.scatter(np.array([1]), np.array(Q[1]), marker='o', c=blue, edgecolors=blue)

plt.subplots_adjust(wspace=0, hspace=0)
plt.xticks([0, 1, 2, 3], [r'$1$', r'$2$', r'$3$', r'$4$'])
plt.xlabel(r'$\rm nstates$', **textp)
plt.xlim([-0.5, 3.5])
ax2.tick_params(axis='both', which='major', **labelp)
plt.tight_layout(pad=30, rect=aspect)
fig.savefig("./new_plots/ga_23_late_tsep_101214_taucut.pdf", transparent=True)

######################gV
plt.rcParams['figure.figsize'] = figsize
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)     
ax1.set_ylabel(ov00_label, fontsize=20)
ax1.set_ylim(gV_ylim)
ax1.tick_params(axis='both', which='major', **labelp)

for i in range(4):
    ax1.errorbar(np.array([i]), np.array(V4[i]), yerr=np.array(V4_err[i]), marker='o', color=blue, **errorp, label='tau cut')
    ax1.errorbar(np.array([i-0.2]), np.array(V4_m1[i]), yerr=np.array(V4_err_m1[i]), marker='^', color=green, **errorp, label='tau cut-1')
    ax1.errorbar(np.array([i+0.2]), np.array(V4_p1[i]), yerr=np.array(V4_err_p1[i]), marker='x', color=peach, **errorp, label='tau cut+1')
    
ax1.errorbar(np.array([1]), np.array(V4[1]), yerr=np.array(V4_err[1]), marker='o', mfc=blue, color=blue, **errorb)
    
ax1.fill_between(x = [-1, 4], y1 = [y1v, y1v], y2 = [y2v, y2v], alpha = 0.3, color=yellow)
ax1.fill_between(x = [-1, 4], y1 = V4[1]+V4_err[1], y2 = V4[1]-V4_err[1], alpha = 0.3, color=blue)

print(V4[1])
print(V4_err[1])
print(Q[1])
    
ax2.set_ylabel(q_label, **textp)
ax2.set_ylim([0, 1.1])
ax2.plot(np.arange(0 - 0.5, 4 + 0.5, 1), 0.1 * np.ones([5]), 'r--')

for i in range(4):
    ax2.scatter(np.array([i]), np.array(Q[i]), marker='o', c='', edgecolors=blue)
    ax2.scatter(np.array([i-0.2]), np.array(Q_m1[i]), marker='o', c='', edgecolors=green)
    ax2.scatter(np.array([i+0.2]), np.array(Q_p1[i]), marker='o', c='', edgecolors=peach)
    
ax2.scatter(np.array([1]), np.array(Q[1]), marker='o', c=blue, edgecolors=blue)

plt.subplots_adjust(wspace=0, hspace=0)
plt.xticks([0, 1, 2, 3], [r'$1\ state$', r'$2\ states$', r'$3\ states$', r'$4\ states$'])
plt.xlim([-0.5, 3.5])
ax2.tick_params(axis='both', which='major', **labelp)
plt.tight_layout(pad=30, rect=aspect)
fig.savefig("./new_plots/gv_23_late_tsep_101214_taucut.pdf", transparent=True)


# %%
### 2pt+3pt [10, 12, 14]
# stability plot of late tsep [10, 12, 14] # varying 2pt_tmin

#############################################################
################ set parameters #############################

pt3_A3_t = []
pt3_A3_tau = []
pt3_V4_t = []
pt3_V4_tau = []

pt3_A3_t = [10, 12, 12, 14, 14] # [10, 12, 12, 14, 14]
pt3_A3_tau = [5, 5, 6, 6, 7] # [5, 5, 6, 6, 7]
pt3_V4_t = pt3_A3_t
pt3_V4_tau = pt3_A3_tau

pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]

pt2_nstates = 2
pt3_nstates = pt2_nstates

sum_A3 = np.array([])
sum_V4 = sum_A3 
sum_nstates = 0
sum_tau_cut = 1

use_p0 = False
include_2pt = True
include_3pt = True
include_sum = False

# for tmin = [6, 12]
pt2_t_list = [
    [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    [8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    [9, 10, 11, 12, 13, 14, 15, 16, 17],
    [10, 11, 12, 13, 14, 15, 16, 17],
    [11, 12, 13, 14, 15, 16, 17]
] 

##############################################################
##############################################################

tA3 = []
tA3_err = []
tV4 = []
tV4_err = []
tQ = []

for tmin in range(6, 12):

    pt2_t = np.array(pt2_t_list[tmin-6])

    data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)


    if use_p0 == False:
        best_p0 = 0

    best_p0 = {'E0': 0.4926861952896794, 'log(dE1)': -1.273594898116278, 'z0': 0.0003260892031046884, 'z1': 0.000263341199546721, 'z0_ps': 0.002999999999999449, 'z1_ps': 8.742252813171513e-17,
               'log(dE2)': -1.2923954972884248, 'z2': -0.00041038719837353235, 'z2_ps': 1.3527493699158273e-19, 'log(dE3)': -1.2363538425685927, 'z3': -0.00021858370411089373, 'z3_ps': -2.0481201222876496e-19,
               'log(dE4)': -0.8211404557969483, 'z4': 0.00022859544792623758, 'z4_ps': -2.269364858601806e-19, 'log(E_fh)': -1.2500000000000864, 'z_fh_ss': -2.1292822629149773e-20, 'z_fh_ps': 1.7147835093259457e-19,
               'A3_00': 1.2661290957780242, 'V4_00': 1.0264985362386028, 'A3_01': -0.3516778615231407, 'V4_01': -0.04468534689899525, 'A3_11': 0.5133635952637498, 'V4_11': 1.0495939630022617, 'A3_02': -0.16492505182541714,
               'V4_02': -0.10206490798498477, 'A3_12': -0.42836059937108667, 'V4_12': -0.03879580934696547, 'A3_22': 0.490651712437167, 'V4_22': 0.9461073860857595, 'A3_03': 0.3915796384039534, 'V4_03': 0.40932562996043054,
               'A3_13': 0.1994825622768483, 'V4_13': 0.005796290661118874, 'A3_23': 0.3126615943412394, 'V4_23': -0.3775397962601796, 'A3_33': 0.09835455689160501, 'V4_33': 0.9985879807579602, 'A3_04': 1.0032935731651338,
               'V4_04': 0.8348906017782508, 'A3_14': 0.01354959867275492, 'V4_14': 0.21364847516473137, 'A3_24': -0.10285387912945705, 'V4_24': -0.243947645810037, 'A3_34': -0.07604449549348145, 'V4_34': -0.04355265360537515,
               'A3_44': 0.018277411532817922, 'V4_44': 1.0002794679730667, 'fh_A3_0': 0.0, 'fh_V4_0': 0.0, 'd0_ss_A3': -1.5e-06, 'd0_ss_V4': 1.3e-06, 'd0_ps_A3': -9e-06, 'd0_ps_V4': 7.5e-06, 'fh_A3_1': 0.0, 'fh_V4_1': 0.0,
               'd1_ss_A3': 0.0, 'd1_ss_V4': 0.0, 'd1_ps_A3': 0.0, 'd1_ps_V4': 0.0, 'fh_A3_2': 0.0, 'fh_V4_2': 0.0, 'd2_ss_A3': 0.0, 'd2_ss_V4': 0.0, 'd2_ps_A3': 0.0, 'd2_ps_V4': 0.0, 'fh_A3_3': 0.0, 'fh_V4_3': 0.0,
               'd3_ss_A3': 0.0, 'd3_ss_V4': 0.0, 'd3_ps_A3': 0.0, 'd3_ps_V4': 0.0, 'd4_ss_A3': 0.0, 'd4_ss_V4': 0.0, 'd4_ps_A3': 0.0, 'd4_ps_V4': 0.0, 'fh_A3_4': 0.0, 'fh_V4_4': 1.0}

    fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

    fit_result = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0)[0]

    best_p0 = {key: fit_result.p[key].mean for key in fit_result.p}   

    print(fit_result.format(100)) 
    
    tA3.append(fit_result.p['A3_00'].mean)
    tA3_err.append(fit_result.p['A3_00'].sdev)
    tV4.append(fit_result.p['V4_00'].mean)
    tV4_err.append(fit_result.p['V4_00'].sdev)
    tQ.append(fit_result.Q)


# %%
# stability plot of late tsep [10, 12, 14] # varying 2pt_tmin
    
# tmin = 6, 7, 8, 9, 10, 11

#####################gA
plt.rcParams['figure.figsize'] = figsize
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)     
ax1.set_ylabel(oa00_label, fontsize=20)
ax1.set_ylim(gA_ylim)
ax1.tick_params(axis='both', which='major', **labelp)

for i in range(6):
    ax1.errorbar(np.array([i+6]), np.array(tA3[i]), yerr=np.array(tA3_err[i]), marker='o', color=blue, **errorp)
    
ax1.errorbar(np.array([7]), np.array(tA3[1]), yerr=np.array(tA3_err[1]), marker='o', mfc=blue, color=blue, **errorb)
    
ax1.fill_between(x = [5, 12], y1 = [y1a, y1a], y2 = [y2a, y2a], alpha = 0.3, color=yellow) # 23s best fit
ax1.fill_between(x = [5, 12], y1 = tA3[1]+tA3_err[1], y2 = tA3[1]-tA3_err[1], alpha = 0.3, color=blue) # late tsep 23 best fit 

ax2.set_ylabel(q_label, **textp)
ax2.set_ylim([0, 1.1])
ax2.plot(np.arange(6 - 0.5, 12 + 0.5, 1), 0.1 * np.ones([7]), 'r--')

for i in range(6):
    ax2.scatter(np.array([i+6]), np.array(tQ[i]), marker='o', c='', edgecolors=blue)
    
ax2.scatter(np.array([7]), np.array(tQ[1]), marker='o', c=blue, edgecolors=blue)

plt.subplots_adjust(wspace=0, hspace=0)
#plt.xlabel('2pt tmin')
plt.xlim([5.5, 11.5])
plt.xlabel(c2pt_tmin, **textp)
ax2.tick_params(axis='both', which='major', **labelp)
plt.tight_layout(pad=30, rect=aspect)
fig.savefig("./new_plots/ga_23_late_tsep_101214_2pttmin.pdf", transparent=True)

######################gV
plt.rcParams['figure.figsize'] = figsize
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)     
ax1.set_ylabel(ov00_label, **textp)
ax1.set_ylim(gV_ylim)
ax1.tick_params(axis='both', which='major', **labelp)

for i in range(6):
    ax1.errorbar(np.array([i+6]), np.array(tV4[i]), yerr=np.array(tV4_err[i]), marker='o', color=blue, **errorp)
    
ax1.errorbar(np.array([7]), np.array(tV4[1]), yerr=np.array(tV4_err[1]), marker='o', mfc=blue, color=blue, **errorb)
    
ax1.fill_between(x = [5, 12], y1 = [y1v, y1v], y2 = [y2v, y2v], alpha = 0.3, color=yellow)
ax1.fill_between(x = [5, 12], y1 = tV4[1]+tV4_err[1], y2 = tV4[1]-tV4_err[1], alpha = 0.3, color=blue) # late tsep 23 best fit 
    
ax2.set_ylabel(q_label, **textp)
ax2.set_ylim([0, 1.1])
ax2.plot(np.arange(6 - 0.5, 12 + 0.5, 1), 0.1 * np.ones([7]), 'r--')

for i in range(6):
    ax2.scatter(np.array([i+6]), np.array(tQ[i]), marker='o', c='', edgecolors=blue)

ax2.scatter(np.array([7]), np.array(tQ[1]), marker='o', c=blue, edgecolors=blue)
    
plt.subplots_adjust(wspace=0, hspace=0)
#plt.xlabel('2pt tmin')
plt.xlim([5.5, 11.5])
ax2.tick_params(axis='both', which='major', **labelp)
plt.tight_layout(pad=30, rect=aspect)
fig.savefig("./new_plots/gv_23_late_tsep_101214_2pttmin.pdf", transparent=True)

# %%
# comparison of best fits

n_ga = [1.266, 0.011]
n_gv = [1.024, 0.001]
y1a = n_ga[0] - n_ga[1]
y2a = n_ga[0] + n_ga[1]
y1v = n_gv[0] - n_gv[1]
y2v = n_gv[0] + n_gv[1]
 
A3 = [1.25343903, 1.2602802, 1.2164349200469218]                 
A3_err = [0.01930155, 0.02283146, 0.012686761251489917]
V4 = [1.02243618, 1.02293516, 1.0139706273890496]                         
V4_err = [0.00086969, 0.00207524, 0.0044685616627316935]
Q = [0.88, 0.59, 0.8738588083167814]


c_list = [5, 2, 5]

#####################gA
plt.rcParams['figure.figsize'] = figsize
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)    
ax1.set_ylabel(oa00_label, **textp)
ax1.set_ylim(gA_ylim)
ax1.tick_params(axis='both', which='major', **labelp)

for i in range(3):
    ax1.errorbar(np.array([i]), np.array(A3[i]), yerr=np.array(A3_err[i]), marker='o', color="k", **errorp)

ax1.fill_between(x = [-1, 3], y1 = [y1a, y1a], y2 = [y2a, y2a], alpha = 0.2, color="k")

ax2.set_ylabel(q_label, **textp)
ax2.set_ylim([0, 1.1])
ax2.plot(np.arange(0 - 0.5, 3 + 0.5, 1), 0.1 * np.ones([4]), 'r--')

for i in range(3):
    ax2.scatter(np.array([i]), np.array(Q[i]), marker='o', c='', edgecolors="k")

plt.subplots_adjust(wspace=0, hspace=0)
#plt.xlabel('$tau\ cut$', **textp)
plt.xticks([0, 1, 2], [r'$23s$', r'$23$', r'$23\ late\ tsep$'])
plt.xlim([-0.5, 2.5])
ax2.tick_params(axis='both', which='major', **labelp)
plt.tight_layout(pad=30, rect=aspect)
fig.savefig("./new_plots/gv_23_late_tsep_101214_compare_ga.pdf", transparent=True)

######################gV
plt.rcParams['figure.figsize'] = figsize
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)    
ax1.set_ylabel(ov00_label, **textp)
ax1.set_ylim(gV_ylim)
ax1.tick_params(axis='both', which='major', **labelp)

for i in range(3):
    ax1.errorbar(np.array([i]), np.array(V4[i]), yerr=np.array(V4_err[i]), marker='o', color="k", **errorp)

ax1.fill_between(x = [-1, 3], y1 = [y1v, y1v], y2 = [y2v, y2v], alpha = 0.2, color="k")
    
ax2.set_ylabel(q_label, **textp)
ax2.set_ylim([0, 1.1])
ax2.plot(np.arange(0 - 0.5, 3 + 0.5, 1), 0.1 * np.ones([4]), 'r--')

for i in range(3):
    ax2.scatter(np.array([i]), np.array(Q[i]), marker='o', c='', edgecolors="k")

plt.subplots_adjust(wspace=0, hspace=0)
#plt.xlabel('$tau\ cut$', **textp)
plt.xticks([0, 1, 2], [r'$23s$', r'$23$', r'$23\ late\ tsep$'])
plt.xlim([-0.5, 2.5])
ax2.tick_params(axis='both', which='major', **labelp)
plt.tight_layout(pad=30, rect=aspect)
fig.savefig("./new_plots/gv_23_late_tsep_101214_compare_gv.pdf", transparent=True)
# %%

# %%

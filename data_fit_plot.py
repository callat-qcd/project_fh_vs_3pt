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
#%matplotlib inline

# %%
from module.prepare_data import Prepare_data
from module.fit import Fit
from module.plot import plot_pt2
from module.plot import plot_pt3
from module.plot import plot_sum
from module.plot import plot_pt3_no_tra
from module.plot import plot_pt3_no_sca
from module.plot import plot_sum_no_tra
from module.plot import plot_sum_no_sca
from module.prior_setting import prior_ho_width_1 ## notice here
prior = prior_ho_width_1

# %% ########## best fit: ############
# 2pt + 3pt + sum

# 2pt_range = [3, 18], nstates=5
# 3pt_range = [3, 15], tau_cut=1 for [2, 11], tau_cut=2 for [11, 15]
# sum_range = [3, 14], nstates=5

# 2pt + sum

# 2pt_range = [5, 18], nstates=2
# sum_range = [5, 14]

# 2pt + 3pt

# 2pt_range = [3, 18], nstates=5
# 3pt_range = [3, 15], tau_cut=1

# %%
#############################################
#################### 2pt+3pt+sum best fit
#############################################
file_name = 'a09m310_e_gA_srcs0-15.h5'
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

use_p0 = False
include_2pt = True
include_3pt = True
include_sum = True

save = False

#######################################################

data_avg_dict_completed = prepare_data.add_sum_data_scattered(data_avg_dict, sum_tau_cut, sum_A3)


if use_p0 == False:
    best_p0 = 0
    

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

best_p0 = {key: fit_result.p[key].mean for key in fit_result.p}   

#print(fit_result.format(100)) 
print(fit_result)

# %%
########## plot fit and data ###########

print('plot of 2pt+3pt+sum best fit')

# plot_pt2(pt2_data_range, data_avg_dict_completed) # if you just want to plot data, cut the last two parameters: fit_result and fitter

# plot_pt3(pt3_data_range, data_avg_dict_completed, 1) # 1 is tau cut

# plot_sum(pt3_data_range, data_avg_dict_completed)

#plot_pt2(pt2_data_range, data_avg_dict_completed, fit_result, fitter, "_23s", True) # if you just want to plot data, cut the last two parameters: fit_result and fitter

#plot_pt3(pt3_data_range, data_avg_dict_completed, 1, fit_result, fitter, "_23s", True)

#plot_pt3_no_tra(pt3_data_range, data_avg_dict_completed, 1, fit_result, fitter, "_23s", True)

#plot_pt3_no_sca(pt3_data_range, data_avg_dict_completed, 1, fit_result, fitter, "_23s", True)

plot_sum(pt3_data_range, data_avg_dict_completed, fit_result, fitter, pt2_nstates, sum_nstates, "_23s", True) #if you just want to plot data, cut the last three parameters

#plot_sum_no_tra(pt3_data_range, data_avg_dict_completed, fit_result, fitter, pt2_nstates, sum_nstates, "_23s", True)

#plot_sum_no_sca(pt3_data_range, data_avg_dict_completed, fit_result, fitter, pt2_nstates, sum_nstates, "_23s", True)

# %%
# #############################
# ### 2pt+sum best fit
# #############################
# file_name = 'a09m310_e_gA_srcs0-15.h5'
# file_path = os.getcwd() + '/' + file_name # just put the data file inside the 'my_project' folder, if no database, data file should be put in the same path as this file.

# pt2_data_range = [0, 96]
# pt3_data_range = [2, 15]

# prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

# data_avg_dict = prepare_data.read_data_with_average()

# # do the fit

# ################ set parameters #############################

# pt2_t = np.array([5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])
# pt2_nstates = 2

# pt3_A3_t = []
# pt3_A3_tau = []
# pt3_V4_t = []
# pt3_V4_tau = []

# pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
# pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]
# pt3_nstates = pt2_nstates

# sum_A3 = np.array([5, 6, 7, 8, 9, 10, 11, 12, 13])
# sum_V4 = np.array([5, 6, 7, 8, 9, 10, 11, 12, 13])
# sum_nstates = 2
# sum_tau_cut = 1

# use_p0 = False
# include_2pt = True
# include_3pt = False
# include_sum = True

# #########################################################

# data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)


# if use_p0 == False:
#     best_p0 = 0

# fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

# fit_result = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0)[0]

# best_p0 = {key: fit_result.p[key].mean for key in fit_result.p}   

# #print(fit_result.format(maxline=True)) 


# ########## plot fit and data ###########

# print('plot of 2pt+sum best fit')

# plot_pt2(pt2_data_range, data_avg_dict_completed) # if you just want to plot data, cut the last two parameters: fit_result and fitter

# plot_sum(pt3_data_range, data_avg_dict_completed) #if you just want to plot data, cut the last three parameters

# plot_pt2(pt2_data_range, data_avg_dict_completed, fit_result, fitter, "_2s") # if you just want to plot data, cut the last two parameters: fit_result and fitter

# plot_sum(pt3_data_range, data_avg_dict_completed, fit_result, fitter, pt2_nstates, sum_nstates, "_2s") #if you just want to plot data, cut the last three parameters

# # %%
# ############################################
# ### 2pt+3pt best fit
# ##############################################
# file_name = 'a09m310_e_gA_srcs0-15.h5'
# file_path = os.getcwd() + '/' + file_name # just put the data file inside the 'my_project' folder, if no database, data file should be put in the same path as this file.

# pt2_data_range = [0, 96]
# pt3_data_range = [2, 15]

# prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

# data_avg_dict = prepare_data.read_data_with_average()

# # do the fit

# ################ set parameters #############################

# pt2_t = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])
# pt2_nstates = 5

# pt3_A3_t = []
# pt3_A3_tau = []
# pt3_V4_t = []
# pt3_V4_tau = []

# for t in range(3, 15):
#     for tau in range(1, int(t/2)+1):
#         pt3_A3_t.append(t)
#         pt3_A3_tau.append(tau)
#         pt3_V4_t.append(t)
#         pt3_V4_tau.append(tau)

# pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
# pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]
# pt3_nstates = pt2_nstates

# sum_A3 = np.array([])
# sum_V4 = sum_A3 
# sum_nstates = 0
# sum_tau_cut = 1

# use_p0 = False
# include_2pt = True
# include_3pt = True
# include_sum = False

# ######################################################

# data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)


# if use_p0 == False:
#     best_p0 = 0
    
# best_p0 = {'E0': 0.4926861952896794, 'log(dE1)': -1.273594898116278, 'z0': 0.0003260892031046884, 'z1': 0.000263341199546721, 'z0_ps': 0.002999999999999449, 'z1_ps': 8.742252813171513e-17,
#            'log(dE2)': -1.2923954972884248, 'z2': -0.00041038719837353235, 'z2_ps': 1.3527493699158273e-19, 'log(dE3)': -1.2363538425685927, 'z3': -0.00021858370411089373, 'z3_ps': -2.0481201222876496e-19,
#            'log(dE4)': -0.8211404557969483, 'z4': 0.00022859544792623758, 'z4_ps': -2.269364858601806e-19, 'log(E_fh)': -1.2500000000000864, 'z_fh_ss': -2.1292822629149773e-20, 'z_fh_ps': 1.7147835093259457e-19,
#            'A3_00': 1.2661290957780242, 'V4_00': 1.0264985362386028, 'A3_01': -0.3516778615231407, 'V4_01': -0.04468534689899525, 'A3_11': 0.5133635952637498, 'V4_11': 1.0495939630022617, 'A3_02': -0.16492505182541714,
#            'V4_02': -0.10206490798498477, 'A3_12': -0.42836059937108667, 'V4_12': -0.03879580934696547, 'A3_22': 0.490651712437167, 'V4_22': 0.9461073860857595, 'A3_03': 0.3915796384039534, 'V4_03': 0.40932562996043054,
#            'A3_13': 0.1994825622768483, 'V4_13': 0.005796290661118874, 'A3_23': 0.3126615943412394, 'V4_23': -0.3775397962601796, 'A3_33': 0.09835455689160501, 'V4_33': 0.9985879807579602, 'A3_04': 1.0032935731651338,
#            'V4_04': 0.8348906017782508, 'A3_14': 0.01354959867275492, 'V4_14': 0.21364847516473137, 'A3_24': -0.10285387912945705, 'V4_24': -0.243947645810037, 'A3_34': -0.07604449549348145, 'V4_34': -0.04355265360537515,
#            'A3_44': 0.018277411532817922, 'V4_44': 1.0002794679730667, 'fh_A3_0': 0.0, 'fh_V4_0': 0.0, 'd0_ss_A3': -1.5e-06, 'd0_ss_V4': 1.3e-06, 'd0_ps_A3': -9e-06, 'd0_ps_V4': 7.5e-06, 'fh_A3_1': 0.0, 'fh_V4_1': 0.0,
#            'd1_ss_A3': 0.0, 'd1_ss_V4': 0.0, 'd1_ps_A3': 0.0, 'd1_ps_V4': 0.0, 'fh_A3_2': 0.0, 'fh_V4_2': 0.0, 'd2_ss_A3': 0.0, 'd2_ss_V4': 0.0, 'd2_ps_A3': 0.0, 'd2_ps_V4': 0.0, 'fh_A3_3': 0.0, 'fh_V4_3': 0.0,
#            'd3_ss_A3': 0.0, 'd3_ss_V4': 0.0, 'd3_ps_A3': 0.0, 'd3_ps_V4': 0.0, 'd4_ss_A3': 0.0, 'd4_ss_V4': 0.0, 'd4_ps_A3': 0.0, 'd4_ps_V4': 0.0, 'fh_A3_4': 0.0, 'fh_V4_4': 1.0}

# fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

# fit_result = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0)[0]

# best_p0 = {key: fit_result.p[key].mean for key in fit_result.p}   

# #print(fit_result.format(100)) 


# ########## plot fit and data ###########

# print('plot of 2pt+3pt best fit')

# plot_pt2(pt2_data_range, data_avg_dict_completed) # if you just want to plot data, cut the last two parameters: fit_result and fitter

# plot_pt3(pt3_data_range, data_avg_dict_completed, 1)

# plot_pt2(pt2_data_range, data_avg_dict_completed, fit_result, fitter, "_23") # if you just want to plot data, cut the last two parameters: fit_result and fitter

# plot_pt3(pt3_data_range, data_avg_dict_completed, 1, fit_result, fitter, "_23")


# %%
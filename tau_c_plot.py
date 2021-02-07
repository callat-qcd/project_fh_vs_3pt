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

figsize = (7, 4)
aspect=[0.15, 0.15, 0.8, 0.8]
plt.rcParams['figure.figsize'] = figsize

errorp = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1}

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

# pt2_t = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])
# pt2_nstates = 5

# pt3_A3_t = [3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14]
# pt3_A3_tau = [1, 1, 2, 1, 2, 1, 2, 3, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 7]
# pt3_V4_t = pt3_A3_t
# pt3_V4_tau = pt3_A3_tau

# pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
# pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]
# pt3_nstates = pt2_nstates

# sum_A3 = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
# sum_V4 = sum_A3 
# sum_nstates = 5
# sum_tau_cut = 1

# use_p0 = False
# include_2pt = True
# include_3pt = True
# include_sum = True

# save = False

# #######################################################

plt.figure(figsize=figsize)
ax = plt.axes(aspect)
for sum_tau_cut in range(0, 5):

    data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)

    sum_fit_list = []

    tmin = max(2, 2*sum_tau_cut)

    for t in range(tmin, 14):
        sum_fit_list.append(data_avg_dict_completed['sum_A3_fit_'+str(t)])

    temp_mean = np.array([val.mean for val in sum_fit_list])
    temp_sdev = np.array([val.sdev for val in sum_fit_list])

    ax.errorbar(np.arange(tmin, 14), temp_mean, yerr=temp_sdev, label='tau c='+str(sum_tau_cut), marker='o', **errorp)

ax.set_ylim([1, 1.4])
ax.set_xlim([1, 16])
ax.legend(loc='upper right')
plt.savefig(f"./new_plots/tau_c_sum_plot.pdf", transparent=True)
plt.show()
# %%

# %%
import matplotlib.pyplot as plt
import numpy as np 
import matplotlib as mpl
mpl.pyplot.ion()
#%matplotlib inline

# %%
from module.prepare_data import Prepare_data
from module.fit import Fit
from module.plot import moose_plot
from module.plot import moose_23_late_plot
from module.plot import tau_c_plot
from module.plot import m_z_plot
from module.plot import sum_gA_plot
from best_fits import combined_best_fit
from best_fits import late_23_fit

plt.rcParams.update({"text.usetex": True})

file_path = './a09m310_e_gA_srcs0-15.h5'

# %% moose plot
data_avg_dict_completed, fit_result, fitter = combined_best_fit(file_path) 

pt3_data_range = range(2, 15)
tau_c_2 = 1.5 # 2-0.5 to make the color change at half way
divide_n = 11

moose_plot(data_avg_dict_completed, fit_result, fitter, pt3_data_range, tau_c_2, divide_n)

# %% moose plot for 23 late tsep
data_avg_dict_completed, fit_result, fitter = late_23_fit(file_path, True) 

pt3_data_range = [10, 12, 14]

moose_23_late_plot(data_avg_dict_completed, fit_result, fitter, pt3_data_range)

# %% tau_c plot
data_avg_dict_completed, fit_result, fitter = combined_best_fit('./a09m310_e_gA_srcs0-15.h5') 
f_range = [1, 6] # tau_c range of fit part
d_range = [1, 6] # tau_c range of data part

tau_c_plot(file_path, f_range, d_range, fit_result)

# %% m_eff, z_eff and sum_gA
data_avg_dict_completed, fit_result, fitter = combined_best_fit(file_path) 

m_z_plot(data_avg_dict_completed, fit_result, fitter)
sum_gA_plot(data_avg_dict_completed, fit_result, fitter)
# %%

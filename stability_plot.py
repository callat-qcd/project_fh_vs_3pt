# %%
import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  
import math
import os 
import hashlib

from fh_db.ff_fit.models import FF_fit as ff_fit_FF_fit
os.environ["DJANGO_ALLOW_ASYNC_UNSAFE"] = "true" # for jupyter to use database

plt.rcParams.update({"text.usetex": True})

# %%
from module.plot import tmin_plot
from module.plot import tmin_div_plot
from module.plot import tmin_late_plot
from module.plot import tmax_plot
from module.plot import tau_cut_plot

# labels
c2pt_tmin = r"$C_{\textrm{2pt}}\ t_{\textrm{min}}$"
c2pt_tmax = r"$C_{\textrm{2pt}}\ t_{\textrm{max}}$"
c3pt_tmin = r"$C_{\textrm{3pt}}\ t_{\textrm{min}}$"
c3pt_tmax = r"$C_{\textrm{3pt}}\ t_{\textrm{max}}$"
csum_tmin = r"$C_{\textrm{sub}}\ t_{\textrm{min}}$"
csum_tmax = r"$C_{\textrm{sub}}\ t_{\textrm{max}}$"
q_label = r"$Q$"
w_label = r"$w$"
oa00_label = r"$O^A_{00}$"
ov00_label = r"$O^V_{00}$"
e0_label = r"$E_{0}$"
z0_label = r"$z_{0}$"
nstate_label = r"$n_{\textrm{states}}$"
t_label = r"$t$"
tau_label = r"$\tau$"
meff_label = r"$m_{\textrm{eff}}$"
zeff_label = r"$z_{\textrm{eff}}$"
oaeff_label = r"$O^A_{\textrm{eff}}$"
oveff_label = r"$O^V_{\textrm{eff}}$"

gA_ylim=[1.205, 1.32]
gV_ylim=[1.0051, 1.03]
E0_ylim=[0.484, 0.497]

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
####################################
########## 23s 2pt tmin #############
#####################################
n_range=[4, 8]
t_range=[3, 8]
best_n=5
best_t=3

nstate_name='2pt'
tmin_name='2pt'

fit_name='23s'
xlabel=c2pt_tmin

situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=True, 
pt2_tmax=18,
pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=1,
sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5)

tmin_plot(n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, gA_ylim, gV_ylim, E0_ylim, fit_name, xlabel)

# %%
####################################
########## 23s 2pt tmax #############
#####################################
t_range = [14, 20]
best_n = 5
best_t = 18

tmax_name = '2pt'
fit_name = '23s'
xlabel = c2pt_tmax

situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=True,
pt2_tmin=3, pt2_nstates=5,
pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=1,
sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5)

tmax_plot(t_range, best_n, best_t, tmax_name, situation_list, gA_ylim, gV_ylim, E0_ylim, fit_name, xlabel)

# %%
####################################
########## 23s 3pt tmin #############
#####################################
n_range=[4, 8]
t_range=[3, 8]
best_n=5
best_t=3

nstate_name='3pt'
tmin_name='3pt_gA' # here because gA and gV performs differently in 3pt, so you may change paras of gA/ gV separately

fit_name='23s'
xlabel=c3pt_tmin

situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=True, 
pt2_tmin=3, pt2_tmax=18,
pt3_A3_tsep_max=15, id_num=1,
sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5)

tmin_plot(n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, gA_ylim, gV_ylim, E0_ylim, fit_name, xlabel)
# %%
####################################
########## 23s 3pt tmax #############
#####################################
t_range = [8, 16]
best_n = 5
best_t = 15

tmax_name = '3pt' # varying tmax is not so interesting as tmin 
fit_name = '23s'
xlabel = c3pt_tmax

situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=True,
pt2_tmin=3, pt2_tmax=18, pt2_nstates=5,
pt3_A3_tsep_min=3, id_num=1,
sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5)

tmax_plot(t_range, best_n, best_t, tmax_name, situation_list, gA_ylim, gV_ylim, E0_ylim, fit_name, xlabel)

# %%
####################################
########## 23s sum tmin #############
#####################################
n_range=[4, 8]
t_range=[3, 8]
best_n=5
best_t=3

nstate_name='sum'
tmin_name='sum_gA'

fit_name='23s'
xlabel=csum_tmin

situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=True, 
pt2_tmin=3, pt2_tmax=18, pt2_nstates=5,
pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=1,
sum_A3_tsep_max=14, sum_tau_cut=1)

tmin_plot(n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, gA_ylim, gV_ylim, E0_ylim, fit_name, xlabel)

# %%
####################################
########## 23s tau cut #############
#####################################
# taucut-1, taucut, taucut+1, taucut+2, taucut+3
E0 = [0.49225, 0.4904, 0.4897, 0.4896, 0.4894]
E0_err = [0.00049, 0.0016, 0.0018, 0.0020, 0.0022]
A3 = [1.326, 1.253, 1.258, 1.253, 1.253]
A3_err = [0.012, 0.019, 0.022, 0.021, 0.022]
V4 = [1.0224, 1.02244, 1.0230, 1.0227, 1.0227]
V4_err = [0.0011, 0.00087, 0.0020, 0.0020, 0.0021]
Q = [0.076, 0.88, 0.61, 0.45, 0.21]

n = 5

fit_name = '23s'

tau_cut_plot(E0, E0_err, A3, A3_err, V4, V4_err, Q, n, fit_name, gA_ylim, gV_ylim, E0_ylim)
# %%

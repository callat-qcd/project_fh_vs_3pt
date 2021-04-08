# %%
import numpy as np 
import gvar as gv  
import hashlib
import os

from fh_db.ff_fit.models import FF_fit as ff_fit_FF_fit
os.environ["DJANGO_ALLOW_ASYNC_UNSAFE"] = "true" # for jupyter to use database

# %%
from module.prepare_data import Prepare_data
from module.fit import Fit
from module.p0 import best_p0
from module.prior_setting import prior_ho_width_1
prior = prior_ho_width_1

# %%
def fit_and_save(data_file_name, fit_type, save, include_2pt, include_3pt, include_sum, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut, id_num=None, pt2_range=[None, None], pt3_A3_range=[None, None], pt3_V4_range=[None, None], sum_A3_range=[None, None], sum_V4_range=[None, None], pt3_tau_dict=None):
    if fit_type == "scattered":

        pt2_t = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])

        pt2_range = [pt2_t[0], pt2_t[-1]]

        pt3_A3_t = [3, 5, 5, 7, 7, 7]
        pt3_A3_tau = [1, 1, 2, 1, 2, 3]

        pt3_A3_range = [pt3_A3_t[0], pt3_A3_t[-1]]
        pt3_V4_range = pt3_A3_range

        pt3_V4_t = pt3_A3_t
        pt3_V4_tau = pt3_A3_tau
        
        pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
        pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]

        sum_A3 = np.array([3, 5])
        sum_V4 = sum_A3 

        sum_A3_range = [sum_A3[0], sum_A3[-1]]
        sum_V4_range = sum_A3_range

        data_avg_dict_completed = prepare_data.add_sum_data_scattered(data_avg_dict, sum_tau_cut, sum_A3)

    elif fit_type == "continuous":
        pt2_t = np.arange(pt2_range[0], pt2_range[1])

        pt3_A3_t = []
        pt3_A3_tau = []

        for t in range(pt3_A3_range[0], pt3_A3_range[1]):
            for tau in pt3_tau_dict['A3_tsep'+str(t)]:
                pt3_A3_t.append(t)
                pt3_A3_tau.append(tau)

        pt3_V4_t = []
        pt3_V4_tau = []

        for t in range(pt3_V4_range[0], pt3_V4_range[1]):
            for tau in pt3_tau_dict['V4_tsep'+str(t)]:
                pt3_V4_t.append(t)
                pt3_V4_tau.append(tau)

        pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
        pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]

        sum_A3 = np.arange(sum_A3_range[0], sum_A3_range[1])
        sum_V4 = np.arange(sum_V4_range[0], sum_V4_range[1])

        print(pt2_t)
        print(pt3_A3)
        print(pt3_V4)
        print(sum_A3)
        print(sum_V4)

        data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)

    print(pt2_nstates, pt3_nstates, sum_nstates)

    fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

    fit_result, hexcode, dr_hex = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0, save)  

    print(fit_result.format(100))
    print(fit_result.p['A3_00'].mean)
    print(fit_result.p['A3_00'].sdev) 

    #gv.dump(fit_result, 'n2_result'+str(pt2_range[0]))

    if save == True:
        tau_dict = gv.dumps(pt3_tau_dict)

        tau_dict_hex = tau_dict.hex()

        ff_fit_ff_fit, created = ff_fit_FF_fit.objects.get_or_create(
        data_file_name=data_file_name, # Store the name of data file
        include_3pt=include_3pt, # Whether include 3pt part in fitting
        include_sum=include_sum, # Whether include sum part in fitting
        data_and_results=dr_hex, # Use gv.dump to package data, prama, chi2 and dof, and turn into…
        prior_hexcode=hexcode, # put all priors into a string with fixed order and turn it…
        pt2_nstates=pt2_nstates, # Number of states in 2pt fit
        pt3_nstates=pt3_nstates, # Number of states in 3pt fit
        sum_nstates=sum_nstates, # Number of states in sum tau fit
        sum_tau_cut=sum_tau_cut, # Start of sum in sum tau fit
        E0=fit_result.p['E0'].mean, # fit result of E0
        E0_err=fit_result.p['E0'].sdev, # fit result of E0 err
        z0=fit_result.p['z0'].mean, # fit result of z0
        z0_err=fit_result.p['z0'].sdev, # fit result of z0 err
        A300=fit_result.p['A3_00'].mean, # fit result of A300
        A300_err=fit_result.p['A3_00'].sdev, # fit result of A300 err
        V400=fit_result.p['V4_00'].mean, # fit result of V400
        V400_err=fit_result.p['V4_00'].sdev, # fit result of V400 err
        Q_value=fit_result.Q, # Q value of fitting
        log_GBF=fit_result.logGBF, # logGBF of fitting
        fit_type=fit_type, # "scattered" or "continuous" 
        include_2pt=include_2pt, # (Optional) Whether include 2pt part in fitting
        pt2_tmin=pt2_range[0], # (Optional) Minimum of t in 2pt fit
        pt2_tmax=pt2_range[1], # (Optional) Maximum of t in 2pt fit
        pt3_A3_tsep_min=pt3_A3_range[0], # (Optional) Minimum of tsep in 3pt fit
        pt3_A3_tsep_max=pt3_A3_range[1], # (Optional) Maximum of tsep in 3pt fit
        pt3_V4_tsep_min=pt3_V4_range[0], # (Optional) Minimum of tsep in 3pt fit
        pt3_V4_tsep_max=pt3_V4_range[1], # (Optional) Maximum of tsep in 3pt fit
        pt3_tau_dict=tau_dict_hex, # (Optional) Use gv.dump to package tau dict of 3pt, and turn into hexcode
        sum_A3_tsep_min=sum_A3_range[0], # (Optional) Minimum of tsep in sum tau fit
        sum_A3_tsep_max=sum_A3_range[1], # (Optional) Maximum of tsep in sum tau fit
        sum_V4_tsep_min=sum_V4_range[0], # (Optional) Minimum of tsep in sum tau fit
        sum_V4_tsep_max=sum_V4_range[1], # (Optional) Maximum of tsep in sum tau fit
        id_num=id_num, # (Optional) used to divide and select
        )



# %%
file_name = 'a09m310_e_gA_srcs0-15.h5'
file_path = './' + file_name # just put the data file in the same path as this code.

pt2_data_range = [0, 96] # range in data file
pt3_data_range = [2, 15] # range in data file

prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

data_avg_dict = prepare_data.read_data_with_average()

#print([i.mean for i in data_avg_dict['pt3_A3_tsep_10']])

include_2pt = True
include_3pt = True
include_sum = True

sum_tau_cut = 1
id_num = 1

fit_type = "continuous" 
save = False

##################################################################################################################
#########################################################

if fit_type == "scattered":

    pt2_nstates = 5
    pt3_nstates = 5
    sum_nstates = 5

    fit_and_save(file_name, fit_type, save, include_2pt, include_3pt, include_sum, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut, id_num)


############################################################

elif fit_type == "continuous":
    pt3_tau_dict = dict()
    for t in range(2, 11):
        pt3_tau_dict['A3_tsep'+str(t)] = np.arange(1, int(t/2)+1)
        pt3_tau_dict['V4_tsep'+str(t)] = np.arange(1, int(t/2)+1)

    for t in range(11, 15):
        pt3_tau_dict['A3_tsep'+str(t)] = np.arange(2, int(t/2)+1) # for 23 fit, start from 1
        pt3_tau_dict['V4_tsep'+str(t)] = np.arange(2, int(t/2)+1)

    pt2_range_list = [[tmin, tmax] for tmin in range(3, 4) for tmax in range(18, 19)]
    pt2_nstates_list = [nstates for nstates in range(5, 6)]

    pt3_A3_range_list = [[tmin, tmax] for tmin in range(3, 4) for tmax in range(15, 16)]
    #pt3_V4_range_list = [[tmin, tmax] for tmin in range(3, 4) for tmax in range(15, 16)]
    #pt3_nstates_list = [nstates for nstates in range(5, 6)] # pt2_nstates = pt3_nstates

    sum_A3_range_list = [[tmin, tmax] for tmin in range(3, 4) for tmax in range(14, 15)]
    #sum_V4_range_list = [[tmin, tmax] for tmin in range(3, 4) for tmax in range(14, 15)]
    sum_nstates_list = [nstates for nstates in range(5, 6)]

    times = 0

    for pt2_range in pt2_range_list:
        for pt2_nstates in pt2_nstates_list:
            for pt3_A3_range in pt3_A3_range_list:
                for sum_A3_range in sum_A3_range_list:
                    for sum_nstates in sum_nstates_list:
                        print('#################################################')
                        times += 1
                        pt3_nstates = pt2_nstates 

                        #sum_A3_range = [pt3_A3_range[0], pt3_A3_range[1]-1]
                        print(sum_A3_range)

                        pt3_V4_range = pt3_A3_range ##
                        sum_V4_range = sum_A3_range ##

                        fit_and_save(file_name, fit_type, save, include_2pt, include_3pt, include_sum, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut, id_num, pt2_range=pt2_range, pt3_A3_range=pt3_A3_range, pt3_V4_range=pt3_V4_range, sum_A3_range=sum_A3_range, sum_V4_range=sum_V4_range, pt3_tau_dict=pt3_tau_dict)
    print('total '+str(times)+' fits')
#########################################################


# %%
t_list = []
tau_list = []

for t in range(3, 15): # 15 -> 6
    if t % 2 == 1:
        if t in range(2, 11):
            for tau in range(1, int(t/2)+1):
                t_list.append(t)
                tau_list.append(tau)

        elif t in range(11, 15):
            for tau in range(2, int(t/2)+1):
                t_list.append(t)
                tau_list.append(tau)

print(t_list)
print(tau_list)
# %%

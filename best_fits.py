# %% ########## best fits: ############
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
import numpy as np 

from module.prepare_data import Prepare_data
from module.fit import Fit
from module.p0 import best_p0
from module.prior_setting import prior_ho_width_1
prior = prior_ho_width_1

# %%
pt2_data_range = [0, 96]
pt3_data_range = [2, 15]

# %%
def combined_best_fit(file_path):
    file_name = 'a09m310_e_gA_srcs0-15.h5'

    prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

    data_avg_dict = prepare_data.read_data_with_average()

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

    data_avg_dict_completed = prepare_data.add_sum_data_scattered(data_avg_dict, sum_tau_cut, sum_A3)

    fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

    fit_result = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0)[0]

    return data_avg_dict_completed, fit_result, fitter

def late_23_fit(file_path, best_fit=False):
    # stability plot of late tsep [10, 12, 14] # varying nstates and tau_cut
    file_name = 'a09m310_e_gA_srcs0-15.h5'

    prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

    data_avg_dict = prepare_data.read_data_with_average()

    ######################################################
    ################ set parameters ######################

    pt3_A3_t = []
    pt3_A3_tau = []

    t_dict = {}
    tau_dict = {}
    A3 = {}
    A3_err = {}
    Q = {}

    t_dict['opt-1'] = [10, 12, 14]
    tau_dict['opt-1'] = [5, 6, 7]
    t_dict['opt'] = [10, 12, 12, 14, 14]
    tau_dict['opt'] = [5, 5, 6, 6, 7]
    t_dict['opt+1'] = [10, 10, 12, 12, 12, 14, 14, 14]
    tau_dict['opt+1'] = [4, 5, 4, 5, 6, 5, 6, 7]

    if best_fit == False:
        situation_list = ['opt+1', 'opt', 'opt-1']
        nstate_range = range(1, 5)

    elif best_fit == True:
        situation_list = ['opt']
        nstate_range = range(2, 3)

    for situation in situation_list:
        A3[situation] = []
        A3_err[situation] = []
        Q[situation] = []

        pt3_A3_t = t_dict[situation]
        pt3_A3_tau = tau_dict[situation]
        pt3_V4_t = pt3_A3_t
        pt3_V4_tau = pt3_A3_tau

        pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
        pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]

        sum_A3 = np.array([])
        sum_V4 = sum_A3 
        sum_nstates = 0
        sum_tau_cut = 1

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

        data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)

        for nstates in nstate_range:
            pt2_nstates = nstates
            pt3_nstates = pt2_nstates
            pt2_t = np.array(pt2_t_list[nstates])

            fitter = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

            fit_result = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0)[0]  

            A3[situation].append(fit_result.p['A3_00'].mean)
            A3_err[situation].append(fit_result.p['A3_00'].sdev)
            Q[situation].append(fit_result.Q)

    if best_fit == False:
        return A3, A3_err, Q

    elif best_fit == True:
        return data_avg_dict_completed, fit_result, fitter
    
    
# %%
if __name__ == '__main__':
    data_avg_dict_completed, fit_result, fitter = combined_best_fit('/home/greyyy/Desktop/qcd/fh_vs_3pt/a09m310_e_gA_srcs0-15.h5') # path of data file

    # %%
    data_avg_dict_completed, fit_result, fitter = late_23_fit('/home/greyyy/Desktop/qcd/fh_vs_3pt/a09m310_e_gA_srcs0-15.h5', True)

# %%
# import numpy as np
# from module.prior_setting import prior_ho_width_1
# from module.prior_setting import prior_sw_width_1
# from module.prior_setting import prior_id1_width_1
# from module.prior_setting import prior_id2_width_1

# prior = prior_ho_width_1(5, 5, 5)
# prior['E1'] = prior['E0'] + prior['dE1']
# prior['E2'] = prior['E0'] + prior['dE1'] + prior['dE2']
# prior['E3'] = prior['E0'] + prior['dE1'] + prior['dE2'] + prior['dE3']
# prior['E4'] = prior['E0'] + prior['dE1'] + prior['dE2'] + prior['dE3'] + prior['dE4']

# print(prior['E0'])
# print(prior['dE1'])

# print(prior['E1'], prior['E2'], prior['E3'], prior['E4'])

# %%
# %%
import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  

import matplotlib as mpl
mpl.pyplot.ion()
#%matplotlib inline

# %%
from module.prepare_data import Prepare_data
from module.fit import Fit
from module.plot import excited_states_plot
from module.prior_setting import prior_ho_width_1
prior = prior_ho_width_1

from best_fits import combined_best_fit


plt.rcParams.update({"text.usetex": True})

file_path = './a09m310_e_gA_srcs0-15.h5'
file_name = 'a09m310_e_gA_srcs0-15.h5'
pt3_data_start = 2
pt3_data_end = 15 # sum_end = pt3_end - 1
linspace_num = 100

pt2_nstates = 5
pt3_nstates = pt2_nstates
sum_nstates = 5
include_2pt = True
include_3pt = True
include_sum = True
sum_tau_cut_plot = 1 # this is for plot

# %%
def pt3_gA_ratio(data_avg_dict_completed, fit_result, fitter):
    # tra data: data_avg_dict - fit_result with p_ntra
    # sca data: data_avg_dict - fit_result with p_nsca, here sca does not include g.s.
    # gs data: data_avg_dict - fit_result with p_ngs

    # tra fit: fit_result with p_tra
    # sca fit: fit_result with p_sca
    # gs fit: fit_result with p_gs

    gA_tsep_data = [] # errorbar do not need plot density
    gA_tau_data = [] 
    pt2_data = []  
    pt3_gA_data = []  

    

    gA_tsep_fit = np.linspace(0, 35, linspace_num)
    gA_tau_fit = gA_tsep_fit / 2 

    for t in range(pt3_data_start, pt3_data_end): # tsep and tau to calc fit values of no transition 
        if t % 2 == 0: # only when tsep is even, tau = tsep/2
            gA_tsep_data.append(t)
            gA_tau_data.append(t/2) # tau = tsep/2
            pt2_data.append(data_avg_dict_completed['pt2_tsep_'+str(t)])
            pt3_gA_data.append(data_avg_dict_completed['pt3_A3_tsep_'+str(t)][int(t/2)])

    pt2_data = np.array(pt2_data)
    pt3_gA_data = np.array(pt3_gA_data)


    pt2_fit = fitter.pt2_fit_function(gA_tsep_fit, fit_result.p)['pt2']

    ########################
    ####### tra data #######
    p_ntra = gv.BufferDict(fit_result.p) # no transition
    for key in p_ntra:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] != key.split("_")[1][1]):
            p_ntra[key] = 0


    # fit function values of no transition
    pt3_gA_fit_ntra = fitter.pt3_fit_function(np.array(gA_tsep_data), np.array(gA_tsep_data), np.array(gA_tau_data), np.array(gA_tau_data), p_ntra)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt3_gA_tra_data = pt3_gA_data - pt3_gA_fit_ntra
    pt3_gA_tra_data = pt3_gA_tra_data / pt2_data

    ########################
    ####### sca data #######
    p_nsca = gv.BufferDict(fit_result.p) # no scattering, but have g.s.
    for key in p_nsca:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] == key.split("_")[1][1]) and (key.split("_")[1][1] != '0'):
            p_nsca[key] = 0

    # fit function values of no transition
    pt3_gA_fit_nsca = fitter.pt3_fit_function(np.array(gA_tsep_data), np.array(gA_tsep_data), np.array(gA_tau_data), np.array(gA_tau_data), p_nsca)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt3_gA_sca_data = pt3_gA_data - pt3_gA_fit_nsca
    pt3_gA_sca_data = pt3_gA_sca_data / pt2_data

    #####################
    ###### gs data ######
    p_ngs = gv.BufferDict(fit_result.p) # no g.s.
    for key in p_ngs:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][1] == '0'):
            p_ngs[key] = 0
    p_ngs['z0'] = 0

    # fit function values of no transition
    pt2_fit_ngs = fitter.pt2_fit_function(np.array(gA_tsep_data), p_ngs)['pt2']

    pt3_gA_fit_ngs = fitter.pt3_fit_function(np.array(gA_tsep_data), np.array(gA_tsep_data), np.array(gA_tau_data), np.array(gA_tau_data), p_ngs)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt2_gs_data = pt2_data - pt2_fit_ngs
    pt3_gA_gs_data = pt3_gA_data - pt3_gA_fit_ngs
    pt3_gA_gs_data = pt3_gA_gs_data / pt2_gs_data

    #####################
    ###### all data #####

    # fit function values of no transition
    pt3_gA_fit_ngs = fitter.pt3_fit_function(np.array(gA_tsep_data), np.array(gA_tsep_data), np.array(gA_tau_data), np.array(gA_tau_data), p_ngs)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt3_gA_all_data = pt3_gA_data - pt3_gA_fit_ngs

    pt3_gA_all_data = pt3_gA_all_data / pt2_data

    pt3_gA_all_data += pt3_gA_sca_data
    pt3_gA_all_data += pt3_gA_tra_data


    ########################
    ####### tra fit ########
    p_tra = gv.BufferDict(fit_result.p) 
    for key in p_tra:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] == key.split("_")[1][1]):
            p_tra[key] = 0

    # fit function values of only transition
    pt3_gA_tra_fit = fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_tra)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt3_gA_tra_fit = pt3_gA_tra_fit / pt2_fit

    ########################
    ####### sca fit ########
    p_sca = gv.BufferDict(fit_result.p) # no g.s.
    for key in p_sca:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] != key.split("_")[1][1]):
            p_sca[key] = 0
    p_sca['A3_00'] = 0


    # fit function values of only transition
    pt3_gA_sca_fit = fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_sca)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt3_gA_sca_fit = pt3_gA_sca_fit / pt2_fit

    ########################
    ####### gs fit #########
    p_gs = gv.BufferDict(fit_result.p) # no g.s.
    for key in p_gs:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][1] != '0'):
            p_gs[key] = 0

        if ('z' in key and key != 'z0'):
            p_gs[key] = 0

    # fit function values of only g.s.
    pt2_gs_fit = fitter.pt2_fit_function(np.array(gA_tsep_fit),p_gs)['pt2']

    pt3_gA_gs_fit = fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_gs)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt3_gA_gs_fit = pt3_gA_gs_fit / pt2_gs_fit


    ########################
    ####### all fit ########

    # fit function values of only g.s.
    pt3_gA_all_fit = fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_gs)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt3_gA_all_fit += fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_sca)['pt3_A3']

    pt3_gA_all_fit += fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_tra)['pt3_A3']

    pt3_gA_all_fit = pt3_gA_all_fit / pt2_fit

    pt3_list = [pt3_gA_tra_fit, pt3_gA_sca_fit, pt3_gA_gs_fit, pt3_gA_all_fit, pt3_gA_tra_data, pt3_gA_sca_data, pt3_gA_gs_data, pt3_gA_all_data]

    return pt3_list

def sum_gA_ratio(data_avg_dict_completed, fit_result, fitter):
    # tra data: ( data_avg_dict[sum] - fit_result[sum] with p_ntra ) / data_avg_dict[pt2]
    # sca data: ( data_avg_dict[sum] - fit_result[sum] with p_nsca ) / data_avg_dict[pt2], here sca does not include g.s.
    # gs data: ( data_avg_dict[sum] - fit_result[sum] with p_ngs ) / ( data_avg_dict[pt2] - fit_result[2pt] with p_ngs )

    # tra fit: fit_result[sum] / fit_result[pt2] with p_tra then subtract
    # sca fit: fit_result[sum] / fit_result[pt2] with p_sca then subtract
    # gs fit: fit_result[sum] / fit_result[pt2] with p_gs for both then subtract

    sum_data_range = [pt3_data_start, pt3_data_end-1]

    gA_tsep_data = [] # errorbar do not need plot density
    pt2_data = []  
    sum_gA_data = []  

    gA_tsep_fit = np.linspace(0, 35, linspace_num)

    plot_gap = (35 - 0) / (linspace_num - 1) 

    for t in range(pt3_data_start, pt3_data_end): # tsep to calc fit values of no transition 
        gA_tsep_data.append(t)
        pt2_data.append(data_avg_dict_completed['pt2_tsep_'+str(t)])
        sum_gA_data.append(data_avg_dict_completed['sum_A3_tsep_'+str(t)])

    pt2_data = np.array(pt2_data)
    sum_gA_data = np.array(sum_gA_data)

    pt2_fit = fitter.pt2_fit_function(gA_tsep_fit, fit_result.p)['pt2']

    ########################
    ####### tra data #######
    p_ntra = gv.BufferDict(fit_result.p) # no transition
    for key in p_ntra:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] != key.split("_")[1][1]):
            p_ntra[key] = 0

    # fit function values of no transition
    sum_gA_fit_ntra = fitter.summation_same_can(np.array(gA_tsep_data), np.array(gA_tsep_data), p_ntra)['sum_A3'] # here gV_tsep same as gA

    ratio_gA_tra_data = (sum_gA_data - sum_gA_fit_ntra) / pt2_data

    sum_gA_tra_data = []

    for t in range(sum_data_range[1] - sum_data_range[0]):
        sum_gA_tra_data.append(ratio_gA_tra_data[t+1] - ratio_gA_tra_data[t])

    sum_gA_tra_data = np.array(sum_gA_tra_data)

    ########################
    ####### sca data #######
    p_nsca = gv.BufferDict(fit_result.p) # no scattering, but have g.s.
    for key in p_nsca:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] == key.split("_")[1][1]) and (key.split("_")[1][1] != '0'):
            p_nsca[key] = 0

    # fit function values of no transition
    sum_gA_fit_nsca = fitter.summation_same_can(np.array(gA_tsep_data), np.array(gA_tsep_data), p_nsca)['sum_A3'] # here gV_tsep same as gA

    ratio_gA_sca_data = (sum_gA_data - sum_gA_fit_nsca) / pt2_data

    sum_gA_sca_data = []

    for t in range(sum_data_range[1] - sum_data_range[0]):
        sum_gA_sca_data.append(ratio_gA_sca_data[t+1] - ratio_gA_sca_data[t])

    sum_gA_sca_data = np.array(sum_gA_sca_data)

    #####################
    ###### gs data ######      summation of C3_gs / C2_gs 
    p_ngs = gv.BufferDict(fit_result.p) # no g.s.
    for key in p_ngs:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][1] == '0'):
            p_ngs[key] = 0
    p_ngs['z0'] = 0

    # fit function values of no transition
    pt2_fit_ngs = fitter.pt2_fit_function(np.array(gA_tsep_data), p_ngs)['pt2']

    sum_gA_fit_ngs = fitter.summation_same_can(np.array(gA_tsep_data), np.array(gA_tsep_data), p_ngs)['sum_A3'] # here gV_tsep same as gA

    ratio_gA_gs_data = (sum_gA_data - sum_gA_fit_ngs) / (pt2_data - pt2_fit_ngs)

    sum_gA_gs_data = []

    for t in range(sum_data_range[1] - sum_data_range[0]):
        sum_gA_gs_data.append(ratio_gA_gs_data[t+1] - ratio_gA_gs_data[t])

    sum_gA_gs_data = np.array(sum_gA_gs_data)

    #####################
    ###### all data #####      

    # fit function values of no transition

    sum_gA_fit_ngs = fitter.summation_same_can(np.array(gA_tsep_data), np.array(gA_tsep_data), p_ngs)['sum_A3'] # here gV_tsep same as gA

    ratio_gA_all_data = (sum_gA_data - sum_gA_fit_ngs) / (pt2_data)
    ratio_gA_all_data += ratio_gA_sca_data
    ratio_gA_all_data += ratio_gA_tra_data

    sum_gA_all_data = []

    for t in range(sum_data_range[1] - sum_data_range[0]):
        sum_gA_all_data.append(ratio_gA_all_data[t+1] - ratio_gA_all_data[t])

    sum_gA_all_data = np.array(sum_gA_all_data)

    ########################
    ####### tra fit ########
    p_tra = gv.BufferDict(fit_result.p) 
    for key in p_tra:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] == key.split("_")[1][1]):
            p_tra[key] = 0

    # fit function values of only transition
    ratio_gA_tra_fit = fitter.summation_same_can(np.array(gA_tsep_fit), np.array(gA_tsep_fit), p_tra)['sum_A3'] / pt2_fit # here gV_tsep same as gA

    sum_gA_tra_fit = []

    for t in range(len(gA_tsep_fit) - int(1/plot_gap)):
        temp = (ratio_gA_tra_fit[t+int(1/plot_gap)] - ratio_gA_tra_fit[t]) # here gA[t+1] - gA[t] != 1
        sum_gA_tra_fit.append(temp)

    sum_gA_tra_fit = np.array(sum_gA_tra_fit)

    ########################
    ####### sca fit ########
    p_sca = gv.BufferDict(fit_result.p) # sca only
    for key in p_sca:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] != key.split("_")[1][1]):
            p_sca[key] = 0
    p_sca['A3_00'] = 0

    # fit function values of only transition
    ratio_gA_sca_fit = fitter.summation_same_can(np.array(gA_tsep_fit), np.array(gA_tsep_fit), p_sca)['sum_A3'] / pt2_fit # here gV_tsep same as gA

    sum_gA_sca_fit = []

    for t in range(len(gA_tsep_fit) - int(1/plot_gap)):
        temp = (ratio_gA_sca_fit[t+int(1/plot_gap)] - ratio_gA_sca_fit[t]) # here gA[t+1] - gA[t] != 1
        sum_gA_sca_fit.append(temp)

    sum_gA_sca_fit = np.array(sum_gA_sca_fit)

    ########################
    ####### gs fit #########
    p_gs = gv.BufferDict(fit_result.p) # g.s. only
    for key in p_gs:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][1] != '0'):
            p_gs[key] = 0
        
        if ('z' in key and key != 'z0'):
            p_gs[key] = 0


    # fit function values of only transition
    pt2_gs_fit = fitter.pt2_fit_function(np.array(gA_tsep_fit), p_gs)['pt2']
    ratio_gA_gs_fit = fitter.summation_same_can(np.array(gA_tsep_fit), np.array(gA_tsep_fit), p_gs)['sum_A3'] / pt2_gs_fit # here gV_tsep same as gA

    sum_gA_gs_fit = []

    for t in range(len(gA_tsep_fit) - int(1/plot_gap)):
        temp = (ratio_gA_gs_fit[t+int(1/plot_gap)] - ratio_gA_gs_fit[t]) # here gA[t+1] - gA[t] != 1
        sum_gA_gs_fit.append(temp)

    sum_gA_gs_fit = np.array(sum_gA_gs_fit)

    ########################
    ####### all fit ########

    # fit function values of only transition
    ratio_gA_all_fit = fitter.summation_same_can(np.array(gA_tsep_fit), np.array(gA_tsep_fit), p_gs)['sum_A3'] / pt2_fit # here gV_tsep same as gA

    ratio_gA_all_fit += fitter.summation_same_can(np.array(gA_tsep_fit), np.array(gA_tsep_fit), p_sca)['sum_A3'] / pt2_fit

    ratio_gA_all_fit += fitter.summation_same_can(np.array(gA_tsep_fit), np.array(gA_tsep_fit), p_tra)['sum_A3'] / pt2_fit

    sum_gA_all_fit = []

    for t in range(len(gA_tsep_fit) - int(1/plot_gap)):
        temp = (ratio_gA_all_fit[t+int(1/plot_gap)] - ratio_gA_all_fit[t]) # here gA[t+1] - gA[t] != 1
        sum_gA_all_fit.append(temp)

    sum_gA_all_fit = np.array(sum_gA_all_fit)

    ratio_list = [ratio_gA_tra_fit, ratio_gA_sca_fit, ratio_gA_gs_fit, ratio_gA_all_fit, ratio_gA_tra_data, ratio_gA_sca_data, ratio_gA_gs_data, ratio_gA_all_data]

    return ratio_list

# %%
data_avg_dict_completed, fit_result, fitter = combined_best_fit(file_path) 
fitter_plot = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut_plot, include_2pt, include_3pt, include_sum)

pt3_list = pt3_gA_ratio(data_avg_dict_completed, fit_result, fitter_plot)
ratio_list = sum_gA_ratio(data_avg_dict_completed, fit_result, fitter_plot)

excited_states_plot(pt3_list, ratio_list)
# %%

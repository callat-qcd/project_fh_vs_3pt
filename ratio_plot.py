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

figsize = (7, 4)
aspect=[0.15, 0.15, 0.8, 0.8]
plt.rcParams['figure.figsize'] = figsize

errorp = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1}
labelp = {"labelsize": 14}
textp = {"fontsize": 15}

# %%
def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))

def pt3_gA_ratio(pt3_data_range, data_avg_dict_completed, fit_result, fitter, div_2pt):
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

    linspace_num = 2000

    gA_tsep_fit = np.linspace(0, 35, linspace_num)
    gA_tau_fit = gA_tsep_fit / 2 

    for t in range(pt3_data_range[0], pt3_data_range[1]): # tsep and tau to calc fit values of no transition 
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

    # print('\n p_ntra: \n')
    # print(p_ntra)

    # fit function values of no transition
    pt3_gA_fit_ntra = fitter.pt3_fit_function(np.array(gA_tsep_data), np.array(gA_tsep_data), np.array(gA_tau_data), np.array(gA_tau_data), p_ntra)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt3_gA_tra_data = pt3_gA_data - pt3_gA_fit_ntra
    if div_2pt == True:
        pt3_gA_tra_data = pt3_gA_tra_data / pt2_data

    ########################
    ####### sca data #######
    p_nsca = gv.BufferDict(fit_result.p) # no scattering, but have g.s.
    for key in p_nsca:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] == key.split("_")[1][1]) and (key.split("_")[1][1] != '0'):
            p_nsca[key] = 0

    # print('\n p_nsca: \n')
    # print(p_nsca)

    # fit function values of no transition
    pt3_gA_fit_nsca = fitter.pt3_fit_function(np.array(gA_tsep_data), np.array(gA_tsep_data), np.array(gA_tau_data), np.array(gA_tau_data), p_nsca)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt3_gA_sca_data = pt3_gA_data - pt3_gA_fit_nsca
    if div_2pt == True:
        pt3_gA_sca_data = pt3_gA_sca_data / pt2_data

    #####################
    ###### gs data ######
    p_ngs = gv.BufferDict(fit_result.p) # no g.s.
    for key in p_ngs:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][1] == '0'):
            p_ngs[key] = 0
    p_ngs['z0'] = 0

    # print('\n p_ngs: \n')
    # print(p_ngs)

    # fit function values of no transition
    pt2_fit_ngs = fitter.pt2_fit_function(np.array(gA_tsep_data), p_ngs)['pt2']

    pt3_gA_fit_ngs = fitter.pt3_fit_function(np.array(gA_tsep_data), np.array(gA_tsep_data), np.array(gA_tau_data), np.array(gA_tau_data), p_ngs)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt2_gs_data = pt2_data - pt2_fit_ngs
    pt3_gA_gs_data = pt3_gA_data - pt3_gA_fit_ngs
    if div_2pt == True:
        pt3_gA_gs_data = pt3_gA_gs_data / pt2_gs_data

    #####################
    ###### all data #####

    # fit function values of no transition
    pt3_gA_fit_ngs = fitter.pt3_fit_function(np.array(gA_tsep_data), np.array(gA_tsep_data), np.array(gA_tau_data), np.array(gA_tau_data), p_ngs)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt3_gA_all_data = pt3_gA_data - pt3_gA_fit_ngs

    if div_2pt == True:
        pt3_gA_all_data = pt3_gA_all_data / pt2_data

    pt3_gA_all_data += pt3_gA_sca_data
    pt3_gA_all_data += pt3_gA_tra_data


    ########################
    ####### tra fit ########
    p_tra = gv.BufferDict(fit_result.p) 
    for key in p_tra:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] == key.split("_")[1][1]):
            p_tra[key] = 0

    # print('\n p_tra: \n')
    # print(p_tra)

    # fit function values of only transition
    pt3_gA_tra_fit = fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_tra)['pt3_A3'] # here gV_tsep and gV_tau same as gA
    if div_2pt == True:
        pt3_gA_tra_fit = pt3_gA_tra_fit / pt2_fit

    ########################
    ####### sca fit ########
    p_sca = gv.BufferDict(fit_result.p) # no g.s.
    for key in p_sca:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] != key.split("_")[1][1]):
            p_sca[key] = 0
    p_sca['A3_00'] = 0

    # print('\n p_sca: \n')
    # print(p_sca)

    # fit function values of only transition
    pt3_gA_sca_fit = fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_sca)['pt3_A3'] # here gV_tsep and gV_tau same as gA
    if div_2pt == True:
        pt3_gA_sca_fit = pt3_gA_sca_fit / pt2_fit

    ########################
    ####### gs fit #########
    p_gs = gv.BufferDict(fit_result.p) # no g.s.
    for key in p_gs:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][1] != '0'):
            p_gs[key] = 0

        if ('z' in key and key != 'z0'):
            p_gs[key] = 0

    # print('\n p_gs: \n')
    # print(p_gs)

    # fit function values of only g.s.
    pt2_gs_fit = fitter.pt2_fit_function(np.array(gA_tsep_fit),p_gs)['pt2']

    pt3_gA_gs_fit = fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_gs)['pt3_A3'] # here gV_tsep and gV_tau same as gA
    if div_2pt == True:
        pt3_gA_gs_fit = pt3_gA_gs_fit / pt2_gs_fit


    ########################
    ####### all fit ########

    # fit function values of only g.s.
    pt3_gA_all_fit = fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_gs)['pt3_A3'] # here gV_tsep and gV_tau same as gA

    pt3_gA_all_fit += fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_sca)['pt3_A3']

    pt3_gA_all_fit += fitter.pt3_fit_function(np.array(gA_tsep_fit), np.array(gA_tsep_fit), np.array(gA_tau_fit), np.array(gA_tau_fit), p_tra)['pt3_A3']

    if div_2pt == True:
        pt3_gA_all_fit = pt3_gA_all_fit / pt2_fit

    return pt3_gA_tra_fit, pt3_gA_sca_fit, pt3_gA_gs_fit, pt3_gA_all_fit, pt3_gA_tra_data, pt3_gA_sca_data, pt3_gA_gs_data, pt3_gA_all_data


def sum_gA_ratio(pt3_data_range, data_avg_dict_completed, fit_result, fitter):
    # tra data: ( data_avg_dict[sum] - fit_result[sum] with p_ntra ) / data_avg_dict[pt2]
    # sca data: ( data_avg_dict[sum] - fit_result[sum] with p_nsca ) / data_avg_dict[pt2], here sca does not include g.s.
    # gs data: ( data_avg_dict[sum] - fit_result[sum] with p_ngs ) / ( data_avg_dict[pt2] - fit_result[2pt] with p_ngs )

    # tra fit: fit_result[sum] / fit_result[pt2] with p_tra then subtract
    # sca fit: fit_result[sum] / fit_result[pt2] with p_sca then subtract
    # gs fit: fit_result[sum] / fit_result[pt2] with p_gs for both then subtract

    sum_data_range = [pt3_data_range[0], pt3_data_range[1]-1]

    gA_tsep_data = [] # errorbar do not need plot density
    pt2_data = []  
    sum_gA_data = []  

    linspace_num = 2000

    gA_tsep_fit = np.linspace(0, 35, linspace_num)

    plot_gap = (35 - 0) / (linspace_num - 1) 

    for t in range(pt3_data_range[0], pt3_data_range[1]): # tsep to calc fit values of no transition 
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

    # print('\n p_ntra: \n')
    # print(p_ntra)

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

    # print('\n p_nsca: \n')
    # print(p_nsca)

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

    # print('\n p_ngs: \n')
    # print(p_ngs)

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

    # print('\n p_tra: \n')
    # print(p_tra)

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

    # print('\n p_sca: \n')
    # print(p_sca)

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

    # print('\n p_gs: \n')
    # print(p_gs)

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

    # print('\n p_gs: \n')
    # print(p_gs)

    # print('\n p_sca: \n')
    # print(p_sca)

    # print('\n p_tra: \n')
    # print(p_tra)

    # fit function values of only transition
    ratio_gA_all_fit = fitter.summation_same_can(np.array(gA_tsep_fit), np.array(gA_tsep_fit), p_gs)['sum_A3'] / pt2_fit # here gV_tsep same as gA

    ratio_gA_all_fit += fitter.summation_same_can(np.array(gA_tsep_fit), np.array(gA_tsep_fit), p_sca)['sum_A3'] / pt2_fit

    ratio_gA_all_fit += fitter.summation_same_can(np.array(gA_tsep_fit), np.array(gA_tsep_fit), p_tra)['sum_A3'] / pt2_fit

    sum_gA_all_fit = []

    for t in range(len(gA_tsep_fit) - int(1/plot_gap)):
        temp = (ratio_gA_all_fit[t+int(1/plot_gap)] - ratio_gA_all_fit[t]) # here gA[t+1] - gA[t] != 1
        sum_gA_all_fit.append(temp)

    sum_gA_all_fit = np.array(sum_gA_all_fit)


    return plot_gap, ratio_gA_tra_fit, ratio_gA_sca_fit, ratio_gA_gs_fit, ratio_gA_all_fit, ratio_gA_tra_data, ratio_gA_sca_data, ratio_gA_gs_data, ratio_gA_all_data

def ratio_plot(pt3_data_range, plot_gap, pt3_gA_tra_fit, pt3_gA_sca_fit, pt3_gA_gs_fit, pt3_gA_all_fit, pt3_gA_tra_data, pt3_gA_sca_data, pt3_gA_gs_data, pt3_gA_all_data, ratio_gA_tra_fit, ratio_gA_sca_fit, ratio_gA_gs_fit, ratio_gA_all_fit, ratio_gA_tra_data, ratio_gA_sca_data, ratio_gA_gs_data, ratio_gA_all_data, div_2pt, sum_tau_cut_plot, plot_in_fm):

    if plot_in_fm == True:
        omega_imp_a09 = 0.08730 # converse lattice to fm
        plot_xmax = 3
    elif plot_in_fm == False:
        omega_imp_a09 = 1
        plot_xmax = 3 / 0.08730

    linspace_num = 2000

    pt3_tsep_data = np.arange(pt3_data_range[0], pt3_data_range[1], 2) # errorbar do not need plot density

    pt3_tsep_fit = np.linspace(0, 35, linspace_num) # till 3fm

    plot_gap = (35 - 0) / (linspace_num - 1) 

    sum_tsep_data = np.arange(pt3_data_range[0], pt3_data_range[1]-1)

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



    ###### tra plot ######
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)

    temp_mean = np.array([val.mean for val in (pt3_gA_tra_fit / pt3_gA_gs_fit)]) 
    temp_sdev = np.array([val.sdev for val in (pt3_gA_tra_fit / pt3_gA_gs_fit)])
    ax.fill_between(pt3_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, color=orange, alpha=0.3, label='3pt tra')

    temp_mean = np.array([val.mean for val in (pt3_gA_tra_data / pt3_gA_gs_data)])
    temp_sdev = np.array([val.sdev for val in (pt3_gA_tra_data / pt3_gA_gs_data)])
    ax.errorbar(pt3_tsep_data, temp_mean, yerr=temp_sdev, marker='o', color=orange, **errorp)

    temp_mean = np.array([val.mean for val in (sum_gA_tra_fit / sum_gA_gs_fit)])
    temp_sdev = np.array([val.sdev for val in (sum_gA_tra_fit / sum_gA_gs_fit)])
    ax.fill_between(sum_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, color=blue, alpha=0.3, label='sum tra')

    temp_mean = np.array([val.mean for val in (sum_gA_tra_data / sum_gA_gs_data)])
    temp_sdev = np.array([val.sdev for val in (sum_gA_tra_data / sum_gA_gs_data)])
    ax.errorbar(sum_tsep_data, temp_mean, yerr=temp_sdev, marker='o', color=blue, **errorp)
    

    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * 0., color='k', linestyle='-', linewidth=1)
    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * 0.01, color=red, linestyle='--', linewidth=1)
    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * -0.01, color=red, linestyle='--', linewidth=1)

    if plot_in_fm == False:
        ax.set_xlabel(r"$t_{\rm sep}$", **textp)
    elif plot_in_fm == True:
        ax.set_xlabel(r"$t_{\rm sep} / {\rm fm}$", **textp)

    ax.tick_params(axis='both', which='major', **labelp)
    if div_2pt == True:
        ax.text(20 * omega_imp_a09, 0.10, 'tra', fontsize=15)

    ax.set_ylim([-0.15, 0.15])
    ax.set_xlim([0, plot_xmax + 0.8])

    ax.legend(loc='upper right')

    plt.savefig(f"./new_plots/ratio_tra.pdf", transparent=True)
    plt.show()


    ###### sca plot ######
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)

    temp_mean = np.array([val.mean for val in (pt3_gA_sca_fit / pt3_gA_gs_fit)]) 
    temp_sdev = np.array([val.sdev for val in (pt3_gA_sca_fit / pt3_gA_gs_fit)])
    ax.fill_between(pt3_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, color=orange, alpha=0.3, label='3pt sca')

    temp_mean = np.array([val.mean for val in (pt3_gA_sca_data / pt3_gA_gs_data)])
    temp_sdev = np.array([val.sdev for val in (pt3_gA_sca_data / pt3_gA_gs_data)])
    ax.errorbar(pt3_tsep_data, temp_mean, yerr=temp_sdev, marker='o', color=orange, **errorp)

    temp_mean = np.array([val.mean for val in (sum_gA_sca_fit / sum_gA_gs_fit)])
    temp_sdev = np.array([val.sdev for val in (sum_gA_sca_fit / sum_gA_gs_fit)])
    ax.fill_between(sum_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, color=blue, alpha=0.3, label='sum sca')

    temp_mean = np.array([val.mean for val in (sum_gA_sca_data / sum_gA_gs_data)])
    temp_sdev = np.array([val.sdev for val in (sum_gA_sca_data / sum_gA_gs_data)])
    ax.errorbar(sum_tsep_data, temp_mean, yerr=temp_sdev, marker='o', color=blue, **errorp)
    

    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * 0., color='k', linestyle='-', linewidth=1)
    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * 0.01, color=red, linestyle='--', linewidth=1)
    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * -0.01, color=red, linestyle='--', linewidth=1)

    if plot_in_fm == False:
        ax.set_xlabel(r"$t_{\rm sep}$", **textp)
    elif plot_in_fm == True:
        ax.set_xlabel(r"$t_{\rm sep} / {\rm fm}$", **textp)

    ax.tick_params(axis='both', which='major', **labelp)
    if div_2pt == True:
        ax.text(20 * omega_imp_a09, 0.10, 'sca', fontsize=15)

    ax.set_ylim([-0.15, 0.15])
    ax.set_xlim([0, plot_xmax + 0.8])

    ax.legend(loc='upper right')

    plt.savefig(f"./new_plots/ratio_sca.pdf", transparent=True)
    plt.show()




    ###### 3pt/2pt plot ######
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)

    pt2_es_fit = pt3_gA_all_fit - pt3_gA_gs_fit - pt3_gA_tra_fit - pt3_gA_sca_fit
    pt2_es_data = pt3_gA_all_data - pt3_gA_gs_data - pt3_gA_tra_data - pt3_gA_sca_data

    plot_fit_list = [pt3_gA_tra_fit, pt3_gA_sca_fit, pt2_es_fit]
    plot_data_list = [pt3_gA_tra_data, pt3_gA_sca_data, pt2_es_data]
    label_list = ['3pt tra ', '3pt sca ', '2pt es ']
    color_list = [blue, red, grape]

    for i in range(len(plot_fit_list)):
        temp_mean = np.array([val.mean for val in plot_fit_list[i] / pt3_gA_gs_fit])
        temp_sdev = np.array([val.sdev for val in plot_fit_list[i] / pt3_gA_gs_fit])
        ax.fill_between(pt3_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, color=color_list[i], alpha=0.3, label=label_list[i])

        temp_mean = np.array([val.mean for val in plot_data_list[i] / pt3_gA_gs_data])
        temp_sdev = np.array([val.sdev for val in plot_data_list[i] / pt3_gA_gs_data])
        ax.errorbar(pt3_tsep_data, temp_mean, yerr=temp_sdev, marker='o', color=color_list[i], **errorp)

    temp_mean = np.array([val.mean for val in (pt3_gA_all_fit - pt3_gA_gs_fit) / pt3_gA_gs_fit])
    temp_sdev = np.array([val.sdev for val in (pt3_gA_all_fit - pt3_gA_gs_fit) / pt3_gA_gs_fit])
    ax.fill_between(pt3_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, facecolor='None',edgecolor=grey, hatch='/', label='3pt')

    temp_mean = np.array([val.mean for val in (pt3_gA_sca_fit + pt2_es_fit) / pt3_gA_gs_fit])
    temp_sdev = np.array([val.sdev for val in (pt3_gA_sca_fit + pt2_es_fit) / pt3_gA_gs_fit])
    ax.fill_between(pt3_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, facecolor=red,edgecolor=grape, hatch='\\', alpha=0.5, label='2pt+sca')

    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * 0., color='k', linestyle='-', linewidth=1)
    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * 0.01, color=red, linestyle='--', linewidth=1)
    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * -0.01, color=red, linestyle='--', linewidth=1)

    if plot_in_fm == False:
        ax.set_xlabel(r"$t_{\rm sep}$", **textp)
    elif plot_in_fm == True:
        ax.set_xlabel(r"$t_{\rm sep} / {\rm fm}$", **textp)

    ax.tick_params(axis='both', which='major', **labelp)
    if div_2pt == True:
        ax.text(20 * omega_imp_a09, 0.10, '3pt / 2pt', fontsize=15)
    elif div_2pt == False:
        ax.text(20 * omega_imp_a09, 0.10, '3pt', fontsize=15)

    ax.set_ylim([-0.15, 0.15])
    ax.set_xlim([0, plot_xmax + 0.8])

    ax.legend(loc='upper right')

    plt.savefig(f"./new_plots/ratio_pt3_pt2.pdf", transparent=True)
    plt.show()





    ###### sum-sub plot ######
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)

    pt2_es_fit = sum_gA_all_fit - sum_gA_gs_fit - sum_gA_tra_fit - sum_gA_sca_fit
    pt2_es_data = sum_gA_all_data - sum_gA_gs_data - sum_gA_tra_data - sum_gA_sca_data

    plot_fit_list = [sum_gA_tra_fit, sum_gA_sca_fit, pt2_es_fit]
    plot_data_list = [sum_gA_tra_data, sum_gA_sca_data, pt2_es_data]
    label_list = ['sum tra ', 'sum sca ', '2pt es ']
    color_list = [blue, red, grape]

    for i in range(len(plot_fit_list)):
        temp_mean = np.array([val.mean for val in plot_fit_list[i] / sum_gA_gs_fit])
        temp_sdev = np.array([val.sdev for val in plot_fit_list[i] / sum_gA_gs_fit])
        ax.fill_between(sum_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, color=color_list[i], alpha=0.3, label=label_list[i])

        temp_mean = np.array([val.mean for val in plot_data_list[i] / sum_gA_gs_data])
        temp_sdev = np.array([val.sdev for val in plot_data_list[i] / sum_gA_gs_data])
        ax.errorbar(sum_tsep_data, temp_mean, yerr=temp_sdev, marker='o', color=color_list[i], **errorp)

    temp_mean = np.array([val.mean for val in (sum_gA_all_fit - sum_gA_gs_fit) / sum_gA_gs_fit])
    temp_sdev = np.array([val.sdev for val in (sum_gA_all_fit - sum_gA_gs_fit) / sum_gA_gs_fit])
    ax.fill_between(sum_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, facecolor='None',edgecolor=grey, hatch='/', label='sum')

    temp_mean = np.array([val.mean for val in (sum_gA_sca_fit + pt2_es_fit) / sum_gA_gs_fit])
    temp_sdev = np.array([val.sdev for val in (sum_gA_sca_fit + pt2_es_fit) / sum_gA_gs_fit])
    ax.fill_between(sum_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, facecolor=red,edgecolor=grape, hatch='\\', alpha=0.5, label='2pt+sca')

    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * 0., color='k', linestyle='-', linewidth=1)
    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * 0.01, color=red, linestyle='--', linewidth=1)
    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * -0.01, color=red, linestyle='--', linewidth=1)

    if plot_in_fm == False:
        ax.set_xlabel(r"$t_{\rm sep}$", **textp)
    elif plot_in_fm == True:
        ax.set_xlabel(r"$t_{\rm sep} / {\rm fm}$", **textp)

    ax.tick_params(axis='both', which='major', **labelp)
    ax.text(20 * omega_imp_a09, 0.10, 'sum sub', fontsize=15)
    ax.set_ylim([-0.15, 0.15])
    ax.set_xlim([0, plot_xmax + 0.8])


    ax.legend(loc='upper right')

    plt.savefig(f"./new_plots/ratio_sum_sub.pdf", transparent=True)
    plt.show()


    

    ###### compare plot ######
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)


    temp_mean = np.array([val.mean for val in (pt3_gA_all_fit - pt3_gA_gs_fit) / pt3_gA_gs_fit]) 
    temp_sdev = np.array([val.sdev for val in (pt3_gA_all_fit - pt3_gA_gs_fit) / pt3_gA_gs_fit])
    ax.fill_between(pt3_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, color=orange, alpha=0.3, label='3pt')

    temp_mean = np.array([val.mean for val in (pt3_gA_all_data - pt3_gA_gs_data) / pt3_gA_gs_data])
    temp_sdev = np.array([val.sdev for val in (pt3_gA_all_data - pt3_gA_gs_data) / pt3_gA_gs_data])
    ax.errorbar(pt3_tsep_data, temp_mean, yerr=temp_sdev, marker='o', color=orange, **errorp)

    temp_mean = np.array([val.mean for val in (sum_gA_all_fit - sum_gA_gs_fit) / sum_gA_gs_fit])
    temp_sdev = np.array([val.sdev for val in (sum_gA_all_fit - sum_gA_gs_fit) / sum_gA_gs_fit])
    ax.fill_between(sum_tsep_fit, temp_mean + temp_sdev, temp_mean - temp_sdev, color=blue, alpha=0.3, label='sum')

    temp_mean = np.array([val.mean for val in (sum_gA_all_data - sum_gA_gs_data) / sum_gA_gs_data])
    temp_sdev = np.array([val.sdev for val in (sum_gA_all_data - sum_gA_gs_data) / sum_gA_gs_data])
    ax.errorbar(sum_tsep_data, temp_mean, yerr=temp_sdev, marker='o', color=blue, **errorp)

    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * 0., color='k', linestyle='-', linewidth=1)
    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * 0.01, color=red, linestyle='--', linewidth=1)
    ax.plot(pt3_tsep_fit, np.ones(linspace_num) * -0.01, color=red, linestyle='--', linewidth=1)

    if plot_in_fm == False:
        ax.set_xlabel(r"$t_{\rm sep}$", **textp)
    elif plot_in_fm == True:
        ax.set_xlabel(r"$t_{\rm sep} / {\rm fm}$", **textp)

    if div_2pt == True:
        ax.text(20 * omega_imp_a09, 0.08, '3pt / 2pt', fontsize=15)
    elif div_2pt == False:
        ax.text(20 * omega_imp_a09, 0.08, '3pt', fontsize=15)

    ax.tick_params(axis='both', which='major', **labelp)
    ax.text(20 * omega_imp_a09, 0.1, 'tau cut: '+str(sum_tau_cut_plot), fontsize=15)
    # ax.set_ylim([1, 1.4])
    ax.set_ylim([-0.15, 0.15])
    ax.set_xlim([0, plot_xmax + 0.8])

    legend_without_duplicate_labels(ax)

    plt.savefig(f"./new_plots/ratio_compare.pdf", transparent=True)
    plt.show()

    print(sum_gA_gs_fit)
    print(pt3_gA_gs_fit)

# %%
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

#######################################################

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
sum_tau_cut_plot = 1 # this is for plot
fitter_plot = Fit(file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut_plot,  include_2pt, include_3pt, include_sum)

div_2pt = True
plot_in_fm = True

pt3_gA_tra_fit, pt3_gA_sca_fit, pt3_gA_gs_fit, pt3_gA_all_fit, pt3_gA_tra_data, pt3_gA_sca_data, pt3_gA_gs_data, pt3_gA_all_data = pt3_gA_ratio(pt3_data_range, data_avg_dict_completed, fit_result, fitter_plot, div_2pt)

plot_gap, ratio_gA_tra_fit, ratio_gA_sca_fit, ratio_gA_gs_fit, ratio_gA_all_fit, ratio_gA_tra_data, ratio_gA_sca_data, ratio_gA_gs_data, ratio_gA_all_data = sum_gA_ratio(pt3_data_range, data_avg_dict_completed, fit_result, fitter_plot)

ratio_plot(pt3_data_range, plot_gap, pt3_gA_tra_fit, pt3_gA_sca_fit, pt3_gA_gs_fit, pt3_gA_all_fit, pt3_gA_tra_data, pt3_gA_sca_data, pt3_gA_gs_data, pt3_gA_all_data, ratio_gA_tra_fit, ratio_gA_sca_fit, ratio_gA_gs_fit, ratio_gA_all_fit, ratio_gA_tra_data, ratio_gA_sca_data, ratio_gA_gs_data, ratio_gA_all_data, div_2pt, sum_tau_cut_plot, plot_in_fm)







# %% demonstrate that gs+tra+sca = original sum
# p_gs = gv.BufferDict(fit_result.p) # g.s. only
# for key in p_gs:
#     if key.split("_")[0] in ["A3"] and (key.split("_")[1][1] != '0'):
#         p_gs[key] = 0
    
#     if ('z' in key and key != 'z0'):
#         p_gs[key] = 0


# p_sca = gv.BufferDict(fit_result.p) # sca only
# for key in p_sca:
#     if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] != key.split("_")[1][1]):
#         p_sca[key] = 0
# p_sca['A3_00'] = 0

# p_tra = gv.BufferDict(fit_result.p) 
# for key in p_tra:
#     if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] == key.split("_")[1][1]):
#         p_tra[key] = 0

# num = 700

# result_a = fitter_plot.summation_same_can(np.linspace(0, 35, num), np.linspace(0, 35, num), fit_result.p)['sum_A3']

# result_b = fitter_plot.summation_same_can(np.linspace(0, 35, num), np.linspace(0, 35, num), p_gs)['sum_A3'] + fitter_plot.summation_same_can(np.linspace(0, 35, num), np.linspace(0, 35, num), p_sca)['sum_A3'] + fitter_plot.summation_same_can(np.linspace(0, 35, num), np.linspace(0, 35, num), p_tra)['sum_A3']

# for val in (result_a-result_b):
#     if val > 1e-18:
#         print(val)

# gap = int((num-1)/35)+1
# omega_imp_a09 = 0.08730

# pt2_fitter = fitter_plot.pt2_fit_function(np.linspace(0, 35, num), fit_result.p, sum_nstates)['pt2']
# mean = []
# sdev = []
# for i in range(num-gap):
#     temp = result_a[i+gap] / pt2_fitter[i+gap] - result_a[i] / pt2_fitter[i]
#     mean.append(temp.mean)
#     sdev.append(temp.sdev)

# y1 = np.array(mean) + np.array(sdev)
# y2 = np.array(mean) - np.array(sdev)

# plt.figure()
# ax = plt.axes()
# ax.fill_between(np.linspace(0, 35, num)[:-gap]*omega_imp_a09, y1, y2, color=blue, alpha=0.3)
# ax.set_xlim([0, 5])
# ax.set_ylim([1.0, 1.4])
# plt.show()

# %%

# %%
import h5py as h5 
import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  
import lsqfit as lsf
import math
import hashlib

# %%
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

color_list = [grey, fuschia, violet, grape, blue, turquoise, green, lime, yellow, sunkist, orange, peach, red]

marker_list = ["o" for _ in range(10)] #['v', 's', 'o', 'D', '^', 'X', 'P']

figsize = (7, 4)
aspect=[0.15, 0.15, 0.8, 0.8]
#gridspec = {'height_ratios': [3, 1], 'left': 0.05, 'right': 0.95, 'bottom': 0.05, 'top': 0.95}
gridspec_tmin = {'height_ratios': [3, 1, 1], 'left': 0.12, 'right': 0.95, 'bottom': 0.15, 'top': 0.95}
gridspec_tmin_div = {'height_ratios': [3, 1], 'left': 0.12, 'right': 0.95, 'bottom': 0.15, 'top': 0.95}
gridspec_tmax = {'height_ratios': [3, 1], 'left': 0.12, 'right': 0.95, 'bottom': 0.15, 'top': 0.95}
gridspec_prior_width = {'height_ratios': [3, 1], 'left': 0.12, 'right': 0.95, 'bottom': 0.15, 'top': 0.95}
textp = {"fontsize": 14}
labelp = {"labelsize": 14}
errorp = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1}
errorb = {"markersize": 5, "linestyle": "none", "capsize": 3, "elinewidth": 1}

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
t_label = r"$t_{\rm sep}$"
tau_label = r"$\tau - t_{\rm sep}/2$"
fm_t_label = r"$t_{\rm sep} / {\rm fm}$"
fm_tau_label = r"$(\tau - t_{\rm sep}/2) / {\rm fm}$"
meff_label = r"$m_{\textrm{eff}}$"
zeff_label = r"$z_{\textrm{eff}}$"
oaeff_label = r"$O^A_{\textrm{eff}}$"
oveff_label = r"$O^V_{\textrm{eff}}$"

plt.rcParams['figure.figsize'] = figsize

omega_imp_a09 = 0.08730 # converse lattice to fm

# %%
def plot_pt2(pt2_data_range, data_avg_dict, fit_result=None, fitter=None, plot_type=None, plot_in_fm=False):
    '''plot effective mass with 2pt data, you can also plot fit on data'''
    m0_eff = []
    m0_eff_err = []
    m0_eff_fit = []
    m0_eff_fit_err = []
    
    plot_space = 0.05 # smaller, more smoothy

    if plot_in_fm == False:
        x_errorbar = np.arange(pt2_data_range[0], pt2_data_range[1]-1)
        x_fill = np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)[:int(-1/plot_space)]

    elif plot_in_fm == True:
        x_errorbar = np.arange(pt2_data_range[0], pt2_data_range[1]-1) * omega_imp_a09 
        x_fill = np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)[:int(-1/plot_space)] * omega_imp_a09 

    else:
        x_errorbar = np.arange(pt2_data_range[0], pt2_data_range[1]-1)
        x_fill = np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)[:int(-1/plot_space)]
        print("Input error: plot_in_fm")

    for i in range(pt2_data_range[0], pt2_data_range[1]-1):
        temp = gv.log(data_avg_dict['pt2_tsep_'+str(i)] / data_avg_dict['pt2_tsep_'+str(i+1)])
        m0_eff.append(temp.mean)
        m0_eff_err.append(temp.sdev)

    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    ax.errorbar(x_errorbar, np.array(m0_eff), yerr=np.array(m0_eff_err), marker='o', color="k", **errorp)
    
    if fit_result != None and fitter != None:
        pt2_fitter = fitter.pt2_fit_function(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space), fit_result.p)['pt2']
        
        for i in range(len(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)) - int(1/plot_space)):
            temp = gv.log(pt2_fitter[i] / pt2_fitter[i+ int(1/plot_space)])
            m0_eff_fit.append(temp.mean)
            m0_eff_fit_err.append(temp.sdev)
            
        m0_eff_fit_y1 = []
        m0_eff_fit_y2 = []

        for i in range(len(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)) - int(1/plot_space)):
            m0_eff_fit_y1.append(m0_eff_fit[i] + m0_eff_fit_err[i])
            m0_eff_fit_y2.append(m0_eff_fit[i] - m0_eff_fit_err[i])
        ax.fill_between(x_fill, np.array(m0_eff_fit_y1), np.array(m0_eff_fit_y2), color=blue, alpha=0.3, label='fit')
            
    x_lim = [2, 25]
    if plot_in_fm == False:
        ax.set_xlim(x_lim)
        ax.set_xlabel(t_label, **textp)
    elif plot_in_fm == True:
        ax.set_xlim([num*omega_imp_a09 for num in x_lim])
        ax.set_xlabel(fm_t_label, **textp)

    ax.set_ylim([0.45, 0.62])
    ax.set_ylabel(meff_label, **textp)        
    
    ax.tick_params(axis='both', which='major', **labelp)

    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/meff{plot_type}.pdf", transparent=True)
    plt.show()
    
    
    #zeff
    zeff = []

    for i in range(pt2_data_range[0], pt2_data_range[1]-1):
        meff = gv.log(data_avg_dict['pt2_tsep_'+str(i)] / data_avg_dict['pt2_tsep_'+str(i+1)])
        zeff.append(data_avg_dict['pt2_tsep_'+str(i)]*np.exp(meff*i))
    
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    ax.errorbar(x_errorbar, [i.mean for i in zeff], yerr=[i.sdev for i in zeff], marker='o', color="k", **errorp)
 
    if fit_result != None and fitter != None:
        z0_eff_fit = []
        z0_eff_fit_err = []
    
        x = np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)[:int(-1/plot_space)]
        pt2_fitter = fitter.pt2_fit_function(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space), fit_result.p)['pt2']
        
        for idx, i in enumerate(range(len(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)) - int(1/plot_space))):
            meff = gv.log(pt2_fitter[i] / pt2_fitter[i+ int(1/plot_space)])
            zeff = pt2_fitter[i] * np.exp(meff*x[idx])
            z0_eff_fit.append(zeff.mean)
            z0_eff_fit_err.append(zeff.sdev)
        
        z0_eff_fit_y1 = []
        z0_eff_fit_y2 = []

        for i in range(len(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)) - int(1/plot_space)):
            z0_eff_fit_y1.append(z0_eff_fit[i] + z0_eff_fit_err[i])
            z0_eff_fit_y2.append(z0_eff_fit[i] - z0_eff_fit_err[i])
        
        ax.fill_between(x_fill, np.array(z0_eff_fit_y1), np.array(z0_eff_fit_y2), color=blue, alpha=0.3, label='fit')

    x_lim = [2, 25]
    if plot_in_fm == False:
        ax.set_xlim(x_lim)
        ax.set_xlabel(t_label, **textp)
    elif plot_in_fm == True:
        ax.set_xlim([num*omega_imp_a09 for num in x_lim])
        ax.set_xlabel(fm_t_label, **textp)

    ax.set_ylim([0, 3E-7])
    ax.set_ylabel(zeff_label, **textp)


    ax.tick_params(axis='both', which='major', **labelp)
    
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/zeff{plot_type}.pdf", transparent=True)
    plt.show()

    # return [np.arange(pt2_data_range[0], pt2_data_range[1]-1), np.array(m0_eff), np.array(m0_eff_err)] # save this for data plot of half stat

def plot_pt3(pt3_data_range, data_avg_dict_completed, tau_cut, fit_result=None, fitter=None, plot_type=None, plot_in_fm=False):
    '''plot form factor with 3pt data, you can also plot fit on data'''
    gA = {}
    gA_err = {}
    gV = {}
    gV_err = {}
    
    gA_fit = {}
    gA_fit_err = {}
    gV_fit = {}
    gV_fit_err = {}
    
    gA_tsep = []
    gA_tau = []
    gV_tsep = []
    gV_tau = []
    
    plot_density = 100 # bigger, more smoothy
    
    for i in range(pt3_data_range[0], pt3_data_range[1]): ##
            for j in np.linspace(tau_cut, i-tau_cut, plot_density):
                gA_tsep.append(i)
                gA_tau.append(j)
                
    for i in range(pt3_data_range[0], pt3_data_range[1]): ##
            for j in np.linspace(tau_cut, i-tau_cut, plot_density):
                gV_tsep.append(i)
                gV_tau.append(j)
                
    for i in range(pt3_data_range[0], pt3_data_range[1]): ##
        gA['tsep_'+str(i)] = []
        gA_err['tsep_'+str(i)] = []
        gA_fit['tsep_'+str(i)] = []
        gA_fit_err['tsep_'+str(i)] = []
        gV['tsep_'+str(i)] = []
        gV_err['tsep_'+str(i)] = []
        gV_fit['tsep_'+str(i)] = []
        gV_fit_err['tsep_'+str(i)] = []
        
    for i in range(pt3_data_range[0], pt3_data_range[1]):
        for j in range(tau_cut, i-tau_cut+1):
            temp1 = ( data_avg_dict_completed['pt3_A3_tsep_'+str(i)][j] + data_avg_dict_completed['pt3_A3_tsep_'+str(i)][i-j] ) / (2 * data_avg_dict_completed['pt2_tsep_'+str(i)])
            
            gA['tsep_'+str(i)].append(temp1.mean)
            gA_err['tsep_'+str(i)].append(temp1.sdev)
            
            temp2 = ( data_avg_dict_completed['pt3_V4_tsep_'+str(i)][j] + data_avg_dict_completed['pt3_V4_tsep_'+str(i)][i-j] ) / (2 * data_avg_dict_completed['pt2_tsep_'+str(i)])
            
            gV['tsep_'+str(i)].append(temp2.mean)
            gV_err['tsep_'+str(i)].append(temp2.sdev)
          
    if fit_result != None and fitter != None:
        pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1]), fit_result.p)['pt2']
        pt3_gA_fitter = fitter.pt3_fit_function(np.array(gA_tsep), np.array(gV_tsep), np.array(gA_tau), np.array(gV_tau), fit_result.p)['pt3_A3']
        pt3_gV_fitter = fitter.pt3_fit_function(np.array(gA_tsep), np.array(gV_tsep), np.array(gA_tau), np.array(gV_tau), fit_result.p)['pt3_V4']
        index = 0
        for i in range(pt3_data_range[0], pt3_data_range[1]):
            for j in range(plot_density):
                index = int((i-pt3_data_range[0])*plot_density + j)
                
                temp1 = pt3_gA_fitter[index] / pt2_fitter[i - pt3_data_range[0]]
                gA_fit['tsep_'+str(i)].append(temp1.mean)
                gA_fit_err['tsep_'+str(i)].append(temp1.sdev)
                
                temp2 = pt3_gV_fitter[index] / pt2_fitter[i - pt3_data_range[0]]
                gV_fit['tsep_'+str(i)].append(temp2.mean)
                gV_fit_err['tsep_'+str(i)].append(temp2.sdev)
            
        
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    for idx, i in enumerate(range(pt3_data_range[0], pt3_data_range[1])):
        if plot_in_fm == False:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1)
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density)

        elif plot_in_fm == True:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1) * omega_imp_a09
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density) * omega_imp_a09

        else:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1)
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density)
            print("Input error: plot_in_fm")

        color = color_list[idx]
        ax.errorbar(x_errorbar, np.array(gA['tsep_' + str(i)]), yerr=np.array(gA_err['tsep_' + str(i)]), marker='o', color=color, **errorp)# tau cut = 1

        if fit_result != None and fitter != None:
            gA_fit_y1 = np.array(gA_fit['tsep_'+str(i)]) + np.array(gA_fit_err['tsep_'+str(i)])
            gA_fit_y2 = np.array(gA_fit['tsep_'+str(i)]) - np.array(gA_fit_err['tsep_'+str(i)])
            ax.fill_between(x_fill, gA_fit_y1, gA_fit_y2, color=color, alpha=0.3) # tau_cut=1

    x_lim = [-6.5, 6.5]
    if plot_in_fm == False:
        ax.set_xlim(x_lim)
        ax.set_xlabel(tau_label, **textp)
    elif plot_in_fm == True:
        ax.set_xlim([num*omega_imp_a09 for num in x_lim])
        ax.set_xlabel(fm_tau_label, **textp)

    ax.set_ylim([1, 1.3])
    ax.set_ylabel(oaeff_label, **textp)
    #
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/oaeff{plot_type}.pdf", transparent=True)
    plt.show()
    
    
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    for idx, i in enumerate(range(pt3_data_range[0], pt3_data_range[1])):
        if plot_in_fm == False:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1)
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density)

        elif plot_in_fm == True:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1) * omega_imp_a09
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density) * omega_imp_a09

        else:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1)
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density)
            print("Input error: plot_in_fm")

        color = color_list[idx]
        ax.errorbar(x_errorbar, np.array(gV['tsep_' + str(i)]), yerr=np.array(gV_err['tsep_' + str(i)]), marker='o', color=color, **errorp)# tau cut = 1

        if fit_result != None and fitter != None:
            gV_fit_y1 = np.array(gV_fit['tsep_'+str(i)]) + np.array(gV_fit_err['tsep_'+str(i)])
            gV_fit_y2 = np.array(gV_fit['tsep_'+str(i)]) - np.array(gV_fit_err['tsep_'+str(i)])
            ax.fill_between(x_fill, gV_fit_y1, gV_fit_y2, color=color, alpha=0.3) # tau_cut=1

    x_lim = [-6.5, 6.5]
    if plot_in_fm == False:
        ax.set_xlim(x_lim)
        ax.set_xlabel(tau_label, **textp)
    elif plot_in_fm == True:
        ax.set_xlim([num*omega_imp_a09 for num in x_lim])
        ax.set_xlabel(fm_tau_label, **textp)
    
    ax.set_ylim([1.0, 1.15])
    ax.set_ylabel(oveff_label, **textp)
    #
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/oveff{plot_type}.pdf", transparent=True)
    plt.show()

def plot_pt3_no_tra(pt3_data_range, data_avg_dict_completed, tau_cut, fit_result, fitter, plot_type, plot_in_fm=False):
    '''plot form factor with 3pt data, you can also plot fit on data'''    
    gA_fit = {} # no transition
    gA_fit_err = {}

    gA_nsca = {}
    gA_ntra = {} # data - no scattering
    
    gA_tsep = []
    gA_tau = []
    gV_tsep = []
    gV_tau = []
    
    plot_density = 100 # bigger, more smoothy
    
    for i in range(pt3_data_range[0], pt3_data_range[1]): 
            for j in np.linspace(tau_cut, i-tau_cut, plot_density):
                gA_tsep.append(i)
                gA_tau.append(j)
                gV_tsep.append(i)
                gV_tau.append(j)
                
                
    for i in range(pt3_data_range[0], pt3_data_range[1]): 
        gA_fit['tsep_'+str(i)] = []
        gA_fit_err['tsep_'+str(i)] = []

        gA_nsca['tsep_'+str(i)] = []
          
    p_ntra = gv.BufferDict(fit_result.p) # no transition
    for key in p_ntra:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] != key.split("_")[1][1]):
            p_ntra[key] = 0

    p_nsca = gv.BufferDict(fit_result.p) # no scattering
    for key in p_nsca:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] == key.split("_")[1][1]):
            p_nsca[key] = 0
                
######################################## no transition
    pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1]), fit_result.p)['pt2']
    pt3_gA_fitter = fitter.pt3_fit_function(np.array(gA_tsep), np.array(gV_tsep), np.array(gA_tau), np.array(gV_tau), p_ntra)['pt3_A3']

    index = 0
    for i in range(pt3_data_range[0], pt3_data_range[1]):
        for j in range(plot_density):
            index = int((i-pt3_data_range[0])*plot_density + j)
            
            temp1 = pt3_gA_fitter[index] / pt2_fitter[i - pt3_data_range[0]]
            gA_fit['tsep_'+str(i)].append(temp1.mean)
            gA_fit_err['tsep_'+str(i)].append(temp1.sdev)
        
######################################### data - no scattering
    gA_tsep_nsca = [] # errorbar do not need plot density
    gA_tau_nsca = []      

    for i in range(pt3_data_range[0], pt3_data_range[1]): # tsep and tau to generating fitter of no scattering 
        for j in range(tau_cut, i-tau_cut+1):
            gA_tsep_nsca.append(i)
            gA_tau_nsca.append(j)

    pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1]), fit_result.p)['pt2']
    pt3_gA_fitter = fitter.pt3_fit_function(np.array(gA_tsep_nsca), np.array(gV_tsep), np.array(gA_tau_nsca), np.array(gV_tau), p_nsca)['pt3_A3'] # fitter of no scattering

    index = 0
    for i in range(pt3_data_range[0], pt3_data_range[1]):
        for j in range(tau_cut, i-tau_cut+1):
            index = int((pt3_data_range[0]+i+1-4*tau_cut)*(i-pt3_data_range[0])/2 + j - 1)
            
            temp1 = pt3_gA_fitter[index] / pt2_fitter[i - pt3_data_range[0]]
            gA_nsca['tsep_'+str(i)].append(temp1)

        gA_ntra['tsep_'+str(i)] = np.array( ( data_avg_dict_completed['pt3_A3_tsep_'+str(i)][1:-1] ) / (data_avg_dict_completed['pt2_tsep_'+str(i)] ) ) - np.array(gA_nsca['tsep_'+str(i)]) # ntra = data - nsca

    mean_demo = []
    sdev_demo = []
    for i in range(pt3_data_range[0], pt3_data_range[1], 2): # only even tsep has tau = tsep/2
        mean_demo.append(gA_ntra['tsep_'+str(i)][int(i/2-1)].mean)
        sdev_demo.append(gA_ntra['tsep_'+str(i)][int(i/2-1)].sdev)

    print("########################### data - no scattering ########################"+plot_type)
    print(mean_demo)
    print("###")
    print(sdev_demo)
#########################################

    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    for idx, i in enumerate(range(pt3_data_range[0], pt3_data_range[1])):
        if plot_in_fm == False:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1)
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density)

        elif plot_in_fm == True:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1) * omega_imp_a09
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density) * omega_imp_a09

        else:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1)
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density)
            print("Input error: plot_in_fm")

        color = color_list[idx]
        gA_fit_y1 = np.array(gA_fit['tsep_'+str(i)]) + np.array(gA_fit_err['tsep_'+str(i)])
        gA_fit_y2 = np.array(gA_fit['tsep_'+str(i)]) - np.array(gA_fit_err['tsep_'+str(i)])
        # no transition
        ax.fill_between(x_fill, gA_fit_y1, gA_fit_y2, color=color, alpha=0.3) # tau_cut=1              

        gA_ntra_mean = np.array([value.mean for value in gA_ntra['tsep_'+str(i)]])
        gA_ntra_sdev = np.array([value.sdev for value in gA_ntra['tsep_'+str(i)]])
        # data - no scattering
        ax.errorbar(x_errorbar, gA_ntra_mean, yerr=gA_ntra_sdev, marker='x', color=color, **errorp)
    
    x_lim = [-6.5, 6.5]
    if plot_in_fm == False:
        ax.set_xlim(x_lim)
        ax.set_xlabel(tau_label, **textp)
    elif plot_in_fm == True:
        ax.set_xlim([num*omega_imp_a09 for num in x_lim])
        ax.set_xlabel(fm_tau_label, **textp)

    ax.set_ylim([1, 1.3])
    ax.set_ylabel(oaeff_label, **textp)
    ax.tick_params(axis='both', which='major', **labelp)
    #
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/oaeff{plot_type}_ntra.pdf", transparent=True)
    plt.show()

def plot_pt3_no_sca(pt3_data_range, data_avg_dict_completed, tau_cut, fit_result, fitter, plot_type, plot_in_fm=False):
    '''plot form factor with 3pt data, you can also plot fit on data'''    
    gA_fit = {} 
    gA_fit_err = {} # no scattering except for g.s.

    gA_nsca = {}
    gA_ntra = {} # data - no transition and no g.s.
    
    gA_tsep = []
    gA_tau = []
    gV_tsep = []
    gV_tau = []
    
    plot_density = 100 # bigger, more smoothy
    
    for i in range(pt3_data_range[0], pt3_data_range[1]): 
            for j in np.linspace(tau_cut, i-tau_cut, plot_density):
                gA_tsep.append(i)
                gA_tau.append(j)
                gV_tsep.append(i)
                gV_tau.append(j)
                
                
    for i in range(pt3_data_range[0], pt3_data_range[1]): 
        gA_fit['tsep_'+str(i)] = []
        gA_fit_err['tsep_'+str(i)] = []

        gA_ntra['tsep_'+str(i)] = []

    p_nsca = gv.BufferDict(fit_result.p) # no scattering except for g.s.
    for key in p_nsca:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] == key.split("_")[1][1]) and key != 'A3_00':
            p_nsca[key] = 0

    p_ntra = gv.BufferDict(fit_result.p) # no transition and no g.s.
    for key in p_ntra:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] != key.split("_")[1][1]):
            p_ntra[key] = 0

        p_ntra['A3_00'] = 0
                
######################################## no scattering except for g.s.
    pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1]), fit_result.p)['pt2']
    pt3_gA_fitter = fitter.pt3_fit_function(np.array(gA_tsep), np.array(gV_tsep), np.array(gA_tau), np.array(gV_tau), p_nsca)['pt3_A3']

    index = 0
    for i in range(pt3_data_range[0], pt3_data_range[1]):
        for j in range(plot_density):
            index = int((i-pt3_data_range[0])*plot_density + j)
            
            temp1 = pt3_gA_fitter[index] / pt2_fitter[i - pt3_data_range[0]]
            gA_fit['tsep_'+str(i)].append(temp1.mean)
            gA_fit_err['tsep_'+str(i)].append(temp1.sdev)
        
######################################### data - no transition and no g.s.
    gA_tsep_ntra = [] # errorbar do not need plot density
    gA_tau_ntra = []      

    for i in range(pt3_data_range[0], pt3_data_range[1]): # tsep and tau to generating fitter of no transition and no g.s.
        for j in range(tau_cut, i-tau_cut+1):
            gA_tsep_ntra.append(i)
            gA_tau_ntra.append(j)

    pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1]), fit_result.p)['pt2']
    pt3_gA_fitter = fitter.pt3_fit_function(np.array(gA_tsep_ntra), np.array(gV_tsep), np.array(gA_tau_ntra), np.array(gV_tau), p_ntra)['pt3_A3'] # fitter of no transition and no g.s.

    index = 0
    for i in range(pt3_data_range[0], pt3_data_range[1]):
        for j in range(tau_cut, i-tau_cut+1):
            index = int((pt3_data_range[0]+i+1-4*tau_cut)*(i-pt3_data_range[0])/2 + j - 1)
            
            temp1 = pt3_gA_fitter[index] / pt2_fitter[i - pt3_data_range[0]]
            gA_ntra['tsep_'+str(i)].append(temp1)

        gA_nsca['tsep_'+str(i)] = np.array( ( data_avg_dict_completed['pt3_A3_tsep_'+str(i)][1:-1] ) / (data_avg_dict_completed['pt2_tsep_'+str(i)] ) ) - np.array(gA_ntra['tsep_'+str(i)]) # nsca = data - ntra

    mean_demo = []
    sdev_demo = []
    for i in range(pt3_data_range[0], pt3_data_range[1], 2): # only even tsep has tau = tsep/2
        mean_demo.append(gA_nsca['tsep_'+str(i)][int(i/2-1)].mean)
        sdev_demo.append(gA_nsca['tsep_'+str(i)][int(i/2-1)].sdev)

    print("########################## data - no transition and no g.s. ################################"+plot_type)
    print(mean_demo)
    print("###")
    print(sdev_demo)
#########################################

    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    for idx, i in enumerate(range(pt3_data_range[0], pt3_data_range[1])):
        if plot_in_fm == False:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1)
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density)

        elif plot_in_fm == True:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1) * omega_imp_a09
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density) * omega_imp_a09

        else:
            x_errorbar = np.arange(tau_cut - i/2, i/2 - tau_cut + 1)
            x_fill = np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density)
            print("Input error: plot_in_fm")

        color = color_list[idx]
        gA_fit_y1 = np.array(gA_fit['tsep_'+str(i)]) + np.array(gA_fit_err['tsep_'+str(i)])
        gA_fit_y2 = np.array(gA_fit['tsep_'+str(i)]) - np.array(gA_fit_err['tsep_'+str(i)])
        # no transition
        ax.fill_between(x_fill, gA_fit_y1, gA_fit_y2, color=color, alpha=0.3) # tau_cut=1              

        gA_nsca_mean = np.array([value.mean for value in gA_nsca['tsep_'+str(i)]])
        gA_nsca_sdev = np.array([value.sdev for value in gA_nsca['tsep_'+str(i)]])
        # data - no scattering
        ax.errorbar(x_errorbar, gA_nsca_mean, yerr=gA_nsca_sdev, marker='x', color=color, **errorp)

    x_lim = [-6.5, 6.5]
    if plot_in_fm == False:
        ax.set_xlim(x_lim)
        ax.set_xlabel(tau_label, **textp)
    elif plot_in_fm == True:
        ax.set_xlim([num*omega_imp_a09 for num in x_lim])
        ax.set_xlabel(fm_tau_label, **textp)

    ax.set_ylim([1, 1.3])
    ax.set_ylabel(oaeff_label, **textp)
    ax.tick_params(axis='both', which='major', **labelp)
    #
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/oaeff{plot_type}_nsca.pdf", transparent=True)
    plt.show()

def plot_sum(pt3_data_range, data_avg_dict_completed, fit_result=None, fitter=None, pt2_nstates=None, sum_nstates=None, plot_type=None, plot_in_fm=False):
    '''plot form factor with sum data, you can also plot fit on data'''
    gA = []
    gA_err = []
    gV = []
    gV_err = []
    
    gA_fit = []
    gA_fit_err = []
    gV_fit = []
    gV_fit_err = []
    
    plot_space = 0.05

    if plot_in_fm == False:
        x_errorbar = np.arange(pt3_data_range[0], pt3_data_range[1]-1)
        x_fill = np.arange(pt3_data_range[0], pt3_data_range[1], plot_space)[:int(-1 / plot_space)]

    elif plot_in_fm == True:
        x_errorbar = np.arange(pt3_data_range[0], pt3_data_range[1]-1) * omega_imp_a09
        x_fill = np.arange(pt3_data_range[0], pt3_data_range[1], plot_space)[:int(-1 / plot_space)] * omega_imp_a09

    else:
        x_errorbar = np.arange(pt3_data_range[0], pt3_data_range[1]-1)
        x_fill = np.arange(pt3_data_range[0], pt3_data_range[1], plot_space)[:int(-1 / plot_space)]
        print("Input error: plot_in_fm")

    for i in range(pt3_data_range[0], pt3_data_range[1]-1):
        temp1 = data_avg_dict_completed['sum_A3_fit_'+str(i)]
        gA.append(temp1.mean)
        gA_err.append(temp1.sdev)
        
        temp2 = data_avg_dict_completed['sum_V4_fit_'+str(i)]
        gV.append(temp2.mean)
        gV_err.append(temp2.sdev)

    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    ax.errorbar(x_errorbar, np.array(gA), yerr=np.array(gA_err), marker='o', color="k", **errorp)
    # print(gA)
    # print(gA_err)
    if fit_result != None and fitter != None and sum_nstates != None:
        if sum_nstates == pt2_nstates:
            pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), fit_result.p)['pt2']
            sum_A3_fitter = fitter.summation_same_can(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), fit_result.p)['sum_A3']
            sum_V4_fitter = fitter.summation_same_can(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), fit_result.p)['sum_V4']
            
        else:
            pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), fit_result.p, sum_nstates)['pt2']
            sum_A3_fitter = fitter.summation(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), fit_result.p)['sum_A3']
            sum_V4_fitter = fitter.summation(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), np.arange(pt3_data_range[0], pt3_data_range[1], plot_space), fit_result.p)['sum_V4']
        
        for i in range(len(np.arange(pt3_data_range[0], pt3_data_range[1], plot_space)) - int(1 / plot_space)):
            temp1 = sum_A3_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i]
            gA_fit.append(temp1.mean)
            gA_fit_err.append(temp1.sdev)

            temp2 = sum_V4_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_V4_fitter[i] / pt2_fitter[i]
            gV_fit.append(temp2.mean)
            gV_fit_err.append(temp2.sdev)
            
        gA_fit_y1 = np.array(gA_fit) + np.array(gA_fit_err)
        gA_fit_y2 = np.array(gA_fit) - np.array(gA_fit_err)
        
        ax.fill_between(x_fill, gA_fit_y1, gA_fit_y2, color=blue, alpha=0.3, label='fit')
    

    x_lim = [1.5, 13.5]
    if plot_in_fm == False:
        ax.set_xlim(x_lim)
        ax.set_xlabel(t_label, **textp)
    elif plot_in_fm == True:
        ax.set_xlim([num*omega_imp_a09 for num in x_lim])
        ax.set_xlabel(fm_t_label, **textp)

    ax.set_ylim([1, 1.4])
    ax.set_ylabel(oaeff_label, **textp)
    
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/soaeff{plot_type}.pdf", transparent=True)
    plt.show()
    
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    ax.errorbar(x_errorbar, np.array(gV), yerr=np.array(gV_err), marker='o', color="k", **errorp)
    # print(gV)
    # print(gV_err)
    if fit_result != None and fitter != None and sum_nstates != None:
        gV_fit_y1 = np.array(gV_fit) + np.array(gV_fit_err)
        gV_fit_y2 = np.array(gV_fit) - np.array(gV_fit_err)
        
        ax.fill_between(x_fill, gV_fit_y1, gV_fit_y2, color=blue, alpha=0.3)
    
    
    x_lim = [1.5, 13.5]
    if plot_in_fm == False:
        ax.set_xlim(x_lim)
        ax.set_xlabel(t_label, **textp)
    elif plot_in_fm == True:
        ax.set_xlim([num*omega_imp_a09 for num in x_lim])
        ax.set_xlabel(fm_t_label, **textp)

    ax.set_ylim([1.0, 1.1])
    ax.set_ylabel(oveff_label, **textp)
    
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/soveff{plot_type}.pdf", transparent=True)
    plt.show()

    #print(gA[8])

    # return [np.arange(pt3_data_range[0], pt3_data_range[1]-1), np.array(gA), np.array(gA_err)] # save this for data plot of half stat
        
def plot_sum_no_tra(pt3_data_range, data_avg_dict_completed, fit_result, fitter, pt2_nstates, sum_nstates, plot_type, plot_in_fm=False):
    '''plot form factor with sum data, you can also plot fit on data'''
    gA = []
    gA_err = []
    
    gA_fit = []
    gA_fit_err = []
    gV_fit = []
    gV_fit_err = []
    
    plot_space = 0.05

    plt_min = 0
    plt_max = 60 
    fillx = np.arange(plt_min, plt_max, plot_space)[:int(-1 / plot_space)] # as tsep input

    if plot_in_fm == False:
        x_errorbar = np.arange(pt3_data_range[0], pt3_data_range[1]-1)
        x_fill = np.arange(plt_min, plt_max, plot_space)[:int(-1 / plot_space)]

    elif plot_in_fm == True:
        x_errorbar = np.arange(pt3_data_range[0], pt3_data_range[1]-1) * omega_imp_a09
        x_fill = np.arange(plt_min, plt_max, plot_space)[:int(-1 / plot_space)] * omega_imp_a09

    else:
        x_errorbar = np.arange(pt3_data_range[0], pt3_data_range[1]-1)
        x_fill = np.arange(plt_min, plt_max, plot_space)[:int(-1 / plot_space)]
        print("Input error: plot_in_fm")

    for i in range(pt3_data_range[0], pt3_data_range[1]-1):
        temp1 = data_avg_dict_completed['sum_A3_fit_'+str(i)]
        gA.append(temp1.mean)
        gA_err.append(temp1.sdev)

##################### fit on data
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    ax.errorbar(x_errorbar, np.array(gA), yerr=np.array(gA_err), marker='o', color=blue, **errorp)

    if sum_nstates == pt2_nstates:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), fit_result.p)['pt2']
        sum_A3_fitter = fitter.summation_same_can(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_A3']
        sum_V4_fitter = fitter.summation_same_can(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_V4']
        
    else:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), fit_result.p, sum_nstates)['pt2']
        sum_A3_fitter = fitter.summation(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_A3']
        sum_V4_fitter = fitter.summation(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_V4']
    
    for i in range(len(np.arange(plt_min, plt_max, plot_space)) - int(1 / plot_space)):
        temp1 = sum_A3_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i]
        gA_fit.append(temp1.mean)
        gA_fit_err.append(temp1.sdev)

        temp2 = sum_V4_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_V4_fitter[i] / pt2_fitter[i]
        gV_fit.append(temp2.mean)
        gV_fit_err.append(temp2.sdev)
        
    gA_fit_y1 = np.array(gA_fit) + np.array(gA_fit_err)
    gA_fit_y2 = np.array(gA_fit) - np.array(gA_fit_err)
    
    ax.fill_between(x_fill, gA_fit_y1, gA_fit_y2, color=blue, alpha=0.3, label='fit')

######################################## 
    p_ntra = gv.BufferDict(fit_result.p) # no transition
    for key in p_ntra:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] != key.split("_")[1][1]):
            p_ntra[key] = 0

    p_nsca = gv.BufferDict(fit_result.p) # no scattering
    for key in p_nsca:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] == key.split("_")[1][1]):
            p_nsca[key] = 0

    p_gs = gv.BufferDict(fit_result.p) # g.s. only
    for key in p_gs:
        if key.split("_")[0] in ["A3"] and key != 'A3_00':
            p_gs[key] = 0

    p_gs_pt2 = gv.BufferDict(fit_result.p) # g.s. only for 2pt
    for key in p_gs_pt2:
        if 'z'in key and key != 'z0' and key != 'log(z2)':
            p_gs_pt2[key] = 0

    del p_gs_pt2['log(z2)']
    p_gs_pt2['z2'] = 0

##################### original gA from 3pt with tau = tsep/2
    y3 = np.array([1.137246887069871, 1.0397951114352422, 1.0581101894217029, 1.0944768415815544, 1.129653822067308, 1.1607318494702763, 1.1883613569687608])
    dy3 = np.array([0.0006817408659827335, 0.0008354900714438219, 0.0011941822664530966, 0.0019041344115398975, 0.003132857777381124, 0.005059734516954951, 0.008446467991949032])
    x3 = np.arange(2, 15, 2)

    if plot_in_fm == True:
        x3 = x3 * omega_imp_a09

    ax.errorbar(x3, y3, yerr=dy3, ls='none', marker='o', capsize=3, color=orange) 

    pt2_fitter = fitter.pt2_fit_function(fillx, fit_result.p)['pt2']
    pt3_gA_fitter = fitter.pt3_fit_function(fillx, fillx, fillx/2, fillx/2, fit_result.p)['pt3_A3'] # here fillx/2 means tau = tsep/2

    demo_ori = pt3_gA_fitter/pt2_fitter

    ax.fill_between(x_fill, np.array([y.mean-y.sdev for y in demo_ori]), np.array([y.mean+y.sdev for y in demo_ori]), alpha=0.2, color=orange)


##################### data - no scattering gA from 3pt with tau = tsep/2
    y3_ntra = np.array([1.0168142241834535, 1.1330664811954303, 1.189630793941224, 1.218697582268415, 1.233856395733631, 1.2438491219561556, 1.2531474983332445])

    dy3_ntra = np.array([0.12563164492523324, 0.07344246748838909, 0.04554307000176309, 0.03209866967751453, 0.025578449043411906, 0.022593798888500174, 0.021829703982068273])

    x3_ntra = np.arange(2, 15, 2)

    if plot_in_fm == True:
        x3_ntra = x3_ntra * omega_imp_a09

    #ax.errorbar(x3_ntra, y3_ntra, yerr=dy3_ntra, ls='none', marker='x', capsize=3, color=yellow) 

    pt2_fitter = fitter.pt2_fit_function(fillx, fit_result.p)['pt2']
    pt3_gA_fitter = fitter.pt3_fit_function(fillx, fillx, fillx/2, fillx/2, p_ntra)['pt3_A3']  # here fillx/2 means tau = tsep/2

    demo_ntra = pt3_gA_fitter/pt2_fitter

    #ax.fill_between(x_fill, np.array([y.mean-y.sdev for y in demo_ntra]), np.array([y.mean+y.sdev for y in demo_ntra]), alpha=0.2, color=yellow)


###################### no transition sum-sub
    gA_fit_ntra = []
    gA_fit_err_ntra = []

    if sum_nstates == pt2_nstates:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2)['pt2']
        sum_A3_fitter = fitter.summation_same_can(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), p_ntra)['sum_A3']
        
    else:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2, sum_nstates)['pt2']
        sum_A3_fitter = fitter.summation(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), p_ntra)['sum_A3']
    
    for i in range(len(np.arange(plt_min, plt_max, plot_space)) - int(1 / plot_space)):
        temp1 = sum_A3_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i]
        gA_fit_ntra.append(temp1.mean)
        gA_fit_err_ntra.append(temp1.sdev)

    ###
    print(p_ntra)
        
    gA_fit_y1_ntra = np.array(gA_fit_ntra) + np.array(gA_fit_err_ntra)
    gA_fit_y2_ntra = np.array(gA_fit_ntra) - np.array(gA_fit_err_ntra)
    
    #ax.fill_between(x_fill, gA_fit_y1_ntra, gA_fit_y2_ntra, color=green, alpha=0.3)

######################## data - no scattering sum-sub
    gA_fit_ntra = []

    if sum_nstates == pt2_nstates:
        pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1]), p_gs_pt2)['pt2']
        sum_A3_fitter = fitter.summation_same_can(np.arange(pt3_data_range[0], pt3_data_range[1]), np.arange(pt3_data_range[0], pt3_data_range[1]), p_nsca)['sum_A3']
        
    else:
        pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1]), p_gs_pt2, sum_nstates)['pt2']
        sum_A3_fitter = fitter.summation(np.arange(pt3_data_range[0], pt3_data_range[1]), np.arange(pt3_data_range[0], pt3_data_range[1]), p_nsca)['sum_A3']
    
    for i in range(pt3_data_range[0], pt3_data_range[1]-1):
        index = i - pt3_data_range[0]
        temp1 = data_avg_dict_completed['sum_A3_fit_'+str(i)] - (sum_A3_fitter[index+1] / pt2_fitter[index+1] - sum_A3_fitter[index] / pt2_fitter[index])
        gA_fit_ntra.append(temp1) # ntra = data - nsca

    #ax.errorbar(x_errorbar, np.array([val.mean for val in gA_fit_ntra]), yerr=np.array([val.sdev for val in gA_fit_ntra]), ls='none', marker='x', capsize=3, color=green)

######################
###################### g.s. only
    gA_fit_gs = []
    gA_fit_err_gs = []

    if sum_nstates == pt2_nstates:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2)['pt2']
        sum_A3_fitter = fitter.summation_same_can(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), p_gs)['sum_A3']
        
    else:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2, sum_nstates)['pt2']
        sum_A3_fitter = fitter.summation(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), p_gs)['sum_A3']
    
    for i in range(len(np.arange(plt_min, plt_max, plot_space)) - int(1 / plot_space)):
        temp1 = sum_A3_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i]
        gA_fit_gs.append(temp1.mean)
        gA_fit_err_gs.append(temp1.sdev)

        
    gA_fit_y1_gs = np.array(gA_fit_gs) + np.array(gA_fit_err_gs)
    gA_fit_y2_gs = np.array(gA_fit_gs) - np.array(gA_fit_err_gs)
    
    ax.fill_between(x_fill, gA_fit_y1_gs, gA_fit_y2_gs, color=red, alpha=0.3)

######################
###################### sum-sub / g.s. only of 2pt
    gA_fit_gs = []
    gA_fit_err_gs = []

    if sum_nstates == pt2_nstates:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2)['pt2']
        sum_A3_fitter = fitter.summation_same_can(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_A3']
        
    else:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2, sum_nstates)['pt2']
        sum_A3_fitter = fitter.summation(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_A3']
    
    for i in range(len(np.arange(plt_min, plt_max, plot_space)) - int(1 / plot_space)):
        temp1 = sum_A3_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i]
        gA_fit_gs.append(temp1.mean)
        gA_fit_err_gs.append(temp1.sdev)

        
    gA_fit_y1_gs = np.array(gA_fit_gs) + np.array(gA_fit_err_gs)
    gA_fit_y2_gs = np.array(gA_fit_gs) - np.array(gA_fit_err_gs)
    
    #ax.fill_between(x_fill, gA_fit_y1_gs, gA_fit_y2_gs, color=grey, alpha=0.3)

##################

    x_lim = [plt_min-0.5, plt_max+0.5]
    if plot_in_fm == False:
        ax.set_xlim(x_lim)
        ax.set_xlabel(t_label, **textp)
    elif plot_in_fm == True:
        ax.set_xlim([num*omega_imp_a09 for num in x_lim])
        ax.set_xlabel(fm_t_label, **textp)
    
    ax.set_xlim([0, 5])
    ax.set_ylim([1, 1.4])
    ax.set_ylabel(oaeff_label, **textp)
    ax.tick_params(axis='both', which='major', **labelp)
    ax.text(2.7, 1.35, r'$ntra$', fontsize=15)
    
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/soaeff{plot_type}_ntra.pdf", transparent=True)
    plt.show()


def plot_sum_no_sca(pt3_data_range, data_avg_dict_completed, fit_result, fitter, pt2_nstates, sum_nstates, plot_type, plot_in_fm=False):
    '''plot form factor with sum data, you can also plot fit on data'''
    gA = []
    gA_err = []
    
    gA_fit = []
    gA_fit_err = []
    gV_fit = []
    gV_fit_err = []
    
    plot_space = 0.05

    plt_min = 2
    plt_max = 60 
    fillx = np.arange(plt_min, plt_max, plot_space)[:int(-1 / plot_space)] # as tsep input

    if plot_in_fm == False:
        x_errorbar = np.arange(pt3_data_range[0], pt3_data_range[1]-1)
        x_fill = np.arange(plt_min, plt_max, plot_space)[:int(-1 / plot_space)]

    elif plot_in_fm == True:
        x_errorbar = np.arange(pt3_data_range[0], pt3_data_range[1]-1) * omega_imp_a09
        x_fill = np.arange(plt_min, plt_max, plot_space)[:int(-1 / plot_space)] * omega_imp_a09

    else:
        x_errorbar = np.arange(pt3_data_range[0], pt3_data_range[1]-1)
        x_fill = np.arange(plt_min, plt_max, plot_space)[:int(-1 / plot_space)]
        print("Input error: plot_in_fm")

    for i in range(pt3_data_range[0], pt3_data_range[1]-1):
        temp1 = data_avg_dict_completed['sum_A3_fit_'+str(i)]
        gA.append(temp1.mean)
        gA_err.append(temp1.sdev)

##################### fit on data
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    ax.errorbar(x_errorbar, np.array(gA), yerr=np.array(gA_err), marker='o', color=blue, **errorp)

    if sum_nstates == pt2_nstates:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), fit_result.p)['pt2']
        sum_A3_fitter = fitter.summation_same_can(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_A3']
        sum_V4_fitter = fitter.summation_same_can(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_V4']
        
    else:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), fit_result.p, sum_nstates)['pt2']
        sum_A3_fitter = fitter.summation(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_A3']
        sum_V4_fitter = fitter.summation(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_V4']
    
    for i in range(len(np.arange(plt_min, plt_max, plot_space)) - int(1 / plot_space)):
        temp1 = sum_A3_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i]
        gA_fit.append(temp1.mean)
        gA_fit_err.append(temp1.sdev)

        temp2 = sum_V4_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_V4_fitter[i] / pt2_fitter[i]
        gV_fit.append(temp2.mean)
        gV_fit_err.append(temp2.sdev)
        
    gA_fit_y1 = np.array(gA_fit) + np.array(gA_fit_err)
    gA_fit_y2 = np.array(gA_fit) - np.array(gA_fit_err)
    
    
    ax.fill_between(x_fill, gA_fit_y1, gA_fit_y2, color=blue, alpha=0.3, label='fit')

######################################## 
    p_nsca = gv.BufferDict(fit_result.p) # no scattering except for g.s.
    for key in p_nsca:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] == key.split("_")[1][1]) and key != 'A3_00':
            p_nsca[key] = 0

    p_ntra = gv.BufferDict(fit_result.p) # no transition and no g.s.
    for key in p_ntra:
        if key.split("_")[0] in ["A3"] and (key.split("_")[1][0] != key.split("_")[1][1]):
            p_ntra[key] = 0

    p_ntra['A3_00'] = 0

    p_gs = gv.BufferDict(fit_result.p) # g.s. only
    for key in p_gs:
        if key.split("_")[0] in ["A3"] and key != 'A3_00':
            p_gs[key] = 0

    p_gs_pt2 = gv.BufferDict(fit_result.p) # g.s. only for 2pt
    for key in p_gs_pt2:
        if 'z'in key and key != 'z0' and key != 'log(z2)':
            p_gs_pt2[key] = 0

    del p_gs_pt2['log(z2)']
    p_gs_pt2['z2'] = 0

##################### original gA from 3pt with tau = tsep/2
    y3 = np.array([1.137246887069871, 1.0397951114352422, 1.0581101894217029, 1.0944768415815544, 1.129653822067308, 1.1607318494702763, 1.1883613569687608])
    dy3 = np.array([0.0006817408659827335, 0.0008354900714438219, 0.0011941822664530966, 0.0019041344115398975, 0.003132857777381124, 0.005059734516954951, 0.008446467991949032])
    x3 = np.arange(2, 15, 2)

    if plot_in_fm == True:
        x3 = x3 * omega_imp_a09

    ax.errorbar(x3, y3, yerr=dy3, ls='none', marker='o', capsize=3, color=orange) 

    pt2_fitter = fitter.pt2_fit_function(fillx, fit_result.p)['pt2']
    pt3_gA_fitter = fitter.pt3_fit_function(fillx, fillx, fillx/2, fillx/2, fit_result.p)['pt3_A3'] # here fillx/2 means tau = tsep/2

    demo_ori = pt3_gA_fitter/pt2_fitter

    ax.fill_between(x_fill, np.array([y.mean-y.sdev for y in demo_ori]), np.array([y.mean+y.sdev for y in demo_ori]), alpha=0.2, color=orange)


##################### data - no transition and no g.s. gA from 3pt with tau = tsep/2
    y3_nsca = np.array([0.6727799364966929, 0.7774506761828093, 0.9104118135085515, 1.0112160578126128, 1.0830704000959943, 1.1351676308811318, 1.1748843233329838])

    dy3_nsca = np.array([0.12391746097652938, 0.06866160207895389, 0.03807136583051402, 0.022605124172479476, 0.0138835672428265, 0.008723160838174095, 0.00821622163788732])

    x3_nsca = np.arange(2, 15, 2)

    if plot_in_fm == True:
        x3_nsca = x3_nsca * omega_imp_a09

    ax.errorbar(x3_nsca, y3_nsca, yerr=dy3_nsca, ls='none', marker='x', capsize=3, color=yellow) 

    pt2_fitter = fitter.pt2_fit_function(fillx, fit_result.p)['pt2']
    pt3_gA_fitter = fitter.pt3_fit_function(fillx, fillx, fillx/2, fillx/2, p_nsca)['pt3_A3']  # here fillx/2 means tau = tsep/2

    demo_nsca = pt3_gA_fitter/pt2_fitter

    ax.fill_between(x_fill, np.array([y.mean-y.sdev for y in demo_nsca]), np.array([y.mean+y.sdev for y in demo_nsca]), alpha=0.2, color=yellow)


###################### no scattering except for g.s. sum-sub
    gA_fit_nsca = []
    gA_fit_err_nsca = []

    if sum_nstates == pt2_nstates:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2)['pt2']
        sum_A3_fitter = fitter.summation_same_can(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), p_nsca)['sum_A3']
        
    else:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2, sum_nstates)['pt2']
        sum_A3_fitter = fitter.summation(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), p_nsca)['sum_A3']
    
    for i in range(len(np.arange(plt_min, plt_max, plot_space)) - int(1 / plot_space)):
        temp1 = sum_A3_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i]
        gA_fit_nsca.append(temp1.mean)
        gA_fit_err_nsca.append(temp1.sdev)

    ###
    print(p_nsca)
    print(p_gs_pt2)

        
    gA_fit_y1_nsca = np.array(gA_fit_nsca) + np.array(gA_fit_err_nsca)
    gA_fit_y2_nsca = np.array(gA_fit_nsca) - np.array(gA_fit_err_nsca)
    
    ax.fill_between(x_fill, gA_fit_y1_nsca, gA_fit_y2_nsca, color=green, alpha=0.3)


######################## data - no transition and no g.s. sum-sub
    gA_fit_nsca = []

    if sum_nstates == pt2_nstates:
        pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1]), p_gs_pt2)['pt2']
        sum_A3_fitter = fitter.summation_same_can(np.arange(pt3_data_range[0], pt3_data_range[1]), np.arange(pt3_data_range[0], pt3_data_range[1]), p_ntra)['sum_A3']
        
    else:
        pt2_fitter = fitter.pt2_fit_function(np.arange(pt3_data_range[0], pt3_data_range[1]), p_gs_pt2, sum_nstates)['pt2']
        sum_A3_fitter = fitter.summation(np.arange(pt3_data_range[0], pt3_data_range[1]), np.arange(pt3_data_range[0], pt3_data_range[1]), p_ntra)['sum_A3']
    
    for i in range(pt3_data_range[0], pt3_data_range[1]-1):
        index = i - pt3_data_range[0]
        temp1 = data_avg_dict_completed['sum_A3_fit_'+str(i)] - (sum_A3_fitter[index+1] / pt2_fitter[index+1] - sum_A3_fitter[index] / pt2_fitter[index])
        gA_fit_nsca.append(temp1) # nsca = data - ntra

    #ax.errorbar(x_errorbar, np.array([val.mean for val in gA_fit_nsca]), yerr=np.array([val.sdev for val in gA_fit_nsca]), ls='none', marker='x', capsize=3, color=green)

######################
###################### g.s. only
    gA_fit_gs = []
    gA_fit_err_gs = []

    if sum_nstates == pt2_nstates:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2)['pt2']
        sum_A3_fitter = fitter.summation_same_can(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), p_gs)['sum_A3']
        
    else:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2, sum_nstates)['pt2']
        sum_A3_fitter = fitter.summation(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), p_gs)['sum_A3']
    
    for i in range(len(np.arange(plt_min, plt_max, plot_space)) - int(1 / plot_space)):
        temp1 = sum_A3_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i]
        gA_fit_gs.append(temp1.mean)
        gA_fit_err_gs.append(temp1.sdev)

        
    gA_fit_y1_gs = np.array(gA_fit_gs) + np.array(gA_fit_err_gs)
    gA_fit_y2_gs = np.array(gA_fit_gs) - np.array(gA_fit_err_gs)
    
    ax.fill_between(x_fill, gA_fit_y1_gs, gA_fit_y2_gs, color=red, alpha=0.3)

######################
###################### g.s. only for 2pt
    gA_fit_gs = []
    gA_fit_err_gs = []

    if sum_nstates == pt2_nstates:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2)['pt2']
        sum_A3_fitter = fitter.summation_same_can(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_A3']
        
    else:
        pt2_fitter = fitter.pt2_fit_function(np.arange(plt_min, plt_max, plot_space), p_gs_pt2, sum_nstates)['pt2']
        sum_A3_fitter = fitter.summation(np.arange(plt_min, plt_max, plot_space), np.arange(plt_min, plt_max, plot_space), fit_result.p)['sum_A3']
    
    for i in range(len(np.arange(plt_min, plt_max, plot_space)) - int(1 / plot_space)):
        temp1 = sum_A3_fitter[i+int(1 / plot_space)] / pt2_fitter[i+int(1 / plot_space)] - sum_A3_fitter[i] / pt2_fitter[i]
        gA_fit_gs.append(temp1.mean)
        gA_fit_err_gs.append(temp1.sdev)

        
    gA_fit_y1_gs = np.array(gA_fit_gs) + np.array(gA_fit_err_gs)
    gA_fit_y2_gs = np.array(gA_fit_gs) - np.array(gA_fit_err_gs)
    
    ax.fill_between(x_fill, gA_fit_y1_gs, gA_fit_y2_gs, color=grey, alpha=0.3)

##################
    
    x_lim = [plt_min-0.5, plt_max+0.5]
    if plot_in_fm == False:
        ax.set_xlim(x_lim)
        ax.set_xlabel(t_label, **textp)
    elif plot_in_fm == True:
        ax.set_xlim([num*omega_imp_a09 for num in x_lim])
        ax.set_xlabel(fm_t_label, **textp)

    ax.set_ylim([1, 1.4])
    ax.set_ylabel(oaeff_label, **textp)
    ax.tick_params(axis='both', which='major', **labelp)
    ax.text(2.7, 1.35, r'$nsca$', fontsize=15)
    
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/soaeff{plot_type}_nsca.pdf", transparent=True)
    plt.show()

# %%
def tmin_plot(n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, gA_ylim, gV_ylim, E0_ylim, fit_name, xlabel):
    value={}
    value['Q']=[]
    value['logGBF']=[]
    value['E0']=[]
    value['E0_err']=[]
    value['gA']=[]
    value['gA_err']=[]
    value['gV']=[]
    value['gV_err']=[]
    x = []

    for n_ in range(n_range[1]): # n_ represents index of n in n_range
        value['Q'].append([])
        value['logGBF'].append([])
        value['E0'].append([])
        value['E0_err'].append([])
        value['gA'].append([])
        value['gA_err'].append([])
        value['gV'].append([])
        value['gV_err'].append([])
        x.append([])

    for n in range(n_range[0], n_range[1]): # n represents nstates value
        for situation in situation_list:         
            nstate_dict = {}
            nstate_dict['2pt'] = situation.pt2_nstates
            nstate_dict['3pt'] = situation.pt3_nstates
            nstate_dict['sum'] = situation.sum_nstates
            
            if nstate_dict[nstate_name] == n:
                tmin_dict = {}
                tmin_dict['2pt'] = situation.pt2_tmin
                tmin_dict['3pt_gA'] = situation.pt3_A3_tsep_min
                tmin_dict['3pt_gV'] = situation.pt3_V4_tsep_min
                tmin_dict['sum_gA'] = situation.sum_A3_tsep_min
                tmin_dict['sum_gV'] = situation.sum_V4_tsep_min
                
                x[n].append(tmin_dict[tmin_name])  # here is the varying parameter

        x[n].sort() # fix n, search for all situations to complete x[n], then sort x[n]

        for i in range(len(x[n])):
            for situation in situation_list:
                nstate_dict = {}
                nstate_dict['2pt'] = situation.pt2_nstates
                nstate_dict['3pt'] = situation.pt3_nstates
                nstate_dict['sum'] = situation.sum_nstates

                tmin_dict = {}
                tmin_dict['2pt'] = situation.pt2_tmin
                tmin_dict['3pt_gA'] = situation.pt3_A3_tsep_min
                tmin_dict['3pt_gV'] = situation.pt3_V4_tsep_min
                tmin_dict['sum_gA'] = situation.sum_A3_tsep_min
                tmin_dict['sum_gV'] = situation.sum_V4_tsep_min
                
                if nstate_dict[nstate_name] == n and tmin_dict[tmin_name] == x[n][i]:
                    value['Q'][n].append(situation.Q_value)
                    value['logGBF'][n].append(situation.log_GBF)
                    value['E0'][n].append(situation.E0)
                    value['E0_err'][n].append(situation.E0_err)
                    value['gA'][n].append(situation.A300)
                    value['gA_err'][n].append(situation.A300_err)
                    value['gV'][n].append(situation.V400)
                    value['gV_err'][n].append(situation.V400_err)
                
            
    print(x)
    
    best_n_ = best_n - n_range[0] # n_ represents index, n represents value
    best_t_ = best_t - t_range[0] # t_ represents index, t represents value

    #####################################################################################
    #####################################################################################
    # gA - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      


    # ax1
    ax1.set_ylabel(oa00_label, **textp)
    ax1.set_ylim(gA_ylim)


    for n in range(n_range[0], n_range[1]):
        ax1.errorbar(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['gA'][n]), yerr=np.array(value['gA_err'][n]), marker='o', color=color_list[n-2], **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), (value['gA'][best_n][best_t_]+value['gA_err'][best_n][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['gA'][best_n][best_t - t_range[0]]-value['gA_err'][best_n][best_t - t_range[0]])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n - 2], alpha=0.2)

    # best fit
    ax1.errorbar(np.array([best_t + (best_n-4)*0.1]), np.array([value['gA'][best_n][best_t_]]), yerr=np.array([value['gA_err'][best_n][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)
    print("Best fit")
    print(np.array([value['gA'][best_n][best_t_]]), np.array([value['gA_err'][best_n][best_t_]]))
    
    
    # ax2
    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])

    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        ax2.scatter(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['Q'][n]), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([value['Q'][best_n][best_t_]]), marker='o', c=color_list[best_n-2])

    
    # ax3
    ax3.set_ylabel(w_label, **textp)
    ax3.set_ylim([0, 1.1])
    
    ax3.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.3 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')
    
    log_max = {}
    
    for t in range(t_range[0], t_range[1]):
        t_ = t - t_range[0]
        logGBF_list = []
        for n in range(n_range[0], n_range[1]):
            logGBF_list.append(value['logGBF'][n][t_])
            
        log_max['t='+str(t)] = max(logGBF_list)
        
        w_list = []
        for n in range(n_range[0], n_range[1]):
            w = np.exp(value['logGBF'][n][t_] - log_max['t='+str(t)])
            ax3.scatter(np.array([t + (n-4)*0.1]), np.array([w]), marker='o', c='', edgecolors=color_list[n-2])
        
    # best fit
    ax3.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([np.exp( value['logGBF'][best_n][best_t_] - log_max['t='+str(best_t)] )]), marker='o', c=color_list[best_n-2])

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([t_range[0] - 0.5, t_range[1] - 0.5])
    ax1.tick_params(axis='both', which='major', **labelp)

    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_gA-'+tmin_name+'_tmin.pdf', transparent=True)
    
    #####################################################################################
    #####################################################################################
    # gV - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)     


    # ax1
    ax1.set_ylabel(ov00_label, **textp)
    ax1.set_ylim(gV_ylim)


    for n in range(n_range[0], n_range[1]):
        ax1.errorbar(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['gV'][n]), yerr=np.array(value['gV_err'][n]), marker='o', color=color_list[n-2], **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), (value['gV'][best_n][best_t_]+value['gV_err'][best_n][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['gV'][best_n][best_t - t_range[0]]-value['gV_err'][best_n][best_t - t_range[0]])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n - 2], alpha=0.2)

    # best fit
    ax1.errorbar(np.array([best_t + (best_n-4)*0.1]), np.array([value['gV'][best_n][best_t_]]), yerr=np.array([value['gV_err'][best_n][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)
    print("Best fit")
    print(np.array([value['gV'][best_n][best_t_]]), np.array([value['gV_err'][best_n][best_t_]]))
    
    
    # ax2
    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])

    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        ax2.scatter(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['Q'][n]), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([value['Q'][best_n][best_t_]]), marker='o', c=color_list[best_n-2])

    
    # ax3
    ax3.set_ylabel(w_label, **textp)
    ax3.set_ylim([0, 1.1])
    
    ax3.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.3 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')
    
    log_max = {}
    
    for t in range(t_range[0], t_range[1]):
        t_ = t - t_range[0]
        logGBF_list = []
        for n in range(n_range[0], n_range[1]):
            logGBF_list.append(value['logGBF'][n][t_])
            
        log_max['t='+str(t)] = max(logGBF_list)
        
        w_list = []
        for n in range(n_range[0], n_range[1]):
            w = np.exp(value['logGBF'][n][t_] - log_max['t='+str(t)])
            ax3.scatter(np.array([t + (n-4)*0.1]), np.array([w]), marker='o', c='', edgecolors=color_list[n-2])
        
    # best fit
    ax3.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([np.exp( value['logGBF'][best_n][best_t_] - log_max['t='+str(best_t)] )]), marker='o', c=color_list[best_n-2])

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([t_range[0] - 0.5, t_range[1] - 0.5])
    ax1.tick_params(axis='both', which='major', **labelp)
    ax2.tick_params(axis='both', which='major', **labelp)
    ax2.tick_params(axis='both', which='major', **labelp)

    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_gV-'+tmin_name+'_tmin.pdf', transparent=True)

    #####################################################################################
    #####################################################################################
    # E0 - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      

    # ax1
    ax1.set_ylabel(e0_label, **textp)
    ax1.set_ylim(E0_ylim)


    for n in range(n_range[0], n_range[1]):
        ax1.errorbar(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['E0'][n]), yerr=np.array(value['E0_err'][n]), marker='o', color=color_list[n-2], **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), (value['E0'][best_n][best_t_]+value['E0_err'][best_n][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['E0'][best_n][best_t - t_range[0]]-value['E0_err'][best_n][best_t - t_range[0]])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n - 2], alpha=0.2)

    # best fit
    ax1.errorbar(np.array([best_t + (best_n-4)*0.1]), np.array([value['E0'][best_n][best_t_]]), yerr=np.array([value['E0_err'][best_n][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)

    # ax2
    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])


    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        ax2.scatter(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['Q'][n]), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([value['Q'][best_n][best_t_]]), marker='o', c=color_list[best_n-2])

    
    # ax3
    ax3.set_ylabel(w_label, **textp)
    ax3.set_ylim([0, 1.1])
    
    ax3.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.3 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')
    
    log_max = {}
    
    for t in range(t_range[0], t_range[1]):
        t_ = t - t_range[0]
        logGBF_list = []
        for n in range(n_range[0], n_range[1]):
            logGBF_list.append(value['logGBF'][n][t_])
            
        log_max['t='+str(t)] = max(logGBF_list)
        
        w_list = []
        for n in range(n_range[0], n_range[1]):
            w = np.exp(value['logGBF'][n][t_] - log_max['t='+str(t)])
            ax3.scatter(np.array([t + (n-4)*0.1]), np.array([w]), marker='o', c='', edgecolors=color_list[n-2])
        
    # best fit
    ax3.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([np.exp( value['logGBF'][best_n][best_t_] - log_max['t='+str(best_t)] )]), marker='o', c=color_list[best_n-2])
    

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([t_range[0] - 0.5, t_range[1] - 0.5])
    ax1.tick_params(axis='both', which='major', **labelp)
    ax2.tick_params(axis='both', which='major', **labelp)
    ax3.tick_params(axis='both', which='major', **labelp)

    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_E0-'+tmin_name+'_tmin.pdf', transparent=True)

# %%
def tmin_div_plot(n_range, t_range, best_n, best_t, tmin_name, situation_list, gA_ylim, gV_ylim, E0_ylim, fit_name, xlabel):
    value={}
    value['Q']=[]
    value['E0']=[]
    value['E0_err']=[]
    value['gA']=[]
    value['gA_err']=[]
    value['gV']=[]
    value['gV_err']=[]
    x = []

    for n_ in range(n_range[1]):
        value['Q'].append([])
        value['E0'].append([])
        value['E0_err'].append([])
        value['gA'].append([])
        value['gA_err'].append([])
        value['gV'].append([])
        value['gV_err'].append([])
        x.append([])


    for n in range(n_range[0], n_range[1]):
        n_ = n - n_range[0]
        for situation in situation_list[n_]:         
            value['Q'][n].append(situation.Q_value)
            value['E0'][n].append(situation.E0)
            value['E0_err'][n].append(situation.E0_err)
            value['gA'][n].append(situation.A300)
            value['gA_err'][n].append(situation.A300_err)
            value['gV'][n].append(situation.V400)
            value['gV_err'][n].append(situation.V400_err)

            tmin_dict = {}
            tmin_dict['2pt'] = situation.pt2_tmin
            tmin_dict['3pt_gA'] = situation.pt3_A3_tsep_min
            tmin_dict['3pt_gV'] = situation.pt3_V4_tsep_min
            tmin_dict['sum_gA'] = situation.sum_A3_tsep_min
            tmin_dict['sum_gV'] = situation.sum_V4_tsep_min

            x[n].append(tmin_dict[tmin_name])  # here is the varying parameter

    print(x)
    
    best_n_ = best_n - n_range[0]
    best_t_ = []
    for n_ in range(len(best_t)):
        best_t_.append(best_t[n_] - t_range[n_][0])
    
    plot_tmin = t_range[n_range[1] - n_range[0] - 1][0]
    plot_tmax = t_range[0][1]

    #####################################################################################
    #####################################################################################
    # gA - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmin_div)      


    # ax1
    ax1.set_ylabel(oa00_label, **textp)
    ax1.set_ylim(gA_ylim)


    for n in range(n_range[0], n_range[1]):
        n_ = n - n_range[0]
        ax1.errorbar(np.arange(t_range[n_][0], t_range[n_][1]) + (n-2) * 0.2, np.array(value['gA'][n]), yerr=np.array(value['gA_err'][n]), marker='o', color=color_list[n-2], **errorp)

    # best fit
    ax1.fill_between(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['gA'][best_n][best_t_[best_n_]]+value['gA_err'][best_n][best_t_[best_n_]])*np.ones([plot_tmax - plot_tmin + 1]), (value['gA'][best_n][best_t_[best_n_]]-value['gA_err'][best_n][best_t_[best_n_]])*np.ones([plot_tmax - plot_tmin + 1]), color=color_list[best_n - 2], alpha=0.2)
    for n_ in range(n_range[1] - n_range[0]):
        n = n_ + n_range[0]
        if n != best_n:
            ax1.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['gA'][n][best_t_[n_]]+value['gA_err'][n][best_t_[n_]])*np.ones([plot_tmax - plot_tmin + 1]), color = color_list[n-2], linestyle='--')
            ax1.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['gA'][n][best_t_[n_]]-value['gA_err'][n][best_t_[n_]])*np.ones([plot_tmax - plot_tmin + 1]), color = color_list[n-2], linestyle='--')
        
        ax1.errorbar(np.array([best_t[n_] + (n - 2)*0.2]), np.array([value['gA'][n][best_t_[n_]]]), yerr=np.array([value['gA_err'][n][best_t_[n_]]]), marker='o', mfc=color_list[n-2], color=color_list[n-2], **errorb)
    
    ax1.plot()
    
    print("Best fit")
    print(np.array([value['gA'][best_n][best_t_[best_n_]]]), np.array([value['gA_err'][best_n][best_t_[best_n_]]]))
    
    
    # ax2
    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])

    ax2.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), 0.1 * np.ones([plot_tmax - plot_tmin + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        n_ = n - n_range[0]
        ax2.scatter(np.arange(t_range[n_][0], t_range[n_][1]) + (n-2) * 0.2, np.array(value['Q'][n]), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    for n_ in range(n_range[1] - n_range[0]):
        n = n_ + n_range[0]    
        ax2.scatter(np.array([best_t[n_] + (n-2)*0.2]), np.array([value['Q'][n][best_t_[n_]]]), marker='o', c=color_list[n-2])
    

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([plot_tmin - 0.5, plot_tmax - 0.5])
    ax1.tick_params(axis='both', which='major', **labelp)
    ax2.tick_params(axis='both', which='major', **labelp)

    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_gA-'+tmin_name+'_tmin.pdf', transparent=True)
    
    #####################################################################################
    #####################################################################################
    # gV - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmin_div)      


    # ax1
    ax1.set_ylabel(ov00_label, **textp)
    ax1.set_ylim(gV_ylim)


    for n in range(n_range[0], n_range[1]):
        n_ = n - n_range[0]
        ax1.errorbar(np.arange(t_range[n_][0], t_range[n_][1]) + (n-2) * 0.2, np.array(value['gV'][n]), yerr=np.array(value['gV_err'][n]), marker='o', color=color_list[n-2], **errorp)

    # best fit
    ax1.fill_between(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['gV'][best_n][best_t_[best_n_]]+value['gV_err'][best_n][best_t_[best_n_]])*np.ones([plot_tmax - plot_tmin + 1]), (value['gV'][best_n][best_t_[best_n_]]-value['gV_err'][best_n][best_t_[best_n_]])*np.ones([plot_tmax - plot_tmin + 1]), color=color_list[best_n - 2], alpha=0.2)
    for n_ in range(n_range[1] - n_range[0]):
        n = n_ + n_range[0]
        if n != best_n:
            ax1.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['gV'][n][best_t_[n_]]+value['gV_err'][n][best_t_[n_]])*np.ones([plot_tmax - plot_tmin + 1]), color = color_list[n-2], linestyle='--')
            ax1.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['gV'][n][best_t_[n_]]-value['gV_err'][n][best_t_[n_]])*np.ones([plot_tmax - plot_tmin + 1]), color = color_list[n-2], linestyle='--')

        ax1.errorbar(np.array([best_t[n_] + (n - 2)*0.2]), np.array([value['gV'][n][best_t_[n_]]]), yerr=np.array([value['gV_err'][n][best_t_[n_]]]), marker='o', mfc=color_list[n-2], color=color_list[n-2], **errorb)
        
    print("Best fit")
    print(np.array([value['gV'][best_n][best_t_[best_n_]]]), np.array([value['gV_err'][best_n][best_t_[best_n_]]]))
    
    
    # ax2
    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])

    ax2.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), 0.1 * np.ones([plot_tmax - plot_tmin + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        n_ = n - n_range[0]
        ax2.scatter(np.arange(t_range[n_][0], t_range[n_][1]) + (n-2) * 0.2, np.array(value['Q'][n]), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    for n_ in range(n_range[1] - n_range[0]):
        n = n_ + n_range[0]    
        ax2.scatter(np.array([best_t[n_] + (n-2)*0.2]), np.array([value['Q'][n][best_t_[n_]]]), marker='o', c=color_list[n-2])


    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([plot_tmin - 0.5, plot_tmax - 0.5])
    ax1.tick_params(axis='both', which='major', **labelp)
    ax2.tick_params(axis='both', which='major', **labelp)

    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_gV-'+tmin_name+'_tmin.pdf', transparent=True)

    #####################################################################################
    #####################################################################################
    # E0 - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmin_div)    

    # ax1
    ax1.set_ylabel(e0_label, **textp)
    ax1.set_ylim(E0_ylim)


    for n in range(n_range[0], n_range[1]):
        n_ = n - n_range[0]
        ax1.errorbar(np.arange(t_range[n_][0], t_range[n_][1]) + (n-2) * 0.2, np.array(value['E0'][n]), yerr=np.array(value['E0_err'][n]), marker='o', color=color_list[n-2], **errorp)

    # best fit
    ax1.fill_between(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['E0'][best_n][best_t_[best_n_]]+value['E0_err'][best_n][best_t_[best_n_]])*np.ones([plot_tmax - plot_tmin + 1]), (value['E0'][best_n][best_t_[best_n_]]-value['E0_err'][best_n][best_t_[best_n_]])*np.ones([plot_tmax - plot_tmin + 1]), color=color_list[best_n - 2], alpha=0.2)
    for n_ in range(n_range[1] - n_range[0]):
        n = n_ + n_range[0]
        if n != best_n:
            ax1.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['E0'][n][best_t_[n_]]+value['E0_err'][n][best_t_[n_]])*np.ones([plot_tmax - plot_tmin + 1]), color = color_list[n-2], linestyle='--')
            ax1.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), (value['E0'][n][best_t_[n_]]-value['E0_err'][n][best_t_[n_]])*np.ones([plot_tmax - plot_tmin + 1]), color = color_list[n-2], linestyle='--')

        ax1.errorbar(np.array([best_t[n_] + (n - 2)*0.2]), np.array([value['E0'][n][best_t_[n_]]]), yerr=np.array([value['E0_err'][n][best_t_[n_]]]), marker='o', mfc=color_list[n-2], color=color_list[n-2], **errorb)
    
    
    # ax2
    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])

    ax2.plot(np.arange(plot_tmin - 0.5, plot_tmax + 0.5, 1), 0.1 * np.ones([plot_tmax - plot_tmin + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        n_ = n - n_range[0]
        ax2.scatter(np.arange(t_range[n_][0], t_range[n_][1]) + (n-2) * 0.2, np.array(value['Q'][n]), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    for n_ in range(n_range[1] - n_range[0]):
        n = n_ + n_range[0]    
        ax2.scatter(np.array([best_t[n_] + (n-2)*0.2]), np.array([value['Q'][n][best_t_[n_]]]), marker='o', c=color_list[n-2])
    

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([plot_tmin - 0.5, plot_tmax - 0.5])
    ax1.tick_params(axis='both', which='major', **labelp)
    ax2.tick_params(axis='both', which='major', **labelp)

    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_E0-'+tmin_name+'_tmin.pdf', transparent=True)

# %%
def tmin_late_plot(n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, gA_ylim, gV_ylim, E0_ylim, z0_ylim, fit_name, xlabel):
    value={}
    value['Q']=[]
    value['logGBF']=[]
    value['E0']=[]
    value['E0_err']=[]
    value['z0']=[]
    value['z0_err']=[]
    value['gA']=[]
    value['gA_err']=[]
    value['gV']=[]
    value['gV_err']=[]
    x = []

    for n_ in range(n_range[1]):
        value['Q'].append([])
        value['logGBF'].append([])
        value['E0'].append([])
        value['E0_err'].append([])
        value['z0'].append([])
        value['z0_err'].append([])
        value['gA'].append([])
        value['gA_err'].append([])
        value['gV'].append([])
        value['gV_err'].append([])
        x.append([])


    for n in range(n_range[0], n_range[1]): # n represents nstates value
        for situation in situation_list:         
            nstate_dict = {}
            nstate_dict['2pt'] = situation.pt2_nstates
            nstate_dict['3pt'] = situation.pt3_nstates
            nstate_dict['sum'] = situation.sum_nstates
            
            if nstate_dict[nstate_name] == n:
                tmin_dict = {}
                tmin_dict['2pt'] = situation.pt2_tmin
                tmin_dict['3pt_gA'] = situation.pt3_A3_tsep_min
                tmin_dict['3pt_gV'] = situation.pt3_V4_tsep_min
                tmin_dict['sum_gA'] = situation.sum_A3_tsep_min
                tmin_dict['sum_gV'] = situation.sum_V4_tsep_min
                
                x[n].append(tmin_dict[tmin_name])  # here is the varying parameter

        x[n].sort() # fix n, search for all situations to complete x[n], then sort x[n]

        for i in range(len(x[n])):
            for situation in situation_list:
                nstate_dict = {}
                nstate_dict['2pt'] = situation.pt2_nstates
                nstate_dict['3pt'] = situation.pt3_nstates
                nstate_dict['sum'] = situation.sum_nstates

                tmin_dict = {}
                tmin_dict['2pt'] = situation.pt2_tmin
                tmin_dict['3pt_gA'] = situation.pt3_A3_tsep_min
                tmin_dict['3pt_gV'] = situation.pt3_V4_tsep_min
                tmin_dict['sum_gA'] = situation.sum_A3_tsep_min
                tmin_dict['sum_gV'] = situation.sum_V4_tsep_min
                
                if nstate_dict[nstate_name] == n and tmin_dict[tmin_name] == x[n][i]:
                    value['Q'][n].append(situation.Q_value)
                    value['logGBF'][n].append(situation.log_GBF)
                    value['E0'][n].append(situation.E0)
                    value['E0_err'][n].append(situation.E0_err)
                    value['gA'][n].append(situation.A300)
                    value['gA_err'][n].append(situation.A300_err)
                    value['gV'][n].append(situation.V400)
                    value['gV_err'][n].append(situation.V400_err)

    print(x)
    
    best_n_ = best_n - n_range[0]
    best_t_ = best_t - t_range[0]

    #####################################################################################
    #####################################################################################
    # gA - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      


    # ax1
    ax1.set_ylabel(oa00_label, **textp)
    ax1.set_ylim(gA_ylim)


    for n in range(n_range[0], n_range[1]):
        ax1.errorbar(np.array(x[n]) + (n-4) * 0.1, np.array(value['gA'][n]), yerr=np.array(value['gA_err'][n]), marker='o', color=color_list[n-2], **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), (value['gA'][best_n][best_t_]+value['gA_err'][best_n][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['gA'][best_n][best_t - t_range[0]]-value['gA_err'][best_n][best_t - t_range[0]])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n - 2], alpha=0.2)

    # best fit
    ax1.errorbar(np.array([best_t + (best_n-4)*0.1]), np.array([value['gA'][best_n][best_t_]]), yerr=np.array([value['gA_err'][best_n][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)
    print("Best fit")
    print(np.array([value['gA'][best_n][best_t_]]), np.array([value['gA_err'][best_n][best_t_]]))
    
    
    # ax2
    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])

    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        ax2.scatter(np.array(x[n]) + (n-4) * 0.1, np.array(value['Q'][n]), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([value['Q'][best_n][best_t_]]), marker='o', c=color_list[best_n-2])

    
    # ax3
    ax3.set_ylabel(w_label, **textp)
    ax3.set_ylim([0, 1.1])
    
    ax3.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.3 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')
    
    log_max = {}
    
    for t in range(t_range[0], t_range[1]):
        t_ = t - t_range[0]
        logGBF_list = []
        for n in range(n_range[0], n_range[1]):
            logGBF_list.append(value['logGBF'][n][t_])
            
        log_max['t='+str(t)] = max(logGBF_list)
        
        w_list = []
        for n in range(n_range[0], n_range[1]):
            w = np.exp(value['logGBF'][n][t_] - log_max['t='+str(t)])
            ax3.scatter(np.array([t + (n-4)*0.1]), np.array([w]), marker='o', c='', edgecolors=color_list[n-2])
        
    # best fit
    ax3.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([np.exp( value['logGBF'][best_n][best_t_] - log_max['t='+str(best_t)] )]), marker='o', c=color_list[best_n-2])

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([t_range[0] - 0.5, t_range[1] - 0.5])
    ax1.tick_params(axis='both', which='major', **labelp)
    ax2.tick_params(axis='both', which='major', **labelp)
    ax3.tick_params(axis='both', which='major', **labelp)

    plt.tight_layout(pad=30, rect=aspect)

    #plt.savefig('./new_plots/'+fit_name+'_gA-'+tmin_name+'_tmin.pdf', transparent=True)
    
    #####################################################################################
    #####################################################################################
    # gV - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      


    # ax1
    ax1.set_ylabel(ov00_label, **textp)
    ax1.set_ylim(gV_ylim)


    for n in range(n_range[0], n_range[1]):
        ax1.errorbar(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['gV'][n]), yerr=np.array(value['gV_err'][n]), marker='o', color=color_list[n-2], **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), (value['gV'][best_n][best_t_]+value['gV_err'][best_n][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['gV'][best_n][best_t - t_range[0]]-value['gV_err'][best_n][best_t - t_range[0]])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n - 2], alpha=0.2)

    # best fit
    ax1.errorbar(np.array([best_t + (best_n-4)*0.1]), np.array([value['gV'][best_n][best_t_]]), yerr=np.array([value['gV_err'][best_n][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)
    print("Best fit")
    print(np.array([value['gV'][best_n][best_t_]]), np.array([value['gV_err'][best_n][best_t_]]))
    
    
    # ax2
    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])

    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        ax2.scatter(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['Q'][n]), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([value['Q'][best_n][best_t_]]), marker='o', c=color_list[best_n-2])

    
    # ax3
    ax3.set_ylabel(w_label, **textp)
    ax3.set_ylim([0, 1.1])
    
    ax3.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.3 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')
    
    log_max = {}
    
    for t in range(t_range[0], t_range[1]):
        t_ = t - t_range[0]
        logGBF_list = []
        for n in range(n_range[0], n_range[1]):
            logGBF_list.append(value['logGBF'][n][t_])
            
        log_max['t='+str(t)] = max(logGBF_list)
        
        w_list = []
        for n in range(n_range[0], n_range[1]):
            w = np.exp(value['logGBF'][n][t_] - log_max['t='+str(t)])
            ax3.scatter(np.array([t + (n-4)*0.1]), np.array([w]), marker='o', c='', edgecolors=color_list[n-2])
        
    # best fit
    ax3.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([np.exp( value['logGBF'][best_n][best_t_] - log_max['t='+str(best_t)] )]), marker='o', c=color_list[best_n-2])

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([t_range[0] - 0.5, t_range[1] - 0.5])
    ax1.tick_params(axis='both', which='major', **labelp)
    ax2.tick_params(axis='both', which='major', **labelp)
    ax3.tick_params(axis='both', which='major', **labelp)

    plt.tight_layout(pad=30, rect=aspect)

    #plt.savefig('./new_plots/'+fit_name+'_gV-'+tmin_name+'_tmin.pdf', transparent=True)

    #####################################################################################
    #####################################################################################
    # E0 - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      

    # ax1
    ax1.set_ylabel(e0_label, **textp)
    ax1.set_ylim(E0_ylim)


    for n in range(n_range[0], n_range[1]):
        ax1.errorbar(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['E0'][n]), yerr=np.array(value['E0_err'][n]), marker='o', color=color_list[n-2], **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), (value['E0'][best_n][best_t_]+value['E0_err'][best_n][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['E0'][best_n][best_t - t_range[0]]-value['E0_err'][best_n][best_t - t_range[0]])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n - 2], alpha=0.2)

    # best fit
    ax1.errorbar(np.array([best_t + (best_n-4)*0.1]), np.array([value['E0'][best_n][best_t_]]), yerr=np.array([value['E0_err'][best_n][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)

    # ax2
    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])


    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        ax2.scatter(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['Q'][n]), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([value['Q'][best_n][best_t_]]), marker='o', c=color_list[best_n-2])

    
    # ax3
    ax3.set_ylabel(w_label, **textp)
    ax3.set_ylim([0, 1.1])
    
    ax3.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.3 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')
    
    log_max = {}
    
    for t in range(t_range[0], t_range[1]):
        t_ = t - t_range[0]
        logGBF_list = []
        for n in range(n_range[0], n_range[1]):
            logGBF_list.append(value['logGBF'][n][t_])
            
        log_max['t='+str(t)] = max(logGBF_list)
        
        w_list = []
        for n in range(n_range[0], n_range[1]):
            w = np.exp(value['logGBF'][n][t_] - log_max['t='+str(t)])
            ax3.scatter(np.array([t + (n-4)*0.1]), np.array([w]), marker='o', c='', edgecolors=color_list[n-2])
        
    # best fit
    ax3.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([np.exp( value['logGBF'][best_n][best_t_] - log_max['t='+str(best_t)] )]), marker='o', c=color_list[best_n-2])
    

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([t_range[0] - 0.5, t_range[1] - 0.5])
    ax1.tick_params(axis='both', which='major', **labelp)
    ax2.tick_params(axis='both', which='major', **labelp)
    ax3.tick_params(axis='both', which='major', **labelp)

    plt.tight_layout(pad=30, rect=aspect)

    #plt.savefig('./new_plots/'+fit_name+'_E0-'+tmin_name+'_tmin.pdf', transparent=True)
    plt.savefig('./new_plots/'+fit_name+'_E0-'+tmin_name+'_tmin.pdf', transparent=True)
    
    #####################################################################################
    #####################################################################################
    # z0 - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      

    # ax1
    ax1.set_ylabel(z0_label, **textp)
    ax1.set_ylim(z0_ylim)


    for n in range(n_range[0], n_range[1]):
        ax1.errorbar(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['z0'][n]), yerr=np.array(value['z0_err'][n]), marker='o', color=color_list[n-2], **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), (value['z0'][best_n][best_t_]+value['z0_err'][best_n][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['z0'][best_n][best_t - t_range[0]]-value['z0_err'][best_n][best_t - t_range[0]])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n - 2], alpha=0.2)

    # best fit
    ax1.errorbar(np.array([best_t + (best_n-4)*0.1]), np.array([value['z0'][best_n][best_t_]]), yerr=np.array([value['z0_err'][best_n][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)

    # ax2
    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])


    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    for n in range(n_range[0], n_range[1]):
        ax2.scatter(np.arange(t_range[0], t_range[1]) + (n-4) * 0.1, np.array(value['Q'][n]), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([value['Q'][best_n][best_t_]]), marker='o', c=color_list[best_n-2])

    
    # ax3
    ax3.set_ylabel(w_label, **textp)
    ax3.set_ylim([0, 1.1])
    
    ax3.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.3 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')
    
    log_max = {}
    
    for t in range(t_range[0], t_range[1]):
        t_ = t - t_range[0]
        logGBF_list = []
        for n in range(n_range[0], n_range[1]):
            logGBF_list.append(value['logGBF'][n][t_])
            
        log_max['t='+str(t)] = max(logGBF_list)
        
        w_list = []
        for n in range(n_range[0], n_range[1]):
            w = np.exp(value['logGBF'][n][t_] - log_max['t='+str(t)])
            ax3.scatter(np.array([t + (n-4)*0.1]), np.array([w]), marker='o', c='', edgecolors=color_list[n-2])
        
    # best fit
    ax3.scatter(np.array([best_t + (best_n-4)*0.1]), np.array([np.exp( value['logGBF'][best_n][best_t_] - log_max['t='+str(best_t)] )]), marker='o', c=color_list[best_n-2])
    

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([t_range[0] - 0.5, t_range[1] - 0.5])
    ax1.tick_params(axis='both', which='major', **labelp)
    ax2.tick_params(axis='both', which='major', **labelp)
    ax3.tick_params(axis='both', which='major', **labelp)

    plt.tight_layout(pad=30, rect=aspect)

    #plt.savefig('./new_plots/'+fit_name+'_E0-'+tmin_name+'_tmin.pdf', transparent=True)

# %%
def tmax_plot(t_range, best_n, best_t, tmax_name, situation_list, gA_ylim, gV_ylim, E0_ylim, fit_name, xlabel):

    value={}
    value['Q']=[]
    value['E0']=[]
    value['E0_err']=[]
    value['gA']=[]
    value['gA_err']=[]
    value['gV']=[]
    value['gV_err']=[]

    x=[]

    for situation in situation_list:
        tmax_dict = {}
        tmax_dict['2pt'] = situation.pt2_tmax
        tmax_dict['3pt'] = situation.pt3_A3_tsep_max
        tmax_dict['sum'] = situation.sum_A3_tsep_max

        x.append(tmax_dict[tmax_name])  # here is the varying parameter

    x.sort()
    for i in range(len(x)):
        for situation in situation_list:
            tmax_dict = {}
            tmax_dict['2pt'] = situation.pt2_tmax
            tmax_dict['3pt'] = situation.pt3_A3_tsep_max
            tmax_dict['sum'] = situation.sum_A3_tsep_max
            if tmax_dict[tmax_name] == x[i]:
                value['Q'].append(situation.Q_value)
                value['E0'].append(situation.E0)
                value['E0_err'].append(situation.E0_err)
                value['gA'].append(situation.A300)
                value['gA_err'].append(situation.A300_err)
                value['gV'].append(situation.V400)
                value['gV_err'].append(situation.V400_err)

    print(x)
    
    best_t_ = best_t - t_range[0]

    #####################################################################################
    #####################################################################################

    # gA - tmax
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)      
    ax1.set_ylabel(oa00_label, **textp)
    ax1.set_ylim(gA_ylim)
    ax1.tick_params(axis='both', which='major', **labelp)

    ax1.errorbar(np.array(x)-1, np.array(value['gA']), yerr=np.array(value['gA_err']), marker='o', color=color_list[best_n-2], **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1)-1, (value['gA'][best_t_]+value['gA_err'][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['gA'][best_t_]-value['gA_err'][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n-2], alpha=0.2)
    # best fit
    ax1.errorbar(np.array([best_t])-1, np.array([value['gA'][best_t_]]), yerr=np.array([value['gA_err'][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)


    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1)-1, 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    ax2.scatter(np.array(x)-1, np.array(value['Q']), marker='o', c='', edgecolors=color_list[best_n-2])

    # best fit
    ax2.scatter(np.array([best_t])-1, np.array([value['Q'][best_t_]]), marker='o', c=color_list[best_n-2])


    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([t_range[0] - 1.5, t_range[1] - 1.5])
    ax2.tick_params(axis='both', which='major', **labelp)
    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_gA-'+tmax_name+'_tmax.pdf', transparent=True)
    
    #####################################################################################
    #####################################################################################

    # gV - tmax
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)     
    ax1.set_ylabel(ov00_label, **textp)
    ax1.set_ylim(gV_ylim)
    ax1.tick_params(axis='both', which='major', **labelp)

    ax1.errorbar(np.array(x)-1, np.array(value['gV']), yerr=np.array(value['gV_err']), marker='o', color=color_list[best_n-2], **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1)-1, (value['gV'][best_t_]+value['gV_err'][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['gV'][best_t_]-value['gV_err'][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n-2], alpha=0.2)
    # best fit
    ax1.errorbar(np.array([best_t])-1, np.array([value['gV'][best_t_]]), yerr=np.array([value['gV_err'][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)


    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1)-1, 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    ax2.scatter(np.array(x)-1, np.array(value['Q']), marker='o', c='', edgecolors=color_list[best_n-2])

    # best fit
    ax2.scatter(np.array([best_t])-1, np.array([value['Q'][best_t_]]), marker='o', c=color_list[best_n-2])


    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([t_range[0] - 1.5, t_range[1] - 1.5])
    ax2.tick_params(axis='both', which='major', **labelp)
    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_gV-'+tmax_name+'_tmax.pdf', transparent=True)

    #####################################################################################
    #####################################################################################

    # E0 - tmax
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)     
    ax1.set_ylabel(e0_label, **textp)
    ax1.set_ylim(E0_ylim)


    ax1.errorbar(np.array(x)-1, np.array(value['E0']), yerr=np.array(value['E0_err']), marker='o', color=color_list[best_n-2], **errorp)

    ax1.fill_between(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1)-1, (value['E0'][best_t_]+value['E0_err'][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), (value['E0'][best_t_]-value['E0_err'][best_t_])*np.ones([t_range[1] - t_range[0] + 1]), color=color_list[best_n-2], alpha=0.2)
    # best fit
    ax1.errorbar(np.array([best_t])-1, np.array([value['E0'][best_t_]]), yerr=np.array([value['E0_err'][best_t_]]), marker='o', mfc=color_list[best_n-2], color=color_list[best_n-2], **errorb)


    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1)-1, 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

    ax2.scatter(np.array(x)-1, np.array(value['Q']), marker='o', c='', edgecolors=color_list[best_n-2])

    # best fit
    ax2.scatter(np.array([best_t])-1, np.array([value['Q'][best_t_]]), marker='o', c=color_list[best_n-2])


    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel(xlabel, **textp)
    plt.xlim([t_range[0] - 1.5, t_range[1] - 1.5])
    ax1.tick_params(axis='both', which='major', **labelp)
    ax2.tick_params(axis='both', which='major', **labelp)
    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_E0-'+tmax_name+'_tmax.pdf', transparent=True)

# %%
def tau_cut_plot(E0, E0_err, A3, A3_err, V4, V4_err, Q, n, fit_name, gA_ylim, gV_ylim, E0_ylim):
    #####################gA
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)      
    ax1.set_ylabel(oa00_label, **textp)
    ax1.set_ylim(gA_ylim)

    ax1.errorbar(np.arange(0, 5), np.array(A3), yerr=np.array(A3_err), marker='o', color=color_list[n-2], **errorp)

    ax1.fill_between(np.arange(0 - 0.5, 5 + 0.5, 1), (A3[2] + A3_err[2])*np.ones([6]), (A3[2] - A3_err[2])*np.ones([6]), color=color_list[n-2], alpha=0.2)
    # best fit
    ax1.errorbar(np.array([1]), np.array([A3[1]]), yerr=np.array([A3_err[1]]), marker='o', mfc=color_list[n-2], color=color_list[n-2], **errorb)


    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange(0 - 0.5, 5 + 0.5, 1), 0.1 * np.ones([6]), 'r--')

    ax2.scatter(np.arange(0, 5), np.array(Q), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([1]), np.array([Q[1]]), marker='o', c=color_list[n-2])


    plt.subplots_adjust(wspace=0, hspace=0)
    #plt.xlabel('$tau\ cut$', **textp)
    plt.xticks([0, 1, 2, 3, 4], [r'$-1$', r'$tau\ cut$', r'$+1$', r'$+2$', r'$+3$'])
    plt.xlim([-0.5, 4.5])
    ax2.tick_params(axis='both', which='major', **labelp)
    plt.tight_layout(pad=30, rect=aspect)
    plt.savefig("./new_plots/"+fit_name+"_gA_tins.pdf", transparent=True)

    #####################gV
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)      
    ax1.set_ylabel(ov00_label, **textp)
    ax1.set_ylim(gV_ylim)


    ax1.errorbar(np.arange(0, 5), np.array(V4), yerr=np.array(V4_err), marker='o', color=color_list[n-2], **errorp)

    ax1.fill_between(np.arange(0 - 0.5, 5 + 0.5, 1), (V4[2] + V4_err[2])*np.ones([6]), (V4[2] - V4_err[2])*np.ones([6]), color=color_list[n-2], alpha=0.2)
    # best fit
    ax1.errorbar(np.array([1]), np.array([V4[1]]), yerr=np.array([V4_err[1]]), marker='o', mfc=color_list[n-2], color=color_list[n-2], **errorb)


    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange(0 - 0.5, 5 + 0.5, 1), 0.1 * np.ones([6]), 'r--')

    ax2.scatter(np.arange(0, 5), np.array(Q), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([1]), np.array([Q[1]]), marker='o', c=color_list[n-2])


    plt.subplots_adjust(wspace=0, hspace=0)
    #plt.xlabel('$tau\ cut$', **textp)
    plt.xticks([0, 1, 2, 3, 4], [r'$-1$', r'$tau\ cut$', r'$+1$', r'$+2$', r'$+3$'])
    plt.xlim([-0.5, 4.5])
    ax2.tick_params(axis='both', which='major', **labelp)
    plt.tight_layout(pad=30, rect=aspect)
    plt.savefig("./new_plots/"+fit_name+"_gV_tins.pdf", transparent=True)


    ######################E0
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)      
    ax1.set_ylabel(e0_label, **textp)
    ax1.set_ylim(E0_ylim)

    ax1.errorbar(np.arange(0, 5), np.array(E0), yerr=np.array(E0_err), marker='o', color=color_list[n-2], **errorp)

    ax1.fill_between(np.arange(0 - 0.5, 5 + 0.5, 1), (E0[2] + E0_err[2])*np.ones([6]), (E0[2] - E0_err[2])*np.ones([6]), color=color_list[n-2], alpha=0.2)
    # best fit
    ax1.errorbar(np.array([1]), np.array([E0[1]]), yerr=np.array([E0_err[1]]), marker='o', mfc=color_list[n-2], color=color_list[n-2], **errorb)


    ax2.set_ylabel(q_label, **textp)
    ax2.set_ylim([0, 1.1])
    ax2.plot(np.arange(0 - 0.5, 5 + 0.5, 1), 0.1 * np.ones([6]), 'r--')

    ax2.scatter(np.arange(0, 5), np.array(Q), marker='o', c='', edgecolors=color_list[n-2])

    # best fit
    ax2.scatter(np.array([1]), np.array([Q[1]]), marker='o', c=color_list[n-2])


    plt.subplots_adjust(wspace=0, hspace=0)
    #plt.xlabel('$tau\ cut$', **textp)
    plt.xticks([0, 1, 2, 3, 4], [r'$-1$', r'$tau\ cut$', r'$+1$', r'$+2$', r'$+3$'])
    plt.xlim([-0.5, 4.5])
    ax2.tick_params(axis='both', which='major', **labelp)
    plt.tight_layout(pad=30, rect=aspect)
    plt.savefig("./new_plots/"+fit_name+"_E0_tins.pdf", transparent=True)
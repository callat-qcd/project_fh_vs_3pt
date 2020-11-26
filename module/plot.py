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

color_list = [red, orange, yellow, green, blue, violet, grey, grape, peach, fuschia]
marker_list = ["o" for _ in range(10)] #['v', 's', 'o', 'D', '^', 'X', 'P']

figsize = (7, 4)
aspect=[0.15, 0.15, 0.8, 0.8]
gridspec = {'height_ratios': [3, 1], 'left': 0.05, 'right': 0.95, 'bottom': 0.05, 'top': 0.95}
gridspec_tmin = {'height_ratios': [3, 1, 1], 'left': 0.1, 'right': 0.95, 'bottom': 0.15, 'top': 0.95}
gridspec_tmin_div = {'height_ratios': [3, 1], 'left': 0.1, 'right': 0.95, 'bottom': 0.15, 'top': 0.95}
gridspec_tmax = {'height_ratios': [3, 1], 'left': 0.1, 'right': 0.95, 'bottom': 0.15, 'top': 0.95}
gridspec_prior_width = {'height_ratios': [3, 1], 'left': 0.1, 'right': 0.95, 'bottom': 0.15, 'top': 0.95}
textp = {"fontsize": 12}
labelp = {"labelsize": 12}
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
t_label = r"$t$"
tau_label = r"$\tau$"
meff_label = r"$m_{\textrm{eff}}$"
zeff_label = r"$z_{\textrm{eff}}$"
oaeff_label = r"$O^A_{\textrm{eff}}$"
oveff_label = r"$O^V_{\textrm{eff}}$"

plt.rcParams['figure.figsize'] = figsize

# %%
def plot_pt2(pt2_data_range, data_avg_dict, fit_result=None, fitter=None, plot_type=None):
    '''plot effective mass with 2pt data, you can also plot fit on data'''
    m0_eff = []
    m0_eff_err = []
    m0_eff_fit = []
    m0_eff_fit_err = []
    
    plot_space = 0.05 # smaller, more smoothy

    for i in range(pt2_data_range[0], pt2_data_range[1]-1):
        temp = gv.log(data_avg_dict['pt2_tsep_'+str(i)] / data_avg_dict['pt2_tsep_'+str(i+1)])
        m0_eff.append(temp.mean)
        m0_eff_err.append(temp.sdev)

    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    ax.errorbar(np.arange(pt2_data_range[0], pt2_data_range[1]-1), np.array(m0_eff), yerr=np.array(m0_eff_err), marker='o', color="k", **errorp)
    
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
        ax.fill_between(np.arange(pt2_data_range[0], pt2_data_range[1], plot_space)[:int(-1/plot_space)], np.array(m0_eff_fit_y1), np.array(m0_eff_fit_y2), color=blue, alpha=0.3, label='fit')
            
    
    ax.set_xlim([2, 25])
    ax.set_ylim([0.45, 0.62])
    ax.set_xlabel(t_label, **textp)
    ax.set_ylabel(meff_label, **textp)
    
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
    ax.errorbar(np.arange(pt2_data_range[0], pt2_data_range[1]-1), [i.mean for i in zeff], yerr=[i.sdev for i in zeff], marker='o', color="k", **errorp)
 
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
        
        ax.fill_between(x, np.array(z0_eff_fit_y1), np.array(z0_eff_fit_y2), color=blue, alpha=0.3, label='fit')

   
    ax.set_xlim([2, 25])
    ax.set_ylim([0, 3E-7])
    ax.set_xlabel(t_label, **textp)
    ax.set_ylabel(zeff_label, **textp)
    
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/zeff{plot_type}.pdf", transparent=True)
    plt.show()

    return [np.arange(pt2_data_range[0], pt2_data_range[1]-1), np.array(m0_eff), np.array(m0_eff_err)] # save this for data plot of half stat

def plot_pt3(pt3_data_range, data_avg_dict_completed, tau_cut, fit_result=None, fitter=None, plot_type=None):
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
        color = color_list[idx]
        ax.errorbar(np.arange(-(i-2)/2, (i)/2, 1), np.array(gA['tsep_' + str(i)]), yerr=np.array(gA_err['tsep_' + str(i)]), marker='o', color=color, **errorp)# tau cut = 1

        if fit_result != None and fitter != None:
            gA_fit_y1 = np.array(gA_fit['tsep_'+str(i)]) + np.array(gA_fit_err['tsep_'+str(i)])
            gA_fit_y2 = np.array(gA_fit['tsep_'+str(i)]) - np.array(gA_fit_err['tsep_'+str(i)])
            ax.fill_between(np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density), gA_fit_y1, gA_fit_y2, color=color, alpha=0.3) # tau_cut=1
    
    ax.set_xlim([-6.5, 6.5])
    ax.set_ylim([1, 1.3])
    ax.set_xlabel(t_label, **textp)
    ax.set_ylabel(oaeff_label, **textp)
    #
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/oaeff{plot_type}.pdf", transparent=True)
    plt.show()
    
    
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    for idx, i in enumerate(range(pt3_data_range[0], pt3_data_range[1])):
        color = color_list[idx]
        ax.errorbar(np.arange(-(i-2)/2, (i)/2, 1), np.array(gV['tsep_' + str(i)]), yerr=np.array(gV_err['tsep_' + str(i)]), marker='o', color=color, **errorp)# tau cut = 1

        if fit_result != None and fitter != None:
            gV_fit_y1 = np.array(gV_fit['tsep_'+str(i)]) + np.array(gV_fit_err['tsep_'+str(i)])
            gV_fit_y2 = np.array(gV_fit['tsep_'+str(i)]) - np.array(gV_fit_err['tsep_'+str(i)])
            ax.fill_between(np.linspace(tau_cut - i/2, i/2 - tau_cut, plot_density), gV_fit_y1, gV_fit_y2, color=color, alpha=0.3) # tau_cut=1
    
    ax.set_xlim([-6.5, 6.5])
    ax.set_ylim([1.0, 1.15])
    ax.set_xlabel(t_label, **textp)
    ax.set_ylabel(oveff_label, **textp)
    #
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/oveff{plot_type}.pdf", transparent=True)
    plt.show()
        
def plot_sum(pt3_data_range, data_avg_dict_completed, fit_result=None, fitter=None, pt2_nstates=None, sum_nstates=None, plot_type=None):
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

    for i in range(pt3_data_range[0], pt3_data_range[1]-1):
        temp1 = data_avg_dict_completed['sum_A3_fit_'+str(i)]
        gA.append(temp1.mean)
        gA_err.append(temp1.sdev)
        
        temp2 = data_avg_dict_completed['sum_V4_fit_'+str(i)]
        gV.append(temp2.mean)
        gV_err.append(temp2.sdev)

    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    ax.errorbar(np.arange(pt3_data_range[0], pt3_data_range[1]-1), np.array(gA), yerr=np.array(gA_err), marker='o', color="k", **errorp)
    print(gA)
    print(gA_err)
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
        
        fillx = np.arange(pt3_data_range[0], pt3_data_range[1], plot_space)[:int(-1 / plot_space)]
        ax.fill_between(fillx, gA_fit_y1, gA_fit_y2, color=blue, alpha=0.3, label='fit')
        print(fillx)
        print(gA_fit_y1)
        print(gA_fit_y2)
    
    
    ax.set_xlim([1.5, 13.5])
    ax.set_ylim([1, 1.4])
    ax.set_xlabel(t_label, **textp)
    ax.set_ylabel(oaeff_label, **textp)
    
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/soaeff{plot_type}.pdf", transparent=True)
    plt.show()
    
    plt.figure(figsize=figsize)
    ax = plt.axes(aspect)
    ax.errorbar(np.arange(pt3_data_range[0], pt3_data_range[1]-1), np.array(gV), yerr=np.array(gV_err), marker='o', color="k", **errorp)
    print(gV)
    print(gV_err)
    if fit_result != None and fitter != None and sum_nstates != None:
        gV_fit_y1 = np.array(gV_fit) + np.array(gV_fit_err)
        gV_fit_y2 = np.array(gV_fit) - np.array(gV_fit_err)
        
        fillx = np.arange(pt3_data_range[0], pt3_data_range[1], plot_space)[:int(-1 / plot_space)]
        ax.fill_between(fillx, gV_fit_y1, gV_fit_y2, color=blue, alpha=0.3)
        print(fillx)
        print(gV_fit_y1)
        print(gV_fit_y2)
    
    
    ax.set_xlim([1.5, 13.5])
    ax.set_ylim([1.0, 1.1])
    ax.set_xlabel(t_label, **textp)
    ax.set_ylabel(oveff_label, **textp)
    
    #plt.tight_layout(pad=0, rect=aspect)
    plt.savefig(f"./new_plots/soveff{plot_type}.pdf", transparent=True)
    plt.show()

    #print(gA[8])

    return [np.arange(pt3_data_range[0], pt3_data_range[1]-1), np.array(gA), np.array(gA_err)] # save this for data plot of half stat


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

    for n_ in range(n_range[1]):
        value['Q'].append([])
        value['logGBF'].append([])
        value['E0'].append([])
        value['E0_err'].append([])
        value['gA'].append([])
        value['gA_err'].append([])
        value['gV'].append([])
        value['gV_err'].append([])
        x.append([])


    for n in range(n_range[0], n_range[1]):
        for situation in situation_list:         
            nstate_dict = {}
            nstate_dict['2pt'] = situation.pt2_nstates
            nstate_dict['3pt'] = situation.pt3_nstates
            nstate_dict['sum'] = situation.sum_nstates
            
            if nstate_dict[nstate_name] == n:
                value['Q'][n].append(situation.Q_value)
                value['logGBF'][n].append(situation.log_GBF)
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
    best_t_ = best_t - t_range[0]

    #####################################################################################
    #####################################################################################
    # gA - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      #不同subplot共享x轴


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

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      #不同subplot共享x轴


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

    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_gV-'+tmin_name+'_tmin.pdf', transparent=True)

    #####################################################################################
    #####################################################################################
    # E0 - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      #不同subplot共享x，y轴

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

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmin_div)      #不同subplot共享x轴


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

    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_gA-'+tmin_name+'_tmin.pdf', transparent=True)
    
    #####################################################################################
    #####################################################################################
    # gV - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmin_div)      #不同subplot共享x轴


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

    plt.tight_layout(pad=30, rect=aspect)

    plt.savefig('./new_plots/'+fit_name+'_gV-'+tmin_name+'_tmin.pdf', transparent=True)

    #####################################################################################
    #####################################################################################
    # E0 - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmin_div)      #不同subplot共享x，y轴

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


    for n in range(n_range[0], n_range[1]):
        for situation in situation_list:         
            nstate_dict = {}
            nstate_dict['2pt'] = situation.pt2_nstates
            nstate_dict['3pt'] = situation.pt3_nstates
            nstate_dict['sum'] = situation.sum_nstates
            
            if nstate_dict[nstate_name] == n:
                value['Q'][n].append(situation.Q_value)
                value['logGBF'][n].append(situation.log_GBF)
                value['E0'][n].append(situation.E0)
                value['E0_err'][n].append(situation.E0_err)
                value['z0'][n].append(situation.z0)
                value['z0_err'][n].append(situation.z0_err)
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

    plt.tight_layout(pad=30, rect=aspect)

    #plt.savefig('./new_plots/'+fit_name+'_gA-'+tmin_name+'_tmin.pdf', transparent=True)
    
    #####################################################################################
    #####################################################################################
    # gV - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      #不同subplot共享x轴


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

    plt.tight_layout(pad=30, rect=aspect)

    #plt.savefig('./new_plots/'+fit_name+'_gV-'+tmin_name+'_tmin.pdf', transparent=True)

    #####################################################################################
    #####################################################################################
    # E0 - tmin
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      #不同subplot共享x，y轴

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

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin)      #不同subplot共享x，y轴

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
        value['Q'].append(situation.Q_value)
        value['E0'].append(situation.E0)
        value['E0_err'].append(situation.E0_err)
        value['gA'].append(situation.A300)
        value['gA_err'].append(situation.A300_err)
        value['gV'].append(situation.V400)
        value['gV_err'].append(situation.V400_err)

        tmax_dict = {}
        tmax_dict['2pt'] = situation.pt2_tmax
        tmax_dict['3pt'] = situation.pt3_A3_tsep_max
        tmax_dict['sum'] = situation.sum_A3_tsep_max

        x.append(tmax_dict[tmax_name])  # here is the varying parameter

    print(x)
    
    best_t_ = best_t - t_range[0]

    #####################################################################################
    #####################################################################################

    # gA - tmax
    fig=plt.figure()
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)      #不同subplot共享x轴
    ax1.set_ylabel(oa00_label, **textp)
    ax1.set_ylim(gA_ylim)

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

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)      #不同subplot共享x轴
    ax1.set_ylabel(ov00_label, **textp)
    ax1.set_ylim(gV_ylim)

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

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)      #不同subplot共享x，y轴
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

import h5py

import matplotlib.pyplot as plt
import numpy as np
import gvar as gv

possible_FHdatas = ['all', 'gA', 'gV']
possible_gAsets = ['A1', 'A2', 'A3', 'A4']
possible_gVsets = ['V1', 'V2', 'V3', 'V4']

# For the 3-point correlation functions, the datasets which switch between 2-d array data (for plotting and interacting with the data source) and 1-d data (for running the model-functions and fitting)
# Universal precondition for the below functions: x_1d, x_2d, y_1d, and y_2d, in the end, must have the
# same number of elements. Also, x_1d and y_1d have to have the same dimensions, as must x_2d and y_2d
def convert_1dto2d(x_1d, x_2d, y_1d):
    y_2d = np.array([np.array([None for j in range(len(x_2d[i]))]) for i in range(len(x_2d))])
    i = 0
    j = 0
    for k in range(len(x_1d)):
        y_2d[i][j] = y_1d[k]
        j += 1
        if j >= len(x_2d[i]):
            i +=1
            j = 0
            
    return y_2d

def convert_2dto1d(x_1d, x_2d, y_2d):
    y_1d = np.array([None for k in range(len(x_1d))])
    k = 0
    for i in range(len(x_2d)):
        for j in range(len(x_2d[i])):
            y_1d[k] = y_2d[i][j]
            k += 1
    
    return y_1d
    

# Gets a NumPy array of the 2-point data 
def get_2pt(fname):
    file = h5py.File(fname, "r")
    return_array = np.array(file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['spec']['ml0p0126']['proton']['px0_py0_pz0']['spin_up']) + np.array(file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['spec']['ml0p0126']['proton']['px0_py0_pz0']['spin_dn']) + np.array(file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['spec']['ml0p0126']['proton_np']['px0_py0_pz0']['spin_up']) + np.array(file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['spec']['ml0p0126']['proton_np']['px0_py0_pz0']['spin_dn'])
    
    return np.array([[[np.real(return_array[i][j][k][0]) for k in range(2)] for j in range(64)] for i in range(1053)])/4.

# Gets a 3-point data array with a specified gA or gV dataset and a tsep range
def get_3pt(fname, dataset, tsep_min, tsep_max):
    file = h5py.File(fname, "r")
    if dataset in possible_gAsets:
        spindn_coeff = -1
    elif dataset in possible_gVsets:
        spindn_coeff = 1
    else:
        raise ValueError("the entry for the g-value dataset must be A1, A2, A3, A4, V1, V2, V3, or V4")

    return_array = (np.array([file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['formfac']['ml0p0126']['proton_UU_up_up_tsep_{}_sink_mom_px0_py0_pz0'.format(i)][dataset]['px0_py0_pz0']['local_curr'] for i in range(tsep_min, tsep_max)]) - np.array([file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['formfac']['ml0p0126']['proton_DD_up_up_tsep_{}_sink_mom_px0_py0_pz0'.format(i)][dataset]['px0_py0_pz0']['local_curr'] for i in range(tsep_min, tsep_max)])) + spindn_coeff*(np.array([file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['formfac']['ml0p0126']['proton_UU_dn_dn_tsep_{}_sink_mom_px0_py0_pz0'.format(i)][dataset]['px0_py0_pz0']['local_curr'] for i in range(tsep_min, tsep_max)]) - np.array([file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['formfac']['ml0p0126']['proton_DD_dn_dn_tsep_{}_sink_mom_px0_py0_pz0'.format(i)][dataset]['px0_py0_pz0']['local_curr'] for i in range(tsep_min, tsep_max)])) - (np.array([file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['formfac']['ml0p0126']['proton_np_UU_up_up_tsep_{}_sink_mom_px0_py0_pz0'.format(-i)][dataset]['px0_py0_pz0']['local_curr'] for i in range(tsep_min, tsep_max)]) - np.array([file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['formfac']['ml0p0126']['proton_np_DD_up_up_tsep_{}_sink_mom_px0_py0_pz0'.format(-i)][dataset]['px0_py0_pz0']['local_curr'] for i in range(tsep_min, tsep_max)])) - spindn_coeff*(np.array([file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['formfac']['ml0p0126']['proton_np_UU_dn_dn_tsep_{}_sink_mom_px0_py0_pz0'.format(-i)][dataset]['px0_py0_pz0']['local_curr'] for i in range(tsep_min, tsep_max)]) - np.array([file['gf1p0_w3p0_n30_M51p2_L58_a1p5']['formfac']['ml0p0126']['proton_np_DD_dn_dn_tsep_{}_sink_mom_px0_py0_pz0'.format(-i)][dataset]['px0_py0_pz0']['local_curr'] for i in range(tsep_min, tsep_max)]))
    
    if dataset in possible_gAsets:
        return np.array([[[np.imag(return_array[i][j][k]) for k in range(64)] for i in range(10)] for j in range(1053)])/4.
        
    else:
        return np.array([[[np.real(return_array[i][j][k]) for k in range(64)] for i in range(10)] for j in range(1053)])/4.
        
# Gets a FH dataset
def get_FH(fname, dataset):
    file = h5py.File(fname, "r")
    if dataset not in possible_FHdatas:
        raise ValueError("the entry for the FH data set must be either gA or gV")
        
    elif dataset=='all':
        return np.transpose(np.array([file['gA_ss'], file['gA_ps'], file['gV_ss'], file['gV_ps']]), axes=[1,2,0])/4.
    else:
        return np.transpose(np.array([file['{}_ss'.format(dataset)], file['{}_ps'.format(dataset)]]), axes=[1,2,0])/4.
    
# Designed to create a plot of raw data, without fit functions
# Precondition: xdat is the x data, ydat is the y-data, xtitle is the x-axis title, ytitle is the y-axis title, graph title is the graph's title, nplots is the number of plots, legendlabels is the labels that go into the legend, logy is used to toggle the log-scale y-axis option
# Precondition: every parameter except the axis and graph titles must be an array-like type with nplots elements or its default value.
def raw_data_plot(graphtitle, xtitle, ytitle, xdat, ydat, ylim = None, nplots = 1, yerror = None, legendlabels = None, logy = False, filename = None):
    fig = plt.figure(figsize = (10,7))
    ax = plt.axes()
    if logy:
        ax.semilogy()
    
    if legendlabels is None:
        if yerror is None:
            for i in range(nplots):
                if nplots != 1:
                    ax.errorbar(xdat[i], ydat[i], fmt='.')
                    if ylim is not None:
                        plt.ylim(ylim[0], ylim[1])
                        
                else:
                    ax.errorbar(xdat, ydat, fmt='.')
                    if ylim is not None:
                        plt.ylim(ylim[0], ylim[1])
        
        else:
            for i in range(nplots):
                if nplots != 1:
                    ax.errorbar(xdat[i], ydat[i], yerr=yerror[i], fmt='.')
                    if ylim is not None:
                        plt.ylim(ylim[0], ylim[1])
                else:
                    ax.errorbar(xdat, ydat, yerr=yerror, fmt='.')
                    if ylim is not None:
                        plt.ylim(ylim[0], ylim[1])
                    
    else:
        if yerror is None:
            for i in range(nplots):
                if nplots != 1:
                    ax.errorbar(xdat[i], ydat[i], fmt='.', label=legendlabels[i])
                    if ylim is not None:
                        plt.ylim(ylim[0], ylim[1])
                else:
                    ax.errorbar(xdat, ydat, fmt='.', label=legendlabels[i])
                    if ylim is not None:
                        plt.ylim(ylim[0], ylim[1])
        
        else:
            for i in range(nplots):
                if nplots != 1:
                    ax.errorbar(xdat[i], ydat[i], yerr=yerror[i], fmt='.', label=legendlabels[i])
                    if ylim is not None:
                        plt.ylim(ylim[0], ylim[1])
                    
                else:
                    ax.errorbar(xdat, ydat, yerr=yerror, fmt='.', label=legendlabels[i])
                    if ylim is not None:
                        plt.ylim(ylim[0], ylim[1])
                
        plt.legend()
        
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(graphtitle)
    
    if filename is not None:
        plt.savefig(filename)
    
    plt.show()
    
    plt.clf()

#ax.errorbar(xdat[i], ydat[i], ls='-', )
#ax.fill_between(xfit[i], yfit_upperbound[i], yfit_lowerbound[i], alpha=0.3)
    
#Plots out the fit-function with the data
# Precondition: each fit-result must correspond to a data result
def fit_data_plot(graphtitle, xtitle, ytitle, xdat, ydat, xfit, yfit, yfit_upperbound, yfit_lowerbound, colors, g00=None, ylim = None, nplots = 1, ydat_error = None, datalabels = None, fitlabels = None, logy = False, filename = None):
    fig = plt.figure(figsize = (10,7))
    ax = plt.axes()
    if logy:
        ax.semilogy()
    
    if fitlabels is None:
        if datalabels is None:
            if ydat_error is None:
                for i in range(nplots):
                    if nplots != 1:
                        ax.errorbar(xdat[i], ydat[i], fmt='.', color = colors[i])
                        ax.errorbar(xfit[i], yfit[i], ls='-', color = colors[i])
                        ax.fill_between(xfit[i], yfit_upperbound[i], yfit_lowerbound[i], color=colors[i], alpha=0.3)
                        
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                        
                    else:
                        ax.errorbar(xdat, ydat, fmt='.', color = colors)
                        ax.errorbar(xfit, yfit, ls='-', color = colors)
                        ax.fill_between(xfit, yfit_upperbound, yfit_lowerbound, color=colors, alpha=0.3)
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                    
                    
        
            else:
                for i in range(nplots):
                    if nplots != 1:
                        ax.errorbar(xdat[i], ydat[i], yerr=ydat_error[i], fmt='.', color = colors[i])
                        ax.errorbar(xfit[i], yfit[i], ls='-', color = colors[i])
                        ax.fill_between(xfit[i], yfit_upperbound[i], yfit_lowerbound[i], color=colors[i], alpha=0.3)
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                    else:
                        ax.errorbar(xdat, ydat, yerr=ydat_error, fmt='.', color = colors)
                        ax.errorbar(xfit, yfit, ls='-', color = colors)
                        ax.fill_between(xfit, yfit_upperbound, yfit_lowerbound, color=colors, alpha=0.3)
                        if g00 is not None:
                            ax.axhspan(g00.mean - g00.sdev, g00.mean + g00.sdev, color='cyan', alpha=0.2)
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                    
        else:
            if ydat_error is None:
                for i in range(nplots):
                    if nplots != 1:
                        ax.errorbar(xdat[i], ydat[i], fmt='.', label=datalabels[i], color = colors[i])
                        ax.errorbar(xfit[i], yfit[i], ls='-', color = colors[i])
                        ax.fill_between(xfit[i], yfit_upperbound[i], yfit_lowerbound[i], color=colors[i], alpha=0.3)
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                    else:
                        ax.errorbar(xdat, ydat, fmt='.', label=datalabels, color = colors)
                        ax.errorbar(xfit, yfit, ls='-', color = colors)
                        ax.fill_between(xfit, yfit_upperbound, yfit_lowerbound, color=colors, alpha=0.3)
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
        
            else:
                for i in range(nplots):
                    if nplots != 1:
                        ax.errorbar(xdat[i], ydat[i], yerr=ydat_error[i], fmt='.', color = colors[i], label=datalabels[i])
                        ax.errorbar(xfit[i], yfit[i], ls='-', color = colors[i])
                        ax.fill_between(xfit[i], yfit_upperbound[i], yfit_lowerbound[i], color=colors[i], alpha=0.3)
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                    
                    else:
                        ax.errorbar(xdat, ydat, yerr=ydat_error, fmt='.', color = colors, label=datalabels)
                        ax.errorbar(xfit, yfit, ls='-', color = colors)
                        ax.fill_between(xfit, yfit_upperbound, yfit_lowerbound, color=colors, alpha=0.3)
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                
            plt.legend()
        
    else:
        if datalabels is None:
            if ydat_error is None:
                for i in range(nplots):
                    if nplots != 1:
                        ax.errorbar(xdat[i], ydat[i], fmt='.', color = colors[i])
                        ax.errorbar(xfit[i], yfit[i], ls='-', color = colors[i])
                        ax.fill_between(xfit[i], yfit_upperbound[i], yfit_lowerbound[i], color=colors[i], alpha=0.3, label = fitlabels[i])
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                        
                    else:
                        ax.errorbar(xdat, ydat, fmt='.', color = colors)
                        ax.errorbar(xfit, yfit, ls='-', color = colors)
                        ax.fill_between(xfit, yfit_upperbound, yfit_lowerbound, color=colors, alpha=0.3, label = fitlabels)
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
        
            else:
                for i in range(nplots):
                    if nplots != 1:
                        ax.errorbar(xdat[i], ydat[i], yerr=ydat_error[i], fmt='.', color = colors[i])
                        ax.errorbar(xfit[i], yfit[i], ls='-', color = colors[i])
                        ax.fill_between(xfit[i], yfit_upperbound[i], yfit_lowerbound[i], color=colors[i], alpha=0.3, label = fitlabels[i])
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                    else:
                        ax.errorbar(xdat, ydat, yerr=ydat_error, fmt='.', color = colors)
                        ax.errorbar(xfit, yfit, ls='-', color = colors)
                        ax.fill_between(xfit, yfit_upperbound, yfit_lowerbound, color=colors, alpha=0.3, label = fitlabels)
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                    
            plt.legend()
            
        else:
            if ydat_error is None:
                for i in range(nplots):
                    if nplots != 1:
                        ax.errorbar(xdat[i], ydat[i], fmt='.', label=datalabels[i], color = colors[i])
                        ax.errorbar(xfit[i], yfit[i], ls='-', color = colors[i])
                        ax.fill_between(xfit[i], yfit_upperbound[i], yfit_lowerbound[i], color=colors[i], alpha=0.3, label = fitlabels[i])
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                    else:
                        ax.errorbar(xdat, ydat, fmt='.', label=datalabels, color = colors)
                        ax.errorbar(xfit, yfit, ls='-', color = colors)
                        ax.fill_between(xfit, yfit_upperbound, yfit_lowerbound, color=colors, alpha=0.3, label = fitlabels)
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
        
            else:
                for i in range(nplots):
                    if nplots != 1:
                        ax.errorbar(xdat[i], ydat[i], yerr=ydat_error[i], fmt='.', color = colors[i], label=datalabels[i])
                        ax.errorbar(xfit[i], yfit[i], ls='-', color = colors[i])
                        ax.fill_between(xfit[i], yfit_upperbound[i], yfit_lowerbound[i], color=colors[i], alpha=0.3, label = fitlabels[i])
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                    
                    else:
                        ax.errorbar(xdat, ydat, yerr=ydat_error, fmt='.', color = colors, label=datalabels)
                        ax.errorbar(xfit, yfit, ls='-', color = colors)
                        ax.fill_between(xfit, yfit_upperbound, yfit_lowerbound, color=colors, alpha=0.3, label = fitlabels)
                        if ylim is not None:
                            plt.ylim(ylim[0], ylim[1])
                
            plt.legend()
    
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(graphtitle)
    
    if filename is not None:
        plt.savefig(filename)
    
    plt.show()
    
    plt.clf()

# Like the regular raw_data_plot, but only plots one set of data
#def raw_data_plot(graphtitle, xtitle, ytitle, xdat, ydat, yerror = None, legendlabels = None, logy = False):
    #fig = plt.figure(figsize = (10,7))
    #ax = plt.axes()
    #if logy:
    #    ax.semilogy()

    #if legendlabels is None:
    #    if yerror is None:
    #        ax.errorbar(xdat, ydat, fmt='.')
    
    #    else:
    #        ax.errorbar(xdat, ydat, yerr=yerror, fmt='.')

    #else:
     #   if yerror is None:
      #      ax.errorbar(xdat, ydat, fmt='.', label=legendlabels[i])
    
      #  else:
      #      ax.errorbar(xdat, ydat, yerr=yerror, fmt='.', label=legendlabels[i])
            
      #  plt.legend()
    
    #plt.xlabel(xtitle)
    #plt.ylabel(ytitle)
    #plt.title(graphtitle)

    #plt.show()

    #plt.clf()

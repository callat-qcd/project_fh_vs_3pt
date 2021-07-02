#!/usr/bin/env python3
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import gvar as gv
if gv.__version__ != '11.5.2':
    print('your gvar is version ', gv.__version__)
    sys.exit('For now, you must have gvar version 11.5.2 to be able to read the results')

s_mev = 197.3 / 0.08730
s_gev = s_mev / 1000


# set the enery levels
mpi = 0.14073
mN  = 0.4904

fit_types = ['spec_results','spec_results_ga','spec_results_gv','spec_results_pt2']

models = ['SqW', 'HO', '1/n', '1/n2']
models = ['HO', '1/n', '1/n2']
file_m = {'SqW':'sw_result', 'HO':'ho_result', '1/n':'n_result', '1/n2':'n2_result'}
fits = dict()
for m in fit_types:
    fits[m] = dict()
    for k in ['logGBF','w','Q','E0','E1','E2','E3','E4','gA','z0','pdE1', 'pdE2', 'pdE3', 'pdE4']:
        fits[m][k] = []
for t in ['3','4','5','6','7']:
    for m in fit_types:
        f = gv.load('data/'+m+'/'+file_m['HO']+t) # [prior, posterior, Q, logGBF]
        fits[m]['logGBF'].append(f[3])
        fits[m]['Q'].append(f[2])
        fits[m]['E0'].append(f[1]['E0'])
        for n in [1,2,3,4]:
            tmp = f[1]['E0']
            for l in range(1,n+1):
                tmp += f[1]['dE'+str(l)]
            fits[m]['E'+str(n)].append(tmp)
            fits[m]['pdE'+str(n)].append(f[0]['dE'+str(n)])
        #print(f[0]['dE2'])
        fits[m]['gA'].append(f[1]['A3_00'])
        fits[m]['z0'].append(f[1]['z0'])
'''
# make weights
for t in range(5):
    lGBF_rel = []
    for m in models:
        lGBF_rel.append(fits[m]['logGBF'][t])
    #print(lGBF_rel)
    lGBF_rel = np.array(lGBF_rel) - max(lGBF_rel)
    #print(lGBF_rel)
    w = np.exp(lGBF_rel)
    w = w / w.sum()
    #print(w)
    for i_m, m in enumerate(models):
        fits[m]['w'].append(w[i_m])
'''

# excited state energy gaps
def dE_mod(n,mod):
    p = gv.BufferDict()
    if mod == 'SqW':
        p['log(dE_%d)' %n] = gv.gvar(np.log(2*mpi * (2*n-1)), 0.5 / (2*n-1) )
    elif mod == 'HO':
        p['log(dE_%d)' %n] = gv.gvar(np.log(2*mpi), .5)
    elif mod == '1/n':
        p['log(dE_%d)' %n] = gv.gvar(np.log(2*mpi / n), .5*n)
    elif mod == '1/n2':
        p['log(dE_%d)' %n] = gv.gvar(np.log(2*mpi / n**2), .5*n**2)
    dEn = p['dE_%d' %n]
    #print(mod, n, dEn)
    return dEn

def Npi(mN,mpi,L):
    E = dict()
    n_max = 20
    for px in range(n_max):
        for py in range(n_max):
            for pz in range(n_max):
                p = np.array([px,py,pz])
                psq = int(np.sum(p**2))
                if not np.all(p==0) and psq not in E:
                    E[psq] = np.sqrt(mN**2 +(2*np.pi/L)**2 * psq) + np.sqrt(mpi**2 +(2*np.pi/L)**2 * psq)
    return E

# plot fit results of En
fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
figsize  = (2*fig_width, fig_width )
aspect=[0.1, 0.1, 0.895, 0.895]
textp = {"fontsize": 18}
labelp = {"labelsize": 14}
labelsize = {'size': 16}

plt.rcParams.update({"text.usetex": True})

red = "#FF6F6F"
peach = "#FF9E6F"
orange = "#FFBC6F"
sunkist = "#FFDF6F"
yellow = "#FFEE6F"
lime = "#CBF169"
green = "#5CD25C"
turquoise = "#4AAB89"
blue = "#508EAD"
grape = "#635BB1"
violet = "#7C5AB8"
fuschia = "#C3559F"

color_list = [violet, green, fuschia, yellow, blue, orange, turquoise, red]

plt.ion()
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

gridspec_tmin = {'height_ratios': [10, 2, 2],
    'left': 0.06, 'bottom': 0.1, 'right': 0.945, 'top': 0.99, 'hspace':0}

fig, (ax_es, ax_e0, ax_Q) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin, figsize=figsize)

# add N-pi levels
E_Npi = Npi(mN=mN, mpi=mpi, L=32)
ax_es.axhline(E_Npi[1]+10, linestyle=':', color='k', label=r'$N(\mathbf{q})\pi(\mathbf{-q})$')
for psq in E_Npi:
    #ax_es.axhline(E_Npi[psq], linestyle=(0, (5, 10)), color='k',alpha=.3, label=lbl)
    ax_es.axhline(E_Npi[psq], linestyle=':', color='k',alpha=.3)
# add lowest N-pipi level
ax_es.axhline(mN+2*mpi, linestyle=(0, (3, 5, 1, 5)), color='k',alpha=.2)

m_lbl = {'spec_results': '2pt+gA+gV', 'spec_results_pt2':'2pt', 'spec_results_ga':'2pt+gA','spec_results_gv':'2pt+gV'}

shift = {'SqW':-.3, 'HO':-.1, '1/n':.1, '1/n2':.3}
n_clr = {'SqW':red, 'HO':orange, '1/n':green, '1/n2':blue}
mrkr  = {'SqW':'s', 'HO':'o', '1/n':'X', '1/n2':'D'}
p_width = 0.025

shift = {'spec_results': -.3, 'spec_results_pt2':.3, 'spec_results_ga':-.1,'spec_results_gv':.1}
n_clr = {'spec_results': red, 'spec_results_pt2':orange, 'spec_results_ga':green,'spec_results_gv':blue}
mrkr  = {'spec_results': 's', 'spec_results_pt2':'o', 'spec_results_ga':'X','spec_results_gv':'D'}
p_width = 0.025


for i_t,t in enumerate([3,4,5,6,7]):
    for m in fit_types:
        mfc = n_clr[m]
        clr = 'k'
        e0 = fits[m]['E0'][i_t]

        if t==3 and m == 'spec_results':
            mfc = 'k'
            clr = n_clr[m]

        t_0  = t + shift[m]
        ax_e0.errorbar(t_0, e0.mean, yerr=e0.sdev,
                    color=clr, marker=mrkr[m], mfc=mfc)
        #gA = fits[m]['gA'][i_t]
        #ax_gA.errorbar(t_0, gA.mean, yerr=gA.sdev,
        #            color=clr, marker=mrkr[m], mfc=mfc)
        ax_Q.plot(t_0, fits[m]['Q'][i_t], color=n_clr[m], marker=mrkr[m], mfc=mfc)
        #ax_w.plot(t_0, fits[m]['w'][i_t], color=clr, marker=mrkr[m], mfc=mfc)

        if t==3 and m == 'spec_results':
            ax_e0.axhspan(e0.mean-e0.sdev, e0.mean+e0.sdev, color=n_clr[m], alpha=.3)
            #ax_gA.axhspan(gA.mean-gA.sdev, gA.mean+gA.sdev, color=n_clr[m], alpha=.3)

        for n in [1,2,3,4]:
            # plot prior
            #En = fits[m]['E'+str(n-1)][i_t] + dE_mod(n, m)
            if n == 2 and False:
                print(En)
                print(fits[m]['E1'][i_t])
                print(fits[m]['pdE2'][i_t])
            En = fits[m]['E'+str(n-1)][i_t] + fits[m]['pdE'+str(n)][i_t]
            n_base  = t_0 + 2*p_width*(n-1)
            p_range = np.arange(n_base-p_width, n_base+p_width+.001, .001)
            if n == 2 and False:
                print(En)
            ax_es.fill_between(p_range, En.mean-En.sdev, En.mean+En.sdev,
                    color=n_clr[m], alpha=.5)

            # plot posterior
            lbl = ""
            if i_t == 1 and n == 1:
                lbl = m_lbl[m]
            en = fits[m]['E'+str(n)][i_t]
            if t==3 and m == 'spec_results':
                ax_es.axhspan(en.mean-en.sdev, en.mean+en.sdev, color=n_clr[m], alpha=.3)
            ax_es.errorbar(t_0+2*p_width*(n-1), en.mean, yerr=en.sdev,
                    color=clr, marker=mrkr[m], mfc=mfc, label=lbl, linestyle='None')

            if t in [5,6] and n==1 and m in ['spec_results','spec_results_pt2']:
                print(m,t,en)

ax_e0.set_ylim(0.481, 0.499)
ax_e0.set_ylabel(r'$a_{09} E_0$',**textp)
ax_e0.tick_params(axis='both', which='major', **labelp)
ax_e0r = ax_e0.twinx()
ax_e0r.set_ylim(ax_e0.get_ylim()[0]*s_gev, ax_e0.get_ylim()[1]*s_gev)
ax_e0r.set_yticks([s_gev*t for t in ax_e0.get_yticks()[1:-1]])
ax_e0r.tick_params(axis='both', which='major', **labelp)
ax_e0r.set_yticklabels(["%.2f" %t for t in ax_e0r.get_yticks()])
ax_e0r.set_ylabel(r'$E_0 / {\rm GeV}$', **textp)

#ax_gA.set_ylim(1.16, 1.34)
#ax_gA.set_ylabel(r'$\mathring{g}_A$', **textp)
#ax_gA.tick_params(axis='both', which='major', **labelp)
#ax_gA.set_yticks([1.19, fits['HO']['gA'][0].mean, 1.31])
#ax_gA.set_yticklabels(["%.2f" %e for e in ax_gA.get_yticks()])

#ax_w.set_ylim(0,1)
#ax_w.set_yticks([0.1,0.75])
#ax_w.set_ylabel(r'$w$', **textp)
ax_Q.set_xlabel(r'$t_{\rm sep}^{\rm min} : C_2$', **textp)
#ax_w.tick_params(axis='both', which='major', **labelp)
ax_Q.set_ylim(0,1)
ax_Q.set_yticks([0.1,0.75])
ax_Q.set_ylabel(r'$Q$', **textp)
ax_Q.tick_params(axis='both', which='major', **labelp)


ax_es.set_ylim(0.51, 3.8)
ax_es.legend(loc=1,ncol=5,**textp)
ax_es.tick_params(axis='both', which='major', **labelp)
ax_es.set_ylabel(r'$a_{09}E_n$', **textp)
ax_es.set_yticks([fits['spec_results']['E%d' %n][0].mean for n in range(1,5)])
ax_es.set_yticklabels(["%.3f" %e for e in ax_es.get_yticks()])
ax_esr = ax_es.twinx()
ax_esr.set_ylim(ax_es.get_ylim()[0]*s_gev, ax_es.get_ylim()[1]*s_gev)
ax_esr.tick_params(axis='both', which='major', **labelp)
ax_esr.set_yticks([s_gev*t for t in ax_es.get_yticks()])
ax_esr.set_yticklabels(["%.2f" %t for t in ax_esr.get_yticks()])
ax_esr.set_ylabel(r'$E_n / {\rm GeV}$', **textp)

if not os.path.exists('new_plots'):
    os.makedirs('new_plots')
plt.savefig('new_plots/es_model_sensitivity_2pt_g.pdf', transparent=True)

plt.ioff()
plt.show()

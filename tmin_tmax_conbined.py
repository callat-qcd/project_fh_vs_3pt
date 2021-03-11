# %%
import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  
import os

from fh_db.ff_fit.models import FF_fit as ff_fit_FF_fit
os.environ["DJANGO_ALLOW_ASYNC_UNSAFE"] = "true" # for jupyter to use database

plt.rcParams.update({"text.usetex": True})

# %%
grey = "#808080" 
red = "#FF6F6F" # nstates = 8
peach = "#FF9E6F" # nstates = 7
orange = "#FFBC6F" # nstates = 6
sunkist = "#FFDF6F"
yellow = "#FFEE6F" 
lime = "#CBF169" 
green = "#5CD25C" 
turquoise = "#4AAB89" # nstates = 5
blue = "#508EAD" # nstates = 4
grape = "#635BB1" # nstates = 3
violet = "#7C5AB8"
fuschia = "#C3559F" # nstates = 2

color_list = [fuschia, grape, blue, turquoise, orange, peach, red]

plt.rcParams.update({"text.usetex": True})
fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
fig_size  = (fig_width, fig_width / gr)
plt_axes  = [0.15,0.15,0.845,0.845]
fs_text   = 16 # font size of text
fs_leg    = 16 # legend font size
tick_size = 16 # tick size
plt.rcParams['figure.figsize'] = fig_size

gridspec_tmax = {'height_ratios': [3, 1], 'left': 0.12, 'right': 0.95, 'bottom': 0.15, 'top': 0.95}
textp = {"fontsize": 14}
labelp = {"labelsize": 14}
errorp = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1}
errorb = {"markersize": 5, "linestyle": "none", "capsize": 3, "elinewidth": 1}


oa00_label = r"$\mathring{g}_A$"
c_tmax = r"$t_{\textrm{sep}}^{\textrm{max}}$"
c_tmin = r"$t_{\textrm{sep}}^{\textrm{min}}$"
q_label = r"$Q$"

gA_ylim=[1.01, 1.349]




# %%
####################################
########## 23s 3pt/sum tmax #############
#####################################
fit_name = '23s'

value={}
value['Q']=[]
value['gA']=[]
value['gA_err']=[]

x=[]


situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', fit_type='continuous', include_2pt=True, include_3pt=True, include_sum=True,
pt2_tmin=3, pt2_tmax=18, pt2_nstates=5,
pt3_A3_tsep_max=15, id_num=3,
sum_A3_tsep_max=14, sum_tau_cut=1)

for situation in situation_list:
    x.append(situation.pt3_A3_tsep_min)

x.sort()
print(x)
l1 = len(x)

for i in range(l1):
    for situation in situation_list:
        if situation.pt3_A3_tsep_min == x[i]:
            value['Q'].append(situation.Q_value)
            value['gA'].append(situation.A300)
            value['gA_err'].append(situation.A300_err)

situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', fit_type='continuous', include_2pt=True, include_3pt=True, include_sum=True,
pt2_tmin=3, pt2_tmax=18, pt2_nstates=5,
pt3_A3_tsep_min=3, id_num=3,
sum_A3_tsep_min=3, sum_tau_cut=1, sum_nstates=5)

for situation in situation_list:
    x.append(situation.pt3_A3_tsep_max)

x.sort()
print(x)
l2 = len(x)

for i in range(l1, l2):
    for situation in situation_list:
        if situation.pt3_A3_tsep_max == x[i]:
            value['Q'].append(situation.Q_value)
            value['gA'].append(situation.A300)
            value['gA_err'].append(situation.A300_err)

    x[i] = x[i] - 1 # tmax is not included in fits

print(value['gA'][1])
print(value['gA'][12])

best_gA = value['gA'][1]
best_gA_err = value['gA_err'][1]
best_Q = value['Q'][1]

t_range = [2, 15]
nstate = 5

# %%
fig=plt.figure(figsize=fig_size)
plt.rcParams['figure.figsize'] = fig_size
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_tmax)      
ax1.set_ylabel(oa00_label, fontsize=fs_text)
ax1.set_ylim(gA_ylim)

ax1.tick_params(direction='in', labelsize=tick_size)

ax1.vlines(7.5, gA_ylim[0], gA_ylim[1], colors = "k", linestyles = "dashed")

ax1.errorbar(np.array(x), np.array(value['gA']), yerr=np.array(value['gA_err']), marker='o', color=color_list[nstate-2], **errorp)

#best fit
ax1.fill_between(np.arange(1.5, 15.5), (best_gA + best_gA_err)*np.ones([14]), (best_gA - best_gA_err)*np.ones([14]), color=color_list[nstate-2], alpha=0.2)

ax1.errorbar(np.array([3]), np.array([best_gA]), yerr=np.array([best_gA_err]), marker='o', mfc=color_list[nstate-2], color=color_list[nstate-2], **errorb)

ax1.errorbar(np.array([14]), np.array([best_gA]), yerr=np.array([best_gA_err]), marker='o', mfc=color_list[nstate-2], color=color_list[nstate-2], **errorb)


ax2.set_ylabel(q_label, fontsize=fs_text)
ax2.set_ylim([0, 1.1])

ax2.vlines(7.5, 0, 1.1, colors = "k", linestyles = "dashed")

ax2.plot(np.arange(t_range[0] - 0.5, t_range[1] + 0.5, 1), 0.1 * np.ones([t_range[1] - t_range[0] + 1]), 'r--')

ax2.scatter(np.array(x), np.array(value['Q']), marker='o', c='', edgecolors=color_list[nstate-2])

#best fit
ax2.scatter(np.array([3]), np.array([best_Q]), marker='o', c=color_list[nstate-2])

ax2.scatter(np.array([14]), np.array([best_Q]), marker='o', c=color_list[nstate-2])

ax2.tick_params(direction='in', labelsize=tick_size)

plt.subplots_adjust(wspace=0, hspace=0)
plt.xlabel(c_tmin+'\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ '+c_tmax, fontsize=fs_text)
plt.xlim([1.5, 14.5])
plt.tight_layout(pad=30, rect=plt_axes)


plt.savefig('./new_plots/'+fit_name+'_gA-both_tmin_tmax_conbined.pdf', transparent=True)

# %%

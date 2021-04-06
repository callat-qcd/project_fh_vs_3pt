# %%
import h5py as h5 
import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  

from fh_db.ff_fit.models import FF_fit as ff_fit_FF_fit
os.environ["DJANGO_ALLOW_ASYNC_UNSAFE"] = "true" # for jupyter to use database

plt.rcParams.update({"text.usetex": True})

# %%
grey = "#808080" # nstates = 1
red = "#FF6F6F" # nstates = 2
peach = "#FF9E6F" # nstates = 3
orange = "#FFBC6F" # nstates = 4
sunkist = "#FFDF6F"
yellow = "#FFEE6F" 
lime = "#CBF169"
green = "#5CD25C" # nstates = 5
turquoise = "#4AAB89"
blue = "#508EAD" # nstates = 6
grape = "#635BB1"
violet = "#7C5AB8" # nstates = 7
fuschia = "#C3559F"

plt.rcParams.update({"text.usetex": True})
fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
fig_size  = (fig_width, fig_width / gr)
plt_axes  = [0.12,0.15,0.87,0.83]
fs_text   = 18 # font size of text
fs_text_gA = 20 # font size of ga label
fs_leg    = 16 # legend font size
tick_size = 16 # tick size
plt.rcParams['figure.figsize'] = fig_size


gridspec_prior_width = {'height_ratios': [3, 1], 'left': 0.12, 'right': 0.99, 'bottom': 0.15, 'top': 0.98}
textp = {"fontsize": 18}
labelp = {"labelsize": 18}
errorp = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1}
errorb = {"markersize": 5, "linestyle": "none", "capsize": 3, "elinewidth": 1}

# labels
q_label = r"$Q$"
ga_label = r"$\mathring{g}_A$"
gA_ylim = [1.01, 1.349]

# %%
# 2pt+3pt+sum varying g.s. and 1 ex prior width

value={}
value['Q']=[]
value['E0']=[]
value['E0_err']=[]
value['A3']=[]
value['A3_err']=[]
value['V4']=[]
value['V4_err']=[]

# gs01, gs05, gs1, gs2, gs10
for hexcode in ['a6d7c8abd3ab957e15af222a35f1f618', '70101230395214272d3ea3a896359b68', '8e08f23bc983bf0fa9778157733d8235', 'a68192fce72d0097536b9f00721ce302', '722949c46e4fd76d4b718d12631a1653']:
    for situation in ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode=hexcode,
    include_2pt=True, include_3pt=True, include_sum=True, 
    pt2_tmax=18, pt2_nstates=5,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=9,
    sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5):
                value['Q'].append(situation.Q_value)
                value['E0'].append(situation.E0)
                value['E0_err'].append(situation.E0_err)
                value['A3'].append(situation.A300)
                value['A3_err'].append(situation.A300_err)
                value['V4'].append(situation.V400)
                value['V4_err'].append(situation.V400_err)
                print(situation.prior_hexcode)
            

value['Q_']=[]
value['E0_']=[]
value['E0_err_']=[]
value['A3_']=[]
value['A3_err_']=[]
value['V4_']=[]
value['V4_err_']=[]

# gs and 1 ex
for hexcode in ['8b7225f794495bc718ed3a8c0ce5b26c', 'eb82329129cbfa975b81b08755b0d4dd', '8e08f23bc983bf0fa9778157733d8235', '1e8538003f62d1f766972446f1b08c27', '73ae4b7fa602c032c3e5f2c476d919f4']:
    for situation in ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode=hexcode,
    include_2pt=True, include_3pt=True, include_sum=True, 
    pt2_tmax=18, pt2_nstates=5,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=9,
    sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5):
                value['Q_'].append(situation.Q_value)
                value['E0_'].append(situation.E0)
                value['E0_err_'].append(situation.E0_err)
                value['A3_'].append(situation.A300)
                value['A3_err_'].append(situation.A300_err)
                value['V4_'].append(situation.V400)
                value['V4_err_'].append(situation.V400_err)
                print(situation.prior_hexcode)
            

#####################################################################################
#####################################################################################
# gA - prior width
fig=plt.figure()
plt.rcParams['figure.figsize'] = fig_size
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_prior_width)      #不同subplot共享x，y轴
ax1.set_ylabel(ga_label, fontsize=20)
ax1.set_ylim(gA_ylim)
ax1.tick_params(axis='both', which='major', **labelp)


mask = [0, 1, 3, 4]
ax1.errorbar([-0.05, 0.95, 2.95, 3.95], np.array(value['A3'])[mask], yerr=np.array(value['A3_err'])[mask],label='g.s.', marker="o", color=red, **errorp)
ax1.errorbar([0.05, 1.05, 3.05, 4.05], np.array(value['A3_'])[mask], yerr=np.array(value['A3_err_'])[mask], label='g.s. and 1 e.s.', marker="o", color=blue, **errorp)

ax1.fill_between(np.arange(- 0.5, 5.5, 1), (value['A3'][2]+value['A3_err'][2])*np.ones([6]), (value['A3'][2]-value['A3_err'][2])*np.ones([6]), color=green, alpha=0.3)
# highest logGBF
ax1.errorbar(np.array([2]), np.array([value['A3'][2]]), yerr=np.array([value['A3_err'][2]]), color=green, marker="o", label="best fit", **errorb)
#ax1.errorbar(np.array([2])+0.1, np.array([value['A3_'][2]]), yerr=np.array([value['A3_err_'][2]]), marker='o', markersize=5, mfc='r', linestyle='none', capsize=3, color='r', elinewidth=0.8)

ax1.legend(loc='lower center', ncol=3, fontsize=14)

ax2.set_ylabel(q_label, **textp)
ax2.set_ylim([0, 1.1])
ax2.plot(np.arange( - 0.5, 5 + 0.5, 1), 0.1 * np.ones([6]), 'r--')


ax2.scatter([-0.05, 0.95, 2.95, 3.95], np.array(value['Q'])[mask], marker='o', c='', edgecolors=red)
ax2.scatter([0.05, 1.05, 3.05, 4.05], np.array(value['Q_'])[mask], marker='o', c='', edgecolors=blue)

# highest logGBF
ax2.scatter(np.array([2]), np.array([value['Q'][2]]), marker="o", c=green)

ax2.tick_params(axis='both', which='major', **labelp)

plt.subplots_adjust(wspace=0, hspace=0)
plt.xlabel(r'$ \rm prior\ width$', **textp)
plt.xticks([0, 1, 2, 3, 4], [r'$0.1\sigma$', r'$0.5\sigma$', r'$1\sigma$', r'$2\sigma$', r'$10\sigma$'], **textp) 
plt.xlim([- 0.5, 4.5])
plt.savefig('./new_plots/'+'23s_gA-prior_width.pdf', transparent=True)
plt.tight_layout(pad=30, rect=plt_axes)

plt.show()
# %%

# %%

# %%
import h5py as h5 
import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  
import lsqfit as lsf
import os 
import hashlib

from fh_db.ff_fit.models import FF_fit as ff_fit_FF_fit
os.environ["DJANGO_ALLOW_ASYNC_UNSAFE"] = "true" # for jupyter to use database

# %%
file_path = os.getcwd() + '/FK_Fpi_data.h5' # put the data file in the same path as this file.
myfile = h5.File(file_path, 'r')

print([key for key in myfile]) # see what inside the data file

mpi_data = myfile['a09m310']['mpi'] # a09m310 as an example
mpi_list = [data for data in mpi_data]
mpi_array = np.array(mpi_list)

mpi = mpi_array.mean() # read the mass of pion from the data file
mN = 0.49 # this is the mass of proton, you can use best fit results of E0 as the mass of proton here
L = 32 #this is the size of the box

print(mpi)

# %%
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

errorp = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1}
errorb = {"markersize": 5, "linestyle": "none", "capsize": 3, "elinewidth": 1}
gridspec_tmin = {'height_ratios': [2, 2, 2, 2, 2, 1.5, 1.5], 'left': 0.15, 'right': 0.88, 'bottom': 0.1, 'top': 0.95, 'hspace':0}
fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
figsize  = (fig_width, fig_width * 1.2)
aspect=[0.15, 0.15, 0.845, 0.845]
textp = {"fontsize": 18}
labelp = {"labelsize": 18}

plt.rcParams.update({"text.usetex": True})

# %%
# from the database read the fit results and plot to check whether they are stable

tmin_list = [] # here check the stability when varying the tmin of 2pt fit as an example

E0_list = []
E0err_list = []
E1_list = []
E1err_list = []
E2_list = []
E2err_list = []
E3_list = []
E3err_list = []
E4_list = []
E4err_list = []

Q_list = []
log_list = []

# ho prior
for situation in ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=True, 
pt2_tmax=18, pt2_nstates=5,
pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=1,
sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5):
    dump_temp = bytes.fromhex(situation.data_and_results)
    data_and_results = gv.loads(dump_temp)
    
    tmin_list.append(situation.pt2_tmin)
    
    E0_list.append(data_and_results['params']['E0'].mean)
    E0err_list.append(data_and_results['params']['E0'].sdev)
    
    E1_list.append((data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    E1err_list.append((data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    E2_list.append((data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    E2err_list.append((data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    E3_list.append((data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    E3err_list.append((data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)

    E4_list.append((data_and_results['params']['dE4'] + data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    E4err_list.append((data_and_results['params']['dE4'] + data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    Q_list.append(situation.Q_value)
    log_list.append(situation.log_GBF)
    
    print('hologGBF: ' + str(situation.log_GBF))

print(tmin_list)
    
tmin_list = []    
    
swE0_list = []
swE0err_list = []
swE1_list = []
swE1err_list = []
swE2_list = []
swE2err_list = []
swE3_list = []
swE3err_list = []
swE4_list = []
swE4err_list = []

swQ_list = []
swlog_list = []
    
# sw prior    
for situation in ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='ce6f28b84ebcdf0de58d5fa15396687a', include_2pt=True, include_3pt=True, include_sum=True, 
pt2_tmax=18, pt2_nstates=5,
pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=1,
sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5):
    dump_temp = bytes.fromhex(situation.data_and_results)
    data_and_results = gv.loads(dump_temp)
    
    tmin_list.append(situation.pt2_tmin)
    
    swE0_list.append(data_and_results['params']['E0'].mean)
    swE0err_list.append(data_and_results['params']['E0'].sdev)
    
    swE1_list.append((data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    swE1err_list.append((data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    swE2_list.append((data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    swE2err_list.append((data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    swE3_list.append((data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    swE3err_list.append((data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)

    swE4_list.append((data_and_results['params']['dE4'] + data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    swE4err_list.append((data_and_results['params']['dE4'] + data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    swQ_list.append(situation.Q_value)
    swlog_list.append(situation.log_GBF)
    
    print('swlogGBF: ' + str(situation.log_GBF))

print(tmin_list)
    
tmin_list = []    
    
id1E0_list = []
id1E0err_list = []
id1E1_list = []
id1E1err_list = []
id1E2_list = []
id1E2err_list = []
id1E3_list = []
id1E3err_list = []
id1E4_list = []
id1E4err_list = []

id1Q_list = []
id1log_list = []
    
# id1
for situation in ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='6d93adeadf0937e2ab03d441887082ca', include_2pt=True, include_3pt=True, include_sum=True, 
pt2_tmax=18, pt2_nstates=5,
pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=1,
sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5):
    dump_temp = bytes.fromhex(situation.data_and_results)
    data_and_results = gv.loads(dump_temp)
    
    tmin_list.append(situation.pt2_tmin)
    
    id1E0_list.append(data_and_results['params']['E0'].mean)
    id1E0err_list.append(data_and_results['params']['E0'].sdev)
    
    id1E1_list.append((data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    id1E1err_list.append((data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    id1E2_list.append((data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    id1E2err_list.append((data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    id1E3_list.append((data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    id1E3err_list.append((data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)

    id1E4_list.append((data_and_results['params']['dE4'] + data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    id1E4err_list.append((data_and_results['params']['dE4'] + data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    id1Q_list.append(situation.Q_value)
    id1log_list.append(situation.log_GBF)
    
    print('id1logGBF: ' + str(situation.log_GBF))

print(tmin_list)

tmin_list = []    
    
id2E0_list = []
id2E0err_list = []
id2E1_list = []
id2E1err_list = []
id2E2_list = []
id2E2err_list = []
id2E3_list = []
id2E3err_list = []
id2E4_list = []
id2E4err_list = []

id2Q_list = []
id2log_list = []
    
# id2
for situation in ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8251565b53438fc6cf6fee9df7b98c09', include_2pt=True, include_3pt=True, include_sum=True, 
pt2_tmax=18, pt2_nstates=5,
pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=1,
sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5):
    dump_temp = bytes.fromhex(situation.data_and_results)
    data_and_results = gv.loads(dump_temp)
    
    tmin_list.append(situation.pt2_tmin)
    
    id2E0_list.append(data_and_results['params']['E0'].mean)
    id2E0err_list.append(data_and_results['params']['E0'].sdev)
    
    id2E1_list.append((data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    id2E1err_list.append((data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    id2E2_list.append((data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    id2E2err_list.append((data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    id2E3_list.append((data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    id2E3err_list.append((data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)

    id2E4_list.append((data_and_results['params']['dE4'] + data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).mean)
    id2E4err_list.append((data_and_results['params']['dE4'] + data_and_results['params']['dE3'] + data_and_results['params']['dE2'] + data_and_results['params']['dE1'] + data_and_results['params']['E0']).sdev)
    
    id2Q_list.append(situation.Q_value)
    
    id2log_list.append(situation.log_GBF)
    
    print('id2logGBF: ' + str(situation.log_GBF))

print(tmin_list)

##############################################################################    
# plot fit results of En
fig=plt.figure()
plt.rcParams['figure.figsize'] = figsize
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, 1, sharex = True, gridspec_kw=gridspec_tmin) 

ax1.set_ylabel(r"$a_{09}E_4$", **textp)
ax2.set_ylabel(r"$a_{09}E_3$", **textp)
ax3.set_ylabel(r"$a_{09}E_2$", **textp)
ax4.set_ylabel(r"$a_{09}E_1$", **textp)
ax5.set_ylabel(r"$a_{09}E_0$", **textp)

# square well
ax5.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE0_list), yerr=np.array(swE0err_list),marker='o', color=peach, **errorp)
ax4.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE1_list), yerr=np.array(swE1err_list),marker='^', color=peach, **errorp)
ax3.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE2_list), yerr=np.array(swE2err_list),marker='s', color=peach, **errorp)
ax2.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE3_list), yerr=np.array(swE3err_list),marker='P', color=peach, **errorp)
ax1.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE4_list), yerr=np.array(swE4err_list),marker='o', color=peach, label='SqW', **errorp)

# harmonic oscillator
ax5.errorbar(np.array(tmin_list), np.array(E0_list), yerr=np.array(E0err_list), marker='o', color=red, **errorp)
ax4.errorbar(np.array(tmin_list), np.array(E1_list), yerr=np.array(E1err_list), marker='^', color=red, **errorp)
ax3.errorbar(np.array(tmin_list), np.array(E2_list), yerr=np.array(E2err_list), marker='s', color=red, **errorp)
ax2.errorbar(np.array(tmin_list), np.array(E3_list), yerr=np.array(E3err_list), marker='P', color=red, **errorp)
ax1.errorbar(np.array(tmin_list), np.array(E4_list), yerr=np.array(E4err_list), marker='P', color=red, label='HO', **errorp)

# id 1
ax5.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E0_list), yerr=np.array(id1E0err_list),marker='o', color=green, **errorp)
ax4.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E1_list), yerr=np.array(id1E1err_list),marker='^', color=green, **errorp)
ax3.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E2_list), yerr=np.array(id1E2err_list),marker='s', color=green, **errorp)
ax2.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E3_list), yerr=np.array(id1E3err_list),marker='P', color=green, **errorp)
ax1.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E4_list), yerr=np.array(id1E4err_list),marker='P', color=green, label='$1/n$', **errorp)

# id 2
ax5.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E0_list), yerr=np.array(id2E0err_list),marker='o', color=blue, **errorp)
ax4.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E1_list), yerr=np.array(id2E1err_list),marker='^', color=blue, **errorp)
ax3.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E2_list), yerr=np.array(id2E2err_list),marker='s', color=blue, **errorp)
ax2.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E3_list), yerr=np.array(id2E3err_list),marker='P', color=blue, **errorp)
ax1.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E4_list), yerr=np.array(id2E4err_list),marker='P', color=blue, label='$1/n^2$', **errorp)

# best fit
ax5.errorbar(np.array(tmin_list[4]), np.array(E0_list[4]), yerr=np.array(E0err_list[4]), marker='o', mfc=red, color=red, **errorb)
ax4.errorbar(np.array(tmin_list[4]), np.array(E1_list[4]), yerr=np.array(E1err_list[4]), marker='^', mfc=red, color=red, **errorb)
ax3.errorbar(np.array(tmin_list[4]), np.array(E2_list[4]), yerr=np.array(E2err_list[4]), marker='s', mfc=red, color=red, **errorb)
ax2.errorbar(np.array(tmin_list[4]), np.array(E3_list[4]), yerr=np.array(E3err_list[4]), marker='P', mfc=red, color=red, **errorb)
ax1.errorbar(np.array(tmin_list[4]), np.array(E4_list[4]), yerr=np.array(E4err_list[4]), marker='P', mfc=red, color=red, **errorb)

ax5.fill_between(np.arange(2.5, 8.5, 1), (E0_list[4]+E0err_list[4])*np.ones([6]), (E0_list[4]-E0err_list[4])*np.ones([6]), color=red, alpha=0.2)
ax4.fill_between(np.arange(2.5, 8.5, 1), (E1_list[4]+E1err_list[4])*np.ones([6]), (E1_list[4]-E1err_list[4])*np.ones([6]), color=red, alpha=0.2)
ax3.fill_between(np.arange(2.5, 8.5, 1), (E2_list[4]+E2err_list[4])*np.ones([6]), (E2_list[4]-E2err_list[4])*np.ones([6]), color=red, alpha=0.2)
ax2.fill_between(np.arange(2.5, 8.5, 1), (E3_list[4]+E3err_list[4])*np.ones([6]), (E3_list[4]-E3err_list[4])*np.ones([6]), color=red, alpha=0.2)
ax1.fill_between(np.arange(2.5, 8.5, 1), (E4_list[4]+E4err_list[4])*np.ones([6]), (E4_list[4]-E4err_list[4])*np.ones([6]), color=red, alpha=0.2)

const = 0.197/0.09 # lattice to GeV

def set_yaxis(axi, axj, Ei_list, Eierr_list):   
    axi.set_ylim([Ei_list[0]-3*Eierr_list[0], Ei_list[0]+3*Eierr_list[0]])
    axi.set_yticks([ round(Ei_list[0]-2*Eierr_list[0], 3), round(Ei_list[0], 3), round(Ei_list[0]+2*Eierr_list[0], 3) ])
    axj.set_ylim([(Ei_list[0]-3*Eierr_list[0])*const, (Ei_list[0]+3*Eierr_list[0])*const])
    axj.set_yticks([ round((Ei_list[0]-2*Eierr_list[0])*const-0.003, 2), round(Ei_list[0]*const, 2), round((Ei_list[0]+2*Eierr_list[0])*const, 2) ])
    axi.tick_params(axis='both', which='major', **labelp)
    axj.tick_params(axis='both', which='major', **labelp)

ax8 = ax1.twinx()
ax9 = ax2.twinx()
ax10 = ax3.twinx()
ax11 = ax4.twinx()
ax12 = ax5.twinx()


set_yaxis(ax2, ax9, E3_list, E3err_list)
set_yaxis(ax3, ax10, E2_list, E2err_list)
set_yaxis(ax4, ax11, E1_list, E1err_list)
set_yaxis(ax5, ax12, E0_list, E0err_list)

ax1.legend(ncol=4, loc='upper center')

ax1.set_ylim([E4_list[0]-3*E4err_list[0], E4_list[0]+5*E4err_list[0]])
ax1.set_yticks([ round(E4_list[0]-2*E4err_list[0], 3), round(E4_list[0], 3), round(E4_list[0]+2*E4err_list[0], 3) ])
ax8.set_ylim([(E4_list[0]-3*E4err_list[0])*const, (E4_list[0]+5*E4err_list[0])*const])
ax8.set_yticks([ round((E4_list[0]-2*E4err_list[0])*const, 2), round(E4_list[0]*const, 2), round((E4_list[0]+2*E4err_list[0])*const, 2) ])
ax1.tick_params(axis='both', which='major', **labelp)
ax8.tick_params(axis='both', which='major', **labelp)

ax10.set_ylabel(r'$E / {\rm GeV}$', fontsize=16)


# ax6
ax6.set_ylabel(r"$Q$", **textp)
ax6.set_ylim([0, 1.2])
ax6.set_yticks([0, 0.75])
ax6.plot(np.arange(2.5, 8.5, 1), 0.1 * np.ones([6]), 'r--')

ax6.scatter(np.array(tmin_list), np.array(Q_list), marker='o', c='', edgecolors=red)
ax6.scatter(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swQ_list), marker='o', c='', edgecolors=peach)
ax6.scatter(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1Q_list), marker='o', c='', edgecolors=green)
ax6.scatter(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2Q_list), marker='o', c='', edgecolors=blue)

# best fit
ax6.scatter(np.array(tmin_list[4]), np.array(Q_list[4]), marker='o', c=red, edgecolors=red)

# ax7
ax7.set_ylabel(r"$w$", **textp)
ax7.set_ylim([0, 1.2])
ax7.set_yticks([0, 0.75])
ax7.plot(np.arange(2.5, 8.5, 1), 0.3 * np.ones([6]), 'r--')

log_list_3 = [ li[4] for li in [log_list, swlog_list, id1log_list, id2log_list]]
log_list_4 = [ li[3] for li in [log_list, swlog_list, id1log_list, id2log_list]]
log_list_5 = [ li[2] for li in [log_list, swlog_list, id1log_list, id2log_list]]
log_list_6 = [ li[1] for li in [log_list, swlog_list, id1log_list, id2log_list]]
log_list_7 = [ li[0] for li in [log_list, swlog_list, id1log_list, id2log_list]]

log_max_list = [max(log_list_7), max(log_list_6), max(log_list_5), max(log_list_4), max(log_list_3)]

w_list = [ np.exp(log_list[i] - log_max_list[i]) for i in range(5) ]
sww_list = [ np.exp(swlog_list[i] - log_max_list[i]) for i in range(5) ]
id1w_list = [ np.exp(id1log_list[i] - log_max_list[i]) for i in range(5) ]
id2w_list = [ np.exp(id2log_list[i] - log_max_list[i]) for i in range(5) ]

ax7.scatter(np.array(tmin_list), np.array(w_list), marker='o', c='', edgecolors=red)
ax7.scatter(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(sww_list), marker='o', c='', edgecolors=peach)
ax7.scatter(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1w_list), marker='o', c='', edgecolors=green)
ax7.scatter(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2w_list), marker='o', c='', edgecolors=blue)

# best fit
ax7.scatter(np.array(tmin_list[4]), np.array(w_list[4]), marker='o', c=red, edgecolors=red)


ax6.tick_params(axis='both', which='major', **labelp)
ax7.tick_params(axis='both', which='major', **labelp)

ax7.set_xlabel(r"$t_{\rm sep}^{\rm min}:C_2$", **textp)

plt.subplots_adjust(wspace=0, hspace=0)
plt.xlabel(r"$t_{\rm sep}^{\rm min}:C_2$", **textp)
plt.xlim([2.5, 7.5])
plt.tight_layout(pad=30, rect=aspect)

plt.savefig('./new_plots/spec_23s_2pttmin.pdf', transparent=True)

#################################################################################

##### dump all posteriors into a dict ######
###### note that tmin_list = [7, 6, 5, 4, 3], so append from i=4 to i=0
spectrum_posterior = dict()
spectrum_posterior['SqW'] = dict()
spectrum_posterior['SqW']['E0'] = [gv.gvar(swE0_list[4-i], swE0err_list[4-i]) for i in range(5)]
spectrum_posterior['SqW']['E1'] = [gv.gvar(swE1_list[4-i], swE1err_list[4-i]) for i in range(5)]
spectrum_posterior['SqW']['E2'] = [gv.gvar(swE2_list[4-i], swE2err_list[4-i]) for i in range(5)]
spectrum_posterior['SqW']['E3'] = [gv.gvar(swE3_list[4-i], swE3err_list[4-i]) for i in range(5)]
spectrum_posterior['SqW']['E4'] = [gv.gvar(swE4_list[4-i], swE4err_list[4-i]) for i in range(5)]
spectrum_posterior['SqW']['Q'] = [swQ_list[4-i] for i in range(5)]
spectrum_posterior['SqW']['w'] = [sww_list[4-i] for i in range(5)]

spectrum_posterior['HO'] = dict()
spectrum_posterior['HO']['E0'] = [gv.gvar(E0_list[4-i], E0err_list[4-i]) for i in range(5)]
spectrum_posterior['HO']['E1'] = [gv.gvar(E1_list[4-i], E1err_list[4-i]) for i in range(5)]
spectrum_posterior['HO']['E2'] = [gv.gvar(E2_list[4-i], E2err_list[4-i]) for i in range(5)]
spectrum_posterior['HO']['E3'] = [gv.gvar(E3_list[4-i], E3err_list[4-i]) for i in range(5)]
spectrum_posterior['HO']['E4'] = [gv.gvar(E4_list[4-i], E4err_list[4-i]) for i in range(5)]
spectrum_posterior['HO']['Q'] = [Q_list[4-i] for i in range(5)]
spectrum_posterior['HO']['w'] = [w_list[4-i] for i in range(5)]

spectrum_posterior['1/n'] = dict()
spectrum_posterior['1/n']['E0'] = [gv.gvar(id1E0_list[4-i], id1E0err_list[4-i]) for i in range(5)]
spectrum_posterior['1/n']['E1'] = [gv.gvar(id1E1_list[4-i], id1E1err_list[4-i]) for i in range(5)]
spectrum_posterior['1/n']['E2'] = [gv.gvar(id1E2_list[4-i], id1E2err_list[4-i]) for i in range(5)]
spectrum_posterior['1/n']['E3'] = [gv.gvar(id1E3_list[4-i], id1E3err_list[4-i]) for i in range(5)]
spectrum_posterior['1/n']['E4'] = [gv.gvar(id1E4_list[4-i], id1E4err_list[4-i]) for i in range(5)]
spectrum_posterior['1/n']['Q'] = [id1Q_list[4-i] for i in range(5)]
spectrum_posterior['1/n']['w'] = [id1w_list[4-i] for i in range(5)]

spectrum_posterior['1/n2'] = dict()
spectrum_posterior['1/n2']['E0'] = [gv.gvar(id2E0_list[4-i], id2E0err_list[4-i]) for i in range(5)]
spectrum_posterior['1/n2']['E1'] = [gv.gvar(id2E1_list[4-i], id2E1err_list[4-i]) for i in range(5)]
spectrum_posterior['1/n2']['E2'] = [gv.gvar(id2E2_list[4-i], id2E2err_list[4-i]) for i in range(5)]
spectrum_posterior['1/n2']['E3'] = [gv.gvar(id2E3_list[4-i], id2E3err_list[4-i]) for i in range(5)]
spectrum_posterior['1/n2']['E4'] = [gv.gvar(id2E4_list[4-i], id2E4err_list[4-i]) for i in range(5)]
spectrum_posterior['1/n2']['Q'] = [id2Q_list[4-i] for i in range(5)]
spectrum_posterior['1/n2']['w'] = [id2w_list[4-i] for i in range(5)]

print(spectrum_posterior)

gv.dump(spectrum_posterior, 'spectrum_posterior')

"""

import gvar as gv
dic = gv.load('spectrum_posterior')
# dic['HO']['E0'] = [tmin=3, tmin=4, tmin=5, tmin=6, tmin=7]

"""
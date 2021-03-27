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

# %% out put the momentum distribution
from itertools import combinations

def balls_in_boxes(n, m):
    """Generate combinations of n balls in m boxes.

    >>> list(balls_in_boxes(4, 2))
    [(0, 4), (1, 3), (2, 2), (3, 1), (4, 0)]
    >>> list(balls_in_boxes(3, 3))
    [(0, 0, 3), (0, 1, 2), (0, 2, 1), (0, 3, 0), (1, 0, 2), (1, 1, 1), (1, 2, 0), (2, 0, 1), (2, 1, 0), (3, 0, 0)]

    """
    for c in combinations(range(n + m - 1), m - 1):
        yield tuple(b - a - 1 for a, b in zip((-1,) + c, c + (n + m - 1,)))
        
        
def gen_momentum(npi, mom):
    def non_zero_position(array):
        return tuple(i for i, val in enumerate(array) if val > 0)
    
    for array in balls_in_boxes(mom, 3*npi):    # array:tuple
        if sum(array)%2==(npi)%2:
            non_zero_pos = non_zero_position(array)  # zero_pos:tuple
            for j in range(len(non_zero_pos)+1):
                for c in combinations(non_zero_pos, j):
                    tarray = np.array(array)
                    for i in c:
                        tarray[i] *= -1
                    yield tarray

def permuteUnique(nums): # Full Permutation
    # write your code here
    def search(nums,lists,index):
        if(len(nums)==len(lists)):
            result.append(lists)
            return
        for i in range(0,len(nums)):
            if(flag[i]==1):continue

            elif(i!=0 and nums[i]==nums[i-1] and flag[i-1]==0):
                continue

            else :
                flag[i]=1
                search(nums,lists+[nums[i]],index+1)
                flag[i]=0

    result=[]

    flag=[]
    for i in range(len(nums)):
        flag.append(0)
    search(sorted(nums),[],0)
    return result

def all_mom_array(npi, mom): # output all possible array, npi is num of pions, npi*3 is num of momentum composante, mom is total momentum number like (2, 1, 1) has mom = 4
    a = np.zeros([mom]) # a "0" represent 1 unit momentum
    b = np.ones([3*npi-1]) # "1" are used to divide zeros
    l = list(a)+list(b)
    
    all_list = permuteUnique(l)
    
    mom_list = [] 
    
    for listx in all_list:
        count = 0
        for i in range(len(listx)):
            if listx[i] == 1.0:
                mom_list.append(count)
                count = 0
            elif listx[i] == 0.0:
                count += 1
                
            if i == len(listx)-1: # last momentum
                mom_list.append(count) # every digit in mom_list represent momentum of a pion in one direction
                
    all_mom_array = np.array(mom_list).reshape([-1, 3*npi]) # each array in all_mom_array represent a kind of distribution of momentum
    
    return all_mom_array # all possiblility

# %% # calc total energy with momentum distribution
def tot_e(mom_array, mpi, mN):
    npi = int(len(mom_array) / 3)
    ep = np.zeros([npi])
    sx = sy = sz = 0
    for i in range(npi):
        s = mom_array[3*i+0]**2 + mom_array[3*i+1]**2 + mom_array[3*i+2]**2
        ep[i] = np.sqrt(mpi**2 + s*(2*np.pi/L)**2)
        
        sx += mom_array[3*i+0]
        sy += mom_array[3*i+1]
        sz += mom_array[3*i+2]
    
    sN = sx**2 + sy**2 + sz**2
    eN = np.sqrt(mN**2 + sN*(2*np.pi/L)**2)
    
    tot_e = eN
    for i in range(npi):
        tot_e += ep[i]
    
    return tot_e

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
labelsize = {'size': 18}

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
ax5.errorbar(np.array(tmin_list[0]), np.array(E0_list[0]), yerr=np.array(E0err_list[0]), marker='o', mfc=red, color=red, **errorb)
ax4.errorbar(np.array(tmin_list[0]), np.array(E1_list[0]), yerr=np.array(E1err_list[0]), marker='^', mfc=red, color=red, **errorb)
ax3.errorbar(np.array(tmin_list[0]), np.array(E2_list[0]), yerr=np.array(E2err_list[0]), marker='s', mfc=red, color=red, **errorb)
ax2.errorbar(np.array(tmin_list[0]), np.array(E3_list[0]), yerr=np.array(E3err_list[0]), marker='P', mfc=red, color=red, **errorb)
ax1.errorbar(np.array(tmin_list[0]), np.array(E4_list[0]), yerr=np.array(E4err_list[0]), marker='P', mfc=red, color=red, **errorb)

ax5.fill_between(np.arange(2.5, 8.5, 1), (E0_list[0]+E0err_list[0])*np.ones([6]), (E0_list[0]-E0err_list[0])*np.ones([6]), color=red, alpha=0.2)
ax4.fill_between(np.arange(2.5, 8.5, 1), (E1_list[0]+E1err_list[0])*np.ones([6]), (E1_list[0]-E1err_list[0])*np.ones([6]), color=red, alpha=0.2)
ax3.fill_between(np.arange(2.5, 8.5, 1), (E2_list[0]+E2err_list[0])*np.ones([6]), (E2_list[0]-E2err_list[0])*np.ones([6]), color=red, alpha=0.2)
ax2.fill_between(np.arange(2.5, 8.5, 1), (E3_list[0]+E3err_list[0])*np.ones([6]), (E3_list[0]-E3err_list[0])*np.ones([6]), color=red, alpha=0.2)
ax1.fill_between(np.arange(2.5, 8.5, 1), (E4_list[0]+E4err_list[0])*np.ones([6]), (E4_list[0]-E4err_list[0])*np.ones([6]), color=red, alpha=0.2)

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

#set_yaxis(ax1, ax8, E4_list, E4err_list)
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
ax6.scatter(np.array(tmin_list[0]), np.array(Q_list[0]), marker='o', c=red, edgecolors=red)

# ax7
ax7.set_ylabel(r"$w$", **textp)
ax7.set_ylim([0, 1.2])
ax7.set_yticks([0, 0.75])
ax7.plot(np.arange(2.5, 8.5, 1), 0.3 * np.ones([6]), 'r--')

w_list = [0.7591742989696858, 0.0017540544417473593, 1, 0.008397830121424019, 0.03593942605030466]
sww_list = [1, 1, 0.6273322515434925, 1, 1]
id1w_list = [4.4597223650273686e-08, 4.1640568963641826e-08, 1.0574320691685268e-07, 1.6622175012779632e-07, 2.0924988449570086e-07]
id2w_list = [3.21556125844409e-15, 2.69591162005402e-15, 3.8414543732742373e-16, 4.1399535039242913e-16, 1.157588012901151e-15]

ax7.scatter(np.array(tmin_list), np.array(w_list), marker='o', c='', edgecolors=red)
ax7.scatter(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(sww_list), marker='o', c='', edgecolors=peach)
ax7.scatter(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1w_list), marker='o', c='', edgecolors=green)
ax7.scatter(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2w_list), marker='o', c='', edgecolors=blue)


ax6.tick_params(axis='both', which='major', **labelp)
ax7.tick_params(axis='both', which='major', **labelp)

ax7.set_xlabel(r"$t_{\rm sep}^{\rm min}:C_2$", **textp)

plt.subplots_adjust(wspace=0, hspace=0)
plt.xlabel(r"$t_{\rm sep}^{\rm min}:C_2$", **textp)
plt.xlim([2.5, 7.5])
plt.tight_layout(pad=30, rect=aspect)

plt.savefig('./new_plots/spec_23s_2pttmin.pdf', transparent=True)

#################################################################################




# %%
# read the prior of En from the database
E_prior_list = []
Eerr_prior_list = []

for situation in ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=True, 
pt2_tmin=3, pt2_tmax=18, pt2_nstates=5,
pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=1,
sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5):
    dump_temp = bytes.fromhex(situation.data_and_results)
    data_and_results = gv.loads(dump_temp)
    
    E0 = data_and_results['prior']['E0']
    E1 = data_and_results['prior']['E0'] + np.exp(data_and_results['prior']['log(dE1)'])
    E2 = data_and_results['prior']['E0'] + np.exp(data_and_results['prior']['log(dE1)']) + np.exp(data_and_results['prior']['log(dE2)'])
    E3 = data_and_results['prior']['E0'] + np.exp(data_and_results['prior']['log(dE1)']) + np.exp(data_and_results['prior']['log(dE2)']) + np.exp(data_and_results['prior']['log(dE3)'])
    E4 = data_and_results['prior']['E0'] + np.exp(data_and_results['prior']['log(dE1)']) + np.exp(data_and_results['prior']['log(dE2)']) + np.exp(data_and_results['prior']['log(dE3)']) + np.exp(data_and_results['prior']['log(dE4)'])

print(E4.mean)
print(E4.sdev)
    
for En in [E0, E1, E2, E3, E4]:#, E5, E6]:
    E_prior_list.append(En.mean)
    Eerr_prior_list.append(En.sdev)
    
E_fit_list = []
Eerr_fit_list = []

for E_list in [E0_list, E1_list, E2_list, E3_list, E4_list]:#, E5_list, E6_list]:
    E_fit_list.append(E_list[0]) # here [0] means choosing tmin=3 situation for plotting
    
for Eerr_list in [E0err_list, E1err_list, E2err_list, E3err_list, E4err_list]:#, E5err_list, E6err_list]:
    Eerr_fit_list.append(Eerr_list[0]) # here [0] means choosing tmin=3 situation for plotting
    
    
plt.figure(figsize=(8, 10))
# fit results 
plt.fill_between(np.array([1.625, 1.875]), E_fit_list[3]+Eerr_fit_list[3], E_fit_list[3]-Eerr_fit_list[3], color='k', alpha=0.5, zorder=2)
plt.fill_between(np.array([1.125, 1.375]), E_fit_list[2]+Eerr_fit_list[2], E_fit_list[2]-Eerr_fit_list[2], color='k', alpha=0.5, zorder=2)
plt.fill_between(np.array([0.625, 0.875]), E_fit_list[1]+Eerr_fit_list[1], E_fit_list[1]-Eerr_fit_list[1], color='k', alpha=0.5, zorder=2)
plt.fill_between(np.array([0.125, 0.375]), E_fit_list[0]+Eerr_fit_list[0], E_fit_list[0]-Eerr_fit_list[0], color='k', alpha=0.5, label='fit results', zorder=2)
# prior
plt.fill_between(np.array([0.0, 0.5]), E_prior_list[0]+Eerr_prior_list[0], E_prior_list[0]-Eerr_prior_list[0], color=blue, alpha=0.6, label='prior', zorder=2)

plt.fill_between(np.array([0.5, 1.0]), E_prior_list[1]+Eerr_prior_list[1], E_prior_list[1]-Eerr_prior_list[1], color=blue, alpha=0.5, zorder=2)

plt.fill_between(np.array([1.0, 1.5]), E_prior_list[2]+Eerr_prior_list[2], E_prior_list[2]-Eerr_prior_list[2], color=blue, alpha=0.6, zorder=2)

plt.fill_between(np.array([1.5, 2.0]), E_prior_list[3]+Eerr_prior_list[3], E_prior_list[3]-Eerr_prior_list[3], color=blue, alpha=0.5, zorder=2)

print(E_fit_list)
print(Eerr_fit_list)
print(E_prior_list)
print(Eerr_prior_list)

print(E4_list[0])
print(E4err_list[0])

# theoretical prediction
plot_range = np.array([-0.5, 2.5])

for npi in range(1, 9):
    for mom in range(0, 6):
        e_array = np.array(tuple(map(
            lambda array: tot_e(array,mpi,mN),
            gen_momentum(npi, mom))))
        e_array_uni = np.unique(e_array)
        for e in e_array_uni:
            if e < 2:
                plt.plot(plot_range, e*np.ones([2]), color='k', linestyle='--', linewidth=0.3, zorder=0)
                

n_label = -1                    
for npi in range(1, 9):
    for mom in range(0, 2):
        if npi%2 == mom%2:
            all_possible_array = all_mom_array(npi, mom)
            for array in all_possible_array:
                e = tot_e(array, mpi, mN)
                if e < 2:
                    if n_label != npi:
                        plt.plot(np.linspace(-0.5, 2.5, npi+4), e*np.ones([npi+4]), color=color_list[npi-1], linestyle='--', marker='x', linewidth=1.5, zorder=1, label=str(npi)+'pions')
                        n_label = npi
                        print(e)

plt.xlim([-0.5, 3.2])
plt.ylim([0.4, 2])
plt.xticks([0.25, 0.75, 1.25, 1.75], ['E0', 'E1', 'E2', 'E3'], size=18)
plt.legend(loc='upper right')
#plt.title('a09_23s_E1_best_fit', size=18)
plt.savefig('./new_plots/'+'23s_spec_oh.pdf', transparent=True)
plt.show()


# %%
'''
fig=plt.figure()
plt.rcParams['figure.figsize'] = figsize
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, gridspec_kw=gridspec_tmin) 

ax1.set_ylabel(r"$E_n$", **textp)

# harmonic oscillator
ax1.errorbar(np.array(tmin_list), np.array(E0_list), yerr=np.array(E0err_list), marker='o', color=red, label='ho', **errorp)
ax1.errorbar(np.array(tmin_list), np.array(E1_list), yerr=np.array(E1err_list), marker='^', color=red, **errorp)
ax1.errorbar(np.array(tmin_list), np.array(E2_list), yerr=np.array(E2err_list), marker='s', color=red, **errorp)
ax1.errorbar(np.array(tmin_list), np.array(E3_list), yerr=np.array(E3err_list), marker='P', color=red, **errorp)

# square well
ax1.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE0_list), yerr=np.array(swE0err_list),marker='o', color=peach, label='sw', **errorp)
ax1.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE1_list), yerr=np.array(swE1err_list),marker='^', color=peach, **errorp)
ax1.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE2_list), yerr=np.array(swE2err_list),marker='s', color=peach, **errorp)
ax1.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE3_list), yerr=np.array(swE3err_list),marker='P', color=peach, **errorp)

# id 1
ax1.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E0_list), yerr=np.array(id1E0err_list),marker='o', color=green, label='id1', **errorp)
ax1.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E1_list), yerr=np.array(id1E1err_list),marker='^', color=green, **errorp)
ax1.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E2_list), yerr=np.array(id1E2err_list),marker='s', color=green, **errorp)
ax1.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E3_list), yerr=np.array(id1E3err_list),marker='P', color=green, **errorp)

# id 2
ax1.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E0_list), yerr=np.array(id2E0err_list),marker='o', color=blue, label='id2', **errorp)
ax1.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E1_list), yerr=np.array(id2E1err_list),marker='^', color=blue, **errorp)
ax1.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E2_list), yerr=np.array(id2E2err_list),marker='s', color=blue, **errorp)
ax1.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E3_list), yerr=np.array(id2E3err_list),marker='P', color=blue, **errorp)

# best fit
ax1.errorbar(np.array(tmin_list[0]), np.array(E0_list[0]), yerr=np.array(E0err_list[0]), marker='o', mfc=red, color=red, **errorb)
ax1.errorbar(np.array(tmin_list[0]), np.array(E1_list[0]), yerr=np.array(E1err_list[0]), marker='^', mfc=red, color=red, **errorb)
ax1.errorbar(np.array(tmin_list[0]), np.array(E2_list[0]), yerr=np.array(E2err_list[0]), marker='s', mfc=red, color=red, **errorb)
ax1.errorbar(np.array(tmin_list[0]), np.array(E3_list[0]), yerr=np.array(E3err_list[0]), marker='P', mfc=red, color=red, **errorb)

ax1.fill_between(np.arange(2.5, 8.5, 1), (E0_list[0]+E0err_list[0])*np.ones([6]), (E0_list[0]-E0err_list[0])*np.ones([6]), color=red, alpha=0.2)
ax1.fill_between(np.arange(2.5, 8.5, 1), (E1_list[0]+E1err_list[0])*np.ones([6]), (E1_list[0]-E1err_list[0])*np.ones([6]), color=red, alpha=0.2)
ax1.fill_between(np.arange(2.5, 8.5, 1), (E2_list[0]+E2err_list[0])*np.ones([6]), (E2_list[0]-E2err_list[0])*np.ones([6]), color=red, alpha=0.2)
ax1.fill_between(np.arange(2.5, 8.5, 1), (E3_list[0]+E3err_list[0])*np.ones([6]), (E3_list[0]-E3err_list[0])*np.ones([6]), color=red, alpha=0.2)


# ax2
ax2.set_ylabel(r"$Q$", **textp)
ax2.set_ylim([0, 1.1])
ax2.plot(np.arange(2.5, 8.5, 1), 0.1 * np.ones([6]), 'r--')

ax2.scatter(np.array(tmin_list), np.array(Q_list), marker='o', c='', edgecolors=red)
ax2.scatter(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swQ_list), marker='o', c='', edgecolors=peach)
ax2.scatter(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1Q_list), marker='o', c='', edgecolors=green)
ax2.scatter(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2Q_list), marker='o', c='', edgecolors=blue)

# best fit
ax2.scatter(np.array(tmin_list[0]), np.array(Q_list[0]), marker='o', c=red, edgecolors=red)

# ax3
ax3.set_ylabel(r"$w$", **textp)
ax3.set_ylim([0, 1.1])

ax3.plot(np.arange(2.5, 8.5, 1), 0.3 * np.ones([6]), 'r--')

w_list = [0.7591742989696858, 0.0017540544417473593, 1, 0.008397830121424019, 0.03593942605030466]
sww_list = [1, 1, 0.6273322515434925, 1, 1]
id1w_list = [4.4597223650273686e-08, 4.1640568963641826e-08, 1.0574320691685268e-07, 1.6622175012779632e-07, 2.0924988449570086e-07]
id2w_list = [3.21556125844409e-15, 2.69591162005402e-15, 3.8414543732742373e-16, 4.1399535039242913e-16, 1.157588012901151e-15]

ax3.scatter(np.array(tmin_list), np.array(w_list), marker='o', c='', edgecolors=red)
ax3.scatter(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(sww_list), marker='o', c='', edgecolors=peach)
ax3.scatter(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1w_list), marker='o', c='', edgecolors=green)
ax3.scatter(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2w_list), marker='o', c='', edgecolors=blue)


plt.subplots_adjust(wspace=0, hspace=0)
plt.xlabel(r"$C_{\textrm{2pt}}\ t_{\textrm{min}}$", **textp)
plt.xlim([2.5, 7.5])
ax1.tick_params(axis='both', which='major', **labelp)
ax2.tick_params(axis='both', which='major', **labelp)
ax3.tick_params(axis='both', which='major', **labelp)

plt.tight_layout(pad=30, rect=aspect)

plt.savefig('./new_plots/spec_23s_2pttmin.pdf', transparent=True)
'''

'''
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(ncols=3, nrows=2, gridspec_kw=gridspec_tmin) 

ax1.set_ylabel(r"$E_0$", **textp)
ax1.yaxis.set_label_coords(0, 1)
ax2.set_ylabel(r"$E_1$", **textp)
ax2.yaxis.set_label_coords(0, 1)
ax3.set_ylabel(r"$E_2$", **textp)
ax3.yaxis.set_label_coords(0, 1)
ax4.set_ylabel(r"$E_3$", **textp)
ax4.yaxis.set_label_coords(0, 1)

# harmonic oscillator
ax1.errorbar(np.array(tmin_list), np.array(E0_list), yerr=np.array(E0err_list), marker='o', color=red, label='ho', **errorp)
ax2.errorbar(np.array(tmin_list), np.array(E1_list), yerr=np.array(E1err_list), marker='^', color=red, **errorp)
ax3.errorbar(np.array(tmin_list), np.array(E2_list), yerr=np.array(E2err_list), marker='s', color=red, **errorp)
ax4.errorbar(np.array(tmin_list), np.array(E3_list), yerr=np.array(E3err_list), marker='P', color=red, **errorp)

# square well
ax1.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE0_list), yerr=np.array(swE0err_list),marker='o', color=peach, label='sw', **errorp)
ax2.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE1_list), yerr=np.array(swE1err_list),marker='^', color=peach, **errorp)
ax3.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE2_list), yerr=np.array(swE2err_list),marker='s', color=peach, **errorp)
ax4.errorbar(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swE3_list), yerr=np.array(swE3err_list),marker='P', color=peach, **errorp)

# id 1
ax1.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E0_list), yerr=np.array(id1E0err_list),marker='o', color=green, label='id1', **errorp)
ax2.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E1_list), yerr=np.array(id1E1err_list),marker='^', color=green, **errorp)
ax3.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E2_list), yerr=np.array(id1E2err_list),marker='s', color=green, **errorp)
ax4.errorbar(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1E3_list), yerr=np.array(id1E3err_list),marker='P', color=green, **errorp)

# id 2
ax1.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E0_list), yerr=np.array(id2E0err_list),marker='o', color=blue, label='id2', **errorp)
ax2.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E1_list), yerr=np.array(id2E1err_list),marker='^', color=blue, **errorp)
ax3.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E2_list), yerr=np.array(id2E2err_list),marker='s', color=blue, **errorp)
ax4.errorbar(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2E3_list), yerr=np.array(id2E3err_list),marker='P', color=blue, **errorp)

# best fit
ax1.errorbar(np.array(tmin_list[0]), np.array(E0_list[0]), yerr=np.array(E0err_list[0]), marker='o', mfc=red, color=red, **errorb)
ax2.errorbar(np.array(tmin_list[0]), np.array(E1_list[0]), yerr=np.array(E1err_list[0]), marker='^', mfc=red, color=red, **errorb)
ax3.errorbar(np.array(tmin_list[0]), np.array(E2_list[0]), yerr=np.array(E2err_list[0]), marker='s', mfc=red, color=red, **errorb)
ax4.errorbar(np.array(tmin_list[0]), np.array(E3_list[0]), yerr=np.array(E3err_list[0]), marker='P', mfc=red, color=red, **errorb)

ax1.fill_between(np.arange(2.5, 8.5, 1), (E0_list[0]+E0err_list[0])*np.ones([6]), (E0_list[0]-E0err_list[0])*np.ones([6]), color=red, alpha=0.2)
ax2.fill_between(np.arange(2.5, 8.5, 1), (E1_list[0]+E1err_list[0])*np.ones([6]), (E1_list[0]-E1err_list[0])*np.ones([6]), color=red, alpha=0.2)
ax3.fill_between(np.arange(2.5, 8.5, 1), (E2_list[0]+E2err_list[0])*np.ones([6]), (E2_list[0]-E2err_list[0])*np.ones([6]), color=red, alpha=0.2)
ax4.fill_between(np.arange(2.5, 8.5, 1), (E3_list[0]+E3err_list[0])*np.ones([6]), (E3_list[0]-E3err_list[0])*np.ones([6]), color=red, alpha=0.2)


# ax5
ax5.set_ylabel(r"$Q$", **textp)
ax5.yaxis.set_label_coords(0, 1)
ax5.set_ylim([0.01, 1.2])
ax5.plot(np.arange(2.5, 8.5, 1), 0.1 * np.ones([6]), 'r--')

ax5.scatter(np.array(tmin_list), np.array(Q_list), marker='o', c='', edgecolors=red)
ax5.scatter(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(swQ_list), marker='o', c='', edgecolors=peach)
ax5.scatter(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1Q_list), marker='o', c='', edgecolors=green)
ax5.scatter(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2Q_list), marker='o', c='', edgecolors=blue)

# best fit
ax5.scatter(np.array(tmin_list[0]), np.array(Q_list[0]), marker='o', c=red, edgecolors=red)

# ax6
ax6.set_ylabel(r"$w$", **textp)
ax6.yaxis.set_label_coords(0, 1)
ax6.set_ylim([0.01, 1.2])

ax6.plot(np.arange(2.5, 8.5, 1), 0.3 * np.ones([6]), 'r--')

w_list = [0.7591742989696858, 0.0017540544417473593, 1, 0.008397830121424019, 0.03593942605030466]
sww_list = [1, 1, 0.6273322515434925, 1, 1]
id1w_list = [4.4597223650273686e-08, 4.1640568963641826e-08, 1.0574320691685268e-07, 1.6622175012779632e-07, 2.0924988449570086e-07]
id2w_list = [3.21556125844409e-15, 2.69591162005402e-15, 3.8414543732742373e-16, 4.1399535039242913e-16, 1.157588012901151e-15]

ax6.scatter(np.array(tmin_list), np.array(w_list), marker='o', c='', edgecolors=red)
ax6.scatter(np.array(tmin_list)-0.1 * np.ones(len(tmin_list)), np.array(sww_list), marker='o', c='', edgecolors=peach)
ax6.scatter(np.array(tmin_list)+0.1 * np.ones(len(tmin_list)), np.array(id1w_list), marker='o', c='', edgecolors=green)
ax6.scatter(np.array(tmin_list)+0.2 * np.ones(len(tmin_list)), np.array(id2w_list), marker='o', c='', edgecolors=blue)

const = 0.197/0.09 # lattice to GeV
#ax1.set_xlabel(r"$C_{\textrm{2pt}}\ t_{\textrm{min}}$", **textp)
ax1.set_xlim([2.5, 7.5])
ax1.set_ylim([0.483, 0.4949])
ax7 = ax1.twinx()
ax7.set_ylim([0.483*const, 0.4949*const])
ax7.set_ylabel(r"$GeV$", fontsize=14)
ax7.yaxis.set_label_coords(1.05, 1)
#ax2.set_xlabel(r"$C_{\textrm{2pt}}\ t_{\textrm{min}}$", **textp)
ax2.set_xlim([2.5, 7.5])
ax2.set_ylim([0.63, 0.88])
ax8 = ax2.twinx()
ax8.set_ylim([0.63*const, 0.88*const])
ax8.set_ylabel(r"$GeV$", fontsize=14)
ax8.yaxis.set_label_coords(1.05, 1)
#ax3.set_xlabel(r"$C_{\textrm{2pt}}\ t_{\textrm{min}}$", **textp)
ax3.set_xlim([2.5, 7.5])
ax3.set_ylim([0.85, 1.6])
ax9 = ax3.twinx()
ax9.set_ylim([0.85*const, 1.6*const])
ax9.set_ylabel(r"$GeV$", fontsize=14)
ax9.yaxis.set_label_coords(1.05, 1)
#ax4.set_xlabel(r"$C_{\textrm{2pt}}\ t_{\textrm{min}}$", **textp)
ax4.set_xlim([2.5, 7.5])
ax4.set_ylim([1.3, 1.95])
ax10 = ax4.twinx()
ax10.set_ylim([1.3*const, 1.95*const])
ax10.set_ylabel(r"$GeV$", fontsize=14)
ax10.yaxis.set_label_coords(1.05, 1)

ax5.set_xlabel(r"$t_{\rm sep}^{\rm min}:C_2$", **textp)
ax5.set_xlim([2.5, 7.5])
#ax6.set_xlabel(r"$C_{\textrm{2pt}}\ t_{\textrm{min}}$", **textp)
ax6.set_xlim([2.5, 7.5])

ax1.tick_params(axis='both', which='major', **labelp)
ax2.tick_params(axis='both', which='major', **labelp)
ax3.tick_params(axis='both', which='major', **labelp)
ax4.tick_params(axis='both', which='major', **labelp)
ax5.tick_params(axis='both', which='major', **labelp)
ax6.tick_params(axis='both', which='major', **labelp)
ax7.tick_params(axis='both', which='major', **labelp)
ax8.tick_params(axis='both', which='major', **labelp)
ax9.tick_params(axis='both', which='major', **labelp)
ax10.tick_params(axis='both', which='major', **labelp)

plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout(pad=30, rect=aspect)

plt.savefig('./new_plots/spec_23s_2pttmin.pdf', transparent=True)
'''
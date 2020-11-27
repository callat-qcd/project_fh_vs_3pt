# %%
import numpy as np 
import gvar as gv  
import hashlib
import os

from fh_db.ff_fit.models import FF_fit as ff_fit_FF_fit
os.environ["DJANGO_ALLOW_ASYNC_UNSAFE"] = "true" # for jupyter to use database

# %%
from module.prepare_data import Prepare_data
from module.fit import Fit

# %%
def fit_and_save(data_file_name, fit_type, save, include_2pt, include_3pt, include_sum, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut, id_num=None, pt2_range=[None, None], pt3_A3_range=[None, None], pt3_V4_range=[None, None], sum_A3_range=[None, None], sum_V4_range=[None, None], pt3_tau_dict=None):
    if fit_type == "scattered":

        pt2_t = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])

        pt3_A3_t = [3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14]
        pt3_A3_tau = [1, 1, 2, 1, 2, 1, 2, 3, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 7]

        pt3_V4_t = pt3_A3_t
        pt3_V4_tau = pt3_A3_tau
        
        pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
        pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]

        sum_A3 = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
        sum_V4 = sum_A3 

    elif fit_type == "continuous":
        pt2_t = np.arange(pt2_range[0], pt2_range[1])

        pt3_A3_t = []
        pt3_A3_tau = []

        for t in range(pt3_A3_range[0], pt3_A3_range[1]):
            for tau in pt3_tau_dict['A3_tsep'+str(t)]:
                pt3_A3_t.append(t)
                pt3_A3_tau.append(tau)

        pt3_V4_t = []
        pt3_V4_tau = []

        for t in range(pt3_V4_range[0], pt3_V4_range[1]):
            for tau in pt3_tau_dict['V4_tsep'+str(t)]:
                pt3_V4_t.append(t)
                pt3_V4_tau.append(tau)

        pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
        pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]

        sum_A3 = np.arange(sum_A3_range[0], sum_A3_range[1])
        sum_V4 = np.arange(sum_V4_range[0], sum_V4_range[1])

        print(pt2_t)
        print(pt3_A3)
        print(pt3_V4)
        print(sum_A3)
        print(sum_V4)

    data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)

    if use_p0 == False:
        best_p0 = 0

    best_p0 = {'E0': 0.49007432827585923, 'log(dE1)': -1.2182574657830274,  'z0': 0.0003202622719326246, 'z1': 0.00033086411929253656, 'z0_ps': 0.003, 'z1_ps': 1.4290420366321432e-21,'log(dE2)': -1.0203037715038503, 'z2': -0.0003420981842067054, 'z2_ps': 4.77757294807751e-19, 'log(dE3)': -0.6763116611503818, 'z3': 0.0006114436301814257, 'z3_ps': 6.18600096063414e-20, 'log(dE4)': 0.1096154276143707, 'z4': 0.00030912500545967415, 'z4_ps': -2.1747716630250984e-21, 'log(dE5)': -1.25, 'z5': -1.5592681273181436e-22, 'z5_ps': 1.3797601616142584e-19, 
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

    fitter = Fit(file_name, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

    fit_result, hexcode, dr_hex = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0, save)

    best_p0 = {key: fit_result.p[key].mean for key in fit_result.p}   

    print(fit_result.format(100)) 

    if save == True:
        tau_dict = gv.dumps(pt3_tau_dict)

        tau_dict_hex = tau_dict.hex()

        ff_fit_ff_fit, created = ff_fit_FF_fit.objects.get_or_create(
        data_file_name=data_file_name, # Store the name of data file
        include_3pt=include_3pt, # Whether include 3pt part in fitting
        include_sum=include_sum, # Whether include sum part in fitting
        data_and_results=dr_hex, # Use gv.dump to package data, prama, chi2 and dof, and turn into…
        prior_hexcode=hexcode, # put all priors into a string with fixed order and turn it…
        pt2_nstates=pt2_nstates, # Number of states in 2pt fit
        pt3_nstates=pt3_nstates, # Number of states in 3pt fit
        sum_nstates=sum_nstates, # Number of states in sum tau fit
        sum_tau_cut=sum_tau_cut, # Start of sum in sum tau fit
        E0=fit_result.p['E0'].mean, # fit result of E0
        E0_err=fit_result.p['E0'].sdev, # fit result of E0 err
        z0=fit_result.p['z0'].mean, # fit result of z0
        z0_err=fit_result.p['z0'].sdev, # fit result of z0 err
        A300=fit_result.p['A3_00'].mean, # fit result of A300
        A300_err=fit_result.p['A3_00'].sdev, # fit result of A300 err
        V400=fit_result.p['V4_00'].mean, # fit result of V400
        V400_err=fit_result.p['V4_00'].sdev, # fit result of V400 err
        Q_value=fit_result.Q, # Q value of fitting
        log_GBF=fit_result.logGBF, # logGBF of fitting
        fit_type=fit_type, # "scattered" or "continuous" 
        include_2pt=include_2pt, # (Optional) Whether include 2pt part in fitting
        pt2_tmin=pt2_range[0], # (Optional) Minimum of t in 2pt fit
        pt2_tmax=pt2_range[1], # (Optional) Maximum of t in 2pt fit
        pt3_A3_tsep_min=pt3_A3_range[0], # (Optional) Minimum of tsep in 3pt fit
        pt3_A3_tsep_max=pt3_A3_range[1], # (Optional) Maximum of tsep in 3pt fit
        pt3_V4_tsep_min=pt3_V4_range[0], # (Optional) Minimum of tsep in 3pt fit
        pt3_V4_tsep_max=pt3_V4_range[1], # (Optional) Maximum of tsep in 3pt fit
        pt3_tau_dict=tau_dict_hex, # (Optional) Use gv.dump to package tau dict of 3pt, and turn into hexcode
        sum_A3_tsep_min=sum_A3_range[0], # (Optional) Minimum of tsep in sum tau fit
        sum_A3_tsep_max=sum_A3_range[1], # (Optional) Maximum of tsep in sum tau fit
        sum_V4_tsep_min=sum_V4_range[0], # (Optional) Minimum of tsep in sum tau fit
        sum_V4_tsep_max=sum_V4_range[1], # (Optional) Maximum of tsep in sum tau fit
        id_num=id_num, # (Optional) used to divide and select
        )



# %%

file_name = 'a09m310_e_gA_srcs0-15.h5'
file_path = './' + file_name # just put the data file in the same path as this code.

pt2_data_range = [0, 96] # range in data file
pt3_data_range = [2, 15] # range in data file

prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

data_avg_dict = prepare_data.read_data_with_average()

use_p0 = False
include_2pt = True
include_3pt = True
include_sum = True

sum_tau_cut = 1
id_num = 1

fit_type = "continuous" 
save = False

##################################################################################################################
#########################################################

if fit_type == "scattered":

    pt2_nstates = 5
    pt3_nstates = 5
    sum_nstates = 5

    fit_and_save(file_name, fit_type, save, include_2pt, include_3pt, include_sum, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut, id_num)


############################################################

elif fit_type == "continuous":
    pt3_tau_dict = dict()
    for t in range(2, 11):
        pt3_tau_dict['A3_tsep'+str(t)] = np.arange(1, int(t/2)+1)
        pt3_tau_dict['V4_tsep'+str(t)] = np.arange(1, int(t/2)+1)

    for t in range(11, 15):
        pt3_tau_dict['A3_tsep'+str(t)] = np.arange(2, int(t/2)+1)
        pt3_tau_dict['V4_tsep'+str(t)] = np.arange(2, int(t/2)+1)

    pt2_range_list = [[tmin, tmax] for tmin in range(4, 5) for tmax in range(18, 19)]
    pt2_nstates_list = [nstates for nstates in range(5, 6)]

    pt3_A3_range_list = [[tmin, tmax] for tmin in range(3, 4) for tmax in range(15, 16)]
    pt3_V4_range_list = [[tmin, tmax] for tmin in range(3, 4) for tmax in range(15, 16)]
    #pt3_nstates_list = [nstates for nstates in range(5, 6)] # pt2_nstates = pt3_nstates

    sum_A3_range_list = [[tmin, tmax] for tmin in range(3, 4) for tmax in range(14, 15)]
    sum_V4_range_list = [[tmin, tmax] for tmin in range(3, 4) for tmax in range(14, 15)]
    sum_nstates_list = [nstates for nstates in range(5, 6)]

    times = 0

    for pt2_range in pt2_range_list:
        for pt2_nstates in pt2_nstates_list:
            for pt3_A3_range in pt3_A3_range_list:
                for pt3_V4_range in pt3_V4_range_list:
                    for sum_A3_range in sum_A3_range_list:
                        for sum_V4_range in sum_V4_range_list:
                            for sum_nstates in sum_nstates_list:
                                print('#################################################')
                                times += 1
                                pt3_nstates = pt2_nstates ##
                                fit_and_save(file_name, fit_type, save, include_2pt, include_3pt, include_sum, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut, id_num, pt2_range=pt2_range, pt3_A3_range=pt3_A3_range, pt3_V4_range=pt3_V4_range, sum_A3_range=sum_A3_range, sum_V4_range=sum_V4_range, pt3_tau_dict=pt3_tau_dict)
    print('total '+str(times)+' fits')
#########################################################


# %%

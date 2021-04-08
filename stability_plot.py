# %%
import os
import matplotlib.pyplot as plt

from fh_db.ff_fit.models import FF_fit as ff_fit_FF_fit
from module.plot import tmin_tmax_combined_plot
from module.plot import tmax_scattered_plot
from module.plot import late_23_tau_inc_plot
from module.plot import tmin_plot
from module.plot import tmax_plot
from module.plot import tmin_div_plot
from module.plot import tau_cut_plot
from module.plot import late_23_tmin_plot
from module.plot import prior_width_plot
from module.plot import best_fit_comparison_plot
from best_fits import late_23_fit


os.environ["DJANGO_ALLOW_ASYNC_UNSAFE"] = "true" # for jupyter to use database
plt.rcParams.update({"text.usetex": True})
file_path = './a09m310_e_gA_srcs0-15.h5'

c2pt_tmin = r"$t_{\rm sep}^{\rm min}:C_2$"
c2pt_tmax = r"$t_{\rm sep}^{\rm max}:C_2$"
c3pt_tmin = r"$t_{\rm sep}^{\rm min}:C_3$"
csum_tmin = r"$t_{\rm sep}^{\rm min}:\rm FH$"

# %%
def tmin_tmax_combined():
    #####################################
    ####### 23s tmin tmax combined ######
    #####################################
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

    tmin_tmax_combined_plot(x, value)

def both_even_tmax():
    #####################################
    ######## 23s 3pt/sum even tmax ######
    #####################################
    t_list = [8, 10, 12, 14]
    best_n = 5

    tmax_name = '3pt' # varying tmax is not so interesting as tmin 
    save_name = 'both_even'

    situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', fit_type='scattered', include_2pt=True, include_3pt=True, include_sum=True,
    pt2_nstates=5,
    id_num=2,
    sum_tau_cut=1, sum_nstates=5)

    tmax_scattered_plot(best_n, tmax_name, situation_list, save_name)

def both_odd_tmax():
    #####################################
    ######## 23s 3pt/sum odd tmax #######
    #####################################
    t_list = [7, 9, 11, 13]
    best_n = 5

    tmax_name = '3pt' # varying tmax is not so interesting as tmin 
    save_name = 'both_odd'

    situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', fit_type='scattered', include_2pt=True, include_3pt=True, include_sum=True,
    pt2_nstates=5,
    id_num=3,
    sum_tau_cut=1, sum_nstates=5)

    tmax_scattered_plot(best_n, tmax_name, situation_list, save_name)

def late_tsep_23_tau_inc():
    A3, A3_err, Q = late_23_fit(file_path)

    late_23_tau_inc_plot(A3, A3_err, Q)

def pt2_tmin_23s():
    ####################################
    ########## 23s 2pt tmin #############
    #####################################
    n_range=[4, 8]
    t_range=[3, 8]
    best_n=5
    best_t=3

    nstate_name='2pt'
    tmin_name='2pt'

    fit_name='23s'
    xlabel=c2pt_tmin

    situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=True, 
    pt2_tmax=18,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=1,
    sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5)

    tmin_plot(n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, fit_name, xlabel)

def pt3_tmin_23s():
    ####################################
    ########## 23s 3pt tmin #############
    #####################################
    n_range=[4, 8]
    t_range=[3, 8]
    best_n=5
    best_t=3

    nstate_name='3pt'
    tmin_name='3pt_gA' # here because gA and gV performs differently in 3pt, so you may change paras of gA/ gV separately

    fit_name='23s'
    xlabel=c3pt_tmin

    situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=True, 
    pt2_tmin=3, pt2_tmax=18,
    pt3_A3_tsep_max=15, id_num=1,
    sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5)

    tmin_plot(n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, fit_name, xlabel)

def sum_tmin_23s():
    ####################################
    ########## 23s sum tmin #############
    #####################################
    n_range=[4, 8]
    t_range=[3, 8]
    best_n=5
    best_t=3

    nstate_name='sum'
    tmin_name='sum_gA'

    fit_name='23s'
    xlabel=csum_tmin

    situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=True, 
    pt2_tmin=3, pt2_tmax=18, pt2_nstates=5,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=1,
    sum_A3_tsep_max=14, sum_tau_cut=1)

    tmin_plot(n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, fit_name, xlabel)

def tau_cut_23s():
    ####################################
    ########## 23s tau cut #############
    #####################################
    # taucut-1, taucut, taucut+1, taucut+2, taucut+3

    A3 = [1.326, 1.253, 1.258, 1.253, 1.253]
    A3_err = [0.012, 0.019, 0.022, 0.021, 0.022]
    Q = [0.076, 0.88, 0.61, 0.45, 0.21]

    n = 5

    fit_name = '23s'

    tau_cut_plot(A3, A3_err, Q, n, fit_name)

def pt2_tmin_23():
    ####################################
    ########## 23 2pt tmin #############
    #####################################
    n_range=[4, 8]
    t_range=[3, 8]
    best_n=5
    best_t=3

    nstate_name='2pt'
    tmin_name='2pt'

    fit_name='23'
    xlabel=c2pt_tmin

    situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=False, 
    pt2_tmax=18,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=5,
    sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5)

    tmin_plot(n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, fit_name, xlabel)

def pt3_tmin_23():
    ####################################
    ########## 23 3pt tmin #############
    #####################################
    n_range=[4, 8]
    t_range=[3, 8]
    best_n=5
    best_t=3

    nstate_name='3pt'
    tmin_name='3pt_gA'

    fit_name='23'
    xlabel=c3pt_tmin

    situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=True, include_sum=False, 
    pt2_tmin=3, pt2_tmax=18,
    pt3_A3_tsep_max=15, id_num=4,
    sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5)

    tmin_plot(n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, fit_name, xlabel)

def tau_cut_23():
    ####################################
    ########## 23 tau cut #############
    #####################################
    # taucut, taucut+1, taucut+2, taucut+3, taucut+4

    A3 = [1.260, 1.259, 1.254, 1.251, 1.265]
    A3_err = [0.023, 0.022, 0.023, 0.025, 0.030]
    Q = [0.74, 0.55, 0.7, 0.34, 0.27]

    n = 5

    fit_name = '23'

    tau_cut_plot(A3, A3_err, Q, n, fit_name)

def pt2_tmax_2s():
    t_range = [7, 19] 
    best_n = 2
    best_t = 18
    tmax_name = '2pt'

    fit_name = '2s'
    xlabel = c2pt_tmax

    situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=False, include_sum=True, 
    pt2_tmin=5, pt2_nstates=2,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=4,
    sum_A3_tsep_min=5, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=2)

    tmax_plot(t_range, best_n, best_t, tmax_name, situation_list, fit_name, xlabel)

def pt2_tmin_2s():
    n_range=[1, 4]
    t_range=[[9, 15], [4, 9], [2, 8]]
    best_n=2
    best_t=[12, 5, 5]
    tmin_name='2pt'

    fit_name='2s'
    xlabel=c2pt_tmin

    situation_1 = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='a99c26dd6000f43a5cce4158bd42c3f7', include_2pt=True, include_3pt=False, include_sum=True, 
    pt2_tmax=18, pt2_nstates=1,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=4,
    sum_A3_tsep_min=11, sum_V4_tsep_min=7, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=1)

    situation_2 = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=False, include_sum=True, 
    pt2_tmax=18, pt2_nstates=2,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=5,
    sum_A3_tsep_min=5, sum_V4_tsep_min=5, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=2)

    situation_3 = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=False, include_sum=True, 
    pt2_tmax=18, pt2_nstates=3,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=4,
    sum_A3_tsep_min=4, sum_V4_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=3)

    situation_list = [situation_1, situation_2, situation_3]

    tmin_div_plot(n_range, t_range, best_n, best_t, tmin_name, situation_list, fit_name, xlabel)

def sum_tmin_2s():
    n_range=[1, 4]
    t_range=[[10, 13], [4, 9], [2, 8]]
    best_n=2
    best_t=[11, 5, 3]
    tmin_name='sum_gA'

    fit_name='2s'
    xlabel=csum_tmin

    situation_1 = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='a99c26dd6000f43a5cce4158bd42c3f7', include_2pt=True, include_3pt=False, include_sum=True, 
    pt2_tmin=12, pt2_tmax=18, pt2_nstates=1,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=6,
    sum_V4_tsep_min=7, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=1)

    situation_2 = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=False, include_sum=True, 
    pt2_tmin=5, pt2_tmax=18, pt2_nstates=2,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=6,
    sum_V4_tsep_min=5, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=2)

    situation_3 = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=False, include_sum=True, 
    pt2_tmin=5, pt2_tmax=18, pt2_nstates=3,
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=6,
    sum_V4_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=3)

    situation_list = [situation_1, situation_2, situation_3]

    tmin_div_plot(n_range, t_range, best_n, best_t, tmin_name, situation_list, fit_name, xlabel)

def ga_summary():
    # 23s, 2s, 23
    A3 = [1.253440575215748, 1.2486801489322228, 1.2602800946082546]                 
    A3_err = [0.019302130858843348, 0.013323406775549013, 0.022831409526096327]
    Q = [0.88, 0.43, 0.59]

    best_fit_comparison_plot(A3, A3_err, Q)

def late_tsep_23_E0():
    n_range=[1, 5]
    t_range=[2, 17]
    best_n=2
    best_t=7

    E0_ylim = [0.475, 0.505]

    nstate_name='2pt'
    tmin_name='2pt'

    fit_name='2'
    xlabel=c2pt_tmin

    situation_list = ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode='8e08f23bc983bf0fa9778157733d8235', include_2pt=True, include_3pt=False, include_sum=False, 
    pt2_tmax=18, 
    pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=3,
    sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5)

    late_23_tmin_plot(E0_ylim, n_range, t_range, best_n, best_t, nstate_name, tmin_name, situation_list, fit_name, xlabel)

def prior_width_23s():
    # 2pt+3pt+sum varying g.s. and 1 ex prior width
    value={}
    value['Q']=[]
    value['A3']=[]
    value['A3_err']=[]

    # gs01, gs05, gs1, gs2, gs10
    for hexcode in ['a6d7c8abd3ab957e15af222a35f1f618', '70101230395214272d3ea3a896359b68', '8e08f23bc983bf0fa9778157733d8235', 'a68192fce72d0097536b9f00721ce302', '722949c46e4fd76d4b718d12631a1653']:
        for situation in ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode=hexcode,
        include_2pt=True, include_3pt=True, include_sum=True, 
        pt2_tmax=18, pt2_nstates=5,
        pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=9,
        sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5):
                    value['Q'].append(situation.Q_value)
                    value['A3'].append(situation.A300)
                    value['A3_err'].append(situation.A300_err)
                    print(situation.prior_hexcode)
                

    value['Q_']=[]
    value['A3_']=[]
    value['A3_err_']=[]

    # gs and 1 ex
    for hexcode in ['8b7225f794495bc718ed3a8c0ce5b26c', 'eb82329129cbfa975b81b08755b0d4dd', '8e08f23bc983bf0fa9778157733d8235', '1e8538003f62d1f766972446f1b08c27', '73ae4b7fa602c032c3e5f2c476d919f4']:
        for situation in ff_fit_FF_fit.objects.filter(data_file_name='a09m310_e_gA_srcs0-15.h5', prior_hexcode=hexcode,
        include_2pt=True, include_3pt=True, include_sum=True, 
        pt2_tmax=18, pt2_nstates=5,
        pt3_A3_tsep_min=3, pt3_A3_tsep_max=15, id_num=9,
        sum_A3_tsep_min=3, sum_A3_tsep_max=14, sum_tau_cut=1, sum_nstates=5):
                    value['Q_'].append(situation.Q_value)
                    value['A3_'].append(situation.A300)
                    value['A3_err_'].append(situation.A300_err)
                    print(situation.prior_hexcode)

    prior_width_plot(value)

# %%
if __name__ == "__main__":
    sum_tmin_2s()

# %%

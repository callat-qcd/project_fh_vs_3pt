# %%
import h5py as h5 
import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  
import lsqfit as lsf
import math
import os 
import hashlib

from prior_setting import prior_con
prior = prior_con
# %%
class Prepare_data():
    def __init__(self, file_name, file_path, pt2_data_range, pt3_data_range):
        self.file_name = file_name
        self.file_path = file_path
        self.pt2_data_range = pt2_data_range
        self.pt3_data_range = pt3_data_range

    def rotate(self, lis, n): # used to rotate the list
        return lis[n:] + lis[:n]

    def find_key(self, part_of_file, name_of_key): # find the keys that named with specific string in the dictionary 
        i = 0
        list_of_key = []
        for key in part_of_file:
            if name_of_key in key:
                list_of_key.append(key)
                i += 1

        if self.pt3_data_range[1] > 10:
            list_of_key = self.rotate(list_of_key, self.pt3_data_range[1]-10)

        return list_of_key

    def read_data_with_average(self):
        myfile = h5.File(self.file_path,'r')

        two_point_ml0p00951 = myfile['gf1p0_w3p5_n45_M51p1_L56_a1p5']['spec']['ml0p00951'] # different data file need different path here

        two_point_proton = np.array([two_point_ml0p00951['proton']['px0_py0_pz0']['spin_dn'][:],
        two_point_ml0p00951['proton']['px0_py0_pz0']['spin_up'][:]])
        two_point_proton_np = np.array([two_point_ml0p00951['proton_np']['px0_py0_pz0']['spin_dn'][:],
        two_point_ml0p00951['proton_np']['px0_py0_pz0']['spin_up'][:]])
        
        two_point_data = (two_point_proton[0] + two_point_proton[1] + two_point_proton_np[0] + two_point_proton_np[1])/4
        two_point_data = np.squeeze(two_point_data)

        three_point_ml0p00951 = myfile['gf1p0_w3p5_n45_M51p1_L56_a1p5']['formfac']['ml0p00951']

        proton_DD_dn_dn = self.find_key(three_point_ml0p00951, 'proton_DD_dn_dn')
        proton_DD_up_up = self.find_key(three_point_ml0p00951, 'proton_DD_up_up')
        proton_np_DD_dn_dn = self.find_key(three_point_ml0p00951, 'proton_np_DD_dn_dn')
        proton_np_DD_up_up = self.find_key(three_point_ml0p00951, 'proton_np_DD_up_up')

        proton_UU_dn_dn = self.find_key(three_point_ml0p00951, 'proton_UU_dn_dn')
        proton_UU_up_up = self.find_key(three_point_ml0p00951, 'proton_UU_up_up')
        proton_np_UU_dn_dn = self.find_key(three_point_ml0p00951, 'proton_np_UU_dn_dn')
        proton_np_UU_up_up = self.find_key(three_point_ml0p00951, 'proton_np_UU_up_up')
        
    
        data_dict = {}

        start = 0
        end = 392

        #temp = two_point_data[::2]

        for i in range(len(two_point_data[0])):
            data_dict['pt2_tsep_' + str(i)] = two_point_data[start:end, i, 0] # 2pt ss data    
            #data_dict['pt2_tsep_' + str(i)] = temp[:, i, 0]
            ## shape = (784,)


        pro_up = {}
        pro_dn = {}
        pnp_dn = {}
        pnp_up = {}

        for i in range(self.pt3_data_range[0], self.pt3_data_range[1]): # average
            for key in ['A3', 'V4']:
                pro_up[key] = three_point_ml0p00951[proton_UU_up_up[i-self.pt3_data_range[0]]][key]['px0_py0_pz0']['local_curr'][start:end] - three_point_ml0p00951[proton_DD_up_up[i-self.pt3_data_range[0]]][key]['px0_py0_pz0']['local_curr'][start:end]

                pro_dn[key] = three_point_ml0p00951[proton_UU_dn_dn[i-self.pt3_data_range[0]]][key]['px0_py0_pz0']['local_curr'][start:end] - three_point_ml0p00951[proton_DD_dn_dn[i-self.pt3_data_range[0]]][key]['px0_py0_pz0']['local_curr'][start:end]

                pnp_dn[key] = three_point_ml0p00951[proton_np_UU_dn_dn[i-self.pt3_data_range[0]]][key]['px0_py0_pz0']['local_curr'][start:end] - three_point_ml0p00951[proton_np_DD_dn_dn[i-self.pt3_data_range[0]]][key]['px0_py0_pz0']['local_curr'][start:end]

                pnp_up[key] = three_point_ml0p00951[proton_np_UU_up_up[i-self.pt3_data_range[0]]][key]['px0_py0_pz0']['local_curr'][start:end] - three_point_ml0p00951[proton_np_DD_up_up[i-self.pt3_data_range[0]]][key]['px0_py0_pz0']['local_curr'][start:end]

            
            data_dict['pt3_A3_tsep_' + str(i)] = ( np.imag(pro_up['A3']) - np.imag(pro_dn['A3']) + np.imag(pnp_dn['A3']) - np.imag(pnp_up['A3']) )/4

            data_dict['pt3_V4_tsep_' + str(i)] = ( np.real(pro_up['V4']) + np.real(pro_dn['V4']) - np.real(pnp_dn['V4']) - np.real(pnp_up['V4']) )/4


            ## shape = (784, tsep+1)

        print(np.shape(data_dict['pt2_tsep_1']))
        print(np.shape(data_dict['pt3_A3_tsep_3']))

        data_avg_dict = gv.dataset.avg_data(data_dict)

        return data_avg_dict


    def add_sum_data(self, data_avg_dict, sum_tau_cut):
        for i in range(self.pt3_data_range[0], self.pt3_data_range[1]):
            data_avg_dict['sum_A3_tsep_' + str(i)] = gv.gvar(0, 0)
            data_avg_dict['sum_V4_tsep_' + str(i)] = gv.gvar(0, 0)

            for j in range(sum_tau_cut, i-sum_tau_cut+1): # do the summation
                data_avg_dict['sum_A3_tsep_' + str(i)] += data_avg_dict['pt3_A3_tsep_' + str(i)][j]
                data_avg_dict['sum_V4_tsep_' + str(i)] += data_avg_dict['pt3_V4_tsep_' + str(i)][j]

        for i in range(self.pt3_data_range[0], self.pt3_data_range[1]-1): # use ratio form to do the fit
            # fit_14 will use tsep_15 here, so the tsep range of sum part should not be bigger than 13
            data_avg_dict['sum_A3_fit_'+str(i)] = (data_avg_dict['sum_A3_tsep_' + str(i+1)]/data_avg_dict['pt2_tsep_'+str(i+1)]) - (data_avg_dict['sum_A3_tsep_' + str(i)]/data_avg_dict['pt2_tsep_'+str(i)])

            data_avg_dict['sum_V4_fit_'+str(i)] = (data_avg_dict['sum_V4_tsep_' + str(i+1)]/data_avg_dict['pt2_tsep_'+str(i+1)]) - (data_avg_dict['sum_V4_tsep_' + str(i)]/data_avg_dict['pt2_tsep_'+str(i)])

        return data_avg_dict
        
# %%
class Fit():
    def __init__(self, file_name, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum):
        self.file_name = file_name
        
        self.pt2_nstates = pt2_nstates
        self.pt3_nstates = pt3_nstates
        self.sum_nstates = sum_nstates
        self.sum_tau_cut = sum_tau_cut
        
        self.include_2pt = include_2pt
        self.include_3pt = include_3pt
        self.include_sum = include_sum
        
        self.prior = prior(self.pt2_nstates, self.pt3_nstates, self.sum_nstates)

    def pt2_fit_function(self, pt2_t, p, pt2_nstates=None):
        E_list = {}
        
        if pt2_nstates==None: # not None when divide sum by 2pt
            pt2_nstates = self.pt2_nstates
            not_sum = True
            
        else:
            not_sum = False
        
        for i in range(pt2_nstates): #initialize       
            E_list['E'+str(i)] = p['E0']

        for i in range(1, pt2_nstates): #define Ei      
            for j in range(1, i+1):
                    E_list['E'+str(i)] += p['dE'+str(j)]

        val = {}

        val['pt2'] = p['z0']*p['z0']*np.exp(-E_list['E0']*pt2_t)

        for i in range(1, pt2_nstates-1):
            val['pt2'] += p['z'+str(i)]*p['z'+str(i)]*np.exp(-E_list['E'+str(i)]*pt2_t)
            
        if not_sum==True:
            val['pt2'] += p['z'+str(pt2_nstates-1)]*p['z'+str(pt2_nstates-1)]*np.exp(-E_list['E'+str(pt2_nstates-1)]*pt2_t)
            
        else:
            if self.sum_nstates == 1:
                return val
            
            E_list['E_sum'] = E_list['E'+ str(self.sum_nstates - 2)] + p['E_sum'] # use garbage can of sum
            
            val['pt2'] += p['z_sum']*p['z_sum']*np.exp(-E_list['E_sum']*pt2_t)

        return val

    def pt3_fit_function(self, pt3_tsep_A3, pt3_tsep_V4, pt3_tau_A3, pt3_tau_V4, p):
        E_list = {}
        for i in range(self.pt3_nstates): #initialize       
            E_list['E'+str(i)] = p['E0']

        for i in range(1, self.pt3_nstates): #define Ei      
            for j in range(1, i+1):
                    E_list['E'+str(i)] += p['dE'+str(j)]

        val = {}

        val['pt3_A3'] = p['A3_00']*p['z0']*p['z0']*np.exp(-E_list['E0']*pt3_tsep_A3) 

        val['pt3_V4'] = p['V4_00']*p['z0']*p['z0']*np.exp(-E_list['E0']*pt3_tsep_V4) 

        for i in range(self.pt3_nstates):    
            for j in range(self.pt3_nstates):
                if i+j >= 1:
                    if j == i:
                        val['pt3_A3'] += p['A3_'+str(j)+str(i)]*p['z'+str(j)]*p['z'+str(i)]*np.exp(-E_list['E'+str(j)]*pt3_tsep_A3)

                        val['pt3_V4'] += p['V4_'+str(j)+str(i)]*p['z'+str(j)]*p['z'+str(i)]*np.exp(-E_list['E'+str(j)]*pt3_tsep_V4)

                    else:
                        mi = np.minimum(j, i)
                        ma = np.maximum(j, i)
                        val['pt3_A3'] += p['A3_'+str(mi)+str(ma)]*p['z'+str(j)]*p['z'+str(i)]*np.exp(-E_list['E'+str(j)]*pt3_tsep_A3)*np.exp((E_list['E'+str(j)]-E_list['E'+str(i)])*pt3_tau_A3)

                        val['pt3_V4'] += p['V4_'+str(mi)+str(ma)]*p['z'+str(j)]*p['z'+str(i)]*np.exp(-E_list['E'+str(j)]*pt3_tsep_V4)*np.exp((E_list['E'+str(j)]-E_list['E'+str(i)])*pt3_tau_V4)

        return val

    def summation(self, A3_t, V4_t, p):
        E_list = {}
        for i in range(self.sum_nstates): #initialize       
            E_list['E'+str(i)] = p['E0']

        for i in range(1, self.sum_nstates): #define Ei      
            for j in range(1, i+1):
                    E_list['E'+str(i)] += p['dE'+str(j)]

        if self.sum_nstates == 1:
            cut = self.sum_tau_cut

            val = {}

            val['sum_A3'] = p['z0'] * p['A3_00'] * p['z0'] * np.exp(-E_list['E0'] * A3_t) * (A3_t - 2*cut + 1)

            val['sum_V4'] = p['z0'] * p['V4_00'] * p['z0'] * np.exp(-E_list['E0'] * V4_t) * (V4_t - 2*cut + 1)

            return val

        E_list['E_sum'] = E_list['E'+ str(self.sum_nstates - 2)] + p['E_sum'] 

        E = np.zeros([self.sum_nstates-1], dtype=gv.GVar)
        z = np.zeros([self.sum_nstates-1], dtype=gv.GVar)
        for i in range(self.sum_nstates-1):
            E[i] = E_list['E'+str(i)]
            z[i] = p['z'+str(i)]

        E_sum = E_list['E_sum']
        z_sum = p['z_sum']

        D = np.zeros([self.sum_nstates-1, self.sum_nstates-1], dtype=gv.GVar)
        for i in range(self.sum_nstates-1):
            for j in range(self.sum_nstates-1):
                D[i][j] = E[i] - E[j]
        
        A3 = np.zeros([self.sum_nstates-1, self.sum_nstates-1], dtype=gv.GVar)
        V4 = np.zeros([self.sum_nstates-1, self.sum_nstates-1], dtype=gv.GVar)
        sumA3 = np.zeros([self.sum_nstates], dtype=gv.GVar)
        sumV4 = np.zeros([self.sum_nstates], dtype=gv.GVar)
        for i in range(self.sum_nstates-1):
            for j in range(self.sum_nstates-1):
                mi = np.minimum(j, i)
                ma = np.maximum(j, i)
                A3[i][j] = p['A3_'+str(mi)+str(ma)]
                V4[i][j] = p['V4_'+str(mi)+str(ma)]

        for i in range(self.sum_nstates):
            sumA3[i] = p['sum_A3_'+str(i)]
            sumV4[i] = p['sum_V4_'+str(i)]

        cut = self.sum_tau_cut

        val = {}

        val['sum_A3'] = z[0] * A3[0][0] * z[0] * np.exp(-E[0] * A3_t) * (A3_t - 2*cut + 1)

        val['sum_V4'] = z[0] * V4[0][0] * z[0] * np.exp(-E[0] * V4_t) * (V4_t - 2*cut + 1)


        for i in range(self.sum_nstates-1):
            for j in range(self.sum_nstates-1):
                if i+j >= 1:
                    if j == i: 
                        val['sum_A3'] += z[j] * A3[j][i] * z[i] * np.exp(-E[j] * A3_t) * (A3_t - 2*cut + 1)

                        val['sum_V4'] += z[j] * V4[j][i] * z[i] * np.exp(-E[j] * V4_t) * (V4_t - 2*cut + 1) 

                    else:
                        val['sum_A3'] += z[j] * A3[j][i] * z[i] * np.exp(-E[j] * A3_t) * ((np.exp(cut * D[j][i]) * (1 - np.exp((A3_t - 2*cut + 1) * D[j][i]) )) / (1 - np.exp(D[j][i]) ))

                        val['sum_V4'] += z[j] * V4[j][i] * z[i] * np.exp(-E[j] * V4_t) * ((np.exp(cut * D[j][i]) * (1 - np.exp((V4_t - 2*cut + 1) * D[j][i]) )) / (1 - np.exp(D[j][i]) ))
                        

        for i in range(self.sum_nstates-1):
            val['sum_A3'] += z_sum * sumA3[i] * z[i] * np.exp(-E_sum * A3_t) * (((np.exp(cut * (E_sum - E[i])) ) * (1 - np.exp((A3_t - 2*cut + 1) * (E_sum - E[i])) )) / (1 - np.exp(E_sum - E[i]) ))

            val['sum_A3'] += z[i] * sumA3[i] * z_sum * np.exp(-E[i] * A3_t) * (((np.exp(cut * (E[i] - E_sum)) ) * (1 - np.exp((A3_t - 2*cut + 1) * (E[i] - E_sum)) )) / (1 - np.exp(E[i] - E_sum) ))

            val['sum_V4'] += z_sum * sumV4[i] * z[i] * np.exp(-E_sum * V4_t) * (((np.exp(cut * (E_sum - E[i])) ) * (1 - np.exp((V4_t - 2*cut + 1) * (E_sum - E[i])) )) / (1 - np.exp(E_sum - E[i]) ))

            val['sum_V4'] += z[i] * sumV4[i] * z_sum * np.exp(-E[i] * V4_t) * (((np.exp(cut * (E[i] - E_sum)) ) * (1 - np.exp((V4_t - 2*cut + 1) * (E[i] - E_sum)) )) / (1 - np.exp(E[i] - E_sum) ))
            

        val['sum_A3'] += z_sum * sumA3[self.sum_nstates-1] * z_sum * np.exp(-E_sum * A3_t) * (A3_t - 2*cut + 1)

        val['sum_V4'] += z_sum * sumV4[self.sum_nstates-1] * z_sum * np.exp(-E_sum * V4_t) * (V4_t - 2*cut + 1)

        return val 


    def fcn(self, x, p):
        val = {}
        
        if self.include_2pt == True:
            pt2_t = x['pt2']

            val['pt2'] = self.pt2_fit_function(pt2_t, p)['pt2']

        if self.include_3pt == True:
            pt3_tsep_A3 = x['pt3_A3'][0]
            pt3_tau_A3 = x['pt3_A3'][1]
            pt3_tsep_V4 = x['pt3_V4'][0]
            pt3_tau_V4 = x['pt3_V4'][1]

            val['pt3_A3'] = self.pt3_fit_function(pt3_tsep_A3, pt3_tsep_V4, pt3_tau_A3, pt3_tau_V4, p)['pt3_A3']

            val['pt3_V4'] = self.pt3_fit_function(pt3_tsep_A3, pt3_tsep_V4, pt3_tau_A3, pt3_tau_V4, p)['pt3_V4']

        if self.include_sum == True:
            sum_tsep_A3 = x['sum_A3']
            sum_tsep_V4 = x['sum_V4']

            sum_tsep_A3_fit_1 = sum_tsep_A3
            sum_tsep_V4_fit_1 = sum_tsep_V4

            sum_tsep_A3_fit_2 = sum_tsep_A3 + 1
            sum_tsep_V4_fit_2 = sum_tsep_V4 + 1

            pt2_t_A3_fit_1 = sum_tsep_A3
            pt2_t_V4_fit_1 = sum_tsep_V4

            pt2_t_A3_fit_2 = sum_tsep_A3 + 1
            pt2_t_V4_fit_2 = sum_tsep_V4 + 1

            val['sum_A3'] = (self.summation(sum_tsep_A3_fit_2, sum_tsep_V4_fit_2, p)['sum_A3']/self.pt2_fit_function(pt2_t_A3_fit_2, p, self.sum_nstates)['pt2']) - (self.summation(sum_tsep_A3_fit_1, sum_tsep_V4_fit_1, p)['sum_A3']/self.pt2_fit_function(pt2_t_A3_fit_1, p, self.sum_nstates)['pt2']) # using same nstates in the ratio

            val['sum_V4'] = (self.summation(sum_tsep_A3_fit_2, sum_tsep_V4_fit_2, p)['sum_V4']/self.pt2_fit_function(pt2_t_V4_fit_2, p, self.sum_nstates)['pt2']) - (self.summation(sum_tsep_A3_fit_1, sum_tsep_V4_fit_1, p)['sum_V4']/self.pt2_fit_function(pt2_t_V4_fit_1, p, self.sum_nstates)['pt2'])
  
        return val

    def fit(self, data_avg_dict, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0, priors=None):
        if priors == None:
            priors = self.prior

        t_tsep_tau = {}
        Amp = {}

        if self.include_2pt == True:
            pt2_amp = []
            for t in pt2_t:
                pt2_amp.append(data_avg_dict['pt2_tsep_'+str(t)])

            t_tsep_tau['pt2'] = pt2_t
            Amp['pt2'] = np.array(pt2_amp)
        
        if self.include_3pt == True:
            pt3_A3_tsep = pt3_A3[0]
            pt3_A3_tau = pt3_A3[1]
            pt3_V4_tsep = pt3_V4[0]
            pt3_V4_tau = pt3_V4[1]

            pt3_A3_amp = []
            pt3_V4_amp = []

            for i in range(len(pt3_A3[0])):
                t = pt3_A3[0][i]
                tau = pt3_A3[1][i]

                pt3_A3_amp.append((data_avg_dict['pt3_A3_tsep_' + str(t)][tau] + data_avg_dict['pt3_A3_tsep_' + str(t)][t - tau])/2) # average tau=i and tau=tsep-i to make data symmetric

            for i in range(len(pt3_V4[0])):
                t = pt3_V4[0][i]
                tau = pt3_V4[1][i]

                pt3_V4_amp.append((data_avg_dict['pt3_V4_tsep_' + str(t)][tau] + data_avg_dict['pt3_V4_tsep_' + str(t)][t - tau])/2) # average tau=i and tau=tsep-i to make data symmetric


            t_tsep_tau['pt3_A3'] = [np.array(pt3_A3_tsep), np.array(pt3_A3_tau)]
            t_tsep_tau['pt3_V4'] = [np.array(pt3_V4_tsep), np.array(pt3_V4_tau)]
            Amp['pt3_A3'] = np.array(pt3_A3_amp)
            Amp['pt3_V4'] = np.array(pt3_V4_amp)

        if self.include_sum == True:
            sum_A3_factor=[]
            sum_V4_factor=[]

            for t in sum_A3:
                sum_A3_factor.append(data_avg_dict['sum_A3_fit_'+str(t)])

            for t in sum_V4:
                sum_V4_factor.append(data_avg_dict['sum_V4_fit_'+str(t)])

            t_tsep_tau['sum_A3'] = sum_A3
            t_tsep_tau['sum_V4'] = sum_V4

            Amp['sum_A3'] = np.array(sum_A3_factor)
            Amp['sum_V4'] = np.array(sum_V4_factor)

        #print(t_tsep_tau)
        #print(Amp)

        if best_p0 == 0:
            fit_result = lsf.nonlinear_fit(data=(t_tsep_tau, Amp), prior=priors, fcn=self.fcn, maxit=100000, fitter='scipy_least_squares') # scipy_least_squares   # gsl_multifit

        else :
            fit_result = lsf.nonlinear_fit(data=(t_tsep_tau, Amp), prior=priors, fcn=self.fcn, maxit=100000, fitter='scipy_least_squares', p0=best_p0)

        return fit_result

# %%
def data_plot_pt2(pt2_data_range, data_avg_dict):
    '''plot effective mass with 2pt data'''
    m0_eff = []
    m0_eff_err = []

    for i in range(pt2_data_range[0], pt2_data_range[1]-1):
        temp = gv.log(data_avg_dict['pt2_tsep_'+str(i)] / data_avg_dict['pt2_tsep_'+str(i+1)])
        m0_eff.append(temp.mean)
        m0_eff_err.append(temp.sdev)

    plt.figure(figsize=(10, 7))
    plt.errorbar(np.arange(pt2_data_range[0], pt2_data_range[1]-1), np.array(m0_eff), yerr=np.array(m0_eff_err), fmt='x', label='data')
    plt.grid()
    plt.minorticks_on()
    plt.xlim([2, 25])
    plt.ylim([0.45, 0.65])
    plt.xlabel('t', size=20)
    plt.ylabel('meff', size=20)
    plt.legend()
    plt.show()

    print(m0_eff[8])

    return [np.arange(pt2_data_range[0], pt2_data_range[1]-1), np.array(m0_eff), np.array(m0_eff_err)]

def data_plot_sum(pt3_data_range, data_avg_dict_completed):
    '''plot form factor with sum data'''
    gA = []
    gA_err = []

    for i in range(pt3_data_range[0], pt3_data_range[1]-1):
        temp = data_avg_dict_completed['sum_A3_fit_'+str(i)]
        gA.append(temp.mean)
        gA_err.append(temp.sdev)

    plt.figure(figsize=(10, 7))
    plt.errorbar(np.arange(pt3_data_range[0], pt3_data_range[1]-1), np.array(gA), yerr=np.array(gA_err), fmt='x', label='data')
    plt.grid()
    plt.minorticks_on()
    plt.xlim([2, 15])
    #plt.ylim([0.45, 0.65])
    plt.xlabel('tsep', size=20)
    plt.ylabel('gA', size=20)
    plt.legend()
    plt.show()

    print(gA[8])

    return [np.arange(pt3_data_range[0], pt3_data_range[1]-1), np.array(gA), np.array(gA_err)]

# %%
if __name__ == '__main__':
    file_name = 'a09m310_e_gA_srcs0-15.h5'
    file_path = os.getcwd() + '/' + file_name # just put the data file inside the 'my_project' folder, if no database, data file should be put in the same path as this file.

    pt2_data_range = [0, 96]
    pt3_data_range = [2, 15]

    prepare_data = Prepare_data(file_name, file_path, pt2_data_range, pt3_data_range)

    data_avg_dict = prepare_data.read_data_with_average()

    x_y_yerr_02 = data_plot_pt2(pt2_data_range, data_avg_dict)

# %%
# do the fit

#############################################################
################ set parameters #############################

pt2_t = np.array([4, 5, 6, 7, 8, 9])
pt2_nstates = 5

pt3_A3_t = [4, 4, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9]
pt3_A3_tau = [1, 2, 1, 2, 1, 2, 3, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4]
pt3_V4_t = pt3_A3_t
pt3_V4_tau = pt3_A3_tau

pt3_A3 = [np.array(pt3_A3_t), np.array(pt3_A3_tau)]
pt3_V4 = [np.array(pt3_V4_t), np.array(pt3_V4_tau)]
pt3_nstates = pt2_nstates

sum_A3 = np.array([7, 8, 9, 10, 11, 12])
sum_V4 = sum_A3 
sum_nstates = 1
sum_tau_cut = 1

use_p0 = False
include_2pt = True
include_3pt = False
include_sum = False

##############################################################
##############################################################

data_avg_dict_completed = prepare_data.add_sum_data(data_avg_dict, sum_tau_cut)


######### plot sum data ###########
temp = data_plot_sum(pt3_data_range, data_avg_dict_completed)
###################################


if use_p0 == False:
    best_p0 = 0

fitter = Fit(file_name, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum)

fit_result = fitter.fit(data_avg_dict_completed, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0)

best_p0 = {key: fit_result.p[key].mean for key in fit_result.p}   

#print(fit_result)#.format(100)) 

# %%
plt.figure(figsize=(10, 7))
plt.errorbar(x_y_yerr_02[0]-0.2, x_y_yerr_02[1], yerr=x_y_yerr_02[2], fmt='x', label='[0 : Ncfg/2]')
plt.errorbar(x_y_yerr_13[0], x_y_yerr_13[1], yerr=x_y_yerr_13[2], fmt='x', label='[Ncfg/4 : 3*Ncfg/4]')
plt.errorbar(x_y_yerr_24[0]+0.2, x_y_yerr_24[1], yerr=x_y_yerr_24[2], fmt='x', label='[Ncfg/2 : Ncfg]')
plt.grid()
plt.minorticks_on()
plt.xlim([0, 25])
plt.ylim([0.425, 0.625])
plt.xlabel('t', size=20)
plt.ylabel('meff', size=20)
plt.legend()
plt.show()


#%% 2pt data
# x_y_yerr_02 = [np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
#        17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
#        34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
#        51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
#        68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
#        85, 86, 87, 88, 89, 90, 91, 92, 93, 94]), np.array([-3.67684845,  1.07723362,  0.77885928,  0.65924149,  0.59895217,
#         0.56430859,  0.54205121,  0.52755443,  0.51810395,  0.51101949,
#         0.5059069 ,  0.50144945,  0.4980635 ,  0.49636072,  0.49560785,
#         0.49376855,  0.4914562 ,  0.49520114,  0.49956237,  0.49592251,
#         0.48830321,  0.48761648,  0.49505285,  0.51243485,  0.54529019,
#         0.59427716,  0.7143718 ,  1.03345144,         None, -1.73203927,
#        -0.29677057, -0.1101338 ,  0.12071698,  0.26767568,  0.29360272,
#         0.31561214,  0.3628247 ,  0.32981601,  0.24070773,  0.17893353,
#         0.14854797,  0.22649988,  0.38462836,  0.50633066,  0.53484998,
#         0.52472882,  0.16288408,  0.23988648,  0.48658953,  0.41233572,
#        -0.37068323,  0.88011729,  2.4716216 , -3.3622326 , -0.79165739,
#        -0.06104182,  3.38510837,         None, -0.67849918, -0.20383616,
#         0.03387748, -0.08992951, -0.53118861, -0.52340753, -0.29957456,
#        -0.10194997,  0.12404368,  0.50096342,  0.63397353,  0.81585827,
#                None, -1.66446424, -0.48382026,  1.19051196,         None,
#        -1.16280734, -0.58497046, -0.07553343,  0.29178016,  0.19497202,
#        -0.79550904, -1.05315846, -0.83147703, -0.71260434, -0.72480144,
#        -0.71946489, -0.73953997, -0.77783562, -0.81186969, -0.86040539,
#        -0.92786698, -1.01448793, -1.13058564, -1.2999641 , -1.6597967 ]), np.array([7.92838616e-03, 9.78157337e-04, 7.82939341e-04, 7.62453072e-04,
#        8.11891923e-04, 9.00349875e-04, 1.01153970e-03, 1.16056317e-03,
#        1.29578985e-03, 1.40792940e-03, 1.55063082e-03, 1.81919680e-03,
#        2.20230649e-03, 2.64857848e-03, 3.32753412e-03, 4.01360222e-03,
#        4.73345098e-03, 6.05869807e-03, 8.30822908e-03, 1.16813897e-02,
#        1.58843703e-02, 2.14384602e-02, 2.85506941e-02, 3.79147459e-02,
#        5.58425657e-02, 8.71089623e-02, 1.58224848e-01, 4.61284961e-01,
#        7.02419891e+00, 5.42265012e+00, 4.15778433e-01, 2.01521893e-01,
#        1.12945096e-01, 1.06165512e-01, 1.10972573e-01, 1.27397675e-01,
#        1.64475069e-01, 1.92206249e-01, 1.87062048e-01, 1.82833424e-01,
#        1.85598513e-01, 1.85326498e-01, 2.53474560e-01, 4.08199623e-01,
#        6.32757982e-01, 9.32877132e-01, 7.62184396e-01, 8.95316886e-01,
#        1.70873458e+00, 2.30189435e+00, 2.54780256e+00, 5.66513541e+00,
#        1.02832213e+02, 1.06962694e+02, 2.31289518e+00, 8.23207187e-01,
#        8.56943888e+01, 8.99511731e+01, 9.43288842e-01, 3.83054138e-01,
#        4.67267019e-01, 5.08998688e-01, 6.06301571e-01, 4.01391489e-01,
#        2.69046169e-01, 2.78646996e-01, 4.25995585e-01, 1.13273798e+00,
#        2.81746437e+00, 8.74246350e+00, 2.11017551e+01, 6.72753593e+00,
#        7.94790458e-01, 4.73369043e+00, 8.30387569e+00, 1.49548242e+00,
#        4.04525865e-01, 2.68155504e-01, 5.80971225e-01, 8.36295026e-01,
#        9.22612195e-01, 6.29448292e-01, 2.39480222e-01, 1.19725560e-01,
#        7.76222192e-02, 4.77398076e-02, 2.88909186e-02, 1.92753454e-02,
#        1.22003270e-02, 7.57268487e-03, 5.17212745e-03, 3.70535652e-03,
#        2.86175743e-03, 2.02692254e-03, 1.43330295e-03])]

# x_y_yerr_13 = [np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
#        17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
#        34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
#        51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
#        68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
#        85, 86, 87, 88, 89, 90, 91, 92, 93, 94]), np.array([-3.67661887,  1.07824788,  0.77976848,  0.65954264,  0.59860357,
#         0.56318155,  0.5403244 ,  0.52497289,  0.51495362,  0.50807042,
#         0.50371206,  0.50150977,  0.49974566,  0.4982125 ,  0.49752893,
#         0.49827102,  0.49811171,  0.4964816 ,  0.49192757,  0.48839   ,
#         0.49079145,  0.48196505,  0.46450398,  0.45340407,  0.44927058,
#         0.45968435,  0.53157294,  0.6390065 ,  0.7604725 ,  0.75219788,
#         1.15451225,         None, -0.73514649,  0.04848123,  0.32440039,
#         0.60707723,  0.88354902,  0.14591195, -0.22439623,  0.23397238,
#         1.11892371,         None, -0.67186828,  0.01405469,  0.39947616,
#         0.34154419, -0.05013246, -0.2115552 , -0.09754741,  0.04803926,
#         0.31300955,  0.84927593,         None, -0.97604536, -0.09583016,
#         0.16605254,  0.704027  ,         None, -0.8905558 , -0.00740986,
#         0.74603459,         None, -0.3016409 ,  0.175948  ,         None,
#        -0.44745168,         None, -1.1810101 ,  0.03006279,         None,
#        -1.67753119,  0.75551258,         None, -0.37173313,  0.2126276 ,
#         2.85103934,         None,  1.05661296,         None, -0.48618324,
#         0.01888143,         None, -2.23425052, -1.14966603, -0.89424967,
#        -0.78882707, -0.75101806, -0.78036047, -0.82568615, -0.8744083 ,
#        -0.93949889, -1.02323691, -1.13514412, -1.30179499, -1.65926648]), np.array([1.08371812e-02, 9.64807037e-04, 7.84021253e-04, 7.76723269e-04,
#        7.91987868e-04, 8.60570391e-04, 9.79128554e-04, 1.11891420e-03,
#        1.27950374e-03, 1.45145821e-03, 1.60232314e-03, 1.83558912e-03,
#        2.17778384e-03, 2.52686148e-03, 3.24225185e-03, 4.25830961e-03,
#        5.36264238e-03, 6.81988612e-03, 8.86661094e-03, 1.16149378e-02,
#        1.54141304e-02, 1.96585577e-02, 2.50804215e-02, 3.03263108e-02,
#        3.79430240e-02, 5.07356646e-02, 7.74230019e-02, 1.33164553e-01,
#        2.63797262e-01, 4.52772346e-01, 1.56086231e+00, 4.97102110e+00,
#        1.64030444e+00, 3.00699341e-01, 2.61804218e-01, 5.53298360e-01,
#        1.45633954e+00, 8.06538148e-01, 1.10380384e+00, 5.02044016e-01,
#        2.74423679e+00, 9.98296630e+00, 3.49415137e+00, 7.74535454e-01,
#        9.12510067e-01, 1.11022221e+00, 1.09115250e+00, 1.06390381e+00,
#        7.22621600e-01, 5.94852347e-01, 8.98458507e-01, 3.28115006e+00,
#        8.98049504e+00, 2.46248701e+00, 5.52152421e-01, 8.41436646e-01,
#        3.16818030e+00, 7.88724351e+00, 1.68423765e+00, 5.70107891e-01,
#        2.58896021e+00, 2.38240047e+01, 6.95284575e+00, 1.13560277e+01,
#        3.85533405e+01, 5.17075806e+00, 1.35742949e+01, 2.83390076e+00,
#        8.38092576e-01, 1.05703414e+01, 6.83054478e+00, 3.58800929e+00,
#        7.76788870e+00, 9.27640653e-01, 1.50794294e+00, 6.69631504e+01,
#        7.29154746e+01, 8.61127133e+00, 1.38835734e+01, 9.69139534e-01,
#        9.48769733e-01, 6.49465642e+00, 3.56015137e+00, 3.62226504e-01,
#        1.26368222e-01, 5.94803698e-02, 3.20035521e-02, 2.00587650e-02,
#        1.25476327e-02, 7.72655168e-03, 4.92114610e-03, 3.38960577e-03,
#        2.78856406e-03, 1.97644769e-03, 1.35049856e-03])]

# x_y_yerr_24 = [np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
#        17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
#        34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
#        51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
#        68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
#        85, 86, 87, 88, 89, 90, 91, 92, 93, 94]), np.array([-3.68574404,  1.07971123,  0.7804216 ,  0.65893883,  0.5968307 ,
#         0.56118943,  0.53939589,  0.5249932 ,  0.51517045,  0.50755545,
#         0.5019771 ,  0.49924372,  0.4983276 ,  0.49940397,  0.50024855,
#         0.50072991,  0.50180989,  0.49712188,  0.48973929,  0.48811239,
#         0.48795376,  0.47427167,  0.45297091,  0.42964772,  0.43089586,
#         0.46427435,  0.51744391,  0.55141324,  0.4754617 ,  0.29436391,
#         0.24111136,  0.3019304 ,  0.39990933,  0.39389532,  0.30704483,
#         0.33538677,  0.37084553,  0.65546909,  2.21163188,         None,
#         0.28780877,  2.37235767,         None, -0.28012711, -0.05423558,
#        -0.25597821, -0.48096263, -0.35225594, -0.07666802,  0.17455587,
#         0.45299775,  0.94289765,         None, -0.93557695,         None,
#        -1.33979768,  0.44505352,  1.26059793, -1.7852241 , -0.18894697,
#         3.88124602,         None,  0.83307056,         None, -2.61805128,
#        -0.85489108, -0.32508527, -0.03524922, -0.23206244, -1.00797742,
#        -0.66370733, -0.02324355,  1.81691375,         None, -0.43352135,
#        -0.05862518,  0.16864161, -0.48736331, -0.50971633, -0.05541892,
#         0.20557539,         None, -2.57796249, -1.1924514 , -0.90373059,
#        -0.80860564, -0.77359211, -0.79489586, -0.83170868, -0.8773691 ,
#        -0.94025889, -1.02324491, -1.13597782, -1.30245018, -1.65906219]), np.array([1.09984292e-02, 9.21243229e-04, 7.55733444e-04, 7.55626147e-04,
#        7.74741359e-04, 8.29195134e-04, 9.56518376e-04, 1.07241658e-03,
#        1.20962809e-03, 1.36688273e-03, 1.51180913e-03, 1.72757881e-03,
#        2.08326227e-03, 2.55150381e-03, 3.25802214e-03, 4.30659133e-03,
#        5.61529559e-03, 7.12339655e-03, 8.92196130e-03, 1.19117519e-02,
#        1.59322106e-02, 1.97742274e-02, 2.31779118e-02, 2.62611980e-02,
#        3.26188082e-02, 4.55895440e-02, 7.14530967e-02, 1.04188187e-01,
#        1.15186829e-01, 1.02636027e-01, 1.04421432e-01, 1.17993742e-01,
#        1.61058067e-01, 1.91199544e-01, 2.00890518e-01, 2.31494392e-01,
#        2.66928605e-01, 5.88820581e-01, 8.90715611e+00, 1.27250901e+01,
#        9.10359007e-01, 2.31462236e+01, 3.08495801e+01, 2.28855333e+00,
#        1.15950531e+00, 1.15055618e+00, 9.37724155e-01, 4.43616450e-01,
#        2.18855217e-01, 2.18984494e-01, 4.34519183e-01, 1.58378422e+00,
#        1.30287670e+01, 6.65910351e+00, 1.13641632e+01, 5.10639451e+00,
#        2.01846889e+00, 1.23848183e+01, 1.30379900e+01, 1.16181655e+00,
#        1.71595392e+02, 1.81405243e+02, 1.45136933e+01, 7.40015487e+01,
#        4.93320301e+01, 2.67240404e+00, 9.72982742e-01, 1.02183653e+00,
#        9.90242635e-01, 1.48778354e+00, 5.34427150e-01, 3.50048567e-01,
#        6.20973045e+00, 8.63916644e+00, 6.50812187e-01, 5.63596823e-01,
#        9.75168188e-01, 8.79917094e-01, 6.69863279e-01, 5.55102049e-01,
#        1.06533685e+00, 8.73262679e+00, 5.88517522e+00, 4.15174671e-01,
#        1.34090579e-01, 6.22421964e-02, 3.30561411e-02, 1.98838552e-02,
#        1.25001036e-02, 7.78283240e-03, 5.09108191e-03, 3.40368458e-03,
#        2.67411390e-03, 1.89205305e-03, 1.30915702e-03])]

# %%  sum data
#x_y_yerr_02 = [np.array([ 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13]), np.array([1.06353772, 1.08438903, 1.11961391, 1.15138211, 1.17700582, 1.2002507 , 1.21912177, 1.23286619, 1.24894384, 1.27104163, 1.29518782, 1.29280979]), np.array([0.00119628, 0.00158653, 0.00224862, 0.00320557, 0.00462822, 0.00674897, 0.00966203, 0.01323506, 0.01721373, 0.02345746, 0.03266905, 0.04620311])]

#x_y_yerr_13 = [np.array([ 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13]), np.array([1.0624236 , 1.0833115 , 1.11797503, 1.1469344 , 1.17280742, 1.19588988, 1.21522291, 1.22442263, 1.2323056 , 1.25085608, 1.27411004, 1.29440225]), np.array([0.00119577, 0.00163956, 0.00227257, 0.00321268, 0.00455229, 0.00643394, 0.00917139, 0.01276198, 0.01682244, 0.02300662, 0.03241269, 0.04480563])]

#x_y_yerr_24 = [np.array([ 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13]), np.array([1.06119903, 1.08221049, 1.11683807, 1.14674461, 1.17247396, 1.19356657, 1.2125429 , 1.21730622, 1.22227836, 1.23578339, 1.2452322 , 1.28339873]), np.array([0.00118336, 0.00163121, 0.00235201, 0.00340018, 0.00480871, 0.00652394, 0.0086798 , 0.01186616, 0.01606915, 0.02185578, 0.03104191, 0.04290947])]
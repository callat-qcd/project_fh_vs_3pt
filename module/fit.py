import h5py as h5 
import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  
import lsqfit as lsf
import math
import hashlib

class Fit():
    def __init__(self, file_name, prior, pt2_nstates, pt3_nstates, sum_nstates, sum_tau_cut,  include_2pt, include_3pt, include_sum):
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
        
        if pt2_nstates == 1:
            return val

        for i in range(1, pt2_nstates-1):
            val['pt2'] += p['z'+str(i)]*p['z'+str(i)]*np.exp(-E_list['E'+str(i)]*pt2_t)
            
        if not_sum==True and pt2_nstates > 1:
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
    
    def summation_same_can(self, A3_t, V4_t, p):
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

        E = np.zeros([self.sum_nstates], dtype=gv.GVar)
        z = np.zeros([self.sum_nstates], dtype=gv.GVar)
        for i in range(self.sum_nstates):
            E[i] = E_list['E'+str(i)]
            z[i] = p['z'+str(i)]

        D = np.zeros([self.sum_nstates, self.sum_nstates], dtype=gv.GVar)
        for i in range(self.sum_nstates):
            for j in range(self.sum_nstates):
                D[i][j] = E[i] - E[j]
        
        A3 = np.zeros([self.sum_nstates, self.sum_nstates], dtype=gv.GVar)
        V4 = np.zeros([self.sum_nstates, self.sum_nstates], dtype=gv.GVar)

        for i in range(self.sum_nstates):
            for j in range(self.sum_nstates):
                mi = np.minimum(j, i)
                ma = np.maximum(j, i)
                A3[i][j] = p['A3_'+str(mi)+str(ma)]
                V4[i][j] = p['V4_'+str(mi)+str(ma)]

        cut = self.sum_tau_cut

        val = {}

        val['sum_A3'] = z[0] * A3[0][0] * z[0] * np.exp(-E[0] * A3_t) * (A3_t - 2*cut + 1)

        val['sum_V4'] = z[0] * V4[0][0] * z[0] * np.exp(-E[0] * V4_t) * (V4_t - 2*cut + 1)


        for i in range(self.sum_nstates):
            for j in range(self.sum_nstates):
                if i+j >= 1:
                    if j == i: 
                        val['sum_A3'] += z[j] * A3[j][i] * z[i] * np.exp(-E[j] * A3_t) * (A3_t - 2*cut + 1)

                        val['sum_V4'] += z[j] * V4[j][i] * z[i] * np.exp(-E[j] * V4_t) * (V4_t - 2*cut + 1) 

                    else:  ### transition terms
                        val['sum_A3'] += z[j] * A3[j][i] * z[i] * np.exp(-E[j] * A3_t) * ((np.exp(cut * D[j][i]) * (1 - np.exp((A3_t - 2*cut + 1) * D[j][i]) )) / (1 - np.exp(D[j][i]) ))

                        val['sum_V4'] += z[j] * V4[j][i] * z[i] * np.exp(-E[j] * V4_t) * ((np.exp(cut * D[j][i]) * (1 - np.exp((V4_t - 2*cut + 1) * D[j][i]) )) / (1 - np.exp(D[j][i]) ))

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

            if self.pt2_nstates != self.sum_nstates:
                val['sum_A3'] = (self.summation(sum_tsep_A3_fit_2, sum_tsep_V4_fit_2, p)['sum_A3']/self.pt2_fit_function(pt2_t_A3_fit_2, p, self.sum_nstates)['pt2']) - (self.summation(sum_tsep_A3_fit_1, sum_tsep_V4_fit_1, p)['sum_A3']/self.pt2_fit_function(pt2_t_A3_fit_1, p, self.sum_nstates)['pt2']) # using same nstates in the ratio

                val['sum_V4'] = (self.summation(sum_tsep_A3_fit_2, sum_tsep_V4_fit_2, p)['sum_V4']/self.pt2_fit_function(pt2_t_V4_fit_2, p, self.sum_nstates)['pt2']) - (self.summation(sum_tsep_A3_fit_1, sum_tsep_V4_fit_1, p)['sum_V4']/self.pt2_fit_function(pt2_t_V4_fit_1, p, self.sum_nstates)['pt2'])
  
            if self.pt2_nstates == self.sum_nstates:
                val['sum_A3'] = (self.summation_same_can(sum_tsep_A3_fit_2, sum_tsep_V4_fit_2, p)['sum_A3']/self.pt2_fit_function(pt2_t_A3_fit_2, p)['pt2']) - (self.summation_same_can(sum_tsep_A3_fit_1, sum_tsep_V4_fit_1, p)['sum_A3']/self.pt2_fit_function(pt2_t_A3_fit_1, p)['pt2'])

                val['sum_V4'] = (self.summation_same_can(sum_tsep_A3_fit_2, sum_tsep_V4_fit_2, p)['sum_V4']/self.pt2_fit_function(pt2_t_V4_fit_2, p)['pt2']) - (self.summation_same_can(sum_tsep_A3_fit_1, sum_tsep_V4_fit_1, p)['sum_V4']/self.pt2_fit_function(pt2_t_V4_fit_1, p)['pt2'])

        return val

    def fit(self, data_avg_dict, pt2_t, pt3_A3, pt3_V4, sum_A3, sum_V4, best_p0, save=False, priors=None):
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

        hexcode = None
        dr_hex = None

        if save == True:
            str_of_prior = ''
            for key in self.prior:
                if key in ['E0', 'log(dE1)', 'log(dE2)', 'log(E_sum)', 'z0', 'z1', 'z2', 'A3_00', 'V4_00', 'A3_01', 'V4_01', 'A3_11', 'V4_11']:
                    str_of_prior += key
                    str_of_prior += str(self.prior[key])
            
            hextype = hashlib.md5()
            str_of_prior = bytes(str_of_prior, encoding= 'utf-8')
            hextype.update(str_of_prior)
            hexcode = hextype.hexdigest()

            g = dict(data=(t_tsep_tau, Amp), prior=priors, params=fit_result.p, chi2=fit_result.chi2, dof=fit_result.dof)

            data_and_results = gv.dumps(g)

            dr_hex = data_and_results.hex()

            print(hexcode)

        return fit_result, hexcode, dr_hex
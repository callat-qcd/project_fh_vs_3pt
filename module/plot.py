# functions named with 'data_plot' plots data to get estimation
# functions named with 'fit_plot' plots fit results on the data

import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv 
class Plot():
    def __init__(self, data_avg_dict, data_avg_dict_completed, pt2_data_range, pt3_data_range, pt3_tau_cut):
        # here the tau_cut is just about plot range
        self.data_avg_dict = data_avg_dict
        self.data_avg_dict_completed = data_avg_dict_completed
        self.pt2_data_range = pt2_data_range
        self.pt3_data_range = pt3_data_range
        self.pt3_tau_cut = pt3_tau_cut
        return

    def data_plot_E0(self, E0):
        '''plot 2pt data to estimate E0 and z0'''
        m0_eff = []
        m0_eff_err = []
    
        for i in range(self.pt2_data_range[0], self.pt2_data_range[1]-1):
            temp = gv.log(self.data_avg_dict['pt2_tsep_'+str(i)] / self.data_avg_dict['pt2_tsep_'+str(i+1)])
            m0_eff.append(temp.mean)
            m0_eff_err.append(temp.sdev)

        plt.figure()
        plt.title('m0_eff')
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), np.array(m0_eff), yerr=np.array(m0_eff_err), fmt='x', label='data')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), (0.665-0.02)*np.ones([len(m0_eff)]), (0.665+0.02)*np.ones([len(m0_eff)]), color='g', alpha = 0.3, label='prior')
        plt.plot(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), 0.665*np.ones([len(m0_eff)]), 'r-')
        plt.xlim([0, 15])
        plt.ylim([0.5, 0.8])
        plt.legend()
        plt.show()

        # plot z0
        z0_mean = []
        z0_err = []

        for i in range(self.pt2_data_range[0], self.pt2_data_range[1]):
            z0_mean.append(gv.sqrt(self.data_avg_dict['pt2_tsep_'+str(i)] * gv.exp(E0 * i)).mean)
            z0_err.append(gv.sqrt(self.data_avg_dict['pt2_tsep_'+str(i)] * gv.exp(E0 * i)).sdev)

        plt.figure()
        plt.title('z0')
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), np.array(z0_mean), yerr=np.array(z0_err), fmt='x', label='data')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), (0.000775-0.0000775)*np.ones([len(z0_mean)]), (0.000775+0.0000775)*np.ones([len(z0_mean)]), color='g', alpha = 0.3, label='prior')
        plt.plot(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), 0.000775*np.ones([len(z0_mean)]), 'r-')
        plt.xlim([0, 15])
        plt.ylim([0.0005, 0.001])
        plt.legend()
        plt.show()
        
    def data_plot_E0_ps(self, E0):
        '''plot 2pt data to estimate z0_ps'''
        m0_ps_eff = []
        m0_ps_eff_err = []
    
        for i in range(self.pt2_data_range[0], self.pt2_data_range[1]-1):
            temp = gv.log(self.data_avg_dict['pt2_ps_tsep_'+str(i)] / self.data_avg_dict['pt2_ps_tsep_'+str(i+1)])
            m0_ps_eff.append(temp.mean)
            m0_ps_eff_err.append(temp.sdev)

        plt.figure()
        plt.title('m0_ps_eff')
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), np.array(m0_ps_eff), yerr=np.array(m0_ps_eff_err), fmt='x', label='data')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), (0.665-0.02)*np.ones([len(m0_ps_eff)]), (0.665+0.02)*np.ones([len(m0_ps_eff)]), color='g', alpha = 0.3, label='prior')
        plt.plot(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), 0.665*np.ones([len(m0_ps_eff)]), 'r-')
        plt.xlim([0, 15])
        plt.ylim([0.5, 0.8])
        plt.legend()
        plt.show()

        # plot z0_ps
        z0_ps_mean = []
        z0_ps_err = []

        for i in range(self.pt2_data_range[0], self.pt2_data_range[1]):
            z0_ps_mean.append(gv.sqrt(self.data_avg_dict['pt2_ps_tsep_'+str(i)] * gv.exp(E0 * i)).mean)
            z0_ps_err.append(gv.sqrt(self.data_avg_dict['pt2_ps_tsep_'+str(i)] * gv.exp(E0 * i)).sdev)

        plt.figure()
        plt.title('z0_ps')
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), np.array(z0_ps_mean), yerr=np.array(z0_ps_err), fmt='x', label='data')
        #plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), (0.000775-0.0000775)*np.ones([len(z0_ps_mean)]), (0.000775+0.0000775)*np.ones([len(z0_ps_mean)]), color='g', alpha = 0.3, label='prior')
        #plt.plot(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), 0.000775*np.ones([len(z0_ps_mean)]), 'r-')
        plt.xlim([0, 15])
        plt.ylim([0, 0.005])
        plt.legend()
        plt.show()

    def data_plot_E1(self, E0, z0, E1):
        '''plot 2pt data to estimate E1 and z1'''
        m1_eff = []
        m1_eff_err = []

        E1_data = np.zeros([self.pt2_data_range[1] - self.pt2_data_range[0]], dtype=gv.GVar)

        for i in range(self.pt2_data_range[0], self.pt2_data_range[1]):
            E1_data[i-self.pt2_data_range[0]] = data_avg_dict['pt2_tsep_'+str(i)] - z0*z0*gv.exp(-E0*i)

        for i in range(self.pt2_data_range[1] - self.pt2_data_range[0] - 1):
            m1_eff.append(gv.log(E1_data[i]/E1_data[i+1]).mean)
            m1_eff_err.append(gv.log(E1_data[i]/E1_data[i+1]).sdev)

        plt.figure()
        plt.title('m1_eff')
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), np.array(m1_eff), yerr=np.array(m1_eff_err), fmt='x', label='data')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), (1.045-0.2)*np.ones([len(m1_eff)]), (1.045+0.2)*np.ones([len(m1_eff)]), color='g', alpha = 0.3, label='prior')
        plt.plot(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), 1.0*np.ones([len(m1_eff)]), 'r-')
        plt.xlim([0, 14])
        plt.ylim([0.5, 2.0])
        plt.legend()
        plt.show()

        # plot z1
        z1_mean = []
        z1_err = []

        for i in range(self.pt2_data_range[1] - self.pt2_data_range[0]):
            temp = gv.sqrt(E1_data[i] * gv.exp(E1 * (i + self.pt2_data_range[0])))
            z1_mean.append(temp.mean)
            z1_err.append(temp.sdev)

        plt.figure()
        plt.title('z1')
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), np.array(z1_mean), yerr=np.array(z1_err), fmt='x', label='data')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), (0.0005-0.0003)*np.ones([len(z1_mean)]), (0.0005+0.0003)*np.ones([len(z1_mean)]), color='g', alpha = 0.3, label='prior')
        plt.plot(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), 0.0005*np.ones([len(z1_mean)]), 'r-')
        plt.xlim([0, 14])
        plt.ylim([-0.00025, 0.001])
        plt.legend(loc='upper right')
        plt.show()

    def data_plot_E2(self, E0, z0, E1, z1, E2):
        '''plot 2pt data to estimate E2 and z2'''
        m2_eff = []
        m2_eff_err = []

        E2_data = np.zeros([self.pt2_data_range[1] - self.pt2_data_range[0]], dtype=gv.GVar)

        for i in range(self.pt2_data_range[0], self.pt2_data_range[1]):
            E2_data[i-self.pt2_data_range[0]] = data_avg_dict['pt2_tsep_'+str(i)] - z0*z0*gv.exp(-E0*i) - z1*z1*gv.exp(-E1*i)

        for i in range(self.pt2_data_range[1] - self.pt2_data_range[0] - 1):
            m2_eff.append(gv.log(E2_data[i]/E2_data[i+1]).mean)
            m2_eff_err.append(gv.log(E2_data[i]/E2_data[i+1]).sdev)

        plt.figure()
        plt.title('m2_eff')
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), np.array(m2_eff), yerr=np.array(m2_eff_err), fmt='x', label='data')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), (1.2-0.2)*np.ones([len(m2_eff)]), (1.2+0.2)*np.ones([len(m2_eff)]), color='g', alpha = 0.3, label='prior')
        plt.plot(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), 1.2*np.ones([len(m2_eff)]), 'r-')
        plt.xlim([0, 14])
        plt.ylim([-0.5, 3.0])
        plt.legend()
        plt.show()

        # plot z2
        z2_mean = []
        z2_err = []

        for i in range(self.pt2_data_range[1] - self.pt2_data_range[0]):
            temp = gv.sqrt(E2_data[i] * gv.exp(E2 * (i + self.pt2_data_range[0])))
            z2_mean.append(temp.mean)
            z2_err.append(temp.sdev)

        plt.figure()
        plt.title('z2')
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), np.array(z2_mean), yerr=np.array(z2_err), fmt='x', label='data')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), (0.0002-0.0001)*np.ones([len(z2_mean)]), (0.0002+0.0001)*np.ones([len(z2_mean)]), color='g', alpha = 0.3, label='prior')
        plt.plot(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), 0.0002*np.ones([len(z2_mean)]), 'r-')
        plt.xlim([0, 14])
        plt.ylim([-0.0005, 0.0015])
        plt.legend(loc='upper right')
        plt.show()

    
    def data_plot_E3(self, E0, z0, E1, z1, E2, z2):
        '''plot 2pt data to estimate E3 and z3'''
        m3_eff = []
        m3_eff_err = []

        E3_data = np.zeros([self.pt2_data_range[1] - self.pt2_data_range[0]], dtype=gv.GVar)

        for i in range(self.pt2_data_range[0], self.pt2_data_range[1]):
            E3_data[i-self.pt2_data_range[0]] = data_avg_dict['pt2_tsep_'+str(i)] - z0*z0*gv.exp(-E0*i) - z1*z1*gv.exp(-E1*i) - z2*z2*gv.exp(-E2*i)

        for i in range(self.pt2_data_range[1] - self.pt2_data_range[0] - 1):
            m3_eff.append(gv.log(E3_data[i]/E3_data[i+1]).mean)
            m3_eff_err.append(gv.log(E3_data[i]/E3_data[i+1]).sdev)

        plt.figure()
        plt.title('m3_eff')
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), np.array(m3_eff), yerr=np.array(m3_eff_err), fmt='x', label='data')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), (1.2-0.2)*np.ones([len(m3_eff)]), (1.2+0.2)*np.ones([len(m3_eff)]), color='g', alpha = 0.3, label='prior')
        plt.plot(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), 1.2*np.ones([len(m3_eff)]), 'r-')
        plt.xlim([0, 14])
        plt.ylim([-0.5, 3.0])
        plt.legend()
        plt.show()

    def data_plot_fh(self, E0_ss=None, E0_ps=None):
        '''plot fh data to estimate d0, maybe dn can also be plotted'''
        d0_gA_ss = []
        d0_gA_ss_err = []
        d0_gA_ps = []
        d0_gA_ps_err = []
        d0_gV_ss = []
        d0_gV_ss_err = []
        d0_gV_ps = []
        d0_gV_ps_err = []
        
        if E0_ss == None:
            for i in range(self.pt2_data_range[0], self.pt2_data_range[1]-1):
                temp1 = (self.data_avg_dict['pt2_tsep_'+str(i)] / self.data_avg_dict['pt2_tsep_'+str(i+1)]) * self.data_avg_dict['fh_gA_ss_tsep_1']
                d0_gA_ss.append(temp1.mean)
                d0_gA_ss_err.append(temp1.sdev)

                temp3 = (self.data_avg_dict['pt2_tsep_'+str(i)] / self.data_avg_dict['pt2_tsep_'+str(i+1)]) * self.data_avg_dict['fh_gV_ss_tsep_1']
                d0_gV_ss.append(temp3.mean)
                d0_gV_ss_err.append(temp3.sdev)


        else:
            for i in range(self.pt2_data_range[0], self.pt2_data_range[1]-1):
                temp1 = np.exp(E0_ss) * self.data_avg_dict['fh_gA_ss_tsep_1']
                d0_gA_ss.append(temp1.mean)
                d0_gA_ss_err.append(temp1.sdev)

                temp3 = np.exp(E0_ss) * self.data_avg_dict['fh_gV_ss_tsep_1']
                d0_gV_ss.append(temp3.mean)
                d0_gV_ss_err.append(temp3.sdev)

                
        if E0_ps == None:
            for i in range(self.pt2_data_range[0], self.pt2_data_range[1]-1):
                temp2 = (self.data_avg_dict['pt2_ps_tsep_'+str(i)] / self.data_avg_dict['pt2_ps_tsep_'+str(i+1)]) * self.data_avg_dict['fh_gA_ps_tsep_1']
                d0_gA_ps.append(temp2.mean)
                d0_gA_ps_err.append(temp2.sdev)

                temp4 = (self.data_avg_dict['pt2_ps_tsep_'+str(i)] / self.data_avg_dict['pt2_ps_tsep_'+str(i+1)]) * self.data_avg_dict['fh_gV_ps_tsep_1']
                d0_gV_ps.append(temp4.mean)
                d0_gV_ps_err.append(temp4.sdev)
                
        else:
            for i in range(self.pt2_data_range[0], self.pt2_data_range[1]-1):
                temp2 = np.exp(E0_ps) * self.data_avg_dict['fh_gA_ps_tsep_1']
                d0_gA_ps.append(temp2.mean)
                d0_gA_ps_err.append(temp2.sdev)

                temp4 = np.exp(E0_ps) * self.data_avg_dict['fh_gV_ps_tsep_1']
                d0_gV_ps.append(temp4.mean)
                d0_gV_ps_err.append(temp4.sdev)
                
                
        plt.figure()
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), d0_gA_ss, yerr=d0_gA_ss_err, fmt='x', label='gA_ss')
        plt.xlim([0, 25])
        plt.ylim([-0.000005, 0])
        plt.legend()
        plt.title('d0_gA_ss')
        plt.show()
        
        plt.figure()
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), d0_gA_ps, yerr=d0_gA_ps_err, fmt='x', label='gA_ps')
        plt.xlim([0, 25])
        plt.ylim([-0.00002, 0])
        plt.legend()
        plt.title('d0_gA_ps')
        plt.show()
        
        plt.figure()
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), d0_gV_ss, yerr=d0_gV_ss_err, fmt='x', label='gV_ss')
        plt.xlim([0, 25])
        plt.ylim([0, 0.000003])
        plt.legend()
        plt.title('d0_gV_ss')
        plt.show()
        
        plt.figure()
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), d0_gV_ps, yerr=d0_gV_ps_err, fmt='x', label='gV_ps')
        plt.xlim([0, 25])
        plt.ylim([0, 0.00002])
        plt.legend()
        plt.title('d0_gV_ps')
        plt.show()
        
            
    def fit_plot_E0(self, best_fitter, best_fit):
        '''plot E0 with fit over data'''
        m0_eff = []
        m0_eff_err = []
        m0_eff_fit = []
        m0_eff_fit_err = []

        pt2_fitter = best_fitter.pt2_fit_function(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), best_fit.p)['pt2']
    
        for i in range(self.pt2_data_range[0], self.pt2_data_range[1]-1):
            temp = gv.log(self.data_avg_dict['pt2_tsep_'+str(i)] / self.data_avg_dict['pt2_tsep_'+str(i+1)])
            m0_eff.append(temp.mean)
            m0_eff_err.append(temp.sdev)

        for i in range(self.pt2_data_range[1] - self.pt2_data_range[0]-1):
            temp = gv.log(pt2_fitter[i] / pt2_fitter[i+1])
            m0_eff_fit.append(temp.mean)
            m0_eff_fit_err.append(temp.sdev)

        m0_eff_y1 = []
        m0_eff_y2 = []

        for i in range(self.pt2_data_range[1] - self.pt2_data_range[0]-1):
            m0_eff_y1.append(m0_eff_fit[i] + m0_eff_fit_err[i])
            m0_eff_y2.append(m0_eff_fit[i] - m0_eff_fit_err[i])

        plt.figure()
        plt.title('m0_eff')
        plt.xlabel('t')
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), np.array(m0_eff), yerr=np.array(m0_eff_err), fmt='x', label='data')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), np.array(m0_eff_y1), np.array(m0_eff_y2), color='g', alpha = 0.3, label='fit')
        plt.xlim([3, 25])
        plt.ylim([0.5, 0.8])
        plt.legend()
        plt.show()
        
    def fit_plot_E0_ps(self, best_fitter, best_fit):
        '''plot E0_ps with fit over data'''
        m0_eff = []
        m0_eff_err = []
        m0_eff_fit = []
        m0_eff_fit_err = []

        pt2_fitter = best_fitter.pt2_ps_fit_function(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), best_fit.p)['pt2_ps']
    
        for i in range(self.pt2_data_range[0], self.pt2_data_range[1]-1):
            temp = gv.log(self.data_avg_dict['pt2_ps_tsep_'+str(i)] / self.data_avg_dict['pt2_ps_tsep_'+str(i+1)])
            m0_eff.append(temp.mean)
            m0_eff_err.append(temp.sdev)

        for i in range(self.pt2_data_range[1] - self.pt2_data_range[0]-1):
            temp = gv.log(pt2_fitter[i] / pt2_fitter[i+1])
            m0_eff_fit.append(temp.mean)
            m0_eff_fit_err.append(temp.sdev)

        m0_eff_y1 = []
        m0_eff_y2 = []

        for i in range(self.pt2_data_range[1] - self.pt2_data_range[0]-1):
            m0_eff_y1.append(m0_eff_fit[i] + m0_eff_fit_err[i])
            m0_eff_y2.append(m0_eff_fit[i] - m0_eff_fit_err[i])

        plt.figure()
        plt.title('m0_ps_eff')
        plt.xlabel('t')
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), np.array(m0_eff), yerr=np.array(m0_eff_err), fmt='x', label='data')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), np.array(m0_eff_y1), np.array(m0_eff_y2), color='g', alpha = 0.3, label='fit')
        plt.xlim([3, 25])
        plt.ylim([0.5, 0.8])
        plt.legend()
        plt.show()

    def fit_plot_pt3(self, best_fitter, best_fit):
        '''plot form factors with fit over pt3_data'''
        A3_data = {}
        A3_data_err = {}
        A3_fit = {}
        A3_fit_err = {}

        V4_data = {}
        V4_data_err = {}
        V4_fit = {}
        V4_fit_err = {}

        pt3_tsep_A3 = []
        pt3_tau_A3 = []
        pt3_tsep_V4 = []
        pt3_tau_V4 = []

        for i in range(self.pt3_data_range[0]+1, self.pt3_data_range[1]-1):
            if i-self.pt3_tau_cut+1 <= self.pt3_tau_cut:
                print('cut error 1')
            for j in range(self.pt3_tau_cut, i-self.pt3_tau_cut+1):
                pt3_tsep_A3.append(i)
                pt3_tau_A3.append(j)

        for i in range(self.pt3_data_range[0]+1, self.pt3_data_range[1]-1):
            if i-self.pt3_tau_cut+1 <= self.pt3_tau_cut:
                print('cut error 2')
            for j in range(self.pt3_tau_cut, i-self.pt3_tau_cut+1):   
                pt3_tsep_V4.append(i)
                pt3_tau_V4.append(j)

        pt2_fitter = best_fitter.pt2_fit_function(np.arange(self.pt3_data_range[0], self.pt3_data_range[1]), best_fit.p)['pt2']

        pt3_A3_fitter = best_fitter.pt3_fit_function(np.array(pt3_tsep_A3), np.array(pt3_tsep_V4), np.array(pt3_tau_A3), np.array(pt3_tau_V4), best_fit.p)['pt3_A3']

        pt3_V4_fitter = best_fitter.pt3_fit_function(np.array(pt3_tsep_A3), np.array(pt3_tsep_V4), np.array(pt3_tau_A3), np.array(pt3_tau_V4), best_fit.p)['pt3_V4']

        for i in range(self.pt3_data_range[0]+1, self.pt3_data_range[1]-1):
            A3_data['tsep_'+str(i)] = []
            A3_data_err['tsep_'+str(i)] = []
            A3_fit['tsep_'+str(i)] = []
            A3_fit_err['tsep_'+str(i)] = []
            V4_data['tsep_'+str(i)] = []
            V4_data_err['tsep_'+str(i)] = []
            V4_fit['tsep_'+str(i)] = []
            V4_fit_err['tsep_'+str(i)] = []

        index_A3 = 0
        for i in range(self.pt3_data_range[0]+1, self.pt3_data_range[1]-1): #cut tsep=3 and 12
            if i-self.pt3_tau_cut+1 <= self.pt3_tau_cut:
                print('cut error 3')
            for j in range(self.pt3_tau_cut, i-self.pt3_tau_cut+1):
                temp = self.data_avg_dict['pt3_A3_tsep_'+str(i)][j] / self.data_avg_dict['pt2_tsep_'+str(i)]
                A3_data['tsep_'+str(i)].append(temp.mean)
                A3_data_err['tsep_'+str(i)].append(temp.sdev)
                
                index_A3 = int((i+1)*(i-4)/2 + j -self.pt3_tau_cut) # from t=4, tau_cut=1
                #index_A3 = int((i-3)*(i-4)/2 + j -self.pt3_tau_cut) # from t=4, tau_cut=2
                temp = pt3_A3_fitter[index_A3] / pt2_fitter[i - self.pt3_data_range[0]]
                A3_fit['tsep_'+str(i)].append(temp.mean)
                A3_fit_err['tsep_'+str(i)].append(temp.sdev)

        index_V4 = 0
        for i in range(self.pt3_data_range[0]+1, self.pt3_data_range[1]-1): #cut tsep=3 and 12
            if i-self.pt3_tau_cut+1 <= self.pt3_tau_cut:
                print('cut error 4')
            for j in range(self.pt3_tau_cut, i-self.pt3_tau_cut+1):
                temp = self.data_avg_dict['pt3_V4_tsep_'+str(i)][j] / self.data_avg_dict['pt2_tsep_'+str(i)]
                V4_data['tsep_'+str(i)].append(temp.mean)
                V4_data_err['tsep_'+str(i)].append(temp.sdev)

                index_V4 = int((i+1)*(i-4)/2 + j -self.pt3_tau_cut) # from t=4, tau_cut=1
                #index_V4 = int((i-3)*(i-4)/2 + j -self.pt3_tau_cut) # from t=4, tau_cut=2
                temp = pt3_V4_fitter[index_V4] / pt2_fitter[i - self.pt3_data_range[0]]
                V4_fit['tsep_'+str(i)].append(temp.mean)
                V4_fit_err['tsep_'+str(i)].append(temp.sdev)

        plt.figure()
        for i in range(self.pt3_data_range[0]+1, self.pt3_data_range[1]-1): #cut tsep=3 and 12
            plt.errorbar(np.arange(-(i-2)/2, (i)/2, 1), np.array(A3_data['tsep_' + str(i)]), yerr=np.array(A3_data_err['tsep_' + str(i)]), fmt='x', ecolor='r', label = 'data')# tau cut = 1
            
            A3_fit_y1 = np.array(A3_fit['tsep_'+str(i)]) + np.array(A3_fit_err['tsep_'+str(i)])
            A3_fit_y2 = np.array(A3_fit['tsep_'+str(i)]) - np.array(A3_fit_err['tsep_'+str(i)])
            plt.fill_between(np.arange(-(i-2)/2, (i)/2, 1), A3_fit_y1, A3_fit_y2, color='g', alpha=0.3, label = 'fit') # tau_cut=1
        plt.title('pt3_A300')
        plt.ylim([1.1, 1.35])
        plt.xlabel('centered tau')
        #plt.legend()
        plt.show()


        plt.figure()
        for i in range(self.pt3_data_range[0]+1, self.pt3_data_range[1]-1): #cut tsep=3 and 12
            plt.errorbar(np.arange(-(i-2)/2, (i)/2, 1), np.array(V4_data['tsep_' + str(i)]), yerr=np.array(V4_data_err['tsep_' + str(i)]), fmt='x', ecolor='r', label = 'data')# tau cut = 1
            
            V4_fit_y1 = np.array(V4_fit['tsep_'+str(i)]) + np.array(V4_fit_err['tsep_'+str(i)])
            V4_fit_y2 = np.array(V4_fit['tsep_'+str(i)]) - np.array(V4_fit_err['tsep_'+str(i)])
            plt.fill_between(np.arange(-(i-2)/2, (i)/2, 1), V4_fit_y1, V4_fit_y2, color='g', alpha=0.3, label = 'fit') # tau_cut=1
        plt.title('pt3_V400')
        plt.ylim([0.95, 1.15])
        plt.xlabel('centered tau')
        #plt.legend()
        plt.show()

    def fit_plot_sum(self, best_fitter, best_fit):
        '''plot form factors with fit over sum_data'''
        sum_A3_data = []
        sum_A3_data_err = []
        sum_A3_fit = []
        sum_A3_fit_err = []

        sum_V4_data = []
        sum_V4_data_err = []
        sum_V4_fit = []
        sum_V4_fit_err = []

        pt2_fitter = best_fitter.pt2_fit_function(np.arange(self.pt3_data_range[0], self.pt3_data_range[1]), best_fit.p)['pt2']

        sum_A3_fitter = best_fitter.summation(np.arange(self.pt3_data_range[0], self.pt3_data_range[1]), np.arange(self.pt3_data_range[0], self.pt3_data_range[1]), best_fit.p)['sum_A3']

        sum_V4_fitter = best_fitter.summation(np.arange(self.pt3_data_range[0], self.pt3_data_range[1]), np.arange(self.pt3_data_range[0], self.pt3_data_range[1]), best_fit.p)['sum_V4']

        for i in range(self.pt3_data_range[1] - self.pt3_data_range[0] - 1):
            temp = sum_A3_fitter[i+1] / pt2_fitter[i+1] - sum_A3_fitter[i] / pt2_fitter[i]
            sum_A3_fit.append(temp.mean)
            sum_A3_fit_err.append(temp.sdev)

            temp = sum_V4_fitter[i+1] / pt2_fitter[i+1] - sum_V4_fitter[i] / pt2_fitter[i]
            sum_V4_fit.append(temp.mean)
            sum_V4_fit_err.append(temp.sdev)

        for i in range(self.pt3_data_range[0], self.pt3_data_range[1]-1):
            temp = self.data_avg_dict_completed['sum_A3_fit_'+str(i)]
            sum_A3_data.append(temp.mean)
            sum_A3_data_err.append(temp.sdev)

            temp = self.data_avg_dict_completed['sum_V4_fit_'+str(i)]
            sum_V4_data.append(temp.mean)
            sum_V4_data_err.append(temp.sdev)

        sum_A3_fit_y1 = np.array(sum_A3_fit) + np.array(sum_A3_fit_err)
        sum_A3_fit_y2 = np.array(sum_A3_fit) - np.array(sum_A3_fit_err)

        plt.figure()
        plt.title('sum_A3_00')
        plt.xlabel('tsep')
        plt.fill_between(np.arange(self.pt3_data_range[0], self.pt3_data_range[1]-1), sum_A3_fit_y1, sum_A3_fit_y2, color='g', alpha=0.3, label='fit')
        
        plt.errorbar(np.arange(self.pt3_data_range[0], self.pt3_data_range[1]-1),np.array(sum_A3_data), yerr= np.array(sum_A3_data_err), fmt='x', ecolor='r', label='data')
        plt.ylim([1.1, 1.35])
        plt.legend()
        plt.show()

        sum_V4_fit_y1 = np.array(sum_V4_fit) + np.array(sum_V4_fit_err)
        sum_V4_fit_y2 = np.array(sum_V4_fit) - np.array(sum_V4_fit_err)

        plt.figure()
        plt.title('sum_V4_00')
        plt.xlabel('tsep')
        plt.fill_between(np.arange(self.pt3_data_range[0], self.pt3_data_range[1]-1), sum_V4_fit_y1, sum_V4_fit_y2, color='g', alpha=0.3, label='fit')
        
        plt.errorbar(np.arange(self.pt3_data_range[0], self.pt3_data_range[1]-1),np.array(sum_V4_data), yerr= np.array(sum_V4_data_err), fmt='x', ecolor='r', label='data')
        plt.ylim([0.95, 1.15])
        plt.legend()
        plt.show()
        
    def fit_plot_fh_ss(self, best_fitter, best_fit):
        '''plot form factors with fit over sum_data'''
        fh_ss_A3_data = []
        fh_ss_A3_data_err = []
        fh_ss_A3_fit = []
        fh_ss_A3_fit_err = []

        fh_ss_V4_data = []
        fh_ss_V4_data_err = []
        fh_ss_V4_fit = []
        fh_ss_V4_fit_err = []

        pt2_fitter = best_fitter.pt2_fit_function(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), best_fit.p)['pt2']

        fh_ss_A3_fitter = best_fitter.fh_ss_fit_function(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), best_fit.p)['fh_ss_A3']

        fh_ss_V4_fitter = best_fitter.fh_ss_fit_function(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), best_fit.p)['fh_ss_V4']

        for i in range(self.pt2_data_range[1] - self.pt2_data_range[0] - 1):
            temp = fh_ss_A3_fitter[i+1] / pt2_fitter[i+1] - fh_ss_A3_fitter[i] / pt2_fitter[i]
            fh_ss_A3_fit.append(temp.mean)
            fh_ss_A3_fit_err.append(temp.sdev)

            temp = fh_ss_V4_fitter[i+1] / pt2_fitter[i+1] - fh_ss_V4_fitter[i] / pt2_fitter[i]
            fh_ss_V4_fit.append(temp.mean)
            fh_ss_V4_fit_err.append(temp.sdev)

        for i in range(self.pt2_data_range[0], self.pt2_data_range[1]-1):
            temp = self.data_avg_dict_completed['fh_ss_A3_fit_'+str(i)]
            fh_ss_A3_data.append(temp.mean)
            fh_ss_A3_data_err.append(temp.sdev)

            temp = self.data_avg_dict_completed['fh_ss_V4_fit_'+str(i)]
            fh_ss_V4_data.append(temp.mean)
            fh_ss_V4_data_err.append(temp.sdev)

        fh_ss_A3_fit_y1 = np.array(fh_ss_A3_fit) + np.array(fh_ss_A3_fit_err)
        fh_ss_A3_fit_y2 = np.array(fh_ss_A3_fit) - np.array(fh_ss_A3_fit_err)

        plt.figure()
        plt.title('fh_ss_A3_00')
        plt.xlabel('tsep')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), fh_ss_A3_fit_y1, fh_ss_A3_fit_y2, color='g', alpha=0.3, label='fit')
        
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1),np.array(fh_ss_A3_data), yerr= np.array(fh_ss_A3_data_err), fmt='x', ecolor='r', label='data')
        plt.ylim([1.1, 1.35])
        plt.xlim([0, 20])
        plt.legend()
        plt.show()

        fh_ss_V4_fit_y1 = np.array(fh_ss_V4_fit) + np.array(fh_ss_V4_fit_err)
        fh_ss_V4_fit_y2 = np.array(fh_ss_V4_fit) - np.array(fh_ss_V4_fit_err)

        plt.figure()
        plt.title('fh_ss_V4_00')
        plt.xlabel('tsep')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), fh_ss_V4_fit_y1, fh_ss_V4_fit_y2, color='g', alpha=0.3, label='fit')
        
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1),np.array(fh_ss_V4_data), yerr= np.array(fh_ss_V4_data_err), fmt='x', ecolor='r', label='data')
        plt.ylim([0.95, 1.15])
        plt.xlim([0, 20])
        plt.legend()
        plt.show()
        
    def fit_plot_fh_ps(self, best_fitter, best_fit):
        '''plot form factors with fit over sum_data'''
        fh_ps_A3_data = []
        fh_ps_A3_data_err = []
        fh_ps_A3_fit = []
        fh_ps_A3_fit_err = []

        fh_ps_V4_data = []
        fh_ps_V4_data_err = []
        fh_ps_V4_fit = []
        fh_ps_V4_fit_err = []

        pt2_ps_fitter = best_fitter.pt2_ps_fit_function(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), best_fit.p)['pt2_ps']

        fh_ps_A3_fitter = best_fitter.fh_ps_fit_function(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), best_fit.p)['fh_ps_A3']

        fh_ps_V4_fitter = best_fitter.fh_ps_fit_function(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), np.arange(self.pt2_data_range[0], self.pt2_data_range[1]), best_fit.p)['fh_ps_V4']

        for i in range(self.pt2_data_range[1] - self.pt2_data_range[0] - 1):
            temp = fh_ps_A3_fitter[i+1] / pt2_ps_fitter[i+1] - fh_ps_A3_fitter[i] / pt2_ps_fitter[i]
            fh_ps_A3_fit.append(temp.mean)
            fh_ps_A3_fit_err.append(temp.sdev)

            temp = fh_ps_V4_fitter[i+1] / pt2_ps_fitter[i+1] - fh_ps_V4_fitter[i] / pt2_ps_fitter[i]
            fh_ps_V4_fit.append(temp.mean)
            fh_ps_V4_fit_err.append(temp.sdev)

        for i in range(self.pt2_data_range[0], self.pt2_data_range[1]-1):
            temp = self.data_avg_dict_completed['fh_ps_A3_fit_'+str(i)]
            fh_ps_A3_data.append(temp.mean)
            fh_ps_A3_data_err.append(temp.sdev)

            temp = self.data_avg_dict_completed['fh_ps_V4_fit_'+str(i)]
            fh_ps_V4_data.append(temp.mean)
            fh_ps_V4_data_err.append(temp.sdev)

        fh_ps_A3_fit_y1 = np.array(fh_ps_A3_fit) + np.array(fh_ps_A3_fit_err)
        fh_ps_A3_fit_y2 = np.array(fh_ps_A3_fit) - np.array(fh_ps_A3_fit_err)

        plt.figure()
        plt.title('fh_ps_A3_00')
        plt.xlabel('tsep')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), fh_ps_A3_fit_y1, fh_ps_A3_fit_y2, color='g', alpha=0.3, label='fit')
        
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1),np.array(fh_ps_A3_data), yerr= np.array(fh_ps_A3_data_err), fmt='x', ecolor='r', label='data')
        plt.ylim([1.1, 1.35])
        plt.xlim([0, 20])
        plt.legend()
        plt.show()

        fh_ps_V4_fit_y1 = np.array(fh_ps_V4_fit) + np.array(fh_ps_V4_fit_err)
        fh_ps_V4_fit_y2 = np.array(fh_ps_V4_fit) - np.array(fh_ps_V4_fit_err)

        plt.figure()
        plt.title('fh_ps_V4_00')
        plt.xlabel('tsep')
        plt.fill_between(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1), fh_ps_V4_fit_y1, fh_ps_V4_fit_y2, color='g', alpha=0.3, label='fit')
        
        plt.errorbar(np.arange(self.pt2_data_range[0], self.pt2_data_range[1]-1),np.array(fh_ps_V4_data), yerr= np.array(fh_ps_V4_data_err), fmt='x', ecolor='r', label='data')
        plt.ylim([0.95, 1.15])
        plt.xlim([0, 20])
        plt.legend()
        plt.show()
   
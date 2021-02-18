import h5py as h5 
import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  
import lsqfit as lsf
import math

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

        start = 0 ### used to adjust stat range for half stat
        end = 784

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
            ## full tau: shape = (784, 96)

        print(np.shape(data_dict['pt2_tsep_1']))
        print(np.shape(data_dict['pt3_A3_tsep_3']))

        data_avg_dict = gv.dataset.avg_data(data_dict)

        return data_avg_dict

    def add_sum_data(self, data_avg_dict, sum_tau_cut):
        if sum_tau_cut == -1: # sum full tau
            for i in range(self.pt3_data_range[0], self.pt3_data_range[1]):
                data_avg_dict['sum_A3_tsep_' + str(i)] = gv.gvar(0, 0)
                data_avg_dict['sum_V4_tsep_' + str(i)] = gv.gvar(0, 0)

                for j in range(len(data_avg_dict['pt3_A3_tsep_3'])): # do the summation
                    data_avg_dict['sum_A3_tsep_' + str(i)] += data_avg_dict['pt3_A3_tsep_' + str(i)][j]
                    data_avg_dict['sum_V4_tsep_' + str(i)] += data_avg_dict['pt3_V4_tsep_' + str(i)][j]

            for i in range(self.pt3_data_range[0], self.pt3_data_range[1]-1): # use ratio form to do the fit
                # fit_14 will use tsep_15 here, so the tsep range of sum part should not be bigger than 13
                data_avg_dict['sum_A3_fit_'+str(i)] = (data_avg_dict['sum_A3_tsep_' + str(i+1)]/data_avg_dict['pt2_tsep_'+str(i+1)]) - (data_avg_dict['sum_A3_tsep_' + str(i)]/data_avg_dict['pt2_tsep_'+str(i)])

                data_avg_dict['sum_V4_fit_'+str(i)] = (data_avg_dict['sum_V4_tsep_' + str(i+1)]/data_avg_dict['pt2_tsep_'+str(i+1)]) - (data_avg_dict['sum_V4_tsep_' + str(i)]/data_avg_dict['pt2_tsep_'+str(i)])

            return data_avg_dict


        else:
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

    def add_sum_data_scattered(self, data_avg_dict, sum_tau_cut, sum_tsep):
        sum_tsep = [t for t in sum_tsep]
        sum_tsep.append(2*sum_tsep[-1] - sum_tsep[-2]) # the last element 
        if sum_tau_cut == -1: # sum full tau
            for i in range(self.pt3_data_range[0], self.pt3_data_range[1]):
                data_avg_dict['sum_A3_tsep_' + str(i)] = gv.gvar(0, 0)
                data_avg_dict['sum_V4_tsep_' + str(i)] = gv.gvar(0, 0)

                for j in range(len(data_avg_dict['pt3_A3_tsep_3'])): # do the summation
                    data_avg_dict['sum_A3_tsep_' + str(i)] += data_avg_dict['pt3_A3_tsep_' + str(i)][j]
                    data_avg_dict['sum_V4_tsep_' + str(i)] += data_avg_dict['pt3_V4_tsep_' + str(i)][j]

            gap = sum_tsep[1] - sum_tsep[0]
            for t in range(self.pt3_data_range[0], self.pt3_data_range[1]-gap):
                data_avg_dict['sum_A3_fit_'+str(t)] = (
                    (data_avg_dict['sum_A3_tsep_' + str(t+gap)]/data_avg_dict['pt2_tsep_'+str(t+gap)]) - (data_avg_dict['sum_A3_tsep_' + str(t)]/data_avg_dict['pt2_tsep_'+str(t)])
                    ) / gap

                data_avg_dict['sum_V4_fit_'+str(t)] = (
                    (data_avg_dict['sum_V4_tsep_' + str(t+gap)]/data_avg_dict['pt2_tsep_'+str(t+gap)]) - (data_avg_dict['sum_V4_tsep_' + str(t)]/data_avg_dict['pt2_tsep_'+str(t)])
                    ) / gap

            return data_avg_dict


        else:
            for i in range(self.pt3_data_range[0], self.pt3_data_range[1]):
                data_avg_dict['sum_A3_tsep_' + str(i)] = gv.gvar(0, 0)
                data_avg_dict['sum_V4_tsep_' + str(i)] = gv.gvar(0, 0)

                for j in range(sum_tau_cut, i-sum_tau_cut+1): # do the summation
                    data_avg_dict['sum_A3_tsep_' + str(i)] += data_avg_dict['pt3_A3_tsep_' + str(i)][j]
                    data_avg_dict['sum_V4_tsep_' + str(i)] += data_avg_dict['pt3_V4_tsep_' + str(i)][j]

            gap = sum_tsep[1] - sum_tsep[0]
            for t in range(self.pt3_data_range[0], self.pt3_data_range[1]-gap):
                data_avg_dict['sum_A3_fit_'+str(t)] = (
                    (data_avg_dict['sum_A3_tsep_' + str(t+gap)]/data_avg_dict['pt2_tsep_'+str(t+gap)]) - (data_avg_dict['sum_A3_tsep_' + str(t)]/data_avg_dict['pt2_tsep_'+str(t)])
                    ) / gap

                data_avg_dict['sum_V4_fit_'+str(t)] = (
                    (data_avg_dict['sum_V4_tsep_' + str(t+gap)]/data_avg_dict['pt2_tsep_'+str(t+gap)]) - (data_avg_dict['sum_V4_tsep_' + str(t)]/data_avg_dict['pt2_tsep_'+str(t)])
                    ) / gap


            return data_avg_dict
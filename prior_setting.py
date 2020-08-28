import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  


def prior_con(pt2_nstates, pt3_nstates, sum_nstates):
    prior = {}

    prior['E0'] = gv.gvar(0.50, 0.02)
    prior['z0'] = gv.gvar(0.00034, 0.00034)

    prior['log(dE1)'] = gv.gvar(-1, 0.69)
    prior['z1'] = gv.gvar(0, 0.00025)

    prior['log(dE2)'] = gv.gvar(-1, 0.69)
    prior['z2'] = gv.gvar(0, 0.00015)

    max_nstate = max(np.array([pt2_nstates, pt3_nstates, sum_nstates]))

    if max_nstate > 3:
        for i in range(3, max_nstate):
            prior['log(dE'+str(i)+')'] = gv.gvar(-1, 0.69)
            prior['z'+str(i)] = gv.gvar(0, 0.00015)

    # same garbage for sum and fh, so below are shared 
    prior['log(E_sum)'] = gv.gvar(-1, 0.69) 
    prior['z_sum'] = gv.gvar(0, 0.00015) 

    prior['A3_00'] = gv.gvar(1.25, 0.05)
    prior['V4_00'] = gv.gvar(1.0, 0.2)

    ff_nstates = max(np.array([pt3_nstates, sum_nstates]))

    for i in range(ff_nstates):
        for j in range(ff_nstates):
            if i+j >= 1:
                if j < i:  
                    prior['A3_'+str(j)+str(i)] = gv.gvar(0, 1)
                    prior['V4_'+str(j)+str(i)] = gv.gvar(0, 1)

                elif j == i:
                    prior['A3_'+str(j)+str(i)] = gv.gvar(0, 1)
                    prior['V4_'+str(j)+str(i)] = gv.gvar(1, 0.2)

    for i in range(sum_nstates-1):
        prior['sum_A3_'+str(i)] = gv.gvar(0, 1)
        prior['sum_V4_'+str(i)] = gv.gvar(0, 1)

    prior['sum_A3_'+str(sum_nstates-1)] = gv.gvar(0, 1)
    prior['sum_V4_'+str(sum_nstates-1)] = gv.gvar(1, 0.2)

    return prior

def prior_bigger(max_nstate, pt3_nstates, sum_nstates, fh_nstates):
    prior = {}

    prior['E0'] = gv.gvar(0.50, 0.02)
    prior['z0'] = gv.gvar(0.00034, 0.00034)

    prior['log(dE1)'] = gv.gvar(-1, 0.69)
    prior['z1'] = gv.gvar(0, 0.00025)

    prior['log(dE2)'] = gv.gvar(-1 + np.log(2), 0.69)
    prior['z2'] = gv.gvar(0, 0.00015)

    max_nstate = max(np.array([pt2_nstates, pt3_nstates, sum_nstates]))

    if max_nstate > 3:
        for i in range(3, max_nstate):
            prior['log(dE'+str(i)+')'] = gv.gvar(-1 + np.log(i), 0.69)
            prior['z'+str(i)] = gv.gvar(0, 0.00015)

    # same garbage for sum and fh, so below are shared 
    prior['log(E_sum)'] = gv.gvar(-1 + np.log(sum_nstates-1), 0.69) 
    prior['z_sum'] = gv.gvar(0, 0.00015) 

    prior['A3_00'] = gv.gvar(1.25, 0.05)
    prior['V4_00'] = gv.gvar(1.0, 0.2)

    ff_nstates = max(np.array([pt3_nstates, sum_nstates]))

    for i in range(ff_nstates):
        for j in range(ff_nstates):
            if i+j >= 1:
                if j < i:  
                    prior['A3_'+str(j)+str(i)] = gv.gvar(0, 1)
                    prior['V4_'+str(j)+str(i)] = gv.gvar(0, 1)

                elif j == i:
                    prior['A3_'+str(j)+str(i)] = gv.gvar(0, 1)
                    prior['V4_'+str(j)+str(i)] = gv.gvar(1, 0.2)

    for i in range(sum_nstates-1):
        prior['sum_A3_'+str(i)] = gv.gvar(0, 1)
        prior['sum_V4_'+str(i)] = gv.gvar(0, 1)

    prior['sum_A3_'+str(sum_nstates-1)] = gv.gvar(0, 1)
    prior['sum_V4_'+str(sum_nstates-1)] = gv.gvar(1, 0.2)

    return prior


def prior_smaller(max_nstate, pt3_nstates, sum_nstates, fh_nstates):
    prior = {}

    prior['E0'] = gv.gvar(0.50, 0.02)
    prior['z0'] = gv.gvar(0.00034, 0.00034)

    prior['log(dE1)'] = gv.gvar(-1, 0.69)
    prior['z1'] = gv.gvar(0, 0.00025)

    prior['log(dE2)'] = gv.gvar(-1 - np.log(2), 0.69)
    prior['z2'] = gv.gvar(0, 0.00015)

    max_nstate = max(np.array([pt2_nstates, pt3_nstates, sum_nstates]))

    if max_nstate > 3:
        for i in range(3, max_nstate):
            prior['log(dE'+str(i)+')'] = gv.gvar(-1 - np.log(i), 0.69)
            prior['z'+str(i)] = gv.gvar(0, 0.00015)

    # same garbage for sum and fh, so below are shared 
    prior['log(E_sum)'] = gv.gvar(-1 - np.log(sum_nstates-1), 0.69) 
    prior['z_sum'] = gv.gvar(0, 0.00015) 

    prior['A3_00'] = gv.gvar(1.25, 0.05)
    prior['V4_00'] = gv.gvar(1.0, 0.2)

    ff_nstates = max(np.array([pt3_nstates, sum_nstates]))

    for i in range(ff_nstates):
        for j in range(ff_nstates):
            if i+j >= 1:
                if j < i:  
                    prior['A3_'+str(j)+str(i)] = gv.gvar(0, 1)
                    prior['V4_'+str(j)+str(i)] = gv.gvar(0, 1)

                elif j == i:
                    prior['A3_'+str(j)+str(i)] = gv.gvar(0, 1)
                    prior['V4_'+str(j)+str(i)] = gv.gvar(1, 0.2)

    for i in range(sum_nstates-1):
        prior['sum_A3_'+str(i)] = gv.gvar(0, 1)
        prior['sum_V4_'+str(i)] = gv.gvar(0, 1)

    prior['sum_A3_'+str(sum_nstates-1)] = gv.gvar(0, 1)
    prior['sum_V4_'+str(sum_nstates-1)] = gv.gvar(1, 0.2)

    return prior
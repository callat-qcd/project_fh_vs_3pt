import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  

def prior_ho_width_1(pt2_nstates, pt3_nstates, sum_nstates):
    prior = gv.BufferDict()

    prior['E0'] = gv.gvar(0.50, 0.05)
    prior['z0'] = gv.gvar(0.00034, 0.00034)

    prior['log(dE1)'] = gv.gvar(-1.25, 0.5)
    prior['z1'] = gv.gvar(0, 0.00025)

    prior['log(dE2)'] = gv.gvar(-1.25, 0.5)
    prior['z2'] = gv.gvar(0, 0.00015)

    max_nstate = max(np.array([pt2_nstates, pt3_nstates, sum_nstates]))

    if max_nstate > 3:
        for i in range(3, max_nstate):
            prior['log(dE'+str(i)+')'] = gv.gvar(-1.25, 0.5)
            prior['z'+str(i)] = gv.gvar(0, 0.00015)

    prior['log(dE4)'] = gv.gvar(-1.25, 0.5*5)

    prior['log(E_sum)'] = gv.gvar(-1.25, 0.5) 
    prior['z_sum'] = gv.gvar(0, 0.00015) 

    prior['A3_00'] = gv.gvar(1.2, 0.2)
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

def prior_ho_width_gs_times(pt2_nstates, pt3_nstates, sum_nstates, times): # change prior width of g.s.
    prior = gv.BufferDict()

    prior['E0'] = gv.gvar(0.50, 0.05*times)
    prior['z0'] = gv.gvar(0.00034, 0.00034*times)

    prior['log(dE1)'] = gv.gvar(-1.25, 0.5)
    prior['z1'] = gv.gvar(0, 0.00025)

    prior['log(dE2)'] = gv.gvar(-1.25, 0.5)
    prior['z2'] = gv.gvar(0, 0.00015)

    max_nstate = max(np.array([pt2_nstates, pt3_nstates, sum_nstates]))

    if max_nstate > 3:
        for i in range(3, max_nstate):
            prior['log(dE'+str(i)+')'] = gv.gvar(-1.25, 0.5)
            prior['z'+str(i)] = gv.gvar(0, 0.00015)

    prior['log(E_sum)'] = gv.gvar(-1.25, 0.5) 
    prior['z_sum'] = gv.gvar(0, 0.00015) 

    prior['A3_00'] = gv.gvar(1.2, 0.2*times)
    prior['V4_00'] = gv.gvar(1.0, 0.2*times)

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

def prior_ho_width_ge_times(pt2_nstates, pt3_nstates, sum_nstates, times): # change prior width of g.s. and 1 e.s.
    prior = gv.BufferDict()

    prior['E0'] = gv.gvar(0.50, 0.05*times)
    prior['z0'] = gv.gvar(0.00034, 0.00034*times)

    prior['log(dE1)'] = gv.gvar(-1.25, 0.5*times)
    prior['z1'] = gv.gvar(0, 0.00025*times)

    prior['log(dE2)'] = gv.gvar(-1.25, 0.5)
    prior['z2'] = gv.gvar(0, 0.00015)

    max_nstate = max(np.array([pt2_nstates, pt3_nstates, sum_nstates]))

    if max_nstate > 3:
        for i in range(3, max_nstate):
            prior['log(dE'+str(i)+')'] = gv.gvar(-1.25, 0.5)
            prior['z'+str(i)] = gv.gvar(0, 0.00015)

    # same garbage for sum and fh, so below are shared 
    prior['log(E_sum)'] = gv.gvar(-1.25, 0.5) 
    prior['z_sum'] = gv.gvar(0, 0.00015) 

    prior['A3_00'] = gv.gvar(1.2, 0.2*times)
    prior['V4_00'] = gv.gvar(1.0, 0.2*times)

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

    prior['A3_01'] = gv.gvar(0, 1*times)
    prior['A3_11'] = gv.gvar(0, 1*times)
    prior['V4_01'] = gv.gvar(0, 1*times)
    prior['V4_11'] = gv.gvar(1, 0.2*times)

    for i in range(sum_nstates-1):
        prior['sum_A3_'+str(i)] = gv.gvar(0, 1)
        prior['sum_V4_'+str(i)] = gv.gvar(0, 1)

    prior['sum_A3_'+str(sum_nstates-1)] = gv.gvar(0, 1)
    prior['sum_V4_'+str(sum_nstates-1)] = gv.gvar(1, 0.2)

    return prior

def prior_sw_width_1(pt2_nstates, pt3_nstates, sum_nstates): # square well
    prior = gv.BufferDict()

    prior['E0'] = gv.gvar(0.50, 0.05)
    prior['z0'] = gv.gvar(0.00034, 0.00034)

    prior['log(dE1)'] = gv.gvar(-1.25, 0.5)
    prior['z1'] = gv.gvar(0, 0.00025)

    prior['log(dE2)'] = gv.gvar(-1.25+np.log(3), 0.5/3)
    prior['z2'] = gv.gvar(0, 0.00015)

    max_nstate = max(np.array([pt2_nstates, pt3_nstates, sum_nstates]))

    if max_nstate > 3:
        for i in range(3, max_nstate):
            prior['log(dE'+str(i)+')'] = gv.gvar(-1.25+np.log(2*i-1), 0.5/(2*i-1))
            prior['z'+str(i)] = gv.gvar(0, 0.00015)

    prior['log(dE4)'] = gv.gvar(-1.25+np.log(2*4-1), 0.5/(2*4-1)*5)

    prior['log(E_sum)'] = gv.gvar(-1.25, 0.5) 
    prior['z_sum'] = gv.gvar(0, 0.00015) 

    prior['A3_00'] = gv.gvar(1.2, 0.2)
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

def prior_id1_width_1(pt2_nstates, pt3_nstates, sum_nstates): # 1/n
    prior = gv.BufferDict()

    prior['E0'] = gv.gvar(0.50, 0.05)
    prior['z0'] = gv.gvar(0.00034, 0.00034)

    prior['log(dE1)'] = gv.gvar(-1.25, 0.5)
    prior['z1'] = gv.gvar(0, 0.00025)

    prior['log(dE2)'] = gv.gvar(-1.25-np.log(2), 0.5*2)
    prior['z2'] = gv.gvar(0, 0.00015)

    max_nstate = max(np.array([pt2_nstates, pt3_nstates, sum_nstates]))

    if max_nstate > 3:
        for i in range(3, max_nstate):
            prior['log(dE'+str(i)+')'] = gv.gvar(-1.25-np.log(i), 0.5*i)
            prior['z'+str(i)] = gv.gvar(0, 0.00015)

    prior['log(E_sum)'] = gv.gvar(-1.25, 0.5) 
    prior['z_sum'] = gv.gvar(0, 0.00015) 

    prior['A3_00'] = gv.gvar(1.2, 0.2)
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

def prior_id2_width_1(pt2_nstates, pt3_nstates, sum_nstates): # 1/n**2
    prior = gv.BufferDict()

    prior['E0'] = gv.gvar(0.50, 0.05)
    prior['z0'] = gv.gvar(0.00034, 0.00034)

    prior['log(dE1)'] = gv.gvar(-1.25, 0.5)
    prior['z1'] = gv.gvar(0, 0.00025)

    prior['log(dE2)'] = gv.gvar(-1.25-2*np.log(2), 0.5*4)
    prior['z2'] = gv.gvar(0, 0.00015)

    max_nstate = max(np.array([pt2_nstates, pt3_nstates, sum_nstates]))

    if max_nstate > 3:
        for i in range(3, max_nstate):
            prior['log(dE'+str(i)+')'] = gv.gvar(-1.25-2*np.log(i), 0.5*(i**2))
            prior['z'+str(i)] = gv.gvar(0, 0.00015)

    prior['log(E_sum)'] = gv.gvar(-1.25, 0.5) 
    prior['z_sum'] = gv.gvar(0, 0.00015) 

    prior['A3_00'] = gv.gvar(1.2, 0.2)
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
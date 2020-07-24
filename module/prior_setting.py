import matplotlib.pyplot as plt
import numpy as np 
import gvar as gv  
def prior(max_nstate, pt3_nstates, sum_nstates, fh_nstates):
    prior = {}

    prior['E0'] = gv.gvar(0.665, 0.02)
    prior['log(dE1)'] = gv.gvar(-0.97, 0.69)

    prior['z0'] = gv.gvar(0.000775, 0.0000775)
    prior['z1'] = gv.gvar(0.0005, 0.0003)

    prior['z0_ps'] = gv.gvar(0.0015, 0.0005)
    prior['z1_ps'] = gv.gvar(0.0015, 0.001)

    prior['log(dE2)'] = gv.gvar(-0.97, 0.69)
    prior['z2'] = gv.gvar(0.0005, 0.0003)
    prior['z2_ps'] = gv.gvar(0.001, 0.001)

    if max_nstate > 3:
        for i in range(3, max_nstate):
            prior['log(dE'+str(i)+')'] = gv.gvar(-0.97, 0.69)
            prior['z'+str(i)] = gv.gvar(0.0003, 0.0003)
            prior['z'+str(i)+'_ps'] = gv.gvar(0.0005, 0.0005)

    # same garbage for sum and fh, so below are shared 
    prior['log(E_fh)'] = gv.gvar(-0.97, 0.69) 
    prior['z_fh_ss'] = gv.gvar(0.0005, 0.0003) 
    prior['z_fh_ps'] = gv.gvar(0.001, 0.001)

    prior['A3_00'] = gv.gvar(1.23, 0.05)
    prior['V4_00'] = gv.gvar(1.0, 0.2)

    nstate_array=np.array([pt3_nstates, sum_nstates, fh_nstates])
    ff_nstates = max(nstate_array)

    for i in range(ff_nstates):
        for j in range(ff_nstates):
            if i+j >= 1:
                if j < i:  
                    prior['A3_'+str(j)+str(i)] = gv.gvar(0, 1)
                    prior['V4_'+str(j)+str(i)] = gv.gvar(0, 1)

                elif j == i:
                    prior['A3_'+str(j)+str(i)] = gv.gvar(0, 1)
                    prior['V4_'+str(j)+str(i)] = gv.gvar(1, 0.2)

    for i in range(fh_nstates-1):
        prior['fh_A3_'+str(i)] = gv.gvar(0, 1)
        prior['fh_V4_'+str(i)] = gv.gvar(0, 1)
        prior['d'+str(i)+'_ss_A3'] = gv.gvar(-0.0000015, 0.0000015)
        prior['d'+str(i)+'_ss_V4'] = gv.gvar(0.0000013, 0.0000013)
        prior['d'+str(i)+'_ps_A3'] = gv.gvar(-0.000009, 0.000009)
        prior['d'+str(i)+'_ps_V4'] = gv.gvar(0.0000075, 0.0000075)

    prior['d0_ss_A3'] = gv.gvar(-0.0000015, 0.0000015) # rewrite the d0
    prior['d0_ss_V4'] = gv.gvar(0.0000013, 0.0000013) 

    prior['d'+str(fh_nstates-1)+'_ss_A3'] = gv.gvar(-0.0000015, 0.0000015) # add the d(nstates-1)
    prior['d'+str(fh_nstates-1)+'_ss_V4'] = gv.gvar(0.0000013, 0.0000013) 

    prior['d0_ps_A3'] = gv.gvar(-0.000009, 0.000009) # rewrite the d0
    prior['d0_ps_V4'] = gv.gvar(0.0000075, 0.0000075) 

    prior['d'+str(fh_nstates-1)+'_ps_A3'] = gv.gvar(-0.000009, 0.000009) # add the d(nstates-1)
    prior['d'+str(fh_nstates-1)+'_ps_V4'] = gv.gvar(0.0000075, 0.0000075)

    prior['fh_A3_'+str(fh_nstates-1)] = gv.gvar(0, 1)
    prior['fh_V4_'+str(fh_nstates-1)] = gv.gvar(1, 0.2)

    return prior
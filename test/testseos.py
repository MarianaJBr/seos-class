from classyy import Class
import numpy as np
cosmo = Class()
params={'h':'0.6775',
        'T_cmb':'2.7255', 
        'omega_b':'0.0230',
        'N_ur':'3.046',
        'Omega_cdm':'0.12038',
        'Omega_dcdmdr':'0.0',
        'Gamma_dcdm':'0.0',
        'N_ncdm':'0',
        'Omega_k':'0.0',
        'Omega_fld':'0.7',
	'w0_fld':'-0.92',
	'wa_fld':'0.07',
	'q_fld':'2',
	'zt_fld':'0.8',
	'cs2_fld':'1.0',
        'P_k_max_h/Mpc':'100','output':'mPk','k_pivot':'0.05',
        'A_s':'2.215e-9',
        'n_s':'0.9619',
        'nonlinear_verbose':'1',
        'lensing_verbose':'1',
        'z_reio':'11.357',
        'gauge':'synchronous',
        'alpha_s':'0.',
        'z_pk':'0.',
        'P_k_ini type': 'analytic_Pk'
        }
cosmo.set(params)
cosmo.compute()
print cosmo.sigma8()
print cosmo.q_fld()
print cosmo.zt_fld()

#universal constants definitions

import numpy as np

epsilon = 1e-10

hPl = 6.582119569*1e-25 #Planck constant
alpha_EM = 1./137 #EM fine structure, maybe to more precise?
g_EM = np.sqrt(4*np.pi*alpha_EM)
g_W = 0.652905 #Weak coupling (2018PDG and physics.nist.gov)
alpha_W = g_W**2/(4*np.pi)
G_F = 1.1663787e-5 #Fermi coupling constant (in GeV^-2) 2021PDG
N_8 = 1.53e-7 # from 2110.10698 (without phase)
m_W = 80.379 # W mass in GeV 2018PDG
v = 2*m_W/g_W # Higgs vev
m_Z = 91.1876
theta_w = np.arccos(m_W/m_Z)
m_H = 125.1
m_el = 0.000511
m_mu = 0.105658
m_tau = 1.77686
m_q = [0.00216, 0.00467, 0.093, 1.275, 4.18, 173.0] #quark mass 2019PDG u,d,s,c,b,t
m_pi = 0.13957
m_pi0 = 0.13498
m_eta = 0.54786
m_etap = 0.95778
m_rho = 0.77526
m_omega = 0.78265
m_phi = 1.01946
m_B0 = 5.27966 #B0 mass in GeV 2022PDG
m_B = 5.27934 #B^+- mass in GeV 2022PDG
m_Bs = 5.4154 #Bs mass in GeV 2023PDG
m_Bc = 6.27447 #Bc mass in GeV 2023PDG
m_Bstar = 5.32471 #B*^+- mass in GeV 2023PDG
m_Bsstar = 5.4154 #Bs* mass in GeV 2023PDG
m_Bcstar = 6.331 #Bc*^+- value in GeV Lattice2016
#K masses in GeV 2024PDG
m_K = {"K" :    0.493677,
       "K_0" :  0.497611,
       "K0star_700" :   0.845,
       "K0star_1430" :  1.425,
       "Kstar_892" :    0.89176,
       "Kstar_892_0" :  0.89555,
       "Kstar_1410" :   1.414,
       "Kstar_1680" :   1.718,
       "K1_1270" :      1.253,
       "K1_1400" :      1.403,
       "K2star_1430" :  1.4273,
       "K2star_1430_0" :1.4324,
       }

m_D0 = 1.86483; #D0 mass in GeV 2020PDG
m_D = 1.86968; #D^+- mass in GeV 2020PDG
m_Dstar = 2.01026 #D*^+- mass in GeV 2023PDG
m_D0star = 2.00685 #D*0 mass in GeV 2023PDG
m_Ds = 1.96834; #Ds^+- mass in GeV 2020PDG
m_Dsstar = 2.1122; #Ds*^+- mass in GeV 2023PDG
m_p = 0.93827


BR_Pi0_2Gamma = 0.98823
BR_Eta_2Gamma = 0.3941
BR_EtaPrim_2Gamma = 0.0222
BR_EtaPrim_2PiEta = 0.427 #EtaPrim -> 2PiEta PDG 2020
BR_EtaPrim_2Pi0Eta = 0.228 #EtaPrim -> 2Pi0Eta PDG 2020
BR_Rho_PiGamma = 4.5*1e-4
BR_Rho0_Pi0Gamma = 4.7*1e-4
BR_Omega_Pi0Gamma = 0.084
BR_Phi_Pi0Gamma = 1.3*1e-3
BR_Rho0_EtaGamma = 3*1e-4
BR_Omega_EtaGamma = 4.5*1e-4
BR_Phi_EtaGamma = 0.0103

ratio_EtaPrim_2PiEta_neutral = BR_EtaPrim_2Pi0Eta/(BR_EtaPrim_2PiEta + BR_EtaPrim_2Pi0Eta) #portion of neutral pipieta decay from the whole (should be 1/3 based on the symmetry)
ratio_EtaPrim_2PiEta_charged = BR_EtaPrim_2PiEta/(BR_EtaPrim_2PiEta + BR_EtaPrim_2Pi0Eta) #portion of charged pipieta decay from the whole (should be 2/3 based on the symmetry)

BR_D_NuTau = 1.20e-3 # D  -> nu_tau tau PDG 2022
BR_Ds_NuTau = 0.0532 # Ds -> nu_tau tau PDG 2022

#CKM 2020PDG
V_ub = 0.00361
V_us = 0.2265
V_ud = 0.974
V_cb = 0.0405
V_cs = 0.9732
V_cd = 0.2263
V_tb = 0.99917
V_ts = 0.03978
V_td = 0.0085

B0 = m_pi0**2 /(m_q[0] + m_q[1])
fpi = 0.093
feta = 0.093
fetap = 0.073
th_eta_etap = np.arcsin(-1/3.)

#Hadronic decay constants 
f_DV_cd  = 0.04550 #GeV PDG2018
f_DsV_cs = 0.2453 #GeV PDG2022
f_BV_ub  = 0.00077 #GeV PDG2018
f_pi = 0.1302 #GeV FLAG Review 2019
f_K  = 0.1557 #GeV FLAG Review 2019
f_B0 = 0.1905  #GeV FLAG Review 2019
f_Bs = 0.2303  #GeV FLAG Review 2019

#.. to be validated: (taken from JHEP 0710)
f_eta = 1.2*fpi
f_etaprime = 0.45*fpi

g_rho = 0.162 #GeV^2 Phys. Lett. B635 (2006) 93–99, This value is very inconsistent in literature

#Life times PDG2022 in 
tau_K  = 1.238e-8 / hPl #GeV^-1
tau_KS  = 0.8954e-10 / hPl #GeV^-1
tau_KL  = 5.116e-8 / hPl #GeV^-1
tau_D  = 1033.1e-15 / hPl #GeV^-1
tau_D0 = 410.3e-15 / hPl #GeV^-1
tau_Ds = 504e-15 / hPl #GeV^-1
tau_tau = 290.3e-15 / hPl #GeV^-1
tau_B  = 1638e-15 / hPl #GeV^-1
tau_B0 = 1519e-15 / hPl #GeV^-1
tau_Bs = 1520e-15 / hPl #GeV^-1

#resulting widths
Gamma_K = 1/tau_K
Gamma_KS = 1/tau_KS
Gamma_KL = 1/tau_KL
Gamma_D = 1/tau_D
Gamma_D0 = 1/tau_D0
Gamma_Ds = 1/tau_Ds
Gamma_B = 1/tau_B
Gamma_B0 =1/tau_B0
Gamma_Bs = 1/tau_Bs

#scalars
Gamma_f0 = [0.55,0.055,0.35] #GeV
m_f0 = [0.475,0.98,1.35] #GeV
f_f0 = [0.28,1.8,-0.99] #GeV

#vectors
Gamma_rho0 = [0.1474,0.4,0.25] #GeV
m_rho0 = [0.77526,1.465,1.72] #GeV
f_rho0 = [0.616,0.223,-0.339] #GeV
Gamma_omega0 = [0.00868,0.29,0.315] #GeV
m_omega0 = [0.78266,1.41,1.67] #GeV
f_omega0 = [1.011,-0.881,0.369] #GeV

CF = 4/3 # Casimir fundamental rep for SU(3)


# effective dark scalar couplings

C_sd_S = 3*m_q[2]*g_W**2/(v*64*np.pi**2)*(V_us*V_ud*m_q[0]**2+V_cs*V_cd*m_q[3]**2+V_ts*V_td*m_q[5]**2)/(m_W**2)
C_bd_S = 3*m_q[4]*g_W**2/(v*64*np.pi**2)*(V_ub*V_ud*m_q[0]**2+V_cb*V_cd*m_q[3]**2+V_tb*V_td*m_q[5]**2)/(m_W**2)
# C_bs_S = 3*V_tb*V_ts*m_q[5]**2/(v*v*v*16*np.pi**2) / 1e3 # for recasting ALP B meson production to DS (add 1e-3 for MeV^-1)
C_bs_S = 3*m_q[4]*g_W**2/(v*64*np.pi**2)*(V_ub*V_us*m_q[0]**2+V_cb*V_cs*m_q[3]**2+V_tb*V_ts*m_q[5]**2)/(m_W**2)
C_bs_2S = 3*m_q[4]*g_W**2/(m_H**2*64*np.pi**2)*(V_ub*V_us*m_q[0]**2+V_cb*V_cs*m_q[3]**2+V_tb*V_ts*m_q[5]**2)/(m_W**2)

# production cross sections 
sigma_pp = { # based on https://pdg.lbl.gov/2015/hadronic-xsections/rpp2014-hadronic_xsec_fits.pdf
  30  : 38.38*1E9,
  70  : 38.38*1E9,
  120 : 38.54*1E9,
  400 : 39.85*1E9,
  800 : 41.14*1E9
}

sigma_cc = { # based on [2407.08673]
  30  : 0.,
  70  : 2.2*1E6,  # ± 1.2
  120 : 4.6*1E6,  # ± 1.8
  400 : 20.*1E6,  # ± 4.5 
  800 : 39.5*1E6  # ± 7.7
}

sigma_bb = { # based on [2407.08673]
  30  : 0.,
  70  : 0.,
  120 : 20., #± 20
  400 : 2.63*1E3, #±0.76
  800 : 18.6*1E3 # ± 4.6
}


# photon-proton sigma_abs from https://www.nist.gov/pml/xcom-photon-cross-sections-database
sigma_pgamma = {            
              'NuTeV':   1.45*1e11,
              'NA62':    6.*1e12,
              'CHARM':   6.*1e12, 
              'NuCal':   5.1*1e12, 
              'SHiP':    12.66*1e12, 
              'DarkQuest':5.1*1e12, 
              'DUNE':    0.35*1e12,
              'SHADOWS': 6.*1e12, 
              'KOTO':    38.1*1e12, 
              'KOTO2':   38.1*1e12,
              'BEBC':    6.*1e12,
              'ORCA':    0.35*1e12,
}

# target material volume averaged proton count per nucleus
Z_target =  {'NuTeV':   6.55872,
             'NA62':    29,
             'CHARM':   29, 
             'NuCal':   26, 
             'SHiP':    42, 
             'DarkQuest':26, 
             'DUNE':    6,
             'SHADOWS': 29, 
             'KOTO':    79, 
             'KOTO2':   79,
             'BEBC':    29,
             'ORCA':    6,
}

# target material volume averaged nucleon count per nucleus
A_target =  {'NuTeV':   12.5056,
             'NA62':    63.546,
             'CHARM':   63.546, 
             'NuCal':   56, 
             'SHiP':    95, 
             'DarkQuest':56, 
             'DUNE':    12,
             'SHADOWS': 63.546, 
             'KOTO':    197, 
             'KOTO2':   197,
             'BEBC':    63.546,
             'ORCA':    12,
}
#universal constants definitions

import numpy as np

epsilon = 1e-10

hPl = 6.582119569*1e-25 #red. Planck constant in GeV.s
alpha_EM = 1./137 #EM fine structure, maybe to more precise?
g_EM = np.sqrt(4*np.pi*alpha_EM)
g_W = 0.652905 #Weak coupling (2018PDG and physics.nist.gov)
alpha_W = g_W**2/(4*np.pi)
G_F = 1.1663787e-5 #Fermi coupling constant (in GeV^-2) 2021PDG
N_8 = 1.53e-7 # from 2110.10698 (without phase)
m_W = 80.379 # W mass in GeV 2018PDG
v = 2*m_W/g_W # Higgs vev
m_Z = 91.1876
m_H = 125.2
theta_w = np.arccos(m_W/m_Z)
sin_w2 = 1-(m_W/m_Z)**2
#Lepton masses in GEV PDG 2024
m_el = 5.1099895e-4
m_mu = 0.1056583755
m_tau = 1.77693
#Quark masses in GEV PDG 2024
m_u = 2.16e-3
m_c = 1.273
m_t = 172.57
m_d = 4.7e-3
m_s = 93.5e-3
m_b = 4.183
# fermion mass collections
m_q = [m_u, m_d, m_s, m_c, m_b, m_t] #lightest to heaviest
m_qt = [m_u,m_c,m_t]
m_qb = [m_d,m_s,m_b]
m_l = [m_el,m_mu,m_tau]

m_pi = 0.13957
m_pi0 = 0.13498
m_eta = 0.54786
m_etap = 0.95778
m_rho = 0.77526 # pdg rpp2019-list-rho-770
m_B0 = 5.27963 #B0 mass in GeV 2018PDG
m_B = 5.27932 #B^+- mass in GeV 2018PDG
m_K0 = 0.49761 #K0 mass in GeV 2018PDG
m_K = 0.49367 #K^+- mass in GeV 2018PDG
m_K0star = 0.89555 #K0* mass in GeV 2018PDG
m_Kstar = 0.89176 #K^+-* mass in GeV 2018PDG
m_D0 = 1.86483; #D0 mass in GeV 2020PDG
m_D = 1.86968; #D^+- mass in GeV 2020PDG
m_Ds = 1.96834; #Ds^+- mass in GeV 2020PDG
m_Dsstar = 2.1122 #Ds^+-* mass in GeV 2020PDG
m_omega  = 0.783 #omega0 mass in GeV 2020PDG
m_phi = 1.0195   #Phi0 mass in GeV 2020PDG
m_Jpsi = 3.0969  #J/Psi mass in GeV 2020PDG

BR_Pi0_2Gamma = 0.98823
BR_Eta_2Gamma = 0.3941
BR_EtaPrim_2Gamma = 0.0222
BR_EtaPrim_2PiEta = 0.427 #EtaPrim -> 2PiEta PDG 2020
BR_EtaPrim_2Pi0Eta = 0.228 #EtaPrim -> 2Pi0Eta PDG 2020
ratio_EtaPrim_2PiEta_neutral = BR_EtaPrim_2Pi0Eta/(BR_EtaPrim_2PiEta + BR_EtaPrim_2Pi0Eta) #portion of neutral pipieta decay from the whole (should be 1/3 based on the symmetry)
ratio_EtaPrim_2PiEta_charged = BR_EtaPrim_2PiEta/(BR_EtaPrim_2PiEta + BR_EtaPrim_2Pi0Eta) #portion of charged pipieta decay from the whole (should be 2/3 based on the symmetry)

BR_D_NuTau = 1.20e-3 # D  -> nu_tau tau PDG 2022
BR_Ds_NuTau = 0.0532 # Ds -> nu_tau tau PDG 2022

#.. this information is subject to production energy these are values for 400GeV production
ratio_D_cc  = 0.31
ratio_Ds_cc = 0.09

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

V_CKM = np.array([[0.974, 0.2265, 0.00361],
                  [0.2263, 0.9732, 0.0405],
                  [0.0085, 0.03978, 0.99917]])

B0 = m_pi0**2 /(m_q[0] + m_q[1])
delta_I = (m_q[1] - m_q[0])/(m_q[0] + m_q[1])

fpi = 0.093
feta = 0.093
fetap = 0.073
th_eta_etap = np.arcsin(-1/3.)

#Hadronic decay constants 
f_DV_cd  = 0.04550 #GeV PDG2018
f_DsV_cs = 0.2453  #GeV PDG2022
f_BV_ub  = 0.00077 #GeV PDG2018
f_piV_ud = 0.12713 #GeV PDG 2018
f_KV_us  = 0.03509 #GeV PDG 2018
f_pi = 0.1302  #GeV FLAG Review 2019
f_K  = 0.1557  #GeV FLAG Review 2019
f_Bs = 0.2303  #GeV FLAG Review 2019

#.. to be validated:
f_eta = 0.0817   # as presented in INR-TH-2018-014
f_etap = -0.0947 # as presented in INR-TH-2018-014
g_rho = 0.162 #GeV^2 Phys. Lett. B635 (2006) 93–99 ..This value is very inconsistent in literature
g_Dsstar = 0.650 #GeV^2 Eur. Phys. J. Plus 133, 134 (2018). 
g_omega  = 0.195*m_omega #GeV^2 hep-ph/0602110
g_phi    = f_pi*m_pi  
g_phi0   = f_pi*m_pi0
g_Jpsi   = 0.406*m_Jpsi  #GeV^2 hep-ph/9703364v

#Life times PDG2022 in 
tau_K  = 1.238e-8 / hPl #GeV^-1
tau_D  = 1033.1e-15 / hPl #GeV^-1
tau_D0 = 410.3e-15 / hPl #GeV^-1
tau_Ds = 504e-15 / hPl #GeV^-1
tau_tau = 290.3e-15 / hPl #GeV^-1
tau_B  = 1638e-15 / hPl #GeV^-1
tau_B0 = 1519e-15 / hPl #GeV^-1
tau_Bs = 1520e-15 / hPl #GeV^-1

#resulting widths
Gamma_K =1/tau_K
Gamma_B0 =1/tau_B0
Gamma_B = 1/tau_B
Gamma_D0 = 1/tau_D0
Gamma_D = 1/tau_D

CF = 4/3 # Casimir fundamental rep for SU(3)

C_bs_S = 3*V_tb*V_ts*m_q[5]**2/(v*v*v*16*np.pi**2) / 1e3 # for recasting ALP B meson production to DS (add 1e-3 for MeV^-1)

prod_xsec_ref = { # reference production cross section and standard deviation values from [2407.08673]
    'Dmeson' : { # in nb
        'E137':     [0, 0],
        'E141':     [0, 0],
        "KOTO":     [0, 0],
        "KOTO2":    [0, 0],
        "NuCal":    [2.2e3, 1.2e3],
        "DUNE" :    [4.6e3, 1.8e3],
        "DarkQuest":[4.6e3, 1.8e3],
        "NA62"  :   [20.0e3, 4.5e3],
        "BEBC"  :   [20.0e3, 4.5e3],
        "CHARM" :   [20.0e3, 4.5e3],
        'SHADOWS':  [20.0e3, 4.5e3],
        "SHiP"  :   [20.0e3, 4.5e3],
        "ORCA"  :   [20.0e3, 4.5e3],
        "NuTeV" :   [39.5e3, 7.7e3]
    },
    'Bmeson' : { # in nb
        'E137':     [0, 0],
        'E141':     [0, 0],
        "KOTO":     [0, 0],
        "KOTO2":    [0, 0],
        "NuCal":    [0, 0],
        "DUNE" :    [0.02, 0.02],
        "DarkQuest":[0.02, 0.02],
        "NA62"  :   [2.63, 0.76],
        "BEBC"  :   [2.63, 0.76],
        "CHARM" :   [2.63, 0.76],
        'SHADOWS':  [2.63, 0.76],
        "SHiP"  :   [2.63, 0.76],
        "ORCA"  :   [2.63, 0.76],
        "NuTeV" :   [18.6, 4.6]
    }
}


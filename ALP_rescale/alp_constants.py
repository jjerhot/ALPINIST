#universal constants definitions

import numpy as np

epsilon = 1e-10

hPl = 6.582119569*1e-25 #Planck constant
alpha_EM = 1./137 #EM fine structure, maybe to more precise?
g_EM = np.sqrt(4*np.pi*alpha_EM)
g_W = 0.652905 #Weak coupling (2018PDG and physics.nist.gov)
alpha_W = g_W**2/(4*np.pi)
m_W = 80.379 # W mass in GeV 2018PDG
v = 2*m_W/g_W # Higgs vev
m_Z = 91.1876
theta_w = np.arccos(m_W/m_Z)
m_el = 0.000511
m_mu = 0.105658
m_tau = 1.77686
m_q = [0.00216, 0.00467, 0.093, 1.275, 4.18, 173.0] #quark mass 2019PDG u,d,s,c,b,t
m_pi = 0.13957
m_pi0 = 0.13498
m_eta = 0.54786
m_etap = 0.95778
m_B0 = 5.27963 #B0 mass in GeV 2018PDG
m_B = 5.27932 #B^+- mass in GeV 2018PDG
Gamma_B0 = hPl/(1.52*1e-12) #B0 total decay width in GeV
Gamma_B = hPl/(1.638*1e-12) #B^+- total decay width in GeV
Gamma_D0 = hPl/(0.41*1e-12) #D0 total decay width in GeV
Gamma_D = hPl/(1.04*1e-12) #D^+- total decay width in GeV
m_K0 = 0.49761 #K0 mass in GeV 2018PDG
m_K = 0.49367 #K^+- mass in GeV 2018PDG
m_K0star = 0.89555 #K0* mass in GeV 2018PDG
m_Kstar = 0.89176 #K^+-* mass in GeV 2018PDG
m_D0 = 1.86483; #D0 mass in GeV 2020PDG
m_D = 1.86968; #D^+- mass in GeV 2020PDG

BR_Pi0_2Gamma = 0.98823
BR_Eta_2Gamma = 0.3941
BR_EtaPrim_2Gamma = 0.0222
BR_EtaPrim_2PiEta = 0.427 #EtaPrim -> 2PiEta PDG 2020
BR_EtaPrim_2Pi0Eta = 0.228 #EtaPrim -> 2Pi0Eta PDG 2020
ratio_EtaPrim_2PiEta_neutral = BR_EtaPrim_2Pi0Eta/(BR_EtaPrim_2PiEta + BR_EtaPrim_2Pi0Eta) #portion of neutral pipieta decay from the whole (should be 1/3 based on the symmetry)
ratio_EtaPrim_2PiEta_charged = BR_EtaPrim_2PiEta/(BR_EtaPrim_2PiEta + BR_EtaPrim_2Pi0Eta) #portion of charged pipieta decay from the whole (should be 2/3 based on the symmetry)


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

CF = 4/3 # Casimir fundamental rep for SU(3)
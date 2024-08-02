import numpy as np
# import mpmath as mp #for polylog
from scipy.interpolate import RectBivariateSpline, interp1d
from os import path
from ALP_rescale.general import constants as c
from ALP_rescale.general import functions as f

class Load_digitized:
    def __init__(self, mode):
        self.digitized_exists = 0
        self.digitized_ratio = 0
        self._decay_mode = self._decay_mode_filename = mode
        if mode == '2Pi0Eta': self._decay_mode_filename = '2PiEta'
        if mode == '2Pi0EtaPrim': self._decay_mode_filename = '2PiEtaPrim'  
        if mode == '2Pi2Pi0' or mode == '4Pi': self._decay_mode_filename = '2Rho' 
        if mode == '2KPi0': self._decay_mode_filename = '2KPi'        
        width_file_name = path.dirname(path.realpath(__file__))+'/../../widths/integrated/digitized_width_'+self._decay_mode_filename+'.dat'
        width_file_name_fermion = path.dirname(path.realpath(__file__))+'/../../widths/fermion_ratios/ratio_'+self._decay_mode_filename+'.dat'
        if path.exists(width_file_name):
            width_digitized = np.loadtxt(width_file_name)
            nlines = int(width_digitized.size/2)
            m_a_list = np.array([width_digitized[i,0] for i in range(nlines)])
            Gamma_a_list = np.array([width_digitized[i,1] for i in range(nlines)])
            self._width_inter = interp1d(m_a_list, Gamma_a_list,fill_value="extrapolate")

            if not mode == ('2Gamma' or '2El' or '2Mu'): self.digitized_exists = 1 #for 2Gamma/2El/2Mu use analytical instead           
        else: 
            print('[Warning:] \t',width_file_name,'not found')
        if path.exists(width_file_name_fermion):
            ratio_digitized = np.loadtxt(width_file_name_fermion)
            nlines = int(ratio_digitized.size/2)
            m_a_list = np.array([ratio_digitized[i,0] for i in range(nlines)])
            ratio_a_list = np.array([ratio_digitized[i,1] for i in range(nlines)])
            self._ratio_inter = interp1d(m_a_list, ratio_a_list,fill_value="extrapolate")
            self.digitized_ratio = 1
        else:
            width_file_name_fermion = path.dirname(path.realpath(__file__))+'/../../widths/fermion_widths/digitized_width_'+self._decay_mode_filename+'.dat'
            if path.exists(width_file_name_fermion):
                width_digitized_ferm = np.loadtxt(width_file_name_fermion)
                nlines = int(width_digitized_ferm.size/2)
                m_a_list_ferm = np.array([width_digitized_ferm[i,0] for i in range(nlines)])
                Gamma_a_ferm_list = np.array([width_digitized_ferm[i,1] for i in range(nlines)])
                self._width_inter_ferm = interp1d(m_a_list_ferm, Gamma_a_ferm_list,fill_value="extrapolate")

            else:
                print('[Warning:] \t',width_file_name_fermion,'not found')

    def _mass_bound(self,m_a):
        return {
            'total'         : 1,
            '2Gamma'        : 1,
            '2El'           : 1 if m_a > (2*c.m_el) else 0,
            '2Mu'           : 1 if m_a > (2*c.m_mu) else 0,
            '3Pi0'          : 1 if m_a > (3*c.m_pi0) and m_a < c.m_etap else 0,
            '3Pi'           : 1 if m_a > (2*c.m_pi+c.m_pi0) and m_a < c.m_etap else 0,
            '2PiGamma'      : 1 if m_a > (2*c.m_pi) else 0,
            '2Pi0Eta'       : 1 if m_a > (2*c.m_pi0+c.m_eta) else 0,
            '2PiEta'        : 1 if m_a > (2*c.m_pi+c.m_eta) else 0,
            '2Pi0EtaPrim'   : 1 if m_a > (2*c.m_pi0+c.m_etap) else 0,
            '2PiEtaPrim'    : 1 if m_a > (2*c.m_pi+c.m_etap) else 0,
            '2KPi0'         : 1 if m_a > (2*c.m_K+c.m_pi0) else 0,
            '2Pi2Pi0'       : 1 if m_a > (2*c.m_pi+2*c.m_pi0) else 0,
            '4Pi'           : 1 if m_a > (4*c.m_pi) else 0,
            '2Omega'        : 1 if m_a > (2*c.m_omega) else 0,
            '2Kstar'        : 1 if m_a > (2*c.m_K+2*c.m_pi) else 0,
            '2Phi'          : 1 if m_a > (4*c.m_K) else 0,
        }.get(self._decay_mode, 0) 

    def _rescale(self):
        return {
            'total'         : 1.,
            '2Gamma'        : 1.,
            '2El'           : 1.,
            '2Mu'           : 1.,
            '3Pi0'          : c.BR_Pi0_2Gamma**3,
            '3Pi'           : c.BR_Pi0_2Gamma,
            '2PiGamma'      : 1.,
            '2Pi0Eta'       : c.BR_Pi0_2Gamma**2 * c.BR_Eta_2Gamma *c.ratio_EtaPrim_2PiEta_neutral,
            '2PiEta'        : c.BR_Eta_2Gamma * c.ratio_EtaPrim_2PiEta_charged,
            '2Pi0EtaPrim'   : c.BR_Pi0_2Gamma**2 * c.BR_EtaPrim_2Gamma * c.ratio_EtaPrim_2PiEta_neutral,
            '2PiEtaPrim'    : c.BR_EtaPrim_2Gamma * c.ratio_EtaPrim_2PiEta_charged, 
            '2KPi0'         : c.BR_Pi0_2Gamma * 1./6,
            '2Pi2Pi0'       : c.BR_Pi0_2Gamma**2 * 2./3,
            '4Pi'           : 1./3,
            '2Omega'        : 1,
            '2Kstar'        : 1.,
            '2Phi'          : 1.,           
        }.get(self._decay_mode, 0) 
    
    def width_from_file(self,m_a,g_A):
        if self._mass_bound(m_a) and g_A != 0: return self._rescale() * g_A**2 * self._width_inter(m_a) * 1e-9 * 1e+6 # was in eV scale and Lambda = 1TeV
        else: return 0.

    def width_from_file_input(self,m_a,g_A):
        if self._mass_bound(m_a) and g_A != 0: return g_A**2 * self._width_inter(m_a) * 1e-9 * 1e+6 # was in eV scale and Lambda = 1TeV
        else: return 0.

    def width_ferm_from_file(self,m_a,g_A):
        if self._mass_bound(m_a) and g_A != 0: 
            if self.digitized_ratio == 1:
                return self._rescale() * g_A**2 * self._ratio_inter(m_a) * self._width_inter(m_a) * 1e-9 * 1e+6 # was in eV scale and Lambda = 1TeV
            else:
                return self._rescale() * g_A**2 * self._width_inter_ferm(m_a) * 1e-9 * 1e+6
        else: return 0.

    def width_ferm_from_file_input(self,m_a,g_A):
        if self._mass_bound(m_a) and g_A != 0:
            if self.digitized_ratio == 1:
                return g_A**2 * self._ratio_inter(m_a) * self._width_inter(m_a) * 1e-9 * 1e+6 # was in eV scale and Lambda = 1TeV
            else:
                return g_A**2 * self._width_inter_ferm(m_a) * 1e-9 * 1e+6
        else: return 0.

class Total_width:  #total decay width
    def __init__(self, is_fermionic = False):
        self._is_fermionic = is_fermionic
        self._total = Load_digitized("total")
        self._digamma = Load_digitized("2Gamma")
        if is_fermionic:
            self._3Pi0 = Load_digitized("3Pi0")
            self._3Pi = Load_digitized("3Pi")
            self._2PiGamma = Load_digitized("2PiGamma")
            self._2PiEta = Load_digitized("2PiEta")
            self._2PiEtaPrim = Load_digitized("2PiEtaPrim")
            self._2KPi0 = Load_digitized("2KPi0")
            self._4Pi = Load_digitized("4Pi")
            self._2Omega = Load_digitized("2Omega")
            self._2Kstar = Load_digitized("2Kstar")
            self._2Phi = Load_digitized("2Phi")

    def width_a(self,m_a,g_gg,g_GG,g_ee=0,g_mumu=0,g_qq=0):
        if self._is_fermionic:
            return (a_2Gamma(m_a, g_gg) + a_2El(m_a,g_ee) + a_2Mu(m_a,g_mumu) 
                    + self._3Pi0.width_ferm_from_file_input(m_a,g_qq) 
                    + self._3Pi.width_ferm_from_file_input(m_a,g_qq) 
                    + self._2PiGamma.width_ferm_from_file_input(m_a,g_qq) 
                    + self._2PiEta.width_ferm_from_file_input(m_a,g_qq) 
                    + self._2PiEtaPrim.width_ferm_from_file_input(m_a,g_qq) 
                    + self._2KPi0.width_ferm_from_file_input(m_a,g_qq) 
                    + self._4Pi.width_ferm_from_file_input(m_a,g_qq) 
                    + self._2Omega.width_ferm_from_file_input(m_a,g_qq) 
                    + self._2Kstar.width_ferm_from_file_input(m_a,g_qq) 
                    + self._2Phi.width_ferm_from_file_input(m_a,g_qq)
                    )
        else:
            if g_GG != 0. and m_a > 0.3 and m_a <= 3.0:
                #reproduce hadronic branching ratio with digitized photon width and then obtain hadronic_width with our analytical photon width
                width_a_hadronic = self._total.width_from_file(m_a,g_GG) - self._digamma.width_from_file(m_a,g_GG)
                if width_a_hadronic < 0.:
                    width_a_hadronic = 0.
            else:
                width_a_hadronic = 0.
            return a_2Gamma(m_a, g_gg) + a_2El(m_a,g_ee) + a_2Mu(m_a,g_mumu) + width_a_hadronic

### analytical widths:

# ALP decay:
def a_2Gamma(m_a, g_gg): #ALP -> 2 gamma width
    return c.g_EM**4. * g_gg**2. * m_a**3./ (4*np.pi) 

def a_2Mu(m_a, g_ll): #ALP -> 2mu width
    if m_a > 2*c.m_mu: return g_ll**2. * c.m_mu**2 * np.sqrt(m_a**2 - 4 * c.m_mu**2)  / (8. * np.pi)
    return 0.

def a_2El(m_a, g_ll): #ALP -> 2mu width
    if m_a > 2*c.m_el: return g_ll**2. * c.m_el**2 * np.sqrt(m_a**2 - 4 * c.m_el**2)  / (8. * np.pi)
    return 0.

analytical_widths ={'2Gamma':   a_2Gamma,
                    '2Mu':      a_2Mu,
                    '2El':      a_2El}

# ALP production
def B_K_a(m_a, g_bs): #B -> K + ALP
    f_0 = 0.330/(1-m_a**2/37.46) #hep-ph/0406232
    if c.m_B > m_a + c.m_K: return np.abs(g_bs)**2 /(64*np.pi*c.m_B**3)*(c.m_B**2 - c.m_K**2)**2 * np.sqrt(f.lambda_Kallen(c.m_B,c.m_K,m_a)) * f_0**2
    return 0.

def B_Kstar_a(m_a, g_bs): #B -> K* + ALP
    A_0 = 1.364/(1-m_a**2/(c.m_B**2)) - 0.990/(1-m_a**2/36.78) #hep-ph/0412079
    if c.m_B > m_a + c.m_Kstar: return np.abs(g_bs)**2 /(64*np.pi*c.m_B**3) * np.sqrt(f.lambda_Kallen(c.m_B,c.m_Kstar,m_a)**3) * A_0**2
    return 0.

def D_pi_a(m_a, g_cu): #D -> pi + ALP
    f_0 = 0.612/(1-m_a**2/6.46)
    if c.m_D0 > m_a + c.m_pi: return np.abs(g_cu)**2 /(64*np.pi*c.m_D**3)*(c.m_D**2 - c.m_pi**2)**2 * np.sqrt(f.lambda_Kallen(c.m_D,c.m_pi,m_a)) * f_0**2
    return 0.

# Yukawa
def B_K_a_Y(m_a, g_bs): #B -> K + ALP
    f_0 = 0.330/(1-m_a**2/37.46) #hep-ph/0406232
    if c.m_B > m_a + c.m_K: return np.abs(g_bs)**2 /(16*np.pi*c.m_B**3)/(c.m_q[4]-c.m_q[2])**2 * (c.m_B**2 - c.m_K**2)**2 * np.sqrt(f.lambda_Kallen(c.m_B,c.m_K,m_a)) * f_0**2
    return 0.

def B_Kstar_a_Y(m_a, g_bs): #B -> K* + ALP
    A_0 = 1.364/(1-m_a**2/(c.m_B**2)) - 0.990/(1-m_a**2/36.78) #hep-ph/0412079
    if c.m_B > m_a + c.m_Kstar: return np.abs(g_bs)**2 /(16*np.pi*c.m_B**3)/(c.m_q[4]+c.m_q[2])**2 * np.sqrt(f.lambda_Kallen(c.m_B,c.m_Kstar,m_a)**3) * A_0**2
    return 0.

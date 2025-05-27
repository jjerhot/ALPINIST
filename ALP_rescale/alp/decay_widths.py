import numpy as np
# import mpmath as mp #for polylog
from scipy.interpolate import RectBivariateSpline, interp1d
from os import path
from ..general import constants as c
from ..general import functions as f

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

class Load_digitized_fermions: # digitized widths from 2310.03524
    def __init__(self, mode):
        self.digitized_exists = 0
        self._decay_mode = self._decay_mode_filename = mode
        width_file_name = path.dirname(path.realpath(__file__))+'/../../widths/integrated_fermion_2310.03524/digitized_width_'+self._decay_mode_filename+'.txt'
        if path.exists(width_file_name):
            width_digitized = np.loadtxt(width_file_name)
            nlines = int(width_digitized.size/2)
            m_a_list = np.array([width_digitized[i,0] for i in range(nlines)])
            Gamma_a_list = np.array([width_digitized[i,1] for i in range(nlines)])
            self._width_inter = interp1d(m_a_list, Gamma_a_list,fill_value="extrapolate")

            if not mode == ('2El' or '2Mu'): self.digitized_exists = 1 #for 2El/2Mu use analytical instead           
        else: 
            print('[Warning:] \t',width_file_name,'not found')

    def _mass_bound(self,m_a):
        return {
            'TotalHad'      : 1,
            '2Gamma'        : 1,
            '2El'           : 1 if m_a > (2*c.m_el) else 0,
            '2Mu'           : 1 if m_a > (2*c.m_mu) else 0,
            '3Pi0'          : 1 if m_a > (3*c.m_pi0) and m_a < c.m_etap else 0,
            '3Pi'           : 1 if m_a > (2*c.m_pi+c.m_pi0) and m_a < c.m_etap else 0,
            '3Eta'          : 1 if m_a > (3*c.m_eta) else 0,
            '2PiGamma'      : 1 if m_a > (2*c.m_pi) else 0,
            '2Pi0Eta'       : 1 if m_a > (2*c.m_pi0+c.m_eta) else 0,
            '2PiEta'        : 1 if m_a > (2*c.m_pi+c.m_eta) else 0,
            '2Pi0EtaPrim'   : 1 if m_a > (2*c.m_pi0+c.m_etap) else 0,
            '2PiEtaPrim'    : 1 if m_a > (2*c.m_pi+c.m_etap) else 0,
            '2PiOmega'      : 1 if m_a > (2*c.m_pi+c.m_omega) else 0,
            '2KPi0'         : 1 if m_a > (2*c.m_K+c.m_pi0) else 0,
            '2K0Pi0'        : 1 if m_a > (2*c.m_K0+c.m_pi0) else 0,
            'KK0Pi'         : 1 if m_a > (c.m_K+c.m_K0+c.m_pi) else 0,
            '2Pi2Pi0'       : 1 if m_a > (2*c.m_pi+2*c.m_pi0) else 0,
            '4Pi'           : 1 if m_a > (4*c.m_pi) else 0,
            '2Omega'        : 1 if m_a > (2*c.m_omega) else 0,
            '2Kstar'        : 1 if m_a > (2*c.m_K+2*c.m_pi) else 0,
            '2Phi'          : 1 if m_a > (4*c.m_K) else 0,
        }.get(self._decay_mode, 0) 

    def _rescale(self):
        return {
            'TotalHad'         : 1.,
            '2Gamma'        : 1.,
            '2El'           : 1.,
            '2Mu'           : 1.,
            '3Pi0'          : c.BR_Pi0_2Gamma**3,
            '3Pi'           : c.BR_Pi0_2Gamma,
            '3Eta'          : c.BR_Eta_2Gamma**3,
            '2PiGamma'      : 1.,
            '2Pi0Eta'       : c.BR_Pi0_2Gamma**2 * c.BR_Eta_2Gamma,
            '2PiEta'        : c.BR_Eta_2Gamma,
            '2Pi0EtaPrim'   : c.BR_Pi0_2Gamma**2 * c.BR_EtaPrim_2Gamma,
            '2PiEtaPrim'    : c.BR_EtaPrim_2Gamma, 
            '2PiOmega'      : 1, 
            '2KPi0'         : c.BR_Pi0_2Gamma,
            '2K0Pi0'        : c.BR_Pi0_2Gamma,
            'KK0Pi'         : 1,
            '2Pi2Pi0'       : c.BR_Pi0_2Gamma**2,
            '4Pi'           : 1.,
            '2Omega'        : 1,
            '2Kstar'        : 1.,
            '2Phi'          : 1.,           
        }.get(self._decay_mode, 0) 
    
    def width_from_file(self,m_a,g_A):
        if self._mass_bound(m_a) and g_A != 0: return self._rescale() * g_A**2 * self._width_inter(m_a) * 1e-9 * 1e+9 /4. # was in eV scale and Lambda = 1PeV
        else: return 0.

    def width_from_file_input(self,m_a,g_A):
        if self._mass_bound(m_a) and g_A != 0: return g_A**2 * self._width_inter(m_a) * 1e-9 * 1e+9 /4. # was in eV scale and Lambda = 1PeV, factor 2 difference coupling
        else: return 0.

class Total_width:  #total decay width
    def __init__(self, is_fermionic = False):
        self._is_fermionic = is_fermionic
        if is_fermionic:
            self._total = Load_digitized_fermions("TotalHad")
            self._digamma = Load_digitized_fermions("2Gamma")
        else:
            self._total = Load_digitized("total")
            self._digamma = Load_digitized("2Gamma")
            

    def width_a(self,m_a,g_gg,g_GG,g_ee=0,g_mumu=0,g_qq=0):
        if self._is_fermionic:
            return self._total.width_from_file(m_a,g_qq) + self._digamma.width_from_file(m_a,g_qq) + a_2El(m_a,g_ee) + a_2Mu(m_a,g_mumu) 

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

def a_2Mu(m_a, g_ll): #ALP -> 2mu width; factor of 2 difference in 2310.03524 coupling compared to 2110.10698
    if m_a > 2*c.m_mu: return g_ll**2. * c.m_mu**2 * np.sqrt(m_a**2 - 4 * c.m_mu**2)  / (8. * np.pi)
    return 0.

def a_2El(m_a, g_ll): #ALP -> 2ell width
    if m_a > 2*c.m_el: return g_ll**2. * c.m_el**2 * np.sqrt(m_a**2 - 4 * c.m_el**2)  / (8. * np.pi)
    return 0.

analytical_widths ={'2Gamma':   a_2Gamma,
                    '2Mu':      a_2Mu,
                    '2El':      a_2El}
import numpy as np
from scipy.interpolate import interp1d
from os import path
from ALP_rescale.general import constants as c

class DV_digitized:
    def __init__(self, mode):
        self.digitized_exists = 0
        self._decay_mode = self._decay_mode_filename = mode    
        file_name = path.dirname(path.realpath(__file__))+'/../../widths/vector/R_'+self._decay_mode_filename+'.dat'
        self.m_min = 0
        self.m_max = 0
        self.r_max = 0
        if path.exists(file_name):
            ratio_digitized = np.loadtxt(file_name)
            self.m_min = ratio_digitized[0,0]
            self.m_max = ratio_digitized[-1,0]
            self.r_max = ratio_digitized[-1,1] # use the last value for m>m_max
            nlines = int(ratio_digitized.size/2)
            m_list = np.array([ratio_digitized[i,0] for i in range(nlines)])
            ratio_list = np.array([ratio_digitized[i,1] for i in range(nlines)])
            self._ratio_inter = interp1d(m_list, ratio_list,fill_value="extrapolate")

            if not mode == ('2El' or '2Mu'): self.digitized_exists = 1 #for 2El/2Mu use analytical instead           
        else: 
            print('[Warning:] \t',file_name,'not found')

    def _mass_bound(self,m_dp):
        if m_dp < self.m_min: return 0
        return {
            'total'         : 1,
            '2El'           : 1 if m_dp > (2*c.m_el) else 0,
            '2Mu'           : 1 if m_dp > (2*c.m_mu) else 0,
            '2Pi'           : 1 if m_dp > (2*c.m_pi) else 0,
            '3Pi'           : 1 if m_dp > (2*c.m_pi+c.m_pi0) else 0,
            '4Pi'           : 1 if m_dp > (4*c.m_pi) else 0,
            '2Pi2Pi0'       : 1 if m_dp > (2*c.m_pi+2*c.m_pi0) else 0,
            '2K'            : 1 if m_dp > (2*c.m_K) else 0,
            '2KPi0'          : 1 if m_dp > (2*c.m_K+c.m_pi0) else 0,    
            'hadrons'          : 1 if m_dp > (2*c.m_pi) else 0,    
        }.get(self._decay_mode, 0) 

    def get_r(self,m_dp):
        if self._mass_bound(m_dp):
            if m_dp > self.m_max:
                return self.r_max
            return self._ratio_inter(m_dp)
        else: return 0.

class DP_widths:
    def __init__(self):
        self._ratio_2Pi = DV_digitized('2Pi')
        self._ratio_3Pi = DV_digitized('3Pi')
        self._ratio_4Pi = DV_digitized('4Pi')
        self._ratio_2Pi2Pi0 = DV_digitized('2Pi2Pi0')
        self._ratio_2K = DV_digitized('2K')
        self._ratio_2KPi0 = DV_digitized('2KPi0')
        self._ratio_hadronic = DV_digitized('hadrons')

    def dp_2mu_for_ratio(self, m_dp, epsilon): #A' -> 2mu width
        return epsilon**2. * c.alpha_EM * m_dp * np.sqrt(1 - 4 * c.m_mu**2 / m_dp**2) * (1 + 2*c.m_mu**2/m_dp**2)  / 3.
        
    def dp_2el(self, m_dp, epsilon): #A' -> 2el width
        if m_dp > 2*c.m_el: return epsilon**2. * c.alpha_EM * m_dp * np.sqrt(1 - 4 * c.m_el**2 / m_dp**2) * (1 + 2*c.m_el**2/m_dp**2)  / 3.
        return 0.

    def dp_2mu(self, m_dp, epsilon): #A' -> 2mu width
        if m_dp > 2*c.m_mu: return epsilon**2. * c.alpha_EM * m_dp * np.sqrt(1 - 4 * c.m_mu**2 / m_dp**2) * (1 + 2*c.m_mu**2/m_dp**2)  / 3.
        return 0.

    def dp_2pi(self, m_dp, epsilon): #A' -> 2pi width
        if m_dp > 2*c.m_pi: return self._ratio_2Pi.get_r(m_dp) * self.dp_2mu_for_ratio(m_dp, epsilon)
        return 0.

    def dp_3pi(self, m_dp, epsilon): #A' -> 2pi width
        if m_dp > (2*c.m_pi+c.m_pi0): return self._ratio_3Pi.get_r(m_dp) * self.dp_2mu_for_ratio(m_dp, epsilon)
        return 0.

    def dp_4pi(self, m_dp, epsilon): #A' -> 2pi width
        if m_dp > 4*c.m_pi: return self._ratio_4Pi.get_r(m_dp) * self.dp_2mu_for_ratio(m_dp, epsilon)
        return 0.

    def dp_2pi2pi0(self, m_dp, epsilon): #A' -> 2pi width
        if m_dp > (2*c.m_pi+2*c.m_pi0): return self._ratio_2Pi2Pi0.get_r(m_dp) * self.dp_2mu_for_ratio(m_dp, epsilon)
        return 0.

    def dp_2K(self, m_dp, epsilon): #A' -> 2pi width
        if m_dp > 2*c.m_K: return self._ratio_2K.get_r(m_dp) * self.dp_2mu_for_ratio(m_dp, epsilon) / 2
        return 0.

    def dp_2K_all(self, m_dp, epsilon): #A' -> 2pi width
        if m_dp > 2*c.m_K: return self._ratio_2K.get_r(m_dp) * self.dp_2mu_for_ratio(m_dp, epsilon)
        return 0.
    
    def dp_2Kpi0(self, m_dp, epsilon): #A' -> 2pi width
        if m_dp > (2*c.m_K+c.m_pi0): return self._ratio_2KPi0.get_r(m_dp) * self.dp_2mu_for_ratio(m_dp, epsilon) / 6
        return 0.
 
    def dp_2Kpi0_all(self, m_dp, epsilon): #A' -> 2pi width
        if m_dp > (2*c.m_K+c.m_pi0): return self._ratio_2KPi0.get_r(m_dp) * self.dp_2mu_for_ratio(m_dp, epsilon)
        return 0.

    def dp_hadr(self, m_dp, epsilon): #A' -> 2pi width
        if m_dp > 2*c.m_pi: return self._ratio_hadronic.get_r(m_dp) * self.dp_2mu_for_ratio(m_dp, epsilon)
        return 0.

    def dp_total(self, m_dp, epsilon):
        return self.dp_2el(m_dp, epsilon) + self.dp_2mu(m_dp, epsilon) + self.dp_hadr(m_dp, epsilon)
    
    _DP_widths ={'2El':  dp_2el,
                '2Mu':  dp_2mu,
                '2Pi':  dp_2pi,
                '3Pi':  dp_3pi,
                '4Pi':  dp_4pi,
                '2Pi2Pi0':  dp_2pi2pi0,
                '2K':  dp_2K,
                '2KPi0':  dp_2Kpi0,
                '2K_all':  dp_2K_all,
                '2KPi0_all':  dp_2Kpi0_all,
                'total':  dp_total,
                'hadrons':  dp_hadr,
                }

    def get_width(self,channel,m_dp,epsilon): return self._DP_widths[channel].__get__(self, type(self))(m_dp,epsilon)

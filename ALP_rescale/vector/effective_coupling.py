from os import path
from scipy.interpolate import interp1d
import numpy as np

class mixing:
    def __init__(self,ratio_rho,ratio_omega,ratio_phi):
        self.ratio_rho = ratio_rho
        self.ratio_omega = ratio_omega
        self.ratio_phi = ratio_phi
        
        self.__load_ratio_rho()
        self.__load_ratio_omega()
        self.__load_ratio_phi()
        return
        
    def __load_ratio_rho(self):
        file_name = path.dirname(path.realpath(__file__))+'/../../widths/vector/mixing/R_rho.dat'
        self._m_min_rho = 0
        self._m_max_rho = 0
        self._r_max_rho = 0
        if path.exists(file_name):
            ratio_digitized = np.loadtxt(file_name)
            self._m_min_rho = ratio_digitized[0,0]
            self._m_max_rho = ratio_digitized[-1,0]
            self._r_max_rho = ratio_digitized[-1,1] # use the last value for m>m_max
            nlines = int(ratio_digitized.size/2)
            m_list = np.array([ratio_digitized[i,0] for i in range(nlines)])
            ratio_list = np.array([ratio_digitized[i,1] for i in range(nlines)])
            self._ratio_inter_rho = interp1d(m_list, ratio_list,fill_value="extrapolate")

        else: 
            print('[Warning:] \t',file_name,'not found')

    def get_mixing_rho(self,m_dp):
        if m_dp > self._m_max_rho:
            return self._r_max_rho*self.ratio_rho
        elif m_dp < self._m_min_rho:
            return 0.
        else:
            return self._ratio_inter_rho(m_dp)*self.ratio_rho
        
    def __load_ratio_omega(self):
        file_name = path.dirname(path.realpath(__file__))+'/../../widths/vector/mixing/R_omega.dat'
        self._m_min_omega = 0
        self._m_max_omega = 0
        self._r_max_omega = 0
        if path.exists(file_name):
            ratio_digitized = np.loadtxt(file_name)
            self._m_min_omega = ratio_digitized[0,0]
            self._m_max_omega = ratio_digitized[-1,0]
            self._r_max_omega = ratio_digitized[-1,1] # use the last value for m>m_max
            nlines = int(ratio_digitized.size/2)
            m_list = np.array([ratio_digitized[i,0] for i in range(nlines)])
            ratio_list = np.array([ratio_digitized[i,1] for i in range(nlines)])
            self._ratio_inter_omega = interp1d(m_list, ratio_list,fill_value="extrapolate")

        else: 
            print('[Warning:] \t',file_name,'not found')

    def get_mixing_omega(self,m_dp):
        if m_dp > self._m_max_omega:
            return self._r_max_omega*self.ratio_omega
        elif m_dp < self._m_min_omega:
            return 0.
        else:
            return self._ratio_inter_omega(m_dp)*self.ratio_omega
        
    def __load_ratio_phi(self):
        file_name = path.dirname(path.realpath(__file__))+'/../../widths/vector/mixing/R_phi.dat'
        self._m_min_phi = 0
        self._m_max_phi = 0
        self._r_max_phi = 0
        if path.exists(file_name):
            ratio_digitized = np.loadtxt(file_name)
            self._m_min_phi = ratio_digitized[0,0]
            self._m_max_phi = ratio_digitized[-1,0]
            self._r_max_phi = ratio_digitized[-1,1] # use the last value for m>m_max
            nlines = int(ratio_digitized.size/2)
            m_list = np.array([ratio_digitized[i,0] for i in range(nlines)])
            ratio_list = np.array([ratio_digitized[i,1] for i in range(nlines)])
            self._ratio_inter_phi = interp1d(m_list, ratio_list,fill_value="extrapolate")

        else: 
            print('[Warning:] \t',file_name,'not found')

    def get_mixing_phi(self,m_dp):
        if m_dp > self._m_max_phi:
            return self._r_max_phi*self.ratio_phi
        elif m_dp < self._m_min_phi:
            return 0.
        else:
            return self._ratio_inter_phi(m_dp)*self.ratio_phi
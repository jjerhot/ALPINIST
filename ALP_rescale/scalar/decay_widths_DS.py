import numpy as np
import mpmath as mp #for polylog
from scipy.interpolate import interp1d
from os import path
from ALP_rescale.general import constants as c
from ALP_rescale.general import functions as f

# S production: https://arxiv.org/pdf/1809.01876.pdf

def B_K_S(m_S, y): #B -> K + S
    f_0 = 0.330/(1-m_S**2/37.46) #hep-ph/0406232
    if c.m_B > m_S + c.m_K: return mp.fabs(y)**2 * mp.fabs((3 * c.m_q[5]**2 * c.V_ts * c.V_tb)/(c.v**3 * 16 * np.pi**2 ))**2 /(16*np.pi*c.m_B**3)*(c.m_B**2 - c.m_K**2)**2 * np.sqrt(f.lambda_Kallen(c.m_B,c.m_K,m_S)) * f_0**2
    return 0.

def B_Kstar_S(m_S, y): #B -> K* + S
    A_0 = 1.364/(1-m_S**2/(c.m_B**2)) - 0.990/(1-m_S**2/36.78) #hep-ph/0412079
    if c.m_B > m_S + c.m_Kstar: return mp.fabs(y)**2 * mp.fabs((3 * c.m_q[5]**2 * c.V_ts * c.V_tb)/(c.v**3 * 16 * np.pi**2 ))**2 /(16*np.pi*c.m_B**3) * np.sqrt(f.lambda_Kallen(c.m_B,c.m_Kstar,m_S)**3) * A_0**2
    return 0.

# S decay: https://arxiv.org/pdf/1809.01876.pdf

class DS_digitized:
    def __init__(self, mode):
        self.digitized_exists = 0
        self._decay_mode = self._decay_mode_filename = mode
        if mode == ('2Pi' or '2K'):
            self.digitized_exists = 1 #for 2Gamma/2El/2Mu use analytical instead      
        gamma_mod_file_name = path.dirname(path.realpath(__file__))+'/DS_formfactor/gamma_'+self._decay_mode_filename+'_Mod.csv'
        self._gamma_mod_inter = 0
        if path.exists(gamma_mod_file_name):
            ff_digitized = np.loadtxt(gamma_mod_file_name)
            nlines = int(ff_digitized.size/2)
            m_S_list = np.array([ff_digitized[i,0] for i in range(nlines)])
            val_list = np.array([ff_digitized[i,1] for i in range(nlines)])
            self._gamma_mod_inter = interp1d(m_S_list, val_list,fill_value="extrapolate")
        else: 
            print('[Warning:] \t',gamma_mod_file_name,'not found')
        delta_mod_file_name = path.dirname(path.realpath(__file__))+'/DS_formfactor/delta_'+self._decay_mode_filename+'_Mod.csv'
        self._delta_mod_inter = 0
        if path.exists(delta_mod_file_name):
            ff_digitized = np.loadtxt(delta_mod_file_name)
            nlines = int(ff_digitized.size/2)
            m_S_list = np.array([ff_digitized[i,0] for i in range(nlines)])
            val_list = np.array([ff_digitized[i,1] for i in range(nlines)])
            self._delta_mod_inter = interp1d(m_S_list, val_list,fill_value="extrapolate")
        else: 
            print('[Warning:] \t',delta_mod_file_name,'not found')
        theta_mod_file_name = path.dirname(path.realpath(__file__))+'/DS_formfactor/theta_'+self._decay_mode_filename+'_Mod.csv'
        self._theta_mod_inter = 0
        if path.exists(theta_mod_file_name):
            ff_digitized = np.loadtxt(theta_mod_file_name)
            nlines = int(ff_digitized.size/2)
            m_S_list = np.array([ff_digitized[i,0] for i in range(nlines)])
            val_list = np.array([ff_digitized[i,1] for i in range(nlines)])
            self._theta_mod_inter = interp1d(m_S_list, val_list,fill_value="extrapolate")
        else: 
            print('[Warning:] \t',theta_mod_file_name,'not found')
        gamma_ang_file_name = path.dirname(path.realpath(__file__))+'/DS_formfactor/gamma_'+self._decay_mode_filename+'_Ang.csv'
        self._gamma_ang_inter = 0
        if path.exists(gamma_ang_file_name):
            ff_digitized = np.loadtxt(gamma_ang_file_name)
            nlines = int(ff_digitized.size/2)
            m_S_list = np.array([ff_digitized[i,0] for i in range(nlines)])
            val_list = np.array([ff_digitized[i,1] for i in range(nlines)])
            self._gamma_ang_inter = interp1d(m_S_list, val_list,fill_value="extrapolate")
        else: 
            print('[Warning:] \t',gamma_ang_file_name,'not found')
        delta_ang_file_name = path.dirname(path.realpath(__file__))+'/DS_formfactor/delta_'+self._decay_mode_filename+'_Ang.csv'
        self._delta_ang_inter = 0
        if path.exists(delta_ang_file_name):
            ff_digitized = np.loadtxt(delta_ang_file_name)
            nlines = int(ff_digitized.size/2)
            m_S_list = np.array([ff_digitized[i,0] for i in range(nlines)])
            val_list = np.array([ff_digitized[i,1] for i in range(nlines)])
            self._delta_ang_inter = interp1d(m_S_list, val_list,fill_value="extrapolate")
        else: 
            print('[Warning:] \t',delta_ang_file_name,'not found')
        theta_ang_file_name = path.dirname(path.realpath(__file__))+'/DS_formfactor/theta_'+self._decay_mode_filename+'_Ang.csv'
        self._theta_ang_inter = 0
        if path.exists(theta_ang_file_name):
            ff_digitized = np.loadtxt(theta_ang_file_name)
            nlines = int(ff_digitized.size/2)
            m_S_list = np.array([ff_digitized[i,0] for i in range(nlines)])
            val_list = np.array([ff_digitized[i,1] for i in range(nlines)])
            self._theta_ang_inter = interp1d(m_S_list, val_list,fill_value="extrapolate")
        else: 
            print('[Warning:] \t',theta_ang_file_name,'not found')

    # def gamma(self,m_S): return self._gamma_ang_inter(m_S)
    # def delta(self,m_S): return self._delta_ang_inter(m_S)
    # def theta(self,m_S): return self._theta_ang_inter(m_S)
    def gamma(self,m_S): return self._gamma_mod_inter(m_S)*np.exp(1.j*self._gamma_ang_inter(m_S))
    def delta(self,m_S): return  self._delta_mod_inter(m_S)*np.exp(1.j*self._delta_ang_inter(m_S))
    def theta(self,m_S): return  self._theta_mod_inter(m_S)*np.exp(1.j*self._theta_ang_inter(m_S))

    def S_2Pi(self, m_S, y): #S -> 2pi width
        gamma = self._gamma_mod_inter(m_S)*np.exp(1.j*self._gamma_ang_inter(m_S))
        delta = self._delta_mod_inter(m_S)*np.exp(1.j*self._delta_ang_inter(m_S))
        theta = self._theta_mod_inter(m_S)*np.exp(1.j*self._theta_ang_inter(m_S))
        if m_S > 2*c.m_pi and m_S < 2.: return y**2. * 3 * np.sqrt(1 - 4 * c.m_pi**2 / m_S**2) * np.absolute((7*gamma + 7*delta + 2*theta)/9)**2 / (32. * np.pi * m_S * c.v**2)
        return 0.

    def S_2K(self, m_S, y): #S -> 2K width
        gamma = self._gamma_mod_inter(m_S)*np.exp(1.j*self._gamma_ang_inter(m_S))
        delta = self._delta_mod_inter(m_S)*np.exp(1.j*self._delta_ang_inter(m_S))
        theta = self._theta_mod_inter(m_S)*np.exp(1.j*self._theta_ang_inter(m_S))
        if m_S > 2*c.m_K and m_S < 2.: return y**2. * np.sqrt(1 - 4 * c.m_K**2 / m_S**2) * np.absolute((7*gamma + 7*delta + 2*theta)/9)**2 / (8. * np.pi * m_S * c.v**2)
        return 0.


    # def S_2Pi(m_S, y): #S -> 2pi width
    #     p = 0.73
    #     q = 0.52
    #     omii = 1
    #     omij = 0
    #     gamma = c.m_pi**2*(omii+omij/np.sqrt(3))
    #     delta = 2/np.sqrt(3)*(c.m_K**2 - c.m_pi**2/2)*omij
    #     theta = (2*c.m_pi**2 + p*m_S**2)*omii + 2/np.sqrt(3)*(2*c.m_K**2 + q * m_S**2)*omij
    #     if m_S > 2*c.m_pi: return y**2. * 3 * np.sqrt(1 - 4 * c.m_pi**2 / m_S**2) * ((7*gamma + 7*delta + 2*theta)/9)**2 / (32. * np.pi * m_S * c.v**2)
    #     return 0.

    # def S_2K(m_S, y): #S -> 2K width
    #     p = 0.73
    #     q = 0.52
    #     omii = 1
    #     omij = 0
    #     gamma = c.m_pi**2*(np.sqrt(3)*omij+omii)/2
    #     delta = (c.m_K**2 - c.m_pi**2/2)*omii
    #     theta = np.sqrt(3)/2*(2*c.m_pi**2 + p*m_S**2)*omij + (2*c.m_K**2 + q * m_S**2)*omii
    #     if m_S > 2*c.m_K: return y**2. * np.sqrt(1 - 4 * c.m_K**2 / m_S**2) * ((7*gamma + 7*delta + 2*theta)/9)**2 / (8. * np.pi * m_S * c.v**2)
    #     return 0.


class DS_widths:
    def __init__(self):
        self._width_2Pi = DS_digitized('2Pi')
        self._width_2K = DS_digitized('2K')

    def S_2Pi(self, m_S, y):
        return self._width_2Pi.S_2Pi(m_S,y)

    def S_2K(self, m_S, y):
        return self._width_2K.S_2K(m_S,y)

    def S_below2El(self, m_S, y): #no decays possible below 2el production
        if m_S <= 0.00104713: return 10**20
        return 0.

    def S_2El(self, m_S, y): #S -> 2el width
        #if m_S > 2*c.m_el: 
        # if m_S > 0.00104713: 
        if m_S > 2*c.m_el: return y**2. * c.m_el**2 * m_S * np.sqrt(1 - 4 * c.m_el**2 / m_S**2)**3  / (8. * np.pi * c.v**2)
        return 0.

    def S_2Mu(self, m_S, y): #S -> 2mu width
        # if m_S > 2*c.m_mu:
        # if m_S > 0.217985:
        if m_S > 2*c.m_mu: return y**2. * c.m_mu**2 * m_S * np.sqrt(1 - 4 * c.m_mu**2 / m_S**2)**3  / (8. * np.pi * c.v**2)
        return 0.

    def S_2hadrest(self, m_S, y): #S -> 2eta, 2rho, 4pi, ..
        C = 5.1*1E-9
        if m_S > 4*c.m_pi and m_S < 2.: return y**2. * C * m_S**3 * np.sqrt(1 - 16*c.m_pi**2 / m_S**2)
        return 0.
    
    def S_2Pi2Pi0(self, m_S, y): #S -> 2pi+2pi0 assuming it's the only hadronic channel + accounting for symmetry factor
        return self.S_2hadrest(m_S, y)*2/3

    def S_4Pi(self, m_S, y): #S -> 4pi assuming it's the only hadronic channel + accounting for symmetry factor
        return self.S_2hadrest(m_S, y)*1/3

    def S_ss(self, m_S, y):
        m_s = 0.095
        if m_S >= 2.: return y**2. * 3. * m_s**2 * m_S * np.sqrt(1 - 4 * c.m_K**2 / m_S**2)**3  / (8. * np.pi * c.v**2)
        return 0.

    def S_cc(self, m_S, y):
        m_c = 1.3
        if m_S >= 2. and m_S > 2*c.m_D: return y**2. * 3. * m_c**2 * m_S * np.sqrt(1 - 4 * c.m_D**2 / m_S**2)**3  / (8. * np.pi * c.v**2)
        return 0.

    def S_bb(self, m_S, y):
        m_b = 4.18
        if m_S >= 2. and m_S > 2*c.m_B: return y**2. * m_b**2 * m_S * np.sqrt(1 - 4 * c.m_B**2 / m_S**2)**3  / (8. * np.pi * c.v**2)
        return 0.

    def S_2Tau(self, m_S, y): #S -> 2tau width
        # if m_S > 2*c.m_mu:
        # if m_S > 0.217985:
        if m_S > 2*c.m_tau: return y**2. * c.m_tau**2 * m_S * np.sqrt(1 - 4 * c.m_tau**2 / m_S**2)**3  / (8. * np.pi * c.v**2)
        return 0.

    def S_2gluon(self, m_S, y): #S -> 2gluon width - active above 1GeV
        def _tau_i(m_i):
            return m_S**2/(4*m_i**2)
        def _f(tau):
            if tau <= 1: return mp.asin(np.emath.sqrt(tau))**2
            return (-(np.log((1+np.sqrt(1-1/tau))/(1-np.sqrt(1-1/tau))) - 1j*np.pi)**2)/4
        def _r(m_i):
            return (_tau_i(m_i) + (_tau_i(m_i) -1)*_f(_tau_i(m_i)))/_tau_i(m_i)**2
        if m_S > 2: return (y**2 * (f.alpha_s(m_S)/1.6)**2 * m_S**3)/(32 * np.pi**3 * c.v**2) * np.absolute(sum([_r(m_q) for m_q in c.m_q]) )**2
        return 0.

    _DS_widths ={'2El':     S_2El,
                '2Mu':      S_2Mu,
                '2Pi':      S_2Pi,
                '2K':       S_2K,
                '2Pi2Pi0':  S_2Pi2Pi0,
                '4Pi':      S_4Pi,
                }

    def get_width(self,channel,m_S,y): return self._DS_widths[channel].__get__(self, type(self))(m_S,y)

    def S_total(self, m_S, y):
        return self.S_below2El(m_S, y) + self.S_2El(m_S, y) + self.S_2Mu(m_S, y) + self.S_2Pi(m_S, y) + self.S_2K(m_S, y) + self.S_2hadrest(m_S, y) + self.S_ss(m_S, y) + self.S_cc(m_S, y) + self.S_bb(m_S, y) + self.S_2Tau(m_S, y) + self.S_2gluon(m_S, y)

    def S_total_hadrons(self, m_S, y):
        return self.S_2Pi(m_S, y) + self.S_2K(m_S, y) + self.S_2hadrest(m_S, y) + self.S_ss(m_S, y) + self.S_cc(m_S, y) + self.S_bb(m_S, y) + self.S_2gluon(m_S, y)

import numpy as np
# from scipy.interpolate import RectBivariateSpline, interp1d
from os import path
import mpmath as mp #for polylog
from ALP_rescale.general import constants as c
from ALP_rescale.general import functions as f

class mixing:

    def __init__(self, Lambda):
        self._lambda_a = Lambda
        #define kinetic mixing parameters

        self._inv_mass_sum = 1./c.m_q[0] + 1./c.m_q[1] + 1./c.m_q[2]
        self._kappa_u = (1./c.m_q[0])*1./self._inv_mass_sum
        self._kappa_d = (1./c.m_q[1])*1./self._inv_mass_sum
        self._kappa_s = (1./c.m_q[2])*1./self._inv_mass_sum

    def _c_q(self, kappa, cg, cq):
        return cq + 64 * np.pi**2 * kappa * cg

    def _K_a_pi0(self,cg,cu,cd,cs):
        return - 0.5 * (self._c_q(self._kappa_u,cg,cu) - self._c_q(self._kappa_d,cg,cd))
    def _K_a_eta(self,cg,cu,cd,cs):
        return - 0.5 * ( (self._c_q(self._kappa_u,cg,cu) + self._c_q(self._kappa_d,cg,cd)) * (1/np.sqrt(3.) * np.cos(c.th_eta_etap) - np.sqrt(2/3.) * np.sin(c.th_eta_etap)) - self._c_q(self._kappa_s,cg,cs) * (2/np.sqrt(3.) * np.cos(c.th_eta_etap) + np.sqrt(2/3.) * np.sin(c.th_eta_etap))  )
    def _K_a_etap(self,cg,cu,cd,cs):
        return - 0.5 * ( (self._c_q(self._kappa_u,cg,cu) + self._c_q(self._kappa_d,cg,cd)) * (1/np.sqrt(3.) * np.sin(c.th_eta_etap) + np.sqrt(2/3.) * np.cos(c.th_eta_etap)) - self._c_q(self._kappa_s,cg,cs) * (2/np.sqrt(3.) * np.sin(c.th_eta_etap) - np.sqrt(2/3.) * np.cos(c.th_eta_etap))  )

    #define mass mixing parameters
    def _m_a_eta_Squared(self,cg):
        return - 64 * np.pi**2 * cg * np.sqrt(6.)*c.B0*np.sin(c.th_eta_etap)/self._inv_mass_sum
    def _m_a_etap_Squared(self,cg):
        return 64 * np.pi**2 * cg * np.sqrt(6.)*c.B0*np.cos(c.th_eta_etap)/self._inv_mass_sum

    def _FVMD(self,m_a): # implementation of VMD-F-function. Intermediate regime for 1.4GeV < m_a < 2GeV is based on a cubic fit
        if m_a <= 1.4:
            return 1

        elif m_a > 1.4 and m_a <= 2:
            return 5.7022 * m_a**3 - 29.4815 * m_a**2 + 49.0191 * m_a - 25.4899

        else:
            return (1.4/m_a)**4

    #define ALP-meson mixing parameters

    def pi(self, m_a, cg, cu, cd, cs):
        return self._FVMD(m_a) * f.alpha_s(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * self._K_a_pi0(cg,cu,cd,cs) * m_a**2/(m_a**2 - c.m_pi0**2)
    def eta(self, m_a, cg, cu, cd, cs):
        return self._FVMD(m_a) * f.alpha_s(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * (self._K_a_eta(cg,cu,cd,cs) * m_a**2 + self._m_a_eta_Squared(cg))/(m_a**2 - c.m_eta**2)
    def etap(self, m_a, cg, cu, cd, cs):
        return self._FVMD(m_a) * f.alpha_s(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * (self._K_a_etap(cg,cu,cd,cs) * m_a**2 + self._m_a_etap_Squared(cg))/(m_a**2 - c.m_etap**2)

class gluon_coupling:
    def __init__(self, Lambda):
        self._lambda = Lambda
    def _f(self,tau):
        if tau >= 1: return mp.asin(1/np.sqrt(tau))
        else: return np.pi/2 + 1j/2*mp.log((1+np.sqrt(1-tau))/(1-np.sqrt(1-tau)))
    def _B1(self,tau): return 1 - tau * np.abs(self._f(tau))**2
    def g_GG_eff(self, m_a, C_GG, C_qq=0): #for decays to mesons
        return (C_GG + C_qq/2*(self._B1(4*c.m_q[0]**2/m_a**2) + self._B1(4*c.m_q[1]**2/m_a**2) + self._B1(4*c.m_q[2]**2/m_a**2) + self._B1(4*c.m_q[3]**2/m_a**2) + self._B1(4*c.m_q[4]**2/m_a**2)))/self._lambda

class photon_coupling: #using (22) [1708.00443]
    def __init__(self, Lambda):
        self._lambda = Lambda
        self._th = mixing(self._lambda)
    def _f(self,tau):
        if tau >= 1: return mp.asin(1/np.sqrt(tau))
        else: return np.pi/2 + 1j/2*mp.log((1+np.sqrt(1-tau))/(1-np.sqrt(1-tau)))
    def _B1(self,tau): return 1 - tau * np.abs(self._f(tau))**2
    def _B2(self,tau): return 1 - (tau - 1) * np.abs(self._f(tau))**2
    def _chiC(self,m_a,C_GG): #chiral contribution
        if m_a > c.m_etap: return 0.
        return - 1.92 * C_GG
    def g_gg_eff(self, m_a, C_GG, C_WW, C_BB, C_ll = 0, C_qq = 0):
        C_gg_eff = C_BB + C_WW + C_WW*(2* c.alpha_EM * self._B2(4*c.m_W**2/m_a**2))/(np.pi*mp.sin(c.theta_w)**2) + self._chiC(m_a, C_GG) #contributions of gauge boson couplings to photon coupling (up to 1 loop)
        C_gg_eff = C_gg_eff + C_ll/(16 * np.pi**2)*(self._B1(4*c.m_el**2/m_a**2)+self._B1(4*c.m_mu**2/m_a**2)+self._B1(4*c.m_tau**2/m_a**2)) # + leptonic 1-loop contribution
        C_gg_eff = C_gg_eff + C_qq/(16 * np.pi**2)*(4/3*self._B1(4*c.m_q[3]**2/m_a**2)+1/3*self._B1(4*c.m_q[4]**2/m_a**2)) # + quark 1-loop contribution, b and c quarks
        C_gg_eff = C_gg_eff + self._lambda*(self._th.pi(m_a,C_GG,C_qq,C_qq,C_qq)/(16*np.pi**2 * c.fpi) + self._th.eta(m_a,C_GG,C_qq,C_qq,C_qq)/(16*np.pi**2 * c.feta) + self._th.etap(m_a,C_GG,C_qq,C_qq,C_qq)/(16*np.pi**2 * c.fetap)) # + meson mixing contributions

        return C_gg_eff / self._lambda

class lepton_coupling:
    def __init__(self, Lambda):
        self._lambda = Lambda
    def g_ee_eff(self, C_WW, C_BB, C_ll=0):
        return (C_ll - 6*c.alpha_EM**2*(3./4 * C_WW/np.sin(c.theta_w)**4 + 5./4 * C_BB/np.cos(c.theta_w)**4) * np.log(self._lambda**2/c.m_W**2) - 12*c.alpha_EM**2 * (C_WW + C_BB) * np.log(c.m_W**2/c.m_el) )/self._lambda
    def g_mumu_eff(self, C_WW, C_BB, C_ll=0):
        return (C_ll - 6*c.alpha_EM**2*(3./4 * C_WW/np.sin(c.theta_w)**4 + 5./4 * C_BB/np.cos(c.theta_w)**4) * np.log(self._lambda**2/c.m_W**2) - 12*c.alpha_EM**2 * (C_WW + C_BB) * np.log(c.m_W**2/c.m_mu) )/self._lambda

def F_loop(x):
    return (x*(1+x*(np.log(x)-1))/(1-x)**2)     

class bs_coupling:
    def __init__(self, Lambda, AA, BB):
        self._lambda = Lambda
        self._A = AA
        self._B = BB
        #define B decay branching fraction
        self._V_qb = [c.V_ub, c.V_cb, c.V_tb]
        self._V_qs = [c.V_us, c.V_cs, c.V_ts]
        self._xi_q = [c.m_q[0]**2/c.m_W**2, c.m_q[3]**2/c.m_W**2, c.m_q[5]**2/c.m_W**2]

        self._f_q_loop = [F_loop(c.m_q[0]**2/c.m_W**2), F_loop(c.m_q[3]**2/c.m_W**2), F_loop(c.m_q[5]**2/c.m_W**2)]

        self._h_B = sum([np.prod(q) for q in zip(self._V_qb, self._V_qs, self._xi_q)]) 

    def _f_B(self,mu):
        V_qb = [c.V_ub, c.V_cb, c.V_tb]
        V_qs = [c.V_us, c.V_cs, c.V_ts]
        xi_q = [c.m_q[0]**2/c.m_W**2, c.m_q[3]**2/c.m_W**2, c.m_q[5]**2/c.m_W**2]
        def bracket(xx):
            return (- 3*xx/2*mp.log(c.m_W**2/mu**2)**2
                    + (xx*(3*xx-2)*(3*xx+4)/(2 * (xx-1)**2) - (xx-2)*(3*xx+1)/(xx-1)*mp.log(xx-1))*mp.log(xx)
                    + (3*xx**3-14*xx**2-8*xx+4)/(2*(xx-1)**2)*mp.log(xx)**2
                    + (9*xx*(xx+1)/(2*(xx-1)) - xx*(3*xx*2-2*xx+8)/(xx-1)**2*mp.log(xx))*mp.log(c.m_W**2/mu**2)
                    + (np.pi**2*(4+11*xx-7*xx**2)+3*xx*(13*xx-3))/(12*(xx-1))
                    + (xx-2)*(3*xx+1)/(xx-1)*mp.polylog(2,1/xx)
                    - (xx+2)*(xx**2+2*xx-1)/(xx-1)**2*mp.polylog(2,(xx-1)/xx))
        bracket_q = map(bracket,xi_q)
        res = sum([np.prod(q) for q in zip(V_qb, V_qs, bracket_q)]) * 3 / 2 * c.CF
        return res #eq B6 2102.04474

    def _g_B(self,mu):
        V_qb = [c.V_ub, c.V_cb, c.V_tb]
        V_qs = [c.V_us, c.V_cs, c.V_ts]
        xi_q = [c.m_q[0]**2/c.m_W**2, c.m_q[3]**2/c.m_W**2, c.m_q[5]**2/c.m_W**2]
        def bracket(xx):
            return ((xx+5)/(xx-1) + 2*(xx**2-2*xx+4)/(xx-1)**2 * mp.log(xx) + 2*mp.log(c.m_W**2 / mu**2))
        bracket_q = map(bracket,xi_q)
        res = sum([np.prod(q) for q in zip(V_qb, V_qs, xi_q, bracket_q)]) / 4
        return res #eq B7 2102.04474

    def g_bs_eff(self, m_a, C_GG, C_WW): # combination of CWW and CGG
        C_aG_bs = f.alpha_s(m_a)**2 * c.alpha_W / (4*np.pi) * (self._f_B(self._lambda) + 2*c.CF*self._A*self._g_B(self._lambda) + 2*c.CF*self._B*self._h_B)

        C_aW_bs = 3 * c.alpha_W**2 * sum([np.prod(q) for q in zip(self._V_qb, self._V_qs, self._f_q_loop)]) # effective a_bs coupling

        return C_aG_bs * C_GG / self._lambda + C_aW_bs * C_WW / self._lambda

    def g_bs_fixed_all(self, C_GG, C_WW, C_BB, C_ll, C_qq): # fixed at Lambda ~ 6.3 TeV (12.6 TeV in 2110.10698 notation)
        return 1E-5 * c.V_tb*c.V_ts*(-6.1*C_GG - 2.8 * C_WW - 0.02 * C_BB + 4.15 * C_ll - 9.2 * C_qq)/self._lambda 

    def g_bs_fixed_cgg(self, C_GG): # fixed at 1TeV (based on BR from 2102.04474)
        m_a_ref = 1 #1GeV
        f_0 = 0.330/(1-m_a_ref**2/37.46)
        g_bs = C_GG/self._lambda * np.sqrt(c.Gamma_B * 2 * 0.11 * (32)**3 * np.pi**5 * c.m_B**3 / ((c.m_B**2-c.m_K**2)**2 * f_0**2 * np.sqrt(f.lambda_Kallen(c.m_B,c.m_K,m_a_ref))))
        return g_bs
    
class cu_coupling:
    def __init__(self, Lambda, AA, BB):
        self._lambda = Lambda
        self._A = AA
        self._B = BB
        #define D decay branching fraction
        self._V_kc = [c.V_cd, c.V_cs, c.V_cb]
        self._V_ku = [c.V_ud, c.V_us, c.V_ub]
        self._xi_k = [c.m_q[1]**2/c.m_W**2, c.m_q[2]**2/c.m_W**2, c.m_q[4]**2/c.m_W**2]

        self._f_k_loop = [F_loop(c.m_q[1]**2/c.m_W**2), F_loop(c.m_q[2]**2/c.m_W**2), F_loop(c.m_q[4]**2/c.m_W**2)]

        self._h_D = sum([np.prod(k) for k in zip(self._V_kc, self._V_ku, self._xi_k)]) 

    def _f_D(self,mu):
        V_kc = [c.V_cd, c.V_cs, c.V_cb]
        V_ku = [c.V_ud, c.V_us, c.V_ub]
        xi_k = [c.m_q[1]**2/c.m_W**2, c.m_q[2]**2/c.m_W**2, c.m_q[4]**2/c.m_W**2]
        def bracket(xx):
            return (- 3*xx/2*mp.log(c.m_W**2/mu**2)**2
                    + (xx*(3*xx-2)*(3*xx+4)/(2 * (xx-1)**2) - (xx-2)*(3*xx+1)/(xx-1)*mp.log(xx-1))*mp.log(xx)
                    + (3*xx**3-14*xx**2-8*xx+4)/(2*(xx-1)**2)*mp.log(xx)**2
                    + (9*xx*(xx+1)/(2*(xx-1)) - xx*(3*xx*2-2*xx+8)/(xx-1)**2*mp.log(xx))*mp.log(c.m_W**2/mu**2)
                    + (np.pi**2*(4+11*xx-7*xx**2)+3*xx*(13*xx-3))/(12*(xx-1))
                    + (xx-2)*(3*xx+1)/(xx-1)*mp.polylog(2,1/xx)
                    - (xx+2)*(xx**2+2*xx-1)/(xx-1)**2*mp.polylog(2,(xx-1)/xx))
        bracket_k = map(bracket,xi_k)
        res = sum([np.prod(k) for k in zip(V_kc, V_ku, bracket_k)]) * 3 / 2 * c.CF
        return res #eq B6 2102.04474

    def _g_D(self,mu):
        V_kc = [c.V_cd, c.V_cs, c.V_cb]
        V_ku = [c.V_ud, c.V_us, c.V_ub]
        xi_k = [c.m_q[1]**2/c.m_W**2, c.m_q[2]**2/c.m_W**2, c.m_q[4]**2/c.m_W**2]
        def bracket(xx):
            return ((xx+5)/(xx-1) + 2*(xx**2-2*xx+4)/(xx-1)**2 * mp.log(xx) + 2*mp.log(c.m_W**2 / mu**2))
        bracket_k = map(bracket,xi_k)
        res = sum([np.prod(k) for k in zip(V_kc, V_ku, xi_k, bracket_k)]) / 4
        return res #eq B7 2102.04474

    def g_cu_eff(self, m_a, C_GG, C_WW):
        C_aG_cu = f.alpha_s(m_a)**2 * c.alpha_W / (4*np.pi) * (self._f_D(self._lambda) + 2*c.CF*self._A*self._g_D(self._lambda) + 2*c.CF*self._B*self._h_D)
        C_aW_cu = 3 * c.alpha_W**2 * sum([np.prod(k) for k in zip(self._V_kc, self._V_ku, self._f_k_loop)]) 

        return C_aG_cu * C_GG / self._lambda + C_aW_cu * C_WW / self._lambda
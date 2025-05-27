import numpy as np
# from scipy.interpolate import RectBivariateSpline, interp1d
from os import path
import mpmath as mp #for polylog
from scipy.integrate import quad
from ALP_rescale.general import constants as c
from ALP_rescale.general import functions as f

class mixing: #couplings for m_a ~ 2GeV

    def __init__(self, Lambda):
        self._lambda_a = Lambda
        #define kinetic mixing parameters

        self._inv_mass_sum = 1./c.m_q[0] + 1./c.m_q[1] + 1./c.m_q[2]
        self._kappa_u = (1./c.m_q[0])*1./self._inv_mass_sum
        self._kappa_d = (1./c.m_q[1])*1./self._inv_mass_sum
        self._kappa_s = (1./c.m_q[2])*1./self._inv_mass_sum
        
        _ew = below_EW() #running from m_t to 2GeV
        self._cu_G = _ew.Cqq(2.,1,0,type="u",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cu_g = _ew.Cqq(2.,0,1,type="u",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cu_t = _ew.Cqq(2.,0,0,type="u",ct=1,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cu_b = _ew.Cqq(2.,0,0,type="u",ct=0,cel=0, cb=1, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cu_u = _ew.Cqq(2.,0,0,type="u",ct=0,cel=0, cb=0, cc=0, cs=0, cu=1, cd=0, cmu=0, ctau=0)

        self._cd_G = _ew.Cqq(2.,1,0,type="d",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cd_g = _ew.Cqq(2.,0,1,type="d",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cd_t = _ew.Cqq(2.,0,0,type="d",ct=1,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cd_b = _ew.Cqq(2.,0,0,type="d",ct=0,cel=0, cb=1, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cd_d = _ew.Cqq(2.,0,0,type="d",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=1, cmu=0, ctau=0)

        self._cs_G = _ew.Cqq(2.,1,0,type="s",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cs_g = _ew.Cqq(2.,0,1,type="s",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cs_t = _ew.Cqq(2.,0,0,type="s",ct=1,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cs_b = _ew.Cqq(2.,0,0,type="s",ct=0,cel=0, cb=1, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cs_s = _ew.Cqq(2.,0,0,type="s",ct=0,cel=0, cb=0, cc=0, cs=1, cu=0, cd=0, cmu=0, ctau=0)
        
    def _c_q(self, kappa, cg, cq):
        return cq + 64 * np.pi**2 * kappa * cg

    #kinetic mixing + term from isospin breaking (proportional to delta_I)
    def _K_a_pi0(self,cg,cu,cd,cs,m_a):
        return - 0.5 * (self._c_q(self._kappa_u,cg,cu) - self._c_q(self._kappa_d,cg,cd) - c.delta_I * c.m_pi0**2 / 3. * ( (self._c_q(self._kappa_u,cg,cu) + self._c_q(self._kappa_d,cg,cd) + 2*self._c_q(self._kappa_s,cg,cs))/(m_a**2 - c.m_etap**2) + 2*(self._c_q(self._kappa_u,cg,cu) + self._c_q(self._kappa_d,cg,cd) - self._c_q(self._kappa_s,cg,cs))/(m_a**2 - c.m_eta**2) ) )        # return - 0.5 * (self._c_q(self._kappa_d,cg,cd) - self._c_q(self._kappa_u,cg,cu) + c.delta_I * c.m_pi0**2 / 3. * ( (self._c_q(self._kappa_u,cg,cu) + self._c_q(self._kappa_d,cg,cd) + 2*self._c_q(self._kappa_s,cg,cs))/(m_a**2 - c.m_etap**2) + 2*(self._c_q(self._kappa_u,cg,cu) + self._c_q(self._kappa_d,cg,cd) - self._c_q(self._kappa_s,cg,cs))/(m_a**2 - c.m_eta**2) ) )
    def _K_a_eta(self,cg,cu,cd,cs,m_a):
        return - 0.5 * ( (self._c_q(self._kappa_u,cg,cu) + self._c_q(self._kappa_d,cg,cd) - c.delta_I * c.m_pi0**2 * (self._c_q(self._kappa_u,cg,cu) - self._c_q(self._kappa_d,cg,cd))/(m_a**2 - c.m_pi0**2) ) * (1/np.sqrt(3.) * np.cos(c.th_eta_etap) - np.sqrt(2/3.) * np.sin(c.th_eta_etap)) 
                        - self._c_q(self._kappa_s,cg,cs) * (2/np.sqrt(3.) * np.cos(c.th_eta_etap) + np.sqrt(2/3.) * np.sin(c.th_eta_etap))  )
    def _K_a_etap(self,cg,cu,cd,cs,m_a): # up to the sign?
        return - 0.5 * ( (self._c_q(self._kappa_u,cg,cu) + self._c_q(self._kappa_d,cg,cd) - c.delta_I * c.m_pi0**2 * (self._c_q(self._kappa_u,cg,cu) - self._c_q(self._kappa_d,cg,cd))/(m_a**2 - c.m_pi0**2) ) * (1/np.sqrt(3.) * np.sin(c.th_eta_etap) + np.sqrt(2/3.) * np.cos(c.th_eta_etap))
                        - self._c_q(self._kappa_s,cg,cs) * (2/np.sqrt(3.) * np.sin(c.th_eta_etap) - np.sqrt(2/3.) * np.cos(c.th_eta_etap))  )

    def _K_a_pi0_old(self,cg,cu,cd,cs,m_a):
        return - 0.5 * (self._c_q(self._kappa_u,cg,cu) - self._c_q(self._kappa_d,cg,cd) )
    def _K_a_eta_old(self,cg,cu,cd,cs,m_a):
        return - 0.5 * ( (self._c_q(self._kappa_u,cg,cu) + self._c_q(self._kappa_d,cg,cd) ) * (1/np.sqrt(3.) * np.cos(c.th_eta_etap) - np.sqrt(2/3.) * np.sin(c.th_eta_etap)) 
                        - self._c_q(self._kappa_s,cg,cs) * (2/np.sqrt(3.) * np.cos(c.th_eta_etap) + np.sqrt(2/3.) * np.sin(c.th_eta_etap))  )
    def _K_a_etap_old(self,cg,cu,cd,cs,m_a): # up to the sign?
        return - 0.5 * ( (self._c_q(self._kappa_u,cg,cu) + self._c_q(self._kappa_d,cg,cd) ) * (1/np.sqrt(3.) * np.sin(c.th_eta_etap) + np.sqrt(2/3.) * np.cos(c.th_eta_etap))
                        - self._c_q(self._kappa_s,cg,cs) * (2/np.sqrt(3.) * np.sin(c.th_eta_etap) - np.sqrt(2/3.) * np.cos(c.th_eta_etap))  )

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

    # def pi(self, m_a, cg, cu, cd, cs):
    #     return self._FVMD(m_a) * f.alpha_s(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * self._K_a_pi0(cg,cu,cd,cs,m_a) * m_a**2/(m_a**2 - c.m_pi0**2)
    # def eta(self, m_a, cg, cu, cd, cs):
    #     return self._FVMD(m_a) * f.alpha_s(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * (self._K_a_eta(cg,cu,cd,cs,m_a) * m_a**2 + self._m_a_eta_Squared(cg))/(m_a**2 - c.m_eta**2)
    # def etap(self, m_a, cg, cu, cd, cs):
    #     return self._FVMD(m_a) * f.alpha_s(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * (self._K_a_etap(cg,cu,cd,cs,m_a) * m_a**2 + self._m_a_etap_Squared(cg))/(m_a**2 - c.m_etap**2)

    def pi_old(self, m_a, cg, cu, cd, cs):
        return self._FVMD(m_a) * f.alpha_s(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * self._K_a_pi0_old(cg,cu,cd,cs,m_a) * m_a**2/(m_a**2 - c.m_pi0**2)
    def eta_old(self, m_a, cg, cu, cd, cs):
        return self._FVMD(m_a) * f.alpha_s(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * (self._K_a_eta_old(cg,cu,cd,cs,m_a) * m_a**2 + self._m_a_eta_Squared(cg))/(m_a**2 - c.m_eta**2)
    def etap_old(self, m_a, cg, cu, cd, cs):
        return self._FVMD(m_a) * f.alpha_s(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * (self._K_a_etap_old(cg,cu,cd,cs,m_a) * m_a**2 + self._m_a_etap_Squared(cg))/(m_a**2 - c.m_etap**2)

    def pi(self, m_a, cG, cg=0, cu=0, cd=0, cs=0, cb=0, ct=0): # change fermion coupling definition by 2
        cu_new = cG*self._cu_G + cg*self._cu_g + ct*self._cu_t/2 + cb*self._cu_b/2 + cu*self._cu_u/2
        cd_new = cG*self._cd_G + cg*self._cd_g + ct*self._cd_t/2 + cb*self._cd_b/2 + cd*self._cd_d/2
        cs_new = cG*self._cs_G + cg*self._cs_g + ct*self._cs_t/2 + cb*self._cs_b/2 + cs*self._cs_s/2
        
        return self._FVMD(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * self._K_a_pi0(cG,cu_new,cd_new,cs_new,m_a) * m_a**2/(m_a**2 - c.m_pi0**2)
    def eta(self, m_a, cG, cg=0, cu=0, cd=0, cs=0, cb=0, ct=0): # change fermion coupling definition by 2
        cu_new = cG*self._cu_G + cg*self._cu_g + ct*self._cu_t/2 + cb*self._cu_b/2 + cu*self._cu_u/2
        cd_new = cG*self._cd_G + cg*self._cd_g + ct*self._cd_t/2 + cb*self._cd_b/2 + cd*self._cd_d/2
        cs_new = cG*self._cs_G + cg*self._cs_g + ct*self._cs_t/2 + cb*self._cs_b/2 + cs*self._cs_s/2
        
        return self._FVMD(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * (self._K_a_eta(cG,cu_new,cd_new,cs_new,m_a) * m_a**2 + self._m_a_eta_Squared(cG))/(m_a**2 - c.m_eta**2)
    def etap(self, m_a, cG, cg=0, cu=0, cd=0, cs=0, cb=0, ct=0): # change fermion coupling definition by 2
        cu_new = cG*self._cu_G + cg*self._cu_g + ct*self._cu_t/2 + cb*self._cu_b/2 + cu*self._cu_u/2
        cd_new = cG*self._cd_G + cg*self._cd_g + ct*self._cd_t/2 + cb*self._cd_b/2 + cd*self._cd_d/2
        cs_new = cG*self._cs_G + cg*self._cs_g + ct*self._cs_t/2 + cb*self._cs_b/2 + cs*self._cs_s/2
        
        return self._FVMD(m_a) * (- c.fpi ) / ( 2 * self._lambda_a) * (self._K_a_etap(cG,cu_new,cd_new,cs_new,m_a) * m_a**2 + self._m_a_etap_Squared(cG))/(m_a**2 - c.m_etap**2)

class below_EW:
    def __init__(self, matching_scale = c.m_q[5]): #default scale = m_t
        self._lambda = matching_scale
        
    def CGG(self,mu,cG,ct=0,cb=None, cc=None, cs=None, cu=None, cd=None):
        '''Gluon coupling at scale mu; by default all top-type and bottom-type quark couplings are equal'''
        if cb is None: cb=ct
        
        if cc is None: cc=ct
        if cu is None: cu=ct
        if cs is None: cs=cb
        if cd is None: cd=cb
        
        sum = 0

        if mu<=c.m_t: sum+=ct
        if mu<=c.m_b: sum+=cb
        if mu<=c.m_c: sum+=cc
        if mu<=c.m_s: sum+=cs
        if mu<=c.m_d: sum+=cd
        if mu<=c.m_u: sum+=cu

        return cG - 1./2 * sum

    def Cgg(self,mu,cg,ct=0,cel=None, cb=None, cc=None, cs=None, cu=None, cd=None, cmu=None, ctau=None):
        '''Photon coupling at scale mu; by default all fermion couplings are equal, lepton, top-type and bottom-type can be also split'''
        if cel is None: cel=ct
        if cmu is None: cmu=cel
        if ctau is None: ctau=cel
        
        if cb is None: cb=ct
        
        if cc is None: cc=ct
        if cu is None: cu=ct
        if cs is None: cs=cb
        if cd is None: cd=cb        

        sum = 0

        if mu<=c.m_t: sum+=3*(2./3)**2*ct
        if mu<=c.m_b: sum+=3*(-1./3)**2*cb
        if mu<=c.m_c: sum+=3*(2./3)**2*cc
        if mu<=c.m_s: sum+=3*(-1./3)**2*cs
        if mu<=c.m_d: sum+=3*(-1./3)**2*cd
        if mu<=c.m_u: sum+=3*(2./3)**2*cu
        if mu<=c.m_el: sum+=cel
        if mu<=c.m_mu: sum+=cmu
        if mu<=c.m_tau: sum+=ctau

        return cg - sum

    def Cqq(self,mu,cG,cg,type,ct=None,cel=None, cb=None, cc=None, cs=None, cu=None, cd=None, cmu=None, ctau=None,nq = 6):
        '''Quark coupling at scale mu; by default all fermion couplings are equal, lepton, top-type and bottom-type can be also split'''
        quark_types = ["t","b","c","s","u","d"]
        if not type in quark_types:
            return ValueError("Quark type", type, "not matching any known value")
        charge = 0
        if type=="t":
            cq = ct
            charge = 2./3
        if type=="b":
            cq = cb
            charge = -1./3
        if type=="c":
            cq = cc
            charge = 2./3
        if type=="s":
            cq = cs
            charge = -1./3
        if type=="u":
            cq = cu
            charge = 2./3
        if type=="d":
            cq = cd
            charge = -1./3
        
        if cq is None:
            return ValueError("Coupling value not set for quark type", type)
        
        if cel is None: cel=ct
        if cmu is None: cmu=cel
        if ctau is None: ctau=cel
        
        if cb is None: cb=ct
        
        if cc is None: cc=ct
        if cu is None: cu=ct
        if cs is None: cs=cb
        if cd is None: cd=cb        
        
        if nq < 1 or nq > len(c.m_q) or mu > self._lambda: return cq

        if nq>1:
            if mu < c.m_q[nq-2]:
                return - 4*self.CGG(c.m_q[nq-1],cG,ct=ct,cb=cb, cc=cc, cs=cs, cu=cu, cd=cd)/f.beta_QCD(c.m_q[nq-1])*(f.alpha_s(c.m_q[nq-2])-f.alpha_s(c.m_q[nq-1]))/np.pi - charge**2*3*self.Cgg(c.m_q[nq-1],cg,ct=ct,cb=cb, cc=cc, cs=cs, cu=cu, cd=cd, cel=cel, cmu=cmu,ctau=ctau)/f.beta_QED(c.m_q[nq-1])*(f.alpha_EM(c.m_q[nq-2])-f.alpha_EM(c.m_q[nq-1]))/np.pi + self.Cqq(mu,cG,cg,type=type,ct=ct,cel=cel, cb=cb, cc=cc, cs=cs, cu=cu, cd=cd, cmu=cmu, ctau=ctau,nq=nq-1)
            else:
                return cq - 4*self.CGG(c.m_q[nq-1],cG,ct=ct,cb=cb, cc=cc, cs=cs, cu=cu, cd=cd)/f.beta_QCD(c.m_q[nq-1])*(f.alpha_s(mu)-f.alpha_s(c.m_q[nq-1]))/np.pi - charge**2*3*self.Cgg(c.m_q[nq-1],cg,ct=ct,cb=cb, cc=cc, cs=cs, cu=cu, cd=cd, cel=cel, cmu=cmu,ctau=ctau)/f.beta_QED(c.m_q[nq-1])*(f.alpha_EM(mu)-f.alpha_EM(c.m_q[nq-1]))/np.pi 
        else:
            return cq - 4*self.CGG(c.m_q[nq-1],cG,ct=ct,cb=cb, cc=cc, cs=cs, cu=cu, cd=cd)/f.beta_QCD(c.m_q[nq-1])*(f.alpha_s(mu)-f.alpha_s(c.m_q[nq-1]))/np.pi - charge**2*3*self.Cgg(c.m_q[nq-1],cg,ct=ct,cb=cb, cc=cc, cs=cs, cu=cu, cd=cd, cel=cel, cmu=cmu,ctau=ctau)/f.beta_QED(c.m_q[nq-1])*(f.alpha_EM(mu)-f.alpha_EM(c.m_q[nq-1]))/np.pi 

class above_EW:
    def __init__(self, matching_scale = 1000.): #default scale = 1TeV
        self._lambda = matching_scale

    def C_GG_tilde(self, C_GG, C_u=[0,0,0], C_d=[0,0,0], C_Q=[0,0,0]): # fermion couplings are ndarrays
        if len(C_u) != 3 or  len(C_d) != 3 or  len(C_Q) != 3:
            raise ValueError("There are 3 generations of fermions, C_u,C_d,C_Q must be an array with 3 elements")
        return C_GG + 1./2 * (np.array(C_u) + np.array(C_d) - 2*np.array(C_Q)).sum()
    
    def C_WW_tilde(self, C_WW, C_Q=[0,0,0], C_L=[0,0,0]): # fermion couplings are ndarrays
        if len(C_Q) != 3 or  len(C_L) != 3:
            raise ValueError("There are 3 generations of fermions, C_Q,C_L must be an array with 3 elements")
        return C_WW - 1./2 * (3*np.array(C_Q) + np.array(C_L)).sum()
    
    def C_BB_tilde(self, C_BB, C_u=[0,0,0], C_d=[0,0,0], C_Q=[0,0,0], C_e=[0,0,0], C_L=[0,0,0]): # fermion couplings are ndarrays
        if len(C_u) != 3 or len(C_d) != 3 or len(C_Q) != 3 or len(C_e) != 3 or len(C_L) != 3:
            raise ValueError("There are 3 generations of fermions, C_u,C_d,C_Q,C_e,C_L must be an array with 3 elements")
        return C_BB + (4./3*np.array(C_u) + 1./3*np.array(C_d) - 1./6*np.array(C_Q) + np.array(C_e) - 1./2*np.array(C_L)).sum()
    
    def __one_minus_e18U(self,mu,lam):
        return 9./2 * f.alpha_t(mu)/f.alpha_s(mu) * (1-(f.alpha_s(lam)/f.alpha_s(mu))**(1/7))
    
    def C_GG_tilde_mu(self, mu, C_GG_tilde, C_tt):
        return C_GG_tilde - 2./9 * self.__one_minus_e18U(mu,self._lambda) * C_tt
    
    def C_WW_tilde_mu(self, mu, C_WW_tilde, C_tt):
        return C_WW_tilde - 1./6 * self.__one_minus_e18U(mu,self._lambda) * C_tt
    
    def C_BB_tilde_mu(self, mu, C_BB_tilde, C_tt):
        return C_BB_tilde - 17./54 * self.__one_minus_e18U(mu,self._lambda) * C_tt
    
    def __I_t(self,mu_w,C_GG_tilde,C_WW_tilde,C_BB_tilde,C_tt):
        def integrand_Lambda_mu(mu):
            return 1./mu * self.__one_minus_e18U(mu_w,mu)\
                * (2.*f.alpha_s(mu)**2/ (np.pi**2) * self.C_GG_tilde_mu(mu, C_GG_tilde, C_tt) \
                 + 9.*f.alpha_ww(mu)**2 / (16.*np.pi**2) * self.C_WW_tilde_mu(mu, C_WW_tilde, C_tt) \
                 + 17.*f.alpha_bb(mu)**2 / (48.*np.pi**2) * self.C_BB_tilde_mu(mu, C_BB_tilde, C_tt))
        integral_Lambda_mu = quad(integrand_Lambda_mu, self._lambda, mu_w)
        return - 2./3 * self.__one_minus_e18U(mu_w,self._lambda) * C_tt - 2./3 * integral_Lambda_mu[0]
    
    def C_uu_mu(self,mu_w,C_GG_tilde,C_WW_tilde,C_BB_tilde,C_tt, C_uu = None, type = None):
        quark_types = ["t","c","u"]
        if not type in quark_types:
            return ValueError("Quark type", type, "not matching any known value")
        delta = 0.
        if type == "t": delta = 1
        if C_uu is None: C_uu = C_tt
        def integrand_Lambda_mu(mu):
            return 1./mu * (1-2./3*(1+delta/2)*self.__one_minus_e18U(mu_w,mu))\
                * (2.*f.alpha_s(mu)**2/ np.pi**2 * self.C_GG_tilde_mu(mu, C_GG_tilde, C_tt) \
                 + 9.*f.alpha_ww(mu)**2 / (16.*np.pi**2) * self.C_WW_tilde_mu(mu, C_WW_tilde, C_tt) \
                 + 17.*f.alpha_bb(mu)**2 / (48.*np.pi**2) * self.C_BB_tilde_mu(mu, C_BB_tilde, C_tt))
        integral_Lambda_mu = quad(integrand_Lambda_mu, self._lambda, mu_w)
        return C_uu - 2./3*(1+delta/2)*self.__one_minus_e18U(mu_w,self._lambda)*C_tt + integral_Lambda_mu[0]

    def C_dd_mu(self,mu_w,C_GG_tilde,C_WW_tilde,C_BB_tilde,C_tt, C_dd = None, type = None):
        quark_types = ["b","s","d"]
        if not type in quark_types:
            return ValueError("Quark type", type, "not matching any known value")
        V_CKM = 0
        if type=="b": V_CKM = c.V_tb
        elif type=="s": V_CKM = c.V_ts
        elif type=="d": V_CKM = c.V_td
        if C_dd is None: C_dd = C_tt
        
        def integrand_Lambda_mu(mu):
            return 1./mu \
                * (2.*f.alpha_s(mu)**2/ np.pi**2 * self.C_GG_tilde_mu(mu, C_GG_tilde, C_tt) \
                 + 9.*f.alpha_ww(mu)**2 / (16.*np.pi**2) * self.C_WW_tilde_mu(mu, C_WW_tilde, C_tt) \
                 + 5.*f.alpha_bb(mu)**2 / (48.*np.pi**2) * self.C_BB_tilde_mu(mu, C_BB_tilde, C_tt))
        integral_Lambda_mu = quad(integrand_Lambda_mu, self._lambda, mu_w)
        return C_dd + (- 1 + V_CKM**2/6)*self.__I_t(mu_w,C_GG_tilde,C_WW_tilde,C_BB_tilde,C_tt) + integral_Lambda_mu[0]

    def C_ee_mu(self,mu_w,C_GG_tilde,C_WW_tilde,C_BB_tilde,C_tt, C_ee = None):
        if C_ee is None: C_ee = C_tt
        def integrand_Lambda_mu(mu):
            return 1./mu * (1-self.__one_minus_e18U(mu_w,mu))\
                * ( 9.*f.alpha_ww(mu)**2 / (16.*np.pi**2) * self.C_WW_tilde_mu(mu, C_WW_tilde, C_tt) \
                 + 15.*f.alpha_bb(mu)**2 / (16.*np.pi**2) * self.C_BB_tilde_mu(mu, C_BB_tilde, C_tt))
        integral_Lambda_mu = quad(integrand_Lambda_mu, self._lambda, mu_w)
        return C_ee - self.__I_t(mu_w,C_GG_tilde,C_WW_tilde,C_BB_tilde,C_tt) + integral_Lambda_mu[0]
    
    def C_qq_ij_mu(self,mu_w,C_GG_tilde,C_WW_tilde,C_BB_tilde,C_tt, C_qq_ij = 0, i = None, j = None):
        top_type = {"t":c.V_tb,"c":c.V_cb,"u":c.V_ub}
        bottom_type = {"b":c.V_tb,"s":c.V_ts,"d":c.V_td}
        if i in top_type.keys() and j in bottom_type.keys() or i in bottom_type.keys() and j in top_type.keys():
            return ValueError("Both quarks have to be either of the top-type or bottom-type")
        types = top_type
        types.update(bottom_type)
        if not i in list(types.keys()):
            return ValueError("Quark type i", i, "not matching any known value")
        if not j in list(types.keys()):
            return ValueError("Quark type j", j, "not matching any known value")
        if i==j:
            return ValueError("Quarks have to be of a different flavor")
        V_i = types[i]
        V_j = types[j]
        x_t = c.m_t**2/c.m_W**2
        return C_qq_ij + V_i*V_j*(- 1./6*self.__I_t(mu_w,C_GG_tilde,C_WW_tilde,C_BB_tilde,C_tt) + f.alpha_t(mu_w)/(4*np.pi)*(self.C_uu_mu(mu_w,C_GG_tilde,C_WW_tilde,C_BB_tilde,C_tt,type="t")*( np.log(mu_w**2/c.m_t**2)/2 - 1./4 - 3./2*( 1 - x_t + np.log(x_t))/(1-x_t)**2) - 3*f.alpha_EM(mu_w)/(2*np.pi*c.sin_w2)*self.C_WW_tilde_mu(mu_w, C_WW_tilde, C_tt)*( 1 - x_t + x_t*np.log(x_t))/(1-x_t)**2 ))
    
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
        
        _ew = below_EW() #running from m_t to 2GeV
        self._cu_G = _ew.Cqq(2.,1,0,type="u",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cu_g = _ew.Cqq(2.,0,1,type="u",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cu_t = _ew.Cqq(2.,0,0,type="u",ct=1,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cu_b = _ew.Cqq(2.,0,0,type="u",ct=0,cel=0, cb=1, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cu_u = _ew.Cqq(2.,0,0,type="u",ct=0,cel=0, cb=0, cc=0, cs=0, cu=1, cd=0, cmu=0, ctau=0)

        self._cd_G = _ew.Cqq(2.,1,0,type="d",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cd_g = _ew.Cqq(2.,0,1,type="d",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cd_t = _ew.Cqq(2.,0,0,type="d",ct=1,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cd_b = _ew.Cqq(2.,0,0,type="d",ct=0,cel=0, cb=1, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cd_d = _ew.Cqq(2.,0,0,type="d",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=1, cmu=0, ctau=0)

        self._cs_G = _ew.Cqq(2.,1,0,type="s",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cs_g = _ew.Cqq(2.,0,1,type="s",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cs_t = _ew.Cqq(2.,0,0,type="s",ct=1,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cs_b = _ew.Cqq(2.,0,0,type="s",ct=0,cel=0, cb=1, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cs_s = _ew.Cqq(2.,0,0,type="s",ct=0,cel=0, cb=0, cc=0, cs=1, cu=0, cd=0, cmu=0, ctau=0)
        
        self._cc_G = _ew.Cqq(2.,1,0,type="c",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cc_g = _ew.Cqq(2.,0,1,type="c",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cc_t = _ew.Cqq(2.,0,0,type="c",ct=1,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cc_b = _ew.Cqq(2.,0,0,type="c",ct=0,cel=0, cb=1, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cc_c = _ew.Cqq(2.,0,0,type="c",ct=0,cel=0, cb=0, cc=1, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        
        self._cb_G = _ew.Cqq(2.,1,0,type="b",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cb_g = _ew.Cqq(2.,0,1,type="b",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cb_t = _ew.Cqq(2.,0,0,type="b",ct=1,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._cb_b = _ew.Cqq(2.,0,0,type="b",ct=0,cel=0, cb=1, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        
        self._ct_G = _ew.Cqq(2.,1,0,type="t",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._ct_g = _ew.Cqq(2.,0,1,type="t",ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._ct_t = _ew.Cqq(2.,0,0,type="t",ct=1,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        self._ct_b = _ew.Cqq(2.,0,0,type="t",ct=0,cel=0, cb=1, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0)
        
    def _f(self,tau):
        if tau >= 1: return mp.asin(1/np.sqrt(tau))
        else: return np.pi/2 + 1j/2*mp.log((1+np.sqrt(1-tau))/(1-np.sqrt(1-tau)))
    def _B1(self,tau): return 1 - tau * np.abs(self._f(tau))**2
    def _B2(self,tau): return 1 - (tau - 1) * np.abs(self._f(tau))**2
    def _chiC(self,m_a,C_GG): #chiral contribution
        if m_a > c.m_etap: return 0.
        return - 1.92 * C_GG
    def g_gg_eff(self, m_a, C_GG, C_WW, C_BB, ct=0,cel=0, cb=0, cc=0, cs=0, cu=0, cd=0, cmu=0, ctau=0):
        cu_new = C_GG*self._cu_G + (C_BB+C_WW)*self._cu_g + ct*self._cu_t + cb*self._cu_b + cu*self._cu_u
        cd_new = C_GG*self._cd_G + (C_BB+C_WW)*self._cd_g + ct*self._cd_t + cb*self._cd_b + cd*self._cd_d
        cs_new = C_GG*self._cs_G + (C_BB+C_WW)*self._cs_g + ct*self._cs_t + cb*self._cs_b + cs*self._cs_s
        cc_new = C_GG*self._cc_G + (C_BB+C_WW)*self._cc_g + ct*self._cc_t + cb*self._cc_b + cc*self._cc_c
        cb_new = C_GG*self._cb_G + (C_BB+C_WW)*self._cb_g + ct*self._cb_t + cb*self._cb_b 
        ct_new = C_GG*self._ct_G + (C_BB+C_WW)*self._ct_g + ct*self._ct_t + cb*self._ct_b 
        
        # print(cu_new,cd_new,cs_new,cc_new,cb_new,ct_new)
        
        C_gg_eff = C_BB + C_WW + C_WW*(2* c.alpha_EM * self._B2(4*c.m_W**2/m_a**2))/(np.pi*mp.sin(c.theta_w)**2) + self._chiC(m_a, C_GG) #contributions of gauge boson couplings to photon coupling (up to 1 loop)
        C_gg_eff = C_gg_eff + 1./(16 * np.pi**2)*(cel*self._B1(4*c.m_el**2/m_a**2) + cmu*self._B1(4*c.m_mu**2/m_a**2) + ctau*self._B1(4*c.m_tau**2/m_a**2)) # + leptonic 1-loop contribution
        C_gg_eff = C_gg_eff + 1./(16 * np.pi**2)*(cu_new*4/3*self._B1(4*c.m_q[0]**2/m_a**2) + cd_new*1/3*self._B1(4*c.m_q[1]**2/m_a**2) + cs_new*1/3*self._B1(4*c.m_q[2]**2/m_a**2) + cc_new*4/3*self._B1(4*c.m_q[3]**2/m_a**2) + cb_new*1/3*self._B1(4*c.m_q[4]**2/m_a**2) + ct_new*4/3*self._B1(4*c.m_q[5]**2/m_a**2)) # + quark 1-loop contribution
        C_gg_eff = C_gg_eff + self._lambda*(self._th.pi(m_a,C_GG,cu,cd,cs)/(16*np.pi**2 * c.fpi) + self._th.eta(m_a,C_GG,cu,cd,cs)/(16*np.pi**2 * c.feta) + self._th.etap(m_a,C_GG,cu,cd,cs)/(16*np.pi**2 * c.fetap)) # + meson mixing contributions

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
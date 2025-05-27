import numpy as np
from scipy.interpolate import interp1d
from os import path
from ALP_rescale.general import constants as c
from ALP_rescale.general import functions as f

s_interp = np.linspace((c.m_pi0+c.m_pi)**2,5.31**2,2000)
Fabs2_pi = interp1d(s_interp, np.abs(f.F_pi_analytic(s_interp))**2, bounds_error=False, fill_value=0.)
del s_interp
class hnl_decay_width:
	"""decay widths of the hnl as taken from JHEP, 0710
	"""
	@staticmethod
	def widths(decay_mode, hnl_mass, U2 = [1.,1.,1.], dirac_HNL = False, only_charged_current = False):
		if dirac_HNL: return hnl_decay_width.dirac_widths(decay_mode, hnl_mass, U2, only_charged_current=only_charged_current)
		return 2. * hnl_decay_width.dirac_widths(decay_mode, hnl_mass, U2, only_charged_current=only_charged_current)
	@staticmethod
	def dirac_widths(decay_mode, hnl_mass, U2 = [1.,1.,1.], only_charged_current=False):
		U2_el, U2_mu, U2_tau = U2

		if   decay_mode == "PiNu":  return sum(U2) * hnl_decay_width.N_NuH(hnl_mass, c.m_pi0, 1., c.f_pi,1.) if not only_charged_current else 0.
		elif decay_mode == "EtaNu": return sum(U2) * hnl_decay_width.N_NuH(hnl_mass, c.m_eta, 1., c.f_eta,   1.) if not only_charged_current else 0.
		elif decay_mode == "EtapNu":return sum(U2) * hnl_decay_width.N_NuH(hnl_mass, c.m_etap,1., c.f_etap,	 1.) if not only_charged_current else 0.

		elif decay_mode == "PiEl":	return hnl_decay_width.N_LH(hnl_mass, c.m_el, c.m_pi, 1., c.f_piV_ud,  U2_el)
		elif decay_mode == "PiMu": 	return hnl_decay_width.N_LH(hnl_mass, c.m_mu, c.m_pi, 1., c.f_piV_ud,  U2_mu)
		elif decay_mode == "PiTau":	return hnl_decay_width.N_LH(hnl_mass, c.m_tau,c.m_pi, 1., c.f_piV_ud,  U2_tau)
		elif decay_mode == "KEl":	return hnl_decay_width.N_LH(hnl_mass, c.m_el, c.m_K,  1., c.f_KV_us,   U2_el)
		elif decay_mode == "KMu": 	return hnl_decay_width.N_LH(hnl_mass, c.m_mu, c.m_K,  1., c.f_KV_us,   U2_mu)
		elif decay_mode == "KTau": 	return hnl_decay_width.N_LH(hnl_mass, c.m_tau,c.m_K,  1., c.f_KV_us,   U2_tau)
		elif decay_mode == "DEl":	return hnl_decay_width.N_LH(hnl_mass, c.m_el, c.m_D,  1.,  c.f_DV_cd,  U2_el)
		elif decay_mode == "DMu": 	return hnl_decay_width.N_LH(hnl_mass, c.m_mu, c.m_D,  1.,  c.f_DV_cd,  U2_mu)
		elif decay_mode == "DTau": 	return hnl_decay_width.N_LH(hnl_mass, c.m_tau,c.m_D,  1.,  c.f_DV_cd,  U2_tau)
		elif decay_mode == "DsEl":	return hnl_decay_width.N_LH(hnl_mass, c.m_el, c.m_Ds, 1.,  c.f_DsV_cs, U2_el)
		elif decay_mode == "DsMu": 	return hnl_decay_width.N_LH(hnl_mass, c.m_mu, c.m_Ds, 1.,  c.f_DsV_cs, U2_mu)
		elif decay_mode == "DsTau":	return hnl_decay_width.N_LH(hnl_mass, c.m_tau,c.m_Ds, 1.,  c.f_DsV_cs, U2_tau)
		elif decay_mode == "DsstarEl":return  hnl_decay_width.N_LHv(hnl_mass, c.m_el, c.m_Dsstar, c.V_cs, c.g_Dsstar, U2_el)
		elif decay_mode == "DsstarMu":return  hnl_decay_width.N_LHv(hnl_mass, c.m_mu, c.m_Dsstar, c.V_cs, c.g_Dsstar, U2_mu)
		elif decay_mode == "DsstarTau":return hnl_decay_width.N_LHv(hnl_mass, c.m_tau,c.m_Dsstar, c.V_cs, c.g_Dsstar, U2_tau)

		elif decay_mode == "OmegaNu":	return sum(U2) * hnl_decay_width.N_NuHv(hnl_mass, c.m_omega, c.g_omega, 2./3*c.sin_w2,    1.) if not only_charged_current else 0.
		elif decay_mode == "PhiNu": 	return sum(U2) * hnl_decay_width.N_NuHv(hnl_mass, c.m_phi,   c.g_phi,   (4./3*c.sin_w2-1)/np.sqrt(2), 1.) if not only_charged_current else 0.
		elif decay_mode == "JpsiNu": 	return sum(U2) * hnl_decay_width.N_NuHv(hnl_mass, c.m_Jpsi,  c.g_Jpsi,  (1.-8./3*c.sin_w2)/np.sqrt(2),  1.) if not only_charged_current else 0.

		elif decay_mode == "NuNuNu": 	return sum(U2) * hnl_decay_width.N_3Nus(hnl_mass, 1) if not only_charged_current else 0.

		elif decay_mode == "NuElEl":	return (hnl_decay_width.N_NuaLaLa(hnl_mass, c.m_el, U2_el)   + (U2_mu+U2_tau)*hnl_decay_width.N_NuaLbLb(hnl_mass, c.m_el, 1.)) if not only_charged_current else hnl_decay_width.N_LaLbNub(hnl_mass, c.m_el, c.m_el,  U2_el)
		elif decay_mode == "NuMuMu":	return (hnl_decay_width.N_NuaLaLa(hnl_mass, c.m_mu, U2_mu)   + (U2_el+U2_tau)*hnl_decay_width.N_NuaLbLb(hnl_mass, c.m_mu, 1.)) if not only_charged_current else hnl_decay_width.N_LaLbNub(hnl_mass, c.m_mu, c.m_mu,  U2_mu)
		elif decay_mode == "NuTauTau":	return (hnl_decay_width.N_NuaLaLa(hnl_mass, c.m_tau, U2_tau) + (U2_el+U2_mu)*hnl_decay_width.N_NuaLbLb(hnl_mass, c.m_tau, 1.)) if not only_charged_current else hnl_decay_width.N_LaLbNub(hnl_mass, c.m_tau, c.m_tau,  U2_tau)

		elif decay_mode == "NuElMu":	return hnl_decay_width.N_LaLbNub(hnl_mass, c.m_el, c.m_mu,  U2_el) + hnl_decay_width.N_LaLbNub(hnl_mass, c.m_mu,  c.m_el, U2_mu)
		elif decay_mode == "NuElTau":	return hnl_decay_width.N_LaLbNub(hnl_mass, c.m_el, c.m_tau, U2_el) + hnl_decay_width.N_LaLbNub(hnl_mass, c.m_tau, c.m_el, U2_tau)
		elif decay_mode == "NuMuTau":	return hnl_decay_width.N_LaLbNub(hnl_mass, c.m_mu, c.m_tau, U2_mu) + hnl_decay_width.N_LaLbNub(hnl_mass, c.m_tau, c.m_mu, U2_tau)

		# pipi states replace rho resonance approximation
		elif decay_mode == "PiPiNu":  return sum(U2) * hnl_decay_width.N_NuPiPi(hnl_mass, 1.) if not only_charged_current else 0.
		elif decay_mode == "PiPiEl":  return hnl_decay_width.N_LPiPi(hnl_mass, c.m_el,  U2_el) 
		elif decay_mode == "PiPiMu":  return hnl_decay_width.N_LPiPi(hnl_mass, c.m_mu,  U2_mu) 
		elif decay_mode == "PiPiTau": return hnl_decay_width.N_LPiPi(hnl_mass, c.m_tau, U2_tau) 

		# elif decay_mode == "RhoNu":		return sum(U2) * hnl_decay_width.N_NuHv(hnl_mass, c.m_rho,  c.g_rho,   1.-2.*c.sin_w2,   1.)  if not only_charged_current else 0.
		# elif decay_mode == "RhoEl":	return hnl_decay_width.N_LHv(hnl_mass, c.m_el, c.m_rho, c.V_ud, c.g_rho, U2_el)
		# elif decay_mode == "RhoMu":	return hnl_decay_width.N_LHv(hnl_mass, c.m_mu, c.m_rho, c.V_ud, c.g_rho, U2_mu)
		# elif decay_mode == "RhoTau":return hnl_decay_width.N_LHv(hnl_mass, c.m_tau,c.m_rho, c.V_ud, c.g_rho, U2_tau)  
		else: print("[Info:] HNL branching for decay mode ",decay_mode,"not found.")
		return 0.
	@staticmethod 
	def total_width(hnl_mass, U2, dirac_HNL = False, only_charged_current=False):
		semileptonic_width = hnl_decay_width.semileptonic_width_single_meson(hnl_mass, U2, dirac_HNL,only_charged_current)
		if 1. < hnl_mass: semileptonic_width = max(hnl_decay_width.semileptonic_width_approx(hnl_mass,U2,dirac_HNL,only_charged_current), semileptonic_width)
		return hnl_decay_width.leptonic_width(hnl_mass,U2,dirac_HNL,only_charged_current) + semileptonic_width
	@staticmethod
	def N_NuH(m_N, m_H, V_H, f_H, U2_l):
		if  m_H < m_N : return U2_l * np.divide(c.G_F**2 * f_H**2 * m_N**3 * V_H**2, 32*np.pi) * (1 - (m_H/m_N)**2)**2 
		return 0.
	@staticmethod
	def N_NuHv(m_N, m_Hv, g_Hv, kappa_Hv, U2_l):
		if  m_Hv < m_N :
			x2 = (m_Hv/m_N)**2
			return U2_l * np.divide(c.G_F**2 * g_Hv**2 * kappa_Hv**2 * m_N**3, 32*np.pi* m_Hv**2) * (1. + 2*x2) * (1. - x2)**2 
		return 0.
	@staticmethod
	def N_LH(m_N, m_l, m_H, V_H, f_H, U2_l):
		if m_l + m_H < m_N : 
			xl = m_l/m_N
			xh = m_H/m_N
			return U2_l * np.divide(c.G_F**2 * f_H**2 * V_H**2 * m_N**3 , 16.*np.pi) * ( (1-xl**2)**2 - xh**2*(1+xl**2))*np.sqrt(f.lambda_Kallen(1,xh,xl))
		return 0.
	@staticmethod
	def  N_LHv(m_N, m_l, m_Hv, V_Hv, g_Hv, U2_l ):
		if m_l + m_Hv < m_N : 
			xl = m_l/m_N
			xh = m_Hv/m_N
			return U2_l * np.divide(c.G_F**2 * g_Hv**2 * m_N**3 * V_Hv**2, 16*np.pi*m_Hv**2) * ((1-xl**2)**2 + xh**2*(1+xl**2) - 2*xh**4)*np.sqrt(f.lambda_Kallen(1.,xh,xl))
		return 0.
	@staticmethod
	def N_3Nus(m_N,  U2):
		return U2 * m_N**5 * np.divide(c.G_F**2, 192*np.pi**3)
	@staticmethod
	def N_LaLbNub(m_N, m_la, m_lb,  U2_la):
		if m_la + m_lb > m_N: return 0.
		x2 = max(m_la, m_lb)**2/m_N**2
		return U2_la * m_N**5 * np.divide(c.G_F**2, 192*np.pi**3) * (1 - 8*x2 + 8*x2**3 - x2**4 - 12*x2**2*np.log(x2))
	@staticmethod
	def N_NuaLbLb(m_N, m_lb, U2_la):
		if 2 * m_lb > m_N: return 0.
		x2 = (m_lb / m_N)**2
		x4 = x2**2
		L =  np.log(np.where(x2>1e-3, np.divide(1. - 3.*x2 - (1.-x2)*np.sqrt(1. - 4.*x2), x2*(1.+np.sqrt(1.-4.*x2))), x4)) #Taylor expansion to prevent log errors
		C1 = 0.25 * (1 - 4* c.sin_w2 + 8*c.sin_w2**2)
		C2 = 0.5 * c.sin_w2 * (2*c.sin_w2 - 1)
		return U2_la * m_N**5 * np.divide(c.G_F**2, 192*np.pi**3) *  (
			C1 *((1-14*x2 - 2*x4-12*x2**3)*np.sqrt(1-4*x2) + 12*x4 * (x4-1)*L) + 
			4 * C2 * (x2*(2+10*x2 - 12*x4)*np.sqrt(1-4*x2) + 6*x4*(1-2*x2+2*x4)*L)   )
	@staticmethod
	def N_NuaLaLa(m_N, m_l, U2_la):
		if 2 * m_l < m_N: 
			x2 = m_l**2 / m_N**2
			x4 = x2**2
			L =  np.log(np.where(x2>1e-3, np.divide(1. - 3.*x2 - (1.-x2)*np.sqrt(1. - 4.*x2), x2*(1.+np.sqrt(1.-4.*x2))), x4)) 
			C3 = 0.25 * (1 + 4* c.sin_w2 + 8*c.sin_w2**2)
			C4 = 0.5 * c.sin_w2 * (2*c.sin_w2 + 1)
			return U2_la * m_N**5 * np.divide(c.G_F**2, 192*np.pi**3) *  (
				C3 *((1-14*x2 - 2*x4-12*x2**3)*np.sqrt(1-4*x2) + 12*x4 * (x4-1)*L) + 4 * C4 * (
				x2*(2+10*x2 - 12*x4)*np.sqrt(1-4*x2) + 6*x4*(1-2*x2+2*x4)*L)
			)
		return 0.
	@staticmethod
	def N_NuPiPi(m_N, U2_la, n_intpoints = 1000):
		if 2*c.m_pi < m_N:
			M2_Zs = np.linspace(4*c.m_pi**2, m_N**2, n_intpoints+2, dtype=np.double)[1:-1]
			x_Z = np.sqrt(M2_Zs) / m_N 
			Integrant =   (1. - x_Z**2)**2 * (1. + x_Z**2) * np.power(1 -4.*c.m_pi**2 / M2_Zs, 1.5) *  Fabs2_pi(M2_Zs)
			return U2_la * np.divide(c.G_F**2*m_N**3 *c.V_ud**2, 768.*np.pi**3) * (1-2*c.sin_w2**2)**2  * np.trapz(Integrant, x=M2_Zs)
		return 0.
	@staticmethod
	def N_LPiPi(m_N, m_l, U2_la, n_intpoints = 1000):
		if c.m_pi + c.m_pi0 + m_l < m_N:
			M2_Ws = np.linspace((c.m_pi0+c.m_pi)**2, (m_N - m_l)**2, n_intpoints+2, dtype=np.double)[1:-1]
			x_l = m_l / m_N
			x_W = np.sqrt(M2_Ws) / m_N 
			Integrant =  ((1-x_l**2)**2 + x_W**2 *(1+x_l**2) - 2*x_W**4) * np.sqrt(f.lambda_Kallen(1, x_W, x_l)) * np.power(1 - (c.m_pi+c.m_pi0)**2 / M2_Ws, 1.5) *  Fabs2_pi(M2_Ws)
			return U2_la * np.divide(c.G_F**2*m_N**3 *c.V_ud**2, 384.*np.pi**3) * np.trapz(Integrant, x=M2_Ws)
		return 0.
	# misc widths for hadronic width approximation:
	@staticmethod
	def leptonic_width(m_N, U2, dirac_HNL=False, only_charged_current=False):
		Gamma = 0.
		decay_channels = ["NuNuNu", "NuElEl", "NuMuMu", "NuTauTau", "NuElMu", "NuElTau", "NuMuTau"]
		for channel in decay_channels: Gamma += hnl_decay_width.dirac_widths(channel, m_N, U2, only_charged_current=only_charged_current)
		return Gamma if dirac_HNL else 2.*Gamma
	@staticmethod
	def semileptonic_width_single_meson(hnl_mass, U2, dirac_HNL=False, only_charged_current = False):
		Gamma = 0.
		decay_channels = ["PiNu", "EtaNu", "EtapNu", "PiPiNu", "OmegaNu", "PhiNu", "PiEl", "PiMu", "PiTau", "KEl", "KMu", "KTau", "DEl", "DMu","DTau", "DsEl", "DsMu", "DsTau", "PiPiEl", "PiPiMu", "PiPiTau", "DsstarEl", "DsstarMu", "DsstarTau"]
		for channel in decay_channels: Gamma += hnl_decay_width.dirac_widths(channel, hnl_mass, U2, only_charged_current=only_charged_current)
		return Gamma if dirac_HNL else 2.*Gamma
	@staticmethod
	def semileptonic_width_approx(m_N, U2, dirac_HNL=False, only_charged_current = False):
		Gamma = 0.
		# m_us = [c.m_u, c.m_c] #PDG quark masses 
		# m_ds = [c.m_d, c.m_s, c.m_b]
		# V_us = [[c.V_ud, c.V_us, c.V_ub], [c.V_cd, c.V_cs, c.V_cb]]
		if not only_charged_current:
			Gamma += hnl_decay_width.N_NuaQQ(m_N,c.m_u,True, 1)
			Gamma += hnl_decay_width.N_NuaQQ(m_N,c.m_d,False, 1)
			Gamma *= sum(U2)
		for m_l, U2_l in zip([c.m_el, c.m_mu, c.m_tau], U2):
			if not U2_l: continue
			Gamma +=  hnl_decay_width.N_LaUD(m_N, m_l, c.m_u, c.m_d, U2_l, c.V_ud**2)
		alpha_pi = f.alpha_s(m_N)/np.pi 
		Gamma *= 1.01907 * (1 - 0.006 + alpha_pi + 1.63982 * alpha_pi**2 + 6.37101*alpha_pi**3 + 49.07570*alpha_pi**4 ) # QCD loop corrections and EW radiative correction
		return Gamma if dirac_HNL else 2.*Gamma
	@staticmethod
	def quarks_width(m_N, U2, dirac_HNL=False, QCD_enhancement_factor = False):
		Gamma = 0.
		m_us = [c.m_u, c.m_c] #PDG quark masses 
		m_ds = [c.m_d, c.m_s, c.m_b]
		Vs = [[c.V_ud, c.V_us, c.V_ub], [c.V_cd, c.V_cs, c.V_cb]]
		
		for m_u in m_us: Gamma += hnl_decay_width.N_NuaQQ(m_N,m_u,True, 1)
		for m_d in m_ds: Gamma += hnl_decay_width.N_NuaQQ(m_N,m_d,False, 1)
		Gamma *= sum(U2)
		for m_l, U2_l in zip([c.m_el, c.m_mu, c.m_tau], U2):
			if not U2_l: continue
			for m_u,V_uds in zip(m_us, Vs):
				for m_d, V_ud in zip(m_ds, V_uds): Gamma +=  hnl_decay_width.N_LaUD(m_N, m_l, m_u, m_d, U2_l, V_ud**2)
		alpha_pi = f.alpha_s(m_N)/np.pi 
		if QCD_enhancement_factor: Gamma *= 1.01907 * (1 - 0.006 + alpha_pi + 1.63982 * alpha_pi**2 + 6.37101*alpha_pi**3 + 49.07570*alpha_pi**4 ) # QCD loop corrections and EW radiative correction
		return Gamma if dirac_HNL else 2.*Gamma
	@staticmethod
	def N_LaUD(m_N, m_la, m_u, m_d , U2_la, V2_ud):
		if m_la + m_u + m_d > m_N: return 0
		xf = np.array([m_la/m_N, m_d/m_N, m_u/m_N])
		mask = xf>1e-2
		xf = np.where(mask,xf,0)
		massless_factor = 3 * V2_ud * U2_la * m_N**5 * np.divide(c.G_F**2, 192*np.pi**3) 
		if np.sum(mask) == 0: return massless_factor
		if np.sum(mask) == 1: 
			x2 = np.sum(xf)**2
			return massless_factor * (1 - 8*x2 + 8*x2**3 - x2**4 - 12*x2**2*np.log(x2))
		x2 = np.linspace((xf[0]+xf[1])**2,(1-xf[2])**2,1002)[1:-1]
		return massless_factor * 12 * np.trapz((x2-xf[0]**2 - xf[1]**2)*(1-x2-xf[2]**2)*np.sqrt(f.lambda_Kallen(np.sqrt(x2),xf[0],xf[1])*f.lambda_Kallen(1,np.sqrt(x2),xf[2])) / x2,x2)
	@staticmethod
	def N_NuaQQ(m_N, m_q, up_like, U2_la):
		if 2 * m_q > m_N: return 0
		x2 = (m_q / m_N)**2
		x4 = x2**2
		L =  np.log(np.where(x2>1e-4, np.divide(1. - 3.*x2 - (1.-x2)*np.sqrt(1. - 4.*x2), x2*(1.+np.sqrt(1.-4.*x2))), x4)) #Taylor expansion to prevent log errors
		C1 = 0.25 * (1 - 8./3 * c.sin_w2 + 32./9*c.sin_w2**2) if up_like  else 0.25 * (1 - 4./3* c.sin_w2 + 8/9*c.sin_w2**2)
		C2 = 0.333333*c.sin_w2*(4./3*c.sin_w2 - 1) if up_like else 0.166667*c.sin_w2*(2./3*c.sin_w2 - 1)
		return 3* U2_la * m_N**5 * np.divide(c.G_F**2, 192*np.pi**3) *  (
			C1 *( (1-14*x2 - 2*x4-12*x2**3)*np.sqrt(1-4*x2) + 12*x4 * (x4-1)*L) + 
			4 * C2 * ( x2*(2+10*x2 - 12*x4)*np.sqrt(1-4*x2) + 6*x4*(1-2*x2+2*x4)*L)
		)

class Total_width:
	def __init__(self) -> None:
		path_ = path.dirname(path.abspath(__file__))+"/../../widths/hnl/total_Dirac_widths/"
		#leptonic dirac widths
		self.el_L_singleMeson_interp  = interp1d(*np.loadtxt(path_+"ElMixing_leptonic.dat", delimiter=' ').T)
		self.mu_L_singleMeson_interp  = interp1d(*np.loadtxt(path_+"MuMixing_leptonic.dat", delimiter=' ').T)
		self.tau_L_singleMeson_interp = interp1d(*np.loadtxt(path_+"TauMixing_leptonic.dat",delimiter=' ').T)
		#semileptonic dirac widths
		self.el_SL_singleMeson_interp  = interp1d(*np.loadtxt(path_+"ElMixing_semileptonic_singleMeson.dat", delimiter=' ').T)
		self.mu_SL_singleMeson_interp  = interp1d(*np.loadtxt(path_+"MuMixing_semileptonic_singleMeson.dat", delimiter=' ').T)
		self.tau_SL_singleMeson_interp = interp1d(*np.loadtxt(path_+"TauMixing_semileptonic_singleMeson.dat",delimiter=' ').T)
		self.el_SL_approx_interp = interp1d(*np.loadtxt(path_+"ElMixing_semileptonic_approx.dat", delimiter=' ').T)
		self.mu_SL_approx_interp = interp1d(*np.loadtxt(path_+"MuMixing_semileptonic_approx.dat", delimiter=' ').T)
		self.tau_SL_approx_interp = interp1d(*np.loadtxt(path_+"TauMixing_semileptonic_approx.dat", delimiter=' ').T)
		return
	def width_N(self, m_N, U2, use_dirac_width = False, use_analytic = False, use_only_single_meson=False, only_charged_current=False):
		"""alias for method hnl_decay_width.total_width with the optino of loading pregenerated totla width

		Args:
			m_N (float_like): hnl mass in GeV
			U2 (length3-float-iterable): Vector of [U2_el, U2_mu, U2_tau]
			use_dirac_width (bool, optional): use lepton number conserving decay modes of the hnl. Defaults to False.
			use_analytic (bool, optional): calculate toatal width rather than using pregenerated interpoalation. Defaults to False.
			use_only_single_meson (bool, optional): disregards contributions with more than one hadron in the final state. Defaults to False.
			only_charged_current (bool, optional): for compatibility with other models: disregard HNL coupling to Z boson. Defaults to False.

		Returns:
			float_like: total hnl decay width in Ge
		"""
		use_analytic = use_analytic or only_charged_current # only charged currents are not available precompiled
		leptonic_dirac_width = hnl_decay_width.leptonic_width(m_N, U2, dirac_HNL=True, only_charged_current=only_charged_current) if use_analytic else U2[0]*self.el_L_singleMeson_interp(m_N) + U2[1]*self.mu_L_singleMeson_interp(m_N) + U2[2]*self.tau_L_singleMeson_interp(m_N)
		semileptonic_dirac_width = hnl_decay_width.semileptonic_width_single_meson(m_N, U2, dirac_HNL=True,only_charged_current=only_charged_current) if use_analytic else U2[0]*self.el_SL_singleMeson_interp(m_N) + U2[1]*self.mu_SL_singleMeson_interp(m_N) + U2[2]*self.tau_SL_singleMeson_interp(m_N)
		if not use_only_single_meson and m_N>1: semileptonic_dirac_width = max(semileptonic_dirac_width,
								  hnl_decay_width.semileptonic_width_approx(m_N, U2, dirac_HNL=True, only_charged_current=only_charged_current) if use_analytic else U2[0]*self.el_SL_approx_interp(m_N) + U2[1]*self.mu_SL_approx_interp(m_N) + U2[2]*self.tau_SL_approx_interp(m_N)
								  )
		
		if use_dirac_width: return (leptonic_dirac_width+semileptonic_dirac_width)
		return 2. * (leptonic_dirac_width+semileptonic_dirac_width)

		
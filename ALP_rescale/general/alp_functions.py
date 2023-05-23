#universal functions definitions

import numpy as np
from . import alp_constants as c

def alpha_s(m_a): #leading order alpha_s running
	nf = 6
	beta_0 = 11. - 2/3.*nf

	Lambda_QCD = 0.34

	if m_a <= 1:
		return 1
	elif m_a > 1 and m_a < 1.5: #cubic fit
		return 5.2377 * m_a**3 - 19.913 * m_a**2 + 24.1129 * m_a - 8.4376

	else:
		return 4*np.pi/(beta_0*np.log(m_a**2 /Lambda_QCD**2 ))

	print("Error: alpha_s function - something went wrong")
	return 0

def alpha_s_part(mu): # running weak isospin
	Lambda = c.m_Z # at Z mass
	alpha_s_ref = 0.1192
	beta_s = 7
	return alpha_s_ref/(1-beta_s/(2*np.pi)*alpha_s_ref*np.log(Lambda/mu))

def alpha_t(mu): #running top
	yt = 1 #at the weak scale
	return yt**2/(4*np.pi) 

def alpha_ww(mu): # running weak isospin
	Lambda = c.m_Z # at Z mass
	alpha_EM_Z = 127.952 # at Z mass
	alpha_ww_ref = alpha_EM_Z/np.sin(c.theta_w)**2
	beta_ww = 19/6
	return alpha_ww_ref/(1-beta_ww/(2*np.pi)*alpha_ww_ref*np.log(Lambda/mu))

def alpha_bb(mu): # running weak hypercharge
	Lambda = c.m_Z # at Z mass
	alpha_EM_Z = 127.952 # at Z mass
	alpha_bb_ref = alpha_EM_Z/np.cos(c.theta_w)**2
	beta_bb = -41/6
	return alpha_bb_ref/(1-beta_bb/(2*np.pi)*alpha_bb_ref*np.log(Lambda/mu))

def tau(Gamma): #lifetime in ps
	return c.hPl*1E12/Gamma
def Gamma(tau): #width in GeV, tau in ps
	return c.hPl*1E12/tau

def lambda_Kallen(xx,yy,zz): #Kallen function
    return (xx**2-(yy+zz)**2)*(xx**2-(yy-zz)**2)
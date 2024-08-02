#universal functions definitions

import numpy as np
from ALP_rescale.general import constants as c

def alpha_s_perturbative(q, order = 0, nf = -1, Lambda_QCD = 0.34):
	beta_0 = 11. - 2/3.*nf
	log_mL = np.log((q/Lambda_QCD)**2)
	alpha_LO = 4*np.pi/(beta_0*log_mL)
	if order == 0: 	return  alpha_LO
	log_log_mL = np.log(log_mL)
	beta_1 = 102. - 38/3*nf
	if order == 1:	return alpha_LO * (1 - beta_1 / beta_0**2 * log_log_mL/log_mL)
	beta_2 = 2857/2 - 5033/18*nf + 325/54*nf**2
	if order == 2: return alpha_s_perturbative(q,1, nf=nf, Lambda_QCD=Lambda_QCD) + alpha_LO*(
		beta_1**2 / (beta_0**4*log_mL**2) * (log_log_mL**2- log_log_mL - 1 + beta_2 * beta_0 / beta_1**2 )
		)
	xi = 1.2026
	beta_3 =(149753/6 + 3564*xi) - (1078361/162 + 6508/27*xi)*nf + (50065/162 + 6472/81*xi)*nf**2 + 1093/729*nf**3
	beta_2
	return alpha_s_perturbative(q, 2, nf=nf, Lambda_QCD=Lambda_QCD) + alpha_LO * ( 
		+ beta_1**3 / (beta_0**6*log_mL**3) * (-log_log_mL**3 + 5/2*log_log_mL**2 +2*log_log_mL - 1/2 - 3*beta_2*beta_0/beta_1**2*log_log_mL + beta_3*beta_0**2/(2*beta_1**3))
		)

def adaptive_nf_LambdaQCD(q):
	nf, Lambda_QCD = [2,0.31] # ensuring Lambda_{QCD}^{nf=3}=0.34GeV
	for m_q in [c.m_s,c.m_c,c.m_b,c.m_t]: 
		if m_q < q: 
			nf+=1
			Lambda_QCD = Lambda_QCD*(Lambda_QCD/m_q)**( 2/(33 - 2*nf))
	return nf, Lambda_QCD

def alpha_s(q, order = 3, q_suspension = [0, 1], nf = -1, Lambda_QCD = 0.34):
	"""Strong coupling constant alpha_s at energy scale q.

	Args:
		q (float): energy scale
		order (str, optional): Perturbative alpha order in beta expansion. Defaults to 3.
		q_suspension (list of floats, optional): limit points of interpolation region between non-perturbative and pertrubative regime. Defaults to [0, 1].
		nf (int, optional): number of flavours. Defaults to -1 (adapting uisng using adaptive_nf_LambdaQCD(q)).
		Lambda_QCD (float, optional): QCD Landau Pole scale. Defaults to 0.34. When nf is -1 pole adapts using adaptive_nf_LambdaQCD(q).

	Returns:
		float: strong coupling constant with cubic interpolation into the non-pertrubative regime.
	"""
	if q < q_suspension[0]: return 1
	if q >= q_suspension[1]: 
		if nf == -1: nf, Lambda_QCD = adaptive_nf_LambdaQCD(q)
		return alpha_s_perturbative(q, order=order, nf=nf, Lambda_QCD=Lambda_QCD)
	# fitting suspended part to 3rd order polynomial
	if nf == -1: nf, Lambda_QCD = adaptive_nf_LambdaQCD(q_suspension[1])
	x1 = q_suspension[0]
	y1 = 1
	m1 = 0
	x2 = q_suspension[1]
	y2 = alpha_s_perturbative(x2, order=order, nf=nf, Lambda_QCD=Lambda_QCD)
	m2 = (alpha_s_perturbative(1.000001*x2, order=order, nf=nf, Lambda_QCD=Lambda_QCD) - y2) * 1000000 / x2

	a_fit = (m1+m2-2*(y2-y1)/(x2-x1))/(x1-x2)**2
	b_fit = (m2-m1)/(2*(x2-x1))-(3/2)*(x1+x2)*a_fit
	c_fit = m1-3*(x1**2)*a_fit -2*x1*b_fit
	d_fit = y1-(x1**3)*a_fit -(x1**2)*b_fit-x1*c_fit
	return a_fit * q**3 + b_fit * q**2 + c_fit * q + d_fit

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
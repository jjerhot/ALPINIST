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

def F_pi_analytic(s):
	"""Resonance dominance model with parameters taken from [1205.2228]
	Args:
		s (float): invariant mass squared of pion pair 
	Returns:
		F_pi (float): structure function
	"""
	Gamma_omega =  8.13e-3 
	m_omega = 781.91e-3
	BW_omega = m_omega**2 / (m_omega**2 - s - 1.j * m_omega*Gamma_omega)

	beta_pi = lambda s: np.sqrt(np.clip(1-4*c.m_pi**2/s,0,None))
	k = lambda s : np.sqrt(np.clip(s/4 - c.m_pi**2,0,None))#0.5*np.sqrt(s)*beta_pi(s)
	h = lambda s : 2./np.pi * k(s)/np.sqrt(s) * np.log((np.sqrt(s)+2.*k(s))/(2*c.m_pi))
	hprime = lambda s: np.divide( s*beta_pi(s) + 16*c.m_pi**2 *np.log( (k(s) + 2*np.sqrt(s)) / (4*c.m_pi)),  8*np.pi*s*k(s)) 
	d = lambda m : 3./np.pi * c.m_pi**2/k(m**2)**2 * np.log((m+2*k(m**2))/(2.*c.m_pi)) + m/np.pi*(0.5/k(m**2) - c.m_pi**2/ k(m**2)**3)
	f = lambda s, m, Gamma : Gamma*m**2/k(m**2)**3 * ( k(s)**2 * (h(s)- h(m**2)) + (m**2 - s) * k(m**2) * hprime(m**2) ) 
	G = lambda s, m, Gamma : Gamma * s / m**2 * ( beta_pi(s)/beta_pi(m**2) )**3

	BW_rho = lambda s, m_rho, Gamma_rho : np.divide(	m_rho**2 * (1 + d(m_rho) * Gamma_rho/m_rho),
														m_rho**2 - s + f(s, m_rho, Gamma_rho) - 1.j * m_rho * G(s, m_rho, Gamma_rho),
											dtype=np.csingle)
	
	c_omega = 1.664e-3 * np.exp(-0.011j)
	c_rhoi	= 0.158	   * np.exp(3.76j)
	c_rhoii = 0.068    * np.exp(1.39j)
	c_rhoiii= 0.0051   * np.exp(0.70j)

	return np.divide(BW_rho(s, 775.02e-3, 149.59e-3 ) * (1 + c_omega * BW_omega)/(1+c_omega)+ c_rhoi * BW_rho(s, 1493e-3, 427e-3) +  c_rhoii * BW_rho(s, 1861e-3, 316e-3) + c_rhoiii * BW_rho(s, 2254e-3, 109e-3),
					1 + c_rhoi + c_rhoii + c_rhoiii)

def alpha_s_part(mu): # running weak isospin
	Lambda = c.m_Z # at Z mass
	alpha_s_ref = 0.1192
	beta_s = 7
	return alpha_s_ref/(1-beta_s/(2*np.pi)*alpha_s_ref*np.log(Lambda/mu))

def alpha_t(mu, order = 3, q_suspension = [0, 1], nf = -1, Lambda_QCD = 0.34):
	"""Top yukawa coupling constant alpha_t at energy scale q. Assuming that y_t(q = m_t) = 1
		Interwoven with strong coupling alpha_s.

	Args:
		q (float): energy scale
		order (str, optional): Perturbative alpha strong order in beta expansion. Defaults to 3.
		q_suspension (list of floats, optional): limit points of interpolation region between non-perturbative and pertrubative regime. Defaults to [0, 1].
		nf (int, optional): number of flavours. Defaults to -1 (adapting using using adaptive_nf_LambdaQCD(q)).
		Lambda_QCD (float, optional): QCD Landau Pole scale. Defaults to 0.34. When nf is -1 pole adapts using adaptive_nf_LambdaQCD(q).

	Returns:
		float: strong coupling constant with cubic interpolation into the non-pertrubative regime.
	"""
	alpha_s_ref = alpha_s(c.m_q[5],order,q_suspension,nf,Lambda_QCD)
	alpha_s_mu = alpha_s(mu,order,q_suspension,nf,Lambda_QCD)
	return 1./(4*np.pi) * (alpha_s_mu/alpha_s_ref)**(8/7)/(1 + 9/2*(4*np.pi*alpha_s_ref)*((alpha_s_mu/alpha_s_ref)**(1/7) - 1))


def alpha_ww(mu): # running weak isospin
	Lambda = c.m_Z # at Z mass
	alpha_EM_Z = 1./127.952 # at Z mass
	alpha_ww_ref = alpha_EM_Z/np.sin(c.theta_w)**2
	beta_ww = 19/6
	return alpha_ww_ref/(1-beta_ww/(2*np.pi)*alpha_ww_ref*np.log(Lambda/mu))

def alpha_bb(mu): # running weak hypercharge
	Lambda = c.m_Z # at Z mass
	alpha_EM_Z = 1./127.952 # at Z mass
	alpha_bb_ref = alpha_EM_Z/np.cos(c.theta_w)**2
	beta_bb = -41/6
	return alpha_bb_ref/(1-beta_bb/(2*np.pi)*alpha_bb_ref*np.log(Lambda/mu))

def alpha_EM(mu): # running EM constant
	Lambda = c.m_Z # at Z mass
	alpha_EM_Z = 1/127.952 # at Z mass
	return alpha_EM_Z/(1-beta_QED(mu)/(2*np.pi)*alpha_EM_Z*np.log(Lambda/mu))

def beta_QCD(mu):
	nq = 0
	for mq in c.m_q:
		if mu > mq: nq += 1.
	return 11.-2./3*nq

def beta_QED(mu):
	beta = 0
	for m in c.m_qt:
		if mu > m: beta += 3 * (2./3)**2
	for m in c.m_qb:
		if mu > m: beta += 3 * (-1./3)**2
	for m in c.m_l:
		if mu > m: beta += 1
	return -4./3 * beta

def tau(Gamma): #lifetime in ps
	return c.hPl*1E12/Gamma

def Gamma(tau): #width in GeV, tau in ps
	return c.hPl*1E12/tau

def lambda_Kallen(xx,yy,zz): #Kallen function
    return (xx**2-(yy+zz)**2)*(xx**2-(yy-zz)**2)
#universal functions definitions

import os
import numpy as np
from . import exotic_constants as c
from . import exotic_production_setup as s
import random
from scipy import integrate
from scipy.optimize import leastsq
from scipy.stats import rv_continuous

rng = np.random.default_rng() # default random generator used for random processes 

def number_to_3sigfigs_str(number):
    assert 0 <= number < 1e12, f"number {number} outside of interpretation range (0, 999e9)"
    translator = ['','k','M','B']
    number_split = f"{number:3.3e}".split("e")
    order_1000 = int(int(number_split[1])/3)
    raise_by = int(number_split[1])%3
    sigfigs = str(float(f"{(float(number_split[0])):.2e}")*10**raise_by)
    while '.' in sigfigs and sigfigs[-1] in ['0','.']: sigfigs = sigfigs[:-1]
    return sigfigs+translator[order_1000]

def export(experiment,fileName,data,header=""):
	output_dir = os.getcwd()+'/../tab_prod/'
	outPath = output_dir + experiment + '/'
	if not os.path.isdir(outPath): os.makedirs(outPath)
	np.savetxt(outPath + fileName,data,fmt='%.3e',header=header)
	print('\n[Info:] \t', 'File ' + fileName + ' saved to ' + outPath)
	return

def lambda_Kallen(xx,yy,zz):
	"""Kallen function for sqrt entries xx, yy, zz """
	return (xx**2-(yy+zz)**2)*(xx**2-(yy-zz)**2)

def lorentz_transf(gam,bet):
	"""returns lorentz transformation matrix for given gamma and beta"""
	return [[gam,			gam*bet[0],										gam*bet[1],										gam*bet[2]],
			[gam*bet[0], 	1+(gam-1)*bet[0]**2/np.linalg.norm(bet)**2, 	(gam-1)*bet[0]*bet[1]/np.linalg.norm(bet)**2,	(gam-1)*bet[0]*bet[2]/np.linalg.norm(bet)**2],
			[gam*bet[1],	(gam-1)*bet[1]*bet[0]/np.linalg.norm(bet)**2,	1+(gam-1)*bet[1]**2/np.linalg.norm(bet)**2,		(gam-1)*bet[1]*bet[2]/np.linalg.norm(bet)**2],
			[gam*bet[2],	(gam-1)*bet[2]*bet[0]/np.linalg.norm(bet)**2,	(gam-1)*bet[2]*bet[1]/np.linalg.norm(bet)**2,	1+(gam-1)*bet[2]**2/np.linalg.norm(bet)**2]]

def lorentz_transf_vectorized(gam,b0,b1,b2,bsq):
	"""optimized implementation of lorentz_transf"""
	return np.array(
			[[gam,      gam*b0,            		gam*b1,             	gam*b2],
			[gam*b0,    1+(gam-1)*b0**2/bsq, 	(gam-1)*b0*b1/bsq,		(gam-1)*b0*b2/bsq],
			[gam*b1,    (gam-1)*b0*b1/bsq,   	1+(gam-1)*b1**2/bsq,	(gam-1)*b1*b2/bsq],
			[gam*b2,    (gam-1)*b0*b2/bsq,   	(gam-1)*b1*b2/bsq,		1+(gam-1)*b2**2/bsq]]).T

def beta_lor(p,E):
	"""Lotrentz beta given momentum p and Energy E"""
	return p/E

def gamma_lor(beta):
	"""Lotrentz gamma given Lorentz beta"""
	return 1/np.sqrt(1-beta**2)
        
def f_kde(kde,norm,e,th):
	"returns probability evaluation of kde for given energy and angle with normalisation norm"
	val = [0]
	if e>0 and th>0:
		val = norm * np.exp(kde.score_samples([[np.log(e),np.log(th)]]))/(e*th)
	return val

def decay_meson(p_cm,m_exo,lor):
	"""returns energy and angle of exotic in random decay m_1 -> m_2 + exotic, given mother lorentz boost lor."""
	th_exo_cm = np.arccos(random.uniform(-1, 1))
	phi_exo_cm = random.uniform(0, 2*np.pi)

	p_lor_cm = [np.sqrt(p_cm**2+m_exo**2),p_cm*np.sin(th_exo_cm)*np.cos(phi_exo_cm),p_cm*np.sin(th_exo_cm)*np.sin(phi_exo_cm),p_cm*np.cos(th_exo_cm)]
	p_lor_exo = np.dot(lor,p_lor_cm)
	#return e_exo and th_exo
	return p_lor_exo[0], np.arccos(p_lor_exo[3]/np.linalg.norm(p_lor_exo[1:]))

def decay_meson_vectorized(p_cm, m_exo,lors):
	"""returns energy and angle of exotic in random decay m_1 -> m_2 + exotic, given mother lorentz boosts lors."""
	if np.ndim(p_cm) == 0: # creating N random mother-restframe momenta
		N = np.shape(lors)[0]
		th_cm = np.arccos(2*np.random.rand(N)-1)
		phi_cm = 2*np.pi*np.random.rand(N)
		p_cm = np.transpose([np.sqrt(p_cm**2+m_exo**2)*np.ones(N),p_cm*np.sin(th_cm)*np.cos(phi_cm),p_cm*np.sin(th_cm)*np.sin(phi_cm),p_cm*np.cos(th_cm)])
		del th_cm, phi_cm
	p_lor = np.einsum('ijk,ik->ij', lors, p_cm) # boosting according to mother boost
	#return e_exo and th_exo
	return np.array((p_lor[:,0], np.arccos(p_lor[:,3]/np.linalg.norm(p_lor[:,1:],axis =1))))

def decay_meson_cartesian_momentum(p_cm, m_exo, lor):
	"""returns carthesian momenta of exotic in random decay m_1 -> m_2 + exotic, given lorentz transform from mother rest to lab frame lor."""
	th_exo_cm = np.arccos(random.uniform(-1, 1))
	phi_exo_cm = random.uniform(0, 2*np.pi)

	p_lor_cm = [np.sqrt(p_cm**2+m_exo**2),p_cm*np.sin(th_exo_cm)*np.cos(phi_exo_cm),p_cm*np.sin(th_exo_cm)*np.sin(phi_exo_cm),p_cm*np.cos(th_exo_cm)]
	p_lor_exo = np.dot(lor,p_lor_cm)

	#return e_exo and th_exo
	return p_lor_exo[1:]

def decay_meson_cartesian_momentum_vectorized(p_cm, m_daughter, lors): #returns energy and angle of exotic in random decay m_1 -> m_2 + exotic
	"""returns carthesian momenta of exotic in random decay m_1 -> m_2 + exotic, given lorentz transforms from mother rest to lab frame lors."""
	if np.ndim(p_cm) == 0: # creating N random mother-restframe momenta
		N = np.shape(lors)[0]
		th_exo_cm = np.arccos(2*np.random.rand(N)-1)
		phi_exo_cm = 2*np.pi*np.random.rand(N)
		#generating random phi and theta in cm frame
		p_lor_daughter_cm = np.array([np.sqrt(p_cm**2+m_daughter**2)*np.ones(N),p_cm*np.sin(th_exo_cm)*np.cos(phi_exo_cm),p_cm*np.sin(th_exo_cm)*np.sin(phi_exo_cm),p_cm*np.cos(th_exo_cm)])
		#boosting according to mother boost
		p_lor_daughter = np.einsum('ijk,ki->ij', lors, p_lor_daughter_cm)
		del p_lor_daughter_cm
	else: p_lor_daughter = np.einsum('ijk,ki->ij', lors, p_cm)
	#return carthesian 3 momentum 
	return p_lor_daughter[:,1:]

def lorentz_boosts_vectorized(event_list, mass, ndecays=1):
	"""returns lorentz transform matrices given event list and mass repeated ndecays amount of times."""
	ps = np.array(event_list)[:,1:]
	ps_squared = np.linalg.norm(ps, axis=1)**2
	es   	   = np.sqrt(mass**2 + ps_squared)
	lors = lorentz_transf_vectorized(es/mass, ps[:,0]/es, ps[:,1]/es, ps[:,2]/es, ps_squared / es**2 )
	return np.repeat(lors, ndecays, axis=0)

def rescale_p_2gamma(m_exo, m_meson):
	return 2*(1-m_exo**2/m_meson**2)**3

def rescale_v_p_gamma(m_exo, m_vector, m_scalar):
	return (m_vector**2-m_exo**2-m_scalar**2)**2*np.sqrt((m_vector**2-m_exo**2+m_scalar**2)**2 - 4*m_vector**2 * m_scalar**2)/(m_vector**2-m_exo**2)**3

class Input_reweight:
	""""Utility class to 
			reshape an external meson momentum distribution to match a reference	->  reshape_input
			draw a meson momentum distribution list based on external reference 	->  generate_plab_distribution"""
	########## external utility functions ##########
	@staticmethod 
	def reshape_input(meson_list, ndecays, normalize_to = "LEBC", meson = "D0", p_beam = 400):
		weights = Input_reweight.get_reshape_weights(meson_list[:,1:], normalize_to, meson, p_beam) * ndecays
		ceils = np.ceil(weights)
		draw = ceils * np.random.rand(weights.size)
		repeats = np.where(weights<draw, weights, ceils).astype(int)
		reshaped_list = np.empty((np.sum(repeats),4),dtype=np.float64)
		curr = 0
		for meson_event, repeat in zip(meson_list, repeats):
			reshaped_list[curr:(curr+repeat),:] = np.repeat([meson_event], repeat, axis=0)
			curr += repeat
		return reshaped_list
	@staticmethod
	def generate_momentum_array(meson, N_Interactions, random_generator = np.random.default_rng(), p_beam = 400, reference = "Alpinist"):
		"""Generate random meson momenta array using empirical distributions for N_Interactions ccbar events according to reference at beam momentum p_beam in GeV. 

		Args:
			meson (string): desired meson options ["D0","D0bar","D","Dbar","Ds","Dsbar"].
			N_Interactions (_type_): number of ccbar events.
			random_generator (_type_, optional): random generator to use complying with properties of np.random.default_rng(), Defaults to np.random.default_rng().
			p_beam (int, optional): beam momentum in GeV. Defaults to 400.
			reference (str, optional): literature reference. Defaults to "Alpinist".

		Returns:
			momentum array (np.array(np.float64)): array with rows [particle ID, px, py, pz] (momenta in GeV)
		"""
		N_mesons = int(N_Interactions*Input_reweight.cc_multiplicities(reference, meson, p_beam))
		b = Input_reweight.get_b(reference, meson, p_beam) 
		n = Input_reweight.get_n(reference, meson, p_beam)
		rho_ = Input_reweight.N_per_NxFpos[reference]
		xF = np.concatenate( (-random_generator.beta(1., (n+1)/(rho_-1), size = int(N_mesons*(1-1/rho_)) ), random_generator.beta(1., n+1, size = int(N_mesons/rho_) )) )
		pT = random_generator.exponential(1/b, size=xF.shape[0])
		if Input_reweight.square_pT[reference]: pT = np.sqrt(pT)
		PID =  [Input_reweight.meson_to_PID[meson]*np.ones_like(xF, dtype=np.float64)]
		return np.concatenate((PID, Input_reweight.xFpT_to_plab(xF, pT, Input_reweight.m_meson[meson], rng, p_beam = p_beam))).T
	
	########## matching functions and dicts ##########
	@staticmethod
	def get_b(reference, meson, p_beam):
		if reference == "CHARM": return 2. # Phys. Lett. 166B (1986), 473-478
		if reference == "BEBC":  return 2. # Phys. Lett. 160B (1985) 207
		if reference == "E653":  return 0.81 # + 0.10 - 0.08, Phys. Lett. B 263, 573 (1991)
		if reference == "LEBC": #LEBC distribution Nucl.Phys.B 199 (1982) 424-440
			if meson == "D0": 	return 1.04 # +- 0.19
			if meson == "D0bar":return 1.49 # +- 0.32
			if meson == "D":	return 0.75 # +- 0.14
			if meson == "Dbar": return 0.93 # +- 0.18
			print(f"[Warning:] \t No pT fit data from LEBC for requested {meson} meson, defaulting to D/Dbar distributions")	
			if "bar" in meson: 	return 0.93 # +- 0.18
			return 0.75 # +- 0.14
		if reference == "LEBC_l": #LEBC distribution Nucl.Phys.B 199 (1982) 424-440
			if meson == "D0": 	return 1.04 - 0.19
			if meson == "D0bar":return 1.49 - 0.32
			if meson == "D":	return 0.75 - 0.14
			if meson == "Dbar": return 0.93 - 0.18
			print(f"[Warning:] \t No pT fit data from LEBC for requested {meson} meson, defaulting to D/Dbar distributions")	
			if "bar" in meson: 	return 0.93 - 0.18
			return 0.75 - 0.14
		if reference == "LEBC_u": #LEBC distribution Nucl.Phys.B 199 (1982) 424-440
			if meson == "D0": 	return 1.04 + 0.19
			if meson == "D0bar":return 1.49 + 0.32
			if meson == "D":	return 0.75 + 0.14
			if meson == "Dbar": return 0.93 + 0.18
			print(f"[Warning:] \t No pT fit data from LEBC for requested {meson} meson, defaulting to D/Dbar distributions")	
			if "bar" in meson: 	return 0.93 + 0.18
			return 0.75 + 0.14
		# else assuming Alpinist distribution fit results
		sqrts_CM = np.sqrt(c.m_p**2 + c.m_p*(c.m_p+2*np.sqrt(c.m_p**2 + p_beam**2)))
		ref_data = np.loadtxt('./tab_mesons/charm/Alpinist_dxcc_fit.dat',usecols=(0,3,4)).T
		if reference == "Alpinist_l": return np.interp(sqrts_CM, ref_data[0], ref_data[1]-ref_data[2])
		if reference == "Alpinist_u": return np.interp(sqrts_CM, ref_data[0], ref_data[1]+ref_data[2]) 
		return np.interp(sqrts_CM, ref_data[0], ref_data[1])


	@staticmethod
	def get_n(reference, meson, p_beam):
		if reference == "CHARM": return 5. # +0-1. (variation between publications) Phys. Lett. 166B (1986), 473-478
		if reference == "BEBC":  return 4. # Physics Letters B 160, 207 (1985)
		if reference == "E653":  return 6.8 # + 2.1 - 1.9,  Phys. Lett. B 263, 573 (1991)
		if reference == "LEBC": # LEBC distribution Nucl.Phys.B 199 (1982) 424-440
			if meson == "D0": 	return 5.4 # +- 1.1
			if meson == "D0bar":return 8.1 # +- 1.9
			if meson == "D":	return 3.1 # +- 0.8
			if meson == "Dbar": return 5.4 # +- 1.2
			print(f"[Warning:] \t No xF fit data from LEBC for requested {meson} meson, defaulting to D/Dbar distributions")
			if "bar" in meson: 5.4 # +- 1.2
			return 3.1 # +- 0.8
		if reference == "LEBC_l": # LEBC distribution Nucl.Phys.B 199 (1982) 424-440
			if meson == "D0": 	return 5.4 + 1.1
			if meson == "D0bar":return 8.1 + 1.9
			if meson == "D":	return 3.1 + 0.8
			if meson == "Dbar": return 5.4 + 1.2
			print(f"[Warning:] \t No xF fit data from LEBC for requested {meson} meson, defaulting to D/Dbar distributions")
			if "bar" in meson: 5.4 + 1.2
			return 3.1 + 0.8
		if reference == "LEBC_u": # LEBC distribution Nucl.Phys.B 199 (1982) 424-440
			if meson == "D0": 	return 5.4 - 1.1
			if meson == "D0bar":return 8.1 - 1.9
			if meson == "D":	return 3.1 - 0.8
			if meson == "Dbar": return 5.4 - 1.2
			print(f"[Warning:] \t No xF fit data from LEBC for requested {meson} meson, defaulting to D/Dbar distributions")
			if "bar" in meson: 5.4 - 1.2
			return 3.1 - 0.8
		# else assuming Alpinist distribution fit results
		sqrts_CM = np.sqrt(c.m_p**2 + c.m_p*(c.m_p+2*np.sqrt(c.m_p**2 + p_beam**2)))
		ref_data = np.loadtxt('./tab_mesons/charm/Alpinist_dxcc_fit.dat',usecols=(0,1,2)).T
		if reference == "Alpinist_l": return np.interp(sqrts_CM, ref_data[0], ref_data[1]+ref_data[2])
		if reference == "Alpinist_u": return np.interp(sqrts_CM, ref_data[0], ref_data[1]-ref_data[2]) 
		return np.interp(sqrts_CM, ref_data[0], ref_data[1])
	@staticmethod
	def get_xF_distribution(xF, reference, meson):
		return np.power(1-np.abs(xF), Input_reweight.get_n(reference, meson))
	@staticmethod
	def get_pT_distribution(pT, reference, meson):
		if Input_reweight.square_pT[reference]: return np.exp(-Input_reweight.get_b(reference, meson)*pT**2)
		return np.exp(-Input_reweight.get_b(reference, meson)*pT)
	m_meson = {	"D0"	: c.m_D0,
				"D0bar" : c.m_D0,
				"D"		: c.m_D,
				"Dbar"	: c.m_D,
				"Ds"	: c.m_Ds,
				"Dsbar"	: c.m_Ds,
	}
	meson_to_PID = {
				"D"		: 411,
				"Dbar"	: -411,
				"D0"	: 421,
				"D0bar" : -421,
				"Ds"	: 431,
				"Dsbar"	: -431,
	}
	square_pT = {
		"E653":		True,
		"E653_l":	True,
		"E653_u":	True,
		"CHARM":	False,
		"BEBC":		False,
		"LEBC":		True,
		"LEBC_l":	True,
		"LEBC_u":	True,
		"Alpinist": True,
		"Alpinist_l": True,
		"Alpinist_u": True
	}
	N_per_NxFpos = { # assymmetry between xF<0 to xF>0, this is 2 for p and 1.6 for pi beams
		"E653":		2,
		"E653_l":	2,
		"E653_u":	2,
		"CHARM":	2,
		"BEBC":		2,
		"LEBC":		2,
		"LEBC_l":	2,
		"LEBC_u":	2,
		"Alpinist": 2,
		"Alpinist_l": 2,
		"Alpinist_u": 2
	}
	@staticmethod
	def cc_multiplicities(reference, meson, p_beam=0) :
		"""meson multiplicities per ccbar event

		Args:
			reference (str): reference 
			meson (str): meson choice
			p_beam (int, optional): beam momentum in GeV. Defaults to 0.

		Returns:
			float: N_meson/N_ccbar
		"""
		if reference in ["LEBC","LEBC_l","LEBC_u"] :
			return {
			"D"		: 0.3904,
			"Dbar"	: 0.4247,
			"D0"	: 0.7192,
			"D0bar"	: 0.5411,
			"Ds"	: 0.1,
			"Dsbar"	: 0.1
			}[meson]
		if reference == "CHARM":
			return { 
			"D"		: 0.333333,
			"Dbar"	: 0.333333,
			"D0"	: 0.666667,
			"D0bar"	: 0.666667,
			"Ds"	: 0, # F mesons considered only in seperate sarch
			"Dsbar"	: 0  # F mesons considered only in seperate sarch
			}[meson],
		if reference == "BEBC":
			return { #Phys. Lett. 160B (1985) 207
			"D"		: 0.333333,
			"Dbar"	: 0.333333,
			"D0"	: 0.666667,
			"D0bar"	: 0.666667,
			"Ds"	: 0, 
			"Dsbar"	: 0,
			}[meson]
		if reference == "E653":
			if "Ds" in meson: return 0.1
			return 0.5
		return {	# assume 
			120:{
			"D"		: 0.22053,
			"Dbar"	: 0.304277,
			"D0"	: 0.40128,
			"D0bar"	: 0.58138,
			"Ds"	: 0.060693,
			"Dsbar"	: 0.083927
			},
			400:{
				"D"		: 0.25717,
				"Dbar"	: 0.299203,
				"D0"	: 0.48242,
				"D0bar"	: 0.571143,
				"Ds"	: 0.07694,
				"Dsbar"	: 0.09039
			},
			800:{

				"D"		: 0.267443,
				"Dbar"	: 0.29886,
				"D0"	: 0.5039,
				"D0bar"	: 0.56981,
				"Ds"	: 0.08118,
				"Dsbar"	: 0.091053
			},
		}[p_beam][meson]

	literature = {"LEBC": "Nucl.Phys.B 199 (1982) 424", "CHARM":"Phys.Lett.B 166 (1986) 473", "BEBC":" Phys. Lett. 160B (1985) 207", "E653":"Phys. Lett. B 263 (1991) 573", "Alpinist":"2407.08673"} #literature entries for refrences, "SMTs": "SHiP-NOTE-2015-009" is not publically available
	########## internal utility functions ##########
	@staticmethod
	def plab_to_xF(p_meson, m_meson, p_beam = 400, beam_particle_mass=c.m_p, target_mass=c.m_p):
		pMom = np.linalg.norm(p_meson, axis = -1)
		e 	 = np.sqrt(pMom**2 + m_meson**2)
		sqrts_CM = np.sqrt(beam_particle_mass**2 + target_mass*(target_mass+2*np.sqrt(beam_particle_mass**2 + p_beam**2)))
		pl_CM_max = np.sqrt(sqrts_CM**2/4 - m_meson**2)
		bg_CM = p_beam/sqrts_CM
		g_CM  = np.sqrt(p_beam*p_beam + sqrts_CM*sqrts_CM)/sqrts_CM
		pl_CM = g_CM*p_meson[:,2] - bg_CM*e
		return pl_CM/pl_CM_max
	@staticmethod
	def plab_to_pT(p_meson):
		return np.sqrt(p_meson[:,0]**2+p_meson[:,1]**2)#np.linalg.norm(p_meson[:,:-1], axis = -1)
	@staticmethod
	def plab_to_pT2(p_meson):
		return p_meson[:,0]**2+p_meson[:,1]**2
	@staticmethod
	def xFpT_to_plab(xF, pT, m_meson, rng = np.random.default_rng(), p_beam = 400, beam_particle_mass=c.m_p, target_mass=c.m_p):
		sqrts_CM = np.sqrt(beam_particle_mass**2 + target_mass*(target_mass+2*np.sqrt(beam_particle_mass**2 + p_beam**2)))
		pl_CM_max = np.sqrt(sqrts_CM**2/4 - m_meson**2)
		bg_CM = p_beam/sqrts_CM
		g_CM  = np.sqrt(p_beam*p_beam + sqrts_CM*sqrts_CM)/sqrts_CM
		pl_CM = xF * pl_CM_max
		delta_g = g_CM**2 - bg_CM**2
		pl  = np.divide(g_CM*pl_CM+bg_CM*np.sqrt(pl_CM**2+delta_g*(m_meson**2 + pT**2)), delta_g)
		alpha_ = np.random.rand(np.size(xF))
		return np.array([rng.choice([-1,1],size = np.size(xF))*np.sqrt(alpha_)*pT, rng.choice([-1,1],size = np.size(xF))*np.sqrt((1-alpha_))*pT, pl])
	@staticmethod
	def fit_xF_dist(xF): 
		""" returns a, b for assumed a*(1-xF)^b to the xF multiplicity for xF in (0, 0.5)"""
		dS_hist, xF_lims = np.histogram(xF, bins = 100, range = (0, 0.5))
		dS_dist = dS_hist / float(np.max(dS_hist))
		xF_bins = (xF_lims[:-1] + xF_lims[1:]) / 2
		min_fun = lambda pars : np.where(dS_dist!=0, np.log(dS_dist) - pars[1]*np.log(pars[0]*(1-np.abs(xF_bins))), 0)#fitting over log to account for orders of magnitude in dS
		fit_params = leastsq(min_fun, x0 = [1,5])[0]
		return fit_params
	@staticmethod
	def get_reshape_weights(p_mesons, normalize_to = "LEBC", meson = "D0", p_beam = 400, xF_validation = np.array([])):
		xF = Input_reweight.plab_to_xF(p_mesons, Input_reweight.m_meson[meson], p_beam)
		a, b = Input_reweight.fit_xF_dist(xF)
		print(a,b)
		if xF_validation.size: 
			xF_intpoints = np.linspace(0, np.max(xF), 1000)	
			normalization = np.trapz(Input_reweight.get_xF_distribution(xF_intpoints,reference=normalize_to,meson=meson) / (a*(1-xF_intpoints)**b), xF_intpoints)
			return Input_reweight.get_xF_distribution(xF_validation,reference=normalize_to,meson=meson) / (a*(1-np.abs(xF_validation))**b) / normalization
		weights = Input_reweight.get_xF_distribution(xF,reference=normalize_to,meson=meson) / (a*(1-np.abs(xF))**b)
		return  weights * (float(weights.size) / np.sum(weights))

class Meson_form_factors:
	"""meson form factor functions for semileptonic decays of pseudoscalar mesons including neutrinos with references as given for individual functions"""
	@staticmethod
	def f_plus(q_sq, H_in, H_out):
		if H_in == "D":
			if H_out == "K":	return Meson_form_factors.f_plus_template(q_sq, 0.224, [0.7877,-0.97,-0.3], c.m_D, c.m_K)#FLAG Review 2021  M_pole = mDs* = 2.1122 GeV, alt: return Meson_form_factors.f_D_H(q_sq,  0.224, 0.7647, 2.084, c.m_D, c.m_pi)# ETM 17D 
			if H_out == "Pi":	return Meson_form_factors.f_D_H(q_sq, 0.1314, 0.6117, 1.985, c.m_D, c.m_pi)# ETM 17D
		elif H_in == "Ds":
			if H_out == "Eta":	return Meson_form_factors.f_Ds_eta(q_sq, 0.224, 0.495, 0.198)# 1805.08567v3
		
		elif H_in == "B":
			if H_out == "D":	return Meson_form_factors.f_plus_template(q_sq, 0., [0.896, -7.94, 51.4], c.m_B, c.m_D) # FLAG Review 2021 M_pole = 0 
			if H_out == "Pi":	return Meson_form_factors.f_plus_template(q_sq, 0.0353, [0.404, -0.68, -0.86], c.m_B, c.m_pi)# FLAG Review 2021  M_pole = mB* = 5.32465GeV

		elif H_in =="Bs":
			if H_out == "Ds":	return Meson_form_factors.f_plus_template(q_sq, 0.02225,[-0.075, 3.24, 0.7], c.m_Bs, c.m_Ds)#FLAG Review 2021 M_pole = Bc0 = 6.704 GeV
			if H_out == "K":	return Meson_form_factors.f_plus_template(q_sq, 0.0353, [0.374, -0.672, 0.07, 1.34], c.m_Bs, c.m_K)#FLAG Review 2021 M_pole =  mB* = 5.32465GeV
		return 0.
	@staticmethod
	def f_0(q_sq, H_in, H_out):
		if H_in == "D":
			if H_out == "K":	return Meson_form_factors.f_0_template(q_sq, 0.1862721522, [0.6959, 0.775, 0.3], c.m_D, c.m_K)#FLAG Review 2021 M_pole = mDs0 = 2.317 GeV ,alt: return Meson_form_factors.f_D_H(q_sq, 0, 0.7647, 0.066, c.m_D, c.m_K) # ETM 17D
			if H_out == "Pi":	return Meson_form_factors.f_D_H(q_sq, 0.0342, 0.6117, 1.188, c.m_D, c.m_pi) # ETM 17D
		elif H_in == "Ds":
			if H_out == "Eta":	return Meson_form_factors.f_Ds_eta(q_sq, 0., 0.495, 0.) # 1805.08567v3
		elif H_in == "B": 
			if H_out == "D":	return Meson_form_factors.f_0_template(q_sq, 0.,[ 0.7821, -3.28], c.m_B, c.m_D)# FLAG Review 2021 M_pole = 0 , 17.608
			if H_out == "Pi":	return Meson_form_factors.f_0_template(q_sq, 0., [0.490, -1.61,-0.5072], c.m_B, c.m_pi)# FLAG Review 2021
		elif H_in =="Bs":
			if H_out == "Ds":	return Meson_form_factors.f_0_template(q_sq, 0.02496,[0.666,-0.26,-0.1], c.m_Bs, c.m_Ds) #FLAG Review 2021 M_pole = Bc* = 6.329 GeV ,-21828
			if H_out == "K":	return Meson_form_factors.f_0_template(q_sq, 0.0310,[0.2203, 0.089, 0.24], c.m_Bs, c.m_K)#FLAG Review 2021 M_pole = B*(0+) = 5.68 GeV ,13.548
		return 0.
	@staticmethod
	def g_HHv(q_sq, H, Hv): 
		if H == "D":
			if Hv=="Kstar": return Meson_form_factors.g_HHv_template(q_sq, c.m_D, c.m_Kstar, 1.03, 0.76, 0.66, 0.49, 0.27, 0.17, 0.30, 0.67, 0., 0., 0.20, 0.16, c.m_Dsstar) #1805.08567v3
		elif H == "B":
			# if Hv=="Dstar":	return Meson_form_factors.g_HHv_template(q_sq, c.m_B, c.m_Dstar, 0.76, 0.69, 0.66, 0.62, 0.57, 0.59, 0.78, 1.40, 0., 0., 0., 0.41, c.m_Bcstar) #1805.08567v3
			if Hv=="Dstar":return Meson_form_factors.V(q_sq,0.895, c.m_B, c.m_Dstar)/(c.m_B + c.m_Dstar)#1711.11013v3 (HPQCD17)
			if Hv=="Rho":	return Meson_form_factors.g_HHv_template(q_sq, c.m_B, c.m_rho, 0.295, 0.231, 0.269, 0.282, 0.875, 0.796, 0.54, 1.34, 0., 0.055, 0., -0.21, c.m_Bstar)
		elif H == "Bs":
			# if Hv=="Dsstar":return Meson_form_factors.g_HHv_template(q_sq, c.m_Bs, c.m_Dsstar, 0.95, 0.67, 0.70, 0.75, 0.372, 0.350, 0.463, 1.04, 0.561, 0.600, 0.510, 0.070, c.m_Bcstar) #1805.08567v3
			if Hv=="Dsstar":return Meson_form_factors.V(q_sq,0.883)/(c.m_Bs + c.m_Dsstar)#1711.11013v3 (HPQCD17)
			if Hv=="Kstar":	return Meson_form_factors.g_HHv_template(q_sq, c.m_Bs, c.m_Kstar, 0.291, 0.289, 0.287, 0.286, -0.516, -0.383, 0., 1.05, 2.10, 1.58, 1.06, -0.074, c.m_Bsstar) #1805.08567v3
		return 0.
	@staticmethod
	def f_HHv(q_sq, H, Hv):#all factors in this function are taken from 1805.08567v3
		if H == "D":
			if Hv=="Kstar": return Meson_form_factors.f_HHv_template(q_sq, c.m_D, c.m_Kstar, 1.03, 0.76, 0.66, 0.49, 0.27, 0.17, 0.30, 0.67, 0., 0., 0.20, 0.16, c.m_Dsstar) # 1805.08567v3
		elif H == "B":
			# if Hv=="Dstar":	return Meson_form_factors.f_HHv_template(q_sq, c.m_B, c.m_Dstar, 0.76, 0.69, 0.66, 0.62, 0.57, 0.59, 0.78, 1.40, 0., 0., 0., 0.41, c.m_Bcstar) #1805.08567v3
			if Hv=="Dstar":return Meson_form_factors.A1(q_sq,0.895, c.m_B, c.m_Dstar)*(c.m_B + c.m_Dstar)#1711.11013v3 (HPQCD17)
			if Hv=="Rho":	return Meson_form_factors.f_HHv_template(q_sq, c.m_B, c.m_rho, 0.295, 0.231, 0.269, 0.282, 0.875, 0.796, 0.54, 1.34, 0., 0.055, 0., -0.21, c.m_Bstar) #1805.08567v3
		elif H == "Bs":
			# if Hv=="Dsstar":return Meson_form_factors.f_HHv_template(q_sq, c.m_Bs, c.m_Dsstar, 0.95, 0.67, 0.70, 0.75, 0.372, 0.350, 0.463, 1.04, 0.561, 0.600, 0.510, 0.070, c.m_Bcstar) #1805.08567v3
			if Hv=="Dsstar":return Meson_form_factors.A1(q_sq,0.883, c.m_B, c.m_Dstar)*(c.m_Bs + c.m_Dsstar)#1711.11013v3 (HPQCD17)
			if Hv=="Kstar":	return Meson_form_factors.f_HHv_template(q_sq, c.m_Bs, c.m_Kstar, 0.291, 0.289, 0.287, 0.286, -0.516, -0.383, 0., 1.05, 2.10, 1.58, 1.06, -0.074, c.m_Bsstar) #1805.08567v3
		return 0.
	@staticmethod
	def a_plus_HHv(q_sq, H, Hv):
		if H == "D":
			if Hv=="Kstar": return Meson_form_factors.a_plus_HHv_template(q_sq, c.m_D, c.m_Kstar, 1.03, 0.76, 0.66, 0.49, 0.27, 0.17, 0.30, 0.67, 0., 0., 0.20, 0.16, c.m_Dsstar) #1805.08567v3
		elif H == "B":
			# if Hv=="Dstar":	return Meson_form_factors.a_plus_HHv_template(q_sq, c.m_B, c.m_Dstar, 0.76, 0.69, 0.66, 0.62, 0.57, 0.59, 0.78, 1.40, 0., 0., 0., 0.41, c.m_Bcstar) #1805.08567v3
			if Hv=="Dstar": return - Meson_form_factors.A2(q_sq,0.895, c.m_B, c.m_Dstar)/(c.m_B + c.m_Dstar)#1711.11013v3 (HPQCD17)
			if Hv=="Rho":	return Meson_form_factors.a_plus_HHv_template(q_sq, c.m_B, c.m_rho, 0.295, 0.231, 0.269, 0.282, 0.875, 0.796, 0.54, 1.34, 0., 0.055, 0., -0.21, c.m_Bstar) #1805.08567v3
		elif H == "Bs":
			# if Hv=="Dsstar":return Meson_form_factors.a_plus_HHv_template(q_sq, c.m_Bs, c.m_Dstar, 0.95, 0.67, 0.70, 0.75, 0.372, 0.350, 0.463, 1.04, 0.561, 0.600, 0.510, 0.070, c.m_Bcstar) #1805.08567v3
			if Hv=="Dsstar":return - Meson_form_factors.A2(q_sq,0.883)/(c.m_Bs + c.m_Dsstar)#1711.11013v3 (HPQCD17)
			if Hv=="Kstar":	return Meson_form_factors.a_plus_HHv_template(q_sq, c.m_Bs, c.m_Kstar, 0.291, 0.289, 0.287, 0.286, -0.516, -0.383, 0., 1.05, 2.10, 1.58, 1.06, -0.074, c.m_Bsstar) #1805.08567v3
		return 0.	
	@staticmethod
	def a_minus_HHv(q_sq, H, Hv):
		if H == "D":
			if Hv=="Kstar": return Meson_form_factors.a_minus_HHv_template(q_sq, c.m_D, c.m_Kstar, 1.03, 0.76, 0.66, 0.49, 0.27, 0.17, 0.30, 0.67, 0., 0., 0.20, 0.16, c.m_Dsstar,c.m_Ds) #1805.08567v3
		elif H == "B":
			# if Hv=="Dstar":	return Meson_form_factors.a_minus_HHv_template(q_sq, c.m_B, c.m_Dstar, 0.76, 0.69, 0.66, 0.62, 0.57, 0.59, 0.78, 1.40, 0., 0., 0., 0.41, c.m_Bcstar,c.m_Bc) #1805.08567v3
			if Hv=="Dstar": return (Meson_form_factors.A0(q_sq,0.895) * 2*c.m_Dstar - Meson_form_factors.f_HHv(q_sq,"B","Dstar") - (c.m_B**2 - c.m_Dstar**2)*Meson_form_factors.a_plus_HHv(q_sq,"B","Dstar"))/q_sq#1711.11013v3 (HPQCD17)
			if Hv=="Rho":	return Meson_form_factors.a_minus_HHv_template(q_sq, c.m_B, c.m_rho, 0.295, 0.231, 0.269, 0.282, 0.875, 0.796, 0.54, 1.34, 0., 0.055, 0., -0.21, c.m_Bstar, c.m_B) #1805.08567v3
		elif H == "Bs":
			# if Hv=="Dsstar":return Meson_form_factors.a_minus_HHv_template(q_sq, c.m_Bs, c.m_Dstar, 0.95, 0.67, 0.70, 0.75, 0.372, 0.350, 0.463, 1.04, 0.561, 0.600, 0.510, 0.070, c.m_Bcstar, c.m_Bc) #1805.08567v3
			if Hv=="Dsstar": return (Meson_form_factors.A0(q_sq,0.883) * 2*c.m_Dsstar - Meson_form_factors.f_HHv(q_sq,"Bs","Dsstar") - (c.m_Bs**2 - c.m_Dsstar**2)*Meson_form_factors.a_plus_HHv(q_sq,"Bs","Dsstar"))/q_sq#1711.11013v3 (HPQCD17)
			if Hv=="Kstar":	return Meson_form_factors.a_minus_HHv_template(q_sq, c.m_Bs, c.m_Kstar, 0.291, 0.289, 0.287, 0.286, -0.516, -0.383, 0., 1.05, 2.10, 1.58, 1.06, -0.074, c.m_Bsstar, c.m_Bs) #1805.08567v3
		return 0.
	
	@staticmethod
	def g_HHv_template(q_sq, m_H, m_Hv, f_V, f_A0, f_A1, f_A2, s_V, s_A0, s_A1, s_A2, xi_V, xi_A0, xi_A1, xi_A2, M_V):
		return f_V /((m_H+m_Hv)  * (1 - q_sq / M_V**2 ) * (1 - s_V * q_sq/M_V**2 - xi_V * q_sq**2 / M_V**4))
	@staticmethod
	def f_HHv_template(q_sq, m_H, m_Hv, f_V, f_A0, f_A1, f_A2, s_V, s_A0, s_A1, s_A2, xi_V, xi_A0, xi_A1, xi_A2, M_V):
		return f_A1 * (m_H+m_Hv) / ( 1 - s_A1 * q_sq / M_V**2 - xi_A1 * q_sq**2 / M_V**4)
	@staticmethod
	def a_plus_HHv_template(q_sq, m_H, m_Hv, f_V, f_A0, f_A1, f_A2, s_V, s_A0, s_A1, s_A2, xi_V, xi_A0, xi_A1, xi_A2, M_V):
		return - f_A2 / ((m_H+m_Hv) * ( 1 - s_A2 * q_sq / M_V**2 - xi_A2 * q_sq**2 / M_V**4))
	@staticmethod
	def a_minus_HHv_template(q_sq, m_H, m_Hv, f_V, f_A0, f_A1, f_A2, s_V, s_A0, s_A1, s_A2, xi_V, xi_A0, xi_A1, xi_A2, M_V, M_P):
		return (2 * m_Hv * f_A0 / ((1 - q_sq / M_P**2 ) * (1 - s_A0 * q_sq/M_V**2 - xi_A0 * q_sq**2 / M_V**4))
	  			- Meson_form_factors.f_HHv_template(q_sq, m_H, m_Hv, f_V, f_A0, f_A1, f_A2, s_V, s_A0, s_A1, s_A2, xi_V, xi_A0, xi_A1, xi_A2, M_V)
				- (m_H**2 - m_Hv**2) * Meson_form_factors.a_plus_HHv_template(q_sq, m_H, m_Hv, f_V, f_A0, f_A1, f_A2, s_V, s_A0, s_A1, s_A2, xi_V, xi_A0, xi_A1, xi_A2, M_V)
				) / q_sq
	@staticmethod
	def f_B(q_sq, P, a_0, a_1, a_2, m_H1, m_H2):
		t_plus = (m_H1+m_H2)**2
		t_0    = (m_H1+m_H2)*(np.sqrt(m_H1)-np.sqrt(m_H2))**2
		z = np.divide(np.sqrt(t_plus - q_sq) - np.sqrt(t_plus-t_0), np.sqrt(t_plus-q_sq) + np.sqrt(t_plus-t_0))
		return 1/(1.-P*q_sq)*( a_0 + a_1 * (z - 0.3333*z**3) + a_2 * (z**2 + 0.66 * z**3) )
	@staticmethod
	def f_D_H(q_sq, P, f_0, a, m_H1, m_H2):
		t_plus = (m_H1+m_H2)**2
		t_0    = (m_H1+m_H2)*(np.sqrt(m_H1)-np.sqrt(m_H2))**2
		z = np.divide(np.sqrt(t_plus - q_sq) - np.sqrt(t_plus-t_0), np.sqrt(t_plus-q_sq) + np.sqrt(t_plus-t_0))
		z0 =  np.divide(np.sqrt(t_plus) - np.sqrt(t_plus-t_0), np.sqrt(t_plus) + np.sqrt(t_plus-t_0))
		return 1/(1.-P*q_sq)*( f_0 - a * (z - z0) * (1 + 0.5 * (z+z0)) )
	@staticmethod
	def f_plus_template(q_sq, P, a_list, m_H1, m_H2):
		t_plus = (m_H1+m_H2)**2
		sqrt_dt = np.power(4*t_plus*m_H1*m_H2,0.25)
		z = np.divide(np.sqrt(t_plus - q_sq) - sqrt_dt, np.sqrt(t_plus-q_sq) + sqrt_dt)
		N = len(a_list); sum_ = 0.
		for n in range(N): sum_ += a_list[n] * (z**n -(-1)**(n-N) * n/N * z**N)
		return sum_/(1.-P*q_sq)
	@staticmethod
	def f_0_template(q_sq, P, a_list, m_H1, m_H2):
		t_plus = (m_H1+m_H2)**2
		sqrt_dt = np.power(4*t_plus*m_H1*m_H2,0.25)
		z = np.divide(np.sqrt(t_plus - q_sq) - sqrt_dt, np.sqrt(t_plus-q_sq) + sqrt_dt)
		N = len(a_list); sum_ = 0.
		for n in range(N): sum_ += a_list[n] * z**n
		return sum_/(1.-P*q_sq)
	@staticmethod
	def f_Ds_eta(q_sq, P, f_0, a):
		return  np.divide(f_0, (1.-P*q_sq) *( 1 - a * P*q_sq))

	@staticmethod
	def hA1(w,h_A1_1):
		z = (np.sqrt(w+1)-np.sqrt(2)) / (np.sqrt(w+1)+np.sqrt(2) )
		rho_sq = 1.23 # 1.17 #
		# hA1_1 = 0.883#0.902
		return h_A1_1 *( 1- 8*rho_sq*z + (53*rho_sq-15)*z**2 - (231*rho_sq - 91)*z**3)
	@staticmethod
	def R0(w):
		R0_1 = 1.057
		return R0_1 -  0.11*(w - 1) + 0.01*(w - 1)**2
	@staticmethod
	def R1(w):
		R1_1 = 1.52#1.40#
		return R1_1 -  0.12*(w - 1) + 0.05*(w - 1)**2
	@staticmethod
	def R2(w):
		R2_1 = 0.93 #0.85#
		return R2_1 -  0.11*(w - 1) + 0.06*(w - 1)**2
	@staticmethod
	def A0(q_sq,h_A1_1, mH = c.m_Bs, mHv = c.m_Dsstar):
		r  = mHv / mH
		RDs = 2*np.sqrt(r)/ (1+r)
		w = (mH**2+ mHv**2 - q_sq)/(2*mH*mHv)
		return Meson_form_factors.hA1(w,h_A1_1)*Meson_form_factors.R0(w)/RDs
	@staticmethod
	def A1(q_sq,h_A1_1, mH = c.m_Bs, mHv = c.m_Dsstar):
		r  = mHv / mH
		RDs = 2*np.sqrt(r)/ (1+r)
		w = (mH**2+ mHv**2 - q_sq)/(2*mH*mHv)
		return Meson_form_factors.hA1(w,h_A1_1)*RDs * 2 / (w+1)
	@staticmethod
	def A2(q_sq,h_A1_1, mH = c.m_Bs, mHv = c.m_Dsstar):
		r  = mHv / mH
		RDs = 2*np.sqrt(r)/ (1+r)
		w = (mH**2+ mHv**2 - q_sq)/(2*mH*mHv)
		return Meson_form_factors.hA1(w, h_A1_1)*Meson_form_factors.R2(w)/RDs
	@staticmethod
	def V(q_sq, h_A1_1,mH = c.m_Bs, mHv = c.m_Dsstar):
		r  = mHv / mH
		RDs = 2*np.sqrt(r)/ (1+r)
		w = (mH**2+ mHv**2 - q_sq)/(2*mH*mHv)
		return Meson_form_factors.hA1(w,h_A1_1)*Meson_form_factors.R1(w)/RDs

class alp_production_BRs:
	"""Branching ratios of pseudoscalar meson decays into axion-like-particles"""
	@staticmethod
	def get_BRs(parents, daughters, alp_mass):
		branchings = []
		for parent, daughter in zip(parents, daughters):
			branchings.append(alp_production_BRs(parent, daughter, alp_mass))
		return np.array(branchings)	
	@staticmethod
	def branching_ratios(parent, daughter, alp_mass, g = 1e-4):
		if parent == "Bmeson":
			if 	 daughter == "K":      return alp_production_BRs.B_K_a(alp_mass, g_bs=g)
			elif daughter == "Kstar" : return alp_production_BRs.B_Kstar_a(alp_mass, g_bs=g)
		if parent == "B0meson":
			if 	 daughter == "K0" :    return alp_production_BRs.B0_K0_a(alp_mass, g_bs=g)
			elif daughter == "K0star": return alp_production_BRs.B0_K0star_a(alp_mass, g_bs=g)
		if parent == "Dmeson":
			if 	 daughter == "Pi":      return alp_production_BRs.D_pi_a(alp_mass, g_cu=g)
		if parent == "D0meson":
			if 	 daughter == "Pi0" :    return alp_production_BRs.D0_pi0_a(alp_mass, g_cu=g)
		return 0.
	@staticmethod
	def B_K_a(m_a, g_bs): #B -> K + ALP
		f_0 = 0.330/(1-m_a**2/37.46) #hep-ph/0406232
		if c.m_B > m_a + c.m_K: return np.abs(g_bs)**2 /(64*np.pi*c.m_B**3)*(c.m_B**2 - c.m_K**2)**2 * np.sqrt(lambda_Kallen(c.m_B,c.m_K,m_a)) * f_0**2 * c.tau_B
		return 0.
	@staticmethod
	def B0_K0_a(m_a, g_bs): #B0 -> K0 + ALP
		f_0 = 0.330/(1-m_a**2/37.46) #hep-ph/0406232
		if c.m_B0 > m_a + c.m_K0: return np.abs(g_bs)**2 /(64*np.pi*c.m_B0**3)*(c.m_B0**2 - c.m_K0**2)**2 * np.sqrt(lambda_Kallen(c.m_B0,c.m_K0,m_a)) * f_0**2 * c.tau_B0
		return 0.
	@staticmethod
	def B_Kstar_a(m_a, g_bs): #B -> K* + ALP
		A_0 = 1.364/(1-m_a**2/(c.m_B**2)) - 0.990/(1-m_a**2/36.78) #hep-ph/0412079
		if c.m_B > m_a + c.m_Kstar: return np.abs(g_bs)**2 /(64*np.pi*c.m_B**3) * np.sqrt(lambda_Kallen(c.m_B,c.m_Kstar,m_a)**3) * A_0**2 * c.tau_B
		return 0.
	@staticmethod
	def B0_K0star_a(m_a, g_bs): #B0 -> K0* + ALP
		A_0 = 1.364/(1-m_a**2/(c.m_B**2)) - 0.990/(1-m_a**2/36.78) #hep-ph/0412079
		if c.m_B0 > m_a + c.m_K0star: return np.abs(g_bs)**2 /(64*np.pi*c.m_B0**3) * np.sqrt(lambda_Kallen(c.m_B0,c.m_K0star,m_a)**3) * A_0**2 * c.tau_B0
		return 0.
	@staticmethod
	def D_pi_a(m_a, g_cu): #D -> pi + ALP
		f_0 = 0.612/(1-m_a**2/6.46)
		if c.m_D > m_a + c.m_pi: return np.abs(g_cu)**2 /(64*np.pi*c.m_D**3)*(c.m_D**2 - c.m_pi**2)**2 * np.sqrt(lambda_Kallen(c.m_D,c.m_pi,m_a)) * f_0**2 * c.tau_D
		return 0.
	@staticmethod
	def D0_pi0_a(m_a, g_cu): #D0 -> pi0 + ALP
		f_0 = 0.612/(1-m_a**2/6.46)
		if c.m_D0 > m_a + c.m_pi0: return np.abs(g_cu)**2 /(64*np.pi*c.m_D0**3)*(c.m_D0**2 - c.m_pi0**2)**2 * np.sqrt(lambda_Kallen(c.m_D0,c.m_pi0,m_a)) * f_0**2 * c.tau_D0
		return 0.

class hnl_production_BRs:
	"""Branching ratios and singly differential weights of pseudoscalar meson decays into heavy neutral leptons taken from JHEP, 0710""" 
	@staticmethod
	def get_BRs(parents, daughters, hnl_mass):
		branchings = []
		for parent, daughter in zip(parents, daughters):
			branchings.append(hnl_production_BRs.branching_ratios(parent, daughter, hnl_mass))
		return np.array(branchings)
	@staticmethod
	def branching_ratios(parent, daughter, hnl_mass, U2 = 1):
		if isinstance(daughter,str):
			if parent == "Dmeson":
				if daughter == "El": return hnl_production_BRs.Meson_to_LN(hnl_mass, c.m_el, c.m_D, 1, c.f_DV_cd, c.tau_D, U_l2 = U2)
				if daughter == "Mu": return hnl_production_BRs.Meson_to_LN(hnl_mass, c.m_mu, c.m_D, 1, c.f_DV_cd, c.tau_D, U_l2 = U2)
			elif parent == "Dsmeson":
				if daughter == "El": return hnl_production_BRs.Meson_to_LN(hnl_mass, c.m_el, c.m_Ds, 1, c.f_DsV_cs, c.tau_Ds, U_l2 = U2)
				if daughter == "Mu": return hnl_production_BRs.Meson_to_LN(hnl_mass, c.m_mu, c.m_Ds, 1, c.f_DsV_cs, c.tau_Ds, U_l2 = U2)
				if daughter == "Tau":return hnl_production_BRs.Meson_to_LN(hnl_mass, c.m_tau, c.m_Ds, 1, c.f_DsV_cs, c.tau_Ds, U_l2 = U2)
			elif parent == "Tau":
				if daughter == "Pi": return hnl_production_BRs.tau_to_HN(hnl_mass, c.m_pi, c.V_ud, c.f_pi, U2_tau=U2)
				if daughter == "K" : return hnl_production_BRs.tau_to_HN(hnl_mass, c.m_K, c.V_us,  c.f_K, U2_tau=U2) 
				if daughter == "Rho":return hnl_production_BRs.tau_to_rhoN(hnl_mass, U2_tau=U2) 
			elif parent == "Bmeson": 
				if daughter == "El": return hnl_production_BRs.Meson_to_LN(hnl_mass, c.m_el, c.m_B, 1, c.f_BV_ub, c.tau_B, U_l2 = U2)
				if daughter == "Mu": return hnl_production_BRs.Meson_to_LN(hnl_mass, c.m_mu, c.m_B, 1, c.f_BV_ub, c.tau_B, U_l2 = U2)
				if daughter == "Tau":return hnl_production_BRs.Meson_to_LN(hnl_mass, c.m_tau,c.m_B, 1, c.f_BV_ub, c.tau_B, U_l2 = U2)
		else: # Here Clebsh-Gordan Coefficients to be respected (factor 1/2 for H->HvLN decays into with H == pi0 || Hv == rho0) 
			if parent == "D0meson":
				if daughter == ["Pi","El"] or daughter == ["El","Pi"]: 		return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_el, "D", "Pi", c.m_D0, c.m_pi, c.V_cd, c.tau_D0, U_l2 = U2)
				if daughter == ["K","El"] or daughter == ["El","K"]: 		return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_el, "D", "K", c.m_D0, c.m_K, c.V_cs, c.tau_D0, U_l2 = U2)
				if daughter == ["Kstar","El"] or daughter == ["El","Kstar"]:return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_el, "D", "Kstar", c.m_D0, c.m_Kstar, c.V_cs, c.tau_D0, U_l2 = U2)
				if daughter == ["Pi","Mu"] or daughter == ["Mu","Pi"]: 		return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_mu, "D", "Pi", c.m_D0, c.m_pi, c.V_cd, c.tau_D0, U_l2 = U2)
				if daughter == ["K","Mu"] or daughter == ["Mu","K"]: 		return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_mu, "D", "K", c.m_D0, c.m_K, c.V_cs, c.tau_D0, U_l2 = U2)
				if daughter == ["Kstar","Mu"] or daughter == ["Mu","Kstar"]:return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_mu, "D", "Kstar", c.m_D0, c.m_Kstar, c.V_cs, c.tau_D0, U_l2 = U2)
			if parent == "Dmeson":
				if daughter == ["Pi0","El"] or daughter == ["El","Pi0"]: 		return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_el, "D", "Pi", c.m_D, c.m_pi0, c.V_cd, c.tau_D, U_l2 = U2)/2
				if daughter == ["K0","El"] or daughter == ["El","K0"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_el, "D", "K", c.m_D, c.m_K0, c.V_cs, c.tau_D, U_l2 = U2)
				if daughter == ["K0star","El"] or daughter == ["El","K0star"]: 	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_el, "D", "Kstar", c.m_D, c.m_K0star, c.V_cs, c.tau_D, U_l2 = U2)
				if daughter == ["Pi0","Mu"] or daughter == ["Mu","Pi0"]: 		return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_mu, "D", "Pi", c.m_D, c.m_pi0, c.V_cd, c.tau_D, U_l2 = U2)/2
				if daughter == ["K0","Mu"] or daughter == ["Mu","K0"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_mu, "D", "K", c.m_D, c.m_K0, c.V_cs, c.tau_D, U_l2 = U2)
				if daughter == ["K0star","Mu"] or daughter == ["Mu","K0star"]: 	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_mu, "D", "Kstar", c.m_D, c.m_K0star, c.V_cs, c.tau_D, U_l2 = U2)
			elif parent == "Dsmeson":
				if daughter == ["Eta","El"] or daughter == ["El","Eta"]: return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_el, "D", "Eta", c.m_Ds, c.m_eta, c.V_cs, c.tau_Ds, U_l2 = U2)
				if daughter == ["Eta","Mu"] or daughter == ["Mu","Eta"]: return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_mu, "D", "Eta", c.m_Ds, c.m_eta, c.V_cs, c.tau_Ds, U_l2 = U2)
			elif parent  == "B0meson":
				if daughter == ["Pi", "El"] or daughter == ["El", "Pi"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_el, "B", "Pi", c.m_B0, c.m_pi, c.V_ub, c.tau_B0, U_l2 = U2)
				elif daughter == ["Rho", "El"] or daughter == ["El", "Rho"]:		return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_el, "B", "Rho", c.m_B0, c.m_rho, c.V_ub, c.tau_B0, U_l2 = U2)
				elif daughter == ["D", "El"] or daughter == ["El", "D"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_el, "B", "D", c.m_B0, c.m_D, c.V_cb, c.tau_B0, U_l2 = U2)
				elif daughter == ["Dstar", "El"] or daughter == ["El", "Dstar"]:	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_el, "B", "Dstar", c.m_B0, c.m_Dstar, c.V_cb, c.tau_B0, U_l2 = U2)
				elif daughter == ["Pi", "Mu"] or daughter == ["Mu", "Pi"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_mu, "B", "Pi", c.m_B0, c.m_pi, c.V_ub, c.tau_B0, U_l2 = U2)
				elif daughter == ["Rho", "Mu"] or daughter == ["Mu", "Rho"]:		return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_mu, "B", "Rho", c.m_B0, c.m_rho, c.V_ub, c.tau_B0, U_l2 = U2)
				elif daughter == ["D", "Mu"] or daughter == ["Mu", "D"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_mu, "B", "D", c.m_B0, c.m_D, c.V_cb, c.tau_B0, U_l2 = U2)
				elif daughter == ["Dstar", "Mu"] or daughter == ["Mu", "Dstar"]:	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_mu, "B", "Dstar", c.m_B0, c.m_Dstar, c.V_cb, c.tau_B0, U_l2 = U2)
				elif daughter == ["Pi", "Tau"] or daughter == ["Tau", "Pi"]: 		return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_tau, "B", "Pi", c.m_B0, c.m_pi, c.V_ub, c.tau_B0, U_l2 = U2)
				elif daughter == ["Rho", "Tau"] or daughter == ["Tau", "Rho"]:		return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_tau, "B", "Rho", c.m_B0, c.m_rho, c.V_ub, c.tau_B0, U_l2 = U2)
				elif daughter == ["D", "Tau"] or daughter == ["Tau", "D"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_tau, "B", "D", c.m_B0, c.m_D, c.V_cb, c.tau_B0, U_l2 = U2)
				elif daughter == ["Dstar", "Tau"] or daughter == ["Tau", "Dstar"]:	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_tau, "B", "Dstar", c.m_B0, c.m_Dstar, c.V_cb, c.tau_B0, U_l2 = U2)
			elif parent == "Bmeson": 
				if daughter == ["Pi0", "El"] or daughter == ["El", "Pi0"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_el, "B", "Pi", c.m_B, c.m_pi0, c.V_ub, c.tau_B, U_l2 = U2)/2
				elif daughter == ["Rho0", "El"] or daughter == ["El", "Rho0"]:		return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_el, "B", "Rho", c.m_B, c.m_rho, c.V_ub, c.tau_B, U_l2 = U2)/2
				elif daughter == ["D0", "El"] or daughter == ["El", "D0"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_el, "B", "D", c.m_B, c.m_D0, c.V_cb, c.tau_B, U_l2 = U2)
				elif daughter == ["D0star", "El"] or daughter == ["El", "D0star"]:	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_el, "B", "Dstar", c.m_B, c.m_D0star, c.V_cb, c.tau_B, U_l2 = U2)
				elif daughter == ["Pi0", "Mu"] or daughter == ["Mu", "Pi0"]: 		return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_mu, "B", "Pi", c.m_B, c.m_pi0, c.V_ub, c.tau_B, U_l2 = U2)/2
				elif daughter == ["Rho0", "Mu"] or daughter == ["Mu", "Rho0"]:		return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_mu, "B", "Rho", c.m_B, c.m_rho, c.V_ub, c.tau_B, U_l2 = U2)/2
				elif daughter == ["D0", "Mu"] or daughter == ["Mu", "D0"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_mu, "B", "D", c.m_B, c.m_D0, c.V_cb, c.tau_B, U_l2 = U2)
				elif daughter == ["D0star", "Mu"] or daughter == ["Mu", "D0star"]:	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_mu, "B", "Dstar", c.m_B, c.m_D0star, c.V_cb, c.tau_B, U_l2 = U2)
				elif daughter == ["Pi0", "Tau"] or daughter == ["Tau", "Pi0"]: 		return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_tau, "B", "Pi", c.m_B, c.m_pi0, c.V_ub, c.tau_B, U_l2 = U2)/2
				elif daughter == ["Rho0", "Tau"] or daughter == ["Tau", "Rho0"]:	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_tau, "B", "Rho", c.m_B, c.m_rho, c.V_ub, c.tau_B, U_l2 = U2)/2
				elif daughter == ["D0", "Tau"] or daughter == ["Tau", "D0"]: 		return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_tau, "B", "D", c.m_B, c.m_D0, c.V_cb, c.tau_B, U_l2 = U2)
				elif daughter == ["D0star", "Tau"] or daughter == ["Tau", "D0star"]:return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_tau, "B", "Dstar", c.m_B, c.m_D0star, c.V_cb, c.tau_B, U_l2 = U2)
			elif parent == "Bsmeson":
				if daughter == ["K", "El"] or daughter == ["El", "K"]: 				return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_el, "Bs", "K", c.m_Bs, c.m_K, c.V_ub, c.tau_Bs, U_l2 = U2)
				elif daughter == ["Kstar", "El"] or daughter == ["El", "Kstar"]:	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_el, "Bs", "Kstar", c.m_Bs, c.m_Kstar, c.V_ub, c.tau_Bs, U_l2 = U2)
				elif daughter == ["Ds", "El"] or daughter == ["El", "Ds"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_el, "Bs", "D", c.m_Bs, c.m_Ds, c.V_cb, c.tau_Bs, U_l2 = U2)
				elif daughter == ["Dsstar", "El"] or daughter == ["El", "Dsstar"]:	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_el, "Bs", "Dsstar", c.m_Bs, c.m_Dsstar, c.V_cb, c.tau_Bs, U_l2 = U2)
				elif daughter == ["K", "Mu"] or daughter == ["Mu", "K"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_mu, "Bs", "K", c.m_Bs, c.m_K, c.V_ub, c.tau_Bs, U_l2 = U2)
				elif daughter == ["Kstar", "Mu"] or daughter == ["Mu", "Kstar"]:	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_mu, "Bs", "Kstar", c.m_Bs, c.m_Kstar, c.V_ub, c.tau_Bs, U_l2 = U2)
				elif daughter == ["Ds", "Mu"] or daughter == ["Mu", "Ds"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_mu, "Bs", "D", c.m_Bs, c.m_Ds, c.V_cb, c.tau_Bs, U_l2 = U2)
				elif daughter == ["Dsstar", "Mu"] or daughter == ["Mu", "Dsstar"]:	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_mu, "Bs", "Dsstar", c.m_Bs, c.m_Dsstar, c.V_cb, c.tau_Bs, U_l2 = U2)
				elif daughter == ["K", "Tau"] or daughter == ["Tau", "K"]: 			return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_tau, "Bs", "K", c.m_Bs, c.m_K, c.V_ub, c.tau_Bs, U_l2 = U2)
				elif daughter == ["Kstar", "Tau"] or daughter == ["Tau", "Kstar"]:	return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_tau, "Bs", "Kstar", c.m_Bs, c.m_Kstar, c.V_ub, c.tau_Bs, U_l2 = U2)
				elif daughter == ["Ds", "Tau"] or daughter == ["Tau", "Ds"]: 		return hnl_production_BRs.meson_to_HLN(hnl_mass, c.m_tau, "Bs", "D", c.m_Bs, c.m_Ds, c.V_cb, c.tau_Bs, U_l2 = U2)
				elif daughter == ["Dsstar", "Tau"] or daughter == ["Tau", "Dsstar"]:return hnl_production_BRs.meson_to_HvLN(hnl_mass, c.m_tau, "Bs", "Dsstar", c.m_Bs, c.m_Dsstar, c.V_cb, c.tau_Bs, U_l2 = U2)
			elif parent == "Tau": # The case of tau->nu_tau l_alpha N is much less likely than 3 body D meson production in the low mass regime 
				if daughter == ["El", "Nu"] or daughter == ["Nu", "El"]: 	return hnl_production_BRs.tau_to_NLaNua(hnl_mass, c.m_el, U2_tau = U2)
				elif daughter == ["Mu", "Nu"] or daughter == ["Nu", "Mu"]: 	return hnl_production_BRs.tau_to_NLaNua(hnl_mass, c.m_mu, U2_tau = U2)
		print(f"[Warning:] \t Branching {parent} to {daughter} not parsed")
		return 0.
	@staticmethod
	def Meson_to_LN(m_N, m_l, m_H, V_H, f_H, tau_H, U_l2=1.):
		if m_l + m_N < m_H : 
			y_l, y_N = [m_l/m_H,m_N/m_H]
			return U_l2 * np.divide( c.G_F**2 * f_H**2 * V_H**2 * m_H**3 * tau_H, 8*np.pi) * (y_l**2 + y_N**2 - (y_N**2-y_l**2)**2)*np.sqrt(lambda_Kallen(1,y_N,y_l))
		return 0.
	@staticmethod
	def meson_to_HLN(m_N, m_l, H1, H2, m_H1, m_H2, V_HH, tau_H1, U_l2=1., n_intpoints = 1000):
		if m_H2 + m_l + m_N < m_H1 : 
			# f = lambda xi: hnl_production_BRs.diff_weight_meson_to_HLN(m_N, m_l, xi*m_H1**2, H1, H2, m_H1, m_H2)
			# Integral = integrate.quad(f, ((m_l+m_N)/m_H1)**2, ((m_H1-m_H2)/m_H1)**2 )[0]
			q2s = np.linspace((m_l+m_N)**2, (m_H1-m_H2)**2, n_intpoints+2, dtype=np.double)[1:-1]
			Integrand = hnl_production_BRs.diff_weight_meson_to_HLN(m_N, m_l, q2s, H1, H2, m_H1, m_H2)
			return U_l2 * np.divide(c.G_F**2 * tau_H1 * m_H1**5* V_HH**2, 64*np.pi**3) * np.trapz(Integrand, dx = (q2s[-2] - q2s[1]) / (n_intpoints*m_H1**2))#* Integral
		return 0.
	@staticmethod
	def meson_to_HvLN(m_N, m_l, H, Hv, m_H, m_Hv, V_HHv, tau_H, U_l2=1., n_intpoints = 1000):
		if m_Hv + m_l + m_N < m_H : 
			# f = lambda xi: hnl_production_BRs.diff_weight_meson_to_HvLN(m_N, m_l, xi*m_H**2, H, Hv, m_H, m_Hv)
			# Integral = integrate.quad(f, ((m_l+m_N)/m_H)**2, ((m_H-m_Hv)/m_H)**2 )[0]
			q2s = np.linspace((m_l+m_N)**2, (m_H-m_Hv)**2, n_intpoints+2, dtype=np.double)[1:-1]
			Integrand = hnl_production_BRs.diff_weight_meson_to_HvLN(m_N, m_l, q2s, H, Hv, m_H, m_Hv)
			return U_l2 * np.divide(c.G_F**2 * m_H**7 * tau_H * V_HHv**2, 64 * np.pi**3* m_Hv**2 ) * np.trapz(Integrand, dx = (q2s[-2] - q2s[1]) / (n_intpoints*m_H**2))#* Integral
		return 0.
	@staticmethod
	def tau_to_HN(m_N, m_H, V_H, f_H, U2_tau=1.):
		if m_H+m_N < c.m_tau: 
				return U2_tau * np.divide(c.tau_tau * c.G_F**2 * f_H**2 * V_H**2 * c.m_tau**3, 16*np.pi) * (
				(1.-m_N**2/c.m_tau**2)**2 - m_H**2/c.m_tau**2 * (1 + m_N**2/c.m_tau**2)) * np.sqrt(
				(1.-(m_H-m_N)**2/c.m_tau**2)*(1.- (m_H+m_N)**2/c.m_tau**2))
		return 0.
	@staticmethod
	def tau_to_rhoN(m_N, U2_tau=1.):
		if c.m_rho+m_N< c.m_tau: 
				return U2_tau * np.divide(c.tau_tau * c.G_F**2 * c.g_rho**2 * c.V_ud**2 * c.m_tau**3, 8*np.pi) * (
				(1.-m_N**2/c.m_tau**2)**2 + c.m_rho**2/c.m_tau**2 * (1 + (m_N**2 - 2*c.m_rho**2)/c.m_tau**2)) * np.sqrt(
				(1.-(c.m_rho-m_N)**2/c.m_tau**2)*(1.- (c.m_rho+m_N)**2/c.m_tau**2))
		return 0.
	@staticmethod
	def tau_to_NLaNua(m_N, m_l, U2_tau=1.):
		if m_l + m_N < c.m_tau:
			y_l = m_l / c.m_tau
			y_N = m_N / c.m_tau
			xis = np.linspace(y_l**2, (1.-y_N)**2, dtype=np.double)
			Integrant  = np.divide(
					(xis-y_l**2)**2 * np.sqrt(lambda_Kallen(1, np.sqrt(xis), y_N)) * ( (xis+2*y_l**2)*(1-y_N**2)**2 + xis*(xis-y_l**2)*(1+y_N**2-y_l**2) - xis*y_l**4 - 2*xis**3), xis**3
			)
			return U2_tau * np.divide(c.G_F**2*c.m_tau**5*c.tau_tau, 96.*np.pi**3) *np.trapz(Integrant, x=xis)
		return 0.
	@staticmethod
	def tau_to_NutauLN(m_N, m_l, U2_l=1.):
		if m_l + m_N < c.m_tau:
			y_l = m_l / c.m_tau
			y_N = m_N / c.m_tau
			xis = np.linspace((y_l+y_N)**2, 1., dtype=np.double)
			Integrant  = np.divide(
					(1-xis)**2 * np.sqrt(lambda_Kallen(np.sqrt(xis), y_N, y_l)) * ( 2*xis**3 + xis**2 + xis*(1-xis)*(y_N**2+y_l**2) - (2+xis)*(y_N**2-y_l**2)**2 ), xis**3
			)
			return U2_l * np.divide(c.G_F**2*c.m_tau**5*c.tau_tau, 96.*np.pi**3) *np.trapz(Integrant, x=xis)
		return 0.
	@staticmethod
	def get_differential_weight_function(parent, daughters):
		return lambda m_N, M12_2, M23_2: hnl_production_BRs.differential_weights(parent, daughters, M12_2, M23_2, m_N)
	@staticmethod
	def differential_weights(parent, daughters, M12_2, M23_2, hnl_mass):
		if parent == "D0meson":
			if daughters == ["El","Pi"]: 	 return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_el, M12_2, "D", "Pi", c.m_D0, c.m_pi)
			if daughters == ["El","K"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_el, M12_2,  "D", "K", c.m_D0, c.m_K)
			if daughters == ["El","Kstar"]: return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_el, M12_2, "D", "Kstar", c.m_D0, c.m_Kstar)
			if daughters == ["El","Pi"]: 	 return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_el, M12_2, "D", "Pi", c.m_D0, c.m_pi)
			if daughters == ["Mu","K"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_mu, M12_2,  "D", "K", c.m_D0, c.m_K)
			if daughters == ["Mu","Kstar"]: return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_mu, M12_2, "D", "Kstar", c.m_D0, c.m_Kstar)
		elif parent == "Dmeson":
			if daughters == ["El","Pi0"]: 	 return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_el, M12_2,  "D", "Pi", c.m_D, c.m_pi0)
			if daughters == ["El","K0",]: 	 return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_el, M12_2,  "D", "K", c.m_D, c.m_K0)
			if daughters == ["El","K0star"]: return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_el, M12_2, "D", "Kstar", c.m_D, c.m_K0star)
			if daughters == ["Mu","Pi0"]: 	 return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_mu, M12_2,  "D", "Pi", c.m_D, c.m_pi0) 
			if daughters == ["Mu","K0"]:	 return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_mu, M12_2,  "D", "K", c.m_D, c.m_K0)
			if daughters == ["Mu","K0star"]: return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_mu, M12_2, "D", "Kstar", c.m_D, c.m_K0star)
		elif parent == "Dsmeson":
			if daughters == ["El","Eta"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_el, M12_2,  "D", "Eta", c.m_Ds, c.m_eta)
			if daughters == ["Mu","Eta"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_mu, M12_2,  "D", "Eta", c.m_Ds, c.m_eta)
		elif parent  == "B0meson":
			if daughters == ["El","D"]:			return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_el,  M12_2, "B", "D", c.m_B0, c.m_D)
			elif daughters == ["El","Dstar"]:	return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_el, M12_2, "B", "Dstar", c.m_B0, c.m_Dstar)
			elif daughters == ["El", "Pi"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_el, M12_2, "B", "Pi", c.m_B0, c.m_pi)
			elif daughters == ["Mu","D"]:		return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_mu,  M12_2, "B", "D", c.m_B0, c.m_D)
			elif daughters == ["Mu","Dstar"]:	return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_mu, M12_2, "B", "Dstar", c.m_B0, c.m_Dstar)
			elif daughters == ["Mu", "Pi"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_mu, M12_2, "B", "Pi", c.m_B0, c.m_pi)
			elif daughters == ["Tau","D"]: 		return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_tau, M12_2, "B", "D", c.m_B0, c.m_D)
			elif daughters == ["Tau","Dstar"]:	return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_tau,M12_2, "B", "Dstar", c.m_B0, c.m_Dstar)
			elif daughters == ["Tau", "Pi"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_tau, M12_2, "B", "Pi", c.m_B0, c.m_pi)
		elif parent == "Bmeson": 
			if daughters == ["El","D0"]: 		return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_el, M12_2, "B", "D", c.m_B, c.m_D0)
			elif daughters == ["El", "D0star"]:	return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_el,M12_2, "B", "Dstar", c.m_B, c.m_D0star)
			elif daughters == ["El", "Pi0"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_el, M12_2, "B", "Pi", c.m_B, c.m_pi0)
			elif daughters == ["Mu", "D0"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_mu, M12_2, "B", "D", c.m_B, c.m_D0)
			elif daughters == ["Mu", "D0star"]:	return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_mu,M12_2, "B", "Dstar", c.m_B, c.m_D0star)
			elif daughters == ["Mu", "Pi0"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_mu, M12_2, "B", "Pi", c.m_B, c.m_pi0)
			elif daughters == ["Tau", "D0"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_tau,M12_2, "B", "D", c.m_B, c.m_D0)
			elif daughters == ["Tau", "D0star"]:return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_tau,M12_2, "B", "Dstar", c.m_B, c.m_D0star)
			elif daughters == ["Tau", "Pi0"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_tau, M12_2, "B", "Pi", c.m_B, c.m_pi0)
		elif parent == "Bsmeson":
			if daughters == ["El","Ds"]: 		return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_el, M12_2, "B", "D", c.m_Bs, c.m_Ds)
			elif daughters == ["El","Dsstar"]:	return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_el,M12_2,  "B", "Dstar", c.m_Bs, c.m_Dsstar)
			elif daughters == ["Mu", "Ds"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_mu, M12_2, "B", "D", c.m_Bs, c.m_Ds)
			elif daughters == ["Mu", "Dsstar"]:	return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_mu,M12_2, "B", "Dstar", c.m_Bs, c.m_Dsstar)
			elif daughters == ["Tau", "Ds"]: 	return hnl_production_BRs.diff_weight_meson_to_HLN(hnl_mass, c.m_tau,M12_2,  "B", "D", c.m_Bs, c.m_Ds)
			elif daughters == ["Tau", "Dsstar"]:return hnl_production_BRs.diff_weight_meson_to_HvLN(hnl_mass, c.m_tau,M12_2,  "B", "Dstar", c.m_Bs, c.m_Dsstar)
		elif parent == "Tau": # The case of tau->nu_tau l_alpha N is much less likely than 3 body D meson production in the low mass regime 
			if daughters == ["El", "Nu"]: 	return hnl_production_BRs.diff_weight_tau_to_NLaNua(hnl_mass, c.m_el, M12_2, M23_2)
			elif daughters == ["Mu", "Nu"]: return hnl_production_BRs.diff_weight_tau_to_NLaNua(hnl_mass, c.m_mu, M12_2, M23_2)
		return 0.
	@staticmethod
	def diff_weight_meson_to_HLN(m_N, m_l, MNl_2, H1, H2, m_H1, m_H2):
		if m_H2 + m_l + m_N < m_H1 : 
			y_H = m_H2/m_H1
			y_l = m_l/m_H1
			y_N = m_N/m_H1
			q2 = MNl_2
			xi = q2 / m_H1**2

			f_plus_sq = np.abs(Meson_form_factors.f_plus(q2, H1, H2 ))**2
			f_0_sq =    np.abs(Meson_form_factors.f_0(q2, H1, H2 ))**2
			lambda_lN = lambda_Kallen(y_N, np.sqrt(xi),y_l)
			lambda_H = lambda_Kallen(1.,y_H, np.sqrt(xi))
			L = np.sqrt(lambda_H*lambda_lN)
			G_minus= xi*(y_N**2+y_l**2) - (y_N**2 - y_l**2)**2

			return np.divide((0.3333*f_plus_sq**2*L**3 + 0.5*f_plus_sq*L*G_minus*lambda_Kallen(1.,y_H, np.sqrt(xi)) + 0.5*f_0_sq*L*G_minus*(1.-y_H**2)**2), xi**3) #
		return 0.
	@staticmethod
	def diff_weight_meson_to_HvLN(m_N, m_l, MNl_2, H, Hv, m_H, m_Hv):
		if m_Hv + m_l + m_N < m_H : 
			y_H = m_Hv/m_H
			y_l = m_l/m_H
			y_N = m_N/m_H
			q2 = MNl_2
			xi = q2 / m_H**2

			f = Meson_form_factors.f_HHv(q2, H, Hv )
			gs = Meson_form_factors.g_HHv(q2, H, Hv )
			a_plus = Meson_form_factors.a_plus_HHv(q2, H, Hv )
			a_minus= Meson_form_factors.a_minus_HHv(q2, H, Hv )

			lambda_lN = lambda_Kallen(y_N, np.sqrt(xi),y_l)
			lambda_H = lambda_Kallen(1.,y_H, np.sqrt(xi))
			L = np.sqrt(lambda_H*lambda_lN)
			G_plus = xi*(y_N**2+y_l**2) + (y_N**2 - y_l**2)**2
			G_minus= xi*(y_N**2+y_l**2) - (y_N**2 - y_l**2)**2
			F =  (1.-xi)**2 - 2*y_H**2*(1.+xi)+y_H**4

			return L * np.divide(
							0.33333333333 * m_H**2 * y_H**2 * xi * gs**2 * F * (2*xi**2 - G_plus)
		 					+ 0.04166666666 * f**2 *  np.divide(3*F*(xi**2 - (y_l**2-y_N**2)**2) - L**2 + 12 * y_H**2 * xi*(2*xi**2 - G_plus), m_H**2)
							+ 0.04166666666 * m_H**2 * a_plus**2 * F * (F*(2*xi**2 - G_plus) + 3*G_minus*(1.-y_H**2)**2)
							+ 0.125 * m_H**2 * xi**2 * a_minus**2 * F * G_minus
							+ 0.08333333333 * f * a_plus * (3 * xi * F * G_minus + (1.-xi-y_H**2) * (3*F * (xi**2-(y_l**2-y_N**2)**2) - L**2 )  )
							+ 0.25 * xi * f * a_minus * F * G_minus
							+ 0.25 * m_H**2 * xi * a_plus * a_minus * F * G_minus * (1.-y_H**2),
						xi**3)
		return 0.
	@staticmethod
	def diff_weight_tau_to_NLaNua(m_N, m_l, MNl_2, M_lnu_2):
		if m_l + m_N < c.m_tau:
			y_l = m_l / c.m_tau
			y_N = m_N / c.m_tau
			xi = M_lnu_2 / c.m_tau**2 
			return np.divide(
					(xi-y_l**2)**2 * np.sqrt(lambda_Kallen(1, np.sqrt(xi), y_N)) * ( (xi+2*y_l**2)*(1-y_N**2)**2 + xi*(xi-y_l**2)*(1+y_N**2-y_l**2) - xi*y_l**4 - 2*xi**3), xi**3
			)
		return 0.

class ds_production_BRs:
	"""Branching ratios of pseudoscalar meson decays into dark scalars taken from hep-ph/0406232""" 
	@staticmethod
	def get_BRs(parents, daughters, ds_mass):
		branchings = []
		for parent, daughter in zip(parents, daughters):
			branchings.append(hnl_production_BRs.branching_ratios(parent, daughter, ds_mass))
		return np.array(branchings)
	@staticmethod
	def branching_ratios(parent, daughter, ds_mass, coupling = s.coupling_ref, pair_prod = False):
		if pair_prod:
			if parent == "Bmeson":
				if 	 daughter == "K":      return ds_production_BRs.B_K_2S(ds_mass, lambda_s=coupling)
				elif daughter == "Kstar" : return ds_production_BRs.B_Kstar_2S(ds_mass, lambda_s=coupling)
			elif parent == "B0meson":
				if 	 daughter == "K0" :    return ds_production_BRs.B0_K0_2S(ds_mass, lambda_s=coupling)
				elif daughter == "K0star": return ds_production_BRs.B0_K0star_2S(ds_mass, lambda_s=coupling)
		if parent == "Bmeson":
			if 	 daughter == "K":      return ds_production_BRs.B_K_S(ds_mass, s_th=coupling)
			elif daughter == "Kstar" : return ds_production_BRs.B_Kstar_S(ds_mass, s_th=coupling)
		elif parent == "B0meson":
			if 	 daughter == "K0" :    return ds_production_BRs.B0_K0_S(ds_mass, s_th=coupling)
			elif daughter == "K0star": return ds_production_BRs.B0_K0star_S(ds_mass, s_th=coupling)
		return 0.

	@staticmethod
	def B_K_S(m_S, s_th = s.coupling_ref): #B -> K + DS
		f_0 = 0.330/(1-m_S**2/37.46) #hep-ph/0406232
		if c.m_B > m_S + c.m_K: return np.abs(c.C_bs_S*s_th)**2 /(16*np.pi*c.m_B**3) * (c.m_B**2 - c.m_K**2)**2/(c.m_q[4]-c.m_q[2])**2 * np.sqrt(lambda_Kallen(c.m_B,c.m_K,m_S)) * f_0**2 / c.Gamma_B
		return 0.

	def B0_K0_S(m_S, s_th = s.coupling_ref): #B0 -> K0 + DS
		f_0 = 0.330/(1-m_S**2/37.46) #hep-ph/0406232
		if c.m_B0 > m_S + c.m_K0: return np.abs(c.C_bs_S*s_th)**2  /(16*np.pi*c.m_B0**3) * (c.m_B0**2 - c.m_K0**2)**2/(c.m_q[4]-c.m_q[2])**2 * np.sqrt(lambda_Kallen(c.m_B0,c.m_K0,m_S)) * f_0**2 / c.Gamma_B0
		return 0.
	@staticmethod
	def B_Kstar_S(m_S, s_th = s.coupling_ref): #B -> K* + DS
		A_0 = 1.364/(1-m_S**2/(c.m_B**2)) - 0.990/(1-m_S**2/36.78) #hep-ph/0412079
		if c.m_B > m_S + c.m_Kstar: return np.abs(c.C_bs_S*s_th)**2  /(16*np.pi*c.m_B**3)/(c.m_q[4]+c.m_q[2])**2 * np.sqrt(lambda_Kallen(c.m_B,c.m_Kstar,m_S)**3) * A_0**2 / c.Gamma_B
		return 0.
	@staticmethod
	def B0_K0star_S(m_S, s_th = s.coupling_ref): #B0 -> K0* + DS
		A_0 = 1.364/(1-m_S**2/(c.m_B**2)) - 0.990/(1-m_S**2/36.78) #hep-ph/0412079
		if c.m_B0 > m_S + c.m_K0star: return np.abs(c.C_bs_S*s_th)**2 /(16*np.pi*c.m_B0**3)/(c.m_q[4]+c.m_q[2])**2 * np.sqrt(lambda_Kallen(c.m_B0,c.m_K0star,m_S)**3) * A_0**2 / c.Gamma_B0
		return 0.
	@staticmethod
	def B_K_2S(m_S,lambda_s = s.coupling_ref): #B -> K + 2S
		return integrate.quad(lambda m12_2: ds_production_BRs.d_BR_B_K_2S(m_S, m12_2, 0, lambda_s),4*m_S**2,(c.m_B-c.m_K)**2)[0]
	@staticmethod
	def B0_K0_2S(m_S,lambda_s = s.coupling_ref): #B0 -> K0 + 2S
		return integrate.quad(lambda m12_2: ds_production_BRs.d_BR_B0_K0_2S(m_S, m12_2, 0, lambda_s),4*m_S**2,(c.m_B0-c.m_K0)**2)[0]
	@staticmethod
	def B_Kstar_2S(m_S,lambda_s = s.coupling_ref): #B -> K + 2S
		return integrate.quad(lambda m12_2: ds_production_BRs.d_BR_B_Kstar_2S(m_S, m12_2, 0, lambda_s),4*m_S**2,(c.m_B-c.m_Kstar)**2)[0]
	@staticmethod
	def B0_K0star_2S(m_S,lambda_s = s.coupling_ref): #B -> K + 2S
		return integrate.quad(lambda m12_2: ds_production_BRs.d_BR_B0_K0star_2S(m_S, m12_2, 0, lambda_s),4*m_S**2,(c.m_B0-c.m_K0star)**2)[0]
	@staticmethod
	def d_BR_B_K_2S(m_S, m12_2, m23_2, lambda_s = s.coupling_ref): #B -> K + 2S
		f_0 = 0.330/(1-m12_2/37.46) #hep-ph/0406232
		if c.m_B > 2*m_S + c.m_K: return np.abs(c.C_bs_2S*lambda_s)**2 /(128*np.pi**3*c.m_B**3)*(c.m_B**2 - c.m_K**2)**2/(c.m_q[4]-c.m_q[2])**2 * np.sqrt(ds_production_BRs.I_3b(m_S,c.m_B,c.m_K,m12_2)) * f_0**2 / c.Gamma_B
		return 0.
	@staticmethod
	def d_BR_B0_K0_2S(m_S, m12_2, m23_2, lambda_s = s.coupling_ref): #B0 -> K0 + 2S
		f_0 = 0.330/(1-m12_2/37.46) #hep-ph/0406232
		if c.m_B0 > 2*m_S + c.m_K0: return np.abs(c.C_bs_2S*lambda_s)**2 /(128*np.pi**3*c.m_B0**3)*(c.m_B0**2 - c.m_K0**2)**2/(c.m_q[4]-c.m_q[2])**2 * np.sqrt(ds_production_BRs.I_3b(m_S,c.m_B0,c.m_K0,m12_2)) * f_0**2 / c.Gamma_B0
		return 0.
	@staticmethod
	def d_BR_B_Kstar_2S(m_S, m12_2, m23_2, lambda_s = s.coupling_ref): #B -> K + 2S
		A_0 = 1.364/(1-m12_2/(c.m_B**2)) - 0.990/(1-m12_2/36.78) #hep-ph/0412079
		if c.m_B > 2*m_S + c.m_Kstar: return np.abs(c.C_bs_2S*lambda_s)**2 /(128*np.pi**3*c.m_B**3)*lambda_Kallen(c.m_B,c.m_Kstar,np.sqrt(m12_2))/(c.m_q[4]+c.m_q[2])**2 * np.sqrt(ds_production_BRs.I_3b(m_S,c.m_B,c.m_Kstar,m12_2)) * A_0**2 / c.Gamma_B
		return 0.
	@staticmethod
	def d_BR_B0_K0star_2S(m_S, m12_2, m23_2, lambda_s = s.coupling_ref): #B0 -> K0 + 2S
		A_0 = 1.364/(1-m12_2/(c.m_B0**2)) - 0.990/(1-m12_2/36.78) #hep-ph/0412079
		if c.m_B0 > 2*m_S + c.m_K0star: return np.abs(c.C_bs_2S*lambda_s)**2 /(128*np.pi**3*c.m_B0**3)*lambda_Kallen(c.m_B0,c.m_K0star,np.sqrt(m12_2))/(c.m_q[4]+c.m_q[2])**2 * np.sqrt(ds_production_BRs.I_3b(m_S,c.m_B0,c.m_K0star,m12_2)) * A_0**2 / c.Gamma_B0
		return 0.
	@staticmethod
	def I_3b(m_S,m_B,m_K,m12_2):
		return (m12_2**2-2*m12_2*(m_B**2+m_K**2)+(m_B**2-m_K**2)**2)*(1-4*m_S**2/m12_2)

def get_branching_ratios(exo, exo_mass, parent, daughter_sm, coupling=s.coupling_ref):
	"""utility class to get branching ratios of decays of pseudoscalar mesons parent into exotic particles exo.

	Args:
		exo (str): exotic decay daughter
		exo_mass (float): exotic mass
		parent (str): pseudoscalar parent
		daughter_sm (str / list[str]): SM decay daughters
		coupling (float, optional): coupling suppression strenght at which to evaluate branching ratio. Defaults to s.coupling_ref.

	Returns:
		float: branching ratio
	"""
	if exo == "alp": 
		return alp_production_BRs.branching_ratios(parent, daughter_sm, exo_mass, g=coupling)
	elif exo == "hnl": 
		return hnl_production_BRs.branching_ratios(parent, daughter_sm, exo_mass, U2=coupling)
	elif exo == "ds":
		return ds_production_BRs.branching_ratios(parent, daughter_sm, exo_mass, coupling=coupling)
	elif exo == "2ds":
		return ds_production_BRs.branching_ratios(parent, daughter_sm, exo_mass, coupling=coupling, pair_prod=True)
	return 0.
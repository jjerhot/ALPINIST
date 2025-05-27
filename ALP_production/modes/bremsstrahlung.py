# modes/bremsstrahlung.py
# dark scalar/vector bremsstrahlung production following 2108.05900

import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
from ALP_production.general import exotic_constants as c
from ALP_production.general import exotic_functions as f
from ALP_production.general import exotic_production_setup as setup
from ALP_production.modes import dp_brems_functions as bf

class brems_production:
	def __init__(self, experiment,daugther_name, Lambda = 1.5, z_min = 0.1, z_max = 0.9):
		self.exp = experiment
		self.daughter_exo = daugther_name
		self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced in Dark Bremsstrahlung | generated with quasi-real approximation following 2108.05900"

		self.Lambda = Lambda # GeV; cutoff in [1,2] GeV
		self.z_min = z_min
		self.z_max = z_max

		self.th_list = np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins).tolist()
		en_halfstep = np.floor(9.5*setup.p_beam.get(self.exp,400)/setup.energy_bins/2)/10
		self.en_list = [round(en_halfstep + 2*i*en_halfstep,4) for i in range(setup.energy_bins)]
		self.m_lists = []
		for iFile in range(len(setup.mass_bins)):
			self.m_lists.append(np.linspace(setup.mass_min[iFile], setup.mass_max[iFile], num=setup.mass_bins[iFile]).tolist())

	def process_pool(self,nthreads,single_mass_point=0):
		"""Wrapper to run production in parallel threads.
		Args:
			nthreads (int): number of parallel threads to run
			single_mass_point (float, optional): Single mass point for which to run. Defaults to 0.
		"""

		export_units_header = "Theta_X[rad] E_X[GeV] m_X[GeV] dY[per(rad GeV N_pN coupling^2)]"
		if self.daughter_exo == "dp":
			export_units_header.replace('coupling','eps')
			brems_process = self.dv_brems_process
		elif self.daughter_exo == "ds":
			export_units_header.replace('coupling','theta')
			brems_process = self.ds_brems_process

		if single_mass_point>0:
			list_brems = brems_process(single_mass_point)
			f.export(setup.experiments[self.exp],self.daughter + "_Brems_beam" + str(setup.p_beam[self.exp]) + "GeV_" + str(int(single_mass_point*1e6)) + "keV_" + setup.experiments[self.exp] + ".dat",np.reshape(list_brems,(len(self.th_list)*len(self.en_list),4)),header = self.export_header+'\n'+export_units_header)
			return

		for iFile in range(len(setup.mass_bins)):
			iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
			if setup.mass_min[iFile]*1000 < 1:
				iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"

			print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")
			with Pool(processes=nthreads) as pool:
				list_brems = list(tqdm(pool.imap(brems_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))
			pool.close()
			pool.join() 

			# export
			f.export(setup.experiments[self.exp],self.daughter_exo + "_Brems_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_brems,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)), header = self.export_header+'\n'+export_units_header)
		return

	def form_factor_scalar(self,m_X):
		return np.abs(sum([f*m**2/(m**2-m_X**2-1j*m_X*Gamma) for (f,m,Gamma) in zip(c.f_f0,c.m_f0,c.Gamma_f0)]))

	def form_factor_vector(self,m_X):
		f_isovector = np.abs(sum([f*m**2/(m**2-m_X**2-1j*m_X*Gamma) for (f,m,Gamma) in zip(c.f_rho0,c.m_rho0,c.Gamma_rho0)]))
		f_isoscalar = np.abs(sum([f*m**2/(m**2-m_X**2-1j*m_X*Gamma) for (f,m,Gamma) in zip(c.f_omega0,c.m_omega0,c.Gamma_omega0)]))
		return f_isovector+f_isoscalar

	def __form_factor_ppstar(self,pp2): # pp2 = (p-p_s)**2
		return self.Lambda**4/(self.Lambda**4+(pp2-c.m_p**2)**2)

	def __H(self,m_ds,z,pt2):
		return pt2 + z**2 * c.m_p**2 + (1-z)*m_ds**2

	def __wS(self,m_X,z,pt2): # mixing = 1
		gSNN = 1.2E-3
		return gSNN**2 * 1/(8*np.pi**2) * (self.form_factor_scalar(m_X) * self.__form_factor_ppstar(c.m_p**2 - self.__H(m_X,z,pt2)/z))**2 * z * (1 + (1-z)*(4*c.m_p**2 - m_X**2)/self.__H(m_X,z,pt2))/(2 * self.__H(m_X,z,pt2))

	def __wV(self,m_X,z,pt2): # mixing = 1
		return c.alpha_EM * 1/(2*np.pi) * (self.form_factor_vector(m_X) * self.__form_factor_ppstar(c.m_p**2 - self.__H(m_X,z,pt2)/z))**2 * z * (1 - (1-z)*(2*c.m_p**2 + m_X**2)/self.__H(m_X,z,pt2) + self.__H(m_X,z,pt2)/(2*z**2*m_X**2))/(self.__H(m_X,z,pt2))

	def ds_brems_process(self,m_X):
		bcsTot = bf.brems_cross_section()
		data_list = []
		p_ds = 0
		for e_X in self.en_list:
			if e_X>m_X: p_ds = np.sqrt(e_X**2-m_X**2)
			z_X = p_ds/setup.p_beam[self.exp]
			sPrime = 2*c.m_p*setup.p_beam[self.exp]*(1-z_X)+2*c.m_p**2
			s = 2*c.m_p*np.sqrt(setup.p_beam[self.exp]**2+c.m_p**2)
			sigmaPPprime = bcsTot.SigmaPP(sPrime)
			sigmaPP = bcsTot.SigmaPP(s)
			data_list_sublist = []
			if e_X > m_X and z_X > self.z_min and z_X < self.z_max:
				for th_X in self.th_list:
					pt2 = p_ds**2 * np.sin(th_X)**2
					jacobian = e_X*z_X*np.sin(2*th_X)
					dsigma = sigmaPPprime * jacobian * self.__wS(m_X,z_X,pt2)
					data_list_sublist.append([th_X, e_X, m_X, dsigma/sigmaPP])
			else:
				for th_X in self.th_list:
					data_list_sublist.append([th_X, e_X, m_X, 0.])

			data_list.append([data_list_sublist])

		return data_list

	def dv_brems_process(self,m_X):
		bcsTot = bf.brems_cross_section()
		data_list = []
		p_ds = 0
		for e_X in self.en_list:
			if e_X>m_X: p_ds = np.sqrt(e_X**2-m_X**2)
			z_X = p_ds/setup.p_beam[self.exp]
			sPrime = 2*c.m_p*setup.p_beam[self.exp]*(1-z_X)+2*c.m_p**2
			s = 2*c.m_p*np.sqrt(setup.p_beam[self.exp]**2+c.m_p**2)
			sigmaPPprime = bcsTot.SigmaPP(sPrime)
			sigmaPP = bcsTot.SigmaPP(s)
			data_list_sublist = []
			if e_X > m_X and z_X > self.z_min and z_X < self.z_max:
				for th_X in self.th_list:
					pt2 = p_ds**2 * np.sin(th_X)**2
					jacobian = e_X*z_X*np.sin(2*th_X)
					dsigma = sigmaPPprime * jacobian * self.__wV(m_X,z_X,pt2)
					data_list_sublist.append([th_X, e_X, m_X, dsigma/sigmaPP])
			else:
				for th_X in self.th_list:
					data_list_sublist.append([th_X, e_X, m_X, 0.])

			data_list.append([data_list_sublist])

		return data_list
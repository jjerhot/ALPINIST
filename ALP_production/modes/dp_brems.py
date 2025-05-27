# modes/dp_brems.py
# dark photon bremsstrahlung production according to 1311.3870

import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
from ALP_production.general import exotic_constants as c
from ALP_production.general import exotic_functions as f
from ALP_production.general import exotic_production_setup as setup
from ALP_production.modes import dp_brems_functions as bf

class brems_production:
	def __init__(self, experiment,daugther_name):
		self.exp = experiment
		self.daughter = daugther_name
		self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced in Dark Bremsstrahlung | generated with modified Weizsacker-Williams-approximation following 1311.3870"

		self.th_list = np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins).tolist()
		en_halfstep = np.floor(9.5*setup.p_beam.get(self.exp,400)/setup.energy_bins/2)/10
		self.en_list = [round(en_halfstep + 2*i*en_halfstep,4) for i in range(setup.energy_bins)]
		self.m_lists = []
		for iFile in range(len(setup.mass_bins)):
			self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())

	def process_pool(self,nthreads,single_mass_point=0):
		"""Wrapper to run production in parallel threads.
		Args:
			nthreads (int): number of parallel threads to run
			single_mass_point (float, optional): Single mass point for which to run. Defaults to 0.
		"""

		export_units_header = "Theta_X[rad] E_X[GeV] m_X[GeV] dY[per(rad GeV N_pN eps^2)]"

		if single_mass_point>0:
			list_brems = self.dp_brems_process(single_mass_point)
			f.export(setup.experiments[self.exp],self.daughter + "_Brems_beam" + str(setup.p_beam[self.exp]) + "GeV_" + str(int(single_mass_point*1e6)) + "keV_" + setup.experiments[self.exp] + ".dat",np.reshape(list_brems,(len(self.th_list)*len(self.en_list),4)),header = self.export_header+'\n'+export_units_header)
			return
		for iFile in range(len(setup.mass_bins)):
			iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
			if setup.mass_min[iFile]*1000 < 1:
				iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
			print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")

			with Pool(processes=nthreads) as pool:
				list_brems = list(tqdm(pool.imap(self.dp_brems_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

			pool.close()
			pool.join() 

			# export
			f.export(setup.experiments[self.exp],self.daughter + "_Brems_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_brems,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)),header = self.export_header+'\n'+export_units_header)
		return

	def dp_brems_process(self,m_dp):
		"""evaluates modified WWA for dark photon with mass m_dp emitted at angle theta to proton axis with energy e

		Args:
			m_dp (float): dark photon mss

		Returns:
			array(float): [theta, e, m_dp, expected yield]
		"""
		bcsTot = bf.brems_cross_section()
		bcsTot.setM(m_dp)
		bcsTot.setMomentum(setup.p_beam[self.exp])
		bcsTot.setEpsilon(1.) # asserting that refrence epsilon = 1
		data_list = []
		for e_dp in self.en_list:
			z_dp = e_dp/setup.p_beam[self.exp]
			data_list_sublist = []
			#skip first energy bin to satisfy WW approx
			if e_dp > m_dp and e_dp != self.en_list[0]:
				for th_dp in self.th_list:
					pt = e_dp*np.tan(th_dp)
					data_list_sublist.append([th_dp, e_dp, m_dp, bcsTot.evaluate2d([z_dp,pt],m_dp)*e_dp/np.sqrt(th_dp**2+1)])
			else:
				for th_dp in self.th_list:
					data_list_sublist.append([th_dp, e_dp, m_dp, 0.])

			data_list.append([data_list_sublist])

		return data_list
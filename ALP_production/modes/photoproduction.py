import numpy as np 

from tqdm import tqdm
from multiprocessing import Pool
from scipy.interpolate import LinearNDInterpolator
from scipy.integrate import dblquad

from ALP_production.general import exotic_constants as c 
from ALP_production.general import exotic_functions as f
from ALP_production.general import exotic_production_setup as setup

class primakoff_production:
    def __init__(self, experiment,daugther_name, use_trapz = False, use_vegas=False): #use_trapz will result in much faster simulation at the cost of accuracy
        self.use_trapz=use_trapz
        self.use_vegas=use_vegas
        self.exp = experiment
        self.p_beam  = setup.p_beam.get(self.exp,400)
        self.sigma_pp = c.sigma_pp[self.p_beam]
        self.daughter = daugther_name
        self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced by Primakoff effect | generated assuming Budnev photon dist. and ALPtraum photon nucleus interactions"

        self.th_list = np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins).tolist()
        self.en_list = np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)), num=setup.energy_bins).tolist()
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())

        method = "VegasFlow" if use_vegas else ("Trapezoid" if use_trapz else "QUADPACK")
        print(f"[Info:] \t Selected integration method: {method}")


    def process_pool(self,nthreads,single_mass_point=0):
        """Wrapper to run production in parallel threads.
        Args:
            nthreads (int): number of parallel threads to run
            single_mass_point (float, optional): Single mass point for which to run. Defaults to 0.
        """

        export_units_header = "Theta[rad], E_x [GeV], m_x [GeV], dY[per(rad GeV N_PoT g_ga**2)]"

        if single_mass_point>0:
            if self.use_vegas: 
                import tensorflow as tf
                tf.config.threading.set_inter_op_parallelism_threads(int(nthreads))
                list_primakoff = self.primakoff_process(single_mass_point)
            else:
                with Pool(processes=nthreads) as pool:
                    list_primakoff = list(tqdm(pool.imap(self.primakoff_process_single_mass, zip(np.full_like(self.en_list, single_mass_point), self.en_list)),total=len(self.en_list)))
                pool.close()
                pool.join() 
            f.export(setup.experiments[self.exp],self.daughter + "_Primakoff_beam" + str(setup.p_beam[self.exp]) + "GeV_" + str(int(single_mass_point*1e6)) + "keV_" + setup.experiments[self.exp] + ".dat",np.reshape(list_primakoff,(len(self.th_list)*len(self.en_list),4)), header = self.export_header+'\n'+export_units_header)
            return

        for iFile in range(len(setup.mass_bins)):
            iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")
            if self.use_vegas: 
                import tensorflow as tf
                tf.config.threading.set_inter_op_parallelism_threads(int(nthreads))
                list_primakoff = [self.primakoff_process(m, nthreads) for m in tqdm(self.m_lists[iFile],total=len(self.m_lists[iFile]) )]
            else:
                with Pool(processes=nthreads) as pool:
                    list_primakoff = list(tqdm(pool.imap(self.primakoff_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

                pool.close()
                pool.join() 

            # export
            f.export(setup.experiments[self.exp],self.daughter + "_Primakoff_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_primakoff,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)),header = self.export_header+'\n'+export_units_header)
        
        return

    def primakoff_process_single_mass(self, x_info):
        m_x, e_x = x_info
        data_list = []
        if e_x > m_x:
            for th_x in self.th_list:
                data_list.append([th_x, e_x, m_x, Photoproduction_functions.dsigma_pN_alt(m_x, e_x, th_x, g_a = 1,  exp_name=self.exp, use_trapz=self.use_trapz, use_vegas=self.use_vegas)/(self.sigma_pp*c.A_target[self.exp]**0.77)])
        else:
            for th_x in self.th_list:
                data_list.append([th_x, e_x, m_x, 0.])
        return data_list
    def primakoff_process(self, m_x):
        """evaluates primakoff process producing an alp with mass m_x emitted at angle theta to beam axis with energy e
        Args:
            m_x (float): exotic mss

        Returns:
            array(float): [theta, e, m_x, expected yield]
        """
        data_list = []
        for e_x in tqdm(self.en_list):
            data_list_sublist = []
            if e_x > m_x:
                for th_x in self.th_list:
                    data_list_sublist.append([th_x, e_x, m_x, Photoproduction_functions.dsigma_pN_alt(m_x, e_x, th_x, g_a = 1,  exp_name=self.exp, use_trapz=self.use_trapz, use_vegas=self.use_vegas)/(self.sigma_pp*c.A_target[self.exp]**0.77)])
            else:
                for th_x in self.th_list:
                    data_list_sublist.append([th_x, e_x, m_x, 0.])

            data_list.append([data_list_sublist])

        return data_list
	
class photon_from_meson_production:
    def __init__(self, experiment,daugther_name, use_trapz=False, use_vegas=False): #use_trapz will result in much faster simulation at the cost of accuracy
        self.use_trapz=use_trapz
        self.use_vegas=use_vegas
        self.exp = experiment
        self.daughter = daugther_name
        self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced from mesondecay photons | generated using photons from Pythia 8.2s SoftQCD with PDF:pSet=2 and ALPtraum photon nucleus interactions"
        self.p_beam  = setup.p_beam.get(self.exp,400)
        self.th_list = np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins).tolist()
        self.en_list = np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)), num=setup.energy_bins).tolist()
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())

        method = "VegasFlow" if use_vegas else ("Trapezoid" if use_trapz else "QUADPACK")
        print(f"[Info:] \t Selected integration method: {method}")
        self.meson_photon_kde = self._init_photon_distribution(experiment)
    @staticmethod
    def _init_photon_distribution(exp_name):
        if exp_name == "KOTO" or exp_name == "KOTO2":
            log10_p, log10_pt, density = np.loadtxt('/'.join(__file__.split('/')[:-1])+f'/../tab_gammas/{setup.p_beam[exp_name]}GeV/EvsPTLoggamma_200kEvts_{setup.p_beam[exp_name]}GeV_full.dat').T	
        else:
            log10_p, log10_pt, density = np.loadtxt('/'.join(__file__.split('/')[:-1])+f'/../tab_gammas/{setup.p_beam[exp_name]}GeV/EvsPTLoggamma_100kEvts_{setup.p_beam[exp_name]}GeV_full.dat').T		
        xy = np.c_[log10_p, log10_pt]
        return LinearNDInterpolator(xy, density, fill_value=0)
    
    def gamma_density_fun(self, x, pt2):
        """Gamma density as produced from the decay of light mesons produced in soft QCD processes.
        Args:
            x (float): energy share of proton momentum
            pt2 (float): squared momentum transverse to beam axis

        Returns:
            float or array(float): density estimate 
        """
        e_beam = np.sqrt(self.p_beam**2 + c.m_p**2)
        p, pt = x*e_beam, np.sqrt(pt2)
        if hasattr(pt, '__iter__'): return e_beam/(2.*p*pt**2*np.log(10)**2)* self.meson_photon_kde(np.full_like(pt2,np.log10(p)), np.log10(pt))
        return e_beam/(2.*p*pt**2*np.log(10)**2)* self.meson_photon_kde(np.log10(p),np.log10(pt))

    def process_pool(self,nthreads,single_mass_point=0):
        """Wrapper to run production in parallel threads.
        Args:
            nthreads (int): number of parallel threads to run
            single_mass_point (float, optional): Single mass point for which to run. Defaults to 0.
        """

        export_units_header = "Theta[rad], E_x [GeV], m_x [GeV], dY[per(rad GeV N_PoT g_ga**2)]"

        if single_mass_point>0:
            if self.use_vegas: 
                import tensorflow as tf
                tf.config.threading.set_inter_op_parallelism_threads(int(nthreads))
                list_photonfrommeson = self.photon_from_meson_process(single_mass_point)
            else:
                with Pool(processes=nthreads) as pool:
                    list_photonfrommeson = list(tqdm(pool.imap(self.photon_from_meson_process_single_mass, zip(np.full_like(self.en_list, single_mass_point), self.en_list)),total=len(self.en_list)))
                pool.close()
                pool.join() 
            f.export(setup.experiments[self.exp],self.daughter + "_PhotonFromMeson_beam" + str(setup.p_beam[self.exp]) + "GeV_" + str(int(single_mass_point*1e6)) + "keV_" + setup.experiments[self.exp] + ".dat",np.reshape(list_photonfrommeson,(len(self.th_list)*len(self.en_list),4)), header = self.export_header+'\n'+export_units_header)
            return
        
        for iFile in range(len(setup.mass_bins)):
            iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")

            if self.use_vegas: 
                import tensorflow as tf
                tf.config.threading.set_inter_op_parallelism_threads(int(nthreads))
                list_photonfrommeson = [self.photon_from_meson_process(m, nthreads) for m in tqdm(self.m_lists[iFile],total=len(self.m_lists[iFile]) )]
            else:
                with Pool(processes=nthreads) as pool:
                    list_photonfrommeson = list(tqdm(pool.imap(self.photon_from_meson_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

                pool.close()
                pool.join() 

            # export
            f.export(setup.experiments[self.exp],self.daughter + "_PhotonFromMeson_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_photonfrommeson,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)),header = self.export_header+'\n'+export_units_header)
        return

    def photon_from_meson_process_single_mass(self,x_info):
        m_x, e_x = x_info
        data_list = []
        if e_x > m_x:
            for th_x in self.th_list:
                data_list.append([th_x, e_x, m_x, Photoproduction_functions.dsigma_pN_alt(m_x, e_x, th_x, photon_density=self.gamma_density_fun, g_a = 1, exp_name=self.exp, use_trapz= self.use_trapz, use_vegas=self.use_vegas)/c.sigma_pgamma[self.exp]])
        else:
            for th_x in self.th_list:
                data_list.append([th_x, e_x, m_x, 0.])
        return data_list
    def photon_from_meson_process(self,m_x):
        """evaluates photon scattering process producing an alp with mass m_x emitted at angle theta to beam axis with energy e
        Args:
            m_x (float): exotic mss

        Returns:
            array(float): [theta, e, m_x, expected yield]
        """
        data_list = []
        for e_x in self.en_list:
            data_list_sublist = []
            if e_x > m_x:
                for th_x in self.th_list:
                    data_list_sublist.append([th_x, e_x, m_x, Photoproduction_functions.dsigma_pN_alt(m_x, e_x, th_x, photon_density=self.gamma_density_fun, g_a = 1, exp_name=self.exp, use_trapz= self.use_trapz, use_vegas=self.use_vegas)/c.sigma_pgamma[self.exp]])
            else:
                for th_x in self.th_list:
                    data_list_sublist.append([th_x, e_x, m_x, 0.])

            data_list.append([data_list_sublist])

        return data_list
	
class Photoproduction_functions:
    """Utility function to evaluate photo production of axions based on JHEP02(2016)018
    """
    def __init__(self) -> None:
        pass  
    q2_0 = 0.71 # scaling for momentum [GeV^2]
    mu2_p = 7.78 # proton mag. moment
    s = 0.9 # Helm form factor constant [fm]
    hbarCFm = 0.1973  # reduced planck constant in [GeV.fm] 
    q2min = 1e-10
    @staticmethod
    def tt(mx, Ex, pt, theta, phi): 
        """"Momentum transfer after eq 3.17"""
        return -mx**4/(4.*Ex**2) - pt**2 + 2.*Ex*pt*theta*np.cos(phi) - Ex**2*theta**2
    @staticmethod
    def R1(A_target): 
        return np.sqrt((1.23*A_target**(1./3.) - 0.6)**2 + 7./3.*(0.52*np.pi)**2 - 5.*Photoproduction_functions.s**2)/Photoproduction_functions.hbarCFm
    @staticmethod
    def custom_3spherical_bessel_j1_x(x): 
        """Custom implementation of spherical bessel j1. 
        approximation error is < 10^{-8}/280."""
        return np.where(x<0.01, 1. - 0.1*np.power(x,2), 3*np.divide(np.divide(np.sin(x),x) - np.cos(x), np.power(x,2)))
    @staticmethod
    def FF(q2, R1):
        """Helm Form factor"""
        return  Photoproduction_functions.custom_3spherical_bessel_j1_x(np.sqrt(q2)*R1)*np.exp(-q2*(Photoproduction_functions.s / Photoproduction_functions.hbarCFm)**2/2.)
    @staticmethod
    def gamma_p_Budnev(x, qt2):
        """photon from proton distribution - Budnev """

        Etg2 = qt2 + x**2*c.m_p**2
        qi2 = Etg2/(1-x)
        G2 = 1 / (1+qi2/Photoproduction_functions.q2_0)**4
        C = Photoproduction_functions.mu2_p * G2
        D = np.divide(4*c.m_p**2*G2 + qi2*Photoproduction_functions.mu2_p*G2, 4*c.m_p**2+qi2)

        if np.shape(Etg2): return np.where((Etg2 + x) > 1, 0, c.alpha_EM/np.pi * np.divide(qt2 * (1-x) * D / Etg2 + x**2 * C / 2, x * Etg2))
        if (Etg2 + x) > 1.: return 0.
        return np.divide( qt2 * (1-x) * D / Etg2 + x**2 * C / 2, x * np.pi* Etg2)  * c.alpha_EM

    @staticmethod
    def dsigma_gammaN_alt(mx2_4gax2, Ex, pt, theta, phi, R1, Z_target):
        """eq 3.16 in JHEP02(2016)018"""
        pE_x = (pt - Ex*theta)**2 + 2.*Ex*pt*theta*(1. - np.cos(phi))
        return np.divide(c.alpha_EM * Ex**2 * pE_x, 4*(pE_x+mx2_4gax2)**2) * np.power( Z_target * Photoproduction_functions.FF(pE_x+mx2_4gax2, R1 = R1), 2) 

    @staticmethod
    def dsigma_pN_alt(mx, Ex, theta, photon_density = None, g_a = 1., exp_name="NA62", use_trapz = False, use_vegas=False, eps_rel = 1e-4):
        """Integrate photon photon differential cross sections to evaluate pp level
        Args:
            mx (float): alp mass
            Ex (float): alp energy
            theta (float): alp angle with beam axis
            photon_density (function(x, qt2), optional): _description_. Defaults to Budnev distribution.
            g_a (float, optional): axion-photon-coupling. Defaults to 1.
            exp_name (str, optional): experiment facility. Defaults to "NA62".
            use_trapz (bool, optional): use simplified numerical integration (handle with great care). Defaults to False.
            use_vegas(bool, optional):  use vegas integration. Defaults to False.
            eps_rel: relative error to tolerate
        Returns:
            float: proton proton level differential scattering cross section
        """
        if photon_density == None: photon_density = Photoproduction_functions.gamma_p_Budnev
        R1_ = Photoproduction_functions.R1(c.A_target[exp_name])
        E_beam = np.sqrt(setup.p_beam[exp_name]**2 + c.m_p**2)
        q2max = np.power(4.49/R1_, 2 )
        mx2_4gax2 = np.power(mx**2/(2*Ex), 2) # common calculation blocks
        delta_qE_m2 = np.sqrt(np.clip(q2max - mx2_4gax2, 0, None))

        prefactor = (16*g_a**2*c.g_EM**4*np.sin(theta))/(np.pi*E_beam) * 3.89379*1e8 # conversion to pb 
        pt2low = np.power(np.clip(Ex*theta - delta_qE_m2, 0, None), 2)
        pt2up = np.power(Ex*theta + delta_qE_m2, 2)
        if pt2up == 0 or pt2up<=pt2low: return 0.
        if use_vegas:
            from ALP_production.modes.photoproduction_tf import compute_vegas_integral
            integral, error_est = compute_vegas_integral(theta, Ex, Ex/E_beam, mx2_4gax2, [pt2low, pt2up], R1_, c.Z_target[exp_name], photon_density)
            if np.isnan(integral): integral = 0
        else:
            phi_lower = lambda pt2: np.arccos(np.clip(( mx2_4gax2 + pt2 + (Ex*theta)**2 - Photoproduction_functions.q2min)/(2.*Ex*np.sqrt(pt2)*theta),-1,1))
            phi_upper = lambda pt2: np.arccos(np.clip(( mx2_4gax2 + pt2 + (Ex*theta)**2 - q2max                          )/(2.*Ex*np.sqrt(pt2)*theta),-1,1))
            f_integrand = lambda phi, pt2: photon_density(Ex/E_beam,pt2) * Photoproduction_functions.dsigma_gammaN_alt(mx2_4gax2, Ex, np.sqrt(pt2), theta, phi, R1_, c.Z_target[exp_name])
            if use_trapz: return prefactor * f.integrate_2d_trapz(f_integrand, pt2low, pt2up, phi_lower, phi_upper)
            integral, error_est = dblquad(f_integrand, pt2low, pt2up, phi_lower, phi_upper, epsrel=eps_rel)
        return prefactor * integral
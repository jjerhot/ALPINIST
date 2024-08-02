import numpy as np 

from tqdm import tqdm
from multiprocessing import Pool
from scipy.interpolate import LinearNDInterpolator
from scipy.special import spherical_jn
from scipy.integrate import dblquad

from ALP_production.general import exotic_constants as c 
from ALP_production.general import exotic_functions as f
from ALP_production.general import exotic_production_setup as setup


def integrate_2d_trapz(f, a, b, x_lower_func, x_upper_func, y_samples=1000, x_samples=1000):
    """
    np.trapz based implementation of scipy.integrate.dblquad
    """
    y_values = np.linspace(a, b, y_samples)
    x_lower_values = x_lower_func(y_values)
    x_upper_values = x_upper_func(y_values)
    delta_x = x_upper_values - x_lower_values
    
    y_mesh, x_mesh = np.meshgrid(y_values, np.linspace(0, 1, x_samples))
    x_mesh = x_lower_values + delta_x * x_mesh
    
    f_values = np.nan_to_num(f(x_mesh, y_mesh))
   
    integrated_over_x = np.trapz(f_values, dx = delta_x / y_samples, axis=0) # Integrate first over x using trapz
    result = np.trapz(integrated_over_x, y_values)
    
    return result    


class primakoff_production:
    def __init__(self, experiment,daugther_name, use_trapz = False): #use_trapz will result in much faster simulation at the cost of accuracy
        self.use_trapz=use_trapz
        self.exp = experiment
        self.daughter = daugther_name
        self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced by Primakoff effect | generated assuming Budnev photon dist. and ALPtraum photon nucleus interactions"

        self.th_list = np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins).tolist()
        en_halfstep = np.floor(9.5*setup.p_beam.get(self.exp,400)/setup.energy_bins/2)/10
        self.en_list = [round(en_halfstep + 2*i*en_halfstep,4) for i in range(setup.energy_bins)]
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())

        # uncomment for 400GeV original settings
        # for lower, upper, steps in zip([1e-4, 1e-2],[1e-2,3.01],[101,301]):
        #     self.m_lists.append(np.linspace(lower, upper, num=steps).tolist())
        # self.th_list = np.arange(0.00018, 0.01098, 0.00036).tolist()
        # self.en_list = np.arange(5.5, 324.5, 11).tolist()


    def process_pool(self,nthreads):
        """Wrapper to run production in parallel threads.
        Args:
            nthreads (int): number of parallel threads to run
        """
        for iFile in range(len(setup.mass_bins)):
            iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")

            with Pool(processes=nthreads) as pool:
                list_primakoff = list(tqdm(pool.imap(self.primakoff_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

            pool.close()
            pool.join() 

            # export
            export_units_header = "Theta[rad] E_x[GeV] dY[per(rad GeV PoT ga**2)]"
            f.export(setup.experiments[self.exp],self.daughter + "_primakoff_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_primakoff,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)),header = self.export_header+'\n'+export_units_header)
        return

    def primakoff_process(self,m_x):
        """evaluates primakoff process producing an alp with mass m_x emitted at angle theta to beam axis with energy e
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
                    data_list_sublist.append([th_x, e_x, m_x, Photoproduction_functions.dsigma_pN_alt(m_x, e_x, th_x, g_a = 1,  exp_name=self.exp, use_trapz=self.use_trapz)])
            else:
                for th_x in self.th_list:
                    data_list_sublist.append([th_x, e_x, m_x, 0.])

            data_list.append([data_list_sublist])

        return data_list
	
class photon_from_meson_production:
    def __init__(self, experiment,daugther_name, use_trapz=False): #use_trapz will result in much faster simulation at the cost of accuracy
        self.use_trapz=use_trapz

        self.exp = experiment
        self.daughter = daugther_name
        self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced from mesondecay photons | generated using photons from Pythia 8.2s SoftQCD with PDF:pSet=2 and ALPtraum photon nucleus interactions"
        self.p_beam  = setup.p_beam.get(self.exp,400)
        self.th_list = np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins).tolist()
        en_halfstep = np.floor(9.5*self.p_beam/setup.energy_bins/2)/10
        self.en_list = [round(en_halfstep + 2*i*en_halfstep,4) for i in range(setup.energy_bins)]
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())

        # uncomment for 400GeV original settings
        # self.th_list = np.arange(0.00018, 0.01098, 0.00036).tolist()
        # self.en_list = np.arange(5.5, 324.5, 11).tolist()
        # self.m_lists = []
        # for lower, upper, steps in zip([1e-4, 1e-2],[1e-2,3.01],[101,301]):
        #     self.m_lists.append(np.linspace(lower, upper, num=steps).tolist())

        self.meson_photon_kde = self._init_photon_distribution(experiment)
    @staticmethod
    def _init_photon_distribution(exp_name):
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

    def process_pool(self,nthreads):
        """Wrapper to run production in parallel threads.
        Args:
            nthreads (int): number of parallel threads to run
        """
        for iFile in range(len(setup.mass_bins)):
            iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")

            with Pool(processes=nthreads) as pool:
                list_photonfrommeson = list(tqdm(pool.imap(self.photon_from_meson_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

            pool.close()
            pool.join() 

            # export
            export_units_header = "Theta[rad] E_x[GeV] dY[per(rad GeV L_p ga**2)]"
            f.export(setup.experiments[self.exp],self.daughter + "_photonfrommeson_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_photonfrommeson,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)),header = self.export_header+'\n'+export_units_header)
        return

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
                    data_list_sublist.append([th_x, e_x, m_x, Photoproduction_functions.dsigma_pN_alt(m_x, e_x, th_x, photon_density=self.gamma_density_fun, g_a = 1, exp_name=self.exp, use_trapz= self.use_trapz)])
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
    q2_0 = 0.71 # scaling for momentum
    mu2_p = 7.78 # proton mag. moment
    s = 0.9 # Helm form factor constant times    
    hbarCFm = 0.1973  # reduced planck constant in fm 
    q2min = 1e-10
    @staticmethod
    def tt(mx, Ex, pt, theta, phi): 
        """"Momentum transfer after eq 3.17"""
        return -mx**4/(4.*Ex**2) - pt**2 + 2.*Ex*pt*theta*np.cos(phi) - Ex**2*theta**2
    @staticmethod
    def R1(A_material): 
        return np.sqrt((1.23*A_material**(1./3.) - 0.6)**2 + 7./3.*np.pi**2*0.52**2 -  5.*Photoproduction_functions.s**2)/Photoproduction_functions.hbarCFm; 
    @staticmethod
    def FF(q2, exp_name = "NA62"):
        """Helm Form factor"""
        q = np.sqrt(q2)
        return  3.*spherical_jn(1, q*Photoproduction_functions.R1(setup.A_target[exp_name]))/( q*Photoproduction_functions.R1(setup.A_target[exp_name]))*np.exp(-(q * Photoproduction_functions.s * Photoproduction_functions.hbarCFm)**2/2.);
    @staticmethod
    def gamma_p_Budnev(x, qt2):
        """photon from proton distribution - Budnev"""
        if np.shape(x) or np.shape(qt2): return np.where((qt2 + x**2*c.m_p**2 + x) > 1, 0, c.alpha_EM/np.pi*(1.-x)/x*((qt2*(1.-x))/(qt2+x**2*c.m_p**2)**2  *(4.*c.m_p**2*(1.-x)+Photoproduction_functions.mu2_p*(qt2+x**2*c.m_p**2))/(4.*c.m_p**2*(1.-x)+qt2+x**2*c.m_p**2)+x**2/(2.*(qt2+x**2*c.m_p**2))*Photoproduction_functions.mu2_p) / (1.+(qt2+x**2*c.m_p**2)/(Photoproduction_functions.q2_0*(1.-x)))**4)
        if (qt2 + x**2*c.m_p**2 + x) > 1.: return 0.
        return c.alpha_EM/np.pi*(1.-x)/x*((qt2*(1.-x))/(qt2+x**2*c.m_p**2)**2  *(4.*c.m_p**2*(1.-x)+Photoproduction_functions.mu2_p*(qt2+x**2*c.m_p**2))/(4.*c.m_p**2*(1.-x)+qt2+x**2*c.m_p**2)+x**2/(2.*(qt2+x**2*c.m_p**2))*Photoproduction_functions.mu2_p) / (1.+(qt2+x**2*c.m_p**2)/(Photoproduction_functions.q2_0*(1.-x)))**4
    @staticmethod
    def dsigma_gammaN_alt(mx, Ex, pt, theta, phi, exp_name="NA62"):
        """eq 3.16 in JHEP02(2016)018"""
        return np.divide((c.alpha_EM*Ex**2*((pt - Ex*theta)**2 + 2.*Ex*pt*theta*(1. - np.cos(phi))))*setup.Z_target[exp_name]**2*Photoproduction_functions.FF((mx**2/(2*Ex))**2 + (pt - Ex*theta)**2 +  2.*Ex*pt*theta*(1. - np.cos(phi)),exp_name=exp_name)**2,
                (4.*(mx**4/(4.*Ex**2) + (pt - Ex*theta)**2 + 2.*Ex*pt*theta*(1. - np.cos(phi)))**2))
    @staticmethod
    def dsigma_pN_alt(mx, Ex, theta, photon_density = None, g_a = 1., exp_name="NA62", use_trapz = False):
        """Integrate photon photon differential cross sections to evaluate pp level

        Args:
            mx (float): alp mass
            Ex (float): alp energy
            theta (float): alp angle with beam axis
            photon_density (function(x, qt2), optional): _description_. Defaults to Budnev distribution.
            g_a (float, optional): axion-photon-coupling. Defaults to 1.
            exp_name (str, optional): experiment facility. Defaults to "NA62".
            use_trapz (bool, optional): use simplified numerical integration (handle with great care). Defaults to False.

        Returns:
            float: proton proton level differential scattering cross section
        """
        if photon_density == None: photon_density = Photoproduction_functions.gamma_p_Budnev
        q2max = (4.49/Photoproduction_functions.R1(setup.A_target[exp_name]))**2
        E_beam = np.sqrt(setup.p_beam[exp_name]**2 + c.m_p**2)
        philow = lambda pt2: np.arccos(np.clip(((mx**2/(2*Ex))**2 + pt2 + (Ex*theta)**2 - Photoproduction_functions.q2min)/(2.*Ex*np.sqrt(pt2)*theta),-1,1))
        phiup = lambda pt2: np.arccos(np.clip(((mx**2/(2*Ex))**2 + pt2 + (Ex*theta)**2 - q2max)/(2.*Ex*np.sqrt(pt2)*theta),-1,1))
        delta_qE_m2 = np.sqrt(np.clip(q2max - (mx**2/(2*Ex))**2,0,None))
        pt2low = np.clip(Ex*theta - delta_qE_m2, 1e-20, None)**2 # setting minimal pt2 to avoid division by 0
        pt2up = (Ex*theta + delta_qE_m2)**2
        
        prefactor = (16*g_a**2*c.g_EM**4*np.sin(theta))/(np.pi*E_beam) * 3.89379*1e8 # conversion to pb 
        f_integrand = lambda phi, pt2: photon_density(Ex/E_beam,pt2) *  Photoproduction_functions.dsigma_gammaN_alt(mx, Ex, np.sqrt(pt2), theta, phi) 
        if use_trapz: return prefactor * integrate_2d_trapz(f_integrand, pt2low, pt2up, philow, phiup)
        integral, error_est = dblquad(f_integrand, pt2low, pt2up, philow, phiup, epsrel=1e-5)
        # if error_est*1e2 > integral: print(f"[Warning:] \t Folding integral relative error estimate ({error_est/integral:.3f}) above 1% at (m_x [GeV], e_x [GeV], th_x)={(mx, Ex, theta)}")

        return prefactor * integral
# modes/mixing.py

import os
import numpy as np
from ALP_production.general import exotic_production_setup as setup
from ALP_production.general import exotic_constants as c
from ALP_production.general import exotic_functions as f
from sklearn.neighbors import KernelDensity
# from scipy import integrate
# import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from os import path
from multiprocessing import Pool
from tqdm import tqdm
import sys

class mixing_production:
    def __init__(self, experiment, n_events,daugther_name,use_external=False,fix_energy_in_cm=False):
        self.exp = experiment
        self.nevts = n_events
        self.daughter_exo = daugther_name
        self.pp_beta = f.beta_lor(setup.p_beam[self.exp],np.sqrt(setup.p_beam[self.exp]**2+c.m_p**2)+c.m_p)
        self.pp_gamma = f.gamma_lor(self.pp_beta)
        self.fix_energy_in_cm = fix_energy_in_cm
        if fix_energy_in_cm:
            print("[Info:] \t Fixing energy in meson CM frame")
            self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced in mixing processes, fixing energy in CM frame of original meson | generated with statistics of {n_events} events using "
        else:
            print("[Info:] \t Fixing momentum in meson CM frame")
            self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced in mixing processes, fixing momentum in CM frame of original meson | generated with statistics of {n_events} events using "
        self.fail_counter = 0
        
        self.th_a_list = [ th_a for th_a in np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins)]
        self.en_a_list = [ en_a for en_a in np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)),  num=setup.energy_bins)]
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())

        if self.daughter_exo == "alp":
            self.pi0_list = []
            self.eta_list = []
            self.etap_list = []
        elif self.daughter_exo == "dp":
            self.rho_list = []
            self.omega_list = []
            self.phi_list = []
        else:
            print("[Error:] \t Mixing production for "+self.daughter_exo+" not implemented")
            sys.exit(1)

        if use_external:
            print("[Info:] \t Using external meson source")
            filename = os.path.dirname(os.path.realpath(__file__))+"/../tab_mesons/softQCD/softQCD_"+str(self.nevts//1000)+"kevts_pp_"+str(setup.p_beam[self.exp])+"GeV_8.2_pSet2_pom_ok.txt"
            if setup.p_beam[self.exp] == 30: #different file for KOTO
                filename = os.path.dirname(os.path.realpath(__file__))+"/../tab_mesons/softQCD/softQCD_200kevts_pp_30GeV_8.3_pSet2_pom_ok.txt"
            elif setup.p_beam[self.exp] == 800:
                filename = os.path.dirname(os.path.realpath(__file__))+"/../tab_mesons/softQCD/softQCD_100kevts_pp_800GeV_8.3_pSet2_pom_ok.txt"

            if os.path.exists(filename):
                external_list = np.loadtxt(filename, usecols=(2,3,4,5))
                with open(filename) as meson_source_file: self.export_header += meson_source_file.readline()[1:-1]
                for iEvent in external_list:
                    event_mom = np.sqrt(iEvent[1]**2+iEvent[2]**2+iEvent[3]**2)
                    if event_mom <= 0: continue
                    event_th = np.arccos(iEvent[3]/event_mom)
                    if event_th <= 0: event_th = 1E-10 #set min value for log
                    if self.daughter_exo == "alp":
                        if np.abs(iEvent[0]) == 111:#pi0
                            self.pi0_list.append([np.log(c.m_pi0**2+event_mom**2)/2, np.log(event_th)])
                        if np.abs(iEvent[0]) == 221:#eta
                            self.eta_list.append([np.log(c.m_eta*2+event_mom**2)/2, np.log(event_th)])
                        if np.abs(iEvent[0]) == 331:#eta'
                            self.etap_list.append([np.log(c.m_etap**2+event_mom**2)/2, np.log(event_th)])
                    elif self.daughter_exo == "dp":                    
                        if np.abs(iEvent[0]) == 113:#rho
                            self.rho_list.append([np.log(c.m_rho**2+event_mom**2)/2, np.log(event_th)])
                        if np.abs(iEvent[0]) == 223: #omega
                            self.omega_list.append([np.log(c.m_omega*2+event_mom**2)/2, np.log(event_th)])
                        if np.abs(iEvent[0]) == 333: #phi
                            self.phi_list.append([np.log(c.m_phi**2+event_mom**2)/2, np.log(event_th)])

            else:
                print("[Error:] \t Path:",filename,"does not exist. Exiting")
                exit()

        else:
            print("[Info:] \t Starting Pythia")
            if setup.p_beam[self.exp] == 30: #different normalization for KOTO
                self.nevts = 2*1E5
            self._init_pythia_()

            print("[Info:] \t Pythia production finished")

        if setup.p_beam[self.exp] == 30: #different normalization for KOTO
            self.nevts = 0.956*2*1E5

        if self.daughter_exo == "alp":
            self.mult_pi0 = len(self.pi0_list)/self.nevts
            self.mult_eta = len(self.eta_list)/self.nevts
            self.mult_etap = len(self.etap_list)/self.nevts
            print("[Info:] \t Meson multiplicity Pi0:",self.mult_pi0,"Eta:",self.mult_eta,"Eta':",self.mult_etap)

            print("[Info:] \t Interpolating")
            #use the same settings as Mathematica default
            self.kde_pi0 = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.pi0_list))
            self.kde_eta = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.eta_list))
            self.kde_etap = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.etap_list))

        elif self.daughter_exo == "dp":
            self.mult_rho = len(self.rho_list)/self.nevts
            self.mult_omega = len(self.omega_list)/self.nevts
            self.mult_phi = len(self.phi_list)/self.nevts
            print("[Info:] \t Meson multiplicity Rho:",self.mult_rho,"Omega:",self.mult_omega,"Phi:",self.mult_phi)

            print("[Info:] \t Interpolating")
            #use the same settings as Mathematica default
            self.kde_rho = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.rho_list))
            self.kde_omega = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.omega_list))
            self.kde_phi = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.phi_list))


    def _init_pythia_(self):
        #import pythia
        if not os.path.exists(setup.pythia_dir):
            print("[Error:] \t Set correct path to Pythia in general/exotic_production_setup.py")
            sys.exit(1)
        cfg = open(setup.pythia_dir+"/examples/Makefile.inc")
        lib = setup.pythia_dir+"/lib"
        for line in cfg:
            if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
        sys.path.insert(0, lib)
        import pythia8

        pythia = pythia8.Pythia("", False)
        pythia_modifiers = ["SoftQCD:all = on", "111:mayDecay = no", "Beams:idA = 2212","Beams:idB = 2212","Beams:eA = " + str(setup.p_beam[self.exp]),"Beams:eB = 0","Beams:frameType = 2"]
        for flag in pythia_modifiers+["Print:quiet=on","Print:errors=off","Next:numberShowEvent = 0"]: pythia.readString(flag)
        self.export_header += f"Pythia{pythia.parm('Pythia:versionNumber')} with settings [" +" ".join(pythia_modifiers)+"]"

        pythia.init() # initialize

        for _ in range(self.nevts):
            pythia.next()
            for isubEvent in range(pythia.event.size()):
                event_en = pythia.event[isubEvent].e()
                if event_en <= 0: continue
                event_th = pythia.event[isubEvent].theta()
                if event_th <= 0: event_th = 1E-10 #set min value for log
                evtID = np.abs(pythia.event[isubEvent].id())
                if self.daughter_exo == "alp":
                    if evtID == 111:#pi0
                        self.pi0_list.append([np.log(event_en), np.log(event_th)])
                    if evtID == 221:#eta
                        self.eta_list.append([np.log(event_en), np.log(event_th)])
                    if evtID == 331:#eta'
                        self.etap_list.append([np.log(event_en), np.log(event_th)])
                elif self.daughter_exo == "dp":
                    if evtID == 113:#rho
                        self.rho_list.append([np.log(event_en), np.log(event_th)])
                    if evtID == 223:#omega
                        self.omega_list.append([np.log(event_en), np.log(event_th)])
                    if evtID == 333:#phi
                        self.phi_list.append([np.log(event_en), np.log(event_th)])

        return

    def is_mom_valid(self,m,p_orig,p_new = 0, is_theta = True, is_energy = False, is_lab = False):
        if np.isnan(p_new):
            if self.fail_counter < 3:
                if is_theta:
                    print("[Warning:] \t Meson momentum estimation failed. Energy:",p_orig,"Mass:",m, "Setting theta=pi")
            
                else:
                    if is_lab:
                        orig_frame = "Lab."
                        new_frame = "CM"
                    else:
                        orig_frame = "CM"
                        new_frame = "Lab"
                        
                    if is_energy:
                        print("[Warning:] \t "+orig_frame+" energy estimation failed:", p_new,"setting 0. "+new_frame+" energy:",p_orig, "mass:",m,"gamma:", self.pp_gamma,"beta:",self.pp_beta)
                    else:
                        print("[Warning:] \t "+orig_frame+" momentum estimation failed:", p_new,"setting 0. "+new_frame+" momentum:",p_orig, "mass:",m,"gamma:", self.pp_gamma,"beta:",self.pp_beta)
            
            if self.fail_counter == 3:
                print("[Warning:] \t Meson momentum or energy estimation failed (4), suppressing further warnings")
            self.fail_counter+=1
            return False
        
        return True

    def e_lab(self,e_cm,m):
        if e_cm<m:
            self.is_mom_valid(m,e_cm,p_new = np.nan, is_theta = False, is_energy = True, is_lab = True)
            return 0
        return self.pp_gamma*e_cm + self.pp_gamma * self.pp_beta * np.sqrt(e_cm**2-m**2)

    def e_cm(self,e_lab,m):
        if e_lab<m:
            self.is_mom_valid(m,e_lab,p_new = np.nan, is_theta = False, is_energy = True, is_lab = False)
            return 0
        return self.pp_gamma*e_lab - self.pp_gamma * self.pp_beta * np.sqrt(e_lab**2-m**2)

    def p_lab(self,p_cm,m):
        p_out = self.pp_gamma*p_cm + self.pp_gamma * self.pp_beta * np.sqrt(p_cm**2+m**2)
        if not self.is_mom_valid(m,p_cm,p_out,is_theta=False,is_lab=True): p_out = 0
        return p_out

    def p_cm(self,p_lab,m):
        p_out = self.pp_gamma*p_lab - self.pp_gamma * self.pp_beta * np.sqrt(p_lab**2+m**2)
        if not self.is_mom_valid(m,p_lab,p_out,is_theta=False,is_lab=False): p_out = 0
        return p_out

    def theta_meson_p(self,theta_a,e_a,m_a,m_meson):
        if e_a < m_a:
            self.is_mom_valid(m_a,e_a,p_new = np.nan)
            return np.pi 
        p_a = np.sqrt(e_a**2 - m_a**2)
        mom_a_cm = self.p_cm(p_a,m_a)
        mom_meson_lab = self.p_lab(mom_a_cm,m_meson)
        if mom_meson_lab == 0: mom_meson_lab = 1E-6
        s_th_out = abs(np.sin(theta_a) * p_a/mom_meson_lab) #positive angle only
        theta_out = 0
        if s_th_out <= 1:
             theta_out = np.arcsin(s_th_out)
        else:
            theta_out = np.pi
        return theta_out 
    
    def theta_meson_e(self,theta_a,e_a,m_a,m_meson): # theta for fixed E in cm frame
        if e_a < m_a:
            self.is_mom_valid(m_a,e_a,p_new = np.nan)
            return np.pi 
        p_a = np.sqrt(e_a**2 - m_a**2)
        mom_meson_lab_2 = self.e_lab(self.e_cm(e_a,m_a),m_meson)**2-m_meson**2
        mom_meson_lab = 0
        if mom_meson_lab_2 >= 0:
            mom_meson_lab = np.sqrt(mom_meson_lab_2)
        if mom_meson_lab == 0: mom_meson_lab = 1E-6
        s_th_out = abs(np.sin(theta_a) * p_a/mom_meson_lab) #positive angle only
        theta_out = 0
        if s_th_out <= 1:
             theta_out = np.arcsin(s_th_out)
        else:
            theta_out = np.pi
        return theta_out 
    
    def process_pool(self,nthreads,single_mass_point=0):
        """Wrapper to run production in parallel threads.
        Args:
            nthreads (int): number of parallel threads to run
            single_mass_point (float, optional): Single mass point for which to run. Defaults to 0.
        """
        export_units_header = "Theta_X[rad] E_X [GeV] m_X [GeV] dY[per(rad GeV N_pN theta^2)]"
        if single_mass_point>0:
            massName = str(int(single_mass_point*1e6))+'keV'
            if self.daughter_exo == "alp":
                list_pi0, list_eta, list_etap = self.ALP_mixing_process_e(single_mass_point) if self.fix_energy_in_cm else self.ALP_mixing_process_p(single_mass_point)
                f.export(setup.experiments[self.exp],self.daughter_exo+"_mixingPi0_beam" + str(setup.p_beam[self.exp]) + "GeV_" + massName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_pi0,(len(self.th_a_list)*len(self.en_a_list),4)), header = self.export_header.replace('mixing processes','Pi0 mixing')+'\n'+export_units_header.replace('theta','th_pi'))
                f.export(setup.experiments[self.exp],self.daughter_exo+"_mixingEta_beam" + str(setup.p_beam[self.exp]) + "GeV_" + massName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_eta,(len(self.th_a_list)*len(self.en_a_list),4)), header = self.export_header.replace('mixing processes','Eta mixing')+'\n'+export_units_header.replace('theta','th_eta'))
                f.export(setup.experiments[self.exp],self.daughter_exo+"_mixingEtaPrim_beam" + str(setup.p_beam[self.exp]) + "GeV_" + massName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_etap,(len(self.th_a_list)*len(self.en_a_list),4)), header = self.export_header.replace('mixing processes','EtaPrime mixing')+'\n'+export_units_header.replace('theta','th_etap'))
            elif self.daughter_exo == "dp":
                list_rho0, list_omega, list_phi = self.DP_mixing_process_e(single_mass_point) if self.fix_energy_in_cm else self.DP_mixing_process_p(single_mass_point)
                f.export(setup.experiments[self.exp],self.daughter_exo+"_MixingRho_beam" + str(setup.p_beam[self.exp]) + "GeV_" + massName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_rho0,(len(self.th_a_list)*len(self.en_a_list),4)), header = self.export_header.replace('mixing processes','Rho mixing')+'\n'+export_units_header.replace('theta','th_rho'))
                f.export(setup.experiments[self.exp],self.daughter_exo+"_MixingOmega_beam" + str(setup.p_beam[self.exp]) + "GeV_" + massName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_omega,(len(self.th_a_list)*len(self.en_a_list),4)), header = self.export_header.replace('mixing processes','Omega mixing')+'\n'+export_units_header.replace('theta','th_omega'))
                f.export(setup.experiments[self.exp],self.daughter_exo+"_MixingPhi_beam" + str(setup.p_beam[self.exp]) + "GeV_" + massName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_phi,(len(self.th_a_list)*len(self.en_a_list),4)), header = self.export_header.replace('mixing processes','Phi mixing')+'\n'+export_units_header.replace('theta','th_phi'))
            return

        for iFile in range(len(setup.mass_bins)):
            iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")

            if self.daughter_exo == "alp":
                with Pool(processes=nthreads) as pool:
                    if self.fix_energy_in_cm:
                        list_pi0, list_eta, list_etap = zip(*tqdm(pool.imap(self.ALP_mixing_process_e,self.m_lists[iFile]),total=len(self.m_lists[iFile])))
                    else:
                        list_pi0, list_eta, list_etap = zip(*tqdm(pool.imap(self.ALP_mixing_process_p,self.m_lists[iFile]),total=len(self.m_lists[iFile])))
                        

                pool.close()
                pool.join() 

                # export
                f.export(setup.experiments[self.exp],self.daughter_exo+"_mixingPi0_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_pi0,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)), header = self.export_header.replace('mixing processes','Pi0 mixing')+'\n'+export_units_header.replace('theta','th_pi'))
                f.export(setup.experiments[self.exp],self.daughter_exo+"_mixingEta_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_eta,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)), header = self.export_header.replace('mixing processes','Eta mixing')+'\n'+export_units_header.replace('theta','th_eta'))
                f.export(setup.experiments[self.exp],self.daughter_exo+"_mixingEtaPrim_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_etap,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)), header = self.export_header.replace('mixing processes','EtaPrime mixing')+'\n'+export_units_header.replace('theta','th_etap'))

            elif self.daughter_exo == "dp":
                with Pool(processes=nthreads) as pool:
                    if self.fix_energy_in_cm:
                        list_rho0, list_omega, list_phi = zip(*tqdm(pool.imap(self.DP_mixing_process_e,self.m_lists[iFile]),total=len(self.m_lists[iFile])))
                    else:
                        list_rho0, list_omega, list_phi = zip(*tqdm(pool.imap(self.DP_mixing_process_p,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

                pool.close()
                pool.join() 
                
                # export
                f.export(setup.experiments[self.exp],self.daughter_exo+"_MixingRho_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_rho0,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)), header = self.export_header.replace('mixing processes','Rho mixing')+'\n'+export_units_header.replace('theta','th_rho'))
                f.export(setup.experiments[self.exp],self.daughter_exo+"_MixingOmega_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_omega,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)), header = self.export_header.replace('mixing processes','Omega mixing')+'\n'+export_units_header.replace('theta','th_omega'))
                f.export(setup.experiments[self.exp],self.daughter_exo+"_MixingPhi_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_phi,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)), header = self.export_header.replace('mixing processes','Phi mixing')+'\n'+export_units_header.replace('theta','th_phi'))

        return

    def ALP_mixing_process_p(self,m_a):
        """evaluates axion-like-particles with mass m_a mixing with light pseudoscalar mesons emitted at angle theta to beam axis with energy e
        Args:
            m_a (float): alp mass

        Returns:
            array(float): [theta, e, m_a, expected yield]
        """

        data_list_pi0 = []
        data_list_eta = []
        data_list_etap = []
        for e_a in self.en_a_list:
            if e_a > m_a:
                p_a = np.sqrt(e_a**2 - m_a**2)
                e_pi0 = np.sqrt(self.p_lab(self.p_cm(p_a, m_a), c.m_pi0)**2 + c.m_pi0**2)
                e_eta = np.sqrt(self.p_lab(self.p_cm(p_a, m_a), c.m_eta)**2 + c.m_eta**2)
                e_etap = np.sqrt(self.p_lab(self.p_cm(p_a, m_a), c.m_etap)**2 + c.m_etap**2)

                data_sublist_pi0 = []
                data_sublist_eta = []
                data_sublist_etap = []

                for th_a in self.th_a_list:
                    data_sublist_pi0.append([th_a, e_a, m_a, f.f_kde(self.kde_pi0,self.mult_pi0,e_pi0,self.theta_meson_p(th_a,e_a,m_a,c.m_pi0))[0]])
                    data_sublist_eta.append([th_a, e_a, m_a, f.f_kde(self.kde_eta,self.mult_eta,e_eta,self.theta_meson_p(th_a,e_a,m_a,c.m_eta))[0]])
                    data_sublist_etap.append([th_a, e_a, m_a, f.f_kde(self.kde_etap,self.mult_etap,e_etap,self.theta_meson_p(th_a,e_a,m_a,c.m_etap))[0]])

                data_list_pi0.append([data_sublist_pi0])
                data_list_eta.append([data_sublist_eta])
                data_list_etap.append([data_sublist_etap])
            else:
                data_sublist_pi0 = []
                data_sublist_eta = []
                data_sublist_etap = []

                for th_a in self.th_a_list:
                    data_sublist_pi0.append([th_a, e_a, m_a, 0])
                    data_sublist_eta.append([th_a, e_a, m_a, 0])
                    data_sublist_etap.append([th_a, e_a, m_a, 0])

                data_list_pi0.append([data_sublist_pi0])
                data_list_eta.append([data_sublist_eta])
                data_list_etap.append([data_sublist_etap])

        return data_list_pi0, data_list_eta, data_list_etap

    def ALP_mixing_process_e(self,m_a): # process for fixed E in cm frame

        data_list_pi0 = []
        data_list_eta = []
        data_list_etap = []
        for e_a in self.en_a_list:
            if e_a > m_a:
                e_pi0 = self.e_lab(self.e_cm(e_a,m_a),c.m_pi0)
                e_eta = self.e_lab(self.e_cm(e_a,m_a),c.m_eta)
                e_etap = self.e_lab(self.e_cm(e_a,m_a),c.m_etap)

                data_sublist_pi0 = []
                data_sublist_eta = []
                data_sublist_etap = []

                for th_a in self.th_a_list:
                    data_sublist_pi0.append([th_a, e_a, m_a, f.f_kde(self.kde_pi0,self.mult_pi0,e_pi0,self.theta_meson_e(th_a,e_a,m_a,c.m_pi0))[0]])
                    data_sublist_eta.append([th_a, e_a, m_a, f.f_kde(self.kde_eta,self.mult_eta,e_eta,self.theta_meson_e(th_a,e_a,m_a,c.m_eta))[0]])
                    data_sublist_etap.append([th_a, e_a, m_a, f.f_kde(self.kde_etap,self.mult_etap,e_etap,self.theta_meson_e(th_a,e_a,m_a,c.m_etap))[0]])

                data_list_pi0.append([data_sublist_pi0])
                data_list_eta.append([data_sublist_eta])
                data_list_etap.append([data_sublist_etap])
            else:
                data_sublist_pi0 = []
                data_sublist_eta = []
                data_sublist_etap = []

                for th_a in self.th_a_list:
                    data_sublist_pi0.append([th_a, e_a, m_a, 0])
                    data_sublist_eta.append([th_a, e_a, m_a, 0])
                    data_sublist_etap.append([th_a, e_a, m_a, 0])

                data_list_pi0.append([data_sublist_pi0])
                data_list_eta.append([data_sublist_eta])
                data_list_etap.append([data_sublist_etap])

        return data_list_pi0, data_list_eta, data_list_etap

    def DP_mixing_process_p(self,m_X):
        """evaluates dark photons with mass m_X mixing with light vector mesons emitted at angle theta_X to beam axis with energy e_X
        Args:
            m_X (float): dark photon mass

        Returns:
            array(float): [theta_X, e_X, m_X, expected yield]
        """

        data_list_rho = []
        data_list_omega = []
        data_list_phi = []
        
        for e_X in self.en_a_list:
            if e_X > m_X:
                p_X = np.sqrt(e_X**2 - m_X**2)
                e_rho = np.sqrt(self.p_lab(self.p_cm(p_X, m_X), c.m_rho)**2 + c.m_rho**2)
                e_omega = np.sqrt(self.p_lab(self.p_cm(p_X, m_X), c.m_omega)**2 + c.m_omega**2)
                e_phi = np.sqrt(self.p_lab(self.p_cm(p_X, m_X), c.m_phi)**2 + c.m_phi**2)

                data_sublist_rho = []
                data_sublist_omega = []
                data_sublist_phi = []

                for th_X in self.th_a_list:
                    data_sublist_rho.append([th_X, e_X, m_X, f.f_kde(self.kde_rho,self.mult_rho,e_rho,self.theta_meson_p(th_X,e_X,m_X,c.m_rho))[0]])
                    data_sublist_omega.append([th_X, e_X, m_X, f.f_kde(self.kde_omega,self.mult_omega,e_omega,self.theta_meson_p(th_X,e_X,m_X,c.m_omega))[0]])
                    data_sublist_phi.append([th_X, e_X, m_X, f.f_kde(self.kde_phi,self.mult_phi,e_phi,self.theta_meson_p(th_X,e_X,m_X,c.m_phi))[0]])

                data_list_rho.append([data_sublist_rho])
                data_list_omega.append([data_sublist_omega])
                data_list_phi.append([data_sublist_phi])
            else:
                data_sublist_rho = []
                data_sublist_omega = []
                data_sublist_phi = []

                for th_X in self.th_a_list:
                    data_sublist_rho.append([th_X, e_X, m_X, 0])
                    data_sublist_omega.append([th_X, e_X, m_X, 0])
                    data_sublist_phi.append([th_X, e_X, m_X, 0])

                data_list_rho.append([data_sublist_rho])
                data_list_omega.append([data_sublist_omega])
                data_list_phi.append([data_sublist_phi])

        return data_list_rho, data_list_omega, data_list_phi

    def DP_mixing_process_e(self,m_X): # process for fixed E in cm frame

        data_list_rho = []
        data_list_omega = []
        data_list_phi = []
        for e_X in self.en_a_list:
            if e_X > m_X:
                e_rho = self.e_lab(self.e_cm(e_X,m_X),c.m_rho)
                e_omega = self.e_lab(self.e_cm(e_X,m_X),c.m_omega)
                e_phi = self.e_lab(self.e_cm(e_X,m_X),c.m_phi)

                data_sublist_rho = []
                data_sublist_omega = []
                data_sublist_phi = []

                for th_X in self.th_a_list:
                    data_sublist_rho.append([th_X, e_X, m_X, f.f_kde(self.kde_rho,self.mult_rho,e_rho,self.theta_meson_e(th_X,e_X,m_X,c.m_rho))[0]])
                    data_sublist_omega.append([th_X, e_X, m_X, f.f_kde(self.kde_omega,self.mult_omega,e_omega,self.theta_meson_e(th_X,e_X,m_X,c.m_omega))[0]])
                    data_sublist_phi.append([th_X, e_X, m_X, f.f_kde(self.kde_phi,self.mult_phi,e_phi,self.theta_meson_e(th_X,e_X,m_X,c.m_phi))[0]])

                data_list_rho.append([data_sublist_rho])
                data_list_omega.append([data_sublist_omega])
                data_list_phi.append([data_sublist_phi])
            else:
                data_sublist_rho = []
                data_sublist_omega = []
                data_sublist_phi = []

                for th_X in self.th_a_list:
                    data_sublist_rho.append([th_X, e_X, m_X, 0])
                    data_sublist_omega.append([th_X, e_X, m_X, 0])
                    data_sublist_phi.append([th_X, e_X, m_X, 0])

                data_list_rho.append([data_sublist_rho])
                data_list_omega.append([data_sublist_omega])
                data_list_phi.append([data_sublist_phi])

        return data_list_rho, data_list_omega, data_list_phi

class mixing_production_var_mass: #to be finalized
    def __init__(self,experiment, n_events):
        self.exp = experiment
        self.nevts = n_events
        # self.beam_beta = f.beta_lor(setup.p_beam[self.exp],np.sqrt(setup.p_beam[self.exp]**2+c.m_p**2))
        # self.beam_gamma = f.gamma_lor(self.beam_beta)
        self.pp_beta = f.beta_lor(setup.p_beam[self.exp],np.sqrt(setup.p_beam[self.exp]**2+c.m_p**2)+c.m_p)
        self.pp_gamma = f.gamma_lor(self.pp_beta)

        self.th_a_list = [ th_a for th_a in np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins)]
        self.en_a_list = [ en_a for en_a in np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)),  num=setup.energy_bins)]
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())

    # def mute(self):
    #     sys.stdout = open(os.devnull, 'w') 

    def e_lab(self,e_cm,m):
        if e_cm<m:
            return 0
        return self.pp_gamma*e_cm + self.pp_gamma * self.pp_beta * np.sqrt(e_cm**2-m**2)

    def e_cm(self,e_lab,m):
        if e_lab<m:
            return 0
        return self.pp_gamma*e_lab - self.pp_gamma * self.pp_beta * np.sqrt(e_lab**2-m**2)

    def p_lab(self,p_cm,m):
        p_out = self.pp_gamma*p_cm + self.pp_gamma * self.pp_beta * np.sqrt(p_cm**2+m**2)
        if np.isnan(p_out):
            print("[Warning:] \t Lab. momentum estimation failed:", p_out,"setting 0. CM momentum:",p_cm, "mass:",m,"gamma:", self.pp_gamma,"beta:",self.pp_beta)
            p_out = 0
        return p_out

    def p_cm(self,p_lab,m):
        p_out = self.pp_gamma*p_lab - self.pp_gamma * self.pp_beta * np.sqrt(p_lab**2+m**2)
        if np.isnan(p_out):
            print("[Warning:] \t CM momentum estimation failed:", p_out,"setting 0. Lab momentum:",p_lab, "mass:",m,"gamma:", self.pp_gamma,"beta:",self.pp_beta)
            p_out = 0
        return p_out

    def process_pool(self,nthreads,single_mass_point=0):
        """Wrapper to run production in parallel threads.
        Args:
            nthreads (int): number of parallel threads to run
            single_mass_point (float, optional): Single mass point for which to run. Defaults to 0.
        """
        if single_mass_point>0:
            list_pi0, list_eta, list_etap = self.ALP_mixing_process(single_mass_point)
            fileName = + str(int(single_mass_point*1e6)) + "keV_" 
            f.export(setup.experiments[self.exp],"alp_mixingPi0_varMass_beam" + str(setup.p_beam[self.exp]) + "GeV_" + fileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_pi0,(len(self.th_a_list)*len(self.en_a_list),4)))
            f.export(setup.experiments[self.exp],"alp_mixingEta_varMass_beam" + str(setup.p_beam[self.exp]) + "GeV_" + fileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_eta,(len(self.th_a_list)*len(self.en_a_list),4)))
            f.export(setup.experiments[self.exp],"alp_mixingEtaPrim_varMass_beam" + str(setup.p_beam[self.exp]) + "GeV_" + fileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_etap,(len(self.th_a_list)*len(self.en_a_list),4)))
            return

        for iFile in range(len(setup.mass_bins)):
            iFileName = str(setup.mass_min[iFile]*1000) + "to" + str(setup.mass_max[iFile]*1000) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(setup.mass_max[iFile]*1000) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")

            with Pool(processes=nthreads) as pool:
                list_pi0, list_eta, list_etap = zip(*tqdm(pool.imap(self.ALP_mixing_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

            pool.close() #check utilization
            pool.join() #check utilization
            
            # export        
            f.export(setup.experiments[self.exp],"alp_mixingPi0_varMass_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_pi0,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)))
            f.export(setup.experiments[self.exp],"alp_mixingEta_varMass_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_eta,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)))
            f.export(setup.experiments[self.exp],"alp_mixingEtaPrim_varMass_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_etap,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)))

        return
  
    def ALP_mixing_process(self,m_a):
        #import pythia
        if not os.path.exists(setup.pythia_dir):
            print("[Error:] \t Set correct path to Pythia in general/exotic_production_setup.py")
            sys.exit(1)
        cfg = open(setup.pythia_dir+"/examples/Makefile.inc")
        lib = setup.pythia_dir+"/lib"
        for line in cfg:
            if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
        sys.path.insert(0, lib)
        import pythia8

        pythia = pythia8.Pythia("", False)
        #verbosity
        # pythia.readString("Check:levelParticleData  = 0")
        # pythia.readString("Check:event = off")
        # pythia.readString("Check:history  = off")
        #setup
        pythia.readString("SoftQCD:all = on")
        pythia.readString("pdf:pHardSet = 2")
        pythia.readString("Next:numberShowEvent = 5")
        pythia.readString("111:mayDecay = no")
        pythia.readString("Beams:idA = 2212")
        pythia.readString("Beams:idB = 2212")
        pythia.readString("Beams:eA = " + str(setup.p_beam[self.exp]))
        pythia.readString("Beams:eB = 0")
        pythia.readString("Beams:frameType = 2")

        pythia.readString("111:m0 = " + str(m_a))
        pythia.readString("221:m0 = " + str(m_a))
        pythia.readString("331:m0 = " + str(m_a))

        pythia.init() # initialize

        counter=0

        data_pythia_pi0 = []
        data_pythia_eta = []
        data_pythia_etap = []

        for iEvent in range(self.nevts):
            pythia.next()
            for isubEvent in range(pythia.event.size()):
                event_en = pythia.event[isubEvent].e()
                if event_en <= 0: continue
                event_th = pythia.event[isubEvent].theta()
                if event_th <= 0: event_th = 1E-10 #set min value for log
                evtID = np.abs(pythia.event[isubEvent].id())
                if evtID == 111:#pi0
                    data_pythia_pi0.append([np.log(event_en), np.log(event_th)])
                if evtID == 221:#eta
                    data_pythia_eta.append([np.log(event_en), np.log(event_th)])
                if evtID == 331:#eta'
                    data_pythia_etap.append([np.log(event_en), np.log(event_th)])
                counter+=1

        mult_pi0 = len(data_pythia_pi0)/self.nevts
        mult_eta = len(data_pythia_eta)/self.nevts
        mult_etap = len(data_pythia_etap)/self.nevts

        #use the same settings as Mathematica default
        kde_pi0 = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(data_pythia_pi0))
        kde_eta = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(data_pythia_eta))
        kde_etap = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(data_pythia_etap))

        data_list_pi0 = []
        data_list_eta = []
        data_list_etap = []
        for th in self.th_a_list:
            data_sublist_pi0 = []
            data_sublist_eta = []
            data_sublist_etap = []
            for e_a in self.en_a_list:
                data_sublist_pi0.append([th, e_a, m_a, f.f_kde(kde_pi0,mult_pi0,e_a,th)[0]])
                data_sublist_eta.append([th, e_a, m_a, f.f_kde(kde_eta,mult_eta,e_a,th)[0]])
                data_sublist_etap.append([th, e_a, m_a, f.f_kde(kde_etap,mult_etap,e_a,th)[0]])
            data_list_pi0.append([data_sublist_pi0])
            data_list_eta.append([data_sublist_eta])
            data_list_etap.append([data_sublist_etap])

        return data_list_pi0, data_list_eta, data_list_etap
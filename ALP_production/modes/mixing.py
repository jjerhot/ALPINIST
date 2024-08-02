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
    def __init__(self, experiment, n_events,daugther_name,use_external=False):
        self.exp = experiment
        self.nevts = n_events
        self.daughter_exo = daugther_name
        self.pp_beta = f.beta_lor(setup.p_beam[self.exp],np.sqrt(setup.p_beam[self.exp]**2+c.m_p**2)+c.m_p)
        self.pp_gamma = f.gamma_lor(self.pp_beta)
        self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced in mixing processes | generated with statistics of {n_events} events using "

        self.th_a_list = [ th_a for th_a in np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins)]
        self.en_a_list = [ en_a for en_a in np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)),  num=setup.energy_bins)]
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())

        self.pi0_list = []
        self.eta_list = []
        self.etap_list = []

        if use_external:
            print("[Info:] \t Using external meson source")
            filename = os.getcwd()+"/tab_mesons/softQCD/softQCD_"+str(self.nevts//1000)+"kevts_pp_"+str(setup.p_beam[self.exp])+"GeV_8.2_pSet2_pom_ok.txt"
            if setup.p_beam[self.exp] == 30: #different file for KOTO
                filename = os.getcwd()+"/tab_mesons/softQCD/softQCD_200kevts_pp_30GeV_8.3_pSet2_pom_ok.txt"
            elif setup.p_beam[self.exp] == 800:
                filename = os.getcwd()+"/tab_mesons/softQCD/softQCD_100kevts_pp_800GeV_8.3_pSet2_pom_ok.txt"

            if os.path.exists(filename):
                external_list = np.loadtxt(filename)
                with open(filename) as meson_source_file: self.export_header += meson_source_file.readline()[1:-1]
                external_list = np.delete(external_list, 6, 1)
                external_list = np.delete(external_list, 0, 1)
                external_list = np.delete(external_list, 0, 1)
                for iEvent in external_list:
                    event_mom = np.sqrt(iEvent[1]**2+iEvent[2]**2+iEvent[3]**2)
                    if event_mom <= 0: continue
                    event_th = np.arccos(iEvent[3]/event_mom)
                    if event_th <= 0: event_th = 1E-10 #set min value for log
                    if np.abs(iEvent[0]) == 111:#pi0
                        self.pi0_list.append([np.log(c.m_pi0**2+event_mom**2)/2, np.log(event_th)])
                    if np.abs(iEvent[0]) == 221:#eta
                        self.eta_list.append([np.log(c.m_eta*2+event_mom**2)/2, np.log(event_th)])
                    if np.abs(iEvent[0]) == 331:#eta'
                        self.etap_list.append([np.log(c.m_etap**2+event_mom**2)/2, np.log(event_th)])
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

        self.mult_pi0 = len(self.pi0_list)/self.nevts
        self.mult_eta = len(self.eta_list)/self.nevts
        self.mult_etap = len(self.etap_list)/self.nevts
        print("[Info:] \t Meson multiplicity Pi0:",self.mult_pi0,"Eta:",self.mult_eta,"Eta':",self.mult_etap)

        print("[Info:] \t Interpolating")
        #use the same settings as Mathematica default
        self.kde_pi0 = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.pi0_list))
        self.kde_eta = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.eta_list))
        self.kde_etap = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.etap_list))

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
                if evtID == 111:#pi0
                    self.pi0_list.append([np.log(event_en), np.log(event_th)])
                if evtID == 221:#eta
                    self.eta_list.append([np.log(event_en), np.log(event_th)])
                if evtID == 331:#eta'
                    self.etap_list.append([np.log(event_en), np.log(event_th)])

        return
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

    def theta_meson(self,theta_a,e_a,m_a,m_meson):
        if e_a < m_a:
            print("[Warning:] \t Meson momentum estimation failed. Energy:",e_a,"Mass:",m_a, "Setting theta=pi")
            return np.pi 
        p_a = np.sqrt(e_a**2 - m_a**2)
        mom_a_cm = self.p_cm(p_a,m_a)
        mom_meson_lab = self.p_lab(mom_a_cm,m_meson)
        if mom_meson_lab == 0:
            print("[Warning:] \t Meson momentum estimation failed:",mom_meson_lab, "Setting minimal non-zero value (1keV)")
            mom_meson_lab = 1E-6
        s_th_out = abs(np.sin(theta_a) * p_a/mom_meson_lab) #positive angle only
        theta_out = 0
        if s_th_out <= 1:
             theta_out = np.arcsin(s_th_out)
        else:
            theta_out = np.pi
        return theta_out 
    
    def theta_meson_e(self,theta_a,e_a,m_a,m_meson): # theta for fixed E in cm frame
        if e_a < m_a:
            print("[Warning:] \t Meson momentum estimation failed. Energy:",e_a,"Mass:",m_a, "Setting theta=pi")
            return np.pi 
        p_a = np.sqrt(e_a**2 - m_a**2)
        mom_meson_lab_2 = self.e_lab(self.e_cm(e_a,m_a),m_meson)**2-m_meson**2
        mom_meson_lab = 0
        if mom_meson_lab_2 >= 0:
            mom_meson_lab = np.sqrt(mom_meson_lab_2)
        if mom_meson_lab == 0:
            # print("[Warning:] \t Meson momentum estimation failed, mom^2 = ",mom_meson_lab_2, "Setting minimal non-zero value (1keV)")
            mom_meson_lab = 1E-6
        s_th_out = abs(np.sin(theta_a) * p_a/mom_meson_lab) #positive angle only
        theta_out = 0
        if s_th_out <= 1:
             theta_out = np.arcsin(s_th_out)
        else:
            theta_out = np.pi
        return theta_out 
    
    def process_pool(self,nthreads):

        for iFile in range(len(setup.mass_bins)):
            iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")

            with Pool(processes=nthreads) as pool:
                # list_pi0, list_eta, list_etap = zip(*tqdm(pool.imap(self.ALP_mixing_process_e,self.m_lists[iFile]),total=len(self.m_lists[iFile])))
                list_pi0, list_eta, list_etap = zip(*tqdm(pool.imap(self.ALP_mixing_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

            pool.close()
            pool.join() 

            # export
            export_units_header = "Theta[rad] E_x [GeV] dY[per(rad GeV N_pN theta^2)]"
            f.export(setup.experiments[self.exp],self.daughter_exo+"_mixingPi0_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_pi0,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)),self.export_header.replace('mixing processes','Pi0 mixing')+'\n'+export_units_header.replace('theta','th_pi'))
            f.export(setup.experiments[self.exp],self.daughter_exo+"_mixingEta_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_eta,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)),self.export_header.replace('mixing processes','Eta mixing')+'\n'+export_units_header.replace('theta','th_eta'))
            f.export(setup.experiments[self.exp],self.daughter_exo+"_mixingEtaPrim_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_etap,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)),self.export_header.replace('mixing processes','EtaPrime mixing')+'\n'+export_units_header.replace('theta','th_etap'))

        return

    def ALP_mixing_process(self,m_a):
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
                    data_sublist_pi0.append([th_a, e_a, m_a, f.f_kde(self.kde_pi0,self.mult_pi0,e_pi0,self.theta_meson(th_a,e_a,m_a,c.m_pi0))[0]])
                    data_sublist_eta.append([th_a, e_a, m_a, f.f_kde(self.kde_eta,self.mult_eta,e_eta,self.theta_meson(th_a,e_a,m_a,c.m_eta))[0]])
                    data_sublist_etap.append([th_a, e_a, m_a, f.f_kde(self.kde_etap,self.mult_etap,e_etap,self.theta_meson(th_a,e_a,m_a,c.m_etap))[0]])

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

    def process_pool(self,nthreads):
        """Wrapper to run production in parallel threads.
        Args:
            nthreads (int): number of parallel threads to run
        """

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


class mixing_production_vector:
    def __init__(self, experiment, n_events, daugther_name, use_external=False):
        self.exp = experiment
        self.nevts = n_events
        self.daughter_exo = daugther_name
        # self.beam_beta = f.beta_lor(setup.p_beam[self.exp],np.sqrt(setup.p_beam[self.exp]**2+c.m_p**2))
        # self.beam_gamma = f.gamma_lor(self.beam_beta)
        self.pp_beta = f.beta_lor(setup.p_beam[self.exp],np.sqrt(setup.p_beam[self.exp]**2+c.m_p**2)+c.m_p)
        self.pp_gamma = f.gamma_lor(self.pp_beta)
        self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced in light meson mixing | generated with statistics of {n_events} events using "

        self.th_a_list = [ th_a for th_a in np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins)]
        self.en_a_list = [ en_a for en_a in np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)),  num=setup.energy_bins)]
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())


        self.ratio_rho = 0
        self.ratio_omega = 0
        self.ratio_phi = 0
        
        #<V-exo> eq 2.18 from 1801.04847
        if self.daughter_exo == "dp": # dark photon
            self.ratio_rho = 1
            self.ratio_omega = 1
            self.ratio_phi = 1

        elif self.daughter_exo == "dv": # B and B-L vectors
            self.ratio_omega = 4
            self.ratio_phi = 1

        self.load_ratio_rho()
        self.load_ratio_omega()
        self.load_ratio_phi()

        self.rho_list = []
        self.omega_list = []
        self.phi_list = []

        if use_external:
            print("[Info:] \t Using external meson source")
            filename = os.getcwd()+"/tab_mesons/softQCD/softQCD_"+str(self.nevts//1000)+"kevts_pp_"+str(setup.p_beam[self.exp])+"GeV_8.2_pSet2_pom_ok.txt"
            if setup.p_beam[self.exp] == 30: #different file for KOTO
                filename = os.getcwd()+"/tab_mesons/softQCD/softQCD_200kevts_pp_30GeV_8.3_pSet2_pom_ok.txt"
            elif setup.p_beam[self.exp] == 800:
                filename = os.getcwd()+"/tab_mesons/softQCD/softQCD_100kevts_pp_800GeV_8.3_pSet2_pom_ok.txt"

            if os.path.exists(filename):
                external_list = np.loadtxt(filename)
                with open(filename) as meson_source_file: self.export_header += meson_source_file.readline()[1:-1]
                external_list = np.delete(external_list, 6, 1)
                external_list = np.delete(external_list, 0, 1)
                external_list = np.delete(external_list, 0, 1)
                for iEvent in external_list:
                    event_mom = np.sqrt(iEvent[1]**2+iEvent[2]**2+iEvent[3]**2)
                    if event_mom <= 0: continue
                    event_th = np.arccos(iEvent[3]/event_mom)
                    if event_th <= 0: event_th = 1E-10 #set min value for log
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

        self.mult_rho = len(self.rho_list)/self.nevts
        self.mult_omega = len(self.omega_list)/self.nevts
        self.mult_phi = len(self.phi_list)/self.nevts
        print("[Info:] \t Meson multiplicity rho:",self.mult_rho,"omega:",self.mult_omega,"phi:",self.mult_phi)

        print("[Info:] \t Interpolating")
        #use the same settings as Mathematica default
        if self.ratio_rho != 0:
            self.kde_rho = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.rho_list))
        if self.ratio_omega != 0:
            self.kde_omega = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.omega_list))
        if self.ratio_phi != 0:
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
        pythia_modifiers = ["SoftQCD:all = on", "111:mayDecay = no", "111:mayDecay = no","Beams:idA = 2212","Beams:idB = 2212","Beams:eA = " + str(setup.p_beam[self.exp]),"Beams:eB = 0","Beams:frameType = 2"]
        for flag in pythia_modifiers+["Print:quiet=on","Print:errors=off", "Next:numberShowEvent = 0"]: pythia.readString(flag)
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
                if evtID == 113:#rho
                    self.rho_list.append([np.log(event_en), np.log(event_th)])
                if evtID == 223:#omega
                    self.omega_list.append([np.log(event_en), np.log(event_th)])
                if evtID == 333:#phi
                    self.phi_list.append([np.log(event_en), np.log(event_th)])

        return
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

    def theta_meson(self,theta_a,e_a,m_a,m_meson):
        if e_a < m_a:
            print("[Warning:] \t Meson momentum estimation failed. Energy:",e_a,"Mass:",m_a, "Setting theta=pi")
            return np.pi 
        p_a = np.sqrt(e_a**2 - m_a**2)
        mom_a_cm = self.p_cm(p_a,m_a)
        mom_meson_lab = self.p_lab(mom_a_cm,m_meson)
        if mom_meson_lab == 0:
            print("[Warning:] \t Meson momentum estimation failed:",mom_meson_lab, "Setting minimal non-zero value (1keV)")
            mom_meson_lab = 1E-6
        s_th_out = abs(np.sin(theta_a) * p_a/mom_meson_lab) #positive angle only
        theta_out = 0
        if s_th_out <= 1:
             theta_out = np.arcsin(s_th_out)
        else:
            theta_out = np.pi
        return theta_out 
    
    def theta_meson_e(self,theta_a,e_a,m_a,m_meson): # theta for fixed E in cm frame
        if e_a < m_a:
            print("[Warning:] \t Meson momentum estimation failed. Energy:",e_a,"Mass:",m_a, "Setting theta=pi")
            return np.pi 
        p_a = np.sqrt(e_a**2 - m_a**2)
        mom_meson_lab_2 = self.e_lab(self.e_cm(e_a,m_a),m_meson)**2-m_meson**2
        mom_meson_lab = 0
        if mom_meson_lab_2 >= 0:
            mom_meson_lab = np.sqrt(mom_meson_lab_2)
        if mom_meson_lab == 0:
            print("[Warning:] \t Meson momentum estimation failed:",mom_meson_lab, "Setting minimal non-zero value (1keV)")
            mom_meson_lab = 1E-6
        s_th_out = abs(np.sin(theta_a) * p_a/mom_meson_lab) #positive angle only
        theta_out = 0
        if s_th_out <= 1:
             theta_out = np.arcsin(s_th_out)
        else:
            theta_out = np.pi
        return theta_out 
    
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
                # list_mixing = list(tqdm(pool.imap(self.DP_mixing_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))
                list_mixing = list(tqdm(pool.imap(self.DP_mixing_process_e,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

            pool.close()
            pool.join() 
            
            # export
        
            export_units_header = "Theta[rad] E_x [GeV] dY[per(rad GeV N_pN eps^2)]"
            f.export(setup.experiments[self.exp],self.daughter_exo+"_Mixing_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_mixing,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)), header = self.export_info_header+'\n'+export_units_header)
            # f.export(setup.experiments[self.exp],"alp_mixingrho_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_rho,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)))
            # f.export(setup.experiments[self.exp],"alp_mixingomega_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_omega,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)))
            # f.export(setup.experiments[self.exp],"alp_mixingphi_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_phi,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)))

        return

    def DP_mixing_process(self,m_dp):
        """evaluates dark photons with mass m_dp mixing with light vector mesons emitted at angle theta to beam axis with energy e
        Args:
            m_dp (float): dark photon mass

        Returns:
            array(float): [theta, e, m_dp, expected yield]
        """

        #mixing angles (use digitized instead)
        # theta_rho = setup.reference_couplings["mixing"]*m_dp**2/(m_dp**2-c.m_rho**2)
        # theta_omega = setup.reference_couplings["mixing"]*m_dp**2/(m_dp**2-c.m_omega**2)
        # theta_phi = setup.reference_couplings["mixing"]*m_dp**2/(m_dp**2-c.m_phi**2)
        theta_rho = self.get_mixing_rho(m_dp)
        theta_omega = self.get_mixing_omega(m_dp)
        theta_phi = self.get_mixing_phi(m_dp)

        data_list = []
        for e_a in self.en_a_list:
            if e_a > m_dp:
                e_rho = np.sqrt(self.p_lab(self.p_cm(np.sqrt(e_a**2 - m_dp**2), m_dp), c.m_rho)**2 + c.m_rho**2)
                e_omega = np.sqrt(self.p_lab(self.p_cm(np.sqrt(e_a**2 - m_dp**2), m_dp), c.m_omega)**2 + c.m_omega**2)
                e_phi = np.sqrt(self.p_lab(self.p_cm(np.sqrt(e_a**2 - m_dp**2), m_dp), c.m_phi)**2 + c.m_phi**2)

                data_sublist = []

                for th_a in self.th_a_list:
                    yi = 0
                    if theta_rho != 0:
                        yi += f.f_kde(self.kde_rho,theta_rho*self.mult_rho,e_rho,self.theta_meson(th_a,e_a,m_dp,c.m_rho))[0]
                    if theta_omega != 0:
                        yi += f.f_kde(self.kde_omega,theta_omega*self.mult_omega,e_omega,self.theta_meson(th_a,e_a,m_dp,c.m_omega))[0]
                    if theta_phi != 0:
                        yi += f.f_kde(self.kde_phi,theta_phi*self.mult_phi,e_phi,self.theta_meson(th_a,e_a,m_dp,c.m_phi))[0]
                    data_sublist.append([th_a, e_a, m_dp, yi])
 
                data_list.append([data_sublist])
            else:
                data_sublist = []
                for th_a in self.th_a_list:
                    data_sublist.append([th_a, e_a, m_dp, 0])

                data_list.append([data_sublist])

        return data_list

    def DP_mixing_process_e(self,m_dp): # process for fixed E in cm frame

        #mixing angles (use digitized instead)
        # theta_rho = setup.reference_couplings["mixing"]*m_dp**2/(m_dp**2-c.m_rho**2)
        # theta_omega = setup.reference_couplings["mixing"]*m_dp**2/(m_dp**2-c.m_omega**2)
        # theta_phi = setup.reference_couplings["mixing"]*m_dp**2/(m_dp**2-c.m_phi**2)
        theta_rho = self.get_mixing_rho(m_dp)
        theta_omega = self.get_mixing_omega(m_dp)
        theta_phi = self.get_mixing_phi(m_dp)

        data_list = []
        for e_a in self.en_a_list:
            if e_a > m_dp:
                e_rho = self.e_lab(self.e_cm(e_a,m_dp),c.m_rho)
                e_omega = self.e_lab(self.e_cm(e_a,m_dp),c.m_omega)
                e_phi = self.e_lab(self.e_cm(e_a,m_dp),c.m_phi)

                data_sublist = []

                for th_a in self.th_a_list:
                    yi = 0
                    if theta_rho != 0:
                        yi += f.f_kde(self.kde_rho,theta_rho*self.mult_rho,e_rho,self.theta_meson_e(th_a,e_a,m_dp,c.m_rho))[0]
                    if theta_omega != 0:
                        yi += f.f_kde(self.kde_omega,theta_omega*self.mult_omega,e_omega,self.theta_meson_e(th_a,e_a,m_dp,c.m_omega))[0]
                    if theta_phi != 0:
                        yi += f.f_kde(self.kde_phi,theta_phi*self.mult_phi,e_phi,self.theta_meson_e(th_a,e_a,m_dp,c.m_phi))[0]
                    data_sublist.append([th_a, e_a, m_dp, yi])
 
                data_list.append([data_sublist])
            else:
                data_sublist = []
                for th_a in self.th_a_list:
                    data_sublist.append([th_a, e_a, m_dp, 0])

                data_list.append([data_sublist])

        return data_list
    
    def load_ratio_rho(self):
        file_name = path.dirname(path.realpath(__file__))+'/../../widths/vector/mixing/R_rho.dat'
        self._m_min_rho = 0
        self._m_max_rho = 0
        self._r_max_rho = 0
        if path.exists(file_name):
            ratio_digitized = np.loadtxt(file_name)
            self._m_min_rho = ratio_digitized[0,0]
            self._m_max_rho = ratio_digitized[-1,0]
            self._r_max_rho = ratio_digitized[-1,1] # use the last value for m>m_max
            nlines = int(ratio_digitized.size/2)
            m_list = np.array([ratio_digitized[i,0] for i in range(nlines)])
            ratio_list = np.array([ratio_digitized[i,1] for i in range(nlines)])
            self._ratio_inter_rho = interp1d(m_list, ratio_list,fill_value="extrapolate")

        else: 
            print('[Warning:] \t',file_name,'not found')

    def get_mixing_rho(self,m_dp):
        if m_dp > self._m_max_rho:
            return self._r_max_rho*self.ratio_rho
        elif m_dp < self._m_min_rho:
            return 0.
        else:
            return self._ratio_inter_rho(m_dp)*self.ratio_rho
        
    def load_ratio_omega(self):
        file_name = path.dirname(path.realpath(__file__))+'/../../widths/vector/mixing/R_omega.dat'
        self._m_min_omega = 0
        self._m_max_omega = 0
        self._r_max_omega = 0
        if path.exists(file_name):
            ratio_digitized = np.loadtxt(file_name)
            self._m_min_omega = ratio_digitized[0,0]
            self._m_max_omega = ratio_digitized[-1,0]
            self._r_max_omega = ratio_digitized[-1,1] # use the last value for m>m_max
            nlines = int(ratio_digitized.size/2)
            m_list = np.array([ratio_digitized[i,0] for i in range(nlines)])
            ratio_list = np.array([ratio_digitized[i,1] for i in range(nlines)])
            self._ratio_inter_omega = interp1d(m_list, ratio_list,fill_value="extrapolate")

        else: 
            print('[Warning:] \t',file_name,'not found')

    def get_mixing_omega(self,m_dp):
        if m_dp > self._m_max_omega:
            return self._r_max_omega*self.ratio_omega
        elif m_dp < self._m_min_omega:
            return 0.
        else:
            return self._ratio_inter_omega(m_dp)*self.ratio_omega
        
    def load_ratio_phi(self):
        file_name = path.dirname(path.realpath(__file__))+'/../../widths/vector/mixing/R_phi.dat'
        self._m_min_phi = 0
        self._m_max_phi = 0
        self._r_max_phi = 0
        if path.exists(file_name):
            ratio_digitized = np.loadtxt(file_name)
            self._m_min_phi = ratio_digitized[0,0]
            self._m_max_phi = ratio_digitized[-1,0]
            self._r_max_phi = ratio_digitized[-1,1] # use the last value for m>m_max
            nlines = int(ratio_digitized.size/2)
            m_list = np.array([ratio_digitized[i,0] for i in range(nlines)])
            ratio_list = np.array([ratio_digitized[i,1] for i in range(nlines)])
            self._ratio_inter_phi = interp1d(m_list, ratio_list,fill_value="extrapolate")

        else: 
            print('[Warning:] \t',file_name,'not found')

    def get_mixing_phi(self,m_dp):
        if m_dp > self._m_max_phi:
            return self._r_max_phi*self.ratio_phi
        elif m_dp < self._m_min_phi:
            return 0.
        else:
            return self._ratio_inter_phi(m_dp)*self.ratio_phi
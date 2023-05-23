# modes/mixing.py

import os
import numpy as np
from general import exotic_production_setup as setup
from general import exotic_constants as c
from general import exotic_functions as f
from sklearn.neighbors import KernelDensity
# from scipy import integrate
# import matplotlib.pyplot as plt
from multiprocessing import Pool
from tqdm import tqdm
import sys

class mixing_production:
    def __init__(self, experiment, n_events,use_external=False):
        self.exp = experiment
        self.nevts = n_events
        self.beam_beta = f.beta_lor(setup.p_beam[self.exp],np.sqrt(setup.p_beam[self.exp]**2+c.m_p**2))
        self.beam_gamma = f.gamma_lor(self.beam_beta)

        self.th_a_list = [ th_a for th_a in np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins)]
        self.en_a_list = [ en_a for en_a in np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)),  num=setup.energy_bins)]
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append([ ma for ma in np.linspace(setup.mass_min[iFile], setup.mass_max[iFile], num=setup.mass_bins[iFile])])

        self.pi0_list = []
        self.eta_list = []
        self.etap_list = []

        if use_external:
            print("[Info:] \t Using external meson source")
            filename = os.getcwd()+"/tab_mesons/softQCD/softQCD_"+str(self.nevts//1000)+"kevts_pp_"+str(setup.p_beam[self.exp])+"GeV_8.2_pSet2_pom_ok.txt"
            if setup.p_beam[self.exp] == 30: #different file for KOTO
                filename = os.getcwd()+"/tab_mesons/softQCD/softQCD_200kevts_pp_30GeV_8.3_pSet2_pom_ok.txt"

            if os.path.exists(filename):
                external_list = np.loadtxt(filename)
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
        pythia.readString("SoftQCD:all = on")
        pythia.readString("pdf:pHardSet = 2")
        pythia.readString("Next:numberShowEvent = 5")
        pythia.readString("111:mayDecay = no")
        pythia.readString("Beams:idA = 2212")
        pythia.readString("Beams:idB = 2212")
        pythia.readString("Beams:eA = " + str(setup.p_beam[self.exp]))
        pythia.readString("Beams:eB = 0")
        pythia.readString("Beams:frameType = 2")

        pythia.init() # initialize

        for iEvent in range(self.nevts):
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
        return self.beam_gamma*e_cm + self.beam_gamma * self.beam_beta * np.sqrt(e_cm**2-m**2)

    def e_cm(self,e_lab,m):
        if e_lab<m:
            return 0
        return self.beam_gamma*e_lab - self.beam_gamma * self.beam_beta * np.sqrt(e_lab**2-m**2)

    def p_lab(self,p_cm,m):
        p_out = self.beam_gamma*p_cm + self.beam_gamma * self.beam_beta * np.sqrt(p_cm**2+m**2)
        if np.isnan(p_out):
            print("[Warning:] \t Lab. momentum estimation failed:", p_out,"setting 0. CM momentum:",p_cm, "mass:",m,"gamma:", self.beam_gamma,"beta:",self.beam_beta)
            p_out = 0
        return p_out

    def p_cm(self,p_lab,m):
        p_out = self.beam_gamma*p_lab - self.beam_gamma * self.beam_beta * np.sqrt(p_lab**2+m**2)
        if np.isnan(p_out):
            print("[Warning:] \t CM momentum estimation failed:", p_out,"setting 0. Lab momentum:",p_lab, "mass:",m,"gamma:", self.beam_gamma,"beta:",self.beam_beta)
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
    
    def process_pool(self,nthreads):

        for iFile in range(len(setup.mass_bins)):
            iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")

            with Pool(processes=nthreads) as pool:
                list_pi0, list_eta, list_etap = zip(*tqdm(pool.imap(self.ALP_mixing_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

            pool.close()
            pool.join() 
            
            # export
        
            f.export(setup.experiments[self.exp],"alp_mixingPi0_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_pi0,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)))
            f.export(setup.experiments[self.exp],"alp_mixingEta_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_eta,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)))
            f.export(setup.experiments[self.exp],"alp_mixingEtaPrim_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_etap,(len(self.th_a_list)*len(self.en_a_list)*len(self.m_lists[iFile]),4)))

        return

    def ALP_mixing_process(self,m_a):

        data_list_pi0 = []
        data_list_eta = []
        data_list_etap = []
        for e_a in self.en_a_list:
            if e_a > m_a:
                e_pi0 = np.sqrt(self.p_lab(self.p_cm(np.sqrt(e_a**2 - m_a**2), m_a), c.m_pi0)**2 + c.m_pi0**2)
                e_eta = np.sqrt(self.p_lab(self.p_cm(np.sqrt(e_a**2 - m_a**2), m_a), c.m_eta)**2 + c.m_eta**2)
                e_etap = np.sqrt(self.p_lab(self.p_cm(np.sqrt(e_a**2 - m_a**2), m_a), c.m_etap)**2 + c.m_etap**2)

                data_sublist_pi0 = []
                data_sublist_eta = []
                data_sublist_etap = []

                for th_a in self.th_a_list:
                    data_sublist_pi0.append([th_a, e_a, m_a, f.f_kde(self.kde_pi0,setup.reference_couplings["mixing"]**setup.scaling_exponent["mixing"]*self.mult_pi0,e_pi0,self.theta_meson(th_a,e_a,m_a,c.m_pi0))[0]])
                    data_sublist_eta.append([th_a, e_a, m_a, f.f_kde(self.kde_eta,setup.reference_couplings["mixing"]**setup.scaling_exponent["mixing"]*self.mult_eta,e_eta,self.theta_meson(th_a,e_a,m_a,c.m_eta))[0]])
                    data_sublist_etap.append([th_a, e_a, m_a, f.f_kde(self.kde_etap,setup.reference_couplings["mixing"]**setup.scaling_exponent["mixing"]*self.mult_etap,e_etap,self.theta_meson(th_a,e_a,m_a,c.m_etap))[0]])

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
        self.beam_beta = f.beta_lor(setup.p_beam[self.exp],np.sqrt(setup.p_beam[self.exp]**2+c.m_p**2))
        self.beam_gamma = f.gamma_lor(self.beam_beta)

        self.th_a_list = [ th_a for th_a in np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins)]
        self.en_a_list = [ en_a for en_a in np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)),  num=setup.energy_bins)]
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append([ ma for ma in np.linspace(setup.mass_min[iFile], setup.mass_max[iFile], num=setup.mass_bins[iFile])])

    # def mute(self):
    #     sys.stdout = open(os.devnull, 'w') 

    def e_lab(self,e_cm,m):
        if e_cm<m:
            return 0
        return self.beam_gamma*e_cm + self.beam_gamma * self.beam_beta * np.sqrt(e_cm**2-m**2)

    def e_cm(self,e_lab,m):
        if e_lab<m:
            return 0
        return self.beam_gamma*e_lab - self.beam_gamma * self.beam_beta * np.sqrt(e_lab**2-m**2)

    def p_lab(self,p_cm,m):
        p_out = self.beam_gamma*p_cm + self.beam_gamma * self.beam_beta * np.sqrt(p_cm**2+m**2)
        if np.isnan(p_out):
            print("[Warning:] \t Lab. momentum estimation failed:", p_out,"setting 0. CM momentum:",p_cm, "mass:",m,"gamma:", self.beam_gamma,"beta:",self.beam_beta)
            p_out = 0
        return p_out

    def p_cm(self,p_lab,m):
        p_out = self.beam_gamma*p_lab - self.beam_gamma * self.beam_beta * np.sqrt(p_lab**2+m**2)
        if np.isnan(p_out):
            print("[Warning:] \t CM momentum estimation failed:", p_out,"setting 0. Lab momentum:",p_lab, "mass:",m,"gamma:", self.beam_gamma,"beta:",self.beam_beta)
            p_out = 0
        return p_out

    def process_pool(self,nthreads):

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
                data_sublist_pi0.append([th, e_a, m_a, f.f_kde(kde_pi0,setup.reference_couplings["mixing"]**setup.scaling_exponent["mixing"]*mult_pi0,e_a,th)[0]])
                data_sublist_eta.append([th, e_a, m_a, f.f_kde(kde_eta,setup.reference_couplings["mixing"]**setup.scaling_exponent["mixing"]*mult_eta,e_a,th)[0]])
                data_sublist_etap.append([th, e_a, m_a, f.f_kde(kde_etap,setup.reference_couplings["mixing"]**setup.scaling_exponent["mixing"]*mult_etap,e_a,th)[0]])
            data_list_pi0.append([data_sublist_pi0])
            data_list_eta.append([data_sublist_eta])
            data_list_etap.append([data_sublist_etap])

        return data_list_pi0, data_list_eta, data_list_etap
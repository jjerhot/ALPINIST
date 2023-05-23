# modes/meson_decay.py

import os
import numpy as np
from general import exotic_production_setup as setup
from general import exotic_constants as c
from general import exotic_functions as f
from sklearn.neighbors import KernelDensity
from multiprocessing import Pool
from tqdm import tqdm
import sys

class bmeson_decay_production:
    def __init__(self, experiment, n_production,n_decay,daugther_name,use_external=False):
        self.exp = experiment
        self.nprod = n_production
        self.ndecays = n_decay
        self.daughter = daugther_name

        self.th_list = [ th for th in np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins)]
        self.en_list = [ en for en in np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)), num=setup.energy_bins)]
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append([ ma for ma in np.linspace(setup.mass_min[iFile], setup.mass_max[iFile], num=setup.mass_bins[iFile])])

        self.b_meson_list = []

        if use_external:
            print("[Info:] \t Using external B meson source")
            filename = os.getcwd()+"/tab_mesons/beauty/beauty_"+str(self.nprod//1000)+"kEvts_pp_8.2_"+str(setup.p_beam[self.exp])+"GeV_ptHat300MeV_new_ok.txt"

            if os.path.exists(filename):
                external_b_list = np.loadtxt(filename)
                if setup.p_beam[self.exp] == 400: #drop last column for 400GeV
                    external_b_list = np.delete(external_b_list, 6, 1)
                external_b_list = np.delete(external_b_list, 0, 1)
                self.b_meson_list = np.delete(external_b_list, 0, 1)
            else:
                print("[Error:] \t Path:",filename,"does not exist. Exiting")
                sys.exit(1)

        else:
            print("[Info:] \t Starting Pythia")
            self._init_pythia_()

            print("[Info:] \t Pythia production finished")

        self.b_meson_mult = len(self.b_meson_list)/self.nprod
        print("[Info:] \t B meson multiplicity",self.b_meson_mult)

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
        # pythia.readString("HardQCD:gg2bbbar = on")
        # pythia.readString("HardQCD:qqbar2bbbar = on")
        pythia.readString("HardQCD:hardbbbar = on")
        pythia.readString("Next:numberShowEvent = 5")
        pythia.readString("Beams:idA = 2212")
        pythia.readString("Beams:idB = 2212")
        pythia.readString("Beams:eA = " + str(setup.p_beam[self.exp]))
        pythia.readString("Beams:eB = 0")
        pythia.readString("Beams:frameType = 2")

        pythia.init() # initialize

        for iEvent in range(self.nprod):
            pythia.next()
            for isubEvent in range(pythia.event.size()):
                evtID = np.abs(pythia.event[isubEvent].id())
                if evtID == 511 or evtID == 521:#b meson
                    self.b_meson_list.append([evtID,pythia.event[isubEvent].px(),pythia.event[isubEvent].py(),pythia.event[isubEvent].pz()])

        return

    def b_to_k(self,m_exo):
        en_exo = 0
        th_exo = 0

        data_en_th_exo = []
        for event in self.b_meson_list:
            eid = np.abs(event[0])
            m_parent = 0
            m_daughter = 0
            if eid == 511: #B0
                m_parent = c.m_B0
                m_daughter = c.m_K0
            elif eid == 521: #B
                m_parent = c.m_B
                m_daughter = c.m_K
            else:
                continue
            if m_parent < (m_exo+m_daughter):
                continue
            
            p_parent = [event[1],event[2],event[3]]
            e_parent = np.sqrt(m_parent**2+np.linalg.norm(p_parent)**2)
            p_exo_cm = np.sqrt(f.lambda_Kallen(m_parent,m_daughter,m_exo))/(2*m_parent)
            lor = f.lorentz_transf(e_parent/m_parent,p_parent/e_parent)

            for iDecay in range(self.ndecays):
                en_exo,th_exo = f.decay_meson(p_exo_cm,m_exo,lor)
                if en_exo <= 0: continue
                if th_exo <= 0: th_exo = 1E-10 #set min value for log
                data_en_th_exo.append([np.log(en_exo), np.log(th_exo)])

        if len(data_en_th_exo) > 0:
            self.kde_b_to_k = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(data_en_th_exo))
            return True
        
        return False

    def b_to_k_star(self,m_exo):
        en_exo = 0
        th_exo = 0

        data_en_th_exo = []
        for event in self.b_meson_list:
            eid = np.abs(event[0])
            m_parent = 0
            m_daughter = 0
            if eid == 511: #B0
                m_parent = c.m_B0
                m_daughter = c.m_K0star
            elif eid == 521: #B
                m_parent = c.m_B
                m_daughter = c.m_Kstar
            else:
                continue
            if m_parent < (m_exo+m_daughter):
                continue
            
            p_parent = [event[1],event[2],event[3]]
            e_parent = np.sqrt(m_parent**2+np.linalg.norm(p_parent)**2)
            p_exo_cm = np.sqrt(f.lambda_Kallen(m_parent,m_daughter,m_exo))/(2*m_parent)
            lor = f.lorentz_transf(e_parent/m_parent,p_parent/e_parent)

            for iDecay in range(self.ndecays):
                en_exo,th_exo = f.decay_meson(p_exo_cm,m_exo,lor)
                if en_exo <= 0: continue
                if th_exo <= 0: th_exo = 1E-10 #set min value for log
                data_en_th_exo.append([np.log(en_exo), np.log(th_exo)])

        if len(data_en_th_exo) > 0:
            self.kde_b_to_kstar = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(data_en_th_exo))
            return True
        
        return False

    def process_pool(self,nthreads):

        for iFile in range(len(setup.mass_bins)):
            iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")

            with Pool(processes=nthreads) as pool:
                list_bk, list_bkstar = zip(*tqdm(pool.imap(self.b_meson_decay_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

            pool.close()
            pool.join()
            
            # export
        
            f.export(setup.experiments[self.exp],self.daughter + "_BmesonK_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_bk,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)))
            f.export(setup.experiments[self.exp],self.daughter + "_BmesonKstar_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_bkstar,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)))

        return
  
    def b_meson_decay_process(self,m_exo):

        success_bk = self.b_to_k(m_exo)
        success_bkstar = self.b_to_k_star(m_exo)

        data_list_bk_exo = []
        data_list_bkstar_exo = []

        for en_exo in self.en_list:
            data_sublist_bk = []
            data_sublist_bkstar = []
            if en_exo > m_exo:
                if success_bk:
                    for th_exo in self.th_list:
                        data_sublist_bk.append([th_exo, en_exo, m_exo, f.f_kde(self.kde_b_to_k,setup.reference_couplings["Bmeson"]**setup.scaling_exponent["Bmeson"]*self.b_meson_mult,en_exo,th_exo)[0]])
                else:
                    for th_exo in self.th_list:
                        data_sublist_bk.append([th_exo, en_exo, m_exo, 0])
                if success_bkstar:
                    for th_exo in self.th_list:
                        data_sublist_bkstar.append([th_exo, en_exo, m_exo, f.f_kde(self.kde_b_to_kstar,setup.reference_couplings["Bmeson"]**setup.scaling_exponent["Bmeson"]*self.b_meson_mult,en_exo,th_exo)[0]])
                else:
                    for th_exo in self.th_list:
                        data_sublist_bkstar.append([th_exo, en_exo, m_exo, 0])
            else:
                for th_exo in self.th_list:
                    data_sublist_bk.append([th_exo, en_exo, m_exo, 0])
                    data_sublist_bkstar.append([th_exo, en_exo, m_exo, 0])

            data_list_bk_exo.append([data_sublist_bk])
            data_list_bkstar_exo.append([data_sublist_bkstar])

        return data_list_bk_exo, data_list_bkstar_exo
    
class dmeson_decay_production:
    def __init__(self, experiment, n_production,n_decay,daugther_name,use_external=False):
        self.exp = experiment
        self.nprod = n_production
        self.ndecays = n_decay
        self.daughter = daugther_name

        self.th_list = [ th for th in np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins)]
        self.en_list = [ en for en in np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)), num=setup.energy_bins)]
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append([ ma for ma in np.linspace(setup.mass_min[iFile], setup.mass_max[iFile], num=setup.mass_bins[iFile])])

        self.d_meson_list = []

        if use_external:
            print("[Info:] \t Using external D meson source")
            filename = os.getcwd()+"/tab_mesons/charm/hsccbar_"+str(setup.p_beam[self.exp])+"GeV_ok.txt"

            if os.path.exists(filename):
                external_d_list = np.loadtxt(filename)
                if setup.p_beam[self.exp] == 400: #drop last column for 400GeV
                    external_d_list = np.delete(external_d_list, 6, 1)
                external_d_list = np.delete(external_d_list, 0, 1)
                self.d_meson_list = np.delete(external_d_list, 0, 1)
            else:
                print("[Error:] \t Path:",filename,"does not exist. Exiting")
                sys.exit(1)

        else:
            print("[Info:] \t Starting Pythia")
            self._init_pythia_()

            print("[Info:] \t Pythia production finished")

        self.d_meson_mult = len(self.d_meson_list)/self.nprod
        print("[Info:] \t D meson multiplicity",self.d_meson_mult)

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
        pythia.readString("HardQCD:hardccbar = on")
        pythia.readString("Next:numberShowEvent = 5")
        pythia.readString("Beams:idA = 2212")
        pythia.readString("Beams:idB = 2212")
        pythia.readString("Beams:eA = " + str(setup.p_beam[self.exp]))
        pythia.readString("Beams:eB = 0")
        pythia.readString("Beams:frameType = 2")

        pythia.init() # initialize

        for iEvent in range(self.nprod):
            pythia.next()
            for isubEvent in range(pythia.event.size()):
                evtID = np.abs(pythia.event[isubEvent].id())
                if evtID == 411 or evtID == 421:#d meson
                    self.d_meson_list.append([evtID,pythia.event[isubEvent].px(),pythia.event[isubEvent].py(),pythia.event[isubEvent].pz()])

        return

    def d_to_pi(self,m_exo):
        en_exo = 0
        th_exo = 0

        data_en_th_exo = []
        for event in self.d_meson_list:
            eid = np.abs(event[0])
            m_parent = 0
            m_daughter = 0
            if eid == 411: #D
                m_parent = c.m_D
                m_daughter = c.m_pi
            elif eid == 421: #D0
                m_parent = c.m_D0
                m_daughter = c.m_pi0
            else:
                continue
            if m_parent < (m_exo+m_daughter):
                continue
            
            p_parent = [event[1],event[2],event[3]]
            e_parent = np.sqrt(m_parent**2+np.linalg.norm(p_parent)**2)
            p_exo_cm = np.sqrt(f.lambda_Kallen(m_parent,m_daughter,m_exo))/(2*m_parent)
            lor = f.lorentz_transf(e_parent/m_parent,p_parent/e_parent)

            for iDecay in range(self.ndecays):
                en_exo,th_exo = f.decay_meson(p_exo_cm,m_exo,lor)
                if en_exo <= 0: continue
                if th_exo <= 0: th_exo = 1E-10 #set min value for log
                data_en_th_exo.append([np.log(en_exo), np.log(th_exo)])

        if len(data_en_th_exo) > 0:
            self.kde_d_to_pi = KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(data_en_th_exo))
            return True

        return False

    def process_pool(self,nthreads):

        for iFile in range(len(setup.mass_bins)):
            iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")

            with Pool(processes=nthreads) as pool:
                list_dpi = list(tqdm(pool.imap(self.d_meson_decay_process,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

            pool.close()
            pool.join()
            
            # export
        
            f.export(setup.experiments[self.exp],self.daughter + "_DmesonPi_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_dpi,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)))

        return
  
    def d_meson_decay_process(self,m_exo):

        success = self.d_to_pi(m_exo)

        data_list_dpi_exo = []

        for en_exo in self.en_list:
            if en_exo > m_exo and success:
                data_sublist_dpi = []

                for th_exo in self.th_list:
                    data_sublist_dpi.append([th_exo, en_exo, m_exo, f.f_kde(self.kde_d_to_pi,setup.reference_couplings["Dmeson"]**setup.scaling_exponent["Dmeson"]*self.d_meson_mult,en_exo,th_exo)[0]])

                data_list_dpi_exo.append([data_sublist_dpi])
            else:
                data_sublist_dpi = []

                for th_exo in self.th_list:
                    data_sublist_dpi.append([th_exo, en_exo, m_exo, 0])

                data_list_dpi_exo.append([data_sublist_dpi])

        return data_list_dpi_exo
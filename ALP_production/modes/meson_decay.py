# modes/meson_decay.py

import os
import numpy as np
from ALP_production.general import exotic_production_setup as setup
from ALP_production.general import exotic_constants as c
from ALP_production.general import exotic_functions as f
from sklearn.neighbors import KernelDensity
from multiprocessing import Pool
from ROOT import TGenPhaseSpace, TLorentzVector
from tqdm import tqdm
import sys

rng = np.random.default_rng()

class decay_to_2_body:
    """Generic class to decay parent (m_parent) in a two body decay (m_exo, m_daughter)
    """
    def __init__(self, parent_frame_lorentz_boosts,  m_exo, m_parent, m_daughter):
        self.data_en_th_exo = self.exo_en_th_array(parent_frame_lorentz_boosts, m_exo, m_parent, m_daughter)
        self.success = 0<len(self.data_en_th_exo)

    def __add__(self, other):
        if self.success and other.success: self.data_en_th_exo = np.append(self.data_en_th_exo, other.data_en_th_exo, 0)
        elif other.success: self.data_en_th_exo = other.data_en_th_exo
        self.success = self.success or other.success
        return self

    def kde(self):
        return KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.data_en_th_exo))

    @staticmethod
    def exo_en_th_array(parent_frame_lorentz_boosts, m_exo, m_parent, m_daughter):  
        p_exo_cm = np.sqrt(f.lambda_Kallen(m_parent,m_daughter,m_exo))/(2*m_parent)
        data_en_th_exo = f.decay_meson_vectorized(p_exo_cm, m_exo, parent_frame_lorentz_boosts) 
        if not data_en_th_exo.size: return np.zeros((2,1))
        data_en_th_clean = data_en_th_exo[np.logical_and(np.all(np.isfinite(data_en_th_exo),axis=1), 0<=data_en_th_exo[:,0])]
        data_en_th_clean[:,1] = np.clip(data_en_th_clean[:,1],1E-10,None) #set min value for log theta
        return np.log(data_en_th_clean).T
    
    @staticmethod
    def secondary_generator(parent_boosts, m_parent, m_daughter1=c.m_tau, m_daughter2=0, daughter_pid = 15):
        p_daughter1_cm = np.sqrt(f.lambda_Kallen(m_parent,m_daughter1,m_daughter2))/(2*m_parent)
        daughter_momenta  = f.decay_meson_cartesian_momentum_vectorized(p_daughter1_cm, m_daughter1, parent_boosts)
        daughter_event_list = np.append(daughter_pid*np.ones((np.shape(daughter_momenta)[0],1)), daughter_momenta,axis=1)
        return daughter_event_list

class decay_to_n_body: 
    """Generic class to decay parent (m_parent) in a n body decay (m_exo, *m_daughters) using (reweighted) TGenPhaseSpace
    """
    def __init__(self, parent_frame_lorentz_boosts, m_exo, m_parent, m_daughters, n_evts_decay=10000, n_exos=1, weight_function = None):
        exo_weights, exo_parent_frame_momenta = self.generate_restframe_phasespace(m_parent, n_exos*[m_exo]+m_daughters, n_exos = n_exos, n_evts = n_evts_decay, weight_function=weight_function)
        self.data_en_th_exo = self.exo_en_th_array(parent_frame_lorentz_boosts, exo_weights, exo_parent_frame_momenta)
        self.success = 0<len(self.data_en_th_exo)

    def __add__(self, other):
        if self.success and other.success: self.data_en_th_exo = np.append(self.data_en_th_exo, other.data_en_th_exo, 0)
        elif other.success: self.data_en_th_exo = other.data_en_th_exo
        self.success = self.success or other.success
        return self
    
    def kde(self):
        return KernelDensity(kernel='gaussian', bandwidth="silverman").fit(np.array(self.data_en_th_exo))
    
    @staticmethod
    def exo_en_th_array(parent_frame_lorentz_boosts, exo_weights, exo_parent_frame_momenta):
        p_exo_cms = rng.choice(exo_parent_frame_momenta, size = parent_frame_lorentz_boosts.shape[0], p = exo_weights)
        data_en_th_exo = f.decay_meson_vectorized(p_exo_cms, np.nan, parent_frame_lorentz_boosts) 
        data_en_th_clean = data_en_th_exo[np.logical_and(np.all(np.isfinite(data_en_th_exo),axis=1), 0<=data_en_th_exo[:,0])]
        data_en_th_clean[:,1] = np.clip(data_en_th_clean[:,1], 1E-10, None) #set min value for log theta
        return np.log(data_en_th_clean).T

    # 3-body, first daughter is also exotic, weight function returns differential x-sec given m12^2, m23^2
    @staticmethod
    def generate_restframe_phasespace(m_parent, m_daughters, n_exos=1, n_evts = 100000, weight_function=None):
        """
        Generate a 4 momentum array with associated weights based on TGenPhaseSpace
        Args:
            m_parent (float):                   mass of parent particle
            m_daughters (list(float)):          list of masses starting with the exotics 
            n_exos (int, optional):             number of exotics among the daughters. Defaults to 1.
            n_evts (int, optional):             number of events to generate with TGenPhaseSpace. Defaults to 100000.
            weight_function (func, optional):   optional weight function with params
                                                    (m_X, M12^2, M23^2)         for 3 body decays
                                                    (list(daughter_4_momenta))  else.
                                                Defaults to None.
        Returns:
            weights:    resulting weight for event
            p4Exos:     4 momenta of exotics generated in the event
        """        
        events = TGenPhaseSpace()
        parent_at_rest = TLorentzVector(0., 0., 0., m_parent)
        decay_masses = np.array(m_daughters)
        events.SetDecay(parent_at_rest, decay_masses.size, decay_masses)
        weights, pExos = [[],[]]
        for _ in range(n_evts):
            weight = events.Generate()
            daughter_Vectors = []
            for i in range(decay_masses.size): daughter_Vectors.append(events.GetDecay(i))
            if weight_function != None:
                if decay_masses.size == 3:
                    m12_2 = decay_masses[0]**2+decay_masses[1]**2+2.*daughter_Vectors[0].Dot(daughter_Vectors[1]) # m12^2 = (p1+p2)^2
                    m23_2 = decay_masses[1]**2+decay_masses[2]**2+2.*daughter_Vectors[1].Dot(daughter_Vectors[2]) # m23^2 = (p2+p3)^2
                    weight_external = weight_function(decay_masses[0], m12_2, m23_2) #trust tgps to generate m12, m23 in correct bounds
                else: weight_external = weight_function(daughter_Vectors)
                weight *= weight_external
            for i in range(n_exos):
                weights.append(weight)
                pExos.append([daughter_Vectors[i].E(), daughter_Vectors[i].Px(), daughter_Vectors[i].Py(), daughter_Vectors[i].Pz()])
        return np.array(weights)/np.sum(weights), np.array(pExos)

class bmeson_decay_production:
    def __init__(self, experiment, n_production, n_decay, daugther_name, use_external=False, active_coupling="", use_3part_prod = False, use_massless = False, use_empirical = "", delete_intermed_lists=True, target_code = "2212", beam_particle_code = "2212", pythia_Ppdf_code=20, pythia_extra_args=[]):
        self.exp = experiment
        self.nprod = n_production
        self.ndecays = n_decay
        self.daughter_exo = daugther_name 
        self.active_coupling = active_coupling
        self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced in Bmeson decays | generated with statistics of {n_production*n_decay} events using "

        p_beam = setup.p_beam.get(self.exp,400)
        self.th_list = np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins).tolist()
        self.en_list = np.linspace(setup.energy_min.get(p_beam), setup.energy_max.get(p_beam), num=setup.energy_bins).tolist()
        self.norm_xsec = c.sigma_bb[p_beam] / c.sigma_pp[p_beam] * c.A_target[experiment]**(1./3)

        self.th_list = np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins).tolist()
        self.en_list = np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)), num=setup.energy_bins).tolist()
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())

        self.b_meson_list = []
        self.b0_meson_list = []
        self.bs_meson_list = []
        if use_empirical == "SMTs":
            sys.exit("[Error:] \t SHiP meson tables as generated according to SHiP-NOTE-2015-009 are not publically available ")
            # tablepath = './../../SHiPGeneratedFiles/'
            # self.b_meson_list       = rng.choice(np.loadtxt(tablepath+ '521.txt'),size = int(self.nprod * 0.417))
            # self.b0_meson_list      = rng.choice(np.loadtxt(tablepath+ '511.txt'),size = int(self.nprod * 0.418))
            # self.bbar_meson_list    = rng.choice(np.loadtxt(tablepath+'-521.txt'),size = int(self.nprod * 0.417))
            # self.b0bar_meson_list   = rng.choice(np.loadtxt(tablepath+'-511.txt'),size = int(self.nprod * 0.418))
            # self.bs_meson_list      = np.array([[531,0,0,1],[531,0,0,1]])
            # self.bsbar_meson_list   = np.array([[-531,0,0,1],[-531,0,0,1]])
        elif use_external:
            print("[Info:] \t Using external B meson source")
            filename = os.path.dirname(os.path.realpath(__file__))+"/../tab_mesons/beauty/"
            filename += "hsbbbar_"+("massless_" if use_massless else "")+("statcomb23partprod_" if use_3part_prod else "")+str(setup.p_beam[self.exp])+'GeV_'+f.number_to_3sigfigs_str(n_production)+".dat"
            if target_code != "2212": 
                PID_to_extension =  {"2112":"pn", "100822080":"pPb"}
                assert target_code in PID_to_extension.keys(), "[Error:] \t target code " + target_code + " not recognised."
                filename = filename.replace(".dat", "_" + PID_to_extension[target_code ] + ".dat")

            if os.path.exists(filename):
                external_b_list = np.loadtxt(filename,skiprows=1)
                with open(filename) as meson_source_file: self.export_header += meson_source_file.readline()[1:-1]
                self.b_meson_list       = external_b_list[np.where(external_b_list[:,0]== 521)]
                self.bbar_meson_list    = external_b_list[np.where(external_b_list[:,0]==-521)]
                self.b0_meson_list      = external_b_list[np.where(external_b_list[:,0]== 511)]
                self.b0bar_meson_list   = external_b_list[np.where(external_b_list[:,0]==-511)]
                self.bs_meson_list      = external_b_list[np.where(external_b_list[:,0]== 531)]
                self.bsbar_meson_list   = external_b_list[np.where(external_b_list[:,0]==-531)]
            else:
                print("[Error:] \t Path:",filename,"does not exist. Exiting")
                sys.exit(1)
        else:
            print("[Info:] \t Starting Pythia")
            self._init_pythia_(beam_particle_code, target_code, pythia_Ppdf_code, pythia_extra_args)

            print("[Info:] \t Pythia production finished")

        self.b_meson_mult = (len(self.b_meson_list)+len(self.bbar_meson_list))/self.nprod
        print("[Info:] \t B meson multiplicity",self.b_meson_mult)
        self.b0_meson_mult = (len(self.b0_meson_list)+len(self.b0bar_meson_list))/self.nprod
        print("[Info:] \t B0 meson multiplicity",self.b0_meson_mult)
        self.bs_meson_mult = (len(self.bs_meson_list)+len(self.bsbar_meson_list))/self.nprod
        print("[Info:] \t Bs meson multiplicity",self.bs_meson_mult)

        self.multiplicities = {
            "Bmeson"  : self.b_meson_mult,
            "B0meson" : self.b0_meson_mult,
            "Bsmeson" : self.bs_meson_mult
        }

        self.b0_boosts = f.lorentz_boosts_vectorized(event_list=np.concatenate((self.b0_meson_list,self.b0bar_meson_list)), mass=c.m_B0, ndecays=self.ndecays)
        self.b_boosts =  f.lorentz_boosts_vectorized(event_list=np.concatenate((self.b_meson_list,self.bbar_meson_list)),  mass=c.m_B,  ndecays=self.ndecays)
        if self.daughter_exo =="hnl" and (self.bsbar_meson_list.size  or self.bs_meson_list.size): self.bs_boosts = f.lorentz_boosts_vectorized(event_list=np.concatenate((self.bs_meson_list,self.bsbar_meson_list)),  mass=c.m_Bs, ndecays=self.ndecays)
        if delete_intermed_lists: 
            del self.b_meson_list, self.bbar_meson_list, self.b0_meson_list, self.b0bar_meson_list
            if self.daughter_exo =="hnl": del self.bs_meson_list, self.bsbar_meson_list

    def _init_pythia_(self, beam_particle_code, target_code, pythia_Ppdf_code, pythia_extra_args, save_pythia_file=False):
        #import pythia
        if not os.path.exists(setup.pythia_dir):
            print("[Error:] \t Set correct path to Pythia in general/exotic_production_setup.py")
            sys.exit(1)
        pythia_cfg_dir = setup.pythia_dir+"/examples/"
        if not os.path.exists(pythia_cfg_dir):
            pythia_cfg_dir = setup.pythia_dir+"/share/Pythia8/examples/"
        if not os.path.exists(pythia_cfg_dir):
            print("[Error:] \t Set correct path to Pythia")
            sys.exit(1)

        # cfg = open(setup.pythia_dir+"/examples/Makefile.inc")
        cfg = open(pythia_cfg_dir+"Makefile.inc")
        lib = setup.pythia_dir+"/lib"
        for line in cfg:
            if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
        sys.path.insert(0, lib)
        import pythia8

        pythia = pythia8.Pythia("", False) 
        pythia_modifiers = ["HardQCD:hardbbbar = on",f"Beams:idA = {beam_particle_code}",f"Beams:idB = {target_code}","Beams:eA = " + str(setup.p_beam[self.exp]),"Beams:eB = 0","Beams:frameType = 2", f"PDF:pSet={pythia_Ppdf_code}"]#"PhaseSpace:mHatMin=0",
        if isinstance(pythia_extra_args,str): pythia_modifiers.append(pythia_extra_args)
        elif len(pythia_extra_args): pythia_modifiers += pythia_extra_args
        for modifier in pythia_modifiers+["Print:quiet=on","Print:errors=off","Next:numberShowEvent = 0"]: pythia.readString(modifier)
        self.export_header += f"Pythia{pythia.parm('Pythia:versionNumber')} with settings [" +" ".join(pythia_modifiers)+"]"

        pythia.init() # initialize

        evtID_to_list = {
            511: [],
            -511:[],
            521: [],
            -521:[],
            531: [],
            -531:[]
        }
        for _ in range(self.nprod):
            pythia.next()
            for isubEvent in range(pythia.event.size()):
                evt = pythia.event[isubEvent]
                if evt.id() not in evtID_to_list.keys(): continue
                evtID_to_list[evt.id()].append([evt.id(),evt.px(),evt.py(),evt.pz()])

        self.b0_meson_list = np.array(evtID_to_list[521])
        self.b_meson_list  = np.array(evtID_to_list[511])
        self.bs_meson_list = np.array(evtID_to_list[531])
        self.b0bar_meson_list = np.array(evtID_to_list[-521])
        self.bbar_meson_list  = np.array(evtID_to_list[-511])
        self.bsbar_meson_list = np.array(evtID_to_list[-531])
        if save_pythia_file:
            file_info  =  "Beauty meson momenta as generated for " + self.nprod +" impining protons, with PYTHIA (version "+str(pythia.parm("Pythia:versionNumber"))+") and modifiers " + ', '.join(pythia_modifiers)
            np.savetxt(f".ALP_production/tab_mesons/beauty/hsbbbar_{setup.p_beam[self.exp]}GeV_{f.number_to_3sigfigs_str(self.nprod)}.dat", np.concatenate(evtID_to_list.values()), header=file_info, fmt='%i %4f %4f %3f')
        return

    def process_pool(self, nthreads, active_couplings = [""], single_mass_point=0):
        """Wrapper to run production in parallel threads.
        Args:
            nthreads (int): number of parallel threads to run
            single_mass_point (float, optional): Single mass point for which to run. Defaults to 0.
        """

        for active_coupling in active_couplings:   
            export_units_header = "Theta[rad], E_x [GeV], m_x [GeV], dY[per(rad GeV " + {"hnl":f"U2_{active_coupling}", "alp":"g_bs_eff^2/GeV^2", "ds":"Y^2" }[self.daughter_exo] + ')]'

            mixing_info  = "-" + active_coupling + "Mixing" if active_coupling else ""
            if active_coupling: print("[Info:] \t Starting "+self.daughter_exo+" production with active coupling U2_" + active_coupling)
            export_info_header = self.export_header if not active_coupling else self.export_header.replace(' yield', f'({active_coupling} mixing) yield')
            self.active_coupling = active_coupling

            if single_mass_point>0:
                list_b = self.b_meson_decay_process(single_mass_point)
                f.export(setup.experiments[self.exp],self.daughter + "_Bmeson"+mixing_info+"_beam" + str(setup.p_beam[self.exp]) + "GeV_" + str(int(single_mass_point*1e6)) + "keV_" + setup.experiments[self.exp] + ".dat",np.reshape(list_b,(len(self.th_list)*len(self.en_list),4)),header = self.export_header+'\n'+export_units_header)
            else:
                for iFile in range(len(setup.mass_bins)):
                    iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
                    if setup.mass_min[iFile]*1000 < 1:
                        iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
                    print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," GeV with",nthreads,"threads")

                    with Pool(processes=nthreads) as pool:
                        list_b = list(tqdm(pool.imap(self.b_meson_decay_process, self.m_lists[iFile]), total=len(self.m_lists[iFile])))

                    pool.close()
                    pool.join()

                    # export   
                    f.export(setup.experiments[self.exp],self.daughter_exo + "_Bmeson"+mixing_info+"_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_b,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)),header = export_info_header+'\n'+export_units_header)
        return
  
    def b_meson_decay_process(self, m_exo):
        """evaluates open beauty meson decays into feebly interacting particle with mass m_exo emitted at angle theta to beam axis with energy e
        Args:
            m_exo (float): exotic mss

        Returns:
            array(float): [theta, e, m_exo, expected yield]
        """

        kdes = []
        weights = []
        if self.daughter_exo in ["alp", "ds"]:
            for fin_state in setup.kaons_list:
                #charged
                m_K = c.m_K[fin_state]
                if(m_exo+m_K<c.m_B):
                    b_decay = decay_to_2_body(self.b_boosts, m_exo, c.m_B, m_K)
                    if b_decay.success: 
                        kdes.append(b_decay.kde())
                        weights.append(self.multiplicities["Bmeson"]*f.get_branching_ratios(self.daughter_exo, m_exo, "Bmeson", fin_state, 1))
                    del b_decay
                
                #neutral
                if fin_state == "K" or fin_state == "Kstar_892" or fin_state == "K2star_1430":
                    m_K = c.m_K[fin_state+"_0"]
                else: m_K = c.m_K[fin_state]
                if(m_exo+m_K<c.m_B0):
                    b0_decay = decay_to_2_body(self.b0_boosts, m_exo, c.m_B0, m_K)
                    if b0_decay.success: 
                        kdes.append(b0_decay.kde())
                        weights.append(self.multiplicities["B0meson"]*f.get_branching_ratios(self.daughter_exo, m_exo, "B0meson", fin_state, 1))
                    del b0_decay
                
        elif self.daughter_exo == "hnl":
            if self.active_coupling == "El" and m_exo+c.m_el<c.m_B:
                b_decay = decay_to_2_body(self.b_boosts, m_exo, c.m_B, c.m_el)
                if b_decay.success: 
                    kdes.append(b_decay.kde())
                    weights.append(self.multiplicities["Bmeson"]*f.get_branching_ratios("hnl", m_exo, "Bmeson", "El", 1))
                del b_decay

                curr_weight = 0.
                successfull_decays = []
                if m_exo+c.m_el+c.m_pi<c.m_B:
                    b_pi0El_decay = decay_to_n_body(self.b_boosts, m_exo, c.m_B, [c.m_el, c.m_pi0], weight_function=f.hnl_production_BRs.get_differential_weight_function("Bmeson", ["El","Pi0"]))
                    if b_pi0El_decay.success: 
                        successfull_decays.append(b_pi0El_decay)
                        curr_weight += self.multiplicities["Bmeson"]*f.get_branching_ratios("hnl", m_exo, "Bmeson", ["El","Pi0"], 1)
                    del b_pi0El_decay
                if m_exo+c.m_el+c.m_pi<c.m_B0:
                    b0_piEl_decay = decay_to_n_body(self.b0_boosts, m_exo, c.m_B0, [c.m_el, c.m_pi], weight_function=f.hnl_production_BRs.get_differential_weight_function("B0meson", ["El","Pi"]))
                    if b0_piEl_decay.success:
                        successfull_decays.append(b0_piEl_decay)
                        curr_weight += self.multiplicities["B0meson"]*f.get_branching_ratios("hnl", m_exo, "B0meson", ["El","Pi"], 1)
                    del b0_piEl_decay
                if successfull_decays:
                    b_piEl_decay = np.sum(successfull_decays)
                    kdes.append(b_piEl_decay.kde())
                    weights.append(curr_weight)
                    del b_piEl_decay
                    
                curr_weight = 0.
                successfull_decays = []
                if m_exo+c.m_el+c.m_D0<c.m_B:
                    b_D0El_decay = decay_to_n_body(self.b_boosts, m_exo, c.m_B, [c.m_el, c.m_D0], weight_function=f.hnl_production_BRs.get_differential_weight_function("Bmeson", ["El","D0"]))
                    if b_D0El_decay.success: 
                        successfull_decays.append(b_D0El_decay)
                        curr_weight += self.multiplicities["Bmeson"]*f.get_branching_ratios("hnl", m_exo, "Bmeson", ["El","D0"], 1)
                    del b_D0El_decay
                if m_exo+c.m_el+c.m_D<c.m_B0:
                    b0_DEl_decay = decay_to_n_body(self.b0_boosts, m_exo, c.m_B0, [c.m_el, c.m_D], weight_function=f.hnl_production_BRs.get_differential_weight_function("B0meson", ["El","D"]))
                    if b0_DEl_decay.success:
                        successfull_decays.append(b0_DEl_decay)
                        curr_weight += self.multiplicities["B0meson"]*f.get_branching_ratios("hnl", m_exo, "B0meson", ["El","D"], 1)
                    del b0_DEl_decay
                if successfull_decays:
                    b_DEl_decay = np.sum(successfull_decays)
                    kdes.append(b_DEl_decay.kde())
                    weights.append(curr_weight)
                    del b_DEl_decay

                curr_weight = 0.
                successfull_decays = []
                if m_exo+c.m_el+c.m_D0star<c.m_B:
                    b_D0starEl_decay = decay_to_n_body(self.b_boosts, m_exo, c.m_B, [c.m_el, c.m_D0star], weight_function=f.hnl_production_BRs.get_differential_weight_function("Bmeson", ["El","D0star"]))
                    if b_D0starEl_decay.success: 
                        successfull_decays.append(b_D0starEl_decay)
                        curr_weight += self.multiplicities["Bmeson"]*f.get_branching_ratios("hnl", m_exo, "Bmeson", ["El","D0star"], 1)
                    del b_D0starEl_decay
                if m_exo+c.m_el+c.m_Dstar<c.m_B0:
                    b0_DstarEl_decay = decay_to_n_body(self.b0_boosts, m_exo, c.m_B0,[c.m_el, c.m_Dstar], weight_function=f.hnl_production_BRs.get_differential_weight_function("B0meson", ["El","Dstar"]))
                    if b0_DstarEl_decay.success: 
                        successfull_decays.append(b0_DstarEl_decay)
                        curr_weight += self.multiplicities["B0meson"]*f.get_branching_ratios("hnl", m_exo, "B0meson", ["El","Dstar"], 1)
                    del b0_DstarEl_decay
                if successfull_decays:
                    b_DstarEl_decay = np.sum(successfull_decays)
                    kdes.append(b_DstarEl_decay.kde())
                    weights.append(curr_weight)
                    del b_DstarEl_decay

                if m_exo+c.m_el+c.m_Ds<c.m_Bs:
                    bs_DsEl_decay = decay_to_n_body(self.bs_boosts, m_exo, c.m_Bs, [c.m_el, c.m_Ds], weight_function=f.hnl_production_BRs.get_differential_weight_function("Bsmeson", ["El","Ds"]))
                    if bs_DsEl_decay.success: 
                        kdes.append(bs_DsEl_decay.kde())
                        weights.append((self.multiplicities["Bsmeson"])*f.get_branching_ratios("hnl", m_exo, "Bsmeson", ["El","Ds"], 1))
                    del bs_DsEl_decay
                if m_exo+c.m_el+c.m_Dsstar<c.m_Bs:
                    bs_DsstarEl_decay = decay_to_n_body(self.bs_boosts, m_exo, c.m_Bs, [c.m_el, c.m_Dsstar], weight_function=f.hnl_production_BRs.get_differential_weight_function("Bsmeson", ["El","Dsstar"]))
                    if bs_DsstarEl_decay.success: 
                        kdes.append(bs_DsstarEl_decay.kde())
                        weights.append((self.multiplicities["Bsmeson"])*f.get_branching_ratios("hnl", m_exo, "Bsmeson", ["El","Dsstar"], 1))
                    del bs_DsstarEl_decay

            elif self.active_coupling == "Mu" and (m_exo+c.m_mu<c.m_B):
                b_decay = decay_to_2_body(self.b_boosts, m_exo, c.m_B, c.m_mu)
                if b_decay.success: 
                    kdes.append(b_decay.kde())
                    weights.append(self.multiplicities["Bmeson"]*f.get_branching_ratios("hnl", m_exo, "Bmeson", "Mu", 1))
                del b_decay

                curr_weight = 0.
                successfull_decays = []
                if m_exo+c.m_mu+c.m_pi<c.m_B:
                    b_pi0Mu_decay = decay_to_n_body(self.b_boosts, m_exo, c.m_B, [c.m_mu, c.m_pi0], weight_function=f.hnl_production_BRs.get_differential_weight_function("Bmeson", ["Mu","Pi0"]))
                    if b_pi0Mu_decay.success: 
                        successfull_decays.append(b_pi0Mu_decay)
                        curr_weight += self.multiplicities["Bmeson"]*f.get_branching_ratios("hnl", m_exo, "Bmeson", ["Mu","Pi0"], 1)
                    del b_pi0Mu_decay
                if m_exo+c.m_mu+c.m_pi<c.m_B0:
                    b0_piMu_decay = decay_to_n_body(self.b0_boosts, m_exo, c.m_B0, [c.m_mu, c.m_pi], weight_function=f.hnl_production_BRs.get_differential_weight_function("B0meson", ["Mu","Pi"]))
                    if b0_piMu_decay.success:
                        successfull_decays.append(b0_piMu_decay)
                        curr_weight += self.multiplicities["B0meson"]*f.get_branching_ratios("hnl", m_exo, "B0meson", ["Mu","Pi"], 1)
                    del b0_piMu_decay
                if successfull_decays:
                    b_piMu_decay = np.sum(successfull_decays)
                    kdes.append(b_piMu_decay.kde())
                    weights.append(curr_weight)
                    del b_piMu_decay

                curr_weight = 0.
                successfull_decays = []
                if m_exo+c.m_mu+c.m_D0<c.m_B:
                    b_D0Mu_decay = decay_to_n_body(self.b_boosts, m_exo, c.m_B, [c.m_mu, c.m_D0], weight_function=f.hnl_production_BRs.get_differential_weight_function("Bmeson", ["Mu","D0"]))
                    if b_D0Mu_decay.success: 
                        successfull_decays.append(b_D0Mu_decay)
                        curr_weight += self.multiplicities["Bmeson"]*f.get_branching_ratios("hnl", m_exo, "Bmeson", ["Mu","D0"], 1)
                    del b_D0Mu_decay
                if m_exo+c.m_mu+c.m_D<c.m_B0:
                    b0_DMu_decay = decay_to_n_body(self.b0_boosts, m_exo, c.m_B0, [c.m_mu, c.m_D], weight_function=f.hnl_production_BRs.get_differential_weight_function("B0meson", ["Mu","D"]))
                    if b0_DMu_decay.success: 
                        successfull_decays.append(b0_DMu_decay)
                        curr_weight += self.multiplicities["B0meson"]*f.get_branching_ratios("hnl", m_exo, "B0meson", ["Mu","D"], 1)
                    del b0_DMu_decay
                if successfull_decays:
                    b_DMu_decay = np.sum(successfull_decays)
                    kdes.append(b_DMu_decay.kde())
                    weights.append(curr_weight)
                    del b_DMu_decay

                curr_weight = 0.
                successfull_decays = []
                if m_exo+c.m_mu+c.m_D0star<c.m_B:
                    b_D0starMu_decay = decay_to_n_body(self.b_boosts, m_exo, c.m_B, [c.m_mu, c.m_D0star], weight_function=f.hnl_production_BRs.get_differential_weight_function("Bmeson", ["Mu","D0star"]))
                    if b_D0starMu_decay.success: 
                        successfull_decays.append(b_D0starMu_decay)
                        curr_weight += self.multiplicities["Bmeson"]*f.get_branching_ratios("hnl", m_exo, "Bmeson", ["Mu","D0star"], 1)
                    del b_D0starMu_decay
                if m_exo+c.m_mu+c.m_Dstar<c.m_B0:
                    b0_DstarMu_decay = decay_to_n_body(self.b0_boosts, m_exo, c.m_B0,[c.m_mu, c.m_Dstar], weight_function=f.hnl_production_BRs.get_differential_weight_function("B0meson", ["Mu","Dstar"]))
                    if b0_DstarMu_decay.success: 
                        successfull_decays.append(b0_DstarMu_decay)
                        curr_weight += self.multiplicities["B0meson"]*f.get_branching_ratios("hnl", m_exo, "B0meson", ["Mu","Dstar"], 1)
                    del b0_DstarMu_decay
                if successfull_decays:
                    b_DstarMu_decay = np.sum(successfull_decays)
                    kdes.append(b_DstarMu_decay.kde())
                    weights.append(curr_weight)
                    del b_DstarMu_decay

                if m_exo+c.m_mu+c.m_Ds<c.m_Bs:
                    bs_DsMu_decay = decay_to_n_body(self.bs_boosts, m_exo, c.m_Bs, [c.m_mu, c.m_Ds], weight_function=f.hnl_production_BRs.get_differential_weight_function("Bsmeson", ["Mu","Ds"]))
                    if bs_DsMu_decay.success: 
                        kdes.append(bs_DsMu_decay.kde())
                        weights.append((self.multiplicities["Bsmeson"])*f.get_branching_ratios("hnl", m_exo, "Bsmeson", ["Mu","Ds"], 1))
                    del bs_DsMu_decay
                if m_exo+c.m_mu+c.m_Dsstar<c.m_Bs:
                    bs_DsstarMu_decay = decay_to_n_body(self.bs_boosts, m_exo, c.m_Bs, [c.m_mu, c.m_Dsstar], weight_function=f.hnl_production_BRs.get_differential_weight_function("Bsmeson", ["Mu","Dsstar"]))
                    if bs_DsstarMu_decay.success: 
                        kdes.append(bs_DsstarMu_decay.kde())
                        weights.append((self.multiplicities["Bsmeson"])*f.get_branching_ratios("hnl", m_exo, "Bsmeson", ["Mu","Dsstar"], 1))
                    del bs_DsstarMu_decay
                    
            elif self.active_coupling =="Tau" and (m_exo+c.m_tau<c.m_B):
                b_decay = decay_to_2_body(self.b_boosts, m_exo, c.m_B, c.m_tau)
                if b_decay.success: 
                    kdes.append(b_decay.kde())
                    weights.append(self.multiplicities["Bmeson"]*f.get_branching_ratios("hnl", m_exo, "Bmeson", "Tau", 1))
                del b_decay

                curr_weight = 0.
                successfull_decays = []
                if m_exo+c.m_tau+c.m_pi<c.m_B:
                    b_pi0Tau_decay = decay_to_n_body(self.b_boosts, m_exo, c.m_B, [c.m_tau, c.m_pi0], weight_function=f.hnl_production_BRs.get_differential_weight_function("Bmeson", ["Tau","Pi0"]))
                    if b_pi0Tau_decay.success: 
                        successfull_decays.append(b_pi0Tau_decay)
                        curr_weight += self.multiplicities["Bmeson"]*f.get_branching_ratios("hnl", m_exo, "Bmeson", ["Tau","Pi0"], 1)
                    del b_pi0Tau_decay
                if m_exo+c.m_tau+c.m_pi<c.m_B0:
                    b0_piTau_decay = decay_to_n_body(self.b0_boosts, m_exo, c.m_B0, [c.m_tau, c.m_pi], weight_function=f.hnl_production_BRs.get_differential_weight_function("B0meson", ["Tau","Pi"]))
                    if b0_piTau_decay.success:
                        successfull_decays.append(b0_piTau_decay)
                        curr_weight += self.multiplicities["B0meson"]*f.get_branching_ratios("hnl", m_exo, "B0meson", ["Tau","Pi"], 1)
                    del b0_piTau_decay
                if successfull_decays:
                    b_piTau_decay = np.sum(successfull_decays)
                    kdes.append(b_piTau_decay.kde())
                    weights.append(curr_weight)
                    del b_piTau_decay

        data_list_bk_exo = []
        for en_exo in self.en_list:
            data_sublist_bk = []
            if en_exo > m_exo:
                for th_exo in self.th_list:
                    norm_yield = 0.
                    for weight, kde in zip(weights, kdes):
                        if weight==0.: continue
                        norm_yield += f.f_kde(kde, self.norm_xsec*weight, en_exo, th_exo)[0]
                    data_sublist_bk.append([th_exo, en_exo, m_exo, norm_yield])
            else:
                for th_exo in self.th_list:
                    data_sublist_bk.append([th_exo, en_exo, m_exo, 0.])
            data_list_bk_exo.append([data_sublist_bk])  
        return data_list_bk_exo

    def process_pool_2s(self,nthreads,single_mass_point=0):
        """Wrapper to run production in parallel threads.
        Args:
            nthreads (int): number of parallel threads to run
            single_mass_point (float, optional): Single mass point for which to run. Defaults to 0.
        """

        export_units_header = "Theta[rad], E_x [GeV], m_x [GeV], dY per [rad GeV LambdaS]"
        if single_mass_point>0:
            list_bk = self.b_meson_decay_process_2s(single_mass_point)
            f.export(setup.experiments[self.exp],self.daughter + "_Bmeson2S_beam" + str(setup.p_beam[self.exp]) + "GeV_" + str(int(single_mass_point*1e6)) + "keV_" + setup.experiments[self.exp] + ".dat",np.reshape(list_bk,(len(self.th_list)*len(self.en_list),4)),header = self.export_header+'\n'+export_units_header)
            return
        for iFile in range(len(setup.mass_bins)):
            iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," with",nthreads,"threads")

            with Pool(processes=nthreads) as pool:
                list_bk = list(tqdm(pool.imap(self.b_meson_decay_process_2s,self.m_lists[iFile]),total=len(self.m_lists[iFile])))

            pool.close()
            pool.join()
            
            # export
            f.export(setup.experiments[self.exp],self.daughter_exo + "_Bmeson2S_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_bk,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)),header = self.export_header+'\n'+export_units_header)

        return
  
    def b_meson_decay_process_2s(self,m_exo):
        """evaluates open beauty meson decays into two dark scalar particles with mass m_exo emitted at angle theta to beam axis with energy e
        Args:
            m_exo (float): exotic mss

        Returns:
            array(float): [theta, e, m_exo, expected yield]
        """

        kdes = []
        weights = []

        for fin_state in setup.kaons_list:
            #charged
            m_K = c.m_K[fin_state]
            if(2*m_exo+m_K<c.m_B):
                diff_class = f.diff_BR_B_K_2S("Bmeson",fin_state)
                bk2s_decay = decay_to_n_body(self.b_boosts, m_exo, c.m_B, [m_K], n_exos=2, weight_function = diff_class.d_BR_B_K_2S)
                if bk2s_decay.success: 
                    kdes.append(bk2s_decay.kde())
                    weights.append(self.multiplicities["Bmeson"]*f.get_branching_ratios("2ds", m_exo, "Bmeson",fin_state,1)*2)
                del bk2s_decay
            
            #neutral
            if fin_state == "K" or fin_state == "Kstar_892" or fin_state == "K2star_1430":
                m_K = c.m_K[fin_state+"_0"]
            else: m_K = c.m_K[fin_state]
            if(2*m_exo+m_K<c.m_B0):
                diff_class = f.diff_BR_B_K_2S("B0meson",fin_state)
                b0k02s_decay = decay_to_n_body(self.b0_boosts, m_exo, c.m_B0, [m_K], n_exos=2, weight_function = diff_class.d_BR_B_K_2S)
                if b0k02s_decay.success: 
                    kdes.append(b0k02s_decay.kde())
                    weights.append(self.multiplicities["B0meson"]*f.get_branching_ratios("2ds", m_exo, "B0meson",fin_state,1)*2)
                del b0k02s_decay
        if(2*m_exo<c.m_B0): #B0->2S
            b0_decay = decay_to_2_body(self.b0_boosts, m_exo, c.m_B0, m_exo)
            if b0_decay.success: 
                kdes.append(b0_decay.kde())
                weights.append(self.multiplicities["B0meson"]*f.get_branching_ratios("2ds", m_exo, "B0meson", "2S", 1)*2)
            del b0_decay
        if(2*m_exo<c.m_Bs): #Bs->2S
            bs_decay = decay_to_2_body(self.b0_boosts, m_exo, c.m_Bs, m_exo)
            if bs_decay.success: 
                kdes.append(bs_decay.kde())
                weights.append(self.multiplicities["Bsmeson"]*f.get_branching_ratios("2ds", m_exo, "Bsmeson", "2S", 1)*2)
            del bs_decay

        data_list_bk_exo = []
        for en_exo in self.en_list:
            data_sublist_bk = []
            if en_exo > m_exo:
                for th_exo in self.th_list:
                    norm_yield = 0.
                    for weight, kde in zip(weights, kdes):
                        if weight==0.: continue
                        norm_yield += f.f_kde(kde, self.norm_xsec*weight, en_exo, th_exo)[0]
                    data_sublist_bk.append([th_exo, en_exo, m_exo, norm_yield])
            else:
                for th_exo in self.th_list:
                    data_sublist_bk.append([th_exo, en_exo, m_exo, 0.])
            data_list_bk_exo.append([data_sublist_bk])
        return data_list_bk_exo
    
class dmeson_decay_production:
    def __init__(self, experiment, n_production, n_decay, daugther_name, use_external=False, active_coupling="", use_3part_prod = False, use_massless = False, use_reweight="", use_empirical = "", delete_intermed_lists=True, target_code = "2212", beam_particle_code = "2212", pythia_Ppdf_code=20, pythia_extra_args=[]):
        self.exp = experiment
        self.nprod = n_production
        self.ndecays = n_decay
        self.daughter_exo = daugther_name 
        self.active_coupling = active_coupling

        p_beam = setup.p_beam.get(self.exp,400)
        self.th_list = np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins).tolist()
        self.en_list = np.linspace(setup.energy_min.get(p_beam), setup.energy_max.get(p_beam), num=setup.energy_bins).tolist()
        self.norm_xsec = c.sigma_cc[p_beam] / c.sigma_pp[p_beam] * c.A_target[experiment]**(1./3)

        self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced in Dmeson decays | generated with statistics of {n_production*n_decay} events using "
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())

        self.d_meson_list  = np.array([])
        self.d0_meson_list = np.array([])
        self.ds_meson_list = np.array([])
        self.dbar_meson_list  = np.array([])
        self.d0bar_meson_list = np.array([])
        self.dsbar_meson_list = np.array([])
        self.empirical = ""
        self.external_mod = ""
        if use_external and not use_empirical:
            print("[Info:] \t Using external D meson source")
            filename = os.path.dirname(os.path.realpath(__file__))+"/../tab_mesons/charm/"
            filename += "hsccbar_"+("massless_" if use_massless else "")+("statcomb23partprod_" if use_3part_prod else "")+str(setup.p_beam[self.exp])+'GeV_'+f.number_to_3sigfigs_str(n_production)+".dat"
            if use_3part_prod: self.external_mod = "23parton"
            if target_code != "2212": 
                PID_to_extension =  {"2112":"pn", "100822080":"pPb"}
                assert target_code in PID_to_extension.keys(), "[Error:] \t target code " + target_code + " not recognised."
                filename = filename.replace(".dat", "_" + PID_to_extension[target_code ] + ".dat")
                self.external_mod += PID_to_extension[target_code ]
            if os.path.exists(filename):
                external_d_list = np.loadtxt(filename,skiprows=1)
                with open(filename) as meson_source_file: self.export_header += meson_source_file.readline()[1:-1]
                self.d_meson_list =     external_d_list[np.where(external_d_list[:,0]== 411)]
                self.dbar_meson_list =  external_d_list[np.where(external_d_list[:,0]==-411)]
                self.d0_meson_list =    external_d_list[np.where(external_d_list[:,0]== 421)]
                self.d0bar_meson_list = external_d_list[np.where(external_d_list[:,0]==-421)]
                self.ds_meson_list =    external_d_list[np.where(external_d_list[:,0]== 431)]
                self.dsbar_meson_list = external_d_list[np.where(external_d_list[:,0]==-431)]
            else:
                print("[Error:] \t Path:",filename,"does not exist. Exiting")
                sys.exit(1)
        elif use_empirical:
            self.empirical = use_empirical
            self._init_empirical_(use_empirical)
        else:
            print("[Info:] \t Starting Pythia")
            self._init_pythia_(beam_particle_code, target_code, pythia_Ppdf_code, pythia_extra_args)
            print("[Info:] \t Pythia production finished")

        # d_meson_mult = (len(self.d_meson_list)+len(self.dbar_meson_list))/self.nprod
        # print("[Info:] \t D meson multiplicity",d_meson_mult)
        # d0_meson_mult = (len(self.d0_meson_list)+len(self.d0bar_meson_list))/self.nprod
        # print("[Info:] \t D0 meson multiplicity",d0_meson_mult)
        # ds_meson_mult = (len(self.ds_meson_list)+len(self.dsbar_meson_list))/self.nprod
        # print("[Info:] \t Ds meson multiplicity",ds_meson_mult)
        print("[Info:] \t D+ meson multiplicity",len(self.d_meson_list)/self.nprod)
        print("[Info:] \t D- meson multiplicity",len(self.dbar_meson_list)/self.nprod)
        print("[Info:] \t D0 meson multiplicity",len(self.d0_meson_list)/self.nprod)
        print("[Info:] \t D0bar meson multiplicity",len(self.d0bar_meson_list)/self.nprod)
        print("[Info:] \t Ds meson multiplicity",len(self.ds_meson_list)/self.nprod)
        print("[Info:] \t Dsbar meson multiplicity",len(self.dsbar_meson_list)/self.nprod)
            
        d_meson_mult = (len(self.d_meson_list)+len(self.dbar_meson_list))/self.nprod
        # print("[Info:] \t D meson multiplicity",d_meson_mult)
        d0_meson_mult = (len(self.d0_meson_list)+len(self.d0bar_meson_list))/self.nprod
        # print("[Info:] \t D0 meson multiplicity",d0_meson_mult)
        ds_meson_mult = (len(self.ds_meson_list)+len(self.dsbar_meson_list))/self.nprod
        # print("[Info:] \t Ds meson multiplicity",ds_meson_mult)
        
        self.multiplicities = {
            "Dmeson"  : d_meson_mult,
            "D0meson" : d0_meson_mult,
            "Dsmeson" : ds_meson_mult
        }

        # readjusting input to predefined xF distributions 
        if use_reweight:
            if use_reweight not in f.Input_reweight.literature.keys(): 
                print("[Info:] \t Rescaling reference",use_reweight,"not recognised. Exiting")
                sys.exit(1)
            print("[Info:] \t Sampling meson distributions to match xF differential cross sections presented in",f.Input_reweight.literature[use_reweight])
            self.d0_meson_list = f.Input_reweight.reshape_input(self.d0_meson_list, self.ndecays, normalize_to=use_reweight, meson = "D0", p_beam = setup.p_beam[self.exp]) 
            self.d_meson_list  = f.Input_reweight.reshape_input(self.d_meson_list,  self.ndecays, normalize_to=use_reweight, meson = "D",  p_beam = setup.p_beam[self.exp]) 
            self.ds_meson_list = f.Input_reweight.reshape_input(self.ds_meson_list, self.ndecays, normalize_to=use_reweight, meson = "Ds", p_beam = setup.p_beam[self.exp]) 
            self.d0bar_meson_list = f.Input_reweight.reshape_input(self.d0bar_meson_list, self.ndecays, normalize_to=use_reweight, meson = "D0bar", p_beam = setup.p_beam[self.exp]) 
            self.dbar_meson_list  = f.Input_reweight.reshape_input(self.dbar_meson_list,  self.ndecays, normalize_to=use_reweight, meson = "Dbar",  p_beam = setup.p_beam[self.exp]) 
            self.dsbar_meson_list = f.Input_reweight.reshape_input(self.dsbar_meson_list, self.ndecays, normalize_to=use_reweight, meson = "Dsbar", p_beam = setup.p_beam[self.exp]) 
        self.reweight_ref = use_reweight

        # generating lorentz boost lists
        ndecays = self.ndecays if not use_reweight else 1
        if self.daughter_exo == "alp":
            if delete_intermed_lists: del self.ds_meson_list
            self.d0_boosts = f.lorentz_boosts_vectorized(event_list=np.concatenate((self.d0_meson_list, self.d0bar_meson_list)), mass=c.m_D0, ndecays=ndecays)
            if delete_intermed_lists: del self.d0_meson_list, self.d0bar_meson_list
            self.d_boosts =  f.lorentz_boosts_vectorized(event_list=np.concatenate((self.d_meson_list, self.dbar_meson_list)),  mass=c.m_D,  ndecays=ndecays)
            if delete_intermed_lists: del self.d_meson_list, self.dbar_meson_list

        elif self.daughter_exo == "hnl":
            self.d0_boosts = f.lorentz_boosts_vectorized(event_list=np.concatenate((self.d0_meson_list,self.d0bar_meson_list)), mass=c.m_D0, ndecays=ndecays)
            if delete_intermed_lists: del self.d0_meson_list, self.d0bar_meson_list
            self.d_boosts =  f.lorentz_boosts_vectorized(event_list=np.concatenate((self.d_meson_list,self.dbar_meson_list)),  mass=c.m_D,  ndecays=ndecays)
            if delete_intermed_lists: del self.d_meson_list, self.dbar_meson_list
            self.ds_boosts = f.lorentz_boosts_vectorized(event_list=np.concatenate((self.ds_meson_list,self.dsbar_meson_list)), mass=c.m_Ds, ndecays=ndecays)
            if delete_intermed_lists: del self.ds_meson_list, self.dsbar_meson_list

            tau_lept_list = np.append(decay_to_2_body.secondary_generator(rng.choice(self.d_boosts, size = int(c.BR_D_NuTau/c.BR_Ds_NuTau * np.size(self.d_boosts)), replace = False), c.m_D), decay_to_2_body.secondary_generator(self.ds_boosts, c.m_Ds), axis = 0)
            self.tau_lept_mult = c.BR_Ds_NuTau * ds_meson_mult + c.BR_D_NuTau * d_meson_mult
            tau_MC_lept_mult = ds_meson_mult + c.BR_D_NuTau/c.BR_Ds_NuTau*d_meson_mult
            self.multiplicities["Tau"] = self.tau_lept_mult
            print("[Info:] \t Tau lepton effective multiplicity", self.tau_lept_mult, "and statistical multiplicity", tau_MC_lept_mult )
            self.tau_boosts = f.lorentz_boosts_vectorized(event_list=tau_lept_list,  mass=c.m_tau,  ndecays=1)

    def _init_pythia_(self, beam_particle_code, target_code, pythia_Ppdf_code, pythia_extra_args, save_pythia_file=False):
        #import pythia
        if not os.path.exists(setup.pythia_dir):
            print("[Error:] \t Set correct path to Pythia in general/exotic_production_setup.py")
            sys.exit(1)
        pythia_cfg_dir = setup.pythia_dir+"/examples/"
        if not os.path.exists(pythia_cfg_dir):
            pythia_cfg_dir = setup.pythia_dir+"/share/Pythia8/examples/"
        if not os.path.exists(pythia_cfg_dir):
            print("[Error:] \t Set correct path to Pythia")
            sys.exit(1)

        # cfg = open(setup.pythia_dir+"/examples/Makefile.inc")
        cfg = open(pythia_cfg_dir+"Makefile.inc")
        lib = setup.pythia_dir+"/lib"
        for line in cfg:
            if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
        sys.path.insert(0, lib)
        import pythia8
        pythia = pythia8.Pythia("-v", False)
        print("[Info:] \t Found Pythia version "+str(pythia.parm("Pythia:versionNumber")))
        pythia_modifiers = ["HardQCD:hardccbar = on",f"Beams:idA = {beam_particle_code}",f"Beams:idB = {target_code}","Beams:eA = " + str(setup.p_beam[self.exp]),"Beams:eB = 0","Beams:frameType = 2", f"PDF:pSet={pythia_Ppdf_code}"]#"PhaseSpace:mHatMin=0",
        if isinstance(pythia_extra_args,str): pythia_modifiers.append(pythia_extra_args)
        elif len(pythia_extra_args): pythia_modifiers += pythia_extra_args
        for modifier in pythia_modifiers+["Print:quiet=on","Print:errors=off","Next:numberShowEvent = 0"]: pythia.readString(modifier)
        self.export_header += f"Pythia{pythia.parm('Pythia:versionNumber')} with settings [" +" ".join(pythia_modifiers)+"]"
        pythia.init() # initialize
        evtID_to_list = {
            411: [], -411:[],
            421: [], -421:[],
            431: [], -431:[]
        }
        for _ in range(self.nprod):
            pythia.next()
            for isubEvent in range(pythia.event.size()):
                evt = pythia.event[isubEvent]
                if evt.id() not in evtID_to_list.keys(): continue
                evtID_to_list[evt.id()].append([evt.id(),evt.px(),evt.py(),evt.pz()])
        self.d_meson_list     = np.array(evtID_to_list[411])
        self.dbar_meson_list  = np.array(evtID_to_list[-411])
        self.d0_meson_list    = np.array(evtID_to_list[421])
        self.d0bar_meson_list = np.array(evtID_to_list[-421])
        self.ds_meson_list    = np.array(evtID_to_list[431])
        self.dsbar_meson_list = np.array(evtID_to_list[-431])
        if save_pythia_file:
            file_info  =  "Charmed meson momenta as generated for " + self.nprod +" impining protons, with PYTHIA (version "+str(pythia.parm("Pythia:versionNumber"))+") and modifiers " + ', '.join(pythia_modifiers)
            np.savetxt(f".ALP_production/tab_mesons/charm/hsccbar_{setup.p_beam[self.exp]}GeV_{f.number_to_3sigfigs_str(self.nprod)}.dat", np.concatenate(evtID_to_list.values()), header=file_info, fmt='%i %4f %4f %3f')
        return

    def _init_empirical_(self, reference = "LEBC"):
        if reference == "SMTs": sys.exit("[Error:] \t SHiP meson tables as generated according to SHiP-NOTE-2015-009 are not publically available ")
        print(f"[Info:] \t Generating meson arrays from {reference} distributions ("+f.Input_reweight.literature.get(reference,"no lit. entry found")+")" )
        # if reference == "SMTs": Ship Meson files not publically available 
        #     tablepath = './../../SHiPGeneratedFiles/'
        #     self.d_meson_list       = rng.choice(np.loadtxt(tablepath+ '411.txt'),size = int(self.nprod * 0.207))
        #     self.d0_meson_list      = rng.choice(np.loadtxt(tablepath+ '421.txt'),size = int(self.nprod * 0.632))
        #     self.ds_meson_list      = rng.choice(np.loadtxt(tablepath+ '431.txt'),size = int(self.nprod * 0.088))
        #     self.dbar_meson_list    = rng.choice(np.loadtxt(tablepath+'-411.txt'),size = int(self.nprod * 0.207))
        #     self.d0bar_meson_list   = rng.choice(np.loadtxt(tablepath+'-421.txt'),size = int(self.nprod * 0.632))
        #     self.dsbar_meson_list   = rng.choice(np.loadtxt(tablepath+'-431.txt'),size = int(self.nprod * 0.088))
        # else:
        self.d_meson_list       = f.Input_reweight.generate_momentum_array("D",     N_Interactions = self.nprod, random_generator = rng, p_beam = setup.p_beam[self.exp], reference=reference)
        self.d0_meson_list      = f.Input_reweight.generate_momentum_array("D0",    N_Interactions = self.nprod, random_generator = rng, p_beam = setup.p_beam[self.exp], reference=reference)
        self.ds_meson_list      = f.Input_reweight.generate_momentum_array("Ds",    N_Interactions = self.nprod, random_generator = rng, p_beam = setup.p_beam[self.exp], reference=reference)
        self.dbar_meson_list    = f.Input_reweight.generate_momentum_array("Dbar",  N_Interactions = self.nprod, random_generator = rng, p_beam = setup.p_beam[self.exp], reference=reference)
        self.d0bar_meson_list   = f.Input_reweight.generate_momentum_array("D0bar", N_Interactions = self.nprod, random_generator = rng, p_beam = setup.p_beam[self.exp], reference=reference)
        self.dsbar_meson_list   = f.Input_reweight.generate_momentum_array("Dsbar", N_Interactions = self.nprod, random_generator = rng, p_beam = setup.p_beam[self.exp], reference=reference)
        self.export_header += reference + f"({f.Input_reweight.literature.get(reference,'no lit. entry found')}) empirical meson distributions."
        return

    def process_pool(self, nthreads, active_couplings = [""], single_mass_point=0):
        """Wrapper to run production in parallel threads.
        Args:
            nthreads (int): number of parallel threads to run
            single_mass_point (float, optional): Single mass point for which to run. Defaults to 0.
        """

        for active_coupling in active_couplings:
            export_units_header = "Theta[rad] E_x[GeV], m_x [GeV], dY[per(rad GeV " + {"hnl":f"U2_{active_coupling}", "alp":"g_cu_eff^2/GeV^2", "ds":"Y^2" }[self.daughter_exo] + ')]'

            mixing_info  = "-" + active_coupling + "Mixing" if active_coupling else ""
            if active_coupling: print("[Info:] \t Starting "+self.daughter_exo+" production with active coupling U2_" + active_coupling)
            export_info_header = self.export_header if not active_coupling else self.export_header.replace(' yield', f'({active_coupling} mixing) yield')
            self.active_coupling = active_coupling

            if single_mass_point>0:
                list_d = self.d_meson_decay_process(single_mass_point)
                f.export(setup.experiments[self.exp],self.daughter + "_Dmeson"+mixing_info+"_beam" + str(setup.p_beam[self.exp]) + "GeV_" + str(int(single_mass_point*1e6)) + "keV_" + setup.experiments[self.exp] + ".dat",np.reshape(list_d,(len(self.th_list)*len(self.en_list),4)),header = self.export_header+'\n'+export_units_header)
            else:
                for iFile in range(len(setup.mass_bins)):
                    iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
                    if setup.mass_min[iFile]*1000 < 1:
                        iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
                    print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile]," GeV with",nthreads,"threads")

                    with Pool(processes=nthreads) as pool:
                        list_d = list(tqdm(pool.imap(self.d_meson_decay_process, self.m_lists[iFile]), total=len(self.m_lists[iFile])))

                    pool.close()
                    pool.join()

                    # export   
                    f.export(setup.experiments[self.exp],self.daughter_exo + "_Dmeson"+mixing_info+"_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_d,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)),header = export_info_header+'\n'+export_units_header)
        return

    def d_meson_decay_process(self,m_exo):
        """evaluates open charmed meson decays into feebly interacting particle with mass m_exo emitted at angle theta to beam axis with energy e
        Args:
            m_exo (float): exotic mss

        Returns:
            array(float): [theta, e, m_exo, expected yield]
        """
        kdes = []
        weights = []
        if self.daughter_exo == "alp" or self.daughter_exo == "ds":
            successfull_decays  = []
            curr_weight = 0.
            if(m_exo+c.m_pi <c.m_D): 
                d_decay = decay_to_2_body(self.d_boosts, m_exo, c.m_D, c.m_pi)
                if d_decay.success: 
                    successfull_decays.append(d_decay)
                    curr_weight += self.multiplicities["Dmeson"]*f.get_branching_ratios(self.daughter_exo, m_exo, "Dmeson", "Pi", 1)
                del d_decay
            if(m_exo+c.m_pi0<c.m_D0):
                d0_decay = decay_to_2_body(self.d_boosts, m_exo, c.m_D0, c.m_pi0)
                if d0_decay.success: 
                    successfull_decays.append(d0_decay)
                    curr_weight += self.multiplicities["D0meson"]*f.get_branching_ratios(self.daughter_exo, m_exo, "D0meson", "Pi0", 1)
                del d0_decay
            if successfull_decays:
                d_to_pi_decay = np.sum(successfull_decays)
                kdes.append(d_to_pi_decay.kde())
                weights.append(curr_weight)
                del  d_to_pi_decay

        elif self.daughter_exo == "hnl":
            if self.active_coupling == "El":
                if(m_exo+c.m_el< c.m_D): 
                    d_decay = decay_to_2_body(self.d_boosts,  m_exo, c.m_D, c.m_el)  
                    if d_decay.success: 
                        kdes.append(d_decay.kde())
                        weights.append(self.multiplicities["Dmeson"]*f.get_branching_ratios("hnl", m_exo, "Dmeson", "El", 1))
                    del d_decay
                if(m_exo+c.m_el< c.m_Ds): 
                    ds_decay = decay_to_2_body(self.ds_boosts, m_exo, c.m_Ds, c.m_el)  
                    if ds_decay.success:
                        kdes.append(ds_decay.kde())
                        weights.append(self.multiplicities["Dsmeson"]*f.get_branching_ratios("hnl", m_exo, "Dsmeson", "El", 1))
                    del ds_decay
                if(m_exo+c.m_el+c.m_K["K_0"]< c.m_D): 
                    d_decay = decay_to_n_body(np.concatenate((self.d_boosts,self.d0_boosts)),  m_exo, c.m_D, [c.m_el, c.m_K["K_0"]], weight_function=f.hnl_production_BRs.get_differential_weight_function("Dmeson", ["El","K_0"]))  
                    if d_decay.success: 
                        kdes.append(d_decay.kde())
                        weights.append(self.multiplicities["Dmeson"]*f.get_branching_ratios("hnl", m_exo, "Dmeson", ["El","K_0"], 1)
                                        +self.multiplicities["D0meson"]*f.get_branching_ratios("hnl", m_exo, "D0meson", ["El","K"], 1))
                    del d_decay
                if(m_exo+c.m_el+c.m_K["Kstar_892_0"]< c.m_D):
                    d_elkstar_decay = decay_to_n_body(np.concatenate((self.d_boosts,self.d0_boosts)),  m_exo, c.m_D, [c.m_el, c.m_K["Kstar_892_0"]], weight_function=f.hnl_production_BRs.get_differential_weight_function("Dmeson", ["El","Kstar_892_0"]))  
                    if d_elkstar_decay.success: 
                        kdes.append(d_elkstar_decay.kde())
                        weights.append(self.multiplicities["Dmeson"]*f.get_branching_ratios("hnl", m_exo, "Dmeson", ["El","Kstar_892_0"], 1)
                                        +self.multiplicities["D0meson"]*f.get_branching_ratios("hnl", m_exo, "D0meson", ["El","Kstar_892"], 1))
                    del d_elkstar_decay
                if(m_exo+c.m_el+c.m_pi0<c.m_D): 
                    d_elpi_decay = decay_to_n_body(np.concatenate((self.d_boosts,self.d0_boosts)),  m_exo, c.m_D, [c.m_el, c.m_pi0], weight_function=f.hnl_production_BRs.get_differential_weight_function("Dmeson", ["El","Pi0"]))
                    if d_elpi_decay.success: 
                        kdes.append(d_elpi_decay.kde())
                        weights.append(self.multiplicities["Dmeson"]*f.get_branching_ratios("hnl", m_exo, "Dmeson", ["El","Pi0"], 1)+
                                        self.multiplicities["D0meson"]*f.get_branching_ratios("hnl", m_exo, "D0meson", ["El","Pi"], 1))
                    del d_elpi_decay
            
            elif self.active_coupling == "Mu":
                if(m_exo+c.m_mu< c.m_D): 
                    d_decay = decay_to_2_body(self.d_boosts, m_exo, c.m_D, c.m_mu)  
                    if d_decay.success: 
                        kdes.append(d_decay.kde())
                        weights.append(self.multiplicities["Dmeson"]*f.get_branching_ratios("hnl", m_exo, "Dmeson", "Mu", 1))
                    del d_decay
                if(m_exo+c.m_mu< c.m_Ds): 
                    ds_decay = decay_to_2_body(self.ds_boosts, m_exo, c.m_Ds, c.m_mu)
                    if ds_decay.success: 
                        kdes.append(ds_decay.kde())
                        weights.append(self.multiplicities["Dsmeson"]*f.get_branching_ratios("hnl", m_exo, "Dsmeson", "Mu", 1))
                    del ds_decay
                if(m_exo+c.m_mu+c.m_K["K_0"]< c.m_D): 
                    d_decay = decay_to_n_body(np.concatenate((self.d_boosts,self.d0_boosts)),  m_exo, c.m_D, [c.m_mu, c.m_K["K_0"]], weight_function=f.hnl_production_BRs.get_differential_weight_function("Dmeson", ["Mu","K_0"]))
                    if d_decay.success: 
                        kdes.append(d_decay.kde())
                        weights.append(self.multiplicities["Dmeson"]*f.get_branching_ratios("hnl", m_exo, "Dmeson", ["Mu","K_0"], 1)
                                        +self.multiplicities["D0meson"]*f.get_branching_ratios("hnl", m_exo, "D0meson", ["Mu","K"], 1))
                    del d_decay
                if(m_exo+c.m_mu+c.m_K["Kstar_892_0"]< c.m_D): 
                    d_mukstar_decay = decay_to_n_body(np.concatenate((self.d_boosts,self.d0_boosts)),  m_exo, c.m_D, [c.m_mu, c.m_K["Kstar_892_0"]], weight_function=f.hnl_production_BRs.get_differential_weight_function("Dmeson", ["Mu","Kstar_892_0"]))
                    if d_mukstar_decay.success: 
                        kdes.append(d_mukstar_decay.kde())
                        weights.append(self.multiplicities["Dmeson"]*f.get_branching_ratios("hnl", m_exo, "Dmeson", ["Mu","Kstar_892_0"], 1)
                                        +self.multiplicities["D0meson"]*f.get_branching_ratios("hnl", m_exo, "D0meson", ["Mu","Kstar_892"], 1))
                    del d_mukstar_decay
                if(m_exo+c.m_mu+c.m_pi0<c.m_D): 
                    d_mupi_decay = decay_to_n_body(np.concatenate((self.d_boosts,self.d0_boosts)),  m_exo, c.m_D, [c.m_mu, c.m_pi0], weight_function=f.hnl_production_BRs.get_differential_weight_function("Dmeson", ["Mu","Pi0"]))
                    if d_mupi_decay.success: 
                        kdes.append(d_mupi_decay.kde())
                        weights.append(self.multiplicities["Dmeson"]*f.get_branching_ratios("hnl", m_exo, "Dmeson", ["Mu","Pi0"], 1)+
                                       self.multiplicities["D0meson"]*f.get_branching_ratios("hnl", m_exo, "D0meson", ["Mu","Pi"], 1))
                    del d_mupi_decay

            elif self.active_coupling == "Tau":
                if(m_exo+c.m_tau< c.m_Ds): 
                    ds_decay = decay_to_2_body(self.ds_boosts, m_exo, c.m_Ds, c.m_tau)  
                    if ds_decay.success: 
                        kdes.append(ds_decay.kde())
                        weights.append(self.multiplicities["Dsmeson"]*f.get_branching_ratios("hnl", m_exo, "Dsmeson", "Tau", 1))
                    del ds_decay
                if(m_exo+c.m_el  < c.m_tau): 
                    tau_elnu_decay = decay_to_n_body(self.tau_boosts, m_exo, c.m_tau, [c.m_el, 0.], weight_function=f.hnl_production_BRs.get_differential_weight_function("Tau", ["El","Nu"]))
                    if tau_elnu_decay.success:
                        kdes.append(tau_elnu_decay.kde())
                        weights.append(self.multiplicities["Tau"]*f.get_branching_ratios("hnl", m_exo, "Tau",["El","Nu"], 1))
                    del tau_elnu_decay
                if(m_exo+c.m_mu  < c.m_tau): 
                    tau_munu_decay = decay_to_n_body(self.tau_boosts, m_exo, c.m_tau, [c.m_mu, 0.], weight_function=f.hnl_production_BRs.get_differential_weight_function("Tau", ["Mu","Nu"]))
                    if tau_munu_decay.success:
                        kdes.append(tau_munu_decay.kde())
                        weights.append(self.multiplicities["Tau"]*f.get_branching_ratios("hnl", m_exo, "Tau",["Mu","Nu"], 1))
                    del tau_munu_decay
                if(m_exo+c.m_pi  < c.m_tau): 
                    tau_pi_decay = decay_to_2_body(self.tau_boosts, m_exo, c.m_tau, c.m_pi)
                    if tau_pi_decay.success:
                        kdes.append(tau_pi_decay.kde())
                        weights.append(self.multiplicities["Tau"]*f.get_branching_ratios("hnl", m_exo, "Tau", "Pi", 1))
                    del tau_pi_decay
                # if(m_exo+c.m_rho < c.m_tau): 
                    # tau_rho_decay = decay_to_2_body(self.tau_boosts, m_exo, c.m_tau, c.m_rho)
                    # if tau_rho_decay.success: 
                    #     kdes.append(tau_rho_decay.kde())
                    #     weights.append(self.multiplicities["Tau"]*f.get_branching_ratios("hnl", m_exo, "Tau", "Rho", 1))
                    # del tau_rho_decay
                if(m_exo+c.m_pi+c.m_pi0<c.m_tau): #alternative to Nrho decay using full 3 body decay for future updates
                    tau_pi0pi_decay = decay_to_n_body(self.tau_boosts, m_exo, c.m_tau, [c.m_pi0, c.m_pi], weight_function=f.hnl_production_BRs.get_differential_weight_function("Tau", ["Pi0","Pi"]))
                    if tau_pi0pi_decay.success: 
                        kdes.append(tau_pi0pi_decay.kde())
                        weights.append(self.multiplicities["Tau"]*f.get_branching_ratios("hnl", m_exo, "Tau",["Pi0","Pi"], 1))
                    del tau_pi0pi_decay
        data_list_dpi_exo = []
        for en_exo in self.en_list:
            data_sublist_dpi = []
            if en_exo > m_exo:
                for th_exo in self.th_list:
                    norm_yield = 0.
                    for weight, kde in zip(weights, kdes):
                        if weight == 0.: continue
                        norm_yield += f.f_kde(kde, self.norm_xsec*weight, en_exo,th_exo)[0]
                    data_sublist_dpi.append([th_exo, en_exo, m_exo, norm_yield])
            else:
                for th_exo in self.th_list:
                    data_sublist_dpi.append([th_exo, en_exo, m_exo, 0.])

            data_list_dpi_exo.append(data_sublist_dpi)

        return data_list_dpi_exo

class meson_to_dp_decay_production:
    def __init__(self, experiment, n_production, n_decay, daugther_name, use_external=False):
        self.exp = experiment
        self.nprod = n_production
        self.ndecays = n_decay
        self.daughter_exo = daugther_name
        self.export_header = f"Differential {daugther_name} yield at {experiment} experiment produced in light meson decays | generated with statistics of {n_production*n_decay} events using "

        self.th_list = np.linspace(setup.theta_min.get(self.exp,0.00018), setup.theta_max.get(self.exp,0.01098), num=setup.theta_bins).tolist()
        self.en_list = np.linspace(setup.energy_min.get(setup.p_beam.get(self.exp,400)), setup.energy_max.get(setup.p_beam.get(self.exp,400)), num=setup.energy_bins).tolist()
        self.m_lists = []
        for iFile in range(len(setup.mass_bins)):
            self.m_lists.append(np.logspace(np.log10(setup.mass_min[iFile]), np.log10(setup.mass_max[iFile]), num=setup.mass_bins[iFile]).tolist())
            # self.m_lists.append(np.linspace(setup.mass_min[iFile], setup.mass_max[iFile], num=setup.mass_bins[iFile]).tolist())

        self.pi0_list  = np.empty((0, 4), np.float64)
        self.eta_list = np.empty((0, 4), np.float64)
        self.etap_list = np.empty((0, 4), np.float64)
        self.rho0_list  = np.empty((0, 4), np.float64)
        self.rho_list  = np.empty((0, 4), np.float64)
        self.omega_list = np.empty((0, 4), np.float64)
        self.phi_list = np.empty((0, 4), np.float64)

        if use_external:
            print("[Info:] \t Using external meson source")
            filename = os.path.dirname(os.path.realpath(__file__))+"/../tab_mesons/softQCD/softQCD_"+f.number_to_3sigfigs_str(n_production)+"evts_pp_"+str(setup.p_beam[self.exp])+'GeV_'+("8.3" if setup.p_beam[self.exp] in [30,800] else "8.2" )+"_pSet2_pom_ok.txt"

            if os.path.exists(filename):
                external_list = np.loadtxt(filename,usecols=(2,3,4,5))
                with open(filename) as meson_source_file: self.export_header += meson_source_file.readline()[1:-1]
                external_list[:,0] = np.abs(external_list[:,0])
                self.pi0_list = external_list[np.where(external_list[:,0]==111)]
                self.eta_list = external_list[np.where(external_list[:,0]==221)]
                self.etap_list = external_list[np.where(external_list[:,0]==331)]
                self.rho0_list = external_list[np.where(external_list[:,0]==113)]
                self.rho_list = external_list[np.where(external_list[:,0]==213)]
                self.omega_list = external_list[np.where(external_list[:,0]==223)]
                self.phi_list = external_list[np.where(external_list[:,0]==333)]

            else:
                print("[Error:] \t Path:",filename,"does not exist. Exiting")
                sys.exit(1)

        else:
            print("[Info:] \t Starting Pythia")
            self._init_pythia_()

            print("[Info:] \t Pythia production finished")

        pi0_mult = len(self.pi0_list)/self.nprod
        print("[Info:] \t pi0 multiplicity",pi0_mult)
        eta_mult = len(self.eta_list)/self.nprod
        print("[Info:] \t eta multiplicity",eta_mult)
        etap_mult = len(self.etap_list)/self.nprod
        print("[Info:] \t eta' multiplicity",etap_mult)
        rho0_mult = len(self.rho0_list)/self.nprod
        print("[Info:] \t rho0 multiplicity",rho0_mult)
        rho_mult = len(self.rho_list)/self.nprod
        print("[Info:] \t rho multiplicity",rho_mult)
        omega_mult = len(self.omega_list)/self.nprod
        print("[Info:] \t omega multiplicity",omega_mult)
        phi_mult = len(self.phi_list)/self.nprod
        print("[Info:] \t phi multiplicity",phi_mult)
        
        self.multiplicities = {
            "Pi0"  : pi0_mult,
            "Eta" : eta_mult,
            "EtaP" : etap_mult,
            "Rho0"  : rho0_mult,
            "Rho"  : rho_mult,
            "Omega" : omega_mult,
            "Phi" : phi_mult
        }

        # generating lorentz boost lists
        self.pi0_boosts = f.lorentz_boosts_vectorized(event_list=self.pi0_list, mass=c.m_pi0, ndecays=self.ndecays)
        del self.pi0_list
        self.eta_boosts =  f.lorentz_boosts_vectorized(event_list=self.eta_list,  mass=c.m_eta,  ndecays=self.ndecays)
        del self.eta_list
        self.etap_boosts = f.lorentz_boosts_vectorized(event_list=self.etap_list, mass=c.m_etap, ndecays=self.ndecays)
        del self.etap_list
        self.rho0_boosts = f.lorentz_boosts_vectorized(event_list=self.rho0_list, mass=c.m_rho, ndecays=self.ndecays)
        del self.rho0_list
        self.rho_boosts =  f.lorentz_boosts_vectorized(event_list=self.rho_list,  mass=c.m_rho,  ndecays=self.ndecays)
        del self.rho_list
        self.omega_boosts = f.lorentz_boosts_vectorized(event_list=self.omega_list, mass=c.m_omega, ndecays=self.ndecays)
        del self.omega_list
        self.phi_boosts = f.lorentz_boosts_vectorized(event_list=self.phi_list, mass=c.m_phi, ndecays=self.ndecays)
        del self.phi_list

    def _init_pythia_(self):
        #import pythia
        if not os.path.exists(setup.pythia_dir):
            print("[Error:] \t Set correct path to Pythia in general/exotic_production_setup.py")
            sys.exit(1)
        pythia_cfg_dir = setup.pythia_dir+"/examples/"
        if not os.path.exists(pythia_cfg_dir):
            pythia_cfg_dir = setup.pythia_dir+"/share/Pythia8/examples/"
        if not os.path.exists(pythia_cfg_dir):
            print("[Error:] \t Set correct path to Pythia")
            sys.exit(1)

        # cfg = open(setup.pythia_dir+"/examples/Makefile.inc")
        cfg = open(pythia_cfg_dir+"Makefile.inc")
        lib = setup.pythia_dir+"/lib"
        for line in cfg:
            if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
        sys.path.insert(0, lib)
        import pythia8
        pythia = pythia8.Pythia("", False)
        print("[Info:] \t Found Pythia version "+str(pythia.parm("Pythia:versionNumber")))
        pythia_modifiers = ["SoftQCD:all = on", "pdf:pHardSet = 2","Beams:idA = 2212","Beams:idB = 2212","Beams:eA = " + str(np.sqrt(setup.p_beam[self.exp]**2+0.93827**2)), "Beams:eB = "+str(0.93827), "Beams:frameType = 2"]
        for modifier in pythia_modifiers+["Print:quiet=on","Print:errors=off","Next:numberShowEvent = 0"]: pythia.readString(modifier)
        self.export_header += f"Pythia{pythia.parm('Pythia:versionNumber')} with settings [" +" ".join(pythia_modifiers)+"]"
        pythia.init() # initialize
        evtID_to_list = {}
        for id_ in [111,221,331,113,213,223,333]: evtID_to_list[id_] = []
        for _ in range(self.nprod):
            pythia.next()
            for isubEvent in range(pythia.event.size()):
                evtID = pythia.event[isubEvent].idAbs()
                if evtID not in [111,221,331,113,213,223,333]: continue
                evtID_to_list[evtID].append([evtID,pythia.event[isubEvent].px(),pythia.event[isubEvent].py(),pythia.event[isubEvent].pz()])
        self.pi0_list  = np.array(evtID_to_list[111])
        self.rho_list = np.array(evtID_to_list[213])
        self.rho0_list = np.array(evtID_to_list[113])
        self.eta_list  = np.array(evtID_to_list[221])
        self.omega_list= np.array(evtID_to_list[223])
        self.etap_list = np.array(evtID_to_list[331])
        self.phi_list  = np.array(evtID_to_list[333])
        return
    
    def process_pool(self,nthreads,single_mass_point=0):
        """Wrapper to run production in parallel threads.
        Args:
            nthreads (int): number of parallel threads to run
            single_mass_point (float, optional): Single mass point for which to run. Defaults to 0.
        """
        export_units_header = "Theta[rad] E_x [GeV] dY[per(rad GeV N_pN eps^2)]"
        if single_mass_point>0:
            list_out = self.meson_decay_process(single_mass_point)
            f.export(setup.experiments[self.exp],self.daughter + "_MesonDecay_beam" + str(setup.p_beam[self.exp]) + "GeV_" + str(int(single_mass_point*1e6)) + "keV_" + setup.experiments[self.exp] + ".dat",np.reshape(list_out,(len(self.th_list)*len(self.en_list),4)),header = self.export_header+'\n'+export_units_header)
            return
        for iFile in range(len(setup.mass_bins)):
            iFileName = str(int(setup.mass_min[iFile]*1000)) + "to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            if setup.mass_min[iFile]*1000 < 1:
                iFileName = "01to" + str(int(setup.mass_max[iFile]*1000)) + "MeV"
            print("[Info:] \t Starting exotic production for mass range",setup.mass_min[iFile],"to",setup.mass_max[iFile],"GeV with",nthreads,"threads")

            with Pool(processes=nthreads) as pool:
                list_out = list(tqdm(pool.imap(self.meson_decay_process, self.m_lists[iFile]), total=len(self.m_lists[iFile])))

            pool.close()
            pool.join()
            # export
            f.export(setup.experiments[self.exp],self.daughter_exo+"_MesonDecay_beam" + str(setup.p_beam[self.exp]) + "GeV_" + iFileName + "_" + setup.experiments[self.exp] + ".dat",np.reshape(list_out,(len(self.th_list)*len(self.en_list)*len(self.m_lists[iFile]),4)), header = self.export_header+'\n'+export_units_header)
        return

    def meson_decay_process(self,m_exo):
        """evaluates light meson decays into dark photon with mass m_dp emitted at angle theta to beam axis with energy e
        Args:
            m_dp (float): dark photon mass

        Returns:
            array(float): [theta, e, m_dp, expected yield]
        """
        kdes = []
        weights = []

        #pseudoscalar to gamma decays
        if(m_exo<c.m_pi0):
            pi0_decay = decay_to_2_body(self.pi0_boosts, m_exo, c.m_pi0, 0)
            if pi0_decay.success: 
                kdes.append(pi0_decay.kde())
                weights.append(self.multiplicities["Pi0"]*f.rescale_p_2gamma(m_exo,c.m_pi0)*c.BR_Pi0_2Gamma)
        if(m_exo<c.m_eta):
            eta_decay = decay_to_2_body(self.eta_boosts, m_exo, c.m_eta, 0)
            if eta_decay.success: 
                kdes.append(eta_decay.kde())
                weights.append(self.multiplicities["Eta"]*f.rescale_p_2gamma(m_exo,c.m_eta)*c.BR_Eta_2Gamma)
        if(m_exo<c.m_etap):
            etap_decay = decay_to_2_body(self.etap_boosts, m_exo, c.m_etap, 0)
            if etap_decay.success: 
                kdes.append(etap_decay.kde())
                weights.append(self.multiplicities["EtaP"]*f.rescale_p_2gamma(m_exo,c.m_etap)*c.BR_EtaPrim_2Gamma)

        #vector decays
        if(m_exo+c.m_pi0<c.m_rho):
            rho0_decay = decay_to_2_body(self.rho0_boosts, m_exo, c.m_rho, c.m_pi0)
            if rho0_decay.success: 
                kdes.append(rho0_decay.kde())
                weights.append(self.multiplicities["Rho0"]*f.rescale_v_p_gamma(m_exo,c.m_rho,c.m_pi0)*c.BR_Rho0_Pi0Gamma)
        if(m_exo+c.m_pi<c.m_rho):
            rho_decay = decay_to_2_body(self.rho_boosts, m_exo, c.m_rho, c.m_pi)
            if rho_decay.success: 
                kdes.append(rho_decay.kde())
                weights.append(self.multiplicities["Rho"]*f.rescale_v_p_gamma(m_exo,c.m_rho,c.m_pi)*c.BR_Rho_PiGamma)
        if(m_exo+c.m_pi0<c.m_omega):
            omega_decay = decay_to_2_body(self.omega_boosts, m_exo, c.m_omega, c.m_pi0)
            if omega_decay.success: 
                kdes.append(omega_decay.kde())
                weights.append(self.multiplicities["Omega"]*f.rescale_v_p_gamma(m_exo,c.m_omega,c.m_pi0)*c.BR_Omega_Pi0Gamma)
        if(m_exo+c.m_pi0<c.m_phi):
            phi_decay = decay_to_2_body(self.phi_boosts, m_exo, c.m_phi, c.m_pi0)
            if phi_decay.success: 
                kdes.append(phi_decay.kde())
                weights.append(self.multiplicities["Phi"]*f.rescale_v_p_gamma(m_exo,c.m_phi,c.m_pi0)*c.BR_Phi_Pi0Gamma)
        if(m_exo+c.m_eta<c.m_rho):
            rho0_decay = decay_to_2_body(self.rho0_boosts, m_exo, c.m_rho, c.m_eta)
            if rho0_decay.success: 
                kdes.append(rho0_decay.kde())
                weights.append(self.multiplicities["Rho0"]*f.rescale_v_p_gamma(m_exo,c.m_rho,c.m_eta)*c.BR_Rho0_EtaGamma)
        if(m_exo+c.m_eta<c.m_omega):
            omega_decay = decay_to_2_body(self.omega_boosts, m_exo, c.m_omega, c.m_eta)
            if omega_decay.success: 
                kdes.append(omega_decay.kde())
                weights.append(self.multiplicities["Omega"]*f.rescale_v_p_gamma(m_exo,c.m_omega,c.m_eta)*c.BR_Omega_EtaGamma)
        if(m_exo+c.m_eta<c.m_phi):
            phi_decay = decay_to_2_body(self.phi_boosts, m_exo, c.m_phi, c.m_eta)
            if phi_decay.success: 
                kdes.append(phi_decay.kde())
                weights.append(self.multiplicities["Phi"]*f.rescale_v_p_gamma(m_exo,c.m_phi,c.m_eta)*c.BR_Phi_EtaGamma)

        data_list = []
        for en_exo in self.en_list:
            data_sublist = []
            if en_exo > m_exo:
                for th_exo in self.th_list:
                    norm_yield = 0.
                    for weight, kde in zip(weights, kdes):
                        if weight==0.: continue
                        norm_yield += f.f_kde(kde, weight, en_exo, th_exo)[0]
                    data_sublist.append([th_exo, en_exo, m_exo, norm_yield])
            else:
                for th_exo in self.th_list:
                    data_sublist.append([th_exo, en_exo, m_exo, 0.])
            data_list.append([data_sublist])
        return data_list
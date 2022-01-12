import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d
from os import path
import mpmath as mp #for polylog
import alp_setup as setup
import alp_constants as c
import alp_functions as f
import decay_widths as width
import effective_coupling as eff

class Load_data:
    constraint_dictionary = {}
    boundary_dictionary = {}

    def __init__(self,exp,channels_decay,channels_production):
        for chan_dec in channels_decay:
            for chan_prod in channels_production:

                filename_dat = path.dirname(path.realpath(__file__))+"/../tab_decay/"+exp+"/"+exp+'_'+chan_prod+'_'+chan_dec+'.dat'

                if path.exists(filename_dat):
                    
                    experimental_constraint_data_dat = np.loadtxt(filename_dat)
        
                    experimental_constraint_data = np.delete(experimental_constraint_data_dat.reshape((201,101,3)),100,0)
        
                    # Extract the boundaries of the tabulated grid
                    self.boundary_dictionary[exp+'_'+chan_prod+'_'+chan_dec] = np.array([[experimental_constraint_data[0,0,0],experimental_constraint_data[-1,0,0]],[experimental_constraint_data[0,0,1],experimental_constraint_data[0,-1,1]]])
        
                    # Add a small number to avoid taking the logarithm of zero 
                    experimental_constraint_data = experimental_constraint_data[:,:,:] + [0,0,c.epsilon]
                    # Take logarithm to make interpolation easier
                    experimental_constraint_data = np.log(experimental_constraint_data)
                    # Fast interpolation on rectangular grid
                    experimental_constraint_data_inter = RectBivariateSpline(experimental_constraint_data[:,0,0],experimental_constraint_data[0,:,1],experimental_constraint_data[:,:,2], kx=1, ky=1)
        
                    self.constraint_dictionary[exp+'_'+chan_prod+'_'+chan_dec] = experimental_constraint_data_inter

                else:

                    print(filename_dat,' not found')
        
                    # If no file exists, we define the boundaries in such a way that the channel will be skipped in the calculations below
                    self.boundary_dictionary[exp+'_'+chan_prod+'_'+chan_dec] = np.array([[0, -1],[0,-1]])

class Process_data:
    def __init__(self,experiment,channels_decay,channels_production):
        self._channels_decay = channels_decay
        self._channels_production = channels_production
        self._exp = experiment        

        self._data = Load_data(self._exp,self._channels_decay,self._channels_production)

        self._digi_widths = {}
        for channel in self._channels_decay:
            self._digi_widths[channel] = width.Load_digitized(channel)
 
        self._tot = width.Total_width()

        self._processed = 0.

    def _ALP_decays_single_channel(self, production_channel, decay_channel, Gamma_a, m_a):

        boundary = self._data.boundary_dictionary[self._exp+'_'+production_channel+'_'+decay_channel]

        g_a = self._coupling_decay[decay_channel]
        BR_a = 0
        if self._digi_widths[decay_channel].digitized_exists:
            BR_a = self._digi_widths[decay_channel].width_from_file(m_a,g_a) / Gamma_a
        else:
            BR_a = width.analytical_widths[decay_channel](m_a,g_a) / Gamma_a

        # Check if the requested value of m_a and Gamma_a lie within the tabulated range. Otherwise return zero.
        if boundary[0,0] <= m_a <= boundary[0,1] and boundary[1,0] <= Gamma_a <= boundary[1,1]:
            return BR_a * (self._coupling_production[production_channel] / setup.reference_couplings[production_channel])**setup.scaling_exponent[production_channel] * (np.exp(self._data.constraint_dictionary[self._exp+'_'+production_channel+'_'+decay_channel](np.log(float(m_a)),np.log(float(Gamma_a)))[0,0]) - c.epsilon)
        else:
            return 0

    def _ALP_events(self, m_a, Gamma_a):

        number_of_decays =  np.sum([
                                np.sum([
                                    self._ALP_decays_single_channel(self._channels_production[i], self._channels_decay[j], Gamma_a, m_a)
                                    for i in range(len(self._channels_production))
                                ])
                                for j in range(len(self._channels_decay))
                            ])


        if number_of_decays < 0: number_of_decays = 0.

        return number_of_decays

    def ALP_events_EFT(self, m_a, C_GG, C_WW, C_BB, Lambda, AA, BB, C_ll):
#        self._processed += 1./(201*896)
        self._processed += 1./22725
        print("\r" + " processed for " + self._exp + ": " + "{:.2f}".format(self._processed*100) + "%", end="             ")

        #eff. photon coupling:
        ph = eff.photon_coupling(Lambda)
        g_gg = ph.g_gg_eff(m_a, C_GG, C_WW, C_BB, C_ll)

        #eff. lepton coupling
        lp = eff.lepton_coupling(Lambda)
        g_ee = lp.g_ee_eff(C_WW, C_BB, C_ll)
        g_mumu = lp.g_mumu_eff(C_WW, C_BB, C_ll)

        #ALP-meson mixing:
        th = eff.mixing(Lambda)
        th_pi = th.pi(C_GG,m_a)
        th_eta = th.eta(C_GG,m_a)
        th_etap = th.etap(C_GG,m_a)

        #B decay branching fraction:
        bs = eff.bs_coupling(Lambda, AA, BB)

        BR_B_K_a = width.B_K_a(m_a,bs.g_bs_eff(m_a,C_GG,C_WW)) / c.Gamma_B
        BR_B_Kstar_a = width.B_Kstar_a(m_a,bs.g_bs_eff(m_a,C_GG,C_WW)) / c.Gamma_B

        #D decay branching fraction:

        cu = eff.cu_coupling(Lambda, AA, BB)

        BR_D_Pi_a = width.D_pi_a(m_a,cu.g_cu_eff(m_a,C_GG,C_WW)) / c.Gamma_D

        g_GG = C_GG/Lambda

        self._coupling_production={ 'primakoff':         g_gg,
                                    'photonfrommeson':   g_gg,
                                    'mixingPi0':         th_pi,
                                    'mixingEta':         th_eta,
                                    'mixingEtaPrim':     th_etap,
                                    'BmesonK':           BR_B_K_a,
                                    'BmesonKstar':       BR_B_Kstar_a,
                                    'DmesonPi':          BR_D_Pi_a}        

        self._coupling_decay = {'2Gamma':       g_gg,
                                '2El':          g_ee,
                                '2Mu':          g_mumu,
                                '3Pi0':         g_GG,
                                '3Pi':          g_GG,
                                '2PiGamma':     g_GG,
                                '2Pi0Eta':      g_GG,
                                '2PiEta':       g_GG,
                                '2Pi0EtaPrim':  g_GG,
                                '2PiEtaPrim':   g_GG}
        

        return self._ALP_events(m_a, self._tot.width_a(m_a,g_gg, g_GG, g_ee, g_mumu))
import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d
from os import path
import mpmath as mp #for polylog
from . import alp_setup as setup
from scalar import scalar_setup as setds
from scalar import decay_widths_DS as wds
from general import alp_constants as c
from general import alp_functions as f
from . import decay_widths as width
from . import effective_coupling as eff

class Load_data:
    constraint_dictionary = {}
    boundary_dictionary = {}

    def __init__(self,exp,channels_decay,channels_production):
        for chan_dec in channels_decay:
            for chan_prod in channels_production:

                filename_dat = path.dirname(path.realpath(__file__))+"/../../tab_decay/"+setup.experiments[exp]+"/"+exp+'_'+chan_prod+'_'+chan_dec+'.dat'

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
                    print('[Warning:] \t',filename_dat,'not found')
        
                    # If no file exists, we define the boundaries in such a way that the channel will be skipped in the calculations below
                    self.boundary_dictionary[exp+'_'+chan_prod+'_'+chan_dec] = np.array([[0, -1],[0,-1]])

class Process_data:
    def __init__(self,experiment,channels_decay,channels_production, nbins, exo_particle):
        self._channels_decay = channels_decay
        self._channels_production = channels_production
        self._exp = experiment
        self._exotic = exo_particle

        self._data = Load_data(self._exp,self._channels_decay,self._channels_production)

        self._coupling_production = {}
        self._digi_widths = {}
        self._tot = None
        if self._exotic == "DS":
            all_prod_channels = setds.channels_production
            self._ds_widths = wds.DS_widths()
        else:
            all_prod_channels = setup.channels_production
            for channel in self._channels_decay:
                self._digi_widths[channel] = width.Load_digitized(channel)
            self._tot = width.Total_width()

        for mode in all_prod_channels: #initialize with 1
            self._coupling_production[mode] = self._coupling_production.get(mode,0) + 1
  
        #for progress bar
        self._nbins = nbins
        self._processed = 0.
        self._progress = 0

        print("[Info:] \t Running for experiment:", self._exp)

    def Set_Coupling_Production(self,channel,value): #useful e.g. for scanning over production coupling or fixing the value
        self._coupling_production[channel] = value

    def n_events_single_channel_BR(self, production_channel, decay_channel, BR, Gamma, mX):
        boundary = self._data.boundary_dictionary[self._exp+'_'+production_channel+'_'+decay_channel]
        # Check if the requested value of mX and Gamma lie within the tabulated range. Otherwise return zero.
        if boundary[0,0] <= mX <= boundary[0,1] and boundary[1,0] <= Gamma <= boundary[1,1]:
            return BR * (self._coupling_production[production_channel] / setup.reference_couplings[production_channel])**setup.scaling_exponent[production_channel] * (np.exp(self._data.constraint_dictionary[self._exp+'_'+production_channel+'_'+decay_channel](np.log(float(mX)),np.log(float(Gamma)))[0,0]) - c.epsilon)
        else:
            return 0

    def n_events_Bmeson_BR(self, production_channels, decay_channels, BRprod, BRdec, Gamma, mX):

        number_of_decays =  0.
        for chan in production_channels:
            if "meson" not in chan:
                raise ValueError('Unexpected production channel ' + chan)
            BRKexo = BRprod
            if len(production_channels) != 1 and chan == "BmesonKstar": #when summing over both BmesonK and BmesonKstar modify correctly the B->K*a BR
                BRKexo = BRprod * 1 / (c.m_B**2 - c.m_K**2)**2 * ((1.364/(1-mX**2/(c.m_B**2)) - 0.990/(1-mX**2/36.78))/(0.330/(1-mX**2/37.46)))**2 * (f.lambda_Kallen(c.m_B,c.m_Kstar,mX)**(3/2))/f.lambda_Kallen(c.m_B,c.m_K,mX)**(1/2)
            self.Set_Coupling_Production(chan,BRKexo)
            number_of_decays += self.n_events_single_channel_BR(chan, decay_channels[0], BRdec, Gamma, mX)

        if number_of_decays < 0: number_of_decays = 0.

        return number_of_decays

    def _n_decays_single_channel(self, production_channel, decay_channel, Gamma, mX):

        # boundary = self._data.boundary_dictionary[self._exp+'_'+production_channel+'_'+decay_channel]

        BR = 0
        if self._exotic == "DS":
            BR = self._ds_widths.get_width(decay_channel,mX,self._coupling_decay[decay_channel]) / Gamma
        else:
            if self._digi_widths[decay_channel].digitized_exists:
                BR = self._digi_widths[decay_channel].width_from_file(mX,self._coupling_decay[decay_channel]) / Gamma
            else:
                BR = width.analytical_widths[decay_channel](mX,self._coupling_decay[decay_channel]) / Gamma

        return self.n_events_single_channel_BR(production_channel, decay_channel, BR, Gamma, mX)

    def _n_events(self, mX, Gamma):

        number_of_decays =  np.sum([
                                np.sum([
                                    self._n_decays_single_channel(self._channels_production[i], self._channels_decay[j], Gamma, mX)
                                    for i in range(len(self._channels_production))
                                ])
                                for j in range(len(self._channels_decay))
                            ])


        if number_of_decays < 0: number_of_decays = 0.

        return number_of_decays

    def ALP_events_EFT(self, m_a, C_GG, C_WW, C_BB, Lambda, AA, BB, C_ll, C_qq):
        #progress bar
        self._processed += 1./self._nbins
        if self._processed*20 >= self._progress:
            print("\r" + "[Info:] \t Processing [" + "-"*self._progress + " "*(20-self._progress) + "]", end="             ")
            self._progress += 1

        #eff. photon coupling:
        ph = eff.photon_coupling(Lambda)
        g_gg = ph.g_gg_eff(m_a, C_GG, C_WW, C_BB, C_ll, C_qq)

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

        # BR_B_K_a = width.B_K_a(m_a,bs.g_bs_eff(m_a,C_GG,C_WW)) / c.Gamma_B
        # BR_B_Kstar_a = width.B_Kstar_a(m_a,bs.g_bs_eff(m_a,C_GG,C_WW)) / c.Gamma_B

        #D decay branching fraction:

        cu = eff.cu_coupling(Lambda, AA, BB)

        # BR_D_Pi_a = width.D_pi_a(m_a,cu.g_cu_eff(m_a,C_GG,C_WW)) / c.Gamma_D
        BR_D_Pi_a = 0

        #effective gluon coupling for hadronic decays
        gl = eff.gluon_coupling(Lambda)
        g_GG = gl.g_GG_eff(m_a, C_GG, C_qq)

        g_bs = 0
        if(C_GG != 0 and C_WW != 0):
            g_bs = bs.g_bs_fixed_all(C_GG, C_WW, C_BB, C_ll, C_qq)
        elif(C_GG != 0):
            g_bs = bs.g_bs_fixed_cgg(C_GG)
        elif(C_WW != 0):
            g_bs = bs.g_bs_eff(m_a,C_GG,C_WW)

        BR_B_K_a = width.B_K_a(m_a,g_bs) / c.Gamma_B
        BR_B_Kstar_a = width.B_Kstar_a(m_a,g_bs) / c.Gamma_B

        self._coupling_production={ 'primakoff':         g_gg,
                                    'photonfrommeson':   g_gg,
                                    'mixingPi0':         th_pi,
                                    'mixingEta':         th_eta,
                                    'mixingEtaPrim':     th_etap,
                                    'BmesonK':           BR_B_K_a,
                                    'BmesonKstar':       BR_B_Kstar_a,
                                    'DmesonPi':          BR_D_Pi_a,
                                    'recast':            g_gg}        

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
        

        return self._n_events(m_a, self._tot.width_a(m_a,g_gg, g_GG, g_ee, g_mumu))


    def DS_events_EFT(self, m_S, Y): 
        #progress bar
        self._processed += 1./self._nbins
        if self._processed*20 >= self._progress:
            print("\r" + "[Info:] \t Processing [" + "-"*self._progress + " "*(20-self._progress) + "]", end="             ")
            self._progress += 1

        #define B decay branching fraction
        self._coupling_production = {   'BmesonK':      wds.B_K_S(m_S,Y) / c.Gamma_B,
                                        'BmesonKstar':  wds.B_Kstar_S(m_S,Y) / c.Gamma_B}

        self._coupling_decay = {'2El':  Y,
                                '2Mu':  Y,
                                '2Pi':  Y,
                                '2K':   Y
                                }
        return self._n_events(m_S, self._ds_widths.S_total(m_S,Y))

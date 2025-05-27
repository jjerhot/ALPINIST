import numpy as np
from scipy.interpolate import RectBivariateSpline
from os import path
import mpmath as mp #for polylog

from ALP_rescale.scalar import scalar_setup as setds
from ALP_rescale.scalar import decay_widths_DS as wds

from ALP_rescale.alp import alp_setup as setalp
from ALP_rescale.alp import decay_widths as wda
from ALP_rescale.alp import effective_coupling as eff

from ALP_rescale.vector import vector_setup as setdp
from ALP_rescale.vector import decay_widths_DP as wdp
from ALP_rescale.vector import effective_coupling as eff_dp

from ALP_rescale.hnl import hnl_setup as sethnl
from ALP_rescale.hnl import decay_widths as wdh

from ALP_rescale.general import constants as c
from ALP_rescale.general import functions as f
from ALP_rescale.general import setup as s

from ALP_production.modes.bremsstrahlung import brems_production

class Load_data:
    constraint_dictionary = {}
    boundary_dictionary = {}

    def __init__(self,exo,exp,channels_decay,channels_production,prod_scale_sigmas=0):
        for chan_dec in channels_decay:
            for chan_prod in channels_production:
                
                filename_dat = path.dirname(path.realpath(__file__))+"/../../tab_decay/"+s.experiments[exp]+"/"+exo+"/"+exp+'_'+chan_prod+'_'+chan_dec+'.dat'

                if path.exists(filename_dat):
                    
                    experimental_constraint_data_dat = np.loadtxt(filename_dat)
                    if prod_scale_sigmas: 
                        exp_key = s.experiments[exp]
                        naked_chan_prod = chan_prod.split('-')[0]
                        scaling = max(1 + prod_scale_sigmas/np.divide(*c.prod_xsec_ref[naked_chan_prod][exp_key]),0) if c.prod_xsec_ref[naked_chan_prod][exp_key][0] else 0
                        experimental_constraint_data_dat[:,2] *= scaling

                    reshape_dim_0 =  int(experimental_constraint_data_dat.size / 3 / 101)
                    if reshape_dim_0 != 201: print(f'[Warning:] \texptected data at {filename_dat} of shape (20301,3), found ({101*reshape_dim_0},3) instead. \n\t\t\tMake sure to use current reference grids.')

                    experimental_constraint_data = np.delete(experimental_constraint_data_dat.reshape((reshape_dim_0,101,3)),100,0)
        
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
    def __init__(self,experiment,channels_decay,channels_production, nbins, exo_particle,is_quark_coupling=False, run_standalone = False, dirac_width = False, alp_lambda_UV = 1000, prod_scale_sigmas=0):
        self.ds_brems = brems_production(experiment,exo_particle)
        self._channels_decay = channels_decay
        self._channels_production = channels_production
        self._exp = experiment
        self._exotic = exo_particle
        self._is_quark_coupling = is_quark_coupling
        self.dirac_width = dirac_width
        self._run_standalone = run_standalone

        self._data = Load_data(self._exotic,self._exp,self._channels_decay,self._channels_production,prod_scale_sigmas)

        self._coupling_production = {}
        self._reference_couplings = {}
        self._scaling_exponent = {}
        self._digi_widths = {}
        self._tot = None
        self._mixing_dp = None
        if self._exotic == "ds":
            all_prod_channels = setds.channels_production
            self._ds_widths = wds.DS_widths()
            self._scaling_exponent = 2
        elif self._exotic == "hnl":
            all_prod_channels = sethnl.channels_production
            self._tot = wdh.Total_width()
            self._scaling_exponent = 1
        elif self._exotic == "dp":
            all_prod_channels = setdp.channels_production
            self._dp_widths = wdp.DP_widths()
            self._scaling_exponent = setdp.scaling_exponent
            self._reference_couplings = setdp.reference_couplings
            self._mixing_dp = eff_dp.mixing(1.,1.,1.)
        else:
            all_prod_channels = setalp.channels_production
            for channel in self._channels_decay:
                if is_quark_coupling:
                    self._digi_widths[channel] = wda.Load_digitized_fermions(channel)
                else:
                    self._digi_widths[channel] = wda.Load_digitized(channel)
            self._tot = wda.Total_width(is_quark_coupling)
            self._scaling_exponent = setalp.scaling_exponent
            self._reference_couplings = setalp.reference_couplings
            
            self._lambda = alp_lambda_UV
            self._lp = eff.lepton_coupling(self._lambda)
            self._ph = eff.photon_coupling(self._lambda)
            self._gl = eff.gluon_coupling(self._lambda)
            self._th = eff.mixing(self._lambda)
            
            C_f=[0.5,0.5,0.5]
            C_F=[-0.5,-0.5,-0.5]
            above = eff.above_EW(self._lambda)
            self._cu = above.C_uu_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1,1,"u")
            self._cd = above.C_dd_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1,1,"d")
            self._cs = above.C_dd_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1,1,"s")
            self._cc = above.C_uu_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1,1,"c")
            self._cb = above.C_dd_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1,1,"b")
            self._ct = above.C_uu_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1,1,"t")
            self._cel = above.C_ee_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1,1)
            self._cmu = above.C_ee_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1,1)
            self._ctau = above.C_ee_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1,1)
            
            self._cbs_G = above.C_qq_ij_mu(c.m_t,1,0,0,0, C_qq_ij = 0, i = "b", j = "s")
            self._cbs_W = above.C_qq_ij_mu(c.m_t,0,1,0,0, C_qq_ij = 0, i = "b", j = "s")
            self._cbs_B = above.C_qq_ij_mu(c.m_t,0,0,1,0, C_qq_ij = 0, i = "b", j = "s")
            self._cbs_t = above.C_qq_ij_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1, C_qq_ij = 0, i = "b", j = "s")
            
            self._csd_G = above.C_qq_ij_mu(c.m_t,1,0,0,0, C_qq_ij = 0, i = "s", j = "d")
            self._csd_W = above.C_qq_ij_mu(c.m_t,0,1,0,0, C_qq_ij = 0, i = "s", j = "d")
            self._csd_B = above.C_qq_ij_mu(c.m_t,0,0,1,0, C_qq_ij = 0, i = "s", j = "d")
            self._csd_t = above.C_qq_ij_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1, C_qq_ij = 0, i = "s", j = "d")
            
            self._ccu_G = above.C_qq_ij_mu(c.m_t,1,0,0,0, C_qq_ij = 0, i = "c", j = "u")
            self._ccu_W = above.C_qq_ij_mu(c.m_t,0,1,0,0, C_qq_ij = 0, i = "c", j = "u")
            self._ccu_B = above.C_qq_ij_mu(c.m_t,0,0,1,0, C_qq_ij = 0, i = "c", j = "u")
            self._ccu_t = above.C_qq_ij_mu(c.m_t,above.C_GG_tilde(0,C_u=C_f,C_d=C_f,C_Q=C_F),above.C_WW_tilde(0,C_Q=C_F,C_L=C_F),above.C_BB_tilde(0,C_Q=C_F,C_L=C_F,C_e=C_f,C_u=C_f,C_d=C_f),1, C_qq_ij = 0, i = "c", j = "u")

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
        scaling_exponent = self._scaling_exponent[production_channel] if isinstance( self._scaling_exponent,dict) else self._scaling_exponent
        # Check if the requested value of mX and Gamma lie within the tabulated range. Otherwise return zero.
        if boundary[0,0] <= mX <= boundary[0,1] and boundary[1,0] <= Gamma <= boundary[1,1]:
            # return BR * (np.exp(self._data.constraint_dictionary[self._exp+'_'+production_channel+'_'+decay_channel](np.log(float(mX)),np.log(float(Gamma)))[0,0]) - c.epsilon)
            if self._run_standalone:
               return np.exp(self._data.constraint_dictionary[self._exp+'_'+production_channel+'_'+decay_channel](np.log(float(mX)),np.log(float(Gamma)))[0,0]) - c.epsilon
            #.. to be unified later ( remove self.rescale_exponent[production_channel])
            if self._exotic == "alp": return BR * (self._coupling_production[production_channel]/self._reference_couplings[production_channel]) ** self._scaling_exponent[production_channel] * (np.exp(self._data.constraint_dictionary[self._exp+'_'+production_channel+'_'+decay_channel](np.log(float(mX)),np.log(float(Gamma)))[0,0]) - c.epsilon)

            if self._exotic == "dp": 
                if production_channel == "Brems": #apply VMD form factor for BB brems production
                    return BR * self.ds_brems.form_factor_vector(float(mX))**2 * (self._coupling_production[production_channel]) ** scaling_exponent * (np.exp(self._data.constraint_dictionary[self._exp+'_'+production_channel+'_'+decay_channel](np.log(float(mX)),np.log(float(Gamma)))[0,0]) - c.epsilon) 
                return BR * (self._coupling_production[production_channel]) ** scaling_exponent * (np.exp(self._data.constraint_dictionary[self._exp+'_'+production_channel+'_'+decay_channel](np.log(float(mX)),np.log(float(Gamma)))[0,0]) - c.epsilon) 

            return BR * (self._coupling_production[production_channel]) ** scaling_exponent * (np.exp(self._data.constraint_dictionary[self._exp+'_'+production_channel+'_'+decay_channel](np.log(float(mX)),np.log(float(Gamma)))[0,0]) - c.epsilon) 
        else:
            return 0

    def n_events_Bmeson_BR(self, production_channels, decay_channels, BRprod, BRdec, Gamma, mX):

        number_of_decays =  0.
        for chan in production_channels:
            if "meson" not in chan:
                raise ValueError('Unexpected production channel ' + chan)
            BRKexo = BRprod
            if len(production_channels) != 1 and chan == "BmesonKstar": #when summing over both BmesonK and BmesonKstar modify correctly the B->K*a BR
                BRKexo = BRprod / (c.m_B**2 - c.m_K**2)**2  * ((1.364/(1-mX**2/(c.m_B**2)) - 0.990/(1-mX**2/36.78))/(0.330/(1-mX**2/37.46)))**2 * (f.lambda_Kallen(c.m_B,c.m_Kstar,mX)**(3/2))/f.lambda_Kallen(c.m_B,c.m_K,mX)**(1/2)
            self.Set_Coupling_Production(chan,BRKexo)
            number_of_decays += self.n_events_single_channel_BR(chan, decay_channels[0], BRdec, Gamma, mX)

        if number_of_decays < 0: number_of_decays = 0.

        return number_of_decays

    def _n_decays_single_channel(self, production_channel, decay_channel, Gamma, mX):
        BR = 0
        if Gamma != 0:
            if self._exotic == "ds":
                BR = self._ds_widths.get_width(decay_channel,mX,self._coupling_decay) / Gamma
            elif self._exotic == "dp":
                BR = self._dp_widths.get_width(decay_channel,mX,self._coupling_decay) / Gamma
            elif self._exotic=="hnl":
                BR = wdh.hnl_decay_width.widths(decay_channel,mX,self._coupling_decay, self.dirac_width) / Gamma
            else:
                if self._digi_widths[decay_channel].digitized_exists:
                    BR = self._digi_widths[decay_channel].width_from_file(mX,self._coupling_decay[decay_channel]) / Gamma
                else:
                    BR = wda.analytical_widths[decay_channel](mX,self._coupling_decay[decay_channel]) / Gamma
        if BR > 1: BR = 1.
        if BR < 0 or mp.isnan(BR): BR = 0.
        return self.n_events_single_channel_BR(production_channel, decay_channel, BR, Gamma, mX)
    ## this could profit from pool maybe?
    def _n_events(self, mX, Gamma):
        number_of_decays = 0.
        for decay_channel in self._channels_decay:
            for production_channel in self._channels_production:
                number_of_decays += self._n_decays_single_channel(production_channel, decay_channel, Gamma, mX)

        if number_of_decays < 0: number_of_decays = 0.
        return number_of_decays

    def ALP_events_EFT(self, m_a, C_GG, C_WW, C_BB, C_ll, C_qq):
        #progress bar
        self._processed += 1./self._nbins
        if self._processed*20 >= self._progress:
            print("\r" + "[Info:] \t Processing [" + "-"*self._progress + " "*(20-self._progress) + "]", end="             ")
            self._progress += 1

        #eff. photon coupling:
        g_gg = self._ph.g_gg_eff(m_a, C_GG, C_WW, C_BB, ct=self._ct*C_qq, cel=self._cel*C_ll, cb=self._cb*C_qq, cc=self._cc*C_qq, cs=self._cs*C_qq, cu=self._cu*C_qq, cd=self._cd*C_qq, cmu=self._cmu*C_ll, ctau=self._ctau*C_ll)

        #eff. lepton coupling
        g_ee = self._lp.g_ee_eff(C_WW, C_BB, self._cel*C_ll)
        g_mumu = self._lp.g_mumu_eff(C_WW, C_BB, self._cmu*C_ll)

        #ALP-meson mixing:
        th_pi = self._th.pi(m_a,C_GG, cg=C_BB+C_WW, cu=self._cu*C_qq, cd=self._cd*C_qq, cs=self._cs*C_qq, cb=self._cb*C_qq, ct=self._ct*C_qq)
        th_eta = self._th.eta(m_a,C_GG, cg=C_BB+C_WW, cu=self._cu*C_qq, cd=self._cd*C_qq, cs=self._cs*C_qq, cb=self._cb*C_qq, ct=self._ct*C_qq)
        th_etap = self._th.etap(m_a,C_GG, cg=C_BB+C_WW, cu=self._cu*C_qq, cd=self._cd*C_qq, cs=self._cs*C_qq, cb=self._cb*C_qq, ct=self._ct*C_qq)

        # FCNC couplings
        g_bs = (self._cbs_G*C_GG + self._cbs_W*C_WW + self._cbs_B*C_BB + self._cbs_t*C_qq)/self._lambda
        g_cu = (self._ccu_G*C_GG + self._ccu_W*C_WW + self._ccu_B*C_BB + self._ccu_t*C_qq)/self._lambda
        # g_sd = (self._csd_G*C_GG + self._csd_W*C_WW + self._csd_B*C_BB + self._csd_t*C_qq)/self._lambda
        
        #effective gluon coupling for hadronic decays
        g_GG = self._gl.g_GG_eff(m_a, C_GG, C_qq)

        #effective quark coupling for hadronic decays
        g_qq = C_qq/self._lambda

        self._coupling_production={ 'Primakoff':         g_gg,
                                    'PhotonFromMeson':   g_gg,
                                    'MixingPi0':         th_pi,
                                    'MixingEta':         th_eta,
                                    'MixingEtaPrim':     th_etap,
                                    'Bmeson':            np.abs(g_bs),
                                    'Dmeson':            np.abs(g_cu),
                                    'recast':            g_gg}        
        if self._is_quark_coupling:
            self._coupling_decay = {'2Gamma':       g_gg,
                                    '2El':          g_ee,
                                    '2Mu':          g_mumu,
                                    '3Pi0':         g_qq,
                                    '3Pi':          g_qq,
                                    '3Eta':         g_qq,
                                    '2PiGamma':     g_qq,
                                    '2Pi0Eta':      g_qq,
                                    '2PiEta':       g_qq,
                                    '2Pi0EtaPrim':  g_qq,
                                    '2PiEtaPrim':   g_qq,
                                    '2PiOmega':     g_qq,
                                    '2KPi0':        g_qq,
                                    '2Pi2Pi0':      g_qq,
                                    '4Pi':          g_qq,
                                    '2Omega':       g_qq,
                                    '2Kstar':       g_qq,
                                    '2Phi':         g_qq
                                    }
        else:
            self._coupling_decay = {'2Gamma':       g_gg,
                                    '2El':          g_ee,
                                    '2Mu':          g_mumu,
                                    '3Pi0':         g_GG,
                                    '3Pi':          g_GG,
                                    '2PiGamma':     g_GG,
                                    '2Pi0Eta':      g_GG,
                                    '2PiEta':       g_GG,
                                    '2Pi0EtaPrim':  g_GG,
                                    '2PiEtaPrim':   g_GG,
                                    '2KPi0':        g_GG,
                                    '2Pi2Pi0':      g_GG,
                                    '4Pi':          g_GG,
                                    '2Omega':       g_GG,
                                    '2Kstar':       g_GG,
                                    '2Phi':         g_GG
                                    }
        return self._n_events(m_a, self._tot.width_a(m_a, g_gg, g_GG, g_ee, g_mumu, g_qq))


    def DS_events_EFT(self, m_S, Y, LambdaS): 
        #progress bar
        self._processed += 1./self._nbins
        if self._processed*20 >= self._progress:
            print("\r" + "[Info:] \t Processing [" + "-"*self._progress + " "*(20-self._progress) + "]", end="             ")
            self._progress += 1

        #define B decay branching fraction
        self._coupling_production = {   'Brems':    Y,
                                        'Bmeson':    Y,
                                        'Bmeson2S':  LambdaS}

        self._coupling_decay = Y
        return self._n_events(m_S, self._ds_widths.S_total(m_S,Y))

    def HNL_events(self, m_N, U2el, U2mu, U2tau):
        #progress bar
        self._processed += 1./self._nbins
        if self._processed*20 >= self._progress:
            print("\r" + "[Info:] \t Processing [" + "-"*self._progress + " "*(20-self._progress) + "]", end="             ")
            self._progress += 1

        self._coupling_production={ 'Bmeson-ElMixing':      U2el,
                                    'Bmeson-MuMixing':      U2mu,
                                    'Bmeson-TauMixing':     U2tau,
                                    
                                    'Dmeson-ElMixing':      U2el,
                                    'Dmeson-MuMixing':      U2mu,
                                    'Dmeson-TauMixing':     U2tau,
                                    } #consider separating in the future
        self._coupling_decay = [U2el, U2mu, U2tau]
        return self._n_events(m_N, self._tot.width_N(m_N,[U2el, U2mu, U2tau], self.dirac_width, only_charged_current=False))
    
    def DP_events_EFT(self, m_dp, epsilon): 
        #progress bar
        self._processed += 1./self._nbins
        if self._processed*20 >= self._progress:
            print("\r" + "[Info:] \t Processing [" + "-"*self._progress + " "*(20-self._progress) + "]", end="             ")
            self._progress += 1
            
        th_rho = epsilon**2*self._mixing_dp.get_mixing_rho(m_dp)
        th_omega = epsilon**2*self._mixing_dp.get_mixing_omega(m_dp)
        th_phi = epsilon**2*self._mixing_dp.get_mixing_phi(m_dp)
        
        self._coupling_production = {   'Brems':            epsilon,
                                        'MesonDecay':       epsilon,
                                        'MixingRho':        th_rho,
                                        'MixingOmega':      th_omega,
                                        'MixingPhi':        th_phi}

        self._coupling_decay = epsilon
        return self._n_events(m_dp, self._dp_widths.dp_total(m_dp,epsilon))

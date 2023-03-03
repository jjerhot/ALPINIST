import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d
from os import path
import mpmath as mp
import alp_setup as setup
import alp_constants as c
import decay_widths as width
import argparse

# Derived from load_data.py for cross-check of B meson mode with 2mu decay

parser = argparse.ArgumentParser(description='ALP MC rescaling for C_ff coupling. \n Select the experiment')
parser.add_argument("-e","--exp",  default="", type=str, nargs='*', help="Experiments available (case sensitive): exp = NA62 | CHARM | nuCAL | SHiP | DarkQuest | DUNE | SHADOWS. If not specified, running over all experiments available.")
parser.add_argument("-p","--prod",  default="", type=str, nargs='*', help="Production modes available (case sensitive): prod = BmesonK | BmesonKstar. If not specified, running over all production modes available.")
parser.add_argument("-d","--decay",  default="", type=str, nargs='*', help="Decay modes available (case sensitive): decay = 2El | 2Mu. If not specified, running over all decay modes available.")
parser.add_argument("-l","--lambda", dest="lam", default=1000, type=float, help="The \u039B [GeV] energy scale. Default value is \u039B = 1000 GeV")

args = parser.parse_args()

#experiment:
experiments = []
if args.exp == "":
    experiments = setup.experiments.keys()
    print("[Info:] \t Selected all experiments available")
else:
    for exp in args.exp:
        if exp in setup.experiments.keys():
            experiments.append(exp)
        else:
            parser.error("[Error:] \t Experiment " + exp + " not available. Experiment modes available: exp = " + ' | '.join(setup.experiments.keys()) + ". If not specified, running over all experiments available.")
    print("[Info:] \t Selected experiments:", ', '.join(experiments))

#production mode:
channels_production_setup = ['BmesonK','BmesonKstar'] #unique for this case
channels_production = []
if args.prod == "":
    channels_production = channels_production_setup
    print("[Info:] \t Selected all production modes available")
else:
    for prod in args.prod:
        if prod in channels_production_setup:
            channels_production.append(prod)
        else:
            parser.error("[Error:] \t Production mode " + prod + " not available. Production modes available: prod = " + ' | '.join(channels_production_setup) + ". If not specified, running over all production modes available.")
    print("[Info:] \t Selected production modes:", ', '.join(channels_production))

#decay mode:
channels_decay_setup = ['2El','2Mu'] #unique for this case
channels_decay = []
if args.decay == "":
    channels_decay = channels_decay_setup
    print("[Info:] \t Selected all decay modes available.")
else:
    for dec in args.decay:
        if dec in channels_decay_setup:
            channels_decay.append(dec)
        else:
            parser.error("[Error:] \t Decay mode " + dec + " not available. Decay modes available: decay = " + ' | '.join(channels_decay_setup) + ". If not specified, running over all decay modes available.")
    print("[Info:] \t Selected decay modes:", ', '.join(channels_decay))

#scale
if args.lam <= 0:
    parser.error("\u039B has to be a positive number.")

#directory
this_dir = path.dirname(path.realpath(__file__))

reference_couplings = {'BmesonK':           1e-10,
                       'BmesonKstar':       1e-10}
scaling_exponent =    {'BmesonK':           1,
                       'BmesonKstar':       1}
coupling_production = {}
coupling_decay = {}

processed = 0.

constraint_dictionary = {}
boundary_dictionary = {}

#total decay width interpolation - only used for hadronic channels. Turned on for m_a > 300 MeV
total_width_digitized = np.loadtxt(path.dirname(path.realpath(__file__))+'/../widths/2mu_integrated/TotalWidth_gY1e-4.dat')
m_a_tot_steps = 1416
m_a_tot_list = np.array([total_width_digitized[i,0] for i in range(m_a_tot_steps)])
Gamma_a_tot_list = np.array([total_width_digitized[i,1] for i in range(m_a_tot_steps)])
Gamma_a_tot_inter = interp1d(m_a_tot_list, Gamma_a_tot_list)

for exp in experiments:
    for chan_dec in channels_decay:
        for chan_prod in channels_production:

            filename_dat = path.dirname(path.realpath(__file__))+"/../tab_decay/"+setup.experiments[exp]+"/"+exp+'_'+chan_prod+'_'+chan_dec+'.dat'

        # filename_dat = path.dirname(path.realpath(__file__))+"/../tab_decay/"+exp+"/"+exp+'_'+chan_prod+'_2mu'+'.dat'

            if path.exists(filename_dat):
                    
                experimental_constraint_data_dat = np.loadtxt(filename_dat)

                experimental_constraint_data = np.delete(experimental_constraint_data_dat.reshape((201,101,3)),100,0)

                # Extract the boundaries of the tabulated grid
                boundary_dictionary[exp+'_'+chan_prod+'_'+chan_dec] = np.array([[experimental_constraint_data[0,0,0],experimental_constraint_data[-1,0,0]],[experimental_constraint_data[0,0,1],experimental_constraint_data[0,-1,1]]])

                # Add a small number to avoid taking the logarithm of zero 
                experimental_constraint_data = experimental_constraint_data[:,:,:] + [0,0,c.epsilon]
                # Take logarithm to make interpolation easier
                experimental_constraint_data = np.log(experimental_constraint_data)
                # Fast interpolation on rectangular grid
                experimental_constraint_data_inter = RectBivariateSpline(experimental_constraint_data[:,0,0],experimental_constraint_data[0,:,1],experimental_constraint_data[:,:,2])

                constraint_dictionary[exp+'_'+chan_prod+'_'+chan_dec] = experimental_constraint_data_inter

            else:
                print('[Warning:] \t',filename_dat,'not found')

                # If no file exists, we define the boundaries in such a way that the channel will be skipped in the calculations below
                boundary_dictionary[exp+'_'+chan_prod+'_'+chan_dec] = np.array([[0, -1],[0,-1]])

def ALP_decays_single_channel(experiment, production_channel, decay_channel, Gamma_a, m_a):

    boundary = boundary_dictionary[experiment+'_'+production_channel+'_'+decay_channel]

    g_a = coupling_decay[decay_channel]
    BR_a = width.analytical_widths[decay_channel](m_a,g_a) / Gamma_a
    if BR_a > 1: BR_a = 1.

    # Check if the requested value of m_a and Gamma_a lie within the tabulated range. Otherwise return zero.
    if boundary[0,0] <= m_a <= boundary[0,1] and boundary[1,0] <= Gamma_a <= boundary[1,1]:

        # return (coupling_production[production_channel] / reference_couplings[production_channel])**scaling_exponent[production_channel] * (np.exp(constraint_dictionary[experiment+'_'+production_channel](np.log(m_a),np.log(Gamma_a))[0,0]) - c.epsilon)
        return BR_a * (coupling_production[production_channel] / setup.reference_couplings[production_channel])**setup.scaling_exponent[production_channel] * (np.exp(constraint_dictionary[experiment+'_'+production_channel+'_'+decay_channel](np.log(float(m_a)),np.log(float(Gamma_a)))[0,0]) - c.epsilon)

    else:

        return 0

# Model-independent part

def ALP_events(experiment, m_a, Gamma_a):

    number_of_decays =  np.sum([
                            np.sum([
                                ALP_decays_single_channel(experiment, channels_production[i], channels_decay[j], Gamma_a, m_a)
                                for i in range(len(channels_production))
                            ])
                            for j in range(len(channels_decay))
                        ])

    if number_of_decays < 0: number_of_decays = 0.

    return number_of_decays

def ALP_events_EFT(experiment, m_a, g_Y, Lambda): 
    global processed 
    processed += 1./144000
    print("\r" + " processed: " + "{:.2f}".format(processed*100) + "%", end="             ")

    #define B decay branching fraction
    V_qb = [c.V_ub, c.V_cb, c.V_tb]
    V_qs = [c.V_us, c.V_cs, c.V_ts]

    h_bs = c.alpha_EM*g_Y*c.m_q[5]**2/(4*np.pi*c.m_W**2*mp.sin(c.theta_w)**2*c.v) * np.log(Lambda**2/c.m_q[5]**2) * sum([np.prod(q) for q in zip(V_qb, V_qs)])

    BR_B_K_a = width.B_K_a(m_a,h_bs) / c.Gamma_B
    BR_B_Kstar_a = width.B_Kstar_a(m_a,h_bs) / c.Gamma_B

    global coupling_production
    coupling_production = { 'BmesonK':      BR_B_K_a,
                            'BmesonKstar':  BR_B_Kstar_a}

    global coupling_decay
    coupling_decay = {  '2El':          g_Y/c.v,
                        '2Mu':          g_Y/c.v}

    return ALP_events(experiment, m_a, (g_Y*np.power(10,4))**2*np.power(10,Gamma_a_tot_inter(np.log10(m_a))))

def ALP_events_exp(expName, Lambda):
    # make lists of masses (2*E-1 to ~2*E+0) and couplings (E-6 to E-2)
    global processed
    processed = 0.
    g_a_list = [ 10**(exponent/100-1) for exponent in range(-500,-100)] 
    m_a_list = [ 2*10**(exponent/300-1) for exponent in range(0,360)] 

    data_list_gY = [[ [m_a, g_Y, ALP_events_EFT(expName, m_a, g_Y, Lambda)] for g_Y in g_a_list] for m_a in m_a_list]
    data_gY = np.reshape(data_list_gY,(len(m_a_list)*len(g_a_list),3))

    # export
    output_dir = this_dir+'/../tab_toPlot/'
    outPath = output_dir + setup.experiments[expName] + '/'

    if args.prod == args.decay == "": modes = ""
    elif args.prod != "" and args.decay == "": modes = "_" + '-'.join(channels_production)
    elif args.decay != "" and args.prod == "": modes = "_" + '-'.join(channels_decay)
    else: modes = "_" + '-'.join(channels_production) + "_" + '-'.join(channels_decay)

    outfileName = expName + '_gY' + modes + '_new.dat'
    np.savetxt(outPath + outfileName,data_gY)
    print('\n[Info:] \t', 'File ' + outfileName + ' saved to ' + outPath)

    return

for exp in experiments:
    ALP_events_exp(exp,args.lam)
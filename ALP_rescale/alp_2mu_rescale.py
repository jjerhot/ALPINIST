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
parser.add_argument("-e","--exp",  required="", type=str, help="Experiments available (case sensitive): exp = NA62 | CHARM | nuCAL | SHiP | DarkQuest | DUNE | SHADOWS. If not specified, running over all experiments available.")
parser.add_argument("-l","--lambda", dest="lam", default=1000, type=float, help="The \u039B [GeV] energy scale. Default value is \u039B = 1000 GeV")

args = parser.parse_args()

#experiment:
if args.exp == "":
    print("Running for all experiments available")
    experiments = setup.experiments
elif args.exp in setup.experiments:
    print("Running for ", args.exp, " experiment")
    experiments = [args.exp]
else:
    parser.error("Experiment " + args.exp + " not available. Experiments available: exp = NA62 | CHARM | nuCAL | SHiP | DarkQuest | DUNE | SHADOWS. If not specified, running over all experiments available.")
#scale
if args.lam <= 0:
    parser.error("\u039B has to be a positive number.")

channels_decay = ['2Mu']
channels_production = ['BmesonK','BmesonKstar']
reference_couplings = {'BmesonK':           1e-10,
                       'BmesonKstar':       1e-10}
scaling_exponent =    {'BmesonK':           1,
                       'BmesonKstar':       1}
coupling_production = {}

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
    for chan_prod in channels_production:

        filename_dat = path.dirname(path.realpath(__file__))+"/../tab_decay/"+exp+"/"+exp+'_'+chan_prod+'_2mu'+'.dat'

        if path.exists(filename_dat):
                
            experimental_constraint_data_dat = np.loadtxt(filename_dat)

            experimental_constraint_data = np.delete(experimental_constraint_data_dat.reshape((201,101,3)),100,0)

            # Extract the boundaries of the tabulated grid
            boundary_dictionary[exp+'_'+chan_prod] = np.array([[experimental_constraint_data[0,0,0],experimental_constraint_data[-1,0,0]],[experimental_constraint_data[0,0,1],experimental_constraint_data[0,-1,1]]])

            # Add a small number to avoid taking the logarithm of zero 
            experimental_constraint_data = experimental_constraint_data[:,:,:] + [0,0,c.epsilon]
            # Take logarithm to make interpolation easier
            experimental_constraint_data = np.log(experimental_constraint_data)
            # Fast interpolation on rectangular grid
            experimental_constraint_data_inter = RectBivariateSpline(experimental_constraint_data[:,0,0],experimental_constraint_data[0,:,1],experimental_constraint_data[:,:,2])

            constraint_dictionary[exp+'_'+chan_prod] = experimental_constraint_data_inter

        else:

            print(filename_dat,' not found')

            # If no file exists, we define the boundaries in such a way that the channel will be skipped in the calculations below
            boundary_dictionary[exp+'_'+chan_prod] = np.array([[0, -1],[0,-1]])

def ALP_decays_single_channel(experiment, production_channel, m_a, Gamma_a):

    boundary = boundary_dictionary[experiment+'_'+production_channel]

    # Check if the requested value of m_a and Gamma_a lie within the tabulated range. Otherwise return zero.
    if boundary[0,0] <= m_a <= boundary[0,1] and boundary[1,0] <= Gamma_a <= boundary[1,1]:

        return (coupling_production[production_channel] / reference_couplings[production_channel])**scaling_exponent[production_channel] * (np.exp(constraint_dictionary[experiment+'_'+production_channel](np.log(m_a),np.log(Gamma_a))[0,0]) - c.epsilon)

    else:

        return 0

# Model-independent part

def ALP_events(experiment, m_a, g_Y):

    Gamma_a = (g_Y*np.power(10,4))**2*np.power(10,Gamma_a_tot_inter(np.log10(m_a)))
    number_of_decays = np.sum([ALP_decays_single_channel(experiment, channels_production[i], m_a, Gamma_a) for i in range(len(channels_production))])

    Gamma_mumu = width.a_2Mu(m_a, g_Y/c.v)
    BR_mumu = Gamma_mumu/Gamma_a

    if BR_mumu > 1.:
        BR_mumu = 1

    if BR_mumu < 0.:
        BR_mumu = 0

    return number_of_decays * BR_mumu

def ALP_events_EFT(experiment, m_a, g_Y, Lambda): 
    global processed 
    processed += 1./48000
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

    return ALP_events(experiment, m_a, g_Y)

def ALP_events_exp(expName, Lambda):
    # make lists of masses (2*E-1 to ~2*E+0) and couplings (E-6 to E-2)
    global processed
    processed = 0.
    g_a_list = [ 10**(exponent/100-1) for exponent in range(-500,-100)] 
    m_a_list = [ 2*10**(exponent/100-1) for exponent in range(0,120)] 

    data_list_gY = [[ [m_a, g_Y, ALP_events_EFT(expName, m_a, g_Y, Lambda)] for g_Y in g_a_list] for m_a in m_a_list]
    data_gY = np.reshape(data_list_gY,(len(m_a_list)*len(g_a_list),3))

    # export
    output_dir = path.dirname(path.realpath(__file__))+'/../tab_toPlot/'
    outPath = output_dir + expName + '/'
    outfileName = expName + '_gY.dat'
    np.savetxt(outPath + outfileName,data_gY)
    print('file ' + outfileName + ' saved to ' + outPath)

    return

for exp in experiments:
    ALP_events_exp(exp,args.lam)
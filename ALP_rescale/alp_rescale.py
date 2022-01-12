from os import path
import numpy as np
import alp_setup as setup
import load_data as ld
import argparse

parser = argparse.ArgumentParser(description='ALP MC rescaling. \n Select the experiment, production and decay modes and other parameters or use default')
parser.add_argument("-e","--exp",  default="", type=str, help="Experiments available (case sensitive): exp = NA62 | CHARM | nuCAL | SHiP | DarkQuest | DUNE | SHADOWS. If not specified, running over all experiments available.")
parser.add_argument("-p","--prod",  default="", type=str, help="Production modes available (case sensitive): prod = primakoff | photonfrommeson | mixingPi0 | mixingEta | mixingEtaPrim | BmesonK | BmesonKstar | DmesonPi. If not specified, running over all production modes available.")
parser.add_argument("-d","--decay",  default="", type=str, help="Decay modes available (case sensitive): decay = 2Gamma | 2El | 2Mu | 3Pi0 | 3Pi | 2PiGamma | 2Pi0Eta | 2PiEta. If not specified, running over all decay modes available.")
parser.add_argument("-l","--lambda", dest="lam", default=1000, type=float, help="The \u039B [GeV] energy scale. Default value is \u039B = 1000 GeV")
parser.add_argument("-cbb", default=1, type=float, help="The C_BB coupling. Default value is C_BB = 1")
parser.add_argument("-cww", default=1, type=float, help="The C_WW coupling. Default value is C_WW = 1")
parser.add_argument("-cgg", default=1, type=float, help="The C_GG coupling. Default value is C_GG = 1")
parser.add_argument("-a", default=0, type=float, help="The model-dependent A parameter in [2102.04474]. Default value is A = 0")
parser.add_argument("-b", default=0, type=float, help="The model-dependent B parameter in [2102.04474]. Default value is B = 0")
parser.add_argument("-cll", default=0, type=float, help="The C_ll coupling. Default value is C_ll = 0")
#parser.add_argument("-cff", default=0, type=float, help="The C_ff coupling. Default value is C_ff = 0")

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
#production mode:
if args.prod == "":
    print(" for all production modes available")
    channels_production = setup.channels_production
elif args.prod in setup.channels_production:
    print(" for ", args.prod, " production mode")
    channels_production = [args.prod]
else:
    parser.error("Production mode " + args.prod + " not available. Production modes available: prod = primakoff | photonfrommeson | mixingPi0 | mixingEta | mixingEtaPrim | BmesonK | BmesonKstar | DmesonPi. If not specified, running over all production modes available.")
#decay mode:
if args.decay == "":
    print(" for all decay modes available.")
    channels_decay = setup.channels_decay
elif args.decay in setup.channels_decay:
    print(" for ", args.decay, " decay mode.")
    channels_decay = [args.decay]
else:
    parser.error("Decay mode " + args.decay + " not available. Decay modes available: decay = 2Gamma | 2El | 2Mu | 3Pi0 | 3Pi | 2PiGamma | 2Pi0Eta | 2PiEta. If not specified, running over all decay modes available.")
#scale
if args.lam <= 0:
    parser.error("\u039B has to be a positive number.")
#couplings
if args.cbb == args.cww == args.cgg == args.cll == 0 :
    parser.error("At least one coupling has to be != 0")

this_dir = path.dirname(path.realpath(__file__))

def ALP_events_exp(expName, C_GG, C_WW, C_BB, Lambda, AA, BB, C_ll):

    #make logarithmic lists of masses (E-4 to 3 GeV) and couplings (E-11 to E-1) -> note: small stepsize right now
#    C_a_list_log = [ c_a_log for c_a_log in np.arange(-11+np.log10(Lambda),-0.999+np.log10(Lambda),0.05)] # 201 values
#    m_a_list_log = [ m_a_log for m_a_log in np.arange(-4,np.log10(3.05),0.005)] # 896 values (only for log10(3))
    C_a_list_log = [ c_a_log for c_a_log in np.arange(-11+np.log10(Lambda),-0.999+np.log10(Lambda),0.1)] # 101 values
    m_a_list_log = [ m_a_log for m_a_log in np.arange(-4,np.log10(3.05),0.02)] # 224(225) values (only for log10(3))

    process = ld.Process_data(expName,channels_decay,channels_production)
    data_list = [[ [10**m_a, 10**C_a/Lambda, process.ALP_events_EFT(10**m_a, 10**C_a *C_GG, 10**C_a *C_WW, 10**C_a *C_BB, Lambda, AA, BB, 10**C_a *C_ll)] for C_a in C_a_list_log] for m_a in m_a_list_log]
    data = np.reshape(data_list,(len(m_a_list_log)*len(C_a_list_log),3))

    # export
    output_dir = this_dir+'/../tab_toPlot/'
    outPath = output_dir + expName + '/'

    if args.prod == args.decay == "": modes = ""
    elif args.prod != "": modes = "_" + args.prod
    elif args.decay != "": modes = "_" + args.decay
    else: modes = "_" + args.prod + "_" + args.decay

    if C_ll == 0: outfileName = expName + '_CBB' + str(int(C_BB)) + '_CWW' + str(int(C_WW)) + '_CGG' + str(int(C_GG)) + modes + '.dat'
    else: outfileName = expName + '_CBB' + str(int(C_BB)) + '_CWW' + str(int(C_WW)) + '_CGG' + str(int(C_GG)) + '_Cll' + str(int(C_ll)) + modes + '.dat'
    np.savetxt(outPath + outfileName,data)
    print('file ' + outfileName + ' saved to ' + outPath)

    return

for exp in experiments:
    ALP_events_exp(exp,args.cgg,args.cww,args.cbb,args.lam, args.a, args.b, args.cll)
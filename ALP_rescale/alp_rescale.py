from os import path
import numpy as np
import alp_setup as setup
import load_data as ld
from alp_mergeSigRegions import MergeInput
import argparse

parser = argparse.ArgumentParser(description='ALP MC rescaling. \n Select the experiment, production and decay modes and other parameters or use default')
parser.add_argument("-varX", type=str, help="X-axis variable. Options are: ma | CBB | CWW | CGG | Cll | Cqq. Standard choice is ma.")
parser.add_argument("-varY", type=str, help="Y-axis variable. Options are: ma | CBB | CWW | CGG | Cll | Cqq. Standard choice is one of the couplings. For setting other options to zero (couplings only), use -only option as for example: -varY CBB-only")
parser.add_argument("-e","--exp",  default="", type=str, nargs='*', help="Experiments available (case sensitive): exp = NA62 | CHARM | nuCAL | SHiP | DarkQuest | DUNE | SHADOWS | KOTO. If not specified, running over all experiments available.")
parser.add_argument("-p","--prod",  default="", type=str, nargs='*', help="Production modes available (case sensitive): prod = primakoff | photonfrommeson | mixingPi0 | mixingEta | mixingEtaPrim | BmesonK | BmesonKstar | DmesonPi | KSmesonPi0. If not specified, running over all production modes available.")
parser.add_argument("-d","--decay",  default="", type=str, nargs='*', help="Decay modes available (case sensitive): decay = 2Gamma | 2El | 2Mu | 3Pi0 | 3Pi | 2PiGamma | 2Pi0Eta | 2PiEta. If not specified, running over all decay modes available.")
parser.add_argument("-l","--lambda", dest="lam", default=1000, type=float, help="The \u039B [GeV] energy scale. Default value is \u039B = 1000 GeV")
parser.add_argument("-a", default=0, type=float, help="The model-dependent A parameter in [2102.04474]. Default value is A = 0")
parser.add_argument("-b", default=0, type=float, help="The model-dependent B parameter in [2102.04474]. Default value is B = 0")
parser.add_argument("-reg", default=1, type=int, help="Number of signal regions. Default value is reg = 1. ")

args = parser.parse_args()
#X- and Y-axis:
oneCouplingOnly = 0
if '-only' in args.varX:
    parser.error("[Error:] \t Cxx-only option available only for Y-axis.")
if '-only' in args.varY:
    oneCouplingOnly = 1
    args.varY = args.varY.replace('-only', '')
    if args.varY == "ma":
        parser.error("[Error:] \t -only option available only for couplings.")

if args.varX in setup.variables:
    if args.varY in setup.variables:
        if args.varX != args.varY:
            print("[Info:] \t Running with X-axis varible " + args.varX + " and Y-axis variable " + args.varY)
            if oneCouplingOnly:
                print("[Info:] \t Running for " + args.varY + " coupling only, other couplings set to 0")
        else:
            parser.error("[Error:] \t X-axis and Y-axis variable must be different.")
    else:
        parser.error("[Error:] \t Y-axis variable " + args.varY + " not available. Available options: ma | CBB | CWW | CGG | Cll | Cqq.")
else:
    parser.error("[Error:] \t X-axis variable " + args.varX + " not available. Available options: ma | CBB | CWW | CGG | Cll | Cqq.")
if args.varX != "ma" and args.varY != "ma":
    alpMassFixed = True
    alpMass = float(input(" - Enter fixed ALP mass (in GeV): "))
    if alpMass < 0:
        raise ValueError("[Error:] \t Please enter only positive numbers")
else:
    alpMassFixed = False

#couplings:
fixedValues = {}
scaledCouplingsX = {}
scaledCouplingsY = {}
for var in setup.variables:
    if var != args.varX and var != args.varY and var != "ma":
        if oneCouplingOnly:
            fixedValues[var] = 0
        else:
            response_fixed = input(" - Fixed value for coupling " + var + "? [Y/N]: ")
            if response_fixed == "Y" or response_fixed == "y":
                fixedValues[var] = float(input(" - Enter fixed value for " + var + " coupling: "))
            elif response_fixed == "N" or response_fixed == "n":
                response_scaleX = input(" - Scale coupling " + var + " with x-axis " + args.varX + "? [Y/N]: ")
                if response_scaleX == "Y" or response_scaleX == "y":
                    scaledCouplingsX[var] = float(input(" - Enter ratio between " + var + " coupling and x-axis " + args.varX + ": "))
                elif response_scaleX == "N" or response_scaleX == "n":
                    print("[Info:] \t Scaling with y-axis " + args.varY)
                    scaledCouplingsY[var] = float(input(" - Enter ratio between " + var + " coupling and y-axis " + args.varY + ": "))
                else:
                    raise KeyError("[Error:] \t Invalid answer: " + response_scaleX + ". Type [Y/N] or [y/n]")
            else:
                raise KeyError("[Error:] \t Invalid answer: " + response_fixed + ". Type [Y/N] or [y/n]")

if alpMassFixed:
    print("[Info:] \t Mass ma fixed value: " + str(alpMass))
    fixedValues["ma"] = float(alpMass)
for fixed in fixedValues.keys():
    print("[Info:] \t Variable " + fixed + " fixed value: " + str(fixedValues[fixed]))
for couplings in scaledCouplingsX.keys():
    print("[Info:] \t Ratio between " + couplings + " and x-axis " + args.varX + ": " + str(scaledCouplingsX[couplings]))
for couplings in scaledCouplingsY.keys():
    print("[Info:] \t Ratio between " + couplings + " and y-axis " + args.varY + ": " + str(scaledCouplingsY[couplings]))

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
channels_production = []
if args.prod == "":
    channels_production = setup.channels_production
    print("[Info:] \t Selected all production modes available")
else:
    for prod in args.prod:
        if prod in setup.channels_production:
            channels_production.append(prod)
        else:
            parser.error("[Error:] \t Production mode " + prod + " not available. Production modes available: prod = " + ' | '.join(setup.channels_production) + ". If not specified, running over all production modes available.")
    print("[Info:] \t Selected production modes:", ', '.join(channels_production))

#decay mode:
channels_decay = []
if args.decay == "":
    channels_decay = setup.channels_decay
    print("[Info:] \t Selected all decay modes available.")
else:
    for dec in args.decay:
        if dec in setup.channels_decay:
            channels_decay.append(dec)
        else:
            parser.error("[Error:] \t Decay mode " + dec + " not available. Decay modes available: decay = " + ' | '.join(setup.channels_decay) + ". If not specified, running over all decay modes available.")
    print("[Info:] \t Selected decay modes:", ', '.join(channels_decay))
#scale
if args.lam <= 0:
    parser.error("[Error:] \t \u039B has to be a positive number.")

#signal regions:
if args.reg > 0:
    if args.reg == 1:
        regions = [""]
    else:
        if len(experiments) == 1:
            regions = ["_reg" + str(sigReg+1) for sigReg in range(args.reg)]
        else:
            parser.error("[Error:] \t Option for multiple signal regions is available only when running for one experiment")
else:
    parser.error("[Error:] \t At least one signal region needed")

#directory
this_dir = path.dirname(path.realpath(__file__))

variables_values = {'ma' : 0,
                    'CBB' : 0,
                    'CWW' : 0,
                    'CGG' : 0,
                    'Cll' : 0,
                    'Cqq' : 0,                  
                    }

for fixed in fixedValues.keys():
    variables_values[fixed] = fixedValues[fixed]


def ALP_events_exp(expName, Lambda, AA, BB):
    #x- and y-axis grid
    c_a_bin_width = 0.025 #0.1: 101 values, 0.05: 201 values, 0.025: 401 values
    m_a_bin_width = 0.005 #0.02: 225 values, 0.005: 896 values
    #make logarithmic lists of masses (E-4 to 3 GeV) and couplings (E-11 to E-1) -> note: small stepsize right now
    c_a_list_log = [ c_a_log for c_a_log in np.arange(-11+np.log10(Lambda),-0.999+np.log10(Lambda),c_a_bin_width)]
    m_a_list_log = [ m_a_log for m_a_log in np.arange(-4,np.log10(3.05),m_a_bin_width)]
    if args.varX == "ma":
        x_axis_list_log = m_a_list_log
        y_axis_list_log = c_a_list_log
    elif args.varY == "ma":
        y_axis_list_log = m_a_list_log
        x_axis_list_log = c_a_list_log
    else:
        x_axis_list_log = c_a_list_log
        y_axis_list_log = c_a_list_log

    process = ld.Process_data(expName,channels_decay,channels_production, len(x_axis_list_log)*len(y_axis_list_log))

    #fill the table
    data_list = []
    for x in x_axis_list_log:
        variables_values[args.varX] = 10**x
        for x_var in scaledCouplingsX.keys():
            variables_values[x_var] = scaledCouplingsX[x_var]*10**x
        data_sublist = []
        for y in y_axis_list_log:
            variables_values[args.varY] = 10**y
            for y_var in scaledCouplingsY.keys():
                variables_values[y_var] = scaledCouplingsY[y_var]*10**y
            if args.varX == "ma":
                data_sublist.append([10**x, 10**y/Lambda, process.ALP_events_EFT(variables_values['ma'], variables_values['CGG'], variables_values['CWW'], variables_values['CBB'], Lambda, AA, BB, variables_values['Cll'], variables_values['Cqq'])])
            elif args.varY == "ma":
                data_sublist.append([10**x/Lambda, 10**y, process.ALP_events_EFT(variables_values['ma'], variables_values['CGG'], variables_values['CWW'], variables_values['CBB'], Lambda, AA, BB, variables_values['Cll'], variables_values['Cqq'])])
            else:
                data_sublist.append([10**x, 10**y, process.ALP_events_EFT(variables_values['ma'], variables_values['CGG'], variables_values['CWW'], variables_values['CBB'], Lambda, AA, BB, variables_values['Cll'], variables_values['Cqq'])])
        data_list.append([data_sublist])
    data = np.reshape(data_list,(len(x_axis_list_log)*len(y_axis_list_log),3))

    # export
    output_dir = this_dir+'/../tab_toPlot/'
    outPath = output_dir + setup.experiments[expName] + '/'

    if args.prod == args.decay == "": modes = ""
    elif args.prod != "" and args.decay == "": modes = "_" + '-'.join(channels_production)
    elif args.decay != "" and args.prod == "": modes = "_" + '-'.join(channels_decay)
    else: modes = "_" + '-'.join(channels_production) + "_" + '-'.join(channels_decay)

    fixedNames = ""
    for fixed in fixedValues.keys():
        if fixed == "ma":
            fixedNames = fixedNames + "-" + fixed + str(int(fixedValues[fixed]*1000)) + "MeV"
        else:
            fixedNames = fixedNames + "-" + fixed + str(int(fixedValues[fixed]))
    if fixedNames != "":
        fixedNames = "_fixed" + fixedNames

    xscaledNames = ""
    for couplings in scaledCouplingsX.keys():
        xscaledNames = xscaledNames + "-" + couplings + str(int(scaledCouplingsX[couplings]))
    if xscaledNames != "":
        xscaledNames = "_scaleWith" + args.varX + xscaledNames

    yscaledNames = ""
    for couplings in scaledCouplingsY.keys():
        yscaledNames = yscaledNames + "-" + couplings + str(int(scaledCouplingsY[couplings]))
    if yscaledNames != "":
        yscaledNames = "_scaleWith" + args.varY + yscaledNames

    outfileName = expName + "_" + args.varX + "_" + args.varY + fixedNames + xscaledNames + yscaledNames + modes + '.dat'
    np.savetxt(outPath + outfileName,data)
    print('\n[Info:] \t', 'File ' + outfileName + ' saved to ' + outPath)

    return

for exp in experiments:
    if not regions == [""]:
        mergeSigReg = MergeInput(exp,regions,channels_decay,channels_production)
    ALP_events_exp(exp,args.lam, args.a, args.b)
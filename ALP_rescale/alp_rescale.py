from os import path
import numpy as np
import alp_setup as setup
import load_data as ld
import argparse

parser = argparse.ArgumentParser(description='ALP MC rescaling. \n Select the experiment, production and decay modes and other parameters or use default')
parser.add_argument("-varX", type=str, help="X-axis variable. Options are: ma | CBB | CWW | CGG | Cll | Cqq. Standard choice is ma.")
parser.add_argument("-varY", type=str, help="Y-axis variable. Options are: ma | CBB | CWW | CGG | Cll | Cqq. Standard choice is one of the couplings.")
parser.add_argument("-e","--exp",  default="", type=str, help="Experiments available (case sensitive): exp = NA62 | CHARM | nuCAL | SHiP | DarkQuest | DUNE | SHADOWS | KOTO. If not specified, running over all experiments available.")
parser.add_argument("-p","--prod",  default="", type=str, help="Production modes available (case sensitive): prod = primakoff | photonfrommeson | mixingPi0 | mixingEta | mixingEtaPrim | BmesonK | BmesonKstar | DmesonPi | KSmesonPi0. If not specified, running over all production modes available.")
parser.add_argument("-d","--decay",  default="", type=str, help="Decay modes available (case sensitive): decay = 2Gamma | 2El | 2Mu | 3Pi0 | 3Pi | 2PiGamma | 2Pi0Eta | 2PiEta. If not specified, running over all decay modes available.")
parser.add_argument("-l","--lambda", dest="lam", default=1000, type=float, help="The \u039B [GeV] energy scale. Default value is \u039B = 1000 GeV")
parser.add_argument("-a", default=0, type=float, help="The model-dependent A parameter in [2102.04474]. Default value is A = 0")
parser.add_argument("-b", default=0, type=float, help="The model-dependent B parameter in [2102.04474]. Default value is B = 0")
parser.add_argument("-reg", default=1, type=int, help="Number of signal regions. Default value is reg = 1. ")

args = parser.parse_args()
#X- and Y-axis:
if args.varX in setup.variables:
    if args.varY in setup.variables:
        if args.varX != args.varY:
            print("[Info:] \t Running with X-axis varible " + args.varX + " and Y-axis variable " + args.varY)
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
if args.exp == "":
    print("[Info:] \t Running for all experiments available")
    experiments = setup.experiments
elif args.exp in setup.experiments:
    print("[Info:] \t Running for ", args.exp, " experiment")
    experiments = [args.exp]
else:
    parser.error("[Error:] \t Experiment " + args.exp + " not available. Experiments available: exp = NA62 | CHARM | nuCAL | SHiP | DarkQuest | DUNE | SHADOWS. If not specified, running over all experiments available.")
#production mode:
if args.prod == "":
    print("[Info:] \t Running for all production modes available")
    channels_production = setup.channels_production
elif args.prod in setup.channels_production:
    print("[Info:] \t Running for ", args.prod, " production mode")
    channels_production = [args.prod]
else:
    parser.error("[Error:] \t Production mode " + args.prod + " not available. Production modes available: prod = primakoff | photonfrommeson | mixingPi0 | mixingEta | mixingEtaPrim | BmesonK | BmesonKstar | DmesonPi. If not specified, running over all production modes available.")
#decay mode:
if args.decay == "":
    print("[Info:] \t Running for all decay modes available.")
    channels_decay = setup.channels_decay
elif args.decay in setup.channels_decay:
    print("[Info:] \t Running for ", args.decay, " decay mode.")
    channels_decay = [args.decay]
else:
    parser.error("[Error:] \t Decay mode " + args.decay + " not available. Decay modes available: decay = 2Gamma | 2El | 2Mu | 3Pi0 | 3Pi | 2PiGamma | 2Pi0Eta | 2PiEta. If not specified, running over all decay modes available.")
#scale
if args.lam <= 0:
    parser.error("[Error:] \t \u039B has to be a positive number.")

#signal regions:
if args.reg > 0:
    if args.reg == 1:
        regions = [""]
    else:
        regions = ["_reg" + str(sigReg+1) for sigReg in range(args.reg)]
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


def ALP_events_exp(expName, sigRegion, Lambda, AA, BB):
    #x- and y-axis grid
    #make logarithmic lists of masses (E-4 to 3 GeV) and couplings (E-11 to E-1) -> note: small stepsize right now
    c_a_list_log = [ c_a_log for c_a_log in np.arange(-11+np.log10(Lambda),-0.999+np.log10(Lambda),0.1)] # 101 values
    m_a_list_log = [ m_a_log for m_a_log in np.arange(-4,np.log10(3.05),0.02)] # 224(225) values (only for log10(3))
#    c_a_list_log = [ c_a_log for c_a_log in np.arange(-11+np.log10(Lambda),-0.999+np.log10(Lambda),0.05)] # 201 values
#    m_a_list_log = [ m_a_log for m_a_log in np.arange(-4,np.log10(3.05),0.005)] # 896 values (only for log10(3))
    if args.varX == "ma":
        x_axis_list_log = m_a_list_log
        y_axis_list_log = c_a_list_log
    elif args.varY == "ma":
        y_axis_list_log = m_a_list_log
        x_axis_list_log = c_a_list_log
    else:
        x_axis_list_log = c_a_list_log
        y_axis_list_log = c_a_list_log

    process = ld.Process_data(expName,sigRegion,channels_decay,channels_production)

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
    outPath = output_dir + expName + '/'

    if args.prod[0] == args.decay[0] == "": modes = ""
    elif args.prod[0] == "": modes = "_" + "-".join(args.decay)
    elif args.decay[0] == "": modes = "_" + "-".join(args.prod)
    else: modes = "_" + "-".join(args.prod) + "_" + "-".join(args.decay)

    fixedNames = ""
    for fixed in fixedValues.keys():
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

    outfileName = expName + "_" + args.varX + "_" + args.varY + fixedNames + xscaledNames + yscaledNames + modes + sigRegion + '.dat'
    np.savetxt(outPath + outfileName,data)
    print('file ' + outfileName + ' saved to ' + outPath)

    return

for exp in experiments:
    for sigReg in regions:
        ALP_events_exp(exp,sigReg,args.lam, args.a, args.b)

#!/usr/bin/env python3

from os import path
import numpy as np
from alp import alp_setup as setup
from scalar import scalar_setup as setds
from alp import load_data as ld
from general import alp_functions as f
from general.alp_mergeSigRegions import MergeInput
import argparse
import sys


def main(argv=None):
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    parser = argparse.ArgumentParser(description='ALP MC rescaling. \n Select the experiment, production and decay modes and other parameters or use default')
    parser.add_argument("-varX", type=str, help="X-axis variable. Options are: mX | tauX | GammaX | BRdecay | BRprod | coupling (for ALP: CBB | CWW | CGG | Cll | Cqq ; for DS: Y). Standard choice is mX.")
    parser.add_argument("-varY", type=str, help="Y-axis variable. Options are: mX | tauX | GammaX | BRdecay | BRprod | coupling (for ALP: CBB | CWW | CGG | Cll | Cqq ; for DS: Y). Standard choice is one of the couplings. For setting other couplings to zero, use -only option as for example: -varY CBB-only")
    parser.add_argument("-e","--exp",  default="", type=str, nargs='*', help="Experiments available (case sensitive): exp = NA62 | CHARM | NuCal | SHiP | DarkQuest | DUNE | SHADOWS | KOTO. If not specified, running over all experiments available.")
    parser.add_argument("-p","--prod",  default="", type=str, nargs='*', help="Production modes available (case sensitive): prod = primakoff | photonfrommeson | mixingPi0 | mixingEta | mixingEtaPrim | BmesonK | BmesonKstar | DmesonPi | KSmesonPi0. If not specified, running over all production modes available.")
    parser.add_argument("-d","--decay",  default="", type=str, nargs='*', help="Decay modes available (case sensitive): decay = 2Gamma | 2El | 2Mu | 3Pi0 | 3Pi | 2PiGamma | 2Pi0Eta | 2PiEta (for DS 2Pi and 2K). If not specified, running over all decay modes available.")
    parser.add_argument("-l","--lambda", dest="lam", default=1000, type=float, help="The \u039B [GeV] energy scale. Default value is \u039B = 1000 GeV")
    parser.add_argument("-a", default=0, type=float, help="The model-dependent A parameter in [2102.04474]. Default value is A = 0")
    parser.add_argument("-b", default=0, type=float, help="The model-dependent B parameter in [2102.04474]. Default value is B = 0")
    parser.add_argument("-reg", default=1, type=int, help="Number of signal regions. Default value is reg = 1. ")
    parser.add_argument("-exo","--exotic",  default="ALP", type=str, help="Exotic particles available: exo = ALP | DS . Default is ALP")
    parser.add_argument("-comb","--combine", dest='comb', action='store_true', help="Combine experiments")
    parser.add_argument("-no-comb","--no-combine", dest='comb', action='store_false', help="Rescale experiments separately. Default")
    parser.set_defaults(comb=False)

    args = parser.parse_args()
    #particle
    if args.exotic == "DS":
        vars = setds.variables
        coups = setds.couplings
        exps = setds.experiments
        prods = setds.channels_production
        decs = setds.channels_decay
    else:
        vars = setup.variables
        coups = setup.couplings
        exps = setup.experiments
        prods = setup.channels_production
        decs = setup.channels_decay

    #X- and Y-axis:
    oneCouplingOnly = 0
    if '-only' in args.varX:
        parser.error("[Error:] \t Cxx-only option available only for Y-axis.")
    if '-only' in args.varY:
        oneCouplingOnly = 1
        args.varY = args.varY.replace('-only', '')
        if args.varY not in coups:
            parser.error("[Error:] \t -only option available only for ALP couplings.")
    varAxes = [] #X and Y axis variable names

    if args.varX in vars or args.varX in coups:
        if args.varY in vars or args.varY in coups:
            if args.varX != args.varY:
                print("[Info:] \t Running with X-axis varible " + args.varX + " and Y-axis variable " + args.varY)
                varAxes.append(args.varX)
                varAxes.append(args.varY)
                if oneCouplingOnly:
                    print("[Info:] \t Running for " + args.varY + " coupling only, other couplings set to 0")
            else:
                parser.error("[Error:] \t X-axis and Y-axis variable must be different.")
        else:
            parser.error("[Error:] \t Y-axis variable " + args.varY + " not available. Available options: mX | tauX | GammaX | BRdecay | BRprod | coupling (for ALP: CBB | CWW | CGG | Cll | Cqq ; for DS: Y).")
    else:
        parser.error("[Error:] \t X-axis variable " + args.varX + " not available. Available options: mX | tauX | GammaX | BRdecay | BRprod | coupling (for ALP: CBB | CWW | CGG | Cll | Cqq ; for DS: Y).")
    if "mX" not in varAxes:
        massFixed = True
        exoMass = float(input(" - Enter fixed exotic particle mass (in GeV): "))
        if exoMass < 0:
            raise ValueError("[Error:] \t Please enter only positive numbers")
    else:
        massFixed = False

    #couplings:
    fixedValues = {}
    scaledCouplingsX = {}
    scaledCouplingsY = {}
    if "C" in '\t'.join(varAxes) or "Y" in '\t'.join(varAxes):
        for var in coups:
            if var not in varAxes:
            # if var != args.varX and var != args.varY:
            # if var != args.varX and var != args.varY and var != "mX":
                # if "BR" not in var and var != "GammaX" and var != "tauX":
                if oneCouplingOnly:
                    fixedValues[var] = 0
                else:
                    response_fixed = input(" - Fixed value for " + var + "? [Y/N]: ")
                    if response_fixed == "Y" or response_fixed == "y":
                        fixedValues[var] = float(input(" - Enter fixed value for " + var + ": "))
                        # if ("BR" in var) and (response_fixed > 1 or response_fixed < 0):
                        #     raise ValueError("[Error:] \t BR has to be from interval [0,1]")
                    elif response_fixed == "N" or response_fixed == "n":
                        response_scaleX = input(" - Scale " + var + " with x-axis " + args.varX + "? [Y/N]: ")
                        if response_scaleX == "Y" or response_scaleX == "y":
                            scaledCouplingsX[var] = float(input(" - Enter ratio between " + var + " and x-axis " + args.varX + ": "))
                        elif response_scaleX == "N" or response_scaleX == "n":
                            print("[Info:] \t Scaling with y-axis " + args.varY)
                            scaledCouplingsY[var] = float(input(" - Enter ratio between " + var + " and y-axis " + args.varY + ": "))
                        else:
                            raise KeyError("[Error:] \t Invalid answer: " + response_scaleX + ". Type [Y/N] or [y/n]")
                    else:
                        raise KeyError("[Error:] \t Invalid answer: " + response_fixed + ". Type [Y/N] or [y/n]")
    else:
        for var in vars:
            if var in varAxes: continue # ask only variables that are not 
            if var == "tauX" or var == "GammaX": continue
            if var == "mX": continue # fixed mass asked separately
            response_fixed = input(" - Fixed value for " + var + "? [Y/N]: ")
            if response_fixed == "Y" or response_fixed == "y":
                fixedValues[var] = float(input(" - Enter fixed value for " + var + ": "))
                if ("BR" in var) and (fixedValues[var] > 1 or fixedValues[var] < 0):
                    raise ValueError("[Error:] \t BR has to be from interval [0,1]")
                if var == "GammaX":
                    fixedValues["tauX"] = f.tau(fixedValues["GammaX"])
            elif response_fixed == "N" or response_fixed == "n":
                response_scaleX = input(" - Scale " + var + " with x-axis " + args.varX + "? [Y/N]: ")
                if response_scaleX == "Y" or response_scaleX == "y":
                    scaledCouplingsX[var] = float(input(" - Enter ratio between " + var + " and x-axis " + args.varX + ": "))
                elif response_scaleX == "N" or response_scaleX == "n":
                    print("[Info:] \t Scaling with y-axis " + args.varY)
                    scaledCouplingsY[var] = float(input(" - Enter ratio between " + var + " and y-axis " + args.varY + ": "))
                else:
                    raise KeyError("[Error:] \t Invalid answer: " + response_scaleX + ". Type [Y/N] or [y/n]")
            else:
                raise KeyError("[Error:] \t Invalid answer: " + response_fixed + ". Type [Y/N] or [y/n]")
            # else:
    # elif args.varX == "GammaX" or args.varX == "tauX" or args.varY == "GammaX" or args.varY == "tauX":

    if massFixed:
        print("[Info:] \t Mass mX fixed value: " + str(exoMass))
        fixedValues["mX"] = float(exoMass)
    for fixed in fixedValues.keys():
        print("[Info:] \t Variable " + fixed + " fixed value: " + str(fixedValues[fixed]))
    for couplings in scaledCouplingsX.keys():
        print("[Info:] \t Ratio between " + couplings + " and x-axis " + args.varX + ": " + str(scaledCouplingsX[couplings]))
    for couplings in scaledCouplingsY.keys():
        print("[Info:] \t Ratio between " + couplings + " and y-axis " + args.varY + ": " + str(scaledCouplingsY[couplings]))

    #experiment:
    experiments = []
    if args.exp == "":
        experiments = exps.keys()
        print("[Info:] \t Selected all experiments available")
    else:
        for exp in args.exp:
            if exp in exps.keys():
                experiments.append(exp)
            else:
                parser.error("[Error:] \t Experiment " + exp + " not available. Experiment modes available: exp = " + ' | '.join(exps.keys()) + ". If not specified, running over all experiments available.")
        print("[Info:] \t Selected experiments:", ', '.join(experiments))
    if args.comb:
        print("[Info:] \t Datasets for selected experiments will be combined")
    #production mode:
    channels_production = []
    if args.prod == "":
        channels_production = prods
        print("[Info:] \t Selected all production modes available")
    else:
        for prod in args.prod:
            if prod in prods:
                channels_production.append(prod)
                if "BRprod" in varAxes:
                    if "Bmeson" not in prod and "Dmeson" not in prod:
                        parser.error("[Error:] \t Scanning over BRprod possible only for production in meson decay")
            else:
                parser.error("[Error:] \t Production mode " + prod + " not available. Production modes available: prod = " + ' | '.join(prods) + ". If not specified, running over all production modes available.")
        print("[Info:] \t Selected production modes:", ', '.join(channels_production))

    #decay mode:
    channels_decay = []
    if args.decay == "":
        channels_decay = decs
        print("[Info:] \t Selected all decay modes available.")
    else:
        for dec in args.decay:
            if dec in decs:
                channels_decay.append(dec)
            else:
                parser.error("[Error:] \t Decay mode " + dec + " not available. Decay modes available: decay = " + ' | '.join(decs) + ". If not specified, running over all decay modes available.")
        print("[Info:] \t Selected decay modes:", ', '.join(channels_decay))
    if "BRdecay" in varAxes:
        if len(channels_decay) != 1:
            parser.error("[Error:] \t Only one decay mode has to be selected when scanning over BRdecay")
        else:
            print("[Info:] \t Will scan over BRdecay for mode " + channels_decay[0])

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

    variables_values = {}
    for var in vars + coups: #initialize with zeros
        variables_values[var] = variables_values.get(var,0)

    for fixed in fixedValues.keys():
        variables_values[fixed] = fixedValues[fixed]

    #generate x- and y- tables
    x_axis_list_log = generate_log_list(args.varX, args.lam)
    y_axis_list_log = generate_log_list(args.varY, args.lam)

    #rescale
    comb_data = []

    for exp in experiments:
        if not regions == [""]:
            mergeSigReg = MergeInput(exp,regions,channels_decay,channels_production)

        process = ld.Process_data(exp,channels_decay,channels_production, len(x_axis_list_log)*len(y_axis_list_log),args.exotic)

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
                if args.exotic == "ALP" and "C" in '\t'.join(varAxes):
                    if args.varX == "mX":
                        data_sublist.append([10**x, 10**y/args.lam, process.ALP_events_EFT(variables_values['mX'], variables_values['CGG'], variables_values['CWW'], variables_values['CBB'], args.lam, args.a, args.b, variables_values['Cll'], variables_values['Cqq'])])
                    elif args.varY == "mX":
                        data_sublist.append([10**x/args.lam, 10**y, process.ALP_events_EFT(variables_values['mX'], variables_values['CGG'], variables_values['CWW'], variables_values['CBB'], args.lam, args.a, args.b, variables_values['Cll'], variables_values['Cqq'])])
                    else:
                        data_sublist.append([10**x, 10**y, process.ALP_events_EFT(variables_values['mX'], variables_values['CGG'], variables_values['CWW'], variables_values['CBB'], args.lam, args.a, args.b, variables_values['Cll'], variables_values['Cqq'])])
                elif args.exotic == "DS" and "Y" in '\t'.join(varAxes):
                    data_sublist.append([10**x, 10**y, process.DS_events_EFT(variables_values['mX'], np.sqrt(variables_values['Y']))])
                
                elif "BR" in '\t'.join(varAxes):
                    # process.Set_Coupling_Production("BRprod",variables_values['BRprod'])
                    if "GammaX" in varAxes:
                        data_sublist.append([f'{10**x:.4f}', 10**y, process.n_events_Bmeson_BR(channels_production,channels_decay,variables_values['BRprod'],variables_values['BRdecay'], variables_values['GammaX'], variables_values['mX'])])
                    else:
                        data_sublist.append([10**x, 10**y, process.n_events_Bmeson_BR(channels_production,channels_decay,variables_values['BRprod'],variables_values['BRdecay'], f.Gamma(variables_values['tauX']), variables_values['mX'])])
                elif "mX" in varAxes:
                    # process.Set_Coupling_Production('BRprod',variables_values['BRprod'])
                    if "GammaX" in varAxes:
                        data_sublist.append([10**x, 10**y, process.n_events_Bmeson_BR(channels_production,channels_decay,variables_values['BRprod'],variables_values['BRdecay'], variables_values['GammaX'], variables_values['mX'])])
                    else:
                        data_sublist.append([10**x, 10**y, process.n_events_Bmeson_BR(channels_production,channels_decay,variables_values['BRprod'],variables_values['BRdecay'], f.Gamma(variables_values['tauX']), variables_values['mX'])])
                else: parser.error("[Error:] \t Case not resolved")
            data_list.append([data_sublist])
        data = np.reshape(data_list,(len(x_axis_list_log)*len(y_axis_list_log),3))
        if not args.comb: #store output for experiments separately

            # export
            output_dir = this_dir+'/../tab_toPlot/'
            outPath = output_dir + exps[exp] + '/'

            if args.prod == args.decay == "": modes = ""
            elif args.prod != "" and args.decay == "": modes = "_" + '-'.join(channels_production)
            elif args.decay != "" and args.prod == "": modes = "_" + '-'.join(channels_decay)
            else: modes = "_" + '-'.join(channels_production) + "_" + '-'.join(channels_decay)

            fixedNames = ""
            for fixed in fixedValues.keys():
                if fixed == "mX":
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

            outfileName = exp + "_" + args.varX + "_" + args.varY + fixedNames + xscaledNames + yscaledNames + modes + '.dat'
            np.savetxt(outPath + outfileName,data,fmt='%.4e')
            print('\n[Info:] \t', 'File ' + outfileName + ' saved to ' + outPath)
        
        else:
            if len(comb_data) == 0:
                comb_data = data
            else:
                if len(data) == len(comb_data) and (data[:, 0]==comb_data[:, 0]).all() and (data[:, 1]==comb_data[:, 1]).all():
                    for line in range(len(data)):
                        comb_data[line,2] += data[line,2]
                else:
                    print("[Error:] \t X and Y axes differ for the experiments. Exiting.")
                    sys.exit(1)
            
    if args.comb: #export output for combined experiments
        outPath = this_dir+'/../tab_toPlot/combined/'

        if args.prod == args.decay == "": modes = ""
        elif args.prod != "" and args.decay == "": modes = "_" + '-'.join(channels_production)
        elif args.decay != "" and args.prod == "": modes = "_" + '-'.join(channels_decay)
        else: modes = "_" + '-'.join(channels_production) + "_" + '-'.join(channels_decay)

        fixedNames = ""
        for fixed in fixedValues.keys():
            if fixed == "mX":
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

        outfileName = "-".join(experiments) + "_" + args.varX + "_" + args.varY + fixedNames + xscaledNames + yscaledNames + modes + '.dat'
        np.savetxt(outPath + outfileName,comb_data,fmt='%.4e')
        print('\n[Info:] \t', 'File ' + outfileName + ' saved to ' + outPath)
    
        


def generate_log_list(var,Lambda):
    n_bins = 601
    list = []
    if var == "mX": #GeV
        # list = [ m_a_log for m_a_log in np.linspace(np.log10(0.0001),np.log10(3.01),n_bins)]
        list = [ m_a_log for m_a_log in np.linspace(np.log10(0.01),np.log10(3.01),n_bins)]
    elif var == "tauX": #ps
        tau_bin_width = 0.08
        list = [ tau_log for tau_log in np.arange(-2,6,tau_bin_width)]
    elif var == "GammaX": #GeV
        Gamma_bin_width = 0.15
        list = [ Gamma_log for Gamma_log in np.arange(-25,-10,Gamma_bin_width)]
    elif "C" in var: #GeV^-1
        n_bins = 201
        list = [ c_a_log for c_a_log in np.linspace(-11+np.log10(Lambda),-0.999+np.log10(Lambda),n_bins)]
    elif "BR" in var:
        n_bins = 101
        list = [ BR_log for BR_log in np.linspace(-13,-1,n_bins)]
    elif "Y" in var:
        list = [ Y2_log for Y2_log in np.linspace(-13,-1,n_bins)]
    else:
        print("[Error:] \t Variable " + var + " not recognized")
        sys.exit(1)

    return list

def generate_lin_list(min,max,nbins):
    return [ m_a_log for m_a_log in np.linspace(min,max,nbins)]

if __name__ == "__main__":
    sys.exit(main())
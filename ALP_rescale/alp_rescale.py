#!/usr/bin/env python3

 ##
 # @file alp_rescale.py
 # @brief ALP_rescale main module
 # @details Contains main function and parses arguments
 # 
 # Combines datasets for given production and decay channels for selected model-dependent parameters
 # Allows combining sensitivities of various experiments
 # To get available mandatory and optional arguments run:
 # @code
 # python -m ALP_rescale.alp_rescale -h
 # @endcode
 
from os import path, makedirs
import numpy as np
import argparse
import sys

from ALP_rescale.alp import alp_setup as setalp
from ALP_rescale.scalar import scalar_setup as setds
from ALP_rescale.vector import vector_setup as setdp
from ALP_rescale.hnl import hnl_setup as sethnl
from ALP_rescale.general import load_data as ld
from ALP_rescale.general import functions as f
from ALP_rescale.general.setup import experiments as exp_dict
from ALP_rescale.general.mergeSigRegions import MergeInput


def main(argv=None):
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    parser = argparse.ArgumentParser(description='ALP MC rescaling. \n Select the experiment, production and decay modes and other parameters or use default')
    parser.add_argument("-x","--exotic",  default="alp", type=str, help="Exotic particles available: exotic = alp | hnl | ds | dp . Default is ALP")
    parser.add_argument("-varX", type=str, help="X-axis variable. Options are: "+" | ".join(setalp.variables) +" | coupling (for ALP: "+" | ".join(setalp.couplings)+" ; for HNL: "+" | ".join(sethnl.couplings)+" ; for DS: "+" | ".join(setds.couplings)+" ; for DP: "+" | ".join(setdp.couplings)+"). Standard choice is mX.")
    parser.add_argument("-varY", type=str, help="Y-axis variable. Options are: "+" | ".join(setalp.variables) +" | coupling (for ALP: "+" | ".join(setalp.couplings)+" ; for HNL: "+" | ".join(sethnl.couplings)+" ; for DS: "+" | ".join(setds.couplings)+" ; for DP: "+" | ".join(setdp.couplings)+"). Standard choice is one of the couplings. For setting other couplings to zero, use -only option as for example: -varY CBB-only")
    parser.add_argument("-e","--exp",  default="", type=str, nargs='*', help="Experiments available (case sensitive): exp = "+" | ".join(setalp.experiments)+". If not specified, running over all experiments available.")
    parser.add_argument("-p","--prod",  default="", type=str, nargs='*', help="Production modes available (case sensitive): prod = primakoff | photonfrommeson | mixingPi0 | mixingEta | mixingEtaPrim | BmesonK | BmesonKstar | DmesonPi | KSmesonPi0. If not specified, running over all production modes available.")
    parser.add_argument("-d","--decay",  default="", type=str, nargs='*', help="Decay modes available (case sensitive): decay = 2Gamma | 2El | 2Mu | 3Pi0 | 3Pi | 2PiGamma | 2Pi0Eta | 2PiEta (for DS 2Pi and 2K). If not specified, running over all decay modes available.")
    parser.add_argument("-l","--lambda", dest="lam", default=1000, type=float, help="The \u039B [GeV] energy scale. Default value is \u039B = 1000 GeV")
    parser.add_argument("-a", default=0, type=float, help="The model-dependent A parameter in [2102.04474]. Default value is A = 0")
    parser.add_argument("-b", default=0, type=float, help="The model-dependent B parameter in [2102.04474]. Default value is B = 0")
    parser.add_argument("-reg", default=1, type=int, help="Number of signal regions. Default value is reg = 1. ")
    parser.add_argument("-comb","--combine", dest='comb', action='store_true', help="Combine experiments")
    parser.add_argument("-no-comb","--no-combine", dest='comb', action='store_false', help="Rescale experiments separately. Default")
    parser.set_defaults(comb=False)
    parser.add_argument("-dirac","--Dirac-HNL", dest='dirac', action='store_true', help="Consider the decay width of a Dirac rather than Majorana fermion (only for hnl).")
    parser.add_argument("-ternary","--ternary-evaluation", dest='ter', action='store_true', help="Generate ternarily evaluated bounds at a fixed mass (only for hnl).")
    parser.set_defaults(tern=False)
    parser.add_argument("-standalone","--standalone", dest='standalone', action='store_true', help="Run standalone for given channel, no referential coupling reweight")
    parser.set_defaults(standalone=False)

    args = parser.parse_args()

    rescale(args.exotic, args.varX, args.varY, args.exp, args.prod, args.decay, args.lam, args.a, args.b, args.reg, args.comb, ter=args.ter, dirac=args.dirac, standalone=args.standalone)
    
 ##
 # @fn rescale
 # @brief terminal interface for coupling setup
def rescale(exotic, varX, varY, exp, prod, decay, lam, a, b, reg, comb, ter=False, dirac=False,standalone=False):

    if ter: # if ternary eavluation call auxillary module and exit
        if exotic != "hnl":
            print("[Error:] \t Ternary evaluation only available for exotic class hnl. Exiting.")
            sys.exit(1)
        from ALP_rescale.hnl import ternary_rescale
        ternary_rescale.main(exp, prod, decay, comb)
        sys.exit(0)

    for inp_, descript_ in zip([exotic, exp, varX, varY,],['-x (exotic)','-e (experiment)', '-varX (x-axis variable)', '-varY (y-axis variable)']):
        if not inp_: sys.exit(f"[Error:] \t {descript_} is a required input")

    #particle
    if exotic == "ds":
        vars = setds.variables
        coups = setds.couplings
        exps = setds.experiments
        prods = setds.channels_production
        decs = setds.channels_decay
    elif exotic == "dp":
        vars = setdp.variables
        coups = setdp.couplings
        exps = setdp.experiments
        prods = setdp.channels_production
        decs = setdp.channels_decay
    elif exotic == "hnl":
        vars  = sethnl.variables
        coups = sethnl.couplings
        exps  = sethnl.experiments
        prods = sethnl.channels_production
        decs  = sethnl.channels_decay
    else:
        vars = setalp.variables
        coups = setalp.couplings
        exps = setalp.experiments
        prods = setalp.channels_production
        decs = setalp.channels_decay

    #X- and Y-axis:
    oneCouplingOnly = 0
    if '-only' in varX:
        sys.exit("[Error:] \t Cxx-only option available only for Y-axis.")
    if '-only' in varY:
        oneCouplingOnly = 1
        varY = varY.replace('-only', '')
        if varY not in coups:
            sys.exit("[Error:] \t -only option available only for ALP couplings.")
    varAxes = [] #X and Y axis variable names

    if varX in vars or varX in coups or (exotic == "hnl" and varY == "U2"):
        if varY in vars or varY in coups or (exotic == "hnl" and varY == "U2"):
            if varX != varY:
                print("[Info:] \t Running with X-axis varible " + varX + " and Y-axis variable " + varY)
                varAxes.append(varX)
                varAxes.append(varY)
                if oneCouplingOnly:
                    print("[Info:] \t Running for " + varY + " coupling only, other couplings set to 0")
            else:
                sys.exit("[Error:] \t X-axis and Y-axis variable must be different.")            
        else:
            sys.exit("[Error:] \t Y-axis variable " + varY + " not available. Available options: mX | tauX | GammaX | BRdecay | BRprod | coupling (for ALP: CBB | CWW | CGG | Cll | Cqq ; for HNL: U2 | U2el | U2mu | U2tau ; for DS: Y ; for DP: eps).")
    else:
        sys.exit("[Error:] \t X-axis variable " + varX + " not available. Available options: mX | tauX | GammaX | BRdecay | BRprod | coupling (for ALP: CBB | CWW | CGG | Cll | Cqq ; for HNL: U2 | U2el | U2mu | U2tau ; for DS: Y ; for DP: eps).")
    if "mX" not in varAxes:
        massFixed = True
        exoMass = float(input(" - Enter fixed exotic particle mass (in GeV): "))
        if exoMass < 0:
            raise ValueError("[Error:] \t Please enter only positive numbers")
    else:
        massFixed = False

    is_quark_coupling = False
    if "Cqq" in varAxes:
        is_quark_coupling = True
        print("[Warning:] \t Running with quark coupling, omitting gluon coupling")

    #couplings:
    fixedValues = {}
    scaledCouplingsX = {}
    scaledCouplingsY = {}
    if "C" in '\t'.join(varAxes) or "Y" in '\t'.join(varAxes) or "U" in '\t'.join(varAxes):
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
                        unit_input  =  " in GeV" if 'C' in var else ' per GeV' if var == 'Y' else ''
                        fixedValues[var] = float(input(" - Enter fixed value for " + var + unit_input+": "))
                        # if ("BR" in var) and (response_fixed > 1 or response_fixed < 0):
                        #     raise ValueError("[Error:] \t BR has to be from interval [0,1]")
                    elif response_fixed == "N" or response_fixed == "n":
                        response_scaleX = input(" - Scale " + var + " with x-axis " + varX + "? [Y/N]: ")
                        if response_scaleX == "Y" or response_scaleX == "y":
                            scaledCouplingsX[var] = float(input(" - Enter ratio between " + var + " and x-axis " + varX + ": "))
                        elif response_scaleX == "N" or response_scaleX == "n":
                            print("[Info:] \t Scaling with y-axis " + varY)
                            scaledCouplingsY[var] = float(input(" - Enter ratio between " + var + " and y-axis " + varY + ": "))
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
                response_scaleX = input(" - Scale " + var + " with x-axis " + varX + "? [Y/N]: ")
                if response_scaleX == "Y" or response_scaleX == "y":
                    scaledCouplingsX[var] = float(input(" - Enter ratio between " + var + " and x-axis " + varX + ": "))
                elif response_scaleX == "N" or response_scaleX == "n":
                    print("[Info:] \t Scaling with y-axis " + varY)
                    scaledCouplingsY[var] = float(input(" - Enter ratio between " + var + " and y-axis " + varY + ": "))
                else:
                    raise KeyError("[Error:] \t Invalid answer: " + response_scaleX + ". Type [Y/N] or [y/n]")
            else:
                raise KeyError("[Error:] \t Invalid answer: " + response_fixed + ". Type [Y/N] or [y/n]")
            # else:
    # elif varX == "GammaX" or varX == "tauX" or varY == "GammaX" or varY == "tauX":
    xnorm, ynorm = [0., 0.]
    if exotic == "hnl":# in Case U2 is chosen as axis normalize sum of coupling ratios to 1
        if varX == "U2":
            for coupling in scaledCouplingsX.keys():
                if "U" == coupling[0]: xnorm += scaledCouplingsX[coupling]
            for coupling in scaledCouplingsX.keys():
                if "U" == coupling[0]: scaledCouplingsX[coupling] /= xnorm
        elif varY == "U2":
            for coupling in scaledCouplingsY.keys():
                if "U" == coupling[0]: ynorm += scaledCouplingsY[coupling]
            for coupling in scaledCouplingsY.keys():
                if "U" == coupling[0]: scaledCouplingsY[coupling] /= ynorm

    if massFixed:
        print("[Info:] \t Mass mX fixed value: " + str(exoMass))
        fixedValues["mX"] = float(exoMass)
    for fixed in fixedValues.keys():
        print("[Info:] \t Variable " + fixed + " fixed value: " + str(fixedValues[fixed]))
    for couplings in scaledCouplingsX.keys():
        print("[Info:] \t Ratio between " + couplings + " and x-axis " + varX + ": " + str(scaledCouplingsX[couplings]))
    for couplings in scaledCouplingsY.keys():
        print("[Info:] \t Ratio between " + couplings + " and y-axis " + varY + ": " + str(scaledCouplingsY[couplings]))

    #experiment:
    experiments = []
    if exp == "":
        experiments = exps.keys()
        print("[Info:] \t Selected all experiments available")
    else:
        for exp in exp:
            if exp in exps.keys():
                experiments.append(exp)
            else:
                sys.exit("[Error:] \t Experiment " + exp + " not available. Experiment modes available: exp = " + ' | '.join(exps.keys()) + ". If not specified, running over all experiments available.")
        print("[Info:] \t Selected experiments:", ', '.join(experiments))
    if comb:
        print("[Info:] \t Datasets for selected experiments will be combined")
        
    if standalone:
        if len(prod) != 1 or len(decay) != 1:
            sys.exit("[Error:] \t With standalone option the production and decay channel has to be specified")
        else:
            print("[Info:] \t Running in a standalone mode, no referential coupling or BR reweighting")
    #production mode:
    channels_production = []
    if prod == "":
        channels_production = prods
        print("[Info:] \t Selected all production modes available")
    else:
        for prod in prod:
            if prod in prods:
                channels_production.append(prod)
                if "BRprod" in varAxes:
                    if "Bmeson" not in prod and "Dmeson" not in prod:
                        sys.exit("[Error:] \t Scanning over BRprod possible only for production in meson decay")
            else:
                sys.exit("[Error:] \t Production mode " + prod + " not available. Production modes available: prod = " + ' | '.join(prods) + ". If not specified, running over all production modes available.")
        print("[Info:] \t Selected production modes:", ', '.join(channels_production))

    #decay mode:
    channels_decay = []
    if decay == "":
        channels_decay = decs
        print("[Info:] \t Selected all decay modes available.")
    else:
        for dec in decay:
            if dec in decs:
                channels_decay.append(dec)
            else:
                sys.exit("[Error:] \t Decay mode " + dec + " not available. Decay modes available: decay = " + ' | '.join(decs) + ". If not specified, running over all decay modes available.")
        print("[Info:] \t Selected decay modes:", ', '.join(channels_decay))
    if "BRdecay" in varAxes:
        if len(channels_decay) != 1:
            sys.exit("[Error:] \t Only one decay mode has to be selected when scanning over BRdecay")
        else:
            print("[Info:] \t Will scan over BRdecay for mode " + channels_decay[0])

    #scale
    if lam <= 0:
        sys.exit("[Error:] \t \u039B has to be a positive number.")

    #signal regions:
    if reg > 0:
        if reg == 1:
            regions = [""]
        else:
            if len(experiments) == 1:
                regions = ["_reg" + str(sigReg+1) for sigReg in range(reg)]
            else:
                sys.exit("[Error:] \t Option for multiple signal regions is available only when running for one experiment")
    else:
        sys.exit("[Error:] \t At least one signal region needed")


    variables_values = {}
    for var in vars + coups: #initialize with zeros
        variables_values[var] = variables_values.get(var,0)

    for fixed in fixedValues.keys():
        variables_values[fixed] = fixedValues[fixed]    

    rescale_xy(exotic, varAxes, experiments, channels_production, channels_decay, fixedValues, variables_values, scaledCouplingsX, scaledCouplingsY, regions, [lam, a, b], combine_experiments=comb, exp_dict = exps, dirac_width=dirac, standalone=standalone)

    # #directory
    # this_dir = path.dirname(path.realpath(__file__))

    # #generate x- and y- tables
    # x_axis_list_log = generate_log_list(args.varX, args.lam)
    # y_axis_list_log = generate_log_list(args.varY, args.lam)

    #rescale
    
 ##
 # @fn rescale_xy
 # @brief rescale function for custom x- and y-axes
def rescale_xy(exotic, varAxes , experiments, channels_production, channels_decay, fixedValues, variables_values, scaleWithX, scaleWithY,  regions = 1, alp_model_params = [1000.,0.,0.], combine_experiments = False, exp_dict = exp_dict, dirac_width = False, standalone = False):

    varX, varY = varAxes
    #directory
    this_dir = path.dirname(path.realpath(__file__))

    #generate x- and y- tables
    x_axis_list_log = generate_log_list(varX, alp_model_params[0])
    y_axis_list_log = generate_log_list(varY, alp_model_params[0])

    # adjusting for HNL labling scheme
    channels_production_extended = []
    if exotic == "hnl":
        for channel in channels_production:
            for mixing in ["El","Mu","Tau"]:
                channels_production_extended.append(channel+ "-"+mixing+"Mixing")

    comb_data = []

    for exp in experiments:
        if not regions == [""]:
            mergeSigReg = MergeInput(exp,regions,channels_decay,channels_production)

        process = ld.Process_data(exp,channels_decay,channels_production, len(x_axis_list_log)*len(y_axis_list_log),exotic, run_standalone = standalone) if exotic != "hnl" else ld.Process_data(exp,channels_decay,channels_production_extended, len(x_axis_list_log)*len(y_axis_list_log),exotic,dirac_width = dirac_width)

        #fill the table
        data_list = []
        for x in x_axis_list_log:
            variables_values[varX] = 10**x
            for x_var in scaleWithX.keys():
                variables_values[x_var] = scaleWithX[x_var]*10**x
            data_sublist = []
            for y in y_axis_list_log:
                variables_values[varY] = 10**y
                for y_var in scaleWithY.keys():
                    variables_values[y_var] = scaleWithY[y_var]*10**y
                if exotic == "alp" and "C" in '\t'.join(varAxes):
                    if varX == "mX":
                        data_sublist.append([10**x, 10**y/alp_model_params[0], process.ALP_events_EFT(variables_values['mX'], variables_values['CGG'], variables_values['CWW'], variables_values['CBB'], *alp_model_params, variables_values['Cll'], variables_values['Cqq'])])
                    elif varY == "mX":
                        data_sublist.append([10**x/alp_model_params[0], 10**y, process.ALP_events_EFT(variables_values['mX'], variables_values['CGG'], variables_values['CWW'], variables_values['CBB'], *alp_model_params, variables_values['Cll'], variables_values['Cqq'])])
                    else:
                        data_sublist.append([10**x, 10**y, process.ALP_events_EFT(variables_values['mX'], variables_values['CGG'], variables_values['CWW'], variables_values['CBB'], *alp_model_params, variables_values['Cll'], variables_values['Cqq'])])
                elif exotic == "hnl" and "U" in '\t'.join(varAxes):
                    data_sublist.append([10**x, 10**y, process.HNL_events(variables_values['mX'], variables_values['U2el'], variables_values['U2mu'], variables_values['U2tau'])])
                elif exotic == "ds" and "Y" in '\t'.join(varAxes):
                    data_sublist.append([10**x, 10**y, process.DS_events_EFT(variables_values['mX'], np.sqrt(variables_values['Y']), variables_values['Lambda'])])
                # elif exotic == "dp" and "eps2" in '\t'.join(varAxes):
                #     data_sublist.append([10**x, 10**y, process.DP_events_EFT(variables_values['mX'], np.sqrt(variables_values['eps2']))])
                elif exotic == "dp" and "eps" in '\t'.join(varAxes):
                    data_sublist.append([10**x, 10**y, process.DP_events_EFT(variables_values['mX'], variables_values['eps'])])
                
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
                else: raise "[Error:] \t Case not resolved"
            data_list.append([data_sublist])
        data = np.reshape(data_list,(len(x_axis_list_log)*len(y_axis_list_log),3))
        
        if not combine_experiments: #store output for experiments separately

            # export
            output_dir = this_dir+'/../tab_toPlot/'
            outPath = output_dir + exp_dict[exp] + '/' + exotic + '/'
            if not path.isdir(outPath): makedirs(outPath)
            
            # modes = ""
            # if channels_production: modes += "_" + '-'.join(channels_production)
            # if channels_decay: modes += "_" + '-'.join(channels_decay)

            # fixedNames = ""
            # for fixed in fixedValues.keys():
            #     fixed_val = float(fixedValues[fixed])
            #     if fixed == "mX":
            #         fixed_val*=fixed_val*1000
            #         if fixed_val.is_integer(): fixed_val = str(int(fixed_val))
            #         else: fixed_val = str("{:.1E}".format(fixed_val))
            #         fixedNames = fixedNames + "-" + fixed + str(int(fixedValues[fixed]*1000)) + "MeV"
            #     else:
            #         if fixed_val.is_integer(): fixed_val = str(int(fixed_val))
            #         else: fixed_val = str("{:.1E}".format(fixed_val))
            #         fixedNames = fixedNames + "-" + fixed + fixed_val
            # if fixedNames != "":
            #     fixedNames = "_fixed" + fixedNames

            # xscaledNames = ""
            # for couplings in scaleWithX.keys():
            #     coup_val = float(scaleWithX[couplings])
            #     if coup_val.is_integer(): coup_val = str(int(coup_val))
            #     else: coup_val = str("{:.1E}".format(coup_val))
            #     xscaledNames = xscaledNames + "-" + couplings + coup_val
            # if xscaledNames != "":
            #     xscaledNames = "_scaleWith" + varX + xscaledNames

            # yscaledNames = ""
            # for couplings in scaleWithY.keys():
            #     coup_val = float(scaleWithY[couplings])
            #     if coup_val.is_integer(): coup_val = str(int(coup_val))
            #     else: coup_val = str("{:.1E}".format(coup_val))
            #     yscaledNames = yscaledNames + "-" + couplings + coup_val
            # if yscaledNames != "":
            #     yscaledNames = "_scaleWith" + varY + yscaledNames

            # flags = ""
            # if dirac_width: flags += "_diracwidth"

            outfileName = exp + "_" + generate_file_info(channels_production,channels_decay,fixedValues, varX, varY, scaleWithX, scaleWithY, dirac_width, standalone)+'.dat'
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
            
    if combine_experiments: #export output for combined experiments
        outPath = this_dir+'/../tab_toPlot/combined/' + exotic + '/'
        if not path.isdir(outPath): makedirs(outPath)

        outfileName = "-".join(experiments) + "_" + generate_file_info(channels_production,channels_decay,fixedValues, varX, varY, scaleWithX, scaleWithY, dirac_width, standalone)+ '.dat'
        np.savetxt(outPath + outfileName,comb_data,fmt='%.4e')
        print('\n[Info:] \t', 'File ' + outfileName + ' saved to ' + outPath)


 ##
 # @fn generate_file_info
 # @brief prepares an output filename
def generate_file_info(channels_production,channels_decay,fixedValues, varX, varY, scaleWithX, scaleWithY, dirac_width, standalone):
    modes = ""
    if channels_production: modes += "_" + '-'.join(channels_production)
    if channels_decay: modes += "_" + '-'.join(channels_decay)

    fixedNames = ""
    for fixed in fixedValues.keys():
        if fixed == "mX":
            fixedNames = fixedNames + "-" + fixed + str(int(fixedValues[fixed]*1000)) + "MeV"
        else:
            fixedNames = fixedNames + "-" + fixed + str(int(fixedValues[fixed]))
    if fixedNames != "":
        fixedNames = "_fixed" + fixedNames

    xscaledNames = ""
    for couplings in scaleWithX.keys():
        coup_val = float(scaleWithX[couplings])
        if coup_val.is_integer(): coup_val = str(int(coup_val))
        else: coup_val = str("{:.1E}".format(coup_val))
        xscaledNames = xscaledNames + "-" + couplings + coup_val
    if xscaledNames != "":
        xscaledNames = "_scaleWith" + varX + xscaledNames

    yscaledNames = ""
    for couplings in scaleWithY.keys():
        coup_val = float(scaleWithY[couplings])
        if coup_val.is_integer(): coup_val = str(int(coup_val))
        else: coup_val = str("{:.1E}".format(coup_val))
        yscaledNames = yscaledNames + "-" + couplings + coup_val
    if yscaledNames != "":
        yscaledNames = "_scaleWith" + varY + yscaledNames
        
    flags = ""
    if standalone: modes += "_standalone"
    if dirac_width: flags += "_diracwidth"

    return  varX + "_" + varY + fixedNames + xscaledNames + yscaledNames + modes +flags


 ##
 # @fn generate_log_list
 # @brief Get a log list for given variable
def generate_log_list(var,Lambda):
    n_bins = 601
    if var == "mX": #GeV
        # list = [ m_a_log for m_a_log in np.linspace(np.log10(0.0001),np.log10(3.01),n_bins)]
        return np.linspace(np.log10(0.01),np.log10(5.35),n_bins).tolist()
    elif "C" in var: #GeV^-1
        n_bins = 201
        return np.linspace(-11+np.log10(Lambda),-0.999+np.log10(Lambda),n_bins).tolist()
    elif "U" in var: 
        n_bins = 201
        return  np.linspace(-11,-1,n_bins).tolist()
    elif "Y" in var:
        return np.linspace(-13,-1,n_bins).tolist()
    elif "eps" in var:
        n_bins = 101
        return np.linspace(-8,-2,n_bins).tolist()
    elif "Lambda" in var:
        n_bins = 101
        return np.linspace(-13,-1,n_bins).tolist()
    elif var == "tauX": #ps
        tau_bin_width = 0.08
        return np.arange(-2,6,tau_bin_width).tolist()
    elif var == "GammaX": #GeV
        Gamma_bin_width = 0.15
        return np.arange(-25,-10,Gamma_bin_width).tolist()
    elif "BR" in var:
        n_bins = 101
        return np.linspace(-13,-1,n_bins).tolist()
    # elif "eps2" in var:
    #     n_bins = 101
    #     list = [ eps2 for eps2 in np.linspace(-16,-4,n_bins)]
    else:
        print("[Error:] \t Variable " + var + " not recognized")
        sys.exit(1)
    return []

 ##
 # @fn generate_lin_list
 # @brief Get a linear list in range for given bins
def generate_lin_list(min,max,nbins):
    return [ m_a_log for m_a_log in np.linspace(min,max,nbins)]

if __name__ == "__main__":
    sys.exit(main())
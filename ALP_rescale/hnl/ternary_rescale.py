#!/usr/bin/env python3

from os import path, makedirs
import numpy as np
from ALP_rescale.hnl import hnl_setup as setup
from ALP_rescale.general import load_data as ld
#.. from general.mergeSigRegions import MergeInput
# import argparse
import sys


def main(exps, prods=[], decays =[], combine=False ):
    '''Option interpretation for ternary yield evaluation'''

    exo_mass = float(input(" - Enter fixed HNL mass (in GeV): "))
    if exo_mass < 0:
        raise ValueError("[Error:] \t Please enter only positive numbers")
    print("[Info:] \t Mass mX fixed value: " + str(exo_mass) +"GeV")

    #experiment:
    experiments = []
    if exps == "":
        experiments = setup.experiments.keys()
        print("[Info:] \t Selected all experiments available")
    else:
        for exp in exps:
            if exp in  setup.experiments.keys():
                experiments.append(exp)
            else:
                raise("[Error:] \t Experiment " + exp + " not available. Experiment modes available: exp = " + ' | '.join(setup.experiments.keys()) + ". If not specified, running over all experiments available.")
        print("[Info:] \t Selected experiments:", ', '.join(experiments))
    if combine:
        print("[Info:] \t Datasets for selected experiments will be combined")
    #production mode:
    channels_production = []
    if prods == "":
        channels_production = setup.channels_production
        print("[Info:] \t Selected all production modes available")
    else:
        for prod in prods:
            if prod in setup.channels_production:
                channels_production.append(prod)
            else:
                raise("[Error:] \t Production mode " + prod + " not available. Production modes available: prod = " + ' | '.join(setup.channels_production) + ". If not specified, running over all production modes available.")
        print("[Info:] \t Selected production modes:", ', '.join(channels_production))

    #decay mode:
    channels_decay = []
    if decays:
        for dec in decays:
            if dec in setup.channels_decay: channels_decay.append(dec)
            else: raise("[Error:] \t Decay mode " + dec + " not available. Decay modes available: decay = " + ' | '.join(setup.channels_decay) + ". If not specified, running over all decay modes available.")
        print("[Info:] \t Selected decay modes:", ', '.join(channels_decay))
    else:
        channels_decay = setup.channels_decay
        print("[Info:] \t Selected all decay modes available.")

    this_dir = path.dirname(path.realpath(__file__))
    
    variables_values = {}
    for var in setup.variables + setup.couplings: #initialize with zeros
        variables_values[var] = 0 

    #generate x- and y- tables
    axis_base = generate_lin_list(0,1,21).tolist()
    U2_scales = np.power(10, generate_log_list("U2")).tolist()

    #rescale

    channels_production_extended = []
    for channel in channels_production:
        for mixing in ["El","Mu","Tau"]:
            channels_production_extended.append(channel+ "-"+mixing+"Mixing")

    for exp in experiments:
        #.. if regions != [""]: mergeSigReg = MergeInput(exp,regions,channels_decay,channels_production)
        el_axis, mu_axis, tau_axis = uniform_triangles(axis_base)
        process = ld.Process_data(exp, channels_decay, channels_production_extended, len(el_axis)*len(U2_scales), "hnl")
        #fill the table
        data_list = []
        for U2_scale in U2_scales:
            data_sublist = []
            for el_point, mu_point, tau_point in zip(el_axis, mu_axis, tau_axis):
                data_sublist.append([U2_scale, el_point, mu_point, process.HNL_events(exo_mass, U2_scale*el_point, U2_scale*mu_point, U2_scale*tau_point)])
            data_list.append([data_sublist])
        data = np.reshape(data_list, (len(el_axis)*len(U2_scales), 4))

        if not combine: #store output for experiments separately

            # export
            output_dir = this_dir+'/../../tab_toPlot/'
            outPath = output_dir + '/'.join( [setup.experiments[exp], "hnl",  ""] )

            modes = ""
            if prods: modes += "_" + '-'.join(channels_production)
            if decays: modes += "_" + '-'.join(channels_decay)
            
            outfileName = exp + "_U2ternary_mass"+str(int(exo_mass*1000)) + "MeV"+ modes +".dat"
            if not path.isdir(outPath): makedirs(outPath)
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
            
    if combine: #export output for combined experiments
        outPath = this_dir+'/../../tab_toPlot/combined/hnl/'

        modes = ""
        if prods: modes += "_" + '-'.join(channels_production)
        if decays: modes += "_" + '-'.join(channels_decay)

        outfileName = "-".join(experiments) +  "_U2ternary_mass"+str(int(exo_mass*1000)) + "MeV"+ modes + '.dat'
        if not path.isdir(outPath): makedirs(outPath)
        np.savetxt(outPath + outfileName,comb_data,fmt='%.4e')
        print('\n[Info:] \t', 'File ' + outfileName + ' saved to ' + outPath)

def uniform_triangles(basis):
    x = np.array([])
    y = np.array([])
    N = np.size(basis)
    for i in range(N):
        x = np.append(x,np.repeat(basis[i], N-i))
        y = np.append(y, basis[:(N-i)] )
    z = np.clip(1 - x - y, 0, None)
    return x, y, z

def generate_log_list(var):
    n_bins = 201
    list = []
    if var == "U2": 
        n_bins = 201
        list =  np.linspace(-11,-1,n_bins)
    else:
        print("[Error:] \t Variable " + var + " not recognized")
        sys.exit(1)

    return list

def generate_lin_list(min,max,nbins):
    return np.linspace(min,max,nbins)

if __name__ == "__main__":
    sys.exit(main())
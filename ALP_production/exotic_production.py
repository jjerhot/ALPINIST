#!/usr/bin/env python3

from ALP_production.general import exotic_production_setup as setup
from ALP_production.modes import mixing as mix
from ALP_production.modes import meson_decay as meson
from ALP_production.modes import dp_brems as brems
from ALP_production.modes import photoproduction
import argparse
import sys
from multiprocessing import cpu_count
from ALP_production.general.exotic_functions import Input_reweight as ir 
def main(argv=None):
    '''Command line options for exotic production.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    parser = argparse.ArgumentParser(description='Exotic production in fixed target. \n Select exotic particle class, experiment and production mode')
    parser.add_argument("-x","--exo",  default="", type=str, nargs='*', help="Exotic particle classes available (case sensitive): exo = "+ " | ".join(setup.channels_production.keys())+".")
    parser.add_argument("-e","--exp",  default="", type=str, nargs='*', help="Experiments available (case sensitive): exp = " + ' | '.join(setup.experiments.keys()) + ".")
    parser.add_argument("-p","--prod",  default="", type=str, nargs='*', help = "Production modes available (case sensitive, subject to chosen Exotic): prod = Mixing | Primakoff | PhotonFromMeson | Dmeson | Bmeson | Bmeson2S | MesonDecay | Brems.")
    parser.add_argument("-n","--nthreads",  default=1, type=int, help="Number of threads. Default = 1")
    parser.add_argument("-nprod","--num-mesons", dest='nmesons',  default=1E5, type=int, help="Number of generated mesons. Default = 1E5")
    parser.add_argument("-ndec","--num-decays", dest='ndecays',  default=10, type=int, help="Number of meson decay events (relevant only for Bmeson and Dmeson production). Default = 10")
    parser.add_argument("-ext","--external-source", dest='ext', action='store_true', help="Use external source. Default")
    parser.add_argument("-no-ext","--no-external-source", dest='ext', action='store_false', help="Use Pythia as a source")
    parser.add_argument("-3partprod","--include-three-parton-production", dest='three_part', default=False,  action='store_true', help="include meson production through 3 parton processes (only available with external source). WARNING: 3 parton production labled `incomplete` for Pythia 8.3.")
    # parser.add_argument("-reweight","--reweight-external-source", default="", type=str, nargs=1, dest='rew', help="Reweight Pythia charmed meson maps based on measured xF distributions. Available options: " + ' | '.join(ir.literature.keys()) + ".")
    parser.add_argument("-empirical","--empirical-meson-distribution", default="", type=str, nargs=1, dest='empiric', help="Use empirical differential cross sections to generate charmed meson maps. Available options: "+' | '.join(ir.literature.keys()) + ".")
    parser.set_defaults(ext=True)

    args = parser.parse_args()

    #exotic for which to run production
    exotics = []
    if args.exo == "":
        parser.error("-x/--exo is a required parameter\n Exotic particle classes available: exo = "+ " | ".join(setup.channels_production.keys())+".")
    else: 
        for exo in args.exo:
            if exo in setup.channels_production.keys():
                exotics.append(exo)
            else:
                parser.error("Exotic "+exo+" not available.\n Exotic particle classes available: exo = "+ " | ".join(setup.channels_production.keys())+".")

    #experiment:
    experiments = []
    if args.exp == "":
        parser.error("-e/--exp is a required parameter\n Experiment modes available: exp = " + ' | '.join(setup.experiments.keys()) + ".")
    else:
        for exp in args.exp:
            if exp in setup.experiments.keys():
                experiments.append(exp)
            else:
                parser.error("Experiment " + exp + " not available.\n Experiment modes available: exp = " + ' | '.join(setup.experiments.keys()) + ".")
        print("[Info:] \t Selected experiments:", ', '.join(experiments))

    #production mode:
    channels_production = []
    for exo in exotics:
        channels_production_exo = []
        if args.prod == "":
            parser.error("-p/--prod is a required parameter\n Production modes available for "+ exo +": prod = " + ' | '.join(setup.channels_production[exo]) + ".")
        else:
            for prod in args.prod:
                if prod in setup.channels_production[exo]:
                    channels_production_exo.append(prod)
                else:
                    parser.error("Production mode " + prod + " not available.\n Production modes available: prod = " + ' | '.join(setup.channels_production[exo]) + ".")
            print("[Info:] \t Selected production modes:", ', '.join(channels_production_exo))
        channels_production.append(channels_production_exo)

    if args.three_part and not args.ext:
        print("[Error:] \t 3 parton processes only available with external source.")
        sys.exit(1)

    #number of threads:
    n_cores_avail = cpu_count()
    nthreads = max(1, args.nthreads)
    response_overwrite_n = "n"
    while n_cores_avail < nthreads and (response_overwrite_n!="y" or response_overwrite_n!="Y"):
        print("[Info:] \t Selected number of threads (",args.nthreads, ") exceeds number of available CPU cores found (",n_cores_avail,").")
        response_overwrite_n = input("\tProceed with selected number of threads anyway? [Y/N]:\t")
        if response_overwrite_n == "y" or response_overwrite_n =="Y": nthreads = args.nthreads
        elif response_overwrite_n == "n" or response_overwrite_n == "N": 
            nthreads = max(int(input("\tPlease select new desired number of cores:\t")), 1)
        else:
            raise KeyError("[Error:] \t Invalid answer: " + response_overwrite_n + ". Type [Y/N] or [y/n]")

    #additional arguments
    nevents = int(args.nmesons)
    ndecays = int(args.ndecays)
    use_ext = args.ext
    use_3parton_production = args.three_part
    use_external_reweight ="" # args.rew[0] if args.rew else 
    use_empirical_dists = args.empiric[0] if args.empiric else ""


    for iex in range(len(exotics)):
        exo = exotics[iex]
        active_coupling = ["El", "Mu", "Tau"] if exo == "hnl" else [""] # looping over coupling dominance scenarios for hnls
        for exp in experiments:
            for prod in channels_production[iex]: # production channel associated with exotic
                print("[Info:] \t Running production for",exp,"experiment and",prod,"production mode")
                if prod == "Bmeson": 
                    meson.bmeson_decay_production(exp,nevents,ndecays,exo,use_ext, use_3part_prod = use_3parton_production).process_pool(nthreads, active_coupling)
                elif prod == "Bmeson2S":
                    meson.bmeson_decay_production(exp,nevents,ndecays,exo,use_ext).process_pool_2s(nthreads)
                elif prod == "Dmeson": 
                    for use_empirical_dist in use_empirical_dists.split('+'): # convenience for comparing multiple empirical distributions
                        meson.dmeson_decay_production(exp,nevents,ndecays,exo, use_ext, use_3part_prod = use_3parton_production, use_reweight=use_external_reweight, use_empirical=use_empirical_dist).process_pool(nthreads, active_coupling)
                elif prod == "MesonDecay":
                    meson.meson_to_dp_decay_production(exp,nevents,ndecays,exo,use_ext).process_pool(nthreads)
                elif prod == "Mixing":
                    if exo == "dp":  mix.mixing_production_vector(exp,nevents,exo,use_ext).process_pool(nthreads)
                    if exo == "alp": mix.mixing_production(exp,nevents,exo,use_ext).process_pool(nthreads)
                elif prod == "Brems":
                    brems.brems_production(exp,exo).process_pool(nthreads)
                elif prod == "Primakoff":
                    photoproduction.primakoff_production(exp,exo).process_pool(nthreads)
                elif prod == "PhotonFromMeson":
                    photoproduction.photon_from_meson_production(exp,exo).process_pool(nthreads)

if __name__ == "__main__":
    sys.exit(main())
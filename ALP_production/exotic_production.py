#!/usr/bin/env python3

from general import exotic_production_setup as setup
from general import exotic_constants as c
from modes import mixing as mix
from modes import meson_decay as meson
import argparse
import sys

def main(argv=None):
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    parser = argparse.ArgumentParser(description='ALP production in fixed target. \n Select the experiment and production mode')
    parser.add_argument("-e","--exp",  default="", type=str, nargs='*', help="Experiments available (case sensitive): exp = NA62 | CHARM | NuCal | SHiP | DarkQuest | DUNE | SHADOWS | KOTO | KOTO2.")
    parser.add_argument("-p","--prod",  default="", type=str, nargs='*', help="Production modes available (case sensitive): prod = primakoff | photonfrommeson | mixing | Bmeson | Dmeson.")
    parser.add_argument("-n","--nthreads",  default=1, type=int, help="Number of threads. Default = 1")
    parser.add_argument("-nprod","--num-mesons", dest='nmesons',  default=1E5, type=int, help="Number of generated mesons. Default = 1E5")
    parser.add_argument("-ndec","--num-decays", dest='ndecays',  default=10, type=int, help="Number of meson decay events (relevant only for Bmeson and Dmeson production). Default = 10")
    parser.add_argument("-ext","--external-source", dest='ext', action='store_true', help="Use external source. Default")
    parser.add_argument("-no-ext","--no-external-source", dest='ext', action='store_false', help="Use Pythia as a source")
    parser.set_defaults(ext=True)

    args = parser.parse_args()

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
    if args.prod == "":
        parser.error("-p/--prod is a required parameter\n Production modes available: prod = " + ' | '.join(setup.channels_production) + ".")
    else:
        for prod in args.prod:
            if prod in setup.channels_production:
                channels_production.append(prod)
            else:
                parser.error("Production mode " + prod + " not available.\n Production modes available: prod = " + ' | '.join(setup.channels_production) + ".")
        print("[Info:] \t Selected production modes:", ', '.join(channels_production))

    #number of threads:
    nthreads = 1
    if args.nthreads < 1:
        nthreads = 1
    else:
        nthreads = args.nthreads

    #additional arguments
    nevents = int(args.nmesons)
    ndecays = int(args.ndecays)
    use_ext = args.ext

    for exp in experiments:
        for prod in channels_production:
            print("[Info:] \t Running production for",exp,"experiment and",prod,"production mode")
            if prod == "mixing":
                mix.mixing_production(exp,nevents,use_ext).process_pool(nthreads)
                # mix.mixing_production_var_mass(exp,nevents).process_pool(nthreads)
            elif prod == "Bmeson":
                meson.bmeson_decay_production(exp,nevents,ndecays,"alp",use_ext).process_pool(nthreads)
            elif prod == "Dmeson":
                meson.dmeson_decay_production(exp,nevents,ndecays,"alp",use_ext).process_pool(nthreads)

if __name__ == "__main__":
    sys.exit(main())
#!/usr/bin/env python3

from ALP_production.general import exotic_production_setup as setup
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
    parser.add_argument("-trapz","--use-trapz", dest='trapz', default=False, action='store_true', help="Use trapezoid integration (relevant only for photoproduction). Speedup by a factor 100 at the cost of precision.")
    parser.add_argument("-vegas","--use-vegas", dest='vegas', default=False, action='store_true', help="Use vegas integration (relevant only for photoproduction). Preferable when strong GPU available. Requires tensorflow and vegasflow packages.")
    parser.add_argument("-nprod","--num-mesons", dest='nmesons',  default=1E5, type=int, help="Number of generated mesons (relevant only for B/D/B2S and light meson decay and mixing production). Default = 1E5")
    parser.add_argument("-ndec","--num-decays", dest='ndecays',  default=10, type=int, help="Number of meson decay events (relevant only for B/D/B2S and light meson decay production). Default = 10")
    parser.add_argument("-ext","--external-source", dest='ext', action='store_true', help="Use external source (relevant only for B/D/B2S and light meson decay and mixing production). (Default option)")
    parser.add_argument("-no-ext","--no-external-source", dest='ext', action='store_false', help="Use Pythia as a source (relevant only for B/D/B2S and light meson decay and mixing production)")
    parser.add_argument("-3partprod","--include-three-parton-production", dest='three_part', default=False,  action='store_true', help="include meson production through 3 parton processes (only available with external source). WARNING: 3 parton production labled `incomplete` for Pythia 8.3.")
    # parser.add_argument("-reweight","--reweight-external-source", default="", type=str, nargs=1, dest='rew', help="Reweight Pythia charmed meson maps based on measured xF distributions. Available options: " + ' | '.join(ir.literature.keys()) + ".")
    parser.add_argument("-empirical","--empirical-meson-distribution", default="", type=str, nargs=1, dest='empiric', help="Use empirical differential cross sections to generate charmed meson maps. Available options: "+' | '.join(ir.literature.keys()) + ".")
    parser.add_argument("-ac","--active-coupling", default=['El','Mu','Tau'], type=str, nargs='*', dest='active_coupling', help="Active coupling to run simulation for (relevant only for hnl production). Default = 'El Mu Tau'.")
    parser.add_argument("-fix-en","--fix-energy", dest='fix_en', action='store_true', help="Fix energy in CM frame of original meson in mixing production")
    parser.add_argument("-fix-mom","--fix-momentum", dest='fix_en', default=True, action='store_false', help="Fix energy in CM frame of original meson in mixing production")
    parser.add_argument("-pick-mass","--pick-mass", dest='single_mass_point', default=0, type=float, help="Pick a single mass point [GeV] for which to run the desired production")
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
    
    for ac in args.active_coupling:
        if ac not in ['El', 'Mu', 'Tau'] and exo == "hnl":
            parser.error("Active coupling " + ac + " not recognised.\n Active couplings: ac = " + ' | '.join(['El', 'Mu', 'Tau']) + ".")
    active_coupling = args.active_coupling if exo == "hnl" else [""]

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
    print("[Info:] \t Selected number of cores:", str(nthreads))
    if args.trapz:
        print("[Info:] \t Using trapezoid integration for photoproduction.")

    #additional arguments
    nevents = int(args.nmesons)
    ndecays = int(args.ndecays)
    use_ext = args.ext
    use_3parton_production = args.three_part
    use_external_reweight ="" # args.rew[0] if args.rew else 
    use_empirical_dists = args.empiric[0] if args.empiric else ""
    mixing_frame = args.fix_en
    single_mass_point = args.single_mass_point

    for iex in range(len(exotics)):
        exo = exotics[iex]
        # active_coupling = ["El", "Mu", "Tau"] if exo == "hnl" else [""] # looping over coupling dominance scenarios for hnls
        for exp in experiments:
            for prod in channels_production[iex]: # production channel associated with exotic
                print("[Info:] \t Running production for",exp,"experiment and",prod,"production mode")
                if prod == "Bmeson": 
                    from ALP_production.modes import meson_decay as meson
                    meson.bmeson_decay_production(exp,nevents,ndecays,exo,use_ext, use_3part_prod = use_3parton_production).process_pool(nthreads, active_coupling, single_mass_point=single_mass_point)
                elif prod == "Bmeson2S":
                    from ALP_production.modes import meson_decay as meson
                    meson.bmeson_decay_production(exp,nevents,ndecays,exo,use_ext).process_pool_2s(nthreads, single_mass_point=single_mass_point)
                elif prod == "Dmeson": 
                    from ALP_production.modes import meson_decay as meson
                    for use_empirical_dist in use_empirical_dists.split('+'): # convenience for comparing multiple empirical distributions
                        meson.dmeson_decay_production(exp,nevents,ndecays,exo, use_ext, use_3part_prod = use_3parton_production, use_reweight=use_external_reweight, use_empirical=use_empirical_dist).process_pool(nthreads, active_coupling, single_mass_point=single_mass_point)
                elif prod == "MesonDecay":
                    from ALP_production.modes import meson_decay as meson
                    meson.meson_to_dp_decay_production(exp,nevents,ndecays,exo,use_ext).process_pool(nthreads, single_mass_point=single_mass_point)
                elif prod == "Mixing":
                    from ALP_production.modes import mixing as mix
                    mix.mixing_production(exp,nevents,exo,use_ext,fix_energy_in_cm=mixing_frame).process_pool(nthreads, single_mass_point=single_mass_point)
                elif prod == "Brems":
                    if exo == 'dp':
                        from ALP_production.modes import dp_brems as bb_brems
                        bb_brems.brems_production(exp,exo).process_pool(nthreads, single_mass_point=single_mass_point) # Blumlein-Brunner
                    else:
                        from ALP_production.modes import bremsstrahlung as brems
                        brems.brems_production(exp,exo,Lambda = 1.5, z_min = 0.1, z_max = 0.9).process_pool(nthreads, single_mass_point=single_mass_point)
                elif prod == "Primakoff":
                    from ALP_production.modes import photoproduction
                    photoproduction.primakoff_production(exp,exo,use_trapz=args.trapz, use_vegas = args.vegas).process_pool(nthreads, single_mass_point=single_mass_point)
                elif prod == "PhotonFromMeson":
                    from ALP_production.modes import photoproduction
                    photoproduction.photon_from_meson_production(exp,exo,use_trapz=args.trapz, use_vegas = args.vegas).process_pool(nthreads, single_mass_point=single_mass_point)

if __name__ == "__main__":
    sys.exit(main())
import sys, os

pythia_dir = None

if pythia_dir is not None:
  print("[Info:] \t Using custom Pythia directory: "+pythia_dir)
else:
  pythia_dir = os.environ.get('PYTHIA8', None)
  
  if pythia_dir is None:
    print("[Warning:] \t Pythia not found. Please set PYTHIA8 to the corresponding directory if it should be used for exotic production")
  else:
    if not os.path.exists(pythia_dir + "/examples"):
      if os.path.exists(pythia_dir + "/share"): # probably conda environment
        pythia_dir += "/share/Pythia8/"
      else:
        print("[Warning:] \t Pythia folder not containing examples directory")


#currently available experiments and their tags. 
experiments = { 'NA62':         'NA62',
                'CHARM':        'CHARM',
                'NuCal':        'NuCal',
                'SHiP':         'SHiP',
                'DarkQuest':    'DarkQuest',
                'DUNE':         'DUNE',
                'SHADOWS':      'SHADOWS',
                'KOTO':         'KOTO',
                'KOTO2':        'KOTO2',
                'NuTeV':        'NuTeV',
                'BEBC':         'BEBC',
                'ORCA':         'ORCA'}

#available production channels for the exotic particles 
channels_production = {"alp" : [
                        'Primakoff',
                        'PhotonFromMeson',
                        'Mixing',
                        'Bmeson',
                        'Dmeson'
                        ],
                      "hnl" : [
                          'Bmeson',
                          'Dmeson'
                        ],
                      "dp" : [
                          'Brems',
                          'Mixing',
                          'MesonDecay',
                        ],
                      "ds" : [
                          'Brems',
                          'Bmeson',
                          'Bmeson2S'
                        ]
}

coupling_ref = 1
coupling_exp = 2

# beam momenta at the experimental facilities
p_beam =    {'NuTeV':   800,
             'NA62':    400,
             'CHARM':   400, 
             'NuCal':   70, 
             'SHiP':    400, 
             'DarkQuest':120, 
             'DUNE':    120,
             'SHADOWS': 400, 
             'KOTO':    30, 
             'KOTO2':   30,
             'BEBC':    400,
             'ORCA':    400,
            }

# minimum energy for export grid at given beam momentum (in GeV)
energy_min ={ #0.59: 0.00984,
              30: 0.5,
              70: 1.1,
              120: 1.9,
              400: 6.5,
              800: 13.,
            }
# maximum energy for export grid at given beam momentum (in GeV)
energy_max ={ #0.59: 0.58,
              30: 29.5, 
              70: 64.9,
              120: 112.1,
              400: 383.5,
              800: 767.,
            }

energy_bins = 30

# minimum opening angle from proton beam for export grid at given experiment (in rad)
theta_min = { 
              "CHARM": 0.0025,
              "NuCal": 0.00025,
              "SHiP": 0.00078,
              "DarkQuest": 0.0009,
              "SHADOWS": 0.012,
              "KOTO": 0.06, 
              "KOTO2": 0.06,
              "NuTeV": 0.0000295,
              "BEBC": 0.00017,
              "ORCA": 1.5e-6,
            } #default: 0.00018
# maximum opening angle from proton beam for export grid at given experiment (in rad)
theta_max = { 
              "CHARM": 0.0175,
              "NuCal": 0.01525,
              "SHiP": 0.04758,
              "DarkQuest": 0.0549,
              "SHADOWS": 0.12,
              "KOTO": 0.6,
              "KOTO2": 0.12,
              "BEBC": 0.00972,
              "NuTeV": 0.00179,
              "ORCA": 90e-6,
            } #default: 0.01098

kaons_list = ["K","K0star_700","K0star_1430","Kstar_892","Kstar_1410","Kstar_1680","K1_1270","K1_1400","K2star_1430"]

theta_bins = 31

mass_min = [0.0001,0.01]
mass_max = [0.01,5.31]
mass_bins = [101,201]
#.. for linear bins:
# mass_min = [0.0001,0.01]
# mass_max = [0.01,3.01]
# mass_bins = [101,601]

# eStep[beamMom_] := 
#   Switch[beamMom, 30, 1, 70, 2.2, 120, 3.8, 400, 11, _, 
#    Echo["unknown Pbeam, Estep set to 0"]; 0];

# thStep[exp_] := 
#   Switch[exp, "CHARM", 0.0005, "NuCal", 0.0005, "SHiP", 0.00156, 
#    "DarkQuest", 0.0018, "SHADOWS", 0.0036, "KOTO", 0.018, "KOTO2", 
#    0.002, _, 0.00036];


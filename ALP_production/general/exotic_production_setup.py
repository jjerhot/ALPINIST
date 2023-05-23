pythia_dir = "/home/jjan/cernbox/na62/pythia"

experiments = { 'NA62':         'NA62',
                'CHARM':        'CHARM',
                'NuCal':        'NuCal',
                'SHiP':         'SHiP',
                'DarkQuest':    'DarkQuest',
                'DUNE':         'DUNE',
                'SHADOWS':      'SHADOWS',
                'KOTO':         'KOTO',
                'KOTO2':        'KOTO2'}

channels_production = [ 'primakoff',
                        'photonfrommeson',
                        'mixing',
                        'Bmeson',
                        'Dmeson'
                        ]
reference_couplings = {'primakoff':      1e-4,
                       'photonfrommeson':1e-4,
                       'mixing':         1e-4,
                       'Bmeson':         1e-10,
                       'Dmeson':         1e-10}

scaling_exponent =    {'primakoff':      2,
                       'photonfrommeson':2,
                       'mixing':         2,
                       'Bmeson':         1,
                       'Dmeson':         1}

p_beam =    {'NA62':    400,
             'CHARM':   400, 
             'NuCal':   70, 
             'SHiP':    400, 
             'DarkQuest':120, 
             'DUNE':    120,
             'SHADOWS': 400, 
             'KOTO':    30, 
             'KOTO2':   30
            }

energy_min ={ 30: 0.5,
              70: 1.1,
              120: 1.9,
              400: 5.5,
            }

energy_max ={ 30: 29.5, 
              70: 64.9,
              120: 112.1,
              400: 324.5,
            }

energy_bins = 30

theta_min = { "CHARM": 0.0025,
              "NuCal": 0.00025,
              "SHiP": 0.00078,
              "DarkQuest": 0.0009,
              "SHADOWS": 0.012,
              "KOTO": 0.06,
              "KOTO2": 0.06,
            } #default: 0.00018

theta_max = { "CHARM": 0.0175,
              "NuCal": 0.01525,
              "SHiP": 0.04758,
              "DarkQuest": 0.0549,
              "SHADOWS": 0.12,
              "KOTO": 0.6,
              "KOTO2": 0.12
            } #default: 0.01098

theta_bins = 31

mass_min = [0.0001,0.01]
mass_max = [0.01,3.01]
mass_bins = [101,601]

# eStep[beamMom_] := 
#   Switch[beamMom, 30, 1, 70, 2.2, 120, 3.8, 400, 11, _, 
#    Echo["unknown Pbeam, Estep set to 0"]; 0];

# thStep[exp_] := 
#   Switch[exp, "CHARM", 0.0005, "NuCal", 0.0005, "SHiP", 0.00156, 
#    "DarkQuest", 0.0018, "SHADOWS", 0.0036, "KOTO", 0.018, "KOTO2", 
#    0.002, _, 0.00036];
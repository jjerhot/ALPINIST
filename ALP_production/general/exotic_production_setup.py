pythia_dir = None

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
                          'Bmeson',
                          'Bmeson2S'
                        ]
}

coupling_ref = 1
coupling_exp = 2

#legacy reference couplings and scaling exponents
# reference_couplings = {'primakoff':      1e-4,
#                        'photonfrommeson':1e-4,
#                        'mixing':         1e-4,
#                        'Bmeson':         1e-4,
#                        'Dmeson':         1e-4}
# scaling_exponent =    {'primakoff':      2,
#                        'photonfrommeson':2,
#                        'mixing':         2,
#                        'Bmeson':         2,
#                        'Dmeson':         2}

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

# target material volume averaged proton count per nucleus
Z_target =  {'NuTeV':   6.55872,
             'NA62':    29,
             'CHARM':   29, 
             'NuCal':   26, 
             'SHiP':    42, 
             'DarkQuest':26, 
             'DUNE':    6,
             'SHADOWS': 29, 
             'KOTO':    79, 
             'KOTO2':   79,
             'BEBC':    29,
             'ORCA':    6,
            }

# target material volume averaged nucleon count per nucleus
A_target =  {'NuTeV':   12.5056,
             'NA62':    63.546,
             'CHARM':   63.546, 
             'NuCal':   56, 
             'SHiP':    95, 
             'DarkQuest':56, 
             'DUNE':    12,
             'SHADOWS': 63.546, 
             'KOTO':    197, 
             'KOTO2':   197,
             'BEBC':    63.546,
             'ORCA':    12,
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


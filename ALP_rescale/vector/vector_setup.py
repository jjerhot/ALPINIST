from ALP_rescale.general.setup import *

variables = [   'mX',
                'GammaX',
                'tauX',
            ]

couplings = [   'eps'
            ]

experiments = { 'NA62':         'NA62',
                'CHARM':        'CHARM',
                'NuCal':        'NuCal',
                'SHiP':         'SHiP',
                'SHiPecn3':     'SHiP',
                'DarkQuest':    'DarkQuest',
                'DUNE':         'DUNE',
                'SHADOWS':      'SHADOWS',
                'NuTeV':        'NuTeV',
                'ORCA':         'ORCA'}

channels_decay = [  '2El',
                    '2Mu',
                    '2Pi',
                    '3Pi',
                    '4Pi',
                    '2Pi2Pi0',
                    '2K',
                    '2KPi0'
                    ]

channels_production = [ 
                        'Brems',
                        'MesonDecay',
                        'Mixing',
                        ]

reference_couplings = {
                        'Brems':        1,
                        'MesonDecay':   1,
                        'Mixing':       1,
                        }
scaling_exponent =    {
                        'Brems':        2,
                        'MesonDecay':   2,
                        'Mixing':       2,
                        }

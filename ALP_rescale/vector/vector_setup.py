from ALP_rescale.general.setup import *

variables = [   'mX',
                'GammaX',
                'tauX',
            ]

couplings = [   'eps'
            ]

experiments = { 'NA62':         'NA62',
                'HIKE':         'NA62',
                'CHARM':        'CHARM',
                'NuCal':        'NuCal',
                'SHiP':         'SHiP',
                'SHiPecn4':     'SHiP',
                'DarkQuest':    'DarkQuest',
                'DarkQuestPhase2':'DarkQuest',
                'DUNE':         'DUNE',
                'SHADOWS':      'SHADOWS',
                'NuTeV':        'NuTeV',
                'BEBC':         'BEBC',
                'BEBCcuboid':   'BEBC',
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
                        'MixingRho',
                        'MixingOmega',
                        'MixingPhi',
                        ]

reference_couplings = {
                        'Brems':        1,
                        'MesonDecay':   1,
                        'MixingRho':    1,
                        'MixingOmega':  1,
                        'MixingPhi':    1,
                        }
scaling_exponent =    {
                        'Brems':        2,
                        'MesonDecay':   2,
                        'MixingRho':    1,
                        'MixingOmega':  1,
                        'MixingPhi':    1,
                        }

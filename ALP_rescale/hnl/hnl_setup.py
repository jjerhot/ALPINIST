from ALP_rescale.general.setup import *

couplings = [   
            # 'U2',
            'U2el',
            'U2mu',
            'U2tau',
            ]

experiments = { 'NA62':             'NA62',
                'HIKE':             'NA62',
                'CHARM':            'CHARM',
                'NuCal':            'NuCal',
                'SHiP':             'SHiP',
                'SHiPecn3':         'SHiP',
                'DarkQuest':        'DarkQuest',
                'DarkQuestPhase2':  'DarkQuest',
                'DUNE':             'DUNE',
                'SHADOWS':          'SHADOWS',
                'NuTeV':            'NuTeV',
                'BEBC':             'BEBC',
                'BEBCcuboid':       'BEBC'}


channels_decay = [  
                    'PiEl',
                    'PiPiEl',
                    'PiMu',
                    'PiPiMu',
                    'NuElMu',
                    'NuElEl',
                    'NuMuMu',
                    'PiNu',
                    'EtaNu',
                    'PiPiNu'
                ]

channels_production =   [ 
                        'Bmeson',
                        'Dmeson',
                        # 'recast',
                        # 'NoCuts'
                        ]

scaling_exponent = 1
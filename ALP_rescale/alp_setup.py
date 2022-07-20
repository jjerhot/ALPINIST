variables = [   'ma',
                'CBB',
                'CWW',
                'CGG',
                'Cll',
                'Cqq'
            ]

experiments = [ 'NA62',
                'CHARM',
                'NuCal',
                'SHiP',
                'DarkQuest',
                'DUNE',
                'SHADOWS'
                ]
channels_decay = [  '2Gamma',
                    '2El',
                    '2Mu',
                    '3Pi0',
                    '3Pi',
                    '2PiGamma',
                    '2Pi0Eta',
                    '2PiEta',
#                    '2Pi0EtaPrim',
#                    '2PiEtaPrim'
                    ]

channels_production = [ 'primakoff',
                        'photonfrommeson',
                        'mixingPi0',
                        'mixingEta',
                        'mixingEtaPrim',
                        'BmesonK',
                        'BmesonKstar',
                        'DmesonPi',
                        'KSmesonPi0'
                        ]
reference_couplings = {'primakoff':         1e-4,
                       'photonfrommeson':   1e-4,
                       'mixingPi0':         1e-4,
                       'mixingEta':         1e-4,
                       'mixingEtaPrim':     1e-4,
                       'BmesonK':           1e-10,
                       'BmesonKstar':       1e-10,
                       'DmesonPi':          1e-10,
                       'KSmesonPi0':        1e-10}
scaling_exponent =    {'primakoff':         2,
                       'photonfrommeson':   2,
                       'mixingPi0':         2,
                       'mixingEta':         2,
                       'mixingEtaPrim':     2,
                       'BmesonK':           1,
                       'BmesonKstar':       1,
                       'DmesonPi':          1,
                       'KSmesonPi0':        1,}

from ALP_rescale.general.setup import*

couplings = [   'CBB',
                'CWW',
                'CGG',
                'Cll',
                'Cqq'
            ]

channels_decay = ['2Gamma',
                  '2El',
                  '2Mu',
                  '3Pi0',
                  '3Pi',
                  '2PiGamma',
                  '2Pi0Eta',
                  '2PiEta',
                  '2Pi0EtaPrim',
                  '2PiEtaPrim',
                  '2KPi0',
                  '2Pi2Pi0',
                  '4Pi',
                  '2Omega',
                  '2Kstar',
                  '2Phi'
                    ]

channels_production = [ 'Primakoff',
                        'PhotonFromMeson',
                        'MixingPi0',
                        'MixingEta',
                        'MixingEtaPrim',
                        'Bmeson',
                        'Dmeson',
                        'recast',
                        ]
reference_couplings = {'Primakoff':         1,
                       'PhotonFromMeson':   1,
                       'MixingPi0':         1,
                       'MixingEta':         1,
                       'MixingEtaPrim':     1,
                       'Bmeson':            1,
                       'Dmeson':            1,
                       'recast':            1,
                       }
scaling_exponent =    {'Primakoff':         2,
                       'PhotonFromMeson':   2,
                       'MixingPi0':         2,
                       'MixingEta':         2,
                       'MixingEtaPrim':     2,
                       'Bmeson':            2,
                       'Dmeson':            2,
                       'recast':            2,
                       }
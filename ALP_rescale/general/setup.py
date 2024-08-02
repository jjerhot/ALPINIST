variables = [   'mX',
                'GammaX',
                'tauX',
                'BRdecay',
                'BRprod',
            ]

experiments = { 'NA62':         'NA62',
                'HIKE':         'NA62',
                'CHARM':        'CHARM',
                'NuCal':        'NuCal',
                'SHiP':         'SHiP',
                'SHiPecn4':     'SHiP',
                'DarkQuest':    'DarkQuest',
                'DarkQuestPhase2':    'DarkQuest',
                'DUNE':         'DUNE',
                'SHADOWS':      'SHADOWS',
                'KOTOpnn':      'KOTO',
                'KOTOexclPnn':  'KOTO',
                'KOTOdump':     'KOTO',
                'KOTO2pnn':     'KOTO2',
                'KOTO2dump':    'KOTO2',
                'E137':         'E137',
                'E141':         'E141',
                'KLEVER':       'KLEVER',
                'KLEVERext':    'KLEVER',
                'NuTeV':        'NuTeV',
                'BEBC':         'BEBC',
                'BEBCcuboid':   'BEBC',
                'ORCA':         'ORCA'}

PoT = {
    'NA62':         1E18, 
    'CHARM':        2.4E18, #
    'NuCal':        1.7E18,
    'SHiP':         6E20,
    'DarkQuest':    1.44E18,
    'DarkQuestPhase2': 1e20,
    'DUNE':         1.1E22, 
    'SHADOWS':      5E19, 
    'KOTOdump':     2.2E18, #x10 statistics
    'KOTO2dump':    6.*1E21, #x10 statistics
    'BEBC':         2.72E18, #
    'NuTeV':        2.54E18, #
    'ORCA':         8E19
}

n_events_90CL = {
    'alp': {
    'NA62':         2.3, # 1E18
    'CHARM':        2.3,
    'NuCal':        3.6,
    'SHiPecn4':     3*2.3, #2E20
    'SHiP':         2.3, #6E20
    'DarkQuest':    10.,
    'DarkQuestPhase2': 10.,
    'DUNE':         0.23, #x10 statistics
    'SHADOWS':      0.46, #:5E19
    'KOTOpnn':      2.3, #already scaled
    'KOTOexclPnn':  4.2,
    'KOTOdump':     0.23, #x10 statistics
    'KOTO2pnn':     3.14,
    'KOTO2dump':    0.23, #x10 statistics
    'E137':         2.3,
    'E141':         2.3,
    'NuTeV':        2.3,
    'BEBC':         2.3,
    'BEBCcuboid':   2.3,
    'ORCA':         2.3,
    },
    'hnl': {
    'NA62':         1.,# 1E18
    'HIKE':         2.3/15,
    'CHARM':        2.3,#
    'NuCal':        3.6,
    'SHiPecn4':     3*2.3, #2E20
    'SHiP':         2.3/3, #6E20
    'DarkQuest':    10.,#
    'DarkQuestPhase2': 10.,
    'DUNE':         1., # 10 years data taking with 10 events required (1912.07622)
    'SHADOWS':      0.46, #:5E19
    'KOTOpnn':      2.3, #already scaled
    'KOTOexclPnn':  4.2,
    'KOTOdump':     0.23, #x10 statistics
    'KOTO2pnn':     3.14,
    'KOTO2dump':    0.23, #x10 statistics
    'E137':         2.3,
    'E141':         2.3,
    'NuTeV':        2.30259,# #0 obs 0.56+0.005+0.002 exp,
    'BEBC':         3.45, # 1 obs 0.6 exp,
    'BEBCcuboid':   2.3,
    'ORCA':         2.3,
    },
    'dp': {
    'NA62':         2.3,#1E18
    'CHARM':        2.3,
    'NuCal':        3.6,
    'SHiPecn4':     3*2.3, #2E20
    'SHiP':         2.3, #6E20
    'DarkQuest':    10.,
    'DarkQuestPhase2': 10.,
    'DUNE':         0.23, #x10 statistics
    'SHADOWS':      0.46, #:5E19
    'KOTOpnn':      2.3, #already scaled
    'KOTOexclPnn':  4.2,
    'KOTOdump':     0.23, #x10 statistics
    'KOTO2pnn':     3.14,
    'KOTO2dump':    0.23, #x10 statistics
    'E137':         2.3,
    'E141':         2.3,
    'NuTeV':        2.3,
    'BEBC':         2.3,
    'BEBCcuboid':   2.3,
    'ORCA':         2.3 * 8/600 ,
    },
    'ds': {
    'NA62':         2.3,# 1E18
    'CHARM':        2.3,
    'NuCal':        3.6,
    'SHiPecn4':     3*2.3, #2E20
    'SHiP':         2.3, #6E20
    'DarkQuest':    10.,
    'DarkQuestPhase2': 10.,
    'DUNE':         0.23, #x10 statistics
    'SHADOWS':      0.46, #:5E19
    'KOTOpnn':      2.3, #already scaled
    'KOTOexclPnn':  4.2,
    'KOTOdump':     0.23, #x10 statistics
    'KOTO2pnn':     3.14,
    'KOTO2dump':    0.23, #x10 statistics
    'E137':         2.3,
    'E141':         2.3,
    'NuTeV':        2.3,
    'BEBC':         2.3,
    'BEBCcuboid':   2.3,
    'ORCA':         2.3,
    }
}

#.. to be seen if can be removed entirely 
# reference_couplings = {'primakoff':         1e-4,
#                        'photonfrommeson':   1e-4,
#                        'mixingPi0':         1e-4,
#                        'mixingEta':         1e-4,
#                        'mixingEtaPrim':     1e-4,
#                        'Bmeson':            1e-10,
#                        'Dmeson':            1e-10,
#                        'recast':            1,}

# scaling_exponent =    {'primakoff':         2,
#                        'photonfrommeson':   2,
#                        'mixingPi0':         2,
#                        'mixingEta':         2,
#                        'mixingEtaPrim':     2,
#                        'Bmeson':            1,
#                        'Dmeson':            1,
#                        'recast':            2,}

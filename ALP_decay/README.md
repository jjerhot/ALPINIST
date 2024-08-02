ALP_decay module for ALPINIST
=============================

Derived from Pi0MCfromALPSXs.C by Babette Dobrich & Tommaso Spadaro

ALP_decay loads tables from tab_prod/ for given experiment and production mode and simulates chosen decay. Output is written in tab_decay/. For 3-body decays Dalitz plots from widths/Dalitz/ are used.

Table of Contents
-----------------

- [ALP\_decay module for ALPINIST](#alp_decay-module-for-alpinist)
  - [Table of Contents](#table-of-contents)
  - [Description](#description)
  - [Requirements](#requirements)
  - [Usage](#usage)
  - [License](#license)


Description
-----------
ALP_decay loads tables from tab_prod/ for given experiment and production mode and simulates chosen decay. Output table is written in tab_decay/ in format: mass width Nevents. For 3-body decays Dalitz plots from widths/Dalitz/ are used.
The main loops are:
  * Loops over 2 input mass-energy tables (0.01MeV-10MeV-5310MeV), referential coupling value the same as for input tables
  * Loops over masses (201 values): 0.0001MeV-5310MeV
  * Loops over decay widths (101 values): 1E-25 - 1E-10
  * (for HNLs: loops over lepton coupling dominance scenarios for production)

ALP_decay contains the following files:
  * DecayMC.C: Contains main. Parses arguments.
  * DecayMCGlobal.h: Constants and other global parameters.
  * ExpParameters: Contains geometry of experiments and selection cuts.
  * DecayMCProcess.C: Loads input tables, simulates decay, checks if event passes selection, writes output.

Requirements
------------
ALP_decay requires the following to run:

  * [ROOT][ROOT]
  * cmake3 (to compile)


Usage
-----
ALP_decay can be used with both ROOT interpreter or can be compiled.

To compile, execute in ALP_decay/

```sh
cmake3 -B build && cmake3 --build build
```

To run, execute ./DecayMC in ALP_decay/bin/ with -h to show help or with required parameters (-e EXP -p PRODUCTION_MODE -d DECAY_MODE -n NUMBER_OF_EVENTS_PER_BIN). Example for simulating a->2gamma decay in CHARM experiment for primakoff production with 1M events per mass and weight bin:

```sh
./DecayMC -x alp -e CHARM -p primakoff -d 2Gamma -n 1000000
```

Alternatively can be run with ROOT interpreter. Execute root in ALP_decay/src/ and then for running the same as in example above execute:
```sh
.include ../include
.L ExpParameters.C
.x DecayMCProcess.C(0,1,0,0,1000000)
```

The parameter mapping for DecayMCProcess.C is the following:
  * Exotic particles: 0-Axion-Like-Particles, 1-Heavy Neutral Leptons, 2-Dark Photons,3-Dark Scalars
  * Experiment: 0-CHARM, 1-BEBC, 2-NuCal, 3-NuTeV, 4-NA62, 5-DarkQuest, 6-DarkQuestPhase2, 7-DUNE, 8-SHiP, 9-KOTOdump, 10-KOTOpnn, 11-KOTO2dump, 12-KOTO2pnn, 13-KOTOexclPnn, 14-KLEVER, 15-KLEVERext, 16-SHADOWS, 17-SHiPecn4, 18-ORCA, 19-BEBCcuboid
  * Production mode: 
    * for ALP:  0-Bmeson, 1-Bmeson, 2-Dmeson, 3-primakoff, 4-photonfrommeson, 5-mixingPi0, 6-mixingEta, 7-mixingEtaPrim
    * for HNL:  0-Bmeson, 1-Dmeson
    * for DP:   0-Brems, 1-MesonDecay, 2-Mixing
    * for DS:   0-Bmeson, 1-Bmeson2S
  * Decay mode: 
    * for ALP: 0-2Gamma, 1-2El, 2-2Mu, 3-3Pi0, 4-3Pi, 5-2PiGamma, 6-2Pi0Eta, 7-2PiEta, 8-2Pi0EtaPrim, 9-2PiEtaPrim, 10-2Pi, 11-2K, 12-2KPi,13-2Pi2Pi0
    * for HNL: 1-PiEl, 2-PiMu, 3-RhoEl, 4-RhoMu, 5-NuElEl, 6-NuMuMu, 7-NuElMu, 8-PiNu, EtaNu, RhoNu
    * for DP: 1-2El, 2-2Mu, 3-2Pi, 4-3Pi, 5-4Pi, 6-2Pi2Pi0, 7-2K, 8-2KPi
    * for DS: 1-2El, 2-2Mu, 3-2Pi, 4-2K, 5-4Pi, 6-2Pi2Pi0
  * Number of MC events


License
-------

BSD 3-Clause License

Copyright (c) 2022, Jan Jerhot, Babette Dobrich, Fatih Ertas, Felix Kahlhoefer, Tommaso Spadaro
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

[ROOT]: https://root.cern.ch/
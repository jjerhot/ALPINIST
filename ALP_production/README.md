ALP_production module for ALPINIST
==================================

ALP_production simulates the yield of exotic particles produced in a fixed target by a particle beam of given energy. The output is in a form of a table of exotic energy-mass-angle-yield for some fixed value of a model-dependent parameter (BR or coupling) with the angle-energy range being experiment-specific. The output tables are stored in tab_prod/ directory.

ALP_production uses either Pythia for SM meson production or external tables of mesons can be used (located in tab_meson/ directory and given subdirectories).


Table of Contents
-----------------

  * [Description](#description)
  * [Requirements](#requirements)
  * [Usage](#usage)
  * [License](#license)

ALP_production contains the following parts:
  * exotic_production.py - Main module;
  * general/ - contains methods common to all production modes;
  * modes/ - exotic production in specific processes.


Requirements
------------

ALP_rescale requires the following to run:

  * [python][python] 3.1+
  * python packages: numpy, sklearn, multiprocessing, tqdm, argparse
  * [pythia][pythia] 8.2+ 


Usage
-----

To run with default settings, run the exotic_production.py with arguments --exp to select the experiment and --prod to select the production mode. Optional arguments are:
 * -n: number of threads (default 1);
 * -ext/-no-ext: external source or Pythia (default external dataset in tab_mesons/);
 * -nprod: number of generated mesons in Pythia or external source (default 1E5);
 * -ndec: number of decays simulated for each meson produced (default 10).
For listing all options run with --help:

```sh
./exotic_production.py --help
```

License
-------

BSD 3-Clause License

Copyright (c) 2023, Jan Jerhot, Babette Dobrich, Fatih Ertas, Felix Kahlhoefer, Tommaso Spadaro
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

[python]: https://www.python.org/
[pythia]: https://pythia.org/
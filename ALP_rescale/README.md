ALP_rescale module for ALPINIST
===============================

ALP_rescale loads the tables from /tab_decay for the experiment, production and decay mode selected from the lists stored in alp_setup.py and calculates the number of predicted events for the following model-dependent parameters:
C<sub>BB</sub>, C<sub>WW</sub>, C<sub>gg</sub>, C<sub>&#8467;&#8467;</sub>, &#923;, A, B
(all relations used for the calculations are described in appendices of J. Jerhot et al. [arXiv:2201.xxxxx][2201.xxxxx] and references therein).
Output tables in format [m_ALP g_ALP nEvents] are stored in /tab_toPlot for each experiment.

For general Yukawa coupling C<sub>&#8467;&#8467;</sub>=C<sub>qq</sub> a separate module alp_2mu_rescale.py can be used for rescaling for ALP &#8594; 2&#956; decay (using relations from B. Dobrich et al., [arXiv:1810.11336][1810.11336]).

Table of Contents
-----------------

  * [Description](#description)
  * [Requirements](#requirements)
  * [Usage](#usage)
  * [License](#license)

Description
-----------

ALP_rescale contains the following files:
  * alp_rescale.py: Main module
  * load_data.py: Loading and processing of data  
  * effective_coupling.py: Calculation of effective couplings
  * decay_widths.py: Calculation and interpolation of loaded decay widths
  * alp_setup.py: Lists and dictionaries of experiments, production and decay modes covered by ALP_rescale
  * alp_functions.py: Contains several universal functions
  * alp_constants.py: Contains several universal constants
  * alp_2mu_rescale.py: Separate main module for Yukawa coupling

Requirements
------------

ALP_rescale requires the following to run:

  * [python][python] 3.1+


Usage
-----

To run with default settings, run the alp_rescale.py with your installed version of [python][python], e.g.:

```sh
python3 alp_rescale.py
```

The default settings are:
C<sub>BB</sub> = 1, C<sub>WW</sub> = 1, C<sub>gg</sub> = 1, C<sub>&#8467;&#8467;</sub> = 0, &#923; = 1000, A = 0, B = 0 and running over all experiments available and summing over all production and decay modes.

All these parameters can be modified, e.g.:

```sh
python3 alp_rescale.py -e NA62 --prod BmesonK --decay 2Gamma -cbb 0 -cww 0 -cgg 1 -l 10000 -a 3 -b -3 
```

To list all options run:

```sh
python3 alp_rescale.py -h
```

Similarly, alp_2mu_rescale.py can be run as
```sh
python3 alp_2mu_rescale.py
```
and offers several options, listed when run with `-h` or `--help` argument. 

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


[2201.xxxxx]: https://arxiv.org/
[1810.11336]: https://arxiv.org/abs/1810.11336
[python]: https://www.python.org/
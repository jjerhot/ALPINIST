ALP_rescale module for ALPINIST
===============================

ALP_rescale loads the tables from /tab_decay for the experiment, production and decay mode selected from the lists stored in alp_setup.py (scalar_setup.py for Dark Scalar) and calculates the number of predicted events for the following model-independent parameters: m<sub>X</sub>, &#915;<sub>X</sub>, &#964;<sub>X</sub>, BR<sub>prod</sub>, BR<sub>decay</sub> and model-dependent parameters for ALPs:
C<sub>BB</sub>, C<sub>WW</sub>, C<sub>gg</sub>, C<sub>&#8467;&#8467;</sub>, &#923;, A, B or Y for Dark Scalar
(all relations used for the calculations are described in appendices of J. Jerhot et al. [arXiv:2201.05170][2201.05170] and references therein).
Output tables in format [m_ALP g_ALP nEvents] are stored in /tab_toPlot for each experiment.

For general ALP-Yukawa coupling C<sub>&#8467;&#8467;</sub>=C<sub>qq</sub> a separate module alp_2mu_rescale.py can be used for rescaling for ALP &#8594; 2&#956; decay (using relations from B. Dobrich et al., [arXiv:1810.11336][1810.11336]).

Table of Contents
-----------------

  * [Description](#description)
  * [Requirements](#requirements)
  * [Usage](#usage)
  * [License](#license)

Description
-----------

ALP_rescale contains the following parts:
  * alp_rescale.py: Main module
  * load_data.py: Loading and processing of data  
  * general/ - contains methods common to all modes and exotic particles: 
    * alp_constants.py
    * alp_functions.py
    * alp_mergeSigRegions.py: Allows merging results from several signal regions with given weights
  * alp/ - contains ALP-specific modules:
    * effective_coupling.py: Calculation of effective couplings
    * decay_widths.py: Calculation and interpolation of loaded decay widths
    * alp_setup.py: Lists and dictionaries of experiments, production and decay modes covered by ALP_rescale
  * scalar/ - contains Dark Scalar-specific modules:
    * decay_widths_DS.py: Calculation and interpolation of loaded decay widths
    * scalar_setup.py: Lists and dictionaries of experiments, production and decay modes covered by ALP_rescale
    
  * alp_2mu_rescale.py: Separate main module for Yukawa coupling
  * plot_exclusion.py: Separate module for plotting and extracting contours in python

Requirements
------------

ALP_rescale requires the following to run:

  * [python][python] 3.1+
  * python packages: numpy, scipy, mpmath, argparse


Usage
-----

To run with default settings, run the alp_rescale.py with your installed version of [python][python] with arguments -varX and -varY for selected X- and Y-axis variables, e.g. for X-axis mX (exotic mass) and Y-axis CBB (ALP C<sub>BB</sub> coupling) run:

```sh
./alp_rescale.py -varX mX -varY CBB
```

The default settings are:
&#923; = 1000, A = 0, B = 0 and running over all experiments available and summing over all production and decay modes. The setup of couplings is specified later on but can also be skipped by running with -only option for the coupling selected for varY (other couplings are then set to a fixed 0) as e.g.:
```sh
./alp_rescale.py -varX mX -varY CBB-only
```

All the parameters can be modified and multiple experiments, decay and production modes can be ran simultaneously, e.g.:

```sh
./alp_rescale.py -varX mX -varY CWW -e NA62 CHARM --prod BmesonK BmesonKstar --decay 2Gamma 2El --lambda 10000 -a 3 -b -3 
```

Besides model-dependent parameters, also model-independent parameters can be used, such as decay width gammaX, lifetime tauX, production branching ratio BRprod or decay branching ration BRdecay.

To list all options and parameters run:

```sh
./alp_rescale.py -h
```

Similarly, alp_2mu_rescale.py can be run as
```sh
./alp_2mu_rescale.py
```
and offers several options, listed when run with `-h` or `--help` argument. 

To use plot_exclusion.py, run it with selected experiment or choose specific production or decay modes
```sh
./plot_exclusion.py -e CHARM
```

To list all options and parameters run:

```sh
./plot_exclusion.py -h
```

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


[2201.05170]: https://arxiv.org/abs/2201.05170
[1810.11336]: https://arxiv.org/abs/1810.11336
[python]: https://www.python.org/
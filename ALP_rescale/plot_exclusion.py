#!/usr/bin/env python3

 ##
 # @file plot_exclusion.py
 # @brief Reads ALP_rescale output and generates contours
 # @details Provides terminal interface
 # 
 # Allows either just an exclusion contour display or contour or pdf export

import sys
# sys.path.append("../")
import os

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import argparse
import sys

from ALP_rescale.general import setup
from ALP_rescale.alp import alp_setup as setalp
from ALP_rescale.scalar import scalar_setup as setds
from ALP_rescale.vector import vector_setup as setdp
from ALP_rescale.hnl import hnl_setup as sethnl


def main(argv=None):
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    parser = argparse.ArgumentParser(description='ALPINIST plot or extract exclusion contours. \n Select the experiment, production and decay modes and other parameters or use default')
    parser.add_argument("-e","--exp",  default="", type=str, nargs='*', help="Experiments available (case sensitive): exp = "+" | ".join(setup.experiments.keys())+". If not specified, running over all experiments available.")
    parser.add_argument("-x","--exotic",  default="alp", type=str, help="Exotic particles available: exotic = alp | hnl | ds | dp . Default is alp")
    parser.add_argument("-comb","--combine", dest='comb', action='store_true', help="Combine experiments")
    parser.add_argument("-no-comb","--no-combine", dest='comb', action='store_false', help="Rescale experiments separately. Default")
    parser.add_argument("-ternary","--ternary-evaluation", dest='ter', action='store_true', help="Plot U2 ternarily evaluated bounds at 90CL at a fixed mass (only for hnl).")
    parser.set_defaults(comb=False)
    args = parser.parse_args()

    plot_exclusion(args.exotic, args.exp, args.comb, args.ter)

def plot_exclusion(exotic_, exp_, comb_=False, ternary_=False):

    if ternary_: # if ternary eavluation call auxillary module 
        if exotic_ != "hnl":
            print("[Error:] \t Ternary evaluation only available for exotic class HNL. Exiting.")
            sys.exit(1)
        from ALP_rescale.hnl import plot_ternary_exclusion
        plot_ternary_exclusion.main(exp_, comb_)
        sys.exit(0)

    #experiment:
    experiments = []
    if comb_:
        experiments = ["combined"]
        print("[Info:] \t Checking combined datasets")
    else:
        if exp_ == "":
            experiments = setup.experiments.keys()
            print("[Info:] \t Selected all experiments available")
        else:
            for exp in exp_:
                if exp in setup.experiments.keys():
                    experiments.append(exp)
                else:
                    sys.exit("[Error:] \t Experiment " + exp + " not available. Experiment modes available: exp = " + ' | '.join(setup.experiments.keys()) + ". If not specified, running over all experiments available.")
            print("[Info:] \t Selected experiments:", ', '.join(experiments))

    exotic = exotic_
    if not exotic in ['alp' , 'hnl', 'ds', 'dp']:
        sys.exit("[Error:] exotic: "+exotic+" not available.")

    main_loop(experiments,exotic)

    return

def extract_contour(plot):
    contours = []
    mpl_version = tuple(map(int, (sys.modules[plt.__package__].__version__.split("."))))
    if mpl_version < (3,8): # approach deprecated with 3.8
        for col in plot.collections:
            paths = plot.col.get_paths()
            if len(paths) == 0: continue
            for path in paths:
                # if paths[path].vertices[0][0] == paths[path].vertices[-1][0]: continue # skip trivial cases
                contours.append(path.vertices)
    else:
        for path in plot.get_paths():
            contours.append(path.vertices)
    return contours

def main_loop(experiments = [], exotic = "alp"):

    # find files with common x,y-axes and for the same model
    files_common = []
    filenames_common = set()
    n_exp = 0
    for exp in experiments:
        if exp == "combined":
            path = Path(os.path.dirname(__file__)+'/../tab_toPlot/combined') / exotic
        else:
            path = Path(os.path.dirname(__file__)+'/../tab_toPlot') / setup.experiments[exp] / exotic
        if not path.exists() or not path.is_dir():
            raise FileNotFoundError("[Error:] \t directory corresponding to experiment "+exp+" and exotic "+exotic+" not found")
        files_exp = []
        for file in path.iterdir():
            if exp == "combined":
                files_exp.append(input_file(file.name, exotic))
            else:
                if file.name.split("_")[0] == exp:
                    infile = input_file(file.name, exotic)
                    if infile.label == exp:
                        files_exp.append(infile)
        if n_exp == 0:
            files_common = files_exp
            for file in files_exp:
                filenames_common.add(file.name_for_comparing)
        elif len(filenames_common)>0:
            to_compare = [file.name_for_comparing for file in files_exp]
            intersect = set(to_compare).intersection(filenames_common)
            if len(intersect) > 0:
                for i in range(len(files_common)-1,-1,-1):
                    if files_common[i].name_for_comparing not in intersect:
                        files_common.pop(i)
                filenames_common = intersect
                for new_file in files_exp:
                    if new_file.name_for_comparing in intersect:
                        files_common.append(new_file) 
            else:
                files_common = []
                filenames_common = set()
        n_exp+=1

    if len(filenames_common) == 0:
        print("[Info:] \t No files common to",', '.join(experiments),"found, generate them first. Closing.")
        exit()

    #main loop
    while True:

        print("- Options:")
        print("[1] show pdf")
        print("[2] export pdf")
        print("[3] export 90%CL contour")
        print("[4] exit")
        mode_chosen = int(input(" - Enter number 1 to 4 corresponding to one of the options above:"))
        if not 0 < mode_chosen < 4:
            print("[Info:] Exitting")
            exit(0)

        if mode_chosen == 3:
            #select experiment
            exp_selected = ""
            if len(experiments) == 1:
                exp_selected = experiments[0]
            else:
                print("- Experiments available:")
                n_files = 0
                for exp in experiments:
                    n_files+=1
                    print("["+str(n_files)+"] "+exp)
                exp_chosen = int(input(" - Enter number 1 to "+str(len(experiments))+" corresponding to one of the experiments above:"))
                if not 1<= exp_chosen <= len(experiments):
                    raise ValueError("[Error:] \t Invalid answer")
                exp_selected = experiments[exp_chosen-1]
            if experiments == ["combined"]:
                path = Path(os.path.dirname(__file__)+'/../tab_toPlot/combined') / exotic
            else:
                path = Path(os.path.dirname(__file__)+'/../tab_toPlot') / setup.experiments[exp_selected] / exotic
            if not path.exists() or not path.is_dir():
                raise FileNotFoundError("[Error:] \t directory corresponding to",exp_selected,"not found")
            files_exp = []
            for file in path.iterdir():
                if file.name.split("_")[0] == exp_selected:
                    infile = input_file(file.name, exotic)
                    if infile.label == exp_selected:
                        files_exp.append(infile)
                elif exp_selected == "combined":
                    files_exp.append(input_file(file.name, exotic))
            if len(files_exp) == 0:
                print("[Info:] \t No files for",exp_selected,"found, generate them first. Closing.")
                exit()

            #select file
            print("- Files for",exp_selected,"available:")
            n_files = 0
            for file in files_exp:
                n_files+=1
                print("["+str(n_files)+"] "+str(file))
            n_chosen = int(input(" - Enter number 1 to "+str(len(files_exp))+" corresponding to one of the files above:"))
            if not 1<= n_chosen <= len(files_exp):
                raise ValueError("[Error:] \t Invalid answer")
            file_chosen = list(files_exp)[n_chosen-1]
            print("[Info:] \t Chosen file:",str(file_chosen))
            
            #generate contour
            plot = plot_setup(exotic, file_chosen)
            if file_chosen.experiment == "combined":
                style = exp_style[file_chosen.label.split("-")[0]]
                data = np.loadtxt(os.path.dirname(__file__)+'/../tab_toPlot/combined/'+exotic+'/'+str(file_chosen))
                nevts = setup.n_events_90CL[exotic][file.label.split("-")[0]]
            else:
                style = exp_style[file_chosen.label]
                data = np.loadtxt(os.path.dirname(__file__)+'/../tab_toPlot/'+file_chosen.experiment+'/'+exotic+'/'+str(file_chosen))
                nevts = setup.n_events_90CL[exotic][file.label]
            contour = plot.ax.contour(plot.xaxis_list, plot.yaxis_list,np.transpose(data[:,2].reshape((plot.xaxis_steps, plot.yaxis_steps))),[nevts],colors = [(style[2],style[3],style[4])],linewidths=style[1],linestyles=line_styles[0])

            contourNevts = str("{:.2f}".format(nevts))
            #extract and export contour
            outPath = os.path.dirname(__file__)+'/../Figures/contours/'+file_chosen.experiment+'/'+exotic+'/'
            if not os.path.isdir(outPath): os.makedirs(outPath)
            np.savetxt(outPath+'contour-'+contourNevts+'_'+str(file_chosen),np.concatenate(extract_contour(contour)),fmt='%.4e')
            print("[Info:] file contour-"+contourNevts+'_'+str(file_chosen)+' saved to Figures/contours/'+file_chosen.experiment+'/'+exotic+' directory')

            continue

        print("- Files common to",', '.join(experiments),"available:")
        n_files = 0
        for file in filenames_common:
            n_files+=1
            print("["+str(n_files)+"] "+file)
        n_chosen = int(input(" - Enter number 1 to "+str(len(filenames_common))+" corresponding to one of the files above:"))
        if not 1<= n_chosen <= len(filenames_common):
            raise ValueError("[Error:] \t Invalid answer")
        file_chosen = list(filenames_common)[n_chosen-1]
        print("[Info:] \t Chosen file:",file_chosen)

        plt.rcParams['text.usetex'] = True
        common_setup = None

        #categorize files per experiment
        files_exp = {}
        for file in files_common:
            if file.label in files_exp.keys():
                files_exp[file.label].append(file)
            else:
                files_exp[file.label] = [file]

        leg_proxy = []
        leg_label = []
        for exp in files_exp.keys():
            #sort first by number of modes:
            sorted(files_exp[exp], key=lambda file: (len(file.modes_prod)+len(file.modes_prod)))
            n_matched_files = 0
            # for ifile in range(len(files_exp[exp])):
            for file in files_exp[exp]:
                if file.name_for_comparing != file_chosen:
                    continue
                if n_matched_files>4:
                    print("[Warning:] \t Only up to 4 modes per experiment can be plotted. Skipping ",str(file))
                    continue
                if common_setup == None:
                    common_setup = plot_setup(exotic, file)
                if file.experiment == "combined":
                    style = exp_style[file.label.split("-")[0]] #take first experiment as reference
                    data = np.loadtxt(os.path.dirname(__file__)+'/../tab_toPlot/combined/'+exotic+'/'+str(file))
                    nevts = setup.n_events_90CL[exotic][file.label.split("-")[0]]
                else:
                    style = exp_style[file.label]
                    data = np.loadtxt(os.path.dirname(__file__)+'/../tab_toPlot/'+file.experiment+'/'+exotic+'/'+str(file))
                    nevts = setup.n_events_90CL[exotic][file.label]
                if len(data) != len(common_setup):
                    print("[Warning:] \t File",str(file),"has different size than the referential file, skipping")
                    continue
                if style[0] == True and n_matched_files == 0:
                    common_setup.ax.contourf(common_setup.xaxis_list, common_setup.yaxis_list,np.transpose(data[:,2].reshape((common_setup.xaxis_steps, common_setup.yaxis_steps))),[0.,nevts,np.inf],colors = [(style[2],style[3],style[4],0.),(style[2],style[3],style[4],style[5])])
                    leg_proxy.append(plt.Rectangle((0,0),1,1,fc = (style[2],style[3],style[4],style[5])))
                else:
                    print('shape:',np.shape(np.transpose(data[:,2].reshape((common_setup.xaxis_steps, common_setup.yaxis_steps)))))
                    cont = common_setup.ax.contour(common_setup.xaxis_list, common_setup.yaxis_list,np.transpose(data[:,2].reshape((common_setup.xaxis_steps, common_setup.yaxis_steps))),[nevts],colors = [(style[2],style[3],style[4])],linewidths=style[1],linestyles=line_styles[n_matched_files])
                    
                    # h1,l1 = cont.legend_elements()
                    leg_proxy.append(cont.legend_elements()[0][0])
                n_matched_files+=1
                if len(experiments) == 1:
                    leg_label.append(file.name_for_comparing)
                else:
                    leg_label.append(file.experiment)
        # plt.legend(leg_proxy, [r'$\mathrm{CHARM}\, (\mu\mu,\,m_X = 1\,\mathrm{GeV})$',r'$\mathrm{CHARM}\, (ee,\,m_X = 1\,\mathrm{GeV})$'], loc='center left', bbox_to_anchor=(0.05, 0.1), frameon=False)
        plt.legend(leg_proxy, leg_label, loc='center left', bbox_to_anchor=(0.05, 0.1), frameon=False)


        if mode_chosen == 1:
            plt.show()
            continue

        elif mode_chosen == 2:
            contourNevts = str("{:.2f}".format(nevts))
            outPath = os.path.dirname(__file__)+'/../Figures/'+exotic+'/'
            if not os.path.isdir(outPath): os.makedirs(outPath)
            plt.savefig(outPath+file_chosen+'_contour-'+contourNevts+'.pdf')
            print("[Info:] file "+file_chosen+'_contour-'+contourNevts+'.pdf saved to Figures/'+exotic+'/ directory')
            continue

        else:
            print("[Info:] Exitting")
            exit(0)
    return

class setexo:
    def __init__(self, exo):
        ref = {
            "alp": setalp,
            "hnl": sethnl,
            "ds":   setds,
            "dp":   setdp
        }
        self.experiments = {}
        self.channels_production = []
        self.channels_decay = []
        self.variables = []
        self.couplings = []
        if exo in ref.keys():
            self.experiments = ref[exo].experiments
            self.channels_production = ref[exo].channels_production
            self.channels_decay = ref[exo].channels_decay
            self.variables = ref[exo].variables
            self.couplings = ref[exo].couplings
        if exo == "hnl": self.couplings.append("U2")

class input_file:
    '''
    Input file parameters
    '''
    def __init__(self, filename = "", exo = "alp"):
        s = setexo(exo)
        self._filename = filename.split("/")[-1]
        self._path = filename.split(self._filename)[0]
        self.name_for_comparing = self._filename.replace(".dat","")
        self.experiment = ""
        self.label = ""
        self.axes = []
        self.fixed_values = {}
        self.scale_with = ""
        self.scale_values = {}
        self.modes_prod = []
        self.modes_decay = []
                #  experiment = "", label = "", axes = [], fixed_values = [], scale_with = "", scale_values = [], modes_prod = [], modes_decay = []):
        self.label = self.name_for_comparing.split("_")[0]
        if self.label in s.experiments.keys():
            # self.label = self.name_for_comparing.split("_")[0]
            self.experiment = s.experiments[self.label]
        else:
            self.experiment = "combined"
        self.name_for_comparing = self.name_for_comparing.replace(self.label+"_","")
        for mode in reversed(s.channels_production): #check first presence of BmesonKstar
            if mode in self.name_for_comparing:
                self.name_for_comparing = self.name_for_comparing.replace(mode,"")
                self.modes_prod.append(mode)
        for mode in s.channels_decay:
            if mode in self.name_for_comparing:
                self.name_for_comparing = self.name_for_comparing.replace(mode,"")
                self.modes_decay.append(mode)
        # if len(self.modes_prod) == 0:
        #     self.modes_prod.append("")
        # if len(self.modes_decay) == 0:
        #     self.modes_decay.append("")
        if len(self.name_for_comparing.split("_"))>1:
            varX = self.name_for_comparing.split("_")[0]
            if varX in (s.variables + s.couplings):
                self.axes.append(varX)
                # self.name_for_comparing = self.name_for_comparing.replace(varX+"_","")
            else:
                self.label = ""
                self.experiment = ""
            varY = self.name_for_comparing.split("_")[1]
            if varY in (s.variables + s.couplings):
                self.axes.append(varY)
                # self.name_for_comparing = self.name_for_comparing.replace(varY+"_","")
            else:
                self.label = ""
                self.experiment = ""
            for sub in self.name_for_comparing.split("_"):
                if "scaleWith" in sub:
                    for subsub in sub.split("-"):
                        if "scaleWith" in subsub:
                            self.scale_with = subsub.replace("scaleWith","")
                            continue
                        for var in (s.variables + s.couplings):
                            if var in subsub:
                                self.scale_values[var] = float(subsub.replace(var,""))
                                break
                    break
            for sub in self.name_for_comparing.split("_"):
                if "fixed" in sub:
                    for subsub in sub.split("-"):
                        for var in (s.variables + s.couplings):
                            if var in subsub:
                                self.fixed_values[var] = subsub.replace(var,"")
                                break
                    break
            while self.name_for_comparing[-1] == "-" or self.name_for_comparing[-1] == "_":
                self.name_for_comparing = self.name_for_comparing[:-1]
        else:
            self.label = ""
            self.experiment = ""
    
    def __str__(self):
        return self._filename

class plot_setup:
    '''
    Setup of the exclusion plot
    '''
    def __init__(self, exotic, file : input_file):
        if str(file) == "":
            print("[Error:] \t Input file object empty")
            exit()

        if file.experiment in setup.experiments:
            self._data_setup = np.loadtxt(os.path.dirname(__file__)+'/../tab_toPlot/'+file.experiment+'/'+exotic+'/'+str(file))
        else:
            self._data_setup = np.loadtxt(os.path.dirname(__file__)+'/../tab_toPlot/combined/'+exotic+'/'+str(file))

        self.yaxis_steps = self._yaxis_log_steps(self._data_setup)
        self.xaxis_steps = round(len(self._data_setup)/self.yaxis_steps)

        self.xaxis_list = np.array([self._data_setup[i*self.yaxis_steps,0] for i in range(self.xaxis_steps)])
        self.yaxis_list = np.array([self._data_setup[i,1] for i in range(self.yaxis_steps)])

        self.f, self.ax = plt.subplots(1,1)

        self.ax.set_xscale('log')
        self.ax.set_yscale('log')
        self.ax.set_xlabel(label_variable[file.axes[0]])
        self.ax.set_ylabel(label_variable[file.axes[1]])
        self.ax.set_xlim(self.xaxis_list[0], self.xaxis_list[-1])
        self.ax.set_ylim(self.yaxis_list[0], self.yaxis_list[-1])

    def _yaxis_log_steps(self,data):
        steps = round((np.log10(data[-1][1]) - np.log10(data[0][1]))/(np.log10(data[1][1]) - np.log10(data[0][1])) + 1)
        return steps
    
    def __len__(self):
        return len(self._data_setup)

line_styles = ['solid','dashed','dotted','dashdot']


exp_style = { #fill,linewidth,R,G,B,transparency
    'NA62':         [False,0.5,0.5,0.3,0.2,0.3,],
    'HIKE':         [False,0.5,0.5,0.3,0.2,0.3,],
    'CHARM':        [True,0.5,0.1,0.1,0.1,0.3,], 
    'NuCal':        [True,0.5,0.5,0.5,0.5,0.3,],
    'SHiP':         [False,0.5,0.1,0.1,0.1,0.3,],
    'SHiPecn4':     [False,0.5,0.1,0.1,0.1,0.3,],
    'DarkQuest':    [False,0.5,0.1,0.1,0.1,0.3,],
    'DarkQuestPhase2': [False,0.5,0.1,0.5,0.1,0.3,],
    'DUNE':         [False,0.5,0.1,0.1,0.1,0.3,],
    'SHADOWS':      [False,0.5,0.1,0.1,0.1,0.3,],
    'KOTOpnn':      [True,0.5,0.1,0.1,0.1,0.3,],
    'KOTOexclPnn':  [False,0.5,0.1,0.1,0.1,0.3,],
    'KOTOdump':     [False,0.5,0.1,0.1,0.1,0.3,],
    'KOTO2pnn':     [False,0.5,0.1,0.1,0.1,0.3,],
    'KOTO2dump':    [False,0.5,0.1,0.1,0.1,0.3,],
    'E137':         [True,0.5,0.1,0.1,0.1,0.3,],
    'E141':         [True,0.5,0.1,0.1,0.1,0.3,],
    'NuTeV':        [True,0.5,0.1,0.2,0.1,0.3,],
    'BEBC':         [True,0.5,0.1,0.3,0.1,0.3,],
    'BEBCcuboid':   [True,0.5,0.1,0.3,0.1,0.3,],
    'ORCA':         [False,0.5,0.1,0.5,0.1,0.5,]
}



label_variable = {
    'mX':       '$m_X$ [$\mathrm{GeV}/c^2$]',
    'Y':        '$\mathrm{sin}^2 \\theta$',
    'eps':      '$\epsilon$',
    'gY':       '$g_{Y}$',
    'GammaX':   '$\Gamma_X$ [$\mathrm{GeV}$]',
    'tauX':     '$\\tau_X$ [$\mathrm{ps}$]',
    'BRdecay':  '$\mathrm{BR}_\mathrm{decay} \\times \mathrm{BR}_\mathrm{production}$ [-]',
    'BRprod':   '$\mathrm{BR}_\mathrm{production}$ [-]',
    'CBB':      '$c_{BB}/\Lambda$ [$\mathrm{GeV}^{-1}$]',
    'CWW':      '$c_{WW}/\Lambda$ [$\mathrm{GeV}^{-1}$]',
    'CGG':      '$c_{GG}/\Lambda$ [$\mathrm{GeV}^{-1}$]',
    'Cll':      '$c_{\ell\ell}/\Lambda$ [$\mathrm{GeV}^{-1}$]',
    'Cqq':      '$c_{qq}/\Lambda$ [$\mathrm{GeV}^{-1}$]',
    'U2' :      '$U^2$',
    'U2el':     '$U^2_{e}$',
    'U2mu':     '$U^2_{\mu}$',
    'U2tau':    r'$U^2_{\tau}$'
}


if __name__ == "__main__":
    sys.exit(main())
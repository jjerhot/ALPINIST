#!/usr/bin/env python3

import numpy as np
from pathlib import Path
import sys
import os

# import plotly.figure_factory as ff
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mpltern
from scipy.interpolate import InterpolatedUnivariateSpline


from ALP_rescale.hnl import hnl_setup as setup
from ALP_rescale.general.setup import n_events_90CL, PoT

home_dir = Path(os.path.dirname(os.path.realpath(__file__))) / '../../'

def main(input_exps, combine=False ):
    '''Alternative analogue to plot_exclusion for generating ternary plots for hnls.'''

    #experiment:
    experiments = []
    if combine:
        experiments = ["combined"]
        print("[Info:] \t Checking combined datasets")
    else:
        if input_exps == "":
            experiments = setup.experiments.keys()
            print("[Info:] \t Selected all experiments available")
        else:
            for exp in input_exps:
                if exp in setup.experiments.keys():
                    experiments.append(exp)
                else:
                    raise("[Error:] \t Experiment " + exp + " not available. Experiment modes available: exp = " + ' | '.join(setup.experiments.keys()) + ". If not specified, running over all experiments available.")
            print("[Info:] \t Selected experiments:", ', '.join(experiments))

    main_loop(experiments)

    return

def main_loop(experiments = []):
    # find files with common x,y-axes and for the same model
    files_common = []
    filenames_common = set()
    n_exp = 0

    for exp in experiments: # checking for files that selected experiments have in common
        if exp == "combined":
            path = Path(home_dir,'tab_toPlot/combined') / 'hnl'
        else:
            path = Path(home_dir,'tab_toPlot') / setup.experiments[exp] / 'hnl'
        if not path.exists() or not path.is_dir():
            raise FileNotFoundError("[Error:] \t directory corresponding to",exp,"not found")
        files_exp = []
        for file in path.iterdir():
            if not "U2ternary" in file.name: continue
            if exp == "combined":
                files_exp.append(input_file(file.name))
            else:
                if file.name.split("_")[0] == exp:
                    infile = input_file(file.name) 
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

    #main loop evaluating yields in terms of parameter limits
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
                path = Path(home_dir,'tab_toPlot/combined') / 'hnl'
            else:
                path = Path(home_dir,'tab_toPlot') / setup.experiments[exp_selected] / 'hnl'
            if not path.exists() or not path.is_dir():
                raise FileNotFoundError("[Error:] \t directory corresponding to",exp_selected,"not found")
            files_exp = []
            for file in path.iterdir():
                if file.name.split("_")[0] == exp_selected and "U2ternary" in file.name:
                    infile = input_file(file.name)
                    if infile.label == exp_selected:
                        files_exp.append(infile)
                elif exp_selected == "combined":
                    files_exp.append(input_file(file.name))
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
            plot = plot_setup(file_chosen,only_data = True)
            if file_chosen.experiment == "combined":
                data = np.loadtxt(Path(home_dir,'tab_toPlot/combined/hnl/')/str(file_chosen))
                nevts = n_events_90CL['hnl'][file.label.split("-")[0]]
            else:
                
                data = np.loadtxt(Path(home_dir,'tab_toPlot',file_chosen.experiment,'hnl')/str(file_chosen))
                nevts = n_events_90CL['hnl'][file.label]
            
            contour = plot.contour_data

            #extract and export contour
            contourNevts = str("{:.2f}".format(nevts))
            outPath = Path(home_dir,'Figures/contours/'+file_chosen.experiment+'/hnl/')
            if not os.path.isdir(outPath): os.makedirs(outPath)
            np.savetxt(outPath/'contour-'+contourNevts+'_'+str(file_chosen),contour,fmt='%.4e')
            print('[Info:] file contour-'+contourNevts+'_'+str(file_chosen)+' saved to Figures/contours/'+file_chosen.experiment+' directory')

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
        print("[Info:] \t Chosen file:", file_chosen)

        common_setup = None

        #categorize files per experiment
        files_exp = {}
        for file in files_common:
            if file.label in files_exp.keys():
                files_exp[file.label].append(file)
            else:
                files_exp[file.label] = [file]

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
                    common_setup = plot_setup(file)
                if file.experiment == "combined":
                    data = np.loadtxt(Path(home_dir,'tab_toPlot/combined/hnl')/str(file))
                    nevts = n_events_90CL['hnl'][file.label.split("-")[0]]
                else:
                    data = np.loadtxt(Path(home_dir,'tab_toPlot/'+file.experiment+'/hnl/')/str(file))
                    nevts = n_events_90CL['hnl'][file.label]
                if len(data) != len(common_setup):
                    print("[Warning:] \t File",str(file),"has different size than the referential file, skipping")
                    continue
                n_matched_files+=1

        if mode_chosen == 1: #show plot
            common_setup.fig.show()
            continue

        elif mode_chosen == 2: #export plot
            exp_path = Path(home_dir, "Figures/"+file.experiment+"/hnl/")
            if not os.path.exists(exp_path): os.makedirs(exp_path)
            common_setup.fig.write_image(exp_path+file_chosen+'.pdf')
            print("[Info:] file "+file_chosen+'.pdf saved to Figures/'+file.experiment+'/hnl/ directory')
            continue

        else:
            print("[Info:] Exitting")
            exit(0)
        
    return

class input_file:
    '''
    Yield map file parameters 
    '''
    def __init__(self, filename = ""):
        self._filename = filename
        self.name_for_comparing = self._filename.replace(".dat","")
        self.experiment = ""
        self.label = ""
        self.modes_prod = []
        self.modes_decay = []
        self.mass = 0.
        self.label = self.name_for_comparing.split("_")[0]
        if self.label in setup.experiments.keys(): self.experiment = setup.experiments[self.label]
        else: self.experiment = "combined"

        self.name_for_comparing = self.name_for_comparing.replace(self.label+"_","")

        for mode in setup.channels_production:
            if mode in self.name_for_comparing:
                self.name_for_comparing = self.name_for_comparing.replace(mode,"")
                self.modes_prod.append(mode)

        for mode in setup.channels_decay:
            if mode in self.name_for_comparing:
                self.name_for_comparing = self.name_for_comparing.replace(mode,"")
                self.modes_decay.append(mode)

        while self.name_for_comparing[-1] == '_' or self.name_for_comparing[-1] == '-': self.name_for_comparing = self.name_for_comparing[:-1]

        massLabel =  self.name_for_comparing.split("_")[1]
        if "mass" in massLabel: self.mass = self.name_for_comparing.split("_")[1][4:][:-3]

    def __str__(self):
        return self._filename

class plot_setup:
    '''
    Setup of the exclusion plot for ternary plot representation
    '''
    def __init__(self, file : input_file, only_data = False):
        if str(file) == "":
            print("[Error:] \t Input file object empty")
            exit()
        
        self.load_contour(file)

        self.title      = r"$U^2 \text{ sensitivity (} m_N="+ file.mass +r"\,\mathrm{MeV}, \text{"+ file.experiment +"@}"+"{:.2E}".format(PoT[file.experiment]).replace("E+",r"\times 10^{")+r"}\,\mathrm{PoT}\text{)}$"
        self.ncontours  = 50
        self.interpmode = 'cartesian'
        self.colorscale = 'viridis_r'
        min_U2 = np.min(self.contour_data[:,3])
        max_U2 = np.max(self.contour_data[np.where(self.contour_data[:,3]<1.),3])
        self.ticks = np.logspace(np.log10(min_U2), np.log10(max_U2), 5)
        self.tick_labels = [ "{:.1E}".format(tick) for tick in self.ticks]
        if not only_data:
            self.fig = plt.figure(figsize=(5.8, 4.8))
            levels = np.logspace(np.log10(min_U2), np.log10(max_U2), self.ncontours)
            ax = self.fig.add_subplot(1, 1, 1, projection='ternary')
            cs = ax.tricontourf(*self.contour_data.T, levels=levels, norm=colors.LogNorm(vmin=min_U2, vmax=max_U2), cmap=self.colorscale)
            ax.set_title(self.title )
            ax.set_tlabel("$U_e^2$")
            ax.set_llabel(r"$U_\mu^2$")
            ax.set_rlabel(r"$U_\tau^2$")
            cax = ax.inset_axes([1.05, 0.05, 0.1, 0.9], transform=ax.transAxes)
            colorbar = self.fig.colorbar(cs, cax=cax, ax = ax, ticks=self.ticks)
            colorbar.ax.set_title('$U^2$')
            colorbar.ax.set_yticklabels(self.tick_labels)
            plt.show()


    def load_contour(self, file : input_file):
        """Load and interpret yield map

        Args:
            file (input_file): file to load
        """
        if str(file) == "":
            print("[Error:] \t Input file object empty")
            exit()

        if file.experiment in setup.experiments:
            self._data_setup = np.loadtxt(Path(home_dir,'tab_toPlot/'+file.experiment+'/hnl/')/str(file))
        else:
            self._data_setup = np.loadtxt(Path(home_dir,'../tab_toPlot/combined/hnl/')/str(file))
        
        data_points_dict = {(0.,0.): []}
        for line in self._data_setup: 
            if data_points_dict.get((line[1],line[2])) == None: data_points_dict[(line[1], line[2])] = [(line[0],line[3])]
            else: data_points_dict[(line[1],line[2])].append((line[0],line[3]))
        U2_lowlim_90CL = []
        contour_data = []

        for coords in data_points_dict.keys(): # evaluating the expected number of events on a regular triangular grid for the lowest match with n_events_90CL 
            U2s, Ns = np.transpose(data_points_dict[coords])
            roots = InterpolatedUnivariateSpline(U2s, Ns - n_events_90CL['hnl'][file.label]).roots()
            U2_lowlim_90CL= np.min(roots) if roots.size else 1
            contour_data.append([*coords, max(1-coords[0]-coords[1], 0.), U2_lowlim_90CL])

        self.contour_data = np.array(contour_data)

    def __len__(self):
        return len(self._data_setup)

#.. to be commented out 


if __name__ == "__main__":
    sys.exit(main())
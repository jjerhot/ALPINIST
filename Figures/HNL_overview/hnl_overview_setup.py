import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

rc('text',usetex = True)
rc('font',**{'family':'serif'})
rc('xtick',**{'top':True,'direction':'in'})
rc('ytick',**{'right':True,'direction':'in'})
rc('legend',**{'handlelength':1.,'frameon':False})
rc('axes', labelsize=15)
rc('text.latex',preamble= ' '.join([ r'\usepackage{amsmath}',r'\usepackage{siunitx}']))


def plot_patch_ElExcluded(ax, alpha_=0.5, color_="gray"):
    data = np.genfromtxt("other_exclusions/Excluded_1-0-0.csv",delimiter=',', filling_values=np.nan)
    data_Triumf = data[2:,:2]
    data_KEK = data[2:,2:4]
    data_PS191 = data[2:,4:6]
    data_DELPHI = data[1:,8:10]
    data_CMS =data[1:,10:]
    data_CMS_2024 = np.genfromtxt("other_exclusions/CMS2024_U2s.csv",delimiter=',', filling_values=np.nan)[:,:2]
    data_CMS_138invFB = np.genfromtxt("other_exclusions/CMS_displaced_138invFB.csv",delimiter=',', filling_values=np.nan)[:,:2]
    data_NA62 = np.loadtxt("other_exclusions/ExcludedNA62_1-0-0.csv", delimiter=',',skiprows=2)
    data_BESIII = 10**np.loadtxt("other_exclusions/BESIII_el_data.csv", delimiter=',',skiprows=2)
    data_Belle = 10**np.loadtxt("other_exclusions/Belle_el_data.csv", delimiter=',',skiprows=2)
    data_collected = np.array([[1e-1,1]],dtype=np.float64)
    for data_ in [data_Triumf,data_KEK,data_PS191,data_DELPHI,data_BESIII,data_Belle,data_CMS,data_CMS_2024,data_CMS_138invFB,data_NA62]: 
        xs = data_[~np.isnan(data_[:,0]),0]
        xs = np.concatenate(([xs[0]], xs, [xs[-1]]))
        ys = data_[~np.isnan(data_[:,1]),1]
        ys = np.concatenate(([1.], ys, [1.]))
        data_collected=np.append(data_collected,np.array([xs,ys]).T,axis=0)
        
    combined_path = Path(data_collected)
    patch = PathPatch(combined_path, facecolor=color_, lw=0, alpha=alpha_, fill=True)
    ax.add_patch(patch)
    return 
def plot_patch_MuExcluded(ax, alpha_=0.5, color_="gray"):
    data = np.genfromtxt("other_exclusions/Excluded_0-1-0.csv",delimiter=',', filling_values=np.nan)
    data_PS191 = data[2:,2:4]
    data_E949  = data[2:,4:6]
    data_KEK   = data[2:,6:8]
    data_Belle = 10**np.loadtxt("other_exclusions/Belle_mu_data.csv", delimiter=',',skiprows=2)
    data_CMS = np.loadtxt("other_exclusions/ExcludedCMS_0-1-0.csv", delimiter=',',skiprows=2)
    data_NA62 = np.loadtxt("other_exclusions/ExcludedNA62_0-1-0.csv", delimiter=',',skiprows=2)
    data_CMS_2024 = np.genfromtxt("other_exclusions/CMS2024_U2s.csv",delimiter=',', filling_values=np.nan)[:,2:4]
    data_CMS_138invFB = np.genfromtxt("other_exclusions/CMS_displaced_138invFB.csv",delimiter=',', filling_values=np.nan)[:,2:4]
    data_collected = np.array([[1e-1,1],[10,1],[10,1e-2],[1e-1,1e-2]],dtype=np.float64)
    for data_ in [data_PS191,data_E949,data_KEK,data_Belle,data_CMS,data_CMS_2024,data_CMS_138invFB,data_NA62]: #
        xs = data_[~np.isnan(data_[:,0]),0]
        xs = np.concatenate(([xs[0]], xs, [xs[-1]]))
        ys = data_[~np.isnan(data_[:,1]),1]
        ys = np.concatenate(([1.], ys, [1.]))
        data_collected=np.append(data_collected,np.array([xs,ys]).T,axis=0)
        
    combined_path = Path(data_collected)
    patch = PathPatch(combined_path, facecolor=color_, lw=0, alpha=alpha_, fill=True)
    ax.add_patch(patch)
    return 
def plot_patch_TauExcluded(ax, alpha_=0.5, color_="gray"):
    data_T2K = 10**np.loadtxt( "other_exclusions/ExcludedT2K_0-0-1.csv", delimiter=',',skiprows=0)
    data_Delphi = 10**np.loadtxt("other_exclusions/ExcludedDELPHI_0-0-1.csv", delimiter=',',skiprows=0)
    data_Belle = np.loadtxt("other_exclusions/ExcludedBelle_0-0-1.csv", delimiter=',',skiprows=2)
    data_CMS_2024 = np.genfromtxt("other_exclusions/CMS2024_U2s.csv",delimiter=',', filling_values=np.nan)[:,4:]
    data_collected = np.array([[1e-1,1]],dtype=np.float64)
    for data_ in [data_T2K,data_Delphi,data_Belle,data_CMS_2024]:
        xs = data_[~np.isnan(data_[:,0]),0]
        xs = np.concatenate(([xs[0]], xs, [xs[-1]]))
        ys = data_[~np.isnan(data_[:,1]),1]
        ys = np.concatenate(([1.], ys, [1.]))
        data_collected=np.append(data_collected,np.array([xs,ys]).T,axis=0)
        
    combined_path = Path(data_collected)
    patch = PathPatch(combined_path, facecolor=color_, lw=0, alpha=alpha_, fill=True)
    ax.add_patch(patch)
    return 
def plot_patch_ElExcludedBD(ax, alpha_=0.5, color_="darkgray"):
    
    data_CHARM = np.loadtxt("other_exclusions/ExcludedCHARM.csv", delimiter=',',skiprows=2)[:,:2]
    data_BEBC = np.loadtxt("other_exclusions/BEBC_BEBC_Mu.csv", delimiter=',',skiprows=2)
    for data_ in [data_CHARM,data_BEBC]: #
        xs = data_[~np.isnan(data_[:,0]),0]
        ys = data_[~np.isnan(data_[:,1]),1]
        ax.fill_between(xs, ys, np.ones_like(xs), lw =0, alpha=alpha_, color=color_)
    return 
def plot_patch_MuExcludedBD(ax,alpha_=0.5, color_="darkgray"):
    data_CHARM = np.loadtxt("other_exclusions/ExcludedCHARM.csv", delimiter=',',skiprows=2)[:,2:]
    data_NuTeV = np.loadtxt("other_exclusions/ExcludedNuTeV_0-1-0.csv", delimiter=',',skiprows=2)
    data_BEBC = np.loadtxt("other_exclusions/BEBC_BEBC_Mu.csv", delimiter=',',skiprows=2)
    for data_ in [data_CHARM,data_BEBC,data_NuTeV]: #
        xs = data_[~np.isnan(data_[:,0]),0]
        ys = data_[~np.isnan(data_[:,1]),1]
        ax.fill_between(xs, ys, np.ones_like(xs), lw =0, alpha=alpha_, color=color_)
    return 
def plot_patch_TauExcludedBD(ax, alpha_=0.5, color_="darkgray"): 
    return 

def plot_ComparisonPythia83():
    fig = plt.figure(figsize=(10,10),)
    gs = fig.add_gridspec(8,8)
    axs = [fig.add_subplot(gs[:4,:4]),fig.add_subplot(gs[:4,4:]),fig.add_subplot(gs[4:,1:5])]
    ax_leg = fig.add_subplot(gs[5:7,6:-1])
    bd_color='lightgray'
    past_color='gray'
    mesonProduction = 'Pythia83'
    # fig.suptitle(r'HNL sensitivity SHiP$@6\times10^{20}$ PoT')
    for ax, coupling, Bc, tex_coupling, plot_otherExcl,plot_otherExclBD in zip(axs,["El","Mu","Tau"], ["BC6", "BC7","BC8"],["e",r"\mu",r"\tau"],[plot_patch_ElExcluded,plot_patch_MuExcluded,plot_patch_TauExcluded],[plot_patch_ElExcludedBD,plot_patch_MuExcludedBD,plot_patch_TauExcludedBD]):##
        insert = "U2tau_fixed-U2el0-U2mu0" if coupling == "Tau" else  "U2el_fixed-U2mu0-U2tau0" if coupling == "El" else "U2mu_fixed-U2el0-U2tau0" 
        ax.set_ylim(5e-11, 1e-2)
        ax.set_xlim(0.1,6)
        ax.set_xscale("log")
        ax.set_yscale("log")
        plot_otherExcl(ax,color_=past_color)
        plot_otherExclBD(ax, color_=bd_color)
        cols = ["peru","darkred","orchid","teal","royalblue","darkmagenta","indigo","darkkhaki","navy"]
        artists, labels = [], []
        for exp2, nCL, lab_, ls_, col_ in zip(["BEBC","CHARM","NuCal","NuTeV","NA62","DarkQuest","DarkQuestPhase2","DUNE","SHiP"],
                                                ["2.30","2.30","2.30","2.30","2.30","10.00","10.00","1.00","0.77"],
                                                ["BEBC recast","CHARM recast","NuCal recast","NuTeV recast","NA62-bd","DarkQuest Phase-I","DarkQuest Phase-II","DUNE ND","SHiP"],
                                                ["solid","solid","solid","solid","dashed","dashed","dotted","dotted","dotted"],cols):
            if exp2 == "BEBC" and "7" in Bc: nCL = "3.45"
            exp = "SHiP" if exp2 == "SHiPecn3" else "DarkQuest" if exp2 == "DarkQuestPhase2" else exp2
            data_l = np.loadtxt("../../tab_toPlot/"+exp+"/hnl/"+exp2+"_mX_"+insert+f"_{mesonProduction}_Xnorm_l.dat") 
            data_u = np.loadtxt("../../tab_toPlot/"+exp+"/hnl/"+exp2+"_mX_"+insert+f"_{mesonProduction}_Xnorm_u.dat")
            mN, U2 = np.unique(data_l[:,0]), np.unique(data_l[:,1])
            nCL_float = float(nCL)
            invl = np.divide(nCL_float**2,data_l[:,2]+1e-50)
            confidence_band = np.where(data_l[:,2]>nCL_float, invl, np.where(data_u[:,2] < nCL_float, data_u[:,2], np.minimum(data_u[:,2],invl)))
            confidence_band = confidence_band.reshape((mN.size,U2.size)).T
            ax.contourf(mN, U2, confidence_band, (nCL_float,1e20), colors = [col_,col_],alpha=0.33)

            data_c = np.loadtxt("../../tab_toPlot/"+exp+"/hnl/"+exp2+"_mX_"+insert+f"_{mesonProduction}_Xnorm.dat") 
            ax.contour(mN, U2, data_c[:,2].reshape((mN.size,U2.size)).T, levels=[nCL_float] , colors = col_, linestyles = ls_)

            artists.append((Patch(color=col_,alpha= 0.5,lw=0), Line2D([], [], linestyle=ls_, color=col_)))
            labels.append(lab_)
    
        ax.set_xlabel(r"$m_\mathrm{N}$ in GeV")
        ax.set_ylabel(f"$U^2_{tex_coupling}$")
        ax.tick_params("both", which = "both", direction = "in", bottom = True, left = True, top = True, right = True)
        ax.text(0.025,0.025,r"\noindent{\bf "+ Bc + r"} sens. estimates\par\noindent{($90\,\%\,\mathrm{CL}$ excl. bound)}", transform = ax.transAxes, horizontalalignment='left',verticalalignment='bottom', fontsize='large')

    background_artists=[Patch(visible=False), Patch(color=past_color,alpha= 0.5,lw=0), Patch(color=bd_color,alpha= 0.5,lw=0)]
    background_labels =["Past searches", "Precision and Colliders", "Beam Dumps"]
    ax_leg.axis("off")
    leg = ax_leg.legend(artists+background_artists, labels+background_labels, frameon = True,fancybox=False, framealpha=0.8, edgecolor='1.', handlelength = 1., alignment='left', bbox_to_anchor =(2.25,0.5), loc='center right', title =r"\noindent{\large Sens. estimates using \textbf{ALP\normalsize{I}\large{NIST}}}", fontsize='large' )
    for item, label in zip(leg.legend_handles, leg.texts):
        if label._text  == "Past searches":
            width=item.get_window_extent(fig.canvas.get_renderer()).width
            label.set_ha('left')
            label.set_position((-2*width,0))
    ax.set_ybound(1e-9, 1e-2)
    ax.set_xbound(0.1,4)
    fig.tight_layout()
    fig.savefig("PythiaHNLComparisons.pdf")

    plt.show()

    return 
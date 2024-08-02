import sys
import pythia8
import numpy as np
def number_to_3sigfigs_str(number):
    assert 0 <= number < 1e12, f"number {number} outside of interpretation range (0, 999e9)"
    translator = ['','k','M','B']
    number_split = f"{number:3.3e}".split("e")
    order_1000 = int(int(number_split[1])/3)
    raise_by = int(number_split[1])%3
    sigfigs = str(float(f"{(float(number_split[0])):.2e}")*10**raise_by)
    while '.' in sigfigs and sigfigs[-1] in ['0','.']: sigfigs = sigfigs[:-1]
    return sigfigs+translator[order_1000]

pythia_dir = "/Users/Jonathan/misc/Pythia"

cfg = open(pythia_dir+"/examples/Makefile.inc")
lib = pythia_dir+"/lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)
for nprod in [100000,300000,1000000]:
    for p_beam in [30, 70, 120, 400, 800]:
        print("Prdouction for impact on fixed proton by proton with momentum",p_beam)
        pythia = pythia8.Pythia("-v", False)
        
        pythia_modifiers = ["SoftQCD:all = on","pdf:pHardSet = 2", "111:mayDecay = no","Beams:idA = 2212","Beams:idB = 2212","Beams:eA = " + str(np.sqrt(p_beam**2+0.93827**2)), "Beams:eB = "+str(0.93827), "Beams:frameType = 2"]
        for modifier in pythia_modifiers+["Next:numberShowEvent = 5"]:
            pythia.readString(modifier)

        pythia.init() # initialize
        meson_list = []
        for _ in range(nprod):
            pythia.next()
            for isubEvent in range(pythia.event.size()):
                evt = pythia.event[isubEvent]
                if abs(evt.id()) not in [111, 221, 331, 113, 213, 223, 333]: continue
                meson_list.append([evt.id(),evt.px(),evt.py(),evt.pz()])
        
        file_info  =  "Charmed meson momenta as generated for " + str(nprod) +" impining protons, with PYTHIA"+str(pythia.parm("Pythia:versionNumber"))+" and modifiers [" + ', '.join(pythia_modifiers)+']'
        np.savetxt(f"./softQCD_{p_beam}GeV_{number_to_3sigfigs_str(nprod)}.dat", meson_list, header=file_info, fmt='%i %4f %4f %3f')
    
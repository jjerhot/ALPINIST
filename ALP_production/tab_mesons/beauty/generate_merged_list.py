import sys
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

for nprod in [100000, 300000, 1000000]:
    for p_beam, sigma_bb_2parton, sigma_qq_3parton in zip([400],[1.866],[94.01]):

        bb_2parton_path = f'./hsbbbar_{p_beam}GeV_1M.dat'
        bb_3parton_path = f'./3partonbbbar_{p_beam}GeV_10M.dat'

        with open(bb_2parton_path) as f: header_2parton = f.readline().split('modifiers')[-1][1:-1]
        with open(bb_3parton_path) as f: header_3parton = f.readline().split('modifiers')[-1][1:-1]

        n_2parton_file = 1e6
        n_3parton_file = 1e7

        bb_2parton_file = np.loadtxt(bb_2parton_path,skiprows=1)
        bb_3parton_file = np.loadtxt(bb_3parton_path,skiprows=1)
        hadrochance = bb_2parton_file.shape[0 ] / n_2parton_file

        sigma_bb_3parton = sigma_qq_3parton * bb_3parton_file.shape[0] / hadrochance / n_3parton_file
        sigma_bb  = sigma_bb_2parton + sigma_bb_3parton

        N_Draw_mesons_2part = int(sigma_bb_2parton / sigma_bb * nprod * hadrochance)
        N_Draw_mesons_3part = int(sigma_bb_3parton / sigma_bb * nprod * hadrochance)

        idcs_2 = np.random.choice(bb_2parton_file.shape[0], size = N_Draw_mesons_2part ,replace=True)
        idcs_3 = np.random.choice(bb_3parton_file.shape[0], size = N_Draw_mesons_3part ,replace=True)

        meson_list = []
        for idx in idcs_2: meson_list.append(bb_2parton_file[idx])
        for idx in idcs_3: meson_list.append(bb_3parton_file[idx])
        meson_list = np.array(meson_list)
        np.random.shuffle(np.array(meson_list))
        file_info = f"statistical combination of files {bb_2parton_path.split('/')[-1]} and {bb_3parton_path.split('/')[-1]}, assuming the same B meson hadronisation rate with sigma_bb = {sigma_bb_2parton} and {sigma_bb_3parton} for 2 and 3 parton processes respectively. These were generated with modifiers {header_2parton} and {header_3parton}." 

        np.savetxt(f"./hsbbbar_statcomb23partprod_{p_beam}GeV_{number_to_3sigfigs_str(nprod)}.dat", meson_list, header=file_info, fmt='%i %4f %4f %3f')
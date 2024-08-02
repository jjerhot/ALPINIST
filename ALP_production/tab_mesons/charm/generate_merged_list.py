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
    for p_beam, sigma_cc_2parton, sigma_qq_3parton in zip([400],[4.053],[174.7]):

        cc_2parton_path = f'./hsccbar_{p_beam}GeV_1M.dat'
        cc_3parton_path = f'./3partonccbar_{p_beam}GeV_10M.dat'

        with open(cc_2parton_path) as f: header_2parton = f.readline().split('modifiers')[-1][1:-1]
        with open(cc_3parton_path) as f: header_3parton = f.readline().split('modifiers')[-1][1:-1]

        n_2parton_file = 1e6
        n_3parton_file = 1e7

        cc_2parton_file = np.loadtxt(cc_2parton_path,skiprows=1)
        cc_3parton_file = np.loadtxt(cc_3parton_path,skiprows=1)
        hadrochance = cc_2parton_file.shape[0 ] / n_2parton_file

        sigma_cc_3parton = sigma_qq_3parton * cc_3parton_file.shape[0] / hadrochance / n_3parton_file
        sigma_cc  = sigma_cc_2parton + sigma_cc_3parton

        N_Draw_mesons_2part = int(sigma_cc_2parton / sigma_cc * nprod * hadrochance)
        N_Draw_mesons_3part = int(sigma_cc_3parton / sigma_cc * nprod * hadrochance)

        idcs_2 = np.random.choice(cc_2parton_file.shape[0], size = N_Draw_mesons_2part ,replace=True)
        idcs_3 = np.random.choice(cc_3parton_file.shape[0], size = N_Draw_mesons_3part ,replace=True)

        meson_list = []
        for idx in idcs_2: meson_list.append(cc_2parton_file[idx])
        for idx in idcs_3: meson_list.append(cc_3parton_file[idx])
        meson_list = np.array(meson_list)
        np.random.shuffle(np.array(meson_list))
        file_info = f"statistical combination of files {cc_2parton_path.split('/')[-1]} and {cc_3parton_path.split('/')[-1]}, assuming the same D meson hadronisation rate with sigma_cc = {sigma_cc_2parton} and {sigma_cc_3parton} for 2 and 3 parton processes respectively. These were generated with modifiers {header_2parton} and {header_3parton}." 

        np.savetxt(f"./hsccbar_statcomb23partprod_{p_beam}GeV_{number_to_3sigfigs_str(nprod)}.dat", meson_list, header=file_info, fmt='%i %4f %4f %3f')
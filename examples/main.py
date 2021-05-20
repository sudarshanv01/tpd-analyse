
""" Main Script to plot Figure 1 of the Gold paper """

import numpy as np
from glob import glob
import os, sys
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
import mpmath as mp
from ase.io import read
from ase.db import connect
import matplotlib
from ase.units import kB
import csv
from scipy.optimize import curve_fit
from matplotlib import cm
import matplotlib.pyplot as plt
from tpd_analyse.tpd import PlotTPD

if __name__ == '__main__':

    ## This is where all the output is stored
    output = 'output_figures/'
    os.system('mkdir -p ' + output)

    fig, ax = plt.subplots(1, 1, figsize=(8,6), squeeze=False)

    ## Files with the different exposure
    files = glob('input_TPD/Au_211/*.csv')
    ## Order of the reaction
    order = 1

    # For 211 TPD
    T_switch_211 = [170] #[110, 170] #[110, 165, 175] #K converting between 111 and 211 step
    T_max_211 = 250 #K Where the TPD spectra ends
    T_rate_min_211 = [250, 300] #K Where the rate becomes zero - baseline
    beta_211 = 3 #k/s Heating rate
    data_211 = [T_switch_211, T_max_211, T_rate_min_211, beta_211]

    atoms = read('input_TPD/co.traj')
    thermo_gas = 0.00012 * np.array([2121.52, 39.91, 39.45])
    vibration_energies_gas = IdealGasThermo(thermo_gas, atoms = atoms, \
            geometry='linear', symmetrynumber=1, spin=0)
    vibration_energies_ads = HarmonicThermo(0.00012 * np.array([2044.1, 282.2, 201.5, 188.5, 38.3, 11.5]))
    cmap = cm.get_cmap('viridis', 12)
    TPDClass = PlotTPD(exp_data=files,
                        order=order,
                        thermo_ads=vibration_energies_ads,
                        thermo_gas=vibration_energies_gas,
                        plot_temperature=np.linspace(100, 500, 50), 
                        constants=data_211,
                        )
    TPDClass.get_results()
    ## Nested dict of facets and 
    for surface_index, facet_no in enumerate(TPDClass.results):
        # Plot the TPD in normalised form
        for index_exposure, exposure in enumerate(sorted(TPDClass.results[surface_index])):
            if surface_index == len(TPDClass.results)-1:
                theta = TPDClass.theta_rel[surface_index][exposure] * TPDClass.theta_sat[surface_index][exposure] 
                ax[0,0].plot(theta,
                                TPDClass.Ed[surface_index][exposure],
                                '.', \
                                color=cmap(index_exposure), \
                                label=str(exposure)+'L'
                                )
        ax[0,0].set_xlabel(r'$\theta$ / ML')
        ax[0,0].set_ylabel(r'$rate_{norm}$ / arb. units')
        fig.tight_layout()
        fig.savefig(os.path.join(output, 'figure.pdf'))

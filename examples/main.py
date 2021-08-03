
"""
Example to generate the list of plots in:
10.26434/chemrxiv.14525496.v2

The idea is to determine the contribution of 
1. zero-coverage energy
2. configurational entropy
3. adsorbate-adsorbate interaction

to the overall free energy of adsorption for a first order process.

Inputs:
1. the path to the directory containing the input files
2. the path to the directory containing the plots
3. Inputs about each TPD

"""

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

    # Input files for CO on Gold 211
    files = glob('input_TPD/Au_211/*.csv')

    # Output directory
    out_dir = 'plots/'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # plot the figures here
    fig, ax = plt.subplots(4, 1, figsize=(8,6), squeeze=False, constrained_layout=True)
    cmap = cm.get_cmap('viridis', 12)

    # order of the reaction
    order = 1

    # temperature at which the switch occurs if there is more than one peak
    T_switch_211 = [170] 
    # temperature at which the TPD ends
    T_max_211 = 250 
    # baseline is computed at this temperature
    T_rate_min_211 = [250, 300] 
    # Heating rate
    beta_211 = 3 
    # combine everything into one list
    data_211 = [T_switch_211, T_max_211, T_rate_min_211, beta_211]

    # details of the gas phase molecule
    # atoms object for the gas molecule
    atoms = read('input_TPD/co.traj')
    # Vibrational frequencies in eV
    thermo_gas = 0.00012 * np.array([2121.52, 39.91, 39.45])
    # Get the free energies of the gas phase molecule
    vibration_energies_gas = IdealGasThermo(thermo_gas, atoms = atoms,
            geometry='linear', symmetrynumber=1, spin=0)
    # Harmonic frequencies in cm^-1
    vibration_energies_ads = HarmonicThermo(0.00012 * np.array([2044.1, 282.2, 201.5, 188.5, 38.3, 11.5]))

    # TPD class to get the needed data 
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
                theta = TPDClass.results[surface_index][exposure]['theta_rel'] * TPDClass.results[surface_index][exposure]['theta_sat']
                ax[0,0].plot(theta,
                                TPDClass.results[surface_index][exposure]['Ed'],
                                '.', \
                                color=cmap(index_exposure), \
                                label=str(exposure)+'L'
                                )
        ax[0,0].set_xlabel(r'$\theta$ / ML')
        ax[0,0].set_ylabel(r'$rate_{norm}$ / arb. units')
        fig.savefig(os.path.join(out_dir, 'figure.pdf'))

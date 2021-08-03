TPD analysis for first order desorption process
-----------------------------------------------

## Details

Calculate the adsorption energy, configurational entropy and adsorbate-adsorbate interaction parameter for first order desorption from TPD curves. 

Equations used to fit are described in the manuscipt: How to Extract Adsorption Energies, Adsorbate-adsorbate Interaction Parameters, and Saturation Coverages from Temperature Programmed Desorption Experiments ([link](https://chemrxiv.org/engage/chemrxiv/article-details/60c75888567dfe0feeec6887))

## Installation

The package can be installed via

```
pip install tpd-analyse
```

## Usage

All inputs must be passed to the central class `PlotTPD`. An example of some possible inputs are shown below. 

```
    TPDClass = PlotTPD(exp_data=files,
                        order=order,
                        thermo_ads=vibration_energies_ads,
                        thermo_gas=vibration_energies_gas,
                        plot_temperature=np.linspace(100, 500, 50), 
                        T_switch=T_switch_211,
                        T_max=T_max_211,
                        T_rate_min=T_rate_min_211,
                        beta=beta_211)
```

1. `exp_data`: list of filenames. Each filename must be of the format `exposure_<some_value>.csv` where some_value can be the exposure in any units. 
2. `order`: currently only first order desorption reactions are implemented, so `order=1`
3. `thermo_ads`: Here the `ase` class `HarmonicThermo` can be passed with the required inputs for the adsorbate(s) of interest. More information about the class can be found [here](https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html).
4.  `thermo_gas`: Similar to the `thermo_ads`, only for gas molecules. One option is to use `IdealGasThermo`
5.  `plot_temperature`: A temperature range that the equibirium coverage can be determined in. An example range would be `np.linspace(100, 500)
6.  `T_switch`: Is there are more than one processes occuring in one TPD plot, supply a list of temperatures to switch from one to another in K.
7.  `T_max`: Maximum temperature to consider for *all* the TPDs. 
8.  `T_rate_min`: Minimum temperature at which the baseline correction can be done. 
9.  `beta`: Heating rate in K/s
10.  `correct_background` (optional): bool to decide if the background is corrected
11.  `p` (optional): Pressure in the case of *equilibirum*, only useful if equilbirum coverages are needed
12.  `initial_guess_theta` (optional): Initial guess for the coverage in the Newton root solver
13.  `guess_b`(optional): Guess for the ads-ads interaction. 
14.  `calculate_eq_coverage` (oprtional): Decide is the equilibrium coverage is computed.



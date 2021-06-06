TPD analysis for first order desorption process
-----------------------------------------------

* Details

Calculate the adsorption energy for first order desorption from TPD curves. 

Equations used to fit are described in: 10.26434/chemrxiv.14525496

* Prerequisites

1. `numpy`
2. `scipy`
3. `ase`
4. `matplotlib`

* Usage

An example for the case of CO desorption on Au(211) is shown in `examples/main.py`. In short, what you would need is something of this form

```
TPDClass = PlotTPD(exp_data=files,
                    order=order,
                    thermo_ads=vibration_energies_ads,
                    thermo_gas=vibration_energies_gas,
                    plot_temperature=np.linspace(100, 500, 50), 
                    constants=data_211,
                    )
```

Where 
1. `exp_data` (list) is a list of `glob` files
2. `order` (int) is the order of the reaction, only the first order is currently used
3. `thermo_ads` (class) is the ASE thermochemistry class concerning the adsorbate thermochemistry
4. `thermo_gas` (class) is the ASE thermochemistry class concerning the gas molecules thermochemistry
5. `plot_temperature` (list) is the temperature range that you would like the final plot to be over
6. `constants` (list) is a list of start and end temperature and the heating rate as the -1 element


To run the code, call the class like `TPDClass.get_results()`



"""
Main file for TPD analysis.
"""


import numpy as np
from glob import glob
from pprint import pprint
import os, sys
from scipy.optimize import curve_fit, least_squares, minimize
from tpd_analyse.tools.parser_class import experimentalTPD
import mpmath as mp
from mpmath import mpf
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from ase.io import read
from ase.db import connect
import csv
from ase.units import kB
from scipy import optimize
from scipy.optimize import newton

class PlotTPD():

    def __init__(self, exp_data, order, T_switch, T_max, T_rate_min, beta,
                    thermo_gas, thermo_ads=None, correct_background=True, bounds=[], plot_temperature=np.linspace(100,400),
                    p=101325, initial_guess_theta=0.5, guess_b=0.1, calculate_eq_coverage=True, theta_range=None, fit_quad_ads_ads=False):

        """Perform the temperature programmed desorption analysis for a surface 
        based on configurational entropy an interaction parameters and a zero coverage
        energy term.

        Args:
            exp_data (list): globbed files with csv
            order (int): Order of the reaction 
            T_switch (list): Temperatures at which the reaction rate switches (K)
            T_max (float): Temperature at which the TPD is cutoff (K)
            T_rate_min (float): Temperature range where minimum rate of the TPD is expected (K)
            beta (float): Rate of heating (K/s)
            constants (list): Parameters to parse TPD
            thermo_ads (obj): HarmonicThermo for the adsorbate
            thermo_gas (obj): IdealGasThermo for the gas
            bounds (list, optional): Bounds within to fit the coverage of the TPD. Defaults to [].
            plot_temperature (array, optional): Temperature range to plot the equilbirum coverage. Defaults to np.linspace(100,400).
            p (float, optional): Pressure of gas molecule. Defaults to 101325.
            initial_guess_theta (float, optional): Initial guess for theta. Defaults to 0.5.
            guess_b (float, optional): Initial guess for b. Defaults to 0.1.
            theta_range (list, optional): Range of theta to fit. Defaults to None.

        """

        # Define all the __init__ variables
        self.exp_data = exp_data # globbed files with csv
        self.order = order
        self.thermo_ads = thermo_ads # adsorbate HarmonicThermo
        self.thermo_gas = thermo_gas # gas phase IdealGasThermo
        self.plot_temperature = plot_temperature # temperature for plotting equilibrium coverages
        self.p = p # partial pressure of gas
        self.bounds = bounds
        self.correct_background = correct_background
        self.initial_guess_theta = initial_guess_theta
        self.guess_b = guess_b
        self.T_switch = T_switch
        self.T_max = T_max
        self.T_rate_min = T_rate_min
        self.beta = beta
        self.calculate_eq_coverage = calculate_eq_coverage
        self.theta_range = theta_range
        self.fit_quad_ads_ads = fit_quad_ads_ads

        # Results
        self.norm_results = {} # Normalised results 
        self.results = {} # Final results
        self.theta_eq = {} # Equilibrium coverages
        self.theta_eq_p = {} # Equilibrium coverages positive error
        self.theta_eq_n = {} # Equilibrium coverages negative error
        self.dG = {} # dG for adsorption
        self.temperature_range = {} # temperature range to plot

    def _exponential_fit(self, temperature, a, k):
        """Exponential fit to the tail of the TPD to remove pumping 
        related rates

        Args:
            temperature (list): temperature list based on TPD
            a (float): amplitude of exponent 
            k (float): argument of exponent 

            Returns:
            list: rates for each temperature
        """
        rate = a * np.exp(-k * temperature)
        return rate

    def get_results(self):
        """Perform the TPD analysis
        """
        # T_switch, T_max, T_rate_min, beta = self.constants
        # Do some checks on the temperatures
        if isinstance(self.T_switch, (int,float)):
            self.T_switch = [self.T_switch]

        T_max = self.T_max
        T_rate_min = self.T_rate_min
        beta = self.beta
        T_switch = self.T_switch

        assert T_max >= np.max(T_switch); 'The maximum temperature must be greater than when the switch occurs'
        assert np.max(T_rate_min) >= T_max; 'Maximum temperature of the TPD must be lower than that of the flat region'

        # Create the temperature range based on the switch data
        temperature_ranges = []
        for i in range(len(T_switch)+1):
            if i == 0:
                temperature_ranges.append([0, T_switch[i]])
            elif i == len(T_switch):
                if T_switch[i-1] != T_max:
                    temperature_ranges.append([T_switch[i-1], T_max])

        # range of temperatures for different TPD values
        self.temperature_range = temperature_ranges

        # Get the TPD results which includes background subtraction
        # for each exposure
        for index, f in enumerate(sorted(self.exp_data)):
            exposure = float(f.split('/')[-1].split('.')[0].split('_')[1].replace('p', '.'))
            self.norm_results[exposure] = experimentalTPD(tpd_filename=f,
                                temprange=temperature_ranges,
                                tempmin=T_rate_min,
                                beta=beta,
                                order=self.order, 
                                correct_background=self.correct_background,
                                )

            # Iterate over the different facets in temperature ranges
            for surface_index in range(len(temperature_ranges)):
                
                T_range = temperature_ranges[surface_index]
                # Operate only on the right temperature range
                indices = [ a for a in range(len(self.norm_results[exposure].temperature)) \
                        if T_range[0] < self.norm_results[exposure].temperature[a] < T_range[1] ] 

                # Normalise the data
                self.results.setdefault(surface_index, {})[exposure] = {}
                self.results[surface_index][exposure]['temperature'] = self.norm_results[exposure].temperature[indices]
                self.results[surface_index][exposure]['normalized_rate'] = self.norm_results[exposure].normalized_rate[indices]

                # some variables to make it easy to run
                temperatures = self.results[surface_index][exposure]['temperature']
                rates = self.results[surface_index][exposure]['normalized_rate']

                # For each point get the energy of desorption as a function of the coverage
                data = self._Ed_temp_dependent(
                            temperature=temperatures, 
                            rate=rates,
                            beta=beta,)

                # correct for any nans that could be in place
                args_accept = [i for i in range(len(data[0])) \
                                    if np.isfinite(data[0][i]) and \
                                    data[1][i] > 0]

                # Check if there is a coverage range that is chosen 
                if self.theta_range is not None:
                    # If there is a range of theta values to fit then
                    # only use the data that is within that range
                    args_accept = [i for i in args_accept \
                                    if self.theta_range[0] < data[1][i] < self.theta_range[1]]

                self.results[surface_index][exposure]['Ed'] = data[0][args_accept]
                self.results[surface_index][exposure]['theta_rel'] = data[1][args_accept]

                temperature_fit = self.norm_results[exposure].temperature[indices][args_accept]
                self.results[surface_index][exposure]['temperature_fit'] = temperature_fit

                # Fit the Desorption energy curve to the desorption energy equation
                # First get good initial guesses for parameters
                # For E0 we just take the mean of all the values
                guess_E0 = np.mean(self.results[surface_index][exposure]['Ed']) 
                guess_b = self.guess_b 

                if not self.fit_quad_ads_ads:
                    popt, pcov = curve_fit(\
                    lambda temp, E0, b, theta_sat:  self._fit_Ed_theta(temp, E0, b, theta_sat,
                                                    self.results[surface_index][exposure]['theta_rel']),
                                                    xdata = temperature_fit,
                                                    ydata = self.results[surface_index][exposure]['Ed'],
                                                    p0 = [guess_E0, guess_b, self.initial_guess_theta], 
                                                    )
                else:
                    popt, pcov = curve_fit(\
                    lambda temp, E0, b1, b2, theta_sat:  self._fit_quad_Ed_theta(temp, E0, b1, b2, theta_sat,
                                                    self.results[surface_index][exposure]['theta_rel']),
                                                    xdata = temperature_fit,
                                                    ydata = self.results[surface_index][exposure]['Ed'],
                                                    p0 = [guess_E0] + guess_b + [self.initial_guess_theta], 
                                                    )
                if not self.fit_quad_ads_ads:
                    residual = self._least_sq_Ed_theta(popt, temperature=temperature_fit,
                        theta_rel = self.results[surface_index][exposure]['theta_rel'],
                        Ed_real = self.results[surface_index][exposure]['Ed'],)
                else:
                    residual = self._least_sq_quad_Ed_theta(popt, temperature=temperature_fit,
                        theta_rel = self.results[surface_index][exposure]['theta_rel'],
                        Ed_real = self.results[surface_index][exposure]['Ed'],)
                
                self.results[surface_index][exposure]['E0'] = popt[0]
                if not self.fit_quad_ads_ads:
                    self.results[surface_index][exposure]['b'] = popt[1]
                    self.results[surface_index][exposure]['theta_sat'] = popt[2]
                else:
                    self.results[surface_index][exposure]['b'] = popt[1:3]
                    self.results[surface_index][exposure]['theta_sat'] = popt[3]

                self.results[surface_index][exposure]['error'] = residual
                if not self.fit_quad_ads_ads:
                    self.results[surface_index][exposure]['Ed_fitted'] = self._fit_Ed_theta(temperature_fit, \
                                                    *popt, self.results[surface_index][exposure]['theta_rel']) 
                else:
                    self.results[surface_index][exposure]['Ed_fitted'] = self._fit_quad_Ed_theta(temperature_fit, \
                                                    *popt, self.results[surface_index][exposure]['theta_rel'])
                
                self.results[surface_index][exposure]['configurational_entropy'] = \
                    self._get_configuration_entropy_contribution(temperature_fit, \
                                                self.results[surface_index][exposure]['theta_sat'], \
                                                self.results[surface_index][exposure]['theta_rel'], \
                                        )
                self.results[surface_index][exposure]['ads_ads_interaction'] = \
                    self._get_ads_ads_interaction(temperature_fit, \
                                                self.results[surface_index][exposure]['theta_sat'], \
                                                self.results[surface_index][exposure]['theta_rel'], \
                                                self.results[surface_index][exposure]['b'], \
                                        )

                if not self.calculate_eq_coverage:
                    continue

                # Calculate the coverage at equilbirum
                self.theta_eq.setdefault(surface_index, {})[exposure] = {}
                self.theta_eq_p.setdefault(surface_index, {})[exposure] = {}
                self.theta_eq_n.setdefault(surface_index, {})[exposure] = {} 
                self.dG.setdefault(surface_index, {})[exposure] = {}

                self.theta_eq[surface_index][exposure], self.dG[surface_index][exposure]\
                 = self._get_equilibirum_coverage(
                            E0 = self.results[surface_index][exposure]['E0'],
                            b = self.results[surface_index][exposure]['b'][0] if self.fit_quad_ads_ads else self.results[surface_index][exposure]['b'], 
                            )
                self.theta_eq_p[surface_index][exposure], _ \
                 = self._get_equilibirum_coverage(
                            E0 = self.results[surface_index][exposure]['E0'] + self.results[surface_index][exposure]['error'], 
                            b = self.results[surface_index][exposure]['b'][0] if self.fit_quad_ads_ads else self.results[surface_index][exposure]['b'],  
                            )
                self.theta_eq_n[surface_index][exposure], _ \
                 = self._get_equilibirum_coverage(
                            E0 = self.results[surface_index][exposure]['E0'] - self.results[surface_index][exposure]['error'],
                            b = self.results[surface_index][exposure]['b'][0] if self.fit_quad_ads_ads else self.results[surface_index][exposure]['b'],
                            )
                

    def _Ed_temp_dependent(self, temperature, rate, beta):
        """Gets the desorption energy as a function of the temperature
        1. Do trapezoidal integration to get the coverage by integrating over the 
           rate and temperature 
        2. Get the desorption energy by fitting to the form 
                Ed = -kBT log(-dtheta/dt / mu / theta)
        3. Normalise theta by dividing my maximum coverage

        Args:
            temperature (list): temperatures corresponding to the TPD
            rate (list): rate from the TPD
            beta (float): Rate of heating

        Returns:
            list: Desorption energy and coverage
        """
        h = 4.135e-15 # eV.s

        theta = []
        for i in range(len(temperature)):
            cov = np.trapz(rate[i:], temperature[i:])
            theta.append(cov)
        
        theta = np.array(theta)
        rate = np.array(rate)
        dtheta_dT = np.diff(theta) / np.diff(temperature)
        dtheta_dt = beta * dtheta_dT #rate[0:-1]
        temperature = np.array(temperature)
        nu = kB * temperature[0:-1] / h
        Ed = -kB * temperature[0:-1] * np.log( -1 * dtheta_dt / (nu * theta[0:-1]))

        return [Ed, theta[0:-1]/max(theta[0:-1])]

    def _least_sq_Ed_theta(self, x, temperature, theta_rel, Ed_real):
        E_0, b, theta_sat = x
        Ed_fit = []
        for i in range(len(temperature)):
            Ed = E_0 \
                - kB * temperature[i] * np.log(theta_sat * theta_rel[i] / ( 1 - theta_sat * theta_rel[i] ) ) \
                - b * theta_rel[i] * theta_sat
            Ed_fit.append(Ed)
        residual = Ed_real - Ed_fit
        mea = np.mean([np.abs(a) for a in residual])
        return mea 

    def _least_sq_quad_Ed_theta(self, x, temperature, theta_rel, Ed_real):
        E_0, b1, b2, theta_sat = x
        Ed_fit = []
        for i in range(len(temperature)):
            Ed = E_0 \
            -  kB * temperature[i] * np.log(theta_sat*theta_rel[i] / ( 1 - theta_sat*theta_rel[i]))
            -  b1 * theta_rel[i] * theta_sat
            -  b2 * (theta_rel[i] * theta_sat)**2
            Ed_fit.append(Ed)
        residual = Ed_real - Ed_fit
        mea = np.mean([np.abs(a) for a in residual])
        return mea 

    def _get_configuration_entropy_contribution(self, temperature, theta_sat, theta_rel):
        """Get the contribution of the configuration entropy to the total energy of desorption.
        Args:
            temperature (list): temperature range
            theta_sat (float): saturation coverage of TPD
            theta_rel (list): relative coverage of TPD

        Returns:
            E_config (list): contribution of the configuration entropy to the total entropy
        """
        E_config = []
        for i in range(len(temperature)):
            Ec = -  kB * temperature[i] * np.log(theta_sat*theta_rel[i] / ( 1 - theta_sat*theta_rel[i]))
            E_config.append(Ec)
        return E_config
    
    def _get_ads_ads_interaction(self, temperature, theta_sat, theta_rel, b):
        """Get the contribution of the saturation coverage to the total energy of desorption.
        Args:
            temperature (list): temperature range
            theta_sat (float): saturation coverage of TPD
            theta_rel (list): relative coverage of TPD
            b (float): slope of the saturation coverage

        Returns:
            E_b (list) : contribution of the saturation coverage to the total entropy
        """
        E_b = []
        for i in range(len(temperature)):
            E_b.append(-b * theta_rel[i] * theta_sat)
        return E_b

    def _fit_Ed_theta(self, temperature, E_0, b, theta_sat, theta_rel):
        """Fit the desorption energy to the relative coverage
        Fed into scipy curve fit

        Args:
            temperature (list): temperature range
            E_0 (float): energy at zero coverage
            b (float): interaction parameter
            theta_sat (float): saturation coverage of TPD 
            theta_rel (float): relative coverage 

        Returns:
            list: Desorption energy based on fit
        """
        Ed_all = []
        for i in range(len(temperature)):
            Ed = E_0 \
            -  kB * temperature[i] * np.log(theta_sat*theta_rel[i] / ( 1 - theta_sat*theta_rel[i]))
            -  b * theta_rel[i] * theta_sat
            Ed_all.append(Ed)
        return Ed_all

    def _fit_quad_Ed_theta(self, temperature, E_0, b1, b2, theta_sat, theta_rel):
        """ Fit the Ed vs theta relationship using the configurational entropy
        and the ads-ads interaction which has a quadratic dependence.

        Args:
            temperature (list): temperature range
            E_0 (float): energy at zero coverage
            b (float): interaction parameter
            theta_sat (float): saturation coverage of TPD 
            theta_rel (float): relative coverage 

        Returns:
            list: Desorption energy based on fit
        """
        Ed_all = []
        for i in range(len(temperature)):
            Ed = E_0 \
            -  kB * temperature[i] * np.log(theta_sat*theta_rel[i] / ( 1 - theta_sat*theta_rel[i]))
            -  b1 * theta_rel[i] * theta_sat
            -  b2 * (theta_rel[i] * theta_sat)**2
            Ed_all.append(Ed)
        return Ed_all

    def _eq_coverage_function(self, theta, T, G0, b, p):
        """Function to implicitly solve the equilibrium coverage

        Args:
            theta (float): Guessed coverage 
            T (float) : temperature
            G0 (float): Free energy at the half a mono-layer coverage
            b (float): Interaction parameter
            p (float): partial pressure of CO
        """
        kBT = kB * T
        ## start by calculating the equilibirum constant 
        K = np.exp( -1 * ( G0 + b * ( theta - 1./2. ) ) / kBT ) 
        return theta - ( K / ( 1 + K ) )

    def _jacobian(self, theta, T, G0, b, p):
        """Jacobian function for finding the root 

        Args:
            theta (list): Guessed coverage
            T ([type]): [description]
            G0 (float): Free energy at the half a mono-layer coverage
            b (float): Interaction parameter
            p (float): partial pressure of CO
        """
        kBT = kB * T
        ## start by calculating the equilibirum constant 
        K = np.exp( -1 * ( G0 + b * ( theta - 1./2. ) ) / kBT ) 

        return 1 + K / (1+K)**2 * b / kBT 

    def _get_equilibirum_coverage(self, E0, b):
        """Equilibirum coverage based on equilibrium constant that is coverage dependent

        Args:
            E0 (float): Desorption energy at zero coverage
            b (float): Interaction parameter

        Returns:
            list: equilibrium coverage and free energy of CO adsorption
        """

        theta_eq = []
        dG_eq = []

        for index, T in enumerate(self.plot_temperature):
            entropy_ads = self.thermo_ads.get_entropy(temperature=T, verbose=False)
            entropy_gas = self.thermo_gas.get_entropy(temperature=T, \
                                                          pressure=self.p, verbose=False)

            # converting from energies to free energies
            entropy_difference =  entropy_ads - entropy_gas
            partial_co = self.p / 101325.
            # convert desorption energy into adsorption energy
            dG0 = -1 * E0 -1 * T * entropy_difference
            K_guess = np.exp( -1 * dG0 / kB / T )
            theta_guess = K_guess / ( 1 + K_guess ) 

            try:
                theta = newton(
                    func = lambda x: self._eq_coverage_function(x, T, dG0, b, partial_co ), 
                    fprime = lambda x: self._jacobian(x, T, dG0, b, partial_co ),
                    x0=theta_guess, 
                )
            except RuntimeError:
                theta = 0


            dG = ( -dG0 + b * ( theta - 1./2. )  )
            theta_eq.append(theta)
            dG_eq.append(dG)

        return theta_eq, dG_eq

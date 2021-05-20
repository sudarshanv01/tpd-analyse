import numpy as np
from glob import glob
from useful_functions import AutoVivification, get_vibrational_energy
from pprint import pprint
import os, sys
from scipy.optimize import curve_fit, least_squares, minimize
from matplotlib import cm
from tpd_analyse.tools.parser_function import get_stable_site_vibrations, get_gas_vibrations, \
                            get_coverage_details, diff_energies, \
                            get_lowest_absolute_energies,\
                            get_differential_energy,\
                            accept_states, \
                            get_constants, \
                            stylistic_comp, stylistic_exp

from tpd_analyse.tools.parser_class import ParseInfo, experimentalTPD
import mpmath as mp
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from ase.io import read
from ase.db import connect
import matplotlib
import csv
from ase.units import kB
from mpmath import mpf
from scipy import optimize
from matplotlib.ticker import FormatStrFormatter
from matplotlib import cm
from scipy.optimize import newton
import matplotlib.pyplot as plt

class PlotTPD():

    def __init__(self, exp_data, order, constants,
                    thermo_ads, thermo_gas, correct_background=True, bounds=[], plot_temperature=np.linspace(100,400),
                    p=101325, initial_guess_theta=0.5):

        """Perform the temperature programmed desorption analysis for a surface 
        based on configurational entropy an interaction parameters and a zero coverage
        energy term

        Args:
            exp_data (list): globbed files with csv
            order (int): Order of the reaction
            constants (list): Parameters to parse TPD
            color_map (obj): Color map for different exposures
            thermo_ads (obj): HarmonicThermo for the adsorbate
            thermo_gas (obj): IdealGasThermo for the gas
            bounds (list, optional): Bounds within to fit the coverage of the TPD. Defaults to [].
            plot_temperature (array, optional): Temperature range to plot the equilbirum coverage. Defaults to np.linspace(100,400).
            p (float, optional): Pressure of gas molecule. Defaults to 101325.
        """

        # Define all the __init__ variables
        self.constants = constants # Constants for TPD parsing
        self.exp_data = exp_data # globbed files with csv
        self.order = order
        self.thermo_ads = thermo_ads # adsorbate HarmonicThermo
        self.thermo_gas = thermo_gas # gas phase IdealGasThermo
        self.plot_temperature = plot_temperature # temperature for plotting equilibrium coverages
        self.p = p # partial pressure of gas
        self.bounds = bounds
        self.correct_background = correct_background
        self.initial_guess_theta = initial_guess_theta

        # Results
        self.norm_results = AutoVivification()
        self.results = AutoVivification()
        self.Ed = AutoVivification()
        self.theta_rel = AutoVivification()
        self.theta_eq = AutoVivification()
        self.theta_eq_p = AutoVivification()
        self.theta_eq_n = AutoVivification()
        self.E0 = AutoVivification()
        self.Ed_fitted = AutoVivification()
        self.b = AutoVivification()
        self.theta_sat = AutoVivification()
        self.dG = AutoVivification()
        self.error = AutoVivification() # Error in Ed fit
        self.temperature_range = {} # for plotting in main figure

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
        T_switch, T_max, T_rate_min, beta = self.constants

        # Create the temperature range based on the switch data
        assert len(T_switch) < 3 , 'Only implemented switch lower than 2 '
        if len(T_switch) == 1:
            temperature_ranges = [
                                    [0, T_switch[0]], ## Terrace
                                    [T_switch[0], T_max] ## Step
                                ]
        elif len(T_switch) == 2:
            temperature_ranges = [
                [0, T_switch[0]], # Terrace 1
                [T_switch[0], T_switch[1]], # Terrace 2
                [T_switch[1], T_max] # step
            ]

        # range of temperatures for different TPD values
        self.temperature_range = temperature_ranges

        # 1. Get the TPD results which includes background subtraction
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
                indices = [ a for a in range(len(self.norm_results[exposure].temperature)) \
                        if T_range[0] < self.norm_results[exposure].temperature[a] < T_range[1] ] 

                self.results[surface_index][exposure]['temperature'] = self.norm_results[exposure].temperature[indices]
                self.results[surface_index][exposure]['normalized_rate'] = self.norm_results[exposure].normalized_rate[indices]

                # create variable within look for easy calling 
                temperatures = self.results[surface_index][exposure]['temperature']
                rates = self.results[surface_index][exposure]['normalized_rate']

                #2. For each point get the energy of desorption as a function of the coverage
                data = self._Ed_temp_dependent(
                            temperature=temperatures, 
                            rate=rates,
                            beta=beta,
                            )

                ## if there are nans
                args_accept = [i for i in range(len(data[0])) \
                                    if np.isfinite(data[0][i]) and \
                                    data[1][i] > 0]

                self.Ed[surface_index][exposure], self.theta_rel[surface_index][exposure] = data
                self.Ed[surface_index][exposure] = self.Ed[surface_index][exposure][args_accept]
                self.theta_rel[surface_index][exposure] = self.theta_rel[surface_index][exposure][args_accept]
                temperature_fit = self.norm_results[exposure].temperature[indices][args_accept]
                self.results[surface_index][exposure]['temperature_fit'] = temperature_fit

                ## TODO: Is bounds needed in code?
                # # check if relative theta has bounds
                # if self.bounds:
                #     index_bounds = [i for i in range(len(self.theta_rel[surface_index][exposure])) \
                #              if min(self.bounds) <=  self.theta_rel[surface_index][exposure][i] <= max(self.bounds)]
                #     self.Ed[surface_index][exposure] = self.Ed[surface_index][exposure][index_bounds]
                #     self.theta_rel[surface_index][exposure] = self.theta_rel[surface_index][exposure][index_bounds]
                #     temperature_fit = temperature_fit[index_bounds]
                    
                # 3. Fit the Desorption energy curve to the desorption energy equation
                ## First get good initial guesses for parameters
                guess_E0 = self.Ed[surface_index][exposure][0]
                linear_region = [i for i in range(len(self.Ed[surface_index][exposure])) if 0.25 < self.Ed[surface_index][exposure][i] < 0.75  ]
                guess_b = 0.1#-1 * np.polyfit(self.theta_rel[surface_index][exposure][linear_region], self.Ed[surface_index][exposure][linear_region], 1)[0]

                popt, pcov = curve_fit(\
                lambda temp, E0, b, theta_sat: self._fit_Ed_theta(temp, E0, b, theta_sat, \
                                            self.theta_rel[surface_index][exposure]), \
                                            xdata=temperature_fit,
                                            ydata=self.Ed[surface_index][exposure],
                                            p0=[guess_E0, guess_b, self.initial_guess_theta], 
                                                )

                # # Least squares routine
                # TODO better minimisation routine?
                # res = minimize( 
                #     self._least_sq_Ed_theta, 
                #     x0=[guess_E0, guess_b, 0.8], 
                #     args=(temperature_fit, self.theta_rel[surface_index][exposure], self.Ed[surface_index][exposure]),
                #     bounds=bounds,
                #     # options={'ftol': 1e-2 }
                #     tol=1e-6,
                # ) 
                # popt = res.x
                # print(popt, guess_b)
                # print(res.message)
                # print(self._least_sq_Ed_theta(res.x, temperature_fit, self.theta_rel[surface_index][exposure], self.Ed[surface_index][exposure] ))

                residual = self._least_sq_Ed_theta(popt, temperature=temperature_fit, theta_rel=self.theta_rel[surface_index][exposure], \
                    Ed_real=self.Ed[surface_index][exposure])

                self.E0[surface_index][exposure], self.b[surface_index][exposure], self.theta_sat[surface_index][exposure], \
                                = popt
                self.error[surface_index][exposure] = residual#np.sqrt(np.diag(pcov)) #self._least_sq_Ed_theta(res.x, temperature_fit, self.theta_rel[surface_index][exposure], self.Ed[surface_index][exposure] ) #np.sqrt(np.diag(pcov))
                self.Ed_fitted[surface_index][exposure] = self._fit_Ed_theta(temperature_fit, \
                                                *popt, self.theta_rel[surface_index][exposure] )
                # 4. Calculate the coverage at equilbirum
                self.theta_eq[surface_index][exposure], self.dG[surface_index][exposure]\
                 = self._get_equilibirum_coverage(
                            E0=self.E0[surface_index][exposure], 
                            b=self.b[surface_index][exposure],
                            )
                self.theta_eq_p[surface_index][exposure], _ \
                 = self._get_equilibirum_coverage(
                            E0=self.E0[surface_index][exposure]+self.error[surface_index][exposure], 
                            b=self.b[surface_index][exposure],
                            )
                self.theta_eq_n[surface_index][exposure], _ \
                 = self._get_equilibirum_coverage(
                            E0=self.E0[surface_index][exposure]-self.error[surface_index][exposure],
                            b=self.b[surface_index][exposure],
                            )
                
                ######## PLOTS ##########

                # total_coverage = self.theta_sat[surface_index][exposure] * self.theta_rel[surface_index][exposure]
                # interaction_term = - self.b[surface_index][exposure] * self.theta_sat[surface_index][exposure] * self.theta_rel[surface_index][exposure]
                # config_term = - kB * temperature_fit * np.log(self.theta_sat[surface_index][exposure]*self.theta_rel[surface_index][exposure] / ( 1 - self.theta_sat[surface_index][exposure]*self.theta_rel[surface_index][exposure]))




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
        # square_error = np.square([np.abs(a) for a in residual])
        # mean_sq_error = np.mean(square_error)
        # rms = np.sqrt(mean_sq_error) 
        return mea #square_error #mea#Ed_real - Ed_fit

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
            # K = np.exp( - dG / kB / T)
            # theta =  K * partial_co / (1 + K * partial_co)
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

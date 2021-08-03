import numpy as np
import sys
import os
import csv
from scipy.optimize import curve_fit, minimize, least_squares
from copy import deepcopy

class experimentalTPD:
    def __init__(self, tpd_filename, temprange, tempmin, beta, order, correct_background=True):
        """
        Class that gives the normalised rate and temperature data

        ...
        
        Attributes 
        ----------

        tpd_filename : str 
            a filename which has the tpd data stored with comma separated values 
        beta : float 
            heating rate 
        temprange : list 
            where the simulation of the TPD has to be done 
        tempmin : list 
            where the rate is close to zero, gives a baseline 
        order : int 
            order of the reaction 
        
        exp_temperature : list 
            temperatures read in from the comma separated file in units of K 
        exp_rates : list 
            rates read in from the comma separated file in arbitrary units 
        normalized_rate : list 
            rates after removing background pumping effects 
        temperature : list 
            temperatures for the normalised rates
        popt_b : list 
            fit parameters obtained from scipy.optimize for subtraction from pumping function

        Methods 
        -------

        collect_tpd_data()
            reads in tpd data based on np.genfromtxt and stores rates and temperatures from exp
        get_normalized_data()
            calls _normalize_TPD_baseline to get the parameters of the pumping function fit 
        _normalize_TPD_baseline(kwargs)
            fit an exponential function using curvefit to the "tail" of the TPD 
            subtracts all rates with the exponential function 

            
        """ 
        self.tpd_filename = tpd_filename
        self.beta = beta # heating rate
        self.temprange = temprange # Range where the actual TPD simulation needs to be done
        self.tempmin = tempmin # Range where the rate is close to zero
        self.order = order
        self.exp_temperature = []
        self.exp_rate = []
        self.normalized_rate = []
        self.temperature = []
        self.popt_b = [] # fit for the subtraction function
        self.background_correction = [] # background correction
        self.correct_background = correct_background

        # Consistency checks 
        if self.order != 1:
            raise NotImplementedError

        # Progression of commands
        self.collect_tpd_data()
        self.get_normalized_data()

    # Collects the data which has been taken from a publication
    def collect_tpd_data(self):
        # First get the tpd data frm the file created by WebPlotDigitizer
        text_max_exposure = np.genfromtxt(self.tpd_filename,  delimiter=',')
        # Get the temperatures and rates from the experiment
        self.exp_temperature = text_max_exposure[:,0]
        self.exp_rate = text_max_exposure[:,1] 
        # subtract initial 
        self.exp_rate = self.exp_rate 


    # Normalize the TPD data by subtracting the ends
    def get_normalized_data(self):
        ## Pass all temperature ranges to the normalisation function
        kwargs = {'min_range':self.tempmin, 'val_range':self.temprange}
        corrected_rates = self._normalize_TPD_baseline(kwargs)
        # Return the background normalised rate
        self.normalized_rate = corrected_rates #normalized_data['normalized_rate']
        self.temperature = self.exp_temperature#np.array(normalized_data['temperature'])


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

    def _fit_rate_with_constraints(self, x, temperature_tail, rate_tail, j0, T0):
        """Minimise least sqaures of residual to make sure that
        the last point of the signal matches the real signal and 
        the exponent fits the tail of the TPD spectra

        :param temperature: [description]
        :type temperature: [type]
        :param a: [description]
        :type a: [type]
        :param k: [description]
        :type k: [type]
        :param j0: [description]
        :type j0: [type]
        :param T0:
        :type T0:
        """
        a, k = x
        temperature_tail = np.array(temperature_tail)
        ## First fit an exponent to the tail of the TPD
        rate_fit_tail = a * np.exp(-k * temperature_tail)        
        ## Next ensure that the last point is set to the j0 value
        j0_fit = a * np.exp(-k * T0)
        ## determine an error fit 
        error_tail = rate_fit_tail - rate_tail
        error_final = [ 0.01 * ( j0_fit - j0 ) ]
        ## Total error is the sum of lists of error
        error = error_tail #+ error_final

        lsq_residual = np.sum(np.square(error))       
        
        return lsq_residual
    

    def _normalize_TPD_baseline(self, kwargs):
        """Normalises the TPD with the required correction for
        the background signal coming from CO adsorbed on the walls 
        of the TPD apparatus. This effect is modelled using an extra 
        exponential term that decays with a different decay length

        :param kwargs: ranges of the tail and of the TPD
        :type kwargs: dict
        :return: normalised TPD data
        :rtype: dict
        """
        T_min_end, T_max_end = kwargs['min_range']
        rates_end = []
        T_end = []
        for i in range(len(self.exp_temperature)):
            if T_min_end < self.exp_temperature[i] < T_max_end:
                rates_end.append(self.exp_rate[i])
                T_end.append(self.exp_temperature[i])
            
        ## Normalise the TPD by taking the exponent from the tail 
        ## but by making sure that the last point is the same as the
        ## actual TPD graph
        if self.correct_background:
            x0 = [self.exp_rate[0], 0.01]
            min_fit = minimize(self._fit_rate_with_constraints, x0, \
                    args=(T_end, rates_end, self.exp_rate[0], self.exp_temperature[0]),\
                    method='Nelder-Mead',
                    )
        
            background = self._exponential_fit(self.exp_temperature, *min_fit.x)
        else:
            background = len(self.exp_temperature)*[self.exp_rate[0]]
        self.background_correction = background #no_background + background
        self.corrected_TPD = self.exp_rate - self.background_correction


        return self.corrected_TPD


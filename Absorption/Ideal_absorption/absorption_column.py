import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp

"""In this script a simple adsorption column is designed. The script will use Henry's law to calcualte the vapor-liquid equilibrium
is thus in practice only valid for low gas concentrations in water solvent phase."""


class AbsorptionColumn:

    def __init__(self, T, P, standard_henry_constant, henry_temperature_dependence, gas_inlet_concentration, 
                 gas_outlet_threshold_concentration, liquid_inlet_concentration, molecular_weight_gas, gas_inflow):
        self.T = T
        self.P = P
        self.standard_henry_constant = standard_henry_constant
        self.henry_temperature_dependence = henry_temperature_dependence
        self.gas_inlet_concentration = gas_inlet_concentration
        self.gas_outlet_threshold_concentration = gas_outlet_threshold_concentration
        self.liquid_inlet_concentration = liquid_inlet_concentration
        self.molecular_weight_gas = molecular_weight_gas
        self.gas_inflow = gas_inflow
        

    def calculate_henry_constant(T, standard_henry_constant, henry_temperature_dependence):

        """
        Calculate the henry constant for a given temperature

        Parameters:
            T (float) : Temperature in Kelvin
            standard_henry_constant (float) : Standard henry constant
            henry_temperature_dependence (float) : Temperature dependence of the henry constant

        Returns:
            henry_constant (float) : Henry constant at the given temperature
        """

        # Calculate the henry constant for a given temperature 
        henry_constant = standard_henry_constant * np.exp(henry_temperature_dependence * (1 / T - 1 / 298.15))
        return henry_constant

    def calculate_minimum_solvent_flow(henry_constant, gas_inlet_concentration, gas_outlet_threshold_concentration, 
                                       liquid_inlet_concentration, molecular_weight_gas, gas_inflow, P):
        """
        Calculate the minimum solvent flow from the overall mass balance under the assumption 
        that the liquid is saturated at the bottom of the column

        Parameters:
            henry_constant (float) : Henry constant at the given temperature
            gas_inlet_concentration (float) : Concentration of gas in the inlet stream
            gas_outlet_threshold_concentration (float) : Concentration of gas in the outlet stream
            liquid_inlet_concentration (float) : Concentration of gas in the inlet stream
            molecular_weight_gas (float) : Molecular weight of the gas
            gas_inflow (float) : Inflow of the gas in m3/h
            P (float) : Pressure in Pascals

        Returns:
            minimum_solvent_flow (float) : Minimum solvent flow required to reach the target outlet concentration
        """
        # Calculate the saturated liquid outlet concentration using Henry's law 
        liquid_maximum_outlet_concentration = P * henry_constant * molecular_weight_gas * 1000
        
        # Calculate the vapor-liquid equilibrium coefficient
        K = gas_inlet_concentration / liquid_maximum_outlet_concentration

        # Calculate the minimum solvent flow from the overall mass balance under the assumption 
        # that the liquid is saturated at the bottom of the column
        minimum_solvent_flow = gas_inflow*((gas_inlet_concentration - gas_outlet_threshold_concentration) 
                                           / (gas_inlet_concentration / K - liquid_inlet_concentration))
        return minimum_solvent_flow
    
    def calculate_molar_ratio(molar_fraction):
        """
        Calculate the molar ratio of a given molar fraction

        Parameters:
            molar_fraction (float) : Molar fraction of the gas in the liquid

        Returns:
            molar_ratio (float) : Molar ratio of the gas in the liquid
        """
        molar_ratio = molar_fraction / (1 - molar_fraction)
        return molar_ratio
    
    # Rest is going to be added soon. In the coming section the mccabe-thiele method for absorption is going to be elaborated.
    


"""INPUT PARAMETERS FROM USER"""

"""
Parameter                           :   type    :       units       :   Description
T                                   :   float   :       Kelvin      :   Temperature in Kelvin
P                                   :   float   :        bar        :   Pressure in Pascals
standard_henry_constant             :   float   :     mol/kg*bar    :   Standard henry constant
henry_temperature_dependence        :   float   :       Kelvin      :   Temperature dependence of the henry constant
gas_inlet_concentration             :   float   :        ppm        :   Gas inlet concentration
gas_outlet_threshold_concentration  :   float   :        ppm        :   Gas outlet threshold concentration
liquid_inlet_concentration          :   float   :        ppm        :   Liquid inlet concentration
molecular_weight_gas                :   float   :       g/mol       :   Molecular weight of the gas
"""

# Input parameters
T = 298.15
P = 1
standard_henry_constant = 0.083
henry_temperature_dependence = 2100
gas_inlet_concentration = 20
gas_outlet_threshold_concentration = 2 
liquid_inlet_concentration = 0
molecular_weight_gas = 34.082 # for H2S
gas_inflow = 100 # m3/h

# Run the script


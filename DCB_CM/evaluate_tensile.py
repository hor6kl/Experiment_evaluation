
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import math

from common_evaluation import *

class Evaluate_Tensile():
    def __init__(self,force: [], displacement: [], spec_data: 'class', names_par: dict , path_fig: str):
        self.force = force
        self.displacement = displacement
        self.spec_data = spec_data
        self.names_par = names_par
        self.path_fig = path_fig


    def Return_dimensions(self, name_dim: str) -> []:
    # Function returns list of dimensions specified by name in dictionary 
        name_index = self.names_par[name_dim]

        num_index = [self.spec_data.data_name.index(i) for i in name_index ]
        values = [self.spec_data.data_value[i]*self.spec_data.data_dimension[i] for i in num_index ]

        return values


    def Calculate_stress(self) -> []:
    # Function returns stress computed from surface and force
        
        Height_values = self.Return_dimensions('height') 
        Width_values = self.Return_dimensions('width')

        height = min(Height_values)
        width  = min(Width_values)

        surface = height*width

        sigma = [f/surface for f in self.force]

        return sigma

        
    def Calculate_strain(self) -> []:
    # Function returns strain computed from gauge length and displacement 
        
        Length_values = self.Return_dimensions('length') 

        length = min(Length_values)

        strain = [(d)/(length) for d in self.displacement]

        return strain

    def find_close_value_index(self, lst, x, range_val = None):
    # Find the index of a value in the list that is within a small range (range_val) of x.
        for i, y in enumerate(lst):
            if abs(x - y) <= range_val:
                return i
        print(f"Index was not found! Check your input.")

        return -1

    def Evaluate_Young_modulus(self, strain_arg = None) -> float:
    #Funciton computes Young modulus from stress strain curve
        
        # These are strain limits where young modulus should be evaluated according to ISO 527_1
        strain_lim_1 = 0.0005
        strain_lim_2 = 0.0025

        stress = self.Calculate_stress()

        # strain_arg was added to support direct pass of strain from UTM
        # The strain from UTM was displacement functionality was left in code
        if strain_arg == None:
            strain = self.Calculate_strain()
        else:
            strain = strain_arg

        print(f"Maximum stress: {max(stress):.4e}")
        print(f"Maximum strain: {max(strain):.5f}")

        #print(f"strain: {strain}")

        # Range value to find index is done in two iterations to consider different strain steps during experiment
        delta_in = round(len(strain)/3)
        delta = abs(strain[delta_in] - strain[delta_in-1])*2
        index_aprox_1 = self.find_close_value_index(strain, strain_lim_1, delta)
        index_aprox_2 = self.find_close_value_index(strain, strain_lim_1, delta)
        delta_1 = abs(strain[index_aprox_1] - strain[index_aprox_1-1])
        delta_2 = abs(strain[index_aprox_2] - strain[index_aprox_2-1])

        # Finds index of limits
        index_lim_1 = self.find_close_value_index(strain, strain_lim_1, delta_1)
        index_lim_2 = self.find_close_value_index(strain, strain_lim_2, delta_2)
        if index_lim_2 == -1:
            print(f" index_lim_2 was not found and range for Young moduls evaluation is from: {index_lim_1} to end of list")

        # Modification of dimension to correspond with evaluation limits
        stress = stress[index_lim_1:index_lim_2]
        strain = strain[index_lim_1:index_lim_2]

        # This coefficient will not be the same with polyfit. for same resaults domain=[]
        c0, c1 = np.polynomial.polynomial.Polynomial.fit(strain,stress, 1, domain=[])

        print(f"Young modulues E = {c1:.4e}")
    
        delta = -c0/c1

        x_axis = np.linspace(min(strain), max(strain), 200)
        y_fitted_long = c0 + x_axis*c1 

        fig = plt.figure()
        plt.plot(x_axis, y_fitted_long)
        plt.plot(strain, stress, 'o')
    
        plt.legend(["PolyFit", "Experiments"])
        plt.ylabel('$Stress$ [N/m]', fontsize=12)
        plt.xlabel('$Strain$ [m]', fontsize=12)
        plt.title('Young modulus', fontsize=14)
        plt.grid(True)
    
        path_fig_png = self.path_fig + self.spec_data.data_value[0] + "_E_polyfit" + ".png"
    
        fig.savefig(path_fig_png)

        return c1 


    def Evaluate_Poisson(self, delta_x: []) -> float:
    #Funciton computes poisson ratio from biaxial strain in x and y direction 
        
        # These are strain limits where young modulus should be evaluated according to ISO 527_1
        strain_lim_1 = 0.0005
        strain_lim_2 = 0.0025

        strain_y = self.Calculate_strain()

        Width_values = self.Return_dimensions('width') 
        length = min(Width_values)
        strain_x = [(d)/(length) for d in delta_x]

        # Range value to find index
        delta = abs(strain_y[round(len(strain_y)/10)] - strain_y[round(len(strain_y)/10)-1])/2

        # Finds index of limits
        index_lim_1 = self.find_close_value_index(strain_y, strain_lim_1, delta)
        index_lim_2 = self.find_close_value_index(strain_y, strain_lim_2, delta)

        # Modification of dimension to correspond with evaluation limits
        force_y = self.force[index_lim_1:index_lim_2]
        strain_y = strain_y[index_lim_1:index_lim_2]
        strain_x = strain_x[index_lim_1:index_lim_2]

        # This coefficient will not be the same with polyfit. for same results domain=[]
        c0_x, c1_x = np.polynomial.polynomial.Polynomial.fit(force_y, strain_x, 1, domain=[])
        c0_y, c1_y = np.polynomial.polynomial.Polynomial.fit(force_y, strain_y, 1, domain=[])

        poisson = c1_x/c1_y

        print(f"Poisson ratio nu = {poisson:.5f}")
        x_axis = np.linspace(min(force_y), max(force_y), 200)
        y_fitted_long_x = c0_x + x_axis*c1_x
        y_fitted_long_y = c0_y + x_axis*c1_y

        fig = plt.figure()
        plt.plot(x_axis, y_fitted_long_x)
        plt.plot(x_axis, y_fitted_long_y)
        plt.plot(force_y, strain_x, 'o')
        plt.plot(force_y, strain_y, 'x')
    
        plt.legend(["PolyFit", "Experiments"])
        plt.ylabel('$Force$ [N/m]', fontsize=12)
        plt.xlabel('$Strain$ [m]', fontsize=12)
        plt.title('Poisson', fontsize=14)
        plt.grid(True)
    
        path_fig_png = self.path_fig + self.spec_data.data_value[0] + "_nu_polyfit" + ".png"
    
        fig.savefig(path_fig_png)

        return poisson 
        
# ----------------------------------------

# Graphical evaluation
def graph_stress_strain_curves(X: [[float]], Y: [[float]], name_TAR, color_set, marker_set, line_set, path_fig):
    # Function creates graph of stress strain curve with average

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot all curves
    for cons in range(0, len(X)):
        plt.plot(X[cons], Y[cons], color=color_set[cons], marker=marker_set[cons], linestyle=line_set[cons])

    avg_x, avg_y = average_curve(X, Y)

    plt.plot(avg_x, avg_y, color='green')
    
    plt.xlabel(r'Strain $\epsilon$ [-]')
    plt.ylabel(r'Stress $\sigma$ [Pa]')

    plt.legend(name_TAR)

    path_fig_png = path_fig + ".png"
    path_fig_eps = path_fig + ".eps"

    fig.savefig(path_fig_eps)
    fig.savefig(path_fig_png)





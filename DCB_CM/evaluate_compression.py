import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import math

from common_evaluation import *


class Evaluate_Tensile:
    def __init__(
        self,
        force: [],
        displacement: [],
        spec_data: "class",
        names_par: dict,
        path_fig: str,
    ):
        self.force = force
        self.displacement = displacement
        self.spec_data = spec_data
        self.names_par = names_par
        self.path_fig = path_fig

    def Return_dimensions(self, name_dim: str) -> []:
        # Function returns list of dimensions specified by name in dictionary
        name_index = self.names_par[name_dim]

        num_index = [self.spec_data.data_name.index(i) for i in name_index]
        values = [
            self.spec_data.data_value[i] * self.spec_data.data_dimension[i]
            for i in num_index
        ]

        return values

    def return_surface(self) -> float:

        Width_values = self.Return_dimensions("width")
        width = np.mean(Width_values)

        return np.pi * width**2 / 4

    def Calculate_stress(self) -> []:
        # Function returns stress computed from surface and force
        surface = self.return_surface()

        sigma = [f / surface for f in self.force]

        return sigma

    def Calculate_strain(self) -> []:
        # Function returns strain computed from gauge length and displacement

        Height_values = self.Return_dimensions("height")

        height = min(Height_values)

        strain = [(d) / (height) for d in self.displacement]

        return strain

    def return_max_stress(self, range: [] = None) -> float:

        min_index = find_close_value_index(self.displacement, range[0])
        max_index = find_close_value_index(self.displacement, range[1])
        stress = self.Calculate_stress()
        if min_index == -1:
            min_index = 0
        max_stress = max(stress[min_index:max_index])

        return max_stress

    def return_max_strain(self) -> float:
        return max(self.Calculate_strain())

    def Evaluate_Young_modulus(self, strain_arg=None) -> float:
        # Funciton computes Young modulus from stress strain curve

        # These are strain limits where young modulus should be evaluated according to ISO 527_1
        strain_lim_1 = 0.0005
        strain_lim_2 = 0.0025

        stress = self.Calculate_stress()

        # strain_arg was added to support direct pass of strain from UTM
        # The strain from UTM was displacement functionality was left in code
        if strain_arg is None:
            strain = self.Calculate_strain()
        else:
            strain = strain_arg

        #        min_index = find_close_value_index(self.displacement, 0.0000)
        #        max_index = find_close_value_index(self.displacement, 0.0034)
        #        max_stress = max(stress[min_index:max_index])
        #        print(f"Maximum stress: {max_stress:.4e}")

        # print(f"strain: {strain}")

        # Finds index of limits
        index_lim_1 = find_close_value_index(strain, strain_lim_1)
        index_lim_2 = find_close_value_index(strain, strain_lim_2)
        if index_lim_2 == -1:
            print(
                f" index_lim_2 was not found and range for Young moduls evaluation is from: {index_lim_1} to end of list"
            )

        # Modification of dimension to correspond with evaluation limits
        stress = stress[index_lim_1:index_lim_2]
        strain = strain[index_lim_1:index_lim_2]

        # This coefficient will not be the same with polyfit. for same resaults domain=[]
        c0, c1 = np.polynomial.polynomial.Polynomial.fit(strain, stress, 1, domain=[])

        print(f"Young modulues E = {c1:.4e}")

        delta = -c0 / c1

        x_axis = np.linspace(min(strain), max(strain), 200)
        y_fitted_long = c0 + x_axis * c1

        fig = plt.figure()
        plt.plot(x_axis, y_fitted_long)
        plt.plot(strain, stress, "o")

        plt.legend(["PolyFit", "Experiments"])
        plt.ylabel("$Stress$ [N/m]", fontsize=12)
        plt.xlabel("$Strain$ [m]", fontsize=12)
        plt.title("Young modulus", fontsize=14)
        plt.grid(True)

        path_fig_png = (
            self.path_fig + self.spec_data.data_value[0] + "_E_polyfit" + ".png"
        )

        fig.savefig(path_fig_png)

        return c1


# ----------------------------------------


# Graphical evaluation
def graph_stress_strain_curves(
    X: [[float]], Y: [[float]], name_TAR, color_set, marker_set, line_set, path_fig
):
    # Function creates graph of stress strain curve with average

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot all curves
    for cons in range(0, len(X)):
        plt.plot(
            X[cons],
            Y[cons],
            color=color_set[cons],
            marker=marker_set[cons],
            linestyle=line_set[cons],
        )

    avg_x, avg_y = average_curve(X, Y)

    plt.plot(avg_x, avg_y, color="green")

    plt.xlabel(r"Strain $\epsilon$ [-]")
    plt.ylabel(r"Stress $\sigma$ [Pa]")

    plt.legend(name_TAR)
    plt.grid

    path_fig_png = path_fig + ".png"
    path_fig_eps = path_fig + ".eps"

    fig.savefig(path_fig_eps)
    fig.savefig(path_fig_png)

    plt.close

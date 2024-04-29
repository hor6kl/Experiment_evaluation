
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import math

import logging
# Create a logger
logger = logging.getLogger(__name__)
# Configure the logger to write to a file
logging.basicConfig(filename='DCB_evaluation.log', encoding='utf-8', level=logging.INFO)

from common_evaluation import *
common_functions = Common_functions()

class Evaluate_DCB():
    def __init__(self, force: [], displacement: [], time_utm: [], time_delam: [], delam_length: [], spec_data: 'class', names_par: dict, name: str ):
        self.force = force
        self.displacement = displacement
        self.time_utm = time_utm
        self.time_delam = time_delam
        self.delam_length = delam_length
        self.spec_data = spec_data
        self.names_par = names_par
        self.name = name
        # add initial delamination to all values in delam_length vector
        self.delamination_expand()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~
    # possible common functions
    #
    def Return_dimensions(self, name_dim: str) -> []:
    # Function returns list of dimensions specified by name in dictionary 
        name_index = self.names_par[name_dim]

        num_index = [self.spec_data.data_name.index(i) for i in name_index ]
        values = [self.spec_data.data_value[i]*self.spec_data.data_dimension[i] for i in num_index ]

        return values
    #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~
    
    def delamination_expand(self):
        Delam_values = self.Return_dimensions('delam')
        delam_a0 = min(Delam_values)
    
        # add intitial delamination to values
        Del_length_compl = [float(delam_a0) + float(x) for x in self.delam_length]
        self.delam_length = Del_length_compl

    # ----------------------------------------

    def vypocet_energie_DCB(self, F, a, epsilon_n):
        # Wrong not validated delta does not equal 0
    
        B_index = a.data_name.index('B')
        ai_index = a.data_name.index('ai')
    
        # sirka vzorku
        B = float(a.data_value[B_index])*float(a.data_dimension[B_index])
        # non ahezive od osy zatizeni
        a_del = float(a.data_value[ai_index])*float(a.data_dimension[ai_index])
        delta = 0
    
        F_max = max(F)
        index_F_max = (np.argmax(F))
    
        float_lst = [float(x) for x in epsilon_n]
        deltaC = float_lst[index_F_max]  # posunuti hodnoty v bode F_max v [m]
    
        defor_energ = (3*F_max*deltaC)/(2*B*(a_del+delta))
    
        return defor_energ
    
    # ----------------------------------------
    
    def vypocet_energie_DCB_R_MBT(self, Force: float, delam_R: float, Displ: float, delta: float) -> float:
    
        Width_values = self.Return_dimensions('width')
        B  = min(Width_values)
        
        # delamination length
        a_del = delam_R
    
        deltaC = float(Displ)
        
        defor_energ = (3*Force*deltaC)/(2*B*(a_del+delta))
    
        return defor_energ
    
    # ----------------------------------------
    
    def vypocet_energie_DCB_R_CC(self, Force: float, delam_R: float, Displ: float, n: float) -> float:
    
        Width_values = self.Return_dimensions('width')
        B  = min(Width_values)
        
        # delamination length
        a_del = delam_R
    
        deltaC = float(Displ)
    
        defor_energ = (n*Force*deltaC)/(2*B*a_del)
    
        return defor_energ
    
    # ----------------------------------------
    
    def vypocet_energie_DCB_R_MCC(self, Force: float, Displ: float, A1: float) -> float:
    
        Width_values = self.Return_dimensions('width')
        Height_values = self.Return_dimensions('height')
        B = min(Width_values)
        h = min(Height_values)
    
        C = (Displ/Force)
    
        defor_energ = (3 * Force**2 * C**(2/3))/(2*A1*B*h)
    
        return defor_energ
    
    # ----------------------------------------
    
    def found_indexes(self, time_delay_R: float) -> []:
    # Functions finds indexes in time_utm which correspons with time values in time_delam 
    
        # it is assumed that time step is constant for that reason it should not matter where the delta is calculated
        time_step = float(self.time_utm[6]) - float(self.time_utm[5])
    
        index_time = []
    
        for time_R in self.time_delam:
            # In case of time delay 
            time_del1 = float(time_R) - float(time_delay_R)
    
            # cyklus najde index case v datech trhacky
            for index, time_exp in enumerate(self.time_utm):
                if abs(float(time_del1) - float(time_exp)) <= time_step/2:
                    index_time.append(index)
                    break
    
        return index_time
    
    # ----------------------------------------
    
    def compute_compliance_MBT(self, num_iter_r: [], path: str) -> float:
        # Computing delta for R-Curves using Modified Beam Theory (MBT)
    
        defor_energ = np.zeros(len(num_iter_r))
        
        a = []
        C_1_div_3 = []
    
        for counter, j in enumerate(num_iter_r):
            C_1_div_3.append((self.displacement[j] / self.force[j])**(1/3))
            a.append(self.delam_length[counter])
    
        # c0 + c1*x + c2*x**2
        c0, c1 = np.polynomial.polynomial.Polynomial.fit(a, C_1_div_3, 1, domain=[])
        yFitted = np.polyval([c1, c0], a)
    
        delta_y = yFitted[1] - yFitted[0]
        delta_x = a[1] - a[0]
        n = delta_y / delta_x
    
        delta = -c0/c1
    
        # expanding the polyfit curve to be visible in plot
        x_axis = np.linspace(delta*5/4, max(a), 200)
        y_fitted_long = c1*x_axis + c0 
    
        fig = plt.figure()
        plt.plot(x_axis, y_fitted_long)
        plt.plot(a, C_1_div_3, 'o')
    
        plt.legend(["PolyFit", "Experiments"])
        plt.ylabel('$C^{1/3}$ [N/m]', fontsize=12)
        plt.xlabel('a [m]', fontsize=12)
        plt.title('MBT - Compliance', fontsize=14)
        plt.grid(True)
    
        path_fig_png = path + "_MBT" + ".png"
        path_fig_eps = path + "_MBT" + ".eps"
    
    #    fig.savefig(path_fig_eps)
        fig.savefig(path_fig_png)
    
        
    
        return abs(delta) 
    
    # ----------------------------------------
    
    def compute_compliance_CC(self, num_iter_r: [], path: str) -> float:
        # Computing delta for R-Curves using Compliance Calibration (CC) 
    
        defor_energ = np.zeros(len(num_iter_r))
        
        log_a = []
        log_C = []
    
        for counter, j in enumerate(num_iter_r):
            log_C.append(np.log(self.displacement[j] / self.force[j]))
            log_a.append(np.log(self.delam_length[counter]))
    
        # c0 + c1*x + c2*x**2
        c0, c1 = np.polynomial.polynomial.Polynomial.fit(log_a, log_C, 1, domain=[])
        yFitted = np.polyval([c1, c0], log_a)
    
        delta_y = abs(yFitted[1] - yFitted[0])
        delta_x = abs(log_a[1] - log_a[0])
        n = delta_y / delta_x
    
        delta = -c0/c1
    
        # expanding the polyfit curve to be visible in plot
        x_axis = np.linspace(min(log_a), max(log_a), 200)
        y_fitted_long = c1*x_axis + c0 
    
        fig = plt.figure()
        plt.plot(x_axis, y_fitted_long)
        plt.plot(log_a, log_C, 'o')
    
        plt.legend(["PolyFit", "Experiments"])
        plt.ylabel('$/log C$ [N/m]', fontsize=12)
        plt.xlabel('$/log a$ [m]', fontsize=12)
        plt.title('CC - Compliance', fontsize=14)
        plt.grid(True)
    
        path_fig_png = path + "_CC" + ".png"
        path_fig_eps = path + "_CC" + ".eps"
    
    #    fig.savefig(path_fig_eps)
        fig.savefig(path_fig_png)
    
        
    
        return n 
    
    # ----------------------------------------
    
    def compute_compliance_MCC(self, num_iter_r: [], path: str) -> float:
        # Computing delta for R-Curves using Modified Compliance Calibration (MCC) 
    
        Height_values = self.Return_dimensions('height')
        h = min(Height_values)

        defor_energ = np.zeros(len(num_iter_r))
        
        a_div_h = []
        C_pow_1_3 = []
    
        for counter, j in enumerate(num_iter_r):
            a_div_h.append(self.delam_length[counter]/h)
            C_pow_1_3.append((self.displacement[j] / self.force[j])**(1/3))
    
        # c0 + c1*x + c2*x**2
        c0, c1 = np.polynomial.polynomial.Polynomial.fit(C_pow_1_3, a_div_h, 1, domain=[])
        yFitted = np.polyval([c1, c0], C_pow_1_3)

        delta_y = abs(yFitted[1] - yFitted[0])
        delta_x = abs(C_pow_1_3[1] - C_pow_1_3[0])
        n = delta_y / delta_x
    
        delta = -c0/c1
    
        # expanding the polyfit curve to be visible in plot
        x_axis = np.linspace(min(C_pow_1_3), max(C_pow_1_3), 200)
        y_fitted_long = c1*x_axis + c0 
    
        fig = plt.figure()
        plt.plot(x_axis, y_fitted_long)
        plt.plot(C_pow_1_3, a_div_h, 'o')
    
        plt.legend(["PolyFit", "Experiments"])
        plt.ylabel('$a/h$ [N/m]', fontsize=12)
        plt.xlabel('$C^{1/3}$ [m]', fontsize=12)
        plt.title('MCC - Compliance', fontsize=14)
        plt.grid(True)
    
        path_fig_png = path + "_MCC" + ".png"
        path_fig_eps = path + "_MCC" + ".eps"
    
    #    fig.savefig(path_fig_eps)
        fig.savefig(path_fig_png)
    
        return n
    
    # ----------------------------------------
    
    
    
    def create_R_cuve(self, time_delay_R: float, path: str, method: str):
        # parameter method has possible string values 'MBT', 'CC', 'MCC'
    
        print(f'Evaluating DCB Fracture energie with {method} model')
    
        num_iter_r = self.found_indexes(time_delay_R)
    
        DCB_G_I = []
        delam_set = []

        # compute the compliance for one sample
        if method == "MBT":
            delta_MBT = self.compute_compliance_MBT(num_iter_r, path)
            print('delta MBT is: ', delta_MBT)
        elif method == "CC":
            delta_CC = self.compute_compliance_CC(num_iter_r, path)
            print('delta CC is: ', delta_CC)
        elif method == "MCC":
            delta_MCC = self.compute_compliance_MCC(num_iter_r, path)
            print('delta MCC is: ', delta_MCC)
    
        # loop thorough every delamination length value
        for counter, index in enumerate(num_iter_r):
    
            force_del = self.force[index]
            epsilon_del = self.displacement[index]
            Delam_length = self.delam_length[counter]
            
            logger.info("index " + str(index) + " force " + str(force_del) + " epsilon " + str(epsilon_del) + " time " + str(self.time_utm[index]))
    
            if method == "MBT":
                G_I = self.vypocet_energie_DCB_R_MBT(force_del, Delam_length, epsilon_del, delta_MBT)
            elif method == "CC":
                G_I = self.vypocet_energie_DCB_R_CC(force_del, Delam_length, epsilon_del, delta_CC)
            elif method == "MCC":
                G_I = self.vypocet_energie_DCB_R_MCC(force_del, epsilon_del, delta_MCC)
    
            DCB_G_I.append(G_I)
    
        return DCB_G_I, self.delam_length 
    
    # ----------------------------------------

# ----------------------------------------
# Graphical evaluation
def graph_DCB_R_curves(F: [[]], epsilon: [[]], name_TAR, color_set, marker_set, line_set, path_fig, range_min):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # constant average value from range_x1
    avg_x_const, avg_y_const = average_curve(epsilon, F, range_min)
    # average function
    avg_x, avg_y = average_curve(epsilon, F)

    # plot all R curves
    for cons in range(0, len(F)):
        plt.plot(epsilon[cons], F[cons], color=color_set[cons], marker=marker_set[cons], linestyle=line_set[cons])

    plt.plot(avg_x, avg_y, color='green')
    
    plt.xlabel(r'Delamination legth $a$ [m]')
    plt.ylabel(r'Fracture toughness $G_I$ [J$m^{-2}$]')

    plt.legend(name_TAR)

    path_fig_png = path_fig + ".png"
    path_fig_eps = path_fig + ".eps"

    # Plotting the constant curve which is the mean value from set range
    plt.plot(avg_x_const, avg_y_const, linewidth=3, color='k', linestyle='--')
    print("graph_DCB_R was created")

    fig.savefig(path_fig_eps)
    fig.savefig(path_fig_png)







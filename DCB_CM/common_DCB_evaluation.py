#
# Common function for evaluating DCB experiments

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



def graf_Feps(F, epsilon, name_TAR, color_set, marker_set, line_set, path_fig):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for cons in range(0, len(F)):
        for i in range(0, len(epsilon[cons])):
            epsilon[cons][i] = epsilon[cons][i]*1000

    for cons in range(0, len(F)):
        plt.plot(epsilon[cons], F[cons], color=color_set[cons],
                linestyle=line_set[cons])

    # nacteni prumerovane funkce
    F_average, min_count = value_TAR_avereged(F)

    epsi_value = 0
    for cons in range(0, len(F)):
        if min_count == len(F[cons]):
            epsi_value = cons

    plt.plot(epsilon[epsi_value], F_average, color='green')

    plt.xlabel(r'Displacement $\delta$ [mm]')
    plt.ylabel(r'Force $P$ [N]')

# =============================================================================
#     name_TAR.append("DCB výsledná křivka")
# =============================================================================
    plt.legend(name_TAR)

    path_fig_png = path_fig + ".png"
    path_fig_eps = path_fig + ".eps"

    fig.savefig(path_fig_eps)
    fig.savefig(path_fig_png)

    print("graf_Feps vytvoren graf")

# -------------------- graf


def graf_Feps_R_curves(F: [[]], epsilon, name_TAR, color_set, marker_set, line_set, path_fig, Range):

    # Converting from mm to m
    dimensions = 1e-3 
    range_x0 = Range[0]*dimensions
    range_x1 = Range[1]*dimensions
    range_x2 = Range[2]*dimensions
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    y_average = 0
    counter = 0

    print("color set", color_set)
    print(f"Values: {F}")

    # This calucate the average value after stabilization of values
    for cons in range(0, len(F)):
        for i in range(0, len(epsilon[cons])):
#            epsilon[cons][i] = epsilon[cons][i]
            if epsilon[cons][i] > range_x1:
                y_average += F[cons][i]
                counter += 1

    # plot all R curves
    for cons in range(0, len(F)):
        plt.plot(epsilon[cons], F[cons], color=color_set[cons], marker=marker_set[cons], linestyle=line_set[cons])

    # nacteni prumerovane funkce
    mean_x_axis = np.linspace(range_x0, range_x2, 200)
    ys_interp = [np.interp(mean_x_axis, epsilon[i], F[i]) for i in range(len(epsilon))]
    mean_y_axis = np.mean(ys_interp, axis=0)

    plt.plot(mean_x_axis, mean_y_axis, color='green')
    
    plt.xlabel(r'Delamination legth $a$ [m]')
    plt.ylabel(r'Fracture toughness $G_I$ [J$m^{-2}$]')

    plt.legend(name_TAR)

    path_fig_png = path_fig + ".png"
    path_fig_eps = path_fig + ".eps"

    # Plotting the constant curve which is the mean value from set range
    plt.plot([range_x1, range_x2], [y_average/counter, y_average/counter], linewidth=3, color='k', linestyle='--')
    print("graf_Feps_R vytvoren graf")
    print(y_average/counter)

    fig.savefig(path_fig_eps)
    fig.savefig(path_fig_png)



# -------------------- TAR
def value_TAR(path_TAR: str, drop_lines: int) -> tuple[list, list, list, list, list]:
    # Function reads values from .TAR file to variable column_data 

    # Retrieving Name of TAR file from path_TAR
    name_TAR = path_TAR.split('/')
    name_TAR = name_TAR[len(name_TAR)-1]
    name_TAR = name_TAR[0:len(name_TAR)-4]

    f = open(path_TAR, "r", encoding='utf-8', errors='ignore')

    for s in range(0, drop_lines-1):
        f.readline()

    values = f.readline()
    first_row = f.readline()
    num_columns = len(first_row.split())
    
    print(f'TAR is loaded with this column (num. {num_columns} values:')
    print(values)

    column_data = []

    temp_values = first_row.split()

    # Adding the first_row of values which was used to get number of columns
    for counter, value in enumerate(temp_values):
        column_data.append([float(value)])

    for x in f:
        temp_values = x.split()

        for counter, value in enumerate(temp_values):
            column_data[counter].append(float(value))

    print("==============================")
    print(f"Data loaded for {name_TAR}")
    logger.info("==============================")
    logger.info(f"evaluating {name_TAR}")

    f.close()
    return column_data

# ----------------------------------------

def value_TAR_avereged(F):

    F_average = []
    temp1 = []
    F_prumer = 0

    for cons3 in range(0, len(F)):
        temp1.append(len(F[cons3]))
    temp1 = np.array(temp1)
    min_count = min(temp1)

    for cons2 in range(0, min_count):
        F_prumer = 0
        for cons4 in range(0, len(F)):
            F_prumer += F[cons4][cons2]
        F_prumer = F_prumer/len(F)
        F_average.append(F_prumer)

    return F_average, min_count

# ----------------------------------------

def vypocet_energie_DCB(F, a, epsilon_n):
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

def vypocet_energie_DCB_R_MBT(Force: float, a: 'class', delam_R: float, Displ: [], delta: float) -> float:

    B_index = a.data_name.index('B')
    ai_index = a.data_name.index('ai')

    B = float(a.data_value[B_index])*float(a.data_dimension[B_index])  # sirka vzorku
    
    # delamination length
    a_del = delam_R

    deltaC = float(Displ)
    
    defor_energ = (3*Force*deltaC)/(2*B*(a_del+delta))


    return defor_energ

# ----------------------------------------

def vypocet_energie_DCB_R_CC(Force: float, a: 'class', delam_R: float, Displ: [], n: float) -> float:

    B_index = a.data_name.index('B')

    B = float(a.data_value[B_index])*float(a.data_dimension[B_index])  # sirka vzorku
    
    # delamination length
    a_del = delam_R

    deltaC = float(Displ)

    defor_energ = (n*Force*deltaC)/(2*B*a_del)

    return defor_energ

# ----------------------------------------

def vypocet_energie_DCB_R_MCC(Force: float, a: 'class', Displ: [], A1: float) -> float:

    B_index = a.data_name.index('B')
    h_index = a.data_name.index('2hA')

    B = float(a.data_value[B_index])*float(a.data_dimension[B_index])  # sirka vzorku
    h = float(a.data_value[h_index])*float(a.data_dimension[h_index])  # vyska vzorku

    C = (Displ/Force)

    defor_energ = (3 * Force**2 * C**(2/3))/(2*A1*B*h)

    return defor_energ

# ----------------------------------------

def found_indexes(Time_exp_ar: [], Time_R_ar: []) -> []:

    time_delay_R = 0 

    # it is assumed that time step is constant for that reason it should not matter where the delta is calculated
    time_step = float(Time_exp_ar[6]) - float(Time_exp_ar[5])

    index_time = []

    for time_R in Time_R_ar:

        time_del1 = float(time_R) - float(time_delay_R)

        # cyklus najde index casu v datech trhacky
        for index, time_exp in enumerate(Time_exp_ar):
            if float(time_del1) < float(time_exp) + time_step/2 and float(time_del1) > float(time_exp) - time_step/2:
                index_time.append(index)
                break


    return index_time

# ----------------------------------------

def found_indexes_Voltage(Voltage: []) -> []:

    pass_value = 1 # corresponds with 1 V 
    index_voltage = []

    first_index = None 
    for counter, v in enumerate(Voltage):
        if v >= pass_value and first_index == None:
            first_index = counter
        if first_index != None and v <= pass_value:
            index_voltage.append(first_index)
            first_index = None 

    return index_voltage

# ----------------------------------------

def compute_compliance_MBT(Force: [], Displ: [], Del_length_R: [], Time_R_ar: [], Time_exp_ar: [], path: str) -> float:
    # Computing delta for R-Curves using Modified Beam Theory (MBT)

    # Retrieving Name of TAR file from path_TAR
    name = path.split('/')
    name = name[len(name)-1]
    name = name[0:len(name)-4]

    num_iter_r = found_indexes(Time_exp_ar, Time_R_ar)

    defor_energ = np.zeros(len(num_iter_r))
    
    a = []
    C_1_div_3 = []

    for counter, j in enumerate(num_iter_r):
        C_1_div_3.append((Displ[j] / Force[j])**(1/3))
        a.append(Del_length_R[counter])

    coefficients = np.polyfit(a, C_1_div_3, 1)
    yFitted = np.polyval(coefficients, a)

    delta_y = yFitted[1] - yFitted[0]
    delta_x = a[1] - a[0]
    n = delta_y / delta_x

    a1 = coefficients[0]
    b2 = coefficients[1]

    delta = -b2/a1

    # expanding the polyfit curve to be visible in plot
    x_axis = np.linspace(delta*5/4, max(a), 200)
    y_fitted_long = a1 * x_axis + b2

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

def compute_compliance_CC(Force: [], Displ: [], Del_length_R: [], Time_R_ar: [], Time_exp_ar: [], path: str) -> float:
    # Computing delta for R-Curves using Compliance Calibration (CC) 

    # Retrieving Name of TAR file from path_TAR
    name = path.split('/')
    name = name[len(name)-1]
    name = name[0:len(name)-4]

    num_iter_r = found_indexes(Time_exp_ar, Time_R_ar)

    defor_energ = np.zeros(len(num_iter_r))
    
    log_a = []
    log_C = []

    for counter, j in enumerate(num_iter_r):
        log_C.append(np.log(Displ[j] / Force[j]))
        log_a.append(np.log(Del_length_R[counter]))

    coefficients = np.polyfit(log_a, log_C, 1)
    yFitted = np.polyval(coefficients, log_a)

    delta_y = abs(yFitted[1] - yFitted[0])
    delta_x = abs(log_a[1] - log_a[0])
    n = delta_y / delta_x

    a1 = coefficients[0]
    b2 = coefficients[1]

    delta = -b2/a1

    # expanding the polyfit curve to be visible in plot
    x_axis = np.linspace(min(log_a), max(log_a), 200)
    y_fitted_long = a1 * x_axis + b2

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

def compute_compliance_MCC(Force: [], Displ: [], Del_length_R: [], Time_R_ar: [], Time_exp_ar: [], path: str,h: float) -> float:
    # Computing delta for R-Curves using Modified Compliance Calibration (MCC) 

    # Retrieving Name of TAR file from path_TAR
    name = path.split('/')
    name = name[len(name)-1]
    name = name[0:len(name)-4]

    num_iter_r = found_indexes(Time_exp_ar, Time_R_ar)

    defor_energ = np.zeros(len(num_iter_r))
    
    a_div_h = []
    C_pow_1_3 = []

    for counter, j in enumerate(num_iter_r):
        a_div_h.append(Del_length_R[counter]/h)
        C_pow_1_3.append((Displ[j] / Force[j])**(1/3))

    coefficients = np.polyfit(C_pow_1_3, a_div_h, 1)
    yFitted = np.polyval(coefficients, C_pow_1_3)

    delta_y = abs(yFitted[1] - yFitted[0])
    delta_x = abs(C_pow_1_3[1] - C_pow_1_3[0])
    n = delta_y / delta_x

    a1 = coefficients[0]
    b2 = coefficients[1]

    delta = -b2/a1

    # expanding the polyfit curve to be visible in plot
    x_axis = np.linspace(min(C_pow_1_3), max(C_pow_1_3), 200)
    y_fitted_long = a1 * x_axis + b2

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
def compute_compliance_ENF(Force: [], Displ: [], Del_length_R: [], Time_R_ar: [], Time_exp_ar: [], path: str,h: float, num_iter: [int]) -> float:
    # Computing compliance m for ENF 

    # Retrieving Name of TAR file from path_TAR
    name = path.split('/')
    name = name[len(name)-1]
    name = name[0:len(name)-4]

    #num_iter_r = found_indexes(Time_exp_ar, Time_R_ar)

    defor_energ = np.zeros(len(num_iter))
    
    a = []
    C = []

    for counter, j in enumerate(num_iter):
        a.append(Del_length_R[counter]**3)
        C.append(Displ[j] / Force[j])

    coefficients = np.polyfit(C, a, 1)
    yFitted = np.polyval(coefficients, C)

    delta_y = abs(yFitted[1] - yFitted[0])
    delta_x = abs(C[1] - C[0])
    m = delta_x / delta_y

    a1 = coefficients[0]
    b2 = coefficients[1]

    delta = -b2/a1

    # expanding the polyfit curve to be visible in plot
    x_axis = np.linspace(min(C), max(C), 200)
    y_fitted_long = a1 * x_axis + b2

    fig = plt.figure()
    plt.plot(x_axis, y_fitted_long)
    plt.plot(C, a, 'o')

    plt.legend(["PolyFit", "Experiments"])
    plt.ylabel('$a^3$ [N/m]', fontsize=12)
    plt.xlabel('$C$ [m]', fontsize=12)
    plt.title('ENF - Compliance', fontsize=14)
    plt.grid(True)

    path_fig_png = path + "_ENF" + ".png"
    path_fig_eps = path + "_ENF" + ".eps"

#    fig.savefig(path_fig_eps)
    fig.savefig(path_fig_png)

    return m


# ----------------------------------------
def vypocet_energie_ENF(F_max: float, a: 'class', m: float):

    B_index = a.data_name.index('B')
    ai_index = a.data_name.index('ai')

    B = float(a.data_value[B_index]) * float(a.data_dimension[B_index])  # sirka vzorku
    # non ahezive od osy zatizeni
    a_0 = float(a.data_value[ai_index])*float(a.data_dimension[ai_index])

    defor_energ = (3*m*(F_max)**2 *(a_0)**2)/(2*B) # Norm D7905


    # Bernardin Thesis approach Un-known source
    # index_F_max = (np.argmax(F))
    # float_lst = [float(x) for x in epsilon_n]
    # deltaC = float_lst[index_F_max]  # posunuti hodnoty v bode F_max v [m]
    # defor_energ = (9*a_del**2 * F_max * deltaC)/(2*B*(2*L**3 + 3*a_del**3)) # Bernardin Thesis

    return defor_energ

# ----------------------------------------


def vypocet_energie_ENF_R(F_max: float, a: 'class', m: float, del_length: float):

    B_index = a.data_name.index('B')
    ai_index = a.data_name.index('ai')

    B = float(a.data_value[B_index]) * float(a.data_dimension[B_index])  # sirka vzorku
    # non ahezive od osy zatizeni
    a_0 = float(a.data_value[ai_index])*float(a.data_dimension[ai_index])

    defor_energ = (3*m*(F_max)**2 *(del_length)**2)/(2*B) # Norm D7905


    # Bernardin Thesis approach Un-known source
    # index_F_max = (np.argmax(F))
    # float_lst = [float(x) for x in epsilon_n]
    # deltaC = float_lst[index_F_max]  # posunuti hodnoty v bode F_max v [m]
    # defor_energ = (9*a_del**2 * F_max * deltaC)/(2*B*(2*L**3 + 3*a_del**3)) # Bernardin Thesis

    return defor_energ

# ----------------------------------------


def vypocet_energie_MixI_II(F, a):
    # Not validated

    b_prum = ((float(a.data_value[5]) + float(a.data_value[6]) +
              float(a.data_value[7]))/3)*a.Values_dimension[5]  # sirka vzorku
    a_del = float(a.data_value[10])  # delamination length
    F = np.array(F)
    F_max = max(F)
    c = 0  # lever length
    L = 0  # delka mezi osami zatizeni
    Ksi = 0  # crack length correction parametr
    h_prum = ((float(a.data_value[2]) + float(a.data_value[3]) + float(
        a.data_value[4]))/3)*a.Values_dimension[2]  # specimen thisckness
    Ef = 0  # modul of elesticity in fiber direction

    GI = ((12*F ^ 2*(3*c - L) ^ 2)*(a_del + Ksi*h) ^ 2) / \
        (16*b ^ 2 * h ^ 3 * L ^ 2 * Ef)
    GII = (((9*P ^ 2*(c+L) ^ 2))*(a_del + 0.42*Ksi*h) ^ 2) / \
        (16*b ^ 2 * h ^ 3 * L ^ 2 * Ef)
    G = GI + GII

    return G


# ----------------------------------------

def Variacni_koef(F) -> float:

    var_x = 0
    avrg_F = np.average(F)

    for i in F:
        var_x += (1/len(F))*(i-avrg_F)*(i-avrg_F)

    s_x = math.sqrt(var_x)

    v_x = s_x/avrg_F

    return v_x

# ----------------------------------------


def create_R_cuve(Time_exp_ar, Time_R_ar, Del_length_R, Force, a, epsilon, delam_a0_R, time_delay_R, path: str, method: str):
    # parameter method has possible string values 'MBT', 'CC', 'MCC'

    print(f'Evaluating DCB Fracture energie with {method} model')

    h_index = a.data_name.index('2hA')
    h = float(a.data_value[h_index])*float(a.data_dimension[h_index])  # specimen thickness 

    # Retrieving Name of TAR file from path_TAR
    name = path.split('/')
    name = name[len(name)-1]
    name = name[0:len(name)-4]

    num_iter_r = found_indexes(Time_exp_ar, Time_R_ar)

    DCB_G_I = []
    delam_set = []

    # add intitial delamination to values
    Del_length_R_compl = [float(delam_a0_R) + float(x) for x in Del_length_R]

    # compute the compliance for one sample
    if method == "MBT":
        delta_MBT = compute_compliance_MBT(Force, epsilon, Del_length_R_compl, Time_R_ar, Time_exp_ar, path)
        print('delta MBT is: ', delta_MBT)
    elif method == "CC":
        delta_CC = compute_compliance_CC(Force, epsilon, Del_length_R_compl, Time_R_ar, Time_exp_ar, path)
        print('delta CC is: ', delta_CC)
    elif method == "MCC":
        delta_MCC = compute_compliance_MCC(Force, epsilon, Del_length_R_compl, Time_R_ar, Time_exp_ar, path, h)
        print('delta MCC is: ', delta_MCC)

    # loop thorough every delamination length value
    for counter, index in enumerate(num_iter_r):

        force_del = Force[index]
        epsilon_del = epsilon[index]
        Delam_length = Del_length_R_compl[counter]
        
        logger.info("index " + str(index) + " force " + str(force_del) + " epsilon " + str(epsilon_del) + " time " + str(Time_exp_ar[index]))

        if method == "MBT":
            G_I = vypocet_energie_DCB_R_MBT(force_del, a, Delam_length, epsilon_del, delta_MBT)
        elif method == "CC":
            G_I = vypocet_energie_DCB_R_CC(force_del, a, Delam_length, epsilon_del, delta_CC)
        elif method == "MCC":
            G_I = vypocet_energie_DCB_R_MCC(force_del, a, epsilon_del, delta_MCC)

        DCB_G_I.append(G_I)


    return DCB_G_I, Del_length_R_compl

# ----------------------------------------



def create_R_cuve_ENF(Time_exp_ar, Time_R_ar, Del_length_R, Force, a, epsilon, delam_a0_R, time_delay_R, path: str, voltage: list, voltage_list):

    print(f'Evaluating ENF Fracture energie')

    h_index = a.data_name.index('2hA')
    h = float(a.data_value[h_index])*float(a.data_dimension[h_index])  # specimen thickness 

    # Retrieving Name of TAR file from path_TAR
    name = path.split('/')
    name = name[len(name)-1]
    name = name[0:len(name)-4]

    num_iter_r = found_indexes(Time_exp_ar, Time_R_ar)
    num_iter_r_v = found_indexes_Voltage(voltage)

    print(f"num_iter_r_v: {num_iter_r_v}")
    print(f"length: {len(num_iter_r_v)}")

    # indexes are striped to evaluated values from photos
    num_iter_r_v = num_iter_r_v[voltage_list[0]:(voltage_list[0]+voltage_list[1])]

    print(f"num_iter_r_v: {num_iter_r_v}")
    print(f"length: {len(num_iter_r_v)}")
    print(f"del_length length: {len(Del_length_R)}")
    

    ENF_G_II = []
    delam_set = []

    F_max = max(Force)

    # add intitial delamination to values
    Del_length_R_compl = [float(delam_a0_R) + float(x) for x in Del_length_R]

    # compute the compliance for one sample
    m = compute_compliance_ENF(Force, epsilon, Del_length_R_compl, Time_R_ar, Time_exp_ar, path, h, num_iter_r_v)
    print('Compliance m for ENF is: ', m)

    # loop thorough every delamination length value
    for counter, index in enumerate(num_iter_r):

        force_del = Force[index]
        epsilon_del = epsilon[index]
        Delam_length = Del_length_R_compl[counter]
        
        logger.info("index " + str(index) + " force " + str(force_del) + " epsilon " + str(epsilon_del) + " time " + str(Time_exp_ar[index]))

        #G_II = vypocet_energie_ENF(F_max, a, m)
        G_II = vypocet_energie_ENF_R(force_del, a, m, Delam_length)

        ENF_G_II.append(G_II)


    return ENF_G_II, Del_length_R_compl

# ----------------------------------------


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

    def find_close_value_index(self, lst, x, range_val):
    # Find the index of a value in the list that is within a small range of x.
        for i, y in enumerate(lst):
            if abs(x - y) <= range_val:
                return i
        return -1

    def Evaluate_Young_modulus(self) -> float:
    #Funciton computes Young modulus from stress strain curve
        
        # These are strain limits where young modulus should be evaluated according to ISO 527_1
        strain_lim_1 = 0.0005
        strain_lim_2 = 0.0025

        stress = self.Calculate_stress()
        strain = self.Calculate_strain()

        # Range value to find index
        delta = abs(strain[round(len(strain)/10)] - strain[round(len(strain)/10)-1])/2

        # Finds index of limits
        index_lim_1 = self.find_close_value_index(strain, strain_lim_1, delta)
        index_lim_2 = self.find_close_value_index(strain, strain_lim_2, delta)

        # Modification of dimension to correspond with evaluation limits
        stress = stress[index_lim_1:index_lim_2]
        strain = strain[index_lim_1:index_lim_2]

        # This coefficient will not be the same with polyfit. for same resaults domain=[]
        c0, c1 = np.polynomial.polynomial.Polynomial.fit(strain,stress, 1, domain=[])

        print(f"Young modulues E = {c1}")
    
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

        print(f"Poisson ratio nu = {poisson}")
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






















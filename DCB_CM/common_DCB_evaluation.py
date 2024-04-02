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


def graf_Feps_R_curves(F: [[]], epsilon, name_TAR, color_set, marker_set, line_set, path_fig):

    # Converting from mm to m
    dimensions = 1e-3 
    range_x0 = 30*dimensions
    range_x1 = 70*dimensions
    range_x2 = 150*dimensions
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    y_average = 0
    counter = 0

    print("color set", color_set)

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
    
# =============================================================================
#     mean_x_axis = [i for i in range(max(epsilon))]
# =============================================================================
    mean_x_axis = np.linspace(range_x0, range_x2, 200)
    ys_interp = [np.interp(mean_x_axis, epsilon[i], F[i]) for i in range(len(epsilon))]
    mean_y_axis = np.mean(ys_interp, axis=0)

    plt.plot(mean_x_axis, mean_y_axis, color='green')
    
# =============================================================================
#     F_average, min_count = value_TAR_avereged(F)
# 
#     epsi_value = 0
#     for cons in range(0, len(F)):
#         if min_count == len(F[cons]):
#             epsi_value = cons
# 
#     pl.plot(epsilon[epsi_value], F_average, color='green')
# =============================================================================

    plt.xlabel(r'Delamination legth $a$ [m]')
    plt.ylabel(r'Fracture toughness $G_I$ [J$m^{-2}$]')

# =============================================================================
#     name_TAR.append("DCB výsledná křivka")
# =============================================================================
    plt.legend(name_TAR)

#    path_fig_png = id_directory_fig + "\\" + 'graf_R_' + list + ".png"
#    path_fig_eps = id_directory_fig + "\\" + 'graf_R_' + list + ".eps"

    path_fig_png = path_fig + ".png"
    path_fig_eps = path_fig + ".eps"

    plt.plot([range_x1, range_x2], [y_average/counter, y_average/counter], linewidth=3, color='k', linestyle='--')
    print("graf_Feps_R vytvoren graf")
    print(y_average/counter)

    fig.savefig(path_fig_eps)
    fig.savefig(path_fig_png)



# -------------------- TAR
def value_TAR(path_TAR: str):

    number_drop_lines = 7

    # Retrieving Name of TAR file from path_TAR
    name_TAR = path_TAR.split('/')
    name_TAR = name_TAR[len(name_TAR)-1]
    name_TAR = name_TAR[0:len(name_TAR)-4]

    f = open(path_TAR, "r", encoding='utf-8', errors='ignore')

    # for s in range(0,5):
    for s in range(0, number_drop_lines):
        f.readline()

    F = []
    epsilon = []
    epsilon_n = []
    time = []

    for x in f:
        temp_values = x.split()

        F.append(temp_values[0])
        # print(F)
        epsilon.append(temp_values[2])
        # print(epsilon)
        epsilon_n.append(temp_values[2])
        #
        time.append(temp_values[1])

    for i6 in range(0, len(F)):
        epsilon[i6] = float(epsilon[i6])
        # pouze pro vykresleni v mmm
        # epsilon[i6] = float(epsilon[i6])*1000
        F[i6] = float(F[i6])

# pro grafy po jednom
        # graf_Feps(F, epsilon, name_TAR)


    print("==============================")
    print(f"Data loaded for {name_TAR}")
    logger.info("==============================")
    logger.info(f"evaluating {name_TAR}")

    return F, epsilon, name_TAR, epsilon_n, time
    f.close()

# ----------------------------------------


def prumer(x):

    for temp3 in range(0, len(x)):
        x_prumer += x[temp3]
    x_prumer = x_prumer/len(x)

    return x_prumer


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




# ----------------------------------------

def vypocet_energie_DCB(F, a, epsilon_n):

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

    # kontorola

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

def vypocet_energie_DCB_R_MCC(Force: float, a: 'class', delam_R: float, Displ: [], A1: float) -> float:

    B_index = a.data_name.index('B')
    h_index = a.data_name.index('2hA')

    B = float(a.data_value[B_index])*float(a.data_dimension[B_index])  # sirka vzorku
    h = float(a.data_value[h_index])*float(a.data_dimension[h_index])  # vyska vzorku

    C = (Displ/Force)

    defor_energ = (3 * Force**2 * C**(2/3))/(2*A1*B*h)

    return defor_energ

# ----------------------------------------

def found_indexes(Time_exp_ar: [], Time_R_ar: []) -> []:

    time_delay_R = 0.0
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

    #x_axis = np.linspace(delta*5/4, abs(delta*5/4), 200)
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

    #x_axis = np.linspace(delta*5/4, abs(delta*5/4), 200)
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

    #x_axis = np.linspace(delta*5/4, abs(delta*5/4), 200)
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
def vypocet_energie_ENF(F, a, epsilon_n):

    B = float(a.data_value[B_index]) * \
        float(a.data_dimension[B_index])  # sirka vzorku
    # non ahezive od osy zatizeni
    a_del = float(a.data_value[ai_index])*float(a.data_dimension[ai_index])
    L = (float(a.data_value[lj_index]))*(float(a.data_dimension[lj_index]))
    # m = int(a.data_value[m_index])*a.data_dimension[m_index]
    F_max = max(F)
    index_F_max = (np.argmax(F))
    float_lst = [float(x) for x in epsilon_n]
    deltaC = float_lst[index_F_max]  # posunuti hodnoty v bode F_max v [m]

    # defor_energ = (3*m*(F_max)**2 *(a_del)**2)/(2*b_prum) #podle normy
    defor_energ = (9*a_del**2 * F_max * deltaC)/(2*B*(2*L**3 + 3*a_del**3))

    return defor_energ

# ----------------------------------------


def vypocet_energie_ENF_R(F, a, epsilon_n, delam_R, delta_R):

    B = float(a.data_value[B_index]) * \
        float(a.data_dimension[B_index])  # sirka vzorku
    # non ahezive od osy zatizeni
    a_del = float(delam_R)*float(a.data_dimension[ai_index])
    L = (float(a.data_value[lj_index]))*(float(a.data_dimension[lj_index]))
    h0 = float(a.data_value[h0_index])*(float(a.data_dimension[h0_index]))
    hA = float(a.data_value[hA_index])*(float(a.data_dimension[hA_index]))
    hB = float(a.data_value[hB_index])*(float(a.data_dimension[hB_index]))
    # m = int(a.data_value[m_index])*a.data_dimension[m_index]

    E1 = 51.10e9
    h2 = (h0+hA+hB)/3
    h = h2/2

    print(a_del)

    deltaC = float(delta_R)

    # defor_energ = (3*m*(F_max)**2 *(a_del)**2)/(2*b_prum) #podle normy
    # defor_energ = (9*a_del**2 *F *deltaC)/(2*B*(2*L**3 +3*a_del**3))  #podle Bernardina

    defor_energ = (9*F**2*a_del**2)/(16*B**2*h**3*E1)

    return defor_energ

# ----------------------------------------


def vypocet_energie_MixI_II(F, a):

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
# ----------------------------------------

def Variacni_koef(F):

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

    print(f'Evaluating DCB Fracture energiee with {method} model')

    h_index = a.data_name.index('2hA')
    h = float(a.data_value[h_index])*float(a.data_dimension[h_index])  # specimen thickness 

    # Retrieving Name of TAR file from path_TAR
    name = path.split('/')
    name = name[len(name)-1]
    name = name[0:len(name)-4]

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
    for i in range(0, len(Time_R_ar)):

        time_del1 = float(Time_R_ar.iloc[i]) - float(time_delay_R)
        delta = float(Time_exp_ar[6]) - float(Time_exp_ar[5])
        # cyklus najde index casu v datech trhacky
        for j in range(0, len(Time_exp_ar)):
            if float(time_del1) < float(Time_exp_ar[j]) + delta/2 and float(time_del1) > float(Time_exp_ar[j]) - delta/2:
                break

        force_del = Force[j]
        epsilon_del = epsilon[j]
        Delam_length = Del_length_R_compl[i]
        
        logger.info("index " + str(j) + " force " + str(force_del) + " epsilon " + str(epsilon_del) + " time " + str(time_del1))

        if method == "MBT":
            G_I = vypocet_energie_DCB_R_MBT(force_del, a, Delam_length, epsilon_del, delta_MBT)
        elif method == "CC":
            G_I = vypocet_energie_DCB_R_CC(force_del, a, Delam_length, epsilon_del, delta_CC)
        elif method == "MCC":
            G_I = vypocet_energie_DCB_R_MCC(force_del, a, Delam_length, epsilon_del, delta_MCC)

        DCB_G_I.append(G_I)


    return DCB_G_I, Del_length_R_compl

# ----------------------------------------







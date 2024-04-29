
# import load_workbook
from openpyxl import load_workbook
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import os
import math

# Importing common python functions and classes 
import sys
sys.path.append('/home/horakl/School/FAV/Doktorske_studium/Prace_projekty/Mechanical_testing/Python_files/DCB_CM')
from evaluate_dcb import *
from common_evaluation import *

# ----------------------------------------

id_directory_python = os.path.dirname(os.path.abspath(__file__))

if id_directory_python.endswith("python") or id_directory_python.endswith("Python"):
    id_directory_main = id_directory_python[:-6]

id_directory_data       = id_directory_main + "data"
id_directory_fig        = id_directory_main + "fig"
id_directory_fig_check  = os.path.join(id_directory_fig, "fig_check")
id_directory_raw_data   = id_directory_main + "raw_data"

spec_data, ex_list = import_excel_values(id_directory_python)
delam_length, delam_time, strip_range = import_delam_length_values(id_directory_python, 1e-3)

# ----------------------------------------

# Names in excel 
width_names     = ["B"]
height_names    = ["2hA"]
length_names    = ["lj"]
delam_names     = ["ai"]

names_param = {'height': height_names, 'width': width_names, 'length': length_names, 'delam': delam_names}

color_index     = spec_data[0].data_name.index('color')
marker_index    = spec_data[0].data_name.index('marker')
line_index      = spec_data[0].data_name.index('line')

# ----------------------------------------

F = []
F_max = []
deltaY = []
name_TAR = []
time = []
delam_set_array = []

DCB_G_I_set_MBT = []
DCB_G_I_set_CC = []
DCB_G_I_set_MCC = []
delam_set_array = []

marker_set = []
line_set = []
color_set = []

for counter, specimen in enumerate(spec_data):
    if specimen.data_value[1] == 1 or specimen.data_value[1] == str(1):

        # path where to access data and store data
        name_TAR_str = specimen.data_value[0]
        path_TAR = os.path.join(id_directory_raw_data, (name_TAR_str + ".TRA"))
        path_fig_check = os.path.join(id_directory_fig_check, name_TAR_str)

        # reads values from .TAR file 
        drop_lines = 7
        column_data = value_TAR(path_TAR, drop_lines)
        # Output should be manualy changed according columns in TAR file
        Force_i = column_data[0]
        time_i = column_data[1]
        deltaY_i = column_data[2]

        delam_R = delam_length[counter]
        delam_t = delam_time[counter]
        
        # Creating list of list.
        F.append(Force_i)
        deltaY.append(deltaY_i)
        name_TAR.append(name_TAR_str)
        time.append(time_i)

        color_set.append(specimen.data_value[color_index])
        line_set.append(specimen.data_value[line_index])
        marker_set.append(specimen.data_value[marker_index])

        time_delay_R = 0

        # initiaze class
        evaluate_DCB = Evaluate_DCB(Force_i, deltaY_i, time_i, delam_t, delam_R, specimen, names_param, name_TAR_str)

        DCB_G_I_MBT, delam_set = evaluate_DCB.create_R_cuve(time_delay_R, path_fig_check, "MBT")
        DCB_G_I_CC, delam_set = evaluate_DCB.create_R_cuve(time_delay_R, path_fig_check, "CC")
        DCB_G_I_MCC, delam_set = evaluate_DCB.create_R_cuve(time_delay_R, path_fig_check, "MCC")

        F_max.append(max(Force_i))

        DCB_G_I_set_MBT.append(DCB_G_I_MBT)
        DCB_G_I_set_CC.append(DCB_G_I_CC)
        DCB_G_I_set_MCC.append(DCB_G_I_MCC)
        delam_set_array.append(delam_set)


avrg_F_max = np.average(F_max)
print("Prumer maximalnich sil: " + ex_list + " " + str(avrg_F_max))

path_fig        = os.path.join(id_directory_fig, ('graf_' + ex_list))
path_fig_R_MBT  = os.path.join(id_directory_fig, ('MBT_graf_R_' + ex_list))
path_fig_R_CC   = os.path.join(id_directory_fig, ('CC_graf_R_' + ex_list))
path_fig_R_MCC  = os.path.join(id_directory_fig, ('MCC_graf_R_' + ex_list))

graf_Feps(F, deltaY, name_TAR, color_set, marker_set, line_set, path_fig)

# Converting from mm to m
dimensions = 1e-3 
range_min = 70
range_min = range_min*dimensions

graph_DCB_R_curves(DCB_G_I_set_MBT, delam_set_array, name_TAR, color_set, marker_set, line_set, path_fig_R_MBT, range_min)
graph_DCB_R_curves(DCB_G_I_set_CC, delam_set_array, name_TAR, color_set, marker_set, line_set, path_fig_R_CC, range_min)
graph_DCB_R_curves(DCB_G_I_set_MCC, delam_set_array, name_TAR, color_set, marker_set, line_set, path_fig_R_MCC, range_min)


# ----------------------------------------
# Evaluation
# ----------------------------------------

# constant average value from range_x1
avg_x_MBT, avg_y_MBT = average_curve(delam_set_array, DCB_G_I_set_MBT, range_min)
print(f"Average value was computed on range: {avg_x_MBT} ")
print(f"Average value for MBT: {avg_y_MBT[0]} ")

avg_x_CC, avg_y_CC = average_curve(delam_set_array, DCB_G_I_set_CC, range_min)
print(f"Average value for CC: {avg_y_CC[0]} ")

avg_x_MCC, avg_y_MCC = average_curve(delam_set_array, DCB_G_I_set_MCC, range_min)
print(f"Average value for MCC: {avg_y_MCC[0]} ")

avg_G_list = avg_range_list(delam_set_array, DCB_G_I_set_MBT, range_min)
print(f"Average values for single sample: {avg_G_list}")
print(f"Average value from single samples: {np.mean(avg_G_list)}")

stat_eval = Statistical_evaluation(avg_G_list)
stat_eval.Check_for_outliers()
B = stat_eval.B_basis_normal()

print(f"B-basis: {B}")

# ----------------------------------------




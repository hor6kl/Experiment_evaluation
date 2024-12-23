# import load_workbook
from openpyxl import load_workbook
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import os
import math

# Importing common python functions and classes
import sys

sys.path.append(
    "/home/horakl/School/FAV/Doktorske_studium/Prace_projekty/Mechanical_testing/Python_files/DCB_CM"
)
from evaluate_enf import *
from common_evaluation import *

# ----------------------------------------

id_directory_python = os.path.dirname(os.path.abspath(__file__))

if id_directory_python.endswith("python") or id_directory_python.endswith("Python"):
    id_directory_main = id_directory_python[:-6]

id_directory_data = id_directory_main + "data"
id_directory_fig = id_directory_main + "fig"
id_directory_fig_check = os.path.join(id_directory_fig, "fig_check")
id_directory_raw_data = id_directory_main + "raw_data"

spec_data, ex_list = import_excel_values(id_directory_python)
delam_length, delam_time, strip_range = import_delam_length_values(
    id_directory_python, 1e-3
)

# ----------------------------------------

# Names in excel
width_names = ["b1", "b2", "b3"]
height_names = ["h1", "h2", "h3"]
length_names = ["lj"]
delam_names = ["ai"]

names_param = {
    "height": height_names,
    "width": width_names,
    "length": length_names,
    "delam": delam_names,
}

color_index = spec_data[0].data_name.index("color")
marker_index = spec_data[0].data_name.index("marker")
line_index = spec_data[0].data_name.index("line")

# ----------------------------------------

F = []
F_max = []
deltaY = []
name_TAR = []
time = []
voltage = []
ENF_G_II_set = []
delam_set_array = []


marker_set = []
line_set = []
color_set = []

G_II = []
indexes = []

for counter, specimen in enumerate(spec_data):
    if specimen.data_value[1] == 1 or specimen.data_value[1] == str(1):

        # path where to access data and store data
        name_TAR_str = specimen.data_value[0]
        path_TAR = os.path.join(id_directory_raw_data, (name_TAR_str + ".TRA"))
        path_fig_check = os.path.join(id_directory_fig_check, name_TAR_str)

        # reads values from .TAR file
        drop_lines = 5
        column_data = value_TAR(path_TAR, drop_lines)
        # Output should be manualy changed according columns in TAR file
        Force_i = column_data[0]
        time_i = column_data[1]
        deltaY_i = column_data[2]
        voltage_i = column_data[4]

        delam_R = delam_length[counter]
        delam_t = delam_time[counter]

        # Creating list of list.
        F.append(Force_i)
        deltaY.append(deltaY_i)
        name_TAR.append(name_TAR_str)
        time.append(time_i)
        voltage.append(voltage_i)

        color_set.append(specimen.data_value[color_index])
        line_set.append(specimen.data_value[line_index])
        marker_set.append(specimen.data_value[marker_index])

        time_delay_R = 0

        # initiaze class
        evaluate_ENF = Evaluate_ENF(
            Force_i,
            deltaY_i,
            time_i,
            delam_t,
            delam_R,
            specimen,
            names_param,
            name_TAR_str,
            voltage_i,
            strip_range[counter],
        )

        ENF_G_II, delam_set = evaluate_ENF.create_R_cuve_ENF(
            time_delay_R, path_fig_check
        )

        F_max.append(max(Force_i))
        G_II.append(evaluate_ENF.return_GII(path_fig_check))
        indexes.append(evaluate_ENF.return_voltage_indexes())

        ENF_G_II_set.append(ENF_G_II)
        delam_set_array.append(delam_set)

avrg_F_max = np.average(F_max)
print("Prumer maximalnich sil: " + ex_list + " " + str(avrg_F_max))

path_fig = os.path.join(id_directory_fig, ("graf_" + ex_list))
path_fig_R_ENF = os.path.join(id_directory_fig, ("ENF_graf_R_" + ex_list))

graf_Feps(F, deltaY, name_TAR, color_set, marker_set, line_set, path_fig, indexes)

# Converting from mm to m
dimensions = 1e-3
range_min = 60
range_min = range_min * dimensions

graph_ENF_R_curves(
    ENF_G_II_set,
    delam_set_array,
    name_TAR,
    color_set,
    marker_set,
    line_set,
    path_fig_R_ENF,
    range_min,
)

# constant average value from range_x1
avg_x, avg_y = average_curve(delam_set_array, ENF_G_II_set, range_min)
print(f"Average value was computed on range: {avg_x} ")
print(f"Average value for MBT: {avg_y[0]} ")


for G, delam in zip(ENF_G_II_set, delam_set_array):
    avg_G_list = average_on_range(delam, G, range_min)
    print(f"Average values for single sample: {avg_G_list}")


# ----------------------------------------
# Evaluation
# ----------------------------------------

# constant average value from range_x1
stat_eval = Statistical_evaluation(G_II)
stat_eval.Check_all()
# ----------------------------------------
# ----------------------------------------

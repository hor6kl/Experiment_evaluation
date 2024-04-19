# ----------------------------------------
# 
# ----------------------------------------
# import load_workbook
from openpyxl import load_workbook
import numpy as np
import matplotlib.pyplot as pl
import os
import math

# Load common library with evaluation funcitons
import sys
sys.path.append('/home/horakl/School/FAV/Doktorske_studium/Prace_projekty/Mechanical_testing/Python_files/DCB_CM')
from common_DCB_evaluation import *

import pandas as pd
# ----------------------------------------

# ---------------------------------
# class to strucutre data from experiments excel where dimensions of specimens are saved
class Specimen_data:
    def __init__(self, data_name, data_value, data_dimension):
        range = len(data_name)
        self.data_name = data_name
        self.data_value = data_value
        self.data_dimension = data_dimension
# ---------------------------------

separator = "/"  # linux "/", Windows "\\"

# Retruns path of this python file
id_directory_python = os.path.dirname(os.path.abspath(__file__))

# Returns path to folder structure
if id_directory_python.endswith("python") or id_directory_python.endswith("Python"):
    id_directory_main = id_directory_python[:-6]
id_directory_data = id_directory_main + "data"
id_directory_fig = id_directory_main + "fig"
id_directory_raw_data = id_directory_main + "raw_data"


# set file path
file_name = "experiments.xlsx"
file_name_R = 'experiments_R_curves.csv'
file_path_name = id_directory_data + separator + file_name
file_path_name = id_directory_data + separator + file_name
file_path_name_R = id_directory_data + separator + file_name_R
# nacte xlxs
wb = load_workbook(file_path_name)

# nacte listy
sheets_all = wb.sheetnames

print("seznam listu v xlxs: " + str(sheets_all))
print("-------------------------------------------")

# cyklus pro listy NEAKTIVNI
for list in sheets_all:
    print("prace s: " + str(list))
    sheet = wb[list]


sheet = wb.active

rows = [[cell.value for cell in row] for row in sheet.rows]
cols = [[cell.value for cell in column] for column in sheet.columns]

Quantity_TRA_row = rows[cols[0].index('Quantity_TRA')]
Quantity_TRA = Quantity_TRA_row[1:Quantity_TRA_row.index(None)]

Quantity_row = rows[cols[0].index('Quantity')]
Values_name = Quantity_row[0:Quantity_row.index(None) - 1]
print(Values_name)

# pocet sloupcu pro data
j1 = len(Values_name) + 1

SI_row = rows[cols[0].index('SI')]
SI = SI_row[0:SI_row.index(None)]
print(SI)


#  najde radek s daty -1
i2 = 1
while sheet.cell(row=i2, column=1).value != '#':
    i2 += 1
Values_row = i2 + 1



# ----------------------------------------

# Names in excel 
width_names = ["W1", "W2", "W3"]
height_names = ["H1", "H2", "H3"]
length_names = ["lj"]

names_param = {'height': height_names, 'width': width_names, 'length': length_names}

color_index = Values_name.index('color')
marker_index = Values_name.index('marker')
line_index = Values_name.index('line')

# reads lines in excel 
spec_data = []
pocet = 0
while (sheet.cell(row=Values_row + pocet, column=1).value) != None:

    # nacte hodnoty
    Values = []
    for i3 in range(1, j1):
        hodnota_temp = sheet.cell(row=Values_row + pocet, column=i3).value
        Values.append(hodnota_temp)

    # pole jmen a hodnot z xlsx
    spec_data.append(Specimen_data(Values_name, Values, SI))

    pocet += 1

# ----------------------------------------
F = []
delta_y = []
delta_x = []
time = []

name_TAR = []


marker_set = []
line_set = []
color_set = []

for specimen in spec_data:
    if specimen.data_value[1] == 1 or specimen_data[cd].data_value[1] == str(1):


        # path where to access data and store data
        name_TAR_str = specimen.data_value[0]
        path_TAR = id_directory_raw_data + separator + name_TAR_str + ".TRA"
        path_fig = id_directory_fig + separator + 'graf_' + list

        # reads values from .TAR file 
        drop_lines = 6
        column_data = value_TAR(path_TAR, drop_lines)
        # Output should be manualy changed according columns in TAR file
        F1 = column_data[0]
        time1 = column_data[1]
        delta_y1 = column_data[4]
        delta_x1 = column_data[5]
        
        # Creating list of list.
        F.append(F1)
        delta_y.append(delta_y1)
        name_TAR.append(name_TAR_str)
        time.append(time1)

        color_set.append(specimen.data_value[color_index])
        line_set.append(specimen.data_value[line_index])
        marker_set.append(specimen.data_value[marker_index])

        # Evaluating Young modulus
        Eval_tensile = Evaluate_Tensile(F1, delta_y1, specimen, names_param, path_fig)
        Eval_tensile.Evaluate_Young_modulus()
        Eval_tensile.Evaluate_Poisson(delta_x1)



avrg_F_max = np.average(F_max)
print("Prumer maximalnich sil: " + list + " " + str(avrg_F_max))

path_fig = id_directory_fig + separator + 'graf_' + list
path_fig_R_ENF = id_directory_fig + separator + 'ENF_graf_R_' + list

graf_Feps(F, delta_y, name_TAR, color_set, marker_set, line_set, path_fig)



# ----------------------------------------

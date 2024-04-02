
# import load_workbook
from openpyxl import load_workbook
import numpy as np
import matplotlib.pyplot as pl
import os
import math

import sys
sys.path.append('/home/horakl/School/FAV/Doktorske_studium/Prace_projekty/Mechanical_testing/Python_files/DCB_CM')
from common_DCB_evaluation import *

import pandas as pd

# if url.endswith('.com'):

separator = "/"  # linux "/", Windows "\\"

id_directory_python = os.path.dirname(os.path.abspath(__file__))
# print(id_directory_python)

if id_directory_python.endswith("python") or id_directory_python.endswith("Python"):
    id_directory_main = id_directory_python[:-6]
# print(id_directory_main)
id_directory_data = id_directory_main + "data"
id_directory_fig = id_directory_main + "fig"
id_directory_raw_data = id_directory_main + "raw_data"


# set file path
file_name = "experiments.xlsx"
file_name_R = 'experiments_R_curves.csv'
file_path_name = id_directory_data + "\\" + file_name
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

# najde radku s jmeny velicin
# i1 = 1
# while sheet.cell(row=i1, column=1).value != 'Quantity':
#     i1 += 1


#  najde radek s daty -1
i2 = 1
while sheet.cell(row=i2, column=1).value != '#':
    i2 += 1
Values_row = i2 + 1


# #  najde radek s nasobky velicin
# i4 = 1
# while sheet.cell(row=i4, column = 1).value != 'SI':
#     i4 += 1

# # ---------------------------------
class Data1:
    def __init__(self, data_name, data_value, data_dimension):
        range = len(data_name)
        self.data_name = data_name
        self.data_value = data_value
        self.data_dimension = data_dimension
# # ---------------------------------


B_index = Values_name.index("B")

#h0_index = Values_name.index("2h0")
hA_index = Values_name.index("2hA")
#hB_index = Values_name.index("2hB")

ai_index = Values_name.index("ai")

lj_index = Values_name.index("lj")

color_index = Values_name.index('color')
marker_index = Values_name.index('marker')
line_index = Values_name.index('line')

# cte radky
a = []
pocet = 0
while (sheet.cell(row=Values_row + pocet, column=1).value) != None:

    # nacte hodnoty
    Values = []
    for i3 in range(1, j1):
        hodnota_temp = sheet.cell(row=Values_row + pocet, column=i3).value
        Values.append(hodnota_temp)

# pole jmen a hodnot z xlsx
    a.append(Data1(Values_name, Values, SI))

    pocet += 1

# -------------------- graf


# ----------------------------------------


# parametry pro R-curves
R_sheet = pd.read_csv(file_path_name_R)

time_R_all = R_sheet.iloc[6, 7:]
delam_a0_R = R_sheet.iloc[7:, 6]
time_delay_R = R_sheet.iloc[7:, 2]

# Convering from mm to m
dimensions = 1e-3 
delam_a0_R = [float(x)*dimensions for x in delam_a0_R]

#
F = []
epsilon = []
name_TAR = []
epsilon_n = []
time = []
F_max = []
Energ_crit_array = []
cons_1 = 0


marker_set = []
line_set = []
color_set = []
DCB_G_I_set_MBT = []
DCB_G_I_set_CC = []
DCB_G_I_set_MCC = []
delam_set_array = []

for cd in range(0, len(a)):
    if a[cd].data_value[1] == 1:

        delam_R = R_sheet.iloc[7+cd, 7:]        

        # cleaning list from 'nan' values
        delam_R = [float(x) for x in delam_R if not isinstance(x, float) or not math.isnan(x)]
        time_R = time_R_all[0:len(delam_R)]

        # cleaning list from zero values at the begining exceprt for one
        for counter in range(0, len(delam_R)-1):
            if not (delam_R[counter] == 0.0 and delam_R[counter+1] == 0.0):
                break
        delam_R = delam_R[counter:len(delam_R)]
        time_R = time_R[counter:len(time_R)]

        # Convering from mm to m
        delam_R = [x*dimensions for x in delam_R]


        # path where to access data and store data
        name_TAR_str = a[cd].data_value[0]
        path_TAR = id_directory_raw_data + separator + name_TAR_str + ".TRA"
        path_fig = id_directory_fig + separator + name_TAR_str + ".TRA"

        F1, epsilon1, name_TAR1, epsilon_n1, time1 = value_TAR(path_TAR)
        
        # Creating list of list.
        F.append(F1)
        epsilon.append(epsilon1)
        name_TAR.append(name_TAR1)
        epsilon_n.append(epsilon_n1)
        time.append(time1)
        color_set.append(a[cd].data_value[color_index])
        line_set.append(a[cd].data_value[line_index])
        marker_set.append(a[cd].data_value[marker_index])

#        F_np = np.array(F)
        F_np = F

        if "DCB" in a[cd].data_value[0]:
            Energ_crit = vypocet_energie_DCB(
                F_np[cons_1], a[cd], epsilon[cons_1])
            DCB_G_I_MBT, delam_set = create_R_cuve(time[cons_1], time_R, delam_R, F_np[cons_1], a[cd], epsilon[cons_1], delam_a0_R[cd], time_delay_R.iloc[cd], path_fig, "MBT")
            DCB_G_I_CC, delam_set = create_R_cuve(time[cons_1], time_R, delam_R, F_np[cons_1], a[cd], epsilon[cons_1], delam_a0_R[cd], time_delay_R.iloc[cd], path_fig, "CC")
            DCB_G_I_MCC, delam_set = create_R_cuve(time[cons_1], time_R, delam_R, F_np[cons_1], a[cd], epsilon[cons_1], delam_a0_R[cd], time_delay_R.iloc[cd], path_fig, "MCC")

        elif "ENF" in a[cd].data_value[0]:
            Energ_crit = vypocet_energie_ENF(
                F_np[cons_1], a[cd], epsilon[cons_1])
            DCB_G_I, delam_set = create_R_cuve(
                time[cons_1], time_R, delam_R, F_np[cons_1], a[cd], epsilon[cons_1], delam_a0_R.iloc[cd], time_delay_R.iloc[cd])
        elif "MIX" in a[cd].data_value[0]:
            Energ_crit = vypocet_energie_MixI_II(F_np[cons_1], a[cd])
        print("crit energie pro " + name_TAR[cons_1] + " " + str(Energ_crit))
        Energ_crit_array.append(Energ_crit)
        F_max.append(max(F_np[cons_1]))

        DCB_G_I_set_MBT.append(DCB_G_I_MBT)
        DCB_G_I_set_CC.append(DCB_G_I_CC)
        DCB_G_I_set_MCC.append(DCB_G_I_MCC)
        delam_set_array.append(delam_set)

        cons_1 += 1


avrg_F_max = np.average(F_max)
print("Prumer maximalnich sil: " + list + " " + str(avrg_F_max))

avrg_energ_crit = np.average(Energ_crit_array)
print("Prumer krit energie: " + list + " " + str(avrg_energ_crit))

var_koef_F_max = Variacni_koef(Energ_crit_array)
print("Variacni koeficient pro " + list + " " + str(var_koef_F_max))

path_fig = id_directory_fig + separator + 'graf_' + list
path_fig_R_MBT = id_directory_fig + separator + 'MBT_graf_R_' + list
path_fig_R_CC = id_directory_fig + separator + 'CC_graf_R_' + list
path_fig_R_MCC = id_directory_fig + separator + 'MCC_graf_R_' + list

graf_Feps(F, epsilon, name_TAR, color_set, marker_set, line_set, path_fig)

graf_Feps_R_curves(DCB_G_I_set_MBT, delam_set_array, name_TAR, color_set, marker_set, line_set, path_fig_R_MBT)
graf_Feps_R_curves(DCB_G_I_set_CC, delam_set_array, name_TAR, color_set, marker_set, line_set, path_fig_R_CC)
graf_Feps_R_curves(DCB_G_I_set_MCC, delam_set_array, name_TAR, color_set, marker_set, line_set, path_fig_R_MCC)


# F_average, temp2 = value_TAR_avereged(F)
# F_max = vypocet_energie(F_average)

# print("kriticka deformacni energie pro prumer:")
# print(F_max)


# poznamka displacement, green vysledna krivka, red MKP analyza

#


# ----------------------------------------


# ----------------------------------------

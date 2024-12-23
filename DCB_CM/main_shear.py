import os

# Importing common python functions and classes
import sys

sys.path.append(
    "/home/horakl/School/FAV/Doktorske_studium/Prace_projekty/Mechanical_testing/Python_files/DCB_CM"
)
from evaluate_shear import *
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

# ----------------------------------------

# Names in excel
width_names = ["W1"]
height_names = ["H1"]
thickness_names = ["t1", "t2"]

names_param = {
    "height": height_names,
    "width": width_names,
    "thickness": thickness_names,
}

color_index = spec_data[0].data_name.index("color")
marker_index = spec_data[0].data_name.index("marker")
line_index = spec_data[0].data_name.index("line")

# ----------------------------------------
F = []
deltaY = []
deltaX = []
time = []

name_TAR = []

marker_set = []
line_set = []
color_set = []

stress = []
strain = []

max_stress_list = []

for specimen in spec_data:
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

        # Creating list of list.
        F.append(Force_i)
        deltaY.append(deltaY_i)
        name_TAR.append(name_TAR_str)
        time.append(time_i)

        color_set.append(specimen.data_value[color_index])
        line_set.append(specimen.data_value[line_index])
        marker_set.append(specimen.data_value[marker_index])

        # Evaluating Young modulus
        Eval_tensile = Evaluate_Shear(
            Force_i, deltaY_i, specimen, names_param, path_fig_check
        )
        Eval_tensile.Evaluate_Young_modulus()

        stress.append(Eval_tensile.Calculate_stress())
        strain.append(Eval_tensile.Calculate_strain())

        max_stress_list.append(Eval_tensile.return_max_stress())

path_fig = os.path.join(id_directory_fig, ("graf_" + ex_list))
path_fig_ss = os.path.join(id_directory_fig, ("graf_SS_" + ex_list))

graf_Feps(F, deltaY, name_TAR, color_set, marker_set, line_set, path_fig)

graph_stress_strain_curves(
    strain, stress, name_TAR, color_set, marker_set, line_set, path_fig_ss
)

# ----------------------------------------
# Evaluation
# ----------------------------------------
stat_eval = Statistical_evaluation(max_stress_list)
stat_eval.Check_all()

# ----------------------------------------

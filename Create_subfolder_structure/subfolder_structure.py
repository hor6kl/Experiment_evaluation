#
# Function creates subfolder structure defined in subfolder_name tuple.
# As a input argument in console, main folder path should be passed
#

import sys
import os

def create_subfolder(folder_path: str, subfolder_name: str) -> None:
    # Create the subfolder inside the main folder
    subfolder_path = os.path.join(folder_path, subfolder_name)
    if not os.path.exists(subfolder_path):
        #os.makedirs(subfolder_path)
        print(f"Main folder: '{folder_path}' does not exist!")

    print(f"Subfolder '{subfolder_name}' created successfully inside '{folder_path}'!")

# Check if the correct number of parameters are passed
if len(sys.argv) != 2:
    print("Please input one argument with main folder path")
    print("Usage: python subfolder_structure.py <folder_path>")
    sys.exit(1)

# Get the folder path name from the command line arguments
folder_path = sys.argv[1]
subfolder_name = ['python', 'photo', 'raw_data', ['fig', 'fig_check'], 'data']


# Loop to iterate through every subfolder name
for subfolder in subfolder_name:

    # logic to create subsub folders
    if type(subfolder) == list:
       for counter, subsub in enumerate(subfolder):
           # creates sub folder on first postition in sub list
            if counter == 0:
                subfolder_path = folder_path + "//" + subsub
                create_subfolder(folder_path, subsub)
            else:
                create_subfolder(subfolder_path, subsub)

    else:
        # Call the function with the folder path and subfolder name
        create_subfolder(folder_path, subfolder)




import cv2
import numpy as np
import tkinter as tk
from tkinter import messagebox
import threading  # Import the threading module

import csv

# import OS module
import os



class PixelSelector:
    def __init__(self):
        self.selected_sx = []
        self.selected_sy = [] 
        self.selected_x = []
        self.selected_y = [] 
        self.output_x = []
        self.output_y = [] 

        self.first_img = False

    def get_pixel_coordinates(self, event, x, y, flags, param):

        # setting scale 
        if event == cv2.EVENT_LBUTTONDOWN and len(self.selected_sx) < 2:
            self.selected_sx.append(x)
            self.selected_sy.append(y)

            image[y, x] = [255,0,0]
            cv2.imshow('Select Pixel', image)

            print(f"Selected pixel scale coordinates: [{self.selected_sx}, {self.selected_sy}]")

        # logic for finding midlle pixel  
        elif event == cv2.EVENT_LBUTTONDOWN and len(self.output_x) == 0:
            self.selected_x.append(x)
            self.selected_y.append(y)

            image[y, x] = [255,0,0]
            cv2.imshow('Select Pixel', image)

            print(f"Selected pixel coordinates: [{self.selected_x}, {self.selected_y}]")
            print(f"value of first img {self.first_img}")

        elif event == cv2.EVENT_LBUTTONDOWN:
            self.output_x.append(x)
            self.output_y.append(y)

            image[y, x] = [255,0,0]
            cv2.imshow('Select Pixel', image)

            print(f"Selected pixel output coordinates: [{self.selected_x}, {self.selected_y}]")

    def get_middle_pixel(self):
        if len(self.selected_x) == 2:

            x1 = self.selected_x[0]
            x2 = self.selected_x[1]
            y1 = self.selected_y[0]
            y2 = self.selected_y[1]
        
            x = (x2-x1)/2 + x1
            y = (y2-y1)/2 + y1

            self.output_x.append(x) 
            self.output_y.append(y)

            image[int(y), int(x)] = [0, 255, 0]  # Color the pixel green

            cv2.imshow('Select Pixel', image)
            print(f"middle coardinates: {x} {y}")


    def select_pixel(self, image_path):
        global image
        image = cv2.imread(image_path)
        cv2.namedWindow('Select Pixel', cv2.WINDOW_NORMAL)  # Create a resizable window
        cv2.resizeWindow('Select Pixel', 800, 600)  # Resize the window
        cv2.imshow('Select Pixel', image)
        cv2.setMouseCallback('Select Pixel', self.get_pixel_coordinates)


        if self.first_img == False:
            #This function with first image and setting scale
            while len(self.selected_x) < 2: 
                if cv2.waitKey(1) & 0xFF == ord('q'):
                    break
                elif cv2.waitKey(1) & 0xFF == ord('p') and self.first_img == False:
                    return self.selected_sx, self.selected_sy, None, None
                elif len(self.output_x) == 2:
                    return self.selected_sx, self.selected_sy, self.output_x, self.output_y                

        else:
            #This while loop functions with update of middle pixel
            while len(self.selected_x) < 2 or self.first_img == False:
                if cv2.waitKey(1) & 0xFF == ord('q'):
                    break
                elif cv2.waitKey(1) & 0xFF == ord('p') and self.first_img == False:
                    return self.selected_sx, self.selected_sy, None, None
                elif cv2.waitKey(1) & 0xFF == ord('p') and self.first_img == True:
                    return None, None, None, None                
                elif len(self.output_x) == 2:
                    return self.selected_sx, self.selected_sy, self.output_x, self.output_y                


        self.get_middle_pixel()
        cv2.imshow('Select Pixel', image)
        cv2.waitKey(0)

        cv2.destroyAllWindows()
        return self.selected_sx, self.selected_sy, self.output_x, self.output_y


    def set_first_img(self):
        self.first_img = True 

    def reset_values(self):
        self.selected_x = []
        self.selected_y = [] 
        self.output_x = []
        self.output_y = [] 


class Instruction:
    
    def __init__(self):
        self.selected_x = False 
        self.selected_y = False 

    def show_instructions(self):
        root = tk.Tk()
        root.title("Instructions")
    
        tk.Label(root, text="Instructions:").pack()
    
        check_var_1 = tk.IntVar()

        tk.Checkbutton(root, text="someting", variable = check_var_1, state='disabled').pack()
        tk.Checkbutton(root, text="else", state='disabled').pack()


        root.mainloop()
    
    
    def update_checkboxes(self):
        step1_text = "Step 1: Select the first pixel"
        step2_text = "Step 2: Select the second pixel"
    
        if len(selected_x) >= 1:
            x1 = selected_x[0]
            y1 = selected_y[0]
            step1_text = f"Step 1: Selected Pixel - ({x1}, {y1})"
    
        if len(selected_x) == 2:
            x2 = selected_x[1]
            y2 = selected_y[1]
            step2_text = f"Step 2: Selected Pixel - ({x2}, {y2})"


class MyTkApp(threading.Thread):
    def __init__(self):
        self.root=tk.Tk()
        self.s = tk.StringVar()
        self.s.set('Foo')
        l = tk.Label(self.root,textvariable=self.s)
        l.pack()
        threading.Thread.__init__(self)

    def run(self):
        self.root.mainloop()




# driver function 
if __name__=="__main__": 

    ####
    # Input parameters
    ####
    path_photo = '/home/horakl/School/FAV/Doktorske_studium/Prace_projekty/Mechanical_testing/Python_files/Distance_in_image/Input_folder'

    # Get the list of all files and directories
    dir_list = os.listdir(path_photo)
    print(dir_list)


    # Create an instance of Instructions 
    instructions = Instruction()

    # Show instructions pop-up window in a separate thread
    instructions_thread = threading.Thread(target=instructions.show_instructions)
    instructions_thread.start()

    # Create an instance of PixelSelector
    pixel_selector = PixelSelector()


    for counter, file_name in enumerate(dir_list):

        file_path_name = path_photo + '/' + file_name

        # Call select_pixel method to select a pixel on the image
        selected_sx, selected_sy, output_x, output_y = pixel_selector.select_pixel(file_path_name)
        
        pixel_selector.reset_values()
        if counter == 0:
            pixel_selector.set_first_img()
            # scale value
            scale_value = 10 # in mm
        
            scale_x1 = selected_sx[0]
            scale_x2 = selected_sx[1]
            scale_y1 = selected_sy[0]
            scale_y2 = selected_sy[1]
            
            scale_pixel_dist = np.sqrt((scale_x1 - scale_x2) ** 2 + (scale_y1 - scale_y2) ** 2)

            print(f"Scale value was set as {scale_value} mm")
            print(f"Pixel scale was computed as {scale_pixel_dist} divided by scale value")
    
    
    
        if output_x == None:
            distance_ = 0
        else:
        
            x1 = output_x[0]
            x2 = output_x[1]
            y1 = output_y[0]
            y2 = output_y[1]
        
            pixel_distance = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        
            distance_ = pixel_distance* (scale_value/scale_pixel_dist)

            
        print(f"Distance between selected points: {distance_} mm. Evaluated {file_name}.")
    
        path = '/home/horakl/School/FAV/Doktorske_studium/Prace_projekty/Mechanical_testing/Python_files/Distance_in_image'
        name_file = '/output.csv'
    
        path_name = path + name_file
        header = ['name', 'distance']
    
        # open the file in the write mode
        with open(path_name, 'a') as f:
            # create the csv writer
            writer = csv.writer(f)
    
            if counter == 0:
                # write the header
                writer.writerow(header)
    
            # write a row to the csv file
            writer.writerow([distance_, file_name])

    


    





















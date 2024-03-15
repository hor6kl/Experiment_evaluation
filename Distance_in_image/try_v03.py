import cv2
import numpy as np
import tkinter as tk
from tkinter import messagebox
import threading  # Import the threading module

import csv

# import OS module
import os

import queue


# Shared queue
shared_q_pass = queue.Queue()
shared_q_quit = queue.Queue()
shared_q_next = queue.Queue()
shared_q_back = queue.Queue()


class PixelSelector(threading.Thread):
    def __init__(self, dirt_list, path_photo, app):
        super(PixelSelector, self).__init__()
        self.selected_sx = []
        self.selected_sy = [] 
        self.selected_x = []
        self.selected_y = [] 
        self.output_x = []
        self.output_y = [] 


        self.first_img = False
        
        self.app = app
        self.dirt_list = dirt_list
        self.path_photo = path_photo

        #call function

    def run(self):
        while shared_q_quit.empty():
            self.evaluate_distance(self.dirt_list, self.path_photo)

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

#            print(f"Selected pixel coordinates: [{self.selected_x}, {self.selected_y}]")
#            print(f"value of first img {self.first_img}")

        elif event == cv2.EVENT_LBUTTONDOWN:
            self.output_x.append(x)
            self.output_y.append(y)

            image[y, x] = [255,0,0]
            cv2.imshow('Select Pixel', image)

#            print(f"Selected pixel output coordinates: [{self.selected_x}, {self.selected_y}]")



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
                if cv2.waitKey(1) & 0xFF == ord('q') or (not shared_q_quit.empty()):
                    break
                elif (cv2.waitKey(1) & 0xFF == ord('n') and self.first_img == False) or ( not shared_q_next.empty()):
                    if not shared_q_next.empty():
                        shared_q_next.get()
                    return self.selected_sx, self.selected_sy, None, None
                elif len(self.output_x) == 2:
                    return self.selected_sx, self.selected_sy, self.output_x, self.output_y                
                elif not shared_q_next.empty():
                    print('rajska')
                    shared_q_next.get()


        else:
            #This while loop functions with update of middle pixel
            while len(self.selected_x) < 2 or self.first_img == False:
                if cv2.waitKey(1) & 0xFF == ord('q'):
                    break
                elif (cv2.waitKey(1) & 0xFF == ord('n') and self.first_img == False) or ( not shared_q_next.empty()):
                    if not shared_q_next.empty():
                        shared_q_next.get()
                    return self.selected_sx, self.selected_sy, None, None
                elif (cv2.waitKey(1) & 0xFF == ord('p') and self.first_img == True) or ( not shared_q_pass.empty()):
                    if not shared_q_pass.empty():
                        shared_q_pass.get()
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


    def evaluate_distance(self, dirt_list, path_photo):

        for counter, file_name in enumerate(dir_list):
    
            file_path_name = path_photo + '/' + file_name
    
            # Call select_pixel method to select a pixel on the image
            selected_sx, selected_sy, output_x, output_y = self.select_pixel(file_path_name)
    
    
            self.reset_values()
            if counter == 0:
                self.set_first_img()
            
                scale_x1 = selected_sx[0]
                scale_x2 = selected_sx[1]
                scale_y1 = selected_sy[0]
                scale_y2 = selected_sy[1]
                
                scale_pixel_dist = np.sqrt((scale_x1 - scale_x2) ** 2 + (scale_y1 - scale_y2) ** 2)
    
                print(f"Scale value was set as {scale_value} mm")
                print(f"Pixel scale was computed as {scale_pixel_dist} divided by scale value")
        
                self.app.pixel_scale.set(self.app.pixel_scale.get() + ' = ' + str(scale_pixel_dist) + ' px')
                self.app.check_var_1.set(1)
        
            if output_x == None:
                distance_ = 0
            else:
            
                x1 = output_x[0]
                x2 = output_x[1]
                y1 = output_y[0]
                y2 = output_y[1]
            
                pixel_distance = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
            
                distance_ = pixel_distance* (scale_value/scale_pixel_dist)
    
                
            print(f'Distance between selected points: {distance_} mm. Evaluated {file_name}. {counter}{len(dir_list)}')
        
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





class Application(tk.Frame):
    def __init__(self, root=()):
        super().__init__(root)
        self.root = root 
        self.pack()
        self.check_var_1 = tk.IntVar()
        self.check_var_2 = tk.IntVar()
        self.check_var_3 = tk.IntVar()

        self.pixel_scale = tk.StringVar() 
        self.create_widgets()

    def create_widgets(self):
        self.root.title('Pixel Selector Instructions')

        self.pixel_scale.set('Scale value was set: 10 mm')

        tk.Checkbutton(self.root, text="1# Select pixels for scale", variable = self.check_var_1, state='disabled').pack()

        var = tk.StringVar()
        label = tk.Label( self.root, textvariable=self.pixel_scale)
        label.pack()

        tk.Checkbutton(self.root, text="2# Selct pixels on adherends where initial delamination was marked", variable = self.check_var_2, state='disabled').pack()
        tk.Checkbutton(self.root, text="3# Selct pixel in the end of the crack", variable = self.check_var_3, state='disabled').pack()

        self.quit_button = tk.Button(self.root, text="Pass (p)", fg="black",
                              command=self.on_pass_button_click)
        self.quit_button.pack(side="left")
        self.quit_button = tk.Button(self.root, text="Next (n)", fg="blue",
                              command=self.on_next_button_click)
        self.quit_button.pack(side="left")
        self.quit_button = tk.Button(self.root, text="Back (b)", fg="blue",
                              command=self.on_back_button_click)
        self.quit_button.pack(side="left")
        self.quit_button = tk.Button(self.root, text="Quit (q)", fg="red",
                              command=self.on_quit_button_click)
        self.quit_button.pack(side="left")

    def on_pass_button_click(self):
        shared_q_pass.put(True)
    def on_next_button_click(self):
        shared_q_next.put(True)
    def on_back_button_click(self):
        shared_q_back.put(True)

    def on_quit_button_click(self):
        shared_q_quit.put("Quit")
        self.root.quit()
        self.root.destroy()




# driver function 
if __name__=="__main__": 

    ####
    # Input parameters
    ####
    path_photo = '/home/horakl/School/FAV/Doktorske_studium/Prace_projekty/Mechanical_testing/Python_files/Distance_in_image/Input_folder'

    # scale value
    scale_value = 10 # in mm



    # Get the list of all files and directories
    dir_list = os.listdir(path_photo)
    print(dir_list)


    root = tk.Tk()
    instruction_widget = Application(root=root)
    thread_pixel = PixelSelector(dir_list, path_photo, instruction_widget) 

    # First is started the pixel selector in thread. Then instruction widget is started in main loop.
    thread_pixel.start() 
    instruction_widget.mainloop()  # Start the GUI thread

    







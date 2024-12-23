# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: horakl
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Inputs:
#           path_photo  - path to folder where all images are loaded
#           scale_value - float value which corresponds with scale which is choosen in first step
#           output_path - path where to output csv file with all values
#
#   Outputs:
#           output.csv - file with all measurements and names of images
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Description:
#           Script loads all images in specified folder one by one images are loaded.
#           1. Step is to choose scale. Two pixels are selected for this reason. After selection it is possilbe to Pass (P) or continue with selection and then press Next (N).
#           2. Next steps are without selecting scale. Select two pixels between them middle pixel is created and then select last pixel. Distance is computed between midlle pixel and last selected pixel.
#           3. Step 2. is repeated for every image in folder.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Improvements:
#           Add functionality to choose if middle pixel should be created or not
#           Check button to choose between ENF/DCB
#           Pop up help where functionality of every key in every scenario
#           Proper close of model when all images are evaluated and proper quit without erro
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    def __init__(self, dirt_list, path_photo, app, output_path):
        super(PixelSelector, self).__init__()
        self.selected_sx = []
        self.selected_sy = []
        self.selected_x = []
        self.selected_y = []
        self.output_x = []
        self.output_y = []

        self.history = {"BGR": [], "pos": []}

        self.first_img = False

        self.app = app
        self.dirt_list = dirt_list
        self.output_path = output_path
        self.path_photo = path_photo

    # call function
    def run(self):
        self.evaluate_distance(self.dirt_list, self.path_photo)

    def get_color_array(
        self, y: "position in y direction", x: "position in x direction"
    ) -> list[int, int, int]:
        B = image[y, x, 0]
        G = image[y, x, 1]
        R = image[y, x, 2]
        color = [B, G, R]

        return color

    def get_pixel_coordinates(self, event, x, y, flags, param) -> None:

        # setting scale
        if event == cv2.EVENT_LBUTTONDOWN and len(self.selected_sx) < 2:
            self.selected_sx.append(x)
            self.selected_sy.append(y)

            color = self.get_color_array(y, x)
            self.history["BGR"].append(color)
            self.history["pos"].append([x, y])

            print(self.history)

            image[y, x] = [255, 0, 0]
            cv2.imshow("Select Pixel", image)

            print(
                f"Selected pixel scale coordinates: [{self.selected_sx}, {self.selected_sy}]"
            )

        # logic for finding midlle pixel
        elif event == cv2.EVENT_LBUTTONDOWN and len(self.output_x) == 0:
            self.selected_x.append(x)
            self.selected_y.append(y)

            color = self.get_color_array(y, x)
            self.history["BGR"].append(color)
            self.history["pos"].append([x, y])

            image[y, x] = [255, 0, 0]
            cv2.imshow("Select Pixel", image)

        #            print(f"Selected pixel coordinates: [{self.selected_x}, {self.selected_y}]")
        #            print(f"value of first img {self.first_img}")

        elif event == cv2.EVENT_LBUTTONDOWN:
            self.output_x.append(x)
            self.output_y.append(y)

            color = self.get_color_array(y, x)
            self.history["BGR"].append(color)
            self.history["pos"].append([x, y])

            image[y, x] = [255, 0, 0]
            cv2.imshow("Select Pixel", image)

            self.app.check_var_3.set(1)

    def select_pixel(self, image_path):
        global image
        image = cv2.imread(image_path)
        cv2.namedWindow("Select Pixel", cv2.WINDOW_NORMAL)  # Create a resizable window
        cv2.resizeWindow("Select Pixel", 1680, 1050)  # Resize the window
        cv2.imshow("Select Pixel", image)
        cv2.setMouseCallback("Select Pixel", self.get_pixel_coordinates)

        # This part of code would be good to re-do....
        if self.first_img == False:
            # This function with first image and setting scale
            while len(self.selected_x) < 2:
                keypress = cv2.waitKey(1)
                if keypress & 0xFF == ord("q") or (not shared_q_quit.empty()):
                    break
                elif keypress & 0xFF == ord("p") or (not shared_q_pass.empty()):
                    if not shared_q_pass.empty():
                        shared_q_pass.get()
                    return self.selected_sx, self.selected_sy, None, None
                elif keypress & 0xFF == ord("n") or (not shared_q_next.empty()):
                    print("next1")
                    if not shared_q_next.empty():
                        shared_q_next.get()
                    return (
                        self.selected_sx,
                        self.selected_sy,
                        self.output_x,
                        self.output_y,
                    )
                elif keypress & 0xFF == ord("b") or (not shared_q_back.empty()):
                    print("back1")
                    self.history_recall()
                    if not shared_q_back.empty():
                        shared_q_back.get()
                    continue
                if len(self.selected_x) == 2 and len(self.output_x) == 0:
                    self.get_middle_pixel()
                    cv2.imshow("Select Pixel", image)

        else:
            # This while loop functions with update of middle pixel
            while True:
                keypress = cv2.waitKey(1)
                if keypress & 0xFF == ord("q") or (not shared_q_quit.empty()):
                    break
                elif keypress & 0xFF == ord("n") or (not shared_q_next.empty()):
                    if not shared_q_next.empty():
                        shared_q_next.get()
                    return (
                        self.selected_sx,
                        self.selected_sy,
                        self.output_x,
                        self.output_y,
                    )
                elif keypress & 0xFF == ord("p") or (not shared_q_pass.empty()):
                    if not shared_q_pass.empty():
                        shared_q_pass.get()
                    return None, None, None, None
                elif keypress & 0xFF == ord("b") or (not shared_q_back.empty()):
                    self.history_recall()
                    if not shared_q_back.empty():
                        shared_q_back.get()
                    continue
                if len(self.selected_x) == 2 and len(self.output_x) == 0:
                    self.get_middle_pixel()
                    cv2.imshow("Select Pixel", image)

    def history_recall(self) -> None:
        # Function is retrieving values from history dictionary and deleting values from list with cooarinates

        try:
            position = self.history["pos"].pop()
            color = self.history["BGR"].pop()
            image[position[1], position[0]] = color
            cv2.imshow("Select Pixel", image)

            if self.first_img:
                if len(self.selected_x) == 2 and len(self.output_x) < 2:
                    position = self.history["pos"].pop()
                    color = self.history["BGR"].pop()
                    image[position[1], position[0]] = color

                    self.selected_x.pop()
                    self.selected_y.pop()
                    self.output_x.pop()
                    self.output_y.pop()

                    # Unchecking check box for step 2
                    self.app.check_var_2.set(0)
                elif len(self.selected_x) < 2:
                    self.selected_x.pop()
                    self.selected_y.pop()
                elif len(self.output_x) == 2:
                    self.output_x.pop()
                    self.output_y.pop()
                    # Unchecking check box for step 3
                    self.app.check_var_3.set(0)

            else:
                self.selected_sx.pop()
                self.selected_sy.pop()

            print("Back function executed.")

        except:
            print("Maximum number of back opperations reached!")

        cv2.imshow("Select Pixel", image)

    def set_first_img(self):
        self.first_img = True

    def reset_values(self):
        self.selected_x.clear()
        self.selected_y.clear()
        self.output_x.clear()
        self.output_y.clear()
        self.app.check_var_2.set(0)
        self.app.check_var_3.set(0)

    def get_middle_pixel(self):
        if len(self.selected_x) == 2:

            x1 = self.selected_x[0]
            x2 = self.selected_x[1]
            y1 = self.selected_y[0]
            y2 = self.selected_y[1]

            x = (x2 - x1) / 2 + x1
            y = (y2 - y1) / 2 + y1

            x = int(x)
            y = int(y)

            self.output_x.append(x)
            self.output_y.append(y)

            color = self.get_color_array(y, x)
            self.history["BGR"].append(color)
            self.history["pos"].append([x, y])

            image[y, x] = [0, 255, 0]  # Color the pixel green

            cv2.imshow("Select Pixel", image)
            print(f"middle coardinates: {x} {y}")
            self.app.check_var_2.set(1)

    def evaluate_distance(self, dirt_list, path_photo):

        for counter, file_name in enumerate(dir_list):

            file_path_name = path_photo + "/" + file_name

            # Call select_pixel method to select a pixel on the image
            selected_sx, selected_sy, output_x, output_y = self.select_pixel(
                file_path_name
            )

            if counter == 0:
                self.set_first_img()

                scale_x1 = selected_sx[0]
                scale_x2 = selected_sx[1]
                scale_y1 = selected_sy[0]
                scale_y2 = selected_sy[1]

                scale_pixel_dist = np.sqrt(
                    (scale_x1 - scale_x2) ** 2 + (scale_y1 - scale_y2) ** 2
                )

                print(f"Scale value was set as {scale_value} mm")
                print(
                    f"Pixel scale was computed as {scale_pixel_dist} divided by scale value"
                )

                self.app.pixel_scale.set(
                    self.app.pixel_scale.get() + " = " + str(scale_pixel_dist) + " px"
                )
                self.app.check_var_1.set(1)

            print("gkjsfkljs", output_x)
            if output_x is None:
                distance_ = 0
            else:

                x1 = output_x[0]
                x2 = output_x[1]
                y1 = output_y[0]
                y2 = output_y[1]

                pixel_distance = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

                distance_ = pixel_distance * (scale_value / scale_pixel_dist)

            print(
                f"Distance between selected points: {distance_} mm. Evaluated {file_name}. {counter+1}/{len(dir_list)}"
            )

            name_file = "/output.csv"

            path_name = self.output_path + name_file
            header = ["name", "distance"]

            # open the file in the write mode
            with open(path_name, "a") as f:
                # create the csv writer
                writer = csv.writer(f)

                if counter == 0:
                    # write the header
                    writer.writerow(header)

                # write a row to the csv file
                writer.writerow([distance_, file_name])

            self.reset_values()


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
        self.root.title("Pixel Selector Instructions")

        self.pixel_scale.set("Scale value was set: 10 mm")

        tk.Checkbutton(
            self.root,
            text="1# Select pixels for scale",
            variable=self.check_var_1,
            state="disabled",
        ).pack()

        var = tk.StringVar()
        label = tk.Label(self.root, textvariable=self.pixel_scale)
        label.pack()

        tk.Checkbutton(
            self.root,
            text="2# Selct pixels on adherends where initial delamination was marked",
            variable=self.check_var_2,
            state="disabled",
        ).pack()
        tk.Checkbutton(
            self.root,
            text="3# Selct pixel in the end of the crack",
            variable=self.check_var_3,
            state="disabled",
        ).pack()

        self.quit_button = tk.Button(
            self.root, text="Pass (p)", fg="black", command=self.on_pass_button_click
        )
        self.quit_button.pack(side="left")
        self.quit_button = tk.Button(
            self.root, text="Next (n)", fg="blue", command=self.on_next_button_click
        )
        self.quit_button.pack(side="left")
        self.quit_button = tk.Button(
            self.root, text="Back (b)", fg="blue", command=self.on_back_button_click
        )
        self.quit_button.pack(side="left")
        self.quit_button = tk.Button(
            self.root, text="Quit (q)", fg="red", command=self.on_quit_button_click
        )
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
if __name__ == "__main__":

    ####
    # Input parameters
    ####
    path_photo = "/home/horakl/School/FAV/Doktorske_studium/Prace_projekty/Mechanical_testing/Python_files/Distance_in_image/Input_folder"

    output_path = "/home/horakl/School/FAV/Doktorske_studium/Prace_projekty/Mechanical_testing/Python_files/Distance_in_image/Output_folder"
    # scale value
    scale_value = 10  # in mm

    # Get the list of all files and directories
    dir_list = os.listdir(path_photo)
    dir_list = sorted(dir_list)
    print(dir_list)

    root = tk.Tk()
    instruction_widget = Application(root=root)
    thread_pixel = PixelSelector(dir_list, path_photo, instruction_widget, output_path)

    # First is started the pixel selector in thread. Then instruction widget is started in main loop.
    thread_pixel.start()
    instruction_widget.mainloop()  # Start the GUI thread

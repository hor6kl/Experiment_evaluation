import threading
import tkinter as tk
import time
import queue


# Shared queue
shared_q_pass = queue.Queue()
shared_q_quit = queue.Queue()
shared_q_next = queue.Queue()
shared_q_back = queue.Queue()

class App(threading.Thread):
    def __init__(self, arg=()):
        super().__init__()
        self.arg = arg

    def run(self):
        while shared_q_quit.empty():
            # Check the queue for messages
            if not shared_q_pass.empty():
                message = shared_q_pass.get()
                print("Received message:", message)
                app.check_var_1.set(1)
            if not shared_q_next.empty():
                massege = shared_q_next.get()
                print("Received message:", message)
            if not shared_q_back.empty():
                massege = shared_q_back.get()
                print("Received message:", message)

class Application(tk.Frame):
    def __init__(self, root=()):
        super().__init__(root)
        self.root = root 
        self.pack()
        self.check_var_1 = tk.IntVar()
        self.check_var_2 = tk.IntVar()
        self.check_var_3 = tk.IntVar()
        self.create_widgets()

    def create_widgets(self):
        self.root.title('Pixel Selector Instructions')


        tk.Checkbutton(self.root, text="1# Select pixels for scale", variable = self.check_var_1, state='disabled').pack()

        var = tk.StringVar()
        label = tk.Label( self.root, textvariable=var)
        var.set(f'Scale value was set: 10 mm')
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
        shared_q_pass.put("Button Pass Pressed")
    def on_next_button_click(self):
        shared_q_pass.put("Button Next Pressed")
    def on_back_button_click(self):
        shared_q_pass.put("Button Back Pressed")

    def on_quit_button_click(self):
        shared_q_quit.put("Quit")
        self.root.quit()
        self.root.destroy()



if __name__ == "__main__":
    root = tk.Tk()
    app = Application(root=root)
    app2 = App(arg=app)

    app2.start()  # Start the App thread
    app.mainloop()  # Start the GUI thread





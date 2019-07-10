#!/usr/bin/env python3
import tkinter as tk
import pickle,os
from tkinter import *
import tkinter as tk
import numpy as np
from PIL import ImageTk,Image

root=tk.Tk()

grid= tk.Frame(root)
grid.pack()

#img0=PhotoImage(file="0.png")
#img1=PhotoImage(file="1.png")
images={0:ImageTk.PhotoImage(Image.open("0.png")),1:ImageTk.PhotoImage(Image.open("1.png"))}

fill = 1

class button:
    def __init__(self, x, y):
        self.type=0
        self.but=Button(grid,command=self.change, image=images[0], borderwidth=0)
        self.but.grid(row=y, column=x)
        #Changed
        self.already_changed = False

    def change(self):
        if self.type==fill:
            self.but.config(image=img0)
            self.type=0
        else:
            self.but.config(image=images[fill]) #I left this in here, but you should NEVER use eval(). It's unsafe.
            self.type=fill

    #Changed
    def mouse_entered(self):
        if not self.already_changed:
            self.change()
            self.already_changed = True

    def mouse_up(self):
        self.already_changed = False

#Changed
class Container:
    def __init__(self, x, y):
        grid_buttons = []
        
        for Y in range(y):
            grid_buttons.append([])
            for X in range(x):
                grid_buttons[Y].append(button(X, Y))
        
        self.buttons = grid_buttons
        grid.bind_all("<Button-1>", self.mouse_down)
        grid.bind_all("<ButtonRelease-1>", self.mouse_up)
        grid.bind_all("<B1-Motion>", self.mouse_motion)
        self.mouse_pressed = False

    def mouse_down(self, e):
        self.mouse_pressed = True

    def mouse_up(self, e):
        self.mouse_pressed = False
        for row in self.buttons:
            for but in row:
                but.mouse_up()
    
    def mouse_motion(self, e):
        for row in self.buttons:
            for but in row:
                if grid.winfo_containing(e.x_root, e.y_root) is but.but:
                    but.mouse_entered()
    
    def return_pressed_button_locations(self):
        pressedMatrix = np.zeros([len(self.buttons[0]),len(self.buttons)])
        for buttonrow,row in zip(self.buttons,range(len(self.buttons))):
            for button,col in zip(buttonrow,range(len(buttonrow))):
                if button.type == fill:
                    pressedMatrix[col,row] = 1
        print(pressedMatrix)

container = Container(15,15)
OKbutton = Button(root,text="OK", fg="green",command=lambda:container.return_pressed_button_locations())
quitbutton = Button(root,text="Quit", fg="red",command=quit)
OKbutton.pack()
quitbutton.pack()
root.mainloop()

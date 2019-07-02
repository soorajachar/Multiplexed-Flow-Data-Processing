#!/usr/bin/env python3
import numpy as np
import pickle,os,string
import tkinter as tk
from PIL import ImageTk,Image
from tkinter import *

root = Tk()
frame = Frame(root)
frame.pack()

#callbacks
def enableEntry():
    entry.configure(state="normal")
    entry.update()

def disableEntry():
    entry.configure(state="disabled")
    entry.update()

#GUI widgets
entry = Entry(frame, width=80)
entry.pack(side='right')

var = StringVar()
disableEntryRadioButton = Radiobutton(frame, text="Disable", variable=var, value="0", command=disableEntry)
disableEntryRadioButton.pack(anchor=W)
enableEntryRadioButton = Radiobutton(frame, text="Enable", variable=var, value="1", command=enableEntry)
enableEntryRadioButton.pack(anchor=W)

root.mainloop()

#!/usr/bin/env python3
from tkinter import *
import tkinter.ttk as ttk
import tkinter as tk
import pickle 

hexcolors = pickle.load(open('usedColors.pkl','rb'))

root = tk.Tk()

v = tk.IntVar()
v.set(1)  # initializing the choice, i.e. Python

languages = [
    ("Python",1),
    ("Perl",2),
    ("Java",3),
    ("C++",4),
    ("C",5)
]

def ShowChoice():
    print(v.get())

tk.Label(root,
         text="""Choose your favourite
programming language:""",
         justify = tk.LEFT,
         padx = 20).pack()
i=0
for val, language in enumerate(languages):
    tk.Radiobutton(root, text=language,indicatoron = 0,width = 20,padx = 20, variable=v, command=ShowChoice,bg=hexcolors[i],value=val).pack(anchor=tk.W)
    i+=1

root.mainloop()

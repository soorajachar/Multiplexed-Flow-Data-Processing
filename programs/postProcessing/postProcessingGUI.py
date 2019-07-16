#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thursday July 5 13:52:27 2018

@author: acharsr
"""
import pickle,sys,os,json,math,subprocess,string
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
from pathlib import Path
import pickle
from itertools import product,combinations
import numpy as np
import tkinter as tk
import pandas as pd
import facetPlotLibrary as fpl
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
sys.path.insert(0, '../dataProcessing/')
from miscFunctions import sortSINumerically,reindexDataFrame
from umap import UMAP
from sklearn.manifold import Isomap
from dimensionalityReductionGUI import DimensionReductionHomePage

class PostProcessingStartPage(tk.Frame):
    def __init__(self, master,folderName):
        tk.Frame.__init__(self, master)
        
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l1 = tk.Label(mainWindow, text="""Postprocessing Type:""").grid(row=0,column=0)
        postProcessTypeVar = tk.StringVar(value='dr')
        rb1a = tk.Radiobutton(mainWindow,text='Dimensionality Reduction',variable=postProcessTypeVar,value='dr')
        rb1a.grid(row=1,column=0)

        def collectInputs():
            postProcessType = postProcessTypeVar.get()
            if postProcessType == 'dr':
                master.switch_frame(DimensionReductionHomePage,folderName,PostProcessingStartPage)
        
        def gotoHomePage():
            os.chdir(master.homedirectory)
            master.switch_frame(master.homepage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).pack(in_=buttonWindow,side=tk.LEFT)
        tk.Button(buttonWindow, text="Back",command=lambda: gotoHomePage()).pack(in_=buttonWindow,side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(in_=buttonWindow,side=tk.LEFT)

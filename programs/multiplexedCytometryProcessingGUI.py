#!/usr/bin/env python3
import sys,os,argparse
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import tkinter as tk
from tkinter import *
import pandas as pd
import pickle,warnings
sys.path.insert(0, 'experimentSetup/')
from experimentSetupGUI import ExperimentSetupStartPage
sys.path.insert(0, 'dataProcessing/')
from dataProcessingGUI import DataProcessingStartPage
from miscFunctions import parseCommandLineNNString,exitEnabledSubprocessRun
sys.path.insert(0, 'figuresPipeline/')
from facetPlottingGUI import FacetPlottingStartPage
sys.path.insert(0, 'postProcessing/')
from postProcessingGUI import PostProcessingStartPage 

pathToExperimentSpreadsheet = '../experiments/'
secondPath = '../outputData'

concUnit = 1e9
unitPrefixDictionary = {1e12:'pM',1e9:'nM',1e6:'uM',1e3:'mM',1e0:'M'}
concUnitPrefix = unitPrefixDictionary[concUnit]
        
def switchDirectories(researcherName):
    cwd = os.getcwd()
    #Local; only used by me
    if 'Volumes' not in cwd:
        os.chdir('../experiments/')
    #NIH Server; used by other people
    else:
        os.chdir('../../'+researcherName+'/experiments')

class MainApp(tk.Tk):
    def __init__(self):
        self.root = tk.Tk.__init__(self)
        self._frame = None
        self.homepage = StartPage
        self.homedirectory = os.getcwd()
        self.switch_frame(StartPage)

    def switch_frame(self, frame_class,*args):
        """Destroys current frame and replaces it with a new one."""
        new_frame = frame_class(self,*args)
        if self._frame is not None:
            self._frame.destroy()
        self._frame = new_frame
        self._frame.pack()

class StartPage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)

        os.chdir(master.homedirectory) 
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        l1 = tk.Label(mainWindow, text="""Action: """)
        v = tk.StringVar(value='es')
        rb1a = tk.Radiobutton(mainWindow, text="Setup experiment",padx = 20, variable=v, value='es')
        rb1b = tk.Radiobutton(mainWindow,text="Process data",padx = 20, variable=v, value='pd')
        rb1c = tk.Radiobutton(mainWindow,text="Plot data",padx = 20, variable=v, value='plt')
        rb1d = tk.Radiobutton(mainWindow,text="Postprocess data",padx = 20, variable=v, value='pp')

        l3 = tk.Label(mainWindow, text="""Experiment number: """)
        t = tk.Entry(mainWindow) 
        
        l1.grid(row=0,column=0)
        rb1a.grid(row=1,column=0,sticky=tk.W)
        rb1b.grid(row=2,column=0,sticky=tk.W)
        rb1c.grid(row=3,column=0,sticky=tk.W)
        rb1d.grid(row=4,column=0,sticky=tk.W)

        l3.grid(row=0,column=1)
        t.grid(row=1,column=1)
        def collectInputs():
            action = v.get()
            inputString = str(t.get())
            ex_data = pd.read_excel(pathToExperimentSpreadsheet+'masterExperimentSpreadsheet.xlsx')
            experimentsToRun = parseCommandLineNNString(inputString)
            cwd = os.getcwd()
            for expNum in experimentsToRun:
                expIndex = expNum - 1
                folderName = ex_data['Full Name'][expIndex]
                researcherName = ex_data['Researcher'][expIndex]
                switchDirectories(researcherName)
                if action == 'es':
                    master.switch_frame(ExperimentSetupStartPage,folderName,cwd)
                elif action == 'pd':
                    os.chdir(folderName)
                    master.switch_frame(DataProcessingStartPage,folderName,expNum,ex_data)
                elif action == 'plt':
                    os.chdir(folderName)
                    master.switch_frame(FacetPlottingStartPage,folderName)
                elif action == 'pp':
                    os.chdir(folderName)
                    master.switch_frame(PostProcessingStartPage,folderName)
            
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=1)

def main():
    app = MainApp()
    app.mainloop()

main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thursday July 5 13:52:27 2018

@author: acharsr
"""
from pathlib import Path
import pickle
from itertools import product,combinations
import numpy as np
import pickle,sys,os,json,math,subprocess,string
from tkinter import *
import tkinter as tk
from createExperimentLayout import produceGUIBasedIndexingCoordinates

experimentParameters = {}
parametersUpdatedByGridGUI = {}

class ExperimentSetupStartPage(tk.Frame):
    def __init__(self, master,fName,path):
        global folderName
        global pathToButtonImages
        folderName = fName
        pathToButtonImages = path+'/experimentSetup'
        tk.Frame.__init__(self, master)
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        #Make experiment folder and subfolders
        def createExperimentFolders(folderName):
            subprocess.run(['mkdir',folderName])
            experimentFolderNames = ['inputFiles','semiProcessedData','fullyProcessedFigures','postProcessedData','postProcessedData/dimensionalReductions','semiProcessedData/singleCellData']
            for experimentFolderName in experimentFolderNames:
                subprocess.run(['mkdir',folderName+'/'+experimentFolderName])
            print('Experiment Created')
            master.switch_frame(ExperimentSetupStartPage,fName,path)
        
        ceButton = tk.Button(mainWindow, text="Create experiment",command=lambda: createExperimentFolders(folderName))
        ceButton.pack(fill=X)
        
        epButton = tk.Button(mainWindow, text="Input experiment parameters",command=lambda: master.switch_frame(ExperimentParameterPage,folderName))
        epButton.pack(fill=X)
        
        def plateLayout():
            os.chdir(folderName)
            experimentParameters = json.load(open('inputFiles/experimentParameters-'+folderName+'.json','r'))
            parametersUpdatedByGridGUI = pickle.load(open('inputFiles/gui-parametersUpdatedByGridGUI.pkl','rb'))
            produceGUIBasedIndexingCoordinates(folderName,pathToButtonImages,experimentParameters,parametersUpdatedByGridGUI,master)
        
        elButton = tk.Button(mainWindow, text="Create experiment layout",command=lambda: plateLayout())
        elButton.pack(fill=X)
        
        for i,button in enumerate([ceButton,epButton,elButton]):
            if i == 1:
                requiredFiles = [folderName]
                for requiredFile in requiredFiles:
                    if requiredFile not in os.listdir():
                        button.config(state=tk.DISABLED)
            elif i == 2:
                if folderName in os.listdir():
                    os.chdir(folderName)
                    requiredFiles = ['experimentParameters-'+folderName+'.json','gui-parametersUpdatedByGridGUI.pkl']
                    for requiredFile in requiredFiles:
                        if requiredFile not in os.listdir('inputFiles')+os.listdir('semiProcessedData'):
                            button.config(state=tk.DISABLED)
                    os.chdir('..')
                else:
                    button.config(state=tk.DISABLED)

        def gotoHomePage():
            os.chdir(master.homedirectory)
            master.switch_frame(master.homepage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: gotoHomePage()).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class ExperimentParameterPage(tk.Frame):
    def __init__(self, master,folderName):
        tk.Frame.__init__(self, master)
        os.chdir(folderName)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        v = tk.IntVar()
        v2 = tk.IntVar()
        v3 = tk.IntVar()

        l1 = tk.Label(mainWindow, text="""Was this experiment performed with A/B plates?:""")
        rb1a = tk.Radiobutton(mainWindow, text="Yes",padx = 20, variable=v, value=1)
        rb1b = tk.Radiobutton(mainWindow,text="No",padx = 20, variable=v, value=0)
        l1.grid(row=0,column=0)
        rb1a.grid(row=0,column=1)
        rb1b.grid(row=0,column=2)
        
        l2 = tk.Label(mainWindow, text="""Was this experiment performed replicatewise?:""")
        rb2a = tk.Radiobutton(mainWindow, text="Yes",padx = 20, variable=v2, value=1)
        rb2b = tk.Radiobutton(mainWindow,text="No",padx = 20, variable=v2, value=0)
        l2.grid(row=1,column=0)
        rb2a.grid(row=1,column=1)
        rb2b.grid(row=1,column=2)
        
        l3 = tk.Label(mainWindow, text="""Was this experiment performed with 96 or 384 well plates?:""")
        rb3a = tk.Radiobutton(mainWindow, text="96",padx = 20, variable=v3, value=96)
        rb3b = tk.Radiobutton(mainWindow,text="384",padx = 20, variable=v3, value=384)
        l3.grid(row=2,column=0)
        rb3a.grid(row=2,column=1)
        rb3b.grid(row=2,column=2)
        
        l4 = tk.Label(mainWindow, text="Enter the total number of plates used in this experiment")
        e1 = tk.Entry(mainWindow)
        l4.grid(row=3,column=0)
        e1.grid(row=3,column=1)

        l5 = tk.Label(mainWindow, text="Enter the number of condition levels (including Time): ")
        e2 = tk.Entry(mainWindow)
        l5.grid(row=4,column=0)
        e2.grid(row=4,column=1)
         
        def collectInputs():
            experimentParameters['numPlates'] = int(e1.get())
            experimentParameters['numAllLevels'] = int(e2.get())
            if v.get() == 1:
                experimentParameters['paired'] = True
                if v2.get() == 1:
                    experimentParameters['replicateWise'] = True
                else:
                    experimentParameters['replicateWise'] = False
            else:
                experimentParameters['paired'] = False
                experimentParameters['replicateWise'] = False
            if v3.get() == 384:
                experimentParameters['overallPlateDimensions'] = [16,24]
                parametersUpdatedByGridGUI['currentPlateDimensions'] = [16,24]
            else:
                experimentParameters['overallPlateDimensions'] = [8,12]
                parametersUpdatedByGridGUI['currentPlateDimensions'] = [8,12]
            master.switch_frame(allLevelNamePage,folderName)
        
        buttonWindow = Frame(self)
        buttonWindow.pack(side=TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=1)
        
class allLevelNamePage(tk.Frame):
    def __init__(self, master,folderName):
        tk.Frame.__init__(self, master)
        numAllLevels = experimentParameters['numAllLevels']
        entryList1 = []
        entryList2 = []
        v = tk.IntVar()
        columnVariableRadioButtons = []
        
        mainWindow = Frame(self)
        mainWindow.pack(side=TOP,padx=10,pady=10)
        
        lt1 = tk.Label(mainWindow, text="Level Name").grid(row=0,column=1)
        lt2 = tk.Label(mainWindow, text="Number of Level Values").grid(row=0,column=2)
        lt3 = tk.Label(mainWindow, text="Column Variable?").grid(row=0,column=3)
        for conditionLevelNumber in range(1,numAllLevels+1):
            l1 = tk.Label(mainWindow, text="Condition "+str(conditionLevelNumber))
            e1 = tk.Entry(mainWindow)
            e2 = tk.Entry(mainWindow)
            rb1 = tk.Radiobutton(mainWindow, text="",variable=v, value=conditionLevelNumber-1)
            l1.grid(row=conditionLevelNumber,column=0)
            e1.grid(row=conditionLevelNumber,column=1)
            e2.grid(row=conditionLevelNumber,column=2)
            rb1.grid(row=conditionLevelNumber,column=3)
            entryList1.append(e1)
            entryList2.append(e2)
            columnVariableRadioButtons.append(rb1)
        
        def collectInputs():
            conditionNames = []
            numConditionLevelValues = []
            tiledLevels = []
            for allLevelNumber in range(numAllLevels):
                #Remove column variable from condition name list
                if v.get() == allLevelNumber:
                    experimentParameters['columnVariableName'] = str(entryList1[allLevelNumber].get())
                    experimentParameters['numColumnLevelValues'] = int(entryList2[allLevelNumber].get())
                else:
                    conditionNames.append(str(entryList1[allLevelNumber].get()))
                    numConditionLevelValues.append(int(entryList2[allLevelNumber].get()))

            experimentParameters['numConditionLevels'] = numAllLevels - 1
            experimentParameters['conditionLevelNames'] = conditionNames
            experimentParameters['allLevelNames'] = [experimentParameters['columnVariableName']]+conditionNames
            experimentParameters['numConditionLevelValues'] = numConditionLevelValues
            parametersUpdatedByGridGUI['numLevelsUnparsed'] = numAllLevels
            experimentParameters[''] = tiledLevels
            with open('inputFiles/gui-parametersUpdatedByGridGUI.pkl','wb') as f:
                pickle.dump(parametersUpdatedByGridGUI,f)
            master.switch_frame(columnLevelValuesPage,folderName)
        
        buttonWindow = Frame(self)
        buttonWindow.pack(side=TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=numAllLevels+1,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(ExperimentParameterPage)).grid(row=numAllLevels+1,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=numAllLevels+1,column=2)

class columnLevelValuesPage(tk.Frame):
    def __init__(self, master,folderName):
        numColumnLevelValues = experimentParameters['numColumnLevelValues']
        
        tk.Frame.__init__(self, master)
        
        mainWindow = Frame(self)
        mainWindow.pack(side=TOP,padx=10,pady=10)
        
        lt = tk.Label(mainWindow,text=experimentParameters['columnVariableName']+':').grid(row=0,column=0)
        col_wrap = 12
        for col in range(1,numColumnLevelValues+1):
            lt1 = tk.Label(mainWindow, text='Level Value '+str(col),width=10).grid(row=int((col-1)/col_wrap)*2,column=((col-1)%col_wrap)+1)
        entryList = []
        for columnLevelValueNumber in range(numColumnLevelValues):
            e1 = tk.Entry(mainWindow,width=10)
            e1.grid(row=2*int(columnLevelValueNumber/col_wrap)+1,column=(columnLevelValueNumber%col_wrap)+1)
            entryList.append(e1)

        def collectInputs():
            columnLevelValues = []
            for entry in entryList:
                columnLevelValues.append(float(entry.get()))
            experimentParameters['columnLevelValues'] = columnLevelValues
            master.switch_frame(conditionLevelValuesPage,folderName)
        
        buttonWindow = Frame(self)
        buttonWindow.pack(side=TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=2*int(numColumnLevelValues/col_wrap)+2,column=5)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(allLevelNamePage)).grid(row=2*int(numColumnLevelValues/col_wrap)+2,column=6)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=2*int(numColumnLevelValues/col_wrap)+2,column=7)

class conditionLevelValuesPage(tk.Frame):
    def __init__(self, master,folderName):
        numConditionLevels = experimentParameters['numConditionLevels']
        maxLevelValues = max(experimentParameters['numConditionLevelValues'])
        
        tk.Frame.__init__(self, master)
        
        mainWindow = Frame(self)
        mainWindow.pack(side=TOP,padx=10,pady=10)
        
        for col in range(1,maxLevelValues+1):
            lt1 = tk.Label(mainWindow, text="Level Value "+str(col)).grid(row=0,column=col)
        fullEntryList = []
        for conditionLevelNumber in range(numConditionLevels):
            l1 = tk.Label(mainWindow, text="Level values for \""+experimentParameters['conditionLevelNames'][conditionLevelNumber]+"\":").grid(row=conditionLevelNumber+1,column=0)
            levelEntryList = []
            for col in range(1,maxLevelValues+1):
                if col < experimentParameters['numConditionLevelValues'][conditionLevelNumber]+1:
                    e1 = tk.Entry(mainWindow)
                    e1.grid(row=conditionLevelNumber+1,column=col)
                    levelEntryList.append(e1)
            fullEntryList.append(levelEntryList)

        def collectInputs():
            conditionLevels = {}
            for lvlentrylist,i in zip(fullEntryList,range(numConditionLevels)):
                tempLevels = []
                for entry in lvlentrylist:
                    tempLevels.append(str(entry.get()))
                conditionLevels[experimentParameters['conditionLevelNames'][i]] = tempLevels

            experimentParameters['conditionLevelValues'] = conditionLevels
            experimentParameters['allLevelValues'] = conditionLevels
            experimentParameters['allLevelValues'][experimentParameters['columnVariableName']] = experimentParameters['columnLevelValues']
            with open('inputFiles/experimentParameters-'+folderName+'.json', 'w') as fp:
                json.dump(experimentParameters, fp)
            master.switch_frame(ExperimentSetupStartPage,folderName,pathToButtonImages)

        buttonWindow = Frame(self)
        buttonWindow.pack(side=TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=numConditionLevels+1,column=int(maxLevelValues/2))
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(columnLevelValuesPage)).grid(row=numConditionLevels+1,column=int(maxLevelValues/2)+1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=numConditionLevels+1,column=int(maxLevelValues/2)+2)

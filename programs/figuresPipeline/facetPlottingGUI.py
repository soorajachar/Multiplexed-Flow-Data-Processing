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
sys.path.insert(0, '../dataProcessing/')
import initialDataProcessing 

#Get level names and values into an easily accessible dictionary
def createLabelDict(df,dataType):
    fulldf = df.stack()
    labelDict = {}
    for i in range(fulldf.index.nlevels):
        levelName = fulldf.index.levels[i].name
        if levelName != 'Event':
            labelDict[levelName] = list(pd.unique(fulldf.index.get_level_values(levelName)))
    return labelDict

class checkUncheckAllButton(tk.Button):
    def __init__(self,parent,checkButtonList,**kwargs):
        tk.Button.__init__(self,parent,**kwargs)
        self.checkButtonList = checkButtonList
        self.parent = parent

    def checkAll(self):
        for checkButton in self.checkButtonList:
            checkButton.select()
    
    def uncheckAll(self):
        for checkButton in self.checkButtonList:
            checkButton.deselect()

class FacetPlottingStartPage(tk.Frame):
    def __init__(self, master,fName):
        global folderName
             
        folderName = fName
        
        tk.Frame.__init__(self, master)
        
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l2 = tk.Label(mainWindow, text="""Datatype: """)
        v2 = tk.StringVar(value='cyt')
        rb2a = tk.Radiobutton(mainWindow, text="Cytokine",padx = 20, variable=v2, value='cyt')
        rb2b = tk.Radiobutton(mainWindow,text="Cell",padx = 20, variable=v2, value='cell')
        rb2c = tk.Radiobutton(mainWindow,text="Proliferation",padx = 20, variable=v2, value='prolif')
        rb2d = tk.Radiobutton(mainWindow,text="Single Cell",padx = 20, variable=v2, value='singlecell')
        
        v3 = tk.BooleanVar(value=False)
        cb = tk.Checkbutton(mainWindow,padx = 20, variable=v3)

        l2.grid(row=0,column=0)
        rb2a.grid(row=1,column=0,sticky=tk.W)
        rb2b.grid(row=2,column=0,sticky=tk.W)
        rb2c.grid(row=3,column=0,sticky=tk.W)
        rb2d.grid(row=4,column=0,sticky=tk.W)
        
        l3 = tk.Label(mainWindow,text='Modified?').grid(row=0,column=1,sticky=tk.W)
        cb.grid(row=1,column=1,sticky=tk.W)
        
        def collectInputs():
            global useModifiedDf
            modified = v3.get()
            useModifiedDf = modified
            if modified:
                modifiedString = '-modified'
            else:
                modifiedString = ''
            global dataType
            dataType = str(v2.get())
            global experimentDf
            experimentDf = pickle.load(open('semiProcessedData/'+initialDataProcessing.dataTypeDataFrameFileNames[dataType]+'-'+folderName+modifiedString+'.pkl','rb'))
            global trueLabelDict
            trueLabelDict = {}
            trueLabelDict = createLabelDict(experimentDf,dataType)
            master.switch_frame(PlotTypePage)
        
        def gotoHomePage():
            os.chdir(master.homedirectory)
            master.switch_frame(master.homepage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: gotoHomePage()).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class PlotTypePage(tk.Frame):
    def __init__(self, master):
        plottableFigureDict = {'1d':['histogram','kde'],'categorical':['bar','violin','box','point'],'2d':['line','scatter'],'3d':['heatmap']}
        
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10)
        
        l1 = tk.Label(mainWindow, text="What type of figure do you want to plot?",pady=10).grid(row=0,column = 0,columnspan = len(plottableFigureDict),sticky=tk.N)
         
        plotTypeRadioButtons = []
        plotSelectionString = tk.StringVar(value='1d/histogram')
        maxNumSubplots = 0
        for plotTypeIndex,plotTypeTitle in enumerate(plottableFigureDict):
            plotTypeTitleLabel = tk.Label(mainWindow,text=plotTypeTitle).grid(row=1,column=plotTypeIndex,sticky=tk.NW)
            temprblist = []
            tempselectionstring = []
            for subPlotTypeIndex,subPlotTitle in enumerate(plottableFigureDict[plotTypeTitle]):
                rb = tk.Radiobutton(mainWindow, text=subPlotTitle,padx = 20, variable=plotSelectionString, value=plotTypeTitle+'/'+subPlotTitle)
                rb.grid(row=subPlotTypeIndex+2,column=plotTypeIndex,sticky=tk.NW)
                temprblist.append(rb)
            plotTypeRadioButtons.append(temprblist)
            if len(plottableFigureDict[plotTypeTitle]) > maxNumSubplots:
                maxNumSubplots = len(plottableFigureDict[plotTypeTitle])
        
        stripSwarmBool = tk.BooleanVar()
        cb = tk.Checkbutton(mainWindow,text='Add strip/swarm points to categorical plot',variable=stripSwarmBool,pady=20)
        cb.grid(row=maxNumSubplots+2,column=0,columnspan=len(plottableFigureDict))
        
        def collectInputs():
            global plotType
            global subPlotType
            global addDistributionPoints
            addDistributionPoints = stripSwarmBool.get()
            plotType,subPlotType = plotSelectionString.get().split('/')
            master.switch_frame(selectLevelsPage)
        
        def quitCommand():
            exitBoolean = True
            quit()

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP)
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(FacetPlottingStartPage,folderName)).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quitCommand()).pack(side=tk.LEFT)

class selectLevelsPage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        labelWindow = tk.Frame(self)
        labelWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        global figureLevelList,fullFigureLevelBooleanList
        fullFigureLevelBooleanList = []
        figureLevelList = []
        
        l1 = tk.Label(labelWindow, text="""Which levels names do you want to be included within this figure??:""").pack()
        mainWindow = tk.Frame(self)
        levelNameCheckButtons = []
        checkButtonVariableList = []
        for levelName,i in zip(trueLabelDict.keys(),range(len(trueLabelDict.keys()))):
            includeLevelBool = tk.BooleanVar(value=True)
            cb = tk.Checkbutton(mainWindow, text=levelName,padx = 20, variable=includeLevelBool,onvalue=True)
            cb.grid(row=i+3,column=1,sticky=tk.W)
            cb.select()
            levelNameCheckButtons.append(cb)
            checkButtonVariableList.append(includeLevelBool)
        
        checkButtonWindow = tk.Frame(self)
        checkAllButton1 = checkUncheckAllButton(checkButtonWindow,levelNameCheckButtons, text='Check All')
        checkAllButton1.configure(command=checkAllButton1.checkAll)
        checkAllButton1.pack(side=tk.LEFT)
        
        uncheckAllButton1 = checkUncheckAllButton(checkButtonWindow,levelNameCheckButtons, text='Uncheck All')
        uncheckAllButton1.configure(command=checkAllButton1.uncheckAll)
        uncheckAllButton1.pack(side=tk.LEFT)
        
        checkButtonWindow.pack(side=tk.TOP)
        mainWindow.pack(side=tk.TOP,padx=10)

        def collectInputs():
            includeLevelList = []
            for checkButtonVariable in checkButtonVariableList:
                includeLevelList.append(checkButtonVariable.get())
            for figureLevelBool,levelName in zip(includeLevelList,trueLabelDict):
                if figureLevelBool:
                    figureLevelList.append(levelName)
                    fullFigureLevelBooleanList.append(True)
                else:
                    fullFigureLevelBooleanList.append(False)
            master.switch_frame(selectLevelValuesPage)
        
        def quitCommand():
            exitBoolean = True
            quit()

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).pack(in_=buttonWindow,side=tk.LEFT)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(PlotTypePage)).pack(in_=buttonWindow,side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quitCommand()).pack(in_=buttonWindow,side=tk.LEFT)

class selectLevelValuesPage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        
        global includeLevelValueList
        includeLevelValueList = []
        
        labelWindow = tk.Frame(self)
        labelWindow.pack(side=tk.TOP,padx=10,fill=tk.X,expand=True)
        
        l1 = tk.Label(labelWindow, text='Which specific level values do you want to include in the figure?',pady=10).grid(row=0,column = 0,columnspan=len(trueLabelDict)*6)
        levelValueCheckButtonList = []
        overallCheckButtonVariableList = []
        checkAllButtonList = []
        uncheckAllButtonList = []
        i=0
        maxNumLevelValues = 0
        for levelName in trueLabelDict:
            j=0
            levelCheckButtonList = []
            levelCheckButtonVariableList = []
            levelLabel = tk.Label(labelWindow, text=levelName+':')
            levelLabel.grid(row=1,column = i*6,sticky=tk.N,columnspan=5)
            for levelValue in trueLabelDict[levelName]:
                includeLevelValueBool = tk.BooleanVar()
                cb = tk.Checkbutton(labelWindow, text=levelValue, variable=includeLevelValueBool)
                cb.grid(row=j+4,column=i*6+2,columnspan=2,sticky=tk.W)
                labelWindow.grid_columnconfigure(i*6+3,weight=1)
                cb.select()
                levelCheckButtonList.append(cb)
                levelCheckButtonVariableList.append(includeLevelValueBool)
                j+=1
            
            checkAllButton1 = checkUncheckAllButton(labelWindow,levelCheckButtonList, text='Check All')
            checkAllButton1.configure(command=checkAllButton1.checkAll)
            checkAllButton1.grid(row=2,column=i*6,sticky=tk.N,columnspan=3)
            checkAllButtonList.append(checkAllButton1)
            
            uncheckAllButton1 = checkUncheckAllButton(labelWindow,levelCheckButtonList, text='Uncheck All')
            uncheckAllButton1.configure(command=checkAllButton1.uncheckAll)
            uncheckAllButton1.grid(row=2,column=i*6+3,sticky=tk.N,columnspan=3)
            uncheckAllButtonList.append(checkAllButton1)

            levelValueCheckButtonList.append(levelCheckButtonList)
            overallCheckButtonVariableList.append(levelCheckButtonVariableList)
            if len(trueLabelDict[levelName]) > maxNumLevelValues:
                maxNumLevelValues = len(trueLabelDict[levelName])
            i+=1

        def collectInputs():
            for checkButtonVariableList in overallCheckButtonVariableList:
                tempLevelValueList = []
                for checkButtonVariable in checkButtonVariableList:
                    tempLevelValueList.append(checkButtonVariable.get())
                includeLevelValueList.append(tempLevelValueList)
            master.switch_frame(assignLevelsToParametersPage)
        
        def quitCommand():
            exitBoolean = True
            quit()
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=maxNumLevelValues+4,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(selectLevelsPage)).grid(row=maxNumLevelValues+4,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quitCommand()).grid(row=maxNumLevelValues+4,column=2)

class assignLevelsToParametersPage(tk.Frame):
    
    def __init__(self, master):
        parameterTypeDict = {
                'categorical':['Color','Order', 'Row', 'Column','None'],
                '1d':['Color','Row','Column','None'],
                '2d':['Marker','Color','Size','Row','Column','X Axis Values','None'],
                '3d':['Row','Column','X Axis Values','Y Axis Values']}
        
        tk.Frame.__init__(self, master)
            
        global parametersSelected
        parametersSelected = {}
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l1 = tk.Label(mainWindow, text='Which plotting parameter do you want to assign to each of your figure levels?',pady=10).grid(row=0,column = 0,columnspan = len(figureLevelList))
        rblist = []
        parameterVarList = []
        for figureLevel,figureLevelIndex in zip(figureLevelList,range(len(figureLevelList))):
            v = tk.IntVar()
            temprblist = []
            levelLabel = tk.Label(mainWindow, text=figureLevel+':')
            levelLabel.grid(row=1,column=figureLevelIndex,sticky=tk.NW)
            for plottingParameter,parameterIndex in zip(parameterTypeDict[plotType],range(len(parameterTypeDict[plotType]))):
                rb = tk.Radiobutton(mainWindow, text=plottingParameter,padx = 20, variable=v, value=parameterIndex)
                rb.grid(row=parameterIndex+2,column=figureLevelIndex,sticky=tk.NW)
                temprblist.append(rb)
            rblist.append(temprblist)
            parameterVarList.append(v)
        
        def collectInputs():
            for parameterVar,levelName in zip(parameterVarList,figureLevelList):
                if parameterTypeDict[plotType][parameterVar.get()] not in parametersSelected.keys():
                    parametersSelected[parameterTypeDict[plotType][parameterVar.get()]] = levelName
                else:
                    if not isinstance(parametersSelected[parameterTypeDict[plotType][parameterVar.get()]],list):
                        parametersSelected[parameterTypeDict[plotType][parameterVar.get()]] = [parametersSelected[parameterTypeDict[plotType][parameterVar.get()]]]+[levelName]
                    else:
                        parametersSelected[parameterTypeDict[plotType][parameterVar.get()]] = parametersSelected[parameterTypeDict[plotType][parameterVar.get()]].append(levelName)
            master.switch_frame(plotElementsGUIPage)
        
        def quitCommand():
            exitBoolean = True
            quit()
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=len(parameterTypeDict[plotType])+2,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(selectLevelValuesPage)).grid(row=len(parameterTypeDict[plotType])+2,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quitCommand()).grid(row=len(parameterTypeDict[plotType])+2,column=2)

def getDefaultAxisTitles():
        experimentParameters = json.load(open('inputFiles/experimentParameters-'+folderName+'.json','r'))
        xaxistitle = ''
        yaxistitle = ''
        cbartitle = ''
        cytokineDefault = 'Concentration (nM)'
        if plotType == '1d':
            xaxistitle = 'GFI'
            yaxistitle = 'Frequency'
        else:
            if plotType == 'categorical':
                xaxistitle = parametersSelected['Order']
                if dataType == 'cyt':
                    yaxistitle = cytokineDefault
            else:
                xaxistitle = parametersSelected['X Axis Values']
                if dataType == 'cyt':
                    if plotType == '2d':
                        yaxistitle = cytokineDefault
                    else:
                        cbartitle = cytokineDefault
        if xaxistitle == 'Time':
            xaxistitle += ' (hours)'
        return [xaxistitle,yaxistitle,cbartitle]

class plotElementsGUIPage(tk.Frame):
    def __init__(self, master):
        axisDict = {'categorical':['X','Y'],'1d':['Y'],'2d':['X','Y'],'3d':['X','Y','Colorbar']}
        scalingList = ['Linear','Logarithmic','Biexponential']
        axisSharingList = ['row','col','']
        axisTitleDefaults = getDefaultAxisTitles() 

        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)

        tk.Label(mainWindow, text='Title: ').grid(row=1,column=0,sticky=tk.W)
        for scaling,scalingIndex in zip(scalingList,range(len(scalingList))):
            tk.Label(mainWindow, text=scaling+' Scaling: ').grid(row=scalingIndex+2,column=0,sticky=tk.W)
        tk.Label(mainWindow, text='Linear Range (Biexponential Scaling): ').grid(row=len(scalingList)+2,column=0,sticky=tk.W)
        tk.Label(mainWindow, text='Convert to numeric: ').grid(row=len(scalingList)+3,column=0,sticky=tk.W)
        tk.Label(mainWindow, text='Share axis across: ').grid(row=len(scalingList)+4,column=0,sticky=tk.W)
        tk.Label(mainWindow, text='Axis limits: ').grid(row=len(scalingList)+5,column=0,sticky=tk.W)

        entryList = []
        scalingVariableList = []
        radioButtonList = []
        checkButtonList = []
        checkButtonVarList = []
        radioButtonList2 = []
        radioButtonVarList2 = []
        linearRangeScalingList = []
        limitEntryList = []
        for axis,axisIndex in zip(axisDict[plotType],range(len(axisDict[plotType]))):
            tk.Label(mainWindow, text=axis+ ' Axis').grid(row=0,column=axisIndex+1)
            
            e1 = tk.Entry(mainWindow)
            e1.grid(row=1,column=axisIndex+1)
            e1.insert(0, axisTitleDefaults[axisIndex])
            entryList.append(e1)
            
            axisRadioButtonList = []
            v = tk.StringVar(value='Linear')
            for scaling,scalingIndex in zip(scalingList,range(len(scalingList))):
                rb = tk.Radiobutton(mainWindow,variable=v,value=scaling)
                rb.grid(row=scalingIndex+2,column=axisIndex+1)
                axisRadioButtonList.append(rb)
            radioButtonList.append(axisRadioButtonList)
            scalingVariableList.append(v)
            
            e2 = tk.Entry(mainWindow)
            e2.grid(row=len(scalingList)+2,column=axisIndex+1)
            linearRangeScalingList.append(e2)
            
            b = tk.BooleanVar(value=False)
            cb = tk.Checkbutton(mainWindow,variable=b)
            cb.grid(row=len(scalingList)+3,column=axisIndex+1)
            checkButtonList.append(cb)
            checkButtonVarList.append(b)
       
            shareWindow = tk.Frame(mainWindow)
            shareWindow.grid(row=len(scalingList)+4,column=axisIndex+1)
            shareString = tk.StringVar(value='None')
            rb2a = tk.Radiobutton(shareWindow,variable=shareString,text='None',value='None')
            rb2b = tk.Radiobutton(shareWindow,variable=shareString,text=axisSharingList[axisIndex],value=axisSharingList[axisIndex])
            rb2a.grid(row=0,column=1)
            rb2b.grid(row=0,column=3)
            radioButtonList2.append([rb2a,rb2b])
            radioButtonVarList2.append(shareString)

            limitWindow = tk.Frame(mainWindow)
            limitWindow.grid(row=len(scalingList)+5,column=axisIndex+1)
            #ll = tk.Label(limitWindow,text='Lower:').grid(row=0,column=0)
            e3 = tk.Entry(limitWindow,width=5)
            e3.grid(row=0,column=1)
            #ul = tk.Label(limitWindow,text='Upper:').grid(row=0,column=2)
            e4 = tk.Entry(limitWindow,width=5)
            e4.grid(row=0,column=3)
            limitEntryList.append([e3,e4])

        def collectInputs():
            plotOptions = {}
            for axis,axisIndex in zip(axisDict[plotType],range(len(axisDict[plotType]))):
                if radioButtonVarList2[axisIndex].get() == 'None':
                    share = False
                else:
                    share = radioButtonVarList2[axisIndex].get()
                plotOptions[axis] = {'axisTitle':entryList[axisIndex].get(),
                        'axisScaling':scalingVariableList[axisIndex].get(),
                        'linThreshold':linearRangeScalingList[axisIndex].get(),
                        'numeric':checkButtonVarList[axisIndex].get(),
                        'share':share,
                        'limit':[limitEntryList[axisIndex][0].get(),limitEntryList[axisIndex][1].get()]}
            subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually = fpl.produceSubsettedDataFrames(experimentDf.stack().to_frame('temp'),fullFigureLevelBooleanList,includeLevelValueList)
            fpl.plotFacetedFigures(folderName,plotType,subPlotType,dataType,subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually,useModifiedDf,experimentDf,plotOptions,parametersSelected,addDistributionPoints)
            trueLabelDict = {} 
            master.switch_frame(FacetPlottingStartPage,folderName)
            

        def quitCommand():
            exitBoolean = True
            quit()
         
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=len(scalingList)+4,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(assignLevelsToParametersPage)).grid(row=len(scalingList)+4,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quitCommand()).grid(row=len(scalingList)+4,column=2)

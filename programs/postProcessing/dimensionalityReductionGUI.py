#!/usr/bin/env python3 
import json,pickle,math,matplotlib,sys,os,string
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import numpy as np
import pandas as pd
import seaborn as sns
import tkinter as tk
import itertools
from matplotlib import pyplot as plt
from scipy import stats
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
sys.path.insert(0, '../dataProcessing/')
from miscFunctions import sortSINumerically,reindexDataFrame,setMaxWidth
sys.path.insert(0, '../figuresPipeline/')
from facetPlottingGUI import checkUncheckAllButton
import facetPlotLibrary as fpl
from umap import UMAP
from sklearn.manifold import Isomap

idx = pd.IndexSlice
dataTypeObservableRangeDict = {'cyt':1,'cell':3,'prolif':1}
realDataTypeNameDict = {'cyt':'Supernatant','cell':'Surface/Intracellular Marker','prolif':'Proliferation'}

#Get level names and values into an easily accessible dictionary
def createLabelDict(df,levelRange):
    fulldf = df.stack()
    labelDict = {}
    for i in range(levelRange[0],levelRange[1]):
        levelName = fulldf.index.levels[i].name
        if levelName != 'Event':
            labelDict[levelName] = list(pd.unique(fulldf.index.get_level_values(levelName)))
    return labelDict

class DimensionReductionHomePage(tk.Frame):
    def __init__(self, master,fName):
        global folderName
        folderName = fName
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        l1 = tk.Label(mainWindow, text="""Action: """).grid(row=0,column=0,sticky=tk.W)
        v = tk.StringVar(value='pd')
        rb1a = tk.Radiobutton(mainWindow,text="Reduce dimensions",padx = 20, variable=v, value='pd')
        rb1b = tk.Radiobutton(mainWindow,text="Plot dimensionally reduced data",padx = 20, variable=v, value='plt')
        
        rb1a.grid(row=1,column=0,sticky=tk.W)
        rb1b.grid(row=2,column=0,sticky=tk.W)
        
        def collectInputs():
            action = v.get()
            if action == 'pd':
                master.switch_frame(DimensionReductionProcessingPage)
            else:
                master.switch_frame(DimensionPlottingPage)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=1)

class DimensionReductionProcessingPage(tk.Frame):
    def __init__(self, master):
        dataTypeFileDict = {'cyt':'cytokineConcentrationPickleFile-'+folderName+'-modified.pkl','cell':'cellStatisticPickleFile-'+folderName+'-modified.pkl','prolif':'proliferationStatisticPickleFile-'+folderName+'-modified.pkl'}
        
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l2 = tk.Label(mainWindow, text="DataTypes to postprocess:").grid(row=0,column=1)
        dataTypeCheckButtons = []
        checkButtonVariableList = []
        for row,dataType in enumerate(dataTypeFileDict):
            dataTypeBool = tk.BooleanVar()
            cb = tk.Checkbutton(mainWindow, text=realDataTypeNameDict[dataType],padx = 20, variable=dataTypeBool)
            cb.grid(row=row+1,column=1,sticky=tk.W)
            if dataTypeFileDict[dataType] not in os.listdir('semiProcessedData'):
                cb.config(state=tk.DISABLED)
            else:
                cb.select()
            dataTypeCheckButtons.append(cb)
            checkButtonVariableList.append(dataTypeBool)
        
        def collectInputs():
            global dataTypeDfDict
            dataTypeDfDict = {}
            for checkButtonVariable,dataType in zip(checkButtonVariableList,dataTypeFileDict):
                if checkButtonVariable.get():
                    dataTypeDfDict[dataType] = pickle.load(open('semiProcessedData/'+dataTypeFileDict[dataType],'rb'))
            if 'cyt' in dataTypeDfDict.keys():
                dataTypeDfDict['cyt'] = np.log10(dataTypeDfDict['cyt'])
            master.switch_frame(SelectDimensionsPage,folderName,dataTypeDfDict)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(DimensionReductionHomePage)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=1)

class SelectDimensionsPage(tk.Frame):
    def __init__(self, master,fName,dataTypeDfDct):
        tk.Frame.__init__(self, master)
        
        labelWindow = tk.Frame(self)
        labelWindow.pack(side=tk.TOP,padx=10,fill=tk.X,expand=True)
        l1 = tk.Label(labelWindow, text='Which specific measurables do you want to include in the dimension reduction?',pady=10).pack()
        dataTypeLevelCheckButtonList = []
        dataTypeLevelCheckButtonVariableList = []
        for dataType in dataTypeDfDict:
            dataTypeDf = dataTypeDfDict[dataType]
            observableLevelDict = createLabelDict(dataTypeDf,[0,dataTypeObservableRangeDict[dataType]])
            levelValueCheckButtonList = []
            overallCheckButtonVariableList = []
            checkAllButtonList = []
            uncheckAllButtonList = []
            dataTypeWindow = tk.Frame(labelWindow,borderwidth=1,relief='groove')
            dataTypeLabel = tk.Label(dataTypeWindow,text=realDataTypeNameDict[dataType]+':',font="-weight bold").grid(row=0,column=0,columnspan=len(observableLevelDict)*6)
            i=0
            maxNumLevelValues = 0
            for levelName in observableLevelDict:
                j=0
                levelCheckButtonList = []
                levelCheckButtonVariableList = []
                levelLabel = tk.Label(dataTypeWindow, text=levelName+':')
                levelLabel.grid(row=1,column = i*6,sticky=tk.N,columnspan=5)
                for levelValue in observableLevelDict[levelName]:
                    includeLevelValueBool = tk.BooleanVar()
                    cb = tk.Checkbutton(dataTypeWindow, text=levelValue, variable=includeLevelValueBool)
                    cb.grid(row=j+4,column=i*6+2,columnspan=2,sticky=tk.W)
                    dataTypeWindow.grid_columnconfigure(i*6+3,weight=1)
                    cb.select()
                    levelCheckButtonList.append(cb)
                    levelCheckButtonVariableList.append(includeLevelValueBool)
                    j+=1
                
                checkAllButton1 = checkUncheckAllButton(dataTypeWindow,levelCheckButtonList, text='Check All')
                checkAllButton1.configure(command=checkAllButton1.checkAll)
                checkAllButton1.grid(row=2,column=i*6,sticky=tk.N,columnspan=3)
                checkAllButtonList.append(checkAllButton1)
                
                uncheckAllButton1 = checkUncheckAllButton(dataTypeWindow,levelCheckButtonList, text='Uncheck All')
                uncheckAllButton1.configure(command=checkAllButton1.uncheckAll)
                uncheckAllButton1.grid(row=2,column=i*6+3,sticky=tk.N,columnspan=3)
                uncheckAllButtonList.append(checkAllButton1)

                levelValueCheckButtonList.append(levelCheckButtonList)
                overallCheckButtonVariableList.append(levelCheckButtonVariableList)
                if len(observableLevelDict[levelName]) > maxNumLevelValues:
                    maxNumLevelValues = len(observableLevelDict[levelName])
                i+=1
            dataTypeLevelCheckButtonList.append(levelValueCheckButtonList)
            dataTypeLevelCheckButtonVariableList.append(overallCheckButtonVariableList)
            dataTypeWindow.pack(side=tk.LEFT,fill=tk.Y)
        
        def collectInputs():
            global dimensionDict
            dimensionDict = {}
            for obl,dataType in zip(dataTypeLevelCheckButtonVariableList,dataTypeDfDict):
                includeLevelValueList = []
                for checkButtonVariableList in obl:
                    tempLevelValueList = []
                    for checkButtonVariable in checkButtonVariableList:
                        tempLevelValueList.append(checkButtonVariable.get())
                    includeLevelValueList.append(tempLevelValueList)
                dimensionDict[dataType] = includeLevelValueList
            reduceDimensions(dimensionDict,dataTypeDfDict)
            master.switch_frame(DimensionReductionHomePage,folderName,dataTypeDfDict)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.BOTTOM,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=maxNumLevelValues+4,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(DimensionReductionHomePage)).grid(row=maxNumLevelValues+4,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=maxNumLevelValues+4,column=2)

def reduceDimensions(dimensionDict,dataTypeDfDict):
    dimensionReductionMatrixList = []
    didNotUseList = []
    #Iterate through datatypes
    for dataType in dataTypeDfDict:
        rowList = []
        #Iterate through each row in datatype df; grab dimension names
        for row in range(dataTypeDfDict[dataType].shape[0]):
            names = dataTypeDfDict[dataType].iloc[row,:].name
            dimensionNames = names[:dataTypeObservableRangeDict[dataType]]
            rowList.append(dimensionNames)
        #Go thorugh each level that was selected, add level values of each level
        selectedLevelList = []
        for i,level in enumerate(dataTypeDfDict[dataType].index.names[:len(dimensionDict[dataType])]):
            levelList = []
            for levelValue,includeLevelValue in zip(list(pd.unique(dataTypeDfDict[dataType].index.get_level_values(level))),dimensionDict[dataType][i]):
                if includeLevelValue:
                    levelList.append(levelValue)
                else:
                    didNotUseList.append(levelValue)
            selectedLevelList.append(levelList)
        
        #Get all possible combinations of level values from each level
        allPossibleSelectedLevelCombinations = itertools.product(*selectedLevelList)
        rowindexlist = []
        #From original dataframe; select all rows that appear in the all possible combination list 
        for levelCombination in allPossibleSelectedLevelCombinations:
            if levelCombination in rowList:
                indices = [i for i, x in enumerate(rowList) if x == levelCombination]
                rowindexlist+=indices
        subsettedDf = dataTypeDfDict[dataType].iloc[rowindexlist,:]
        dimensionReductionDf = subsettedDf.stack().unstack(dataTypeDfDict[dataType].index.names[:dataTypeObservableRangeDict[dataType]])
        commonIndex = dimensionReductionDf.index
        scaler = MinMaxScaler()
        standardizedMatrix = scaler.fit_transform(dimensionReductionDf)
        dimensionReductionMatrixList.append(standardizedMatrix)
    
    finalDimensionReductionMatrix = np.hstack(dimensionReductionMatrixList)
    matrixTSNE = TSNE().fit_transform(finalDimensionReductionMatrix)
    matrixUmap = UMAP().fit_transform(finalDimensionReductionMatrix)
    matrixIsomap = Isomap().fit_transform(finalDimensionReductionMatrix)
    indexingTuple = []
    tempDataType = list(dataTypeDfDict.keys())[0]
    reindexingDf = dataTypeDfDict[tempDataType]
    for i in range(dataTypeObservableRangeDict[tempDataType]):
        indexElement = dataTypeDfDict[tempDataType].index.get_level_values(dataTypeDfDict[tempDataType].index.names[i])[0]
        reindexingDf = reindexingDf.xs([indexElement],level=[dataTypeDfDict[tempDataType].index.names[i]])
    dfTSNE = reindexDataFrame(pd.DataFrame(matrixTSNE,index=commonIndex,columns=['Dimension 1','Dimension 2']),reindexingDf,True)
    dfUmap = reindexDataFrame(pd.DataFrame(matrixUmap,index=commonIndex,columns=['Dimension 1','Dimension 2']),reindexingDf,True)
    dfIsomap = reindexDataFrame(pd.DataFrame(matrixIsomap,index=commonIndex,columns=['Dimension 1','Dimension 2']),reindexingDf,True)
    allDimRedDf = pd.concat([dfTSNE,dfUmap,dfIsomap],names=['DimensionalReductionMethod'],keys=['TSNE','UMAP','ISOMAP'])
    if len(didNotUseList) == 0:
        dimensionReductionTitle = '-'.join(dataTypeDfDict.keys())+'-all'
    else:
        dimensionReductionTitle = '-'.join(dataTypeDfDict.keys())+'-no-'+'-'.join(didNotUseList)
    print(allDimRedDf)
    with open('postProcessedData/dimensionalReductions/dimensionalReduction-'+dimensionReductionTitle+'.pkl','wb') as f:
        pickle.dump(allDimRedDf,f)

class DimensionPlottingPage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l1 = tk.Label(mainWindow,text='Select a dimensional reduction to plot: ').pack()
        fileList = ['temp']
        for fileName in os.listdir('postProcessedData/dimensionalReductions/'):
            if 'pkl' in fileName:
                fileList.append(fileName.split('.pk')[0])
        fileList = fileList[1:]
        dropVar = tk.StringVar()
        dropMenu = tk.OptionMenu(mainWindow,dropVar,*fileList)
        dropMenu.pack()
        setMaxWidth(fileList,dropMenu)
        
        def collectInputs():
            global fileName
            fileName = dropVar.get()+'.pkl'
            global fileDf
            fileDf = pickle.load(open('postProcessedData/dimensionalReductions/'+fileName,'rb'))
            global trueLabelDict
            trueLabelDict = createLabelDict(fileDf,[0,len(fileDf.index.names)])
            master.switch_frame(selectLevelsPage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.BOTTOM,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=0,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(DimensionReductionHomePage,folderName)).grid(row=0,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=0,column=2)

class selectLevelsPage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        labelWindow = tk.Frame(self)
        labelWindow.pack(side=tk.TOP,padx=10,pady=10)
        
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
            global figureLevelList,fullFigureLevelBooleanList
            figureLevelList = []
            fullFigureLevelBooleanList = []
            for checkButtonVariable in checkButtonVariableList:
                includeLevelList.append(checkButtonVariable.get())
            for figureLevelBool,levelName in zip(includeLevelList,trueLabelDict):
                if figureLevelBool:
                    figureLevelList.append(levelName)
                    fullFigureLevelBooleanList.append(True)
                else:
                    fullFigureLevelBooleanList.append(False)
            print(figureLevelList)
            print(fullFigureLevelBooleanList)
            master.switch_frame(selectConditionsPage)
        
        def quitCommand():
            exitBoolean = True
            quit()

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).pack(in_=buttonWindow,side=tk.LEFT)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(DimensionPlottingPage)).pack(in_=buttonWindow,side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quitCommand()).pack(in_=buttonWindow,side=tk.LEFT)

class selectConditionsPage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        
        global includeConditionList
        includeConditionList = []
        
        labelWindow = tk.Frame(self)
        labelWindow.pack(side=tk.TOP,padx=10,fill=tk.X,expand=True)
        
        #l1 = tk.Label(labelWindow, text='Which specific conditions do you want to include in the figure?',pady=10).grid(row=0,column = 0)
        levelValueCheckButtonList = []
        overallCheckButtonVariableList = []
        checkAllButtonList = []
        uncheckAllButtonList = []

        maxLenIndex = 0
        global conditionLabelDict
        conditionLabelDict = createLabelDict(fileDf,[0,len(fileDf.index.names)])
        maxNumLevelValues = 0
        i=0
        for levelName in conditionLabelDict:
            j=0
            levelCheckButtonList = []
            levelCheckButtonVariableList = []
            levelLabel = tk.Label(labelWindow, text=levelName+':')
            levelLabel.grid(row=1,column = i*6,sticky=tk.N,columnspan=5)
            for levelValue in conditionLabelDict[levelName]:
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
            if len(conditionLabelDict[levelName]) > maxNumLevelValues:
                maxNumLevelValues = len(conditionLabelDict[levelName])
            i+=1

        def collectInputs():
            for checkButtonVariableList in overallCheckButtonVariableList:
                tempLevelValueList = []
                for checkButtonVariable in checkButtonVariableList:
                    tempLevelValueList.append(checkButtonVariable.get())
                includeConditionList.append(tempLevelValueList)
            master.switch_frame(assignLevelsToParametersPage)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=maxNumLevelValues+4,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(selectLevelsParametersPage)).grid(row=maxNumLevelValues+4,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=maxNumLevelValues+4,column=2)

class assignLevelsToParametersPage(tk.Frame):
    def __init__(self, master):
        parameterTypeDict = {
                'categorical':['Color','Order', 'Row', 'Column','None'],
                '1d':['Color','Row','Column'],
                '2d':['Marker','Color','Size','Row','Column','None'],
                '3d':['Row','Column','X Axis Values','Y Axis Values']}
        axisDict = {'categorical':['X','Y'],'1d':['Y'],'2d':['X','Y'],'3d':['X','Y','Colorbar']}

        plotType = '2d'
        subPlotType = 'scatter'
        dataType = 'dr'
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
                parametersSelected[parameterTypeDict[plotType][parameterVar.get()]] = levelName
            parametersSelected['X Axis Values'] = 'Dimension 1'
            parametersSelected['Y Axis Values'] = 'Dimension 2'
            subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually = fpl.produceSubsettedDataFrames(fileDf,fullFigureLevelBooleanList,includeConditionList)
            plotOptions = {}
            for axis,axisIndex in zip(axisDict[plotType],range(len(axisDict[plotType]))):
                plotOptions[axis] = {'axisTitle':fileDf.columns[axisIndex],'axisScaling':'Linear','linThreshold':0,'numeric':False,'share':False,'limit':['','']}
            fpl.plotFacetedFigures(folderName,plotType,subPlotType,dataType,subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually,False,fileDf,plotOptions,parametersSelected,False,alternateTitle=fileName.split('.pk')[0])
            master.switch_frame(DimensionReductionHomePage,folderName) 
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=len(parameterTypeDict[plotType])+2,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(selectConditionsPage)).grid(row=len(parameterTypeDict[plotType])+2,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=len(parameterTypeDict[plotType])+2,column=2)

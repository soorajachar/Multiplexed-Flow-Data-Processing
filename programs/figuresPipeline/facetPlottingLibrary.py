#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on fri sep 7 13:12:56 2018

@author: acharsr
"""

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import markers
import numpy as np
import seaborn as sns
import pandas as pd
import pickle,os,math,sys,itertools,re
from matplotlib.widgets import RadioButtons,Button,CheckButtons
sys.path.insert(0, '../dataProcessing/')
from miscFunctions import sortSINumerically
import subprocess
plt.switch_backend('QT4Agg') #default on my system

def facetPlottingGUI(df,plotType):
    fulldf = df.stack()
    
    buttonWidth = 0.1/2
    buttonLength = 0.075/2
    buttonXStart = 0.5-(0.01+buttonWidth)
    buttonYStart = 0.01

    #Get level names and values into an easily accessible dictionary
    labelDict = {}
    for i in range(fulldf.index.nlevels):
        levelName = fulldf.index.levels[i].name
        labelDict[levelName] = list(pd.unique(fulldf.index.get_level_values(levelName)))

    #Grab within figure boolean data
    fig = plt.figure()
    plt.axis('off')

    #Grab levels that will be used within a figure
    #plt.text(0.5, 1,'Which levels names do you want to\nbe included within this figure?',transform=plt.gca().transAxes,ha='center')
    rax = plt.axes([0.4, 0.5, 0.05*len(labelDict.keys()),0.05*len(labelDict.keys())])
    check = CheckButtons(rax, labelDict.keys(),actives=[True]*len(labelDict.keys()))
    plt.text(0.5, 1.2,'Which levels names do you want to be included within this figure?',ha='center')
    rax.spines['bottom'].set_visible(False)
    rax.spines['left'].set_visible(False)
    rax.spines['right'].set_visible(False)
    rax.spines['top'].set_visible(False)

    class Index(object):
        def OK(self, event):
            withinFigureBoolean = check.get_status()
            plt.close()
            print(withinFigureBoolean)
            with open('semiProcessedData/wfBool.pkl','wb') as f:
                pickle.dump(withinFigureBoolean,f)
        def Quit(self, event):
            sys.exit(0)    

    callback = Index()
    axOK = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
    axQuit = plt.axes([buttonXStart+buttonWidth+0.01,buttonYStart, buttonWidth, buttonLength])
    bOK = Button(axOK, 'OK')
    bOK.on_clicked(callback.OK)
    bQuit = Button(axQuit, 'Quit')
    bQuit.on_clicked(callback.Quit)
    plt.show()

    maxLabelLength = 0
    for labelName in labelDict:
        if len(labelDict[labelName]) > maxLabelLength:
            maxLabelLength = len(labelDict[labelName])
    fig = plt.figure()
    #fig = plt.figure(figsize=(maxLabelLength*0.6,2*len(labelDict.keys())))
    plt.axis('off')
    plt.text(0.5, 1.1,'Which specific level values do you want to include in the figure?',transform=plt.gca().transAxes,ha='center')
    i=0
    checkbuttons = []
    for levelName in labelDict:
        #rax2 = plt.axes([0.2+0.15*i, 0.1, 0.15,0.1*len(labelDict[levelName])])
        rectLength = 0.1*len(labelDict[levelName])
        rectWidth = (1- (0.02+0.01*len(labelDict.keys())))/len(labelDict.keys()) 
        if rectLength > 0.8:
            rectLength = 0.8
        
        rax2 = plt.axes([0.01*(i+1)+rectWidth*i, 0.94-(buttonWidth+rectLength), rectWidth,rectLength])
        rax2.spines['bottom'].set_visible(False)
        rax2.spines['left'].set_visible(False)
        rax2.spines['right'].set_visible(False)
        plt.text(0.5, 1.01,levelName,transform=plt.gca().transAxes,ha='center')
        checkbuttons.append(CheckButtons(rax2,labelDict[levelName],actives=[True]*len(labelDict[levelName])))
        i+=1

    class Index2(object):
        def OK(self, event):
            specificValueBooleanList = []
            for checkbutton in checkbuttons:
                specificValueBooleanList.append(checkbutton.get_status())
            plt.close()
            with open('semiProcessedData/svBoolList.pkl','wb') as f:
                pickle.dump(specificValueBooleanList,f)
        def Quit(self, event):
            sys.exit(0)    

    callback = Index2()
    axOK = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
    axQuit = plt.axes([buttonXStart+buttonWidth+0.01,buttonYStart, buttonWidth, buttonLength])
    bOK = Button(axOK, 'OK')
    bOK.on_clicked(callback.OK)
    bQuit = Button(axQuit, 'Quit')
    bQuit.on_clicked(callback.Quit)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.show()

    #Ask user which levels to assign to which of the four parameters that can be changed in a bar/point plot: color (c), order (o), x position of subplot (x), y position of subplot (y)
    #Order essentially is which level is plotted on the x axis, as bar plots are always categorical
    if plotType == 'categorical':
        parameterTypes = ['Color','Order', 'Row', 'Column']
    elif plotType == 'frequency':
        parameterTypes = ['Color','Row','Column']
        pass
    else:
        parameterTypes = ['Marker','Color','Size','Row','Column','X Axis Values']

    withinFigureBoolean = pickle.load(open('semiProcessedData/wfBool.pkl','rb'))
    selectedLevels = []
    print('wfbool')
    print(withinFigureBoolean)
    for levelName,index in zip(labelDict,range(len(withinFigureBoolean))):
        if withinFigureBoolean[index]:
            selectedLevels.append(levelName)

    fig = plt.figure(figsize=(2*len(selectedLevels),0.75*len(parameterTypes)))
    plt.axis('off')
    plt.text(0.5, 1.1,'Which plotting parameters do you want to assign to your figure levels?',transform=plt.gca().transAxes,ha='center')
    i=0
    radiobuttons = []
    for levelName in selectedLevels:
        rectLength = 0.15*len(parameterTypes)
        rectWidth = (1- (0.02+0.01*len(selectedLevels)))/len(selectedLevels) 
        if rectLength > 0.8:
            rectLength = 0.8
        rax3 = plt.axes([0.01*(i+1)+rectWidth*i, 0.94-(buttonWidth+rectLength), rectWidth,rectLength])
        rax3.spines['bottom'].set_visible(False)
        rax3.spines['left'].set_visible(False)
        rax3.spines['right'].set_visible(False)
        plt.text(0.5, 1.01,levelName,ha='center',transform=plt.gca().transAxes)
        radiobuttons.append(RadioButtons(rax3,parameterTypes,activecolor='black'))
        i+=1
    print(selectedLevels)

    class Index3(object):
        def OK(self, event):
            radioValues = {}
            for radiobutton,levelName in zip(radiobuttons,selectedLevels):
                radioValues[radiobutton.value_selected] = levelName
            plt.close()
            with open('semiProcessedData/radioVals.pkl','wb') as f:
                pickle.dump(radioValues,f)
        def Quit(self, event):
            sys.exit(0)    

    callback = Index3()
    axOK = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
    axQuit = plt.axes([buttonXStart+buttonWidth+0.01,buttonYStart, buttonWidth, buttonLength])
    bOK = Button(axOK, 'OK')
    bOK.on_clicked(callback.OK)
    bQuit = Button(axQuit, 'Quit')
    bQuit.on_clicked(callback.Quit)
    plt.show()

def produceSubsettedDataFrames(folderName,secondPath,df,useModifiedDf):
    
    withinFigureBoolean = pickle.load(open('semiProcessedData/wfBool.pkl','rb'))
    specificValueBooleanList = pickle.load(open('semiProcessedData/svBoolList.pkl','rb'))

    #Load in dataframe for experiment
    #fulldf = np.log10(df)
    #Ask user which levels they want to separate plots by, and which level values within each level they would like to plot (plot individually option)
    print(df)
    fulldf = df.stack()

    allLevels = []

    figureLevels = []
    withinFigureLevels = []
    withinFigureSubsettedLevelValues = []
    figureSubsettedLevelValues = []
    
    levelValuesPlottedIndividually = []
    subsettedDfList = []
    subsettedDfListTitles = []
    for i in range(0,fulldf.index.nlevels):
        currentLevelName = fulldf.index.levels[i].name
        currentLevelValues = pd.unique(fulldf.index.get_level_values(currentLevelName))
        
        if withinFigureBoolean[i]:
            allLevels.append('wf')
            withinFigureLevels.append(currentLevelName)
            numSpecified = 0
            for specificBoolean in specificValueBooleanList[i]:
                if specificBoolean:
                    numSpecified+=1
            if numSpecified < len(list(currentLevelValues)):
                tempSubsettedLevelValues = []
                for j in range(len(currentLevelValues)):
                    if specificValueBooleanList[i][j]:
                        tempSubsettedLevelValues.append(currentLevelValues[j])
                        levelValuesPlottedIndividually.append(str(currentLevelValues[j]))
                withinFigureSubsettedLevelValues.append(tempSubsettedLevelValues)
            else:
                withinFigureSubsettedLevelValues.append(slice(None))
        else:
            print('figure')
            allLevels.append('f')
            figureLevels.append(currentLevelName)
            numSpecified = 0
            for specificBoolean in specificValueBooleanList[i]:
                if specificBoolean:
                    numSpecified+=1
            if numSpecified < len(list(currentLevelValues)):
                tempSubsettedLevelValues = []
                for j in range(len(currentLevelValues)):
                    if specificValueBooleanList[i][j]:
                        tempSubsettedLevelValues.append(currentLevelValues[j])
                        levelValuesPlottedIndividually.append(str(currentLevelValues[j]))
                figureSubsettedLevelValues.append(tempSubsettedLevelValues)
            else:
                figureSubsettedLevelValues.append(currentLevelValues)
    
    figureCombos = itertools.product(*figureSubsettedLevelValues)
    for figureCombo in figureCombos:
        subsettingTuple = []
        subsettingTitle = []
        j=0
        k=0
        for level in allLevels:
            if level == 'f':
                subsettingTuple.append(figureCombo[j])
                subsettingTitle.append(str(figureCombo[j]))
                j+=1
            else:
                subsettingTuple.append(withinFigureSubsettedLevelValues[k])
                k+=1
        try:
            copiedDf = fulldf.loc[tuple(subsettingTuple)].copy()
        except:
            print('key error: ')
            print(subsettingTuple)
        else:
            print(copiedDf)
            subsettedvalues = copiedDf.values
            subsettedindex = copiedDf.index.remove_unused_levels()
            subsetteddf = pd.DataFrame(subsettedvalues,subsettedindex)
            subsettedDfList.append(subsetteddf)
            subsettedDfListTitles.append(subsettingTitle)

    return subsettedDfList,subsettedDfListTitles,withinFigureLevels,figureLevels,levelValuesPlottedIndividually

def assignParametersToLevels(df,levelsToPlot,plotType):

    legendParameterToLevelNameDict = pickle.load(open('semiProcessedData/radioVals.pkl','rb'))
    subprocess.run(['rm','semiProcessedData/radioVals.pkl'])
    subprocess.run(['rm','semiProcessedData/svBoolList.pkl'])
    subprocess.run(['rm','semiProcessedData/wfBool.pkl'])
    return legendParameterToLevelNameDict

def changeLevelNamesForPlotting(plottingDf,plotType,dataType,subsettedDfTitle,legendParameterToLevelNameDict,figureParameters,kwargs):
    if plotType == 'ordered':
        if legendParameterToLevelNameDict['X Axis Values'] == 'Time':
            kwargs['x'] = 'Time (hours)'
            plottingDf.index.set_names('Time (hours)',level='Time',inplace=True)
    elif legendParameterToLevelNameDict['Order'] == 'Time':
        kwargs['x'] = 'Time (hours)'
        #kwargs['order'] = list(pd.unique(subsettedDf.index.get_level_values('Time')))
        plottingDf.index.set_names('Time (hours)',level='Time',inplace=True)
    if dataType == 'cyt':
        plottingDf.columns = ['Concentration (nM)']
        kwargs['y'] = 'Concentration (nM)'
    else:
        if 'Statistic' in figureParameters:
            plottingDf.columns = [subsettedDfTitle[figureParameters.index('Statistic')]]
            kwargs['y'] = subsettedDfTitle[figureParameters.index('Statistic')]
    return plottingDf,kwargs

def createFacePlotName(folderName,dataType,plotType,subPlotType,legendParameterToLevelNameDict,subsettedDfTitle,levelsPlottedIndividually,useModifiedDf):
    delimiter1 = '-'
    delimiter2 = ','
    
    legendParameterString = delimiter2.join(list(legendParameterToLevelNameDict.keys()))
    levelNameString = delimiter2.join(list(legendParameterToLevelNameDict.values()))
    figureLevelNameString = delimiter2.join(subsettedDfTitle)

    if len(levelsPlottedIndividually) == 0:
        individualLevelString = 'all'
    else:
        individualLevelString = delimiter2.join(levelsPlottedIndividually)
    if useModifiedDf:
        modifiedString = '-modified'
    else:
        modifiedString = ''

    if len(subsettedDfTitle) == 0:
        initialString = delimiter1.join([subPlotType,dataType,folderName,legendParameterString,levelNameString,individualLevelString])
    else:
        initialString = delimiter1.join([subPlotType,dataType,folderName,legendParameterString,levelNameString,individualLevelString,figureLevelNameString])
    fullTitleString = initialString+modifiedString
    return fullTitleString

def plotFacetedFigures(folderName,plotType,subPlotType,dataType,subsettedDfList,subsettedDfListTitles,legendParameterToLevelNameDict,figureParameters,levelsPlottedIndividually,useModifiedDf):
    print(legendParameterToLevelNameDict)
    #Has issues with repeated values (aka CD54 shows up in TCells and APCs) 
    for subsettedDf,subsettedDfTitle in zip(subsettedDfList,subsettedDfListTitles):
        kwargs = {}
        for parameter in legendParameterToLevelNameDict:
            currentLevel = legendParameterToLevelNameDict[parameter]
            if parameter == 'Color':
                kwargs['hue'] = currentLevel
            elif parameter == 'Size':
                kwargs['size'] = currentLevel
            elif parameter == 'Marker':
                kwargs['style'] = currentLevel
            elif parameter == 'Order':
                kwargs['x'] = currentLevel
                kwargs['order'] = list(pd.unique(subsettedDf.index.get_level_values(currentLevel)))
            elif parameter == 'Column':
                kwargs['col'] = currentLevel
                unorderedLevelValues = list(pd.unique(subsettedDf.index.get_level_values(currentLevel)))
                if 'Concentration' in currentLevel:
                    levelValues = sortSINumerically(unorderedLevelValues,True,True)[0]
                else:
                    levelValues = unorderedLevelValues
                kwargs['col_order'] = levelValues
            elif parameter == 'Row':
                kwargs['row'] = currentLevel
                unorderedLevelValues = list(pd.unique(subsettedDf.index.get_level_values(currentLevel)))
                if 'Concentration' in currentLevel:
                    levelValues = sortSINumerically(unorderedLevelValues,True,True)[0]
                else:
                    levelValues = unorderedLevelValues
                kwargs['row_order'] = levelValues
            elif parameter == 'X Axis Values':
                #kwargs['x'] = currentLevel
                kwargs['x'] = 'Time (hrs)' 
                #kwargs['x'] = 'IFNg Pulse Concentration (nM)' 

        fullTitleString = createFacePlotName(folderName,dataType,plotType,subPlotType,legendParameterToLevelNameDict,subsettedDfTitle,levelsPlottedIndividually,useModifiedDf)
        plottingDf,kwargs = changeLevelNamesForPlotting(subsettedDf,plotType,dataType,subsettedDfTitle,legendParameterToLevelNameDict,figureParameters,kwargs)
        plottingDf = plottingDf.reset_index()
        if folderName == '20190323-B16OVAIFNgPulsed_OT1_Timeseries_4':
            currentLevelValues = list(subsettedDf.index.get_level_values('IFNgPulseConcentration'))
            sortedOldLevelValues,newLevelValues = sortSINumerically(currentLevelValues,False,True)
            scaledNewLevelValues = [float(i) * float(1e9) for i in newLevelValues]
            concDf = pd.DataFrame({'IFNg Pulse Concentration (nM)':scaledNewLevelValues})
            plottingDf = pd.concat([plottingDf,concDf],axis=1)
        if plotType == 'categorical':
            print(plottingDf['Concentration (nM)'])
            plottingDf['Concentration (nM)'] = np.log10(plottingDf['Concentration (nM)'])
            print(plottingDf['Concentration (nM)'])
            ax = sns.catplot(**kwargs,data=plottingDf,kind=subPlotType)
        elif plotType == 'frequency':
            gfiDf = pd.DataFrame({'GFI':df.values.ravel()})
            plottingDf = pd.concat([plottingDf,gfiDf],axis=1)
            g = sns.FacetGrid(plottingDf,**kwargs,legend_out=True,sharey=False)
            g.map(sns.distplot,'GFI',kde=False,bins=256)
            for ax in g.axes.flat:
                box = ax.get_position()
                ax.set_position([box.x0,box.y0,box.width*0.85,box.height])
            plt.legend(loc='upper left',bbox_to_anchor=(1,0.5))
        else:
            if 'style' not in kwargs.keys():
                styleValues = []
                for i in range(plottingDf.shape[0]):
                    styleValues.append(' ')
                styleDf = pd.DataFrame({'.':styleValues})
                plottingDf = pd.concat([plottingDf,styleDf],axis=1)
                kwargs['style'] = '.'
            ax = sns.relplot(**kwargs,data=plottingDf,markers=True,kind=subPlotType)
        if dataType == 'cyt' or 'GFI' in subsettedDfTitle:
            if 'GFI' in subsettedDfTitle:
                print('holeefuk')
                ax.fig.get_axes()[0].set_yscale('symlog',linthreshx=10)
            else:
                print('wat')
                print(kwargs)
                #plt.semilogy(nonposy='clip')
                #ax.fig.get_axes()[0].set_yscale('log')
                #symlogfloor = 10**(math.floor(np.log10(np.amin(plottingDf['Concentration (nM)']))))
                #print(symlogfloor)
                #ax.fig.get_axes()[0].set_yscale('symlog',linthreshx=symlogfloor)
        if dataType == 'ordered' and 'Concentration' in kwargs['x']:
            print('wat2')
            ax.fig.get_axes()[0].set_xscale('symlog',linthreshx=1e-4)
        #ax.fig.get_axes()[0].set_xscale('symlog',linthreshx=1e-4)
        plt.subplots_adjust(top=0.94)
        plt.suptitle('-'.join(subsettedDfTitle),fontsize = 'large',fontweight='bold')
        plt.savefig('fullyProcessedFigures/'+fullTitleString+'.png',bbox_inches='tight')
        print(fullTitleString+' plot saved')
        plt.clf()

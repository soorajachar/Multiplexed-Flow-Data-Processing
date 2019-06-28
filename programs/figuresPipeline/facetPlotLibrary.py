#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on fri sep 7 13:12:56 2018

@author: acharsr
"""

import matplotlib
from matplotlib import pyplot as plt
from matplotlib import markers
import numpy as np
import seaborn as sns
import pandas as pd
import pickle,os,math,sys,itertools,re
from matplotlib.widgets import RadioButtons,Button,CheckButtons,TextBox
#from facetPlotHeatmaps import heatmapSpecific_GUIWindow,label_index,label_columns,label_headers,draw_borders
from facetPlotHeatmaps import heatmapSpecific_GUIWindow,draw_faceted_heatmap,returnHeatmapAspectRatios
from facetPlotScatterLineBar import scatterLineBarSpecific_GUIWindow,scaleXandY2D
from facetPlotHistogramKDE import setXandYTicks
sys.path.insert(0, '../dataProcessing/')
from miscFunctions import sortSINumerically,reindexDataFrame
import subprocess
from matplotlib.colors import LogNorm,SymLogNorm
from processProliferationData import returnGates,returnTicks

#Button width conserved across gui figures
buttonWidth = 0.1/2
buttonLength = 0.075/2
buttonXStart = 0.5-(0.01+buttonWidth)
buttonYStart = 0.01

#Bandwidth for kde plots; changes per machine for some reason so need to be modified often
#bandwidth = 1.5
bandwidth = 0.15
col_wrap_min = 6
#shareYVar = True
titleBool = True
secondPathBool = False

#Get level names and values into an easily accessible dictionary
def createLabelDict(df,dataType):
    fulldf = df.stack()
    labelDict = {}
    for i in range(fulldf.index.nlevels):
        levelName = fulldf.index.levels[i].name
        if levelName != 'Event':
            labelDict[levelName] = list(pd.unique(fulldf.index.get_level_values(levelName)))
    return labelDict

def levelsInFigure_GUIWindow(labelDict):
    #Grab levels that will be used within a figure
    fig = plt.figure()
    plt.axis('off')

    #plt.text(0.5, 1,'Which levels names do you want to\nbe included within this figure?',transform=plt.gca().transAxes,ha='center')
    rax = plt.axes([0.4, 0.5, 0.05*len(labelDict.keys()),0.05*len(labelDict.keys())])
    check = CheckButtons(rax, labelDict.keys(),actives=[True]*len(labelDict.keys()))
    plt.text(0.5, 1.2,'Which levels names do you want to be included within this figure?',ha='center')
    rax.spines['bottom'].set_visible(False)
    rax.spines['left'].set_visible(False)
    rax.spines['right'].set_visible(False)
    rax.spines['top'].set_visible(False)

    class GUIButtons(object):
        def OKcheck1(self, event):
            withinFigureBoolean = check.get_status()
            plt.close()
            with open('semiProcessedData/gui-wfBool.pkl','wb') as f:
                pickle.dump(withinFigureBoolean,f)
        def Quit(self, event):
            sys.exit(0)

    callback = GUIButtons()
    axOK = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
    axQuit = plt.axes([buttonXStart+buttonWidth+0.01,buttonYStart, buttonWidth, buttonLength])
    bOK = Button(axOK, 'OK')
    bOK.on_clicked(callback.OKcheck1)
    bQuit = Button(axQuit, 'Quit')
    bQuit.on_clicked(callback.Quit)
    plt.show()

def levelValuesInFigure_GUIWindow(labelDict):
    maxLabelLength = 0
    for labelName in labelDict:
        if len(labelDict[labelName]) > maxLabelLength:
            maxLabelLength = len(labelDict[labelName])
    fig = plt.figure()
    plt.axis('off')
    plt.text(0.5, 1.1,'Which specific level values do you want to include in the figure?',transform=plt.gca().transAxes,ha='center')
    i=0
    checkbuttons = []
    for levelName in labelDict:
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

    class GUIButtons2(object):
        def OKcheck2(self, event):
            specificValueBooleanList = []
            for checkbutton in checkbuttons:
                specificValueBooleanList.append(checkbutton.get_status())
            plt.close()
            with open('semiProcessedData/gui-svBoolList.pkl','wb') as f:
                pickle.dump(specificValueBooleanList,f)
        def Quit(self, event):
            sys.exit(0)

    callback = GUIButtons2()
    axOK = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
    axQuit = plt.axes([buttonXStart+buttonWidth+0.01,buttonYStart, buttonWidth, buttonLength])
    bOK = Button(axOK, 'OK')
    bOK.on_clicked(callback.OKcheck2)
    bQuit = Button(axQuit, 'Quit')
    bQuit.on_clicked(callback.Quit)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.show()

#Ask user which levels to assign to which of the four parameters that can be changed in a bar/point plot: color (c), order (o), x position of subplot (x), y position of subplot (y)
#Order essentially is which level is plotted on the x axis, as bar plots are always categorical
def levelsToPlottingParameters_GUIWindow(labelDict,plotType,subPlotType,dataType):
    if plotType == 'categorical':
        parameterTypes = ['Color','Order', 'Row', 'Column','None']
    elif plotType == '1d':
        parameterTypes = ['Color','Row','Column']
    elif plotType == '3d':
        if subPlotType == 'heatmap':
            parameterTypes = ['Row','Column','X Axis Values','Y Axis Values']
    else:
        parameterTypes = ['Marker','Color','Size','Row','Column','X Axis Values','None']

    withinFigureBoolean = pickle.load(open('semiProcessedData/gui-wfBool.pkl','rb'))
    selectedLevels = []
    for levelName,index in zip(labelDict,range(len(withinFigureBoolean))):
        if withinFigureBoolean[index]:
            selectedLevels.append(levelName)

    #2*4*1.5x4*3*0.45 is perfect ratio: 12x6
    fig = plt.figure(figsize=(3*len(parameterTypes),1.2*len(parameterTypes)))
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

    class GUIButtons3(object):
        def OKradio3(self, event):
            radioValues = {}
            for radiobutton,levelName in zip(radiobuttons,selectedLevels):
                if radiobutton.value_selected in radioValues.keys():
                    if isinstance(radioValues[radiobutton.value_selected], (list,)):
                        radioValues[radiobutton.value_selected] = radioValues[radiobutton.value_selected]+[levelName]
                    else:
                        radioValues[radiobutton.value_selected] = [radioValues[radiobutton.value_selected]]+[levelName]
                else:
                    radioValues[radiobutton.value_selected] = levelName
            plt.close()
            with open('semiProcessedData/gui-radioVals.pkl','wb') as f:
                pickle.dump(radioValues,f)
        def Quit(self, event):
            sys.exit(0)

    callback = GUIButtons3()
    axOK = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
    axQuit = plt.axes([buttonXStart+buttonWidth+0.01,buttonYStart, buttonWidth, buttonLength])
    bOK = Button(axOK, 'OK')
    bOK.on_clicked(callback.OKradio3)
    bQuit = Button(axQuit, 'Quit')
    bQuit.on_clicked(callback.Quit)
    plt.show()

def produceSubsettedDataFrames(df):

    withinFigureBoolean = pickle.load(open('semiProcessedData/gui-wfBool.pkl','rb'))
    specificValueBooleanList = pickle.load(open('semiProcessedData/gui-svBoolList.pkl','rb'))

    #Load in dataframe for experiment
    fulldf = df.stack()

    allLevels = []

    figureLevels = []
    withinFigureSubsettedLevelValues = []
    figureSubsettedLevelValues = []

    levelValuesPlottedIndividually = []
    subsettedDfList = []
    subsettedDfListTitles = []
    k=0
    for i in range(0,fulldf.index.nlevels):
        currentLevelName = fulldf.index.levels[i].name
        if currentLevelName == 'Event':
            allLevels.append('wf')
            withinFigureSubsettedLevelValues.append(slice(None))
        else:
            currentLevelValues = pd.unique(fulldf.index.get_level_values(currentLevelName))
            if withinFigureBoolean[k]:
                allLevels.append('wf')
                numSpecified = 0
                for specificBoolean in specificValueBooleanList[k]:
                    if specificBoolean:
                        numSpecified+=1
                if numSpecified < len(list(currentLevelValues)):
                    tempSubsettedLevelValues = []
                    for j in range(len(currentLevelValues)):
                        if specificValueBooleanList[k][j]:
                            tempSubsettedLevelValues.append(currentLevelValues[j])
                            levelValuesPlottedIndividually.append(str(currentLevelValues[j]))
                    withinFigureSubsettedLevelValues.append(tempSubsettedLevelValues)
                else:
                    withinFigureSubsettedLevelValues.append(slice(None))
            else:
                allLevels.append('f')
                figureLevels.append(currentLevelName)
                numSpecified = 0
                for specificBoolean in specificValueBooleanList[k]:
                    if specificBoolean:
                        numSpecified+=1
                if numSpecified < len(list(currentLevelValues)):
                    tempSubsettedLevelValues = []
                    for j in range(len(currentLevelValues)):
                        if specificValueBooleanList[k][j]:
                            tempSubsettedLevelValues.append(currentLevelValues[j])
                            levelValuesPlottedIndividually.append(str(currentLevelValues[j]))
                    figureSubsettedLevelValues.append(tempSubsettedLevelValues)
                else:
                    figureSubsettedLevelValues.append(currentLevelValues)
            k+=1

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
    return subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually

def createFacetPlotName(folderName,dataType,plotType,subPlotType,legendParameterToLevelNameDict,subsettedDfTitle,levelsPlottedIndividually,useModifiedDf,axisScaling):
    delimiter1 = '-'
    delimiter2 = ','

    legendParameterString = delimiter2.join(list(legendParameterToLevelNameDict.keys()))

    flattened_list = []
    for val in legendParameterToLevelNameDict.values():
        if isinstance(val, (list,)):
            for val2 in val:
                flattened_list.append(val2)
        else:
            flattened_list.append(val)

    levelNameString = delimiter2.join(flattened_list)
    figureLevelNameString = delimiter2.join(subsettedDfTitle)

    if len(levelsPlottedIndividually) == 0:
        individualLevelString = 'all'
    else:
        individualLevelString = delimiter2.join(levelsPlottedIndividually)
    if useModifiedDf:
        modifiedString = '-modified'
    else:
        modifiedString = ''

    axisScalingStringList = []
    for axis in axisScaling:
        if 'X' in axis:
            if axisScaling[axis] == 'Logarithmic':
                axisScalingStringList.append('logX')
            elif axisScaling[axis] == 'Biexponential':
                axisScalingStringList.append('biexpX')
            else:
                axisScalingStringList.append('linX')
        else:
            if axisScaling[axis] == 'Logarithmic':
                axisScalingStringList.append('logY')
            elif axisScaling[axis] == 'Biexponential':
                axisScalingStringList.append('biexpY')
            else:
                axisScalingStringList.append('linY')
    axisScalingString = delimiter2.join(axisScalingStringList)

    if len(subsettedDfTitle) == 0:
        initialString = delimiter1.join([subPlotType,dataType,folderName,legendParameterString,levelNameString,axisScalingString])
    else:
        initialString = delimiter1.join([subPlotType,dataType,folderName,legendParameterString,levelNameString,figureLevelNameString,axisScalingString])
    fullTitleString = initialString+modifiedString
    if '/' in fullTitleString:
        fullTitleString = fullTitleString.replace("/", "_")
    if '.' in fullTitleString:
        fullTitleString = fullTitleString.replace(".", "_")
    if ' ' in fullTitleString:
        fullTitleString = fullTitleString.replace(" ", "_")
    return fullTitleString

def plotFacetedFigures(folderName,plotType,subPlotType,dataType,subsettedDfList,subsettedDfListTitles,figureParameters,levelsPlottedIndividually,useModifiedDf,fulldf):
    legendParameterToLevelNameDict = pickle.load(open('semiProcessedData/gui-radioVals.pkl','rb'))
    plotOptions = pickle.load(open('semiProcessedData/gui-plotOptions.pkl','rb'))

    subprocess.run(['rm','semiProcessedData/gui-radioVals.pkl'])
    subprocess.run(['rm','semiProcessedData/gui-plotOptions.pkl'])
    subprocess.run(['rm','semiProcessedData/gui-wfBool.pkl'])
    subprocess.run(['rm','semiProcessedData/gui-svBoolList.pkl'])
    
    #Has issues with repeated values (aka CD54 shows up in TCells and APCs)
    for subsettedDf,subsettedDfTitle in zip(subsettedDfList,subsettedDfListTitles):
        #Assign all levels to plot parameters in catplot/relplot; reassign x/y axis level names to desired x/y axis titles
        kwargs = {}
        facetgridkwargs = {}
        for parameter in legendParameterToLevelNameDict:
            if parameter == 'Y Axis Values':
                if subPlotType == 'heatmap':
                    kwargs['y'] = legendParameterToLevelNameDict[parameter]
            else:
                currentLevel = legendParameterToLevelNameDict[parameter]
            if parameter == 'Color':
                kwargs['hue'] = currentLevel
            elif parameter == 'Size':
                kwargs['size'] = currentLevel
            elif parameter == 'Marker':
                kwargs['style'] = currentLevel
            elif parameter == 'Order':
                subsettedDf.index.set_names(plotOptions['axisTitles']['X Axis'],level=currentLevel,inplace=True)
                kwargs['x'] = plotOptions['axisTitles']['X Axis']
            elif parameter == 'Column':
                kwargs['col'] = currentLevel
                unorderedLevelValues = list(pd.unique(subsettedDf.index.get_level_values(currentLevel)))
                if 'Concentration' in currentLevel:
                    levelValues = sortSINumerically(unorderedLevelValues,True,True)[0]
                else:
                    levelValues = unorderedLevelValues
                kwargs['col_order'] = levelValues
                facetgridkwargs['col'] = kwargs['col']
                facetgridkwargs['col_order'] = kwargs['col_order']
            elif parameter == 'Row':
                kwargs['row'] = currentLevel
                unorderedLevelValues = list(pd.unique(subsettedDf.index.get_level_values(currentLevel)))
                if 'Concentration' in currentLevel:
                    levelValues = sortSINumerically(unorderedLevelValues,True,True)[0]
                else:
                    levelValues = unorderedLevelValues
                kwargs['row_order'] = levelValues
                facetgridkwargs['row'] = kwargs['row']
                facetgridkwargs['row_order'] = kwargs['row_order']
            elif parameter == 'X Axis Values':
                if len(plotOptions['axisTitles']['X Axis']) > 1:
                    subsettedDf.index.set_names(plotOptions['axisTitles']['X Axis'],level=currentLevel,inplace=True)
                    kwargs['x'] = plotOptions['axisTitles']['X Axis']
                else:
                    kwargs['x'] = currentLevel
            else:
                pass
        #Assign y axis (or color bar axis) parameter
        if subPlotType not in ['kde','histogram']:
            if subPlotType != 'heatmap':
                if 'Statistic' in figureParameters:
                    subsettedDf.columns = [subsettedDfTitle[figureParameters.index('Statistic')]]
                    kwargs['y'] = subsettedDfTitle[figureParameters.index('Statistic')]
                else:
                    subsettedDf.columns = [plotOptions['axisTitles']['Y Axis']]
                    kwargs['y'] = plotOptions['axisTitles']['Y Axis']
            else:
                if 'Statistic' in figureParameters:
                    subsettedDf.columns = [subsettedDfTitle[figureParameters.index('Statistic')]]
                    kwargs['z'] = subsettedDfTitle[figureParameters.index('Statistic')]
                else:
                    subsettedDf.columns = [plotOptions['axisTitles']['Colorbar Axis']]
                    kwargs['z'] = plotOptions['axisTitles']['Colorbar Axis']
        else:
            pass

        if dataType == 'singlecell':
            plottingDf = subsettedDf.stack().to_frame('GFI')
        else:
            plottingDf = subsettedDf.copy()

        #Converts wide form dataframe into long form required for cat/relplot
        plottingDf = plottingDf.reset_index()
        #Use plot options file to initialize numeric x axis ordering
        #NEED TO GET WORKING WITH HEATMAPS
        if subPlotType != 'heatmap':
            if plotOptions['numericX'][0]:
                currentLevelValues = list(plottingDf[kwargs['x']])
                sortedOldLevelValues,newLevelValues = sortSINumerically(currentLevelValues,False,True)
                #Need to interpret parenthetical units for x to get 1e9
                if 'M)' in kwargs['x']:
                    s = kwargs['x']
                    units = '1'+s[s.find("(")+1:s.find(")")]
                    scaledSortedUnits,sortedUnits = sortSINumerically([units],False,True)
                else:
                    sortedUnits = [1]
                scaledNewLevelValues = [float(i) / float(sortedUnits[0]) for i in newLevelValues]
                plottingDf[kwargs['x']] = scaledNewLevelValues
                if plotType == 'categorical':
                    kwargs['order'] = scaledNewLevelValues
            else:
                if plotType == 'categorical':
                    kwargs['order'] = list(pd.unique(subsettedDf.index.get_level_values(kwargs['x'])))
        else:
            print('wat')
            pass
            """
            for numericAxis in ['numericX','numericY']:
                axisLetter = numericAxis[-1]
                lowercaseAxis = axisLetter.lower()
                if plotOptions[numericAxis][0]:
                    currentLevelValues = list(plottingDf[kwargs[lowercaseAxis]])
                    sortedOldLevelValues,newLevelValues = sortSINumerically(currentLevelValues,False,True)
                    #Need to interpret parenthetical units for x to get 1e9
                    if 'M)' in kwargs[lowercaseAxis]:
                        s = kwargs[lowercaseAxis]
                        units = '1'+s[s.find("(")+1:s.find(")")]
                        scaledSortedUnits,sortedUnits = sortSINumerically([units],False,True)
                    else:
                        sortedUnits = [1]
                    scaledNewLevelValues = [float(i) / float(sortedUnits[0]) for i in newLevelValues]
                    plottingDf[kwargs[lowercase]] = scaledNewLevelValues
            """
        #Make sure there are not more than 6 plots in a row (if no row facet variable)
        if 'row' not in facetgridkwargs.keys():
            if 'col' in facetgridkwargs.keys():
                colwrap = min([len(kwargs['col_order']),col_wrap_min])
                kwargs['col_wrap'] = colwrap

        print(plottingDf)
        print(kwargs)
        fullTitleString = createFacetPlotName(folderName,dataType,plotType,subPlotType,legendParameterToLevelNameDict,subsettedDfTitle,levelsPlottedIndividually,useModifiedDf,plotOptions['axisScaling'])
        if len(subsettedDf.index) > 0:
            plotSubsettedFigure(subsettedDf,plottingDf,kwargs,facetgridkwargs,plotType,subPlotType,dataType,fullTitleString,plotOptions,subsettedDfTitle)

def plotSubsettedFigure(subsettedDf,plottingDf,kwargs,facetgridkwargs,plotType,subPlotType,dataType,fullTitleString,plotOptions,subsettedDfTitle):
    print(plotOptions)
    #Use appropriate facet plot
    #1d plots: #Histograms,KDEs
    if plotType == '1d':
        #Will need to update to make sure it pulls from y axis variable
        if subPlotType == 'kde':
            fg = sns.FacetGrid(plottingDf,sharey=False,legend_out=True,**kwargs)
            fg.map(sns.kdeplot,'GFI',shade=True,bw=bandwidth)
        elif subPlotType == 'histogram':
            fg = sns.FacetGrid(plottingDf,sharey=False,legend_out=True,**kwargs)
            fg.map(sns.distplot,'GFI',bins=256,kde=False)

        if dataType == 'singlecell':
            setXandYTicks(plottingDf,fg,subPlotType,kwargs)
        fg.add_legend()
    #2d plots: categorical plots (bar, point), ordered plots (line,scatter)
    elif plotType in ['categorical','ordered']:
        if plotType == 'categorical':
            if subPlotType == 'stripbar':
                fg = sns.catplot(**kwargs,data=plottingDf,kind='bar',alpha=0.8)
                secondkwargs = kwargs.copy()
                for key in ['row','col','col_order','row_order']:
                    if key in secondkwargs.keys():
                        secondkwargs.pop(key,None)
                axisIndex  = 0
                if 'row' in kwargs and 'col' in kwargs:
                    for rowVal in pd.unique(plottingDf[kwargs['row']]):
                        for colVal in pd.unique(plottingDf[kwargs['col']]):
                            secondPlottingDf = plottingDf[plottingDf[kwargs['row']] == rowVal]
                            secondPlottingDf = secondPlottingDf[secondPlottingDf[kwargs['col']] == colVal]
                            sns.stripplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                            if rowVal != pd.unique(plottingDf[kwargs['row']])[-1]:
                                fg.fig.axes[axisIndex].set_xlabel('')
                            axisIndex+=1
                else:
                    if 'row' in kwargs:
                        for rowVal in pd.unique(plottingDf[kwargs['row']]):
                            secondPlottingDf = plottingDf[plottingDf[kwargs['row']] == rowVal]
                            sns.stripplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                            if rowVal != pd.unique(plottingDf[kwargs['row']])[-1]:
                                fg.fig.axes[axisIndex].set_xlabel('')
                            axisIndex+=1
                    elif 'col' in kwargs:
                        for colVal in pd.unique(plottingDf[kwargs['col']]):
                            secondPlottingDf = plottingDf[plottingDf[kwargs['col']] == colVal]
                            sns.stripplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                            axisIndex+=1
                    else:
                        sns.stripplot(**secondkwargs,data=plottingDf,ax=fg.fig.axes[axisIndex])
            else:
                fg = sns.catplot(**kwargs,data=plottingDf,kind=subPlotType)
        elif plotType == 'ordered':
            shareYVar = 'row'
            #Make sure there are markers at each column variable
            if 'style' not in kwargs.keys():
                fg = sns.relplot(data=plottingDf,marker='o',kind=subPlotType,facet_kws={'sharey':shareYVar},ci=False,**kwargs)
            else:
                fg = sns.relplot(data=plottingDf,markers=True,kind=subPlotType,facet_kws={'sharey':shareYVar},ci=False,**kwargs)
        scaleXandY2D(fg,plotOptions)
    #3d plots: heatmaps and 3d scatter/lineplots
    elif plotType == '3d':
        if subPlotType == 'heatmap':
            a,h = returnHeatmapAspectRatios(subsettedDf,kwargs)
            print(facetgridkwargs)
            fg = sns.FacetGrid(plottingDf,height=h,aspect=a,gridspec_kws={"wspace":0.4},**facetgridkwargs)
            if plotOptions['axisScaling']['Colorbar Axis'] != 'Linear' and subsettedDf.values.min() > 0:
                logbool = True
                symlogbool = False
                log_norm = LogNorm(vmin=subsettedDf.min().min(), vmax=subsettedDf.max().max())
                sym_log_norm = ''
                cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(subsettedDf.min().min())), 1+math.ceil(math.log10(subsettedDf.max().max())))]
                lin_thresh = ''
            elif plotOptions['axisScaling']['Colorbar Axis'] != 'Linear' and subsettedDf.values.min() <= 0:
                logbool = False 
                symlogbool = True
                lin_thresh =plotOptions['linThreshold']['Colorbar Axis']
                sym_log_norm = SymLogNorm(linthresh=lin_thresh,linscale=1)
                log_norm = ''
                cbar_ticks= ''
            else:
                logbool = False
                symlogbool = False
                log_norm = ''
                sym_log_norm = ''
                cbar_ticks= ''
                lin_thresh = ''
            fg.map_dataframe(draw_faceted_heatmap, indexingdf=subsettedDf,xaxis=kwargs['x'], yaxis=kwargs['y'], zaxis=kwargs['z']
                    ,lognorm=log_norm,cbarticks=cbar_ticks,logarithmic=logbool,symlog=symlogbool,symlognorm=sym_log_norm,linthresh=lin_thresh)
            #STILL NEED TO WORK ON PUTTING COLORBARS ONLY AT THE END OF A ROW OF HEATMAPS
            for i in range(1,len(fg.fig.get_axes()),2):
                cbarax = fg.fig.get_axes()[i]
                cbarax.set_frame_on(True)

    #SupTitle
    #Make room for suptitle
    if 'row' in kwargs:
        if len(kwargs['row_order']) <=2:
            plt.subplots_adjust(top=0.9)
        else:
            plt.subplots_adjust(top=0.9+len(kwargs['row_order'])*0.005)
    else:
        plt.subplots_adjust(top=0.9)
    #Do not include placeholder celltypes in suptitle
    if 'NotApplicable' in subsettedDfTitle:
        subsettedDfTitle.remove('NotApplicable')
    if 'row' not in kwargs and 'col' not in kwargs:
        plt.suptitle('-'.join(subsettedDfTitle),fontsize = 'x-large',fontweight='bold')
    else:
        if titleBool:
            plt.suptitle('-'.join(subsettedDfTitle),fontsize = 'x-large',fontweight='bold')

    #Save figure
    fg.fig.savefig('fullyProcessedFigures/'+fullTitleString+'.png',bbox_inches='tight')
    fg.fig.savefig('../../outputFigures/'+fullTitleString+'.png',bbox_inches='tight')
    print(fullTitleString+' plot saved')
    plt.clf()

def facetPlottingGUI(df,plotType,subPlotType,dataType):
    labelDict = createLabelDict(df,dataType)
    levelsInFigure_GUIWindow(labelDict)
    levelValuesInFigure_GUIWindow(labelDict)
    levelsToPlottingParameters_GUIWindow(labelDict,plotType,subPlotType,dataType)
    if plotType in ['categorical','ordered','1d']:
        scatterLineBarSpecific_GUIWindow(labelDict,plotType,dataType)
    elif subPlotType == 'heatmap':
        heatmapSpecific_GUIWindow(labelDict,plotType,dataType)
    else:
        pass

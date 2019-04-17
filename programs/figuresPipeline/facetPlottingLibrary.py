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
from matplotlib.widgets import RadioButtons,Button,CheckButtons,TextBox
sys.path.insert(0, '../dataProcessing/')
from miscFunctions import sortSINumerically
import subprocess
plt.switch_backend('QT4Agg') #default on my system

def facetPlottingGUI(df,plotType,dataType):
    
    fulldf = df.stack()
    
    #Button width conserved across gui figures
    buttonWidth = 0.1/2
    buttonLength = 0.075/2
    buttonXStart = 0.5-(0.01+buttonWidth)
    buttonYStart = 0.01

    #Get level names and values into an easily accessible dictionary
    labelDict = {}
    for i in range(fulldf.index.nlevels):
        levelName = fulldf.index.levels[i].name
        labelDict[levelName] = list(pd.unique(fulldf.index.get_level_values(levelName)))

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
    if plotType == 'categorical':
        parameterTypes = ['Color','Order', 'Row', 'Column']
    elif plotType == 'frequency':
        parameterTypes = ['Color','Row','Column']
        pass
    else:
        parameterTypes = ['Marker','Color','Size','Row','Column','X Axis Values']

    withinFigureBoolean = pickle.load(open('semiProcessedData/gui-wfBool.pkl','rb'))
    selectedLevels = []
    for levelName,index in zip(labelDict,range(len(withinFigureBoolean))):
        if withinFigureBoolean[index]:
            selectedLevels.append(levelName)

    fig = plt.figure(figsize=(2*len(selectedLevels),0.45*len(parameterTypes)))
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
    
    #Ask scale to use for axes (both y and x if relplot; only y if categorical)
    if plotType == 'categorical' or plotType == 'frequency':
        axes = ['Y Axis']
    else:
        axes = ['X Axis','Y Axis']
    constantAxes = ['X Axis','Y Axis']

    axisScalingOptions = ['Linear','Logarithmic','Biexponential']

    fig = plt.figure(figsize=(8,2.5*len(axisScalingOptions)))
    plt.axis('off')
    radiobuttons = []
    axis_title_text_boxes = {}
    lin_thresh_text_boxes = {}
    for axis,i in zip(constantAxes,range(len(constantAxes))):
        #Add axis scaling radio buttons
        rectLength = 0.15*len(axisScalingOptions)
        rectWidth = (1- (0.02+0.01*len(constantAxes)))/len(constantAxes)
        
        rax3 = plt.axes([0.01*(i+1)+rectWidth*i, 0.94-(buttonWidth+rectLength), rectWidth,rectLength])
        plt.text(0.5, 1.01,axis,ha='center',transform=plt.gca().transAxes)
        rax3.spines['bottom'].set_visible(False)
        rax3.spines['left'].set_visible(False)
        rax3.spines['right'].set_visible(False)
        plt.tick_params(which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False)
        if(axis == 'X Axis'):
            if(len(constantAxes) == len(axes)):
                radiobuttons.append(RadioButtons(rax3,axisScalingOptions,activecolor='black'))
        else:
            radiobuttons.append(RadioButtons(rax3,axisScalingOptions,activecolor='black'))

        #Determine initial axis title for textboxes
        if 'Y' in axis:
            if dataType == 'cyt':
                initial_name = 'Concentration (nM)'
            else:
                if plotType == 'frequency':
                    initial_name = 'Count'
                else:
                    initial_name = ' '
        else:
            if plotType == 'ordered':
                radioValues = pickle.load(open('semiProcessedData/gui-radioVals.pkl','rb'))
                initial_name = radioValues['X Axis Values']
            elif plotType == 'frequency':
                initial_name = 'GFI'
            else:
                radioValues = pickle.load(open('semiProcessedData/gui-radioVals.pkl','rb'))
                initial_name = radioValues['Order']
        
        #Add in axis title text boxes
        axbox = plt.axes([0.1+0.5*i, 0.25, 0.35, 0.075])
        text_box = TextBox(axbox, 'Title: ', initial=initial_name)
        axis_title_text_boxes[axis] = text_box
        
        #Add in linear threshold (for biexponential scaling) textboxes
        if(axis == 'X Axis'):
            if(len(constantAxes) == len(axes)):
                linthreshbox = plt.axes([0.1+0.5*i, 0.35, 0.35, 0.075])
                text_box2 = TextBox(linthreshbox, 'Linear  \nThreshold: ', initial=' ')
                lin_thresh_text_boxes[axis] = text_box2
        else:
            linthreshbox = plt.axes([0.1+0.5*i, 0.35, 0.35, 0.075])
            text_box2 = TextBox(linthreshbox, 'Linear  \nThreshold: ', initial=' ')
            lin_thresh_text_boxes[axis] = text_box2
   
    if plotType == 'ordered':
        checkbuttons = []
        rax2 = plt.axes([0.05,0.09,0.15,0.15])
        rax2.spines['bottom'].set_visible(False)
        rax2.spines['left'].set_visible(False)
        rax2.spines['right'].set_visible(False)
        rax2.spines['top'].set_visible(False)
        checkbuttons.append(CheckButtons(rax2,['Sort X numerically'],actives=[False]))
    
    linThresholdValues = {}
    axisTitleValues = {'X Axis':initial_name,'Y Axis':' '}
    def submitAxisTitleX(text):
        axisTitleValues['X Axis'] = text
    def submitAxisTitleY(text):
        axisTitleValues['Y Axis'] = text
    axis_title_text_boxes['X Axis'].on_submit(submitAxisTitleX)
    axis_title_text_boxes['Y Axis'].on_submit(submitAxisTitleY)
    
    def submitLinThresholdX(text):
        linThresholdValues['X Axis'] = float(text)
    def submitLinThresholdY(text):
        linThresholdValues['Y Axis'] = float(text)
    if(len(constantAxes) == len(axes)):
        lin_thresh_text_boxes['X Axis'].on_submit(submitLinThresholdX)
        lin_thresh_text_boxes['Y Axis'].on_submit(submitLinThresholdY)
    else:
        lin_thresh_text_boxes['Y Axis'].on_submit(submitLinThresholdY)

    class GUIButtons4(object):
        def OKradiotext4(self, event):
            plotOptions = {}
            radioValues = {}
            for radiobutton,axis in zip(radiobuttons,axes):
                radioValues[axis] = radiobutton.value_selected
            if plotType == 'ordered':
                numericXBoolean = checkbuttons[0].get_status()
                plotOptions['numericX'] = numericXBoolean
            else:
                plotOptions['numericX'] = False
            plt.close()
            plotOptions['axisScaling'] = radioValues
            plotOptions['linThreshold'] = linThresholdValues
            plotOptions['axisTitles'] = axisTitleValues
            print(plotOptions)
            print(os.getcwd())
            with open('semiProcessedData/gui-plotOptions.pkl','wb') as f:
                pickle.dump(plotOptions,f)
        def Quit(self, event):
            sys.exit(0)    

    callback = GUIButtons4()
    axOK = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
    axQuit = plt.axes([buttonXStart+buttonWidth+0.01,buttonYStart, buttonWidth, buttonLength])
    bOK = Button(axOK, 'OK')
    bOK.on_clicked(callback.OKradiotext4)
    bQuit = Button(axQuit, 'Quit')
    bQuit.on_clicked(callback.Quit)
    plt.show()

def produceSubsettedDataFrames(folderName,secondPath,df,useModifiedDf):
    
    withinFigureBoolean = pickle.load(open('semiProcessedData/gui-wfBool.pkl','rb'))
    specificValueBooleanList = pickle.load(open('semiProcessedData/gui-svBoolList.pkl','rb'))

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

def createFacetPlotName(folderName,dataType,plotType,subPlotType,legendParameterToLevelNameDict,subsettedDfTitle,levelsPlottedIndividually,useModifiedDf,axisScaling):
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
    return fullTitleString

def plotFacetedFigures(folderName,plotType,subPlotType,dataType,subsettedDfList,subsettedDfListTitles,figureParameters,levelsPlottedIndividually,useModifiedDf):
    legendParameterToLevelNameDict = pickle.load(open('semiProcessedData/gui-radioVals.pkl','rb'))
    plotOptions = pickle.load(open('semiProcessedData/gui-plotOptions.pkl','rb'))
    
    subprocess.run(['rm','semiProcessedData/gui-*.pkl'])
    
    #Has issues with repeated values (aka CD54 shows up in TCells and APCs) 
    for subsettedDf,subsettedDfTitle in zip(subsettedDfList,subsettedDfListTitles):
        #Assign all levels to plot parameters in catplot/relplot; reassign x/y axis level names to desired x/y axis titles
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
            elif parameter == 'Row':
                kwargs['row'] = currentLevel
                unorderedLevelValues = list(pd.unique(subsettedDf.index.get_level_values(currentLevel)))
                if 'Concentration' in currentLevel:
                    levelValues = sortSINumerically(unorderedLevelValues,True,True)[0]
                else:
                    levelValues = unorderedLevelValues
                kwargs['row_order'] = levelValues
            elif parameter == 'X Axis Values':
                subsettedDf.index.set_names(plotOptions['axisTitles'][0],level=currentLevel,inplace=True)
                kwargs['x'] = plotOptions['axisTitles']['X Axis']
        #Assign y axis parameter
        if 'Statistic' in figureParameters:
            subsettedDf.columns = [subsettedDfTitle[figureParameters.index('Statistic')]]
            kwargs['y'] = subsettedDfTitle[figureParameters.index('Statistic')]
        else:
            plotOptions['axisTitles']['Y Axis']
            subsettedDf.columns = [plotOptions['axisTitles']['Y Axis']]
            kwargs['y'] = plotOptions['axisTitles']['Y Axis']
        
        plottingDf = subsettedDf.copy()
        #Converts wide form dataframe into long form required for cat/relplot
        plottingDf = plottingDf.reset_index()
        #Use plot options file to initialize plot parameters: axis scaling, numeric x axis ordering
        #Numeric X Axis Ordering
        if plotOptions['numericX']:
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
            kwargs['order'] = list(pd.unique(subsettedDf.index.get_level_values(kwargs['x'])))
            print(kwargs['order'])

        #Plot data using kwargs; uses catplot for categorical data and relplot for numeric data and a facet grid with a distplot for frequency data (histograms)
        if plotType == 'categorical':
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
        #Axis Scaling
        for axis in plotOptions['axisScaling']:
            if 'Y' in axis:
                if plotOptions['axisScaling'][axis] == 'Logarithmic':
                    ax.fig.get_axes()[0].set_yscale('log')
                elif plotOptions['axisScaling'][axis] == 'Biexponential':
                    ax.fig.get_axes()[0].set_yscale('symlog',linthreshx=plotOptions['linThreshold'][axis])
                else:
                    pass
            else:
                if plotOptions['axisScaling'][axis] == 'Logarithmic':
                    ax.fig.get_axes()[0].set_xscale('log')
                elif plotOptions['axisScaling'][axis] == 'Biexponential':
                    ax.fig.get_axes()[0].set_xscale('symlog',linthreshx=plotOptions['linThreshold'][axis])
                else:
                    pass
        plt.subplots_adjust(top=0.85)
        if 'NotApplicable' in subsettedDfTitle:
            subsettedDfTitle.remove('NotApplicable')

        plt.suptitle('-'.join(subsettedDfTitle),fontsize = 'x-large',fontweight='bold')
        fullTitleString = createFacetPlotName(folderName,dataType,plotType,subPlotType,legendParameterToLevelNameDict,subsettedDfTitle,levelsPlottedIndividually,useModifiedDf,plotOptions['axisScaling'])
        plt.savefig('fullyProcessedFigures/'+fullTitleString+'.png',bbox_inches='tight')
        print(fullTitleString+' plot saved')
        plt.clf()

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
from facetPlotHeatmaps import heatmapSpecific_GUIWindow,label_index,label_columns,label_headers,draw_borders
from facetPlotScatterLineBar import scatterLineBarSpecific_GUIWindow
sys.path.insert(0, '../dataProcessing/')
from miscFunctions import sortSINumerically,reindexDataFrame
import subprocess
from matplotlib.colors import LogNorm
    
#Button width conserved across gui figures
buttonWidth = 0.1/2
buttonLength = 0.075/2
buttonXStart = 0.5-(0.01+buttonWidth)
buttonYStart = 0.01

#Get level names and values into an easily accessible dictionary
def createLabelDict(df):
    fulldf = df.stack()
    labelDict = {}
    for i in range(fulldf.index.nlevels):
        levelName = fulldf.index.levels[i].name
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
def levelsToPlottingParameters_GUIWindow(labelDict,plotType,dataType):
    if plotType == 'categorical':
        parameterTypes = ['Color','Order', 'Row', 'Column']
    elif plotType == 'frequency':
        parameterTypes = ['Color','Row','Column']
    elif plotType == 'heatmap':
        parameterTypes = ['Row','Column','X Axis Values','Y Axis Values']
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

def produceSubsettedDataFrames(folderName,secondPath,df,useModifiedDf):
    
    withinFigureBoolean = pickle.load(open('semiProcessedData/gui-wfBool.pkl','rb'))
    specificValueBooleanList = pickle.load(open('semiProcessedData/gui-svBoolList.pkl','rb'))

    #Load in dataframe for experiment
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
    return fullTitleString

def plotFacetedFigures(folderName,plotType,subPlotType,dataType,subsettedDfList,subsettedDfListTitles,figureParameters,levelsPlottedIndividually,useModifiedDf,fulldf):
    legendParameterToLevelNameDict = pickle.load(open('semiProcessedData/gui-radioVals.pkl','rb'))
    plotOptions = pickle.load(open('semiProcessedData/gui-plotOptions.pkl','rb'))
    
    subprocess.run(['rm','semiProcessedData/gui-radioVals.pkl'])
    subprocess.run(['rm','semiProcessedData/gui-plotOptions.pkl'])
    subprocess.run(['rm','semiProcessedData/gui-wfBool.pkl'])
    subprocess.run(['rm','semiProcessedData/gui-svBoolList.pkl'])
    
    print(legendParameterToLevelNameDict)
    #Has issues with repeated values (aka CD54 shows up in TCells and APCs) 
    for subsettedDf,subsettedDfTitle in zip(subsettedDfList,subsettedDfListTitles):
        #Assign all levels to plot parameters in catplot/relplot; reassign x/y axis level names to desired x/y axis titles
        kwargs = {}
        heatmapkwargs = {}
        for parameter in legendParameterToLevelNameDict:
            if isinstance(legendParameterToLevelNameDict[parameter], (list,)):
                if parameter == 'Y Axis Values' and plotType == 'heatmap':
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
                heatmapkwargs['col'] = kwargs['col']
                heatmapkwargs['col_order'] = kwargs['col_order']
            elif parameter == 'Row':
                kwargs['row'] = currentLevel
                unorderedLevelValues = list(pd.unique(subsettedDf.index.get_level_values(currentLevel)))
                if 'Concentration' in currentLevel:
                    levelValues = sortSINumerically(unorderedLevelValues,True,True)[0]
                else:
                    levelValues = unorderedLevelValues
                kwargs['row_order'] = levelValues
                heatmapkwargs['row'] = kwargs['row']
                heatmapkwargs['row_order'] = kwargs['row_order']
            elif parameter == 'X Axis Values':
                print(currentLevel)
                subsettedDf.index.set_names(plotOptions['axisTitles']['X Axis'],level=currentLevel,inplace=True)
                kwargs['x'] = plotOptions['axisTitles']['X Axis']
        #Assign y axis (or color bar axis) parameter
        if plotType != 'heatmap':
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

        plottingDf = subsettedDf.copy()
        #Converts wide form dataframe into long form required for cat/relplot
        plottingDf = plottingDf.reset_index()
        #Use plot options file to initialize plot parameters: axis scaling, numeric x axis ordering
        
        #Numeric X Axis Ordering; NEED TO GET WORKING WITH HEATMAPS
        if plotType != 'heatmap':
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
        
        print(kwargs)
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
        elif plotType == 'heatmap':
            idx = pd.IndexSlice
            def draw_heatmap(data,xaxis,yaxis,zaxis, **kwargs):
                unsortedPivotedData = data.pivot_table(index=yaxis,columns=xaxis, values=zaxis)
                indexingList = []
                for name in fulldf.index.names:
                    if name not in unsortedPivotedData.index.names:
                        indexingList.append(list(pd.unique(fulldf.index.get_level_values(name)))[0])
                    else:
                        indexingList.append(list(pd.unique(unsortedPivotedData.index.get_level_values(name))))
                indexlist = []
                for levelIndex in range(len(indexingList)):
                    levelVal = indexingList[levelIndex]
                    if isinstance(levelVal, (list,)):
                        pass
                    else:
                        levelVal = [levelVal]
                    tempindexlist = []
                    for i in range(fulldf.shape[0]):
                        if fulldf.iloc[i,:].name[levelIndex] in levelVal:
                            tempindexlist.append(i)
                    indexlist.append(tempindexlist)
                unique_values = set(map(tuple, indexlist))
                indexlist = list(set.intersection(*map(set,unique_values)))
                indexdf = fulldf.iloc[indexlist,:]
                unusedLevelList = []
                for name in indexdf.index.names:
                    if len(pd.unique(indexdf.index.get_level_values(name))) == 1:
                        unusedLevelList.append(name)
                indexdf.reset_index(level=unusedLevelList, drop=True,inplace=True)
                d = reindexDataFrame(unsortedPivotedData,indexdf,False)
                plt.axis('off')

                if plotOptions['axisScaling']['Colorbar Axis'] == 'Logarithmic':
                    data = d.copy()
                    log_norm = LogNorm(vmin=data.min().min(), vmax=data.max().max())
                    if cbarBoolean:
                        cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(data.min().min())), 1+math.ceil(math.log10(data.max().max())))]
                    g = sns.heatmap(d, norm=log_norm,**kwargs)
                else:
                    g = sns.heatmap(d, **kwargs)

                ax1=plt.gca()
                label_index(ax1,d)
                draw_borders(g,d)
                label_columns(ax1,d)
                if plotOptions['IncludeLevelNames']:
                    label_headers(ax1,d)
            #16x16 heatmap
            basedim = 16
            scale=0.4
            hbase = 6
            abase = 1.25
            data = subsettedDf.copy()
            print(data)
            if 'row' in kwargs:
                data = data.xs(kwargs['row_order'][0],level=kwargs['row'])
            if 'col' in kwargs:
                data = data.xs(kwargs['col_order'][0],level=kwargs['col'])
            heatmapdf = data.pivot_table(index=kwargs['y'],columns=kwargs['x'], values=kwargs['z'])
            hstart = max(hbase+scale*hbase*(heatmapdf.index.size-basedim)/basedim,hbase)
            astart = max(abase+scale*abase*(heatmapdf.columns.size-basedim)/basedim,abase)
            if 'row_order' in kwargs:
                h = hstart*len(kwargs['row_order'])*scale
            else:
                h = hstart
            if 'col_order' in kwargs:
                a = astart*len(kwargs['col_order'])*scale
            else:
                a = astart
            print('h: '+str(h))
            print('a: '+str(a))
            fg = sns.FacetGrid(plottingDf,height=h,aspect=a,gridspec_kws={"wspace":0.4},**heatmapkwargs)
            if 'row' not in kwargs.keys() and 'col' not in kwargs.keys():
                cbarBoolean = True
                cbarBoolean = False
            else:
                cbarBoolean = False
            fg.map_dataframe(draw_heatmap, xaxis=kwargs['x'], yaxis=kwargs['y'], zaxis=kwargs['z'],cbar=cbarBoolean)
        else:
            if 'style' not in kwargs.keys():
                styleValues = []
                for i in range(plottingDf.shape[0]):
                    styleValues.append(' ')
                styleDf = pd.DataFrame({'.':styleValues})
                plottingDf = pd.concat([plottingDf,styleDf],axis=1)
                kwargs['style'] = '.'
            ax = sns.relplot(**kwargs,data=plottingDf,markers=True,kind=subPlotType)

        #Axis Scaling; NEED TO GET WORKING WITH HEATMAPS
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
                    if plotType != 'heatmap':
                        ax.fig.get_axes()[0].set_xscale('log')
                elif plotOptions['axisScaling'][axis] == 'Biexponential':
                    ax.fig.get_axes()[0].set_xscale('symlog',linthreshx=plotOptions['linThreshold'][axis])
                else:
                    pass

        #Make room for title
        if 'row' in kwargs:
            if len(kwargs['row_order']) <=2:
                plt.subplots_adjust(top=0.9)
            else:
                plt.subplots_adjust(top=0.95)
        else:
            plt.subplots_adjust(top=0.9)

        if 'NotApplicable' in subsettedDfTitle:
            subsettedDfTitle.remove('NotApplicable')

        plt.suptitle('-'.join(subsettedDfTitle),fontsize = 'x-large',fontweight='bold')
        fullTitleString = createFacetPlotName(folderName,dataType,plotType,subPlotType,legendParameterToLevelNameDict,subsettedDfTitle,levelsPlottedIndividually,useModifiedDf,plotOptions['axisScaling'])
        """
        for axis in plt.gcf().axes:
            xticks = axis.xaxis.get_major_ticks()
            for i in range(len(xticks)):
                if i % 4 != 0:
                    xticks[i].set_visible(False) 
        """
        fg.fig.savefig('fullyProcessedFigures/'+fullTitleString+'.png',bbox_inches='tight')
        #plt.savefig('fullyProcessedFigures/'+fullTitleString+'.png',bbox_inches='tight')
        print(fullTitleString+' plot saved')
        plt.clf()

def facetPlottingGUI(df,plotType,dataType):
    labelDict = createLabelDict(df)
    levelsInFigure_GUIWindow(labelDict)
    levelValuesInFigure_GUIWindow(labelDict)
    levelsToPlottingParameters_GUIWindow(labelDict,plotType,dataType)
    if plotType in ['categorical','ordered','frequency']:
        scatterLineBarSpecific_GUIWindow(labelDict,plotType,dataType)
    elif plotType == 'heatmap':
        heatmapSpecific_GUIWindow(labelDict,plotType,dataType)
    else:
        pass

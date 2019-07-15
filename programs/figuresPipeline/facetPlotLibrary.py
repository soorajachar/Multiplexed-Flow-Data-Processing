#!/usr/bin/env python3
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys,os,json,pickle,math,itertools
from facetPlot3D import draw_faceted_heatmap,returnHeatmapAspectRatios
from miscFunctions import sortSINumerically,reindexDataFrame
import subprocess
from matplotlib.colors import LogNorm,SymLogNorm
sys.path.insert(0, '../dataProcessing/')
from miscFunctions import returnGates,returnTicks
from operator import itemgetter
import facetPlot1D as fp1D
import facetPlotCategorical as fpCategorical
import facetPlot2D as fp2D
import facetPlot3D as fp3D

def produceSubsettedDataFrames(fulldf,withinFigureBoolean,specificValueBooleanList):
    print(withinFigureBoolean)
    print(specificValueBooleanList)
    #Get all possible subsetted indices
    figureSubsettedLevelValues = []
    withinFigureSubsettedLevelValues = []
    figureSubsettingLevels = []
    figureLevelNames = []
    figureLevelIndices = []
    levelValuesPlottedIndividually = []
    for levelIndex,currentLevelName in enumerate(fulldf.index.names):
        #Event level always has every value included (single cell)
        #Otherwise, check which level values in the level were selected by the user
        if currentLevelName != 'Event':
            currentLevelValues = pd.unique(fulldf.index.get_level_values(currentLevelName))
            levelValues = []
            #Go through each level value in the level; regardless of figure selection status, and add based on level value selection status
            for levelValue,specificBoolean in zip(currentLevelValues,specificValueBooleanList[levelIndex]):
                if specificBoolean:
                    levelValues.append(levelValue)
            #If we will include this level within the figure
            if withinFigureBoolean[levelIndex]:
                withinFigureSubsettedLevelValues.append(levelValues)
                if len(levelValues) == len(currentLevelValues):
                    for levelValue in levelValues:
                        levelValuesPlottedIndividually.append(levelValue)
            #Only need to add level values to figure list; will be xs'd out of the full dataframe in the subsetted, within figure dataframes
            else:
                figureSubsettedLevelValues.append(levelValues)
                figureLevelNames.append(currentLevelName)
                figureLevelIndices.append(levelIndex)
    if len(figureLevelIndices) > 0:
        #Get all row level values present in the dataframe
        rowList = []
        for row in range(fulldf.shape[0]):
            allCurrentLevelValues = fulldf.iloc[row,:].name
            currentLevelValues = itemgetter(*figureLevelIndices)(allCurrentLevelValues)
            if not isinstance(currentLevelValues,tuple):
                currentLevelValues = tuple([currentLevelValues])
            rowList.append(currentLevelValues)
        allPossibleSubsettingCombos = itertools.product(*figureSubsettedLevelValues)
        actualSubsettingCombos = []
        #From original dataframe; select all rows that appear in the all possible combination list 
        for subsettingCombo in allPossibleSubsettingCombos:
            if subsettingCombo in rowList:
                actualSubsettingCombos.append(subsettingCombo) 
        #Use these levels to cross section the fulldf, generating a list of subsetted dfs that will each have their own figure
        allPossibleSubsettedDfList = []
        for actualSubsettingCombo in actualSubsettingCombos:
            possibleSubsettedDf = fulldf.xs(actualSubsettingCombo, level=figureLevelNames)
            allPossibleSubsettedDfList.append(possibleSubsettedDf)
    else:
        actualSubsettingCombos = ['All']
        allPossibleSubsettedDfList = [fulldf]
    actualLevelValueDfList = []
    #Go through each subsetteddf, and only grab rows with level values that are selected
    for possibleSubsettedDf in allPossibleSubsettedDfList:
        allPossibleLevelValueCombos = itertools.product(*withinFigureSubsettedLevelValues)
        rowList = []
        for row in range(possibleSubsettedDf.shape[0]):
            if 'Event' in fulldf.index.names:
                allCurrentLevelValues = possibleSubsettedDf.iloc[row,:].name[:-1]
            else:
                allCurrentLevelValues = possibleSubsettedDf.iloc[row,:].name
            rowList.append(allCurrentLevelValues)
        actualLevelValueCombos = []
        #From original dataframe; select all rows that appear in the all possible level value combination list 
        levelValueRowList = []
        for levelValueCombo in allPossibleLevelValueCombos:
            if levelValueCombo in rowList:
                indices = [i for i, x in enumerate(rowList) if x == levelValueCombo]
                levelValueRowList+=indices
        actualLevelValueDf = possibleSubsettedDf.iloc[levelValueRowList,:]
        actualLevelValueDfList.append(actualLevelValueDf)
     
    #Remove columns from non postprocessed data
    if 'Dimension' not in fulldf.columns[0]:
        for i in range(len(actualLevelValueDfList)):
            actualLevelValueDfList[i].columns.name = ''
    print(actualLevelValueDfList) 
    print(actualSubsettingCombos) 
    return actualLevelValueDfList,actualSubsettingCombos,figureLevelNames,levelValuesPlottedIndividually

def createFacetPlotName(folderName,dataType,plotType,subPlotType,legendParameterToLevelNameDict,subsettedDfTitle,levelsPlottedIndividually,useModifiedDf,plotOptions):
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
    figureLevelNameString = delimiter2.join(list(map(str,subsettedDfTitle)))

    if len(levelsPlottedIndividually) == 0:
        individualLevelString = 'all'
    else:
        individualLevelString = delimiter2.join(list(map(str,levelsPlottedIndividually)))
    if useModifiedDf:
        modifiedString = '-modified'
    else:
        modifiedString = ''

    axisScalingStringList = []
    for axis in plotOptions:
        if 'X' in axis:
            if plotOptions[axis]['axisScaling'] == 'Logarithmic':
                axisScalingStringList.append('logX')
            elif plotOptions[axis]['axisScaling'] == 'Biexponential':
                axisScalingStringList.append('biexpX')
            else:
                axisScalingStringList.append('linX')
        else:
            if plotOptions[axis]['axisScaling'] == 'Logarithmic':
                axisScalingStringList.append('logY')
            elif plotOptions[axis]['axisScaling'] == 'Biexponential':
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

def plotFacetedFigures(folderName,plotType,subPlotType,dataType,subsettedDfList,subsettedDfListTitles,figureParameters,levelsPlottedIndividually,useModifiedDf,fulldf,plotOptions,legendParameterToLevelNameDict,addDistributionPoints,alternateTitle=''):

    #Has issues with repeated values (aka CD54 shows up in TCells and APCs)
    for subsettedDf,subsettedDfTitle in zip(subsettedDfList,subsettedDfListTitles):
        print(subsettedDfListTitles)
        print(subsettedDfTitle)
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
                subsettedDf.index.set_names(plotOptions['X']['axisTitle'],level=currentLevel,inplace=True)
                kwargs['x'] = plotOptions['X']['axisTitle']
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
                if len(plotOptions['X']['axisTitle']) > 1 and dataType != 'dr':
                    subsettedDf.index.set_names(plotOptions['X']['axisTitle'],level=currentLevel,inplace=True)
                    kwargs['x'] = plotOptions['X']['axisTitle']
                else:
                    kwargs['x'] = currentLevel
            else:
                pass
        #Assign y axis (or color bar axis) parameter
        #if subPlotType not in ['kde','histogram']:
        if dataType == 'dr':
            kwargs['y'] = plotOptions['Y']['axisTitle']
        else:
            if subPlotType != 'heatmap':
                if 'Statistic' in figureParameters:
                    subsettedDf.columns = [subsettedDfTitle[figureParameters.index('Statistic')]]
                    kwargs['y'] = subsettedDfTitle[figureParameters.index('Statistic')]
                else:
                    subsettedDf.columns = [plotOptions['Y']['axisTitle']]
                    kwargs['y'] = plotOptions['Y']['axisTitle']
            else:
                if 'Statistic' in figureParameters:
                    subsettedDf.columns = [subsettedDfTitle[figureParameters.index('Statistic')]]
                    kwargs['z'] = subsettedDfTitle[figureParameters.index('Statistic')]
                else:
                    subsettedDf.columns = [plotOptions['Colorbar']['axisTitle']]
                    kwargs['z'] = plotOptions['Colorbar']['axisTitle']
        #else:
        #pass
        
        if dataType == 'singlecell':
            plottingDf = subsettedDf.stack().to_frame('GFI')
        else:
            plottingDf = subsettedDf.copy()
        #Converts wide form dataframe into long form required for cat/relplot
        plottingDf = plottingDf.reset_index()
        print(plottingDf)
        #Use plot options file to initialize numeric x axis ordering
        #NEED TO GET WORKING WITH HEATMAPS
        if subPlotType != 'heatmap':
            if plotType != '1d':
                if plotOptions['X']['numeric']:
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
            pass
        col_wrap_min = 6
        #Make sure there are not more than 6 plots in a row (if no row facet variable)
        if 'row' not in facetgridkwargs.keys():
            if 'col' in facetgridkwargs.keys():
                colwrap = min([len(kwargs['col_order']),col_wrap_min])
                kwargs['col_wrap'] = colwrap
        
        fullTitleString = createFacetPlotName(folderName,dataType,plotType,subPlotType,legendParameterToLevelNameDict,subsettedDfTitle,levelsPlottedIndividually,useModifiedDf,plotOptions)
        if len(subsettedDf.index) > 0:
            plotSubsettedFigure(subsettedDf,plottingDf,kwargs,facetgridkwargs,plotType,subPlotType,dataType,fullTitleString,plotOptions,subsettedDfTitle,addDistributionPoints,alternateTitle)

#Will not be needed when seaborn 0.9.1 releases
def sanitizeSameValueLevels(plottingDf,kwargs):
    #First find levels that have values that are the same as other levels:
    #Get numeric axes levels
    numericAxes = []
    for axis in ['x','y','z']:
        if axis in kwargs.keys():
            numericAxes.append(kwargs[axis])
    #Exclude levels used for numeric axes (x,y,color)
    nonNumericNameList = []
    for levelName in plottingDf.columns:
        if levelName not in numericAxes:
            nonNumericNameList.append(levelName)
    #Get all possible pairs of nonnumericlevels
    possiblePairs = list(itertools.combinations(nonNumericNameList,2))
    overlappingPairs = []
    for possiblePair in possiblePairs:
        uniqueLevelValues1 = pd.unique(plottingDf[possiblePair[0]])
        uniqueLevelValues2 = pd.unique(plottingDf[possiblePair[1]])
        #If overlap, add both levelnames to list
        if bool(set(uniqueLevelValues1) & set(uniqueLevelValues2)):
            overlappingPairs.append(possiblePair)
    #Add dummy spaces to make level values different while not affecting the legend
    #NEED TO THINK ABOUT HOW THIS AFFECTS NUMERIC VARIABLES
    #NEED TO THINK ABOUT HOW TO DEAL WITH VARIABLES I NEED TO SORT NUMERICALLY; MOST LIKELY BY DETERMINING ORDER BEFORE SANITIZING AND ADDING SPACES TO ORDER
    for overlappingPair in overlappingPairs:
        newVals = []
        for value in plottingDf[overlappingPair[1]]:
            newval = value + ' '
            newVals.append(newval)
        plottingDf[overlappingPair[1]] = newVals
    return plottingDf

def plotSubsettedFigure(subsettedDf,plottingDf,kwargs,facetgridkwargs,plotType,subPlotType,dataType,fullTitleString,plotOptions,subsettedDfTitle,addDistributionPoints,alternateTitle):

    titleBool = True
    secondPathBool = False
    
    #Sanitize plottingDf to have every level value be unique across levels
    if subPlotType != 'heatmap':
        plottingDf = sanitizeSameValueLevels(plottingDf,kwargs)
    
    #Add in sharex/sharey options
    if plotType != '1d':
        facetKwargs = {'sharex':plotOptions['X']['share'],'sharey':plotOptions['Y']['share']}
    else:
        facetKwargs = {'sharex':False,'sharey':plotOptions['Y']['share']}

    auxillaryKwargs = {}
    auxillaryKwargs['plotType'] = plotType
    auxillaryKwargs['subPlotType'] = subPlotType
    auxillaryKwargs['facetgridkwargs'] = facetgridkwargs
    #Use appropriate facet plot
    #1D plots: (Histograms and KDEs); use facetgrid with appropriate axis level seaborn function
    if plotType == '1d':
        auxillaryKwargs['dataType'] = dataType
        facetPlotType = fp1D
    #1.5D/categorical plots (bar/point/box etc.); use catplot figure level seaborn function
    elif plotType == 'categorical':
        auxillaryKwargs['addDistributionPoints'] = addDistributionPoints
        facetPlotType = fpCategorical
    #2D plots (line and scatter); use relplot figure level seaborn function
    elif plotType == '2d':
        facetPlotType = fp2D
    #3D plots (heatmaps and 3d scatter/lineplots); use facet grid with appropriate axis level seaborn function
    elif plotType == '3d':
        facetPlotType = fp3D
    fg = facetPlotType.plot(plottingDf,subsettedDf,kwargs,facetKwargs,auxillaryKwargs,plotOptions)
    #SupTitle
    #Make room for suptitle
    if 'row' in kwargs:
        if len(kwargs['row_order']) <=2:
            plt.subplots_adjust(top=0.9)
        else:
            plt.subplots_adjust(top=0.9+len(kwargs['row_order'])*0.005)
    else:
        plt.subplots_adjust(top=0.9)
    if subsettedDfTitle != 'All':
        subsettedDfTitle = list(map(str,subsettedDfTitle))
    else:
        subsettedDfTitle = [subsettedDfTitle]
    #Do not include placeholder celltypes in suptitle
    if 'NotApplicable' in subsettedDfTitle:
        subsettedDfTitle.remove('NotApplicable')
    if 'row' not in kwargs and 'col' not in kwargs:
        plt.suptitle('-'.join(subsettedDfTitle),fontsize = 'x-large',fontweight='bold')
    else:
        if titleBool:
            plt.suptitle('-'.join(subsettedDfTitle),fontsize = 'x-large',fontweight='bold')

    #Save figure
    if dataType == 'dr':
        fullTitleString = alternateTitle
    fg.fig.savefig('fullyProcessedFigures/'+fullTitleString+'.png',bbox_inches='tight')
    if secondPathBool:
        fg.fig.savefig('../../outputFigures/'+fullTitleString+'.png',bbox_inches='tight')
    print(fullTitleString+' plot saved')
    plt.clf()

#!/usr/bin/env python3
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys,os,json,pickle,math,itertools
from facetPlotHeatmaps import draw_faceted_heatmap,returnHeatmapAspectRatios
from miscFunctions import sortSINumerically,reindexDataFrame
import subprocess
from matplotlib.colors import LogNorm,SymLogNorm
sys.path.insert(0, '../dataProcessing/')
from miscFunctions import returnGates,returnTicks

def produceSubsettedDataFrames(df,withinFigureBoolean,specificValueBooleanList):

    #sns.set_palette('muted')

    #Load in dataframe for experiment
    if 'Dimension' not in df.columns[0]:
        fulldf = df.stack()
    else:
        fulldf = df.copy()
    print(fulldf)
    print(withinFigureBoolean)
    print(specificValueBooleanList)
    print(fulldf)
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
            if 'Dimension' not in df.columns[0]:
                copiedDf = fulldf.loc[tuple(subsettingTuple)].copy()
            else:
                copiedDf = fulldf.loc[tuple(subsettingTuple),:].copy()
        except:
            print('key error: ')
            print(subsettingTuple)
        else:
            subsettedvalues = copiedDf.values
            subsettedindex = copiedDf.index.remove_unused_levels()
            if 'Dimension' not in df.columns[0]:
                subsetteddf = pd.DataFrame(subsettedvalues,index=subsettedindex)
            else:
                subsetteddf = pd.DataFrame(subsettedvalues,index=subsettedindex,columns=df.columns)
            subsettedDfList.append(subsetteddf)
            subsettedDfListTitles.append(subsettingTitle)
    return subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually

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
    #legendParameterToLevelNameDict = pickle.load(open('semiProcessedData/gui-radioVals.pkl','rb'))
    #plotOptions = pickle.load(open('semiProcessedData/gui-plotOptions.pkl','rb'))
    #legendParameterToLevelNameDict = 

    #Has issues with repeated values (aka CD54 shows up in TCells and APCs)
    for subsettedDf,subsettedDfTitle in zip(subsettedDfList,subsettedDfListTitles):
        print(subsettedDf)
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
        if subPlotType not in ['kde','histogram']:
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
        else:
            pass

        if dataType == 'singlecell':
            plottingDf = subsettedDf.stack().to_frame('GFI')
        else:
            plottingDf = subsettedDf.copy()
        print(plottingDf)
        #Converts wide form dataframe into long form required for cat/relplot
        plottingDf = plottingDf.reset_index()
        print(plottingDf)
        #Use plot options file to initialize numeric x axis ordering
        #NEED TO GET WORKING WITH HEATMAPS
        if subPlotType != 'heatmap':
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

    #Bandwidth for kde plots; changes per machine for some reason so need to be modified often
    #bandwidth = 1.5
    bandwidth = 0.15
    #shareYVar = True
    titleBool = True
    secondPathBool = False
    
    #Sanitize plottingDf to have every level value be unique across levels
    plottingDf = sanitizeSameValueLevels(plottingDf,kwargs)

    #Use appropriate facet plot
    #1d plots: #Histograms,KDEs
    if plotType == '1d':
        #Will need to update to make sure it pulls from y axis variable
        if subPlotType == 'kde':
            fg = sns.FacetGrid(plottingDf,sharey=False,legend_out=True,**kwargs)
            fg.map(sns.kdeplot,'GFI',shade=False,bw=bandwidth)
        elif subPlotType == 'histogram':
            fg = sns.FacetGrid(plottingDf,sharey=False,legend_out=True,**kwargs)
            fg.map(sns.distplot,'GFI',bins=256,kde=False)
        if dataType == 'singlecell':
            #Get GFI xtick values
            xtickValues,xtickLabels = returnTicks([-1000,1000,10000,100000])
            if subPlotType == 'kde':
                #Get count ytick values from histograms
                g = sns.FacetGrid(plottingDf,sharey=False,legend_out=True,**kwargs)
                g.map(sns.distplot,'GFI',bins=256,kde=False)
                ylabels = []
                for axis in plt.gcf().axes:
                    ylabels.append(list(map(int,axis.get_yticks().tolist())))
                plt.clf()
            #Add appropriate xtick values (also ytick values if kde) for each axis in figure
            for axis,i in zip(fg.fig.get_axes(),range(len(fg.fig.get_axes()))):
                axis.set_xticks(xtickValues)
                axis.set_xticklabels(xtickLabels)
                if subPlotType == 'kde':
                    axis.set_yticklabels(ylabels[i])
        fg.add_legend()
    #2d plots: categorical plots (bar, point), ordered plots (line,scatter)
    elif plotType in ['categorical','2d']:
        if plotOptions['Y']['axisScaling'] != 'Linear':
            errorBar = 95
        else:
            errorBar = 'sd'
        if plotType == 'categorical':
            if subPlotType == 'point':
                if addDistributionPoints:
                    fg = sns.catplot(**kwargs,data=plottingDf,kind=subPlotType,ci=errorBar,join=False,color='k',capsize=0.05,markers='_',zorder=3,errwidth=1)
                else:
                    fg = sns.catplot(**kwargs,data=plottingDf,kind=subPlotType,ci=errorBar,join=False,capsize=0.05,errwidth=1)
            elif subPlotType == 'box':
                fg = sns.catplot(**kwargs,data=plottingDf,kind=subPlotType)
            elif subPlotType == 'bar':
                fg = sns.catplot(**kwargs,data=plottingDf,kind=subPlotType,ci=errorBar,alpha=0.8,errwidth=1,capsize=0.05)
            #violin
            else:
                fg = sns.catplot(**kwargs,data=plottingDf,kind=subPlotType,alpha=0.8,scale='width')
            NoneType = type(None)
            if addDistributionPoints:
                secondkwargs = kwargs.copy()
                for key in ['row','col','col_order','row_order','col_wrap']:
                    if key in secondkwargs.keys():
                        secondkwargs.pop(key,None)
                if subPlotType != 'violin':
                    secondkwargs['dodge'] = True
                secondkwargs['edgecolor'] = 'black'
                secondkwargs['linewidth'] = 0.3
                secondkwargs['zorder'] = 1
                axisIndex  = 0
                if 'row' in kwargs and 'col' in kwargs:
                    for rowVal in pd.unique(plottingDf[kwargs['row']]):
                        for colVal in pd.unique(plottingDf[kwargs['col']]):
                            secondPlottingDf = plottingDf[plottingDf[kwargs['row']] == rowVal]
                            secondPlottingDf = secondPlottingDf[secondPlottingDf[kwargs['col']] == colVal]
                            if subPlotType != 'violin':
                                a = sns.stripplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                            else:
                                a = sns.swarmplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                            if rowVal != pd.unique(plottingDf[kwargs['row']])[-1]:
                                fg.fig.axes[axisIndex].set_xlabel('')
                            if not isinstance(a.legend_, NoneType):
                                a.legend_.remove()
                            axisIndex+=1
                else:
                    if 'row' in kwargs:
                        for rowVal in pd.unique(plottingDf[kwargs['row']]):
                            secondPlottingDf = plottingDf[plottingDf[kwargs['row']] == rowVal]
                            if subPlotType != 'violin':
                                a = sns.stripplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                            else:
                                a = sns.swarmplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                            if rowVal != pd.unique(plottingDf[kwargs['row']])[-1]:
                                fg.fig.axes[axisIndex].set_xlabel('')
                            axisIndex+=1
                    elif 'col' in kwargs:
                        for colVal in pd.unique(plottingDf[kwargs['col']]):
                            secondPlottingDf = plottingDf[plottingDf[kwargs['col']] == colVal]
                            if subPlotType != 'violin':
                                a = sns.stripplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                            else:
                                a = sns.swarmplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                            axisIndex+=1
                    else:
                        if subPlotType != 'violin':
                            a = sns.stripplot(**secondkwargs,data=plottingDf,ax=fg.fig.axes[axisIndex])
                        else:
                            a = sns.swarmplot(**secondkwargs,data=plottingDf,ax=fg.fig.axes[axisIndex])
                    if not isinstance(a.legend_, NoneType):
                        a.legend_.remove()
                for ax in fg.axes.flat:
                    ax.set_ylim(min(plottingDf[kwargs['y']])*0.5,max(plottingDf[kwargs['y']])*2)
        else:
            if dataType != 'dr':
                shareXVar = 'col'
                shareYVar = 'row'
            else:
                shareXVar = False
                shareYVar = False
            #Make sure there are markers at each column variable
            if 'style' not in kwargs.keys():
                fg = sns.relplot(data=plottingDf,marker='o',kind=subPlotType,facet_kws={'sharex':shareXVar,'sharey':shareYVar},ci=False,**kwargs)
            else:
                fg = sns.relplot(data=plottingDf,markers=True,kind=subPlotType,facet_kws={'sharex':shareXVar,'sharey':shareYVar},ci=False,**kwargs)
        #X and Y Axis Scaling for 2D plots
        for axis in plotOptions:
            k = len(fg.fig.get_axes())
            if 'Y' in axis:
                if plotOptions[axis]['axisScaling'] == 'Logarithmic':
                    for i in range(k):
                        fg.fig.get_axes()[i].set_yscale('log')
                elif plotOptions[axis]['axisScaling'] == 'Biexponential':
                    for i in range(k):
                        fg.fig.get_axes()[0].set_yscale('symlog',linthreshx=plotOptions[axis]['linThreshold'])
            else:
                if plotOptions[axis]['axisScaling'] == 'Logarithmic':
                    for i in range(k):
                        fg.fig.get_axes()[0].set_xscale('log')
                elif plotOptions[axis]['axisScaling'] == 'Biexponential':
                    for i in range(k):
                        fg.fig.get_axes()[0].set_xscale('symlog',linthreshx=plotOptions[axis]['linThreshold'])
    #3d plots: heatmaps and 3d scatter/lineplots
    elif plotType == '3d':
        if subPlotType == 'heatmap':
            a,h = returnHeatmapAspectRatios(subsettedDf,kwargs)
            fg = sns.FacetGrid(plottingDf,height=h,aspect=a,gridspec_kws={"wspace":0.4},**facetgridkwargs)
            if plotOptions['Colorbar']['axisScaling'] != 'Linear' and subsettedDf.values.min() > 0:
                logbool = True
                symlogbool = False
                log_norm = LogNorm(vmin=subsettedDf.min().min(), vmax=subsettedDf.max().max())
                sym_log_norm = ''
                cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(subsettedDf.min().min())), 1+math.ceil(math.log10(subsettedDf.max().max())))]
                lin_thresh = ''
            elif plotOptions['Colorbar']['axisScaling'] != 'Linear' and subsettedDf.values.min() <= 0:
                logbool = False 
                symlogbool = True
                lin_thresh =plotOptions['Colorbar']['linThreshold']
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
    if dataType == 'dr':
        print(alternateTitle)
        fullTitleString = alternateTitle
    fg.fig.savefig('fullyProcessedFigures/'+fullTitleString+'.png',bbox_inches='tight')
    if secondPathBool:
        fg.fig.savefig('../../outputFigures/'+fullTitleString+'.png',bbox_inches='tight')
    print(fullTitleString+' plot saved')
    plt.clf()

#!/usr/bin/env python3
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import pickle,sys,os
from itertools import groupby
sys.path.insert(0, '../../programs/dataProcessing/')
from miscFunctions import reindexDataFrame
from matplotlib import colors,ticker

def plot(plottingDf,subsettedDf,kwargs,facetKwargs,auxillaryKwargs,plotOptions):
    if plotOptions['Y']['axisScaling'] != 'Linear':
        errorBar = 95
    else:
        errorBar = 'sd'
    print(facetKwargs)
    if auxillaryKwargs['subPlotType'] == 'point':
        if auxillaryKwargs['addDistributionPoints']:
            fg = sns.catplot(**kwargs,**facetKwargs,data=plottingDf,kind=auxillaryKwargs['subPlotType'],ci=errorBar,join=False,color='k',capsize=0.05,markers='_',zorder=3,errwidth=1)
        else:
            fg = sns.catplot(**kwargs,**facetKwargs,data=plottingDf,kind=auxillaryKwargs['subPlotType'],ci=errorBar,join=False,capsize=0.05,errwidth=1)
    elif auxillaryKwargs['subPlotType'] == 'box':
        fg = sns.catplot(**kwargs,**facetKwargs,data=plottingDf,kind=auxillaryKwargs['subPlotType'])
    elif auxillaryKwargs['subPlotType'] == 'bar':
        fg = sns.catplot(**kwargs,**facetKwargs,data=plottingDf,kind=auxillaryKwargs['subPlotType'],ci=errorBar,alpha=0.8,errwidth=1,capsize=0.05)
    #violin
    else:
        fg = sns.catplot(**kwargs,facet_kws=facetKwargs,data=plottingDf,kind=auxillaryKwargs['subPlotType'],alpha=0.8,scale='width')
    NoneType = type(None)
    if auxillaryKwargs['addDistributionPoints']:
        secondkwargs = kwargs.copy()
        for key in ['row','col','col_order','row_order','col_wrap']:
            if key in secondkwargs.keys():
                secondkwargs.pop(key,None)
        if auxillaryKwargs['subPlotType'] != 'violin':
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
                    if auxillaryKwargs['subPlotType'] != 'violin':
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
                    if auxillaryKwargs['subPlotType'] != 'violin':
                        a = sns.stripplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                    else:
                        a = sns.swarmplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                    if rowVal != pd.unique(plottingDf[kwargs['row']])[-1]:
                        fg.fig.axes[axisIndex].set_xlabel('')
                    axisIndex+=1
            elif 'col' in kwargs:
                for colVal in pd.unique(plottingDf[kwargs['col']]):
                    secondPlottingDf = plottingDf[plottingDf[kwargs['col']] == colVal]
                    if auxillaryKwargs['subPlotType'] != 'violin':
                        a = sns.stripplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                    else:
                        a = sns.swarmplot(**secondkwargs,data=secondPlottingDf,ax=fg.fig.axes[axisIndex])
                    axisIndex+=1
            else:
                if auxillaryKwargs['subPlotType'] != 'violin':
                    a = sns.stripplot(**secondkwargs,data=plottingDf,ax=fg.fig.axes[axisIndex])
                else:
                    a = sns.swarmplot(**secondkwargs,data=plottingDf,ax=fg.fig.axes[axisIndex])
            if not isinstance(a.legend_, NoneType):
                a.legend_.remove()
        for ax in fg.axes.flat:
            ax.set_ylim(min(plottingDf[kwargs['y']])*0.5,max(plottingDf[kwargs['y']])*2)
    
    #X and Y Axis Scaling for 2D plots
    for axis in plotOptions:
        k = len(fg.fig.get_axes())
        if 'Y' in axis:
            if plotOptions[axis]['axisScaling'] == 'Logarithmic':
                for i in range(k):
                    fg.fig.get_axes()[i].set_yscale('log')
            elif plotOptions[axis]['axisScaling'] == 'Biexponential':
                for i in range(k):
                    fg.fig.get_axes()[i].set_yscale('symlog',linthreshx=plotOptions[axis]['linThreshold'])
            
            if str(plotOptions[axis]['limit'][0]) != '' or str(plotOptions[axis]['limit'][1]) != '':
                for i in range(k):
                    if str(plotOptions[axis]['limit'][0]) != '' and str(plotOptions[axis]['limit'][1]) != '':
                        fg.fig.get_axes()[i].set_ylim(bottom=float(plotOptions[axis]['limit'][0]),top=float(plotOptions[axis]['limit'][1]))
                    else:
                        if str(plotOptions[axis]['limit'][0]) != '':
                            fg.fig.get_axes()[i].set_ylim(bottom=float(plotOptions[axis]['limit'][0]))
                        else:
                            fg.fig.get_axes()[i].set_ylim(top=float(plotOptions[axis]['limit'][1]))
        else:
            if plotOptions[axis]['axisScaling'] == 'Logarithmic':
                for i in range(k):
                    fg.fig.get_axes()[i].set_xscale('log')
            elif plotOptions[axis]['axisScaling'] == 'Biexponential':
                for i in range(k):
                    fg.fig.get_axes()[i].set_xscale('symlog',linthreshx=plotOptions[axis]['linThreshold']) 
            
            if str(plotOptions[axis]['limit'][0]) != '' or str(plotOptions[axis]['limit'][1]) != '':
                for i in range(k):
                    if str(plotOptions[axis]['limit'][0]) != '' and str(plotOptions[axis]['limit'][1]) != '':
                        fg.fig.get_axes()[i].set_xlim(bottom=float(plotOptions[axis]['limit'][0]),top=float(plotOptions[axis]['limit'][1]))
                    else:
                        if str(plotOptions[axis]['limit'][0]) != '':
                            fg.fig.get_axes()[i].set_xlim(bottom=float(plotOptions[axis]['limit'][0]))
                        else:
                            fg.fig.get_axes()[i].set_xlim(top=float(plotOptions[axis]['limit'][1]))
    return fg

#!/usr/bin/env python3
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import pickle,sys,os
from itertools import groupby
from matplotlib.widgets import RadioButtons,Button,CheckButtons,TextBox
print(os.getcwd())
sys.path.insert(0, '../../programs/dataProcessing/')
print(os.getcwd())
from miscFunctions import reindexDataFrame
from matplotlib import colors,ticker
        
def plot(plottingDf,subsettedDf,kwargs,facetKwargs,auxillaryKwargs,plotOptions):
    #Will need to update to make sure it pulls from y axis variable
    yvar = kwargs.pop('y')
    if auxillaryKwargs['subPlotType'] == 'kde':
        fg = sns.FacetGrid(plottingDf,legend_out=True,**facetKwargs,**kwargs)
        fg.map(sns.kdeplot,yvar,shade=False)
    elif auxillaryKwargs['subPlotType'] == 'histogram':
        fg = sns.FacetGrid(plottingDf,legend_out=True,**facetKwargs,**kwargs)
        fg.map(sns.distplot,yvar,bins=256,kde=False)
    if auxillaryKwargs['dataType'] == 'singlecell':
        #Get GFI xtick values
        xtickValues,xtickLabels = returnTicks([-1000,1000,10000,100000])
        if auxillaryKwargs['subPlotType'] == 'kde':
            #Get count ytick values from histograms
            g = sns.FacetGrid(plottingDf,legend_out=True,**facetKwargs,**kwargs)
            g.map(sns.distplot,'GFI',bins=256,kde=False)
            ylabels = []
            for axis in plt.gcf().axes:
                ylabels.append(list(map(int,axis.get_yticks().tolist())))
            plt.clf()
        #Add appropriate xtick values (also ytick values if kde) for each axis in figure
        for axis,i in zip(fg.fig.get_axes(),range(len(fg.fig.get_axes()))):
            axis.set_xticks(xtickValues)
            axis.set_xticklabels(xtickLabels)
            if auxillaryKwargs['subPlotType'] == 'kde':
                axis.set_yticklabels(ylabels[i])
    fg.add_legend()
    return fg

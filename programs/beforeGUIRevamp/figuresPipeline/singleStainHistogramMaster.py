#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
import json,pickle,math,matplotlib,sys,os
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

def createHistograms(logicleDf,fileName):
    idx = pd.IndexSlice
    logicleDf.columns.name = 'Markers'
    #Remove FSC and SSC populations and convert dataframe into single column
    histogramDf = logicleDf.iloc[:,2:]
    """
    markersToPlot = ['H-2Kb','PD-L1']
    #[46, 44, 41, 38, 34, 30, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 12, 10, 8, 6, 4]
    timePointsToPlot = list(pd.unique(histogramDf.index.get_level_values('Time')))[::4]
    #timePointsToPlot = [46.0,4.0]
    #histogramDf = histogramDf.loc[idx[:,:],idx[markersToPlot]]
    #histogramDf = histogramDf.loc[idx[:,timePointsToPlot],idx[markersToPlot]]
    for marker in markersToPlot:
        subsettedDf = histogramDf.loc[:,marker]
        subsettedDf = subsettedDf.to_frame('GFI')
        plottingDf = subsettedDf.reset_index()
        g = sns.FacetGrid(plottingDf,hue='IFNgPulseConcentration',col='Time',sharey=False,aspect=2,legend_out=True,col_wrap=6)
        #Higher bandwidth; more smooth, lower, less smooth
        g.map(sns.kdeplot,'GFI',shade=True,bw=15)
        g.add_legend()
        plt.savefig('fullyProcessedFigures/histograms-'+fileName+'-'+marker+'.png')
        
    """
    markerSpecificList = []
    #Select only the events in the channel the single stain is testing
    print(histogramDf)
    histogramDf = histogramDf.xs(['Single Cells'],level=['CellType'])
    for marker in pd.unique(histogramDf.index.get_level_values('Antibody')):
        markerSpecificEvents = histogramDf.loc[idx[marker,:,:,:],marker]
        markerSpecificList.append(markerSpecificEvents)

    #Concatenate dataframes containing only eents in the channel tested to a single large dataframe with column name GFI
    markerSpecificDf = pd.concat(markerSpecificList,keys=pd.unique(histogramDf.index.get_level_values('Antibody')),names=logicleDf.index.names)
    markerSpecificDf = markerSpecificDf.to_frame('GFI')
    #Convert wideform into long form datafram required for facetgrid 
    plottingDf = markerSpecificDf.reset_index()

    #Set up a grid of subplot locations to put the single stain histograms in. 
    g = sns.FacetGrid(plottingDf,row='Dilution',col='Antibody',hue='Antibody',sharey=False)
    #Overlay kde and histogram on top of each other
    #Higher bandwidth; more smooth, lower, less smooth
    g.map(sns.kdeplot,'GFI',shade=True,bw=15)
    g.map(sns.distplot,'GFI',bins=256,kde_kws={"color": "k", "alpha": 0})

    #Save figure
    plt.savefig('fullyProcessedFigures/antibodyTitration-'+fileName+'.png')
    plt.close()

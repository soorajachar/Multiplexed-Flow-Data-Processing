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
    histogramDf = logicleDf.iloc[:,2:].stack()
    markerSpecificList = []
    #Select only the events in the channel the single stain is testing
    for marker in pd.unique(histogramDf.index.get_level_values('Antibody')):
        markerSpecificEvents = histogramDf.loc[idx[marker,:,:,:,marker]]
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

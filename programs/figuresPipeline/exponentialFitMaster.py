#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on fri sep 7 13:12:56 2018

@author: acharsr
"""
import matplotlib
import numpy as np
import sys
import pickle
import os
import pandas as pd
from matplotlib import pyplot as plt
from facetPlottingLibrary import createParameterValues,buildLegendHandles,createSubPlotGrid,returnMarkerVals
from miscFunctions import returnNumericTimePoints
from scipy.optimize import curve_fit
from miscFunctions import r_squared,boundedExponential,logisticDoubleExponential

lineplot = False
idx = pd.IndexSlice

def findFits(x,logy,cytokineMin,fitName):
    epsilon = 0.0001
    if fitName == 'boundedExponential':
        lbounds = [(np.max(logy)-cytokineMin)/2,0,cytokineMin-epsilon]
        ubounds = [(np.max(logy)-cytokineMin)*2,10,cytokineMin+epsilon]
        try:
            popt,pcov = curve_fit(boundedExponential,x,logy,sigma=logy,bounds=(lbounds,ubounds))
            rsquared = round(r_squared(x,logy,boundedExponential,popt),3)
        except:
            return [0,0,cytokineMin],0.0
    else:
        lbounds = [(np.max(logy)-cytokineMin)-epsilon,0,0,0,x[logy.tolist().index(np.max(logy))],cytokineMin-epsilon]
        ubounds = [(np.max(logy)-cytokineMin)+epsilon,2,2,x[logy.tolist().index(np.max(logy))],np.max(x),cytokineMin+epsilon]
        try:
            popt,pcov = curve_fit(logisticDoubleExponential,x,logy,sigma=logy,bounds=(lbounds,ubounds))
            rsquared = round(r_squared(x,logy,logisticDoubleExponential,popt),3)
        except:
            return [0,0,0,0,0,cytokineMin],0.0
    return popt,rsquared

def createParameterDataFrame(secondPath,experimentNumber,df):
    x = returnNumericTimePoints(df)
    parameterizedObservables = ['IFNg','TNFa','IL-6','IL-2','IL-17A']
    parameterizedFits = ['boundedExponential','logisticDoubleExponential']
    parameters = [['r2','A','k0','v'],['r2','A','k0','k1','x0','x1','v']]
    dfListOfTuples = []
    for i in range(len(parameterizedFits)):
        fitName = parameterizedFits[i]
        parameterList = parameters[i]
        for parameter in parameterList:
            dfListOfTuples.append((fitName,parameter))
    columnIndex = pd.MultiIndex.from_tuples(dfListOfTuples,names=['Fit','Parameter'])
    emptyParameterMatrix = np.zeros((df.index.size,columnIndex.size))
    parameterDataFrame = pd.DataFrame(emptyParameterMatrix,index=df.index,columns=columnIndex)
    for observable in parameterizedObservables:
        observableDf = df.loc[observable]
        cytokineMin = np.log10(np.min(observableDf.values))
        for row in range(observableDf.shape[0]):
            logy = np.log10(observableDf.iloc[row,:])
            if observable in ['IFNg','TNFa','IL-6']:
                popt,rsquared = findFits(x,logy,cytokineMin,'boundedExponential')
                parameterDataFrame.loc[observable].iloc[row,parameterDataFrame.columns.get_level_values(0)=='boundedExponential'] = [rsquared,*popt]
            else:
                popt,rsquared = findFits(x,logy,cytokineMin,'logisticDoubleExponential')
                parameterDataFrame.loc[observable].iloc[row,parameterDataFrame.columns.get_level_values(0)=='logisticDoubleExponential'] = [rsquared,*popt]
    print(parameterDataFrame)
    with open(secondPath+'fitParameterPickleFiles/fitParameterPickleFile-%d.pkl'%(experimentNumber), "wb") as f:
        pickle.dump(parameterDataFrame, f)

def extractTrainingSetData(df,experimentNumber):
    if(experimentNumber == 68):
        newdf = df.loc[idx[:,'None',:,:],idx[:]]
        newdf.index = newdf.index.droplevel('Antibody')
        return newdf
    elif(experimentNumber == 53):
        newdf = df.loc[idx[:,'OT1',:,:],idx[:]]
        newdf.index = newdf.index.droplevel('TCellType')
        return newdf
    elif(experimentNumber in [69,70]):
        newdf = df.loc[idx[:,'WT',:,:],idx[:]]
        newdf.index = newdf.index.droplevel('Genotype')
        return newdf
    else:
        return df

def updateTrainingSet(secondPath,experimentNumbers):
    trainingSetList = []
    for expNum in experimentNumbers:
        df = pickle.load(open(secondPath+'fitParameterPickleFiles/fitParameterPickleFile-%d.pkl'%(expNum),'rb'))
        trainingSetList.append(extractTrainingSetData(df,expNum))
    fullTrainingSet = pd.concat(trainingSetList,keys=experimentNumbers,names=['DataSet'])
    with open(secondPath+'fitParameterPickleFiles/allParameterDataSets.pkl', "wb") as f:
        pickle.dump(fullTrainingSet, f)
    print(fullTrainingSet)

#Fill out subplots of the dataframe row by row with appropriately crossectioned data sets
def createFacetedPlot(dfc,legendParameterToLevelNameDict,nonLegendParameterToLevelNameDict):
    
    parameterToParameterValsDict, subPlotGridDimensions, subPlotGridLabels,subPlotGridLevels, maxNumConditions = createParameterValues(dfc,legendParameterToLevelNameDict, nonLegendParameterToLevelNameDict,lineplot,True)
    legendHandles = buildLegendHandles(dfc,parameterToParameterValsDict,legendParameterToLevelNameDict,nonLegendParameterToLevelNameDict,maxNumConditions,lineplot)

    fig1 = plt.figure(figsize=(20*subPlotGridDimensions[1],20*subPlotGridDimensions[0]),facecolor='w',edgecolor='k')
    plt.gcf().set_facecolor('white')
    #LODParameters = pickle.load(open('semiProcessedData/LODParameters-20181018-PeptideComparison_OT1_Timeseries_9-nM.pkl', "rb"))

    fig1,axes,ax2,subtitles,subplotdfs = createSubPlotGrid(fig1,dfc,subPlotGridDimensions,subPlotGridLabels,subPlotGridLevels)
    x = returnNumericTimePoints(dfc)
    LODFileLocation = 'semiProcessedData/LODParameters-20181018-PeptideComparison_OT1_Timeseries_9-nM.pkl'
    LODParameters = pickle.load(open(LODFileLocation, "rb"))
    
    for ax,subtitle,subplotdf in zip(axes,subtitles,subplotdfs):
        ax.set_title(subtitle)
        cytokineMin = np.log10(np.min(subplotdf.values))
        #Plot dataframe row by row (to allow timepoints to be reused as x values)
        for row in range(0,subplotdf.shape[0]):
            y = subplotdf.iloc[row,:]
            logy = np.log10(y)
            markerVal,colorVal,alphaVal,sizeVal = returnMarkerVals(subplotdf,row,parameterToParameterValsDict,False)
            ax.scatter(x,np.log10(y),s=sizeVal,marker=markerVal,alpha=alphaVal,color=colorVal)
            try:
                #For boundedExponential
                if(subtitle in ['IFNg','TNFa','IL-6']): #or ('N4' in list(y.name) and '1uM' in list(y.name))):
                    popt,rsquared = findFits(x,logy,cytokineMin,'boundedExponential')
                    if(rsquared > 0.3):
                        curveFitPlotPoints = np.linspace(0.01,np.max(x),201)
                        ax.plot(curveFitPlotPoints,boundedExponential(curveFitPlotPoints,*popt),marker=None,color=colorVal)
                #For logisticDoubleExponential
                elif(subtitle in ['IL-2','IL-17A']):
                    popt,rsquared = findFits(x,logy,cytokineMin,'logisticDoubleExponential')
                    if(rsquared > 0.3):
                        curveFitPlotPoints = np.linspace(0.01,np.max(x),201)
                        ax.plot(curveFitPlotPoints,logisticDoubleExponential(curveFitPlotPoints,*popt),marker=None,color=colorVal)
                #Do not plot expfit
                else:
                    pass
            except Exception as e:
                print(e)
    
    #Add legend at location slightly to the left of the first subplot
    fig1.legend(handles = legendHandles,ncol=len(parameterToParameterValsDict),bbox_to_anchor=(-0.05,1),bbox_transform=ax2.transAxes)
    return fig1

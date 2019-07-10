#!/usr/bin/env python3

# coding: utf-8

# # Crossvalidate autoencoder
# Import dependencies
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from facetPlottingLibrary import createParameterValues,buildLegendHandles,returnMarkerVals,loadInDataFrameNoUI
from mpl_toolkits.mplot3d import Axes3D

lineplot = False
        
def projectNeurons2D(secondPath,crossValidatedDf,modelName,trainingString,crossValidateNumber):
        
    legendParameterToLevelNameDict = loadInDataFrameNoUI(crossValidatedDf) 
    nonlegendParameterToLevelNameDict = {}
    
    parameterToParameterValsDict,_,_,_,maxNumConditions = createParameterValues(crossValidatedDf,legendParameterToLevelNameDict,nonlegendParameterToLevelNameDict,lineplot,True) 
    legendHandles = buildLegendHandles(crossValidatedDf,parameterToParameterValsDict,legendParameterToLevelNameDict,nonlegendParameterToLevelNameDict,maxNumConditions,lineplot)
    
    fig,axes = plt.subplots(2,2,figsize=(15,15))
    axes[1,1].axis('off')
    for ax,(i,j) in zip(axes.flatten(),[(0,1),(0,2),(1,2)]):
        if(i == 1 and j == 2):
            legendax = ax
        #Plot dataframe row by row (to allow timepoints to be reused as x values)
        for row in range(crossValidatedDf.shape[0]):
            markerVal,colorVal,alphaVal,sizeVal = returnMarkerVals(crossValidatedDf,row,parameterToParameterValsDict,lineplot)
            ax.scatter(crossValidatedDf.iloc[row,i],crossValidatedDf.iloc[row,j],s=sizeVal,marker=markerVal,alpha=alphaVal,color=colorVal)

        ax.set_xlabel(crossValidatedDf.columns.values[i])
        ax.set_ylabel(crossValidatedDf.columns.values[j])
        ax.set_xlim([-1,1])
        ax.set_ylim([-1,1])
    
    #Add legend at location slightly to the left of the first subplot
    fig.legend(handles = legendHandles,ncol=len(parameterToParameterValsDict),bbox_to_anchor=(2,1),bbox_transform=legendax.transAxes)
    plt.title('%d-on-%s'%(crossValidateNumber,trainingString),fontweight='bold')
    fig.savefig(secondPath+'projectionsOnNeurons/2d/%s/%s/projectionOnNeurons-%d-on-%s-%s.png'%(modelName,trainingString,crossValidateNumber,modelName,trainingString),bbox_inches='tight',dpi=120)
    plt.close(fig)

def projectNeurons3D(secondPath,crossValidatedDf,modelName,trainingString,crossValidateNumber):
    
    legendParameterToLevelNameDict = loadInDataFrameNoUI(crossValidatedDf) 
    nonlegendParameterToLevelNameDict = {}

    parameterToParameterValsDict,_,_,_,maxNumConditions = createParameterValues(crossValidatedDf,legendParameterToLevelNameDict,nonlegendParameterToLevelNameDict,lineplot,True) 
    legendHandles = buildLegendHandles(crossValidatedDf,parameterToParameterValsDict,legendParameterToLevelNameDict,nonlegendParameterToLevelNameDict,maxNumConditions,lineplot)
    
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel(crossValidatedDf.columns.values[0])
    ax.set_ylabel(crossValidatedDf.columns.values[1])
    ax.set_zlabel(crossValidatedDf.columns.values[2])
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim([-1,1])

    #Plot dataframe row by row (to allow timepoints to be reused as x values)
    for row in range(crossValidatedDf.shape[0]):
        markerVal,colorVal,alphaVal,sizeVal = returnMarkerVals(crossValidatedDf,row,parameterToParameterValsDict,lineplot)
        ax.scatter(crossValidatedDf.iloc[row,0],crossValidatedDf.iloc[row,1],crossValidatedDf.iloc[row,2],s=sizeVal,marker=markerVal,alpha=alphaVal,color=colorVal)
    
    fig.legend(handles = legendHandles,ncol=len(parameterToParameterValsDict),bbox_to_anchor=(-0.05,1),bbox_transform=ax.transAxes)
    plt.title('%d-on-%s'%(crossValidateNumber,trainingString),fontweight='bold')
    fig.savefig(secondPath+'projectionsOnNeurons/3d/%s/%s/projectionOnNeurons3d-%d-on-%s-%s.png'%(modelName,trainingString,crossValidateNumber,modelName,trainingString),bbox_inches='tight',dpi=120)
    plt.close(fig)

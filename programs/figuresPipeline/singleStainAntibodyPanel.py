#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 2 13:12:56 2018

@author: acharsr
"""

import math
import pickle
import numpy as np
import matplotlib
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import os
import sys
import fcsparser
sys.path.insert(0, '../dataprocessing/')
from logicle import logicle,quantile 

#Retrieve fcs data and parse into a dataframe, make histograms for each single stain based on panel excel file
def createHistograms(folderName):
    
    os.chdir('inputFiles')
    #Read in excel file to dataframe
    panelData = pd.read_csv(folderName+'-testingPanel.csv')
    #Grab all subfolders in directory (different single stain conditions)
    dir_list = next(os.walk('.'))[1]
    #Walk through each singleStainCondition
    for singleStainCondition in dir_list:
        print(os.getcwd())
        #Switch to final gated lymphocyte population folder
        os.chdir(singleStainCondition+'/lymphocytes')
        #grab all files in directory (not folders)
        allFiles = list(filter(os.path.isfile, os.listdir()))
        #Create empty list, fill with names of each file in order of collection (eg sample 001 goes in index 0, sample 005 goes in index 5 etc)
        fileNums = []
        for fileName in allFiles:
            if (fileName.find('.fcs')>-1):
                i = int(fileName[fileName.rfind('_')-3:fileName.rfind('_')]) - 1
                fileNums.append(i)
        firstFileList = [None]*(max(fileNums)+1)
        for fileName in allFiles:
            if (fileName.find('.fcs')>-1):
                i = int(fileName[fileName.rfind('_')-3:fileName.rfind('_')]) - 1
                firstFileList[i] = fileName
        fileList = []
        for fname in firstFileList:
            if(fname != None):
                fileList.append(fname)
        dataList = []
        #fscCutoff = 30000
        fscCutoff = 0
        for fileName,index in zip(fileList,range(len(fileList))):
            #NoCorresponding
            #FCSDetectorName
            print(index)
            meta,fcsDf = fcsparser.parse(fileName, reformat_meta=True)
            if(panelData.iloc[index,:]['FCSDetectorName'] != 'NoCorresponding'):
                cellsDf = fcsDf[(fcsDf['FSC-A'] > fscCutoff)]
                dataList.append([panelData.iloc[index,:]['FCSDetectorName'],cellsDf[panelData.iloc[index,:]['FCSDetectorName']]])
        maxNumElements = 0
        for data in dataList:
            numElements = len(data[1])
            if(numElements > maxNumElements):
                maxNumElements = numElements
        dataMatrix = np.empty((len(dataList),maxNumElements,))
        dataMatrix[:] = np.nan
        T = 262143
        d = 4.5
        m = d*np.log(10)
        voltageProblems = []
        wVals = []
        print('WAT')
        for i in range(dataMatrix.shape[0]):
            logicleData0 = dataList[i][1]
            logicleData = logicleData0[logicleData0<T]
            if(len(logicleData[logicleData<0]) > 0):
                r = quantile(logicleData[logicleData<0], 0.05)
                w = (m-np.log(T/abs(r)))/2
            else:
                w=0
            wVals.append(w)
            if(len(logicleData0[logicleData0>=T]) >= 100):
                voltageProblems.append(i+1)
            transformedX = logicle(logicleData,T,m,r)
            dataMatrix[i,:len(logicleData)] = transformedX
        newPanelDataTupleList = []
        for i in range(panelData.shape[0]):
            if(panelData.iloc[i,:]['FCSDetectorName'] != 'NoCorresponding'):
                newPanelDataTupleList.append(tuple(panelData.iloc[i,:].values))
        multiIndexedObject = pd.MultiIndex.from_tuples(newPanelDataTupleList, names=panelData.columns)
        fullSingleStainDf = pd.DataFrame(dataMatrix,index=multiIndexedObject)
        fullSingleStainDf.columns.name = 'Fluorescence'
        idx = pd.IndexSlice
        markers = pd.unique(fullSingleStainDf.index.get_level_values('Marker'))
        sns.set_palette(sns.color_palette("hls", len(markers)))
        markerColors = sns.color_palette()
        numMarkers = fullSingleStainDf.shape[0]
        fig1 = plt.figure(num=1,figsize=(40,30),dpi=120,facecolor='w',edgecolor='k')
        index = 1
        height = 6
        binNum = 256 #whats used in flowjo apparently
        for marker,cindex in zip(markers,range(len(markers))):
            color = markerColors[cindex]
            markerSingleStains = fullSingleStainDf.loc[idx[:,marker,:,:,:,:,:],idx[:]]
            for row in range(markerSingleStains.shape[0]):
                detectorSingleStain = markerSingleStains.iloc[row,:]
                names = list(detectorSingleStain.name)
                logVals = detectorSingleStain.values[~np.isnan(detectorSingleStain.values)]
                cloneTitle = list(markerSingleStains.index.get_level_values('Clone'))[row]
                numberTitle = list(markerSingleStains.index.get_level_values('Dilution'))[row]
                if('Live' in folderName):
                    efficacyTitle = list(markerSingleStains.index.get_level_values('Live efficacy'))[row]
                else:
                    efficacyTitle = list(markerSingleStains.index.get_level_values('F/P efficacy'))[row]
                titleString = marker+'-'+cloneTitle+'-'+str(numberTitle)+'-('+efficacyTitle+')'
                ax = fig1.add_subplot(height,math.ceil(fullSingleStainDf.index.size/height),index)
                if(marker == 'NoStain'):
                    titleString = 'Negative Control'
                    color = 'black'
                plt.hist(logVals,color=color,bins=binNum)
                if(index in voltageProblems):
                    plt.title(titleString,fontweight='bold')
                else:
                    plt.title(titleString)
                
                ax.set_xlabel(list(markerSingleStains.index.get_level_values('FCSDetectorName'))[row])
                #Set xaxis ticks to correct logicle range
                numLinLabels = 3
                numLogLabels = 5
                
                w = wVals[index-1]
                W = w/np.log(10)
                
                linearRange = np.linspace(0,2*w,num=numLinLabels)
                linearRangeActualVals = np.linspace(-1*10**W,10**W,num=numLinLabels)
                
                logarithmicRange = np.linspace(2*w,m,num=numLogLabels)[1:]
                logarithmicRangeActualVals = np.logspace(W,np.log10(T),num=numLogLabels)[1:]
                
                logicleRange = np.concatenate((linearRange,logarithmicRange),axis=0)
                logicleRangeActualVals = np.concatenate((linearRangeActualVals,logarithmicRangeActualVals),axis=0)
                
                plt.xticks(logicleRange)
                labels = [item.get_text() for item in ax.get_xticklabels()]
                for val,index2 in zip(logicleRange,range(len(logicleRange))): 
                    if(index2 < len(linearRange)):
                        labels[index2] = int(logicleRangeActualVals[index2])
                    else:
                        exponent = np.log10(logicleRangeActualVals[index2])
                        expString = '$10^{%.1f}$'%(exponent)
                        labels[index2] = expString
                ax.set_xticklabels(labels)

                index+=1
        
        plt.suptitle(singleStainCondition,fontsize=20)
        plt.tight_layout()
        plt.subplots_adjust(top=0.95)
        fig1.savefig('../../../fullyProcessedFigures/singleStains-%s-logicle.png'%(singleStainCondition))
        print(voltageProblems)
        plt.clf()
        os.chdir('../..')
os.chdir('..')
#folderName = str(sys.argv[1])
#createDataFrames(folderName)

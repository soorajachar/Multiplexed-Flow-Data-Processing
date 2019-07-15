#!/usr/bin/env python3  
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import os,sys,pickle,math,re,subprocess
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import string 
import tkinter as tk
import tkinter.font as tkfont

def r_squared(xdata,ydata,func,popt):
    residuals = ydata- func(xdata, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

def InverseHill(y,parameters):
    Amplitude=parameters[0]
    EC50=parameters[1]
    hill=parameters[2]
    Background=parameters[3]
    return np.power((np.power(10,y)-Background)/(Amplitude-np.power(10,y)),1/hill)*EC50

def Hill(x, Amplitude, EC50, hill,Background):
    return np.log10(Amplitude * np.power(x,hill)/(np.power(EC50,hill)+np.power(x,hill))+Background)

#2 parameter (vshift fixed per cytokine based on lower LOD of cytokine): y = A(1-e^(-tau*x))
def boundedExponential(x, amplitude,tau,vshift):
    return amplitude*(np.subtract(1,np.exp(np.multiply(-1,np.multiply(x,tau)))))+vshift

#5 parameter (vshift fixed per cytokine; based on lower LOD of cytokine): y = A((1/(1+e^(-tau1*(x-td1))))-(1/(1+e^(-tau2*(x-td2)))))
def logisticDoubleExponential(x,amplitude,tau1,tau2,timedelay1,timedelay2,vshift):
    return amplitude*np.subtract(np.divide(1,np.add(1,np.exp(np.multiply(-1*tau1, np.subtract(x,timedelay1))))),np.divide(1,np.add(1,np.exp(np.multiply(-1*tau2,np.subtract(x,timedelay2))))))+vshift

#Sample_A2_A02_002.fcs,
def cleanUpFlowjoCSV(fileArray,folderName,dataType):
    if dataType == 'singlecell':
        dataTypeForCSV = 'cell'
    elif dataType == 'cytcorr':
        dataTypeForCSV = 'cyt'
    else:
        dataTypeForCSV = dataType
    sortedData = []
    sortedFiles = []
    #Samples will be indexed based on well ID (A01, then A02 etc.)
    orderWellID = {}
    #Supports 384 well (16 rows, 24 columns)
    plateColumnList = list(range(1,25))
    plateRowList = list(string.ascii_uppercase)[:16]
    index = 1
    for plateRow in plateRowList:
        for plateColumn in plateColumnList:
            orderWellID[str(plateRow)+str(plateColumn)] = index
            index+=1
    orderWellID['Mean'] = len(orderWellID.keys())+1
    orderWellID['SD'] = len(orderWellID.keys())+2
    for name in fileArray:
        temp = pd.read_csv('semiProcessedData/'+str(name)+'_'+dataTypeForCSV+'.csv')
        temp2 = [] 
        for i in range(0,temp.shape[0]):
            fullfilename = 'semiProcessedData/singleCellData/'+name+'/'+temp.iloc[i,0][:temp.iloc[i,0].find('.')]
            if '_' in temp.iloc[i,0]:
                wellID = temp.iloc[i,0].split('_')[2]
            else:
                wellID = temp.iloc[i,0]
            temp.iloc[i,0] = orderWellID[wellID]
            temp2.append([str(temp.iloc[i,0]).zfill(3),fullfilename])
        
        temp2 = pd.DataFrame(np.matrix(temp2[:-2]),columns=['Unnamed: 0','fileName'])
        temp = temp.sort_values('Unnamed: 0')
        temp2 = temp2.sort_values('Unnamed: 0')
        sortedData.append(temp[:-2])
        sortedFiles.append(temp2)
    
    return sortedData,sortedFiles

#['1uM' '1nM' '100pM' '10pM' '10nM' '100nM']
unitPrefixDictionary = {'fM':1e-15,'pM':1e-12,'nM':1e-9,'uM':1e-6,'mM':1e-3,'M':1e0,'':0,'K':1000}
def sortSINumerically(listSI,sort,descending):
    numericList = []
    for unitString in listSI:
        splitString = re.split('(\d+)',unitString)
        numericList.append(float(splitString[1])*float(unitPrefixDictionary[splitString[2]]))
    originalNumericList = numericList.copy()
    if sort:
        numericList.sort(reverse=descending)
    numericIndices = []
    for elem in numericList:
        numericIndices.append(originalNumericList.index(elem))
    sortedListSI = []
    for elem in numericIndices:
        sortedListSI.append(listSI[elem])
    return sortedListSI,numericList

#Used to interpret experimentnumbers
def parseCommandLineNNString(inputString):
    if(',' in inputString):
        if('-' in inputString): # - and ,
            experimentNumbers = []
            experimentRanges = list(inputString.split(','))
            for experimentRangeString in experimentRanges:
                if('-' in experimentRangeString):
                    experimentNumberRange = list(map(int, experimentRangeString.split('-')))
                    tempExperimentNumbers = list(range(experimentNumberRange[0],experimentNumberRange[1]+1))
                    for eNum in tempExperimentNumbers:
                        experimentNumbers.append(eNum)
                else:
                    experimentNumbers.append(int(experimentRangeString))
        else: #just ,
            experimentNumbers = list(map(int, inputString.split(',')))
    else:
        if('-' in inputString): #just -
            experimentNumberRange = list(map(int, inputString.split('-')))
            experimentNumbers = list(range(experimentNumberRange[0],experimentNumberRange[1]+1))
        else: #just single experiment number
            experimentNumbers = int(inputString)
    if isinstance(experimentNumbers, int):
        return [experimentNumbers]
    else:
        return experimentNumbers

#used for nonlexographic sort reindexing
def reindexDataFrame(dfToReindex,indexdf,singlecellToNonSinglecell):
    if indexdf.index.names[0] in ['Cytokine','Statistic']:
        indexingDf = indexdf.loc[pd.unique(indexdf.index.get_level_values(0))[0]]
    elif indexdf.index.names[0]  in ['CellType']:
        indexingDf = indexdf.loc[pd.unique(indexdf.index.get_level_values(0))[0]].loc[pd.unique(indexdf.index.get_level_values(1))[0]].loc[pd.unique(indexdf.index.get_level_values(2))[0]]
    else:
        indexingDf = indexdf
    if singlecellToNonSinglecell:
        indexingDf = indexingDf.stack()
        indexingDf = indexingDf.to_frame('temp')

    idx = pd.IndexSlice
    reindexedDfMatrix = np.zeros(dfToReindex.shape)
    if not singlecellToNonSinglecell:
        for row in range(indexingDf.shape[0]):
            if isinstance(indexingDf.iloc[row].index.name, (list,)):
                indexingLevelNames = tuple(indexingDf.iloc[row].name)
            else:
                indexingLevelNames = indexingDf.iloc[row].name
            print(dfToReindex)
            dfToReindexValues = dfToReindex.loc[idx[indexingLevelNames],:]
            reindexedDfMatrix[row,:] = dfToReindexValues
        reindexedDf = pd.DataFrame(reindexedDfMatrix,index=indexingDf.index,columns=dfToReindex.columns)
    else:
        k = 0
        row = 0
        reindexedLevelNames = []
        while k < dfToReindex.shape[0]:
            levelNames = tuple(indexingDf.iloc[row].name)
            stackedLevelNames = tuple(list(levelNames)+[slice(None)])
            dfToReindexValues = dfToReindex.loc[idx[stackedLevelNames],:]
            stackedLength = dfToReindexValues.shape[0]
            reindexedDfMatrix[k:k+stackedLength,:] = dfToReindexValues
            if dfToReindex.index.names[-1] == 'Event':
                for eventVal in range(1,1+stackedLength):
                    reindexedLevelNames.append(list(levelNames)+[eventVal])
            else:
                reindexedLevelNames.append(list(levelNames))
            row+=1
            k+=stackedLength
        print(reindexedLevelNames)
        print(dfToReindex.index.names)
        reindexedMultiIndex = pd.MultiIndex.from_tuples(reindexedLevelNames,names=dfToReindex.index.names)
        reindexedDf = pd.DataFrame(reindexedDfMatrix,index=reindexedMultiIndex,columns=dfToReindex.columns)
    return reindexedDf

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

def returnTicks(xticksToUse):
    logxticks = [-1000,-100,-10,0,10,100,1000,10000,100000]
    logicleXTicks = [64, 212, 229, 231, 233, 251, 399, 684, 925]
    xtickValues = []
    xtickLabels = []
    
    for logxtick in xticksToUse:
        if(logxtick < 0):
            xtickLabels.append('$-10^'+str(int(np.log10(-1*logxtick)))+'$')
        elif(logxtick == 0):
            xtickLabels.append('0')
        else:
            xtickLabels.append('$10^'+str(int(np.log10(logxtick)))+'$')
    
    for tickval in xticksToUse:
        xtickValue = logicleXTicks[logxticks.index(tickval)]
        xtickValues.append(xtickValue)
    
    return xtickValues,xtickLabels

def returnGates(logicleData,rawData,generationZeroBoundary):
    parentGenerationPresent = True
    maxGenerationNumber = 6
    newtemplin = logicleData.values.ravel(order='F')
    newtempraw = rawData.values.ravel(order='F')
    generationGatesLinear = [newtemplin[find_nearest(newtempraw,generationZeroBoundary)[1]]]
    #Get CTV GFI means for each generation by dividing initial GFI (raw) by 2 for each division
    generationZeroGFI = 0.75*generationZeroBoundary
    generationMeansLog = [generationZeroGFI]
    for generation in range(maxGenerationNumber):
        generationMeansLog.append(generationMeansLog[generation]/2)
    #Create initial gates at values between each set of division means. Undivided generation and final generation have boundaries at right, left edges of plot respectively
    generationGatesLog = []
    for generation in range(maxGenerationNumber-1):
        generationGatesLog.append((generationMeansLog[generation]+generationMeansLog[generation+1])*0.5)
    for gateval in generationGatesLog:
        generationGatesLinear.append(newtemplin[find_nearest(newtempraw,gateval)[1]])
    return generationGatesLinear

def extractValues(currentLevelLayout,valueToRemove,equalityBoolean):
    if equalityBoolean:
        idx = np.argwhere(np.all(currentLevelLayout[..., :] == valueToRemove, axis=0))
    else:
        idx = np.argwhere(np.all(currentLevelLayout[..., :] != valueToRemove, axis=0))
    currentLevelLayoutExtractedColumns = np.delete(currentLevelLayout, idx, axis=1)
    if equalityBoolean:
        idx2 = np.argwhere(np.all(currentLevelLayoutExtractedColumns[:,...] == valueToRemove, axis=1))
    else:
        idx2 = np.argwhere(np.all(currentLevelLayoutExtractedColumns[:,...] != valueToRemove, axis=1))
    currentLevelLayoutExtracted = np.delete(currentLevelLayoutExtractedColumns,idx2,axis=0)
    return currentLevelLayoutExtracted

def exitEnabledSubprocessRun(command,script,inputVariable,hasInput):
    if hasInput:
        subprocess.run([command,script]+[inputVariable])
    else:
        subprocess.run([command,script])
    if pickle.load(open('inputFiles/gui-exitBoolean.pkl','rb')):
        exitBoolean = False 
        with open('inputFiles/gui-exitBoolean.pkl','wb') as f:
            pickle.dump(exitBoolean,f)
        sys.exit(0)

def setMaxWidth(stringList, element):
    f = tkfont.nametofont(element.cget("font"))
    zerowidth=f.measure("0")
    w=max([f.measure(i) for i in stringList])/zerowidth
    element.config(width=math.ceil(1.2*w))

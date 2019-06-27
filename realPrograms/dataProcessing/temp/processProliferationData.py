#!/usr/bin/env python3  
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import matplotlib,statistics
#matplotlib.use('QT4Agg') 
import os,sys,pickle,json,math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from draggableLine import draggable_lines
from modifyDataFrames import originalIndexOrder
from logicle import logicle,quantile 

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

def createProliferationData(folderName,logicleDataFull,rawDataFull):
    logicleDataUnstacked = logicleDataFull.unstack('Event')
    rawDataUnstacked = rawDataFull.unstack('Event')
    for i in range(logicleDataFull.index.nlevels-1):
        logicleDataUnstacked = logicleDataUnstacked.reindex(logicleDataFull.index.get_level_values(i).unique(),level=i)
        rawDataUnstacked = rawDataUnstacked.reindex(logicleDataFull.index.get_level_values(i).unique(),level=i)
    return logicleDataUnstacked,rawDataUnstacked

def returnGatesAndTicks(logicleData,rawData):
    maxGenerationNumber = 8
    generationZeroGFI = 100000
    logxticks = [-1000,100,1000,10000,100000]
    print(logicleData) 
    #Get CTV GFI means for each generation by dividing initial GFI by 2 for each division
    generationMeansLog = [generationZeroGFI]
    for generation in range(maxGenerationNumber):
        generationMeansLog.append(generationMeansLog[generation]/2)
    #Create initial gates at values between each set of division means. Undivided generation and final generation have boundaries at right, left edges of plot respectively
    generationGatesLog = []
    for generation in range(maxGenerationNumber-1):
        generationGatesLog.append((generationMeansLog[generation]+generationMeansLog[generation+1])*0.5)
    xtickValues = []
    tempraw = []
    templin = []
    generationGatesLinear = []
    for row in range(rawData.shape[0]):
        tempLin = logicleData.iloc[row,:][~np.isnan(logicleData.iloc[row,:])].values
        tempRaw = rawData.iloc[row,:][~np.isnan(rawData.iloc[row,:])].values
        tempraw.append(tempRaw)
        templin.append(tempLin)
    newtempraw = np.concatenate(tempraw)
    newtemplin = np.concatenate(templin)
    xtickLabels = []
    for logxtick in logxticks:
        if(logxtick < 0):
            xtickLabels.append('$-10^'+str(int(np.log10(-1*logxtick)))+'$')
        elif(logxtick == 0):
            xtickLabels.append('0')
        else:
            xtickLabels.append('$10^'+str(int(np.log10(logxtick)))+'$')
    for tickval in logxticks:
        temprawtick,temprawtickindex = find_nearest(newtempraw,tickval)
        templintick = newtemplin[temprawtickindex]
        xtickValues.append(templintick)
    generationGatesLinear = [896.0, 825.0, 752.0, 677.0, 597.0, 510.0, 419.0]
    """
    T = 262143
    d = 4.5
    m = d*np.log(10)
    logicleData0 = np.array(generationMeansLog)
    logicleData = logicleData0[logicleData0<T]
    r = quantile(logicleData[logicleData<0], 0.05)
    w = 0
    #w = (m-np.log(T/abs(r)))/2
    generationGatesLinear = list(logicle(logicleData,T,m,r)*100)
    print(generationGatesLinear)
    #sys.exit(0)
    """
    #for gateval in generationGatesLog:
    #    generationGatesLinear.append(newtemplin[find_nearest(newtempraw,gateval)[1]])
    return generationGatesLinear,xtickValues,xtickLabels

def processProliferationData(folderName,logicleData,rawData,generationGatesLinear,xtickValues,xtickLabels,logicleDataStackedIndex):
    with open('inputFiles/experimentParameters-'+folderName+'.json') as f:
        experimentParameters = json.load(f)
    numConditions = experimentParameters[0][0]
    numTimePoints = experimentParameters[0][1]
    
    binNum = 256
    bandwidth = 0.08
    
    maxLinear = np.nanmax(logicleData)
    minLinear = np.nanmin(logicleData)-1
    print(maxLinear)
    print(minLinear)
    #One at a time
    #for val in range(logicleData.shape[0]):
    #Whole Timepoint at once
    plateGenerationData = []
    sns.set_palette(sns.color_palette("Blues", numTimePoints))
    tempVals = sns.color_palette()
    completeGenerationData = []
    sampleLogicleLengths = []
    sampleLogicleGatedLengths = []
    logicleLengths = []
    for iteration in range(int(logicleData.shape[0]/numTimePoints)):
        fig = plt.figure(figsize=(3*numTimePoints/3, 5*3))
        k=0
        fullLineList = []
        xlabels= []
        plt.suptitle(logicleData.iloc[iteration*numTimePoints,:].name[:-1])
        for val in range(iteration*numTimePoints,(iteration+1)*numTimePoints):
            tempWithNan = logicleData.iloc[val,:].values
            tempNoNan = tempWithNan[~np.isnan(tempWithNan)]
            rawWithNan = rawData.iloc[val,:].values
            rawNoNan = rawWithNan[~np.isnan(rawWithNan)]
            logicleLengths.append(tempNoNan.shape)
            ax = fig.add_subplot(3,math.ceil(numTimePoints/3),k+1)
            sns.distplot(tempNoNan,color=tempVals[k],kde=False,bins=binNum)
            ylabels = ax.get_yticks()
            plt.cla()
            plt.xlim([200,1000])
            plt.title(logicleData.iloc[val,:].name[-1])
            sns.kdeplot(tempNoNan,shade='True',color=tempVals[k],bw=bandwidth)
            sns.distplot(tempNoNan,color=tempVals[k],bins=binNum,kde_kws={"color": "k", "alpha": 0})
            plt.xticks(xtickValues,xtickLabels)
            ax.set_yticklabels(ylabels)
            lineList = []
            for j in range(len(generationGatesLinear)):
                lineList.append(draggable_lines(ax, "v", generationGatesLinear[j]))
            k+=1
            fullLineList.append(lineList)
        print('wat')
        plt.tight_layout()
        plt.show()
	#plt.show(block=False)
        #plt.close()
        sampleVal = iteration*numTimePoints
        conditionGenerationData = []
        for sampleLines in fullLineList:
            sampleLogicleDataNan = logicleData.iloc[sampleVal,:]
            sampleLogicleData = sampleLogicleDataNan[~np.isnan(sampleLogicleDataNan)]
            sampleLogicleLengths.append(sampleLogicleData.shape)
            sampleGenerationData = []
            sampleGenerationGates = [maxLinear]
            for line in sampleLines:
                sampleGenerationGates.append(line.getX())
            sampleGenerationGates.append(minLinear)
            for sampleEvent in sampleLogicleData:
                for generation in range(len(sampleGenerationGates)-1):
                    upperGate = sampleGenerationGates[generation]
                    lowerGate = sampleGenerationGates[generation+1]
                    if(sampleEvent > lowerGate and sampleEvent <= upperGate):
                        sampleGenerationData.append(generation)
            sampleLogicleGatedLengths.append(len(sampleGenerationData))
            sampleVal+=1
            conditionGenerationData+=sampleGenerationData
        names = [logicleData.iloc[iteration,:].name[:-2]]+[slice(None)]+[slice(None)]
        completeGenerationData += conditionGenerationData
    print('logicleLengths')
    for a in logicleLengths:
        print(a)
    print('sampleLogicleLengths')
    for a in sampleLogicleLengths:
        print(a)
    print('sampleLogicleGatedLengths')
    for a in sampleLogicleGatedLengths:
        print(a)
    #logicleDataStackedIndex = pickle.load(open('semiProcessedData/initialSingleCellDf-scale-'+folderName+'.pkl','rb'))['TCell_Gate'].index
    singleCellProliferationDf = pd.DataFrame(completeGenerationData,index=logicleDataStackedIndex,columns=['Generation'])
    print(singleCellProliferationDf)
    with open('semiProcessedData/singleCellDataFrame-proliferation-'+folderName+'.pkl', 'wb') as f:
        pickle.dump(singleCellProliferationDf,f)

def generateBulkProliferationStatistics(folderName):
    singleCellProliferationDfStacked = pickle.load(open('semiProcessedData/singleCellDataFrame-proliferation-'+folderName+'.pkl','rb'))
    numGenerations = len(singleCellProliferationDfStacked.loc[:,'Generation'].unique())
    singleCellProliferationDfUnstacked = singleCellProliferationDfStacked.unstack('Event')
    for i in range(singleCellProliferationDfStacked.index.nlevels-1):
        singleCellProliferationDfUnstacked = singleCellProliferationDfUnstacked.reindex(singleCellProliferationDfStacked.index.get_level_values(i).unique(),level=i)
    #Frequency
    generationFrequencyMatrix = np.zeros([singleCellProliferationDfUnstacked.shape[0],numGenerations])
    for row in range(singleCellProliferationDfUnstacked.shape[0]):
        proliferationDfPerCondition = singleCellProliferationDfUnstacked.iloc[row,:]
        freqTable = pd.value_counts(proliferationDfPerCondition.values.ravel()).sort_index()
        generations = [float(i) for i in list(freqTable.index.get_level_values(0))]
        for generation in generations:
            generationFrequencyMatrix[row,int(generation)] = freqTable[generation]
    generationFrequencyDataFrame = pd.DataFrame(generationFrequencyMatrix,index=singleCellProliferationDfUnstacked.index,columns=range(numGenerations)).astype(int)
    generationFrequencyDataFrame.columns.name = 'Generation'
    
    fractionDilutedMatrix = np.zeros([singleCellProliferationDfUnstacked.shape[0],1])
    precursorFrequencyMatrix = np.zeros([singleCellProliferationDfUnstacked.shape[0],1])
    proliferationIndexMatrix = np.zeros([singleCellProliferationDfUnstacked.shape[0],1])
    divisionIndexMatrix = np.zeros([singleCellProliferationDfUnstacked.shape[0],1])
    for row in range(generationFrequencyDataFrame.shape[0]):
        precursorFrequencyNumerator = 0
        proliferationIndexNumerator = 0
        for generation in range(1,numGenerations):
            precursorFrequencyNumerator+=generationFrequencyDataFrame.iloc[row,generation]/(2**generation)
            proliferationIndexNumerator+=generation*generationFrequencyDataFrame.iloc[row,generation]/(2**generation)
        precursorFrequencyDenominator = 0
        for generation in range(numGenerations):
            precursorFrequencyDenominator+=generationFrequencyDataFrame.iloc[row,generation]/(2**generation)
        
        fractionDiluted = np.sum(generationFrequencyDataFrame.iloc[row,1:])/np.sum(generationFrequencyDataFrame.iloc[row,0:])
        precursorFrequency = precursorFrequencyNumerator/precursorFrequencyDenominator
        proliferationIndex = proliferationIndexNumerator/precursorFrequencyNumerator
        divisionIndex = proliferationIndexNumerator/precursorFrequencyDenominator

        fractionDilutedMatrix[row,0] = fractionDiluted
        precursorFrequencyMatrix[row,0] = precursorFrequency
        proliferationIndexMatrix[row,0] = proliferationIndex
        divisionIndexMatrix[row,0] = divisionIndex
    allProliferationStatisticsList = []
    allProliferationMatricesList = [fractionDilutedMatrix,precursorFrequencyMatrix,proliferationIndexMatrix,divisionIndexMatrix]
    allProliferationStatisticNames = ['fractionDiluted','precursorFrequency','proliferationIndex','divisionIndex']
    for i in range(len(allProliferationMatricesList)):
        statisticDf = pd.DataFrame(allProliferationMatricesList[i],index=singleCellProliferationDfUnstacked.index)#,columns = [allProliferationStatisticNames[i]])
        statisticDfUnstacked = statisticDf.unstack('Time')
        for j in range(statisticDf.index.nlevels-1):
            statisticDfUnstacked = statisticDfUnstacked.reindex(statisticDf.index.get_level_values(j).unique(),level=j)
        allProliferationStatisticsList.append(statisticDfUnstacked)
    completeBulkProliferationDataFrame = pd.concat(allProliferationStatisticsList,keys=allProliferationStatisticNames,names=['Statistic'],axis=0)
    completeBulkProliferationDataFrame.columns = completeBulkProliferationDataFrame.columns.droplevel(0)
    print(completeBulkProliferationDataFrame)
    with open('semiProcessedData/proliferationStatisticPickleFile-'+folderName+'.pkl','wb') as f:
        pickle.dump(completeBulkProliferationDataFrame,f)
    return completeBulkProliferationDataFrame

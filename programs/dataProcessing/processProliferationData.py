#!/usr/bin/env python3  
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import matplotlib,statistics
#matplotlib.use('QT4Agg') 
import os,sys,pickle,json,math,itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('MacOSX') #default on my system
import seaborn as sns
from draggableLine import draggable_lines
from modifyDataFrames import returnModifiedDf
from logicle import logicle,quantile 

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

def createProliferationData(folderName,logicleDataFull,rawDataFull):
    idx = pd.IndexSlice
    logicleDataUnstacked = logicleDataFull.unstack('Event')
    rawDataUnstacked = rawDataFull.unstack('Event')
    
    emptyMatrixLogicle = np.zeros(logicleDataUnstacked.shape)
    emptyMatrixRaw = np.zeros(logicleDataUnstacked.shape)
    
    #Reindex unstacked dataframe by original ordering
    indexDfFull = pickle.load(open('semiProcessedData/cellStatisticPickleFile-'+folderName+'.pkl','rb'))
    indexDf = indexDfFull.loc[list(pd.unique(indexDfFull.index.get_level_values(0)))[0]].loc[list(pd.unique(indexDfFull.index.get_level_values(1)))[0]].loc[list(pd.unique(indexDfFull.index.get_level_values(2)))[0]]
    levelValues = []
    iteration = 0
    for row in range(indexDf.shape[0]):
        levelName = list(indexDf.iloc[row,:].name)
        for timepoint in list(pd.unique(indexDf.columns)):
            tp = tuple(levelName+[timepoint,slice(None)])
            emptyMatrixRaw[iteration,:] = rawDataUnstacked.loc[tp,slice(None)]
            emptyMatrixLogicle[iteration,:] = logicleDataUnstacked.loc[tp,slice(None)]
            levelValues.append(tp[:-1])
            iteration+=1
    multiIndex = pd.MultiIndex.from_tuples(levelValues, names=indexDf.index.names+[indexDf.columns.name])
    
    logicleData = pd.DataFrame(emptyMatrixLogicle,index=multiIndex,columns=logicleDataUnstacked.columns)
    rawData = pd.DataFrame(emptyMatrixRaw,index=multiIndex,columns=rawDataUnstacked.columns)
    logicleData.columns.name = 'Event'
    rawData.columns.name = 'Event'
    
    """
    for row in range(logicleData.shape[0]):
        rowdf = logicleData.iloc[row,:].name
        if rowdf[-1] == 0.5:
            print(rowdf)
            tempdf = logicleData.loc[rowdf,:]
            temp = tempdf[~np.isnan(tempdf)]
            print(temp.shape)
    sys.exit(0)
    """
    return logicleData,rawData

def returnGatesAndTicks(logicleData,rawData):
    maxGenerationNumber = 7
    #Need to enter this programatically
    generationZeroGFI = 60000
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
    #generationGatesLinear = [896.0, 825.0, 752.0, 677.0, 597.0, 510.0, 419.0]
    for gateval in generationGatesLog:
        generationGatesLinear.append(newtemplin[find_nearest(newtempraw,gateval)[1]])
    return generationGatesLinear,xtickValues,xtickLabels

#def processProliferationData(folderName,logicleData,rawData,generationGatesLinear,xtickValues,xtickLabels,logicleDataStackedIndex):
def processProliferationData(folderName,logicleData,rawData,generationGatesLinear,xtickValues,xtickLabels,logicleDataStacked):
    with open('inputFiles/experimentParameters-'+folderName+'.json') as f:
        experimentParameters = json.load(f)
    numConditions = experimentParameters[0][0]
    numTimePoints = experimentParameters[0][1]
    
    binNum = 256
    #Might need to calculate bandwidth manually as there seems to be significan fluctuation in the numbered required for distinct peaks in the kde of the ctv values
    bandwidth = 15
    #bandwidth = 25
    #bandwidth = 'scott'
    
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
    print(logicleData.iloc[0,:])
    print(logicleData.iloc[24*8,:])
    #for iteration in range(2):
    for iteration in range(int(logicleData.shape[0]/numTimePoints)):
        #fig = plt.figure(figsize=(3*numTimePoints/3, 5*3))
        fig = plt.figure()
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
            if math.ceil(numTimePoints/3) > 3:
                ax = fig.add_subplot(3,math.ceil(numTimePoints/3),k+1)
            else:
                ax = fig.add_subplot(math.ceil(numTimePoints/3),3,k+1)
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
        #plt.tight_layout()
        #plt.show()
        plt.close()
            
        print('iter')
        print(iteration)
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
            sampleVal+=1
            #if sampleVal == numTimePoints*0 or sampleVal == numTimePoints*8:
            conditionGenerationData.append(sampleGenerationData)
            #conditionGenerationData+=sampleGenerationData
        completeGenerationData.append(conditionGenerationData)
        #completeGenerationData += conditionGenerationData
    
    temp1 = logicleDataStacked.copy()
    temp2 = pd.concat([logicleDataStacked,temp1],axis=1,keys=['GFI_1','GFI_2'])
    newList = []
    newIndex = []
    sampleLengths = []
    row = 0
    for conditionList in completeGenerationData:
        for timelist in conditionList:
            newList.append(timelist)
            wat = logicleData.iloc[row,:].name
            newIndex.append(wat)
            sampleLengths.append(len(timelist))
            row+=1
            #row+=len(timelist)
    completeGenerationData = list(itertools.chain.from_iterable(newList))
    realIndex = []
    i=0
    for sampleTuple in newList:
        for eventNumber in range(1,sampleLengths[i]+1):
            eventIndex = list(newIndex[i])
            eventIndex.append(eventNumber)
            realIndex.append(eventIndex)
        i+=1
        print(i)
    idx = pd.IndexSlice
    newMultiIndex = pd.MultiIndex.from_tuples(realIndex,names=tuple(logicleDataStacked.index.names))
    singleCellProliferationDf = pd.DataFrame(completeGenerationData,index=newMultiIndex,columns=['Generation'])
    
    #Merge initial two peaks to correct division index
    generationValues = np.subtract(singleCellProliferationDf.values,1).clip(0,7)
    singleCellProliferationDf.values[:] = generationValues
    
    print(singleCellProliferationDf)
    with open('semiProcessedData/singleCellDataFrame-proliferation-'+folderName+'.pkl', 'wb') as f:
        pickle.dump(singleCellProliferationDf,f)

def generateBulkProliferationStatistics(folderName,experimentNumber):
    singleCellProliferationDfStacked = pickle.load(open('semiProcessedData/singleCellDataFrame-proliferation-'+folderName+'.pkl','rb'))
    idx = pd.IndexSlice
    numGenerations = len(singleCellProliferationDfStacked.loc[:,'Generation'].unique())
    singleCellProliferationDfUnstacked = singleCellProliferationDfStacked.unstack('Event')
    #Reindex unstacked dataframe by original ordering
    indexDfFull = pickle.load(open('semiProcessedData/cellStatisticPickleFile-'+folderName+'.pkl','rb'))
    indexDf = indexDfFull.loc[list(pd.unique(indexDfFull.index.get_level_values(0)))[0]].loc[list(pd.unique(indexDfFull.index.get_level_values(1)))[0]].loc[list(pd.unique(indexDfFull.index.get_level_values(2)))[0]]
    levelValues = []
    iteration = 0
    emptyMatrixLogicle = np.zeros(singleCellProliferationDfUnstacked.shape)
    for row in range(indexDf.shape[0]):
        levelName = list(indexDf.iloc[row,:].name)
        for timepoint in list(pd.unique(indexDf.columns)):
            tp = tuple(levelName+[timepoint,slice(None)])
            emptyMatrixLogicle[iteration,:] = singleCellProliferationDfUnstacked.loc[tp,slice(None)]
            levelValues.append(tp[:-1])
            iteration+=1
    multiIndex = pd.MultiIndex.from_tuples(levelValues, names=indexDf.index.names+[indexDf.columns.name])
    
    singleCellProliferationDfUnstacked = pd.DataFrame(emptyMatrixLogicle,index=multiIndex,columns=singleCellProliferationDfUnstacked.columns)
    singleCellProliferationDfUnstacked.columns.name = 'Event'
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
    for i in range(generationFrequencyDataFrame.shape[1]):
        currentGenerationFrequency = np.array(generationFrequencyDataFrame.iloc[:,i].values.ravel())
        allProliferationStatisticNames.append('rawFrequencyGen'+str(i))
        allProliferationMatricesList.append(currentGenerationFrequency)
        allProliferationStatisticNames.append('normalizedFrequencyGen'+str(i))
        allProliferationMatricesList.append(np.divide(currentGenerationFrequency,2**i))
    for i in range(len(allProliferationMatricesList)):
        statisticDf = pd.DataFrame(allProliferationMatricesList[i],index=singleCellProliferationDfUnstacked.index)#,columns = [allProliferationStatisticNames[i]])
        statisticDfUnstacked = statisticDf.unstack('Time')
        indexDfFull = pickle.load(open('semiProcessedData/cellStatisticPickleFile-'+folderName+'.pkl','rb'))
        indexDf = indexDfFull.loc[list(pd.unique(indexDfFull.index.get_level_values(0)))[0]].loc[list(pd.unique(indexDfFull.index.get_level_values(1)))[0]].loc[list(pd.unique(indexDfFull.index.get_level_values(2)))[0]]
        levelValues = []
        iteration = 0
        emptyMatrixLogicle = np.zeros(statisticDfUnstacked.shape)
        for row in range(indexDf.shape[0]):
            levelName = list(indexDf.iloc[row,:].name)
            tp = tuple(levelName)
            emptyMatrixLogicle[iteration,:] = statisticDfUnstacked.loc[tp,:]
            levelValues.append(tp)
            iteration+=1
        multiIndex = pd.MultiIndex.from_tuples(levelValues, names=indexDf.index.names)
        
        statisticDfUnstacked = pd.DataFrame(emptyMatrixLogicle,index=multiIndex,columns=statisticDfUnstacked.columns)
        statisticDfUnstacked.columns.name = 'Time'
        allProliferationStatisticsList.append(statisticDfUnstacked)
    completeBulkProliferationDataFrame = pd.concat(allProliferationStatisticsList,keys=allProliferationStatisticNames,names=['Statistic'],axis=0)
    completeBulkProliferationDataFrame.columns = completeBulkProliferationDataFrame.columns.droplevel(0)
    with open('semiProcessedData/proliferationStatisticPickleFile-'+folderName+'.pkl','wb') as f:
        pickle.dump(completeBulkProliferationDataFrame,f)
    modifiedDataFrame = returnModifiedDf(experimentNumber,completeBulkProliferationDataFrame,'prolif')
    with open('semiProcessedData/proliferationStatisticPickleFile-'+folderName+'-modified.pkl','wb') as f:
        pickle.dump(modifiedDataFrame,f)
    return modifiedDataFrame

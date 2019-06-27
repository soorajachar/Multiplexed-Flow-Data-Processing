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
from scipy.stats import gaussian_kde
plt.switch_backend('MacOSX') #default on my system
import seaborn as sns
from draggableLine import draggable_lines
from multiDragLine import proliferation_gate_lines
from modifyDataFrames import returnModifiedDf
from logicle import logicle,quantile 
from datetime import datetime
import subprocess
    
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
    maxGenerationNumber = 5
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
    #generationGatesLinear = [896.0, 825.0, 752.0, 677.0, 597.0, 510.0, 419.0]
    for gateval in generationGatesLog:
        generationGatesLinear.append(newtemplin[find_nearest(newtempraw,gateval)[1]])
    return generationGatesLinear

def returnRestartIndex(folderName,logicleDataStacked):
    restartdf = pickle.load(open('semiProcessedData/temp-proliferation-'+folderName+'.pkl', 'rb'))
    conditionLevelNames = logicleDataStacked.index.names[:-2] 
    logicleDataUnstacked = logicleDataStacked.groupby(level=conditionLevelNames,sort=False).first().to_frame('GFI')
    if -1 in restartdf.ravel():
        restartrow = 0
        for row in range(restartdf.shape[0]):
            if restartdf[row] == -1:
                restartrow = row
                break
        
        newdf = pd.DataFrame(restartdf,index=logicleDataStacked.index)
        restartrowname = newdf.iloc[restartrow,:].name
        stackedrestartrow = 0
        for row in range(logicleDataUnstacked.shape[0]):
            conditionLevelValues = logicleDataUnstacked.iloc[row,:].name
            conditionDf = tempdf.xs(conditionLevelValues,level=conditionLevelNames)
            if logicleDataUnstacked.iloc[row,:].name == restartrowname[:-2]:
                stackedrestartrow = row
    else:
        stackedrestartrow = logicleDataUnstacked.shape[0]-1
    return stackedrestartrow

def kde_scipy(x, x_grid,bnw,**kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    kde = gaussian_kde(x, bw_method=bnw / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)
    
def automatedProliferationGating(data):
    bnw=10
    x_grid = np.linspace(0, 1000, 1000)
    startingGateList = []
    for timepoint in pd.unique(data.index.get_level_values('Time')):
        x = data.loc[timepoint].values.ravel().astype(float)
        pdf = kde_scipy(x, x_grid,bnw)

        diff = np.gradient(pdf)
        sdiff = np.sign(diff)
        lastlocalmin = np.where(sdiff[:-1] != sdiff[1:])[0][-2]
        startingGateList.append(lastlocalmin)
    return startingGateList

def processProliferationData(folderName,logicleDataStacked,rawDataStacked,manual,plot):

    startingGateLogicle = 725
    offset = 2 
    binNum = 256
    maxLinear = np.nanmax(logicleDataStacked)
    minLinear = np.nanmin(logicleDataStacked)-1
    
    timepoints = pd.unique(logicleDataStacked.index.get_level_values('Time'))
    conditionLevelNames = logicleDataStacked.index.names[:-2] 
    
    logicleDataUnstacked = logicleDataStacked.groupby(level=conditionLevelNames,sort=False).first().to_frame('GFI')
    experimentGenerationData = np.ones((logicleDataStacked.shape[0],1))*-1
    
    #Restart if proliferation processing stopped in the middle for some reason
    fileList = os.listdir('semiProcessedData')
    if 'temp-proliferation-'+folderName+'.pkl' in fileList:
        experimentGenerationData = pickle.load(open('semiProcessedData/temp-proliferation-'+folderName+'.pkl', 'rb'))
        restartindex = returnRestartIndex(folderName,logicleDataStacked)
        restartrow=0
        for row in range(restartindex):
            conditionLevelValues = logicleDataUnstacked.iloc[row,:].name
            conditionDf = logicleDataStacked.xs(conditionLevelValues,level=conditionLevelNames).to_frame('GFI')
            restartrow+=conditionDf.shape[0]
    else:
        experimentGenerationData = np.zeros((logicleDataStacked.shape[0],1))
        restartrow = 0
        restartindex=0
    i=restartrow
    for row in range(restartindex,logicleDataUnstacked.shape[0]):
        startTime = datetime.now()
        sns.set_palette(sns.color_palette("Purples", len(timepoints)))
        #Grab appropriate subsettedlevel of stacked ctv dataframe
        conditionLevelValues = logicleDataUnstacked.iloc[row,:].name
        conditionDf = logicleDataStacked.xs(conditionLevelValues,level=conditionLevelNames).to_frame('GFI')
        rawDf = rawDataStacked.xs(conditionLevelValues,level=conditionLevelNames).to_frame('GFI')
        
        #Manual adjustment of starting gate per timepoint
        if manual:
            #Create long dataframe for use with seaborn
            plottingDf = conditionDf.reset_index()
            
            #Plot facetgrid of histograms of ctv values of all times for current condition
            g = sns.FacetGrid(plottingDf,sharey=False,legend_out=True,col='Time',hue='Time',col_wrap=6)
            g.map(sns.distplot,'GFI',kde=False,bins=256)
            
            #Plot draggable lines in each histogram; used to set boundary between generation zero and rest proliferation generations (varies slightly overtime)
            linelist = []
            #startingGate = rawDf.values[find_nearest(conditionDf,startingGateLogicle)[1]][0]
            #print(startingGate)
            #generationGates = returnGates(conditionDf,rawDf,startingGate)
            for axis in g.axes:
                #linelist.append(proliferation_gate_lines(axis, startingGate,0,1000,generationGates,len(generationGates)))
                linelist.append(draggable_lines(axis, "v", startingGateLogicle,0,1000))
            #g.add_legend()
            plt.suptitle(conditionLevelValues)
            plt.subplots_adjust(top=0.95)
            plt.show()
        else:
            linelist = automatedProliferationGating(conditionDf)

        #Grab ctv gfi values of genzero boundary lines for each timepoint in the current condition
        timepointGenerationDfList = []
        conditionGenerationData = np.zeros((logicleDataStacked.loc[conditionLevelValues].shape[0],1))
        j=0
        gateList = []
        for line,timepoint in zip(linelist,timepoints):
            sampleGenerationGates = [maxLinear]
            if manual:
                startingGate = rawDf.values[find_nearest(conditionDf,line.getX())[1]][0]
            else:
                startingGate = rawDf.values[find_nearest(conditionDf,line+offset)[1]][0]
            generationGatesLinear = returnGates(conditionDf,rawDf,startingGate)
            for generationGate in generationGatesLinear:
                sampleGenerationGates.append(generationGate)
            sampleGenerationGates.append(minLinear)
            sampleLogicleData = conditionDf.xs([timepoint],level=['Time'])
            timepointGenerationDf = np.zeros(sampleLogicleData.shape)
            for sampleEvent,row in zip(sampleLogicleData.values,range(timepointGenerationDf.shape[0])):
                for generation in range(len(sampleGenerationGates)-1):
                    upperGate = sampleGenerationGates[generation]
                    lowerGate = sampleGenerationGates[generation+1]
                    if(sampleEvent > lowerGate and sampleEvent <= upperGate):
                        timepointGenerationDf[row] = generation
            conditionGenerationData[j:j+timepointGenerationDf.shape[0],:] = timepointGenerationDf
            j+=timepointGenerationDf.shape[0]
            gateList.append(sampleGenerationGates)

        #Combine timepoints in condition to single condition generation array
        experimentGenerationData[i:i+conditionGenerationData.shape[0],:] = conditionGenerationData
        i+=conditionGenerationData.shape[0]
        
        #To see gates overlaid on histograms       
        if plot:
            #Create long dataframe for use with seaborn
            plottingDf = conditionDf.reset_index()
            g = sns.FacetGrid(plottingDf,sharey=False,legend_out=True,col='Time',hue='Time',col_wrap=6)
            g.map(sns.distplot,'GFI',kde=False,bins=256)
            linelist = []
            for axis,lines in zip(g.axes,gateList):
                for lineval in lines:
                    axis.axvline(lineval,linestyle=':',color='k')
            g.add_legend()
            plt.suptitle(conditionLevelValues)
            plt.subplots_adjust(top=0.95)
            plt.show()
        
        if manual: 
            with open('semiProcessedData/temp-proliferation-'+folderName+'.pkl', 'wb') as f:
                pickle.dump(experimentGenerationData,f)
        
        print('Condition '+str('-'.join(list(conditionLevelValues)))+' processed')
        print(datetime.now() - startTime)
        

    #Combine conditions in experiment to single experiment generation array
    singleCellProliferationDf = pd.DataFrame(experimentGenerationData,index=logicleDataStacked.index,columns=['Generation'])
    print(singleCellProliferationDf)

    #Remove restart dataframe
    subprocess.run(['rm','semiProcessedData/temp-proliferation-'+folderName+'.pkl'])
     
    with open('semiProcessedData/singleCellDataFrame-proliferation-'+folderName+'.pkl', 'wb') as f:
        pickle.dump(singleCellProliferationDf,f)

def generateBulkProliferationStatistics(folderName,experimentNumber):
    singleCellProliferationDfStacked = pickle.load(open('semiProcessedData/singleCellDataFrame-proliferation-'+folderName+'.pkl','rb'))
    idx = pd.IndexSlice
    numGenerations = len(singleCellProliferationDfStacked.loc[:,'Generation'].unique())
    singleCellProliferationDfUnstacked = singleCellProliferationDfStacked.unstack('Event')
    #Reindex unstacked dataframe by original ordering
    #conditionLevelNames = singleCellProliferationDfStacked.index.names[:-1]
    #indexDf = singleCellProliferationDfStacked.groupby(level=conditionLevelNames,sort=False).first()
    #print(indexDf)
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

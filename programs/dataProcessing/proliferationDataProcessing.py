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
from multiDragLineAnchorKDE import proliferation_gate_lines
from matplotlib.widgets import RadioButtons,Button,CheckButtons,TextBox
from modifyDataFrames import returnModifiedDf
from logicle import logicle,quantile 
from datetime import datetime
from miscFunctions import find_nearest,returnTicks
import subprocess

#Button width conserved across gui figures
buttonWidth = 0.1/2
buttonLength = 0.075/2
buttonXStart = 0.5-(0.01+buttonWidth)
buttonYStart = 0.01
    
def createProliferationSingleCellDataFrame(folderName,secondPath,experimentNumber,useModifiedDf):
    logicleDataProliferation = pickle.load(open('semiProcessedData/logicleProliferationDf.pkl','rb'))
    rawDataProliferation = pickle.load(open('semiProcessedData/rawProliferationDf.pkl','rb'))
    processProliferationData(folderName,logicleDataProliferation,rawDataProliferation)
    bulkprolifdf = generateBulkProliferationStatistics(folderName,experimentNumber)
    return bulkprolifdf

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
    highGFICutoff = 500
    highGFIPercentCutoff = 0.1
    extraneousExtremaCutoff = 100
    x_grid = np.linspace(0, 1000, 1000)
    startingGateList = []
    startingGateLogicle = 725
    for timepoint in pd.unique(data.index.get_level_values('Time')):
        x = data.loc[timepoint].values.ravel().astype(float)
        percentAboveHalf = x[np.where(x>highGFICutoff)].shape[0]/x.shape[0]
        #Parent Generation Present
        if percentAboveHalf >= highGFIPercentCutoff:
            #KDE
            pdf = kde_scipy(x, x_grid,bnw)
            #Derivative of kde
            diff = np.gradient(pdf)
            #sign of derivative for each element (-1 if - etc.)
            sdiff = np.sign(diff)
            #Locations where sign changes ([0] just removes array)
            #lastlocalmin = np.where(sdiff[:-1] != sdiff[1:])[0][-2]
            lastlocalmin = np.where(sdiff[:-1] != sdiff[1:])[0][-1]
            startingGateList.append(lastlocalmin)
        #Parent Generation Not Present
        else:
            startingGateList.append(startingGateLogicle)
    return startingGateList

def conditionOrTimepointCTVGating_GUIWindow():
    #Grab levels that will be used within a figure
    fig = plt.figure()
    plt.axis('off')

    rax = plt.axes([0.4, 0.5, 0.05*5,0.05*5])
    radio = RadioButtons(rax,['Group All','Group Conditions','Group Timepoints'],activecolor='black') 
    plt.text(0.5, 1.2,'Do you want to group all timepoints or all conditions together for initial CTV Gating?',ha='center')
    rax.spines['bottom'].set_visible(False)
    rax.spines['left'].set_visible(False)
    rax.spines['right'].set_visible(False)
    rax.spines['top'].set_visible(False)

    class GUIButtons(object):
        def OKradio1(self, event):
            if radio.value_selected == 'Group All':
                withinFigureBoolean = [True,False,False]
            elif radio.value_selected == 'Group Conditions':
                withinFigureBoolean = [False,True,False]
            else:
                withinFigureBoolean = [False,False,True]
            plt.close()
            with open('semiProcessedData/gui-groupingBool.pkl','wb') as f:
                pickle.dump(withinFigureBoolean,f)
        def Quit(self, event):
            sys.exit(0)

    callback = GUIButtons()
    axOK = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
    axQuit = plt.axes([buttonXStart+buttonWidth+0.01,buttonYStart, buttonWidth, buttonLength])
    bOK = Button(axOK, 'OK')
    bOK.on_clicked(callback.OKradio1)
    bQuit = Button(axQuit, 'Quit')
    bQuit.on_clicked(callback.Quit)
    plt.show()

#Gate manually
def processProliferationData(folderName,logicleDataStacked,rawDataStacked):
    startingGateLogicle = 725
    binNum = 256
    
    conditionOrTimepointCTVGating_GUIWindow()
    conditionOrTimepointFile = pickle.load(open('semiProcessedData/gui-groupingBool.pkl','rb'))
    
    groupedGateList = []

    rowAsAllLevelNames = logicleDataStacked.index.names[:-1]
    rowAsTimePointLevelNames = logicleDataStacked.index.names[-2]
    rowAsConditionLevelNames = logicleDataStacked.index.names[:-2] 
    conditionTuple = []
    for tpl in logicleDataStacked.groupby(level=rowAsConditionLevelNames,sort=False):
        conditionTuple.append(tpl[0])
    timepointTuple = []
    for tpl in logicleDataStacked.groupby(level=rowAsTimePointLevelNames,sort=False):
        timepointTuple.append(tpl[0])
        
    plt.ioff()
    
    #Group all samples together (no iteration at all)
    if conditionOrTimepointFile[0]:
        levelNames = rowAsAllLevelNames
        outerVariableValues = conditionTuple
        innerVariableValues = timepointTuple
        groupVariable = 'Condition-Time'
        titleVariable = 'All'
        iterationRange = 1
    else:
        #Group conditions together (iterate through each timepoint):
        if conditionOrTimepointFile[1]:
            levelNames = rowAsTimePointLevelNames
            outerVariableValues = timepointTuple
            innerVariableValues = conditionTuple
            groupVariable = 'Condition'
            titleVariable = 'Time'
        #Group timepoints together (iterate through each condition):
        else:
            levelNames = rowAsConditionLevelNames
            outerVariableValues = conditionTuple
            innerVariableValues = timepointTuple
            groupVariable = 'Time'
            titleVariable = 'Condition'
        logicleDataUnstacked = logicleDataStacked.groupby(level=levelNames,sort=False).first().to_frame('GFI')
        iterationRange = logicleDataUnstacked.shape[0] 
    
    
    for row in range(iterationRange):
        if conditionOrTimepointFile[0]:
            currentDf = logicleDataStacked.to_frame('GFI')
            rawDf = rawDataStacked.to_frame('GFI')
        else:
            currentLevelValues = logicleDataUnstacked.iloc[row,:].name
            currentDf = logicleDataStacked.xs(currentLevelValues,level=levelNames).to_frame('GFI')
            rawDf = rawDataStacked.xs(currentLevelValues,level=levelNames).to_frame('GFI')
        plottingDf = currentDf.reset_index()
        if conditionOrTimepointFile[0] or conditionOrTimepointFile[1]:
            nonEventLevelNames = currentDf.index.names[:-1]
            groupingColumn = []
            for conditionGroupbyTuple in currentDf.groupby(level=nonEventLevelNames,sort=False):
                conditionGroupbyTupleNew = [str(i) for i in conditionGroupbyTuple[0]]
                conditionName = '-'.join(conditionGroupbyTupleNew)
                for event in range(conditionGroupbyTuple[1].shape[0]):
                    groupingColumn.append(conditionName)
            plottingDf[groupVariable] = groupingColumn
        sns.set_palette(sns.color_palette("Purples", len(pd.unique(plottingDf[groupVariable]))))
        #Plot facetgrid of histograms of ctv values of all times for current condition
        g = sns.FacetGrid(plottingDf,legend_out=True,hue=groupVariable,height=6)
        g.map(sns.kdeplot,'GFI')
        
        #Add gating lines
        startingGateRaw = rawDf.values[find_nearest(currentDf,startingGateLogicle)[1]][0]
        logicleGates = returnGates(currentDf,rawDf,startingGateRaw)
        maxY = list(plt.gca().get_yticks())[-1]
        currentGroupedGates = proliferation_gate_lines(plt.gca(),startingGateLogicle,-1,maxY,logicleGates[1:],len(logicleGates)-1)
        #Add title
        if conditionOrTimepointFile[0]:
            plt.title(titleVariable)
        else:
            if not isinstance(currentLevelValues, (tuple,)):
                currentLevelValues = [str(currentLevelValues)]
            plt.title('-'.join(currentLevelValues)+' ('+titleVariable+' '+str(row+1)+'/'+str(logicleDataUnstacked.shape[0])+')')
         
        #Add continue gating button or subset deeper button
        class GUIButtons(object):
            def NextGateSet(self, event):
                for names in range(len(pd.unique(plottingDf[groupVariable]))):
                    groupedGateList.append(currentGroupedGates)
                plt.close()
                with open('semiProcessedData/gui-splitBool.pkl','wb') as f:
                    pickle.dump(False,f)
            def NextGateSetAfterSplit(self,event):
                plt.close()
            def SplitGateSet(self,event):
                with open('semiProcessedData/gui-splitBool.pkl','wb') as f:
                    pickle.dump(True,f)
                plt.close()
            def Quit(self, event):
                sys.exit(0)
        callback = GUIButtons()
        axNext = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
        axSplit = plt.axes([buttonXStart+buttonWidth+0.01, buttonYStart, buttonWidth, buttonLength])
        axQuit = plt.axes([buttonXStart+2*buttonWidth+0.02,buttonYStart, buttonWidth, buttonLength])
        bNext = Button(axNext, 'OK')
        bNext.on_clicked(callback.NextGateSet)
        bSplit = Button(axSplit, 'Split')
        bSplit.on_clicked(callback.SplitGateSet)
        bQuit = Button(axQuit, 'Quit')
        bQuit.on_clicked(callback.Quit)
        xtickValues,xtickLabels = returnTicks([-1000,100,1000,10000,100000])
        axis = g.fig.get_axes()[0]
        axis.set_xticks(xtickValues)
        axis.set_xticklabels(xtickLabels)
        plt.subplots_adjust(top=0.95)
        plt.show()
        #Create deeper gates
        splitBool = pickle.load(open('semiProcessedData/gui-splitBool.pkl','rb'))
        if splitBool:
            sns.set_palette(sns.color_palette("Purples_r", len(pd.unique(plottingDf[groupVariable]))))
            #Plot facetgrid of histograms of ctv values of all times for current condition
            colwraplength = 8 
            g = sns.FacetGrid(plottingDf,legend_out=True,col=groupVariable,col_wrap=colwraplength)
            g.map(sns.kdeplot,'GFI')
            for axis in g.axes:
                maxY = list(axis.get_yticks())[-1]
                groupedGateList.append(proliferation_gate_lines(axis,startingGateLogicle,-1,maxY,logicleGates[1:],len(logicleGates)-1))
            callback = GUIButtons()
            axNext = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
            axQuit = plt.axes([buttonXStart+buttonWidth+0.01, buttonYStart, buttonWidth, buttonLength])
            bNext = Button(axNext, 'OK')
            bNext.on_clicked(callback.NextGateSet)
            bQuit = Button(axQuit, 'Quit')
            bQuit.on_clicked(callback.Quit)
            xtickValues,xtickLabels = returnTicks([-1000,1000,10000,100000])
            for axis,i in zip(g.fig.get_axes(),range(len(g.fig.get_axes()))):
                axis.set_xticks(xtickValues)
                axis.set_xticklabels(xtickLabels)
            plt.show()
     
    i = 0
    row = 0
    singleCellProliferationDf = pd.DataFrame(np.zeros(logicleDataStacked.shape),logicleDataStacked.index,columns=['Generation'])
    for outerVariable in outerVariableValues:
        for innerVariable in innerVariableValues:
            groupedGate = groupedGateList[i]
            sampleGenerationGates = groupedGate.getAllX()
            if groupVariable == 'Condition':
                indexingTuple = tuple(list(innerVariable)+[outerVariable,slice(None)])
            else:
                indexingTuple = tuple(list(outerVariable)+[innerVariable,slice(None)])
            ctvValues = logicleDataStacked.loc[indexingTuple]
            generationValues = np.zeros(ctvValues.shape)
            for sampleEvent,row in zip(ctvValues,range(ctvValues.shape[0])):
                for generation in range(len(sampleGenerationGates)-1):
                    upperGate = sampleGenerationGates[generation]
                    lowerGate = sampleGenerationGates[generation+1]
                    if(sampleEvent > lowerGate and sampleEvent <= upperGate):
                        generationValues[row] = generation
            singleCellProliferationDf.loc[indexingTuple,:] = generationValues
            print(indexingTuple)

    print(singleCellProliferationDf)
    
    #Remove temporary gui files 
    subprocess.run(['rm','semiProcessedData/gui-groupingBool.pkl'])
    subprocess.run(['rm','semiProcessedData/gui-splitBool.pkl'])
     
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
    temp = indexDfFull.reset_index()
    temp = singleCellProliferationDfUnstacked.reset_index()
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

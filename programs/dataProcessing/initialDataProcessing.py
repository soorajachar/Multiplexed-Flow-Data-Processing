#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import json,pickle,math,matplotlib,sys,os,string
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from itertools import groupby
from miscFunctions import Hill,InverseHill,r_squared,cleanUpFlowjoCSV,extractValues  
from modifyDataFrames import returnModifiedDf
import cytokineDataProcessing,singleCellDataProcessing,cellDataProcessing

dataTypeLevelNames = {'cyt':['Cytokine'],'cell':['CellType','Marker','Statistic'],'prolif':['Statistic'],'singlecell':['CellType']}
dataTypeDataFrameFileNames = {'cyt':'cytokineConcentrationPickleFile','cell':'cellStatisticPickleFile','prolif':'proliferationStatisticPickleFile','singlecell':'initialSingleCellPickleFile'}

def returnMultiIndex(sortedData,sortedFiles,dataType,folderName):
    if(dataType == 'cyt'):
        newMultiIndex = cytokineDataProcessing.parseCytokineCSVHeaders(pd.read_csv('semiProcessedData/A1_'+dataType+'.csv').columns)
    elif(dataType == 'cell'):
        panelData = pd.read_csv('inputFiles/antibodyPanel-'+folderName+'.csv',)
        newMultiIndex = cellDataProcessing.parseCellCSVHeaders(pd.read_csv('semiProcessedData/A1_'+dataType+'.csv').columns,panelData)
    elif(dataType == 'singlecell'):
        #Grabs a file from samples to read marker names off of
        cellTypeList = []
        for fileName in os.listdir('semiProcessedData/singleCellData/A1/'):
            if 'DS' not in fileName:
                cellTypeList.append(fileName)
        newMultiIndex = singleCellDataProcessing.produceSingleCellHeaders(cellTypeList)
    elif(dataType == 'cytcorr'):
        newMultiIndex = []
    if dataType != 'singlecell':
        return sortedData,newMultiIndex
    else:
        return sortedFiles,newMultiIndex

def tilePlateLayouts(experimentParameters,levelLayouts):
    #Retile smaller layouts to fill overall plate dimensions
    print(levelLayouts)
    baseLevelLayout = levelLayouts[0]
    baseLevelLayoutNoZeros = extractValues(baseLevelLayout,0,True)
    emptyLevelLayout = np.zeros(baseLevelLayout.shape)
    fullFilledLevelLayouts = [baseLevelLayout]
    for sparseLevelLayout in levelLayouts[1:]:
        verticalRepetitionsNeeded = int(baseLevelLayoutNoZeros.shape[0]/sparseLevelLayout.shape[0])
        horizontalRepetitionsNeeded = int(baseLevelLayoutNoZeros.shape[1]/sparseLevelLayout.shape[1])
        filledLevelLayout = np.tile(sparseLevelLayout,(verticalRepetitionsNeeded,horizontalRepetitionsNeeded))
        fullFilledLevelLayout = emptyLevelLayout.copy()
        for row in range(filledLevelLayout.shape[0]):
            for col in range(filledLevelLayout.shape[1]):
                fullFilledLevelLayout[row,col] = filledLevelLayout[row,col]
        fullFilledLevelLayouts.append(fullFilledLevelLayout)

    #Tile across extra plates to span all column values
    experimentLevelLayoutDict = {}
    if (experimentParameters['numPlates'] > 2 and experimentParameters['paired']) or (experimentParameters['numPlates'] > 1 and not experimentParameters['paired']):
        if experimentParameters['paired']:
            numReps = int((experimentParameters['numPlates'] - 2)/2)+1
        else:
            numReps = experimentParameters['numPlates']
        for fullFilledLevelLayout,currentConditionLevel in zip(fullFilledLevelLayouts,experimentParameters['allLevelNames']):
            if currentConditionLevel == experimentParameters['columnVariableName']:
                tempList = []
                #Top right value will be the highest level index on the first plate
                levelOffsetMatrix = extractValues(fullFilledLevelLayout,0,True)
                levelOffset = levelOffsetMatrix[0,levelOffsetMatrix.shape[1]-1]
                for rep in range(numReps):
                    baseMatrix = np.zeros(fullFilledLevelLayout.shape)
                    newLevelOffsetMatrix = levelOffsetMatrix+rep*levelOffset
                    fullLevelOffsetMatrix = baseMatrix[:newLevelOffsetMatrix.shape[0],:newLevelOffsetMatrix.shape[1]] = newLevelOffsetMatrix
                    tempList.append(fullLevelOffsetMatrix)
                experimentLevelLayout = np.hstack(tempList)
            else:
                experimentLevelLayout = np.tile(fullFilledLevelLayout,(1,numReps)) 
            experimentLevelLayoutDict[currentConditionLevel] = experimentLevelLayout
    #No tiling across plates needed    
    else:
        for fullFilledLevelLayout,currentConditionLevel in zip(fullFilledLevelLayouts,experimentParameters['allLevelNames']):
            experimentLevelLayoutDict[currentConditionLevel] = fullFilledLevelLayout
            
    return experimentLevelLayoutDict

def arrangeColumVariableBasedCoordinates(experimentParameters,experimentLevelLayoutDict):

    columnVariableLevelLayout = experimentLevelLayoutDict[experimentParameters['columnVariableName']]
    allColumnVariableCoordinates = []
    if experimentParameters['paired'] and experimentParameters['replicateWise']:
        splitColumnVariableLevelLayoutList = np.vsplit(columnVariableLevelLayout,2)
        columnVariableOffset = splitColumnVariableLevelLayoutList[0].shape[0]
        for columVariable in range(1,experimentParameters['numColumnLevelValues']+1):
            currentColumnVariableCoordinates = []
            for index in range(len(splitColumnVariableLevelLayoutList)):
                for col in range(splitColumnVariableLevelLayout[index].shape[1]):
                    for row in range(splitColumnVariableLevelLayout[index].shape[0]):
                        if splitColumnVariableLevelLayout[index][row,col] == columnVariable:
                            currentColumnVariableCoordinates.append([row+(index*columnVariableOffset),col])
            allColumnVariableCoordinates.append(currentColumnVariableCoordinates)

    else:
        for columnVariable in range(1,experimentParameters['numColumnLevelValues']+1):
            currentColumnVariableCoordinates = []
            for col in range(columnVariableLevelLayout.shape[1]):
                for row in range(columnVariableLevelLayout.shape[0]):
                    if columnVariableLevelLayout[row,col] == columnVariable:
                        currentColumnVariableCoordinates.append([row,col])
            allColumnVariableCoordinates.append(currentColumnVariableCoordinates)
    
    return allColumnVariableCoordinates

def createBaseDataFrame(experimentParameters,folderName,experimentNumber,dataType,allColumnVariableCoordinates,experimentLevelLayoutDict):
    if dataType == 'singlecell':
        realDataType = 'singlecell'
        dataType = 'cell'
    else:
        realDataType = dataType

    dim = experimentParameters['overallPlateDimensions']

    plateNames = []
    if experimentParameters['paired']:
        for i in range(0,int(experimentParameters['numPlates']/2)):
            plateNames.append('A'+str(i+1))
            plateNames.append('B'+str(i+1))
    else:
        for i in range(0,experimentParameters['numPlates']):
            plateNames.append('A'+str(i+1))
    
    sortedData,sortedFiles = cleanUpFlowjoCSV(plateNames,folderName,dataType)
    allRawData,newLevelList = returnMultiIndex(sortedData,sortedFiles,realDataType,folderName)
    fullExperimentData = []
    for levelIndex in range(len(newLevelList)):
        dataTypeLevel = newLevelList[levelIndex]
        experimentDataList = []
        numSamples = len(list(allRawData[0].iloc[:,levelIndex+1]))
        #Adds in dummy values (only used for reshaping)
        if numSamples != dim[0]*dim[1]:
            dummyVal = 3.14159
            dummyMatrix = np.ones(dim)*dummyVal
            valueMatrix = experimentLevelLayoutDict[experimentParameters['conditionLevelNames'][0]]
            contiguous = False
            if experimentParameters['paired']:
                for plateIndex in range(0,experimentParameters['numPlates'],2):
                    plateA = dummyMatrix.copy() 
                    plateB = dummyMatrix.copy()
                    index = 0
                    for row in range(valueMatrix.shape[0]):
                        for col in range(valueMatrix.shape[1]):
                            if valueMatrix[row,col] != 0 and index < numSamples:
                                plateA[row,col] = list(allRawData[plateIndex].iloc[:,levelIndex+1])[index]
                                plateB[row,col] = list(allRawData[plateIndex+1].iloc[:,levelIndex+1])[index]
                                index+=1
                    fullPlate = np.vstack([plateA,plateB])
                    experimentDataList.append(fullPlate)
            else:
                for plateIndex in range(experimentParameters['numPlates']):
                    plateA = dummyMatrix.copy() 
                    index = 0
                    for row in range(valueMatrix.shape[0]):
                        for col in range(valueMatrix.shape[1]):
                            if valueMatrix[row,col] != 0 and index < numSamples:
                                plateA[row,col] = list(allRawData[plateIndex].iloc[:,levelIndex+1])[index]
                                index+=1
                    experimentDataList.append(plateA)
        else:
            contiguous = True
            if experimentParameters['paired']:
                for plateIndex in range(0,experimentParameters['numPlates'],2):
                    plateA = np.reshape(list(allRawData[plateIndex].iloc[:,levelIndex+1]),(dim[0],dim[1]))
                    plateB = np.reshape(list(allRawData[plateIndex+1].iloc[:,levelIndex+1]),(dim[0],dim[1]))
                    fullPlate = np.vstack([plateA,plateB])
                    experimentDataList.append(fullPlate)
            else:
                for plateIndex in range(experimentParameters['numPlates']):
                    plateA = np.reshape(list(allRawData[plateIndex].iloc[:,levelIndex+1]),(dim[0],dim[1]))
                    experimentDataList.append(plateA)

        experimentData = np.hstack(experimentDataList)
        newDataLayoutList = []
        for columnVariableCoordinateList in allColumnVariableCoordinates:
            for coordinate in columnVariableCoordinateList:
                newDataLayoutList.append(experimentData[coordinate[0],coordinate[1]])
        newDataLayout = np.reshape(newDataLayoutList,(len(allColumnVariableCoordinates[0]),experimentParameters['numColumnLevelValues']),order='F')
        
        newIndexLayout = []
        for currentIndex in range(len(experimentParameters['conditionLevelNames'])):
            currentConditionLevelName = experimentParameters['conditionLevelNames'][currentIndex]
            currentConditionLevelValues = experimentParameters['conditionLevelValues'][currentConditionLevelName]
            currentLevelLayout = experimentLevelLayoutDict[currentConditionLevelName]
            newLayoutList = []
            columnVariableCoordinateList = allColumnVariableCoordinates[0]
            for coordinate in columnVariableCoordinateList:
                newLayoutList.append(currentLevelLayout[coordinate[0],coordinate[1]])
            newLayout = np.reshape(newLayoutList,(len(columnVariableCoordinateList),1),order='F')
            newIndexLayout.append(newLayout)
        newIndex = np.hstack(newIndexLayout)
        indexingTupleList = []
        for row in range(newIndex.shape[0]):
            indexingTupleNumeric = newIndex[row,:]
            indexingTuple = dataTypeLevel.copy()
            for numericLevelValue,conditionLevelName in zip(indexingTupleNumeric,experimentParameters['conditionLevelNames']):
                conditionLevelValues = experimentParameters['conditionLevelValues'][conditionLevelName]
                levelValue = conditionLevelValues[int(numericLevelValue)-1]
                indexingTuple.append(levelValue)
            indexingTupleList.append(indexingTuple)
        newMultiIndex = pd.MultiIndex.from_tuples(indexingTupleList,names=dataTypeLevelNames[realDataType]+experimentParameters['conditionLevelNames'])
        newLevelDataFrame = pd.DataFrame(newDataLayout,index=newMultiIndex,columns=experimentParameters['columnLevelValues'])
        newLevelDataFrame.columns.name = experimentParameters['columnVariableName']
        fullExperimentData.append(newLevelDataFrame)
    fullExperimentDf = pd.concat(fullExperimentData)
    return fullExperimentDf

def convertDataFramesToExcel(folderName,secondPath,dataType,df,useModifiedDf):
    if useModifiedDf:
        modifiedString = '-modified'
    else:
        modifiedString = ''
    writer = pd.ExcelWriter('semiProcessedData/excelFile-'+folderName+'-'+dataType+modifiedString+'.xlsx')
    if dataType == 'cyt':
        dfg = pickle.load(open('semiProcessedData/cytokineGFIPickleFile-'+folderName+'.pkl','rb'))
        dfc = pickle.load(open('semiProcessedData/'+dataTypeDataFrameFileNames[dataType]+'-'+folderName+modifiedString+'.pkl','rb'))
        dfg.to_excel(writer,'GFI')
        dfc.to_excel(writer,'Concentration')
    else:
        for statistic in list(pd.unique(df.index.get_level_values('Statistic'))):
            statisticDf = df.xs(statistic,level='Statistic')
            statisticDf.to_excel(writer,statistic)
    writer.save()
    print(dataType[0].upper()+dataType[1:]+' Excel file Saved')

def saveFinalDataFrames(folderName,secondPath,experimentNumber,dataType,fullExperimentDf,excel_data):
    modifiedFullExperimentDf = returnModifiedDf(experimentNumber,fullExperimentDf,dataType,excel_data)
    with open('semiProcessedData/'+dataTypeDataFrameFileNames[dataType]+'-'+folderName+'.pkl', "wb") as f:
        pickle.dump(fullExperimentDf, f)
    with open('semiProcessedData/'+dataTypeDataFrameFileNames[dataType]+'-'+folderName+'-modified.pkl', "wb") as f:
        pickle.dump(modifiedFullExperimentDf, f)
    convertDataFramesToExcel(folderName,secondPath,dataType,fullExperimentDf,False)
    convertDataFramesToExcel(folderName,secondPath,dataType,modifiedFullExperimentDf,True)
    print(fullExperimentDf)
    print(modifiedFullExperimentDf)

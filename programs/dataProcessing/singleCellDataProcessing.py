#!/usr/bin/env python3 
import math,pickle,os,sys,fcsparser,json,time,glob
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.preprocessing import StandardScaler
import processProliferationData as proliferationProcesser
from modifyDataFrames import returnModifiedDf
from miscFunctions import reindexDataFrame

def createInitialSingleCellDataFrame(folderName,experimentNumber,fileNameDataFrame):
    with open('inputFiles/experimentParameters-'+folderName+'.json') as f:
        experimentParameters = json.load(f)
    numConditions = experimentParameters[0][0]
    numTimePoints = experimentParameters[0][1]
    parameterNames = []
    #Uses from product if condition labels not entered in manually, otherwise reads tuples for each condition level with all condition labels within it (e.g. (CD28,N4,100pM),(CD28,N4,10pM)...)
    if(experimentParameters[5]):
        multiIndexedObject = pd.MultiIndex.from_tuples(list(zip(*experimentParameters[2])), names=experimentParameters[1])
    else:
        multiIndexedObject = pd.MultiIndex.from_product(experimentParameters[2], names=experimentParameters[1])
    #Adds time units to column headers
    timepoints = experimentParameters[3]
    timePointNames = []
    for i in range(0,len(timepoints)):
        timePointNames.append(float(timepoints[i]))
    
    #Grabs a file from samples to read marker names off of
    for fileName in os.listdir('semiProcessedData/singleCellData/A1/'):
        if 'DS' not in fileName:
            cellType = fileName
    tempFilePath = 'semiProcessedData/singleCellData/A1/'+cellType+'/' 
    fileExtension = '.csv'
    print(os.getcwd())
    tempFileName = glob.glob(tempFilePath+'*'+fileExtension)[0]
    experimentalChannelDf = pd.read_csv(tempFileName, header=0)
    experimentalChannelNames = experimentalChannelDf.columns.tolist()
    experimentalMarkerNames = []
    #Creates column heasdings for all measured parameters
    for i in range(len(experimentalChannelNames)):
        if 'FSC-A' in experimentalChannelNames[i]:
            experimentalMarkerNames.append('Size')
        elif 'SSC-A' in experimentalChannelNames[i]:
            experimentalMarkerNames.append('Granularity')
        elif 'CTV' in experimentalChannelNames[i]:
            experimentalMarkerNames.append('TCell_Gate')
        elif 'CTFR' in experimentalChannelNames[i] or 'FarRed Cell Trace' in experimentalChannelNames[i]:
            experimentalMarkerNames.append('APC_Gate')
        else:
            #"Comp-APC-A :: CTFR"
            if('::' in experimentalChannelNames[i]):
                experimentalMarkerNames.append(experimentalChannelNames[i].split(' :: ')[1])
            else:
                experimentalMarkerNames.append(experimentalChannelNames[i])
    stackedFileFrame = fileNameDataFrame.stack()
    levelNames = list(stackedFileFrame.index.names)
    singleCellLevelNames = levelNames+['Event']
    for scalingType in ['channel','scale']:
        fullFileFrameTemp = stackedFileFrame.copy().to_frame('Temp')
        completeDataFrameList = []
        for row in range(stackedFileFrame.shape[0]):
            levelValues = fullFileFrameTemp.iloc[row].name
            cellType = levelValues[0]
            fileIndex = stackedFileFrame.iloc[row].rfind('/')
            beforeCellType = stackedFileFrame.iloc[row][:fileIndex]
            afterCellType = scalingType+'_'+stackedFileFrame.iloc[row][fileIndex+1:]+'_'+cellType+fileExtension
            fullFileName = beforeCellType+'/'+cellType+'/'+afterCellType
            
            fcsDf = pd.read_csv(fullFileName,header=0)
            eventNumber = fcsDf.shape[0]
            eventList = range(1,eventNumber+1)
            allLevelValues = []
            for event in eventList:
                allLevelValues.append(list(levelValues)+[event])
            newMultiIndex = pd.MultiIndex.from_tuples(allLevelValues,names=singleCellLevelNames)
            newDf = pd.DataFrame(fcsDf.values,index=newMultiIndex,columns=experimentalMarkerNames)
            completeDataFrameList.append(newDf)
            print([scalingType]+list(levelValues))
        completeDataFrame = pd.concat(completeDataFrameList)
        completeDataFrame.columns.name = 'Markers'
        print(completeDataFrame)
        with open('semiProcessedData/initialSingleCellDf-'+scalingType+'-'+folderName+'.pkl','wb') as f:
            pickle.dump(completeDataFrame,f)
        if scalingType == 'channel':
            logicleDataProliferation = completeDataFrame['TCell_Gate'].xs(['TCells'],level=['CellType'])
            with open('semiProcessedData/logicleProliferationDf.pkl','wb') as f:
                pickle.dump(logicleDataProliferation,f)
        else:
            rawDataProliferation = completeDataFrame['TCell_Gate'].xs(['TCells'],level=['CellType'])
            with open('semiProcessedData/rawProliferationDf.pkl','wb') as f:
                pickle.dump(rawDataProliferation,f)

def createProliferationSingleCellDataFrame(folderName,secondPath,experimentNumber,useModifiedDf):
    logicleDataProliferation = pickle.load(open('semiProcessedData/logicleProliferationDf.pkl','rb'))
    rawDataProliferation = pickle.load(open('semiProcessedData/rawProliferationDf.pkl','rb'))
    proliferationProcesser.processProliferationData(folderName,logicleDataProliferation,rawDataProliferation)
    bulkprolifdf = proliferationProcesser.generateBulkProliferationStatistics(folderName,experimentNumber)
    return bulkprolifdf

def createCompleteSingleCellDf(folderName):
   
    idx = pd.IndexSlice
    
    initialSingleCellDf = pickle.load(open('semiProcessedData/initialSingleCellDf-channel-'+folderName+'.pkl','rb'))
    initialSingleCellDf = initialSingleCellDf.drop(['APC_Gate','TCell_Gate'],axis=1)
    bulkCytokineConcentrationDf = pickle.load(open('semiProcessedData/cytokineConcentrationPickleFile-'+folderName+'-modified.pkl','rb'))
    proliferationDf = pickle.load(open('semiProcessedData/singleCellDataFrame-proliferation-'+folderName+'.pkl','rb'))

    bulkCellStatisticDf = pickle.load(open('semiProcessedData/cellStatisticPickleFile-'+folderName+'-modified.pkl','rb'))

    cytokineLevelNamesUM = []
    cytokineLevelNames = []
    cellLevelNames = []

    tempCellLevels = list(bulkCellStatisticDf.iloc[0,:].name)[:3]
    tempCellDf = bulkCellStatisticDf.loc[tempCellLevels[0]].loc[tempCellLevels[1]].loc[tempCellLevels[2]]

    tempCytokineLevels = list(bulkCytokineConcentrationDf.iloc[0,:].name)[:1]
    tempCytokineDf = bulkCytokineConcentrationDf.loc[tempCytokineLevels[0]]
    
    unmodifiedCytokineConcentrationDf = pickle.load(open('semiProcessedData/cytokineConcentrationPickleFile-'+folderName+'.pkl','rb')).loc[tempCytokineLevels[0]]
    
    for row in range(tempCellDf.shape[0]):
        levelNames = list(tempCellDf.iloc[row,:].name)
        for column in tempCellDf.columns:
            cellLevelNames.append(tuple(levelNames+[column]))
    for row in range(tempCytokineDf.shape[0]):
        levelNames = list(tempCytokineDf.iloc[row,:].name)
        for column in tempCytokineDf.columns:
            cytokineLevelNames.append(tuple(levelNames+[column]))
    for row in range(unmodifiedCytokineConcentrationDf.shape[0]):
        levelNames = list(unmodifiedCytokineConcentrationDf.iloc[row,:].name)
        for column in unmodifiedCytokineConcentrationDf.columns:
            cytokineLevelNamesUM.append(tuple(levelNames+[column]))
    
    differences = list((set(tuple(cytokineLevelNamesUM)) | set(tuple(cytokineLevelNames)) | set(tuple(cellLevelNames))) - (set(tuple(cytokineLevelNamesUM)) & set(tuple(cytokineLevelNames)) & set(tuple(cellLevelNames))))
    levelsToKeep = [0]*initialSingleCellDf.shape[0]
    k=0
    row = 0
    while row < initialSingleCellDf.shape[0]:
        levelNames = tuple(initialSingleCellDf.iloc[row,:].name)[:-1]
        stackedLevelNames = tuple(list(levelNames)+[slice(None)])
        stackedLength = initialSingleCellDf.loc[stackedLevelNames,:].shape[0]
        #Check logic here later; should do for now
        if (tuple(levelNames) in differences):
            print(levelNames)
            pass
        else:
            rowVals = range(row,row+stackedLength)
            levelsToKeep[k:k+stackedLength] = rowVals
            k+=stackedLength
        row+=stackedLength
    initialSingleCellDf = initialSingleCellDf.iloc[levelsToKeep[:k],:]
    proliferationDf = proliferationDf.iloc[levelsToKeep[:k],:]
    indexList = []
    numEventsList = []
    for elem in pd.unique(initialSingleCellDf.index):
        indexList.append(elem[:-1])
    indexList = pd.unique(indexList)
    for index in indexList:
        numEventsList.append(initialSingleCellDf.loc[idx[index],:].shape[0])
    completeSingleCellCytokineValues = []
    for cytokine in pd.unique(bulkCytokineConcentrationDf.index.get_level_values(0)):
        individualCytokineSingleCellValues = []
        for index,numEvents in zip(indexList,numEventsList):
            bulkIndex = tuple([cytokine]+list(index[:-1]))
            bulkCytokineValue = bulkCytokineConcentrationDf.loc[idx[bulkIndex],index[-1]]
            singleCellCytokineValues = np.repeat(bulkCytokineValue,numEvents)
            individualCytokineSingleCellValues.append(singleCellCytokineValues)
        completeSingleCellCytokineValues.append(np.concatenate(individualCytokineSingleCellValues))
    singleCellCytokineMatrix = np.stack(completeSingleCellCytokineValues,axis=1)
    singleCellCytokineDf = pd.DataFrame(singleCellCytokineMatrix,index=initialSingleCellDf.index,columns=pd.unique(bulkCytokineConcentrationDf.index.get_level_values(0)))
     
    completeSingleCellDf = pd.concat([initialSingleCellDf,singleCellCytokineDf,proliferationDf],keys=['Markers','Cytokines','Proliferation'],names=['DataType','Parameter'],axis=1)
    print(completeSingleCellDf)
    with open('semiProcessedData/singleCellDataFrame-complete-'+folderName+'.pkl','wb') as f:
        pickle.dump(completeSingleCellDf,f)

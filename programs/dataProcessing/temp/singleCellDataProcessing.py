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

def createInitialSingleCellDataFrame(folderName,experimentNumber):
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
    tempFilePath = 'semiProcessedData/singleCellData/A1/TCells/' 
    fileExtension = '.csv'
    print(os.getcwd())
    tempFileName = glob.glob(tempFilePath+'*'+fileExtension)[0]
    experimentalChannelDf = pd.read_csv(tempFileName, header=0)
    experimentalChannelNames = experimentalChannelDf.columns.tolist()
    experimentalMarkerNames = []
    #Creates column heasings for all measured parameters
    for i in range(len(experimentalChannelNames)):
        if 'FSC-A' in experimentalChannelNames[i]:
            experimentalMarkerNames.append('Size')
        elif 'SSC-A' in experimentalChannelNames[i]:
            experimentalMarkerNames.append('Granularity')
        elif 'CTV' in experimentalChannelNames[i]:
            experimentalMarkerNames.append('TCell_Gate')
        elif 'CTFR' in experimentalChannelNames[i]:
            experimentalMarkerNames.append('APC_Gate')
        else:
            #"Comp-APC-A :: CTFR"
            if('::' in experimentalChannelNames[i]):
                experimentalMarkerNames.append(experimentalChannelNames[i].split(' :: ')[1])
            else:
                experimentalMarkerNames.append(experimentalChannelNames[i])
    
    emptyDf = np.zeros([numConditions,numTimePoints])
    fullDataFrame = pd.DataFrame(emptyDf,index=multiIndexedObject,columns=timePointNames)
    fullDataFrame.columns.name = 'Time'
    #fullDataFrame = returnModifiedDf(experimentNumber,fullDataFrame,'cell') 
    stackedDf = fullDataFrame.stack()
    levelValues = stackedDf.index.get_values()
    levelNames = list(stackedDf.index.names)

    #Grabs plate names (A1,A2 etc.)
    plateNames = experimentParameters[4]
    os.chdir('semiProcessedData/singleCellData/')
    for scalingType in ['channel','scale']:
        completeDataFrameList = []
        index2=0
        for plateName in plateNames:
            os.chdir(plateName)
            os.chdir('TCells')
            allFiles = glob.glob(os.getcwd()+'/'+scalingType+'*'+fileExtension)
            #Create empty list, fill with names of each file in order of collection (eg sample 001 goes in index 0, sample 005 goes in index 5 etc)
            fileNums = []
            for fileName in allFiles:
                if (fileName.find(fileExtension)>-1):
                    i = int(fileName[fileName.rfind('_0')+1:fileName.rfind('_TCells.')]) - 1
                    fileNums.append(i)
            firstFileList = [None]*(max(fileNums)+1)
            for fileName in allFiles:
                if (fileName.find(fileExtension)>-1):
                    i = int(fileName[fileName.rfind('_0')+1:fileName.rfind('_TCells.')]) - 1
                    firstFileList[i] = fileName
            fileList = []
            for fname in firstFileList:
                if(fname != None):
                    fileList.append(fname)
            bigList = []
            for fileName,index in zip(fileList,range(len(fileList))):
                fcsDf = pd.read_csv(fileName,header=0)
                eventNumber = fcsDf.shape[0]
                eventList = range(1,eventNumber+1)
                allLevelValues = []
                for event in eventList:
                    temp = list(levelValues[index2+index])
                    temp.append(event)
                    allLevelValues.append(temp)
                newNames = levelNames.copy()
                newNames.append('Event')
                newMultiIndex = pd.MultiIndex.from_tuples(allLevelValues,names=newNames)
                if fcsDf.shape[1] == 12 and experimentNumber == 87:
                    fcsDf = fcsDf.drop('Comp-APC-A :: CTFR',axis=1)
                newDf = pd.DataFrame(fcsDf.values,index=newMultiIndex,columns=experimentalMarkerNames)
                bigList.append(newDf)
            plateDf = pd.concat(bigList)
            completeDataFrameList.append(plateDf)
            os.chdir('../..')
            index2+=len(fileList)
        completeDf = pd.concat(completeDataFrameList)
        completeDf = returnModifiedDf(experimentNumber,completeDf,'sc')
        with open('../../semiProcessedData/initialSingleCellDf-'+scalingType+'-'+folderName+'.pkl','wb') as f:
            pickle.dump(completeDf,f)
    os.chdir('../..')

def createProliferationSingleCellDataFrame(folderName,secondPath,useModifiedDf):
    logicleDataFull = pickle.load(open('semiProcessedData/initialSingleCellDf-channel-'+folderName+'.pkl','rb'))['TCell_Gate']
    rawDataFull = pickle.load(open('semiProcessedData/initialSingleCellDf-scale-'+folderName+'.pkl','rb'))['TCell_Gate']
    logicleData,rawData = proliferationProcesser.createProliferationData(folderName,logicleDataFull,rawDataFull)
    ggl,tv,tl = proliferationProcesser.returnGatesAndTicks(logicleData,rawData)
    proliferationProcesser.processProliferationData(folderName,logicleData,rawData,ggl,tv,tl,logicleDataFull.index)
    bulkprolifdf = proliferationProcesser.generateBulkProliferationStatistics(folderName)
    return bulkprolifdf

def createCytokineSingleCellDataFrame(folderName):
    idx = pd.IndexSlice
    initialSingleCellDf = pickle.load(open('semiProcessedData/initialSingleCellDf-channel-'+folderName+'.pkl','rb'))
    bulkCytokineConcentrationDf = pickle.load(open('semiProcessedData/cytokineConcentrationPickleFile-'+folderName+'-modified.pkl','rb'))
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
    with open('semiProcessedData/singleCellDataFrame-cytokines-'+folderName+'.pkl', 'wb') as f:
        pickle.dump(singleCellCytokineDf,f)

def createCellSingleCellDf(folderName):
    idx = pd.IndexSlice
    initialSingleCellDf = pickle.load(open('semiProcessedData/initialSingleCellDf-channel-'+folderName+'.pkl','rb'))
    actualMarkers = []
    for marker in pd.unique(initialSingleCellDf.columns.get_level_values(0)):
        if 'APC_' not in marker:
            actualMarkers.append(marker)
    singleCellMarkerDf = initialSingleCellDf.loc[idx[:],idx[tuple(actualMarkers)]]
    with open('semiProcessedData/singleCellDataFrame-markers-'+folderName+'.pkl', 'wb') as f:
        pickle.dump(singleCellMarkerDf,f)

def createCompleteSingleCellDf(folderName):
   markerdf = pickle.load(open('semiProcessedData/singleCellDataFrame-markers-'+folderName+'.pkl','rb'))
   cytokinedf = pickle.load(open('semiProcessedData/singleCellDataFrame-cytokines-'+folderName+'.pkl','rb'))
   #proliferationdf = pickle.load(open('semiProcessedData/singleCellDataFrame-proliferation-'+folderName+'.pkl','rb') )
   completeSingleCellDf = pd.concat([markerdf,cytokinedf],keys=['Markers','Cytokines'],names=['DataType','Parameter'],axis=1)
   print(completeSingleCellDf)
   with open('semiProcessedData/singleCellDataFrame-complete-'+folderName+'.pkl','wb') as f:
       pickle.dump(completeSingleCellDf,f)
   with open('../../output/singleCellDataFrames/singleCellDataFrame-complete-'+folderName+'.pkl','wb') as f:
       pickle.dump(completeSingleCellDf,f)


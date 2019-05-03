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
    cellType = os.listdir('semiProcessedData/singleCellData/A1/')[0]
    tempFilePath = 'semiProcessedData/singleCellData/A1/'+cellType+'/' 
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
    plateLength = 12
    plateWidth = 8
    #If paired
    if(experimentParameters[6]):
        platesPerCondition = int(len(experimentParameters[4])/2)
    else:
        if experimentParameters[8]:
            numPlates = len(experimentParameters[4])
            platesPerTimePoint = math.ceil(numTimePoints/plateLength)
            platesPerCondition = int(int(numPlates/2)/platesPerTimePoint)
        else:
            platesPerCondition = int(len(experimentParameters[4]))
    if not experimentParameters[6] and experimentParameters[8]:
        emptyDf = np.zeros([numConditions,numTimePoints])
        fullDataFrame = pd.DataFrame(emptyDf,index=multiIndexedObject,columns=timePointNames)
        fullDataFrame.columns.name = 'Time'
        levelValues = []
        alltimepoints = fullDataFrame.columns
        newtimepointlist = []
        for alternate in range(platesPerTimePoint):
            alternatingTimepoints = range(alternate,platesPerTimePoint*plateLength,platesPerTimePoint)
            for alternatetimepoint in alternatingTimepoints:
                newtimepointlist.append(alltimepoints[alternatetimepoint])
        for row in range(fullDataFrame.shape[0]):
            conditionLabel = fullDataFrame.iloc[row,:].name
            for col in range(fullDataFrame.shape[1]):
                timeLabel = newtimepointlist[col] 
                levelValues.append(list(conditionLabel)+[timeLabel])
    else:
        timeSplit = int(numTimePoints/platesPerCondition)
        emptyDf = np.zeros([numConditions,numTimePoints])
        fullDataFrame = pd.DataFrame(emptyDf,index=multiIndexedObject,columns=timePointNames)
        fullDataFrame.columns.name = 'Time'
        levelValues = []
        for timeSplitIndex in range(platesPerCondition):
            for row in range(fullDataFrame.shape[0]):
                conditionLabel = fullDataFrame.iloc[row,:].name
                for col in range(timeSplitIndex*timeSplit,(timeSplitIndex+1)*timeSplit):
                    timeLabel = fullDataFrame.iloc[:,col].name
                    levelValues.append(list(conditionLabel)+[timeLabel])
    #fullDataFrame = returnModifiedDf(experimentNumber,fullDataFrame,'cell') 
    stackedDf = fullDataFrame.stack()
    levelNames = list(stackedDf.index.names)
    #Grabs plate names (A1,A2 etc.)
    plateNames = experimentParameters[4]
    os.chdir('semiProcessedData/singleCellData/')
    print(plateNames)
    for scalingType in ['channel','scale']:
        completeDataFrameList = []
        index2=0
        overallTimePoint=0
        print(scalingType)
        fileDict = {}
        for plateName in plateNames:
            print(plateName)
            os.chdir(plateName)
            os.chdir(cellType)
            allFiles = glob.glob(os.getcwd()+'/'+scalingType+'*'+fileExtension)
            #Create empty list, fill with names of each file in order of collection (eg sample 001 goes in index 0, sample 005 goes in index 5 etc)
            fileNums = []
            for fileName in allFiles:
                if (fileName.find(fileExtension)>-1):
                    i = int(fileName[fileName.rfind('_0')+1:fileName.rfind('_'+cellType+'.')]) - 1
                    fileNums.append(i)
            firstFileList = [None]*(max(fileNums)+1)
            for fileName in allFiles:
                if (fileName.find(fileExtension)>-1):
                    i = int(fileName[fileName.rfind('_0')+1:fileName.rfind('_'+cellType+'.')]) - 1
                    firstFileList[i] = fileName
            fileList = []
            for fname in firstFileList:
                if(fname != None):
                    fileList.append(fname)
            bigList = []
            if experimentParameters[6]:
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
                    newDf = pd.DataFrame(fcsDf.values,index=newMultiIndex,columns=experimentalMarkerNames)
                    bigList.append(newDf)
                plateDf = pd.concat(bigList)
                completeDataFrameList.append(plateDf)
                os.chdir('../..')
                index2+=len(fileList)
            else:
                numTimePointsPerPlate = int(numTimePoints/len(plateNames))
                fileIndices = np.reshape(np.reshape(np.arange(96),(8,12)).flatten('F'),(numConditions,numTimePointsPerPlate),order='F')
                for plateTimePoint in range(numTimePointsPerPlate):
                    conditionNumber = 0
                    for fileIndex in fileIndices[:,plateTimePoint]:
                        fileName = fileList[fileIndex]
                        trueLevelValIndex = (numTimePoints*conditionNumber)+overallTimePoint
                        trueLevelVal = levelValues[trueLevelValIndex]
                        stringVals = []
                        for val in trueLevelVal:
                            stringVals.append(str(val))
                        stringTrueLevelVal = '-'.join(stringVals)
                        fileDict[stringTrueLevelVal] = fileName
                        conditionNumber+=1
                    overallTimePoint+=1
                os.chdir('../..')
        if not experimentParameters[6]:
            for levelValue in levelValues:
                stringVals = []
                for val in levelValue:
                    stringVals.append(str(val))
                stringLevelValue = '-'.join(stringVals)
                fileName = fileDict[stringLevelValue]
                fcsDf = pd.read_csv(fileName,header=0)
                eventNumber = fcsDf.shape[0]
                eventList = range(1,eventNumber+1)
                allLevelValues = []
                for event in eventList:
                    allLevelValues.append(list(levelValue)+[event])
                newNames = levelNames.copy()
                newNames.append('Event')
                newMultiIndex = pd.MultiIndex.from_tuples(allLevelValues,names=newNames)
                newDf = pd.DataFrame(fcsDf.values,index=newMultiIndex,columns=experimentalMarkerNames)
                completeDataFrameList.append(newDf)
        completeDf = pd.concat(completeDataFrameList)    
        completeDf = reindexDataFrame(completeDf,fullDataFrame,True)
        completeDf.columns.name = 'Markers'
        print(completeDf)
        with open('../../semiProcessedData/initialSingleCellDf-'+scalingType+'-'+folderName+'.pkl','wb') as f:
            pickle.dump(completeDf,f)
    os.chdir('../..')

def createProliferationSingleCellDataFrame(folderName,secondPath,experimentNumber,useModifiedDf):
    logicleDataFull = pickle.load(open('semiProcessedData/initialSingleCellDf-channel-'+folderName+'.pkl','rb'))['TCell_Gate']
    rawDataFull = pickle.load(open('semiProcessedData/initialSingleCellDf-scale-'+folderName+'.pkl','rb'))['TCell_Gate']
    logicleData,rawData = proliferationProcesser.createProliferationData(folderName,logicleDataFull,rawDataFull)
    ggl,tv,tl = proliferationProcesser.returnGatesAndTicks(logicleData,rawData)
    proliferationProcesser.processProliferationData(folderName,logicleData,rawData,ggl,tv,tl,logicleDataFull)
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

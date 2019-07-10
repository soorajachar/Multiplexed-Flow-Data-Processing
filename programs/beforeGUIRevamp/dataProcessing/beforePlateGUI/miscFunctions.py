#!/usr/bin/env python3  
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import os,sys,pickle,math,re
import numpy as np
from scipy.optimize import curve_fit
from logicle import logicle,quantile 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import string 

listOfCytokines=['IFNg','IL-2','IL-4','IL-6','IL-10','IL-17A','TNFa']
#sampleIDOrder = True
sampleIDOrder = False 

def Hill(x, Amplitude, EC50, hill,Background):
    return np.log10(Amplitude * np.power(x,hill)/(np.power(EC50,hill)+np.power(x,hill))+Background)

#def boundedExponential(x, amplitude,tau):
#    return amplitude*(np.subtract(1,np.exp(np.multiply(-1,np.divide(x,tau)))))-4

#2 parameter (vshift fixed per cytokine based on lower LOD of cytokine): y = A(1-e^(-tau*x))
def boundedExponential(x, amplitude,tau,vshift):
    return amplitude*(np.subtract(1,np.exp(np.multiply(-1,np.multiply(x,tau)))))+vshift

#5 parameter (vshift fixed per cytokine; based on lower LOD of cytokine): y = A((1/(1+e^(-tau1*(x-td1))))-(1/(1+e^(-tau2*(x-td2)))))
def logisticDoubleExponential(x,amplitude,tau1,tau2,timedelay1,timedelay2,vshift):
    return amplitude*np.subtract(np.divide(1,np.add(1,np.exp(np.multiply(-1*tau1, np.subtract(x,timedelay1))))),np.divide(1,np.add(1,np.exp(np.multiply(-1*tau2,np.subtract(x,timedelay2))))))+vshift

def InverseHill(y,parameters):
    Amplitude=parameters[0]
    EC50=parameters[1]
    hill=parameters[2]
    Background=parameters[3]
    return np.power((np.power(10,y)-Background)/(Amplitude-np.power(10,y)),1/hill)*EC50

def r_squared(xdata,ydata,func,popt):
    residuals = ydata- func(xdata, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

def cleanUpFileList(fileArray):
    #Samples will be indexed based on well ID (A01, then A02 etc.)
    if not sampleIDOrder:
        orderWellID = {}
        index = 1
        for plateRow in plateRowList:
            for plateColumn in plateColumnList:
                orderWellID[str(plateRow)+str(plateColumn)] = index
                index+=1
        sortedFileList = [None]*len(fileArray)
        for name in fileArray:
            wellID = name.split('_')[3]
            sortedFileList[orderWellID[wellID]] = name
    #Samples will be indexed based on order of acqusition (sample001 then sample002 etc.)
    else:
        sortedFileList = [None]*len(fileArray)
        for name in fileArray:
            sampleID = int(name.split('_')[5])
            sortedFileList[sampleID-1] = name
    
    return sortedFileList

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
    if not sampleIDOrder:
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
    #Samples will be indexed based on order of acqusition (sample001 then sample002 etc.)
    else:
        for name in fileArray:
            temp = pd.read_csv('semiProcessedData/'+str(name)+'_'+dataTypeForCSV+'.csv')
            temp2 = [] 
            for i in range(0,temp.shape[0]):
                fullfilename = 'semiProcessedData/singleCellData/'+name+'/'+temp.iloc[i,0][:temp.iloc[i,0].find('.')]
                temp.iloc[i,0] = temp.iloc[i,0][temp.iloc[i,0].find('.')-3:temp.iloc[i,0].find('.')]
                temp2.append([temp.iloc[i,0],fullfilename])
            temp2 = pd.DataFrame(np.matrix(temp2[:-2]),columns=['Unnamed: 0','fileName'])
            temp = temp.sort_values('Unnamed: 0')
            temp2 = temp2.sort_values('Unnamed: 0')
            sortedData.append(temp[:-2])
            sortedFiles.append(temp2)
    if(dataType == 'cyt'):
        newMultiIndex = parseCytokineCSVHeaders(pd.read_csv('semiProcessedData/A1_'+dataType+'.csv').columns)
    elif(dataType == 'cell'):
        panelData = pd.read_csv('inputFiles/antibodyPanel-'+folderName+'.csv',)
        newMultiIndex = parseCellCSVHeaders(pd.read_csv('semiProcessedData/A1_'+dataType+'.csv').columns,panelData)
    elif(dataType == 'singlecell'):
        #Grabs a file from samples to read marker names off of
        cellTypeList = []
        for fileName in os.listdir('semiProcessedData/singleCellData/A1/'):
            if 'DS' not in fileName:
                cellTypeList.append(fileName)
        newMultiIndex = produceSingleCellHeaders(cellTypeList)
    elif(dataType == 'cytcorr'):
        newMultiIndex = []

    if dataType != 'singlecell':
        return sortedData,newMultiIndex
    else:
        return sortedFiles,newMultiIndex

def produceSingleCellHeaders(cellTypes):
    newMultiIndexList = []
    for cellType in cellTypes:
        newMultiIndexList.append([cellType])
    return newMultiIndexList

def parseCytokineCSVHeaders(columns):
    #,Beads/IFNg | Geometric Mean (YG586-A),Beads/IL-2 | Geometric Mean (YG586-A),Beads/IL-4 | Geometric Mean (YG586-A),Beads/IL-6 | Geometric Mean (YG586-A),Beads/IL-10 | Geometric Mean (YG586-A),Beads/IL-17A | Geometric Mean (YG586-A),Beads/TNFa | Geometric Mean (YG586-A),
    newMultiIndexList = []
    for column in columns[1:-1]:
        populationNameVsStatisticSplit = column.split(' | ')
        cytokine = populationNameVsStatisticSplit[0].split('/')[-1]
        newMultiIndexList.append([cytokine])
    return newMultiIndexList

def parseCellCSVHeaders(columns,panelData):
    #,Cells/Single Cells/APCs | Geometric Mean (Comp-BV605-A),Cells/Single Cells/APCs | Geometric Mean (Comp-FITC-A),Cells/Single Cells/APCs | Geometric Mean (Comp-PE ( 561 )-A),Cells/Single Cells/APCs | Count,Cells/Single Cells/APCs/CD86+ | Freq. of Parent (%),Cells/Single Cells/APCs/H2Kb+ | Freq. of Parent (%),Cells/Single Cells/APCs/PDL1+ | Freq. of Parent (%),Cells/Single Cells/TCells | Count,Cells/Single Cells/TCells | Geometric Mean (Comp-BUV737-A),Cells/Single Cells/TCells | Geometric Mean (Comp-BV605-A),Cells/Single Cells/TCells | Geometric Mean (Comp-PE-CF594-A),Cells/Single Cells/TCells | Geometric Mean (Comp-PE-Cy7-A),Cells/Single Cells/TCells | Geometric Mean (Comp-PerCP-Cy5-5-A),Cells/Single Cells/TCells/CD27+ | Freq. of Parent (%),Cells/Single Cells/TCells/CD54+ | Freq. of Parent (%),Cells/Single Cells/TCells/CD69+ | Freq. of Parent (%),Cells/Single Cells/TCells/PDL1+ | Freq. of Parent (%),
    newMultiIndexList = []
    for column in columns[1:-1]:
        populationNameVsStatisticSplit = column.split(' | ')
        fullPopulationName = populationNameVsStatisticSplit[0]
        #Statistics can be performed on the whole cell population, in which case the cellType is allEvents
        #GFI and CV need to be specified in terms of a laser channel; count and percent positive do not
        if '/' in fullPopulationName:
            populationDivisionIndices = [i for i,c in enumerate(fullPopulationName) if c=='/']
            #GFI or CV
            if 'Geometric Mean' in populationNameVsStatisticSplit[1] or 'CV' in populationNameVsStatisticSplit[1]:
                cellType = fullPopulationName
                #cellType = fullPopulationName[populationDivisionIndices[-1]+1:]
                if 'Comp-' in populationNameVsStatisticSplit[1]:
                    statisticVsChannelSplit = populationNameVsStatisticSplit[1].split(' (Comp-')
                else:
                    statisticVsChannelSplit = populationNameVsStatisticSplit[1].split(' (')
                statistic = statisticVsChannelSplit[0]
                if 'Geometric' in statistic:
                    statisticName = 'GFI'
                else:
                    statisticName = 'CV'
                
                #Can have IRF4+,TCells,CD45RB+CD25-; first case: cellType = first case:cellType=previousPop,marker=currentPop[:-1],stat=Positive/Negative GFI; 
                #second case: cellType = currentPop,marker=NotApplicable,statistic=GFI; third case: cellType = currentPop,marker=NotApplicable,statistic=GFI
                numPositive = fullPopulationName[populationDivisionIndices[-1]+1:].count('+')
                numNegative = fullPopulationName[populationDivisionIndices[-1]+1:].count('-')
                channel = statisticVsChannelSplit[1][:-1]
                panelIndex = list(panelData['FCSDetectorName']).index(channel)
                marker = panelData['Marker'][panelIndex]
                 
                #if numPositive+numNegative == 0, just use the last gate as the population name:
                if numPositive+numNegative == 0:
                    cellType = fullPopulationName[populationDivisionIndices[-1]+1:]
                    statistic = statisticName
                #If more than 1 +/- signs, include full population gating hiearchy as cell population name
                elif numPositive+numNegative > 1:
                    cellType = fullPopulationName
                    statistic = statisticName
                #If just a single + or - sign, extract marker name out of last gate
                else:
                    positivePop = fullPopulationName[populationDivisionIndices[-1]+1:-1]
                    if positivePop == marker:
                        cellType = fullPopulationName[:populationDivisionIndices[-1]]
                        print([cellType,marker,statistic])
                    else:
                        cellType = fullPopulationName
                    if numNegative == 0:
                        statistic = 'Positive '+statisticName
                    else:
                        statistic = 'Negative '+statisticName
            #% of parent and count
            else:
                numPositive = fullPopulationName[populationDivisionIndices[-1]+1:].count('+')
                numNegative = fullPopulationName[populationDivisionIndices[-1]+1:].count('-')
                if 'Freq' in populationNameVsStatisticSplit[1]:
                    if numPositive+numNegative == 0:
                        marker = 'NotApplicable'
                        statistic ='%'
                        cellType = fullPopulationName[populationDivisionIndices[-1]+1:]
                    elif numPositive+numNegative > 1:
                        marker = 'NotApplicable'
                        statistic = '%'
                        cellType = fullPopulationName
                    else:
                        if numNegative == 0:
                            statistic = '% Positive'
                        else:
                            statistic = '% Negative'
                        marker = fullPopulationName[populationDivisionIndices[-1]+1:len(fullPopulationName)-1]
                        positivePop = fullPopulationName[populationDivisionIndices[-1]+1:-1]
                        if positivePop == marker:
                            cellType = fullPopulationName[:populationDivisionIndices[-1]]
                            print([cellType,marker,statistic])
                        else:
                            cellType = fullPopulationName
                else:
                    marker = 'NotApplicable'
                    statistic = 'Count'
                    if numPositive+numNegative == 0:
                        cellType = fullPopulationName[populationDivisionIndices[-1]+1:] 
                    else:
                        cellType = fullPopulationName
        else:
            cellType = 'allEvents'
            #Statistics can be performed on the whole cell population, in which case the cellType is allEvents
            #DAPI+ | Freq. of Parent (%)
            #Positive cell percentage statistics do not have channel names, so treat differently
            #GFI or CV
            if len(populationNameVsStatisticSplit) == 1:
                populationNameVsStatisticSplit = [' ',populationNameVsStatisticSplit[0]]
            if 'Geometric Mean' in populationNameVsStatisticSplit[1] or 'CV' in populationNameVsStatisticSplit[1]:
                if 'Comp-' in populationNameVsStatisticSplit[1]:
                    statisticVsChannelSplit = populationNameVsStatisticSplit[1].split(' (Comp-')
                else:
                    statisticVsChannelSplit = populationNameVsStatisticSplit[1].split(' (')
                statistic = statisticVsChannelSplit[0]
                if 'Geometric' in statistic:
                    statistic = 'GFI'
                channel = statisticVsChannelSplit[1][:-1]
                panelIndex = list(panelData['FCSDetectorName']).index(channel)
                marker = panelData['Marker'][panelIndex]
            #% of parent and count
            else:
                marker = 'NotApplicable'
                if 'Freq' in populationNameVsStatisticSplit[1]:
                    statistic = '% Positive'
                else:
                    statistic = 'Count'
        newMultiIndexList.append([cellType,marker,statistic])
    commonBranchesIndices = []
    commonBranches = []
    for multiIndexList in newMultiIndexList:
        fullPopulationName = multiIndexList[0]
        populationDivisionIndices = [i for i,c in enumerate(fullPopulationName) if c=='/']
        commonBranches.append(fullPopulationName.split('/'))
        commonBranchesIndices.append(populationDivisionIndices)
    commonIndex = 0
    masterBranchList = []
    for branchlist in commonBranches:
        for branch in branchlist:
            masterBranchList.append(branch)
    uniqueBranches = list(pd.unique(masterBranchList))
    commonToAllStatistics = []
    for uniqueBranch in uniqueBranches:
        isCommonToAllStatistics = True
        for branchlist in commonBranches:
            if uniqueBranch not in branchlist:
                isCommonToAllStatistics = False
        if isCommonToAllStatistics:
            commonToAllStatistics.append(uniqueBranch)
    uncommonStatisticsList = []
    for commonBranch in commonBranches:
        uncommonStatistics = []
        for branch in commonBranch:
            if branch not in commonToAllStatistics:
                uncommonStatistics.append(branch)
        uncommonStatisticsList.append('/'.join(uncommonStatistics))
    for uncommonStatistic,i in zip(uncommonStatisticsList,range(len(uncommonStatisticsList))):
        newMultiIndexList[i] = [uncommonStatistic]+newMultiIndexList[i][1:]
    #for temp in newMultiIndexList:
    #    print(temp)
    return newMultiIndexList

def returnOrderedFiles(allFiles,extension):
    fileNums = []
    for fileName in allFiles:
        if (fileName.find(extension)>-1):
            i = int(fileName[fileName.rfind('_')-3:fileName.rfind('_')]) - 1
            fileNums.append(i)
    firstFileList = [None]*(max(fileNums)+1)
    for fileName in allFiles:
        if (fileName.find(extension)>-1):
            i = int(fileName[fileName.rfind('_')-3:fileName.rfind('_')]) - 1
            firstFileList[i] = fileName
    fileList = []
    for fname in firstFileList:
        if(fname != None):
            fileList.append(fname)
    return fileList

def cleanUpProliferationData(fileArray,folderName):
    for fileName in fileArray:
        
        path = 'semiProcessedData/CTV_Files/'+fileName+'/'
        fileList = returnOrderedFiles(os.listdir(path),'.csv')
        dataList = []
        for sample in fileList:
            sampleDataCSV = pd.read_csv(path+sample)
            raveledData = sampleDataCSV.values.ravel()
            dataList.append(raveledData)
        maxNumElements = 0
        for data in dataList:
            numElements = len(data)
            if(numElements > maxNumElements):
                maxNumElements = numElements
        dataMatrix = np.empty((len(dataList),maxNumElements,))
        dataMatrix[:] = np.nan
        i=0
        for data in dataList:
            dataMatrix[i,:len(data)] = data
            i+=1
        with open('test2.pkl', "wb") as f:
            pickle.dump(dataMatrix, f)
        sys.exit(0)
        
        """
        print(dataMatrix)
        T = 262143
        d = 4.5
        m = d*np.log(10)
        wVals = []
        for i in range(dataMatrix.shape[0]):
            logicleData0 = dataList[i]
            logicleData = logicleData0[logicleData0<T]
            if(len(logicleData[logicleData<0]) > 0):
                r = quantile(logicleData[logicleData<0], 0.05)
                w = (m-np.log(T/abs(r)))/2
            else:
                w=0
            wVals.append(w)
            transformedX = logicle(logicleData,T,m,r)
            dataMatrix[i,:len(logicleData)] = transformedX
        with open('test.pkl', "wb") as f:
            pickle.dump(dataMatrix, f)
        print(os.getcwd())
        print('done')
        """
        binNum=256
        dataMatrix = pickle.load(open('/Volumes/Group05/Altan-Bonnet/Sooraj/experiments/20181128-PeptideComparison_OT1_Timeseries_10/test2.pkl','rb'))
        for val in range(dataMatrix.shape[0]):
            withNan = dataMatrix[val,:]
            noNan = withNan[~np.isnan(withNan)]
            plt.hist(noNan,color='blue',bins=binNum)
            #sns.distplot(noNan,color='blue',bins=binNum)
            #sns.kdeplot(noNan,shade='True',color='blue',bw=0.5)
            plt.show()    

#Grab numeric value of timepoints columns
def returnNumericTimePoints(dfc):
    timePointsString = list(dfc.columns.values)
    timePoints = []
    for currentTimePointsString in timePointsString:
        timePoints.append(float(re.findall(r'\d+', currentTimePointsString)[0]))
    return timePoints

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
    print(sortedListSI)
    print(numericList)
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
            for eventVal in range(1,1+stackedLength):
                reindexedLevelNames.append(list(levelNames)+[eventVal])
            row+=1
            k+=stackedLength
        reindexedMultiIndex = pd.MultiIndex.from_tuples(reindexedLevelNames,names=dfToReindex.index.names)
        reindexedDf = pd.DataFrame(reindexedDfMatrix,index=reindexedMultiIndex,columns=dfToReindex.columns)
    return reindexedDf

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

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


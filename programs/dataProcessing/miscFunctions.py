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

listOfCytokines=['IFNg','IL-2','IL-4','IL-6','IL-10','IL-17A','TNFa']

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

#Sample_A2_A02_002.fcs,
def cleanUpFlowjoCSV(fileArray,folderName,dataType):
    sampleIDOrder = False
    #Samples will be indexed based on well ID (A01, then A02 etc.)
    if not sampleIDOrder:
        orderWellID = {}
        plateColumnList = list(range(1,13))
        plateRowList = ['A','B','C','D','E','F','G','H']
        index = 1
        for plateRow in plateRowList:
            for plateColumn in plateColumnList:
                orderWellID[str(plateRow)+str(plateColumn)] = index
                index+=1
        orderWellID['Mean'] = len(orderWellID.keys())+1
        orderWellID['SD'] = len(orderWellID.keys())+2
        sortedData = []
        for name in fileArray:
            temp = pd.read_csv('semiProcessedData/'+str(name)+'_'+dataType+'.csv')
            for i in range(0,temp.shape[0]):
                if '_' in temp.iloc[i,0]:
                    wellID = temp.iloc[i,0].split('_')[2]
                else:
                    wellID = temp.iloc[i,0]
                print(wellID)
                temp.iloc[i,0] = orderWellID[wellID]
            temp = temp.sort_values('Unnamed: 0')
            sortedData.append(temp[:-2])
    #Samples will be indexed based on order of acqusition (sample001 then sample002 etc.)
    else:
        sortedData = []
        for name in fileArray:
            temp = pd.read_csv('semiProcessedData/'+str(name)+'_'+dataType+'.csv')
            for i in range(0,temp.shape[0]):
                temp.iloc[i,0] = temp.iloc[i,0][temp.iloc[i,0].find('.')-3:temp.iloc[i,0].find('.')]
            temp = temp.sort_values('Unnamed: 0')
            sortedData.append(temp[:-2])
    if(dataType == 'cyt'):
        newMultiIndex = parseCytokineCSVHeaders(pd.read_csv('semiProcessedData/A1_'+dataType+'.csv').columns)
    elif(dataType == 'cell'):
        panelData = pd.read_csv('inputFiles/antibodyPanel-'+folderName+'.csv',)
        newMultiIndex = parseCellCSVHeaders(pd.read_csv('semiProcessedData/A1_'+dataType+'.csv').columns,panelData)
    return sortedData,newMultiIndex

def parseCytokineCSVHeaders(columns):
    #,Beads/IFNg | Geometric Mean (YG586-A),Beads/IL-2 | Geometric Mean (YG586-A),Beads/IL-4 | Geometric Mean (YG586-A),Beads/IL-6 | Geometric Mean (YG586-A),Beads/IL-10 | Geometric Mean (YG586-A),Beads/IL-17A | Geometric Mean (YG586-A),Beads/TNFa | Geometric Mean (YG586-A),
    newMultiIndexList = []
    for column in columns[1:-1]:
        populationNameVsStatisticSplit = column.split(' | ')
        cytokine = populationNameVsStatisticSplit[0].split('/')[1]
        newMultiIndexList.append([cytokine])
    return newMultiIndexList

def parseCellCSVHeaders(columns,panelData):
    #,Cells/Single Cells/APCs | Geometric Mean (Comp-BV605-A),Cells/Single Cells/APCs | Geometric Mean (Comp-FITC-A),Cells/Single Cells/APCs | Geometric Mean (Comp-PE ( 561 )-A),Cells/Single Cells/APCs | Count,Cells/Single Cells/APCs/CD86+ | Freq. of Parent (%),Cells/Single Cells/APCs/H2Kb+ | Freq. of Parent (%),Cells/Single Cells/APCs/PDL1+ | Freq. of Parent (%),Cells/Single Cells/TCells | Count,Cells/Single Cells/TCells | Geometric Mean (Comp-BUV737-A),Cells/Single Cells/TCells | Geometric Mean (Comp-BV605-A),Cells/Single Cells/TCells | Geometric Mean (Comp-PE-CF594-A),Cells/Single Cells/TCells | Geometric Mean (Comp-PE-Cy7-A),Cells/Single Cells/TCells | Geometric Mean (Comp-PerCP-Cy5-5-A),Cells/Single Cells/TCells/CD27+ | Freq. of Parent (%),Cells/Single Cells/TCells/CD54+ | Freq. of Parent (%),Cells/Single Cells/TCells/CD69+ | Freq. of Parent (%),Cells/Single Cells/TCells/PDL1+ | Freq. of Parent (%),
    newMultiIndexList = []
    for column in columns[1:-1]:
        populationNameVsStatisticSplit = column.split(' | ')
        fullPopulationName = populationNameVsStatisticSplit[0]
        populationDivisionIndices = [i for i,c in enumerate(fullPopulationName) if c=='/']
        #Positive cell percentage statistics do not have channel names, so treat differently
        if('Freq.' in populationNameVsStatisticSplit[1]):
            cellType = fullPopulationName[populationDivisionIndices[-2]+1:populationDivisionIndices[-1]]
            marker = fullPopulationName[populationDivisionIndices[-1]+1:len(fullPopulationName)-1]
            statistic = '% Positive'
        elif('Count' in populationNameVsStatisticSplit[1]):
            cellType = fullPopulationName[populationDivisionIndices[-1]+1:]
            #cellType = fullPopulationName[populationDivisionIndices[-2]+1:populationDivisionIndices[-1]]
            marker = 'NotApplicable'
            statistic = populationNameVsStatisticSplit[1]
        else:
            #GFI of positive populations vs overall TCell GFI
            if('+' in populationNameVsStatisticSplit[0]):
                cellType = fullPopulationName[populationDivisionIndices[-2]+1:populationDivisionIndices[-1]]
                marker = fullPopulationName[populationDivisionIndices[-1]+1:len(fullPopulationName)-1]
                statistic = 'Positive GFI'
            else:
                cellType = fullPopulationName[populationDivisionIndices[-1]+1:]
                statisticVsChannelSplit = populationNameVsStatisticSplit[1].split(' (Comp-')
                if('Geometric' in statisticVsChannelSplit[0]):
                    statistic = 'GFI'
                else:
                    statistic = statisticVsChannelSplit[0]
                channel = statisticVsChannelSplit[1][:-1]
                panelIndex = list(panelData['FCSDetectorName']).index(channel)
                marker = panelData['Marker'][panelIndex]

        newMultiIndexList.append([cellType,marker,statistic])
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
unitPrefixDictionary = {'fM':1e-15,'pM':1e-12,'nM':1e-9,'uM':1e-6,'mM':1e-3,'M':1e0,'':0}
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

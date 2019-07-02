#!/usr/bin/env python3
# coding: utf-8
# # Preprocess data Sooraj
# First, we import dependencies.
from sklearn import preprocessing
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from string import ascii_uppercase
from pprint import pprint
import sys
import itertools
from modifyDataFrames import originalIndexOrder

def preprocessDataNoTimeComponent(secondPath,experimentNumber,df,modelName):
    
    # Consider concentrations in logscale, then normalize (mean to 0, std = 1)
    pivotedDf = np.log10(df.stack().unstack(0))
    if('normPerObservable' in modelName):
        for i in range(pivotedDf.columns.size):
            pivotedDf.iloc[:,i] = (pivotedDf.iloc[:,i] - pivotedDf.iloc[:,i].mean()) / pivotedDf.iloc[:,i].std()
        preprocessString = '-normPerObservable'
    else:
        X_mean, X_std = X_log.mean(), X_log.std()
        pivotedDf.iloc[:,:] = (X_log - X_mean)/X_std
        preprocessString = ''
    with open(secondPath+"preprocessedDataFrames/preprocessedDataFrame%s-%d.pkl"%(preprocessString,experimentNumber), 'wb') as f:
        pickle.dump([pivotedDf.values,pivotedDf.index,pivotedDf.columns],f)

def preprocessDataTimePointPartitioning(secondPath,experimentNumber,df,modelName):
    upperLimit = 500
    lowerLimit = 0
    partitions = [(lowerLimit,24),(24,upperLimit)]
    newColumnTuples = []
    for timePoint in map(float,df.columns):
        for timePointPartition,partitionName in zip(partitions[::-1],ascii_uppercase):
            if(timePoint < timePointPartition[1] and timePoint >= timePointPartition[0]):
                newColumnTuples.append((partitionName,timePoint))
    partitionedMultiIndex = pd.MultiIndex.from_tuples(newColumnTuples,names=('Partition',df.columns.name))
    df.columns = partitionedMultiIndex
    pivotedDf = np.log10(df.stack().stack().unstack(0))
    if('normPerObservable' in modelName):
        for i in range(pivotedDf.columns.size):
            pivotedDf.iloc[:,i] = (pivotedDf.iloc[:,i] - pivotedDf.iloc[:,i].mean()) / pivotedDf.iloc[:,i].std()
        preprocessString = '-normPerObservable'
    else:
        X_mean, X_std = pivotedDf.values.mean(), pivotedDf.values.std()
        pivotedDf.iloc[:,:] = (pivotedDf.values - X_mean)/X_std
        preprocessString = ''
    with open(secondPath+"preprocessedDataFrames/preprocessedDataFrame-partitioned%s-%d.pkl"%(preprocessString,experimentNumber), 'wb') as f:
        pickle.dump([pivotedDf.values,pivotedDf.index,pivotedDf.columns],f)

def preprocessDataExpFitParameterizing(secondPath,experimentNumber,modelName):
    idx = pd.IndexSlice
    fullTrainingSet = pickle.load(open(secondPath+'fitParameterPickleFiles/allParameterDataSets.pkl', "rb"))
    test = fullTrainingSet.loc[idx[:,['IFNg','IL-6','TNFa'],:,:],idx[:]]
    pivotedDf = test['boundedExponential'].unstack('Cytokine').loc[idx[experimentNumber,:,:],idx[['A','k0','v'],:]]
    for i in range(pivotedDf.columns.size):
        pivotedDf.iloc[:,i] = (pivotedDf.iloc[:,i] - pivotedDf.iloc[:,i].mean()) / pivotedDf.iloc[:,i].std()
    with open(secondPath+"preprocessedDataFrames/preprocessedDataFrame-parameterized-%d.pkl"%(experimentNumber), 'wb') as f:
        pickle.dump([pivotedDf.values,pivotedDf.index,pivotedDf.columns],f)

def preprocessDataTimePointPartitioningExpFitParameterizing(secondPath,experimentNumber,df,modelName):
    upperLimit = 500
    lowerLimit = 0
    partitions = [(lowerLimit,24),(24,upperLimit)]
    newColumnTuples = []
    for timePoint in map(float,df.columns):
        for timePointPartition,partitionName in zip(partitions[::-1],ascii_uppercase):
            if(timePoint < timePointPartition[1] and timePoint >= timePointPartition[0]):
                newColumnTuples.append((partitionName,timePoint))
    partitionedMultiIndex = pd.MultiIndex.from_tuples(newColumnTuples,names=('Partition',df.columns.name))
    df.columns = partitionedMultiIndex
    pivotedDf = np.log10(df.stack().stack().unstack(0))
    for i in range(pivotedDf.columns.size):
        pivotedDf.iloc[:,i] = (pivotedDf.iloc[:,i] - pivotedDf.iloc[:,i].mean()) / pivotedDf.iloc[:,i].std()
    logPivotedDf = pd.DataFrame(pivotedDf.values,index=pivotedDf.index,columns=pivotedDf.columns)
    idx = pd.IndexSlice
    rawDataSegment = logPivotedDf.xs('B',level='Partition').loc[idx[:,:,:],idx[['IFNg','IL-6','TNFa']]]
    steadyStateMeasurements = pd.unique(rawDataSegment.index.get_level_values('Time'))
    
    parameterizedObservables = [['IFNg','IL-6','TNFa'],['IL-2','IL-17A']]
    fits = ['boundedExponential','logisticDoubleExponential']
    usedParameters = [['k0'],['A','k0','k1','x0','x1']]
    flattenedUsedParameters = list(itertools.chain(*usedParameters))
    dfListOfTuples = []
    featureCategories = ['Raw Data']+fits
    featureNamesList = [['LogConc']]+usedParameters
    observableList = [parameterizedObservables[0],parameterizedObservables[0],parameterizedObservables[1]]
    for i in range(len(featureCategories)):
        featureCategory = featureCategories[i]
        featureNames = featureNamesList[i]
        observableNames = observableList[i]
        for featureName in featureNames:
            for observableName in observableNames:
                dfListOfTuples.append((featureCategory,featureName,observableName))
    columnIndex = pd.MultiIndex.from_tuples(dfListOfTuples,names=['FeatureCategory','FeatureName','Observable'])
    fitParameterDf = pickle.load(open(secondPath+'fitParameterPickleFiles/fitParameterPickleFile-%d.pkl'%(experimentNumber), "rb"))
    for observable in pd.unique(fitParameterDf.index.get_level_values(0)):
        observableFitParameter = fitParameterDf.loc[observable]
        for i in range(fitParameterDf.columns.size):
            observableFitParameter.iloc[:,i] = (observableFitParameter.iloc[:,i] - observableFitParameter.iloc[:,i].mean())/observableFitParameter.iloc[:,i].std()
        fitParameterDf.loc[observable] = observableFitParameter.values
    columns = list(rawDataSegment.columns)+flattenedUsedParameters
    pivotedParameterizedDf = pd.DataFrame(np.zeros((rawDataSegment.index.size,columnIndex.size)),index=rawDataSegment.index,columns=columnIndex)
    for parameterObservableList,fit,parameterList in zip(parameterizedObservables,fits,usedParameters):
        for observable in parameterObservableList:
            observableFitParameterDf = fitParameterDf.loc[observable]
            for row in range(observableFitParameterDf.shape[0]):
                currentName = observableFitParameterDf.iloc[row,:].name
                for time in steadyStateMeasurements:
                    temp = list(currentName).copy()
                    temp.append(time)
                    newTuple = tuple(temp)
                    pivotedParameterizedDf.loc[idx[newTuple],idx[fit,:,observable]] = observableFitParameterDf.loc[idx[currentName],idx[fit,tuple(parameterList)]].values
    for observable in rawDataSegment.columns:
        pivotedParameterizedDf.loc[idx[:],idx[featureCategories[0],featureNamesList[0][0],observable]] = rawDataSegment.loc[idx[:],idx[observable]]
    print(pivotedParameterizedDf)
    with open(secondPath+"preprocessedDataFrames/preprocessedDataFrame-partitioned-parameterized-%d.pkl"%(experimentNumber), 'wb') as f:
        pickle.dump([pivotedParameterizedDf.values,pivotedParameterizedDf.index,pivotedParameterizedDf.columns],f)

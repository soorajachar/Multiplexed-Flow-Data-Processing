#!/usr/bin/env python3 
import json,pickle,math,matplotlib,sys,os,string
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
idx = pd.IndexSlice
import itertools
import os
from operator import itemgetter

def produceSubsettedDataFrames2(fulldf,withinFigureBoolean,specificValueBooleanList):
    
    #Get all possible subsetted indices
    figureSubsettedLevelValues = []
    withinFigureSubsettedLevelValues = []
    figureSubsettingLevels = []
    figureLevelNames = []
    figureLevelIndices = []
    levelValuesPlottedIndividually = []
    for levelIndex,currentLevelName in enumerate(fulldf.index.names):
        #Event level always has every value included (single cell)
        #Otherwise, check which level values in the level were selected by the user
        if currentLevelName != 'Event':
            currentLevelValues = pd.unique(fulldf.index.get_level_values(currentLevelName))
            levelValues = []
            #Go through each level value in the level; regardless of figure selection status, and add based on level value selection status
            for levelValue,specificBoolean in zip(currentLevelValues,specificValueBooleanList[levelIndex]):
                if specificBoolean:
                    levelValues.append(levelValue)
            #If we will include this level within the figure
            if withinFigureBoolean[levelIndex]:
                withinFigureSubsettedLevelValues.append(levelValues)
                if len(levelValues) == len(currentLevelValues):
                    for levelValue in levelValues:
                        levelValuesPlottedIndividually.append(levelValue)
            #Only need to add level values to figure list; will be xs'd out of the full dataframe in the subsetted, within figure dataframes
            else:
                figureSubsettedLevelValues.append(levelValues)
                figureLevelNames.append(currentLevelName)
                figureLevelIndices.append(levelIndex)
    
    #Get all row level values present in the dataframe
    rowList = []
    for row in range(fulldf.shape[0]):
        allCurrentLevelValues = fulldf.iloc[row,:].name
        currentLevelValues = itemgetter(*figureLevelIndices)(allCurrentLevelValues)
        if not isinstance(currentLevelValues,tuple):
            currentLevelValues = tuple([currentLevelValues])
        rowList.append(currentLevelValues)
    allPossibleSubsettingCombos = itertools.product(*figureSubsettedLevelValues)
    actualSubsettingCombos = []
    #From original dataframe; select all rows that appear in the all possible combination list 
    for subsettingCombo in allPossibleSubsettingCombos:
        if subsettingCombo in rowList:
            actualSubsettingCombos.append(subsettingCombo) 
    #Use these levels to cross section the fulldf, generating a list of subsetted dfs that will each have their own figure
    allPossibleSubsettedDfList = []
    for actualSubsettingCombo in actualSubsettingCombos:
        possibleSubsettedDf = fulldf.xs(actualSubsettingCombo, level=figureLevelNames)
        allPossibleSubsettedDfList.append(possibleSubsettedDf)
    actualLevelValueDfList = []
    #Go through each subsetteddf, and only grab rows with level values that are selected
    for possibleSubsettedDf in allPossibleSubsettedDfList:
        allPossibleLevelValueCombos = itertools.product(*withinFigureSubsettedLevelValues)
        rowList = []
        for row in range(possibleSubsettedDf.shape[0]):
            if 'Event' in fulldf.index.names:
                allCurrentLevelValues = possibleSubsettedDf.iloc[row,:].name[:-1]
            else:
                allCurrentLevelValues = possibleSubsettedDf.iloc[row,:].name
            rowList.append(allCurrentLevelValues)
        actualLevelValueCombos = []
        #From original dataframe; select all rows that appear in the all possible level value combination list 
        levelValueRowList = []
        for levelValueCombo in allPossibleLevelValueCombos:
            if levelValueCombo in rowList:
                indices = [i for i, x in enumerate(rowList) if x == levelValueCombo]
                levelValueRowList+=indices
        actualLevelValueDf = possibleSubsettedDf.iloc[levelValueRowList,:]
        actualLevelValueDfList.append(actualLevelValueDf)
    
    print(actualLevelValueDfList)
    print(actualSubsettingCombos)
    return actualLevelValueDfList,actualSubsettingCombos,figureLevelNames,levelValuesPlottedIndividually

#fulldf = pickle.load(open('../experiments/20190701-TCellNumber_CAR_Timeseries_1/semiProcessedData/cytokineConcentrationPickleFile-20190701-TCellNumber_CAR_Timeseries_1.pkl','rb'))
#withinFigureBoolean = [False, True, True, True]
#specificValueBooleanList = [[True, True, True, True, True, True, True], [True, True, True, True], [True, True, True, True], [True, True, True, True, True, True, True, True, True, True, True, True]]
df = pickle.load(open('../experiments/20190701-TCellNumber_CAR_Timeseries_1/postProcessedData/dimensionalReductions/dimensionalReduction-cyt-all.pkl','rb'))
plottingDf = pickle.load(open('../experiments/20190701-TCellNumber_CAR_Timeseries_1/temp.pkl','rb'))
plottingDf2 = df.reset_index()
print(plottingDf)
print(plottingDf2)
kwargs = {'row': 'DimensionalReductionMethod', 'row_order': ['TSNE', 'UMAP', 'ISOMAP'], 'hue': 'TCellNumber', 'col': 'TumorCellNumber', 'col_order': ['100K', '50K', '25K', '10K'], 'size': 'Time', 'x': 'Dimension 1', 'y': 'Dimension 2'}

facetKwargs = {'sharex': False, 'sharey': False} 
#g = sns.relplot(data=plottingDf,row = 'DimensionalReductionMethod', row_order = ['TSNE', 'UMAP', 'ISOMAP'], hue = 'TCellNumber', col = 'TumorCellNumber', col_order = ['100K', '50K', '25K', '10K'], size = 'Time', x = 'Dimension 1', y = 'Dimension 2',facet_kws={'sharex': False, 'sharey': False},)
fg = sns.relplot(data=plottingDf,marker='o',kind='scatter',facet_kws=facetKwargs,ci=False,**kwargs)
plt.show()
fg = sns.relplot(data=plottingDf2,marker='o',kind='scatter',facet_kws=facetKwargs,ci=False,**kwargs)
plt.show()
#a,b,c,d = produceSubsettedDataFrames2(fulldf.stack().to_frame('temp'),withinFigureBoolean,specificValueBooleanList)

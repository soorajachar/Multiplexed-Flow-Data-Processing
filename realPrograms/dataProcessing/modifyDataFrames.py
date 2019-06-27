#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import numpy as np
import pandas as pd
import sys,pickle

#Removes a particular level (and all of its sublevels) from the dataframe
def dropLevel(df,levelValueToDrop,levelOfLevelValueToDrop):
    newLevelValueList = pd.unique(df.index.get_level_values(levelOfLevelValueToDrop)).tolist()
    newLevelValueList.remove(levelValueToDrop)
    subsettingList = []
    for level in df.index.names:
        if level == levelOfLevelValueToDrop:
            subsettingList.append(newLevelValueList)
        else:
            subsettingList.append(slice(None))
    subsettedLevelValues = tuple(subsettingList)
    copiedDf = df.loc[subsettedLevelValues,:].copy()
    subsettedindex = copiedDf.index.remove_unused_levels()
    droppedLevelDf = pd.DataFrame(copiedDf,subsettedindex)
    return droppedLevelDf

#Structure of outlier indices: each element in the list has four parts: outlierRowIndexStart,outlierRowIndexEnd,outlierColumnIndexStart,outlierColumnIndexEnd
def averageNonEdgeOutliers(df,outlierIndices):
    for outlierIndex in outlierIndices:
        outlierRowIndexStart,outlierRowIndexEnd,outlierColumnIndexStart,outlierColumnIndexEnd = (outlierIndex[0],outlierIndex[1],outlierIndex[2],outlierIndex[3])
        for outlierConditionIndex in range(outlierRowIndexStart,outlierRowIndexEnd):
            for cyt in pd.unique(df.index.get_level_values(0)):
                interpolationStart = df.loc[cyt].iloc[outlierConditionIndex,outlierColumnIndexStart-1].copy()
                interpolationIncrement = (df.loc[cyt].iloc[outlierConditionIndex,outlierColumnIndexEnd]-interpolationStart)/((outlierColumnIndexEnd-outlierColumnIndexStart)+1)
                for i,outlierTimePointIndex in zip(range(1,outlierColumnIndexEnd-outlierColumnIndexStart+1),range(outlierColumnIndexStart,outlierColumnIndexEnd)):
                    newVal = interpolationStart+(i*interpolationIncrement)
                    df.loc[cyt].iloc[outlierConditionIndex,outlierTimePointIndex] = newVal
    return df

#Only if first or last timepoint is missing; set equal to previous timepoint or next time point (last and first respectively)
def averageEdgeOutliers(df,outlierIndices):
    for outlierIndex in outlierIndices:
        outlierRowIndexStart,outlierRowIndexEnd,outlierColumnIndexEnd = (outlierIndex[0],outlierIndex[1],outlierIndex[2])
        for outlierConditionIndex in range(outlierRowIndexStart,outlierRowIndexEnd):
            for cyt in pd.unique(df.index.get_level_values(0)):
                if(outlierColumnIndexEnd == df.shape[1]):
                    df.loc[cyt].iloc[outlierRowIndexStart:outlierRowIndexEnd+1,outlierColumnIndexEnd-1] = df.loc[cyt].iloc[outlierRowIndexStart:outlierRowIndexEnd+1,outlierColumnIndexEnd-2]
                else:
                    df.loc[cyt].iloc[outlierRowIndexStart:outlierRowIndexEnd+1,outlierColumnIndexEnd-1] = df.loc[cyt].iloc[outlierRowIndexStart:outlierRowIndexEnd+1,outlierColumnIndexEnd]
    return df

#Replace saturated cytokine measurements with diluted measurements
def replaceSaturatedCytokines(folderName,df,dilutedDf,dilutionFactor):
    LODParameters = pickle.load(open('semiProcessedData/LODParameters-'+folderName+'-nM.pkl', "rb"))
    for row in range(dilutedDf.shape[0]):
        dilutedDfNames = tuple(dilutedDf.iloc[row,:].name)
        upperConcLOD = LODParameters[dilutedDfNames[0]][3]
        dfRowToDesaturate = df.loc[dilutedDfNames,:]
        print(df.loc[dilutedDfNames])
        for i in range(dfRowToDesaturate.shape[0]):
            if df.loc[dilutedDfNames].values[i] == upperConcLOD:
                df.loc[dilutedDfNames].values[i] = dilutionFactor*dilutedDf.iloc[row,i]
        print(df.loc[dilutedDfNames])
    return df

#Perform various modifications of dataframe; specific to each experiment (outlier averaging, erroneous value dropping etc.)
def returnModifiedDf(experimentNumber,df,dataType):
    excel_data = pd.read_excel('../masterExperimentSpreadsheet.xlsx')
    folderName = excel_data['Full Name'][experimentNumber-1]
    onlyPairs = False
    #Drop positive control aCD3/aCD28, interpolate Q4,1nM T+57 and T4,1uM,T+86 (erroneous)
    if(experimentNumber == 28):
        if dataType == 'cyt':
            wrongDataIndices = [[6,7,11,12],[8,9,17,18]]
            df = dropLevel(df,'aCD3/aCD28','Peptide')
            df = averageNonEdgeOutliers(df,wrongDataIndices)
    #Drop positive control aCD3/aCD28 (concentration throws things off)
    elif(experimentNumber == 30):
        if dataType == 'cyt':
            df = df.mean(level=['Cytokine','Peptide','Concentration'])
            df = dropLevel(df,'aCD3/aCD28','Peptide')
    #Drop negative control Null (concentration value throws off things)
    elif(experimentNumber == 53):
        if dataType == 'cyt':
            df = dropLevel(df,'Null','Peptide')
    #Average the two results from the two mice, remove level from index
    elif(experimentNumber == 64):
        if dataType == 'cyt':
            df = df.mean(level=['Cytokine','TCellType','Peptide','Concentration'])
    #Drop positive control aCD3/aCD28, Take difference of timepoints 11 and 13, divide by 2, assign to conditions 7 and 8
    elif(experimentNumber == 67):
        if dataType == 'cyt':
            wrongDataIndices = [[6,8,10,12]]
            df = dropLevel(df,'aCD3/aCD28','Peptide')
            df = averageNonEdgeOutliers(df,wrongDataIndices)
    #Interpolate timepoint 6 and 21 for antigen/concentration measurements 1-8 and 9-16 respectively in IL-2 antibody condition (mistakenly added CD28 to those wells)
    #Interpolate timepoint 21 E1 1uM (erroneously high)
    #Drop +100hr timepoint (cba was incorrect)
    elif(experimentNumber == 68):
        if dataType == 'cyt':
            wrongDataIndices = [[32,40,5,6],[40,48,20,21],[15,16,20,21]]
            df = averageNonEdgeOutliers(df,wrongDataIndices)
            df = df.drop(['100'], axis=1)
    elif(experimentNumber == 71):
        if dataType == 'sc':
            A2B2_Timepoints = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 13.0]
            for tp in A2B2_Timepoints:
                df = dropLevel(df,tp,'Time')
    #Interpolate timepoints 2,4,6,8,10,12,14,16,18,20,22,24 for conditions 1-8 (plate A1 was dropped by robot so only every other timepoint could be taken)
    elif(experimentNumber == 72):
        if dataType == 'cyt':
            wrongDataIndices = [[0,8,1,2],[0,8,3,4],[0,8,5,6],[0,8,7,8],[0,8,9,10],[0,8,11,12],[0,8,13,14],[0,8,15,16],[0,8,17,18],[0,8,19,20],[0,8,21,22]]
            df = averageNonEdgeOutliers(df,wrongDataIndices)
            df = averageEdgeOutliers(df,[[0,8,24]])
    elif(experimentNumber == 69):
        if(onlyPairs):
            if dataType == 'cyt':
                idx = pd.IndexSlice
                df = df.loc[idx[:,:,['N4','Q4','Q7','V4','G4'],:],idx[:]]
                levelsToKeep = []
                for i in range(df.shape[0]):
                    rowlevels = df.iloc[i,:].name
                    if(rowlevels[-1] != '1uM' or rowlevels[-2] != 'N4'):
                        levelsToKeep.append(i)
                df = df.iloc[levelsToKeep,:]
            #Drop first 3 timepoints (outliers because methanol dried out)
            elif dataType == 'cell':
                df = df.iloc[:,3:]
                idx = pd.IndexSlice
                df = df.loc[idx[:,:,:,:,['N4','Q4','Q7','V4','G4'],:],idx[:]]
                levelsToKeep = []
                for i in range(df.shape[0]):
                    rowlevels = df.iloc[i,:].name
                    if(rowlevels[-1] != '1uM' or rowlevels[-2] != 'N4'):
                        levelsToKeep.append(i)
                df = df.iloc[levelsToKeep,:]
        else:
            pass
            #if dataType == 'sc':
            #    df = dropLevel(df,'n/a','Concentration')
    elif(experimentNumber == 78):
        if dataType == 'cyt':
            df = df.drop(index='Blank',level=1) 
            if(onlyPairs):
                idx = pd.IndexSlice
                df = df.loc[idx[:,:,['N4','Q4','T4','V4'],:],idx[:]]
                levelsToKeep = []
                for i in range(df.shape[0]):
                    rowlevels = df.iloc[i,:].name
                    if((rowlevels[-2] == 'N4' and rowlevels[-1] == '100pM') or (rowlevels[-2] == 'Q4' and rowlevels[-1] == '10nM') or (rowlevels[-2] == 'T4' and rowlevels[-1] == '1uM') or (rowlevels[-2] == 'V4' and rowlevels[-1] == '1uM')):
                        levelsToKeep.append(i)
                df = df.iloc[levelsToKeep,:]
        elif dataType == 'cell':
            df = df.drop(index='Blank',level='Genotype') 
            if(onlyPairs):
                idx = pd.IndexSlice
                #df = df.loc[idx[:,:,:,:,['N4','Q4','T4','V4'],:],idx[:]]
                levelsToKeep = []
                for i in range(df.shape[0]):
                    rowlevels = df.iloc[i,:].name
                    if((rowlevels[-2] == 'N4' and rowlevels[-1] == '100pM') or (rowlevels[-2] == 'Q4' and rowlevels[-1] == '10nM') or (rowlevels[-2] == 'T4' and rowlevels[-1] == '1uM') or (rowlevels[-2] == 'V4' and rowlevels[-1] == '1uM')):
                        levelsToKeep.append(i)
                df = df.iloc[levelsToKeep,:]
    elif(experimentNumber == 79):
        if dataType == 'cyt':
            df = df.drop(index='Blank',level=1)
    elif(experimentNumber == 80):
        if dataType == 'cyt':
            df = df.drop(index='Blank',level=1)
            df = df.drop(['46','40'],axis=1)
    elif(experimentNumber == 81):
        if dataType == 'cyt':
            df = df.drop(index='Blank',level=1)
            df = df.drop(['46','40'],axis=1)
    elif(experimentNumber == 82):
        if dataType == 'cyt':
            df = df.drop(index='Blank',level=1)
        elif dataType == 'cell':
            df = df.drop(index='Blank',level='TumorCellNumber')
    elif(experimentNumber == 83):
        print(df)
        if dataType == 'cyt':
            if(onlyPairs):
                df1 = df.iloc[:5,:]
                df2 = df.iloc[8:13,:]
                df3 = df.iloc[16:21,:]
                df4 = df.iloc[24:29,:]
                slicingdf = pd.concat([df1,df2,df3,df4])
                slicingTuples = []
                selectionIndices = []
                for row in range(slicingdf.shape[0]):
                    slicingTuples.append(slicingdf.iloc[row,:].name[-2:])
                for row in range(df.shape[0]):
                    currentrowPepConc = df.iloc[row,:].name[-2:]
                    if(currentrowPepConc in slicingTuples):
                        selectionIndices.append(row)
                df = df.iloc[selectionIndices,:]
        elif dataType == 'cell':
            if(onlyPairs):
                df1 = df.iloc[:5,:]
                df2 = df.iloc[8:13,:]
                df3 = df.iloc[16:21,:]
                df4 = df.iloc[24:29,:]
                slicingdf = pd.concat([df1,df2,df3,df4])
                slicingTuples = []
                selectionIndices = []
                for row in range(slicingdf.shape[0]):
                    slicingTuples.append(slicingdf.iloc[row,:].name[-2:])
                for row in range(df.shape[0]):
                    currentrowPepConc = df.iloc[row,:].name[-2:]
                    if(currentrowPepConc in slicingTuples):
                        selectionIndices.append(row)
                df = df.iloc[selectionIndices,:]
    elif(experimentNumber == 86):
        if dataType == 'cyt':
            df = pd.concat([df.iloc[:,:3],df.iloc[:,5:]],axis=1)
    elif(experimentNumber == 87):
        if dataType == 'cell':
            df = pd.concat([df.iloc[:,:-2],df.iloc[:,-1]],axis=1)
            print(df)
            #df = dropLevel(df,71.0,'Time')
        if dataType == 'sc':
            df = dropLevel(df,71.0,'Time')
    elif experimentNumber == 88:
        if dataType == 'cyt':
            dilutedDf = pickle.load(open('semiProcessedData/correctedCytokineDataFrames-'+folderName+'-Concentration.pkl','rb'))[folderName]
            df = replaceSaturatedCytokines(folderName,df,dilutedDf,10)
    #Mixed N4 10pM and Q4 1uM by accident. Last timepoint did not have enough cells
    elif(experimentNumber == 91):
        levelsToKeep = []
        for i in range(df.shape[0]):
            rowlevels = df.iloc[i,:].name
            if((rowlevels[-2] == 'N4' and rowlevels[-1] == '10pM') or (rowlevels[-2] == 'Q4' and rowlevels[-1] == '1uM')):
                pass
            else:
                levelsToKeep.append(i)
        df = df.iloc[levelsToKeep,:-1]
        if dataType == 'cyt':
            dilutedDf = pickle.load(open('semiProcessedData/correctedCytokineDataFrames-'+folderName+'-Concentration.pkl','rb'))[folderName]
            df = replaceSaturatedCytokines(folderName,df,dilutedDf,5)
        elif dataType == 'cell' or dataType == 'prolif':
            levelsToKeep2 = []
            for i in range(df.shape[0]):
                rowlevels = df.iloc[i,:].name
                if((rowlevels[-2] == 'N4' or rowlevels[-2] == 'Q4')):
                    pass
                else:
                    levelsToKeep2.append(i)
            df = df.iloc[levelsToKeep2,:]
    elif experimentNumber == 92:
        if dataType == 'cell' or dataType == 'prolif':
            levelsToKeep2 = []
            for i in range(df.shape[0]):
                rowlevels = df.iloc[i,:].name
                if((rowlevels[-2] == 'N4' or rowlevels[-2] == 'Q4')):
                    pass
                else:
                    levelsToKeep2.append(i)
            df = df.iloc[levelsToKeep2,:]
    elif experimentNumber == 93:
        if dataType == 'cyt':
            dilutedDf = pickle.load(open('semiProcessedData/correctedCytokineDataFrames-'+folderName+'-Concentration.pkl','rb'))[folderName]
            df = replaceSaturatedCytokines(folderName,df,dilutedDf,10)
    #Last 2-3 timepoints dried out due to low humidity
    elif experimentNumber == 96:
        if dataType == 'cyt':
            dilutedDf = pickle.load(open('semiProcessedData/correctedCytokineDataFrames-'+folderName+'-Concentration.pkl','rb'))[folderName]
            df = replaceSaturatedCytokines(folderName,df,dilutedDf,10)
        df = df.iloc[:,:-3]
    #Need to change original index order to real way of keep index order (with reset_index in single cell processing scripts)
    modifiedDf = df.copy()
    modifiedDf.columns.name = 'Time'
    print(modifiedDf)
    return modifiedDf

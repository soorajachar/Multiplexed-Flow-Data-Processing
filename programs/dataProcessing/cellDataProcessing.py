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
from miscFunctions import Hill,InverseHill,r_squared,cleanUpFlowjoCSV
from modifyDataFrames import returnModifiedDf

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



def createCellDataFrame(folderName,finalDataFrame,concUnitPrefix):
        
        return finalData

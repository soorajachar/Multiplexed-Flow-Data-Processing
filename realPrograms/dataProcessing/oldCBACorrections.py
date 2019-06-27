#!/usr/bin/env python3
from pathlib import Path
import pickle
from itertools import product,combinations
import numpy as np
import os,sys,argparse
import math
import json
import subprocess
import pandas as pd
import string
from initialDataProcessing import calibrateExperiment
from miscFunctions import Hill,InverseHill,r_squared,cleanUpFlowjoCSV
sys.path.insert(0, '../figuresPipeline/')
import facetPlottingLibrary as fpl
import matplotlib.pyplot as plt

listOfCytokines=['IFNg','IL-2','IL-4','IL-6','IL-10','IL-17A','TNFa']
pathToExperimentSpreadsheet = '../../experiments/'

#Create .json file with lists used to construct multiIndex; based entirely on user input. Json file is human readable so mistakes can be corrected directly in the file after the script writes it
def createCorrectedGFI(currentExpNum):
    
    upperCase = string.ascii_uppercase
    plateLength = 12
    plateWidth = 8
    
    numExperiments = int(input('Enter the total number of experiments you ran a cba correction for: '))
    expNums = []
    for exp in range(numExperiments):
        expNum = int(input('Enter the experiment number of experiment '+str(exp)+': '))
        expNums.append(expNum)
    
    names = []
    timepoints = []
    excel_data = pd.read_excel(pathToExperimentSpreadsheet+'masterExperimentSpreadsheet.xlsx')
    folderName = excel_data['Full Name'][currentExpNum-1]
    for expNum in expNums:
        names.append(excel_data['Full Name'][expNum-1])
        timepoints.append(excel_data['NumTimepoints'][expNum-1])
    numConditionsCorrectedList = []
    
    """
    fullLevelList = []
    for name,i in zip(names,range(len(names))):
        numConditionsCorrected = int(input('Enter the number of conditions corrected in ' + str(names[i])+': '))
        numConditionsCorrectedList.append(numConditionsCorrected)
        with open(pathToExperimentSpreadsheet+'/'+names[i]+'/inputFiles/experimentParameters-'+names[i]+'.json') as f:
            experimentParameters = json.load(f)
        levelNames = experimentParameters[1]
        conditionLevelList = []
        for condition in range(numConditionsCorrected):
            levelList = []
            for levelName in levelNames:
                level = input('Enter '+str(levelName)+' for condition ' + str(condition) + ' in experiment ' + names[i]+': ')
                levelList.append(level)
            conditionLevelList.append(levelList)
        fullLevelList.append(conditionLevelList)
    with open('fullLevelList.json','w') as f:
        json.dump(fullLevelList,f)
    """

    fullLevelList = json.load(open('fullLevelList.json','r'))
    os.chdir(pathToExperimentSpreadsheet+'/'+folderName)
    calibrateExperiment(folderName,'',1e9,'nM',excel_data['NumberOfCBAStandardDilutions'][currentExpNum-1],excel_data['CBAStandardDilutedVolume'][currentExpNum-1])
    dflist = []
    for levelList,i in zip(fullLevelList,range(len(fullLevelList))):
        plateLetter = upperCase[i]
        matrixList = []
        cytList = [None]*7
        for conditionList,j in zip(levelList,range(len(levelList))): 
            plateNumber = j+1
            plate = plateLetter+str(plateNumber)
            values,_ = cleanUpFlowjoCSV([plate],names[i],'cytcorr')
            values = values[0]
            for cytIndex in range(1,values.shape[1]-1):
                if j == 0:
                    cytList[cytIndex-1]=values.iloc[:,cytIndex]
                else:
                    cytList[cytIndex-1]=cytList[cytIndex-1].append(values.iloc[:,cytIndex])
        cytMatrixList = []
        for cytl in cytList:
            cytmatrix = np.reshape(cytl.values,(len(levelList),timepoints[i]))
            cytMatrixList.append(cytmatrix)
        fullMatrix = np.vstack(cytMatrixList)
        with open('../'+names[i]+'/inputFiles/experimentParameters-'+names[i]+'.json') as f:
            experimentParameters = json.load(f)
        levelNames = experimentParameters[1]
        columnNames = experimentParameters[3]
        indexList = []
        for cytokine in listOfCytokines:
            for level in fullLevelList[i]:
                indexList.append([cytokine]+level)
        multiIndex = pd.MultiIndex.from_tuples(indexList,names=['Cytokine']+levelNames)
        df = pd.DataFrame(fullMatrix,index=multiIndex,columns=columnNames)
        df.columns.name = 'Time'
        dflist.append(df)
    dfDict = {}
    for name,i in zip(names,range(len(names))):
        dfDict[name] = dflist[i]
    with open('../'+folderName+'/semiProcessedData/correctedCytokineDataFrames-'+folderName+'-GFI.pkl','wb') as f:
        pickle.dump(dfDict,f)

def calibrateDataFrame(folderName,secondPath,experimentNumber,concUnit,concUnitPrefix,dataType,finalDataFrame):
    fittingParameters = pickle.load(open('semiProcessedData/fittingParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "rb"))
    LODParameters = pickle.load(open('semiProcessedData/LODParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "rb"))
    #Begin converting GFI dataframe into corresponding concentration dataframe
    concentrationList = []
    #Step through dataframe one cytokine at a time
    for cytokine in pd.unique(finalDataFrame.index.get_level_values(0)):
        #Retrieve LODs for current cytokine (from constructed calibration curve)
        lowerGFILOD = LODParameters[cytokine][0]
        upperGFILOD = LODParameters[cytokine][1]
        lowerConcLOD = LODParameters[cytokine][2]
        upperConcLOD = LODParameters[cytokine][3]
        cyt = listOfCytokines.index(cytokine)
        smallConcentrationMatrix = np.zeros(finalDataFrame.loc[cytokine].shape)
        #Loop through every value in current cytokine's portion of the dataframe
        for i in range(0,finalDataFrame.loc[cytokine].values.shape[0]):
            for j in range(0,finalDataFrame.loc[cytokine].values.shape[1]):
                currentGFIval = finalDataFrame.loc[cytokine].values[i,j]
                if currentGFIval > upperGFILOD: #If intensity is greater than upper GFI LOD
                    currentConcVal = upperConcLOD #Concentration is equal to upper concentration LOD
                elif currentGFIval <= upperGFILOD and currentGFIval >= lowerGFILOD: #if intensity is between upper and lower GFI LODs
                    currentConcVal = InverseHill(np.log10(currentGFIval),fittingParameters[cyt,:]) #Use previous hill fit parameters for the cytokine to obtain concentration
                else: #If intensity is less than background GFI LOD
                    currentConcVal = lowerConcLOD #Concentration is equal to lower concentration LOD
                smallConcentrationMatrix[i,j] = currentConcVal
        concentrationList.append(smallConcentrationMatrix)

    concentrationMatrix = np.vstack(concentrationList)
    finalDataFrameConcentration = pd.DataFrame(concentrationMatrix,index=finalDataFrame.index,columns=finalDataFrame.columns)
    finalDataFrameConcentration.columns.name = 'Time'
    #Fill in n/a's (caused by undetectable fluorescence) with lower concentration LOD
    finalDataFrameConcentration.fillna(lowerConcLOD, inplace=True)
    return finalDataFrameConcentration

def createCorrectedConcentration(expNum,folderName):
    os.chdir(pathToExperimentSpreadsheet+'/'+folderName)
    dfGFIDict = pickle.load(open('semiProcessedData/correctedCytokineDataFrames-'+folderName+'-GFI.pkl','rb'))
    dfConcDict = {}
    for key in dfGFIDict:
        dfgfi = dfGFIDict[key]
        dfconc = calibrateDataFrame(folderName,'',expNum,1e9,'nM','cyt',dfgfi)
        if '0404' in key:
            dfconc.iloc[:,4:8]/=2
        dfConcDict[key] = dfconc
    with open('semiProcessedData/correctedCytokineDataFrames-'+folderName+'-Concentration.pkl','wb') as f:
        pickle.dump(dfConcDict,f)
    for key in dfConcDict:
        os.chdir('../'+key)
        dfconc = dfConcDict[key]
        with open('semiProcessedData/correctedCytokineDataFrames-'+key+'-Concentration.pkl','wb') as f:
            pickle.dump(dfConcDict,f)

def plotCorrectedConc(folderName):
    dfConcDict = pickle.load(open('semiProcessedData/correctedCytokineDataFrames-'+folderName+'-Concentration.pkl','rb'))
    plotType = 'ordered'
    dataType = 'cyt'
    subPlotType = 'line'
    plt.switch_backend('QT4Agg') #default on my system
    for key in dfConcDict:
        dfconc = dfConcDict[key]
        if '0404' in key:
            #dfconc.iloc[:,4:8]/=3
            fpl.facetPlottingGUI(dfconc,plotType,dataType)
            subsettedDfList,subsettedDfListTitles,levelsToPlot,figureLevels,levelValuesPlottedIndividually = fpl.produceSubsettedDataFrames(key,'',dfconc,False)
            fpl.plotFacetedFigures(key,plotType,subPlotType,dataType,subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually,False)

def main():
    parser = argparse.ArgumentParser(description="Create diluted cytokine concentration dataframes.")
    parser.add_argument("--input", dest='inputString', help ="Run specified dilution dataframes.")
    args = parser.parse_args()
    expNum = int(args.inputString)
    excel_data = pd.read_excel(pathToExperimentSpreadsheet+'masterExperimentSpreadsheet.xlsx')
    folderName = excel_data['Full Name'][expNum-1]
    #createCorrectedGFI(expNum)
    createCorrectedConcentration(expNum,folderName)
    #plotCorrectedConc(folderName)
    
main()

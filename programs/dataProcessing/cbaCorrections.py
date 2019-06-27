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
import seaborn as sns
import string
from initialDataProcessing import calibrateExperiment
from miscFunctions import Hill,InverseHill,r_squared,cleanUpFlowjoCSV
sys.path.insert(0, '../figuresPipeline/')
import facetPlotLibrary as fpl
import matplotlib.pyplot as plt

listOfCytokines=['IFNg','IL-2','IL-4','IL-6','IL-10','IL-17A','TNFa']
pathToExperimentSpreadsheet = '../../experiments/'
upperCase = string.ascii_uppercase

#Create .json file with lists used to construct multiIndex; based entirely on user input. Json file is human readable so mistakes can be corrected directly in the file after the script writes it
def createCBACorrectionParameters(cbaCorrExpNum,folderName,excel_data):
    
    plateLength = 12
    plateWidth = 8
    
    #Get all experiment numbers which had diluted cbas for this correction
    numExperiments = int(input('Enter the total number of experiments you ran a cba correction for: '))
    expNums = []
    dilutionFactors = []
    for exp in range(numExperiments):
        expNum = int(input('Enter the experiment number of experiment '+str(exp)+': '))
        expNums.append(expNum)
        dilutionFactor = int(input('Enter the dilution factor of experiment '+str(exp)+': '))
        dilutionFactors.append(dilutionFactor)

    #Get all experiment names which had diluted cbas for this correction
    names = []
    timepoints = []
    for expNum in expNums:
        names.append(excel_data['Full Name'][expNum-1])
        timepoints.append(excel_data['NumTimepoints'][expNum-1])
    
    numConditionsCorrectedList = []
    fullLevelList = []
    for name,i in zip(names,range(len(names))):
        numConditionsCorrected = int(input('Enter the number of conditions corrected in ' + str(names[i])+': '))
        numConditionsCorrectedList.append(numConditionsCorrected)
        with open('../'+names[i]+'/inputFiles/experimentParameters-'+names[i]+'.json') as f:
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
    with open('inputFiles/cbaCorrectionParameters-'+folderName+'.json','w') as f:
        json.dump([fullLevelList,names,timepoints,dilutionFactors],f)

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

#Need names,timepoints,excel_data, 
def produceCBADilutedData(cbaCorrExpNum,folderName,excel_data):
    fullLevelList = json.load(open('inputFiles/cbaCorrectionParameters-'+folderName+'.json','r'))[0]
    names = json.load(open('inputFiles/cbaCorrectionParameters-'+folderName+'.json','r'))[1]
    timepoints = json.load(open('inputFiles/cbaCorrectionParameters-'+folderName+'.json','r'))[2]
    dilutionFactors = json.load(open('inputFiles/cbaCorrectionParameters-'+folderName+'.json','r'))[3]
    calibrateExperiment(folderName,'',1e9,'nM',excel_data['NumberOfCBAStandardDilutions'][cbaCorrExpNum-1],excel_data['CBAStandardDilutedVolume'][cbaCorrExpNum-1])
    dflist = []
    for levelList,i in zip(fullLevelList,range(len(fullLevelList))):
        plateLetter = upperCase[i]
        matrixList = []
        cytList = [None]*7
        for conditionList,j in zip(levelList,range(len(levelList))): 
            plateNumber = j+1
            print(len(levelList))
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
    dfGFIDict = {}
    dilutionFactorDict = {}
    for name,i in zip(names,range(len(names))):
        dfGFIDict[name] = dflist[i]
        dilutionFactorDict[name] = dilutionFactors[i]

    dfConcDict = {}
    for key in dfGFIDict:
        dfgfi = dfGFIDict[key]
        dfconc = calibrateDataFrame(folderName,'',cbaCorrExpNum,1e9,'nM','cyt',dfgfi)
        if '0404' in key:
            dfconc.iloc[:,4:8]/=2
        dfConcDict[key] = dfconc
    with open('semiProcessedData/correctedCytokineDataFrames-'+folderName+'-Concentration.pkl','wb') as f:
        pickle.dump(dfConcDict,f)
    for key in dfConcDict:
        os.chdir('../'+key)
        dfconc = dfConcDict[key]
        dilutionFactor = dilutionFactorDict[key]
        with open('semiProcessedData/correctedCytokineDataFrameAndDilutionFactor-'+key+'.pkl','wb') as f:
            print('wat')
            pickle.dump([dfconc,dilutionFactor],f)

def plotCorrectedConc(folderName):
    dfConcDict = pickle.load(open('semiProcessedData/correctedCytokineDataFrames-'+folderName+'-Concentration.pkl','rb'))
    plotType = 'ordered'
    dataType = 'cyt'
    subPlotType = 'line'
    plt.switch_backend('QT4Agg') #default on my system
    for key in dfConcDict:
        specificExperimentFolder = '../'+key
        dataFrameAndDilutionFactor = pickle.load(open(specificExperimentFolder+'/semiProcessedData/correctedCytokineDataFrameAndDilutionFactor-'+key+'.pkl','rb'))
        correctedDf = dataFrameAndDilutionFactor[0]
        dilutionFactor = dataFrameAndDilutionFactor[1]
        uncorrectedDf = pickle.load(open(specificExperimentFolder+'/semiProcessedData/cytokineConcentrationPickleFile-'+key+'.pkl','rb'))
        uncorrectedSampleIndices = []
        for correctedSample in range(correctedDf.shape[0]):
            correctedrow = correctedDf.iloc[correctedSample,:]
            for row in range(uncorrectedDf.shape[0]):
                currentrow = uncorrectedDf.iloc[row,:]
                if currentrow.name[-2] == correctedrow.name[-2] and currentrow.name[-1] == correctedrow.name[-1]  and currentrow.name[0] == correctedrow.name[0]:
                    uncorrectedSampleIndices.append(row)
                    break
        uncorrectedSamples = uncorrectedDf.iloc[uncorrectedSampleIndices,:]
        newDf = pd.concat([uncorrectedSamples,correctedDf*dilutionFactor],axis=0,keys=['NotDiluted','Diluted'],names=['DilutionStatus'])
        print(newDf)
        idx=pd.IndexSlice
        plottingDf = newDf.loc[idx[:,['IFNg','IL-2']],:].stack().to_frame('CytokineConcentration').reset_index()
        ax = sns.relplot(data=plottingDf,row='Cytokine',hue='Concentration',col='DilutionStatus',style='Peptide',x='Time',y='CytokineConcentration',facet_kws={'sharey':False})
        k = len(ax.fig.get_axes())
        for i in range(k):
            ax.fig.get_axes()[i].set_yscale('log')
        plt.show()

#Want several things; to be able to create cytokine dilution experiment parameter files (should be specific to cytokine dilution location); 
def main():
    parser = argparse.ArgumentParser(description="Create diluted cytokine concentration dataframes.")
    parser.add_argument("-ep", action='store_true', help = "Create Cytokine Dilution Correction experiment parameters.")
    parser.add_argument("-pd", action='store_true', help = "Process Cytokine Dilution Correction.")
    parser.add_argument("-plt", action='store_true', help = "Plot Cytokine Dilution Corrections vs originals.")
    parser.add_argument("--input", dest='inputString', help ="Run specified dilution dataframes.")
    args = parser.parse_args()

    expNum = int(args.inputString)
    excel_data = pd.read_excel(pathToExperimentSpreadsheet+'/masterExperimentSpreadsheet.xlsx')
    folderName = excel_data['Full Name'][expNum-1]
    
    os.chdir(pathToExperimentSpreadsheet+'/'+folderName)

    if args.ep:
        createCBACorrectionParameters(expNum,folderName,excel_data)
    elif args.pd:
        produceCBADilutedData(expNum,folderName,excel_data)
    elif args.plt:
        plotCorrectedConc(folderName)
    else:
        sys.exit(0)
main()

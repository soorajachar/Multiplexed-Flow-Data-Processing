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
from miscFunctions import Hill,InverseHill,r_squared,cleanUpFlowjoCSV,cleanUpProliferationData
from modifyDataFrames import returnModifiedDf

#Conditions dimension
plateWidth = 8
#Timepoints dimension
plateLength = 12

#Standard BD Biosciences CBA Kit Cytokines
listOfCytokines=['IFNg','IL-2','IL-4','IL-6','IL-10','IL-17A','TNFa']
MWofCytokines=[17200,17200,14000,21900,18900,15500,17500] #g/mol, correct one

cytokineHeaderNames = ['Cytokine']
cellHeaderNames = ['CellType','Marker','Statistic']

#Create Calibration Curves, obtain LODs (limits of detection) of experiment
def calibrateExperiment(folderName,secondPath,concUnit,concUnitPrefix,numberOfCalibrationSamples,initialStandardVolume):
    #Get cytokine calibration curve data
    sortedData,newMultiIndexedObject = cleanUpFlowjoCSV(['Calibration'],folderName,'cyt')
    calibration = sortedData[0]
    data = np.array(calibration.values[:,1:8],dtype=float)
    
    fittingParameters = np.zeros((np.shape(listOfCytokines)[0],4))
    concLOD = np.zeros((np.shape(listOfCytokines)[0],2))
    
    #Initial concentration all cytokine standards is given by CBA kit manual as 5000 pGg/mL: when standards are diluted in 2mL
    conc = 5000 #pg/mL
    serialDilutionFactor = 2 #1:serialDilutionFactor dilution between each standard well
    #Smaller initial dilution (0.5mL instead of 2mL for example) increase the initial concentration of the first calibration sample
    initialConc = (conc*1e-12) /((initialStandardVolume*1e-3)/2) #g/L (pg/mL * 1e-12 g/pg)/(1e-3 L/mL)
    #Calibration samples are always diluted by a factor of serialdilutionFactor (so with 12 calibration samples, the last sample is (serialDilutionFactor^-11) the concentration of the first, which is pure standard (2^0)
    cbaStandardsConcentrations = np.flipud(initialConc*np.power(serialDilutionFactor,np.linspace(-numberOfCalibrationSamples+1,0,numberOfCalibrationSamples)))
    #More x values along the above concentration bounds are sampled to use to construct calibration curve. Plot points are extended slightly at high range to allow visualization of upper LOD (not accessible with experimental dilution)
    cbaStandardsConcentrationsPlotPoints = np.flipud(initialConc*np.power(2,np.linspace(-numberOfCalibrationSamples+1,4,101)))

    fig1=plt.figure(num=1,figsize=(10,10))
    plt.gcf().set_facecolor('white')
    color_list = plt.cm.jet(np.linspace(0,1,7))
    ax=fig1.add_subplot(1,1,1)

    concLOD = {}
    for i in range(len(listOfCytokines)):
        #amplitude bounded from range/2 to range*2, EC50 bounded from minimum to maximum standard concentration tested, Hill coefficient bounded from 0 to 2, Background bounded from 0 to minimum GFI*2
        lowerCurveFitBounds = [(np.max(data[:,i])-np.min(data[:,i]))/2,np.min(cbaStandardsConcentrations),0,0]
        upperCurveFitBounds = [(np.max(data[:,i])-np.min(data[:,i]))*2, np.max(cbaStandardsConcentrations), 2,np.min(data[:,i])*2]
        #use scipy curve fit to determine best hill equation fit for data, searching within the bounds given above
        popt,pcov = curve_fit(Hill, cbaStandardsConcentrations,np.log10(data[:,i]),sigma=np.log10(data[:,i]),bounds=(lowerCurveFitBounds,upperCurveFitBounds))
        rsquared = round(r_squared(cbaStandardsConcentrations,np.log10(data[:,i]),Hill,popt),3)

        for j in range(len(popt)):  
            #Convert just ec50 value to desired units (nM,uM etc)
            if j == 1:
                fittingParameters[i,j] = np.multiply(popt[j],(concUnit/MWofCytokines[i]))
            #other values in 4 parameter logistic equation are tied to intensity y-value, which doesn't change, or are the hill coefficient, which is completely separate, so parameters are kept the same
            else:
                fittingParameters[i,j] = popt[j]
        
        #Convert x values of experimental data points and curve fit points to desired units (nM,uM,etc.)
        convertedCBAStandards = np.multiply(cbaStandardsConcentrations,(concUnit/MWofCytokines[i]))
        convertedCBAStandardsPlotPoints = np.multiply(cbaStandardsConcentrationsPlotPoints,(concUnit/MWofCytokines[i]))
        #Plot on log-log scale the experimental points and the curve fit line with previously determined curve fitting parameters
        plt.loglog(convertedCBAStandards,data[:,i],'o',color=color_list[i,:],label=listOfCytokines[i])
        plt.loglog(convertedCBAStandardsPlotPoints,np.power(10,Hill(convertedCBAStandardsPlotPoints,*fittingParameters[i,:])),color=color_list[i,:],label=listOfCytokines[i]+'_fit; R2 = '+str(rsquared))
        
        #Get LOD for each cytokine calibration curve (aka the linear range of the calibration curve)
        backgroundGFI = fittingParameters[i,3]
        amplitudeGFI = fittingParameters[i,0]
        
        #Approximate LOD by determining concentration values at LOD% and 1-LOD% (3% and 97%) of curve. Must be used on log10(curve), as calibration curve is plotted in logscale
        LODpercent = 0.03
        #LOD% more than background GFI used for lower LOD GFI
        lowerGFILOD = math.log10(10**((1+LODpercent)*math.log10(backgroundGFI)))
        #LOD% less than maximum GFI (Background + amplitude) used for upper LOD GFI
        upperGFILOD = math.log10(10**((1-LODpercent)*math.log10(amplitudeGFI+backgroundGFI)))
        #Log10(upper/lowerGFILOD) converted back into normal GFI by 10 to its power, then fed into inverse hill equation with current cytokine fitting parameters to obtain corresponding concentration values
        lowerConcLOD = InverseHill(lowerGFILOD,fittingParameters[i,:])
        upperConcLOD = InverseHill(upperGFILOD,fittingParameters[i,:])
        listLOD = [10**lowerGFILOD,10**upperGFILOD,lowerConcLOD,upperConcLOD]
        #Create dict with keys as cytokines, values as GFI/conc LODs
        concLOD[listOfCytokines[i]] = listLOD
        #Plot vertical lines at lower and upper concentration limits of detection
        plt.axvline(x=lowerConcLOD,color=color_list[i,:],linestyle=':') 
        plt.axvline(x=upperConcLOD,color=color_list[i,:],linestyle=':') 

    plt.legend(loc=0,numpoints=1)           
    plt.ylabel('GeoMFI') #Geometric mean of fluorescence of each standard
    plt.xlabel('Concentration of Cytokine Standards ('+concUnitPrefix+')')
    plt.title('Calibration of CBA assay \n'+folderName,fontsize=14)
    plt.close()
    fig1.savefig('fullyProcessedFigures/calibrationCurves-'+folderName+'-'+concUnitPrefix+'.png')
    #Save fitting parameters and LOD for curve fit for each cytokine
    with open('semiProcessedData/fittingParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "wb") as f:
        pickle.dump(fittingParameters, f)
    with open('semiProcessedData/LODParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "wb") as f:
        pickle.dump(concLOD, f)

#Reshapes plateData list obtained directly from flowjo produced .csv file into a list of arrays, with each array containing all the conditions at a given timepoint on the plate. 
#PlateData is an n element list where n is the number of samples on a plate. 
def reshapePlate(plateData,numConditions,numTimePoints,columnsPerTimePoint,contiguous):
    
    #The dummy val is used to fill in "missing" conditions when conditions have blanks between them; e.g. a 12 condition single plate experiment (4 blanks between last condition and end of 2nd column, but first column is full)
    dummyVal = 3.14159
    individualSamples = len(plateData)
    #Condition numbers are always conserved across timepoints (not necessarily across plates), so the num samples on a plate divided by the num conditions on that plate will give the number of timepoints on the plate
    numTimePoints = int(individualSamples/numConditions)
    newNumConditions = numConditions
    if not contiguous:
        #Blank conditions in a non-contiguous plate will always be the number of conditions needed to make numConditions a multiple of 8 (plate WIdth)
        numMissingConditions = numConditions % plateWidth
        plateDataCopy = [1]
        #Create ntimepoints number of dummyvals, add this list of dummy vals to the original values; e.g. a 14x6 single plate (with 2 blank conditions) gets a 2x6 array of dummyvals added to it
        dummyVals = [dummyVal]*numTimePoints
        for j in range(numMissingConditions):
            plateDataCopy+=list(plateData[j*numTimePoints:(j+1)*numTimePoints])+dummyVals
        if(len(plateDataCopy) > 1):
            plateData = np.array(plateDataCopy[1:])
        #With dummy vals included; numConditions is updated to include the dummy conditions
        newNumConditions += numMissingConditions
    #Reshape the list of data (could be including dummyvals) into 
    plateLayout = np.reshape(list(plateData),(int(newNumConditions/columnsPerTimePoint),numTimePoints*columnsPerTimePoint))
    flattenedPlateLayout = list(plateLayout.flatten('F'))
    if dummyVal in flattenedPlateLayout:
        for val in flattenedPlateLayout.copy():
            if val == dummyVal:
                flattenedPlateLayout.remove(val)
    miniPlate = np.reshape(np.array(flattenedPlateLayout),(numConditions,numTimePoints),order='F')
    miniPlateList = np.vsplit(miniPlate,columnsPerTimePoint)
    return miniPlateList

def createAndCombineBaseDataFrames(folderName,allRawData,numPlates,numTimePoints,dataTypeLevelList,dataTypeLevelNames,paired,contiguous,replicateWise,alternatingPlatewise):
    
    #Get number of timepoints per plate, which is equal to the total length of the plate (12 columns) divided by the number of conditions per plate (numConditions/2). If numConditions is not a multiple of 16 (for example 14), round up to the highest multiple of 16 (math.ceil)
    totalSamplesA = 0
    totalSamplesB = 0
    baseDataFrameList = []

    if not alternatingPlatewise: 
        for i in range(0,int(numPlates/2)):
            if paired:
                plateIndexA = i*2
                plateIndexB = plateIndexA + 1
                totalSamplesB+=len(allRawData[plateIndexB].iloc[:,0])
            else:
                plateIndexA = i
            totalSamplesA+=len(allRawData[plateIndexA].iloc[:,0])
        
        numConditionsA = int(totalSamplesA/numTimePoints)
        if paired:
            numConditionsB = int(totalSamplesB/numTimePoints)
        columnsPerTimePoint = math.ceil(numConditionsA/plateWidth)

        #Iterates through each plate 2 at a time (robot processes 2 plates at a time: A/B plates; each A/B set combined has all the conditions for a given experiment)
        for levelIndex in range(len(dataTypeLevelList)):
            masterMiniPlateAList = []
            masterMiniPlateBList = []
            for i in range(0,int(numPlates/2)):
                if paired:
                    plateIndexA = i*2
                    plateIndexB = plateIndexA + 1
                    masterMiniPlateAList.append(reshapePlate(allRawData[plateIndexA].iloc[:,levelIndex+1],numConditionsA,numTimePoints,columnsPerTimePoint,contiguous))
                    masterMiniPlateBList.append(reshapePlate(allRawData[plateIndexB].iloc[:,levelIndex+1],numConditionsB,numTimePoints,columnsPerTimePoint,contiguous))
                else:
                    plateIndex = i
                    masterMiniPlateAList.append(reshapePlate(allRawData[plateIndex].iloc[:,levelIndex+1],numConditionsA,numTimePoints,columnsPerTimePoint,contiguous))

            finalPlateMatrix = []
            for i in range(len(masterMiniPlateAList)):
                if paired:
                    if(replicateWise):
                        finalPlateA = []
                        finalPlateB = []
                        for j in range(len(masterMiniPlateAList[0])):
                            finalPlateA.append(masterMiniPlateAList[i][j])
                        for j in range(len(masterMiniPlateBList[0])):
                            finalPlateB.append(masterMiniPlateBList[i][j])
                        finalPlateMatrix.append(np.vstack([np.vstack(finalPlateA),np.vstack(finalPlateB)]))
                    else:
                        tempPlateMatrix = []
                        for j in range(len(masterMiniPlateAList[0])):
                            tempPlateMatrix.append(np.vstack([masterMiniPlateAList[i][j],masterMiniPlateBList[i][j]]))
                        finalPlateMatrix.append(np.vstack(tempPlateMatrix))
                else:
                    finalPlate = []
                    for j in range(len(masterMiniPlateAList[0])):
                        finalPlate.append(masterMiniPlateAList[i][j])
                    finalPlateMatrix.append(np.vstack(finalPlate))
            finalMatrix = np.hstack(finalPlateMatrix)
            multiIndexedObject,timePointNames = createDataFrameLayout(folderName,dataTypeLevelList[levelIndex],dataTypeLevelNames)
            baseDataFrame = pd.DataFrame(finalMatrix,index=multiIndexedObject,columns=timePointNames)
            baseDataFrameList.append(baseDataFrame)
    else:
        for levelIndex in range(len(dataTypeLevelList)):
            maxPlates = 8
            #Number of plates of a single condition required for all timepoints
            platesPerTimePoint = math.ceil(numTimePoints/plateLength)
            platesPerCondition = int(int(numPlates/2)/platesPerTimePoint)
            upperCase = string.ascii_uppercase
            plateIndex = 0
            conditionlist = []
            for i in range(platesPerCondition):
                timepointlist = []
                for j in range(platesPerTimePoint):
                    timepointlist.append(np.reshape(list(allRawData[plateIndex].iloc[:,levelIndex+1]),(plateWidth,plateLength)))
                    plateIndex+=1
                conditionlist.append(np.hstack(timepointlist))
            beforeCorrectingAlternatingMatrix = np.vstack(conditionlist)
            finalMatrix = np.zeros(beforeCorrectingAlternatingMatrix.shape)
            for alternate in range(platesPerTimePoint):
                alternatingColumns = range(alternate,platesPerTimePoint*plateLength,platesPerTimePoint)
                i=0
                for column in range(alternate*plateLength,(alternate+1)*plateLength):
                    finalMatrix[:,alternatingColumns[i]] = beforeCorrectingAlternatingMatrix[:,column]
                    i+=1
            multiIndexedObject,timePointNames = createDataFrameLayout(folderName,dataTypeLevelList[levelIndex],dataTypeLevelNames)
            baseDataFrame = pd.DataFrame(finalMatrix,index=multiIndexedObject,columns=timePointNames)
            baseDataFrameList.append(baseDataFrame)
    
    combinedBaseDataFrame = pd.concat(baseDataFrameList)
    combinedBaseDataFrame.columns.name = 'Time'
    return combinedBaseDataFrame

def createDataFrameLayout(folderName,newLevel,newLevelListNames):

    #Opens experiment parameter file created by setup experiment script and user input with num conditions,timepoints, and condition names
    with open('inputFiles/experimentParameters-'+folderName+'.json') as f:
        experimentParameters = json.load(f)
    numTimePoints = experimentParameters[0][1]
    
    #Uses from product if condition labels not entered in manually, otherwise reads tuples for each condition level with all condition labels within it (e.g. (CD28,N4,100pM),(CD28,N4,10pM)...)
    if(experimentParameters[5]):
        initialConditionList = list(zip(*experimentParameters[2]))
    else:
        multiIndexedObject = pd.MultiIndex.from_product(experimentParameters[2], names=experimentParameters[1])
        initialConditionList = list(multiIndexedObject)
    fullMultiIndexList = []
    for condition in initialConditionList:
        fullMultiIndexList.append(newLevel+list(condition))
    multiIndexedObject = pd.MultiIndex.from_tuples(fullMultiIndexList, names=newLevelListNames+experimentParameters[1])

    #Adds time units to column headers
    timepoints = experimentParameters[3]
    timePointNames = []
    for i in range(0,len(timepoints)):
        timePointNames.append(float(timepoints[i]))
        #timePointNames.append('+'+str(timepoints[i])+'hr')
    
    return multiIndexedObject,timePointNames

#Retrieve data about experiment, create a multiIndex object, and put experiment data into appropriate locations in the multiIndexed dataframe
def createFullDataFrames(folderName,secondPath,experimentNumber,concUnit,concUnitPrefix,dataType):
    
    #Opens experiment parameter file created by setup experiment script and user input with num conditions,timepoints, and condition names
    with open('inputFiles/experimentParameters-'+folderName+'.json') as f:
        experimentParameters = json.load(f)
    numTimePoints = experimentParameters[0][1]
    
    paired = experimentParameters[6]
    contiguous = experimentParameters[7]
    #Grabs plate names (A1,A2 etc.)
    plateNames = experimentParameters[4]
    numPlates = len(plateNames)
    if not paired:
        replicateWise = True
        alternatingPlatewise = experimentParameters[8]
        numPlates*=2
    else:
        replicateWise = experimentParameters[8]
        alternatingPlatewise = False
    
    processedData = []
    if dataType == 'cyt':
        allRawData,newLevelList = cleanUpFlowjoCSV(plateNames,folderName,dataType)
        finalDataFrame = createAndCombineBaseDataFrames(folderName,allRawData,numPlates,numTimePoints,newLevelList,cytokineHeaderNames,paired,contiguous,replicateWise,alternatingPlatewise)

        with open('semiProcessedData/cytokineGFIPickleFile-'+folderName+'.pkl', "wb") as f:
            pickle.dump(finalDataFrame, f)
        
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
        with open('semiProcessedData/cytokineConcentrationPickleFile-'+folderName+'.pkl', "wb") as f:
            pickle.dump(finalDataFrameConcentration, f)
        #Create and save modified Df (original concentration df with various modifications due to experimental error)
        #modifiedDf = returnModifiedDf(experimentNumber,finalDataFrame,dataType)
        #with open('semiProcessedData/cytokineGFIPickleFile-'+folderName+'-modified.pkl', "wb") as f:
        #    pickle.dump(modifiedDf, f)
        modifiedConcentrationDf = returnModifiedDf(experimentNumber,finalDataFrameConcentration,dataType)
        with open('semiProcessedData/cytokineConcentrationPickleFile-'+folderName+'-modified.pkl', "wb") as f:
            pickle.dump(modifiedConcentrationDf, f)
        finalDataFrame = modifiedConcentrationDf
    if dataType == 'cell':
        allRawData,newLevelList = cleanUpFlowjoCSV(plateNames,folderName,dataType)
        finalDataFrame = createAndCombineBaseDataFrames(folderName,allRawData,numPlates,numTimePoints,newLevelList,cellHeaderNames,paired,contiguous,replicateWise,alternatingPlatewise)
        #Create and save modified Df (original concentration df with various modifications due to experimental error)
        with open('semiProcessedData/cellStatisticPickleFile-'+folderName+'.pkl', "wb") as f:
            pickle.dump(finalDataFrame, f)
        modifiedDataFrame = returnModifiedDf(experimentNumber,finalDataFrame,dataType)
        with open('semiProcessedData/cellStatisticPickleFile-'+folderName+'-modified.pkl', "wb") as f:
            pickle.dump(modifiedDataFrame, f)
        finalDataFrame = modifiedDataFrame
    else:
        pass

    return finalDataFrame

def convertDataFramestoExcel(folderName,secondPath,dataType,df,useModifiedDf):
    writer = pd.ExcelWriter('semiProcessedData/excelFile-'+folderName+'-'+dataType+'.xlsx')
    if dataType == 'cyt':
        dfg = pickle.load(open('semiProcessedData/cytokineGFIPickleFile-'+folderName+'-modified.pkl','rb'))
        dfc = pickle.load(open('semiProcessedData/cytokineConcentrationPickleFile-'+folderName+'-modified.pkl','rb'))
        dfg.to_excel(writer,'GFI')
        dfc.to_excel(writer,'Concentration')
    else:
        for statistic in list(pd.unique(df.index.get_level_values('Statistic'))):
            statisticDf = df.xs(statistic,level='Statistic')
            statisticDf.to_excel(writer,statistic)
    writer.save()
    print(df)
    print(dataType[0].upper()+dataType[1:]+' Excel file Saved')

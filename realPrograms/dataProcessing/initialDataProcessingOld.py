#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import json,pickle,math,matplotlib,sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from itertools import groupby
from miscFunctions import Hill,InverseHill,r_squared,cleanUpCSV,cleanUpCellCSV,cleanUpProliferationData
from modifyDataFrames import returnModifiedDf,returnModifiedCellDf

listOfCytokines=['IFNg','IL-2','IL-4','IL-6','IL-10','IL-17A','TNFa']
cellHeaderNames = ['CellType','Marker','Statistic']
MWofCytokines=[17200,17200,14000,21900,18900,15500,17500] #g/mol, correct one
reshapeReplicateWise = True

#Create Calibration Curves, obtain LODs (limits of detection) of experiment
def calibrateExperiment(folderName,secondPath,concUnit,concUnitPrefix,numberOfCalibrationSamples,initialStandardVolume):
    #Get cytokine calibration curve data
    calibration = cleanUpCSV(['Calibration'])[0]
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
    fig1.savefig(secondPath+'calibrationCurves/calibrationCurves-'+folderName+'-'+concUnitPrefix+'.png')
    #Save fitting parameters and LOD for curve fit for each cytokine
    with open('semiProcessedData/fittingParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "wb") as f:
        pickle.dump(fittingParameters, f)
    with open('semiProcessedData/LODParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "wb") as f:
        pickle.dump(concLOD, f)

#Retrieve data about experiment, create a multiIndex object, and put experiment data into appropriate locations in the multiIndexed dataframe
def createCytokineDataFrames(folderName,secondPath,experimentNumber,concUnit,concUnitPrefix):

    #Opens experiment parameter file created by setup experiment script and user input with num conditions,timepoints, and condition names
    with open('inputFiles/experimentParameters-'+folderName+'.json') as f:
        experimentParameters = json.load(f)
    numConditions = experimentParameters[0][0]
    numTimePoints = experimentParameters[0][1]
    
    #Uses from product if condition labels not entered in manually, otherwise reads tuples for each condition level with all condition labels within it (e.g. (CD28,N4,100pM),(CD28,N4,10pM)...)
    if(experimentParameters[5]):
        multiIndexedObject = pd.MultiIndex.from_tuples(list(zip(*experimentParameters[2])), names=experimentParameters[1])
    else:
        multiIndexedObject = pd.MultiIndex.from_product(experimentParameters[2], names=experimentParameters[1])
    
    print(experimentParameters)
    
    #Adds time units to column headers
    timepoints = experimentParameters[3]
    timePointNames = []
    for i in range(0,len(timepoints)):
        timePointNames.append(str(timepoints[i]))
        #timePointNames.append('+'+str(timepoints[i])+'hr')
    
    #Grabs plate names (A1,A2 etc.)
    plateNames = experimentParameters[4]
    numPlates = len(plateNames)
    allRawData = cleanUpCSV(plateNames)
    processedData = []
    #Conditions dimension
    plateWidth = 8
    #Timepoints dimension
    plateLength = 12
    
    #Get number of timepoints per plate, which is equal to the total length of the plate (12 columns) divided by the number of conditions per plate (numConditions/2). If numConditions is not a multiple of 16 (for example 14), round up to the highest multiple of 16 (math.ceil)
    numTimePointsPerPlate = int((2*(plateLength*plateWidth))/((2*plateWidth)*math.ceil(numConditions/(2*plateWidth))))
    print(numTimePointsPerPlate) 
    #Iterate across cytokines
    for cytIndex in range(len(listOfCytokines)):
        cyt=listOfCytokines[cytIndex]
        rawData = []
        numTimePointsRemaining = numTimePoints
        #Iterates through each plate 2 at a time (robot processes 2 plates at a time: A/B plates; each A/B set combined has all the conditions for a given experiment)
        for i in range(0,int(numPlates/2)):
            plateIndexA = i*2
            plateIndexB = plateIndexA + 1
            tempDataA = allRawData[plateIndexA][cyt]
            tempDataB = allRawData[plateIndexB][cyt]
            #Figure out how many columns are on this set of A/B plates. Normally 12 (as the whole plate is used), but the last set of plates can have fewer than 12 if doing an experiment with a timepoint number that isn't a multiple of numTimePointsPerPlate
            if(numTimePointsRemaining < numTimePointsPerPlate):
                numColumnsOnPlate = int(plateLength * numTimePointsRemaining/numTimePointsPerPlate)
            else:
                numColumnsOnPlate = plateLength
            if(not reshapeReplicateWise):
                #Combine two series into a single list (1-96 rowMajor on plate A, then 1-96 rowMajor on plate B)
                combinedList = np.hstack((tempDataA,tempDataB))
                #Reshape into a single combined plate with up to 16 rows and 12 columns, then append to reshaped lists. Row major reshaping, so t1-12 c1, t1-12 c2 etc.
                #USE FOR ITAM-4
                #reshapedCombinedList = np.reshape(combinedList,(13,11))
                #USE FOR EVERYTHING ELSE
                reshapedCombinedList = np.reshape(combinedList,(int(len(combinedList)/numColumnsOnPlate),numColumnsOnPlate))
            else:
                #tempReshapedListA = np.reshape(list(tempDataA),(8,12))
                #tempReshapedListB = np.reshape(list(tempDataB),(8,12))
                tempReshapedListA = np.reshape(list(tempDataA),(int(len(tempDataA)/numColumnsOnPlate),numColumnsOnPlate))
                tempReshapedListB = np.reshape(list(tempDataB),(int(len(tempDataB)/numColumnsOnPlate),numColumnsOnPlate))
                #dataMatrixA = np.reshape(tempReshapedListA.flatten('F'),(32,3),order='F')
                #dataMatrixB = np.reshape(tempReshapedListB.flatten('F'),(32,3),order='F')
                dataMatrixA = np.reshape(tempReshapedListA.flatten('F'),(int(numConditions/2),numTimePointsPerPlate),order='F')
                dataMatrixB = np.reshape(tempReshapedListB.flatten('F'),(int(numConditions/2),numTimePointsPerPlate),order='F')
                reshapedCombinedList = np.vstack((dataMatrixA,dataMatrixB))
                #print(reshapedCombinedList)
            rawData.append(reshapedCombinedList)
            numTimePointsRemaining -= numTimePointsPerPlate
        #start adding combined A/B plates to each other horizontally
        combinedRawData = rawData[0]
        for i in range(1,len(rawData)):
            combinedRawData = np.hstack((combinedRawData,rawData[i]))
        #Flatten combined A/B plates column wise (so t1c1-16, t1c17-32 etc.), then reshape into appropriately sized dataframe based on numConditions (rows) and numTimepoitns (columns)
        #USE FOR ITAM-4
        #dataMatrix = rawData[0]
        #USE FOR EVERYTHING ELSE
        dataMatrix = np.reshape(combinedRawData.flatten('F'),(numConditions,numTimePoints),order='F')
        #Assign previously created multindex to row names of dataframe, then add column level name
        fullDataFrame = pd.DataFrame(dataMatrix,multiIndexedObject,timePointNames)
        fullDataFrame.columns.name = 'Time'
        processedData.append(fullDataFrame)
    #Add Cytokine level to previous multindexed object, then concatenate all dataframes together
    newLevels = multiIndexedObject.names.copy()
    newLevels.insert(0,'Cytokine')
    completeDataSet = pd.concat(processedData, keys=listOfCytokines,names=newLevels)
    with open('semiProcessedData/cytokineGFIPickleFile-'+folderName+'.pkl', "wb") as f:
        pickle.dump(completeDataSet, f)
    
    fittingParameters = pickle.load(open('semiProcessedData/fittingParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "rb"))
    LODParameters = pickle.load(open('semiProcessedData/LODParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "rb"))
    
    #Begin converting GFI dataframe into corresponding concentration dataframe
    completeDataSetConcentration = completeDataSet.copy()
    #Step through dataframe one cytokine at a time
    for cytokine in completeDataSet.index.levels[0]:
        #Retrieve LODs for current cytokine (from constructed calibration curve)
        lowerGFILOD = LODParameters[cytokine][0]
        upperGFILOD = LODParameters[cytokine][1]
        lowerConcLOD = LODParameters[cytokine][2]
        upperConcLOD = LODParameters[cytokine][3]
        cyt = listOfCytokines.index(cytokine)
        #Loop through every value in current cytokine's portion of the dataframe
        for i in range(0,completeDataSet.loc[cytokine].values.shape[0]):
            for j in range(0,completeDataSet.loc[cytokine].values.shape[1]):
                currentGFIval = completeDataSet.loc[cytokine].values[i,j]
                if currentGFIval > upperGFILOD: #If intensity is greater than upper GFI LOD
                    currentConcVal = upperConcLOD #Concentration is equal to upper concentration LOD
                elif currentGFIval < upperGFILOD and currentGFIval > lowerGFILOD: #if intensity is between upper and lower GFI LODs
                    currentConcVal = InverseHill(np.log10(currentGFIval),fittingParameters[cyt,:]) #Use previous hill fit parameters for the cytokine to obtain concentration
                else: #If intensity is less than background GFI LOD
                    currentConcVal = lowerConcLOD #Concentration is equal to lower concentration LOD
                completeDataSetConcentration.loc[cytokine].values[i,j] = currentConcVal
    
    #Fill in n/a's (caused by undetectable fluorescence) with lower concentration LOD
    completeDataSetConcentration.fillna(lowerConcLOD, inplace=True)
    with open('semiProcessedData/cytokineConcentrationPickleFile-'+folderName+'.pkl', "wb") as f:
        pickle.dump(completeDataSetConcentration, f)
    with open(secondPath+'cytokineConcentrationPickleFiles/cytokineConcentrationPickleFile-'+folderName+'.pkl', "wb") as f:
        pickle.dump(completeDataSetConcentration, f)
    #Create and save modified Df (original concentration df with various modifications due to experimental error)
    modifiedDf = returnModifiedDf(experimentNumber,completeDataSet)
    with open('semiProcessedData/modifiedCytokineGFIPickleFile-'+folderName+'.pkl', "wb") as f:
        pickle.dump(modifiedDf, f)
    modifiedConcentrationDf = returnModifiedDf(experimentNumber,completeDataSetConcentration)
    with open('semiProcessedData/modifiedCytokineConcentrationPickleFile-'+folderName+'.pkl', "wb") as f:
        pickle.dump(modifiedConcentrationDf, f)
    with open(secondPath+'modifiedCytokineConcentrationPickleFiles/modifiedCytokineConcentrationPickleFile-'+folderName+'.pkl', "wb") as f:
        pickle.dump(modifiedDf, f)
    print(modifiedDf)

#Retrieve data about experiment, create a multiIndex object, and put experiment data into appropriate locations in the multiIndexed dataframe
def createCellDataFrames(folderName,secondPath,experimentNumber):

    #Opens experiment parameter file created by setup experiment script and user input with num conditions,timepoints, and condition names
    with open('inputFiles/experimentParameters-'+folderName+'.json') as f:
        experimentParameters = json.load(f)
    numConditions = experimentParameters[0][0]
    numTimePoints = experimentParameters[0][1]
    
    #Grabs plate names (A1,A2 etc.)
    plateNames = experimentParameters[4]
    numPlates = len(plateNames)
    allRawData,cellStatisticList = cleanUpCellCSV(plateNames,folderName)
    print(cellStatisticList)
    #Uses from product if condition labels not entered in manually, otherwise reads tuples for each condition level with all condition labels within it (e.g. (CD28,N4,100pM),(CD28,N4,10pM)...)
    if(experimentParameters[5]):
        initialConditionList = list(zip(*experimentParameters[2]))
        fullMultiIndexList = []
        for cellStatistic in cellStatisticList:
            for condition in initialConditionList:
                fullMultiIndexList.append(cellStatistic+list(condition))
        multiIndexedObject = pd.MultiIndex.from_tuples(fullMultiIndexList, names=cellHeaderNames+experimentParameters[1])
    else:
        multiIndexedObject = pd.MultiIndex.from_product(experimentParameters[2], names=experimentParameters[1])
    
    #Adds time units to column headers
    timepoints = experimentParameters[3]
    timePointNames = []
    for i in range(0,len(timepoints)):
        timePointNames.append(str(timepoints[i]))
        #timePointNames.append('+'+str(timepoints[i])+'hr')
    
    processedData = []
    #Conditions dimension
    plateWidth = 8
    #Timepoints dimension
    plateLength = 12
    
    #Get number of timepoints per plate, which is equal to the total length of the plate (12 columns) divided by the number of conditions per plate (numConditions/2). If numConditions is not a multiple of 16 (for example 14), round up to the highest multiple of 16 (math.ceil)
    numTimePointsPerPlate = int((2*(plateLength*plateWidth))/((2*plateWidth)*math.ceil(numConditions/(2*plateWidth))))
    #Iterate across all cell statistics
   for cellIndex in range(len(cellStatisticList)):
        rawData = []
        numTimePointsRemaining = numTimePoints
        print(cellStatisticList[cellIndex])
        #Iterates through each plate 2 at a time (robot processes 2 plates at a time: A/B plates; each A/B set combined has all the conditions for a given experiment)
        for i in range(0,int(numPlates/2)):
            plateIndexA = i*2
            plateIndexB = plateIndexA + 1
            tempDataA = allRawData[plateIndexA].iloc[:,cellIndex+1]
            tempDataB = allRawData[plateIndexB].iloc[:,cellIndex+1]
            #Combine two series into a single list (1-96 rowMajor on plate A, then 1-96 rowMajor on plate B)
            combinedList = np.hstack((tempDataA,tempDataB))
            #Figure out how many columns are on this set of A/B plates. Normally 12 (as the whole plate is used), but the last set of plates can have fewer than 12 if doing an experiment with a timepoint number that isn't a multiple of numTimePointsPerPlate
            if(numTimePointsRemaining < numTimePointsPerPlate):
                numColumnsOnPlate = int(plateLength * numTimePointsRemaining/numTimePointsPerPlate)
            else:
                numColumnsOnPlate = plateLength
            if(not reshapeReplicateWise):
                #Combine two series into a single list (1-96 rowMajor on plate A, then 1-96 rowMajor on plate B)
                combinedList = np.hstack((tempDataA,tempDataB))
                #Reshape into a single combined plate with up to 16 rows and 12 columns, then append to reshaped lists. Row major reshaping, so t1-12 c1, t1-12 c2 etc.
                #USE FOR ITAM-4
                #reshapedCombinedList = np.reshape(combinedList,(13,11))
                #USE FOR EVERYTHING ELSE
                reshapedCombinedList = np.reshape(combinedList,(int(len(combinedList)/numColumnsOnPlate),numColumnsOnPlate))
            else:
                #tempReshapedListA = np.reshape(list(tempDataA),(8,12))
                #tempReshapedListB = np.reshape(list(tempDataB),(8,12))
                tempReshapedListA = np.reshape(list(tempDataA),(int(len(tempDataA)/numColumnsOnPlate),numColumnsOnPlate))
                tempReshapedListB = np.reshape(list(tempDataB),(int(len(tempDataB)/numColumnsOnPlate),numColumnsOnPlate))
                #dataMatrixA = np.reshape(tempReshapedListA.flatten('F'),(32,3),order='F')
                #dataMatrixB = np.reshape(tempReshapedListB.flatten('F'),(32,3),order='F')
                dataMatrixA = np.reshape(tempReshapedListA.flatten('F'),(int(numConditions/2),numTimePointsPerPlate),order='F')
                dataMatrixB = np.reshape(tempReshapedListB.flatten('F'),(int(numConditions/2),numTimePointsPerPlate),order='F')
                reshapedCombinedList = np.vstack((dataMatrixA,dataMatrixB))
                #print(reshapedCombinedList)
            rawData.append(reshapedCombinedList)
            numTimePointsRemaining -= numTimePointsPerPlate
        #start adding combined A/B plates to each other horizontally
        combinedRawData = rawData[0]
        for i in range(1,len(rawData)):
            combinedRawData = np.hstack((combinedRawData,rawData[i]))
        #Flatten combined A/B plates column wise (so t1c1-16, t1c17-32 etc.), then reshape into appropriately sized dataframe based on numConditions (rows) and numTimepoitns (columns)
        #USE FOR ITAM-4
        #dataMatrix = rawData[0]
        #USE FOR EVERYTHING ELSE
        dataMatrix = np.reshape(combinedRawData.flatten('F'),(numConditions,numTimePoints),order='F')
        processedData.append(dataMatrix)
    fullDataMatrix = np.zeros((processedData[0].shape[0]*len(cellStatisticList),processedData[0].shape[1]))
    print(len(processedData))
    if(not reshapeReplicateWise):
        for i in range(len(cellStatisticList)):
            fullDataMatrix[i*processedData[0].shape[0]:(i+1)*processedData[0].shape[0],:] = processedData[i]
    else:
        for i in range(len(cellStatisticList)):
            fullDataMatrix[i*processedData[0].shape[0]:(i+1)*processedData[0].shape[0],:] = processedData[i]
    #print(multiIndexedObject)
    completeDataSet = pd.DataFrame(fullDataMatrix, multiIndexedObject,timePointNames)
    completeDataSet.columns.name = 'Time'
    with open('semiProcessedData/cellStatisticPickleFile-'+folderName+'.pkl', "wb") as f:
        pickle.dump(completeDataSet, f)
    modifiedDf = returnModifiedCellDf(experimentNumber,completeDataSet)
    with open('semiProcessedData/modifiedCellStatisticPickleFile-'+folderName+'.pkl', "wb") as f:
        pickle.dump(modifiedDf, f)

def convertCytokineDataFramestoExcel(folderName,secondPath,useModifiedDf):
    if(useModifiedDf):
        writer = pd.ExcelWriter('semiProcessedData/excelFile-'+folderName+'-modified.xlsx')
        writer2 = pd.ExcelWriter(secondPath+'excelFiles/excelFile-'+folderName+'-modified.xlsx')
        dfg = pickle.load(open('semiProcessedData/modifiedCytokineGFIPickleFile-'+folderName+'.pkl','rb'))
        dfc = pickle.load(open('semiProcessedData/modifiedCytokineConcentrationPickleFile-'+folderName+'.pkl','rb'))
    else:
        writer = pd.ExcelWriter('semiProcessedData/excelFile-'+folderName+'.xlsx')
        writer2 = pd.ExcelWriter(secondPath+'excelFiles/excelFile-'+folderName+'.xlsx')
        dfg = pickle.load(open('semiProcessedData/cytokineGFIPickleFile-'+folderName+'.pkl','rb'))
        dfc = pickle.load(open('semiProcessedData/cytokineConcentrationPickleFile-'+folderName+'.pkl','rb'))
    dfg.to_excel(writer,'GeomFI')
    dfc.to_excel(writer,'Concentration')
    dfg.to_excel(writer2,'GeomFI')
    dfc.to_excel(writer2,'Concentration')
    writer.save()
    writer2.save()
    print(dfc)
    print('Cytokine Excel file Saved')

def convertCellDataFramestoExcel(folderName,secondPath,useModifiedDf):
    if(useModifiedDf):
        df = pickle.load(open('semiProcessedData/modifiedCellStatisticPickleFile-'+folderName+'.pkl','rb'))
        writer = pd.ExcelWriter('semiProcessedData/excelFile-'+folderName+'-cell.xlsx')
        writer2 = pd.ExcelWriter(secondPath+'excelFiles/excelFile-'+folderName+'-cell.xlsx')
    else:
        df = pickle.load(open('semiProcessedData/cellStatisticPickleFile-'+folderName+'.pkl','rb'))
        writer = pd.ExcelWriter('semiProcessedData/excelFile-'+folderName+'-cell.xlsx')
        writer2 = pd.ExcelWriter(secondPath+'excelFiles/excelFile-'+folderName+'-cell.xlsx')
    for statistic in list(pd.unique(df.index.get_level_values('Statistic'))):
        statisticDf = df.xs(statistic,level='Statistic')
        statisticDf.to_excel(writer,statistic)
        statisticDf.to_excel(writer2,statistic)
    writer.save()
    writer2.save()
    print(df)
    print('Cell Excel file Saved')

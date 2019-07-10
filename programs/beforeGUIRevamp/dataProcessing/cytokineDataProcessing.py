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

#Standard BD Biosciences CBA Kit Cytokines
listOfCytokines=['IFNg','IL-2','IL-4','IL-6','IL-10','IL-17A','TNFa']
MWofCytokines=[17200,17200,14000,21900,18900,15500,17500] #g/mol, correct one

def parseCytokineCSVHeaders(columns):
    #,Beads/IFNg | Geometric Mean (YG586-A),Beads/IL-2 | Geometric Mean (YG586-A),Beads/IL-4 | Geometric Mean (YG586-A),Beads/IL-6 | Geometric Mean (YG586-A),Beads/IL-10 | Geometric Mean (YG586-A),Beads/IL-17A | Geometric Mean (YG586-A),Beads/TNFa | Geometric Mean (YG586-A),
    newMultiIndexList = []
    for column in columns[1:-1]:
        populationNameVsStatisticSplit = column.split(' | ')
        cytokine = populationNameVsStatisticSplit[0].split('/')[-1]
        newMultiIndexList.append([cytokine])
    return newMultiIndexList

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

def createCytokineDataFrame(folderName,finalDataFrame,concUnitPrefix):
        
        columnName = finalDataFrame.columns.name
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
        finalDataFrameConcentration.columns.name = columnName 

        return finalDataFrameConcentration

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thursday July 5 13:52:27 2018

@author: acharsr
"""
from pathlib import Path
import pickle
from itertools import product,combinations
import numpy as np
import os,sys
import math
import json
import subprocess

modelNames = ['ae-5-3-5-raw','ae-5-3-5-partitioned','ae-5-3-5-raw-normPerObservable','ae-5-3-5-partitioned-normPerObservable'\
        ,'ae-7-3-7-tanh-raw','ae-7-3-7-tanh-partitioned','ae-7-3-7-tanh-raw-normPerObservable','ae-7-3-7-tanh-partitioned-normPerObservable',\
        'ae-11-3-11-parameterized','ae-16-3-16-parameterized']

#Make experiment folder and subfolders
def createExperimentFolders(folderName):
    subprocess.run(['mkdir',folderName])
    experimentFolderNames = ['inputFiles','semiProcessedData','fullyProcessedFigures','semiProcessedData/singleCellData']
    for experimentFolderName in experimentFolderNames:
        subprocess.run(['mkdir',folderName+'/'+experimentFolderName])

#Make experiment folder and subfolders
def createNNOutputFolders(outputFolderName):
    trainingNumbers = [28,30,67]
    comboLengths = [1,2,3]
    allTrainingCombinations = []
    for comboLength in comboLengths:
        allTrainingCombinations.append(list(combinations(trainingNumbers,comboLength)))
    trainingCombinations = []
    for elem1 in allTrainingCombinations:
        for elem2 in elem1:
            trainingCombinations.append(elem2)
    print(trainingCombinations)
    os.chdir('../../output/'+outputFolderName)
    for modelName in modelNames:
        if not os.path.exists(modelName):
            subprocess.run(['mkdir',modelName])
        os.chdir(modelName)
        for trainingCombination in trainingCombinations:
            trainingCombinationName = ','.join(map(str,trainingCombination))
            if not os.path.exists(trainingCombinationName):
                subprocess.run(['mkdir',trainingCombinationName])
        os.chdir('..')
            
#Create .json file with lists used to construct multiIndex; based entirely on user input. Json file is human readable so mistakes can be corrected directly in the file after the script writes it
def createParameters(folderName,numConditions,numTimePoints,showOldHeatMap):
    
    print('Creating Parameters for: '+str(folderName))    
    print('Number of Conditions: '+str(numConditions))
    print('Number of TimePoints: '+str(numTimePoints))
    
    plateLength = 12
    plateWidth = 8
    
    """
    #Check to see if old timeInfo file exists; if so use that to construct timepoint list
    timeInfoPath = Path('inputFiles/timeInfo-'+folderName+'.txt')
    if(timeInfoPath.is_file()):
        with open('inputFiles/timeInfo-'+folderName+'.txt','r') as f:
            timeArray = f.read().splitlines() 
            totalTimeHours = timeArray[-1]
        #When making multiIndex, open previous heatmaps if they exist to get a quick look at what conditions the experiment had
        if(showOldHeatMap):
            subprocess.run(["open","fullyProcessedFigures/cbaheatmaps-"+folderName+"-tight-nM.png"]) 
    """

    numPlates = int(input('Enter the total number of plates used in the experiment: '))
    paired = str(input('Was this experiment performed with A/B plates (y/n)? '))
    if numConditions % plateWidth == 0:
        contiguous = 'y'
    else:
        contiguous = str(input('Are the only empty wells in a timepoint at the end of the timepoint\'s columns (y/n)? '))
    if paired != 'y':
        replicateWise = 'y'
    else:
        replicateWise = str(input('Do conditions 1-16 in this experiment come from PlateA only or do they also involve PlateB? '))
    #Asks for user input for number of conditions levels, and if they want to fill in levels manually (one tuple at a time). Used mainly for the peptide/concentration levels, as there are different concentrations associated with each peptide (concentrations are not replicated across peptides)
    numConditionLevels = int(input('Enter the number of condition levels: '))
    manualFill = str(input('Do you want to fill in your condition labels array manually (y/n)? '))

    conditionLabelsArray = []
    nestedConditionLabelsArray = []
    conditionLevelsArray = []
    conditionLevelSizeArray = []
    prod = 1
    
    #If manual fill selected, ask which level to start constructing the tuple (can use list product to fill in previous levels)
    if manualFill == 'y':
        manualStart = int(input('At which condition level do you want to start manually filling in condition labels? ')) - 1
        for i in range(0,numConditionLevels):
            conditionlabel = input('Enter name of condition level '+str(i+1)+': ')
            conditionLevelsArray.append(conditionlabel)
        for i in range(0,manualStart):
            conditionLevelSizeArray.append(int(input('Enter the number of labels in the \"'+str(conditionLevelsArray[i])+'\" condition level: ')))
        automaticLabels = []
        for i in range(0,manualStart):
            temp = []
            for j in range(0,conditionLevelSizeArray[i]):
                currentConditionLabel = input('Enter label for condition ' + str(j+1) + ' in condition level \"' + str(conditionLevelsArray[i]) + '\": ')
                temp.append(currentConditionLabel)
            automaticLabels.append(temp)
        for i in range(0,len(automaticLabels)):
            temp = []
            labels = automaticLabels[i]
            prod*=conditionLevelSizeArray[i]
            for j in range(0,prod):
                #0-31 --> 0, 31-63 --> 1; 0-15 --> 0, 16-31 --> 1,32-47 --> 2,48-63 --> 3
                for k in range(0,int(numConditions/prod)):
                    temp.append(labels[j%len(labels)])
            nestedConditionLabelsArray.append(temp)
        for j in range(manualStart,numConditionLevels):
            temp = []
            for i in range(0,numConditions):
                currentConditionLabel = input('Enter label for condition level \"' + str(conditionLevelsArray[j]) + '\" in condition ' + str(i+1) + ': ')
                temp.append(currentConditionLabel)
            nestedConditionLabelsArray.append(temp)
    else:
        for i in range(0,numConditionLevels):
            conditionlabel = input('Enter name of condition level '+str(i+1)+': ')
            conditionLevelsArray.append(conditionlabel)
        for i in range(0,numConditionLevels-1):
            conditionLevelSizeArray.append(int(input('Enter the number of labels in the \"'+str(conditionLevelsArray[i])+'\" condition level: ')))
            prod*=conditionLevelSizeArray[i]
        conditionLevelSizeArray.append(int(numConditions/prod))
        print(conditionLevelSizeArray)
        for i in range(0,numConditionLevels):
            currentConditionLevelLabelArray = []
            for j in range(0,conditionLevelSizeArray[i]):
                currentConditionLabel = input('Enter label for condition ' + str(j+1) + ' in condition level ' + str(i+1) + ': ')
                currentConditionLevelLabelArray.append(currentConditionLabel)
            nestedConditionLabelsArray.append(currentConditionLevelLabelArray) 
    
    manualTimeFill = str(input('Do you  want to enter timepoints in manually (y/n)?'))
    if manualTimeFill == 'y':
        timeArray= []
        for tp in range(1,numTimePoints+1):
            timepoint = int(input('Enter timepoint '+str(tp)+': '))
            timeArray.append(timepoint)
    #Otherwise, make new timepoint list from scratch (starting from 1 - 100 hours)
    else:
        totalTimeHours = 100 
        #Array for timepoint times in hours
        baseTimeArrayHours = [0]
        for i in range(1,numTimePoints+1):
        #Want amount of time between timepoints to go up in a logarithmic fashion so can simply go from 10^0 (starting point) to 10^2 (ending point) in numTimePoints+1 steps (rounding to the nearest whole hour)
            currentTime = round(10**(math.log10(totalTimeHours)*(i/numTimePoints)))
            baseTimeArrayHours.append(currentTime)
        #Find time between each timepoint; this is the time in minutes that we tell the robot to wait after it finishes taking a timepoint
        timeDifferentialArray = []
        for i in range(0,numTimePoints):
            timeDifferentialArray.append(baseTimeArrayHours[i+1] - baseTimeArrayHours[i])
        #Remove initial 0 hour timepoint
        baseTimeArrayHours.pop(0)
        #If there is no time between two timepoints (because of rounding, a timepoint at 1.1 hours and 1.3 hours will both be 1), timePointdifferential array will be 0; we will set to 1 hour (lowest amount of time robot can wait and still finish a timepoint)
        hoursAdded = 0
        for i in range(0,len(baseTimeArrayHours)):
            if(timeDifferentialArray[i] == 0):
                timeDifferentialArray[i] = 1
                hoursAdded+=1
        #Setting timedifferential arrays to 1 throws off total time; to compensate we subtract the number of changed time points from the last timepoint (which is the longest) to keep overall experiment time at totalTimeHours
        timeDifferentialArray[-1] -= hoursAdded
        #Sorting time differential array in ascending order allows for the time between the last two timepoints to be the longest; because of the earlier subtraction the last element in timedifferential array can be shorter than the next to last
        baseTimeArrayHours = sorted(timeDifferentialArray)
        #Time differential array cumulative sum gives the total time elapsed at each timepoint; used for column values in dataframe
        timeArray = np.cumsum(baseTimeArrayHours).tolist()
    plateNameArray = []
    if paired == 'y':
        for i in range(0,int(numPlates/2)):
            plateNameArray.append('A'+str(i+1))
            plateNameArray.append('B'+str(i+1))
    else:
        for i in range(0,numPlates):
            print('wat')
            print(i)
            plateNameArray.append('A'+str(i+1))

    experimentParameters = [[numConditions,numTimePoints],conditionLevelsArray,nestedConditionLabelsArray,timeArray,plateNameArray,manualFill == 'y',paired == 'y',contiguous == 'y',replicateWise =='y']
    with open('inputFiles/experimentParameters-'+folderName+'.json','w') as f:
        json.dump(experimentParameters, f)

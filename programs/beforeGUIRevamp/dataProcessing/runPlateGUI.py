#!/usr/bin/env python3
import pickle,os,sys,subprocess,json
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from miscFunctions import extractValues  

def createWellGridButtonImages(experimentParameters,pathToButtonImages,maxNumberLevelValues):
    #maxNumberLevelValues = max(experimentParameters['numConditionLevelValues']+[experimentParameters['numColumnLevelValues']])
    width = 5
    sizeInPixels = 50
    maxColor = 255

    sns.set()
    current_palette = sns.color_palette('hls',maxNumberLevelValues)
    pal = [tuple([149/maxColor,165/maxColor,166/maxColor])]+current_palette

    rgbTriples = []
    for palTuple in pal:
        rgbTriples.append(tuple([int(i * maxColor) for i in palTuple]))
    hexColors = []
    for rgbTriple in rgbTriples:
        hexColors.append('#%02x%02x%02x' % rgbTriple)

    for buttonColor,i in zip(pal,range(len(pal))):
        fig1 = plt.figure()
        fig1.set_size_inches(sizeInPixels,sizeInPixels)
        ax1 = fig1.add_subplot(111, aspect='equal')
        ax1.add_patch(patches.Rectangle((0, 0), width, width,edgecolor=buttonColor,facecolor=buttonColor))
        plt.ylim(width)
        plt.xlim(width)

        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        plt.axis('off')
        plt.savefig(pathToButtonImages+'/buttonColors/'+str(i)+'.png',bbox_inches='tight', pad_inches=0,dpi=1)
        plt.clf()

    return hexColors

def produceGUIBasedIndexingCoordinates(folderName,pathToButtonImages,experimentParameters,parametersUpdatedByGridGUI):
    levelLayouts = []
    for conditionLevel in range(experimentParameters['numAllLevels']):
        numLevelValues = len(experimentParameters['allLevelValues'][experimentParameters['allLevelNames'][conditionLevel]])
        hexcolors = createWellGridButtonImages(experimentParameters,pathToButtonImages,numLevelValues)
        hexcolorsInput = '_'.join(hexcolors)
        buttonImagesPathAndHexcolorsInput = pathToButtonImages+'__'+hexcolorsInput
        #for i in range(1,numLevelValues+1):
        #subprocess.rm()
        exitEnabledSubprocessRun('python3',pathToButtonImages+'/layoutPlateStructure.py',buttonImagesPathAndHexcolorsInput,True)
        subprocess.run('rm '+pathToButtonImages+'/buttonColors/*.png',shell=True)
        currentFullLevelLayout = pickle.load(open('inputFiles/fullLevelLayout.pkl','rb'))
        levelLayouts.append(currentFullLevelLayout)
        currentLevelLayout = pickle.load(open('inputFiles/tile-levelLayout.pkl','rb'))
        #Strip out all zero columns/rows to make next plate layout
        currentLevelLayoutNoZeros = extractValues(currentLevelLayout,0,True)
        parametersUpdatedByGridGUI['currentPlateDimensions'] = currentLevelLayoutNoZeros.shape
        parametersUpdatedByGridGUI['numLevelsUnparsed'] -= 1
        with open('inputFiles/experimentParameters-'+folderName+'.json','w') as f:
            json.dump(experimentParameters,f)
        with open('inputFiles/gui-parametersUpdatedByGridGUI.pkl','wb') as f:
            pickle.dump(parametersUpdatedByGridGUI,f)
    
    with open('inputFiles/levelLayouts.pkl','wb') as f:
        pickle.dump(levelLayouts,f)
    #Reset parameters updated by GUI back to original; to allow user to rerun just the plating GUI part of the script again
    parametersUpdatedByGridGUI['currentPlateDimensions'] = experimentParameters['overallPlateDimensions']
    parametersUpdatedByGridGUI['numLevelsUnparsed'] = experimentParameters['numAllLevels']
    with open('inputFiles/gui-parametersUpdatedByGridGUI.pkl','wb') as f:
        pickle.dump(parametersUpdatedByGridGUI,f)

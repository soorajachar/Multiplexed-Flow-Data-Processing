#!/usr/bin/env python3
import sys,os,argparse
import pandas as pd
import pickle,warnings,json,subprocess
import setUpExperiment 
import matplotlib
warnings.filterwarnings("ignore", category=DeprecationWarning)
from miscFunctions import parseCommandLineNNString,exitEnabledSubprocessRun
import runPlateGUI as rpgui 
import initialDataProcessing as idp
import cytokineDataProcessing as cydp
#import cellDataProcessing as cdp
import proliferationDataProcessing as pdp
import singleCellDataProcessing as scdp
sys.path.insert(0, '../../figuresPipeline/')
import facetPlotLibrary as fpl
import singleStainHistogramMaster as sshm
import isomapMaster as im
import mutualInformationMaster as mi
import matplotlib.pyplot as plt

pathToExperimentSpreadsheet = '../../../experiments/'
secondPath = '../../outputData'

concUnit = 1e9
unitPrefixDictionary = {1e12:'pM',1e9:'nM',1e6:'uM',1e3:'mM',1e0:'M'}
concUnitPrefix = unitPrefixDictionary[concUnit]

def runPipelinedScript(scriptToRun,inputString,useModifiedDf,cellTypeArray):
    ex_data = pd.read_excel(pathToExperimentSpreadsheet+'masterExperimentSpreadsheet.xlsx')
    experimentsToRun = parseCommandLineNNString(inputString)
    
    for expNum in experimentsToRun:
        expIndex = expNum - 1
        folderName = ex_data['Full Name'][expIndex]
        researcherName = ex_data['Researcher'][expIndex]
        experimentType = ex_data['ExperimentType'][expIndex]
        cwd = os.getcwd()
        #Local; only used by me
        if 'Volumes' not in cwd:
            if scriptToRun != 1:
                os.chdir('experiments/'+folderName)
            else:
                os.chdir('experiments/')
        #NIH Server; used by other people
        else:
            if scriptToRun != 1:
                os.chdir('../../../'+researcherName+'/experiments/'+folderName)
            else:
                os.chdir('../../../'+researcherName+'/experiments')
        print(os.getcwd())
        if(int(scriptToRun/100.) == 0): #Initial data processing
            if(scriptToRun == 1):
                print('Setting up experiment for: '+str(folderName))
                setUpExperiment.createExperimentFolders(folderName)
                print(os.getcwd())
            elif(scriptToRun == 2):
                print('Creating experiment parameters for: '+str(folderName))
                exitEnabledSubprocessRun('python3',cwd+'/setUpExperiment.py','',False)
            elif(scriptToRun == 6):
                print('Creating plate layout for: '+str(folderName))
                experimentParameters = json.load(open('inputFiles/experimentParameters-'+folderName+'.json','r'))
                parametersUpdatedByGridGUI = pickle.load(open('inputFiles/gui-parametersUpdatedByGridGUI.pkl','rb'))
                rpgui.produceGUIBasedIndexingCoordinates(folderName,cwd,experimentParameters,parametersUpdatedByGridGUI)
            elif(scriptToRun == 3):
                print('Creating Dataframes for: '+str(folderName))
                experimentParameters = json.load(open('inputFiles/experimentParameters-'+folderName+'.json','r'))
                levelLayouts = pickle.load(open('inputFiles/levelLayouts.pkl','rb'))
                experimentLevelLayoutDict = idp.tilePlateLayouts(experimentParameters,levelLayouts)
                allColumnVariableCoordinates = idp.arrangeColumVariableBasedCoordinates(experimentParameters,experimentLevelLayoutDict)
                if(cellTypeArray[0]):
                    dataType = 'cyt'
                    numberOfCalibrationSamples = ex_data['NumberOfCBAStandardDilutions'][expIndex]
                    initialStandardVolume = ex_data['CBAStandardDilutedVolume'][expIndex]
                    cydp.calibrateExperiment(folderName,secondPath,concUnit,concUnitPrefix,numberOfCalibrationSamples,initialStandardVolume)
                    basecytdf = idp.createBaseDataFrame(experimentParameters,folderName,expNum,dataType,allColumnVariableCoordinates,experimentLevelLayoutDict)
                    cytdf = cydp.createCytokineDataFrame(folderName,basecytdf,concUnitPrefix)
                    idp.saveFinalDataFrames(folderName,secondPath,expNum,dataType,cytdf,ex_data) 
                if(cellTypeArray[1]):
                    dataType = 'cell'
                    celldf = idp.createBaseDataFrame(experimentParameters,folderName,expNum,dataType,allColumnVariableCoordinates,experimentLevelLayoutDict)
                    idp.saveFinalDataFrames(folderName,secondPath,expNum,dataType,celldf,ex_data) 
                if(cellTypeArray[2]):
                    dataType = 'prolif'
                    fileList = os.listdir('semiProcessedData')
                    if (('initialSingleCellDf-channel-'+folderName+'.pkl' not in fileList) or ('initialSingleCellDf-scale-'+folderName+'.pkl' not in fileList)):
                        fileNameDf = idp.createBaseDataFrame(experimentParameters,folderName,expNum,'cell',allColumnVariableCoordinates,experimentLevelLayoutDict)
                        finaldf = scdp.createInitialSingleCellDataFrame(folderName,expNum,fileNameDf)
                    prolifdf = pdp.createProliferationSingleCellDataFrame(folderName,secondPath,expNum,useModifiedDf)
                    idp.saveFinalDataFrames(folderName,secondPath,expNum,dataType,prolifdf,ex_data) 
                if(cellTypeArray[3]):
                    dataType = 'singlecell'
                    fileList = os.listdir('semiProcessedData')
                    if (('initialSingleCellDf-channel-'+folderName+'.pkl' not in fileList) or ('initialSingleCellDf-scale-'+folderName+'.pkl' not in fileList)):
                        dataType = 'singlecell'
                        fileNameDf = idp.createBaseDataFrame(experimentParameters,folderName,expNum,'cell',allColumnVariableCoordinates,experimentLevelLayoutDict)
                        scdf = scdp.createInitialSingleCellDataFrame(folderName,expNum,fileNameDf)
                    if 'singleCellDataFrame-proliferation-'+folderName+'.pkl' in fileList:
                        if experimentType not in ['AntibodyTest']:
                            scdp.createCompleteSingleCellDf(folderName)
            if(scriptToRun == 1):
                print(os.getcwd())
                os.chdir('../../../../programs/dataProcessing/')
            else:
                os.chdir('../../../../../programs/dataProcessing/')
        
        elif(int(scriptToRun/100.) == 1): #Figures
            
            dfArray = {}
            if(useModifiedDf):
                modifiedstring = '-modified'
            else:
                modifiedstring = ''
            if(cellTypeArray[0]):
                cytokinedf = pickle.load(open('semiProcessedData/cytokineConcentrationPickleFile-'+folderName+modifiedstring+'.pkl','rb'))
                dataType = 'cyt'
                dfArray[dataType] = cytokinedf
            if(cellTypeArray[1]):
                celldf = pickle.load(open('semiProcessedData/cellStatisticPickleFile-'+folderName+modifiedstring+'.pkl','rb'))
                dataType = 'cell'
                dfArray[dataType] = celldf
            if(cellTypeArray[2]):
                proliferationdf = pickle.load(open('semiProcessedData/proliferationStatisticPickleFile-'+folderName+modifiedstring+'.pkl','rb'))
                dataType = 'prolif'
                dfArray[dataType] = proliferationdf
            if(cellTypeArray[3]):
                fileList = os.listdir('semiProcessedData')
                if 'singleCellDataFrame-complete-'+folderName+modifiedstring+'.pkl' in fileList:
                    singlecelldf = pickle.load(open('semiProcessedData/singleCellDataFrame-complete-'+folderName+modifiedstring+'.pkl','rb'))
                else:
                    singlecelldf = pickle.load(open('semiProcessedData/initialSingleCellDf-channel-'+folderName+modifiedstring+'.pkl','rb'))
                print(singlecelldf)
                dataType = 'singlecell'
                dfArray[dataType] = singlecelldf

            if(scriptToRun in [101,102,103,104,106,107,108,109,110,111,112,113,114]): #Facet plots (line,scatter,scatter w/line for expfits)
                for dfKey in dfArray:
                    if scriptToRun == 101:
                        plotType = '3d'
                        subPlotType = 'heatmap'
                    elif(scriptToRun == 102):
                        plotType = 'ordered'
                        subPlotType = 'line'
                    elif(scriptToRun == 103):
                        plotType = 'ordered'
                        subPlotType = 'scatter'
                    elif(scriptToRun == 104):
                        plotType = 'ordered'
                        subPlotType = 'expFit'
                    elif(scriptToRun == 106):
                        plotType = 'categorical'
                        subPlotType = 'bar'
                    elif(scriptToRun == 109):
                        plotType = 'categorical'
                        subPlotType = 'strip'
                    elif(scriptToRun == 110):
                        plotType = 'categorical'
                        subPlotType = 'barstrip'
                    elif(scriptToRun == 111):
                        plotType = 'categorical'
                        subPlotType = 'swarm'
                    elif(scriptToRun == 112):
                        plotType = 'categorical'
                        subPlotType = 'violin'
                    elif(scriptToRun == 113):
                        plotType = 'categorical'
                        subPlotType = 'box'
                    elif(scriptToRun == 114):
                        plotType = 'categorical'
                        subPlotType = 'pointstrip'
                    elif(scriptToRun == 107):
                        plotType = '1d'
                        subPlotType = 'histogram'
                    elif(scriptToRun == 108):
                        plotType = '1d'
                        subPlotType = 'kde'

                    print('Creating '+subPlotType+'plots for: '+str(folderName))
                    if matplotlib.get_backend() != 'Qt4Agg': 
                        plt.switch_backend('QT4Agg') #default on my system
                    fpl.facetPlottingGUI(dfArray[dfKey],plotType,subPlotType,dataType)
                    subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually = fpl.produceSubsettedDataFrames(dfArray[dfKey])
                    fpl.plotFacetedFigures(folderName,plotType,subPlotType,dataType,subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually,useModifiedDf,dfArray[dfKey])
            elif scriptToRun in [105]: #Single Cell Figures
                if scriptToRun == 105:
                    logicleDf = pickle.load(open('semiProcessedData/initialSingleCellDf-channel-'+folderName+modifiedstring+'.pkl','rb'))
                    sshm.createHistograms(logicleDf,folderName)

            os.chdir('../../programs/figuresPipeline/')

        else: #Not supported
            pass
def main():
    #os.chdir('dataProcessing')
    parser = argparse.ArgumentParser(description="Create experiment, process data, run analysis/visualization/neural network scripts on data.")
    
    parser.add_argument("-ce", action='store_true', help="Create experiment folders and subfolders.")
    parser.add_argument("-ep", action='store_true', help="Create input parameters for experiment based on user input.")
    parser.add_argument("-pgui", action='store_true', help="Create plate condition level layout files using a GUI.")
    parser.add_argument("-pd", action='store_true', help = "Create pickle files and multindexed data frame for selected experiments.")
    parser.add_argument("-fp", action='store_true', help = "Create exp fitting parameter dataframes for selected experiments.")
    
    parser.add_argument("-m", action='store_true', help = "Use non-raw dataframe (outliers/replicates removed/averaged).")
    
    parser.add_argument("-cyt", action='store_true', help = "Process cytokine data type.")
    parser.add_argument("-cell", action='store_true', help = "Process cell data type.")
    parser.add_argument("-prolif", action='store_true', help = "Process proliferation data type.")
    parser.add_argument("-sc", action='store_true', help = "Process single cell data type.")
    parser.add_argument("-all", action='store_true', help = "Process all data types.")
    parser.add_argument("-na", action='store_true', help = "Do not process using a specific data type.")
    
    parser.add_argument("-hm", action='store_true', help = "Create heatmaps of the selected experiments.")
    parser.add_argument("-lp", action='store_true', help = "Create lineplots for selected experiments.")
    parser.add_argument("-sp", action='store_true', help = "Create scatterplots for selected experiments.")
    parser.add_argument("-bp", action='store_true', help = "Create barplots for selected experiments.")
    parser.add_argument("-strip", action='store_true', help = "Create striplots for selected experiments.")
    parser.add_argument("-swarm", action='store_true', help = "Create swarmplots for selected experiments.")
    parser.add_argument("-violin", action='store_true', help = "Create violinplots for selected experiments.")
    parser.add_argument("-box", action='store_true', help = "Create boxplots for selected experiments.")
    parser.add_argument("-barstrip", action='store_true', help = "Create barstripplots for selected experiments.")
    parser.add_argument("-pointstrip", action='store_true', help = "Create pointstripplots for selected experiments.")
    parser.add_argument("-hist", action='store_true', help = "Create histograms for selected experiments.")
    parser.add_argument("-kde", action='store_true', help = "Create kde for selected experiments.")
    parser.add_argument("-ef", action='store_true', help = "Create scatterplot with exponential fit for selected experiments.")
    
    parser.add_argument("-sshm", action='store_true', help = "Create single stain logicle histograms for selected experiments.")
    
    parser.add_argument("-imsp", action='store_true', help = "Create isomap scatterplots for selected crossvalidated experiments.")
    parser.add_argument("-imdp", action='store_true', help = "Create isomap density plots for selected crossvalidated experiments.")
    parser.add_argument("-mi", action='store_true', help = "Create mutual information matrix heatmaps for selected crossvalidated experiments.")
    
    parser.add_argument("-ppnnd", action='store_true', help = "Preprocess data for neural network for selected experiments.")
    parser.add_argument("-tnnd", action='store_true', help = "Train data for neural network for selected experiments.")
    parser.add_argument("-cvnnd", action='store_true', help = "Crossvalidate data for neural network for selected experiments.")
    
    parser.add_argument("--input", dest='inputString', help ="Run specified script on these experimental data sets. Separate numbers with , or - for a range, training sets with a /, or a // to separate training from crossval sets.")
    
    args = parser.parse_args()

    if(args.all):
        cellTypeArray = [True,True,True,True]
    elif args.na:
        cellTypeArray = [False,False,False,False]
    else:
        cellTypeArray = [args.cyt,args.cell,args.prolif,args.sc]
    if(len(str(args.inputString)) > 0):
        inputString = str(args.inputString)
    else:
        inputString = ''
    if(args.ce):
        scriptToRun = 1
    elif(args.ep):
        scriptToRun = 2
    elif(args.pd):
        scriptToRun = 3
    elif(args.ppnnd):
        scriptToRun = 4
    elif(args.fp):
        scriptToRun = 5
    elif(args.pgui):
        scriptToRun = 6
    elif(args.hm):
        scriptToRun = 101
    elif(args.lp):
        scriptToRun = 102
    elif(args.sp):
        scriptToRun = 103
    elif(args.ef):
        scriptToRun = 104
    elif(args.sshm):
        scriptToRun = 105
    elif(args.bp):
        scriptToRun = 106
    elif(args.strip):
        scriptToRun = 109
    elif(args.barstrip):
        scriptToRun = 110
    elif(args.swarm):
        scriptToRun = 111
    elif(args.violin):
        scriptToRun = 112
    elif(args.box):
        scriptToRun = 113
    elif(args.pointstrip):
        scriptToRun = 114
    elif(args.hist):
        scriptToRun = 107
    elif(args.kde):
        scriptToRun = 108
    elif(args.tnnd):
        scriptToRun = 201
    elif(args.cvnnd):
        scriptToRun = 202
    elif(args.imsp):
        scriptToRun = 203
    elif(args.imdp):
        scriptToRun = 204
    elif(args.mi):
        scriptToRun = 205
    runPipelinedScript(scriptToRun,inputString,args.m,cellTypeArray)
main()

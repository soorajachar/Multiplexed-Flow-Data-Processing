#!/usr/bin/env python3
import sys,os,argparse
import pandas as pd
import pickle,warnings
import setUpExperiment
warnings.filterwarnings("ignore", category=DeprecationWarning)
import initialDataProcessing as idp
#import preprocessDataMaster as aeppdm
from miscFunctions import parseCommandLineNNString
sys.path.insert(0, '../figuresPipeline/')
import heatMapsMaster as hmm
import facetPlottingLibrary as fpl
import singleStainHistogramMaster as sshm
#import exponentialFitMaster as efpm
import isomapMaster as im
import mutualInformationMaster as mi
import singleCellDataProcessing as scdp
sys.path.insert(0, '../neuralNetwork/')

pathToExperimentSpreadsheet = '../../experiments/'
secondPath = '../../output/'

concUnit = 1e9
unitPrefixDictionary = {1e12:'pM',1e9:'nM',1e6:'uM',1e3:'mM',1e0:'M'}
concUnitPrefix = unitPrefixDictionary[concUnit]

modelNames = ['ae-5-3-5-raw','ae-5-3-5-partitioned','ae-5-3-5-raw-normPerObservable','ae-5-3-5-partitioned-normPerObservable'\
        ,'ae-7-3-7-tanh-raw','ae-7-3-7-tanh-partitioned','ae-7-3-7-tanh-raw-normPerObservable','ae-7-3-7-tanh-partitioned-normPerObservable',\
        'ae-11-3-11-parameterized','ae-16-3-16-parameterized','ae-6-3-6-parameterized']
modelName = modelNames[10]

def runPipelinedScript(scriptToRun,inputString,useModifiedDf,cellTypeArray,hmgroup):
    if(int(scriptToRun/100.) == 2): #Neural network scripts
        if(scriptToRun == 201):
            print('Training network for: '+str(inputString))
            import trainAutoEncoderMaster as aetm
            if '/' in inputString:
                trainingStrings = inputString.split('/')
            else:
                trainingStrings = [inputString]
            trainingExperimentNumbersArray = []
            for trainingString in trainingStrings:
                trainingExperimentNumbersArray.append(parseCommandLineNNString(trainingString))
            if('raw' in modelName):
                aetm.trainNetworkNoTimeComponent(secondPath,modelName,trainingStrings,trainingExperimentNumbersArray)
            elif('partitioned' in modelName):
                aetm.trainNetworkNoTimeComponent(secondPath,modelName,trainingStrings,trainingExperimentNumbersArray)
                #aetm.trainNetworkTimePointPartitioning(secondPath,modelName,trainingStrings,trainingExperimentNumbersArray)
            elif('parameterized' in modelName):
                aetm.trainNetworkParameterized(secondPath,modelName,trainingStrings,trainingExperimentNumbersArray)
                #aetm.trainNetworkTimePointPartitioningParameterized(secondPath,modelName,trainingStrings,trainingExperimentNumbersArray)
        elif(scriptToRun == 202):
            print('Crossvalidating network for: '+str(inputString))
            import crossValidateMaster as aecvm
            unsplitTrainingStrings = inputString.split('//')[0]
            if '/' in unsplitTrainingStrings:
                trainingStrings = unsplitTrainingStrings.split('/')
            else:
                trainingStrings = [unsplitTrainingStrings]            
            twoD = False
            crossValidateExperimentNumbers = parseCommandLineNNString(inputString.split('//')[1])
            if('raw' in modelName):
                aecvm.crossValidateNetworkNoTimeComponent(secondPath,modelName,trainingStrings,crossValidateExperimentNumbers,twoD)
            elif('partitioned' in modelName):
                aecvm.crossValidateNetworkNoTimeComponent(secondPath,modelName,trainingStrings,crossValidateExperimentNumbers,twoD)
            elif('parameterized' in modelName):
                aecvm.crossValidateNetworkParameterized(secondPath,modelName,trainingStrings,crossValidateExperimentNumbers,twoD)
        elif(scriptToRun in [203,204,205]): #1D Isomaps or MutualInfo
            print('Generating 1D isomap plot for: '+str(inputString))
            unsplitTrainingStrings = inputString.split('//')[0]
            if '/' in unsplitTrainingStrings:
                trainingStrings = unsplitTrainingStrings.split('/')
            else:
                trainingStrings = [unsplitTrainingStrings]
            crossValidateExperimentNumbers = parseCommandLineNNString(inputString.split('//')[1])
            if(scriptToRun == 203):
                im.create1DIsomapScatterPlot(secondPath,modelName,trainingStrings,crossValidateExperimentNumbers,inputString.split('//')[1])
            elif(scriptToRun == 204):
                im.create1DIsomapDensityPlot(secondPath,modelName,trainingStrings,crossValidateExperimentNumbers,inputString.split('//')[1])
            elif(scriptToRun == 205):
                print(modelName)
                mi.plotMutualInfoMatrix(secondPath,modelName,trainingStrings,crossValidateExperimentNumbers,inputString.split('//')[1])
    else:
        ex_data = pd.read_excel(pathToExperimentSpreadsheet+'masterExperimentSpreadsheet.xlsx')
        experimentsToRun = parseCommandLineNNString(inputString)
        
        for expNum in experimentsToRun:
            expIndex = expNum - 1
            folderName = ex_data['Full Name'][expIndex]
            researcherName = ex_data['Researcher'][expIndex]
            experimentType = ex_data['ExperimentType'][expIndex]
            """
            if(researcherName != 'Sooraj'):
                if scriptToRun != 1:
                    os.chdir('../../../'+researcherName+'/experiments/'+folderName)
                else:
                    os.chdir('../../../'+researcherName+'/experiments')
            else:
                if scriptToRun != 1:
                    os.chdir('../../experiments/'+folderName)
                else:
                    os.chdir('../../experiments')
            """
            os.chdir('../../experiments/'+folderName)
            print(os.getcwd())
            if(int(scriptToRun/100.) == 0): #Initial data processing
                if(scriptToRun == 1):
                    print('Setting up experiment for: '+str(folderName))
                    setUpExperiment.createExperimentFolders(folderName)
                    print(os.getcwd())
                elif(scriptToRun == 2):
                    print('Creating experiment parameters for: '+str(folderName))
                    openOldHeatMap = False
                    setUpExperiment.createParameters(folderName,ex_data['NumConditions'][expIndex],ex_data['NumTimepoints'][expIndex],openOldHeatMap)
                elif(scriptToRun == 3):
                    print('Creating Dataframes for: '+str(folderName))
                    if(cellTypeArray[0]):
                        numberOfCalibrationSamples = ex_data['NumberOfCBAStandardDilutions'][expIndex]
                        initialStandardVolume = ex_data['CBAStandardDilutedVolume'][expIndex]
                        idp.calibrateExperiment(folderName,secondPath,concUnit,concUnitPrefix,numberOfCalibrationSamples,initialStandardVolume)
                        cytdf = idp.createFullDataFrames(folderName,secondPath,expNum,concUnit,concUnitPrefix,'cyt')
                        idp.convertDataFramestoExcel(folderName,secondPath,'cyt',cytdf,useModifiedDf)
                    if(cellTypeArray[1]):
                        celldf = idp.createFullDataFrames(folderName,secondPath,expNum,concUnit,concUnitPrefix,'cell')
                        idp.convertDataFramestoExcel(folderName,secondPath,'cell',celldf,useModifiedDf)
                    if(cellTypeArray[2]):
                        scdp.createInitialSingleCellDataFrame(folderName,expNum)
                        bulkprolifdf = scdp.createProliferationSingleCellDataFrame(folderName,secondPath,useModifiedDf)
                        idp.convertDataFramestoExcel(folderName,secondPath,'prolif',bulkprolifdf,useModifiedDf)
                    if(cellTypeArray[3]):
                        scdp.createInitialSingleCellDataFrame(folderName,expNum)
                        if experimentType not in ['AntibodyTest']:
                            scdp.createCytokineSingleCellDataFrame(folderName)
                            scdp.createCellSingleCellDf(folderName)
                            scdp.createCompleteSingleCellDf(folderName)
                
                elif(scriptToRun == 4): #Preprocess data for neural network
                    print('Preprocessing data for: '+str(folderName))
                    if('raw' in modelName):
                        aeppdm.preprocessDataNoTimeComponent(secondPath,expNum,pickle.load(open('semiProcessedData/modifiedCytokineConcentrationPickleFile-'+folderName+'.pkl','rb')),modelName)
                    elif('partitioned' in modelName):
                        aeppdm.preprocessDataTimePointPartitioning(secondPath,expNum,pickle.load(open('semiProcessedData/modifiedCytokineConcentrationPickleFile-'+folderName+'.pkl','rb')),modelName)
                    elif('parameterized' in modelName):
                        aeppdm.preprocessDataExpFitParameterizing(secondPath,expNum,modelName)
                        #aeppdm.preprocessDataTimePointPartitioningExpFitParameterizing(secondPath,expNum,pickle.load(open('semiProcessedData/modifiedCytokineConcentrationPickleFile-'+folderName+'.pkl','rb')),modelName)
                elif(scriptToRun == 5):
                    print('Creating Exponential Fitting Parameter Dataframe for: '+str(folderName))
                    efpm.createParameterDataFrame(secondPath,expNum,pickle.load(open('semiProcessedData/modifiedCytokineConcentrationPickleFile-'+folderName+'.pkl','rb')))
                    trainingSetExperiments = []
                    for i in range(len(ex_data['Full Name'])):
                        if ex_data['IncludeInParameterTrainingSet'][i] == '*':
                            trainingSetExperiments.append(i+1)
                    efpm.updateTrainingSet(secondPath,trainingSetExperiments)
                if(scriptToRun == 1):
                    os.chdir('../programs/dataProcessing/')
                else:
                    os.chdir('../../programs/dataProcessing/')
            
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
                    singlecelldf = pickle.load(open('semiProcessedData/singleCellDataFrame-complete-'+folderName+modifiedstring+'.pkl','rb'))
                    dataType = 'singlecell'
                    dfArray[dataType] = singlecelldf

                if(scriptToRun == 101): #Heatmaps
                    print('Creating Heatmaps for: '+str(folderName))
                    for dfKey in dfArray:
                        if hmgroup == 'c':
                            hmm.createCombinedHeatMap(folderName,secondPath,dfArray[dfKey],concUnit,concUnitPrefix,useModifiedDf,dfKey)
                        elif hmgroup == 'i':
                            hmm.createIndividualHeatMaps(folderName,secondPath,dfArray[dfKey],concUnit,concUnitPrefix,useModifiedDf,dfKey)
                        else:
                            hmm.createCombinedHeatMap(folderName,secondPath,dfArray[dfKey],concUnit,concUnitPrefix,useModifiedDf,dfKey)
                            hmm.createIndividualHeatMaps(folderName,secondPath,dfArray[dfKey],concUnit,concUnitPrefix,useModifiedDf,dfKey)
                
                elif(scriptToRun in [102,103,104,106,107]): #Facet plots (line,scatter,scatter w/line for expfits)
                    print(dfArray)
                    for dfKey in dfArray:
                        if(scriptToRun == 102):
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
                        elif(scriptToRun == 107):
                            plotType = 'frequency'
                            subPlotType = 'histogram'
                        elif(scriptToRun == 108):
                            plotType = 'frequency'
                            subPlotType = 'kde'

                        print('Creating '+subPlotType+'plots for: '+str(folderName))
                        fpl.facetPlottingGUI(dfArray[dfKey],plotType,dataType)
                        subsettedDfList,subsettedDfListTitles,levelsToPlot,figureLevels,levelValuesPlottedIndividually = fpl.produceSubsettedDataFrames(folderName,secondPath,dfArray[dfKey],useModifiedDf)
                        fpl.plotFacetedFigures(folderName,plotType,subPlotType,dataType,subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually,useModifiedDf)
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
    parser.add_argument("-ep", action='store_true', help="Create input parameters for experiment based on master experiment spreadsheet input data.")
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
    
    parser.add_argument("--input", dest='inputString', help ="Run specified neural network script on these experimental data sets. Separate numbers with , or - for a range, training sets with a /, or a // to separate training from crossval sets.")
    parser.add_argument("--hmgroup", dest='hmgroup', help ="Enter i for individual, c for combined, and b for both types of heatmaps")
    
    args = parser.parse_args()

    hmgroup = 'c'
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
    elif(args.hm):
        scriptToRun = 101
        hmgroup = str(args.hmgroup)
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
    runPipelinedScript(scriptToRun,inputString,args.m,cellTypeArray,hmgroup)
main()

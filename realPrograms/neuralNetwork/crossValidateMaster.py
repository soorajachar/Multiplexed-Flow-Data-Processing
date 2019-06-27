#!/usr/bin/env python3

# coding: utf-8

# # Crossvalidate autoencoder
# Import dependencies
import re
import sys
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from keras.models import Model, load_model
from sklearn import model_selection
from pprint import pprint
from scipy.stats.stats import pearsonr
sys.path.insert(0, '../figuresPipeline/')
from projectionOnNeuronsMaster import projectNeurons2D,projectNeurons3D
from scatterPlotsMaster import createParameterValues, buildLegendHandles 
from realityCheckMaster import createRealityCheckPlots
from mpl_toolkits.mplot3d import Axes3D

# Set parameters
#Order is training set, prediction set, training start, predicition start
def crossValidateNetworkNoTimeComponent(secondPath,modelName,trainingStrings,crossValidateNumbers,twoD):

    for trainingString in trainingStrings:
        print('trainedModels/trainedModel-%s-%s.h5'%(modelName,trainingString))
        ae = load_model(secondPath+'trainedModels/trainedModel-%s-%s.h5'%(modelName,trainingString))
        for crossValidateNumber in crossValidateNumbers:
            print('Crossvalidating '+modelName+' model trained with datasets '+trainingString+' with experiment number: '+str(crossValidateNumber))
            if 'normPerObservable' in modelName:
                dfToCrossValidateComponents = pickle.load(open(secondPath+'preprocessedDataFrames/preprocessedDataFrame-partitioned-normPerObservable-%d.pkl'%(crossValidateNumber),'rb'))
            else:
                dfToCrossValidateComponents = pickle.load(open(secondPath+'preprocessedDataFrames/preprocessedDataFrame-partitioned-%d.pkl'%(crossValidateNumber),'rb'))
            dfToCrossValidateFull = pd.DataFrame(dfToCrossValidateComponents[0],dfToCrossValidateComponents[1],dfToCrossValidateComponents[2])
            if('raw' not in modelName):
                dfToCrossValidate = dfToCrossValidateFull.xs('A',level='Partition')
            else:
                print('raw')
                dfToCrossValidate = dfToCrossValidateFull
            if('5-3-5' in modelName):
                idx = pd.IndexSlice
                dfToCrossValidate = dfToCrossValidate.loc[idx[:],idx[['IFNg','IL-2','IL-6','IL-17A','TNFa']]]
            X_norm = dfToCrossValidate.values
            print(dfToCrossValidate)

            # # Crossvalidate autoencoders trained on different data
            encoder = Model(ae.input,ae.get_layer('encoder').output)
            Z_enc = encoder.predict(X_norm)
            decoder = ae.get_layer(name='decoder')
            R_enc = ae.predict(X_norm)
            #createRealityCheckPlots(secondPath,crossValidateNumber,modelName,trainingString,dfToCrossValidate,R_enc,X_norm)
            
            #Project dataset on previously trained neural network
            neuronList = ['Node 1','Node 2','Node 3']
            crossValidatedDf = pd.DataFrame(Z_enc,dfToCrossValidate.index,neuronList)
            if(twoD):
                print('\tSaving 2d neuron projection...')
                projectNeurons2D(secondPath,crossValidatedDf,modelName,trainingString,crossValidateNumber)
            else:
                print('\tSaving 3d neuron projection...')
                projectNeurons3D(secondPath,crossValidatedDf,modelName,trainingString,crossValidateNumber)
            print('\tSaving cross validated network...')
            with open(secondPath+'crossValidatedNetworks/crossValidatedNetwork-%d-on-%s-%s.pkl'%(crossValidateNumber,modelName,trainingString), "wb") as f:
                pickle.dump(crossValidatedDf, f)

def crossValidateNetworkPartitioned(secondPath,modelName,trainingStrings,crossValidateNumbers,twoD):

    for trainingString in trainingStrings:
        ae = load_model(secondPath+'trainedModels/trainedModel-%s-%s.h5'%(modelName,trainingString))
        for crossValidateNumber in crossValidateNumbers:
            print('Crossvalidating '+modelName+' model trained with datasets '+trainingString+' with experiment number: '+str(crossValidateNumber))
            dfToCrossValidateComponents = pickle.load(open(secondPath+'preprocessedDataFrames/preprocessedDataFrame-partitioned-%d.pkl'%(crossValidateNumber),'rb'))
            dfToCrossValidateFull = pd.DataFrame(dfToCrossValidateComponents[0],dfToCrossValidateComponents[1],dfToCrossValidateComponents[2])
            partitions = pd.unique(trainingDf.index.get_level_values(dfToCrossValidateFull.index.names[-1]))
            for partition in partitions:
                dfToCrossValidate = dfToCrossValidateFull.xs(partition,level='Partition')
            X_norm = dfToCrossValidate.values
            
            # # Crossvalidate autoencoders trained on different data
            encoder = Model(ae.input,ae.get_layer('encoder').output)
            Z_enc = encoder.predict(X_norm)
            decoder = ae.get_layer(name='decoder')
            R_enc = ae.predict(X_norm)
            createRealityCheckPlots(secondPath,crossValidateNumber,modelName,trainingString,dfToCrossValidate,R_enc,X_norm)
            
            #Project dataset on previously trained neural network
            neuronList = ['Node 1','Node 2','Node 3']
            crossValidatedDf = pd.DataFrame(Z_enc,dfToCrossValidate.index,neuronList)
            if(twoD):
                print('\tSaving 2d neuron projection...')
                projectNeurons2D(secondPath,crossValidatedDf,modelName,trainingString,crossValidateNumber)
            else:
                print('\tSaving 3d neuron projection...')
                projectNeurons3D(secondPath,crossValidatedDf,modelName,trainingString,crossValidateNumber)
            print('\tSaving cross validated network...')
            with open(secondPath+'crossValidatedNetworks/crossValidatedNetwork-%d-on-%s-%s.pkl'%(crossValidateNumber,modelName,trainingString), "wb") as f:
                pickle.dump(crossValidatedDf, f)

def crossValidateNetworkParameterized(secondPath,modelName,trainingStrings,crossValidateNumbers,twoD):

    for trainingString in trainingStrings:
        ae = load_model(secondPath+'trainedModels/trainedModel-%s-%s.h5'%(modelName,trainingString))
        for crossValidateNumber in crossValidateNumbers:
            print('Crossvalidating '+modelName+' model trained with datasets '+trainingString+' with experiment number: '+str(crossValidateNumber))
            dfToCrossValidateComponents = pickle.load(open(secondPath+'preprocessedDataFrames/preprocessedDataFrame-parameterized-%d.pkl'%(crossValidateNumber),'rb'))
            dfToCrossValidateFull = pd.DataFrame(dfToCrossValidateComponents[0],dfToCrossValidateComponents[1],dfToCrossValidateComponents[2])
            dfToCrossValidate = dfToCrossValidateFull
            
            X_norm = dfToCrossValidate.values
            
            # # Crossvalidate autoencoders trained on different data
            encoder = Model(ae.input,ae.get_layer('encoder').output)
            Z_enc = encoder.predict(X_norm)
            decoder = ae.get_layer(name='decoder')
            R_enc = ae.predict(X_norm)
            #createRealityCheckPlots(secondPath,crossValidateNumber,modelName,trainingString,dfToCrossValidate,R_enc,X_norm)
            
            #Project dataset on previously trained neural network
            neuronList = ['Node 1','Node 2','Node 3']
            crossValidatedDf = pd.DataFrame(Z_enc,dfToCrossValidate.index,neuronList)
            if(twoD):
                print('\tSaving 2d neuron projection...')
                projectNeurons2D(secondPath,crossValidatedDf,modelName,trainingString,crossValidateNumber)
            else:
                print('\tSaving 3d neuron projection...')
                projectNeurons3D(secondPath,crossValidatedDf,modelName,trainingString,crossValidateNumber)
            print('\tSaving cross validated network...')
            with open(secondPath+'crossValidatedNetworks/crossValidatedNetwork-%d-on-%s-%s.pkl'%(crossValidateNumber,modelName,trainingString), "wb") as f:
                pickle.dump(crossValidatedDf, f)

def crossValidateNetworkPartitionedParameterized(secondPath,modelName,trainingStrings,crossValidateNumbers,twoD):

    for trainingString in trainingStrings:
        ae = load_model(secondPath+'trainedModels/trainedModel-%s-%s.h5'%(modelName,trainingString))
        for crossValidateNumber in crossValidateNumbers:
            print('Crossvalidating '+modelName+' model trained with datasets '+trainingString+' with experiment number: '+str(crossValidateNumber))
            dfToCrossValidateComponents = pickle.load(open(secondPath+'preprocessedDataFrames/preprocessedDataFrame-partitioned-parameterized-%d.pkl'%(crossValidateNumber),'rb'))
            dfToCrossValidateFull = pd.DataFrame(dfToCrossValidateComponents[0],dfToCrossValidateComponents[1],dfToCrossValidateComponents[2])
            dfToCrossValidate = dfToCrossValidateFull
            
            X_norm = dfToCrossValidate.values
            
            # # Crossvalidate autoencoders trained on different data
            encoder = Model(ae.input,ae.get_layer('encoder').output)
            Z_enc = encoder.predict(X_norm)
            decoder = ae.get_layer(name='decoder')
            R_enc = ae.predict(X_norm)
            #createRealityCheckPlots(secondPath,crossValidateNumber,modelName,trainingString,dfToCrossValidate,R_enc,X_norm)
            
            #Project dataset on previously trained neural network
            neuronList = ['Node 1','Node 2','Node 3']
            crossValidatedDf = pd.DataFrame(Z_enc,dfToCrossValidate.index,neuronList)
            if(twoD):
                print('\tSaving 2d neuron projection...')
                projectNeurons2D(secondPath,crossValidatedDf,modelName,trainingString,crossValidateNumber)
            else:
                print('\tSaving 3d neuron projection...')
                projectNeurons3D(secondPath,crossValidatedDf,modelName,trainingString,crossValidateNumber)
            print('\tSaving cross validated network...')
            with open(secondPath+'crossValidatedNetworks/crossValidatedNetwork-%d-on-%s-%s.pkl'%(crossValidateNumber,modelName,trainingString), "wb") as f:
                pickle.dump(crossValidatedDf, f)

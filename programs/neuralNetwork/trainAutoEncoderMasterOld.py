#!/usr/bin/env python3

# coding: utf-8

# # Train 7-3-7 lin-tanh-lin autoencoder
# First, we import dependencies. Then we read normalized data into $X\_norm$, we split the data in training and test data and train the autoencoder

# In[ ]:
import pickle
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import keras.layers as ks
from keras.utils import plot_model 
from keras.models import Sequential,Model
from keras.layers import Dense,Concatenate,Input
from keras.optimizers import Adam
from sklearn import model_selection

def trainNetworkNoTimeComponent(secondPath,modelName,trainingStrings,trainingNumbersArray):

    for trainingString,trainingNumbers in zip(trainingStrings,trainingNumbersArray):
        print('Training '+modelName+' model with training string: '+str(trainingString))
        trainingSetDict = {}
        trainingDfContainer = []
        for trainingDataSet in trainingNumbers:
            trainingDataFrameComponents = pickle.load(open(secondPath+'preprocessedDataFrames/preprocessedDataFrame-partitioned-%d.pkl'%(trainingDataSet),'rb'))
            trainingDf = pd.DataFrame(trainingDataFrameComponents[0],trainingDataFrameComponents[1],trainingDataFrameComponents[2])
            trainingDataFrameComponents2 = pickle.load(open(secondPath+'preprocessedDataFrames/preprocessedDataFrame-%d.pkl'%(trainingDataSet),'rb'))
            trainingDf2 = pd.DataFrame(trainingDataFrameComponents2[0],trainingDataFrameComponents2[1],trainingDataFrameComponents2[2])
            trainingDfContainer.append(trainingDf)
            partitions = pd.unique(trainingDf.index.get_level_values(trainingDf.index.names[-1]))
            for partition in partitions:
                partitionedDf = trainingDf.xs(partition,level='Partition')
                if(partition not in list(trainingSetDict.keys())):
                    trainingSetDict[partition] = partitionedDf.values
                else:
                    trainingSetDict[partition] = np.concatenate([trainingSetDict[partition],partitionedDf.values])
        if('raw' not in modelName): 
            X_norm = trainingSetDict['A']
        else:
            X_norm = trainingDfContainer[0]

        X_train, X_test = model_selection.train_test_split(X_norm, test_size=0.1)

        ae = Sequential()
        ae.add(Dense(3, activation='tanh', name = "encoder", input_shape = (X_norm.shape[1],)))
        ae.add(Dense(X_norm.shape[1], activation='linear', name = "decoder"))
        ae.compile(loss='mean_squared_error', optimizer = Adam())
        history = ae.fit(X_train, X_train, shuffle = True, batch_size=32, verbose = 0,epochs=2000, validation_data=(X_test, X_test))
        ae.save(secondPath+"trainedModels/trainedModel-%s-%s.h5"%(modelName,trainingString))
        plt.plot(history.history['loss'], label = 'train')
        plt.plot(history.history['val_loss'], label = 'validation')
        plt.xlabel('epoch')
        plt.ylabel('loss')
        plt.legend()
        plt.savefig(secondPath+'lossFigures/lossFigure-'+modelName+'-%s.png'%(trainingString))
        plt.clf()

def trainNetworkTimePointPartitioning(secondPath,trainingStrings,trainingNumbersArray):

    modelName = 'ae-14-3-14-tanh'
    for trainingString,trainingNumbers in zip(trainingStrings,trainingNumbersArray):
        print('Training model with training string: '+str(trainingString))
        trainingSetDict = {}
        for trainingDataSet in trainingNumbers:
            trainingDataFrameComponents = pickle.load(open(secondPath+'preprocessedDataFrames/preprocessedDataFrame-partitioned-%d.pkl'%(trainingDataSet),'rb'))
            trainingDf = pd.DataFrame(trainingDataFrameComponents[0],trainingDataFrameComponents[1],trainingDataFrameComponents[2])
            partitions = pd.unique(trainingDf.index.get_level_values(trainingDf.index.names[-1]))
            for partition in partitions:
                partitionedDf = trainingDf.xs(partition,level='Partition')
                if(partition not in list(trainingSetDict.keys())):
                    trainingSetDict[partition] = partitionedDf.values
                else:
                    trainingSetDict[partition] = np.concatenate([trainingSetDict[partition],partitionedDf.values])
        maxLength = 0
        numPartitions = len(trainingSetDict)
        numObservables = trainingSetDict[list(trainingSetDict.keys())[0]].shape[1] 
        for key in trainingSetDict:
            trainingSet = trainingSetDict[key]
            if maxLength < trainingSet.shape[0]:
                maxLength = trainingSet.shape[0]
        for key in trainingSetDict:
            trainingSet = trainingSetDict[key]
            zeroArray = np.zeros((maxLength,numObservables))
            zeroArray[:trainingSet.shape[0],:trainingSet.shape[1]] = trainingSet
            trainingSetDict[key] = zeroArray
        for key in trainingSetDict:
            fullMatrix = 
        print(np.hstack((trainingSetDict['A'],trainingSetDict['B'])))
        sys.exit(0)
        ae = Sequential()
        ae.add(Dense(3, activation='tanh', name = "encoder", input_shape = (numPartitions*numObservables,)))
        ae.add(Dense(numPartitions*numObservables, activation='linear', name = "decoder"))
        
        """
        partitionedNeuronList = []
        X_train_list = []
        X_test_list = []
        partitionedInputList = []
        
        for key in trainingSetDict:
            layerName = 'inputPartition'+str(key)
            partitionedInput = Input(shape = (numObservables,),name = layerName)
            partitionedInputList.append(partitionedInput)
            X_norm_temp = trainingSetDict[key]
            X_train_temp, X_test_temp = model_selection.train_test_split(X_norm_temp, test_size=0.1)
            X_train_list.append(X_train_temp)
            X_test_list.append(X_test_temp)
        
        combinedInput = ks.concatenate(partitionedInputList)
        x = Dense(3, activation='tanh')(combinedInput)
        partitionedOutputList = []
        for key in trainingSetDict:
            layerName = 'outputPartition'+str(key)
            partitionedOutput = Dense(numObservables,activation='linear',name=layerName)(x)
            partitionedOutputList.append(partitionedOutput)
        ae = Model(inputs=partitionedInputList, outputs=partitionedOutputList)
        """
        """
        plot_model(ae, to_file='model.png')
        ae.summary()
        ae.compile(loss='mean_squared_error', optimizer = Adam())
        history = ae.fit(X_train_list, X_train_list, shuffle = True, batch_size=32, verbose = 0,epochs=2000, validation_data=(X_test_list, X_test_list))
        ae.save(secondPath+"trainedModels/trainedModel-"+modelName+"-%s.h5"%(trainingString))
        plt.plot(history.history['loss'], label = 'train')
        plt.plot(history.history['val_loss'], label = 'validation')
        plt.xlabel('epoch')
        plt.ylabel('loss')
        plt.legend()
        plt.savefig(secondPath+'lossFigures/lossFigure-'+modelName+'-%s.png'%(trainingString))
        plt.clf()
        """

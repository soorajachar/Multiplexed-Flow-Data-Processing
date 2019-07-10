#!/usr/bin/env python3
"""
Created on Mon Nov 19 16:34:25 2018
 
@author: sachar
"""
import sys
import pickle 
from sklearn import manifold
import numpy as np
import matplotlib.pyplot as plt
from itertools import product,combinations
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd
from itertools import chain

#Creates 1D isomap
def create1DIsomap(secondPath,crossValidateNumber,modelName,trainingString):
    #Load in saved dataframe of a crossvalidating dataset projected onto a training set
    crossValidatedDf = pickle.load(open(secondPath+'crossValidatedNetworks/crossValidatedNetwork-%d-on-%s-%s.pkl'%(crossValidateNumber,modelName,trainingString),'rb'))
    #Only compare OT1 peptides in catplot; remove CD28/IL2 added peptides from 68 and other TCellTypes (P14,F5) from 61,64,53
    if(crossValidateNumber == 68):
        crossValidatedDf = crossValidatedDf.xs(['None'],level=['Antibody'])
    elif(crossValidateNumber in [61,64,53]):
        crossValidatedDf = crossValidatedDf.xs(['OT1'],level=['TCellType'])
    #Get projection values across all 3 nodes of autoencodder
    data=np.array([crossValidatedDf.iloc[:,0],crossValidatedDf.iloc[:,1],crossValidatedDf.iloc[:,2]]).T
    #Use the isomap package to map 3d neuron projection of data onto a 1D space
    trans_data = manifold.Isomap(n_neighbors=30, n_components=1, n_jobs=-1).fit_transform(data)
    #Assign 1D isomap values into the same dataframe as the original data
    trans_data_df = pd.DataFrame(trans_data,crossValidatedDf.index)
    #Make sure all isomap projections progress from positive to negative values as peptide affinity decreases (with N4 1uM always being at the highest positive value)
    if(trans_data_df.iloc[0,0] < 0):
        trans_data_df *= -1

    #Unravel dataframe to facetgrid format (include multindex levels as actual values in dataframe) to prepare for facetgrid plots
    peptideVals = trans_data_df.index.get_level_values('Peptide').ravel()
    concentrationVals = trans_data_df.index.get_level_values('Concentration').ravel()
    isomapVals = trans_data_df.values.ravel()
    #New Dataset column allows datasets to be visualized separately
    dataSetVals = np.repeat(crossValidateNumber,len(peptideVals))
    data = {'Dataset':dataSetVals,'Peptide':peptideVals,'Concentration':concentrationVals,'Isomap Value':isomapVals}
    trans_data_df2 = pd.DataFrame(data,columns = list(data.keys()))
    #Save wide and long versions of isomap dataframes
    with open(secondPath+'isomapPickleFiles/isomapPickleFile-%d-on-%s-%s-long.pkl'%(crossValidateNumber,modelName,trainingString), "wb") as f:
        pickle.dump(trans_data_df, f)
    with open(secondPath+'isomapPickleFiles/isomapPickleFile-%d-on-%s-%s-wide.pkl'%(crossValidateNumber,modelName,trainingString), "wb") as f:
        pickle.dump(trans_data_df2, f)

    return trans_data_df2

#Step through each crossvalidating dataset, grab isomap values, convert to facetgrid format, then concatenate and return full dataframe and new color cycle
def concatenateIsomapDataFrames(secondPath,modelName,trainingString,crossValidateNumbers):
    all_trans_data = [] 
    for cv,index in zip(crossValidateNumbers,range(len(crossValidateNumbers))):
        trans_data_df2 = create1DIsomap(secondPath,cv,modelName,trainingString)
        all_trans_data.append(trans_data_df2)
        
    all_trans_data_df = pd.concat(all_trans_data)
    #Set color cycle to standard scale used for other plots
    peptideVals = pd.unique(all_trans_data_df['Peptide'])
    sns.set_palette(sns.color_palette("hls",len(peptideVals)))
    colorVals = sns.color_palette()
    colorVals.reverse()
    
    return all_trans_data_df,colorVals

#Creates scatter plots of isomap values, with different colors representing concentrations
def create1DIsomapScatterPlot(secondPath,modelName,trainingStrings,crossValidateNumbers,crossValidateString):
    for trainingString in trainingStrings:
        all_trans_data_df,colorVals = concatenateIsomapDataFrames(secondPath,modelName,trainingString,crossValidateNumbers)
        
        fig = plt.figure(figsize = (5*len(crossValidateNumbers),15))
        sns.catplot(x='Peptide',y='Isomap Value',hue='Concentration',col='Dataset',data=all_trans_data_df,palette=colorVals)
        plt.savefig(secondPath+'isomapPlots/scatterPlots/%s/%s/isomapPlot-%s-on-%s-%s-scatter.png'%(modelName,trainingString,crossValidateString,modelName,trainingString))
        plt.clf()

#Creates density plots of isomap values, with different colors representing different peptides
def create1DIsomapDensityPlot(secondPath,modelName,trainingStrings,crossValidateNumbers,crossValidateString):
    for trainingString in trainingStrings:
        all_trans_data_df,colorVals = concatenateIsomapDataFrames(secondPath,modelName,trainingString,crossValidateNumbers)

        fig = plt.figure(figsize = (5*len(crossValidateNumbers),15))
        g = sns.FacetGrid(all_trans_data_df, hue='Peptide',row='Dataset',palette=colorVals)
        g.map(sns.kdeplot,'Isomap Value',shade=True)
        g.add_legend()
        plt.savefig(secondPath+'isomapPlots/densityPlots/%s/%s/isomapPlot-row-%s-on-%s-%s-density.png'%(modelName,trainingString,crossValidateString,modelName,trainingString))
        plt.clf()

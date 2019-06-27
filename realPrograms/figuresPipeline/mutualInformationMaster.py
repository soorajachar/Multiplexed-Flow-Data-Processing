#!/usr/bin/env python3
"""
Created on Mon Nov 19 16:34:25 2018
 
@author: sachar
"""
import sys
import pickle 
from sklearn import manifold,metrics
import numpy as np
import matplotlib.pyplot as plt
from itertools import product,combinations
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd
from itertools import chain
from scipy import stats

#https://stackoverflow.com/questions/20491028/optimal-way-to-compute-pairwise-mutual-information-using-numpy (last part of 1st answer)
#Return mutual information between two lists
def calc_MI(x, y, bins):
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = metrics.mutual_info_score(None, None, contingency=c_xy)
    return mi

def compute_mutual_info(pdf_a,pdf_b):
    mutual_info = 0
    for i in range(len(pdf_a)):
        ab = 0.5*(pdf_a[i]+pdf_b[i])
        if(ab == 0):
            a = 0
            b = 0
        else:
            if pdf_a[i] == 0:
                a = 0;
            else:
                a = pdf_a[i] * np.log2(pdf_a[i]/ab)
            if pdf_b[i] == 0:
                b = 0;
            else:
                b = pdf_b[i] * np.log2(pdf_b[i]/ab)
        mutual_info += 0.5*(a+b)
        
    return mutual_info/np.sum(pdf_a)

#Compute kd estimations of pdf for each pairwise comparison dataframe, then use to construct a histogram (more accurate than using histograms initially due to low sample size for certain peptide/concentration measurements)
def computePairwiseMutualInfo(df_X,df_Y,linspaceRange):
    #kdeValues = np.linspace(-1,1,num=200)
    kdeValues = np.linspace(linspaceRange[0],linspaceRange[1],num=(linspaceRange[1]-linspaceRange[0])*100)
    kernelX = stats.gaussian_kde(df_X)
    c_X = kernelX(kdeValues)
    kernelY = stats.gaussian_kde(df_Y)
    c_Y = kernelY(kdeValues)
    #mutualInfo = calc_MI(c_X,c_Y,bins)
    mutualInfo = compute_mutual_info(c_X,c_Y)
    return mutualInfo

#Compute mutual information scores for every peptide/concentration cluster in a crossvalidated dataset from the 1D isomap values
def plotMutualInfoMatrix(secondPath,modelName,trainingStrings,crossValidateNumbers,crossValidateString):
    for trainingString in trainingStrings:
        for crossValidateNumber in crossValidateNumbers:
            #Load in isomap dataframe
            longIsomapDf = pickle.load(open(secondPath+'isomapPickleFiles/isomapPickleFile-%d-on-%s-%s-long.pkl'%(crossValidateNumber,modelName,trainingString),'rb'))
            
            #Unstack the dataframe at the time level in order to access the index of this new dataframe, which gives all peptide/concentration combinations
            peptideConcCombos = longIsomapDf.unstack(-1)
            #Compute mutual information between isomap values of each peptide and concentration combination (self values and lower part of symetric matrix are excluded by initializing with Nans and never filling over them)
            totalMutualInfo = np.empty((peptideConcCombos.shape[0],peptideConcCombos.shape[0],))
            totalMutualInfo[:] = np.nan
            for i in range(totalMutualInfo.shape[0]):
                for j in range(totalMutualInfo.shape[1]):
                    if(i > j):
                        peptideConcCombo = peptideConcCombos.iloc[i,:].name
                        peptideConcCombo2 = peptideConcCombos.iloc[j,:].name
                        levelNames = peptideConcCombos.index.names
                        totalMutualInfo[i,j] = computePairwiseMutualInfo(peptideConcCombos.xs(peptideConcCombo,level=levelNames).values,peptideConcCombos.xs(peptideConcCombo2,level=levelNames).values,[-1,1])
            #Dataframe constructed with peptide/concentration combos along columns and row of dataframes (symetric matrix charting mutual info between different clusters)
            peptideConcentrationSeparatedMutualInfoDf = pd.DataFrame(totalMutualInfo,index=peptideConcCombos.index,columns=peptideConcCombos.index)
            peptideVals = pd.unique(peptideConcentrationSeparatedMutualInfoDf.index.get_level_values(0))
            print(peptideConcentrationSeparatedMutualInfoDf)

            #Compute mutual information between each peptide cluster as a whole (including all concentrations) by averaging the mutual informations of N4-Q4 comparisons with any concentration combination
            #These values reflect how well two peptides of different qualities are differentiated between by the autoencoder's mapping
            peptideSeparatedMutualInfoDf = pd.DataFrame(np.empty((len(peptideVals),len(peptideVals),)),index=peptideVals,columns=peptideVals)
            peptideSeparatedMutualInfoDf[:] = np.nan
            for i in range(peptideSeparatedMutualInfoDf.shape[0]):
                for j in range(peptideSeparatedMutualInfoDf.shape[0]):
                    if(i > j):
                        peptideGroup1 = peptideVals[i]
                        peptideGroup2 = peptideVals[j]
                        peptideGroupInteractions = peptideConcentrationSeparatedMutualInfoDf.loc[peptideGroup1][peptideGroup2]
                        peptideSeparatedMutualInfoDf.iloc[i,j] = np.nanmean(peptideGroupInteractions.values)
            print(peptideSeparatedMutualInfoDf)
            
            #Compute average mutual information among a single peptide's different concentrations
            #These values reflect how well the autoencoder separates the different concentrations of a particular peptide
            multipleConcentrationPeptideVals = []
            for i in range(len(peptideVals)):
                if(peptideConcentrationSeparatedMutualInfoDf[peptideVals[i]].shape[1] > 1):
                    multipleConcentrationPeptideVals.append(peptideVals[i])
            peptideSelfValsDict = {}
            for multipleConcentrationPeptideVal in multipleConcentrationPeptideVals:
                peptideSelfInteractions = peptideConcentrationSeparatedMutualInfoDf.loc[multipleConcentrationPeptideVal][multipleConcentrationPeptideVal]
                peptideSelfValsDict[str(multipleConcentrationPeptideVal)+'-'+str(multipleConcentrationPeptideVal)] = np.nanmean(peptideSelfInteractions.values)
            peptideSeparatedSelfMutualInfoDf = pd.DataFrame.from_dict(peptideSelfValsDict,columns=['Mutual Info'],orient='index')

            #Compute overall mutual information average for peptide separated mutual information values (excluding within peptide group (self) mutual information values for peptides with multiple concentrations)
            #These values provide a metric for us to judge how well the autoencoder differentiates between different quality antigens as a whole (main purpose: separate quality from quanity)
            #Can be used as a performance metric when attempting to understand if a particular architecture choice does better on our goal (closer to 0 mutual information is better: less overlap between distributions)
            overallMutualInfoMetric = np.nanmean(peptideSeparatedMutualInfoDf.values)

            #Dumpy all 3 mutual information dataframes and overall mutual information metric in a 4 element list pickle file
            allMutualInfoData = [peptideConcentrationSeparatedMutualInfoDf,peptideSeparatedMutualInfoDf,peptideSeparatedSelfMutualInfoDf,overallMutualInfoMetric]
            with open(secondPath+'mutualInfoPickleFiles/mutualInfoPickleFile-pc_p_ps_o-%d-on-%s-%s.pkl'%(crossValidateNumber,modelName,trainingString),'wb') as f:
                pickle.dump(allMutualInfoData, f)

            #Save all mutual information dataframes as figures (NOTE: STILL NEED TO ADD IN MULTIINDEX LABELING ON HEATMAP AXES (LIKE CBA HEATMAPS))
            midfNames = ['peptideConcentrationMIHeatMap','peptideMIHeatMap','peptideSelfMIHeatMap']
            for midf,midfName in zip(allMutualInfoData,midfNames):
                #sns.set(font_scale=1.0)
                fig = plt.figure(num=1,dpi=150,facecolor='w',edgecolor='k')
                ax = fig.add_subplot(111)
                cax = sns.heatmap(midf, cmap='coolwarm_r', vmin=0,vmax=1, center=0.5, linewidths=.5)
                cax.collections[0].colorbar.set_label('Mutual Information',labelpad=20)
                plt.title(crossValidateNumber)
                plt.savefig(secondPath+'mutualInfoHeatmaps/%s/%s/%s/%s-%d-on-%s-%s.png'%(midfName+'s',modelName,trainingString,midfName,crossValidateNumber,modelName,trainingString),bbox_inches='tight')
                plt.clf()
            print(allMutualInfoData[-1])

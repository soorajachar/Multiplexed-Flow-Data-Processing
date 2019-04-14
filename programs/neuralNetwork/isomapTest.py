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
sys.path.insert(0, '../dataprocessing/')
sys.path.insert(0, '../figuresPipeline/')
from facetPlottingLibrary import createParameterValues,buildLegendHandles,returnMarkerVals,loadInDataFrameNoUI
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd
from itertools import chain

fig = plt.figure(figsize = (15,15))
cvdfs = [28,30,67,68,61,64,53]
for cv,index in zip(cvdfs,range(len(cvdfs))):
    crossValidatedDf = pickle.load(open('../../output/crossValidatedNetworks/crossValidatedNetwork-%d-on-ae-7-3-7-tanh-partitioned-28,30,67.pkl'%(cv),'rb'))
    if(cv == 68):
        crossValidatedDf = crossValidatedDf.xs(['None'],level=['Antibody'])
    elif(cv in [61,64,53]):
        crossValidatedDf = crossValidatedDf.xs(['OT1'],level=['TCellType'])
    legendParameterToLevelNameDict = loadInDataFrameNoUI(crossValidatedDf) 
    nonlegendParameterToLevelNameDict = {}

    lineplot = False
    parameterToParameterValsDict,_,_,_,maxNumConditions = createParameterValues(crossValidatedDf,legendParameterToLevelNameDict,nonlegendParameterToLevelNameDict,lineplot,True) 
    legendHandles = buildLegendHandles(crossValidatedDf,parameterToParameterValsDict,legendParameterToLevelNameDict,nonlegendParameterToLevelNameDict,maxNumConditions,lineplot)

    colorVals = parameterToParameterValsDict['c']
    data=np.array([crossValidatedDf.iloc[:,0],crossValidatedDf.iloc[:,1],crossValidatedDf.iloc[:,2]]).T
    if(cv == 28):
        trans_data = manifold.Isomap(n_neighbors=30, n_components=1, n_jobs=-1).fit_transform(data)
    else:
        trans_data = manifold.Isomap(n_neighbors=20, n_components=1, n_jobs=-1).fit_transform(data)
    trans_data_df = pd.DataFrame(trans_data,crossValidatedDf.index)
    if(trans_data_df.iloc[0,0] < 0):
        trans_data_df *= -1
    peptideVals = trans_data_df.index.get_level_values('Peptide').ravel()
    concentrationVals = trans_data_df.index.get_level_values('Concentration').ravel()
    #timeVals = trans_data_df.index.get_level_values('Time').ravel()
    isomapVals = trans_data_df.values.ravel()
    data = {'Peptide':peptideVals,'Concentration':concentrationVals,'IsomapVal':isomapVals}
    trans_data_df2 = pd.DataFrame(data,columns = ['Peptide','Concentration','IsomapVal'])
    #Plot dataframe row by row (to allow timepoints to be reused as x values)
    ax = fig.add_subplot(4,4,index+1)
    #sns.set_palette(sns.color_palette("hls",len(pd.unique(trans_data_df.index.get_level_values('Peptide')))))
    #colors = sns.color_palette()
    #colors.reverse()
    #sizes = np.linspace(320,20,num=numCurrentLevelValues) #Default value is 80 in s units (path collections object used for scatter plots)
    sns.stripplot(x='Peptide',y='IsomapVal',hue='Concentration',data=trans_data_df2)#,palette=colors)
    """
    for peptideVal in pd.unique(trans_data_df.index.get_level_values('Peptide')):
        colorVal = colorVals[peptideVal]
        x = list(chain.from_iterable(trans_data_df.loc[peptideVal].values))
        #sns.distplot(x,color=colorVal)
        sns.kdeplot(x,shade=True,color=colorVal)
    """
    plt.title(str(cv))
#fig.legend(handles = legendHandles,ncol=len(parameterToParameterValsDict))
plt.show()

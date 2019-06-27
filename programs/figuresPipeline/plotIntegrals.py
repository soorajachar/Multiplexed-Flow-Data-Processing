#!/usr/bin/env python3 
import pandas as pd
import seaborn as sns
import numpy as np
import pickle,math,sys,re
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.integrate import simps
from scipy.stats import linregress

#['1uM' '1nM' '100pM' '10pM' '10nM' '100nM']
unitPrefixDictionary = {'fM':1e-15,'pM':1e-12,'nM':1e-9,'uM':1e-6,'mM':1e-3,'M':1e0,'':0}
def sortSINumerically(listSI,sort,descending):
    numericList = []
    for unitString in listSI:
        splitString = re.split('(\d+)',unitString)
        numericList.append(float(splitString[1])*float(unitPrefixDictionary[splitString[2]]))
    originalNumericList = numericList.copy()
    if sort:
        numericList.sort(reverse=descending)
    numericIndices = []
    for elem in numericList:
        numericIndices.append(originalNumericList.index(elem))
    sortedListSI = []
    for elem in numericIndices:
        sortedListSI.append(listSI[elem])
    return sortedListSI,numericList

def r_squared(xdata,ydata,func,popt):
    residuals = ydata- func(xdata, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = round(1 - (ss_res / ss_tot),3)
    return r_squared

#2 parameter (vshift fixed per cytokine based on lower LOD of cytokine): y = A(1-e^(-tau*x))
def boundedExponential(x, amplitude,tau,vshift):
    return amplitude*(np.subtract(1,np.exp(np.multiply(-1,np.multiply(x,tau)))))+vshift

#5 parameter (vshift fixed per cytokine; based on lower LOD of cytokine): y = A((1/(1+e^(-tau1*(x-td1))))-(1/(1+e^(-tau2*(x-td2)))))
def logisticDoubleExponential(x,amplitude,tau1,tau2,timedelay1,timedelay2,vshift):
    return amplitude*np.subtract(np.divide(1,np.add(1,np.exp(np.multiply(-1*tau1, np.subtract(x,timedelay1))))),np.divide(1,np.add(1,np.exp(np.multiply(-1*tau2,np.subtract(x,timedelay2))))))+vshift

def findFits(x,logy,cytokineMin,fitName):
    epsilon = 0.0001
    if fitName == 'boundedExponential':
        logy = np.log10(logy)
        cytokineMin = np.log10(cytokineMin)
        lbounds = [(max(logy)-cytokineMin)/2,0,cytokineMin-epsilon]
        ubounds = [(max(logy)-cytokineMin)*2,10,cytokineMin+epsilon]
        try:
            popt,pcov = curve_fit(boundedExponential,x,logy,sigma=logy,bounds=(lbounds,ubounds))
            rsquared = round(r_squared(x,logy,boundedExponential,popt),3)
        except:
            return [0,0,cytokineMin],0.0
    else:
        logy = np.log10(logy)
        cytokineMin = np.log10(cytokineMin)
        #lbounds = [(max(logy)-cytokineMin)-epsilon,0,0,-100,-100,cytokineMin-epsilon]
        #ubounds = [(max(logy)-cytokineMin)+epsilon,5,5,100,100,cytokineMin+epsilon]
        lbounds = [(max(logy)-cytokineMin)-epsilon,0,0,0,x[logy.tolist().index(max(logy))],cytokineMin-epsilon]
        ubounds = [(max(logy)-cytokineMin)+epsilon,2,2,x[logy.tolist().index(max(logy))],max(x),cytokineMin+epsilon]
        trueX = [0]+x
        trueY = [cytokineMin]+list(logy)
        try:
            #popt,pcov = curve_fit(logisticDoubleExponential,trueX,trueY,sigma=trueY,bounds=(lbounds,ubounds))
            popt,pcov = curve_fit(logisticDoubleExponential,x,logy,sigma=logy,bounds=(lbounds,ubounds))
            rsquared = round(r_squared(x,logy,logisticDoubleExponential,popt),3)
        except:
            return [0,0,0,0,0,cytokineMin],0.0
    return popt,rsquared

def createParameterDataFrame(df):
    x = list(df.columns)
    parameterizedObservables = ['IFNg','TNFa','IL-6','IL-2','IL-17A']
    parameterizedFits = ['boundedExponential','logisticDoubleExponential']
    parameters = [['r2','A','k0','v'],['r2','A','k0','k1','x0','x1','v']]
    dfListOfTuples = []
    for i in range(len(parameterizedFits)):
        fitName = parameterizedFits[i]
        parameterList = parameters[i]
        for parameter in parameterList:
            dfListOfTuples.append((fitName,parameter))
    columnIndex = pd.MultiIndex.from_tuples(dfListOfTuples,names=['Fit','Parameter'])
    emptyParameterMatrix = np.zeros((df.index.size,columnIndex.size))
    parameterDataFrame = pd.DataFrame(emptyParameterMatrix,index=df.index,columns=columnIndex)
    overallrow = 0
    for observable in parameterizedObservables:
        observableDf = df.loc[observable]
        cytokineMin = np.min(observableDf.values)
        #cytokineMin = np.log10(np.min(observableDf.values))
        observablestart=0
        for row in range(df.shape[0]):
            if list(df.iloc[row,:].name)[0] == observable:
                observablestart=row
                break
        for row in range(observableDf.shape[0]):
            logy = np.array(observableDf.iloc[row,:])
            #logy = np.array(np.log10(observableDf.iloc[row,:]))
            if observable in ['IFNg','TNFa','IL-6']:
                popt,rsquared = findFits(x,logy,cytokineMin,'boundedExponential')
                parameterDataFrame.iloc[observablestart+row,parameterDataFrame.columns.get_level_values(0)=='boundedExponential'] = [rsquared,*popt]
            else:
                popt,rsquared = findFits(x,logy,cytokineMin,'logisticDoubleExponential')
                parameterDataFrame.iloc[observablestart+row,parameterDataFrame.columns.get_level_values(0)=='logisticDoubleExponential'] = [rsquared,*popt]
    with open('fitParameterPickleFile.pkl', "wb") as f:
        pickle.dump(parameterDataFrame, f)
    return parameterDataFrame

def createSumDataFrame(df,timeStart,timeEnd):
    cytokinelist = []
    sumdflist = []
    """
    for cytokine in ['IFNg','IL-2','IL-6','TNFa']:
        x = list(df.columns)
        
        integraldflist.append(integralDf)
        cytokinelist.append(cytokine)
    fullintegraldf = pd.concat(integraldflist,axis=1)
    """
    return sumdf

def createIntegralDataFrame(df,parameterDf,timeStart,timeEnd):
    cytokinelist = []
    integraldflist = []
    for cytokine in ['IFNg','IL-2','IL-6','TNFa']:
        x = list(df.columns)
        integralResults = []
        cytokineParameterDf = parameterDf.loc[cytokine]
        for row in range(cytokineParameterDf.shape[0]):
            if cytokine in ['IFNg','TNFa','IL-6']:
                equation = 'boundedExponential'
                parameters = list(cytokineParameterDf.iloc[row,:].loc[equation])[1:]
                r2 = list(parameterDf.iloc[row,:].loc[equation])[0]
                result = integrate.quad(boundedExponential, 0, np.max(x),args=tuple(parameters))
            else:
                equation = 'logisticDoubleExponential'
                parameters = list(cytokineParameterDf.iloc[row,:].loc[equation])[1:]
                r2 = list(parameterDf.iloc[row,:].loc[equation])[0]
                result = integrate.quad(logisticDoubleExponential, 0, np.max(x),args=tuple(parameters))
            integralResults.append(result[0])
        integralDf = pd.DataFrame(integralResults,index=cytokineParameterDf.index, columns=['Integral-'+cytokine])
        integraldflist.append(integralDf)
        cytokinelist.append(cytokine)
    fullintegraldf = pd.concat(integraldflist,axis=1)
    return fullintegraldf

def createSlopeDataFrame(df,timeStart,timeEnd):
    timeStartIndex = 0
    timeEndIndex = len(df.columns)
    index = len(df.columns)-1
    for timepoint in df.columns[::-1]:
        if timepoint >= timeStart:
            timeStartIndex = index
        index-=1
    index = 0
    for timepoint in df.columns[::-1]:
        if timepoint <= timeStart:
            timeEndIndex = index
        index+=1
    cytokinelist = []
    slopedflist = []
    observablestart=0
    for cytokine in ['CD45RB']:
        slopeResults = []
        cytokineParameterDf = df.loc[cytokine]
        x = df.columns[timeStartIndex:timeEndIndex+1]
        for row in range(cytokineParameterDf.shape[0]):
            timePartitionedY = cytokineParameterDf.iloc[row,timeStartIndex:timeEndIndex+1]
            result = linregress(x,timePartitionedY.values.ravel())
            slopeResults.append(result[0])
        slopeDf = pd.DataFrame(slopeResults,index=cytokineParameterDf.index, columns=['Slope-'+cytokine])
        slopedflist.append(slopeDf)
    fullslopedf = pd.concat(slopedflist,axis=1)
    return fullslopedf

def plotFits(df,cytokine,equation):
    df = df.loc[cytokine]
    plottingDf2 = df.stack().reset_index()
    plottingDf2.rename(columns={0:cytokine+' Concentration (nM)'},inplace=True)
    idx = pd.IndexSlice
    parameterDf = pickle.load(open('fitParameterPickleFile.pkl','rb')).loc[idx[cytokine],idx[equation,:]]
    curveFitPointsList = []
    numPoints = 101
    plottingMatrix = np.zeros([df.shape[0],numPoints])
    x = list(df.columns)
    curveFitPlotPoints = np.linspace(0.01,np.max(x),numPoints)
    for row in range(parameterDf.shape[0]):
        parameters = list(parameterDf.iloc[row,:])[1:]
        if equation == 'logisticDoubleExponential':
            curveFitPlotPointsY = 10**logisticDoubleExponential(curveFitPlotPoints,*parameters)
        else:
            curveFitPlotPointsY = 10**boundedExponential(curveFitPlotPoints,*parameters)
        plottingMatrix[row,:] = curveFitPlotPointsY 
    dfXVals = np.tile(curveFitPlotPoints,df.shape[0])
    df2y = pd.DataFrame(plottingMatrix,index=df.index,columns=range(1,numPoints+1))
    df2y.columns.name = 'Points'
    plottingDf = df2y.stack().reset_index()
    plottingDf['FitPointX'] = pd.Series(dfXVals,index=plottingDf.index)
    plottingDf.rename(columns={0:'FitPointY'},inplace=True)
    ax = sns.lineplot(x='FitPointX',y='FitPointY',hue='Peptide',style='Concentration',data=plottingDf,ci=False,legend=False)
    #plt.show()
    ax2 = sns.scatterplot(x='Time',y=cytokine+' Concentration (nM)',hue='Peptide',style='Concentration',ci=None,data=plottingDf2,legend=False,alpha = 0.5)
    plt.gcf().get_axes()[0].set_yscale('log')
    #ax2.fig.get_axes()[0].set_yscale('log')
    plt.show()

"""
df = pickle.load(open('../../experiments/20190412-PeptideComparison_OT1_Timeseries_18/semiProcessedData/cytokineConcentrationPickleFile-20190412-PeptideComparison_OT1_Timeseries_18-modified.pkl','rb'))

pdf = createParameterDataFrame(df)
#plotFits(df)

fullintegraldf = createIntegralDataFrame(df,pdf,0,100)

#plottingdf = xyintegraldf.unstack('Cytokine').xs('Integral of Timeseries',axis=1,level='Statistic')
#plottingdf = plottingdf.reset_index()
#sns.relplot(x=xcytokine,size = 'Concentration',y='IL-2',hue='Peptide',data=plottingdf,kind='scatter')
#plt.show()
"""

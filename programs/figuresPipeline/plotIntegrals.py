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
        lbounds = [(max(logy)-cytokineMin)/2,0,cytokineMin-epsilon]
        ubounds = [(max(logy)-cytokineMin)*2,10,cytokineMin+epsilon]
        #popt,pcov = curve_fit(boundedExponential,x,logy,sigma=logy,bounds=(lbounds,ubounds))
        #rsquared = r_squared(x,logy,boundedExponential,popt)
        try:
            popt,pcov = curve_fit(boundedExponential,x,logy,sigma=logy,bounds=(lbounds,ubounds))
            rsquared = round(r_squared(x,logy,boundedExponential,popt),3)
        except:
            return [0,0,cytokineMin],0.0
    else:
        lbounds = [(max(logy)-cytokineMin)-epsilon,0,0,0,x[logy.tolist().index(max(logy))],cytokineMin-epsilon]
        ubounds = [(max(logy)-cytokineMin)+epsilon,2,2,x[logy.tolist().index(max(logy))],max(x),cytokineMin+epsilon]
        #popt,pcov = curve_fit(logisticDoubleExponential,x,logy,sigma=logy,bounds=(lbounds,ubounds))
        #rsquared = round(r_squared(x,logy,logisticDoubleExponential,popt),3)
        try:
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
    for observable in parameterizedObservables:
        observableDf = df.loc[observable]
        cytokineMin = np.min(observableDf.values)
        #cytokineMin = np.log10(np.min(observableDf.values))
        for row in range(observableDf.shape[0]):
            logy = np.array(observableDf.iloc[row,:])
            #logy = np.array(np.log10(observableDf.iloc[row,:]))
            print(logy)
            if observable in ['IFNg','TNFa','IL-6']:
                popt,rsquared = findFits(x,logy,cytokineMin,'boundedExponential')
                print(rsquared)
                parameterDataFrame.loc[observable].iloc[row,parameterDataFrame.columns.get_level_values(0)=='boundedExponential'] = [rsquared,*popt]
            else:
                popt,rsquared = findFits(x,logy,cytokineMin,'logisticDoubleExponential')
                parameterDataFrame.loc[observable].iloc[row,parameterDataFrame.columns.get_level_values(0)=='logisticDoubleExponential'] = [rsquared,*popt]
    print(parameterDataFrame)
    with open('fitParameterPickleFile.pkl', "wb") as f:
        pickle.dump(parameterDataFrame, f)
    return parameterDataFrame

def plotFits(df,cytokine,equation):
    df = df.loc[cytokine]
    plottingDf2 = df.stack().reset_index()
    plottingDf2.rename(columns={0:cytokine+' Concentration (nM)'},inplace=True)
    idx = pd.IndexSlice
    parameterDf = pickle.load(open('fitParameterPickleFile.pkl','rb')).loc[idx[cytokine],idx[equation,:]]
    print(parameterDf)
    curveFitPointsList = []
    numPoints = 101
    plottingMatrix = np.zeros([df.shape[0],numPoints])
    x = list(df.columns)
    curveFitPlotPoints = np.linspace(0.01,np.max(x),numPoints)
    for row in range(parameterDf.shape[0]):
        parameters = list(parameterDf.iloc[row,:])[1:]
        print(parameters)
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
    ax = sns.relplot(x='FitPointX',y='FitPointY',hue='Peptide',style='Concentration',data=plottingDf,kind='line',ci=False)
    ax.fig.get_axes()[0].set_yscale('log')
    plt.show()
    #ax2 = sns.relplot(x='Time',y='IL-2 Concentration (nM)',hue='Peptide',style='Concentration',ci=None,data=plottingDf2)
    #ax2.fig.get_axes()[0].set_yscale('log')
    #plt.show()

def computeAndPlotIntegrals(df,cytokine,equation):
    df = df.loc['IL-2']
    idx = pd.IndexSlice
    parameterDf = pickle.load(open('fitParameterPickleFile.pkl','rb')).loc[idx[cytokine],idx[equation,:]]
    x = list(df.columns)
    integralResults = []
    for row in range(parameterDf.shape[0]):
        parameters = list(parameterDf.iloc[row,:])[1:]
        if equation == 'logisticDoubleExponential':
            result = integrate.quad(logisticDoubleExponential, 0, np.max(list(df.columns)),args=tuple(parameters))
        else:
            result = integrate.quad(boundedExponential, 0, np.max(list(df.columns)),args=tuple(parameters))
        integralResults.append(result[0])
    integralDf = pd.DataFrame(np.array(integralResults),index=df.index,columns=['Integral of Timeseries'])
    
    integralPlottingDf = integralDf.reset_index()
    currentLevelValues = list(integralPlottingDf['Concentration'])
    sortedOldLevelValues,newLevelValues = sortSINumerically(currentLevelValues,False,True)
    #Need to interpret parenthetical units for x to get 1e9
    #units = '1'+s[s.find("(")+1:s.find(")")] 
    units = '1nM'
    scaledSortedUnits,sortedUnits = sortSINumerically([units],False,True)
    scaledNewLevelValues = [float(i) / float(sortedUnits[0]) for i in newLevelValues]
    #print(integralPlottingDf)
    #integralPlottingDf['Concentration'] = scaledNewLevelValues
    #print(integralPlottingDf)
    #sys.exit(0)
    ax = sns.catplot(x = 'Concentration',y='Integral of Timeseries',hue='Peptide',kind='point',ci=None,data=integralPlottingDf,join=False)
    plt.show()

df = pickle.load(open('/Volumes/Group05/Altan-Bonnet/Sooraj/experiments/20190412-PeptideComparison_OT1_Timeseries_18/semiProcessedData/cytokineConcentrationPickleFile-20190412-PeptideComparison_OT1_Timeseries_18-modified.pkl','rb'))
pdf = createParameterDataFrame(df)
print(pdf)
sys.exit(0)
for cytokine in ['IFNg','IL-2','IL-6','TNFa']:
    if cytokine in ['IFNg','TNFa','IL-6']:
        equation = 'boundedExponential'
    else:
        equation = 'logisticDoubleExponential'
    plotFits(df,cytokine,equation)
    computeAndPlotIntegrals(df,cytokine,equation)
#plotFits(df)

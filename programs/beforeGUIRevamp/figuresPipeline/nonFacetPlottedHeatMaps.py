#!/usr/bin/env python3  
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import pickle
import os
import sys
import matplotlib
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
from itertools import groupby

def add_line(ax, xpos, ypos):
    line = plt.Line2D([xpos, xpos+.1], [ypos, ypos], transform=ax.transAxes, color='black')
    line.set_clip_on(False)
    ax.add_line(line)

def add_vline(ax, xpos, ypos,length):
    line = plt.Line2D([xpos, xpos], [ypos, ypos+length],transform=ax.transAxes, color='black')
    line.set_clip_on(False)
    ax.add_line(line)

def add_hline(ax, xpos, ypos,length):
    line = plt.Line2D([xpos, xpos+length], [ypos, ypos],transform=ax.transAxes, color='black')
    line.set_clip_on(False)
    ax.add_line(line)

def label_len(my_index,level):
    labels = my_index.get_level_values(level)
    return [(k, sum(1 for i in g)) for k,g in groupby(labels)]

def label_conditions(ax, df):
    xpos = -.1
    scale = 1./df.index.size
    for level in range(len(list(df.index.names)))[::-1]:
        pos = len(df.index)
        for label, rpos in label_len(df.index,level):
            lypos = (pos - .6 * rpos)*scale
            if (label not in list(list(df.index.names)[-1])):
                t = ax.text(xpos+0.05, lypos, label, ha='center', transform=ax.transAxes)
            else:
                t = ax.text(-0.05, lypos, label, ha='center', transform=ax.transAxes)
            add_line(ax, xpos, pos*scale)
            pos -= rpos
        add_line(ax, xpos, pos*scale)
        xpos -= .1

def label_time(ax,df):
    ypos = -1./df.index.size
    scale = 1./len(df.columns.values)
    xpos=0
    for timepoint in df.columns.values:
        add_vline(ax,xpos,ypos,ypos*-1)
        ax.text(xpos+scale/2, ypos/2, str(timepoint), ha='center', va='center',transform=ax.transAxes)
        xpos+=scale
    ax.text(-0.05*len(df.index.names), -0.5*(1/df.index.size), df.columns.name+' elapsed (hours):', ha='center', va='center',transform=ax.transAxes)
    add_hline(ax,-0.1*len(df.index.names),ypos,1+(0.1*len(df.index.names)))
    add_vline(ax,-0.1*len(df.index.names),ypos,ypos*-1)
    #add_vline(ax,-0.1*len(df.index.names),ypos,1+scale)

def label_headers(ax,df):
    scale = len(df.index.names)
    xpos = -0.1*scale
    lineheight = 1+(1./df.index.size)
    for name in df.index.names:
        add_vline(ax,xpos,0,lineheight)
        ax.text(xpos+0.05,1+(lineheight-1)/2,name,ha='center',va='center',transform=ax.transAxes)
        xpos+=0.1
    add_hline(ax,df.index.nlevels*-0.1,lineheight,1+(df.index.nlevels*0.1))
    add_hline(ax,0,0,1)
    add_hline(ax,0,1,1)
    add_vline(ax,0,0,lineheight)
    add_vline(ax,1,0-1./df.index.size,1+2/df.index.size)
    return 1+(lineheight-1)/2

def getLongestLabel(df,cyt):
    if(cyt):
        temp = df.loc[pd.unique(df.index.get_level_values(0)[0])]
    else:
        temp = df
    maxLabelLength = 0
    allLabelList = list(temp.index.names)
    for level in range(len(list(temp.index.names))):
        for name in temp.index.get_level_values(level):
            allLabelList.append(name)
    for label in allLabelList:
        if len(label) > maxLabelLength:
            maxLabelLength = len(label)
            longestLabel = label
    return longestLabel

def getLongestLabelDimensions(ax,df,longestLabel,r):
    if longestLabel in df.index.names:
        for name in df.index.names:
            t = ax.text(0,0,name,ha='center',va='center',transform=ax.transAxes,fontweight='bold')
            if name == longestLabel:
                bb = t.get_window_extent(renderer=r)
                twidth = bb.width
                theight = bb.height
                break
    else:
        foundLongest = False
        for level in range(len(df.index.names))[::-1]:
            for label, rpos in label_len(df.index,level):
                t = ax.text(0, 0, label, ha='center', transform=ax.transAxes)
                if label == longestLabel:
                    bb = t.get_window_extent(renderer=r)
                    twidth = bb.width
                    theight = bb.height
                    foundLongest = True
                    break
            if(foundLongest):
                break
    return twidth,theight

def returnLabelScalingSize(df):
    observableList = list(df.index.get_level_values(0))
    h = 20
    fig1 = plt.figure(num=1,figsize=(h*1.5,h),dpi=150,facecolor='w',edgecolor='k')
    ax1=fig1.add_subplot(111)
    temp = df.loc[observableList[0]]
    cax1 = sns.heatmap(temp,vmin=np.max([-13+np.log10(1e9),np.nanmin(temp.values)]),vmax=np.log10(1e-9*1e9),yticklabels=False)
    fig1.subplots_adjust(left=.1*len(list(temp.index.names)))
    r = fig1.canvas.get_renderer()
    bbox = ax1.get_window_extent().transformed(fig1.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    scale = 1./temp.index.size
    bwidth = width*fig1.dpi*0.1
    bheight = height*fig1.dpi*scale
    
    size=10
    pad=0.5
    decrement = 0.1
    longestLabel = getLongestLabel(df,True)
    for i in range(0,int(size/decrement)):
        sns.set(font_scale=size)
        twidth,theight = getLongestLabelDimensions(ax1,df,longestLabel,r)
        sns.set(font_scale=10)
        if(twidth < bwidth and theight < bheight):
            size-=decrement*pad
            break
        else:
            size-=decrement 
    plt.clf()
    return size

#Heatmap of Experiment
def createCombinedHeatMap(folderName,secondPath,df,concUnit,concUnitPrefix,useModifiedDf,dataType):
    if dataType == 'cell':
        for statistic in list(pd.unique(df.index.get_level_values('Statistic'))): 
            statisticdf = df.xs(statistic,level='Statistic')
            fig1 = plt.figure(num=1,figsize=(10*len(list(pd.unique(statisticdf.index.get_level_values('Marker')))),10*len(list(pd.unique(statisticdf.index.get_level_values('CellType'))))),dpi=120,facecolor='w',edgecolor='k')
            sns.set(font_scale=1.0)
            indexNumber = 1
            for cellType in list(pd.unique(statisticdf.index.get_level_values('CellType'))):
                celldf = statisticdf.xs(cellType,level='CellType')
                for marker in list(pd.unique(celldf.index.get_level_values('Marker'))):
                    ax1=fig1.add_subplot(len(list(pd.unique(statisticdf.index.get_level_values('CellType')))),len(list(pd.unique(statisticdf.index.get_level_values('Marker')))),indexNumber)
                    currentdf = celldf.loc[marker]
                    if('%' in statistic):
                        cax1 = sns.heatmap(currentdf.iloc[:,:],vmin=0,vmax=100)
                    elif('Mean' in statistic):
                        cax1 = sns.heatmap(currentdf.iloc[:,:],vmin=0)
                    else:
                        cax1 = sns.heatmap(currentdf.iloc[:,:],vmin=0)
                    plt.title(cellType+'-'+marker)
                    cax1.collections[0].colorbar.set_label(statistic,rotation=270,labelpad=20)
                    indexNumber+=1
            plt.subplots_adjust(hspace = 0.2)
            statistic = statistic.replace(' ','')
            if(useModifiedDf):
                fig1.savefig('fullyProcessedFigures/heatMap-'+folderName+'-cell-'+statistic+'-modified.png',bbox_inches='tight')
            else:
                fig1.savefig('fullyProcessedFigures/heatMap-'+folderName+'-cell-'+statistic+'.png',bbox_inches='tight')
            plt.close(fig1)
            print(statistic+' plot saved')
    else:
        fig1 = plt.figure(num=1,figsize=(55,35),dpi=120,facecolor='w',edgecolor='k')
        sns.set(font_scale=1.0)
        observableList = list(pd.unique(df.index.get_level_values(0)))
        print(observableList)
        for observable,indexNumber in zip(observableList,range(len(observableList))):
            if len(observableList) <= 9:
                ax1=fig1.add_subplot(3,3,indexNumber+1)
            else:
                ax1=fig1.add_subplot(4,5,indexNumber+1)
            currentdf = df.loc[observable]
             
            plt.title(observable)
            
            if dataType == 'cyt': 
                #Pick greater of: lower LOD, minimum value in experiment for heatmap min, dont pick smaller of higher LOD, max value in experiment for heatmapm max because signal is never that high
                LODParameters = pickle.load(open('semiProcessedData/LODParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "rb"))
                lowerGFILOD = LODParameters[observable][0]
                upperGFILOD = LODParameters[observable][1]
                lowerConcLOD = LODParameters[observable][2]
                upperConcLOD = LODParameters[observable][3]
                
                cax1 = sns.heatmap(np.log10(currentdf),vmin=np.log10(lowerConcLOD),vmax=np.log10(upperConcLOD))
                cax1.collections[0].colorbar.set_label('Concentration (log$_{10}$(['+observable+' '+concUnitPrefix+']))',rotation=270,labelpad=20)
            else:
                cax1 = sns.heatmap(currentdf)
                cax1.collections[0].colorbar.set_label(observable,rotation=270,labelpad=20)
                plt.yticks(rotation='horizontal')
                observable = observable.replace(' ','')

        plt.subplots_adjust(hspace = 0.2)
        if(useModifiedDf):
            fig1.savefig('fullyProcessedFigures/heatMap-'+folderName+'-'+dataType+'-all-modified.png',bbox_inches='tight')
        else:
            fig1.savefig('fullyProcessedFigures/heatMap-'+folderName+'-'+dataType+'-all.png',bbox_inches='tight')
        plt.close(fig1)
        print('All observable plot saved')

def createIndividualHeatMaps(folderName,secondPath,df,concUnit,concUnitPrefix,useModifiedDf,dataType):
    if dataType == 'cell':
        for cellType in list(pd.unique(df.index.get_level_values('CellType'))):
            celldf = df.xs(cellType,level='CellType')
            for marker in list(pd.unique(celldf.index.get_level_values('Marker'))):
                markerdf = celldf.xs(marker,level='Marker')
                for statistic in list(pd.unique(markerdf.index.get_level_values('Statistic'))):
                    statisticdf = markerdf.xs(statistic,level='Statistic')
                    titlestring = cellType+'-'+marker
                    fileNameString = cellType+'-'+marker+'-'+statistic
                    fileNameString = fileNameString.replace(' ','')
                    h = 20
                    fig1 = plt.figure(num=1,figsize=(h*1.5,h),dpi=150,facecolor='w',edgecolor='k')
                    sns.set(font_scale=returnLabelScalingSize(statisticdf))
                    ax1=fig1.add_subplot(111)
                    if('%' in statistic):
                        cax1 = sns.heatmap(statisticdf,vmax=100)
                    elif('Positive GFI' == statistic):
                        cax1 = sns.heatmap(statisticdf)
                    else:
                        cax1 = sns.heatmap(statisticdf,vmin=0)
                    plt.axis('off')

                    label_conditions(ax1,statisticdf)
                    label_time(ax1,statisticdf)
                    titleheight = label_headers(ax1,statisticdf)
                    
                    fig1.subplots_adjust(left=.1*statisticdf.index.nlevels)
                    ax1.text(0.5,titleheight,titlestring,fontweight='bold',ha='center',va='center',transform=ax1.transAxes)
                    cax1.collections[0].colorbar.set_label(statistic,rotation=270,labelpad=20)
                    if(useModifiedDf):
                        fig1.savefig('fullyProcessedFigures/heatMap-'+folderName+'-'+dataType+'-'+fileNameString+'-modified.png',bbox_inches='tight')
                    else:
                        fig1.savefig('fullyProcessedFigures/heatMap-'+folderName+'-'+dataType+'-'+fileNameString+'.png',bbox_inches='tight') 
                    plt.clf()
                    print(fileNameString+' plot saved')
    else:
        observableList = list(pd.unique(df.index.get_level_values(0)))
        for observable,indexNumber in zip(observableList,range(len(observableList))):
            currentdf = df.loc[observable]
            
            h = 20
            sns.set(font_scale=returnLabelScalingSize(currentdf))
            fig1 = plt.figure(num=1,figsize=(h*1.5,h),dpi=600,facecolor='w',edgecolor='k')
            ax1=fig1.add_subplot(111)

            if dataType == 'cyt': 
                #Pick greater of: lower LOD, minimum value in experiment for heatmap min, dont pick smaller of higher LOD, max value in experiment for heatmapm max because signal is never that high
                LODParameters = pickle.load(open('semiProcessedData/LODParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "rb"))
                lowerGFILOD = LODParameters[observable][0]
                upperGFILOD = LODParameters[observable][1]
                lowerConcLOD = LODParameters[observable][2]
                upperConcLOD = LODParameters[observable][3]
            
                cax1 = sns.heatmap(np.log10(currentdf),vmin=np.log10(lowerConcLOD),vmax=np.log10(upperConcLOD),xticklabels=False,yticklabels=False)
                cax1.collections[0].colorbar.set_label('Concentration (log$_{10}$(['+observable+' '+concUnitPrefix+']))',rotation=270,labelpad=20)
            else:
                cax1 = sns.heatmap(currentdf,xticklabels=False,yticklabels=False)
                cax1.collections[0].colorbar.set_label(observable,rotation=270,labelpad=20)
            ax1.set_xlabel('')
            ax1.set_ylabel('')
            
            label_conditions(ax1,currentdf)
            label_time(ax1,currentdf)
            titleheight = label_headers(ax1,currentdf)
            
            fig1.subplots_adjust(left=.1*currentdf.index.nlevels)
            ax1.text(0.5,titleheight,observableList[indexNumber],fontweight='bold',ha='center',va='center',transform=ax1.transAxes)
            
            if(useModifiedDf):
                fig1.savefig('fullyProcessedFigures/heatMap-'+folderName+'-'+dataType+'-'+str(observable).replace(' ','')+'-modified.png',bbox_inches='tight') 
            else:
                fig1.savefig('fullyProcessedFigures/heatMap-'+folderName+'-'+dataType+'-'+str(observable).replace(' ','')+'.png',bbox_inches='tight') 
            plt.clf()
            print(str(observable)+' plot saved')

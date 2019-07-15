#!/usr/bin/env python3
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import pickle,sys,os
from itertools import groupby
from matplotlib.widgets import RadioButtons,Button,CheckButtons,TextBox
print(os.getcwd())
sys.path.insert(0, '../../programs/dataProcessing/')
print(os.getcwd())
from miscFunctions import reindexDataFrame
from matplotlib import colors,ticker

#Button width conserved across gui figures
buttonWidth = 0.1/2
buttonLength = 0.075/2
buttonXStart = 0.5-(0.01+buttonWidth)
buttonYStart = 0.01

dividerLength = 0.16

def add_vline(ax, xpos, ypos,length):
    line = plt.Line2D([xpos, xpos], [ypos, ypos+length],transform=ax.transAxes, color='black')
    line.set_clip_on(False)
    ax.add_line(line)

def add_hline(ax, xpos, ypos,length):
    line = plt.Line2D([xpos, xpos+length], [ypos, ypos],transform=ax.transAxes, color='black')
    line.set_clip_on(False)
    ax.add_line(line)

def draw_borders(ax,df):
    add_hline(ax,0,0,1)
    add_hline(ax,0,1,1)
    add_vline(ax,1,0,1)
    add_vline(ax,0,0,1)

def label_len(my_index,level):
    labels = my_index.get_level_values(level)
    return [(k, sum(1 for i in g)) for k,g in groupby(labels)]

def label_index(ax, df):
    xpos = -1*dividerLength
    scale = 1./df.index.size
    for level in range(len(list(df.index.names)))[::-1]:
        pos = len(df.index)
        for label, rpos in label_len(df.index,level):
            lypos = (pos - .6 * rpos)*scale
            if (label not in list(list(df.index.names)[-1])):
                t = ax.text(xpos+(dividerLength/2), lypos, label, va='center',ha='center', transform=ax.transAxes,rotation=0)
            else:
                t = ax.text(-1*(dividerLength/2), lypos, label, va='center',ha='center', transform=ax.transAxes,rotation=0)
            add_hline(ax, xpos, pos*scale,dividerLength)
            pos -= rpos
        add_hline(ax, xpos, pos*scale,dividerLength)
        xpos -= dividerLength

def label_columns(ax,df):
    ypos = -1./df.index.size
    scale = 1./len(df.columns.values)
    xpos=0
    for timepoint in df.columns.values:
        add_vline(ax,xpos,ypos*0.8,ypos*-1*0.8)
        #if float(timepoint) % 1 == 0:
        #timepoint = int(timepoint)
        ax.text(xpos+scale/2, ypos/2, str(timepoint), ha='center', va='center',transform=ax.transAxes)
        xpos+=scale
    add_vline(ax,xpos,ypos*0.8,ypos*-1*0.8)
    ax.text(0.5, ypos*1.4, df.columns.name, ha='center', va='center',transform=ax.transAxes)

def label_headers(ax,df):
    scale = len(df.index.names)
    xpos = -1*dividerLength*scale
    lineheight = 1+(1./df.index.size)
    for name in df.index.names:
        ax.text(xpos+(dividerLength/2),1+(lineheight-1)/2,name,ha='center',va='center',transform=ax.transAxes,size='x-small')
        xpos+=dividerLength

def returnHeatmapAspectRatios(data,kwargs):
    #16x16 heatmap
    basedim = 16
    scale=0.4
    hbase = 6
    abase = 1.25
    if 'row' in kwargs:
        data = data.xs(kwargs['row_order'][0],level=kwargs['row'])
    if 'col' in kwargs:
        data = data.xs(kwargs['col_order'][0],level=kwargs['col'])
    heatmapdf = data.pivot_table(index=kwargs['y'],columns=kwargs['x'], values=kwargs['z'])
    hstart = max(hbase+scale*hbase*(heatmapdf.index.size-basedim)/basedim,hbase)
    astart = max(abase+scale*abase*(heatmapdf.columns.size-basedim)/basedim,abase)
    if 'row_order' in kwargs:
        h = hstart*len(kwargs['row_order'])*scale
    else:
        h = hstart
    if 'col_order' in kwargs:
        a = astart*len(kwargs['col_order'])*scale
    else:
        a = astart
    return a,h

def draw_faceted_heatmap(data,indexingdf,xaxis,yaxis,zaxis,lognorm,cbarticks,logarithmic,symlog,symlognorm,linthresh,**kwargs):
    originalColumnOrder = list(pd.unique(indexingdf.index.get_level_values(xaxis)))
    unsortedPivotedData = data.pivot_table(index=yaxis,columns=xaxis, values=zaxis)
    indexdf = indexingdf.groupby(level=yaxis,sort=False).first()
    data = reindexDataFrame(unsortedPivotedData,indexdf,False)
    data = data[originalColumnOrder]
    data.columns.name = xaxis
    plt.axis('off')
    if logarithmic:
        g = sns.heatmap(data, norm=lognorm,**kwargs,cbar=True,cbar_kws={"ticks": cbarticks,'label':zaxis})
    elif symlog:
        linthresh = int(linthresh) 
        maxlog=int(np.ceil(np.log10(data.values.max())))
        minlog=int(np.ceil(np.log10(-1*data.values.min())))
        tick_locations=([-(10**x) for x in range(-linthresh, minlog+1, 1)][::-1]+[0.0]+[(10**x) for x in range(-linthresh,maxlog+1, 1)] )
        #generate logarithmic ticks
        g = sns.heatmap(data, norm=symlognorm,cbar_kws={'label':zaxis,'ticks':tick_locations,'format':ticker.LogFormatterMathtext()},**kwargs)
    else:
        g = sns.heatmap(data, **kwargs,cbar=True,cbar_kws={'label':zaxis})

    #Add hiearchical level names and borders to heatmap
    ax1 = plt.gca()
    label_index(ax1,data)
    draw_borders(g,data)
    label_columns(ax1,data)
     
    label_headers(ax1,data)

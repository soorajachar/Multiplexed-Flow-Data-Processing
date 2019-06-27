#!/usr/bin/env python3
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import pickle,sys,os
from itertools import groupby
from matplotlib.widgets import RadioButtons,Button,CheckButtons,TextBox
sys.path.insert(0, '../../../programs/dataProcessing/')
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
    #heatmapdf = data.pivot_table(index=kwargs['y'],columns=kwargs['x'], values=kwargs['z'])
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
        #tick_locations=([-(10**x) for x in range(-1*int(linthresh), minlog+1, 1)][::-1]+[0.0]+[(10**x) for x in range(-minlog,maxlog+1, 1)])
        tick_locations=([-(10**x) for x in range(-linthresh, minlog+1, 1)][::-1]+[0.0]+[(10**x) for x in range(-linthresh,maxlog+1, 1)] )
        #generate logarithmic ticks
        g = sns.heatmap(data, norm=symlognorm,cbar_kws={'label':zaxis,'ticks':tick_locations,'format':ticker.LogFormatterMathtext()},**kwargs)
        #g = sns.heatmap(data, norm=symlognorm,**kwargs,cbar=True,cbar_kws={'label':zaxis})
    else:
        g = sns.heatmap(data, **kwargs,cbar=True,cbar_kws={'label':zaxis})

    #Add hiearchical level names and borders to heatmap
    ax1 = plt.gca()
    label_index(ax1,data)
    draw_borders(g,data)
    label_columns(ax1,data)
     
    label_headers(ax1,data)

def heatmapSpecific_GUIWindow(labelDict,plotType,dataType):
    #Ask scale to use for colorbar
    axes = ['X Axis','Y Axis','Colorbar Axis']
    axisScalingOptions = ['Linear','Logarithmic','Biexponential']

    fig = plt.figure(figsize=(12,2.5*len(axisScalingOptions)))
    plt.axis('off')
    radiobuttons = []
    axis_title_text_boxes = {}
    lin_thresh_text_boxes = {}
    checkbuttons = []
    checkbuttons2 = []
    for axis,i in zip(axes,range(len(axes))):
        #Add axis scaling radio buttons
        rectLength = 0.15*len(axisScalingOptions)
        rectWidth = (1- (0.02+0.01*len(axes)))/len(axes)
        
        rax3 = plt.axes([0.01*(i+1)+rectWidth*i, 0.94-(buttonWidth+rectLength), rectWidth,rectLength])
        plt.text(0.5, 1.01,axis,ha='center',transform=plt.gca().transAxes)
        rax3.spines['bottom'].set_visible(False)
        rax3.spines['left'].set_visible(False)
        rax3.spines['right'].set_visible(False)
        plt.tick_params(which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False)
        if axis == 'Colorbar Axis':
            radiobuttons.append(RadioButtons(rax3,axisScalingOptions,activecolor='black'))
             
            #Add in linear threshold (for biexponential scaling) textboxes
            print('wat')
            linthreshbox = plt.axes([0.05+0.33*i, 0.35, 0.25, 0.075])
            text_box2 = TextBox(linthreshbox, 'Linear  \nThreshold: ', initial=' ')
            lin_thresh_text_boxes[axis] = text_box2
        
        axbox = plt.axes([0.05+0.33*i, 0.25, 0.25, 0.075])
        #Determine initial axis title for textboxes
        if axis == 'Colorbar Axis':
            if dataType == 'cyt':
                initial_nameC = 'Concentration (nM)'
            else:
                initial_nameC = ' '
            #Add in axis title text boxes
            text_box = TextBox(axbox, 'Title: ', initial=initial_nameC)
            axis_title_text_boxes[axis] = text_box
        elif axis == 'X Axis':
            initial_nameX = ' '
            #Add in axis title text boxes
            text_box = TextBox(axbox, 'Title: ', initial=initial_nameX)
            axis_title_text_boxes[axis] = text_box
        else:
            initial_nameY = ' '
            #Add in axis title text boxes
            text_box = TextBox(axbox, 'Title: ', initial=initial_nameY)
            axis_title_text_boxes[axis] = text_box
        
        
        #Numerical axis sorting check boxes
        if axis is not 'Colorbar Axis':
            rax2 = plt.axes([0.05+i*0.33,0.09,0.15,0.15])
            rax2.spines['bottom'].set_visible(False)
            rax2.spines['left'].set_visible(False)
            rax2.spines['right'].set_visible(False)
            rax2.spines['top'].set_visible(False)
            checkbuttons.append(CheckButtons(rax2,['Sort '+axis[0]+' numerically'],actives=[False]))
        else:
            rax2 = plt.axes([0.05+i*0.33,0.09,0.15,0.15])
            rax2.spines['bottom'].set_visible(False)
            rax2.spines['left'].set_visible(False)
            rax2.spines['right'].set_visible(False)
            rax2.spines['top'].set_visible(False)
            checkbuttons2.append(CheckButtons(rax2,['Include Level Names'],actives=[True]))

    linThresholdValues = {}
    axisTitleValues = {'X Axis':initial_nameX,'Y Axis':initial_nameY,'Colorbar Axis':initial_nameC}
    def submitAxisTitleX(text):
        axisTitleValues['X Axis'] = text
    def submitAxisTitleY(text):
        axisTitleValues['Y Axis'] = text
    def submitAxisTitleC(text):
        axisTitleValues['Colorbar Axis'] = text
    axis_title_text_boxes['X Axis'].on_submit(submitAxisTitleX)
    axis_title_text_boxes['Y Axis'].on_submit(submitAxisTitleY)
    axis_title_text_boxes['Colorbar Axis'].on_submit(submitAxisTitleC)
    
    def submitLinThresholdC(text):
        linThresholdValues['Colorbar Axis'] = float(text)
    lin_thresh_text_boxes['Colorbar Axis'].on_submit(submitLinThresholdC)

    class GUIButtons4(object):
        def OKradiotext4(self, event):
            plotOptions = {}
            radioValues = {}
            radioValues['Colorbar Axis'] = radiobuttons[0].value_selected
            for i in range(len(checkbuttons)):
                numericBoolean = checkbuttons[i].get_status()
                plotOptions['numeric'+axes[i][0]] = numericBoolean
            includeLevelNamesBoolean = checkbuttons2[0].get_status()
            plotOptions['IncludeLevelNames'] = includeLevelNamesBoolean
            plt.close()

            plotOptions['axisScaling'] = radioValues
            plotOptions['linThreshold'] = linThresholdValues
            plotOptions['axisTitles'] = axisTitleValues
            print(plotOptions)
            with open('semiProcessedData/gui-plotOptions.pkl','wb') as f:
                pickle.dump(plotOptions,f)
        
        def Quit(self, event):
            sys.exit(0)    

    callback = GUIButtons4()
    axOK = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
    axQuit = plt.axes([buttonXStart+buttonWidth+0.01,buttonYStart, buttonWidth, buttonLength])
    bOK = Button(axOK, 'OK')
    bOK.on_clicked(callback.OKradiotext4)
    bQuit = Button(axQuit, 'Quit')
    bQuit.on_clicked(callback.Quit)
    plt.show()

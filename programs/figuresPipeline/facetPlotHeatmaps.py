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
                t = ax.text(xpos+(dividerLength/2), lypos, label, ha='center', transform=ax.transAxes,rotation=0)
            else:
                t = ax.text(-1*(dividerLength/2), lypos, label, ha='center', transform=ax.transAxes,rotation=0)
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
        if timepoint%1 ==0:
            timepoint = int(timepoint)
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
            linthreshbox = plt.axes([0.1+0.5*i, 0.35, 0.25, 0.075])
            text_box2 = TextBox(linthreshbox, 'Linear  \nThreshold: ', initial=' ')
            lin_thresh_text_boxes[axis] = text_box2
        
        #Determine initial axis title for textboxes
        if axis == 'Colorbar Axis':
            if dataType == 'cyt':
                initial_name = 'Concentration (nM)'
            else:
                initial_name = ' '
        elif axis == 'X Axis':
            initial_name = 'Time (hours)'
        else:
            initial_name = ' '
        
        #Add in axis title text boxes
        axbox = plt.axes([0.05+0.33*i, 0.25, 0.25, 0.075])
        text_box = TextBox(axbox, 'Title: ', initial=initial_name)
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
    axisTitleValues = {'X Axis':initial_name,'Y Axis':initial_name,'Colorbar Axis':initial_name}
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
                plotOptions[axes[i]] = numericBoolean
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

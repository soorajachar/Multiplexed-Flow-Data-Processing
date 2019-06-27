#!/usr/bin/env python3
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import markers
import numpy as np
import seaborn as sns
import pandas as pd
import pickle,os,math,sys,itertools,re
from matplotlib.widgets import RadioButtons,Button,CheckButtons,TextBox

#Button width conserved across gui figures
buttonWidth = 0.1/2
buttonLength = 0.075/2
buttonXStart = 0.5-(0.01+buttonWidth)
buttonYStart = 0.01

def scaleXandY2D(ax,plotOptions):
    #X and Y Axis Scaling for 2D plots
    for axis in plotOptions['axisScaling']:
        if 'Y' in axis:
            if plotOptions['axisScaling'][axis] == 'Logarithmic':
                ax.fig.get_axes()[0].set_yscale('log')
            elif plotOptions['axisScaling'][axis] == 'Biexponential':
                ax.fig.get_axes()[0].set_yscale('symlog',linthreshx=plotOptions['linThreshold'][axis])
        else:
            if plotOptions['axisScaling'][axis] == 'Logarithmic':
                ax.fig.get_axes()[0].set_xscale('log')
            elif plotOptions['axisScaling'][axis] == 'Biexponential':
                ax.fig.get_axes()[0].set_xscale('symlog',linthreshx=plotOptions['linThreshold'][axis])

def scatterLineBarSpecific_GUIWindow(labelDict,plotType,dataType):
    #Ask scale to use for axes (both y and x if relplot; only y if categorical)
    if plotType == 'categorical' or plotType == '1d':
        axes = ['Y Axis']
    else:
        axes = ['X Axis','Y Axis']
    constantAxes = ['X Axis','Y Axis']

    axisScalingOptions = ['Linear','Logarithmic','Biexponential']

    fig = plt.figure(figsize=(8,2.5*len(axisScalingOptions)))
    plt.axis('off')
    radiobuttons = []
    axis_title_text_boxes = {}
    lin_thresh_text_boxes = {}
    for axis,i in zip(constantAxes,range(len(constantAxes))):
        #Add axis scaling radio buttons
        rectLength = 0.15*len(axisScalingOptions)
        rectWidth = (1- (0.02+0.01*len(constantAxes)))/len(constantAxes)
        
        rax3 = plt.axes([0.01*(i+1)+rectWidth*i, 0.94-(buttonWidth+rectLength), rectWidth,rectLength])
        plt.text(0.5, 1.01,axis,ha='center',transform=plt.gca().transAxes)
        rax3.spines['bottom'].set_visible(False)
        rax3.spines['left'].set_visible(False)
        rax3.spines['right'].set_visible(False)
        plt.tick_params(which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False)
        if(axis == 'X Axis'):
            if(len(constantAxes) == len(axes)):
                radiobuttons.append(RadioButtons(rax3,axisScalingOptions,activecolor='black'))
        else:
            radiobuttons.append(RadioButtons(rax3,axisScalingOptions,activecolor='black'))

        #Determine initial axis title for textboxes
        if 'Y' in axis:
            if dataType == 'cyt':
                initial_name = 'Concentration (nM)'
            else:
                if plotType == '1d':
                    initial_name = 'Count'
                else:
                    initial_name = ' '
        else:
            if plotType == 'ordered':
                radioValues = pickle.load(open('semiProcessedData/gui-radioVals.pkl','rb'))
                initial_name = radioValues['X Axis Values']
            elif plotType == '1d':
                #initial_name = radioValues['Y Axis Values']
                initial_name = 'GFI'
            else:
                radioValues = pickle.load(open('semiProcessedData/gui-radioVals.pkl','rb'))
                initial_name = radioValues['Order']
        
        #Add in axis title text boxes
        axbox = plt.axes([0.1+0.5*i, 0.25, 0.35, 0.075])
        text_box = TextBox(axbox, 'Title: ', initial=initial_name)
        axis_title_text_boxes[axis] = text_box
        
        #Add in linear threshold (for biexponential scaling) textboxes
        if(axis == 'X Axis'):
            if(len(constantAxes) == len(axes)):
                linthreshbox = plt.axes([0.1+0.5*i, 0.35, 0.35, 0.075])
                text_box2 = TextBox(linthreshbox, 'Linear  \nThreshold: ', initial=' ')
                lin_thresh_text_boxes[axis] = text_box2
        else:
            linthreshbox = plt.axes([0.1+0.5*i, 0.35, 0.35, 0.075])
            text_box2 = TextBox(linthreshbox, 'Linear  \nThreshold: ', initial=' ')
            lin_thresh_text_boxes[axis] = text_box2
   
    if plotType == 'ordered':
        checkbuttons = []
        rax2 = plt.axes([0.05,0.09,0.15,0.15])
        rax2.spines['bottom'].set_visible(False)
        rax2.spines['left'].set_visible(False)
        rax2.spines['right'].set_visible(False)
        rax2.spines['top'].set_visible(False)
        checkbuttons.append(CheckButtons(rax2,['Sort X numerically'],actives=[False]))
    
    linThresholdValues = {}
    axisTitleValues = {'X Axis':initial_name,'Y Axis':' '}
    def submitAxisTitleX(text):
        axisTitleValues['X Axis'] = text
    def submitAxisTitleY(text):
        axisTitleValues['Y Axis'] = text
    axis_title_text_boxes['X Axis'].on_submit(submitAxisTitleX)
    axis_title_text_boxes['Y Axis'].on_submit(submitAxisTitleY)
    
    def submitLinThresholdX(text):
        linThresholdValues['X Axis'] = float(text)
    def submitLinThresholdY(text):
        linThresholdValues['Y Axis'] = float(text)
    if(len(constantAxes) == len(axes)):
        lin_thresh_text_boxes['X Axis'].on_submit(submitLinThresholdX)
        lin_thresh_text_boxes['Y Axis'].on_submit(submitLinThresholdY)
    else:
        lin_thresh_text_boxes['Y Axis'].on_submit(submitLinThresholdY)

    class GUIButtons4(object):
        def OKradiotext4(self, event):
            plotOptions = {}
            radioValues = {}
            for radiobutton,axis in zip(radiobuttons,axes):
                radioValues[axis] = radiobutton.value_selected
            if plotType == 'ordered':
                numericXBoolean = checkbuttons[0].get_status()
                plotOptions['numericX'] = numericXBoolean
            else:
                plotOptions['numericX'] = [False]
            plt.close()
            plotOptions['axisScaling'] = radioValues
            plotOptions['linThreshold'] = linThresholdValues
            plotOptions['axisTitles'] = axisTitleValues
            print(plotOptions)
            print(os.getcwd())
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

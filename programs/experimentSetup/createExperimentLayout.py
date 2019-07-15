#!/usr/bin/env python3
import numpy as np
import pickle,os,string,json
from tkinter import *
import tkinter as tk
from functools import reduce
from PIL import ImageTk,Image
import seaborn as sns
import pickle,os,sys,subprocess,json
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def extractValues(currentLevelLayout,valueToRemove,equalityBoolean):
    if equalityBoolean:
        idx = np.argwhere(np.all(currentLevelLayout[..., :] == valueToRemove, axis=0))
    else:
        idx = np.argwhere(np.all(currentLevelLayout[..., :] != valueToRemove, axis=0))
    currentLevelLayoutExtractedColumns = np.delete(currentLevelLayout, idx, axis=1)
    if equalityBoolean:
        idx2 = np.argwhere(np.all(currentLevelLayoutExtractedColumns[:,...] == valueToRemove, axis=1))
    else:
        idx2 = np.argwhere(np.all(currentLevelLayoutExtractedColumns[:,...] != valueToRemove, axis=1))
    currentLevelLayoutExtracted = np.delete(currentLevelLayoutExtractedColumns,idx2,axis=0)
    return currentLevelLayoutExtracted

def createWellGridButtonImages(experimentParameters,pathToButtonImages,maxNumberLevelValues):
    #maxNumberLevelValues = max(experimentParameters['numConditionLevelValues']+[experimentParameters['numColumnLevelValues']])
    width = 5
    sizeInPixels = 50
    maxColor = 255

    sns.set()
    current_palette = sns.color_palette('hls',maxNumberLevelValues)
    pal = [tuple([149/maxColor,165/maxColor,166/maxColor])]+current_palette

    rgbTriples = []
    for palTuple in pal:
        rgbTriples.append(tuple([int(i * maxColor) for i in palTuple]))
    hexColors = []
    for rgbTriple in rgbTriples:
        hexColors.append('#%02x%02x%02x' % rgbTriple)
        
    fig1 = plt.figure()
    fig1.set_size_inches(sizeInPixels,sizeInPixels)
    ax1 = fig1.add_subplot(111, aspect='equal')

    for buttonColor,i in zip(pal,range(len(pal))):
        ax1.add_patch(patches.Rectangle((0, 0), width, width,edgecolor=buttonColor,facecolor=buttonColor))
        plt.ylim(width)
        plt.xlim(width)

        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        plt.axis('off')
        plt.savefig(pathToButtonImages+'/buttonColors/'+str(i)+'.png',bbox_inches='tight', pad_inches=0,dpi=1)
        plt.cla()

    return hexColors

def produceGUIBasedIndexingCoordinates(folderName,pathToButtonImages,experimentParameters,parametersUpdatedByGridGUI,master):
    levelLayouts = []
    with open('inputFiles/levelLayouts-'+folderName+'.pkl','wb') as f:
        pickle.dump(levelLayouts,f)
    conditionLevel = 0
    numLevelValues = len(experimentParameters['allLevelValues'][experimentParameters['allLevelNames'][conditionLevel]])
    hexcolors = createWellGridButtonImages(experimentParameters,pathToButtonImages,numLevelValues)
    hexcolorsInput = '_'.join(hexcolors)
    buttonImagesPathAndHexcolorsInput = pathToButtonImages+'__'+hexcolorsInput
    master.switch_frame(plateGUIFrame,folderName,experimentParameters,parametersUpdatedByGridGUI,buttonImagesPathAndHexcolorsInput,master)

class WellGridButton(Button): #Convention is for class names to start with uppercase letters
    def __init__(self, master, images):
        super(WellGridButton, self).__init__(master, image = images[0], borderwidth = 0)
        self.images = images
        self.type = 0
        self.already_changed = False
        self.fill = 1
        #self.colorOffsetList = colorOffsetList

    def change(self):
        if self.type == self.fill:
            self.type = 0
        else:
            self.type = self.fill
        self.config(image=self.images[self.type])
        #self.config(image=self.images[self.colorOffsetList[self.type]])

    def mouse_entered(self):
        if not self.already_changed:
            self.change()
            self.already_changed = True

    def mouse_up(self):
        self.already_changed = False

class Container(Frame):
    def __init__(self, master, height, width, paired, images,folderName,experimentParameters,parametersUpdatedByGridGUI,buttonImagesPathAndHexcolorsInput,homepage):
        super(Container, self).__init__(master)
        self.folderName = folderName
        self.experimentParameters = experimentParameters
        self.parametersUpdatedByGridGUI = parametersUpdatedByGridGUI
        self.buttonImagesPathAndHexcolorsInput = buttonImagesPathAndHexcolorsInput
        self.homepage = homepage
        buttons = []
        if paired:
            oldHeight=int(height/2)
        else:
            oldHeight = height
        rowLetters = list(string.ascii_uppercase)[:height]
        for row in range(height):
            rowLetter = rowLetters[row % oldHeight]
            buttons.append([])
            for col in range(width):
                button = WellGridButton(self,images)
                button.grid(row = row, column = col)
                l1 = tk.Label(self, text=rowLetter+str(col+1).zfill(2))
                l1.config(font=("Courier", 10))
                l1.grid(row=row,column=col)
                buttons[row].append(button)

        self.buttons = buttons

        self.bind_all("<Button-1>", self.mouse_down)
        self.bind_all("<ButtonRelease-1>", self.mouse_up)
        self.bind_all("<B1-Motion>", self.mouse_motion)

        self.mouse_pressed = False
        self.tiled = False

    def mouse_down(self, e):
        self.update_containing_button(e)
        self.mouse_pressed = True

    def mouse_up(self, e):
        self.mouse_pressed = False
        for row in self.buttons:
            for button in row:
                button.mouse_up()

    def mouse_motion(self, e):
        self.update_containing_button(e)

    def update_containing_button(self, e):
        for row in self.buttons:
            for button in row:
                if self.winfo_containing(e.x_root, e.y_root) is button:
                    button.mouse_entered()
    
    def select_all(self):
        for buttonrow,row in zip(self.buttons,range(len(self.buttons))):
            for button,col in zip(buttonrow,range(len(self.buttons[0]))):
                if button.fill != 1:
                    button.fill+=1
                button.change()
    
    def unselect_all(self):
        for buttonrow,row in zip(self.buttons,range(len(self.buttons))):
            for button,col in zip(buttonrow,range(len(self.buttons[0]))):
                button.fill = 0
                button.change()

    def advance_to_specified_level_value(self,newfill):
        for buttonrow,row in zip(self.buttons,range(len(self.buttons))):
            for button,col in zip(buttonrow,range(len(self.buttons[0]))):
                button.fill = newfill

    def create_full_level_value_layout(self):
        currentlevelLayout = np.zeros([len(self.buttons),len(self.buttons[0])])
        for buttonrow,row in zip(self.buttons,range(len(self.buttons))):
            for button,col in zip(buttonrow,range(len(self.buttons[0]))):
                currentlevelLayout[row,col] = button.type
        with open('inputFiles/fullLevelLayout.pkl','wb') as f:
            pickle.dump(currentlevelLayout,f)
        if not self.tiled:
            with open('inputFiles/tile-levelLayout.pkl','wb') as f:
                pickle.dump(currentlevelLayout,f)
        levelLayouts = pickle.load(open('inputFiles/levelLayouts-'+self.folderName+'.pkl','rb'))
        levelLayouts.append(currentlevelLayout)
        currentLevelLayout = pickle.load(open('inputFiles/tile-levelLayout.pkl','rb'))
        #Strip out all zero columns/rows to make next plate layout
        currentLevelLayoutNoZeros = extractValues(currentLevelLayout,0,True)
        self.parametersUpdatedByGridGUI['currentPlateDimensions'] = currentLevelLayoutNoZeros.shape
        self.parametersUpdatedByGridGUI['numLevelsUnparsed'] -= 1
        with open('inputFiles/gui-parametersUpdatedByGridGUI.pkl','wb') as f:
            pickle.dump(self.parametersUpdatedByGridGUI,f)
        with open('inputFiles/levelLayouts-'+self.folderName+'.pkl','wb') as f:
            pickle.dump(levelLayouts,f)
        print(levelLayouts)
        if self.parametersUpdatedByGridGUI['numLevelsUnparsed'] != 0:
            self.homepage.switch_frame(plateGUIFrame,self.folderName,self.experimentParameters,self.parametersUpdatedByGridGUI,self.buttonImagesPathAndHexcolorsInput,self.homepage)
        else:
            self.parametersUpdatedByGridGUI['currentPlateDimensions'] = self.experimentParameters['overallPlateDimensions']
            self.parametersUpdatedByGridGUI['numLevelsUnparsed'] = self.experimentParameters['numAllLevels']
            with open('inputFiles/gui-parametersUpdatedByGridGUI.pkl','wb') as f:
                pickle.dump(self.parametersUpdatedByGridGUI,f)
            self.homepage.switch_frame(self.homepage.homepage)

    def tile_level(self,tilingDirection):
        self.tiled = True
        levelLayout = np.zeros([len(self.buttons),len(self.buttons[0])])
        for buttonrow,row in zip(self.buttons,range(len(self.buttons))):
            for button,col in zip(buttonrow,range(len(self.buttons[0]))):
                levelLayout[row,col] = button.type
        idx = np.argwhere(np.all(levelLayout[..., :] == 0, axis=0))
        currentLevelLayoutNoZeroColumns = np.delete(levelLayout, idx, axis=1)
        idx2 = np.argwhere(np.all(currentLevelLayoutNoZeroColumns[:,...] == 0, axis=1))
        tile = np.delete(currentLevelLayoutNoZeroColumns,idx2,axis=0)
        tileShape = tile.shape
        tiledLevelLayout = np.zeros([len(self.buttons),len(self.buttons[0])])
        for row in range(tileShape[0]):
            for col in range(tileShape[1]):
                newfill = 1 
                if tilingDirection == 'horizontal':
                    for repetition in range(int(levelLayout.shape[1]/tileShape[1])):
                        tiledLevelLayout[row,col+repetition*tileShape[1]] = newfill
                        button = self.buttons[row][col+repetition*tileShape[1]]
                        if button.fill != newfill:
                            button.fill = newfill
                            button.change()
                        newfill+=1
                else:
                    for repetition in range(int(levelLayout.shape[0]/tileShape[0])):
                        tiledLevelLayout[row+repetition*tileShape[0],col] = newfill
                        button = self.buttons[row+repetition*tileShape[0]][col]
                        if button.fill != newfill:
                            button.fill = newfill
                            button.change()
                        newfill+=1
        print(tiledLevelLayout)
        with open('inputFiles/tile-levelLayout.pkl','wb') as f:
            pickle.dump(tile,f)

class plateGUIFrame(tk.Frame):
    def __init__(self, master,folderName,experimentParameters,parametersUpdatedByGridGUI,buttonImagesPathAndHexcolorsInput,homepage):
        tk.Frame.__init__(self, master)
        hexcolors = buttonImagesPathAndHexcolorsInput.split('__')[1].split('_')[1:]
        pathToButtonImages = buttonImagesPathAndHexcolorsInput.split('__')[0]
        conditionLevel = 0
        numLevelValues = len(experimentParameters['allLevelValues'][experimentParameters['allLevelNames'][experimentParameters['numAllLevels']-parametersUpdatedByGridGUI['numLevelsUnparsed']]])
        hexcolors = createWellGridButtonImages(experimentParameters,pathToButtonImages,numLevelValues)
        hexcolorsInput = '_'.join(hexcolors)
        buttonImagesPathAndHexcolorsInput = pathToButtonImages+'__'+hexcolorsInput
        folderName = os.getcwd().split('/')[-1]
        #experimentParameters = json.load(open('inputFiles/experimentParameters-'+folderName+'.json','r'))
        #parametersUpdatedByGridGUI = pickle.load(open('inputFiles/gui-parametersUpdatedByGridGUI.pkl','rb'))
        dim = parametersUpdatedByGridGUI['currentPlateDimensions']
        currentIndex = experimentParameters['numAllLevels'] - parametersUpdatedByGridGUI['numLevelsUnparsed']
        currentConditionLevelName = experimentParameters['allLevelNames'][currentIndex]
        currentConditionLevelValues = experimentParameters['allLevelValues'][currentConditionLevelName]
        #root = Tk()
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        #Load in all level value button images
        images = {}
        for fileName in os.listdir(pathToButtonImages+'/buttonColors'):
            images[int(fileName.split('.')[0])] = ImageTk.PhotoImage(file=pathToButtonImages+'/buttonColors/'+fileName)
        if experimentParameters['paired'] and currentIndex == 0:
            grid = Container(mainWindow, dim[0]*2, dim[1],experimentParameters['paired'],images,folderName,experimentParameters,parametersUpdatedByGridGUI,buttonImagesPathAndHexcolorsInput,homepage)
        else:
            grid = Container(mainWindow, dim[0], dim[1],experimentParameters['paired'],images,folderName,experimentParameters,parametersUpdatedByGridGUI,buttonImagesPathAndHexcolorsInput,homepage)
        grid.pack()

        tk.Label(mainWindow, text=currentConditionLevelName,font='TkDefaultFont 16 bold').pack()

        levelValueFrame = Frame(mainWindow)
        rblist = []
        v=IntVar(value=1)
        def store():
            grid.advance_to_specified_level_value(v.get())
            subprocess.run('rm '+pathToButtonImages+'/buttonColors/*.png',shell=True)
        def selectAll():
            grid.select_all()
            for rb in rblist:
                v.set(grid.buttons[0][0].fill)
        def unselectAll():
            grid.unselect_all()
            for rb in rblist:
                v.set(0)
        for currentLevelValue,fillVal in zip(currentConditionLevelValues,range(1,len(currentConditionLevelValues)+1)):
            rb = tk.Radiobutton(mainWindow,text=currentLevelValue,bg=hexcolors[fillVal],variable=v,value=fillVal,command=lambda: store())
            rblist.append(rb)
            rb.pack(in_=levelValueFrame,side=LEFT,fill=BOTH,expand=True)
        levelValueFrame.pack(fill=BOTH,expand=True)
        actionButtons = tk.Frame(self)
        actionButtons.pack()
        tileCommands = tk.Frame(actionButtons)
        tileCommands.pack()
        tk.Button(tileCommands, text="Tile Horizontally",command=lambda:grid.tile_level('horizontal')).pack(side=LEFT)
        tk.Button(tileCommands, text="Tile Vertically",command=lambda:grid.tile_level('vertical')).pack(side=LEFT)
        tk.Button(tileCommands, text="Select All",command=lambda:selectAll()).pack(side=LEFT)
        tk.Button(tileCommands, text="Unselect All",command=lambda:unselectAll()).pack(side=LEFT)
        endCommands = tk.Frame(actionButtons)
        endCommands.pack()
        tk.Button(endCommands, text="Next Level",font="-weight bold",command=grid.create_full_level_value_layout).pack(side=LEFT)
        def quitCommand():
            parametersUpdatedByGridGUI['currentPlateDimensions'] = experimentParameters['overallPlateDimensions']
            parametersUpdatedByGridGUI['numLevelsUnparsed'] = experimentParameters['numAllLevels']
            with open('inputFiles/gui-parametersUpdatedByGridGUI.pkl','wb') as f:
                pickle.dump(parametersUpdatedByGridGUI,f)
            quit()
        tk.Button(endCommands, text="Quit",command=lambda: quitCommand()).pack(side=LEFT)

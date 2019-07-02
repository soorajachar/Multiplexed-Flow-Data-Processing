#!/usr/bin/env python3
import numpy as np
import pickle,os,string,json
from tkinter import *
import tkinter as tk
from functools import reduce
from PIL import ImageTk,Image

exitBoolean = False

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
    def __init__(self, master, height, width, paired, images):
        super(Container, self).__init__(master)

        buttons = []
        if paired:
            oldHeight=int(height/2)
        else:
            oldHeight = height
        rowLetters = list(string.ascii_uppercase)[:height]
        for row in range(height):
            rowLetter = rowLetters[row % oldHeight]
            print(rowLetter)
            print(oldHeight)
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

    def advance_to_specified_level_value(self,newfill):
        for buttonrow,row in zip(self.buttons,range(len(self.buttons))):
            for button,col in zip(buttonrow,range(len(self.buttons[0]))):
                button.fill = newfill

    def create_full_level_value_layout(self):
        levelLayout = np.zeros([len(self.buttons),len(self.buttons[0])])
        for buttonrow,row in zip(self.buttons,range(len(self.buttons))):
            for button,col in zip(buttonrow,range(len(self.buttons[0]))):
                levelLayout[row,col] = button.type
        with open('inputFiles/fullLevelLayout.pkl','wb') as f:
            pickle.dump(levelLayout,f)
        if not self.tiled:
            with open('inputFiles/tile-levelLayout.pkl','wb') as f:
                pickle.dump(levelLayout,f)
        print(levelLayout)
        with open('inputFiles/gui-exitBoolean.pkl','wb') as f:
            pickle.dump(exitBoolean,f)
        quit()
    
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
        
def main():
    global folderName
    
    hexcolors = str(sys.argv[1]).split('__')[1].split('_')[1:]
    pathToButtonImages = str(sys.argv[1]).split('__')[0]
    folderName = os.getcwd().split('/')[-1]
    experimentParameters = json.load(open('inputFiles/experimentParameters-'+folderName+'.json','r'))
    parametersUpdatedByGridGUI = pickle.load(open('inputFiles/gui-parametersUpdatedByGridGUI.pkl','rb'))
    dim = parametersUpdatedByGridGUI['currentPlateDimensions']
    currentIndex = experimentParameters['numAllLevels'] - parametersUpdatedByGridGUI['numLevelsUnparsed']
    currentConditionLevelName = experimentParameters['allLevelNames'][currentIndex]
    currentConditionLevelValues = experimentParameters['allLevelValues'][currentConditionLevelName]
    root = Tk()
    
    #Load in all level value button images
    images = {}
    for fileName in os.listdir(pathToButtonImages+'/buttonColors'):
        images[int(fileName.split('.')[0])] = ImageTk.PhotoImage(file=pathToButtonImages+'/buttonColors/'+fileName)
    if experimentParameters['paired'] and currentIndex == 0:
        grid = Container(root, dim[0]*2, dim[1],experimentParameters['paired'],images)
    else:
        grid = Container(root, dim[0], dim[1],experimentParameters['paired'],images)
    grid.pack()

    tk.Label(root, text=currentConditionLevelName,font='TkDefaultFont 16 bold').pack()

    levelValueFrame = Frame(root)
    rblist = []
    v=IntVar()
    def store():
        grid.advance_to_specified_level_value(v.get())
    for currentLevelValue,fillVal in zip(currentConditionLevelValues,range(1,len(currentConditionLevelValues)+1)):
        rb = tk.Radiobutton(root,text=currentLevelValue,bg=hexcolors[fillVal-1],variable=v,value=fillVal,command=lambda: store())
        rblist.append(rb)
        rb.pack(in_=levelValueFrame,side=LEFT,fill=BOTH,expand=True)
    levelValueFrame.pack(fill=BOTH,expand=True)
    
    actionButtons = Frame(root)
    actionButtons.pack()
    tk.Button(root, text="Next Level",command=grid.create_full_level_value_layout).pack(in_=actionButtons,side=LEFT)
    tk.Button(root, text="Tile Horizontally",command=lambda:grid.tile_level('horizontal')).pack(in_=actionButtons,side=LEFT)
    tk.Button(root, text="Tile Vertically",command=lambda:grid.tile_level('vertical')).pack(in_=actionButtons,side=LEFT)
    def quitCommand():
        exitBoolean = True
        with open('inputFiles/gui-exitBoolean.pkl','wb') as f:
            pickle.dump(exitBoolean,f)
        parametersUpdatedByGridGUI['currentPlateDimensions'] = experimentParameters['overallPlateDimensions']
        parametersUpdatedByGridGUI['numLevelsUnparsed'] = experimentParameters['numAllLevels']
        with open('inputFiles/gui-parametersUpdatedByGridGUI.pkl','wb') as f:
            pickle.dump(parametersUpdatedByGridGUI,f)
        quit()
    tk.Button(root, text="Quit",command=lambda: quitCommand()).pack(in_=actionButtons,side=LEFT)

    root.mainloop()

main()

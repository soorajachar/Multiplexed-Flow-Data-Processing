#!/usr/bin/env python3  
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import numpy as np
import seaborn as sns
import pandas as pd
import os,sys,pickle
idx = pd.IndexSlice
sys.path.insert(0, '../../../programs/dataprocessing/')
from miscFunctions import returnGates,find_nearest,reindexDataFrame

class proliferation_gate_lines:
    def __init__(self, ax, X,Ymin,Ymax, gates,numGates):
        self.ax = ax
        self.c = ax.get_figure().canvas
        self.X = X
        self.gates= gates
        self.gateVals = gates


        gateLengths = [y-x for x, y in zip(gates[:-1], gates[1:])][1:]
        self.gateLengths = gateLengths 

        self.numGates = numGates
        
        x = [X, X]
        y = [Ymin, Ymax]
        
        self.line = lines.Line2D(x, y, picker=5,linestyle=':',color='r')
        self.ax.add_line(self.line)
        self.childlines = []
        for gate in range(self.numGates):
            newX = self.gates[gate]
            childline = lines.Line2D([newX,newX], y,linestyle=':',color='k',picker=5)
            self.childlines.append(childline)
        for childline in self.childlines:
            self.ax.add_line(childline)
        self.c.draw_idle()
        self.sid = self.c.mpl_connect('pick_event', self.clickOnParentLine)
        childsids = []
        for childline,i in zip(self.childlines,range(self.numGates)):
            childsids.append(self.c.mpl_connect('pick_event', self.clickOnChildLines))
        self.childsids = childsids

    def clickOnParentLine(self, event):
        if event.artist == self.line:
            print("line selected ", event.artist)
            self.follower = self.c.mpl_connect("motion_notify_event", self.followParentMouse)
            self.releaser = self.c.mpl_connect("button_press_event", self.releaseParentOnClick)
    
    def clickOnChildLines(self, event):
        if event.artist in self.childlines:
            print("line selected ", event.artist)
            self.currentArtist = event.artist
            self.follower = self.c.mpl_connect("motion_notify_event", self.followChildMouse)
            self.releaser = self.c.mpl_connect("button_press_event", self.releaseChildOnClick)
    
    def followChildMouse(self, event):
        for childline,gate in zip(self.childlines,range(self.numGates)):
            if self.childlines[gate] == self.currentArtist:
                self.childlines[gate].set_xdata([event.xdata, event.xdata])
        self.c.draw_idle()
    
    def releaseChildOnClick(self, event):
        newGates = []
        for childline in self.childlines:
            newGates.append(childline.get_xdata()[0])
        self.gates = newGates
        self.c.mpl_disconnect(self.releaser)
        self.c.mpl_disconnect(self.follower)
        print(self.getAllX())

    def followParentMouse(self, event):
        self.line.set_xdata([event.xdata, event.xdata])
        newGates = []
        for childline,gate in zip(self.childlines,range(self.numGates)):
            newX = (event.xdata-self.X)+self.gates[gate]
            childline.set_xdata([newX,newX])
            newGates.append(newX)
        self.gateVals = newGates
        self.c.draw_idle()
    
    def releaseParentOnClick(self, event):
        self.X = self.line.get_xdata()[0]
        newGates = []
        for childline in self.childlines:
            newGates.append(childline.get_xdata()[0])
        self.gates = newGates
        self.c.mpl_disconnect(self.releaser)
        self.c.mpl_disconnect(self.follower)
        print(self.getAllX())

    def getParentX(self):
        return self.X
    
    def getChildX(self):
        childGates = []
        for childline in self.childlines:
            childGates.append(childline.get_xdata()[0])
        return childGates

    def getAllX(self):
        return [self.getParentX()]+self.getChildX()
"""
fig = plt.figure()
ax = fig.add_subplot(111)
Tline = proliferation_gate_lines(ax,0.5,-1,1,[0.4,0.3,0.1],3)
plt.show()
"""

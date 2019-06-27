#!/usr/bin/env python3  
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines

class proliferation_gate_lines:
    def __init__(self, ax, X,Ymin,Ymax, gates,numGates):
        self.ax = ax
        self.c = ax.get_figure().canvas
        self.X = X
        self.gates= gates
        
        gateLengths = [y-x for x, y in zip(gates[:-1], gates[1:])][1:]
        self.gateLengths = gateLengths 

        self.numGates = numGates
        
        x = [X, X]
        y = [Ymin, Ymax]
        
        self.line = lines.Line2D(x, y, picker=5,linestyle=':',color='k')
        self.ax.add_line(self.line)
        self.childlines = []
        for gate in range(self.numGates):
            newX = self.gates[gate]
            childline = lines.Line2D([newX,newX], y,linestyle=':',color='k')
            self.childlines.append(childline)
        for childline in self.childlines:
            self.ax.add_line(childline)
        self.c.draw_idle()
        self.sid = self.c.mpl_connect('pick_event', self.clickonline)

    def clickonline(self, event):
        if event.artist == self.line:
            print("line selected ", event.artist)
            self.follower = self.c.mpl_connect("motion_notify_event", self.followmouse)
            self.releaser = self.c.mpl_connect("button_press_event", self.releaseonclick)

    def followmouse(self, event):
        self.line.set_xdata([event.xdata, event.xdata])
        for childline,gate in zip(self.childlines,range(self.numGates)):
            newX = (event.xdata-self.X)+self.gates[gate]
            childline.set_xdata([newX,newX])
        self.c.draw_idle()
    
    def releaseonclick(self, event):
        self.X = self.line.get_xdata()[0]

        self.c.mpl_disconnect(self.releaser)
        self.c.mpl_disconnect(self.follower)
    
    def getX(self):
        return self.XorY
"""
fig = plt.figure()
ax = fig.add_subplot(111)
Tline = proliferation_gate_lines(ax,0.5,-1,1,[0.4,0.3,0.1],3)
plt.show()
"""

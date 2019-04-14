#!/usr/bin/env python3  
import matplotlib.pyplot as plt
import matplotlib.lines as lines

class draggable_lines:
    def __init__(self, ax, kind, XorY):
        self.ax = ax
        self.c = ax.get_figure().canvas
        self.o = kind
        self.XorY = XorY

        if kind == "h":
            x = [-1, 1]
            y = [XorY, XorY]

        elif kind == "v":
            x = [XorY, XorY]
            y = [-1, 1]
        self.line = lines.Line2D(x, y, picker=5,linestyle=':',color='k')
        self.ax.add_line(self.line)
        self.c.draw_idle()
        self.sid = self.c.mpl_connect('pick_event', self.clickonline)
        self.background = self.c.copy_from_bbox(self.line.axes.bbox)

    def clickonline(self, event):
        if event.artist == self.line:
            print("line selected ", event.artist)
            self.follower = self.c.mpl_connect("motion_notify_event", self.followmouse)
            self.releaser = self.c.mpl_connect("button_press_event", self.releaseonclick)

    def followmouse(self, event):
        if self.o == "h":
            self.line.set_ydata([event.ydata, event.ydata])
        else:
            self.line.set_xdata([event.xdata, event.xdata])
        axes = self.line.axes
        canvas = self.line.figure.canvas
        #print('wat')
        canvas.restore_region(self.background)
        axes.draw_artist(self.line)
        canvas.blit(axes.bbox)

    def releaseonclick(self, event):
        if self.o == "h":
            self.XorY = self.line.get_ydata()[0]
        else:
            self.XorY = self.line.get_xdata()[0]

        print (self.XorY)

        self.c.mpl_disconnect(self.releaser)
        self.c.mpl_disconnect(self.follower)

    def getX(self):
        return self.XorY

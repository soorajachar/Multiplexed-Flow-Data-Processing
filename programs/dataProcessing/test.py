#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
plt.switch_backend('QT4Agg') #default on my system

fig, ax = plt.subplots()
    
if plotType == 'categorical':
    axes = ['Y Axis']
    axbox1 = plt.axes([0.1, 0.55, 0.4, 0.075])
    text_box = TextBox(axbox, axes, initial=axes[0])
else:
    axes = ['X Axis','Y Axis']
    axbox1 = plt.axes([0.1, 0.55, 0.4, 0.075])
    text_box = TextBox(axbox, axes, initial=axes[0])
    axbox2 = plt.axes([0.6, 0.55, 0.4, 0.075])
    text_box2 = TextBox(axbox, axes, initial=axes[1])
"""
class Index5(object):
    def OK(self, event):
        print(text_box2.on_change(submit))
        with open('semiProcessedData/radioVals-axes.pkl','wb') as f:
            pickle.dump(textValues,f)
    def Quit(self, event):
        sys.exit(0)    

callback = Index4()
axOK = plt.axes([buttonXStart, buttonYStart, buttonWidth, buttonLength])
axQuit = plt.axes([buttonXStart+buttonWidth+0.01,buttonYStart, buttonWidth, buttonLength])
bOK = Button(axOK, 'OK')
bOK.on_clicked(callback.OK)
bQuit = Button(axQuit, 'Quit')
bQuit.on_clicked(callback.Quit)

plt.show()
"""

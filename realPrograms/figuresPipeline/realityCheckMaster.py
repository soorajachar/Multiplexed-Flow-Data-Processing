import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

def createRealityCheckPlots(secondPath,crossValidateNumber,modelName,trainingString,dfToCrossValidate,R_enc,X_norm):
    numPlots = int((len(dfToCrossValidate.columns)-1)/4)+1
    fig1,axes = plt.subplots(numPlots,4,figsize = (15,3.75*numPlots))
    for i,ax in enumerate(axes.flatten()):
        if i > 6:
            ax.axis('off')
            continue
        ax.scatter(R_enc[:,i],X_norm[:,i],marker='+')
        R,p = pearsonr(R_enc[:,i],X_norm[:,i])
        ax.set_title(dfToCrossValidate.columns[i]+"; R=%.2f; p = %.2f"%(R,p))
        ylim = ax.get_ylim()
        ax.set_xlim(ylim)
    plt.suptitle("Predicted (x-axis) vs. real (y-axis). As linear as possible")
    fig1.savefig(secondPath+'realityChecks/realityCheck-%d-on-%s-%s.png'%(crossValidateNumber,modelName,trainingString))
    plt.close(fig1)

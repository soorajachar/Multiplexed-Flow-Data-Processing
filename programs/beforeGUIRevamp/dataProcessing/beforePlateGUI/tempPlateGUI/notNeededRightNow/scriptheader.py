#!/usr/bin/env python3 
import json,pickle,math,matplotlib,sys,os,string,subprocess
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

idx = pd.IndexSlice
kwargs = pickle.load(open('gui-kwargs.pkl','rb'))
print(kwargs)
kwargs['paired'] = True
kwargs['numLevelsUnparsed'] = kwargs['numAllLevels']
with open('gui-kwargs.pkl','wb') as f:
    pickle.dump(kwargs,f)
levelLayout = pickle.load(open('levelLayout.pkl','rb'))

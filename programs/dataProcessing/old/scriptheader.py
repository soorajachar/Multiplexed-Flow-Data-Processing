#!/usr/bin/env python3 
import json,pickle,math,matplotlib,sys,os,string
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
idx = pd.IndexSlice
import itertools

df = pickle.load(open('../../experiments/20190701-TCellNumber_CAR_Timeseries_1/postProcessedData/dimensionalReductions/dimensionalReduction-cyt-cell-all.pkl','rb'))
print(df.reset_index())

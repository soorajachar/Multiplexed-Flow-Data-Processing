#!/usr/bin/env python3
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import markers
import numpy as np
import seaborn as sns
import pandas as pd
import pickle,os,math,sys,itertools,re
from matplotlib.widgets import RadioButtons,Button,CheckButtons,TextBox

def scaleXandY2D(ax,plotOptions):

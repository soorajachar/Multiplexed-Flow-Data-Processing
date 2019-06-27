#!/usr/bin/env python3 
import json,pickle,math,matplotlib,sys,os,string
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from processProliferationData import returnTicks
from sklearn.neighbors import KernelDensity
from scipy.stats import gaussian_kde
from statsmodels.nonparametric.kde import KDEUnivariate
from statsmodels.nonparametric.kernel_density import KDEMultivariate
sys.path.insert(0, '../figuresPipeline/')
from facetPlotHeatmaps import draw_faceted_heatmap,returnHeatmapAspectRatios

folderName = '20190404-CD25MutantTimeSeries_OT1_Timeseries_2'
idx=pd.IndexSlice
data = pickle.load(open('../../experiments/'+folderName+'/semiProcessedData/initialSingleCellDf-channel-'+folderName+'.pkl','rb')).loc[idx['TCells','WT','A8','1uM',:],'TCell_Gate']
bnw = 10 
def kde_scipy(x, x_grid,**kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    kde = gaussian_kde(x, bw_method=bnw / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)

kde_funcs = [kde_statsmodels_u, kde_scipy]
kde_funcnames = ['Statsmodels-U', 'Scipy', 'Scikit-learn']

x_grid = np.linspace(0, 1000, 1000)
i=0

zerolist = []
for timepoint in pd.unique(data.index.get_level_values('Time')):
    x = data.loc[idx[:,:,:,:,timepoint]].values.ravel().astype(float)
    pdf = kde_funcs[i](x, x_grid)

    diff = np.gradient(pdf)
    sdiff = np.sign(diff)
    zc = np.where(sdiff[:-1] != sdiff[1:])[0][-2]
    zerolist.append(zc)

plottingDf = data.to_frame('GFI').reset_index()
facetplotkwargs = {'col':'Time','col_wrap':6}
fg = sns.FacetGrid(plottingDf,sharey=False,**facetplotkwargs)
fg.map(sns.kdeplot,'GFI',shade=True,bw=bnw)
for zc,ax in zip(zerolist,fg.fig.get_axes()):
    ax.axvline(x=zc,linestyle=':',color='k')
plt.show()

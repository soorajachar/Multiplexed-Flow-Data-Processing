#!/usr/bin/env python3
"""Logicle transform"""
from special_funcs import productlog
from scipy.optimize import fsolve, brentq
from scipy import interpolate
from numpy import arange, exp, vectorize, log, min, max, sign, concatenate, zeros
from numpy.random import normal, lognormal, shuffle
import pylab
import time
import math
import numpy as np

def EH(x, y, b, d, r):
    e = float(d)/r
    sgn = sign(x)
    return sgn*10**(sgn*e*x) + b*e*x - sgn - y

def hyperlog0(y, b, d, r):
    return brentq(EH, -10**6, 10**6, (y, b, d, r))
hyperlog0 = vectorize(hyperlog0)

def hyperlog(y, b, d, r, order=2, intervals=1000.0):
    ub = log(max(y)+1-min(y))
    xx = exp(arange(0, ub, ub/intervals))-1+min(y)
    yy = hyperlog0(xx, b, d, r)
    t = interpolate.splrep(xx, yy, k=order)
    return interpolate.splev(y, t)

def S(x, y, T, m, w):
    p = w/(2*productlog(0.5*exp(-w/2)*w))
    sgn = sign(x-w)
    xw = sgn*(x-w)
    return sgn*T*exp(-(m-w))*(exp(xw)-p**2*exp(-xw/p)+p**2-1) - y

def logicle0(y, T, m, r):
    if r>=0:
        return m+log((exp(-m)*T+y)/T)
    else:
        w = (m-log(T/abs(r)))/2
        return brentq(S, -100, 100, (y, T, m, w))
logicle0 = vectorize(logicle0)

def logicle(y, T, m, r, order=2, intervals=1000.0):
    ub = log(max(y)+1-min(y)) #upper bound
    xx = exp(arange(0, ub, ub/intervals))-1+min(y)
    yy = logicle0(xx, T, m, r)
    t = interpolate.splrep(xx, yy, k=order)
    return interpolate.splev(y, t)

def quantile(x, n):
    try:
        return sorted(x)[int(n*len(x))]
    except IndexError:
        return 0
if __name__ == '__main__':
    d1 = normal(0, 50, (50000))
    d2 = lognormal(8, 1, (50000))
    d3 = concatenate([d1, d2])

    T = 262144 #Top of datarange (max value; 2^18 (18 bit data range); same as flowjo
    d = 4 #Not sure what flowjo has but I believe it's 5 (number of positive decades displayed in data; 5 means 10^5 is larges dataval; logicle(T) = M
    m = d*log(10) #convert decades to asymptotic decades (e^x); 0 to m is the data range acceptable for the transformation
    r = quantile(d3[d3<0], 0.05) #lowest negative value
    w = (m-log(T/abs(r)))/2 #Width of linear reigion (centered around 0); negative values appear from 0 to w, positive from w to m
#span(1) = T / 10^(W*2+M);

#%% Generate the X tick marks and labels
#options = optimset('Display','off','TolFun',1e-4);
#LINEAR_TICKMARK_SPACING = 50;
#END_LINEAR_SPACING = logicleFun(x0, a, b, c, d, f, x1, 0);
#%linear portion
#sTicks = (floor(span(1)/10)*10):LINEAR_TICKMARK_SPACING: ...
#    floor(END_LINEAR_SPACING/10)*10; 
#%log portion
#decade = floor(log10(max(sTicks)));
#startLogSpace = 10^decade+ max(sTicks);
#sTicks = [sTicks startLogSpace:10^decade:10^(decade+1)];
#decade = decade+1;
#while max(sTicks)<span(2)
#    sTicks = [sTicks 2*10^decade:10^decade:10^(decade+1)];
#    decade = decade+1;
#end
#sTicks(sTicks>span(2))=[]; %remove extra sTicks
#for i=1:length(sTicks)
#    if (sTicks(i)<=100)
#        if(rem(sTicks(i),100)==0)
#            xLabel{i}=num2str(sTicks(i));
#        else
#            xLabel{i}=' ';
#        end
#    elseif floor(log10(sTicks(i)))==log10(sTicks(i))
#        xLabel{i}=num2str(sTicks(i));
#    else
#        xLabel{i}=' ';
#    end
#    S=sTicks(i);
#    xTickGuess = xtickguess(a,b,c,d,f,x1,S);
#    xTicks(i) = fsolve(@(X)logicleFun(X,a,b,c,d,f,x1,S),xTickGuess,options);
#end
#set(gca,'XTick',xTicks);
#set(gca,'XTickLabel',xLabel);
    #To construct appropriate scales, need two regimes; linear from -w to w on x regime (0 to 2w on logicle regime) and w to T on x regime (2w to m on logicle regime)
    linearTickMarkSpacing = 50
    endLinearSpacing = w
    #linearRegimeNegative = np.arrange(0,w,step=linearTickMarkSpacing)
    #linearRegimePositive = np.arrange(w,2*w,step=linearTickMarkSpacking)
    #logatithmicRegime = np.logspace(2*w,T,num = 3)
    pylab.clf()
    pylab.figtext(0.5, 0.94, 'Logicle transform with r=%.2f, d=%d and T=%d\nData is normal(0, 50, 50000) + lognormal(8, 1, 50000)' % (r, d, T),
                  va='center', ha='center', fontsize=12)

    pylab.subplot(3,1,1)
    x = arange(0, m, 0.1)
    pylab.plot(x, S(x, 0, T, m, w))
    locs, labs = pylab.xticks()
    #pylab.xticks([])
    #pylab.yticks([])
    pylab.ylabel('Inverse logicle')

    pylab.subplot(3,1,2)
    pylab.hist(d3, 1250)
    locs, labs = pylab.xticks()
    #pylab.xticks([])
    #pylab.yticks([])
    pylab.ylabel('Raw data')

    pylab.subplot(3,1,3)
    pylab.hist(logicle(d3, T, m, r), 1250)
    locs, labs = pylab.xticks()
    print(locs)
    print(labs)
    
    #pylab.xticks([])
    pylab.yticks([])
    pylab.ylabel('Data after transform')

    # pylab.savefig('logicle.png')
    pylab.show()

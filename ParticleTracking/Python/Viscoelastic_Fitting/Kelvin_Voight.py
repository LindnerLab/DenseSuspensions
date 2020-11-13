# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 14:39:55 2020

@author: Charly
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
import scipy.optimize as opt

def Kelvin_Voight(t, E, sig0, eta):
    return (sig0/E)*(1 - np.exp(-(E/eta)*t))

def Burgers(t, sig0, eta1, eta2, E1):
    return (sig0/eta1)*t + (eta2/E1)*(1 - np.exp(-(E1/eta2)*t))

def func(t, a, b, c, d):
    return a * (1 - np.exp(-b*t)) + c * (1 - np.exp(-d*t))

tracked = pd.read_pickle(r'F:\Lars\Oscillatory Compression\20201103 Creep\Compression_26Start_125End_25_Back\Preprocessed\V1\tracked.pkl')
nParticles = len(tracked.particle.unique())
nFrames = len(tracked.frame.unique())

avg_strain = np.empty((nFrames,1))
for i in range(nFrames):
    avg_strain[i] = np.mean(tracked[tracked.frame == i].Hencky_strain)
    
plt.figure(dpi = 500)
plt.plot(np.arange(0,297,0.1), avg_strain[30:],
         marker='o',
         linestyle='none')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Time (s)')
plt.ylabel('Hencky strain (-)')
plt.show()

xdata = np.arange(0,100.5,0.1).reshape((1005,))
ydata = avg_strain[30:1035].reshape((1005,))
popt, pcov = opt.curve_fit(func, xdata, ydata)
plt.figure(dpi=500)
plt.plot(xdata, ydata,
         marker='o',
         linestyle='none')
yfit = func(xdata, *popt)
plt.plot(xdata, yfit)
# plt.xscale('log')
# plt.yscale('log')
plt.xlabel('Time (s)')
plt.ylabel('Hencky strain (-)')
plt.show()

x = np.exp(np.arange(-1, 5.6, 0.1)).reshape((66,))
y = func(x, 1, 0.01, 0, 0)
plt.figure(dpi=500)
plt.plot(x,y)
plt.xscale('log')
plt.yscale('log')
plt.show()

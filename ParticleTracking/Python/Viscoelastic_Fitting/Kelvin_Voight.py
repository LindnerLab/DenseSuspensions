# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 14:39:55 2020

Contributors : Lars Kool
Affiliations : Laboratoire Physique et Méchanique des Milieux Hétérogènes
(PMMH), ESPCI, Paris, France
               
This project has received funding from the European Union’s Horizon 2020
research and innovation programme under the Marie Skłodowska-Curie grant
agreement No. 813162
More info on this programme: https://caliper-itn.org/
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

def two_timescale_compression(t, a, b, c, d):
    return a * (1 - np.exp(-b*t)) + c * (1 - np.exp(-d*t))

def two_timescale_compression_drift(t, a, b, c, d, e):
    return a * (1 - np.exp(-b*t)) + c * (1 - np.exp(-d*t)) + e*t

def two_timescale_decompression(t, a, b, c, d):
    return a * np.exp(-b*t) + c * np.exp(-d*t)

def two_timescale_decompression_drift(t, a, b, c, d, e, f):
    return a * np.exp(-b*t) + c * np.exp(-d*t) + e*t + f

tracked_compressed = pd.read_pickle(r'E:\Lars\Oscillatory Compression\20201103 Creep\Compression_26Start_125End_25Back\Preprocessed\V1\tracked.pkl')
tracked_decompressed = pd.read_pickle(r'E:\Lars\Oscillatory Compression\20201103 Creep\Decompression_125Start_26End_25Back\Preprocessed\V1\tracked.pkl')
nParticles_compressed = len(tracked_compressed .particle.unique())
nParticles_decompressed = len(tracked_compressed .particle.unique())
nFrames = len(tracked_compressed.frame.unique())

avg_strain_compressed = np.empty((nFrames,1))
avg_strain_decompressed = np.empty((nFrames,1))
for i in range(nFrames):
    avg_strain_compressed [i] = np.mean(tracked_compressed[tracked_compressed.frame == i].Hencky_strain)
    avg_strain_decompressed [i] = np.mean(tracked_decompressed[tracked_decompressed.frame == i].Hencky_strain)
    
plt.figure(dpi = 500)
plt.plot(np.arange(0,300,0.1), avg_strain_compressed,
         marker='o',
         linestyle='none')
plt.xlabel('Time (s)')
plt.ylabel('Hencky strain (-)')
plt.show()

xdata_compressed = np.arange(0,297,0.1).reshape((2970,))
ydata_compressed = avg_strain_compressed[30:].reshape((2970,))
fit_compressed, _ = opt.curve_fit(two_timescale_compression_drift,
                                  xdata_compressed,
                                  ydata_compressed,
                                  p0=[0.0694, 0.1421, 0.1671, 0.0140, 1E-4])
plt.figure(dpi=500)
plt.plot(xdata_compressed, ydata_compressed,
          marker='o',
          markersize=1.5,
          linestyle='none')
yfit_compressed = two_timescale_compression_drift(xdata_compressed, *fit_compressed)
plt.plot(xdata_compressed, yfit_compressed,
         linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Hencky strain (-)')
plt.show()

xdata_decompressed = np.arange(0,299,0.1).reshape((2990,))
ydata_decompressed = avg_strain_decompressed[10:].reshape((2990,))
fit_decompressed, _ = opt.curve_fit(two_timescale_decompression_drift,
                                    xdata_decompressed,
                                    ydata_decompressed,
                                    p0=[0.0694, 0.1421, 0.1671, 0.0140, -1E-4, 0.02])
plt.figure(dpi=500)
plt.plot(xdata_decompressed, ydata_decompressed,
          marker='o',
          markersize=1.5,
          linestyle='none')
yfit_decompressed = two_timescale_decompression_drift(xdata_decompressed, *fit_decompressed)
plt.plot(xdata_decompressed, yfit_decompressed,
          linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Hencky strain (-)')
plt.show()

plt.figure(dpi=500)
plt.plot(xdata_compressed, ydata_compressed,
          marker='o',
          markersize=1.5,
          linestyle='none')
plt.plot(xdata_compressed, yfit_compressed,
         linewidth=1)
plt.plot(xdata_decompressed+300, ydata_decompressed,
          marker='o',
          markersize=1.5,
          linestyle='none')
plt.plot(xdata_decompressed+300, yfit_decompressed,
          linewidth=1)
plt.xlim([-5, 600])
plt.ylim([0, 0.25])
plt.xlabel('Time (s)')
plt.ylabel('Hencky strain (-)')
plt.show()








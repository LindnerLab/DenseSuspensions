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

tracked = pd.read_pickle(r'E:\Lars\Step Stress\20201130 Compression Decompression\Stress175_35_25Back\Preprocessed\V1\tracked2.pkl')
nParticles = len(tracked.particle.unique())
nFrames = len(tracked.frame.unique())

avg_strain = np.empty((nFrames,1))
for i in range(nFrames):
    avg_strain[i] = np.mean(tracked[tracked.frame == i].Hencky_strain)


start_compression = [29, 2987, 7679, 10038]
start_decompression = [1223, 4456, 8861]
plt.figure(dpi=500)
for i in range(3):
    plt.plot(avg_strain[start_compression[i]:start_decompression[i]-1])
plt.xlabel('Time (s)')
plt.ylabel('$\langle \epsilon_{Hencky} \\rangle \\ (-)$')
plt.show()

plt.figure(dpi=500)
for i in range(3):
    plt.plot(avg_strain[start_decompression[i]:start_compression[i+1]-1])
plt.xlabel('Time (s)')
plt.ylabel('$\langle \epsilon_{Hencky} \\rangle \\ (-)$')
plt.show()
    
ydata_compressed = np.array(avg_strain[start_compression[1]:start_decompression[1]-1]).reshape((1468,))
xdata_compressed = np.arange(0,len(ydata_compressed)/10,0.1)
fit_compressed, _ = opt.curve_fit(two_timescale_compression_drift,
                                  xdata_compressed,
                                  ydata_compressed,
                                  p0=[0.12, 0.08, 0.08, 0.0140, 1E-4],
                                  maxfev=10000)
plt.figure(dpi=500)
plt.plot(xdata_compressed, ydata_compressed,
          marker='o',
          markersize=1.5,
          linestyle='none')
yfit_compressed = two_timescale_compression_drift(xdata_compressed, *fit_compressed)
plt.plot(xdata_compressed, yfit_compressed,
          linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('$\langle \epsilon_{Hencky} \\rangle \\ (-)$')
plt.show()

ydata_decompressed = np.array(avg_strain[start_decompression[1]:start_compression[2]-1]).reshape((3222,))
xdata_decompressed = np.arange(0,len(ydata_decompressed)/10,0.1)
fit_decompressed, _ = opt.curve_fit(two_timescale_decompression_drift,
                                  xdata_decompressed,
                                  ydata_decompressed,
                                  p0=[0.12, 0.08, 0.08, 0.0140, 1E-4, 0.02],
                                  maxfev=10000)
plt.figure(dpi=500)
plt.plot(xdata_decompressed, ydata_decompressed,
          marker='o',
          markersize=1.5,
          linestyle='none')
yfit_decompressed = two_timescale_decompression_drift(xdata_decompressed, *fit_decompressed)
plt.plot(xdata_decompressed, yfit_decompressed,
          linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('$\langle \epsilon_{Hencky} \\rangle \\ (-)$')
plt.show()

plt.figure(dpi=500)
plt.plot(xdata_compressed, ydata_compressed,
          marker='o',
          markersize=2,
          c='C0',
          linestyle='none')
plt.plot(xdata_decompressed+146.7, ydata_decompressed,
          marker='o',
          markersize=2,
          c='C0',
          linestyle='none')
# plt.plot(xdata_compressed, yfit_compressed,
#           linewidth=1,
#           color='C1')
# plt.plot(xdata_decompressed+146.7, yfit_decompressed,
#           linewidth=1,
#           color='C1')
plt.xlabel('Time (s)')
plt.ylabel('$\langle \epsilon_{Hencky} \\rangle \\ (-)$')
plt.show()
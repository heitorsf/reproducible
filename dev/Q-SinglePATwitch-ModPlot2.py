#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 02:59:28 2017

@author: heitor
"""

import numpy as np
from matplotlib import pyplot as plt

cells_v = np.loadtxt("vsomaS.txt")
forceS = np.loadtxt("forceS.txt")
forcefr = np.loadtxt("forceFR.txt")
forceff = np.loadtxt("forceFF.txt")
t = np.loadtxt("t.txt")
latS = np.loadtxt("lat"+"S"+".txt")
latFR = np.loadtxt("lat"+"FR"+".txt")
latFF = np.loadtxt("lat"+"FF"+".txt")

plt.figure()
xlim=(0,800)
#ylim = (-10, t.max() + 10)
ylim = (-.3,2.5)

plt.subplot(4,1,4)
plt.plot(t, cells_v, 'k')#, label='Potencial de membrana NM tipo S')
plt.xlim(xlim)
#plt.vlines(spkt, 0, ylim, 'k', linestyle='dashed')
#plt.vlines(spktAx, 0, ylim, 'k', linestyle='dashed')
plt.ylabel('Potencial\nde membrana\n [mV]')
plt.xlabel('tempo [ms]')
plt.legend(fontsize='small')

plt.subplot(4,1,1)
plt.plot(t, forceS[:t.size].transpose(), 'k', label="S ")
#plt.vlines(t[forceS==forceS.max()], -0.2, 3, 'k', linestyle='dashed')
#plt.vlines(23.19, -0.2, 3, 'k', linestyle='dashed')
#plt.vlines(spkt[0], -0.2, 3, 'k', linestyle='dashed')
plt.hlines(-.15,25,t[forceS>0][0], 'k', linestyle='solid')
plt.vlines(25, -.1,-.2, 'k')
plt.vlines(t[forceS>0][0], -.1,-.2, 'k')
#plt.hlines(1.8, t[forceS>0][0], t[forceS.argmax()], 'k')
#plt.vlines(t[forceS>0][0], 1.75,1.85, 'k')
#plt.vlines(t[forceS.argmax()], 1.75,1.85, 'k')
plt.xlim(xlim)
plt.ylim(ylim)
plt.ylabel('Abalo\nMuscular [N]')
plt.text(56, -.1, "%.1f ms" %latS)#, bbox=dict(facecolor='white'))
plt.legend(fontsize='small')
#plt.grid()

plt.subplot(4,1,2)
plt.plot(t, forcefr[:t.size].transpose(), 'k', label="FR")
#plt.vlines(t[forcefr==forcefr.max()], -0.2, 3, 'k', linestyle='dashed')
#plt.vlines(22.13, -0.2, 3, 'k', linestyle='dashed')
#plt.vlines(spkt[0], -0.2, 3, 'k', linestyle='dashed')
plt.hlines(-.15,25,t[forcefr>0][0], 'k', linestyle='solid')
plt.vlines(25, -.1,-.2, 'k')
plt.vlines(t[forcefr>0][0], -.1,-.2, 'k')
plt.xlim(xlim)
plt.ylim(ylim)
plt.ylabel('Abalo\nMuscular [N]')
lat = t[forcefr>0][0]-25
plt.text(56, -.1, "%.1f ms" %latFR)#, bbox=dict(facecolor='white'))
plt.legend(fontsize='small')
#plt.grid()

plt.subplot(4,1,3)
plt.plot(t, forceff[:t.size].transpose(), 'k', label="FF")
#plt.vlines(t[forceff==forceff.max()], -0.2, 3, 'k', linestyle='dashed')
#plt.vlines(21.66, -0.2, 3, 'k', linestyle='dashed')
#plt.vlines(spkt[0], -0.2, 3, 'k', linestyle='dashed')
plt.hlines(-.15,25,t[forceff>0][0], 'k', linestyle='solid')
plt.vlines(25, -.1,-.2, 'k')
plt.vlines(t[forceff>0][0], -.1,-.2, 'k')
plt.xlim(xlim)
plt.ylim(ylim)
plt.ylabel('Abalo\nMuscular [N]')
plt.text(56, -.1, "%.1f ms" %latFF)#, bbox=dict(facecolor='white'))
plt.legend(fontsize='small')
#plt.grid()
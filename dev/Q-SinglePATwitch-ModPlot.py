#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 23:25:24 2017

@author: heitor
"""

import numpy as np
from matplotlib import pyplot as plt

vsomaS = np.loadtxt("vsomaS.txt")
vsomaFR = np.loadtxt("vsomaFR.txt")
vsomaFF = np.loadtxt("vsomaFF.txt")
forceS = np.loadtxt("forceS.txt")
forceFR = np.loadtxt("forceFR.txt")
forceFF = np.loadtxt("forceFF.txt")
t = np.loadtxt("t.txt")

ylim = (-.3,2.5)
xlim = (0,800)

plt.subplot(3,1,1)
plt.plot(t, forceS[:t.size].transpose(), 'k', label="NM tipo S ")
plt.plot(t, vsomaS*.018,'b')
#plt.hlines(-.15,25,t[forceS>0][0], 'k', linestyle='solid')
#plt.vlines(25, -.1,-.2, 'k')
#plt.vlines(t[forceS>0][0], -.1,-.2, 'k')
#plt.hlines(1.8, t[forceS>0][0], t[forceS.argmax()], 'k')
#plt.vlines(t[forceS>0][0], 1.75,1.85, 'k')
#plt.vlines(t[forceS.argmax()], 1.75,1.85, 'k')
plt.xlim(xlim)
plt.ylim((-0.3, 2.3))
plt.ylabel('Twitch [N]')
plt.xlabel('tempo [ms]')
#plt.legend(fontsize='small')
#plt.grid()

plt.subplot(3,1,2)
plt.plot(t, forceFR[:t.size].transpose(), 'k', label="NM tipo S ")
plt.plot(t, vsomaFR*.018,'b')
#plt.hlines(-.15,29.5,t[forceFR>0][0], 'k', linestyle='solid')
#plt.vlines(29.5, -.1,-.2, 'k')
#plt.vlines(t[forceFR>0][0], -.1,-.2, 'k')
#plt.hlines(1.8, t[forceFR>0][0], t[forceFR.argmax()], 'k')
#plt.vlines(t[forceFR>0][0], 1.75,1.85, 'k')
#plt.vlines(t[forceFR.argmax()], 1.75,1.85, 'k')
plt.xlim(xlim)
plt.ylim((-0.3, 2.3))
plt.ylabel('Twitch [N]')
plt.xlabel('tempo [ms]')
#plt.legend(fontsize='small')
#plt.grid()

plt.subplot(3,1,3)
plt.plot(t, forceFF[:t.size].transpose(), 'k', label="NM tipo S ")
plt.plot(t, vsomaFF*.018,'b')
#plt.hlines(-.15,34.5,t[forceFF>0][0], 'k', linestyle='solid')
#plt.vlines(34.5, -.1,-.2, 'k')
#plt.vlines(t[forceFF>0][0], -.1,-.2, 'k')
#plt.hlines(1.8, t[forceFF>0][0], t[forceFF.argmax()], 'k')
#plt.vlines(t[forceFF>0][0], 1.75,1.85, 'k')
#plt.vlines(t[forceFF.argmax()], 1.75,1.85, 'k')
plt.xlim(xlim)
plt.ylim((-0.3, 2.3))
plt.ylabel('Twitch [N]')
plt.xlabel('tempo [ms]')
#plt.legend(fontsize='small')
#plt.grid()

'''


#spkt = np.array([[91.1]])
#vel_cond = 45. # m/s
#vel_cond = vel_cond/1000 # m/ms
#axon_len = 0.6 # m
#atrasoAx = axon_len/vel_cond
#spktAx = spkt + atrasoAx
#tAx = t+atrasoAx

plt.figure()
plt.subplot(2,3,4)
#plt.plot(t, cells_v[1:], 'k')
#plt.plot(tAx, cells_v[1:], '--k')
#plt.plot(spkt, spkv, '|', markeredgewidth=1.5, markersize=10)
#plt.plot(spktAx, spkv, '|', markeredgewidth=1.5,markersize=10)
#plt.plot(spkt, 1, '|k', markeredgewidth=1.5, markersize=10)
plt.plot(t, vsomaS,'-k')
xlim = (0,800)
plt.xlim(xlim)
#plt.plot(spktAx, 1, '|', markeredgewidth=1.5,markersize=10)
#ylim = (-10, cells_v[:].max() + 10)
#ylim = (-10, t.max() + 10)
#plt.ylim((0,2))
#plt.xlim(xlim)
#plt.vlines(spkt, 0, ylim, 'k', linestyle='dashed')
#plt.vlines(spktAx, 0, ylim, 'k', linestyle='dashed')
plt.ylabel('Potencial de membrana')

plt.subplot(2,3,5)
plt.plot(t, vsomaFR,'-k')
plt.xlim(xlim)

plt.subplot(2,3,6)
plt.plot(t, vsomaFF,'-k')
plt.xlim(xlim)

plt.subplot(2,3,1)
plt.plot(t, forceS[:t.size].transpose(), 'k', label="NM tipo S ")
plt.plot(t, vsomaS*.018,'0.5')
#plt.vlines(t[forceS==forceS.max()], -0.2, 3, 'k', linestyle='dashed')
#plt.vlines(25, -0.2, 3, 'k', linestyle='dashed')
plt.hlines(-.15,25,t[forceS>0][0], 'k', linestyle='solid')
plt.vlines(25, -.1,-.2, 'k')
plt.vlines(t[forceS>0][0], -.1,-.2, 'k')
plt.hlines(1.8, t[forceS>0][0], t[forceS.argmax()], 'k')
plt.vlines(t[forceS>0][0], 1.75,1.85, 'k')
plt.vlines(t[forceS.argmax()], 1.75,1.85, 'k')
plt.xlim(xlim)
plt.ylim((-0.3, 2.3))
plt.ylabel('Twitch [N]')
plt.xlabel('tempo [ms]')
#plt.legend(fontsize='small')
#plt.grid()

plt.subplot(2,3,2)
plt.plot(t, forceFR[:t.size].transpose(), 'k', label="NM tipo FR ")
plt.vlines(t[forceFR==forceFR.max()], -0.2, 3, 'k', linestyle='dashed')
plt.vlines(27, -0.2, 3, 'k', linestyle='dashed')
#plt.vlines(spkt[0], -0.2, 3, 'k', linestyle='dashed')
plt.xlim(xlim)
plt.ylim((-0.2, 3))
plt.ylabel('Twitch [N]')
plt.xlabel('tempo [ms]')
#plt.legend(fontsize='small')
#plt.grid()

plt.subplot(2,3,3)
plt.plot(t, forceFF[:t.size].transpose(), 'k', label="NM tipo FF ")
plt.vlines(t[forceFF==forceFF.max()], -0.2, 3, 'k', linestyle='dashed')
plt.vlines(33, -0.2, 3, 'k', linestyle='dashed')
#plt.vlines(spkt[0], -0.2, 3, 'k', linestyle='dashed')
plt.xlim(xlim)
plt.ylim((-0.2, 3))
plt.ylabel('Twitch [N]')
plt.xlabel('tempo [ms]')
#plt.legend(fontsize='small')
#plt.grid()


'''
'''
plt.subplot(2,3,3)
plt.plot(t, forcefr[:t.size].transpose(), 'k', label="NM tipo FR")
plt.vlines(t[forcefr==forcefr.max()], -0.2, 3, 'k', linestyle='dashed')
plt.vlines(22.13, -0.2, 3, 'k', linestyle='dashed')
plt.vlines(spkt[0], -0.2, 3, 'k', linestyle='dashed')
plt.xlim(xlim)
plt.ylim((-0.2, 3))
plt.ylabel('Twitch [N]')
plt.xlabel('tempo [ms]')
plt.legend(fontsize='small')
plt.grid()

plt.subplot(4,1,4)
plt.plot(t, forceff[:t.size].transpose(), 'k', label="NM tipo FF")
plt.vlines(t[forceff==forceff.max()], -0.2, 3, 'k', linestyle='dashed')
plt.vlines(21.66, -0.2, 3, 'k', linestyle='dashed')
plt.vlines(spkt[0], -0.2, 3, 'k', linestyle='dashed')
plt.xlim(xlim)
plt.ylim((-0.2, 3))
plt.ylabel('Twitch [N]')
plt.xlabel('tempo [ms]')
plt.legend(fontsize='small')
plt.grid()
'''
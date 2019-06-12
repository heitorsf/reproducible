#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 23:46:19 2017

@author: heitor
"""
import numpy as np

def forceFromSpikes(spkt, t, mutype):
    dt = t[1]
    force = np.zeros(t.size*2-1) # convolution result size
    for cell_index in range(spkt.shape[0]):
        impulse_train = np.zeros(t.size)
        for instant in spkt[cell_index,:]:
            impulse_train[int(instant/dt)] = 1
        
        if mutype == "S" or mutype =="s":
            Fmaxmu = 1.70    # N
            Tc =  110.       # ms
            Thr = 92.       # ms
        elif mutype == "FR" or mutype == "fr":
            Fmaxmu = 1.9    # N
            Tc =  73.5       # ms
            Thr = 60.       # ms
        elif mutype == "FF" or mutype == "ff":
            Fmaxmu = 2.15    # N
            Tc =  55.       # ms
            Thr = 60.       # ms
        
        k = np.log(2.) / (Thr -Tc * np.log(abs(Thr-Tc) / Tc))
        
        m = k*Tc
        
        p = Fmaxmu * np.exp(-k * Tc * (np.log(Tc) - 1.))
        
        a = p * t**m * np.exp(-k*t)
    #    print('Calculou o a')
        Fmu = np.convolve(a,impulse_train)
        
        force = np.vstack((force, Fmu))
    return force[1:,]

#c = 0.05
#
#print(inst_freq[1:].mean())
#plt.plot(t,Fmu[:t.size])
#plt.plot(t,volt_soma_np)
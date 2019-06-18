#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 21:19:36 2017

@author: heitor
"""

import numpy as np
from scipy.signal import lfilter
from matplotlib import pyplot as plt
import time

def forceSOF(spkt, t, mutype):
    start_tudo = time.time()
    start = time.time()
    dt = t[1]
    forces = np.zeros(t.size)
    if isinstance(mutype,str):
        mutype = [mutype]
    for cell_index in range(spkt.shape[0]): # Para cada celula/spike train:        
        if mutype[cell_index] == "S" or mutype[cell_index] =="s":
            Fmaxmu = 1.70    # N
            Tc =  110.       # ms
            Thr = 92.       # ms
            vel_cond = 45.5 # m/s
            ctet = 0.5
            Ftet = 10.1
        elif mutype[cell_index] == "FR" or mutype[cell_index] == "fr":
            Fmaxmu = 1.9    # N
            Tc =  73.5       # ms
            Thr = 60.       # ms
            vel_cond = 49.5 # m/s
            ctet = 0.2
            Ftet = 15.2
        elif mutype[cell_index] == "FF" or mutype[cell_index] == "ff":
            Fmaxmu = 2.15    # N
            Tc =  55.       # ms
            Thr = 60.       # ms
            vel_cond = 51.5 # m/s
            ctet = 0.8
            Ftet = 19.3
        else:
            raise ValueError("MU type error: "+str(mutype[cell_index]))
        
        print mutype[cell_index]
        vel_cond = vel_cond/1000. # m/ms
        axon_len = 0.6 # m
        atrasoAx = axon_len/vel_cond
        
        Fmu = np.zeros(t.size)
#        plt.figure()
        start = time.time()
        for instant in spkt[cell_index,:]: # Para cada spike
            spike = np.zeros(t.size)
            if instant > 0:
                try:
                    spike[int((instant+atrasoAx)/dt)] = 1/dt
                except ValueError as err:
                    print "Entrou no ValueError:", err.message
                    if err.message!="cannot convert float Nan to integer":
                        raise ValueError(err.message)
                    continue
        
            start = time.time()
            a = [1, -2*np.exp(-(dt/Tc)), np.exp(-2*(dt/Tc))]
            b = [0, Fmaxmu*((dt**2)/Tc)*np.exp(1-dt/Tc)]
        
            start = time.time()
            F = lfilter(b, a, spike)
#            plt.plot(t,F)
            Fmu += F
#        favg = Fmu[int(1400/s):int(1500/s)].mean()
#        Fsat = (1-np.exp(-ctet*Fmu))/(1+np.exp(-ctet*Fmu))
        Fmu = np.clip(Fmu, 0, Ftet)
        forces = np.vstack((forces, Fmu))
    print "Tempo da funcao QSOFilter:", (time.time()-start_tudo)
    return forces[1:,]
#c = 0.05
#
#print(inst_freq[1:].mean())
#plt.plot(t,Fmu[:t.size])
#plt.plot(t,volt_soma_np)


'''
Implementacao original

t = np.arange(0,5000,0.025)

s = 0.025
t_peak = 110.0
A_peak = 1.7

a = [1, -2*np.exp(-(s/t_peak)), np.exp(-2*(s/t_peak))]
b = [0, A_peak*((s**2)/t_peak)*np.exp(1-(s/t_peak))]
force = np.array([])

for freq in np.arange(20,60,1):    
    if freq > 40:
        freq = 40 + 5 - np.exp(-(freq-45))
        
    spkrate = 1000./freq
    
    spkt = np.arange(100,2000,spkrate)/s
    Fmu = np.zeros(np.size(t))
    for i in spkt:
        spike = np.zeros(len(t))
        spike[int(i)] = 1/s
        
        F = lfilter(b, a, spike)
    #    plt.plot(t,F,'0.5')
        Fmu += F
    favg = Fmu[int(1400/s):int(1500/s)].mean()
    force = np.append(force, favg)
    print freq, favg

#    plt.plot(t,Fmu, 'k')
#    plt.xlabel('tempo [ms]')
#    plt.ylabel('Forca [N]')
#    plt.grid()
freq = np.arange(20,60,1)
plt.figure()
plt.plot(freq, force)
'''

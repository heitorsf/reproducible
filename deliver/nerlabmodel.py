#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt; plt.ion()
from scipy.signal import lfilter
from itertools import chain

def getSpikes(t, v, setdvdt, output=False, plotRaster=True, newfigure=True):
    '''
     Expected matrices shapes:
     t: [length] or [length,] or [length,1]
     v: [length, numberOfCells]
     where length is the size of the time vector
    '''
    if not isinstance(t, np.ndarray):
        t = np.array(t)
    if not isinstance(v, np.ndarray):
        v = np.array(v)
    if len(v.shape)<2:
        v.shape = (v.shape[0],1)
    dvdt = np.diff(v, axis=0)
    if t.shape[0]>v.shape[0]:
        t = t[:-1]
    dt = t[1]
    
    ''' Detect Spikes '''
    # creates a boolean mask to find every dvdt element that is >= threshold(setdvdt)
    hasPA = dvdt>=setdvdt*dt
    numCells = len(hasPA[0,:])
    spkt_list = list()
    for cell_id in range(numCells):
        spkt_list_cell = list()
        numInstants = len(hasPA[:,cell_id])
        for instant in range(numInstants):
            if hasPA[instant,cell_id] and not(hasPA[instant-1,cell_id]):
                if v[instant, cell_id] > 5:
                    spkt_list_cell.append(t[instant])
        spkt_list.append(np.array(spkt_list_cell))
    
    ''' Create spiketrain times (spkt) and volts (spkv) matrices '''
    # Set matrices size according to the cell with more spikes detected    
    numPAs = np.array([])
    for i in range(numCells):
        numPAs = np.append(numPAs, len(spkt_list[i]))
    maxPAs = int(numPAs.max())
    
    # Create and fill the matrices
    spkt = np.zeros((maxPAs, numCells))*np.nan
    spkv = np.zeros(spkt.shape)*np.nan
    for i in range(numCells):
        for j in range(len(spkt_list[i])):
            spkt[j,i] = spkt_list[i][j]
            spkv[j,i] = v[:,i][t==spkt_list[i][j]]
#            spkv[j,i] = v[:,i][t==round(float(spkt_list[i][j]), 2)]


    ''' Plot Raster '''
    if plotRaster and newfigure:
        plt.figure()
        plt.xlabel("tempo [ms]")
        plt.ylabel(u'Número do NM')
    if plotRaster:
        if len(spkt[:,:]) > 0:
            #print "Plotting Spikes..."
            for i in range(len(spkt[0,:])):
                plt.plot(spkt[:,i], np.ones(len(spkt[:,i]))*i, '|', markeredgewidth=1.5)
                
            plt.xlim(left=0)
            plt.show()
    if output:
        return spkt.transpose(), spkv.transpose()

#_______________________________________________________________________

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

#_______________________________________________________________________

def forceSOF(spkt, t, mutype, plotForces=False):
    dt = t[1]
    forces = np.zeros(t.size)
    if isinstance(mutype,str):
        mutype = [mutype]
    if plotForces:
        plt.figure()
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
        
        #print mutype[cell_index]
        vel_cond = vel_cond/1000. # m/ms
        axon_len = 0.6 # m
        atrasoAx = axon_len/vel_cond
        
        Fmu = np.zeros(t.size)
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
        
            a = [1, -2*np.exp(-(dt/Tc)), np.exp(-2*(dt/Tc))]
            b = [0, Fmaxmu*((dt**2)/Tc)*np.exp(1-dt/Tc)]
            F = lfilter(b, a, spike)
            Fmu += F
#        favg = Fmu[int(1400/s):int(1500/s)].mean()
#        Fsat = (1-np.exp(-ctet*Fmu))/(1+np.exp(-ctet*Fmu))
        Fmu = np.clip(Fmu, 0, Ftet)
        if plotForces:
            plt.xlabel("tempo [ms]")
            plt.ylabel(u'Força muscular [N]')
            plt.plot(t,Fmu,label=mutype[cell_index])
        forces = np.vstack((forces, Fmu))
        
    if plotForces:
        plt.plot(t,forces.sum(axis=0),'k',label=u'Força total')
        plt.legend()
    return forces[1:,]

#_______________________________________________________________________

def mus(soma,dend):
    # soma parameters
    soma.L = 80
    soma.nseg = 1
    soma.diam = 80
    soma.Ra = 70.0 #citoplasmatic resistivity
    soma.cm = 1.0
    soma.gl_napp = 1/1100.0

    soma.gnabar_napp = 0.05
    soma.gnapbar_napp = .00052
    soma.gkfbar_napp = 0.0028
    soma.gksbar_napp = 0.018
    soma.mact_napp = 13.0 # alpha = .64*vtrap(mact-v2,4) (Na activation)
    soma.rinact_napp = 0.025 # beta = rinact (slow K inactivation)

    soma.ena = 120.0
    soma.ek = -10.0
    soma.el_napp = 0.0
    soma.vtraub_napp = 0.0

    # dendrite parameters
    dend.L = 6150.0
    dend.nseg = 1
    dend.diam = 52
    dend.Ra = 70.0
    dend.cm = 1.0

    for seg in chain(dend):
        seg.pas.g = 1/12550.0
        seg.pas.e = 0.0

def mufr(soma,dend):
    # soma parameters
    soma.L = 85
    soma.nseg = 1
    soma.diam = 85
    soma.Ra = 70.0
    soma.cm = 1.0

    soma.gnabar_napp = 0.07
    soma.gnapbar_napp = 0.0008
    soma.gkfbar_napp = 0.0040
    soma.gksbar_napp = 0.037
    soma.mact_napp= 17.0
    soma.rinact_napp = 0.058
    soma.ena = 120.0
    soma.ek = -10.0
    soma.el_napp = 0.0
    soma.vtraub_napp = 0.0
    soma.gl_napp = 1/1000.0

    # dendrite parameters
    dend.L = 7450.0
    dend.nseg = 1
    dend.diam = 73
    dend.Ra = 70.0
    dend.cm = 1.0

    for seg in chain(dend):
        seg.pas.g = 1/8825.0
        seg.pas.e = 0.0

def muff(soma,dend):
    # soma parameters
    soma.L = 100.25
    soma.nseg = 1
    soma.diam = 100.25
    soma.Ra = 70.0
    soma.cm = 1.0

    soma.gnabar_napp = 0.075
    soma.gnapbar_napp = 0.00065
    soma.gkfbar_napp = 0.00135
    soma.gksbar_napp = 0.016
    soma.mact_napp = 19.2
    soma.rinact_napp = 0.062
    soma.ena = 120.0
    soma.ek = -10.0
    soma.el_napp = 0.0
    soma.vtraub_napp = 0.0
    soma.gl_napp = 1/800.0

    # dendrite parameters
    dend.L = 9350.0
    dend.nseg = 1
    dend.diam = 88
    dend.Ra = 70.0
    dend.cm = 1.0

    for seg in chain(dend):
        seg.pas.g = 1/6500.0
        seg.pas.e = 0.0
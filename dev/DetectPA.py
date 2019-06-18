#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 16:44:58 2017

@author: heitor
"""

import numpy as np
from matplotlib import pyplot as plt; plt.ion()

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
    if plotRaster:
        if len(spkt[:,:]) > 0:
            print "Plotting Raster..."
            for i in range(len(spkt[0,:])):
                plt.plot(spkt[:,i], np.ones(len(spkt[:,i]))*i, '|', markeredgewidth=1.5)
            plt.show()
    if output:
        return spkt.transpose(), spkv.transpose()

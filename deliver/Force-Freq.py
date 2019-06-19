#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 01:24:45 2017

@author: heitor
"""

from netpyne import sim, specs
import numpy as np
from matplotlib import pyplot as plt
import nerlabmodel as ner
import time

times=[]
freqs = np.array([])
forces = np.array([])
cells_vs = np.array([])
summa = np.array([])
#plt.figure()
import sys
mutype = sys.argv[1]
#mutype = 'S'
if mutype == "S" or mutype =="s":
    i_range = np.logspace(np.log10(6),np.log10(25), 15) # S
elif mutype == "FR" or mutype == "fr":
    i_range = np.logspace(np.log10(13),np.log10(56), 15) # FR
elif mutype == "FF" or mutype == "ff":
    i_range = np.logspace(np.log10(22),np.log10(58), 15) # FF
    
for i_stim in i_range:
    start = time.time()
    netParams = specs.NetParams()
    simConfig = specs.SimConfig()
    
    netParams.popParams['Pop1'] = {'cellType': 'MED',
                                   'cellModel': 'NAPP',
                                   'numCells': 1}
    cellRule = {'conds': {'cellType': 'MED', 'cellModel': 'NAPP'},
                'secs': {'soma': {'geom': {'diam': 80, 'L': 80, 'Ra': 70,
                                           'cm': 1, 'nseg': 1},
                                  'mechs': {'napp': {'gnabar': 0.05,
                                                     'gnapbar': 0.00052,
                                                     'gkfbar': 0.0028,
                                                     'gksbar': 0.018,
                                                     'mact': 13.0,
                                                     'rinact': 0.025,
                                                     'el': 0.0,
                                                     'vtraub': 0.0,
                                                     'gl': 1/1100.0}},
                                  'ions': {'na': {'e': 120},
                                           'k': {'e': -10}},
                                  'vinit': 0},
                         'dend': {'geom': {'diam': 52.0, 'L': 6150.0, 'Ra': 70.0,
                                           'cm': 1, 'nseg': 1},
                                  'topol': {'parentSec': 'soma',
                                            'parentX': 0,
                                            'childX': 1},
                                  'mechs': {'pas': {'g': 1/12550.0, 'e': 0.0}},
    #                              'ions': {'caL': {'e': 140}},
                                  'vinit': 0}}}
    netParams.cellParams['HHNapp'] = cellRule
    
    netParams.stimSourceParams['Istep'] = {'type': 'IClamp_X',
                                           'sig_type': 1,
                                           'delay': 20,
                                           'dur': 2000,
                                           'amp': i_stim
                                           }
    netParams.stimTargetParams['Istep->Pop1'] = {'source': 'Istep',
                                                 'conds': {'cellType': 'MED',
                                                           'cellModel': 'NAPP'},
                                                 'sec': 'soma',
                                                 'loc': 0.5}
    
            
    simConfig.duration = 4*1e3
    simConfig.dt = 0.05
    simConfig.seeds = {'conn': 1, 'stim': 1, 'loc': 1}
    simConfig.createNEURONObj = True
    simConfig.createPyStruct = True
    simConfig.verbose = False
    
    
    simConfig.recordCells = ['all']
    simConfig.recordTraces = {'v_soma': {'sec':'soma', 'loc': 0.5, 'var': 'v'}}
    simConfig.recordStim = True
    simConfig.recordStep = simConfig.dt
    
    simConfig.filename = 'Netpyne-cell_output'
    simConfig.saveTxt = True
    
    #simConfig.analysis['plotRaster'] = True
    #simConfig.analysis['plotTraces'] = {'include': [0,1,2], 'oneFigPer': 'trace'}
    #simConfig.analysis['plot2Dnet'] = True  # Plot 2D net cells and connections
    
    data = sim.create(netParams=netParams, simConfig=simConfig, output=True)
    soma = data[1][0].secs.soma.hSec
    dend = data[1][0].secs.dend.hSec
#    mutype = 's'
    if mutype == "S" or mutype =="s":
        ner.mus(soma, dend)
    elif mutype == "FR" or mutype == "fr":
        ner.mufr(soma, dend)
    elif mutype == "FF" or mutype == "ff":
        ner.muff(soma, dend)
        
    sim.simulate()
    sim.analyze()
    
    t = np.arange(0, simConfig.duration+simConfig.dt, simConfig.dt)
    
    v_list = list()
    for i in data[0]['Pop1'].cellGids:
        cell = 'cell_'+str(i)
        v_list.append(np.array(data[4]['v_soma'][cell].to_python()))
    cells_v = np.array(v_list).transpose()
    
    spkt, spkv = ner.getSpikes(t, cells_v, 20, output=True, plotRaster=False, newfigure=False)

    force = ner.forceSOF(spkt, t, mutype)
    
    vel_cond = 45. # m/s
    axon_len = 0.6 # m
    atrasoAx = axon_len/vel_cond
    spktAx = spkt+ atrasoAx
    
    freq = 1000./np.diff(spkt[:,1:]).mean()
    max_force = force[:,int(1000/t[1]):int(1500/t[1])].mean()
    freqs = np.append(freqs, freq)
    forces = np.append(forces, max_force)
    cells_vs = np.append(cells_vs, cells_v)
    summa = np.append(summa, force)
    print "force.shape", force.shape
#    plt.plot(t, force)
    print i_stim, freq#, max_force
    if freq>100:
        break
    
    end = time.time()
    times.append([i_stim, end-start])

print times
np.savetxt("freq"+mutype+".txt",freqs)
np.savetxt("force"+mutype+".txt",forces)
print("Saved "+mutype+" data.")


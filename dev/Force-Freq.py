#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 01:24:45 2017

@author: heitor
"""

import nerlabmodel as ner
from netpyne import sim, specs
import numpy as np
from matplotlib import pyplot as plt
from DetectPA import getSpikes
from mu_type import *
from QForceFunc import forceFromSpikes
from QSOFilter import forceSOF
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
        mus(soma, dend)
    elif mutype == "FR" or mutype == "fr":
        mufr(soma, dend)
    elif mutype == "FF" or mutype == "ff":
        muff(soma, dend)
        
    #data[1][0].secs.soma.hSec.__getattribute__('ena')
    #data[1][0].secs.soma.hSec.__setattr__('ena', 120)
    sim.simulate()
    sim.analyze()
    #sim.createSimulateAnalyze(netParams=netParams, simConfig=simConfig)
    
    #v_soma_cell_0 = np.array(data[4]['v_soma']['cell_0'].to_python())
    
    t = np.arange(0, simConfig.duration+simConfig.dt, simConfig.dt)
#    t = np.array(data[4].t.to_python())

    #plt.plot(t, v_soma_cell_0)
    
    v_list = list()
    for i in data[0]['Pop1'].cellGids:
        cell = 'cell_'+str(i)
        v_list.append(np.array(data[4]['v_soma'][cell].to_python()))
    cells_v = np.array(v_list).transpose()
    
    spkt, spkv = getSpikes(t, cells_v, 20, output=True, plotRaster=False, newfigure=False)
#    spkt,spkv = ner.getSpikes(t, cells_v,output=True)
    
#    force = np.array([1,2])
#    force = forceFromSpikes(spkt, t, mutype)
    force = forceSOF(spkt, t, mutype)
#    plt.figure()
#    plt.plot(t, np.sum(force, axis=0).transpose()[:t.size])
#    plt.plot(t, force[:,:t.size].transpose())
    
    vel_cond = 45. # m/s
    axon_len = 0.6 # m
    atrasoAx = axon_len/vel_cond
    spktAx = spkt+ atrasoAx
#
#    plt.figure()
#    plt.subplot(2,1,1)
#    #plt.plot(t, cells_v[:,0], 'k', label='AMembrane Potential [mV]')
#    plt.plot(spkt, 1, '|k', markeredgewidth=1.2, markersize=10)
#    #plt.plot(spktAx, spkv, '|', markeredgewidth=1.5)
#    xlim = (-10, 8000)
#    ylim = (-10, cells_v[:,0].max() + 10)
#    plt.ylim((-1,3))
#    plt.xlim(xlim)
#    #plt.vlines(spkt, 0, ylim, 'k', linestyle='dashed')
#    #plt.vlines(spktAx, 0, ylim, 'k', linestyle='dashed')
#    plt.ylabel('Spike Train')
#    plt.subplot(2,1,2)
#    plt.plot(t, force[:,:t.size].transpose(), 'k')
#    plt.xlim(xlim)
#    plt.ylabel('Abalo da Unidade Muscular [N]')
#    plt.xlabel('tempo [ms]')
#    plt.legend(fontsize='small')
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
#summa = summa.reshape((force.shape[1], 11))
#plt.plot(freqs, forces)
np.savetxt("freq"+mutype+".txt",freqs)
np.savetxt("force"+mutype+".txt",forces)
print("Saved "+mutype+" data.")
##np.savetxt("FxFcells_vs_FR.txt", cells_vs)
#np.savetxt("FxF-2forcesmax_"+mutype+".txt", forces)
#np.savetxt("FxF-2freqs_"+mutype+".txt", freqs)
'''
import numpy as np
from matplotlib import pyplot as plt

forcesS = 0.1*np.loadtxt("FxF-2forcesmax_S.txt")
forcesFR = 0.1*np.loadtxt("FxF-2forcesmax_FR.txt")
forcesFF = 0.1*np.loadtxt("FxF-2forcesmax_FF.txt")

freqsS = np.loadtxt("FxF-2freqs_S.txt")
freqsFR = np.loadtxt("FxF-2freqs_FR.txt")
freqsFF = np.loadtxt("FxF-2freqs_FF.txt")

matplotlib.rcParams.update({'font.size': 20})
plt.plot(freqsS, forcesS, 'b', linewidth=1.5, label="Tipo S")
plt.plot(freqsFF, forcesFF, 'g', linewidth=1.5, label="Tipo FF")
plt.plot(freqsFR, forcesFR, 'm', linewidth=1.5, label="Tipo FR")
#plt.plot(freqsS, forcesS, 'b')
#plt.plot(freqsFR, forcesFR, 'g')
#plt.plot(freqsFF, forcesFF, 'm')
plt.legend(fontsize='18', loc = 'upper left')
xlim = (0, 100)
ylim = (0, 2.3)
plt.ylim(ylim)
plt.xlim(xlim)
#plt.vlines(spkt, 0, ylim, 'k', linestyle='dashed')
plt.ylabel(u'For\u00e7a [kgf]')
plt.xlabel(u'Frequ\u00eancia de disparo [Hz]')
plt.grid()
'''

#    # Simulations data point:
#plt.plot(hist_times, hist_inst_freq, 'k.',
#         label='%.1f Hz \nr-square: %.5f' %(f_hist, rsquare))
#    # Fitted sine curve (one cycle, with time normalized to sim. duration):
#plt.plot(time_np[time_np<=period] / time_np[time_np>period][0],
#         e_seno[time_np<=period] , 'r', label='Fit')
#plt.xlabel('Normalized with respect to one input cycle')
#plt.ylabel('Instantaneous Frequency [Hz]')
#plt.legend(fontsize='small')


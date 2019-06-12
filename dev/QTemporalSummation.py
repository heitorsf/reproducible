#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 15:48:51 2017

@author: heitor
"""

from netpyne import sim, specs
import numpy as np
from matplotlib import pyplot as plt
from DetectPA import getSpikes
from mu_type import *
from QForceFunc import forceFromSpikes
from QSOFilter import forceSOF

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
                                       'dur': 5000,
                                       'amp': 9
                                       }
netParams.stimTargetParams['Istep->Pop1'] = {'source': 'Istep',
                                             'conds': {'cellType': 'MED',
                                                       'cellModel': 'NAPP'},
                                             'sec': 'soma',
                                             'loc': 0.5}

        
simConfig.duration = 7.2*1e3
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
mutype = 'S'
if mutype == "S" or mutype =="s":
    mus(soma, dend)
elif mutype == "FR" or mutype == "fr":
    mufr(soma, dend)
elif mutype == "FF" or mutype == "ff":
    muff(soma, dend)

#for i in range(len(data[1])):
#    data[1][i].secs.soma.hSec.L =  75 + 5*i
#    data[1][i].secs.soma.hSec.diam = 75 + 5*i
    
#data[1][0].secs.soma.hSec.__getattribute__('ena')
#data[1][0].secs.soma.hSec.__setattr__('ena', 120)
sim.simulate()
sim.analyze()
#sim.createSimulateAnalyze(netParams=netParams, simConfig=simConfig)

#v_soma_cell_0 = np.array(data[4]['v_soma']['cell_0'].to_python())

t = np.arange(0, simConfig.duration, simConfig.dt)

#plt.plot(t, v_soma_cell_0)

v_list = list()
for i in data[0]['Pop1'].cellGids:
    cell = 'cell_'+str(i)
    v_list.append(np.array(data[4]['v_soma'][cell].to_python()))
cells_v = np.array(v_list).transpose()

spkt, spkv = getSpikes(t, cells_v, 20, output=True, plotRaster=False, newfigure=False)

force = forceSOF(spkt, t, mutype)

vel_cond = 45. # m/s
axon_len = 0.6 # m
atrasoAx = axon_len/vel_cond
spktAx = spkt+ atrasoAx

plt.subplot(2,1,1)
#plt.plot(t, cells_v[:,0], 'k', label='AMembrane Potential [mV]')
plt.plot(spkt, 2, '|k', markeredgewidth=1.2, markersize=10)
#plt.plot(spktAx, spkv, '|', markeredgewidth=1.5)
xlim = (-10, 8000)
ylim = (-10, cells_v[:,0].max() + 10)
plt.ylim((-1,3))
plt.xlim(xlim)
#plt.vlines(spkt, 0, ylim, 'k', linestyle='dashed')
#plt.vlines(spktAx, 0, ylim, 'k', linestyle='dashed')
plt.ylabel('Spike Train')
plt.subplot(2,1,2)
plt.plot(t, force[:,:t.size].transpose(), 'k')
plt.xlim(xlim)
plt.ylabel('Forca [N]')
plt.xlabel('tempo [ms]')
plt.legend(fontsize='small')

''' Plot Tips ''''''
plt.plot(hist_times, hist_inst_freq, 'k.',
         label='%.1f Hz \nr-square: %.5f' %(f_hist, rsquare))
plt.plot(time_np[time_np<=period] / time_np[time_np>period][0],
         e_seno[time_np<=period] , 'r', label='Fit')
'''

#np.savetxt('force1', force)
#np.savetxt('spkt1', spkt)
#np.savetxt('t-tetano', t)

#np.savetxt('force2', force)
#np.savetxt('spkt2', spkt)

''' Script pra plotar bonitinho ''''''

import numpy as np
from matplotlib import pyplot as plt

force1 = np.loadtxt("force1")
force2 = np.loadtxt("force2")
spkt1 = np.loadtxt("spkt1")
spkt2 = np.loadtxt("spkt2")
t = np.loadtxt("t-tetano")
freq1 = 1000./np.diff(spkt1[1:]).mean()
freq2 = 1000./np.diff(spkt2[1:]).mean()

plt.figure()
plt.subplot(2,1,1)
plt.plot(spkt1, -1.5*np.ones(spkt1.size), '|k', markeredgewidth=1.5, markersize = 10)
plt.plot(t, force1[:t.size].transpose(), 'k')
plt.text(4000, 3, "%.1f imp/s" %freq1, bbox=dict(facecolor='white'))
plt.legend(fontsize='small', markerscale=0.1)
xlim = (-10, 6000)
ylim = (-3, 11)
plt.ylim(ylim)
plt.xlim(xlim)
#plt.vlines(spkt, 0, ylim, 'k', linestyle='dashed')
plt.hlines(0,-10, 7000, 'grey')
plt.ylabel('Forca [N]')
plt.xlabel('tempo [ms]')
plt.grid()

plt.subplot(2,1,2)
plt.plot(spkt2, -1.5*np.ones(spkt2.size), '|k', markeredgewidth=1.5, markersize = 10)
plt.plot(t, force2[:t.size].transpose(), 'k')
plt.text(4000,3, "%.1f imp/s" %freq2, bbox=dict(facecolor='white'))
plt.legend(fontsize='small')
plt.ylim(ylim)
plt.xlim(xlim)
#plt.vlines(spkt, 0, ylim, 'k', linestyle='dashed')
plt.hlines(0,-10, 7000, 'grey')
plt.ylabel('Forca [N]')
plt.xlabel('tempo [ms]')
plt.grid()
'''
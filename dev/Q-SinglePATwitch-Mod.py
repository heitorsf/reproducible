#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 23:22:46 2017

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
                                       'delay': 15,
                                       'dur': 50,
                                       'amp': 22
                                       }
netParams.stimTargetParams['Istep->Pop1'] = {'source': 'Istep',
                                             'conds': {'cellType': 'MED',
                                                       'cellModel': 'NAPP'},
                                             'sec': 'soma',
                                             'loc': 0.5}

        
simConfig.duration = 4.*1e3
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
mutype = "FF"
if mutype == "S" or mutype =="s":
    mus(soma, dend)
    vel_cond = 45.5 # m/s
elif mutype == "FR" or mutype == "fr":
    mufr(soma, dend)
    vel_cond = 49.5 # m/s
elif mutype == "FF" or mutype == "ff":
    muff(soma, dend)
    vel_cond = 51.5 # m/s
#soma = data[1][1].secs.soma.hSec
#dend = data[1][1].secs.dend.hSec
#muff(soma, dend)
#soma = data[1][2].secs.soma.hSec
#dend = data[1][2].secs.dend.hSec
#mufr(soma, dend)


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
print spkt

vel_cond = vel_cond/1000 # m/ms
axon_len = 0.6 # m
atrasoAx = axon_len/vel_cond
spktAx = spkt + atrasoAx
tAx = t+atrasoAx
'''
plt.figure()
plt.plot(t, cells_v[1:], 'k', label="Vm (potencial de membrana)")
plt.plot(tAx, cells_v[1:], '--k', label="Vpm (potencial na placa motora)")
#plt.plot(spkt, spkv, '|k', markeredgewidth=1.5, markersize=20)
#plt.plot(spktAx, spkv, '|k', markeredgewidth=1.5, markersize=20)
plt.xlim((0,155))
plt.xlabel('tempo [ms]')
plt.ylabel('Potencial [mV]')
plt.legend(fontsize="small")
'''
spkt=np.array([[25.]])
'''force = forceFromSpikes(spkt, t, mutype)'''
force, lat = forceSOF(spkt, t, [mutype])
#plt.figure()
plt.plot(t,force[0], 'k')

#plt.figure()
#plt.plot(t, np.sum(force, axis=0).transpose()[:t.size])
#plt.plot(t, force[:,:t.size].transpose())

'''
plt.subplot(2,1,1)
plt.plot(t, cells_v[:,0], 'k')
#plt.plot(spkt, spkv, '|', markeredgewidth=1.5)
#plt.plot(spktAx, spkv, '|', markeredgewidth=1.5)
xlim = (-10, 2000)
ylim = (-10, cells_v[:,0].max() + 10)
plt.ylim(ylim)
plt.xlim(xlim)
#plt.vlines(spkt, 0, ylim, 'k', linestyle='dashed')
#plt.vlines(spktAx, 0, ylim, 'k', linestyle='dashed')
plt.ylabel('Vm [mV]')

plt.subplot(2,1,2)
plt.plot(t, force[:,:t.size].transpose(), 'k')
plt.xlim(xlim)
#plt.ylim((-0.2, 3))
plt.ylabel('Twitch [N]')
plt.xlabel('tempo [ms]')
plt.legend(fontsize='small')
plt.grid()
'''
''' Plot Tips ''''''
plt.plot(hist_times, hist_inst_freq, 'k.',
         label='%.1f Hz \nr-square: %.5f' %(f_hist, rsquare))
plt.plot(time_np[time_np<=period] / time_np[time_np>period][0],
         e_seno[time_np<=period] , 'r', label='Fit')
'''


''' Script pra plotar bonitinho '''
np.savetxt("force"+mutype+".txt", force)
#np.savetxt("t.txt", t)
#np.savetxt("vsoma"+mutype+".txt", cells_v[1:])
#np.savetxt("lat"+mutype+".txt", np.array([lat]))

'''
import numpy as np
from matplotlib import pyplot as plt

cells_v = np.loadtxt("vsomaS.txt")
forceS = np.loadtxt("forceS.txt")
forcefr = np.loadtxt("forceFR.txt")
forceff = np.loadtxt("forceFF.txt")
t = np.loadtxt("t.txt")

spkt = np.array([[91.1]])
#vel_cond = 45. # m/s
#vel_cond = vel_cond/1000 # m/ms
#axon_len = 0.6 # m
#atrasoAx = axon_len/vel_cond
#spktAx = spkt + atrasoAx
#tAx = t+atrasoAx

plt.figure()
xlim=(0,800)
plt.subplot(4,1,4)
plt.plot(t, cells_v, 'k')
#plt.plot(tAx, cells_v[1:], '--k')
#plt.plot(spkt, spkv, '|', markeredgewidth=1.5, markersize=10)
#plt.plot(spktAx, spkv, '|', markeredgewidth=1.5,markersize=10)
###plt.plot(spkt, 1, '|k', markeredgewidth=1.5, markersize=10)
#plt.plot(spktAx, 1, '|', markeredgewidth=1.5,markersize=10)
xlim = (-10, 1000)
#ylim = (-10, cells_v[:].max() + 10)
ylim = (-10, t.max() + 10)
#plt.ylim((0,2))
plt.xlim(xlim)
#plt.vlines(spkt, 0, ylim, 'k', linestyle='dashed')
#plt.vlines(spktAx, 0, ylim, 'k', linestyle='dashed')
plt.ylabel('Spike')

plt.subplot(4,1,1)
plt.plot(t, forceS[:t.size].transpose(), 'k', label="NM tipo S ")
plt.vlines(t[forceS==forceS.max()], -0.2, 3, 'k', linestyle='dashed')
plt.vlines(23.19, -0.2, 3, 'k', linestyle='dashed')
plt.vlines(spkt[0], -0.2, 3, 'k', linestyle='dashed')
plt.xlim(xlim)
plt.ylim((-0.2, 3))
plt.ylabel('Twitch [N]')
plt.xlabel('tempo [ms]')
plt.legend(fontsize='small')
plt.grid()

plt.subplot(4,1,2)
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

plt.subplot(4,1,3)
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
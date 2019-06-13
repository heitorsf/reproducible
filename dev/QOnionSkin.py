#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 19:59:08 2017

@author: heitor
"""
from netpyne import sim, specs
import numpy as np
from matplotlib import pyplot as plt
import nerlabmodel as ner

netParams = specs.NetParams()
simConfig = specs.SimConfig()

netParams.popParams['Pop1'] = {'cellType': 'MED',
                               'cellModel': 'NAPP',
                               'numCells': 3}
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

netParams.stimSourceParams['Istep'] = {'type': 'IClampX',
                                       'signal_type': 2, #triang
                                       'delay': 20,
                                       'risedur': 250,
                                       'falldur': 250,
                                       'amp': 22
                                       }
netParams.stimTargetParams['Istep->Pop1'] = {'source': 'Istep',
                                             'conds': {'cellType': 'MED',
                                                       'cellModel': 'NAPP'},
                                             'sec': 'soma',
                                             'loc': 0.5}

        
simConfig.duration = 720 #*1e3
simConfig.dt = 0.025
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

(pops, cells, conns, stims, simData) = sim.create(netParams=netParams, simConfig=simConfig, output=True)
mutype = []
for n in range(len(cells)):
    soma = cells[n].secs.soma.hSec
    dend = cells[n].secs.dend.hSec
    if n == 0:
        ner.mus(soma, dend)
        mutype.append('S')
    elif n == 1:
        ner.mufr(soma, dend)
        mutype.append('FR')
    elif n == 2:
        ner.muff(soma, dend)
        mutype.append('FF')

#for i in range(len(data[1])):
#    data[1][i].secs.soma.hSec.L =  75 + 5*i
#    data[1][i].secs.soma.hSec.diam = 75 + 5*i
    
#data[1][0].secs.soma.hSec.__getattribute__('ena')
#data[1][0].secs.soma.hSec.__setattr__('ena', 120)
sim.simulate()
sim.analyze()
#sim.createSimulateAnalyze(netParams=netParams, simConfig=simConfig)

#v_soma_cell_0 = np.array(data[4]['v_soma']['cell_0'].to_python())

t = np.array(simData.t.to_python())
#    np.arange(0, simConfig.duration, simConfig.dt)

#plt.plot(t, v_soma_cell_0)

v_list = list()
for i in pops['Pop1'].cellGids:
    cell = 'cell_'+str(i)
    v_list.append(np.array(simData['v_soma'][cell].to_python()))
cells_v = np.array(v_list).transpose()

spkt, spkv = ner.getSpikes(t, cells_v, 20, output=True, plotRaster=True, newfigure=True)

force = ner.forceSOF(spkt, t, mutype)

ifreq = 1000./np.diff(spkt)
midtimes = spkt[:,1:] - np.diff(spkt)/2.
midtimes[1,51]=0
midtimes[2,28]=0

#vel_cond = 45. # m/s
#axon_len = 0.6 # m
#atrasoAx = axon_len/vel_cond
#spktAx = spkt+ atrasoAx
'''
plt.subplot(2,1,1)
#plt.plot(t, cells_v[:,0], 'k', label='AMembrane Potential [mV]')
plt.plot(spkt, np.arange(min(spkt.shape)), '|k', markeredgewidth=1.2, markersize=10)
#plt.plot(spktAx, spkv, '|', markeredgewidth=1.5)
xlim = (-10, 8000)
ylim = (-10, cells_v[:,0].max() + 10)
plt.ylim((-1,3))
plt.xlim(xlim)
#plt.vlines(spkt, 0, ylim, 'k', linestyle='dashed')
#plt.vlines(spktAx, 0, ylim, 'k', linestyle='dashed')
plt.ylabel('Spike Train')
plt.subplot(2,1,2)'''
#plt.figure()
fig,ax = plt.subplots()
#plt.plot(t, force[:,:t.size].transpose(), 'k')
#plt.plot(t, force.sum(axis=0))
#plt.plot(t, force.transpose())
#plt.plot(spkt[:,1:-1].transpose(),ifreq[:,:-1].transpose(),'.')
for cell in range(min(spkt.shape)):
    if cell == 0:
        ax.plot(spkt[cell,:], -6*np.ones(max(spkt.shape)), '|b',
                markeredgewidth=1.2, markersize=10)
        ax.plot(midtimes[cell,:-1].transpose(),ifreq[cell,:-1].transpose(),
                '.b', label='Tipo S')
    elif cell == 1:
        ax.plot(spkt[cell,:], -4*np.ones(max(spkt.shape)), '|g',
                markeredgewidth=1.2, markersize=10)
        ax.plot(midtimes[cell,:-1].transpose(),ifreq[cell,:-1].transpose(),
                '.g',label='Tipo FR')
    elif cell == 2:
        ax.plot(spkt[cell,:], -2*np.ones(max(spkt.shape)), '|m',
                markeredgewidth=1.2, markersize=10)
        ax.plot(midtimes[cell,:-1].transpose(),ifreq[cell,:-1].transpose(),
                '.m',label='Tipo FF')
ax.plot(0,0,'-k', label=u'For\u00e7a')
ax.set_xlim((20,5000))
ax.set_ylim(ymin=-8)
ax.hlines(0,20,5000)
ax.set_ylabel(u'Frequ\u00eancia de disparo [Hz]')
ax.set_xlabel('tempo [ms]')
ax.legend(fontsize='small', loc='best')
#ax.text(4500,-3, "Spikes")#, bbox=dict(facecolor='white'))
ax2 = ax.twinx()
ax2.plot(t,0.1*forcesum(axis=0), '-k')
ax2.set_xlim((20,5000))
ax2.set_ylim((-.8,5.0/1.2))
ax2.set_ylabel(u'For\u00e7a [kgf]')

fig.tight_layout()
plt.show()

# Plotando o estimulo

s_delay = 20
s_predur = 2500
s_posdur = 2500
s_posamp = 22

i_stim = np.zeros(t.size)
for i in range(t.size):   
    if t[i] < s_delay:
        i_stim[i] = 0
    elif t[i] < s_delay + s_predur:
        i_stim[i] = s_posamp * (t[i]-s_delay)/s_predur
    elif (t[i] < s_delay+s_predur+s_posdur):
        i_stim[i] = s_posamp - (s_posamp * (t[i]-(s_delay+s_predur))/s_posdur)
    else:
        i_stim[i] = 0

#if i_stim.size==t.size:
#    plt.plot(t, i_stim, 'gray')


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
from matplotlib import rcParams

force1 = np.loadtxt("force1")
force2 = np.loadtxt("force2")
spkt1 = np.loadtxt("spkt1")
spkt2 = np.loadtxt("spkt2")
t = np.loadtxt("t-tetano")
freq1 = 1000./np.diff(spkt1[1:]).mean()
freq2 = 1000./np.diff(spkt2[1:]).mean()

matplotlib.rcParams.update({'font.size': 20})
plt.figure()
plt.subplot(2,1,1)
plt.plot(spkt1, -1.5*np.ones(spkt1.size), '|k', markeredgewidth=1.5, markersize = 10)
plt.plot(t, force1[:t.size].transpose(), 'k', linewidth=1.5)
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
plt.plot(t, force2[:t.size].transpose(), 'k', linewidth=1.5)
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

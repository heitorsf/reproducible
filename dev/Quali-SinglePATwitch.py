#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
Created on Fri Mar 10 15:48:51 2017

@author: heitor
'''

from netpyne import sim, specs
import numpy as np
from matplotlib import pyplot as plt
from DetectPA import getSpikes

numCells = 1
netParams = specs.NetParams()
simConfig = specs.SimConfig()

netParams.popParams['Pop1'] = {'cellType': 'MED',
                               'cellModel': 'NAPP',
                               'numCells': numCells}
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


#netParams.synMechParams['AMPA'] = {'mod': 'Exp2Syn', 'tau1': 0.1,
#                                   'tau2': 1.0, 'e': 70}
#
#netParams.stimSourceParams['dtr'] = {'type': 'NetStim', 'rate': 10,
#                                     'noise': 0.5, 'start': 1}
#netParams.stimTargetParams['dtr->Pop1'] = {'source': 'dtr',
#                                           'conds': {'popLabel': 'Pop1'},
#                                           'weight': 0.1,
#                                           'delay': 'uniform(1,5)'}

netParams.stimSourceParams['Istep'] = {'type': 'IClamp_X',
                                       'sig_type': 2,
                                       'posamp': 10,
                                       'predur': 3000,
                                       'posdur': 3000
                                       #'delay': 20,
                                       #'dur': 200,
                                       #'amp': 10
                                       }
netParams.stimTargetParams['Istep->Pop1'] = {'source': 'Istep',
                                             'conds': {'cellType': 'MED',
                                                       'cellModel': 'NAPP'},
                                             'sec': 'soma',
                                             'loc': 0.5}


#netParams.connParams['Pop1->Pop1'] = {
#    'preConds': {'popLabel': 'Pop1'}, 'postConds': {'popLabel': 'Pop1'},
#    'weight': 0.002,                    # weight of each connection
#    'delay': '0.2+gauss(13.0,1.4)',     # delay min=0.2, mean=13.0, var = 1.4
#    'threshold': 10,                    # threshold
#    'convergence': 'uniform(1,15)'}    # convergence (num presyn targeting postsyn) is uniformly distributed between 1 and 15

        
simConfig.duration = 7*1e3
simConfig.dt = 0.01
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
simConfig.analysis['plotTraces'] = {'include': [0], 'oneFigPer': 'trace'}
#simConfig.analysis['plot2Dnet'] = True  # Plot 2D net cells and connections

data = sim.create(netParams=netParams, simConfig=simConfig, output=True)
#for i in range(len(data[1])):
#    data[1][i].secs.soma.hSec.L =  75 + 5*i
#    data[1][i].secs.soma.hSec.diam = 75 + 5*i

#data[1][0].secs.soma.hSec.L = 60
#data[1][0].secs.soma.hSec.diam = 60
#
#data[1][1].secs.soma.hSec.L = 70
#data[1][1].secs.soma.hSec.diam = 70
#
#data[1][2].secs.soma.hSec.L = 90
#data[1][2].secs.soma.hSec.diam = 90
    
#data[1][0].secs.soma.hSec.__getattribute__('ena')
#data[1][0].secs.soma.hSec.__setattr__('ena', 120)
sim.simulate()
sim.analyze()
#sim.createSimulateAnalyze(netParams=netParams, simConfig=simConfig)

#v_soma_cell_0 = np.array(data[4]['v_soma']['cell_0'].to_python())

#t = np.arange(0, simConfig.duration, simConfig.dt)
t = np.array(data[4].t.to_python())

#plt.plot(t, v_soma_cell_0)

v_list = list()
for i in data[0]['Pop1'].cellGids:
    cell = 'cell_'+str(i)
    v_list.append(np.array(data[4]['v_soma'][cell].to_python()))
cells_v = np.array(v_list).transpose()

spkt, spkv = getSpikes(t, cells_v, setdvdt=20, output=True)

fi = 1000*1/np.diff(spkt)
plt.figure()
for i in range(numCells):
    plt.plot(spkt[i,2:], fi[i,1:], '.')

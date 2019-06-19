from itertools import chain

def choose_mu_type(soma,dend, mutype=False):
    if mutype:
        pass
    else:
        mutype = input('S, FR or FF: ')
    if mutype == "S":
#        print ('slow twitch')
        mus(soma,dend)
    elif mutype == "FR":
#        print ('fast twitch, fatigue resistant')
        mufr(soma,dend)
    elif mutype == "FF":
#        print ('fast twitch, fatigable')
        muff(soma,dend)
    return mutype

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

def mus0(soma,dend):
    mus(soma,dend)
    soma.L = 77.5 
    soma.diam = 77.5
    soma.gl_napp = 1/1150.0
    dend.L = 5500.0
    dend.diam = 41.5 
    dend.g_pas = 1/14400.0

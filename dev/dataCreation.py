import numpy as np
from matplotlib import pyplot as plt; plt.ion()
from neuron import h
from scipy.signal import argrelmax
import SaveParams
from interpol import cellDistrib, interpol_lin, interpol_exp, interpol_point, interpol_fit
import time

'''########## Classes ##########'''

class MotorUnitPas(object):
    """Creates a new motor unit with passive dendrite.
    
    Attributes
    ----------
    soma : NEURON Section
        Has the mechanism 'napp' inserted from napp.mod.
    dend : NEURON Section
        Has the mechanism 'pas' inserted and connected to soma.
    axon : Axon (nerlabmodel.Axon)
        Has only length and conduction velocity.
    musc : MuscUnit (nerlabmodel.MuscUnit)
        Muscular unit of the motor unit."""
    def __init__(self):
        self.type = 'S'
        self.soma = h.Section()
        self.soma.insert('napp')

        self.dend = h.Section()
        self.dend.insert('pas')
        self.dend.connect(self.soma,0,1)

        self.axon = Axon()
        self.musc = MuscUnit()

#____________________________________________________________

class MotorUnitCaL(object):
    """Creates a new motor unit with active dendrite.
    
    Attributes
    ----------
    soma : NEURON Section
        The Section has the mechanism 'napp' inserted from napp.mod.
    dend : NEURON Section
        The Section has the mechanism 'caL' inserted and connected to soma.
    axon : Axon (nerlabmodel.Axon)
        Has only length and conduction velocity.
    musc : MuscUnit (nerlabmodel.MuscUnit)
        Muscular unit of the motor unit."""
    def __init__(self):
        self.type = 'S'
        self.soma = h.Section()
        self.soma.insert('napp')

        self.dend = h.Section()
        self.dend.insert('caL')
        self.dend.connect(self.soma,0,1)

        self.axon = Axon()
        self.musc = MuscUnit()

#____________________________________________________________

class MotorUnit(object):
    """Creates a new motor unit with passive dendrite.
    
    Attributes
    ----------
    soma : NEURON Section
        Has the mechanism 'napp' inserted from napp.mod.
    dend : NEURON Section
        Has the mechanism 'pas' inserted and connected to soma.
    axon : Axon (nerlabmodel.Axon)
        Has only length and conduction velocity.
    musc : MuscUnit (nerlabmodel.MuscUnit)
        Muscular unit of the motor unit."""
    def __init__(self):
        self.type = 'S'
        self.soma = h.Section()
        self.soma.insert('napp')

        self.dend = h.Section()
        self.dend.insert('caL')
        self.gama_caL = 0.2
        self.dend.connect(self.soma,0,1)

        self.axon = Axon()
        self.musc = MuscUnit()

#____________________________________________________________

class Axon(object):
    """Creates a new axon.
    
    The object is a simple model of neuron axon based on
    a given length and a conduction velocity.
    
    Attributes
    ----------
    len : float
        The length of the axon in meters.
    velcon : float
        Conduction velocity in m/s."""
    def __init__(self):
        self.len = 0.5 # m, length
        self.velcon = 45.5 # m/s, conduction velocity

    def delay(self):
        """Calculates the time delay caused during axon propagation.
        
        Returns
        -------
        delay : float
            Given by len/velcon."""
        delay = self.len/self.velcon
        return delay 

#____________________________________________________________

class MuscUnit(object):
    """Creates a new muscular unit.
    
    The object contains paramerts of a muscular unit twitch model.

    Attributes
    ----------
    fmax : float
        Maximum twitch force [Newtons].
    tc : float
        Contraction time (or time-to-peak force) [milliseconds].
    thr : float
        Half-decay time [milliseconds].
    ftet : float
        Maximum tetanic force of the muscular unit [Newton]."""
    def __init__(self):
        self.fmax = 1.5 # N, maximum force
        self.tc = 110. # ms, contraction time
        self.thr = 92. # ms, half-decay time
        self.ftet = 10.1 # N, related to cutoff frequency

#____________________________________________________________

'''########## Functions ##########'''

def saveCellParams(fromNetpyne,data1=0,somaSec=0,dendSec=0):
    SaveParams.saveParams(fromNetpyne,data1,somaSec,dendSec)
    print "\nParameters are now saved.\n"

#____________________________________________________________

def loadCellParams(cellnumber,soma,dend):
    SaveParams.loadParams(cellnumber,soma,dend)
    print "\nParameters are now loaded.\n"

#____________________________________________________________

def createData(MotorUnitClass,numCells):
    """Adjusts the parameters of a population of neurons.

    Changes the parameters of a population of motor neurons
    based on the parameters of 3 motor neuron models of
    types S, FR and FF.

    Parameters
    ----------
    MotorUnits : MotorUnitPas, MotorUnitCaL or list of MotorUnit_like
        Motor unit object or list of motor unit objects that will have
        the parameters changed.
    active_dend : boolean or float
        Informs whether MotorUnits is/are MotorUnitPas or MotorUnitCaL.
    outputParams : boolean, optional
        If True, returns a dict with the parameters names and values.
    plotParams : boolean, optional
        If True, calls plotParameters.
    
    Returns
    -------
    params : dict
        A dict with final parameters and values.
    """
    print("\nChanging the parameters of motor neurons...")
    start = time.time()
    import pandas as pd

    MotorUnits = []
    for i in range(numCells):
        MotorUnits.append(MotorUnitClass())

    numS,numFR,numFF = cellDistrib(numCells)

    # Soma parameters
    Diam_soma = interpol_lin(numCells,80,85,100.25)
    L_soma = Diam_soma
    Gnabar = interpol_lin(numCells,0.05,0.07,0.075)
    Gnapbar = interpol_point(numCells,.00052,.0008,.00065)
    Gkfbar = interpol_point(numCells,.0028,.0040,.00135)
    Gksbar = interpol_point(numCells,.018,.037,.016) 
    Gksbar[60:75] = np.linspace(Gksbar[60],.028,15,endpoint=False)
    Gksbar[75:90] = np.linspace(.028,Gksbar[90],15,endpoint=False)
    Mact = interpol_lin(numCells,13,17,19.2)
    Rinact = interpol_lin(numCells,0.025,0.058,0.062,set_firstS=0.019)
    Gls = interpol_lin(numCells,1./1100,1./1000,1./800,set_firstS=1./1110,set_lastFF=1./700)
    
    # Dendrite parameters
    Diam_dend = interpol_lin(numCells,52,73,88,set_firstS=48.5,set_lastFF=90.)#
    L_dend = interpol_lin(numCells,6150,7450,9350)#
    if hasattr(MotorUnitClass(),"gama_caL"):
        active_dend = 1
        GcaLbar = interpol_lin(numCells, 0.00001056, 0.000007, 0.0000062)
        Vtraub_caL = interpol_point(numCells, 35, 35.6, 34)
        LTAU_caL = interpol_point(numCells, 80, 46, 47)
        Gl_caL = interpol_lin(numCells,1./12550,1./8825,1./6500,set_firstS=1./13000,set_lastFF=1./6000)#
    else:
        active_dend = 0
        Gld = interpol_lin(numCells,1./12550,1./8825,1./6500,set_firstS=1./13000,set_lastFF=1./6000)#

    # Axon parameters
    axon_len = 0.6  # meters
    Axon_velcon = interpol_lin(numCells,45.5,49.5,51.5)  # m/s

    # Set the mu.type label for each motor unit
    for i,mu in enumerate(MotorUnits):
        if i in range(0,numS):
            mu.type = 'S'
        elif i in range(numS,numS+numFR):
            mu.type = 'FR'
        elif i in range(numS+numFR,numS+numFR+numFF):
            mu.type = 'FF'

        # Fixed parameters
        mu.soma.ena = 120.0
        mu.soma.ek = -10.0
        mu.soma.el_napp = 0.0
        mu.soma.vtraub_napp = 0.0
        mu.soma.nseg = 1
        mu.soma.Ra = 70.0
        mu.soma.cm = 1.0

        # Soma parameters
        mu.soma.L = Diam_soma[i]
        mu.soma.diam = Diam_soma[i]
        mu.soma.gl_napp = Gls[i]
        mu.soma.gnabar_napp = Gnabar[i]
        mu.soma.gnapbar_napp = Gnapbar[i]
        mu.soma.gkfbar_napp = Gkfbar[i]
        mu.soma.gksbar_napp = Gksbar[i]
        mu.soma.mact_napp = Mact[i]
        mu.soma.rinact_napp = Rinact[i]

        # Dendrite parameters
        mu.dend.nseg = 1
        mu.dend.Ra = 70.0
        mu.dend.cm = 1.0
        mu.dend.L = L_dend[i]
        mu.dend.diam = Diam_dend[i]
        if active_dend>0:
            mu.dend.ecaL = 140
            mu.dend.gama_caL = active_dend
            mu.dend.gcaLbar_caL = GcaLbar[i]
            mu.dend.vtraub_caL = Vtraub_caL[i]
            mu.dend.Ltau_caL = LTAU_caL[i]
            mu.dend.gl_caL = Gl_caL[i]
            mu.dend.el_caL = 0.
        else:
            mu.dend.e_pas = 0.
            mu.dend.g_pas = Gld[i]

        # Axon parameters
        mu.axon.len = axon_len
        mu.axon.velcon = Axon_velcon[i]
    

    MN_type = [mu.type for mu in MotorUnits]

    Soma_ena = np.array([mu.soma.ena for mu in MotorUnits])
    Soma_ek = np.array([mu.soma.ek for mu in MotorUnits])
    Soma_el_napp = np.array([mu.soma.el_napp for mu in MotorUnits])
    Soma_vtraub_napp = np.array([mu.soma.vtraub_napp for mu in MotorUnits])
    Soma_nseg = np.array([mu.soma.nseg for mu in MotorUnits])
    Soma_Ra = np.array([mu.soma.Ra for mu in MotorUnits])
    Soma_cm = np.array([mu.soma.cm for mu in MotorUnits])
    Soma_L = np.array([mu.soma.L for mu in MotorUnits])
    Soma_diam = np.array([mu.soma.diam for mu in MotorUnits])
    Soma_gl_napp = np.array([mu.soma.gl_napp for mu in MotorUnits])
    Soma_gnabar_napp = np.array([mu.soma.gnabar_napp for mu in MotorUnits])
    Soma_gnapbar_napp = np.array([mu.soma.gnapbar_napp for mu in MotorUnits])
    Soma_gkfbar_napp = np.array([mu.soma.gkfbar_napp for mu in MotorUnits])
    Soma_gksbar_napp = np.array([mu.soma.gksbar_napp for mu in MotorUnits])
    Soma_mact_napp = np.array([mu.soma.mact_napp for mu in MotorUnits])
    Soma_rinact_napp = np.array([mu.soma.rinact_napp for mu in MotorUnits])
    
    Dend_nseg = np.array([mu.dend.nseg for mu in MotorUnits])
    Dend_Ra = np.array([mu.dend.Ra for mu in MotorUnits])
    Dend_cm = np.array([mu.dend.cm for mu in MotorUnits])
    Dend_L = np.array([mu.dend.L for mu in MotorUnits])
    Dend_diam = np.array([mu.dend.diam for mu in MotorUnits])
    if active_dend>0:
        Dend_ecaL = np.array([mu.dend.ecaL for mu in MotorUnits])
        Dend_gama_caL = np.array([mu.dend.gama_caL for mu in MotorUnits])
        Dend_gcaLbar_caL = np.array([mu.dend.gcaLbar_caL for mu in MotorUnits])
        Dend_vtraub_caL = np.array([mu.dend.vtraub_caL for mu in MotorUnits])
        Dend_Ltau_caL = np.array([mu.dend.Ltau_caL for mu in MotorUnits])
        Dend_gl_caL = np.array([mu.dend.gl_caL for mu in MotorUnits])
        Dend_el_caL = np.array([mu.dend.el_caL for mu in MotorUnits])
    else:
        Dend_e_pas = np.array([mu.dend.e_pas for mu in MotorUnits])
        Dend_g_pas = np.array([mu.dend.g_pas for mu in MotorUnits])
    
    db = pd.DataFrame({
        "Type" : MN_type,

        "Soma_ena" : Soma_ena,
        "Soma_ek" : Soma_ek,
        "Soma_el_napp" : Soma_el_napp,
        "Soma_vtraub_napp" : Soma_vtraub_napp,
        "Soma_nseg" : Soma_nseg,
        "Soma_Ra" : Soma_Ra,
        "Soma_cm" : Soma_cm,
        "Soma_L" : Soma_L,
        "Soma_diam" : Soma_diam,
        "Soma_gl_napp" : Soma_gl_napp,
        "Soma_gnabar_napp" : Soma_gnabar_napp,
        "Soma_gnapbar_napp" : Soma_gnapbar_napp,
        "Soma_gkfbar_napp" : Soma_gkfbar_napp,
        "Soma_gksbar_napp" : Soma_gksbar_napp,
        "Soma_mact_napp" : Soma_mact_napp,
        "Soma_rinact_napp" : Soma_rinact_napp,

        "Dend_nseg" : Dend_nseg,
        "Dend_Ra" : Dend_Ra,
        "Dend_cm" : Dend_cm,
        "Dend_L" : Dend_L,
        "Dend_diam" : Dend_diam})
    if active_dend>0:
        db.assign(Dend_ecaL = Dend_ecaL)
        db.assign(Dend_gama_caL = Dend_gama_caL)
        db.assign(Dend_gcaLbar_caL = Dend_gcaLbar_caL)
        db.assign(Dend_vtraub_caL = Dend_vtraub_caL)
        db.assign(Dend_Ltau_caL = Dend_Ltau_caL)
        db.assign(Dend_gl_caL = Dend_gl_caL)
        db.assign(Dend_el_caL = Dend_el_caL)
    else:
        db.assign(Dend_e_pas = Dend_e_pas)
        db.assign(Dend_g_pas = Dend_g_pa)


    # Create unique file name based on date and time:
    import subprocess as sp
    command = ['date']
    command.append('+%Y-%m-%d_%H-%M-%S')
    (out,err) = sp.Popen(command,stdout=sp.PIPE).communicate()
    if out.find("\n") != -1:
        dir_suffix = out.rstrip("\n")
    else:
        dir_suffix = out
    filename = "Spinal_Cord_Motor_Neurons_"+dir_suffix+".csv"
    db.to_csv(filename,index=False,na_rep='nan')

    end = time.time()
    elapsed_time = end-start
    print("Done; fitPopParams time: %.2f s" %elapsed_time)

#____________________________________________________________

def loadData(csv_file):
    """Load parameters of motor neuron from csv file.

    Changes the parameters of a population of motor neurons
    based on the parameters of 3 motor neuron models of
    types S, FR and FF.

    Parameters
    ----------
    MotorUnits : MotorUnitPas, MotorUnitCaL or list of MotorUnit_like
        Motor unit object or list of motor unit objects that will have
        the parameters changed.
    csv_file : string
        Name of the file including its extension '.csv'.
    
    Returns
    -------
    params : pandas.DataFrame
        A data frame containing loaded parameters labels and values.
    """
    print("\nLoading the parameters of motor neurons...")

    import pandas as pd
    params = pd.read_csv(csv_file)
    print("\nDone.")
    return params

#____________________________________________________________

def plotParameters(params):
    """Plots the values of parameters along a neuron pool.
    
    Tool to visualize how the parameters change along the pool.
    
    Parameter
    ---------
    params : dict
        Each key is a parameter of the neuron model."""
    for x in params.keys():
        plt.figure()
        plt.title(x)
        plt.plot(params[x])
        plt.plot([0,59,89,119],params[x][[0,59,89,119]],'rx')
        plt.plot([30,75,105],params[x][[30,75,105]],'g.')

#____________________________________________________________

def getSpikes(t, v, axDelay=False, mu=[], output=True, plotRaster=False, newFigure=False):
    """Detects spikes from time courses of membrane potentials.

    Parameters
    ----------
    t : 1-D array_like
        Time array. The shape must be (x), (x,) or (x,1)
    v : array_like
        Membrane potential array. Expected shape: (numberOfCells, timeVectorLength)
    axDelay : boolean, optional
        If True, calls axonDelay(spkt,mu), including the delay caused by axon
        conduction on spikes instants.
    mu : MotorUnit_like list, optional
        Needed when 'axDelay' is True.
    output : boolean, optional
        If False, supress the returning of spike trains and voltages
    plotRaster : boolean, optional
        If True, generates a raster plot of spike trains.
    newFigure : boolean, optional
        If True, the raster plot is generated in a new matplotlib.pyploy.figure().

    Returns
    -------
    spkt : ndarray
        Array of spikes instants.
        The shape is: (number of neurons, max number of spikes).
        Elements with no spikes are numpy.nan.
    spkv : ndarray
        Array of peaks membrane potentials.
        The shape is: (number of neurons, max number of spikes).
        Elements with no spikes are numpy.nan."""
    start = time.time()
    print "  Detecting spikes instants..."
    if not axDelay and len(mu)>0:
        raise ValueError("To get spikes with axons delays please use axDelay=True and pass a munit list as argument (mu).")
    if not isinstance(t, np.ndarray):
        t = np.array(t)
    if not isinstance(v, np.ndarray):
        v = np.array(v)
    if len(v.shape)<2:
        v.shape = (1,v.shape[0])
    #dvdt = np.diff(v, axis=1)
    dt = t[1]
    
    ''' Detect Spikes '''
    numCells = v.shape[0]
    spkt_list = []
    spkv_list = []
    for cell_id in range(numCells):
        maxima_ids = argrelmax(v[cell_id,:])
        spkt_list_cell_ids = [i for i in maxima_ids[0] if (v[cell_id,i]>=50 and v[cell_id,i]<=100)]
        spkt_cell_t = t[spkt_list_cell_ids]
        spkt_cell_v = v[cell_id,spkt_list_cell_ids]

        spkt_list.append(spkt_cell_t)
        spkv_list.append(spkt_cell_v)
    
    ''' Create spiketrain times (spkt) and volts (spkv) matrices '''
    # Set matrices size according to the neuron that fired the biggest number of spikes
    numPAs = np.array([])
    for i in range(numCells):
        numPAs = np.append(numPAs, len(spkt_list[i]))
    maxPAs = int(numPAs.max())
    
    # Create and fill the matrices
    spkt = np.zeros((numCells, maxPAs))*np.nan
    spkv = np.zeros(spkt.shape)*np.nan
    for i in range(numCells):
        for j in range(len(spkt_list[i])):
            spkt[i,j] = spkt_list[i][j]
            spkv[i,j] = v[i,:][t==spkt_list[i][j]]

    if axDelay:
        spkt = axonDelay(spkt,mu)
    elapsed_time = time.time() - start
    print ("  Done; getSpikes time: %.2f s" %elapsed_time)

    ''' Plot Raster '''
    if plotRaster:
        if newFigure:
            plt.figure()
        plotSpikes(spkt)
    if output:
        return spkt, spkv

#____________________________________________________________

def plotSpikes(spkt,spkid=0,newFigure=False,color='b'):
    """Plot spikes instants."""
    start_plot = time.time()
    doPlot = False
    #if len(spkt[:,:]) > 0:
    print "\n  Plotting Raster..."
    if newFigure:
        plt.figure()
    if len(spkt.shape)==2 and spkid==0:
        try:
            for i in range(len(spkt[:,0])):
                plt.plot(spkt[i,:], np.ones(len(spkt[i,:]))*i, color+'|', markeredgewidth=1.5)
        except IndexError:
            print("\n*****\n\nWarning: 'spkt' has no length, aparently no spikes were fired.\n\n*****\n")
        else:
            doPlot = True
    elif len(spkt.shape)==1 and spkid is not 0:
        try:
            plt.plot(spkt,spkid, color+'|', markeredgewidth=1.5)
        except IndexError:
            print("\n*****\n\nWarning: 'spkt' has no length, aparently no spikes were fired.\n\n*****\n")
        else:
            doPlot = True
    else:
        print("\n*****\n\nWarning: check 'spkt'.\n\n*****\n")
    
    if doPlot:
        plt.title("Spike Trains")
        plt.xlabel("Time [ms]")
        plt.ylabel("# Motor neuron")
        plt.show()
        elapsed_time = time.time() - start_plot
        print("  Done; plotRaster time: %.2f s" %elapsed_time)
    else:
        elapsed_time = time.time() - start_plot
        print("  Done with warnings; plotRaster time: %.2f s" %elapsed_time)

#____________________________________________________________

def axonDelay(spkt_soma,mu):
    """Add axon conduction delay to spike instants."""
    if spkt_soma.shape[0] != len(mu):
        raise ValueError("Number of spike trains differs from number of motor units.")
    spkt_endplate = np.zeros(spkt_soma.shape)
    for i in range(spkt_soma.shape[0]):
        spkt_endplate[i] = spkt_soma[i] + mu[i].axon.delay() 
    return spkt_endplate

#____________________________________________________________

def sig(f,c):
    """Saturate data using a sigmoidal function.

    Parameters
    ----------
    f : array_like
        Data to saturate.
    c : float
        Parameter to adjust the shape of the sigmoid.

    Returns
    -------
    out : array
        Saturated data."""
    b = 1.2
    expcf = np.exp(-c*(f-b))
    sat = (1-expcf) / (1+expcf)
    out = (sat - sat[0]) / (1-sat[0])
    return out

def calculateTetTwt(C):
    """Calculate tetanus/twitch ratio.
    
    Calculate the tetanus/twitch ratio for every sigmoid with
    parameter c in vector C.
    
    Parameters
    ----------
    C : 1-D array, float
        Sequence of values to be used as parameter of a sigmoid function.
    
    Returns
    -------
    ttnp : numpy array
        Value of a sigmoid curve at the index where a straight line of the same
        size is at 1, meaning the tetanus/twitch ratio. It it:
        fsat: sigmoid
        fmu: straight line (from 0 to 10)
            ttnp = 1/fsat[fmu==1]"""
    tt = []
    fmu = np.arange(0,10,0.01)
    for c in C:
        fsat = sig(fmu,c)
        tt.append(1./fsat[fmu==1])
    ttnp = np.array(tt)
    return ttnp

def muscularForce(spkt, t, plotForce=False, newFigure=True, fitStyle='exponential'):
    """Generate muscle force signal from spike trains.

    Use digital filter to generate summations of twitches in response
    to spike trains.
    
    Parameters
    ----------
    spkt : array
        Spike train(s). If not 1-D, first dimension is motor unit number.
    t : array
        Time array.
    plotForce : boolean, optional
        Choose whether to plot force signals or not.
    newFigure : boolean, optional
        If True, creates a new matplotlib.pyploy.figure() when plotForce is True.
    fitStyle : string, optional
        Decides how the parameters will vary. Accepts "exponential" (default)
        or "linear".
    
    Returns
    -------
    forces : array
        An array with shape (spkt.shape[0], t.size) containing force signal
        for each spike train. 
        """

    print("\n  Generating muscular force data...")
    print("    Setting muscular force parameters.")
    start = time.time()

    from scipy.signal import lfilter

    # Here we assume that every cell generated a spike train:
    numCells = spkt.shape[0] 

    dt = t[1]
    forces = np.zeros(t.size)
    
    if fitStyle=="linear":
        Fmaxmu = interpol_lin(numCells,1.7,1.9,2.15)  # N
        Tc = interpol_lin(numCells,110.,73.5,55.)  # ms
        Ftet = interpol_lin(numCells,10.1,15.2,19.3)  # N
    elif fitStyle=="exponential":
        fmax = 260.0*1e-3  # N
        fmin = 1.*1e-3  # N
        Fmaxmu = interpol_exp(numCells,fmin,fmax,ascending=True)
        tcmin = 30.  # ms
        tcmax = 100.  # ms
        Tc = interpol_exp(numCells,tcmin,tcmax,ascending=False)
    else:
        raise ValueError("fitStyle not understood, must be \"linear\" or \"exponential\"")

    print("    Calculating auxiliary arrays.")
    A = np.array([])
    B = np.array([])
    for i in range(numCells):
        A = np.concatenate((A, np.array([1, -2*np.exp(-(dt/Tc[i])), np.exp(-2*(dt/Tc[i]))])))
        B = np.concatenate((B, np.array([0, 1*((dt**2)/Tc[i])*np.exp(1-dt/Tc[i])])))
    A = A.reshape((numCells,3))
    B = B.reshape((numCells,2))

    print("    Computing force.")
    C = interpol_exp(numCells,xmin=1.68,xmax=2.25,ascending=False)  # Fsat from 15 to 65 Hz
    TetTwt = calculateTetTwt(C)
    for cell_index in range(spkt.shape[0]):  # For every cell (or for every spike train):
        Fmu = np.zeros(t.size)
        spkt_cell = spkt[cell_index][spkt[cell_index]>0]
        spikes = np.zeros(t.size)
        for instant in spkt_cell:  # For every spike of this cell:
            spikes[int(instant/dt)] = 1./dt
        F = lfilter(B[cell_index],A[cell_index],spikes)
        Fsat = sig(F,C[cell_index]) * TetTwt[cell_index] * Fmaxmu[cell_index]
        forces = np.vstack((forces, Fsat))

    if plotForce:
        print("    Plotting...")
        if newFigure:
            plt.figure()
        plotForces(t,forces)

    elapsed_time = time.time() - start
    print("  Done; muscularForce time: %.2f s\n" %elapsed_time)
    return forces[1:,] #, TetTwt, C,Fmaxmu

#____________________________________________________________

def plotForces(t,forces):
    """Plot force signal."""
    plt.title("FDI Muscular Force")
    plt.xlabel("Time [ms]")
    plt.ylabel("Muscular Force [N]")
    plt.plot(t,forces.sum(axis=0))
    plt.show()

#____________________________________________________________

def instFreq(spkt):
    """Returns instantaneous frequency and midtimes from spike trains."""
    ifreq = 1000./np.diff(spkt) # 1/s (Hz)
    midtimes = spkt[:,1:] - np.diff(spkt)/2. # put the point in the middle of the time samples pairs
    return ifreq,midtimes

#____________________________________________________________

# !/usr/local/bin/python
# -*- coding: utf-8 -*-

"""
My visualisation tools
"""
import pylab
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
from nested_dict import nested_dict
import numpy as np
from numpy import exp
from numpy import linspace
import scipy.optimize as opti
import scipy.fftpack as fftpack
rc('mathtext', default='regular')
N = 100 # number of points per dimension

def mysine(t, f2, amp, off):
    return amp * np.sin(2. * np.pi * f2 * t) + off

def show_ISI_histogram(spiketime_array, simtime, n_bins, deltat, **keywords):
    if len(spiketime_array) == 0:
        print('no spikes recorded')
        exit(1)
    # find number of active neurons
    if deltat >= 1.0: print('ISI may not be representative due to too low recording frequency')
    # First generate a vector with the number of ISI of any binned length
    interspikes = np.diff(spiketime_array)
    stepwidth = simtime / n_bins
    axis_x = np.linspace(interspikes.min(), interspikes.max(), n_bins)
    interspikes_cumcount = np.zeros_like(axis_x)
    for i in axis_x:
        interspikes = interspikes - stepwidth
        interspikes[ interspikes < 0. ] = 0.
        interspikes_cumcount[i] = np.count_nonzero(interspikes[ interspikes < stepwidth ])
    # generate the graph
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = plt.hist(interspikes_cumcount, axis_x, cumulative=True)
#    y_min = offset
#    ax.set_xlim([ t_min, t_max ])
#    ax.set_ylim([ y_min, y_max ])
    ax.set_xlabel('ISI length [ms]', size=12)
    ax.set_ylabel('Occurences in population', size=12)
    if 'output' in keywords:
        plt.savefig(os.path.join(self.output_dir, self.label +
                                 '_Dotplot_' + area + '.' + keywords[ 'output' ]))
    else:
        fig.show()


def show_Rasterplot(figu, spiketime_array, spiketime_senders, t_recstart, t_sim, n_rows=50, nsubplot=111, mingid = 0):
    spike_senders = spiketime_senders
    spike_senders -= mingid
    spike_senders = spike_senders[ spiketime_senders <= n_rows ]
    spike_times2 = spiketime_array[ spiketime_senders <= n_rows ]
    axis = figu.add_subplot(nsubplot)
    axis.set_xlim([ t_recstart, t_sim ])
    axis.set_xlabel('time [ms]', size=9)
    axis.set_ylabel('Neuron ID', size=9)
    #axis.set_ylim([ min(spike_senders)-0.5, max(spike_senders)+0.5 ])
    axis.plot(spike_times2, spike_senders, ".k")  # <----- THIS IS THE RASTER PLOT
    return axis

def show_Rasterplot2(figu, spiketime_array, spiketime_senders, t_recstart, t_sim, n_rows=50, nsubplot=111, mingid=0):
    print('{0} spikes recorded'.format(len(spiketime_array)))
    # spike_senders = np.zeros_like(spiketime_senders)
    # spike_senders = spiketime_senders
    spike_senders = spiketime_senders - mingid
    # spike_senders = spike_senders[ spiketime_senders <= n_rows ]
    # spike_times2 = spiketime_array[ spiketime_senders <= n_rows ]
    axis = figu.add_subplot(nsubplot)
    axis.set_xlim([ t_recstart, t_sim ])
    axis.set_xlabel('time [ms]', size=9)
    axis.set_ylabel('Neuron ID', size=9)
    axis.set_ylim([0.,51.])
    #if not len(spike_senders) == 0:
    #    axis.set_ylim([ min(spike_senders)-0.5, max(spike_senders)+0.5 ])
    axis.plot(spiketime_array, spike_senders, ".k")  # <----- THIS IS THE RASTER PLOT
    return axis, spike_senders, spiketime_array


# ----------------- core routine -------------------------
def predict_str_freq(t1=100., g=6.5, g1=30., C=250., remote=False):
    a = g * t1 / C
    b = g1 * t1 / C
    if (b < np.sqrt((a+1.)*(a+1.) +1.) -1. -a) & (remote == False):
        print('STR will not occur')
    h0 = a + b + 1.
    h1 = a + 1.
    h2 = np.sqrt(h0 * h0 - h1 * h1) - 1.
    if h2 >= 0:
        freq = np.sqrt(h2) / (2 * np.pi * t1)
    else:
        freq = 0
    return freq * 1000.

def predict_str_freq_dep_t(t1=100., t=15, g1=30., C=250., remote=False):
    a = t1 / t
    b = g1 * t1 / C
    if (b < np.sqrt((a+1.)*(a+1.) +1.) -1. -a) & (remote == False):
        print('STR will not occur')
    h0 = a + b + 1.
    h1 = a + 1.
    h2 = np.sqrt(h0 * h0 - h1 * h1) - 1.
    if h2 >= 0:
        freq = np.sqrt(h2) / (2 * np.pi * t1)
    else:
        freq = 0
    return freq * 1000.


# -------------------------------- PREDICTION FOR RESONANCE --------------------------------
def show_freqs_dep_t1C(maxt1=150.0, mint1=1.5, minC=250.0, maxC=550., g=6.5, g1 = 100., N=N):
    time_constants = np.linspace(mint1, maxt1, N)
    capacities = np.linspace(minC, maxC, N)
    time_constants2, capacities2 = np.meshgrid(time_constants, capacities)
    freqs2 = np.zeros((N,N))
    for i in np.arange(0, N):
        for j in np.arange(0, N):
            freqs2[i,j] = predict_str_freq(t1=time_constants[i], g=g, g1=g1, C=capacities[j], remote=True)
    fig2 = plt.figure()
    ax = fig2.gca(projection='3d')
    ax.plot_surface(time_constants2, capacities2, freqs2, cmap=cm.jet, linewidth=0.2, cstride=1, antialiased=False)
    ax.set_xlabel('Time constants [ms]')
    ax.set_ylabel('Membrane Capacity [pF]')
    ax.set_zlabel('Frequency [Hz]')
    plt.ion()
    plt.show()

def show_freqs_dep_t1g1(f_r=10., maxt1=150.0, mint1=1.5, ming1=200.0, maxg1=550., C=250.0, g=6.5,
                        N=N, contours=False, myplane=False, spec_intersect=False, only_fitting=False):
    time_constants = np.linspace(mint1, maxt1, N)
    conductance = np.linspace(ming1, maxg1, N)
    time_constants2, conductance2 = np.meshgrid(time_constants, conductance)
    freqs2 = np.zeros((N,N))
    fitting = np.zeros((2,0))
    for i in np.arange(0, N):
        for j in np.arange(0, N):
            freqs2[i,j] = predict_str_freq(t1=time_constants[i], g=g, g1=conductance[j], C=C, remote=True)
            if abs(freqs2[i,j] - f_r) <= 1e-1:
                fitting = np.concatenate((fitting, np.array([[time_constants[i]], [conductance[j]]])), axis=1)
    fig1 = plt.figure()
    # ax = fig1.gca(projection='3d')
    if not only_fitting:
        ax = fig1.add_subplot(111, projection='3d')
        ax.plot_surface(time_constants2, conductance2, freqs2,
                    cmap=cm.jet, linewidth=0.2, antialiased=False, alpha = 0.5, cstride = 1)
    else:
        ax = fig1.add_subplot(111)
        ax.plot(fitting[0], fitting[1])# color='r')
        # print('fitting: {0}'.format(fitting))
    if contours: # add contours to the side so the graph
        cset = ax.contourf(time_constants2, conductance2, freqs2, zdir='z', offset=-40, cmap=cm.coolwarm)
        cset = ax.contourf(time_constants2, conductance2, freqs2, zdir='x', offset=mint1-40, cmap=cm.coolwarm)
        cset = ax.contourf(time_constants2, conductance2, freqs2, zdir='y', offset=ming1-40, cmap=cm.coolwarm)
    if myplane: # display the 10 Hz plane
        myplane = np.ones_like(freqs2) * 10.0
        ax.plot_surface(time_constants2, conductance2, myplane,
                        cmap=cm.jet, linewidth=0, antialiased=False, alpha = 0.99, color='k')
    if spec_intersect: # show the intersection of the graph with the 10 Hz plane
        #myintersect_bool = np.isclose(freqs2, 10.*np.ones_like(freqs2), 1e-2)
        #myintersect = np.zeros_like(freqs2)
        #myintersect[myintersect_bool == True] += 10.
        myintersectX = np.zeros(0)
        myintersectY = np.zeros(0)
        for i in np.arange(0,N):
            for j in np.arange(0, N):
                if np.isclose(freqs2[i, j], 10., 1e-3):
                    myintersectX = np.append(myintersectX, time_constants[ i ])
                    myintersectY = np.append(myintersectY, conductance[ j ])
        myintersectZ = 10. * np.ones_like(myintersectX)
        # fig2 = plt.figure()
        # ax = fig1.add_subplot(111, projection='3d')
        ax.plot(myintersectX, myintersectY, myintersectZ) # This currently generates errors
    ax.set_xlabel('Time constants [ms]')
    ax.set_ylabel('Conductance g1 [nS]')
    if not only_fitting:
        ax.set_zlabel('Frequency [Hz]')
    # ax.set_zlim([10., 200.])
    plt.ion()
    plt.show()


def show_freqs_dep_gg1(f=10., maxg=200.0, ming=1.5, t1=100., C=250.0, N=N):
    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    conductance1 = np.linspace(ming, maxg, N)
    conductance = np.linspace(ming, maxg, N)
    freqs2 = np.zeros((N,N))
    fitting = np.zeros((2,0))
    for i in np.arange(0, N):
        for j in np.arange(0, N):
            freqs2[i,j] = predict_str_freq(t1=t1, g=conductance[i], g1=conductance1[j], C=C, remote=True)
            if abs(freqs2[i,j] - f) <= 1e-1:
                fitting = np.concatenate((fitting, np.array([[conductance1[i]], [conductance[j]]])), axis=1)
    ax.scatter(fitting[ 0 ], fitting[ 1 ])  # color='r')
    ax.set_xlabel('Conductance g [nS]')
    ax.set_ylabel('Conductance g1 [nS]')
    plt.ion()
    plt.show()


def show_freqs_dep_g1(t1=100., ming1=200.0, maxg1=550., C=250.0, g=10., N=10*N, returning=False):
    conductance = np.linspace(ming1, maxg1, N)
    freqs=np.zeros(N)
    for i in np.arange(0,N):
        freqs[ i ] = predict_str_freq(t1=t1, g=g, g1=conductance[ i ], C=C, remote=True)
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.plot(conductance, freqs)
    ax.set_xlabel('Conductance $g_1$ [nS]')
    ax.set_ylabel('Frequency $f_R$ [Hz]')
    if returning:
        return conductance, freqs
    else:
        plt.ion()
        plt.show()


def show_freqs_dep_t1(mint1=35.0, maxg1=65., g1=32.5, C=250.0, g=25.0, N=10*N, returning=False):
    conductance = np.linspace(ming1, maxg1, N)
    freqs=np.zeros(N)
    for i in np.arange(0,N):
        freqs[ i ] = predict_str_freq(t1=t1, g=g, g1=conductance[ i ], C=C, remote=True)
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.plot(conductance, freqs)
    ax.set_xlabel('Secondary time constant $\tau_1$ [nS]')
    ax.set_ylabel('Spike frequency $f_R$ [Sp/s]')
    if returning:
        return conductance, freqs
    else:
        plt.ion()
        plt.show()



def show_freqs_dep_g(t1=100., ming=10.0, maxg=100., C=250.0, g1=10., N=10*N, returning=False):
    conductance = np.linspace(ming, maxg, N)
    freqs=np.zeros(N)
    for i in np.arange(0,N):
        freqs[ i ] = predict_str_freq(t1=t1, g=conductance[ i ], g1=g1, C=C, remote=True)
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.text(0.01, 0.95, '$C = ${0}pF, $g_1 = ${1}nS, $t_1 = ${2}ms'.format(C, g1, t1),
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes)
    ax.plot(conductance, freqs)
    ax.set_xlabel('Membrane Conductance $g$ [nS]')
    ax.set_ylabel('Frequency $f_R$ [Hz]')
    if returning:
        return conductance, freqs
    else:
        plt.ion()
        plt.show()

def show_freqs_dep_t(t1=100., mint=10.0, maxt=100., C=250.0, g1=10., N=10*N, returning=False):
    timescale = np.linspace(mint, maxt, N)
    freqs=np.zeros(N)
    for i in np.arange(0, N):
        freqs[ i ] = predict_str_freq_dep_t(t1=t1, t=timescale[ i ], g1=g1, C=C, remote=True)
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.text(0.01, 0.95, '$C = ${0}pF, $g_1 = ${1}nS, $t_1 = ${2}ms'.format(C, g1, t1),
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes)
    ax.plot(timescale, freqs)
    ax.set_xlabel(r'Membrane timescale $\tau$ [ms]')
    ax.set_ylabel('Frequency $f_R$ [Hz]')
    ax.set_xlim([ np.floor(mint), np.ceil(maxt) ])
    ax.set_ylim([ 5, 12 ])
    if returning:
        return timescale, freqs
    else:
        plt.ion()
        plt.show()


def show_freqs_dep_C(minC=50.0, maxC=500.0, g1=10., t1=100., g=10.0, N=10*N, returning=False):
    capacity = np.linspace(minC, maxC, N)
    freqs=np.zeros(N)
    for i in np.arange(0,N):
        freqs[ i ] = predict_str_freq(t1=t1, g=g, g1=g1, C=capacity[ i ], remote=True)
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.plot(capacity, freqs)
    ax.text(0.01, 0.95, '$g = ${0}nS, $g_1 = ${1}nS, $t_1 = ${2}ms'.format(g, g1, t1),
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes)
    ax.set_xlabel('Capacity $C_m$ [nS]')
    ax.set_ylabel('Frequency $f_R$ [Hz]')
    if returning:
        return capacity, freqs
    else:
        plt.ion()
        plt.show()

def show_freqs_dep_gC(maxg=200.0, ming=1.5, t1=100., minC=50.0, maxC=500., g1=30., f_r=10., N=N):
    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    capacity = np.linspace(minC, maxC, N)
    conductance = np.linspace(ming, maxg, N)
    freqs2 = np.zeros((N,N))
    fitting = np.zeros((2,0))
    for i in np.arange(0, N):
        for j in np.arange(0, N):
            freqs2[i,j] = predict_str_freq(t1=t1, g=conductance[i], g1=g1, C=capacity[j], remote=True)
            if abs(freqs2[i,j] - f_r) <= 1e-1:
                fitting = np.concatenate((fitting, np.array([[conductance[i]], [capacity[j]]])), axis=1)
    ax.scatter(fitting[ 0 ], fitting[ 1 ])  # color='r')
    ax.set_title('Membrane timescale parameters that yield {0} Hz resonance frequency'.format(f_r))
    ax.text(0.01, 0.95, '$g_1 = ${0}nS, $t_1 = ${1}ms'.format(g1, t1),
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes)
    ax.set_xlabel('Conductance $g$ [nS]')
    ax.set_ylabel('Capacity $C_m$ [pF]')
    plt.ion()
    plt.show()



def show_some_ideas():
    print('good parameter set: predict_str_freq(t1=120., g=6.5, g1=112.24, C=250.)')
    print('good parameter set: predict_str_freq(t1=37.40475, g=10., g1=30., C=250.)')
    print('good parameter set: predict_str_freq(t1=120., g=10., g1=109.04, C=250.)')
    print('good parameter set: predict_str_freq(t1=100., g=10., g1=89.466, C=250.)')

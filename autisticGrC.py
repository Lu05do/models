import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz
from experimental_setup import add_gabazine
from patch import p
from dbbs_models import GranuleCell
from experimental_setup import create_autistic_model, current_injection, voltage_clamp, synapticStim
import copy
from plotly.subplots import make_subplots
import plotly.graph_objects as go
#from sklearn.linear_model import LinearRegression
import scipy.optimize


AutisticGranuleCell = create_autistic_model()


def create_cell(model):
    cell = model()
    for dend in cell.dendrites:
        cell.create_synapse(dend, "AMPA")
        cell.create_synapse(dend, "NMDA")


def analysis(*models, durForFreqCalculation=0, time=0, V=0):
    vm = V[0]['Vm']
    spikes = V[0]['spikes']
    if len(spikes) > 1:
        freq= len(spikes)/durForFreqCalculation * 1e3
    else:
        freq=0
    print("mean freq: ", freq)
    Im = V[0]['Im']
    r_Trelease_n_NMDA=V[0]['difGluNMDA']
    r_Trelease_n_AMPA=V[0]['difGluAMPA']
    r_Trelease_n_NMDA = np.array(list(list(r) for r in r_Trelease_n_NMDA))
    r_Trelease_n_NMDA_sum = np.sum(r_Trelease_n_NMDA.T, axis=1)
    r_Trelease_n_AMPA = np.array(list(list(r) for r in r_Trelease_n_AMPA))
    r_Trelease_n_AMPA_sum = np.sum(r_Trelease_n_AMPA.T, axis=1)
    AMPAcurr = V[0]['AMPAcurr']
    NMDAcurr = V[0]['NMDAcurr']
    AMPAcurr = np.array(list(list(r) for r in AMPAcurr))
    NMDAcurr = np.array(list(list(r) for r in NMDAcurr))
    AMPAcurr_sum = np.sum(AMPAcurr.T, axis=1)
    NMDAcurr_sum = np.sum(NMDAcurr.T, axis=1)

    fig = make_subplots(rows=2, cols=2,subplot_titles=('current [nA]', 'voltage [mV]', 'AMPA', 'NMDA') )
    fig.add_trace(
    go.Scatter(x=time, y=Im),
    row=1, col=1
    )
    fig.add_trace(
    go.Scatter(x=time, y=vm),
    row=1, col=2
    )
    fig.add_trace(
    go.Scatter(x=time, y=AMPAcurr_sum),
    row=2, col=1
    )
    fig.add_trace(
    go.Scatter(x=time, y=NMDAcurr_sum),
    row=2, col=2
    )
    fig['layout']['xaxis']['title']='time [ms]'
    fig.show()


def analysisSyn(*models, time=0, V=0):
    # Current clamp paper Soda et al fig 1B
    vm = V[0]['Vm']
    spikes = V[0]['spikes']
    r_Trelease_n_NMDA=V[0]['difGluNMDA']
    r_Trelease_n_AMPA=V[0]['difGluAMPA']
    r_Trelease_n_NMDA = np.array(list(list(r) for r in r_Trelease_n_NMDA))
    r_Trelease_n_NMDA_sum = np.sum(r_Trelease_n_NMDA.T, axis=1)
    r_Trelease_n_AMPA = np.array(list(list(r) for r in r_Trelease_n_AMPA))
    r_Trelease_n_AMPA_sum = np.sum(r_Trelease_n_AMPA.T, axis=1)
    AMPAcurr = V[0]['AMPAcurr']
    NMDAcurr = V[0]['NMDAcurr']
    AMPAcurr = np.array(list(list(r) for r in AMPAcurr))
    NMDAcurr = np.array(list(list(r) for r in NMDAcurr))
    AMPAcurr_sum = np.sum(AMPAcurr.T, axis=1)
    NMDAcurr_sum = np.sum(NMDAcurr.T, axis=1)

    fig = make_subplots(rows=3, cols=3,subplot_titles=( 'voltage [mV]', 'AMPA curr [nA]','NMDA curr [nA]','  ', 'Glu difWave through AMPA', 'Glu difWave through NMDA', ) )
    fig.add_trace(
    go.Scatter(x=time, y=vm),
    row=1, col=1
    )
    fig.add_trace(
    go.Scatter(x=time, y=AMPAcurr_sum),
    row=1, col=2
    )
    fig.add_trace(
    go.Scatter(x=time, y=NMDAcurr_sum),
    row=1, col=3
    )
    fig.add_trace(
    go.Scatter(x=time, y=r_Trelease_n_AMPA_sum),
    row=2, col=2
    )
    fig.add_trace(
    go.Scatter(x=time, y=r_Trelease_n_NMDA_sum),
    row=2, col=3
    )
    fig['layout']['xaxis']['title']='time [ms]'
    fig.show()


def experiment_1():
    """
    This is where Ludovica will write a detailed description
    of experiment 1.
    """
    gc = create_cell(GranuleCell)
    aut_gc = create_cell(AutisticGranuleCell)

    injectDuration = 800  # ms
    # Inject current into wildtype cell
    time, (wt_results, aut_results) = current_injection(
        gc, aut_gc,
        durationPre=200,
        durationInject=injectDuration,
        durationPost=200,
        inj_current= 2* 1e-3,
        current_for_holding = 0,
        temperature=32,
        v_init=-80,
        dt = 0.025,
    )

    

def experiment_2():
    """
    This is where Ludovica will write a detailed description
    of experiment 2.
    """
    gc = create_cell(GranuleCell)
    aut_gc = create_cell(AutisticGranuleCell)

    time, V = synapticStim(
        gc,
        durationTot=500,
        temperature=32,
        v_init=-65,
        dt = 0.025,
        startTrain=50,
        nPreSynSpikes=10,
        ISI=10
    )
    analysisSyn(
        gc,
        time=time,
        V=V
    )
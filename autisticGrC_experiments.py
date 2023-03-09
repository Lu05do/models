import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz
from patch import p
from dbbs_models import GranuleCell
from experimental_setup import *
import copy
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from sklearn.linear_model import LinearRegression
import scipy.optimize


AutisticGranuleCell = create_autistic_model()

#print(AutisticGranuleCell().soma[0].diam)

def create_cell(model):
    cell = model()
    for dend in cell.dendrites:
        cell.create_synapse(dend, "AMPA")
        cell.create_synapse(dend, "NMDA")
    return cell

def analysis(results, time, durForFreqCalculation, list_freq):
    for result in results:
        vm = result['Vm']
        spikes = result['spikes']
        if len(spikes) > 1:
            freq= len(spikes)/durForFreqCalculation * 1e3
        else:
            freq=0
        print("mean freq: ", freq)
        Im = result['Im']
        r_Trelease_n_NMDA=result['difGluNMDA']
        r_Trelease_n_AMPA=result['difGluAMPA']
        r_Trelease_n_NMDA = np.array(list(list(r) for r in r_Trelease_n_NMDA))
        r_Trelease_n_NMDA_sum = np.sum(r_Trelease_n_NMDA.T, axis=1)
        r_Trelease_n_AMPA = np.array(list(list(r) for r in r_Trelease_n_AMPA))
        r_Trelease_n_AMPA_sum = np.sum(r_Trelease_n_AMPA.T, axis=1)
        AMPAcurr = result['AMPAcurr']
        NMDAcurr = result['NMDAcurr']
        AMPAcurr = np.array(list(list(r) for r in AMPAcurr))
        NMDAcurr = np.array(list(list(r) for r in NMDAcurr))
        AMPAcurr_sum = np.sum(AMPAcurr.T, axis=1)
        NMDAcurr_sum = np.sum(NMDAcurr.T, axis=1)

        fig = make_subplots(rows=2, cols=2, subplot_titles=('current [nA]', 'voltage [mV]', 'AMPA', 'NMDA'))
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
        list_freq.append(freq)


def analysisSyn(results, time=0, V=0):
    # Current clamp paper Soda et al fig 1B
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


currents_gc = list()
freq_wt = list()
def experiment_gc():
    
    for i in np.arange(2, 24, 2):
        inj = i
        currents_gc.append(inj *1e-3)
        gc = create_cell(GranuleCell)
        injectDuration = 800  # ms
        # Inject current into wildtype cell

        time, results = current_injection(   
            gc, 
            durationPre=200,
            durationInject=injectDuration,
            durationPost=200,
            inj_current= inj* 1e-3,
            current_for_holding = 0,
            temperature=32,
            v_init=-80,
            dt = 0.025,
        )
        analysis(results, time, injectDuration, freq_wt)
        
        
currents_aut = list()
freq_aut = list()
def experiment_aut():

    for i in np.arange(2, 24, 2):
        inj = i
        currents_aut.append(inj *1e-3)
        aut_gc = create_cell(AutisticGranuleCell)
        injectDuration = 800  # ms
        # Inject current into autistic cell

        time, results = current_injection(   
            aut_gc, 
            durationPre=200,
            durationInject=injectDuration,
            durationPost=200,
            inj_current= inj* 1e-3,
            current_for_holding = 0,
            temperature=32,
            v_init=-80,
            dt = 0.025,
        )
        analysis(results, time, injectDuration, freq_aut)


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


if __name__ == "__main__":
     experiment_gc()

if __name__ == "__main__":
      experiment_aut()


fig = make_subplots(rows=1, cols=1)
fig.add_trace(
    go.Scatter(x=currents_gc, y=freq_wt, name='WT',
          mode='lines+markers',
         line=go.scatter.Line(color="black"),
         marker=dict(
             color='black',
             size=2)), row=1, col=1
 )
fig.add_trace(
    go.Scatter(x=currents_aut, y=freq_aut, name='KO',
          mode='lines+markers',
          line=go.scatter.Line(color="blue"),
          marker=dict(
              color='blue',
              size=2)),row=1, col=1
  )
fig['layout']['xaxis']['title']='Inj current [pA]'
fig['layout']['yaxis']['title']='Spike freq [Hz]'
fig.show()
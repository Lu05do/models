import matplotlib.pyplot as plt
from arborize.core import flatten_composite
import numpy as np
from numpy import trapz
from patch import p
from dbbs_models import GranuleCell
import copy
import copy
from experimental_setup import *
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import scipy.optimize

def create_wildtype_model():
    # Create template WT
    class GrC(GranuleCell):
        morphologies = copy.deepcopy(GranuleCell.morphologies)
        section_types = copy.deepcopy(GranuleCell.section_types)
        synapse_types = copy.deepcopy(GranuleCell.synapse_types)

    return GrC

grc = create_wildtype_model()
grc= GranuleCell()
for dend in grc.dendrites:
    grc.create_synapse(dend, "AMPA")
    grc.create_synapse(dend, "NMDA")

r_Trelease_n_NMDA= []
r_Trelease_n_AMPA= []
NMDAcurr=[]
AMPAcurr=[]

def create_cell(model):
    cell = model()
    for dend in cell.dendrites:
        cell.create_synapse(dend, "AMPA")
        cell.create_synapse(dend, "NMDA")
    return cell
grc=GranuleCell()

def current_injection(*models, nDend=4, durationPre=100, durationInject=800, durationPost=100, temperature=32, v_init=-80, dt = 0.025, inj_current=4, current_for_holding=0):
    p.time
    p.celsius = temperature
    p.v_init = v_init
    p.dt = dt
    for model in models:
        model.record_soma(grc)
        model._spike_detector = p.NetCon(model.soma[0], None)
        clamp = model.soma[0].iclamp(x=0.5, delay=durationPre, duration=durationInject, amplitude=inj_current)
        model.soma[0].iclamp(x=0.5, delay=0, duration=durationPre+durationInject+durationPost, amplitude=current_for_holding)
        i = p.record(clamp._ref_i)
        model.spike_times = p.Vector()
        model._spike_detector.record(model.spike_times)
        for dend in range(nDend):
            AMPAsyn=model.dendrites[dend]._synapses[0]
            NMDAsyn=model.dendrites[dend]._synapses[1]
            #NMDAsyn.stimulate(start=100, number=2, interval=200)
            NMDAcurr.append(NMDAsyn.record())
            AMPAcurr.append(AMPAsyn.record())
            r_tr_NMDA = p.record(NMDAsyn._point_process._ref_Trelease)
            r_tr_AMPA = p.record(AMPAsyn._point_process._ref_Trelease)
            r_Trelease_n_NMDA.append(r_tr_NMDA)
            r_Trelease_n_AMPA.append(r_tr_AMPA)

    #v = p.record(clamp._ref_v)
    p.finitialize(v_init)
    p.continuerun(durationPre + durationInject + durationPost)

    return list(p.time), [{
        "Vm": list(model.Vm),
        "spikes": list(model.spike_times),
        "model": model.__class__.__name__,
        "Im": list(i),
        "NMDAcurr": list(NMDAcurr),
        "AMPAcurr": list(AMPAcurr),
        "difGluAMPA": list(r_Trelease_n_AMPA),
        "difGluNMDA": list(r_Trelease_n_NMDA)
    } for model in models]

def analysis(results, time, durForFreqCalculation, list_freq=list):
    for result in results:
        vm = result['Vm']
        spikes = result['spikes']
        if len(spikes) > 1:
            freq= len(spikes)/durForFreqCalculation * 1e3
        else:
            freq=0
        print("mean freq: ", freq)
        Im = result['Im']
        r_Trelease_n_NMDA= result['difGluNMDA']
        r_Trelease_n_AMPA= result['difGluAMPA']
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

        fig = make_subplots(rows=2, cols=3, subplot_titles=('current [nA]', 'voltage [mV]', 'AMPA', 'NMDA', 'Glu difWave through AMPA', 'Glu difWave through NMDA'))
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
        row=1, col=3
        )
        fig.add_trace(
        go.Scatter(x=time, y=NMDAcurr_sum),
        row=2, col=1
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
        list_freq.append(freq)

def create_cell(model):
    cell = model()
    for dend in cell.dendrites:
        cell.create_synapse(dend, "AMPA")
        cell.create_synapse(dend, "NMDA")
    return cell
grc=GranuleCell

def experiment1(currents, frequencies, cell):
    for i in np.arange(0, 24, 2):
        inj = i
        currents.append(inj *1e-3)
        model = create_cell(cell)
        injectDuration = 800
        time, results = current_injection(   
            grc, 
            durationPre=100,
            durationInject=injectDuration,
            durationPost=100,
            inj_current= inj* 1e-3,
            current_for_holding = -0.00988, #0.01 WT:25 KO:46.25
            temperature=32,
            v_init=-80,
            dt = 0.025,
        )
        analysis(results, time, injectDuration, frequencies)

currents_gc = list()

freq_wt= list ()

experiment1(currents_gc, freq_wt, grc)

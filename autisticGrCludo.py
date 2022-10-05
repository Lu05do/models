import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz
from patch import p
import arborize.core
from arborize.core import flatten_composite
from dbbs_models import GranuleCell
import copy
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import scipy.optimize

class AutisticGranuleCell(GranuleCell):
    section_types = copy.deepcopy(GranuleCell.section_types)
    synapse_types = copy.deepcopy(GranuleCell.synapse_types)

AutisticGranuleCell.section_types['axon_initial_segment']= flatten_composite(AutisticGranuleCell, AutisticGranuleCell.section_types['axon_initial_segment'])
AutisticGranuleCell.section_types['axon_initial_segment']['mechanisms']['Na','granule_cell_FHF']['gnabar'] *= 1.9
#AutisticGranuleCell.section_types['axon_initial_segment']['attributes']['Kv3_4'][gkbar] *= 1.81  # Potassium [1] - transient
AutisticGranuleCell.section_types['soma']['mechanisms']['Kv4_3']['gkbar'] *= 1.8 # Potassium [2] - persistent
#AutisticGranuleCell.section_types['soma']['mechanisms']['Kv4_3']['gkbar'] *= 1.8 # Potassium [2] - persistent
AutisticGranuleCell.section_types['dendrites']['mechanisms']['Kca1_1']['gbar'] *= 1.8 # Potassium [3] - persistent

AutisticGranuleCell.synapse_types["NMDA"]["point_process"] = ("NMDA", "autistic")
AutisticGranuleCell.synapse_types["NMDA"]["attributes"]["difwave_init"] = 1  #tonic glu on NMDA
AutisticGranuleCell.synapse_types["NMDA"]["attributes"]["gmax"] *= 2.5
AutisticGranuleCell.synapse_types["AMPA"]["point_process"] = ("AMPA", "autistic")
AutisticGranuleCell.synapse_types["AMPA"]["attributes"]["difwave_init"] = 1  #tonic glu on AMPA

def change_soma(model):
     model.soma[0].diam = 5.3 

AutisticGranuleCell.morphologies[0] = (GranuleCell.morphologies[0], change_soma)

gc=GranuleCell()
gc_autistic=AutisticGranuleCell()

nDend=4
for dend in range(nDend):
    syn1 = gc.create_synapse(gc.dendrites[dend], "AMPA")
    syn2 = gc.create_synapse(gc.dendrites[dend], "NMDA")
for dend in range(nDend):
    syn1 = gc_autistic.create_synapse(gc_autistic.dendrites[dend], "AMPA")
    syn2 = gc_autistic.create_synapse(gc_autistic.dendrites[dend], "NMDA")


nDendStim=2
r_Trelease_n_NMDA= []
r_Trelease_n_AMPA= []
NMDAcurr=[]
AMPAcurr=[]

def current_inj(*models, nDend=4, durationTot=1000, durationForIstep=800,temperature=32, v_init=-80, dt = 0.025, inj_current=0, current_for_holding=0):
    p.time
    p.celsius = temperature
    p.v_init = v_init
    p.dt = dt
    for model in models:
        model.record_soma()
        model._spike_detector = p.NetCon(model.soma[0], None)
        clamp = model.soma[0].iclamp(x=0.5, delay=0, duration=durationTot, amplitude=current_for_holding)
        clamp = model.soma[0].iclamp(x=0.5, delay=100, duration=durationForIstep, amplitude=inj_current)
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
    p.continuerun(durationTot)

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

def voltclamp(*models, durationTot=1000, durationVstep=400, temperature=32, v_init=-80, dt = 0.025, vclamped=0, delayStep=100):
    p.time
    p.celsius = temperature
    p.v_init = v_init
    p.dt = dt
    for model in models:
        model.record_soma()
        model._spike_detector = p.NetCon(model.soma[0], None)
        clamp = model.soma[0].vclamp(x=0.5, delay=delayStep, duration=durationVstep, after=0, voltage=vclamped, holding=v_init)
        i = p.record(clamp._ref_i)
        model.spike_times = p.Vector()
        model._spike_detector.record(model.spike_times)
        for dend in range(nDend):
            AMPAsyn=model.dendrites[dend]._synapses[0]
            NMDAsyn=model.dendrites[dend]._synapses[1]
            #NMDAsyn.stimulate(start=100, number=4, interval=200)
            NMDAsyn.record()
            AMPAsyn.record()
            r_tr_NMDA = p.record(NMDAsyn._point_process._ref_Trelease)
            r_tr_AMPA = p.record(AMPAsyn._point_process._ref_Trelease)
            r_Trelease_n_NMDA.append(r_tr_NMDA)
            r_Trelease_n_AMPA.append(r_tr_AMPA)


        #v = p.record(clamp._ref_v)
    p.finitialize(v_init)
    p.continuerun(durationTot)

    return list(p.time), [{
        "Vm": list(model.Vm),
        "spikes": list(model.spike_times),
        "model": model.__class__.__name__,
        "Im": list(i),
        "difGluAMPA": list(r_Trelease_n_AMPA),
        "difGluNMDA": list(r_Trelease_n_NMDA)
    } for model in models]


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
#current clamp paper Soda et al fig 1B
def analysisSyn(*models, time=0, V=0):
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
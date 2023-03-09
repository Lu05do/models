import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz
from patch import p
import arborize.core
from dbbs_models import GranuleCell
import copy
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from sklearn.linear_model import LinearRegression
import scipy.optimize
from arborize.core import flatten_composite


def create_autistic_model():
    # Create template
    class AutisticGranuleCell(GranuleCell):
        morphologies = copy.deepcopy(GranuleCell.morphologies)
        section_types = copy.deepcopy(GranuleCell.section_types)
        synapse_types = copy.deepcopy(GranuleCell.synapse_types)

    alter_currents(AutisticGranuleCell)
    alter_synapses(AutisticGranuleCell)
    alter_soma(AutisticGranuleCell)

    return AutisticGranuleCell

def alter_currents(model):
    stypes = model.section_types
    stypes['axon_initial_segment'] = flatten_composite(model, stypes['axon_initial_segment'])
    # Increase AIS sodium current
    stypes['axon_initial_segment']['mechanisms']['Na','granule_cell_FHF']['gnabar'] *= 1.9
    # Increase AIS transient pota\ssium current
    stypes['axon_initial_segment']['mechanisms']['Kv3_4']['gkbar'] *= 1.81
    # Increase soma persistent potassium current
    stypes['soma']['mechanisms']['Kv4_3']['gkbar'] *= 1.8
    # Increase dendrites persistent potassium current
    stypes['dendrites']['mechanisms']['Kca1_1']['gbar'] *= 1.8

def alter_synapses(model):
    syntypes = model.synapse_types
    syntypes["NMDA"]["point_process"] = ("NMDA", "autistic")
    # Tonic glutamate on NMDA
    syntypes["NMDA"]["attributes"]["difwave_init"] = 0.5 #not considered, to be modified in mod
    syntypes["NMDA"]["attributes"]["gmax"] *= 2.5
    # Tonic glutamate on AMPA
    syntypes["AMPA"]["point_process"] = ("AMPA", "autistic")
    syntypes["AMPA"]["attributes"]["difwave_init"] = 0.5 #not considered, to be modified in mod

def alter_soma(model):
    def change_soma(model):
        # Reduce cell size (Length, diameter -- > increase membrane Resistance)
        model.soma[0].diam = 5.3

    if isinstance(model.morphologies[0], str):
        model.morphologies[0] = (model.morphologies[0], change_soma)

def add_gabazine(model):
    model.section_types["dendrites"]["mechanisms"].pop(("Leak","GABA"))

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

    fig = make_subplots(rows=2, cols=3,subplot_titles=('current [nA]', 'voltage [mV]', 'AMPA', 'NMDA','Glu difWave through AMPA', 'Glu difWave through NMDA' ))
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

    fig = make_subplots(rows=3, cols=3,subplot_titles=( 'voltage [mV]', 'AMPA curr [nA]','NMDA curr [nA]','  ', 'Glu difWave through AMPA', 'Glu difWave through NMDA' ) )
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

AutisticGranuleCell = create_autistic_model()

#gc=GranuleCell()
gc_autistic=AutisticGranuleCell()
for dend in gc_autistic.dendrites:
    gc_autistic.create_synapse(dend, "AMPA")
    gc_autistic.create_synapse(dend, "NMDA")

r_Trelease_n_NMDA= []
r_Trelease_n_AMPA= []
NMDAcurr=[]
AMPAcurr=[]

nDend=4
durationForIstep=800
inj_currents=np.arange(-2, 22, 2) #pA

time, V = synapticStim(gc_autistic, nDend=4, nDendStim=1, durationTot=500, temperature=32, v_init=-65, dt = 0.025, startTrain=50, nPreSynSpikes=1, ISI=10)  # NEURON default is Nano
analysisSyn(gc_autistic, time=time, V=V)

freqWT = np.empty(len(inj_currents), dtype=object)
for count, inj_current in enumerate(inj_currents):
     gc=GranuleCell()
     for dend in range(nDend):
         syn1 = gc.create_synapse(gc.dendrites[dend], "AMPA")
         syn2 = gc.create_synapse(gc.dendrites[dend], "NMDA")
     r_Trelease_n_NMDA= []
     r_Trelease_n_AMPA= []

     time, V = current_inj(gc, nDend=4, durationTot=1000, durationForIstep=durationForIstep,temperature=32, v_init=-80, dt = 0.025, inj_current=inj_current * 1e-3, current_for_holding = -13 * 1e-3)
     spikes = V[0]['spikes']
     if len(spikes) > 1:
         freqWT[count]= len(spikes)/durationForIstep * 1e3
     else:
         freqWT[count]=0
     del spikes
     print(count, inj_current, freqWT)
     analysis(gc, durForFreqCalculation=800, time=time, V=V)
     del gc
lm = LinearRegression().fit(inj_currents.reshape((-1, 1)), freqWT)
# print(lm.intercept_)
# print(lm.coef_)

freqKO = np.empty(len(inj_currents), dtype=object)
for count, inj_current in enumerate(inj_currents):
     gc_autistic=AutisticGranuleCell()
     for dend in range(nDend):
         syn1 = gc_autistic.create_synapse(gc_autistic.dendrites[dend], "AMPA")
         syn2 = gc_autistic.create_synapse(gc_autistic.dendrites[dend], "NMDA")
     r_Trelease_n_NMDA= []
     r_Trelease_n_AMPA= []

     time, V = current_inj(gc_autistic, durationTot=1000, durationForIstep=durationForIstep,temperature=32, v_init=-80, dt = 0.025, inj_current=inj_current * 1e-3, current_for_holding = -13 * 1e-3)
     spikes = V[0]['spikes']
     if len(spikes) > 1:
         freqKO[count]= len(spikes)/durationForIstep * 1e3
     else:
         freqKO[count]=0
     del spikes
     print(count, inj_current, freqKO)
     analysis(gc_autistic, durForFreqCalculation=800, time=time, V=V)
     del gc_autistic
lm = LinearRegression().fit(inj_currents.reshape((-1, 1)), freqKO)
# print(lm.intercept_)
# print("slope: ", lm.coef_)
# print("x-intercept: ", -lm.intercept_/lm.coef_)

fig = make_subplots(rows=1, cols=1)
fig.add_trace(
 go.Scatter(x=inj_currents, y=freqWT, name='WT',
          mode='lines+markers',
         line=go.scatter.Line(color="black"),
         marker=dict(
             color='black',
             size=2)), row=1, col=1
 )
fig.add_trace(
 go.Scatter(x=inj_currents, y=freqKO, name='KO',
         mode='lines+markers',
         line=go.scatter.Line(color="blue"),
         marker=dict(
             color='blue',
             size=2)),row=1, col=1
 )
fig['layout']['xaxis']['title']='Inj current [pA]'
fig['layout']['yaxis']['title']='Spike freq [Hz]'
fig.show()
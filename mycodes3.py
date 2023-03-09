import matplotlib.pyplot as plt
#from arborize.core import flatten_composite
from arborize import define_model
from arborize import bsb_schematic
from bsb.morphologies import Morphology
from arborize import neuron_build
import numpy as np
from numpy import trapz
from patch import p
#from dbbs_models import GranuleCell, AutisticGranuleCell
import copy
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import scipy.optimize


definition = define_model({
    "morphologies" : ["GranuleCell.swc"],
    "synapses" : {
        "AMPA": {
            "point_process": ("AMPA", "granule"),
            "attributes": {
                "tau_facil": 5,
                "tau_rec": 8,
                "tau_1": 1,
                "gmax": 1400,
                "U": 0.43,
            },
        },
        "NMDA": {
            "point_process": ("NMDA", "granule"),
            "attributes": {
                "tau_facil": 5,
                "tau_rec": 8,
                "tau_1": 1,
                "gmax": 23500,
                "U": 0.43,
            },
        },
        "GABA": {"point_process": ("GABA", "granule"), "attributes": {"U": 0.35}}
    },
  "cable_types": {
    
        "soma": {
            "cable": {"Ra": 100, "cm": 2},
            "ions": {"k": {"e": -80.993}, "ca": {"e": 137.5}},
            "mechanisms": {
                "Leak": {"e": -60, "gmax": 0.00029038073716},
                "Kv3_4": {"gkbar": 0.00076192450952},
                "Kv4_3": {"gkbar": 0.00281496839061},
                "Kir2_3": {"gkbar": 0.00074725514702},
                "Ca": {"gcabar": 0.00060938071784},
                "Kv1_1": {"gbar": 0.00569738264555},
                "Kv1_5": {"gKur": 0.00083407556714},
                "Kv2_2": {"gKv2_2bar": 1.203410852e-05},
                ("cdp5", "CR"): {},
            },
        },
        "dendrites": {
            "cable": {"Ra": 100, "cm": 2.5},
            "ions": {"k": {"e": -80.993}, "ca": {"e": 137.5}},
            "mechanisms": {
                "Leak": {"e": -60, "gmax": 0.00025029700737},
                ("Leak", "GABA"): {},
                "Ca": {"gcabar": 0.00500128008459},
                "Kca1_1": {"gbar": 0.01001807454651},
                "Kv1_1": {"gbar": 0.00381819207934},
                ("cdp5", "CR"): {},
            },
            "synapse_type": ["AMPA", "NMDA", "GABA"]
        },
        "axon": {
            "cable": {"Ra": 100, "cm": 1},
            "ions": {"na": {"e": 87.39}, "k": {"e": -80.993}, "ca": {"e": 137.5}},
            "mechanisms": {
                ("Na", "granule_cell"): {"gnabar": 0.02630163681502},
                "Kv3_4": {"gkbar": 0.00237386061632},
                "Leak": {"e": -60, "gmax": 9.364092125e-05},
                "Ca": {"gcabar": 0.00068197420273},
                ("cdp5", "CR"): {},
            },
        },
        "ascending_axon": {
            "cable": {"Ra": 100, "cm": 1},
            "ions": {"na": {"e": 87.39}, "k": {"e": -80.993}, "ca": {"e": 137.5}},
            "mechanisms": {
                ("Na", "granule_cell"): {"gnabar": 0.02630163681502},
                "Kv3_4": {"gkbar": 0.00237386061632},
                "Leak": {"e": -60, "gmax": 9.364092125e-05},
                "Ca": {"gcabar": 0.00068197420273},
                ("cdp5", "CR"): {},
            },
        },
        "axon_initial_segment": {
            "cable": {"Ra": 100, "cm": 1},
                "ions": {"na": {"e": 87.39}, "k": {"e": -80.993}, "ca": {"e": 137.5}},
                "mechanisms": {
                    ("Na", "granule_cell_FHF"): {"gnabar": 1.28725006737226},
                    "Kv3_4": {"gkbar": 0.00649595340654},
                    "Leak": {"e": -60, "gmax": 0.00029276697557},
                    "Ca": {"gcabar": 0.00031198539472},
                    "Km": {"gkbar": 0.00056671971737},
                    ("cdp5", "CR"): {},
                },
        },
        "axon_hillock": {
            "cable": {"Ra": 100, "cm": 2},
                "ions": {"na": {"e": 87.39}, "k": {"e": -80.993}, "ca": {"e": 137.5}},
                "mechanisms": {
                    "Leak": {"e": -60, "gmax": 0.0003695818972},
                    ("Na", "granule_cell_FHF"): {"gnabar": 0.00928805851462},
                    "Kv3_4": {"gkbar": 0.02037346310915},
                    "Ca": {"gcabar": 0.00057726155447},
                    ("cdp5", "CR"): {},
                },
        },
        "parallel_fiber": {
            "cable": {"Ra": 100, "cm": 1},
            "ions": {"na": {"e": 87.39}, "k": {"e": -80.993}, "ca": {"e": 137.5}},
            "mechanisms": {
                ("Na", "granule_cell"): {"gnabar": 0.01771848449261},
                "Kv3_4": {"gkbar": 0.00817568047037},
                "Leak": {"e": -60, "gmax": 3.5301616e-07},
                "Ca": {"gcabar": 0.0002085683353},
                ("cdp5", "CR"): {},
            },
        },
    },
})

print(type(definition))
tags = {
        16: ["axon", "axon_hillock"],
        17: ["axon", "axon_initial_segment"],
        18: ["axon", "ascending_axon"],
        19: ["axon", "parallel_fiber"],
    }

'''grc= GranuleCell()
for dend in grc.dendrites:
    grc.create_synapse(dend, "AMPA")
    grc.create_synapse(dend, "NMDA")'''

r_Trelease_n_NMDA= []
r_Trelease_n_AMPA= []
NMDAcurr=[]
AMPAcurr=[]

#def add_gabazine(AutisticGranuleCell):
#    AutisticGranuleCell.section_types["dendrites"]["mechanisms"].pop(("Leak","GABA"))

def current_injection(*models, nDend=4, durationPre=100, durationInject=800, durationPost=100, temperature=32, v_init=-80, dt = 0.025, inj_current=4, current_for_holding=0):
    p.time 
    p.celsius = temperature
    p.v_init = v_init
    p.dt = dt
    for model in models:
        #model.record_soma()
        model.soma[0].record()
        model._spike_detector = p.NetCon(model.soma[0], None)
        clamp = model.soma[0].iclamp(x=0.5, delay=durationPre, duration=durationInject, amplitude=inj_current)
        model.soma[0].iclamp(x=0.5, delay=0, duration=durationPre+durationInject+durationPost, amplitude=current_for_holding)
        i = p.record(clamp._ref_i)
        model.spike_times = p.Vector()
        model._spike_detector.record(model.spike_times)
        for dend in range(nDend):
            print(dir(model.dendrites[dend].synapse_type))
            print(len(model.dendrites[dend].synapse_type))
            AMPAsyn=model.dendrites[dend].synapse[0]
            NMDAsyn=model.dendrites[dend].synapse[1]
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

def voltage_clamp(*models, durationTot=1000, durationVstep=400, temperature=32, v_init=-80, dt = 0.025, vclamped=0, delayStep=100):
    r_Trelease_n_NMDA= []
    r_Trelease_n_AMPA= []
    
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

def synapticStim(*models, nDend=4, nDendStim=2, durationTot=500, temperature=32, v_init=-65, dt = 0.025, startTrain=0, nPreSynSpikes=0, ISI=0):

    NMDAcurr=[]
    AMPAcurr=[]
    r_Trelease_n_NMDA= []
    r_Trelease_n_AMPA= []

    p.time
    p.celsius = temperature
    p.v_init = v_init
    p.dt = dt
    for model in models:
        model.record_soma()
        model._spike_detector = p.NetCon(model.soma[0], None)
        model.spike_times = p.Vector()
        model._spike_detector.record(model.spike_times)
        clamp = model.soma[0].vclamp(x=0.5, delay=0, duration=durationTot, after=0, voltage=v_init, holding=v_init)
        for dend in range(nDend):
            AMPAsyn=model.dendrites[dend]._synapses[0]
            NMDAsyn=model.dendrites[dend]._synapses[1]
            NMDAcurr.append(NMDAsyn.record())
            AMPAcurr.append(AMPAsyn.record())
            r_tr_NMDA = p.record(NMDAsyn._point_process._ref_Trelease)
            r_tr_AMPA = p.record(AMPAsyn._point_process._ref_Trelease)
            r_Trelease_n_NMDA.append(r_tr_NMDA)
            r_Trelease_n_AMPA.append(r_tr_AMPA)
        for dend in range(nDendStim):
            NMDAsyn.stimulate(start=startTrain, number=nPreSynSpikes, interval=ISI)
            AMPAsyn.stimulate(start=startTrain, number=nPreSynSpikes, interval=ISI)
        #v = p.record(clamp._ref_v)
    p.finitialize(v_init)
    p.continuerun(durationTot)

    return list(p.time), [{
        "Vm": list(model.Vm),
        "spikes": list(model.spike_times),
        "model": model.__class__.__name__,
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

def analysisSyn(results, time=0, V=0):
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

def create_cell(model):
    morpho = Morphology.from_swc("GranuleCell.swc", tags=tags)
    branch_ids = [morpho.branches.index(b) for b in morpho.subtree(["dendrites"]).branches]

    schematic = bsb_schematic(morpho, definition)
    cell = neuron_build(schematic)

    for dend_id in branch_ids:
        cell.insert_synapse("AMPA", (dend_id, 0))
        #cell.synapse(dend, "NMDA")
    return cell

'''def create_cell(model):
    cell = model()
    for dend in cell.dendrites:
        cell.create_synapse(dend, "AMPA")
        cell.create_synapse(dend, "NMDA")
    return cell'''


def experimentwt(currents, frequencies, cell):
    for i in np.arange(0, 24, 2):
        inj = i
        currents.append(inj *1e-3)
        model = create_cell(cell)
        injectDuration = 800
        time, results = current_injection(   
            model, 
            durationPre=100,
            durationInject=injectDuration,
            durationPost=100,
            inj_current= inj* 1e-3,
            current_for_holding = -0.00988, #-0.00988, #0.01 WT:25 KO:46.25
            temperature=32,
            v_init=-80,
            dt = 0.025,
        )
        analysis(results, time, injectDuration, frequencies)

currents_gc = list()
freq_wt= list ()
experimentwt(currents_gc, freq_wt, definition)

def experimentaut(currents, frequencies, cell):
    for i in np.arange(0, 24, 2):
        inj = i
        currents.append(inj *1e-3)
        model = create_cell(cell)
        model.soma[0].diam = 5.3
        injectDuration = 800
        time, results = current_injection(   
            model, 
            durationPre=100,
            durationInject=injectDuration,
            durationPost=100,
            inj_current= inj* 1e-3,
            current_for_holding = -0.00988, #-0.00988, #0.01 WT:25 KO:46.25
            temperature=32,
            v_init=-80,
            dt = 0.025,
        )
        analysis(results, time, injectDuration, frequencies)
currents_aut = list()
freq_aut = list()
experimentaut(currents_aut, freq_aut, AutisticGranuleCell)

#time, V = synapticStim(aut_gc, nDend=4, nDendStim=2, durationTot=1000, temperature=32, v_init=-65, dt = 0.025, startTrain=50, nPreSynSpikes=1, ISI=10)  # NEURON default is Nano
#analysisSyn(aut_gc, time=time, V=V)

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
fig['layout']['xaxis']['title']='Inj current [nA]'
fig['layout']['yaxis']['title']='Spike freq [Hz]'
fig.show()
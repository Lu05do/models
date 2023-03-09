import matplotlib.pyplot as plt
#from arborize.core import flatten_composite
from arborize import define_model
from arborize import bsb_schematic
from bsb.morphologies import Morphology
from arborize import neuron_build
import numpy as np
from numpy import trapz
from patch import p
#from dbbs_models import PurkinjeCell
import copy
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import scipy.optimize
from modelforcodes.GranuleWT import definitionGrCWT, tagsGrC
from modelforcodes.GranuleKO import definitionGrCKO, tagsGrC


r_Trelease_n_NMDA= []
r_Trelease_n_AMPA= []
NMDAcurr=[]
AMPAcurr=[]

def record_spikes(model):
    model._spike_detector = p.NetCon(model.soma[0], None)
    model.spike_times = p.Vector()
    model._spike_detector.record(model.spike_times)

#def add_gabazine(AutisticGranuleCell):
#    AutisticGranuleCell.section_types["dendrites"]["mechanisms"].pop(("Leak","GABA"))


def voltage_clamp(
    *models,
    nDend = 1,
    durationVstep=400,
    vclamped=0,
    delayStep=100,
    dt=0.025,
    durationTot=1000,
    temperature=32,
    v_init=-80,
):
    t = p.time
    p.celsius = temperature
    p.v_init = v_init
    p.dt = dt
    for model in models:
        model.record_soma()
        record_spikes(model)

        clamp = model.soma[0].vclamp(
            x=0.5,
            before=delayStep,
            duration=durationVstep,
            after=0,
            voltage=vclamped,
            holding=v_init
        )
        i = p.record(clamp._ref_i)
        dendrites = []
        for dend in model.dendrites:
            dend_synapses = []
            dendrites.append(dend_synapses)
            for i, synapse in enumerate(dend.synapses):
                synapse.stimulate(start=100, number=1, interval=200)
                i_syn = synapse.record()
                dend_synapses.append(i_syn)
                print("Index", i, "synapse", synapse, dir(synapse))
            # r_tr_NMDA = p.record(NMDAsyn._pp.Trelease)
            # r_tr_AMPA = p.record(AMPAsyn._pp.Trelease)
            # r_Trelease_n_NMDA.append(r_tr_NMDA)
            # r_Trelease_n_AMPA.append(r_tr_AMPA)

        #v = p.record(clamp._ref_v)
    p.finitialize(v_init)
    p.continuerun(durationTot)

    return list(t), [{
        "Vm": list(model.Vm),
        "spikes": list(model.spike_times),
        "model": model.__class__.__name__,
        "synaptic_currents": [
            [list(syn) for syn in dend]
            for dend in dendrites
        ],
    } for model in models]

def analysisVC(time, results):
    for result in results:
        fig = make_subplots(rows = 2, cols = 1, subplot_titles=('voltage [mV]', 'current [nA]'))
        fig.add_trace(
            go.Scatter(
                x=time,y=result["Vm"],
                name=result["model"] + " Vm"
            ),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(
                x =time, y = result["synaptic_currents"][0][0],
                name=result["model"] + " synapse"
            ),
            row=2, col=1      
        )
        fig.update_layout(xaxis_title="time [ms]")
        fig.show()

def create_cell(definition, synapses=["AMPA", "NMDA", "GABA"]):
    morpho = Morphology.from_swc("GranuleCell.swc", tags=tagsGrC)
    branch_ids = [morpho.branches.index(b) for b in morpho.subtree(["dendrites"]).branches]

    schematic = bsb_schematic(morpho, definition)
    cell = neuron_build(schematic)

    for dend_id in branch_ids:
        for synapse in synapses:
            cell.insert_synapse(synapse, (dend_id, 0))

    return cell


def experimentVC():
    model_ampa = create_cell(definitionGrCWT, ["AMPA"])
    #model_gaba = create_cell(definitionGrCWT, ["GABA"])
    #model_nmda = create_cell(definitionGrCWT, ["NMDA"])
    model_ampaKO = create_cell(definitionGrCKO, ["AMPA"])
    #model_nmda = create_cell(definitionGrCKO, ["NMDA"])
    time, results = voltage_clamp(
        # model_gaba,
        # model_nmda,
        model_ampa,
        model_ampaKO,
        durationTot=400, 
        durationVstep=400, 
        temperature=32, 
        v_init=-70, 
        dt = 0.025, 
        vclamped=-60,
        delayStep=100,
    )
    analysisVC(time, results)

experimentVC()

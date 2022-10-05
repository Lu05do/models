from arborize.core import flatten_composite
from dbbs_models import GranuleCell
import copy
from patch import p

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
    # Increase AIS transient potassium current
    stypes['axon_initial_segment']['mechanisms']['Kv3_4']['gkbar'] *= 1.81
    # Increase soma persistent potassium current
    stypes['soma']['mechanisms']['Kv4_3']['gkbar'] *= 1.8
    # Increase dendrites persistent potassium current
    stypes['dendrites']['mechanisms']['Kca1_1']['gbar'] *= 1.8

def alter_synapses(model):
    syntypes = model.synapse_types
    syntypes["NMDA"]["point_process"] = ("NMDA", "autistic")
    # Tonic glutamate on NMDA
    syntypes["NMDA"]["attributes"]["difwave_init"] = 1
    syntypes["NMDA"]["attributes"]["gmax"] *= 2.5
    # Tonic glutamate on AMPA
    syntypes["AMPA"]["point_process"] = ("AMPA", "autistic")
    syntypes["AMPA"]["attributes"]["difwave_init"] = 1

def alter_soma(model):
    def change_soma(model):
        # Reduce cell size (Length, diameter -- > increase membrane Resistance)
        model.soma[0].diam = 5.3
    
    if isinstance(model.morphologies[0], str):
        model.morphologies[0] = (model.morphologies[0], change_soma)

def add_gabazine(model):
    model.section_types["dendrites"]["mechanisms"].pop(("Leak","GABA"))


def current_injection(
    *models,
    nDend=4,
    durationPre=200,
    durationInject=800,
    durationPost=200,
    temperature=32,
    v_init=-80,
    dt = 0.025,
    inj_current=4,
    current_for_holding=0
):

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
        model.soma[0].iclamp(x=0.5, delay=0, duration=durationPre, amplitude=current_for_holding)
        clamp = model.soma[0].iclamp(x=0.5, delay=durationPre, duration=durationInject, amplitude=inj_current)
        model.soma[0].iclamp(x=0.5, delay=durationPre + durationInject, duration=durationPost, amplitude=current_for_holding)
        i = p.record(clamp._ref_i)
        model.spike_times = p.Vector()
        model._spike_detector.record(model.spike_times)
        for dend in range(nDend):
            AMPAsyn=model.dendrites[dend]._synapses[0]
            NMDAsyn=model.dendrites[dend]._synapses[1]
            NMDAcurr.append(NMDAsyn.record())
            AMPAcurr.append(AMPAsyn.record())
            r_tr_NMDA = p.record(NMDAsyn._point_process._ref_Trelease)
            r_tr_AMPA = p.record(AMPAsyn._point_process._ref_Trelease)
            r_Trelease_n_NMDA.append(r_tr_NMDA)
            r_Trelease_n_AMPA.append(r_tr_AMPA)

    p.finitialize(v_init)
    p.continuerun(durationPre + durationInject + durationPost)

    return list(p.time), [
        {
            "Vm": list(model.Vm),
            "spikes": list(model.spike_times),
            "model": model.__class__.__name__,
            "Im": list(i),
            "NMDAcurr": list(NMDAcurr),
            "AMPAcurr": list(AMPAcurr),
            "difGluAMPA": list(r_Trelease_n_AMPA),
            "difGluNMDA": list(r_Trelease_n_NMDA),
        }
        for model in models
    ]


def synapticStim(
    *models,
    durationTot=500,
    temperature=32,
    v_init=-65,
    dt = 0.025,
    startTrain=0,
    nPreSynSpikes=0,
    ISI=0,
    nDendStim=2,
    nDend=4
):

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
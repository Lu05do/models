from patch import p
import plotly.graph_objs as go
import glia
from arborize import define_model
from arborize import bsb_schematic
from arborize import neuron_build
from bsb.morphologies import Morphology
from modelforcodes.Basket import definitionBasket, tagsBk

morpho = Morphology.from_swc("BasketCell.swc", tags=tagsBk)
branch_ids = [morpho.branches.index(b) for b in morpho.subtree(["dendrites"]).branches]
print(branch_ids)

schematic = bsb_schematic(morpho, definitionBasket)
cell = neuron_build(schematic) 
for dend_id in branch_ids:
        cell.insert_synapse("AMPA", (dend_id, 0))
        #cell.insert_synapse("AMPA_AA", (dend_id, 0))
        #cell.insert_synapse("AMPA_MF", (dend_id, 0))
        cell.insert_synapse("NMDA", (dend_id, 0))
        cell.insert_synapse("GABA", (dend_id, 0))

print(dir(cell))

'''p.celsius = 32
s = p.Section()
t = p.time
pp = glia.insert(s, 'AMPA', "autistic")
pp.stimulate(start=0, interval=10, number=10)
r = pp.record()
p.run(100)
go.Figure(go.Scatter(x=list(t), y=list(r))).show()'''


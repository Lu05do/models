from dbbs_models import GranuleCell, quick_test, quick_plot
from neuron import *
import numpy as np
from patch import p
import pickle
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.graph_objects as go

gcs = []
gc = GranuleCell()

for i in np.arange(10, 24, 2):
    inj = i
    gc.soma[0].iclamp(delay = 100, duration= 800, amplitude= inj*1e-3)
    gc.soma[0].iclamp(delay = 0, duration= 800, amplitude= -0.003)
    quick_plot(gc, v_init= -80, duration=1000)


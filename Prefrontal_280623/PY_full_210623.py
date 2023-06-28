#!/usr/bin/env python
# coding: utf-8

# ## Required packages and .py

# In[1]:

from Dendritic_eqs_210623 import *
from Soma_eqs_210623 import *


# In[2]:


from brian2 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


# ## Global parameters

#Number of pyramidal neurons desired
n_PY = 1

#Area
s_soma = 10**-6*cm**2

#somatic and dendritic compartments
R = 10 * Mohm


# ## Axosomatic and dendritic compartment initialization, gap junction and monitors

close('all')
start_scope()
prefs.codegen.target = "cython"
defaultclock.dt = 0.02*ms

#Instantiate neurons
PY_dendrite = NeuronGroup(n_PY,Dendritic_eqs,method='rk4',threshold='v>0*mV',refractory=3*ms)
PY_dendrite.v = -70*mvolt
PY_dendrite.m_na = 0.05
PY_dendrite.h_na = 0.95
PY_dendrite.m_nap = 0.00
PY_dendrite.m_km = 0.05
PY_dendrite.m_kca = 0.05
PY_dendrite.m_hva = 0.05
PY_dendrite.h_hva = 0.95
# PY_dendrite.CA_i = 0.0001*mM #not from the paper
PY_dendrite.CA_i = 3.5E-4*mM
PY_dendrite.g_na = 0.8*msiemens*cm**-2 
PY_dendrite.g_nap = 2.5*msiemens*cm**-2 
PY_dendrite.g_km = 0.02*msiemens*cm**-2
PY_dendrite.rho = 165
#PY_dendrite.rho = 50

PY_soma = NeuronGroup(n_PY,Soma_eqs,method='rk4',threshold='v>0*mV',refractory=3*ms)
PY_soma.h_na = 0.95
PY_soma.m_na = 0.05
PY_soma.m_nap = 0.00
PY_soma.n_k = 0.05
PY_soma.g_na = 3000*msiemens*cm**-2
PY_soma.g_k = 200*msiemens*cm**-2
PY_soma.g_nap = 15*msiemens*cm**-2

# inpext=-20*amp*meter**-2 
# PY_soma.Iext = inpext

#Defining the gap junction between the two compartments
eq_gap_dendrite ='''
    Igap_post = g * (v_post - v_pre) : amp * meter**-2 (summed)
    g : siemens * meter**-2
'''
eq_gap_soma ='''
    vgap_post = v_pre : volt (summed)
'''

#From soma to dendrites
gapPY_som_den = Synapses(PY_soma,PY_dendrite,model=eq_gap_dendrite)
gapPY_som_den.connect(j = 'i')
gapPY_som_den.g = 1 / (R * s_soma * gapPY_som_den.rho) 

#From dendrites to soma
gapPY_den_som = Synapses(PY_dendrite,PY_soma,model=eq_gap_soma)
gapPY_den_som.connect(j = 'i')

#Monitors
V1 = StateMonitor(PY_soma,('v','vgap','Iext','m_na','h_na','m_nap','n_k','alpham_na','alphah_na','betam_na','betah_na','tau_h_na','tau_m_na','m_nainf','h_nainf'),record=True)
V2 = StateMonitor(PY_dendrite,('v','Igap','m_na','h_na','m_nap','m_km','m_kca','m_hva','h_hva','alpham_na','alphah_na','betam_na','betah_na','tau_h_na','tau_m_na','m_nainf','h_nainf','CA_i'),record=True)

#Patch clamp test
runtime = 3100*ms
num_samples = int(runtime/defaultclock.dt)
init_arr = zeros(num_samples)
init_arr[0]=1
init_timedarray = TimedArray(init_arr, dt=defaultclock.dt)

PY_dendrite.Iext = 0*pA*meter**-2
run(50*ms,report='text',report_period=60*second)
store('initial')
# PY_dendrite.Iext = 0.025*amp*meter**-2
PY_dendrite.Iext = 0.05*amp*meter**-2
run(3000*ms,report='text',report_period=60*second)
PY_dendrite.Iext = 0*pA*meter**-2
run(50*ms,report='text',report_period=60*second)


# In[5]:


#Visualize

fig,ax = subplots(5,1, sharex = True)
ax[0].plot(V1.t/ms, V1.v[0]/mV, label='axosomatic')
ax[0].plot(V2.t/ms, V2.v[0]/mV, label='dendrites')
ax[0].legend()
ax[0].set_ylabel ('Voltage (mV)')
ax[1].plot(V1.t/ms, V1.vgap[0]/mV, label='vgap (soma)')
ax[1].legend()
ax[1].set_ylabel ('Voltage (mV)')
ax[2].plot(V2.t/ms, (V2.Igap[0]/pA)*-1, label='- Igap (dendrites)')
ax[2].legend()
ax[2].set_ylabel ('Electric current (pA)')
ax[2].set_xlabel ('Time (ms)')
ax[3].plot(V1.t/ms, V2.m_na[0], label='m_na')
ax[3].plot(V1.t/ms, V2.h_na[0], label='h_na')
#ax[3].plot(V1.t/ms, V1.m_nap[0], label='m_nap')
ax[3].plot(V1.t/ms, V2.m_km[0], label='n_k')
ax[3].legend()
ax[3].set_ylabel ('Activation')
ax[3].set_xlabel ('Time (ms)')
ax[4].plot(V1.t/ms, V2.betam_na[0], label='betam_na')
ax[4].plot(V1.t/ms, V2.alpham_na[0], label='alpham_na')
ax[4].plot(V1.t/ms, V2.betah_na[0], label='betah_na')
ax[4].plot(V1.t/ms, V2.alphah_na[0], label='alphah_na')
ax[4].legend()
ax[4].set_ylabel ('Activation')
ax[4].set_xlabel ('Time (ms)')

figure()
plot(V2.t/ms, V2.CA_i[0])
#!/usr/bin/env python
# coding: utf-8
# Thalamus_matrix.ipynb provides a function to create a proper Thalamic matrix/core with AMPA, NMDA, GABAa synapses and GABAb synapses. 
# ## Required packages and .ipynb

# In[1]:


from brian2 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


# In[ ]:


# get_ipython().run_line_magic('run', "'/home/mathilde/Modèles/Prefrontal_Wei/Cells_and_synapses/TC_eqs.ipynb'")
# get_ipython().run_line_magic('run', "'/home/mathilde/Modèles/Prefrontal_Wei/Cells_and_synapses/RE_eqs.ipynb'")
# get_ipython().run_line_magic('run', "'/home/mathilde/Modèles/Prefrontal_Wei/Cells_and_synapses/Synapses.ipynb'")
from TC_eqs_210623 import *
from RE_eqs_210623 import *
from Synapses_210623 import *


# ## Defining a cortical layer: instantiating neurons, synapses and Poisson inputs

# In[2]:


def create_thalamic_subparts(N):
    
    ###Preferences
    prefs.codegen.target = 'cython'
    defaultclock.dt = 0.02*ms
    
    ###Parameters
    #Channel-specific conductances per unit of surface
    g_syn_gabaa_retc = 0.15*usiemens
    g_syn_gabaa_rere = 0.125*usiemens
    g_syn_ampa_tcre = 0.35*usiemens
    #Reverse potentials
    E_ampa = 0*mV
    E_gabaa = -70*mV
    #Fraction of resources used per AP
    U_gabaa = 0.073
    U_ampa = 0.07
    #Rate constants
    alpha_ampa = 1.1*kHz
    beta_ampa = 0.19*kHz
    alpha_gabaa = 10.5*kHz
    beta_gabaa = 0.166*kHz
    #Number of neurons
    N_RE,N_TC = N,N
    #Fictional & arbitrary position of neurons for easyness of synapses.connect() condition
    neuron_spacing = 1*um
    #Areas of the different neurons
    s_TC = 1.43e-4*cm**2
    s_RE = 2.9E-4*cm**2
    
    ###Instantiate neurons
    #Reticular neurons
    RE = NeuronGroup(N,RE_eqs,method='exponential_euler',threshold='v>0*mV',refractory=3*ms,events={'custom_poisson_PY':'v>1000*mV','custom_poisson_IN':'v>1000*mV','custom_poisson_inter':'v>1000*mV'})
    RE.v = -70*mV
    RE.m_na = 0.01
    RE.h_na = 0.99
    RE.m_t = 0.01
    RE.h_t = 0.01
    RE.n_k = 0.01
    RE.CA_i_RE = 1e-4*mM
    RE.x = '(4*i+1.5)*neuron_spacing'
    
    #Thalamic relay cells
    TC = NeuronGroup(N,TC_eqs,method='exponential_euler',threshold='v>0*mV',refractory=3*ms,events={'custom_poisson_PY':'v>1000*mV','custom_poisson_IN':'v>1000*mV','custom_poisson_inter':'v>1000*mV'})
    TC.v = -70*mV
    TC.m_na = 0.01
    TC.h_na = 0.99
    TC.m_t = 0.01
    TC.h_t = 0.01
    TC.n_k = 0.01
    TC.CA_i_TC = 1e-4*mM
    TC.P1 = 0.0
    TC.Op = 0.5
    TC.Op_L = 0.0
    TC.x = '(4*i+1.5)*neuron_spacing'
    
    ###Instantiate synapses
    S_GABAA_RE_TC = syn_gabaa(RE,TC,'IsynGABAA_RE_TC',s_TC,'abs(i-j)<=6',g_syn_gabaa_retc,E_gabaa,alpha_gabaa,beta_gabaa,U_gabaa) #radii = 6
    S_GABAA_RE_TC.t_last_spike = -100*ms
    S_GABAA_RE_TC.D = 1
    S_GABAB_RE_TC = syn_gabab(RE,TC,'IsynGABAB_RE_TC',s_TC,'abs(i-j)<=6') 
    S_GABAB_RE_TC.t_last_spike = -100*ms
    S_GABAA_RE_RE = syn_gabaa(RE,RE,'IsynGABAA_RE_RE',s_RE,'abs(i-j)<=4 and i!=j',g_syn_gabaa_rere,E_gabaa,alpha_gabaa,beta_gabaa,U_gabaa) 
    S_GABAA_RE_RE.t_last_spike = -100*ms
    S_GABAA_RE_RE.D = 1
    S_AMPA_TC_RE = syn_ampa(TC,RE,'IsynAMPA_TC_RE',s_RE,'abs(i-j)<=4',g_syn_ampa_tcre,g_syn_ampa_tcre,E_ampa,alpha_ampa,beta_ampa,U_ampa,0,'IEPSPs_TC_RE',0) 
    S_AMPA_TC_RE.t_last_spike = -100*ms
    S_AMPA_TC_RE.D = 1
    S_AMPA_TC_RE.A_mEPSP = 0
    S_AMPA_TC_RE.w = 0
    
    ###Define monitors
    #general neuron monitoring
    V1=StateMonitor(RE,('v'),record=True)
    V2=StateMonitor(TC,('v'),record=True)
    #spike monitoring
    R1=SpikeMonitor(RE,record=True)
    R2=SpikeMonitor(TC,record=True)
    #synaptic currents monitoring
    I1=StateMonitor(RE,('Isyn_RE'),record=True)
    I2=StateMonitor(TC,('Isyn_TC'),record=True)
    #synapses monitoring (optionnal)
    #S1=StateMonitor(S_GABAA_RE_TC,('W','g_syn'),record=True)
    #S2=StateMonitor(S_GABAB_RE_TC,('R_gab','G'),record=True)
    #S3=StateMonitor(S_GABAA_RE_RE,('W'),record=True)
    #S4=StateMonitor(S_AMPA_TC_RE,('W'),record=True)
    
    all_neurons=RE,TC
    all_synapses=S_GABAA_RE_TC,S_GABAB_RE_TC,S_GABAA_RE_RE,S_AMPA_TC_RE
    all_monitors=V1,V2,R1,R2,I1,I2
    
    return all_neurons,all_synapses,all_monitors  

if __name__=='__main__':
    # ## Creation of the network, and simulation
    start_scope()
    close('all')
    
    ###Preferences
    #Define the length of the simulation
    runtime=1000*ms
    #Define the number of neurons and radii of connexions
    N=20
    
    ###Creation of the network
    net = Network(collect())
    all_neurons, all_synapses, all_monitors = create_thalamic_subparts(N)
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_monitors)
    V1,V2,R1,R2,I1,I2 = all_monitors
    S_GABAA_RE_TC,S_GABAB_RE_TC,S_GABAA_RE_RE,S_AMPA_TC_RE = all_synapses
    RE,TC = all_neurons
    
    ###Run the simulation
    prefs.codegen.target = 'cython'
    net.run(1000*ms,report='text',report_period=60*second)
    TC.Iext = 20*namp
    RE.Iext = 20*namp
    net.run(20*ms,report='text',report_period=60*second)
    TC.Iext = 0*namp
    RE.Iext = 0*namp
    net.run(1000*ms,report='text',report_period=60*second)
    TC.Iext = 0*namp
    RE.Iext = 0*namp
    net.run(1000*msecond,report='text',report_period=60*second)
    TC.Iext = 0*namp
    RE.Iext = 0*namp
    net.run(1000*ms,report='text',report_period=60*second)
    # ## Plotting results
    
    # In[10]:
    
    
    n = 6
    
    fig,ax = subplots(6,1, sharex = True)
    ax[0].plot(R1.t/ms, R1.i, '.g',alpha=0.5, label='RE',color="royalblue")
    ax[0].set_ylabel('Neuron index')
    ax[0].set_title('Reticular cells',fontsize=10)
    ax[1].plot(R2.t/ms, R2.i, '.g',alpha=0.5, label='TC')
    ax[1].set_ylabel('Neuron index')
    ax[1].set_title('Thalamic relay cells',fontsize=10)
    ax[2].plot(V2.t/ms, V2.v[n]/mV, label='TC cells',color="tab:green")
    ax[2].set_ylabel('Potential (mV)')
    ax[2].set_title('Membrane potential TC neuron n',fontsize=10)
    ax[3].plot(V1.t/ms, V1.v[n]/mV, label='RE cells',color="royalblue")
    ax[3].set_ylabel('Potential (mV)')
    ax[3].set_title('Membrane potential RE neuron n',fontsize=10)
    #ax[4].plot(I1.t/ms, I1.IsynGABAA_RE_RE[n]/nA, label='IsynGABAA_RE_RE',color="orange")
    #ax[4].plot(I1.t/ms, I1.IsynAMPA_TC_RE[n]/nA, label='IsynAMPA_TC_RE',color="yellow")
    ax[4].plot(I2.t/ms, I2.IsynGABAB_RE_TC[n]/nA, label='IsynGABAB_RE_TC',color="pink")
    #ax[4].plot(I2.t/ms, I2.IsynGABAA_RE_TC[n]/nA, label='IsynGABAA_RE_TC',color="crimson")
    ax[4].set_ylabel('Current (nAmp)')
    ax[4].set_title('Synaptic currents TC and RE neuron n',fontsize=10)
    ax[4].legend(loc='upper right')
    ax[5].plot(V1.t/ms, V1.Iext[n]/nA)
    ax[5].set_ylabel('Current (nAmp)')
    ax[5].set_xlabel('Time (ms)')
    ax[5].set_title('External current applied to both TC and RE',fontsize=10)

    fig,ax = subplots(5,1, sharex = True)
    ax[0].plot(V1.t/ms, V1.CA_i_RE[n]/uM, label='RE cells')
    ax[0].set_ylabel('Concentration (uM)')
    ax[0].plot(V2.t/ms, V2.CA_i_TC[n]/uM, label='TC cells')
    ax[0].set_ylabel('Concentration (uM)')
    ax[0].set_title('Calcium concentration in cells',fontsize=10)
    ax[0].legend()
    ax[1].plot(V1.t/ms, V1.I_t[n]/nA, label='RE cells')
    ax[1].plot(V2.t/ms, V1.I_t[n]/nA, label='TC cells')
    ax[1].set_ylabel('Current (nA)')
    ax[1].set_title('I_t currents in cells',fontsize=10)
    ax[1].legend()
    ax[2].plot(V2.t/ms, V2.I_h[n]/nA, color = "tab:blue")
    ax[2].set_ylabel('Current (nA)')
    ax[2].set_title('I_h current in thalamic relay cells',fontsize=10)
    ax[3].plot(V2.t/ms, V2.Op_L[n],label='OpenLocked',color = "tab:orange")
    ax[3].set_ylabel('')
    ax[3].set_title('I_h channels',fontsize=10)
    ax[3].plot(V2.t/ms, V2.Op[n],label='Open',color = "tab:red")
    ax[3].legend()
    ax[4].plot(V2.t/ms, V2.P1[n],label='',color = "tab:orange")
    ax[4].set_ylabel('')
    ax[4].set_title('Regulating factor bounded form',fontsize=10)
    ax[4].set_xlabel('Time (ms)')
    fig,ax = subplots(4,1, sharex = True)
    ax[0].plot(V1.t/ms, V1.v[14]/mV, label='RE cells',color="royalblue")
    ax[0].set_ylabel('Potential (mV)')
    ax[0].set_title('Membrane potential RE neuron 15',fontsize=10)
    ax[1].plot(V1.t/ms, V1.m_t[14], label='m_t')
    ax[1].plot(V1.t/ms, V1.h_t[14], label='h_t')
    ax[1].set_title('RE cells',fontsize=10)
    ax[1].legend()
    ax[2].plot(V2.t/ms, V2.v[14]/mV, label='TC cells',color="tab:green")
    ax[2].set_ylabel('Potential (mV)')
    ax[2].set_title('Membrane potential TC neuron 15',fontsize=10)
    ax[3].plot(V2.t/ms, V2.m_t[14], label='m_t')
    ax[3].plot(V2.t/ms, V2.h_t[14], label='h_t')
    ax[3].legend()
    ax[3].set_title('TC cells',fontsize=10)
    ax[3].set_xlabel('Time (ms)')#connectivity visualization to check out synapses
    
    def visualise_connectivity(S):
        Ns = len(S.source)
        Nt = len(S.target)
        figure(figsize=(10, 4))
        subplot(121)
        plot(zeros(Ns), arange(Ns), 'ok', ms=10)
        plot(ones(Nt), arange(Nt), 'ok', ms=10)
        for i, j in zip(S.i, S.j):
            plot([0, 1], [i, j], '-k')
        xticks([0, 1], ['Source', 'Target'])
        yticks(np.arange(0, Ns+1, 5.0))
        ylabel('Neuron index')
        xlim(-0.1, 1.1)
        ylim(-1, max(Ns, Nt))
        subplot(122)
        plot(S.i, S.j, 'ok')
        xlim(-1, Ns)
        ylim(-1, Nt)
        xticks(np.arange(0, Ns+1, 1.0))
        yticks(np.arange(0, Nt+1, 5.0))
        xlabel('Source neuron index')
        ylabel('Target neuron index')
        suptitle(f'Synapse Connectivity')
    
    visualise_connectivity(S_AMPA_TC_RE)
    Test = create_thalamic_subparts(5)
    Test[0][0] #RE neuronsTest[1][0] #RE neurons
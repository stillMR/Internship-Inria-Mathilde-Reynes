#!/usr/bin/env python
# coding: utf-8

# ## Required packages and .ipynb

# Soma_eqs.ipynb, Dendritic_eqs.ipynb and Synapses.ipynb are required for the proper functioning of this file.

# In[1]:

from Soma_eqs_210623 import *
from Dendritic_eqs_210623 import *
from Synapses_210623 import *


# In[2]:


from brian2 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


# ## Defining a cortical layer: instantiating neurons, synapses and Poisson inputs

# The function create_cortical_layer permits the generation of a single cortical layer of the prefrontal cortex. A cortical layer is composed of pyramidal (PY) and interneurons/inhibitory cells (IN).
# N is considered to be the number of PY cells needed. The number of IN is logically deduced from it (N/4).

# In[3]:


def create_cortical_layer(N,radii_AMPA_PYPY):
    
    ###Preferences
    prefs.codegen.target = 'cython'
    defaultclock.dt = 0.02*ms
    
    ###Parameters
    #Channel-specific conductances per unit of surface
    g_syn_ampa_pypy = 0.08*usiemens
    g_syn_nmda_pypy = 0.006*usiemens #0.006
    g_syn_ampa_pyin = 0.08*usiemens
    g_syn_nmda_pyin = 0.005*usiemens #0.005
    g_syn_gabaa_inpy = 0.25*usiemens
    # g_syn_ampa_pypy = 0.08*usiemens*0
    # g_syn_nmda_pypy = 0.006*usiemens*0 #0.006
    # g_syn_ampa_pyin = 0.08*usiemens*0
    # g_syn_nmda_pyin = 0.005*usiemens*0 #0.005
    # g_syn_gabaa_inpy = 0.25*usiemens*0
    #Reverse potentials
    E_ampa = 0*mV
    E_nmda = 0*mV
    E_gabaa = -70*mV
    #Fraction of resources used per AP
    U_gabaa = 0.073
    U_ampa = 0.07
    #Rate constants
    alpha_ampa = 1.1*kHz
    beta_ampa = 0.19*kHz
    alpha_nmda = 1*kHz
    beta_nmda = 0.0067*kHz
    alpha_gabaa = 10.5*kHz
    beta_gabaa = 0.166*kHz
    #Number of neurons
    N_PY,N_IN = N,N/4
    #Fictional & arbitrary position of neurons for proper synapses.connect() condition
    neuron_spacing = 1*um
    # #Amplitudes for synapses
    # A_PY_PY = 0.2*usiemens
    # A_PY_IN = 0.05*usiemens
    #Amplitudes for minis
    A_PY_PY = 0.2*usiemens*0.6
    A_PY_IN = 0.05*usiemens*0.6
    
    #Areas of the different neurons
    s_Soma_PYIN = 10**-6*cm**2
    s_Dend_PY = 165*s_Soma_PYIN
    s_Dend_IN = 50*s_Soma_PYIN
    
    ###Instantiate neurons
    #Pyramidal
    PY_dendrite = NeuronGroup(N_PY,Dendritic_eqs,method='rk4',threshold='v>-10*mV',refractory=15*ms)
    PY_dendrite.v = -70*mvolt
    PY_dendrite.m_na = 0.05
    PY_dendrite.h_na = 0.95
    PY_dendrite.m_nap = 0.00
    PY_dendrite.m_km = 0.05
    PY_dendrite.m_kca = 0.00
    PY_dendrite.m_hva = 0.05
    PY_dendrite.h_hva = 0.05
    PY_dendrite.CA_i = 1e-4*mM 
    PY_dendrite.g_na = 0.8*msiemens*cm**-2 
    PY_dendrite.g_nap = 2.5*msiemens*cm**-2 
    PY_dendrite.g_km = 0.02*msiemens*cm**-2
    PY_dendrite.rho = 165
    PY_dendrite.x = 'i*neuron_spacing'

    PY_soma = NeuronGroup(N_PY,Soma_eqs,method='rk4',threshold='v>0*mV',refractory=3*ms,events={'custom_poisson_PY':'rand()<mean_rate_PY*dt','custom_poisson_IN':'rand()<mean_rate_IN*dt','custom_poisson_inter':'rand()<mean_rate_inter*dt'})
    PY_soma.h_na = 0.95
    PY_soma.m_na = 0.05
    PY_soma.m_nap = 0.00
    PY_soma.n_k = 0.05
    PY_soma.g_na = 3000*msiemens*cm**-2
    PY_soma.g_nap = 15*msiemens*cm**-2
    PY_soma.g_k = 200*msiemens*cm**-2
    PY_soma.x = 'i*neuron_spacing'
    
    #Interneurons
    IN_dendrite = NeuronGroup(N_IN,Dendritic_eqs,method='rk4',threshold='v>-10*mV',refractory=15*ms)
    IN_dendrite.v = -70*mvolt
    IN_dendrite.m_na = 0.05
    IN_dendrite.h_na = 0.95
    IN_dendrite.m_nap = 0.00
    IN_dendrite.m_km = 0.05
    IN_dendrite.m_kca = 0.00
    IN_dendrite.m_hva = 0.05
    IN_dendrite.h_hva = 0.05
    IN_dendrite.CA_i = 1e-4*mM
    IN_dendrite.g_na = 0.8*msiemens*cm**-2 
    IN_dendrite.g_nap = 0.0*msiemens*cm**-2 #differs from PY_dendrite
    IN_dendrite.g_km = 0.03*msiemens*cm**-2 #differs from PY_dendrite
    IN_dendrite.rho = 50 #differs from PY_dendrite
    IN_dendrite.x = '(4*i+1.5)*neuron_spacing'
    
    IN_soma = NeuronGroup(N_IN,Soma_eqs,method='rk4',threshold='v>0*mV',refractory=3*ms,events={'custom_poisson_PY':'rand()<mean_rate_PY*dt','custom_poisson_IN':'rand()<mean_rate_IN*dt','custom_poisson_inter':'rand()<mean_rate_inter*dt'})
    IN_soma.h_na = 0.95
    IN_soma.m_na = 0.05
    IN_soma.m_nap = 0.00
    IN_soma.n_k = 0.05
    IN_soma.g_na = 3000*msiemens*cm**-2
    IN_soma.g_nap = 0.00*msiemens*cm**-2 #differs from PY_dendrite
    IN_soma.g_k = 200*msiemens*cm**-2
    IN_soma.x = '(4*i+1.5)*neuron_spacing'
    
    ###Define gap junctions between soma and dendritic compartments for both PY and IN
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
    gapPY_som_den.g = 1 / (R * s_Dend_PY)
    
    gapIN_som_den = Synapses(IN_soma,IN_dendrite,model=eq_gap_dendrite)
    gapIN_som_den.connect(j = 'i')
    gapIN_som_den.g = 1 / (R * s_Dend_IN)
    
    #From dendrites to soma
    gapPY_den_som = Synapses(PY_dendrite,PY_soma,model=eq_gap_soma)
    gapPY_den_som.connect(j = 'i')
    
    gapIN_den_som = Synapses(IN_dendrite,IN_soma,model=eq_gap_soma)
    gapIN_den_som.connect(j = 'i')
    
    ###Instantiate synapses
    #from Bazhenov et al. 2002, connections appears to be from somas to dendrites
    S_AMPA_PY_PY = syn_ampa(PY_soma,PY_dendrite,'IsynAMPA_PY_PY',s_Dend_PY,'abs(i-j)<='+str(radii_AMPA_PYPY)+'and i!=j',g_syn_ampa_pypy,g_syn_ampa_pypy,E_ampa,alpha_ampa,beta_ampa,U_ampa,A_PY_PY,'IEPSPs_PY_PY',0)
    S_AMPA_PY_PY.t_last_spike = -100*ms
    S_AMPA_PY_PY.t_last_spike_Poisson_PY = -10*ms
    S_AMPA_PY_PY.t_last_spike_Poisson_IN = -10*ms
    S_AMPA_PY_PY.D = 1
    S_AMPA_PY_PY.A_mEPSP = A_PY_PY
    S_AMPA_PY_PY.w = 0
    S_AMPA_PY_IN = syn_ampa(PY_soma,IN_dendrite,'IsynAMPA_PY_IN',s_Dend_IN,'abs(x_post-x_pre) <= 8.4*um',g_syn_ampa_pyin,g_syn_ampa_pyin,E_ampa,alpha_ampa,beta_ampa,U_ampa,A_PY_IN,'IEPSPs_PY_IN',1)
    S_AMPA_PY_IN.t_last_spike = -100*ms
    S_AMPA_PY_IN.t_last_spike_Poisson_PY = -10*ms
    S_AMPA_PY_IN.t_last_spike_Poisson_IN = -10*ms
    S_AMPA_PY_IN.D = 1
    S_AMPA_PY_IN.A_mEPSP = A_PY_IN
    S_AMPA_PY_IN.w = 0
    S_NMDA_PY_PY = syn_nmda(PY_soma,PY_dendrite,'IsynNMDA_PY_PY',s_Dend_PY,'abs(i-j)<=5 and i!=j',g_syn_nmda_pypy,E_nmda,alpha_nmda,beta_nmda) #not specified in the paper, took the same as AMPA (from info on Bazhenov et al.,2002)
    S_NMDA_PY_PY.t_last_spike = -100*ms
    S_NMDA_PY_IN = syn_nmda(PY_soma,IN_dendrite,'IsynNMDA_PY_IN',s_Dend_IN,'abs(x_post-x_pre) <= 8.4*um',g_syn_nmda_pyin,E_nmda,alpha_nmda,beta_nmda) #same
    S_NMDA_PY_IN.t_last_spike = -100*ms
    S_GABAA_IN_PY = syn_gabaa(IN_soma,PY_dendrite,'IsynGABAA_IN_PY',s_Dend_PY,'abs(x_post-x_pre) <= 4.6*um',g_syn_gabaa_inpy,E_gabaa,alpha_gabaa,beta_gabaa,U_gabaa)
    S_GABAA_IN_PY.t_last_spike = -100*ms
    S_GABAA_IN_PY.D = 1
    
    ###Define monitors
    #general neuron monitoring
    # V1=StateMonitor(PY_dendrite,('v','I_nap','I_kca','I_hva','I_km','CA_i','I_kl','I_l','I_na'),record=True)
    # V2=StateMonitor(PY_soma,('v','I_nap','I_k','I_na','mean_rate_PY','t_last_spike_PY'),record=True)
    # V3=StateMonitor(IN_dendrite,('v','I_nap','I_kca','I_hva','I_km','CA_i','I_kl','I_l','I_na'),record=True)
    # V4=StateMonitor(IN_soma,('v','I_nap','I_k','I_na'),record=True)   
    
    V1=StateMonitor(PY_dendrite,('v'),record=True)
    V2=StateMonitor(PY_soma,('v'),record=True)
    V3=StateMonitor(IN_dendrite,('v'),record=True)
    V4=StateMonitor(IN_soma,('v'),record=True) 
    
    #spike monitoring
    #R1=SpikeMonitor(PY_dendrite,record=True)
    R2=SpikeMonitor(PY_soma,record=True)
    #R3=SpikeMonitor(IN_dendrite,record=True)
    R4=SpikeMonitor(IN_soma,record=True)
    #synaptic currents monitoring
    #I1=StateMonitor(PY_dendrite,('Isyn','IEPSPs_PY_PY'),record=True)
    # I1=StateMonitor(PY_dendrite,('Isyn','IEPSPs_PY_PY','IEPSPs_L_V_L_VI','IEPSPs_L_V_L_II','IEPSPs_L_II_L_V','IEPSPs_L_II_L_VI','IEPSPs_L_VI_L_V','IEPSPs_L_VI_L_II'),record=True)
    # I2=StateMonitor(IN_dendrite,('Isyn','IEPSPs_PY_IN'),record=True)
    I1=StateMonitor(PY_dendrite,'Isyn',record=[0])
    I2=StateMonitor(IN_dendrite,'Isyn',record=[0])
    #synapses monitoring (optionnal)
    #S1=StateMonitor(S_AMPA_PY_PY,('w','g_syn','t_last_spike','apost','apre','D','t_last_spike_Poisson_PY','A_mEPSP'),record=True)
    #S2=StateMonitor(S_AMPA_PY_IN,('W'),record=True)
    #S3=StateMonitor(S_NMDA_PY_PY,('W'),record=True)
    #S4=StateMonitor(S_NMDA_PY_IN,('W'),record=True)
    #S5=StateMonitor(S_GABAA_IN_PY,('W'),record=True)
    #miniEPSPs monitoring
    M0=EventMonitor(PY_soma,'custom_poisson_PY')
    M1=EventMonitor(PY_soma,'custom_poisson_IN')
    M2=EventMonitor(PY_soma,'custom_poisson_inter')
    
    all_neurons=PY_dendrite,PY_soma,IN_dendrite,IN_soma
    all_gap_junctions=gapPY_som_den,gapPY_den_som,gapIN_som_den,gapIN_den_som
    all_synapses=S_AMPA_PY_PY,S_AMPA_PY_IN,S_NMDA_PY_PY,S_NMDA_PY_IN,S_GABAA_IN_PY
    all_monitors=V1,V2,V3,V4,R2,R4,I1,I2,M0,M1,M2
    
    return all_neurons,all_synapses,all_gap_junctions,all_monitors  

if __name__=='__main__':
    # ## Creation of the network, and simulation
    start_scope()
    close('all')
    
    ###Preferences
    #Define the length of the simulation
    runtime=20000*ms
    
    ###Creation of the network
    net = Network(collect())
    all_neurons, all_synapses, all_gap_junctions, all_monitors = create_cortical_layer(100,5)
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
    net.add(all_monitors)
    V1,V2,V3,V4,R2,R4,I1,I2,M0,M1,M2  = all_monitors
    #PY_dendrite,PY_soma,IN_dendrite,IN_soma = all_neurons
    #S_AMPA_PY_PY,S_AMPA_PY_IN,S_NMDA_PY_PY,S_NMDA_PY_IN,S_GABAA_IN_PY = all_synapses
    
    num_samples = int(runtime/defaultclock.dt)
    init_arr = zeros(num_samples)
    init_arr[0]=1
    init_timedarray = TimedArray(init_arr, dt=defaultclock.dt)
    
    ###Run the simulation
    prefs.codegen.target = 'cython'
    net.run(runtime,report='text',report_period=60*second)
    # ## Plotting results
    
    fig,ax = subplots(4,1, sharex = True)
    ax[0].plot(R2.t/ms, R2.i, '.g',alpha=0.5, label='PY_soma',color="royalblue")
    ax[0].set_ylabel('Neuron index')
    ax[0].set_title('Pyramidal cells',fontsize=10)
    #ax[2].plot(V1.t/ms, V1.v[0]/mV, label='PY dendrite',color="royalblue")
    ax[1].plot(V2.t/ms, V2.v[0]/mV, label='PY soma',color="royalblue")
    ax[1].set_ylabel('Potential (mV)')
    ax[1].set_title('Membrane potential PY neuron n',fontsize=10)
    
    # ax[2].plot(R4.t/ms, R4.i, '.b',alpha=0.5, label='PY_soma',color="royalblue")
    # ax[2].set_ylabel('Neuron index')
    # ax[2].set_title('Interneurons',fontsize=10)
    ax[2].plot(I1.t/ms, I1.IEPSPs_PY_PY[0]/nA, label='mEPSPs_PY_PY',color="tab:green")
    ax[2].set_ylabel('Current (Amp)')
    ax[2].set_title('MiniEPSPs PY neuron n',fontsize=10)
    # ax[3].plot(I1.t/ms, I1.IsynAMPA_PY_PY[0]/nA, label='IsynAMPA_PY_PY',color="orange")
    # ax[3].plot(I1.t/ms, I1.IsynNMDA_PY_PY[0]/nA, label='IsynNMDA_PY_PY',color="crimson")
    #ax[4].plot(I1.t/ms, I1.IsynGABAA_IN_PY[20]/amp, label='IsynGABAA_IN_PY',color="orchid")
    # ax[3].set_ylabel('Current (nAmp)')
    # ax[3].set_xlabel('Time (ms)')
    # ax[3].set_title('Synaptic currents PY neuron n',fontsize=10)
    # ax[3].legend()
    ax[3].plot(M0.t/ms, M0.i, '.g',alpha=0.5)
    ax[3].set_title('Poisson Spikes',fontsize=10)
    ax[3].set_ylabel('Neuron index')
    
    figure()
    plot(V2.t/ms, V2.mean_rate_PY[0])
    
    figure()
    plot(V2.t/ms, V2.t_last_spike_PY[0])
    
    # fig,ax = subplots(6,1, sharex = True)
    # fig.suptitle("Poisson inputs through AMPA receptors", fontsize="x-large")
    
    # ax[0].set_yticks(np.arange(0,26,5))
    # ax[0].plot(M1.t/ms, M1.i, '.g',alpha=0.5, label='PY_soma', color="royalblue")
    # ax[0].set_title('Poisson Spikes',fontsize=10)
    # ax[0].set_ylabel('Neuron index')
    # ax[0].grid(True,axis='y', linestyle='--')
    
    # ax[1].plot(R2.t/ms, R2.i, '.g',alpha=0.5, label='PY_soma',color="tab:green")
    # ax[1].set_ylabel('Neuron index')
    # ax[1].set_yticks(np.arange(0,26,5))
    # ax[1].set_title('Pyramidal cells',fontsize=10)
    # ax[1].grid(True,axis='y', linestyle='--')
    
    # ax[2].plot(I1.t/ms, I1.IEPSPs_PYIN[0]/nA, color="royalblue")
    # ax[2].set_title('miniEPSPs from Pyramidal to Pyramidal',fontsize=10)
    # ax[2].set_ylabel('Current (nAmp)')
    
    # ax[3].plot(V2.t/ms, V2.v[1]/mV, color="red", label="pre")
    # ax[3].plot(V1.t/ms, V1.v[0]/mV, color="tab:green", label="post")
    # ax[3].legend()
    # ax[3].set_title('Post and pre neurons',fontsize=10)
    # ax[3].set_ylabel('Membrane potential (mV)')
    
    # ax[4].plot(S1.t/ms, S1.A_mEPSP[5]/usiemens, color="royalblue")
    # ax[4].set_title('mEPSP amplitude',fontsize=10)
    # ax[4].set_ylabel('(usiemens)')
    
    # ax[5].plot(S1.t/ms, S1.g_syn[5]/usiemens, color="royalblue")
    # ax[5].set_title('Synaptic conductance',fontsize=10)
    # ax[5].set_ylabel('(usiemens)')
    
    
    # fig,ax = subplots(7,1, sharex = True)
    # ax[0].plot(V2.t/ms, V2.v[10]/mV, label='PY soma',color="royalblue")
    # ax[0].set_ylabel('Potential (mV)')
    # ax[0].set_title('Membrane potential PY soma index 10',fontsize=10)
    # ax[1].plot(V1.t/ms, V1.v[11]/mV, label='PY dendrite',color="royalblue")
    # ax[1].set_ylabel('Potential (mV)')
    # ax[1].set_title('Membrane potential PY dendrite index 11',fontsize=10)
    # ax[2].plot(S1.t/ms, S1.w[90], label='',color="tab:green")
    # ax[2].set_ylabel('')
    # ax[2].set_title('w',fontsize=10)
    # ax[3].plot(S1.t/ms, S1.t_last_spike[90]/ms, label='',color="tab:green")
    # ax[3].set_ylabel('Time (ms)')
    # ax[3].set_title('Time last spike',fontsize=10)
    # ax[4].plot(S1.t/ms, S1.g_syn[90]/usiemens, label='',color="tab:green")
    # ax[4].set_ylabel('Resistance (usiemens)')
    # ax[4].set_title('g_syn',fontsize=10)
    # ax[5].plot(S1.t/ms, S1.apre[90], label='',color="tab:green")
    # ax[5].set_ylabel('')
    # ax[5].set_title('apre',fontsize=10)
    # ax[6].plot(S1.t/ms, S1.D[90], label='',color="tab:green")
    # ax[6].set_ylabel('')
    # ax[6].set_title('D',fontsize=10)
    # ax[6].set_xlabel('Time (ms)')
    # fig.suptitle('AMPA 90 linking PY soma 10 to PY dendrite 11', fontsize=15)
    
    # fig,ax = subplots(9,1, sharex = True)
    # fig.suptitle("Poisson inputs through AMPA receptors", fontsize="x-large")
    
    # ax[0].set_yticks(np.arange(0,26,5))
    # ax[0].plot(M1.t/ms, M1.i, '.g',alpha=0.5, label='PY_soma', color="royalblue")
    # ax[0].set_title('Poisson Spikes',fontsize=10)
    # ax[0].set_ylabel('Neuron index')
    # ax[0].grid(True,axis='y', linestyle='--')
    
    # ax[1].plot(R2.t/ms, R2.i, '.g',alpha=0.5, label='PY_soma',color="tab:green")
    # ax[1].set_ylabel('Neuron index')
    # ax[1].set_yticks(np.arange(0,26,5))
    # ax[1].set_title('Pyramidal cells',fontsize=10)
    # ax[1].grid(True,axis='y', linestyle='--')
    
    # ax[2].plot(I1.t/ms, I1.IEPSPs_PYIN[0]/nA, color="royalblue")
    # ax[2].set_title('miniEPSPs from Pyramidal to Pyramidal',fontsize=10)
    # ax[2].set_ylabel('Current (nAmp)')
    
    # ax[3].plot(S1.t/ms, S1.t_last_spike_Poisson[5]/ms, color="royalblue")
    # ax[3].set_ylabel('Time (ms)')
    # ax[3].set_title('Time last spike from Poisson',fontsize=10)
    
    # ax[4].plot(V2.t/ms, V2.v[1]/mV, color="red", label="pre")
    # ax[4].plot(V1.t/ms, V1.v[0]/mV, color="tab:green", label="post")
    # ax[4].legend()
    # ax[4].set_title('Post and pre neurons',fontsize=10)
    # ax[4].set_ylabel('Membrane potential (mV)')
    
    # ax[5].plot(S1.t/ms, S1.A_mEPSP[5]/usiemens, color="royalblue")
    # ax[5].set_title('mEPSP amplitude',fontsize=10)
    # ax[5].set_ylabel('(usiemens)')
    
    # ax[6].plot(S1.t/ms, S1.w[5], color="royalblue")
    # ax[6].set_title('w',fontsize=10)
    # ax[6].set_ylabel('STDP (unitless)')
    
    # ax[7].plot(S1.t/ms, S1.apre[5], color="royalblue")
    # ax[7].set_title('apre',fontsize=10)
    # ax[7].set_ylabel('(unitless)')
    
    # ax[8].plot(S1.t/ms, S1.apost[5], color="royalblue")
    # ax[8].set_title('apost',fontsize=10)
    # ax[8].set_ylabel('(unitless)')
    # ax[8].set_xlabel('Time (ms)')
    
    # S_AMPA_PY_PY.j[5]
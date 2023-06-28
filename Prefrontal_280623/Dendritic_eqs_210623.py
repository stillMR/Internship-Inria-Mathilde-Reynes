#!/usr/bin/env python
# coding: utf-8

# ## Required packages

# In[1]:


from brian2 import *
import numpy as np


# ## Global parameters

# In[2]:


#Area
s_soma = 10**-6*cm**2

#Membrane capacitance per 2unit of surface
# Cm = 1*ufarad/cm**2
Cm = 0.75*ufarad/cm**2

#Ion-specific conductances per unit of surface
g_hva = 0.02*msiemens*cm**-2 #0.02
g_kca = 0.3*msiemens*cm**-2 #0.3
g_l = 0.03*msiemens*cm**-2 #0.033
g_kl = 0.003*msiemens*cm**-2 #0.003

#Ion-specific reversal potential
E_na = 50*mV 
E_k = -90*mV
E_ca = 140*mV
E_l = -68*mV
E_kl = -95*mV

#Calcium constants
#A = 5.1819E-5*mM/(ms*uA)*(cm**2/(s_soma*rho))*1000
CAinf = 2.4E-4*mM

#Time constants
tauCA = 165*ms #5ms in paper
tau_m_nap = 0.1991*ms

#Other constants 
Qt = pow(2.3,((36-23)/10))
Tad = pow(2.3,((36-23)/10))

F=96485*coulomb/mole
A_test=(1*10**5/(2*umeter*F))
A_Bazhenov=(5.1819E-5*mM*cm**2)/(ms*uA)


# ## Equations

# In[3]:


Dendritic_eqs = '''

    dv/dt = (- I_kl - I_na - I_nap - I_km - I_kca - I_hva - I_l - Isyn - Igap - IEPSPs + Iext) * (1/Cm) : volt 

    I_kl = g_kl * (v - E_kl) : amp * meter**-2
    
    I_na = g_na * (m_na ** 3) * h_na * (v - E_na) : amp * meter**-2
        dm_na/dt = -(m_na - m_nainf) / tau_m_na : 1
        dh_na/dt = -(h_na - h_nainf) / tau_h_na : 1
        alpham_na = (0.182/ms * (v + 25*mV)/mV / (1 - exp(-(v + 25*mV)/9/mV))) * int(abs((v-10*mV)/(-35*mV)) > 1e-6) + (0.182/ms * 9) * int(abs((v-10*mV)/(-35*mV)) < 1e-6)  : Hz
        betam_na = (-0.124/ms * (v + 25*mV)/mV / (1 - exp((v + 25*mV)/9/mV))) * int(abs((-v+10*mV)/(35*mV)) > 1e-6) + (0.124/ms * 9) * int(abs((-v+10*mV)/(35*mV)) < 1e-6)  : Hz
        alphah_na = (0.024/ms * (v + 40*mV)/mV / (1 - exp(-(v + 40*mV)/5/mV))) * int(abs((v-10*mV)/(-50*mV)) > 1e-6) + (0.024/ms * 5) * int(abs((v-10*mV)/(-50*mV)) < 1e-6)  : Hz
        betah_na = (-0.0091/ms * (v + 65*mV)/mV / (1 - exp((v + 65*mV)/5/mV))) * int(abs((-v+10*mV)/(75*mV)) > 1e-6) + (0.0091/ms * 5) * int(abs((-v+10*mV)/(75*mV)) < 1e-6)  : Hz
        h_nainf = 1 / (1 + exp((v + 55*mV)/6.2/mV)) : 1
        tau_h_na = (1 / (alphah_na + betah_na)) / Qt : second   
        m_nainf = alpham_na / (alpham_na + betam_na) : 1
        tau_m_na = (1 / (alpham_na + betam_na)) / Qt : second
        
    I_nap = g_nap * m_nap * (v - E_na) : amp * meter**-2
        dm_nap/dt = -(m_nap -  m_napinf)/tau_m_nap : 1
        m_napinf = 0.02 / (1 + exp(-(v + 42*mV)/5/mV)) : 1 

    I_km = g_km * m_km * (v - E_k) : amp * meter**-2
        dm_km/dt = -(m_km - m_kminf) / tau_m_km : 1
        alpham_km = 0.001/ms * (v + 30*mV)/mV / (1 - exp(-(v + 30*mV)/9/mV)) : Hz
        betam_km = -0.001/ms * (v + 30*mV)/mV / (1 - exp((v + 30*mV)/9/mV)) : Hz
        m_kminf = alpham_km / (alpham_km + betam_km) : 1
        tau_m_km = (1 / (alpham_km + betam_km)) / Qt: second
        
    I_kca = g_kca * m_kca * (v - E_k) : amp * meter**-2
        dm_kca/dt = -(m_kca - m_kcainf) / tau_m_kca : 1
        alpham_kca = (0.01/ms) * (CA_i/uM) :  Hz
        betam_kca = 0.02/ms : Hz
        m_kcainf = alpham_kca / (alpham_kca + betam_kca) : 1
        tau_m_kca = (1 / (alpham_kca + betam_kca)) / Qt: second 
        drive = -A * I_hva : katal * meter**-3
        A = (5.1819E-5*mM*cm**2)/(ms*uA) : meter**-1 * second**-1 * amp**-1 * mol
        dCA_i/dt = (drive + (CAinf - CA_i) / tauCA) * int(drive > 0*katal*meter**-3) + (0*katal*meter**-3 + (CAinf - CA_i) / tauCA) * int(drive < 0*katal*meter**-3) : mM 

    I_hva = g_hva * (m_hva ** 2) * h_hva * (v - E_ca) : amp * meter **-2
        dm_hva/dt = -(m_hva - m_hvainf) / tau_m_hva : 1
        dh_hva/dt = -(h_hva - h_hvainf) / tau_h_hva : 1
        alpham_hva = 0.055/ms * (v + 27*mV)/mV / (1 - exp(-(v + 27*mV)/3.8/mV)) : Hz
        betam_hva = 0.94 * exp((-v - 75*mV)/17/mV)/ms : Hz
        alphah_hva = 0.000457 * exp(-(v + 13*mV)/50/mV)/ms : Hz
        betah_hva = 0.0065 / (1 + exp(-(v + 15*mV)/28/mV))/ms : Hz
        m_hvainf = alpham_hva / (alpham_hva + betam_hva) : 1
        tau_m_hva = (1 / (alpham_hva + betam_hva)) / Qt : second 
        h_hvainf = 1 / (1 + exp((v + 55*mV)/6.2/mV)) : 1
        tau_h_hva = (1 / (alphah_hva + betah_hva)) / Qt : second   
  
    I_l = g_l * (v - E_l) : amp * meter**-2
    
    Isyn = (IsynAMPA_PY_PY + IsynNMDA_PY_PY + IsynAMPA_PY_IN + IsynNMDA_PY_IN + IsynGABAA_IN_PY + IsynAMPA_TC_PY + IsynAMPA_TC_IN + IsynNMDA_PY_PY_layer + IsynAMPA_PY_PY_layer + IsynAMPA_TC_PY_TCo + IsynAMPA_TC_IN_TCo + IsynAMPA_TC_PY_TM + IsynAMPA_TC_IN_TM ): amp * meter**-2
        IsynAMPA_PY_PY : amp * meter**-2
        IsynNMDA_PY_PY :  amp * meter**-2
        IsynAMPA_PY_IN :  amp * meter**-2
        IsynNMDA_PY_IN : amp * meter**-2
        IsynGABAA_IN_PY :  amp * meter**-2 
        IsynAMPA_TC_PY : amp * meter**-2
        IsynAMPA_TC_IN : amp * meter**-2
        IsynNMDA_PY_PY_layer = IsynAMPA_L_V_L_VI + IsynAMPA_L_V_L_II + IsynAMPA_L_II_L_V + IsynAMPA_L_II_L_VI + IsynAMPA_L_VI_L_II + IsynAMPA_L_VI_L_V :  amp * meter**-2
            IsynAMPA_L_V_L_VI : amp * meter**-2
            IsynAMPA_L_V_L_II : amp * meter**-2
            IsynAMPA_L_II_L_V : amp * meter**-2
            IsynAMPA_L_II_L_VI : amp * meter**-2
            IsynAMPA_L_VI_L_II : amp * meter**-2
            IsynAMPA_L_VI_L_V : amp * meter**-2
        IsynAMPA_PY_PY_layer = IsynNMDA_L_V_L_VI + IsynNMDA_L_V_L_II + IsynNMDA_L_II_L_V + IsynNMDA_L_II_L_VI + IsynNMDA_L_VI_L_II + IsynNMDA_L_VI_L_V :  amp * meter**-2 
            IsynNMDA_L_V_L_VI : amp * meter**-2
            IsynNMDA_L_V_L_II : amp * meter**-2
            IsynNMDA_L_II_L_V : amp * meter**-2
            IsynNMDA_L_II_L_VI : amp * meter**-2
            IsynNMDA_L_VI_L_II : amp * meter**-2
            IsynNMDA_L_VI_L_V : amp * meter**-2
        IsynAMPA_TC_PY_TCo : amp * meter**-2
        IsynAMPA_TC_IN_TCo : amp * meter**-2
        IsynAMPA_TC_PY_TM : amp * meter**-2
        IsynAMPA_TC_IN_TM : amp * meter**-2
        IsynAMPA_TC_INVI_TCo : amp * meter**-2
        
    IEPSPs = IEPSPs_intralayer + IEPSPs_interlayer + IEPSPs_thalcort : amp * meter**-2
        IEPSPs_intralayer = IEPSPs_PY_PY + IEPSPs_PY_IN : amp * meter**-2
            IEPSPs_PY_PY : amp * meter**-2
            IEPSPs_PY_IN : amp * meter**-2
        IEPSPs_interlayer = IEPSPs_L_V_L_VI + IEPSPs_L_V_L_II + IEPSPs_L_II_L_V + IEPSPs_L_II_L_VI + IEPSPs_L_VI_L_V + IEPSPs_L_VI_L_II : amp * meter**-2
            IEPSPs_L_V_L_VI : amp * meter**-2
            IEPSPs_L_V_L_II : amp * meter**-2
            IEPSPs_L_II_L_V : amp * meter**-2
            IEPSPs_L_II_L_VI : amp * meter**-2
            IEPSPs_L_VI_L_II : amp * meter**-2
            IEPSPs_L_VI_L_V : amp * meter**-2
        IEPSPs_thalcort = IEPSPs_TCTCo_PYL_II + IEPSPs_TCTCo_INL_II + IEPSPs_TCM_PYL_V + IEPSPs_TCM_INL_V : amp * meter**-2
            IEPSPs_TCTCo_PYL_II : amp * meter**-2
            IEPSPs_TCTCo_INL_II : amp * meter**-2
            IEPSPs_TCM_PYL_V : amp * meter**-2
            IEPSPs_TCM_INL_V : amp * meter**-2
        
    Igap : amp * meter**-2 
    
    Iext : amp * meter**-2
    
    rho : 1
    
    g_na : siemens * meter**-2
    
    g_nap : siemens * meter**-2
    
    g_km : siemens * meter**-2
    
    x : meter

    '''

if __name__=='__main__':
    
    # ## Plotting results
    close('all')    
    prefs.codegen.target = "numpy"
    defaultclock.dt = 0.02*ms
    
    #Parameters initialization
    PY_dendrite = NeuronGroup(25,Dendritic_eqs,method='rk4',threshold='v>0*mV',refractory=3*ms)
    PY_dendrite.v = -70*mvolt
    PY_dendrite.m_na = 0.05
    PY_dendrite.h_na = 0.95
    PY_dendrite.m_nap = 0.0
    PY_dendrite.m_km = 0.05
    PY_dendrite.m_kca = 0.05
    PY_dendrite.m_hva = 0.05
    PY_dendrite.h_hva = 0.95
    PY_dendrite.CA_i = 3.5E-4*mM #not from the paper
    PY_dendrite.g_na = 0.8*msiemens*cm**-2 
    PY_dendrite.g_nap = 2.5*msiemens*cm**-2 
    PY_dendrite.g_km = 0.02*msiemens*cm**-2
    # PY_dendrite.g_km = 0.8*msiemens*cm**-2
    PY_dendrite.rho = 165
    
    #Recording
    V1 = StateMonitor(PY_dendrite,('v','tau_h_na','tau_m_na','m_na','h_na','m_nap','m_km','m_kca','m_hva','h_hva','alpham_na','alphah_na','betam_na','betah_na','m_nainf','h_nainf','alpham_km','betam_km','Igap','CA_i','I_hva','I_kca','I_kl','I_na','I_nap','I_km','I_l','drive'),record=True)
    
    #Test
    runtime = 700*ms
    num_samples = int(runtime/defaultclock.dt)
    init_arr = zeros(num_samples)
    init_arr[0]=1
    init_timedarray = TimedArray(init_arr, dt=defaultclock.dt)
    
    PY_dendrite.Iext = 0*amp*meter**-2
    run(200*ms)
    store('initial')
    PY_dendrite.Iext = 0.025*amp*meter**-2 #test of a virtual and arbitrary external input
    run(3000*ms)
    PY_dendrite.Iext = 0*amp*meter**-2
    run(200*ms)    #visualize
    
    fig,ax = subplots(3,1, sharex = True)
    ax[0].plot(V1.t/ms, V1.m_na[0], label='m_na')
    ax[0].plot(V1.t/ms, V1.h_na[0], label='h_na')
    ax[0].plot(V1.t/ms, V1.h_hva[0], label='h_hva')
    ax[0].plot(V1.t/ms, V1.m_hva[0], label='m_hva')
    ax[0].plot(V1.t/ms, V1.m_km[0], label='m_km')
    ax[0].plot(V1.t/ms, V1.m_kca[0], label='m_kca')
    ax[0].legend(loc='upper left')
    ax[0].set_ylabel ('Activation')
    ax[1].plot(V1.t/ms, V1.v[0]/mV, label='v')
    ax[1].legend()
    ax[1].set_ylabel ('Voltage (mV)')
    #ax[2].plot(V1.t/ms, V1.alpham_na[0], label='alpham_na')
    #ax[2].plot(V1.t/ms, V1.alphah_na[0], label='alphah_na')
    #ax[2].plot(V1.t/ms, V1.betam_na[0], label='betam_na')
    #ax[2].plot(V1.t/ms, V1.betah_na[0], label='betah_na')
    #ax[2].plot(V1.t/ms, V1.alpham_km[0], label='alpham_km')
    #ax[2].plot(V1.t/ms, V1.betam_km[0], label='betam_km')
    #ax[2].legend()
    #ax[2].plot(V1.t/ms, V1.tau_m_na[0]/ms, label='tau_m_na')
    #ax[2].plot(V1.t/ms, V1.tau_h_na[0]/ms, label='tau_h_na')
    #ax[2].legend()
    #ax[2].set_ylabel ('ms')
    ax[2].plot(V1.t/ms, V1.Igap[0]/pA, label='Igap')
    ax[2].legend(loc='upper left')
    ax[2].set_ylabel ('Current (pA)')
    ax[2].set_xlabel ('Time (ms)')
    
    figure()
    plot(V1.t/ms, V1.CA_i[0]/mM)
    
    figure()
    plot(V1.t/ms, V1.drive[0])
    
    #figure()
    #plot(V1.t/ms,(0.01/ms) * (V1.CA_i[0]/mM),label='alpha_kca')
    
    #figure()
    #plot(V1.t/ms, V1.m_kca[0])

#'I_hva','I_kca','I_kl','I_na','I_nap','I_km','I_kca','I_l'    
    figure()
    plot(V1.t/ms, V1.I_l[0],label='I_l')
    plot(V1.t/ms, V1.I_kl[0],label='I_kl')
    plot(V1.t/ms, V1.I_na[0],label='I_na')
    plot(V1.t/ms, V1.I_nap[0],label='I_nap')
    plot(V1.t/ms, V1.I_km[0],label='I_km')
    plot(V1.t/ms, V1.I_kca[0],label='I_kca')
    plot(V1.t/ms, V1.I_hva[0],label='I_hva')
    legend()
    
    vv=[-100*mV+i*mV for i in range(150)]
    alpham_km_test = [0.001/ms * (v + 30*mV)/mV / (1 - exp(-(v + 30*mV)/9/mV)) for v in vv]
    betam_km_test =  [-0.001/ms * (v + 30*mV)/mV / (1 - exp((v + 30*mV)/9/mV)) for v in vv]
    m_kminf_test = [alpham_km_test[i] / (alpham_km_test[i] + betam_km_test[i]) for i in range(150)]
    tau_m_km_test = [(1 / (alpham_km_test[i] + betam_km_test[i])) / Qt for i in range(150)]
    
    # figure()
    # subplot(121)
    # plot(vv/mV,m_kminf_test)
    # subplot(122)
    # plot(vv/mV,tau_m_km_test/ms)
#!/usr/bin/env python
# coding: utf-8

# ## Required packages

# In[2]:


from brian2 import *
import numpy as np


# ## Global parameters 

# In[3]:


#Membrane capacitance per 2unit of surface
Cm_RE = 1*ufarad/cm**2

#Conductances
g_kl_RE = 0.005*msiemens*cm**-2
g_na_RE = 100*msiemens*cm**-2
g_t = 2.3*msiemens*cm**-2
g_l_RE = 0.05*msiemens*cm**-2
g_k_RE = 10*msiemens*cm**-2
g_t = 2.3*msiemens*cm**-2

#Reversal potentials
E_kl = -95*mV
E_l_TC = -70*mV
E_l_RE = -77*mV
E_na = 50*mV
E_k = -95*mV
E_ca0 = 1000*8.31441*(273.15 + 36)/(2*96489)*mV #13.31*mV

#Calcium parameters
tau_CA_RE = 5*ms
A_RE = 5.1819E-5*(mM*cm**2)/(ms*uA)
CA_inf = 2.4E-4*mM #2.4E-4*mM
CA_0 = 2*mM

#Temperature-dependent variables 
T = 36
Qm_RE = 5**((T-24)/10)
Qh_RE = 3**((T-24)/10)
Qt = 2.9529
Q = 2.3
Qhyp = pow(3,((T-36)/10))

#Parameters to match modelDB
Vtr = -50*mV
VtrK = -50*mV


# ## Equations

# In[5]:


RE_eqs = '''

    dv/dt = (- I_kl - I_na - I_k - I_t - I_l - Isyn_RE + Iext + Idummy) * (1/Cm_RE) : volt
    v2 = v - Vtr : volt
    v2K = v - VtrK : volt
    
    
    I_kl = g_kl_RE * (v - E_kl) : amp * meter**-2
    
    
    I_l = g_l_RE * (v - E_l_RE) : amp * meter**-2


    I_na = g_na_RE * (m_na ** 3) * h_na * (v - E_na) : amp * meter**-2
        dm_na/dt = -(m_na - m_nainf) / tau_m_na : 1
        dh_na/dt = -(h_na - h_nainf) / tau_h_na : 1
        
        alpham_na = 0.32/ms * (13*mV - v2)/mV / (exp((13*mV - v2)/4/mV) - 1) : Hz
        betam_na = 0.28/ms * (v2 - 40*mV)/mV / (exp((v2 - 40*mV)/5/mV) - 1) : Hz
        m_nainf = alpham_na / (alpham_na + betam_na) : 1
        tau_m_na = (1 / (alpham_na + betam_na)) / Qt : second
        
        alphah_na = 0.128 * exp((17*mV - v2)/18/mV)/ms : Hz
        betah_na = 4/(exp((40*mV - v2)/5/mV + 1))/ms  : Hz
        h_nainf = alphah_na / (alphah_na + betah_na) : 1
        tau_h_na = (1 / (alphah_na + betah_na)) / Qt : second                  


    I_k = g_k_RE * (n_k ** 4) * (v - E_k) : amp * meter**-2 
        dn_k/dt = -(n_k - n_kinf) / tau_n_k : 1
        
        alphan_k = 0.032/ms * (15*mV - v2K)/mV / (exp((15*mV - v2K)/5/mV) - 1) : Hz
        betan_k = 0.5/ms * exp((10*mV - v2K)/40/mV) : Hz
        n_kinf = alphan_k / (alphan_k + betan_k) : 1
        tau_n_k = (1 / (alphan_k + betan_k)) / Qhyp : second   


    I_t = g_t * (m_t ** 2) * h_t * (v - E_ca) : amp * meter**-2
        dm_t/dt = -(m_t - m_tinf) / tau_m_t : 1
        dh_t/dt = -(h_t - h_tinf) / tau_h_t : 1
        
        tau_m_t = (3*ms + 1*ms/(exp((v + 27*mV)/10/mV) + exp(-(v + 102*mV)/15/mV))) / Qm_RE : second
        m_tinf = 1 / (1 + exp(-(v + 52*mV)/7.4/mV)) : 1
        
        tau_h_t = (85*ms + 1*ms/(exp((v + 48*mV)/4/mV) + exp(-(v + 407*mV)/50/mV))) / Qh_RE : second
        h_tinf = 1 / (1 + exp((v + 80*mV)/5/mV)) : 1
        
        drive = -A_RE * I_t : katal * meter**-3
        dCA_i_RE/dt = (drive + (CA_inf - CA_i_RE)/tau_CA_RE) * int(drive > 0*katal*meter**-3) + (0*katal*meter**-3 + (CA_inf - CA_i_RE)/tau_CA_RE) * int(drive < 0*katal*meter**-3) : mM 
        ratio = CA_0/CA_i_RE : 1
        E_ca = E_ca0 * log(ratio) : volt
        
    
    Isyn_RE = IsynGABAA_RE_RE + IsynAMPA_TC_RE + IsynAMPA_PY_RE_layer_VI + IsynAMPA_PY_RE_layer_V : amp * meter**-2 
        IsynGABAA_RE_RE : amp * meter**-2
        IsynAMPA_TC_RE : amp * meter**-2
        IsynAMPA_PY_RE_layer_VI :  amp * meter**-2 
        IsynAMPA_PY_RE_layer_V :  amp * meter**-2
    
    
    Iext : amp * meter**-2
    
    
    x : meter
    
    
    Idummy = IEPSPs_L_VI_RETCo + IEPSPs_L_V_RETCo + IEPSPs_TC_RE : amp * meter**-2
        IEPSPs_L_VI_RETCo : amp * meter**-2
        IEPSPs_L_V_RETCo : amp * meter**-2
        IEPSPs_TC_RE : amp * meter**-2
        

    '''

if __name__=='__main__':
    # ## Plotting results
    prefs.codegen.target = "cython"
    defaultclock.dt = 0.02*ms
    
    #Parameters initialization
    RE = NeuronGroup(2,RE_eqs,method='exponential_euler',threshold='v>0*mV',refractory=3*ms)
    RE.v = -70*mV
    RE.m_na = 0.05
    RE.h_na = 0.95
    RE.m_t = 0.01
    RE.h_t = 0.01
    RE.n_k = 0.05
    RE.CA_i_RE = CA_inf
    
    #Recording
    V1 = StateMonitor(RE,('v','drive','ratio','m_na','h_na','n_k','Iext','CA_i_RE','I_k','I_kl','I_l','I_na','I_t','m_t','h_t','E_ca','m_tinf','h_tinf','tau_h_t','tau_m_t'),record=True)
    
    #Test
    RE.Iext = 0*pA
    run(100*ms)
    store('initial')
    RE.Iext = 0*pA
    run(200*ms)
    RE.Iext = 0*nA
    run(500*ms)#visualize

    fig,ax = subplots(6,1, sharex = True)
    ax[0].plot(V1.t/ms, V1.m_t[0], label='m_t')
    ax[0].plot(V1.t/ms, V1.h_t[0], label='h_t')
    ax[0].plot(V1.t/ms, V1.n_k[0], label='n_k')
    ax[0].legend(loc='upper left')
    ax[0].set_ylabel ('Activation')
    ax[1].plot(V1.t/ms, V1.v[0]/mV, label='v')
    ax[1].legend()
    ax[1].set_ylabel ('Voltage (mV)')
    ax[1].set_title('Membrane potential',fontsize=10)
    ax[2].plot(V1.t/ms, V1.Iext[0]/pA, label='Iext')
    ax[2].legend(loc='upper left')
    ax[2].set_ylabel ('Current (pA)')
    ax[2].set_title('External current applied',fontsize=10)
    ax[3].plot(V1.t/ms, V1.CA_i_RE[0]/uM)
    ax[3].set_ylabel ('Concentration (mM)')
    ax[3].set_title('Calcium concentration',fontsize=10)
    ax[4].plot(V1.t/ms, V1.E_ca[0]/mV)
    ax[4].set_title('Calcium reversal potential',fontsize=10)
    ax[4].set_ylabel ('Voltage (mV)')
    ax[5].plot(V1.t/ms, V1.m_tinf[0], label='m_tinf')
    ax[5].plot(V1.t/ms, V1.h_tinf[0], label='h_tinf')
    ax[5].set_xlabel ('Time (ms)')
    ax[5].legend()
    figure()
    plot(V1.t/ms, V1.I_l[0]/nA, label='I_l')
    plot(V1.t/ms, V1.I_kl[0]/nA, label='I_kl')
    plot(V1.t/ms, V1.I_t[0]/nA, label='I_t')
    #plot(V1.t/ms, V1.I_na[0]/nA, label='I_na')
    #plot(V1.t/ms, V1.I_k[0]/nA, label='I_k')
    legend()
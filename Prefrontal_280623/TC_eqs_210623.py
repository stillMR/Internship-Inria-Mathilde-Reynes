#!/usr/bin/env python
# coding: utf-8

# ## Required packages

# In[1]:


from brian2 import *
import numpy as np


# ## Global parameters 

# In[2]:


#Membrane capacitance per 2unit of surface
Cm_TC = 1*ufarad/cm**2

#Conductances
g_kl_TC = 0.03*msiemens*cm**-2
g_na_TC = 90*msiemens*cm**-2
g_t = 2.3*msiemens*cm**-2
g_l_TC = 0.01*msiemens*cm**-2
g_k_TC = 12*msiemens*cm**-2
g_h = 0.017*msiemens*cm**-2
g_a = 0*msiemens*cm**-2 

#Reversal potentials
E_kl = -95*mV
E_l_TC = -70*mV
E_h = -40*mV
E_na = 50*mV
E_k = -95*mV
E_ca0 = 1000*8.31441*(273.15 + 36)/(2*96489)*mV #13.31*mV approx

#Calcium parameters
tau_CA_TC = 5*ms
A_TC = 5.1819E-5*(mM*cm**2)/(ms*uA)
CA_inf = 2.4E-4*mM #2.4E-4
CA_0 = 2*mM #unit was found in Vijayan and Kopell 10.1073/pnas.1215385109

#Temperature-dependent variables 
T = 36
Qm_TC = 3.55**((T-24)/10)
Qh_TC = 3**((T-24)/10)
Qhyp = pow(3,((T-36)/10))
Qt = 2.9529
Q = 2.3
Tad = pow(3,((T-23.5)/10))

#Rates for open and close channels dynamics
k = 2.2
k1 = 7.9012E7*(mM**-4)*ms**-1
k2 = 0.0004*ms**-1 #0.004 in Wei et al.
k3 = 0.1*ms**-1
k4 = 0.001*ms**-1

#k = 2
#k1 = 7.9012E7*mM**-4*usecond**-1
#k2 = 0.004*usecond**-1 
#k3 = 0.1*usecond**-1
#k4 = 0.001*usecond**-1

#ModelDB parameters, Bazhenov et al., 2002
cac=0.002*mM
pc=0.01

#Parameters to match modelDB
Vtr = -50*mV
VtrK = -50*mV


# ## Equations

# In[3]:


TC_eqs = '''
    
    
    dv/dt = (- I_kl - I_na - I_k - I_t - I_l - I_h - I_a - Isyn_TC + Iext + Idummy) * (1/Cm_TC) : volt
    v2 = v - Vtr : volt
    v2K = v - VtrK : volt
    
    I_kl = g_kl_TC * (v - E_kl) : amp * meter**-2
    
    
    I_l = g_l_TC * (v - E_l_TC) : amp * meter**-2
    

    I_na = g_na_TC * (m_na ** 3) * h_na * (v - E_na) : amp * meter**-2 
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


    I_k = g_k_TC * (n_k ** 4) * (v - E_k) : amp * meter**-2
        dn_k/dt = -(n_k - n_kinf) / tau_n_k : 1
        
        alphan_k = 0.032/ms * (15*mV - v2K)/mV / (exp((15*mV - v2K)/5/mV) - 1) : Hz
        betan_k = 0.5/ms * exp((10*mV - v2K)/40/mV) : Hz
        n_kinf = alphan_k / (alphan_k + betan_k) : 1
        tau_n_k = (1 / (alphan_k + betan_k)) / Qhyp : second   


    I_t = g_t * (m_t ** 2) * h_t * (v - E_ca) : amp * meter**-2
        dm_t/dt = -(m_t - m_tinf) / tau_m_t : 1
        dh_t/dt = -(h_t - h_tinf) / tau_h_t : 1
        
        tau_m_t = ((1 / (exp(-(v + 131.6*mV)/16.7/mV) + exp((v + 16.8*mV)/18.2/mV)) + 0.612)/Qm_TC)*ms : second 
        m_tinf = 1 / (1 + exp(-(v + 59*mV)/6.2/mV)) : 1
        
        tau_h_t = ((30.8 + (211.4 + exp((v + 115.2*mV)/5/mV)) / (1 + exp((v + 86*mV)/3.2/mV)))/Qh_TC)*ms : second
        h_tinf = 1 / (1 + exp((v + 83*mV)/4/mV)) : 1
                      
        drive = (-A_TC * I_t)/2 : katal * meter**-3
        dCA_i_TC/dt = (drive + (CA_inf - CA_i_TC)/tau_CA_TC) * int(drive > 0*katal*meter**-3) + (0*katal*meter**-3 + (CA_inf - CA_i_TC)/tau_CA_TC) * int(drive < 0*katal*meter**-3) : mM 
        ratio = CA_0/CA_i_TC : 1
        E_ca = E_ca0 * log(ratio) : volt
   
   
   I_h = g_h * (Op + k*Op_L) * (v - E_h) : amp * meter**-2
        dOp/dt = alpha_Op * (1 - Op - Op_L) - beta_Op*Op : 1
        dP1/dt = k1ca*(1 - P1) - k2*P1 : 1
        dOp_L/dt = k3p*Op*int((Op_L+Op)<1) - k4*Op_L : 1
       
        tau_Op = (20*ms + 1000*ms / (exp((v + 71.5*mV)/(14.2*mV)) + exp(-(v + 89*mV)/(11.6*mV))))/Qhyp : second
        m_Opinf = 1 / (1 + exp((v + 75*mV)/5.5/mV)) : 1
        
        alpha_Op = m_Opinf / tau_Op : Hz
        beta_Op = (1 - m_Opinf) / tau_Op : Hz
        
        k1ca = k2 * (CA_i_TC/cac)**4 : Hz
        k3p = k4 * (P1/pc) : Hz
        
        
    I_a = g_a * (m_a **4) * h_a * (v - E_k) : amp * meter**-2
        dm_a/dt = -(1/tau_m_a)*(m_a - m_ainf) : 1
        dh_a/dt = -(1/tau_h_a)*(h_a - h_ainf) : 1
        
        tau_m_a = (1*ms / (exp((v + 35.82*mV)/19.69/mV) + exp(-(v + 79.69*mV)/12.7/mV)) + 0.37*ms) / Tad : second
        m_ainf = 1 / (1 + exp(-(v + 60*mV)/8.5/mV)) : 1
        
        tau_h_a = (19*ms/Tad) * int(v >= -63*mV) + (1*ms / ((exp((v + 46.05*mV)/5/mV) + exp(-(v + 238.4*mV)/37.45/mV))) / Tad) * int(v < -63*mV) : second
        h_ainf = 1 / (1 + exp((v + 78*mV)/6/mV)) : 1

  
    Isyn_TC = IsynGABAA_RE_TC + IsynGABAB_RE_TC + IsynAMPA_PY_TC_layer_VI + IsynAMPA_PY_TC_layer_V : amp * meter**-2
        IsynGABAA_RE_TC : amp * meter**-2
        IsynGABAB_RE_TC : amp * meter**-2
        IsynAMPA_PY_TC_layer_VI :  amp * meter**-2
        IsynAMPA_PY_TC_layer_V :  amp * meter**-2
    
    
    Iext : amp * meter**-2
    
    Idummy = IEPSPs_L_VI_TCTCo + IEPSPs_L_V_TCTCo : amp * meter**-2
        IEPSPs_L_VI_TCTCo : amp * meter**-2
        IEPSPs_L_V_TCTCo : amp * meter**-2
    
    x : meter
    
    
    '''


if __name__=='__main__':
    # ## Running simulation and plotting results
    prefs.codegen.target = "cython"
    defaultclock.dt = 0.02*ms
    
    #Parameters initialization
    TC = NeuronGroup(1,TC_eqs,method='exponential_euler',threshold='v>0*mV',refractory=3*ms)
    TC.v = -75*mV
    TC.m_na = 0.05
    TC.h_na = 0.95
    TC.m_t = 0.01
    TC.h_t = 0.01
    TC.n_k = 0.01
    TC.CA_i_TC = CA_inf 
    TC.P1 = 0.01
    TC.Op = 0.01
    TC.Op_L = 0.01
    
    #Recording
    V1 = StateMonitor(TC,('v','k1ca','tau_Op','E_ca','m_na','h_na','m_t','n_k','Iext','CA_i_TC','I_k','I_a','I_t','I_kl','I_na','I_h','Op','Op_L','P1','I_l'),record=True)
    
    #Test
    TC.Iext = 0*nA
    run(500*ms)
    store('initial')
    TC.Iext = 0*pA
    run(500*ms)
    TC.Iext = 0*nA
    run(2000*ms)  #visualize
    
    fig,ax = subplots(5,1, sharex = True)
    ax[0].plot(V1.t/ms, V1.v[0]/mV, label='v')
    ax[0].legend()
    ax[0].set_ylabel ('Voltage (mV)')
    ax[0].set_title('Membrane potential',fontsize=10)
    ax[1].plot(V1.t/ms, V1.Iext[0]/pA, label='I_ext')
    ax[1].legend(loc='upper left')
    ax[1].set_ylabel ('Current (pA)')
    ax[1].set_title('External current applied',fontsize=10)
    ax[2].plot(V1.t/ms, V1.Op[0])
    ax[2].set_title('Fully Open channels: Op',fontsize=10)
    ax[3].plot(V1.t/ms, V1.Op_L[0])
    ax[3].set_title('Open Locked channels: Op_L',fontsize=10)
    ax[4].plot(V1.t/ms, V1.P1[0])
    ax[4].set_title('Regulating factor: P1',fontsize=10)
    ax[4].set_xlabel ('Time (ms)')
    #ax[5].plot(V1.t/ms, V1.k1ca[0]/ms**-1)
    #ax[5].set_title('k1ca',fontsize=10)
    #ax[5].set_xlabel ('Time (ms)')
    figure()
    plot(V1.t/ms, V1.I_l[0]/nA, label='I_l')
    plot(V1.t/ms, V1.I_kl[0]/nA, label='I_kl')
    plot(V1.t/ms, V1.I_t[0]/nA, label='I_t')
    #plot(V1.t/ms, V1.I_na[0]/nA, label='I_na')
    #plot(V1.t/ms, V1.I_k[0]/nA, label='I_k')
    plot(V1.t/ms, V1.I_h[0]/nA, label='I_h')
    plot(V1.t/ms, V1.I_a[0]/nA, label='I_a')
    legend()
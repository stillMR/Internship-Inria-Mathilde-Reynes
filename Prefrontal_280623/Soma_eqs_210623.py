#!/usr/bin/env python
# coding: utf-8

# ## Required packages

# In[1]:


from brian2 import *


# ## Global parameters

# In[2]:


#Area
s_soma = 10**-6*cm**2

#Ion-specific conductances per unit of surface
# g_k = 200*msiemens*cm**-2

#Time constants
tau_m_nap = 0.1991*ms

#Ion-specific reversal potential
E_na = 50*mV
E_k = -90*mV

#Other constants 
Qt = pow(2.3,((36-23)/10))
Tad = pow(2.3,((36-23)/10))
Phi_m = pow(2.3,((36-23)/10))

#somatic and dendritic compartments
# R = 10 * Mohm
R = 100 * Mohm
g_ds = 1 / (R * s_soma)


# ## Equations

# In[3]:


Soma_eqs = '''
    v = -70*mV*init_timestep + truev*(1-init_timestep) : volt
    init_timestep = init_timedarray(t) : 1
    
    truev = (vgap + (R * s_soma) * g2) / (1 + (R * s_soma) * g1) : volt    
                            
    g1 = g_na *(m_na ** 3) * h_na + g_k *  n_k + g_nap *  m_nap : siemens * meter**-2
    g2 = g_na *(m_na ** 3) * h_na * E_na + g_k * n_k * E_k + g_nap * m_nap * E_na + Iext : amp*meter**-2
   
    I_na = Phi_m*g_na * (m_na ** 3) * h_na * (v - E_na) : amp * meter**-2
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
    
    I_k = Tad*g_k * n_k * (v - E_k) : amp * meter**-2
        dn_k/dt = -(n_k - n_kinf) / tau_n_k : 1
        alphan_k = (0.02/mV) * (v - 25*mV) / (1 - exp(-(v - 25*mV)/(9*mV))) : 1
        betan_k = (-0.002/mV) * (v - 25*mV) / (1 - exp((v - 25*mV)/(9*mV))) : 1
        n_kinf = alphan_k / (alphan_k + betan_k) : 1
        tau_n_k = (1*msecond / (alphan_k  + betan_k)) / Qt : second
        
    vgap : volt
     
    g_na : siemens * meter**-2
    g_k : siemens * meter**-2
    g_nap : siemens * meter**-2
    
    x : meter
    
    Iext : amp * meter**-2
    
    t_last_spike_PY : second
    
    mean_rate_PY = (2/(1 + exp(-(t-t_last_spike_PY)/(50*ms)))-1)/250*kHz : Hz
    
    mean_rate_IN = (2/(1 + exp(-(t-t_last_spike_PY)/(50*ms)))-1)/250*kHz : Hz
    
    mean_rate_inter = (2/(1 + exp(-(t-t_last_spike_PY)/(80*ms)))-1)/250*kHz : Hz
    

        
    '''

if __name__=='__main__' : 
    close('all')
    prefs.codegen.target = "numpy"
    
    inpext=0*amp*meter**-2 
    vext=10*pA*R
    
    defaultclock.dt = 0.02*ms
    
    #parameters initialization
    PY_soma = NeuronGroup(1,Soma_eqs,method='rk4',threshold='v>0*mV',refractory=3*ms)
    #PY_soma.v = -70*mV
    PY_soma.h_na = 0.92
    PY_soma.m_na = 0.01
    PY_soma.m_nap = 0
    PY_soma.n_k = 0
    PY_soma.g_na = 3000*msiemens*cm**-2
    PY_soma.g_nap = 15*msiemens*cm**-2
    PY_soma.g_k = 200*msiemens*cm**-2
    # PY_soma.g_k = 400*msiemens*cm**-2
    # PY_soma.g_na = 3000*msiemens*cm**-2
    # PY_soma.g_nap = 0*msiemens*cm**-2
    # PY_soma.g_k = 150*msiemens*cm**-2
    
    #recording
    V1 = StateMonitor(PY_soma,('v','h_na','m_na','m_nap','n_k','tau_m_na','tau_h_na','tau_n_k','g1','g2','Iext','alpham_na','betam_na','alphah_na','betah_na','alphan_k','betan_k','I_na','I_nap','I_k'),record=True)
    
    runtime = 400*ms
    num_samples = int(runtime/defaultclock.dt)
    init_arr = zeros(num_samples)
    init_arr[0]=1
    init_timedarray = TimedArray(init_arr, dt=defaultclock.dt)
    
    # PY_soma.Iext = 0*pA*meter**-2
    # run(10*ms)
    # store('initial')
    # # PY_soma.Iext = 10*pA*meter**-2
    # PY_soma.Iext = inpext
    # run(200*ms)
    # PY_soma.Iext = 0*pA*meter**-2
    # run(15*ms)# %matplotlib widget
    
    PY_soma.vgap = 0*mV
    run(100*ms)
    store('initial')
    # PY_soma.Iext = 10*pA*meter**-2
    PY_soma.vgap = 40*mV
    run(200*ms)
    PY_soma.vgap = 0*mV
    run(100*ms)# %matplotlib widget
    
    fig,ax = subplots(4,1, sharex = True)
    ax[0].plot(V1.t/ms, V1.m_na[0], label='m_na')
    ax[0].plot(V1.t/ms, V1.h_na[0], label='h_na')
    ax[0].plot(V1.t/ms, V1.n_k[0], label='n_k')
    ax[0].plot(V1.t/ms, V1.m_nap[0], label='m_nap')
    ax[0].legend(loc='upper right')
    ax[0].set_ylabel ('Activation')
    ax[1].plot(V1.t/ms, V1.v[0]/mV, label='v')
    ax[1].legend(loc='upper right')
    ax[1].set_ylabel ('Voltage (mV)')
    ax[2].plot(V1.t/ms, V1.alpham_na[0], label='alpham_na')
    ax[2].plot(V1.t/ms, V1.alphah_na[0], label='alphah_na')
    ax[2].plot(V1.t/ms, V1.betam_na[0], label='betam_na')
    ax[2].plot(V1.t/ms, V1.betah_na[0], label='betah_na')
    ax[2].plot(V1.t/ms, V1.alphan_k[0], label='alphan_k')
    ax[2].plot(V1.t/ms, V1.betan_k[0], label='betan_k')
    ax[2].legend()
    ax[3].plot(V1.t/ms, V1.Iext[0]/pA, label='Iext')
    ax[3].legend()
    ax[3].set_ylabel ('Current (pA)')
    ax[3].set_xlabel ('Time (ms)')
    
    figure()
    plot(V1.t/ms, V1.I_na[0], label='I_na')
    plot(V1.t/ms, V1.I_nap[0], label='I_nap')
    plot(V1.t/ms, V1.I_k[0], label='I_k')
    legend()
    
    vv=[-100*mV+i*mV for i in range(130)]
    alpham_na = [(0.182/ms * (v + 25*mV)/mV / (1 - exp(-(v + 25*mV)/9/mV))) * int(abs((v-10*mV)/(-35*mV)) > 1e-6) + (0.182/ms * 9) * int(abs((v-10*mV)/(-35*mV)) < 1e-6) for v in vv]
    betam_na = [(-0.124/ms * (v + 25*mV)/mV / (1 - exp((v + 25*mV)/9/mV))) * int(abs((-v+10*mV)/(35*mV)) > 1e-6) + (0.124/ms * 9) * int(abs((-v+10*mV)/(35*mV)) < 1e-6) for v in vv]

    h_nainf = [1 / (1 + exp((v + 55*mV)/6.2/mV)) for v in vv]
    m_nainf = [alpham_na[i] / (alpham_na[i] + betam_na[i]) for i in range(len(vv))]
    
    m_napinf = [0.02 / (1 + exp(-(v + 42*mV)/5/mV)) for v in vv]
    
    alphan_k = [(0.02/mV) * (v - 25*mV) / (1 - exp(-(v - 25*mV)/(9*mV))) for v in vv]
    betan_k = [(-0.002/mV) * (v - 25*mV) / (1 - exp((v - 25*mV)/(9*mV))) for v in vv]
    n_kinf = [alphan_k[i] / (alphan_k[i] + betan_k[i]) for i in range(len(vv))]
    
    g2inf=[PY_soma.g_na*m_nainf[i]**3*h_nainf[i]*E_na+PY_soma.g_k*n_kinf[i]*E_k+PY_soma.g_nap*m_napinf[i]*E_na for i in range(len(vv))]
    g1inf=[PY_soma.g_na*m_nainf[i]**3*h_nainf[i]+PY_soma.g_k*n_kinf[i]+PY_soma.g_nap*m_napinf[i] for i in range(len(vv))]
    
    figure()
    plot(vv,m_nainf,label='m_na')
    plot(vv,h_nainf,label='h_na')
    plot(vv,m_napinf,label='m_nap')
    plot(vv,n_kinf,label='n_k')
    legend()
    
    figure()
    plot(vv/mV,g2inf)
    
    figure()
    plot(vv/mV,g1inf)
    
    figure()
    for vext in [-200*mV+50*k*mV for k in range(10)]:
        plot(vv,[(vext+1/g_ds*g2inf[i])/(1+1/g_ds*g1inf[i]) for i in range(len(vv))],label=str(vext))
    plot(vv,vv,label='vv=vv')
    legend()
  
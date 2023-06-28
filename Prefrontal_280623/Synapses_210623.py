#!/usr/bin/env python
# coding: utf-8


# Synapses.ipynb contains functions to create all the different synapses part of the model: GABAa, GABAb, AMPA and NMDA for both intra and intercortical interactions.
# ## GABAa: with short-term depression

# In[1]:

from brian2 import *

def syn_gabaa(source,target,syntype,surfacetarget,connection_pattern,g_syn,E_syn,alpha,beta,U):
    eq_syn_GABAA = '''_post = D * g_syn * (1/surfacetarget) * W * (v_post - E_syn) : amp * meter**-2 (summed)
        dW/dt = alpha * (1 - W) * T - beta * W : 1 (clock-driven)
        T = A_syn * int((t_last_spike + tmax - t) > 0*ms)* int((t - t_last_spike) > 0*ms) : 1
        D : 1
        g_syn : siemens
        E_syn : volt
        t_last_spike : second
        alpha : Hz
        beta : Hz
        U : 1
        tau = 700*ms : second
        tmax = 0.3*ms : second
        A_syn = 0.5 : 1
        surfacetarget : meter**2
    '''
    pre_code_gabaa = ''' D = 1 - (1 - D*(1 - U)) * exp(-(t - t_last_spike)/tau)
        t_last_spike = t
    '''
    
    S=Synapses(source,target,model=syntype+eq_syn_GABAA,method='rk4',on_pre=pre_code_gabaa)
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(condition=connection_pattern, skip_if_invalid=True)
    S.surfacetarget = surfacetarget
    S.g_syn = g_syn
    S.E_syn = E_syn
    S.alpha = alpha
    S.beta = beta
    S.U = U
    return S


# ## AMPA: with short term depression, Poisson inputs and STDP rule

# In[ ]:


def syn_ampa(source,target,syntype,surfacetarget,connection_pattern,g_syn_max,g_syn,E_syn,alpha,beta,U,A_mEPSP0,poissontype,connection_pattern_mini):
    eq_syn_AMPA = '''_post = D * g_syn * (1/surfacetarget) * W * (v_post - E_syn) : amp * meter**-2 (summed)
        dW/dt = alpha * (1 - W) * T - beta * W : 1 (clock-driven)
        T = A_syn * int((t_last_spike + tmax - t) > 0*ms)* int((t - t_last_spike) > 0*ms) : 1
        D : 1
        g_syn : siemens
        g_syn_max : siemens
        E_syn : volt
        t_last_spike : second
        alpha : Hz
        beta : Hz
        U : 1
        tau = 700*ms : second
        tmax = 0.3*ms : second
        A_syn = 0.5 : 1
        w : 1
        dapre/dt = -apre/taupre : 1 (clock-driven)
        dapost/dt = -apost/taupost : 1 (clock-driven)
        taupre = 20*ms : second
        taupost = 20*ms : second
        Apre : 1
        Apost : 1
        gmax = 3 * g_syn_max : siemens 
        dummy = 0*usiemens : siemens
        connection_pattern_mini : 1
        surfacetarget : meter**2
    '''
    eq_syn_AMPA2 = '''_post = D * A_mEPSP * W_Poisson * (v_post - E_syn) * (1/surfacetarget) : amp * meter**-2 (summed)
        dW_Poisson/dt = alpha * (1 - W_Poisson) * T_Poisson - beta * W_Poisson : 1 (clock-driven)
        T_Poisson = (A_syn * int((t_last_spike_Poisson_PY + tmax - t) > 0*ms)* int((t - t_last_spike_Poisson_PY) > 0*ms))*int(connection_pattern_mini == 0)  + (A_syn * int((t_last_spike_Poisson_IN + tmax - t) > 0*ms)* int((t - t_last_spike_Poisson_IN) > 0*ms))*int(connection_pattern_mini == 1) + (A_syn * int((t_last_spike_Poisson_inter + tmax - t) > 0*ms)* int((t - t_last_spike_Poisson_inter) > 0*ms))*int(connection_pattern_mini == 2) : 1
        t_last_spike_Poisson_PY : second
        t_last_spike_Poisson_IN : second
        t_last_spike_Poisson_inter : second
        f = 0.01 : 1
        A_mEPSP0 : siemens
        A_mEPSP : siemens
    '''
    pre_code_ampa = '''D = 1 - (1 - D*(1 - U)) * exp(-(t - t_last_spike)/tau)
        w = apost
        Apre = 0.0005
        apre += Apre
        t_last_spike = t
        t_last_spike_PY_pre = t
        g_syn = clip((g_syn+w*gmax), dummy, gmax)
        
        A_mEPSP = A_mEPSP + f * w *  A_mEPSP0
    '''
    post_code_ampa = '''
        w = apre
        Apost = -0.0005
        apost += Apost
        g_syn = clip((g_syn+w*gmax), dummy, gmax)
        
        A_mEPSP = A_mEPSP + f * w *  A_mEPSP0
    '''
    poisson_PY = '''
        t_last_spike_Poisson_PY = t
    '''
    poisson_IN = '''
        t_last_spike_Poisson_IN = t
    '''
    poisson_inter = '''
        t_last_spike_Poisson_inter = t
    '''
    S=Synapses(source,target,model=(syntype+eq_syn_AMPA)+(poissontype+eq_syn_AMPA2),method='rk4',on_pre={'pre': pre_code_ampa,'mini_PY': poisson_PY,'mini_IN': poisson_IN,'mini_inter': poisson_inter},on_post=post_code_ampa,on_event={'mini_PY':'custom_poisson_PY','mini_IN':'custom_poisson_IN','mini_inter':'custom_poisson_inter'})
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(condition=connection_pattern, skip_if_invalid=True)
    S.surfacetarget = surfacetarget
    S.g_syn = g_syn
    S.g_syn_max = g_syn_max
    S.E_syn = E_syn
    S.alpha = alpha
    S.beta = beta
    S.U = U
    S.A_mEPSP0 = A_mEPSP0
    S.connection_pattern_mini = connection_pattern_mini
    return S


# ## NMDA: with short-term depression and v_post relationship

# In[3]:


def syn_nmda(source,target,syntype,surfacetarget,connection_pattern,g_syn,E_syn,alpha,beta):
    eq_syn_NMDA = '''_post = g_syn * (1/surfacetarget) * W * f * (v_post - E_syn) : amp * meter**-2 (summed)
        dW/dt = alpha * (1 - W) * T - beta * W : 1 (clock-driven)
        T = A_syn * int((t_last_spike + tmax - t) > 0*ms)* int((t - t_last_spike) > 0*ms) : 1
        f = 1/(1 + exp(-(v_post - v_th)/(sigma))) : 1
        g_syn : siemens
        E_syn : volt
        t_last_spike : second
        alpha : Hz
        beta : Hz
        tmax = 0.3*ms : second
        A_syn = 0.5 : 1
        v_th = -25*mV : volt
        sigma = 12.5*mV : volt
        surfacetarget : meter**2
    '''
    pre_code_nmda = '''
        t_last_spike = t
    '''
    
    S=Synapses(source,target,model=syntype+eq_syn_NMDA,method='rk4',on_pre=pre_code_nmda)
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(condition=connection_pattern, skip_if_invalid=True)
    S.surfacetarget = surfacetarget
    S.g_syn = g_syn
    S.E_syn = E_syn
    S.alpha = alpha
    S.beta = beta
    return S


# ## GABAb with fraction of activated receptors and G-proteins concentration

# In[4]:


def syn_gabab(source,target,syntype,surfacetarget,connection_pattern):
    eq_syn_GABAB = '''_post = g_gabab * (1/surfacetarget) * (G**4 / (G**4 + K)) * (v - E_k) : amp * meter**-2 (summed)
        dR_gab/dt = K1 * (1*uM - R_gab)*T - K2*R_gab : mM (clock-driven)
        dG/dt = K3*R_gab - K4*G : mM (clock-driven)
        T = A_syn * int((t_last_spike + tmax - t) > 0*ms)* int((t - t_last_spike) > 0*ms) : mM
        g_gabab = 0.04*usiemens : siemens
        E_k = -95*mV : volt
        K = 100*uM**4 : mM**4
        K1 = 0.52*ms**-1*mM**-1 : second**-1*mM**-1
        K2 = 0.0013*ms**-1 : Hz
        K3 = 0.098*ms**-1 : Hz
        K4 = 0.033*ms**-1 : Hz
        A_syn = 0.5*mM : mM
        t_last_spike : second
        tmax = 0.3*ms : second
        surfacetarget : meter**2
    '''
    pre_code_gabab = '''
        t_last_spike = t
    '''
    S=Synapses(source,target,model=syntype+eq_syn_GABAB,method='rk4',on_pre=pre_code_gabab)
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(condition=connection_pattern, skip_if_invalid=True)
    S.surfacetarget = surfacetarget
    return S


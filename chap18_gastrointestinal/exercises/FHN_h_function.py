import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.spines.right'] = False  # Me being pedantic and removing the box around the graph
plt.rcParams['axes.spines.top'] = False

def get_v_w(t):   # For any t, find V and W on the stable periodic solution
    tt = np.remainder(t,T)
    dum = t_n-tt
    dum_arg = np.argmax(dum>0)
    tup = t_n[dum_arg]
    tdown = t_n[dum_arg-1]
    vup = v_n[dum_arg]
    vdown = v_n[dum_arg-1]
    v = vdown + ((tt-tdown)/(tup-tdown))*(vup-vdown) # linear interpolation to find V(t)
    wup = w_n[dum_arg]
    wdown = w_n[dum_arg-1]
    w= wdown + ((tt-tdown)/(tup-tdown))*(wup-wdown) # linear interpolation to find W(t)
    return v,w

def get_h(phi):
    vshift = np.empty(n)
    for i in range(n):
        vshift[i],_ = get_v_w(phi+t_n[i])
    integrand = Va*vshift
    h = np.trapz(integrand,t_n)/T
    return h
        


# First, read in all the data from the output of FHN_adjoint.py. The period,
# limit cycle, adjoint solution, parameters, etc. So we don't need to recompute
# all this stuff

[alpha,eps,I,T,t_n,v_n,w_n,Va,Wa]= np.load('FHN_save_variables.npy',allow_pickle=True)
n = len(t_n)

sn = 100
h = np.empty(sn)
hm = h
st = np.linspace(0,T,sn)
for i in range(sn):
    h[i] = get_h(st[i])
    hm = get_h(-st[i])
    
plt.plot(st/T,-h+hm)
plt.xlabel('time/T')
plt.ylabel(r'h$(-\phi)-h(\phi)$')




  

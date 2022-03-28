import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
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
        
def FHN_coupled_deriv(t, y):
    v1,w1,v2,w2 = y 
    v1dot = (1/eps)*( v1*(1-v1)*(v1-alpha) - w1) + I + delta*(v2-v1)
    w1dot = v1 - w1
    v2dot = (1/eps)*( v2*(1-v2)*(v2-alpha) - w2) + I + delta*(v1-v2)
    w2dot = v2 - w2
    return [v1dot,w1dot,v2dot,w2dot]

# First, read in all the data from the output of FHN_adjoint.py. The period,
# limit cycle, adjoint solution, parameters, etc. So we don't need to recompute
# all this stuff

[alpha,eps,I,T,t_n,v_n,w_n,Va,Wa]= np.load('FHN_save_variables.npy',allow_pickle=True)
n = len(t_n)
delta = 0.3   # the coupling strength

sn = 100
h = np.empty(sn)
st = np.linspace(0,T,sn)
for i in range(sn):
    h[i] = get_h(st[i])

plt.figure(1)
plt.plot(st/T,h)
plt.xlabel('$\phi$/T')
plt.ylabel(r'$h(\phi)$')

plt.figure(2)
plt.plot(st/T,np.flip(h) - h)
plt.xlabel('$\phi$/T')
plt.ylabel(r'$h(-\phi)-h(\phi)$')

# Now solve the coupled oscillator equations directly
tend = 10*T
t = np.linspace(0,tend,2000)
y0 = [0.1,0,-0.3,0.1]
soln = solve_ivp(FHN_coupled_deriv, (0, tend), y0, method='Radau',t_eval=t,rtol=1e-8,atol=1e-8)
v1,w1,v2,w2 = soln.y

plt.figure(3)
plt.plot(soln.t,v1,soln.t,v2)
plt.xlim([tend-3*T,tend])

  




#

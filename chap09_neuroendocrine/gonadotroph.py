# -----------------------------------
# 
#  Pituitary cell model using a common core of ionic currents.
#  Adding ER calcium, a cytosolic calcium subspace, and IP3 receptor dynamics.
#  Adapted from the original code (Fletcher, Bertram, Stojilkovic, Mol Cell Endocrin 463 (2018))
#  for Keener and Sneyd, Math Physiology, third edition.
# 
#  ----------------------------------

from numpy import exp, arange
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# The equations
def deRHS(sol, t):
    global ip3b, ip3p, gnav, gca, gk, gkca, gL, gkir, ek, eca, ena, vL, alpha, Cm, sigmaER, sigmaD
    global tauhna, taun, kca, fcyt, Kdiff, kleak, kf, ki, ka, Kd, Aip3, fer, Vpmca, Kpmca, Vserca, Kserca
    
    # There are 7 variables
    v, c, cer, cd, hna, n, hip3r = sol
    
    vKdrive = v - ek
    vnadrive = v - ena
    vcadrive = v - eca
    
    # Sodium
    mnainf = 1 / (1 + exp((-15 - v) / 5))
    hnainf = 1 / (1 + exp((55 + v) / 10))
    ina = gnav * mnainf**3 * hna * vnadrive
    
    # Calcium
    mcainf = 1 / (1 +exp((-15 - v) / 12))
    ica = gca * mcainf * vcadrive
    
    # Potassium
    ninf = 1 / (1 + exp((-v) / 5))
    iKdr = gk * n * vKdrive
    
    c2 = c**2
    nkcainf = c2 / (c2 + kca**2)
    ikca = gkca * nkcainf * vKdrive
    
    nkirinf = 1 / (1 + exp((55 + v) / 10))
    ikir = gkir * nkirinf * vKdrive
    
    ik = ikca + iKdr + ikir
    
    # Non-specific leak
    iL = gL * (v - vL)
    
    # Calcium handling
    jin = -alpha * ica
    
    jpmca = Vpmca * c2 / (c2 + Kpmca**2)
    jserca = Vserca * c2 / (c2 + Kserca**2)
    
    ip3 = ip3b + ip3p * ip3pulse(t)
    
    num = (ip3 * cd * hip3r)
    num3 = num**3
    denom = (ip3 + ki) * (cd + ka)
    denom3 = denom**3
    jip3r = kf * num3 / denom3 * (cer - cd)
    
    jleak = kleak * (cer - cd)
    jrel = jleak + jip3r
    
    jdiff = Kdiff * (cd - c)
    
    # Equations
    vp = -(ina + ica + ik + iL) / Cm
    cp = fcyt * (jin - jpmca + jdiff - jserca)
    cerp = fer * sigmaER * (jserca - jrel)
    cdp = fcyt * sigmaD * (jrel - jdiff)
    hnap = (hnainf - hna) / tauhna
    np = (ninf - n) / taun
    hip3rp = Aip3 * (Kd - (cd + Kd) * hip3r)
    
    return [vp, cp, cerp, cdp, hnap, np, hip3rp]


# IP3 pulse
def ip3pulse(t):
    global tpulse, pnorm, kp
    delpulse = t - tpulse
    pulse = delpulse**2 / pnorm * (delpulse > 0)
    ip3pulse = pulse / (pulse + kp) * (1 + kp)
    return ip3pulse

# Global parameters
ip3b = 0
ip3mx = [0.1, 0.4, 0.7, 2]
tpulse = 5000
gnav = 20
gca = 1.5
gk = 5
gkca = 2.5
gL = 0.1
gkir = 0.9
ek = -75
eca = 60
ena = 75
vL = 0
alpha = 0.0015
Cm = 6
sigmaER = 39
sigmaD = 39
kp = 0.04
tauhna = 2
taun = 20
tau = 175000
pnorm = 4 * tau**2 * exp(-2)
Vserca = 0.16
Kserca = 0.2
Kdiff = 0.6
kca = 0.4
fer = 0.01
fcyt = 0.01
Vpmca = 0.02
Kpmca = 0.1
kleak = 0.0002
kf = 20
ki = 1
ka = 0.8
Kd = 0.8
Aip3 = 0.002

for j in range(4):
    ip3p = ip3mx[j]
    # Initial conditions
    v0 = -60.58
    c0 = 0.0881
    cer0 = 124.49
    cd0 = 0.13
    hna0 = 0.637
    n0 = 0.0
    hip3r0 = 0.862
    init = [v0, c0, cer0, cd0, hna0, n0, hip3r0]
    
    total = 60000
    tstep = 1
    # Specify the output points
    tspan = arange(0, total + tstep, tstep)
    
    # Solve the ODE
    sol = odeint(deRHS, init, tspan)
    
    fig=plt.figure(2 * j + 1)
    plt.plot(tspan / 1000, sol[:, 0],label='c')
    plt.ylabel('V (mV)')
    plt.xlabel('t (s)')
    plt.twinx()
    plt.plot(tspan / 1000, ip3b + ip3p * ip3pulse(tspan),'r--',label=r'IP$_3$')
    plt.ylabel(r'IP$_3$')
    fig.legend()
    
    fig=plt.figure(2 * j + 2)
    plt.plot(tspan / 1000, sol[:, 1], label='c')
    plt.xlabel('t (s)')
    plt.ylabel(r'c ($\mu$M)')
    plt.twinx()
    plt.plot(tspan / 1000, sol[:, 2], 'r--', label=r'$c_e$')
    plt.ylabel(r'$c_e$ ($\mu$M)')
    fig.legend()






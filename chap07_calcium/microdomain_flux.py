import numpy as np
import matplotlib.pyplot as plt


def getc(c):
    global Jc
    Jc = fluxw(c, Dc, cinf)
    e = bisect(gete, c, einf)
    return rho * (e - c) - Jc

def gete(e):
    global Jc
    Je = fluxw(e, De, einf)
    return Je + Jc

def w(c, D):
    wout = D * c + Dbbt * c / (Kb + c)
    wpout = D + Dbbt * Kb / (Kb + c)**2
    return wout, wpout

def fluxw(c, D, cinf):
    winf, wpinf = w(cinf, D)
    w0, wp0 = w(c, D)
    return D * (w0 - winf) / wp0

def bisect(func, a, b, N=25):
    ul = a
    fl = func(ul)
    uu = b
    fu = func(uu)

    for j in range(N):
        u = (ul + uu) / 2
        fc = func(u)
        ftest = (fc * fl >= 0)
        ul = ftest * u + (1 - ftest) * ul
        fl = ftest * fc + (1 - ftest) * fl
        uu = (1 - ftest) * u + ftest * uu
        fu = (1 - ftest) * fc + ftest * fu

    return u

# Set plot parameters
plt.rcParams.update({
    'axes.titlesize': 20,
    'axes.linewidth': 2.0,
    'lines.linewidth': 2.0
})
plt.rcParams['text.usetex'] = True


# Parameters
Dc = 1
De = 1
Dbbt = 5

rho = 2
cinf = 0.5

elist = [1, 5, 10]
keep = []
Kblist = np.arange(0.01, 20.1, 0.1)

for einf in elist:
    J = []
    for Kb in Kblist:
        c = bisect(getc, cinf, einf)
        J.append(fluxw(c, Dc, cinf))

    plt.figure(1)
    plt.plot(Kblist, np.array(J) / (rho * (einf - cinf)))
    keep.append(np.column_stack((Kblist, np.array(J) / (rho * (einf - cinf)))))
    plt.xlabel('$K_b$')
    plt.ylabel('Permeability$/\rho$')

plt.text(2, 0.336, '$e_\infty=1$', fontsize=18)
plt.text(2, 0.362, '$e_\infty=5$', fontsize=18)
plt.text(8, 0.362, '$e_\infty=10$', fontsize=18)



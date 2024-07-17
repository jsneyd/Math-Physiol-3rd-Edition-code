import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def rhs(t, s, p):
    V = s[:N]
    w = s[N:]
    currents = iv(s, p)
    FV, Fw = currents[:, 0], currents[:, 1]
    Fv = dg * (-sc * V + np.concatenate([[0], V[:-1]]) + np.concatenate([V[1:], [0]])) + FV
    return np.concatenate([Fv, Fw])

def iv(s, p):
    u = s[:N]
    w = s[N:]
    Fv = (6 * u * (u - alpha) * (1 - u) - w + Iapp) / eps
    Fw = u - gamma * w
    return np.column_stack([Fv, Fw])


# Define parameters for FHN
Iapp = 0
alpha = 0.2
eps = 0.1
gamma = 0.5
N = 500  # number of grid points

L = 40
h = L / (N - 1)
dg = 1 / (h ** 2)  # Take D = 1, coupling (diffusion) coefficient
X = np.arange(0, N) * h
sc = np.concatenate([[1], 2 * np.ones(N - 2), [1]])

# Initial data
V0 = 1 / np.cosh(X)
u0 = np.concatenate([V0, np.zeros(N)])
tstep = 0.025
t_end = 20

# Specify the output points
tspan = np.arange(0, t_end + tstep, tstep)

sol = solve_ivp(rhs, [0, t_end], u0, t_eval=tspan, method='LSODA')

# Calculate the nullclines
u = np.arange(-0.3, 1.01, 0.01)
w1 = 6 * u * (u - alpha) * (1 - u) + Iapp
w2 = u / gamma

# Plotting
j = 250
plt.figure(1)
plt.plot(u, w1, '--', u, w2, '--', sol.y[j, :N], sol.y[j, N:], linewidth=2)
plt.axis([-0.3, 1, -0.2, 1])
plt.xlabel('v')
plt.ylabel('w')
plt.box(False)
plt.legend(['dv/dt=0', 'dw/dt=0'])

plt.figure(3)
X, T = np.meshgrid(X, sol.t)
plt.pcolormesh(X, T, sol.y[:, :N], shading='auto')
plt.xlabel('x')
plt.ylabel('t')
plt.colorbar(label='v')

# Now find the speed
thresh = 0.5
Tc = np.array([np.argmax(sol.y[:, j] >= thresh) * tstep for j in range(N)])
q = np.polyfit(Tc, X, 1)
speedest = q[0]
spest = q[1] + q[0] * Tc

plt.figure(4)
plt.plot(X, Tc, spest, '--')
plt.xlabel('x')
plt.ylabel('t')
plt.title(f'Speed = {speedest:.2f}')

plt.figure(2)
plt.plot(X, sol.y[j, :N], X, sol.y[j, N:])
plt.xlabel('x')
plt.box(False)
plt.legend(['v', 'w'])
plt.title(f'Speed = {speedest:.2f}')

plt.show()

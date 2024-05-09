
#  ------------------------------------------------
#  Responses of receptive fields to moving bars and steps.
# 
#  For Chapter 19, Section 19.6, as well as a number of exercises.
# 
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  Written by James Keener and James Sneyd
# 
#  ------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.integrate import quad, dblquad


def gaussian(x, s, g):
    return g * s * np.exp(-s**2 * x**2) / np.sqrt(np.pi)

def f(x):    # the receptive field as a sum of Gaussians
    return gaussian(x, s1, g1) - gaussian(x, s2, g2)

def transient(x, t): # transient response assuming a non-instantaneous response
    return f(x) * (t - x / c) * np.exp((x - t) / c)

def imp(t):
    return (1.2 * t * np.exp(-t) - 0.45 * t**2 * np.exp(-t)) * (t >= 0)

def step_response(t):
    return np.array([quad(lambda s: imp(t_ - s), 0, t_)[0]    for t_ in t])

def stim(x, t):    # difference of Heaviside functions
    return (t >= x / c) ^ (t >= (x + w) / c)

def bar_response(t, c):
    return np.array([quad(lambda s, x: stim(x, t_) * imp(t_ - s) * f(x), -np.inf, t_, -np.inf, np.inf)[0] for t_ in t])

# parameters
s1 = 3
s2 = 1
g1 = 3
g2 = 1
c = 1
t = np.linspace(-3, 8, 300)
R = (g1 / 2) * (1 + erf(s1 * c * t)) - (g2 / 2) * (1 + erf(s2 * c * t))

# Calculate moving bar response numerically
steady_f = np.array([quad(f, -np.inf, c * t_) for t_ in t])[:, 0]
transient_f = np.array([quad(lambda x: transient(x, t_), -np.inf, c * t_) for t_ in t])[:, 0]

plt.figure(1)
plt.plot(t, R, linewidth=2)
plt.plot(t, steady_f + transient_f, 'r--', linewidth=2)
plt.xlabel('t')
plt.ylabel('r(t)')


# Now calculate the response to a moving bar, width w
w = 5
c = 1
R1 = (g1/2)*erf(s1*c*t) - (g2/2)*erf(s2*c*t) - (g1/2)*erf(s1*(c*t - w)) + (g2/2)*erf(s2*(c*t - w))
c = 10
R2 = (g1/2)*erf(s1*c*t) - (g2/2)*erf(s2*c*t) - (g1/2)*erf(s1*(c*t - w)) + (g2/2)*erf(s2*(c*t - w))

plt.figure(2)
plt.plot(t, R1, label='c = 1', linewidth=2)
plt.plot(t, R2, label='c = 10', linewidth=2)
plt.xlabel('t')
plt.ylabel('r(t)')
plt.legend()
plt.show()


# Response to a step and a bar when the response is not instantaneous
c = 1
w = 8
n = 200
t = np.linspace(-2, 20, n)
plt.figure(3)
plt.plot(t, imp(t), label='impulse response', linewidth=2)
plt.plot(t, step_response(t), label='step response', linewidth=2)
plt.xlabel('t')
plt.ylabel('r(t)')
plt.legend()
plt.show()


# Response to a moving step when the response is not instantaneous
# You have to be careful with the resolution (the parameter n) here. It's a tradeoff 
# between speed and accuracy. Matlab handles it OK, but Python doesn't even like n=50 much.

w = 8
n = 50
t = np.linspace(-2, 20, n)
clist = [1, 0.5, 3]
for j in range(3):
    c = clist[j]
    bar_response = np.zeros(n)
    for i in range(n):
        integrand = lambda s, x: stim(x, t[i]) * imp(t[i] - s) * f(x)
        bar_response[i] = dblquad(integrand, -20, t[i], -20, 20)[0]   # use 20 instead of np.inf
    plt.plot(t,bar_response,label='c={}'.format(c))
plt.xlabel('t')
plt.ylabel('r(t)')
plt.legend()
plt.show()

# Python code to solve the Longtoin-Milton delay differential model of the pupil
# light reflex.

from pylab import *
from ddeint import ddeint

lam = 30
n = 7
gamma = 5
I = 10
phibar = 1
theta = 10
tau = 1

areafun = lambda x: lam*(theta**n)/(x**n + theta**n)
capF = lambda x: heaviside(x,1)*x
model = lambda x,t : gamma*capF(log(I*areafun(x(t-tau)))/phibar)-x(t)
tt = linspace(0,10,1000) # Time start, time end, num of points/steps
g= lambda t: 10 # solution before the integration interval

yy = ddeint(model,g,tt) # Solving

figure(dpi=600)   # nothing worse than a nasty-looking figure
plot(tt,areafun(yy),c='r')
ylabel('pupil area')
xlabel('time')

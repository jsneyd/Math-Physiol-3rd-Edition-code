% some CO calculations
clearvars
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0);


Thb = 2;

Ko2 = 0.042;

sigmaO = 1.4e-3;
sigmaCO = 3.3e-2;
beta = 200;


Po2 = [40:.01:104];

O = sigmaO*Po2;
U1 = 1.25e-3;
U2 = sigmaCO*7;

Uvals = [0,U1,U2];

for j = 1:length(Uvals)
    U = Uvals(j);
F = O + 4*Thb*O.^4./(Ko2^4+O.^4+beta^4*U.^4);
 

figure(1)
plot(Po2,F);
hold on

end


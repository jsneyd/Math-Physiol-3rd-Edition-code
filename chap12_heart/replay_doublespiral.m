% replay the double spiral simulation

%
%   For Chapter 12, Sections 12.5.3 and 12.5.6 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function main

clear all
close all
clc

set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 2.0, ...
'defaultlinelinewidth', 2.0)

load('doublespiralsimulation.mat')
   
Nsq = N^2;
L=50;
dx = L/N;
X = dx*(1:N)';
Y = X; % a square grid
[X,Y]=meshgrid(dx*(1:N),dx*(1:N));


for j = 1:length(T)
    V = reshape(S(j,1:Nsq),N,N);
    W = reshape(S(j,Nsq+1:2*Nsq),N,N);

    figure(1)
        %mesh(X,Y,V)
        pic = pcolor(X,Y,V);
        shading interp;
        xlabel('x')
        ylabel('y')
        zlabel('v')
        axis([0 L 0 L -.25 1])
        formatSpecF = '%6.2f\n';
        title(strcat('t=',sprintf(formatSpecF,T(j))),'fontsize',18)

          pause(0.1)
    T(j)
end
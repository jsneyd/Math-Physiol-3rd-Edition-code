% 
% For Keener and Sneyd, Mathematical Physiology, Third Edition, Chapter 10

%  -------------------------------------------------------------------
%
%   The repressilator (Elowitz and Leibler (2000) Nature, 403, 335).
%
%   For Chapter 10, Section 10.1.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd.
% 
%  ------------------------------------------------------------------- 
function repressilator 

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
global  alpha n 

% parameters

alpha = 50;
S = 0;
n = 3;
m1=0.7 ;
m2=m1;
m3=0.8;

init=[m1,m2,m3];
total=50;
tstep = 0.1;

tspan = [0:tstep:total];
[T,S] = ode23(@deRHS,tspan, init, odeset('maxstep',10));  

figure(1)
plot(T,S(:,1),T,S(:,2),T,S(:,3))
legend('boxoff')
legend('m_1','m_2','m_3','location','northwest')
box off

%%
function s_prime=deRHS(t,s) 
global  alpha n 

m1=s(1); 
m2=s(2);
m3=s(3);

m1p =  alpha/(1+m2^n) - m1;
m2p =  alpha/(1+m3^n) - m2;
m3p =  alpha/(1+m1^n) - m3;

s_prime =[m1p;m2p;m3p];
 
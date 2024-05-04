%-------------------------------------------------------------------

% Matlab code to compute the flux through a stochastic model of a glucose transporter, 
% using the Gillespie algorithm.

% For Chapter 2,  of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2);
global k kp km  

k = 10; 
km = 1;
kp = 1;  
ci = 0.1; ce = 10; % these concentrations are held fixed

 %states are 
 % S=1, Ce
 % S=2, Pe
 % S=3, Pi
 % S=4, Ci

%the notation here is 
%G_{i,j} = reaction rate for going from state j to state i  

G = [0    km 0 k;
    ce*kp 0  k 0;
    0     k  0 ci*kp;
    k     0 km 0];
CG = cumsum(G)
GS = CG(4,:);
 
 % First do a long simulation using the Gillespie method
% Store in S, with times in T
num = 5000; % number of reaction steps to take
S = zeros(1,num); S(1) = 1;  % S can be 1, 2, 3, or 4
T = zeros(1,num);
count= zeros(4,4); % this is where we count the number of reactions  
trans = zeros(1,num); 
% During the simulation count the number of glucose molecules transferred in each direction.
for i = 1:num-1
        r1 = rand;
        r2 = rand;  % Note that you need two random numbers. One for the time, and one for the choice of state.
        K = GS(S(i));
        T(i+1) = T(i) - (1/K)*log(r1);
        S(i+1) = find(CG(:,S(i))>r2*GS(S(i)),1);  %This determines the next state
   count(S(i+1),S(i)) =  count(S(i+1),S(i)) +1;  % this counts which reaction occurs
  trans(i+1) = (count(2,1)-count(1,2)) ;
end
 
slope = trans(end)/T(end)
  
% now calculate the flux from the deterministic model

A = G-diag(GS);
 
A(4,:) = ones(1,4);
RHS = zeros(4,1);
RHS(4) = 1;
P = A\RHS;
flux =  (G(2,1)*P(1)-G(1,2)*P(2))
 
plot(T,trans,T, flux*T,'--','LineWidth',2)
set(gca,'FontSize',14)
xlabel('time')
ylabel('net number of glucose molecules transferred')
box off












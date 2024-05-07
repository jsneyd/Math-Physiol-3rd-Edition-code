%-------------------------------------------------------------------

% Matlab code for a
% Gillespie simulation of the MM process
% reaction 1  S + E -> C at rate k1
% reaction 2  C -> S + E at rate km1
% reaction 3  C -> S + E at rate k2

% For Chapter 2, Section 2.9.3 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

clear all

set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 2.0, ...
'defaultlinelinewidth', 2.0);


N = 100; %initial substrate number
E = 2; %enzyme numbers

K = 2000;  % number of trials
k1 = 1;
km1 = 1; % without loss of generality
k2  = 0.5;    % without loss of generality

%set the rate constants for two reactions:
r(1) = k1;  % reaction 1  S + E -> C at rate k1 
r(2) = km1; % reaction 2  C -> S + E at rate km1
r(3) = k2;  % reaction 3  C -> P + E at rate k2

%Specify the change matrix
%Ch is the change matrix, a 3 by 2 matrix; three reactions, two state variables
Ch = [-1,1;1,-1;0,-1];  %What happens to s, c  when a reaction occurs:
%initialize the state space

s = N* ones(K,1); % start with N substrate molecules
c = zeros(K,1); %start zero complex molecules

S = s ; %This   keeps track of the trajectories
C = c;
T = zeros(K,1);  %This   tracks the transition times
j = 1; % number of reaction steps
rk = [];

ndx = find(s+c>0);  %update only the s's that are not zero
Nt = length(ndx); %=K at first
Kt = 5000; % max number of time steps to prevent very long runs without termination
while (Nt>0&j<Kt)
j = j+1;
e = E-c;
h(:,1) = r(1)*s.*e;  %reaction rate 1
h(:,2) = r(2)*c; % reaction rate 2
h(:,3) = r(3)*c;  % reaction rate 3 
hc = cumsum(h')'; % the cumulative sum of h
H = sum(h')';

rn = rand(Nt,2); %find a random number for the number of not yet terminated trajectories

T(:,j) = T(:,j-1);  % add the current time to T

for k = 1:Nt
 % update only those for which s+c >0
 
T(ndx(k),j) = - log(rn(k,1))/H(ndx(k)) +T(ndx(k),j-1); % time of next reaction
rk  =  find(rn(k,2) <=hc(ndx(k),:)/H(ndx(k)) ,1); % this determines which reaction occurs
s(ndx(k)) = s(ndx(k)) + Ch(rk,1); % update s 
c(ndx(k)) = c(ndx(k)) + Ch(rk,2); % update c

end 
% save the values of the  trajectories
S(:,j) = s ;
C(:,j) = c ;
ndx = find(s+c>0);  %check to see which trajectories are not extinct yet
if (isempty(ndx)==1)
    Nt = 0;
else
    Nt = length(ndx); %Nt is the number of trajectories not yet extinct

end
end

k = 1; % use the first several trials to plot   sample trajectories

% phase portrait of some trajectories
figure(1)
plot(T(1,:),S(1,:) ,T(1,:),C(1,:))
xlabel('T', 'fontsize', 20)
legend('S','C','fontsize',18)
title('Sample Trajectory','fontsize',20)

figure(2)
plot(S(1,:) , C(1,:))
xlabel('S', 'fontsize', 20)
ylabel( 'C','fontsize',18)
title('Sample Trajectory','fontsize',20)

% now process the data 
% create the pdf for extinction times from the data
[NN,TT]=hist(T(:,end),50);  % histogram of the extinction times with 50 boxes
tt = sort(TT);
dt = mean(TT(2:end)-TT(1:end-1));  %timestep increment on the histogram

% NN/(K*dt) is the approximate pdf for the extinction times from data
% it is normalized to have total integral = 1

figure(3)
NNne0=find(NN>0);
plot(TT(NNne0), NN(NNne0)/(K*dt),'*')
xlabel('Completion time','fontsize',20)
ylabel('Pdf of completion times','fontsize',20)

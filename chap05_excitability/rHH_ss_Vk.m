%rHH steady state as a function of Vk
function rHH_ss
global gnabar gkbar gl Vna Vk0 Vl  rn0 V

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2, ...
   'defaultlinelinewidth', 2.0);
   
%parameters

gnabar = 120.;
gkbar = 36;
gl = 0.3;
Vna = 115;
Vk0 = -12;
Vl = 10.5988;

rn0 = 0.8;
%rn0 = 0.9;

vlist = [0:1.1001:40];
for j = 1:length(vlist)
    V = vlist(j);
    Vk(j) = deRHS(V);
end
figure(1)

plot(Vk,vlist)
xlabel('V_K')
ylabel('V')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the right hand side for ode simulation:
function Vk=deRHS(V)
global gnabar gkbar gl Vna Vk0 Vl  rn0

% given V find Istim
 
Gt= gate_de(V);
AM = Gt(1);
BM = Gt(2);
AN = Gt(3);
BN = Gt(4);
aminf = AM/(AM+BM);

n = AN/(AN+BN);
h = rn0 - n;
Ina = gnabar*aminf.^3.*h.*(V-Vna) ;
 
Icl = gl*(V-Vl);

Vk = (Ina +Icl + gkbar*n.^4*V)./(gkbar*n.^4);
      
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gate = gate_de(V)
 % calculate the gating functions for V                 
AM=.1*(25.-V)./(exp(2.5 -0.1*V)-1.);
BM=4.*exp(-V/18.);
AN=.01*(10.-V)./(exp(1.-0.1*V)-1.);
BN=.125*exp(-V/80.);

gate = [AM BM AN BN];
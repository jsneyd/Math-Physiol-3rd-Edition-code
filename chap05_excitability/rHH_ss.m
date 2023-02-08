%rHH steady state
function rHH_ss
global gnabar gkbar gl Vna Vk Vl  rn0 V

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2, ...
   'defaultlinelinewidth', 2.0);
   
%parameters

gnabar = 120.;
gkbar = 36;
gl = 0.3;
Vna = 115;
Vk = -12;
Vl = 10.5988;

rn0 = 0.8;
%rn0 = 0.9;

vlist = [Vk+.01:1.1001:30];
for j = 1:length(vlist)
    V = vlist(j);
    Ist(j) = deRHS(V);
end
figure(1)

plot(Ist,vlist)
xlabel('I_{app}')
ylabel('V')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the right hand side for ode simulation:
function Istim=deRHS(V)
global gnabar gkbar gl Vna Vk Vl  rn0

% given V find Istim
 
Gt= gate_de(V);
AM=Gt(1);
BM = Gt(2);
AN = Gt(3);
BN = Gt(4);
aminf = AM/(AM+BM);
aminf3 = aminf.^3;
 
n = AN/(AN+BN);
h = rn0 - n;
Ina = gnabar*aminf3.*h.*(V-Vna) ;
Ik = gkbar*n.^4*(V-Vk);
Icl = gl*(V-Vl);
      
Fv = -Ina-Ik-Icl;
Istim = -Fv;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gate = gate_de(V)
 % calculate the gating functions for V                 
AM=.1*(25.-V)./(exp(2.5 -0.1*V)-1.);
BM=4.*exp(-V/18.);
AN=.01*(10.-V)./(exp(1.-0.1*V)-1.);
BN=.125*exp(-V/80.);

AH=.07*exp(-V/20.);
BH=1./(exp(.1*(30.-V))+1.);

gate = [AM BM AN BN AH BH];
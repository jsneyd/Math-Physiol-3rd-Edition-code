% code to solve the microdomain flux equations

%   For Chapter 7, Section 7.6 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd.
% 
%  ------------------------------------------------------------------- 

clear all
close all
clc
set(0,                           ...
    'defaultaxesfontsize', 20,   ...
    'defaultaxeslinewidth', 2.0, ...
    'defaultlinelinewidth', 2.0);
global Kb Dbbt Dc De cinf einf  Jc rho 
% parameters
Dc = 1;
De = 1;
Dbbt = 5;
 
 rho = 2;
cinf = 0.5;
 
 elist = [1,5,10];
 
Kblist = [0.01:0.1:20];
 for je = 1:length(elist)
     einf=elist(je);

for k = 1:length(Kblist)
    Kb = Kblist(k);

    c =  bisect(@getc,cinf,einf);
    J(k)=fluxw(c,Dc,cinf);
 
end
figure(1)
plot(Kblist,J/(rho*(einf-cinf)))
xlabel('K_b')
ylabel('Permeability/\rho')
hold on

 end
text(2,0.336,'e_\infty=1','fontsize',18)
text(2,0.362,'e_\infty=5','fontsize',18)
text(8,0.362,'e_\infty=10','fontsize',18)
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = getc(c)
global Dc cinf einf rho Jc
% first step
% for a given value of c, find e
Jc=fluxw(c,Dc,cinf);
e = bisect2(@gete,c ,einf);

out = rho*(e-c)-Jc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function out = gete(e)
     global   De einf Jc
Je = fluxw(e,De,einf);
out=Je+Jc;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wout,wpout] = w(c,D)
global  Dbbt Kb 

wout = D*c+Dbbt*c./(Kb+c);
wpout = D+Dbbt*Kb./(Kb+c).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = fluxw(c,D,cinf)
 [winf,wpinf]=w(cinf,D);

[w0,wp0] = w(c,D);
out = D*(w0-winf )./wp0;
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%
% this is the bisection algorithm
function root = bisect(feval,a,b)

ul = a;
fl = feval(ul);
uu = b;
fu = feval(uu);

% we make the assumption, without checking, that 
% fu*fl<0

% if not, the algoritm fails to find a root.

N = 25;  % number of iterates
% the main bisection algorithm
for j = 1:N
u = (ul+uu)/2;
fc = feval(u);
ftest = (fc*fl>=0);
ul = ftest*u+(1-ftest)*ul;
fl = ftest*fc + (1-ftest)*fl;

uu = (1-ftest)*u+ ftest*uu;
fu = (1-ftest)*fc + ftest*fu;
end
root = u;
end

% another copy of the bisection algorithm
% this is the bisection algorithm
function root = bisect2(feval,a,b)

ul = a;
fl = feval(ul);
uu = b;
fu = feval(uu);

% we make the assumption, without checking, that 
% fu*fl<0

% if not, the algoritm fails to find a root.

N = 25;  % number of iterates
% the main bisection algorithm
for j = 1:N
u = (ul+uu)/2;
fc = feval(u);
ftest = (fc*fl>=0);
ul = ftest*u+(1-ftest)*ul;
fl = ftest*fc + (1-ftest)*fl;

uu = (1-ftest)*u+ ftest*uu;
fu = (1-ftest)*fc + ftest*fu;
end
root = u;
end



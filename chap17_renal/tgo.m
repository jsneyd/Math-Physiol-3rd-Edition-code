
%    -------------------------------------------------------------------
%
%     Model of tubuloglomerular oscillations.
%
%     For Chapter 17, Section 17.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function tgo
clear all
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
global numx delx delt numtimesteps delay
global c K1 K2 Cop k

numx=100; delx=1/numx; delt=0.001;
numtimesteps=20000; delay=200;
keepend=zeros(numtimesteps+delay,1);
time = [1:numtimesteps+delay]*delt;
c=zeros(numx,1);
c(1)=1;
k=1; Cop=exp(-k);
K1=1.5;


% Now solve the PDE by hand (i.e., without using any packaged routines).
% Take differences in space and time, advance the solution in time, and
% keep track of the solution at the end of the domain, as it is used for
% the delay condition inside phi, as well as for the output.

% The delay is handled in the vector keepend, which we fill up forwards in
% time (i.e., at timestep i + delay), while using the value at timestep i.

gammalist = [3,2];
for jj = 1:length(gammalist)
    gamma=gammalist(jj);;
    K2=gamma/(K1*k*exp(-k));
    for i=1:numtimesteps
        keepend(i+delay)=c(numx);
        c(2:numx) = c(2:numx) - ...
                    delt*k*c(2:numx) - phi(keepend(i))*(delt/delx)*(c(2:numx)-c(1:numx-1));
    end
    figure(jj)
        plot(time,keepend)
        xlabel('t')
        ylabel('macula densa [Cl^-], c')
        save('tgoosc.dat','-ascii','keepend')
        axis([0 20 0 0.7])
        formatSpecF = '%6.2f\n';
        title(strcat('\gamma = ',sprintf(formatSpecF,gamma)))
end

% stability curves
tbar = [0:.001:.3];
for n = 1:4
    w = n*pi./(tbar+1/2);
    g = (-1)^(n+1)*w./(2*sin(w/2));
    ndx = find(g>0);
    figure(3)
        plot(tbar(ndx),g(ndx))
        box off
        hold on
end
axis([0 0.3 0 20])
xlabel('\tau')
ylabel('\gamma')
text(0.1,19,'n=4','fontsize',18)
text(0.17,15,'n=3','fontsize',18)
text(0.2,5.5,'n=2','fontsize',18)
text(.25,3,'n=1','fontsize',18)
text(0.05,17,'unstable','fontsize',18)
text(0.03,3,'stable','fontsize',18)
hold off

end % of main

%%
function out=phi(x)
    global c K1 K2 Cop k
    out = 1 + K1*tanh(K2*(Cop-x));
end



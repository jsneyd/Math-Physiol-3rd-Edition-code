% shooting to find spiral wave profile
 
%   For Chapter 6, Figure 1.14,   of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  ------------------------------------------------------------------- 

 clear all; close all; clc;
 set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);  
 
global c w  
%set parameters 
eps = 0.05;
alp = 0.1;
c=3.0;
mstar = 0.3318;
wlist=[2.8, 9*mstar,3.2];
for j = 1:3

 w=wlist(j);
rstep = -0.001;
r_end = 10;
s0 = [r_end-0.1,0];

tspan = [r_end:rstep:0.0 ];
[r,S] = ode23s(@deRHS,tspan,s0); 
 
ph = S(:,1);
th=S(:,2);

    ps = ph./(sqrt(r.^2-ph.^2)); 
    figure(1)
    ndx = find(ph.^2<r.^2);

    plot(r(ndx),ps(ndx))
    hold on

    if(j==2)
        r0=c/w;
        R =[r0:.01:r_end];
        rho = sqrt(R.^2/r0^2-1);
        th0=2*[0:.01:1]*pi;
        X0=r0.*cos(th0);
        Y0 = r0.*sin(th0);
s=rho-rho(end)+th(1)+pi/2;
X1=r0*cos(s)+r0*rho.*sin(s);
Y1=r0*sin(s)-r0*rho.*cos(s);

        figure(2)
        X=r.*cos(th);
        Y = r.*sin(th);
        plot(X,Y,X0,Y0,'--',X1,Y1,'--')
axis([-10 10 -10 10])

    end
end
 figure(1)
xlabel('r')
ylabel('\psi')
axis([0 5 -4.5 4.5])
legend('boxoff')
legend('\omega/\epsilon=2.8','\omega/\epsilon=2.9864', '\omega/\epsilon=3.2','location','northwest')
box off
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The singular dispersion curve
   w0=[0:.001:(1-2*alp)/2];
   apw = alp+w0;
   w1=1-2*alp-w0;
   sp=(1-2*apw)./sqrt(apw-apw.^2);
   T = log(((1-w0).*w1)./((1-w1).*w0));

   ccrit = sqrt( eps./(T*mstar));

   figure(10)
    plot(2*pi./T,sp,'--',1./T,ccrit)
    box off
    axis([0 10 0 3])
    xlabel('Frequency, 2\pi/T')
    ylabel('Speed, c')
    box off
text(1.8,2.3,'dispersion curve','fontsize',18)
text(7,.9,'critical curve','fontsize',18)
     
   hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function out = deRHS(r,s)
    global c w
ph = s(1);
th=s(2);

     ps = ph./(sqrt(r.^2-ph.^2)); 
    dph = r*(c-w*sqrt(r^2-ph^2));
    dth = ps/r;
out =[dph;dth];
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% event function
function [value,isterminal,direction]=events(t,s)
 
value = (abs(s)-100); 
isterminal = 1;
direction = 0;
end
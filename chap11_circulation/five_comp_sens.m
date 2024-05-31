
%  -------------------------------------------------------------------
%
%   Calculate sensitivies in the five-compartment circulation model
%
%   For Chapter 11, Section 11.5.4 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function five_comp_sens

clear all
close all
clc

% first identify the variables, set nominal values
V(1)=100; %Psa
V(2) = 30; %Ps;
V(3)=2;  % Psv
V(4)=15; %Ppa
V(5) = 5; %Ppv;
V(6) = 1; %Vsa
V(7) = 3.5; %Vsv
V(8) = 0.5; %Vp
V(9) = 5.6% Q
sens=zeros(9,9);

%get the nominal parameters
P=getparams(V);
P(9)=80;
%this checks the calculation
Vn=getpressures(P);


% now calculate sensitivites;
for j=1:9
    % specfiy a perturbation for each parameter
    frac=0.001;
    del=P(j)*frac;
    Pdel=P;
    Pdel(j) = P(j)+del;
    Vdel=getpressures(Pdel);
    sens(:,j)=(Vdel-V)./(frac*V);
end

% the sensitivity matrix is:
sens

end % of main


%%
function params = getparams(V)
% for the input variables find the corresponding parameters
V0s= 0.94;
V0p = 0.26;
F=80;

Psa=V(1);
Ps = V(2);
Psv = V(3);
Ppa= V(4);
Ppv = V(5);
Vsa= V(6);
Vsv=V(7);
Vp=V(8);
Q=V(9);

% now the equations
Rsa = (Psa-Ps)/Q;
Cld = Q/(F*Ppv);;
Csa = 2*(Vsa-V0s)/(Psa+Ps);
Rsv = (Ps-Psv)/Q;
Csv = 2*Vsv/(Psv+Ps);
Rp=(Ppa-Ppv)/Q;
Crd = Q/(Psv*F);
Cp=2*(Vp-V0p)/(Ppa+Ppv);
params = [Csa,Csv,Cp,Rsa,Rsv,Rp,Cld,Crd];

end

%%
function out = getpressures(P)
% input parameters and get pressures and volumes
Csa=P(1);
Csv=P(2);
Cp=P(3);
Rsa=P(4);
Rsv=P(5);
Rp=P(6);
Cld=P(7);
Crd=P(8);
F=P(9);  %nominal rate;

% given these parameters, find the pressures and volumes
V0s= 0.94;
V0p = 0.26;
Vt = 5;
Ve=Vt-V0s-V0p;

alp=Rsv*(Csa+Csv/2)+Rsa*Csa/2+Cp*Rp/2;

Q=Ve/(alp+Cp/(F*Cld)+(Csv+Csa)/(F*Crd));
Psa=Q*(1/(F*Crd)+Rsa+Rsv);
Ps = Q*(1/(F*Crd)+Rsv);
Psv=Q/(F*Crd);
Ppa = Q*(1/(F*Cld)+Rp);
Ppv=Q/(F*Cld);
Vsa=V0s+Csa*(Psa+Ps)/2;
Vsv=Csv*(Psv+Ps)/2;
Vp=V0p+Cp*(Ppa+Ppv)/2;

out=[Psa,Ps,Psv,Ppa,Ppv,Vsa,Vsv,Vp,Q];

end


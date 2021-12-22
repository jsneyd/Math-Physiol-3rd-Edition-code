
function main

clear all;
global diff r c0 P startbigL N0 npts del bigL

diff=1000.0; r=0.05; c0=0.3; P=0.2; startbigL=100.0; N0=0.3;
npts=10000;

bigL=startbigL;
del=bigL/npts;

cinitup=15.0;
cinitdown=0.01; 
[dist c out]=doshoot(cinitup);
plot(dist,c) 
axis([0 100 0 1.3])
hold on
[dist c out]=doshoot(cinitdown);
plot(dist,c)

for i=1:30
    cinit=(cinitup+cinitdown)/2.0d0;
	[dist c out]= doshoot(cinit);
    plot(dist,c)
	if (out<0) 
        cinitdown=(cinitup+cinitdown)/2.0d0;
    end
	if (out>0) 
        cinitup=(cinitup+cinitdown)/2.0d0;
    end
end
    
hold off

%--------------------------------------

function [dist c out]=doshoot(cinit)
global diff r c0 P startbigL N0 npts del bigL

d(1)=0.0; v(1)=0.0; c(1)=cinit;
for i=1:npts-1
    dist(i+1)=(i-1)*del;
    c(i+1) = c(i) + del*d(i);
    d(i+1) = d(i) + del*(1.0/diff)*(v(i)*d(i) + c(i)*2.0*(P/r)*(c(i)-c0) - 2.0*bigN(dist)/r);
    v(i+1) = v(i) + del*2.0*(P/r)*(c(i)-c0);
    if (c(i+1)<0.0)
        out=-100.0;
        return
    end
    if (c(i+1)>100.0)
        out=100.0;
        return
    end
end

out=c(npts)-c0;

%--------------------------------------

function out=bigN(dum)
global diff r c0 P startbigL N0 npts del bigL

out=0.0;
if (dum<bigL/10.0)
    out=N0;
end

% diffusion on a lattice (Chen and Meng)
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
    
%parameters
% numer of lattice points
N = 50;
%number of particles
M = 1000;

%I, J, and K represent the position of each particle
randomly distribute the particles on the cell
I=randi(1:N,M,1);
J=randi(1:N,M,1);
K=randi(1:N,M,1);


plot3d(I,J,K,'*')
% diffusion on a lattice (Chen and Meng)
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
    
%parameters
% numer of lattice points
N = 50;
%number of particles
M = 10;

%I, J, and K represent the position of each particle
%randomly distribute the particles on the cell
pos= [randi(N,M,1),randi(N,M,1),randi(N,M,1)];
 


plot3(pos(:,1),pos(:,2),pos(:,3),'*')
pos
% now allow points to diffuse
r1=randi(3,M,1)
r2=2*randi(2,M,1)-3
for j = 1:M
    pos(j,r1(j)) =  pos(j,r1(j))+r2(j);
end

pos
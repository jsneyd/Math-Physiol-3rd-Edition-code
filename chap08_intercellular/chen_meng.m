% diffusion on a lattice (Chen and Meng)
clear all 
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
    
%parameters
% number of lattice points
N = 50;
%number of particles
M = 1000;
% number of gap junctions
K = 100; %(must be a perfect square)
% number of steps to take
ncase=[2,1];
for nnr = 1:2
gcase=ncase(nnr)

nruns = 5;
for nr = 1:nruns
 
%randomly distribute the particles in the N x N x N cell
 
pos= [randi(N-1,M,1), randi(N-1,M,1),randi(N-1,M,1)];
 
%locate the gap junctions
clear gpos
 if(gcase==1)
     gpos=randi(N,K,2); % random distribution of K gap junctions
 else
        sqK = sqrt(K);  % square array of gap junctions
        strt=fix((N-sqK)/2);

        for j = 1:sqK
            for k = 1:sqK
              ndx = (j-1)*sqK+k;
              gpos(ndx,1)=  strt+j;
              gpos(ndx,2) = strt+k;
            end
        end
    end
 
 
 jstep = 0;
 nk = 0;
 nkthresh = 200;
 clear nstr

while nk<nkthresh % nk is the number of particles in box 2

%  allow points to diffuse
r1=randi(3,M,1); % this is the index that will change
r2=2*randi(2,M,1)-3; % this is the amount of change of the index
 posnew=pos;
for j = 1:M
    posnew(j,r1(j)) =  pos(j,r1(j))+r2(j);
 
 %reflect at the boundaries
 
     if(posnew(j,1)==0) posnew(j,1)=1; end;
     if(posnew(j,2)==0) posnew(j,2)=1; end;
     if(posnew(j,3)==0) posnew(j,3)=1; end;
     if(posnew(j,1)>N) posnew(j,1)=N; end; 
     if(posnew(j,2)>N) posnew(j,2)=N; end;
     if(posnew(j,3)>2*N) posnew(j,3)=2*N;end;
             
      mdx = 0;
            
      if(posnew(j,3)==N)% this  point is at the interface
            
           for jj = 1:K  % check all the gap junctions
                 if((posnew(j,1)==gpos(jj,1))&(posnew(j,2)==gpos(jj,2)))
                     mdx=mdx+1;
              
                 end
            end
                  if (mdx==0) %this is not a gap junction point; reflect it
                 posnew(j,3)=pos(j,3); % set back to previous value
                 else % let it diffuse  
                     nindx = N+2*randi(2,1,1)-3; % either N+1 or N-1
                  posnew(j,3)=  nindx;
                 end
       end
 end
 
                 pos=posnew; %update the points
                 jstep=jstep+1;
%count the number of particles in box 2
 ldx=0;  
for j = 1:M
if(pos(j,3)>N)
    ldx=ldx+1;
end
end
 
if(ldx>nk)
    nstr(nk+1:ldx)=jstep;
end
nk=ldx; % number of particles in box 2
 end
figure(10)
subplot(1,2,gcase)
plot(100*[1:nkthresh]/M,nstr)
hold on
  
              % all done; plot the final point distribution
end

             kdx=0;
             ldx = 0;
             clear lldx
             clear kkdx

           for j=1:M
               if(pos(j,3)>N)
                   ldx=ldx+1;
                   lldx(ldx)=j;
               else
                   kdx=kdx+1;
                  kkdx(kdx)=j;
               end
           end

 
figure(nnr)
           plot3(pos(kkdx,1),pos(kkdx,2),pos(kkdx,3),'b*')
          hold on
           plot3(pos(lldx,1),pos(lldx,2),pos(lldx,3),'r*')
           
           plot3(gpos(:,1),gpos(:,2),N*ones(K,1),'g*')
hold off
end
 figure(10)
 xlabel('Percentage of molecules in target cell')
 ylabel('Steps')
 hold off

 
%  -------------------------------------------------------------------
%
%   Code for diffusion on a lattice and through gap junctions.
%
%   For Chapter 8, Section 8.3.4 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd
%
%  -------------------------------------------------------------------

%warning('off', 'Octave:possible-matlab-short-circuit-operator');

clear all
close all
clc

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
K = 100; % must be a perfect square
% number of steps to take
ncase=[2,1];

for nnr = 1:2
    gcase=ncase(nnr);
    nruns = 10;

    for nr = 1:nruns

        % randomly distribute the particles in the N x N x N cell
        pos= [randi(N-1,M,1), randi(N-1,M,1),randi(N-1,M,1)];

        %locate the gap junctions
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
        nkthresh = 200;  % threshold is 20% of the original number of particles

        while nk<nkthresh % nk is the number of particles in box 2. keep repeating until 20% have transferred to box 2
            %  allow points to diffuse
            r1=randi(3,M,1); % this is the index that will change
            r2=2*randi(2,M,1)-3; % this is the amount of change of the index
            posnew=pos;

            for j = 1:M
                posnew(j,r1(j)) =  pos(j,r1(j))+r2(j);

                %reflect at the boundaries
                if(posnew(j,1)==0) posnew(j,1)=1; end
                if(posnew(j,2)==0) posnew(j,2)=1; end
                if(posnew(j,3)==0) posnew(j,3)=1; end
                if(posnew(j,1)>N) posnew(j,1)=N; end
                if(posnew(j,2)>N) posnew(j,2)=N; end
                if(posnew(j,3)>2*N) posnew(j,3)=2*N;end

                mdx = 0;
                if(posnew(j,3)==N)% this  point is at the interface
                    for jj = 1:K  % check all the gap junctions
                         if((posnew(j,1)==gpos(jj,1))&(posnew(j,2)==gpos(jj,2)))
                            mdx=1;
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

            pos=posnew; % update the points
            jstep=jstep+1;
            %count the number of particles in box 2
            ldx = sum(pos(:,3)>N);

            if(ldx>nk)
                nstr(nk+1:ldx)=jstep;  % This is the step that "fills up" nstr
            end
            nk=ldx; % number of particles in box 2
        end

        figure(10)
        if gcase==1 plot(100*[1:nkthresh]/M,nstr,'r')
        end
        if gcase==2 plot(100*[1:nkthresh]/M,nstr,'b')
        end
        box off
        hold on

         % for external plotting
%          if nnr==1
%              keep1(:,nr) = nstr;
%          else
%              keep2(:,nr) = nstr;
%          end

    end
   
    jstep

    % all done; plot the final point distribution
    up = find(pos(:,3)>N);
    down = find(pos(:,3)<=N);

    figure(nnr)
    plot3(pos(up,1),pos(up,2),pos(up,3),'r*')
    hold on
    plot3(pos(down,1),pos(down,2),pos(down,3),'b*')
    plot3(gpos(:,1),gpos(:,2),N*ones(K,1),'g*')
    hold off
end

% for external plotting
%save('chen_meng_1.dat','keep1','-ascii')
%save('chen_meng_2.dat','keep2','-ascii')

figure(10)
xlabel('Percentage of molecules in target cell')
ylabel('Steps')
hold off


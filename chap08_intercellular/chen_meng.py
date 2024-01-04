#   -------------------------------------------------------------------
# 
#    Code for diffusion on a lattice and through gap junctions.
# 
#    For Chapter 8, Section 8.3.4 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

# Parameters
N = 50  # Number of lattice points
M = 1000  # Number of particles
K = 100  # Number of gap junctions (must be a perfect square)
ncase = [2, 1]

for nnr in range(2):
    gcase = ncase[nnr]
    nruns = 1

    for nr in range(nruns):
        # Randomly distribute the particles in the N x N x N cell
        pos = np.random.randint(1, N, (M, 3))

        # Locate the gap junctions
        if gcase == 1:
            gpos = np.random.randint(1, N+1, (K, 2))  # random distribution
        else:
            sqK = int(np.sqrt(K))
            strt = int((N - sqK) / 2)
            # square array
            gpos = np.array([(strt + j, strt + k) for j in range(sqK) for k in range(sqK)])

        jstep = 0
        nk = 0
        nkthresh = 200 # threshold is 20% of the original number of particles
        nstr = np.zeros(nkthresh)

        # nk is the number of particles in box 2. keep repeating until 20% have transferred to box 2.
        while nk < nkthresh: 
            # allow points to diffuse
            r1 = np.random.randint(0, 3, M)  # Index that will change
            r2 = 2 * np.random.randint(0, 2, M) - 1  # Amount of change of the index
            posnew = pos

            for j in range(M):
                posnew[j, r1[j]] = pos[j, r1[j]] + r2[j]

                # Reflect at the boundaries
                if posnew[j, 0] == 0:
                    posnew[j,0] = 1
                if posnew[j, 1] == 0:
                    posnew[j,1] = 1
                if posnew[j, 2] == 0:
                    posnew[j,2] = 1
                    
                if posnew[j, 0] > N:
                    posnew[j,0] = N
                if posnew[j, 1] > N:
                    posnew[j,1] = N
                if posnew[j, 2] > 2*N:
                    posnew[j,2] = 2*N
                
                mdx = 0
                
                if posnew[j, 2] == N:  # this point is at the interface
                    for jj in range(K): # look through all the gap junctions
                        if (posnew[j,0]==gpos[jj,0]) and (posnew[j,1]==gpos[jj,1]):
                            mdx = 1
                            
                    if mdx==0:                      # not a gap junction point so reflect it
                        posnew[j,2] = pos[j,2]      # set back to previous value
                    else:                           # let it diffuse
                        nindx = N + 2*np.random.randint(1, 3) - 3
                        posnew[j,2] = nindx
            pos = posnew
            jstep += 1

            # Count the number of particles in box 2
            ldx = np.sum(pos[:, 2] > N)

            if ldx > nk:
                nstr[nk:ldx+1] = jstep # This is the step that "fills up" nstr
            nk = ldx

        plt.figure(10)
        if gcase==1:
            plt.plot(100 * np.arange(1, nkthresh + 1) / M, nstr,'r')
        if gcase==2:
            plt.plot(100 * np.arange(1, nkthresh + 1) / M, nstr,'b')
        plt.xlabel('Percentage of molecules in target cell')
        plt.ylabel('Steps')
        plt.title(f'Case {gcase}')
        print(jstep)
    
    
    # All done; plot the final point distribution
    kdx = np.where(pos[:, 2] <= N)[0]
    ldx = np.where(pos[:, 2] > N)[0]

    fig = plt.figure(nnr)
    ax = plt.axes(projection='3d')
    
    ax.scatter3D(pos[kdx, 0], pos[kdx, 1], pos[kdx, 2], c='blue')
    ax.scatter3D(pos[ldx, 0], pos[ldx, 1], pos[ldx, 2], c='red')
    ax.scatter(gpos[:, 0], gpos[:, 1], 50, c='green')

plt.show()
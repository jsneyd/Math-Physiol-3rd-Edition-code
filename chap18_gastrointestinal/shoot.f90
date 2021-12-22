use shootcom
implicit none

diff=1000.0d0; r=0.05d0; c0=0.3; P=0.2d0; startbigL=100.0d0; N0=1.0d-1
npts=10000
!littlel=bigL/r; eps=diff/(r*c0*P); littlen=N0/(c0*c0*P)

open(1,file='shoot.dat')

print*,'enter 1 for single solution'
print*,'enter 2 for efflux conc as function of channel length'
read*,ichoice

if (ichoice .eq. 1) nLstep=1
if (ichoice .ne. 1) nLstep=50

do Lstep=1,nLstep
        if (ichoice .eq. 1) bigL=startbigL
	if (ichoice .ne. 1) bigL= startbigL - (Lstep-1.0d0)*(startbigL-10.0d0)/(float(nLstep-1))
        del=bigL/float(npts)

	cinitup=15.0d0
	cinitdown=0.01d0
	c(1)=cinitdown; call doshoot; if (ichoice .eq. 1) print*, shootout
	c(1)=cinitup; call doshoot; if (ichoice .eq. 1) print*,shootout
	
	shootout=10.0d0
	
	do while (dabs(shootout) > 0.001)
	    c(1)=(cinitup+cinitdown)/2.0d0
	    call doshoot
	    if (shootout<0) cinitdown=(cinitup+cinitdown)/2.0d0
	    if (shootout>0) cinitup=(cinitup+cinitdown)/2.0d0
            if (ichoice == 1) then
	        print*,cinitdown,cinitup,shootout
            endif
	end do
	
        dcdx=(c(npts)-c(npts-1))/del
        outos=c0-diff*dcdx/v(npts)
        if (ichoice .eq. 1) then
		do i=1,npts
		    write(1,*) c(i),v(i)
		end do	
		print*,'solute outflow=',c(npts)*v(npts) - diff*dcdx
		print*,'solute inflow=',2.0d0*(bigL/10.0d0)*N0/r
		print*,'outflow osmolarity =', outos
        endif
	if (ichoice .ne. 1) then
                print*, bigL,outos
	        write(1,*) bigL,outos
	endif
end do

end

!-----------------------------------------------
subroutine doshoot
use shootcom
implicit none
double precision :: N

d(1)=0.0d0; v(1)=0.0d0
do i=1,npts-1
    dist=(i-1)*del
    c(i+1) = c(i) + del*d(i)
    d(i+1) = d(i) + del*(1.0d0/diff)*(v(i)*d(i) + c(i)*2.0d0*(P/r)*(c(i)-c0) - 2.0d0*N(dist)/r)
    v(i+1) = v(i) + del*2.0d0*(P/r)*(c(i)-c0)
    if (c(i+1)<0.0) then
        shootout=-100.0d0
        return
    endif
    if (c(i+1)>100.0d0) then
        shootout=100.0
        return
    endif
end do

shootout=c(npts)-c0

return
end

!-----------------------------
function N(dum)
use shootcom
implicit none
double precision :: N,dum

N=0.0d0
if (dum<bigL/10.0d0) N=N0

return
end
!-------------------------------
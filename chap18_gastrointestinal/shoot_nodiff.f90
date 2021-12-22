
implicit none
double precision :: diff,r,c0,cinit,P,del,bigL,dist,N0,N,c1,c2,c3,solutein,pi
double precision :: littlel,littlen,eps
double precision :: c(10000),d(10000),v(10000)
integer :: i,npts

r=0.05d0; c0=0.3; P=0.2d0; bigL=100.0d0; N0=3.0d-1
npts=10000
littlen=N0/(c0*c0*P)
del=1.0d0/float(npts)

pi=3.14159

open(1,file='shoot_nodiff.dat')

c1=1.0d0; c2=-1.0d0; c3=-N(0.0d0,littlen)
cinit = (-c2 + dsqrt(c2*c2 - 4.0d0*c1*c3))/(2.0d0*c1)
c(1)=cinit; v(1)=0.0d0
c(2)=cinit; v(2)=v(1)+del*2.0d0*(c(1)-1.0d0)
do i=2,npts-1
    dist=(i-1)*del
    c(i+1) = c(i) + del*(1.0d0/v(i))*(-c(i)*2.0d0*(c(i)-1.0d0) + 2.0d0*N(dist,littlen))
    v(i+1) = v(i) + del*2.0d0*(c(i)-1.0d0)
    write(1,*) dist,c(i),v(i)
end do
write(1,*) ' '

print*, 'c(npts)=',c(npts),v(npts)
eps = 0.1d0
do i=1,npts
    dist=i*del
    write(1,*) dist,c(i)+exp(v(npts)*(dist-1.0d0)/eps)*(1.0d0-c(npts))
end do
end

!-----------------------------
function N(dum,littlen)
implicit none
double precision :: N,dum,littlel,littlen

N=0.0d0
if (dum<1.0d0/10.0d0) N=littlen

return
end
!-------------------------------

implicit none
integer :: i
double precision :: addup,tau

open(1,file='plotfdf.dat')
do i=1,100
    tau = 0.1 + (i-1.0d0)*100.0d0/99.0d0
    write(1,*) tau, addup(tau)
end do

end


function addup(tau)
double precision :: tau, addup,pi
integer :: num, i

pi = 3.14159
num = 50000
addup = 0.0
do i=1,num
    addup = addup + dsqrt(1.0d0/(4.0d0*pi*tau*i))*dexp(-i/(4*tau))
end do

end

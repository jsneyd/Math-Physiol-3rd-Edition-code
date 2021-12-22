module shootcom
implicit none
double precision :: diff,r,c0,cinit,P,del,startbigL,bigL,dist,N0,shootout,cinitup,cinitdown,dcdx,outos
double precision :: littlel,eps,littlen
double precision :: c(10000),d(10000),v(10000)
integer :: i,npts,ichoice,Lstep,nLstep
end module shootcom
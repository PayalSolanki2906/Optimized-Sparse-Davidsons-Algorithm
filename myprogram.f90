program davidson_algo
integer,parameter::s=6,k=5,nnz=5050,s2=100
real*8,dimension(0:nnz-1)::h1,h2
integer,dimension(0:nnz-1)::ia,ja
complex*16::HMG(0:nnz-1),evec(0:s2-1,0:k-1)
real*8::tol
real*8::eval(0:k-1),x11
integer::t1,t2,t3,t4,clock_max
real*8::clock_rate
tol=0.000000000005
call system_clock ( t1, clock_rate, clock_max )
do i=0,nnz-1
read(124,*)h1(i),h2(i),ia(i),ja(i)
end do
do i=0,nnz-1
HMG(i)=dcmplx(h1(i),h2(i))
end do
call system_clock ( t2, clock_rate, clock_max )
call system_clock ( t3, clock_rate, clock_max )
call OSDAHC(HMG,nnz,ia,ja,s2,k,tol,eval,evec,20,'l',1000)
call system_clock ( t4, clock_rate, clock_max )
do i=0,k-1
write(201,*)eval(i)
end do
!write(126,*)'time for ham=',  real ( t2 - t1 ) /clock_rate
!write(126,*)'time for diagonalization=',  real ( t4 - t3 ) /clock_rate
end program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

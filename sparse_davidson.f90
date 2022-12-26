
subroutine DAVIDSONSRSPARSE(H,nnz,ia,ja,n,k,tol,eval,evec,b,maxi)
implicit none
integer::k,n,INFO,i,j,i1,m,j1,nnz,m1,b,maxi
real*8,dimension(0:nnz-1)::H
integer,dimension(0:nnz-1)::ia,ja
real*8::tol,olde,t(0:n-1),newe,random_normal,diff,eval(0:k-1)
real*8::temp2(0:n-1),temp3(0:n-1),alpha,beta
real*8::w1(0:n-1),q(0:n-1),evec(0:n-1,0:k-1),const(0:k-1)
real*8, dimension (:,:), allocatable :: A,VP,temp,PPP,V
real*8,dimension (:), allocatable ::WORK
real*8,dimension (:), allocatable ::W,RWORK
character::matdescra(6)
matdescra(1)='S'
matdescra(2)='U'
matdescra(3)='N'
matdescra(4)='F'
alpha=dcmplx(1,0)
beta=dcmplx(0,0)
allocate(V(0:n-1,0:k-1))
do j=0,k-1
call random_number(t)
V(:,j)=t/norm2(t)
end do
newe=tol
do m1=1,maxi,1
do m=k,b,k
allocate(A(0:(m-k)-1,0:m-k-1),W(0:m-k-1),WORK(0:3*(m-k)-1),VP(0:n-1,0:m-k-1))
allocate(temp(0:n-1,0:m-k-1),RWORK(0:3*(m-k)-3),ppp(0:m-k-1,0:m-k-1))
		if(m .gt. k) then
		olde=newe
		VP=V(:,0:m-k-1)
		deallocate(V)
		call QRR(n,m-k,VP)
		temp=0.0d0
                call mkl_dcoomm('n', n, m-k, n, alpha, matdescra, H, ia, ja, nnz, VP, n, beta, temp, n)
		ppp=matmul(transpose(VP),temp)
		A=ppp
		call DSYEV( 'V', 'U', m-k, A,m-k, W, WORK, 3*(m-k)-1, INFO )
		newe=W(k-1)
                write(122,*)W(k-1)
		allocate(V(0:n-1,0:m-1))
		V(:,0:m-k-1)=VP
		do j=0,k-1
		q=0.0d0
		temp3=matmul(VP,A(:,j))
		V(:,j)=temp3
		call mkl_dcoomm('n', n, 1, n, alpha, matdescra, H, ia, ja, nnz, temp3, n, beta, temp2, n)
		w1=temp2-W(j)*temp3
		q=w1/(W(j)-const(j))
 		V(:,m-k+j)=w1
		end do
	end if
diff=olde-newe
if ( abs(diff) .le. tol) then
do i1=0,k-1,1
evec(:,i1)=matmul(VP,A(:,i1))
eval(i1)=W(i1)
end do
end if
deallocate(ppp,temp,RWORK,A,W,WORK,VP)
end do
if ( abs(diff) .le. tol) then
exit
end if
end do
deallocate(V)
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine OSDASR(H,nnz,ia,ja,n,k,tol,eval,evec,b,order,maxi)
implicit none
integer::k,n,INFO,i,j,i1,m,j1,nnz,m1,b,IL,IU,MM,xx,maxi
real*8,dimension(0:nnz-1)::H
integer,dimension(0:nnz-1)::ia,ja
real*8::tol,olde,t(0:n-1),newe,diff,eval(0:k-1)
real*8::temp2(0:n-1),temp3(0:n-1),alpha,beta
real*8::w1(0:n-1),evec(0:n-1,0:k-1)
real*8, dimension (:,:), allocatable :: A,VP,temp,PPP,V,Z
real*8,dimension (:), allocatable ::WORK,IWORK,IFAIL
real*8,dimension (:), allocatable ::W,RWORK
character::matdescra(6)
real*8::VL,VU,tol1
character*1::order
tol1=0.0000000000000001d0
matdescra(1)='S'
matdescra(2)='U'
matdescra(3)='N'
matdescra(4)='F'
alpha=dcmplx(1,0)
beta=dcmplx(0,0)
allocate(V(0:n-1,0:k-1))
do j=0,k-1
call random_number(t)
V(:,j)=t/norm2(t)
end do
call random_number(newe)
do m1=1,maxi,1
do m=k,b,k
allocate(A(0:(m-k)-1,0:m-k-1),W(0:m-k-1),WORK(0:8*(m-k)-1))
allocate(VP(0:n-1,0:m-k-1),ppp(0:m-k-1,0:m-k-1))
allocate(temp(0:n-1,0:m-k-1),RWORK(0:3*(m-k)-3))
allocate(IWORK(0:5*(m-k)-1),IFAIL(0:m-k-1),Z(0:m-k-1,0:k-1))
	if(m .gt. k) then
		olde=newe
		VP=V(:,0:m-k-1)
		deallocate(V)
		call QRR(n,m-k,VP)
		temp=0.0d0
		call mkl_dcoomm('n', n, m-k, n, alpha, &
        matdescra, H, ia, ja, nnz, VP, n, beta, temp, n)
		ppp=matmul(transpose(VP),temp)
		A=ppp
		if (order=='l') then 
			IL=1
			IU=k
			xx=k
		end if
		if (order=='h') then 
			IL=m-k-k+1
			IU=m-k
			xx=1
		end if
		call dsyevx('V','I','U',m-k,A,m-k,VL,VU,IL,IU, &
        tol1,MM,W,Z,m-k,WORK,8*(m-k),IWORK,IFAIL,INFO)
		newe=W(xx-1)
		allocate(V(0:n-1,0:m-1))
		V(:,0:m-k-1)=VP
		do j=0,k-1
			temp3=matmul(VP,Z(:,j))
			V(:,j)=temp3
			call mkl_dcoomm('n', n, 1, n, alpha, &
               matdescra, H, ia, ja, nnz, temp3, &
                       n, beta, temp2, n)
			w1=temp2-W(j)*temp3
			V(:,m-k+j)=w1
		end do
	end if
diff=olde-newe
if (abs(diff) .le. tol) then
do i1=0,k-1,1
evec(:,i1)=matmul(VP,Z(:,i1))
eval(i1)=W(i1)
end do
end if
deallocate(ppp,temp,RWORK,A,W,WORK,VP,Z,IFAIL,IWORK)
end do
if (abs(diff) .le. tol) then
exit
end if
end do
deallocate(V)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine OSDAHC(H,nnz,ia,ja,n,k,tol,eval,evec,b,order,maxi)
implicit none
integer::k,n,INFO,i,j,i1,m,nnz,b,m1,IL,IU,MM,xx,maxi
complex*16,dimension(0:nnz-1)::H
integer,dimension(0:nnz-1)::ia,ja
real*8::tol,olde,t(0:n-1),newe,eval(0:k-1),diff
complex*16::temp2(0:n-1),temp3(0:n-1),alpha,beta
complex*16::w1(0:n-1),evec(0:n-1,0:k-1)
complex*16, dimension (:,:), allocatable :: A,VP,temp,PPP,V,Z
complex*16,dimension (:), allocatable ::WORK
real*8,dimension (:), allocatable ::W,RWORK,IWORK,IFAIL
character::matdescra(6)
real*8::VL,VU,tol1
character*1::order
tol1=0.0000000000000001d0
matdescra(1)='H'
matdescra(2)='U'
matdescra(3)='N'
matdescra(4)='F'
alpha=dcmplx(1,0)
beta=dcmplx(0,0)
allocate(V(0:n-1,0:k-1))
do j=0,k-1
call random_number(t)
V(:,j)=t/norm2(t)
end do
call random_number(newe)
do m1=1,maxi,1
do m=k,b,k
allocate(A(0:(m-k)-1,0:m-k-1),W(0:m-k-1),WORK(0:2*(m-k)-1))
allocate(VP(0:n-1,0:m-k-1),ppp(0:m-k-1,0:m-k-1))
allocate(temp(0:n-1,0:m-k-1),RWORK(0:7*(m-k)-3))
allocate(IWORK(0:5*(m-k)-1),IFAIL(0:m-k-1),Z(0:m-k-1,0:k-1))
		if(m .gt. k) then
		olde=newe
		VP=V(:,0:m-k-1)
		deallocate(V)
		call QRC(n,m-k,VP)
		temp=0.0d0
                call mkl_zcoomm('n', n, m-k, n, alpha, &
        matdescra, H, ia, ja, nnz, VP, n, beta, temp, n)
		ppp=matmul(conjg(transpose(VP)),temp)
		A=ppp
		if (order=='l') then 
			IL=1
			IU=k
			xx=k
		end if
		if (order=='h') then 
			IL=m-k-k+1
			IU=m-k
			xx=1
		end if
		call zheevx( 'V', 'I', 'U', m-k, A, m-k, VL, VU,&
        IL,IU,tol1, MM, W, Z, m-k, WORK, 2*(m-k),&
            RWORK, IWORK, IFAIL, INFO )
		newe=W(xx-1)
		allocate(V(0:n-1,0:m-1))
		V(:,0:m-k-1)=VP
		do j=0,k-1
		temp3=matmul(VP,Z(:,j))
		V(:,j)=temp3
		call mkl_zcoomm('n', n, 1, n, alpha, matdescra,&
         H, ia, ja, nnz, temp3, n, beta, temp2, n)
		w1=temp2-W(j)*temp3
 		V(:,m-k+j)=w1
		end do
	end if
diff=olde-newe
if ( abs(diff) .le. tol) then
do i1=0,k-1,1
evec(:,i1)=matmul(VP,Z(:,i1))
eval(i1)=W(i1)
end do
end if
deallocate(ppp,temp,RWORK,A,W,WORK,VP,Z,IFAIL,IWORK)
end do
if ( abs(diff) .le. tol) then
exit
end if
end do
deallocate(V)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine OPDAVIDSONSR(HN1,n,k,order,b,tol,evec,eval,maxi)
implicit none
integer::k,n,INFO,m1,b,IL,IU,MM,i,j,m,i1,xx,maxi
real*8::tol,diff,olde,II(0:n-1,0:n-1),t(0:n-1),HN1(0:n-1,0:n-1)
real*8::newe,evec(0:n-1,0:k-1),eval(0:k-1)
real*8::temp2(0:n-1,0:n-1),temp3(0:n-1),w1(0:n-1)
real*8, dimension (:,:), allocatable :: A,VP,temp,PPP,V,Z
real*8,dimension (:), allocatable ::WORK,W,IWORK,IFAIL
real*8::VL,VU,tol1
character*1::order
tol1=0.0000000000000001d0
II=0
do i=0,n-1
II(i,i)=1.0d0
end do
allocate(V(0:n-1,0:k-1))
do j=0,k-1
call random_number(t)
V(:,j)=t/norm2(t)
end do
call random_number(newe)
do m1=1,maxi,1
do m=k,b,k
allocate(A(0:(m-k)-1,0:m-k-1),W(0:m-k-1),WORK(0:8*(m-k)-1))
allocate(temp(0:n-1,0:m-k-1),ppp(0:m-k-1,0:m-k-1))
allocate(IWORK(0:5*(m-k)-1),IFAIL(0:m-k-1))
allocate(Z(0:m-k-1,0:k-1),VP(0:n-1,0:m-k-1))
	if(m .gt. k) then
		olde=newe
		VP=V(:,0:m-k-1)
		deallocate(V)
		call QRR(n,m-k,VP)
		temp=matmul(HN1,VP)
		ppp=matmul(transpose(VP),temp)
		A=ppp
		if (order=='l') then 
			IL=1
			IU=k
			xx=k
		end if
		if (order=='h') then 
			IL=m-k-k+1
			IU=m-k
			xx=1
		end if
		call dsyevx('V','I','U',m-k,A,m-k,VL,VU,IL,IU, & 
                   tol1,MM,W,Z,m-k,WORK,8*(m-k),IWORK,IFAIL,INFO)
		newe=W(xx-1)
		allocate(V(0:n-1,0:m-1))
		V(:,k:m-k-1)=VP
		do j=0,k-1
			temp2=HN1-W(j)*II
			temp3=matmul(VP,Z(:,j))
			V(:,j)=temp3
			w1=matmul(temp2,temp3)
			V(:,m-k+j)=w1
		end do
	end if
diff=olde-newe
if ( abs(diff) .le. tol) then
do i1=0,k-1,1
evec(:,i1)=matmul(VP,A(:,i1))
eval(i1)=W(i1)
end do
end if
deallocate(ppp,temp,A,W,WORK,VP,Z,IFAIL,IWORK)
end do
if ( abs(diff) .le. tol) then
exit
end if
end do
deallocate(V)
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine OPDAVIDSONHC(HN1,n,k,order,b,tol,evec,eval,maxi)
implicit none
integer::k,n,INFO,b,i,j,i1,m,m1,IL,IU,MM,xx,maxi
real*8::tol,diff,olde,II(0:n-1,0:n-1),t(0:n-1),newe,eval(0:k-1)
complex*16::HN1(0:n-1,0:n-1),temp2(0:n-1,0:n-1),temp3(0:n-1)
complex*16::w1(0:n-1),evec(0:n-1,0:k-1)
complex*16, dimension (:,:), allocatable :: A,VP,temp,PPP,V,Z
complex*16,dimension (:), allocatable ::WORK
real*8,dimension (:), allocatable ::W,RWORK,IWORK,IFAIL
real*8::VL,VU,tol1
character*1::order
tol1=0.0000000000000001d0
II=0
do i=0,n-1
II(i,i)=1.0d0
end do
allocate(V(0:n-1,0:k-1))
do j=0,k-1
call random_number(t)
V(:,j)=t/norm2(t)
end do
call random_number(newe)
do m1=1,maxi,1
do m=k,b,k
allocate(A(0:(m-k)-1,0:m-k-1),W(0:m-k-1),WORK(0:2*(m-k)-1))
allocate(VP(0:n-1,0:m-k-1),ppp(0:m-k-1,0:m-k-1))
allocate(temp(0:n-1,0:m-k-1),RWORK(0:5*(m-k)-1))
allocate(IWORK(0:5*(m-k)-1),IFAIL(0:m-k-1),Z(0:m-k-1,0:k-1))
	if(m .gt. k) then
		olde=newe
		VP=V(:,0:m-k-1)
		deallocate(V)
		call QRC(n,m-k,VP)
		temp=matmul(HN1,VP)
		ppp=matmul(conjg(transpose(VP)),temp)
		A=ppp
		if (order=='l') then 
			IL=1
			IU=k
			xx=k
		end if
		if (order=='h') then 
			IL=m-k-k+1
			IU=m-k
			xx=1
		end if
		call zheevx( 'V', 'I', 'U', m-k, A, m-k, VL,&
        VU, IL, IU,tol1, MM, W, Z, m-k, WORK, & 
           2*(m-k),RWORK,IWORK, IFAIL, INFO )
		newe=W(xx-1)
		allocate(V(0:n-1,0:m-1))
		V(:,k:m-k-1)=VP
		do j=0,k-1
			temp2=HN1-W(j)*II
			temp3=matmul(VP,Z(:,j))
			V(:,j)=temp3
			w1=matmul(temp2,temp3)
 			V(:,m-k+j)=w1
		end do
	end if
diff=olde-newe
if (abs(diff) .le. tol) then
do i1=0,k-1,1
evec(:,i1)=matmul(VP,Z(:,i1))
eval(i1)=W(i1)
end do
end if
deallocate(ppp,temp,RWORK,A,W,WORK,VP,Z,IFAIL,IWORK)
end do
if (abs(diff) .le. tol) then
exit
end if
end do
deallocate(V)
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine QRR(n,k,wi)
integer::n,k
real*8::wi(0:n-1,0:k-1)
integer::i,i1,INFO
real*8::TAU(0:k-1),WORK(0:k-1)
call DGEQRF( n, k, wi, n, TAU, WORK, k, INFO )
call DORGQR( n, k, k, wi, n, TAU, WORK, n, INFO )
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine QRC(s2,k,wi)
integer::s2
complex*16::wi(0:s2-1,0:k-1)
integer::i,i1,INFO
complex*16::TAU(0:k-1),WORK(0:k-1)
call ZGEQRF( s2, k, wi,s2, TAU, WORK, k, INFO )
call ZUNGQR( s2, k, k, wi, s2, TAU, WORK, s2, INFO )
end subroutine


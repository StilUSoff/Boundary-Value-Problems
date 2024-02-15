module methods_test
contains

character(len=20) function str(k) !function for translating integer to string
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)

end function


subroutine sweep(Im,A,B,C,D,F)
implicit none
integer, intent(in):: Im
real*16, dimension(1:Im), intent(in):: A, B, C, D
real*16, dimension(1:Im), intent(out):: F
real*16, dimension(1:Im):: alpha, beta
real*16 :: k0
integer::i
!Прямой ход
alpha(1) = -A(1) / B(1)
beta(1) = -D(1) / B(1)
do i = 2, (Im-1)
    k0 = (B(i) + C(i)*alpha(i-1))
    alpha(i) = -A(i) / k0
    beta(i) = -(D(i) + C(i)*beta(i-1)) / k0
end do
!Обратный ход
F(Im) = -(D(Im) + C(Im)*beta(Im-1)) / (B(Im) + C(Im)*alpha(Im-1))
do i = (Im-1), 1, -1
    F(i) = alpha(i)*F(i+1) + beta(i)
end do

end subroutine


subroutine finit_differences(X, Fo_final, Bi, p)
implicit none
real*16 Bi, Fo, Fo_final, delt_X, delt_Fo, one, X
integer :: i, j, n,p, x_i, Fo_i
real*16, allocatable, dimension (:) :: AA,CC,BB,DD, O_n

one=1
delt_X=one/(p-1)
delt_Fo=(delt_X**2)/2

allocate(O_n(p+1))

do i=1,p
    O_n(i)=1
end do

allocate(AA(p),BB(p),CC(p),DD(p))

do i=2,p-1
    CC(i)=(1.0)/(delt_X**2)
    AA(i)=CC(i)
    BB(i) = - ((1.0)/(delt_Fo)) - ((2.0)/(delt_X**2))
end do
    CC(1)=0
    AA(1)=(2.0)/(delt_X**2)
    BB(1) = - ((1.0)/(delt_Fo)) - ((2.0)/(delt_X**2))
    CC(p)=(2.0)/(delt_X**2)
    AA(p)=0
    BB(p) = - ((1.0)/(delt_Fo)) - ((2.0)/(delt_X**2)) - ((2.0*Bi)/(delt_X))




Fo=delt_Fo
Fo_i=0
do n=1,100000 !time

    do i=1,p 
        DD(i)=O_n(i)/delt_Fo
    end do

    call sweep(p,AA,BB,CC,DD,O_n)

    Fo=Fo+delt_Fo
    if (Fo>Fo_final) exit
end do
if (X==0) then
    x_i=1
elseif (X==1) then
    x_i=p
end if
write(*,'(A13,$)')'Your theta = '
write(*,*) O_n(x_i)

end subroutine

end module 

program main
use methods_test
implicit none
real*16 X, Fo, Bi
integer p

X=0
p=81
Fo=4.0
Bi=1.6

call finit_differences(X, Fo, Bi, p)

end program
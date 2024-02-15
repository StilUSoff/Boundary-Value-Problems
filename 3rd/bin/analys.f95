module functions_temp_find
contains 

real*16 Function f(Mu, Bi) result(res)
real*16 :: Mu, Bi
res=(cos(Mu))/(sin(Mu))-(Mu/Bi)
end function

real*16 Function A_f(Mu) result(res)
real*16 :: Mu
res=(2*sin(Mu))/(Mu+sin(Mu)*cos(Mu))
end function

real*16 Function O_f(Mu, A, Fo, X) result(res)
real*16 :: Mu, A, Fo, X
res=A*exp( (-1)*Fo*(Mu**2))*cos(Mu*X)
end function

end module

module methods_temp_find
contains

real*16 Function bis(n, E, Bi) result(res)
use functions_temp_find
real*16 :: E, a_Mu, b_Mu, c_Mu, f1, f2, check, pi=3.1415926535897932, Bi, delt, n_r
integer :: n, i

n_r=real(n, 16)
a_Mu=pi*(n_r-1)
b_Mu=pi*(n_r-0.5)
f1=f(b_Mu, Bi)
do i=1,1000
    c_Mu=(a_Mu+b_Mu)/2.0
    f2=f(c_Mu, Bi)
    check=f1*f2
    if (check>0) then 
        b_Mu=c_Mu
        f1=f2
    elseif (check<0) then
        a_Mu=c_Mu
    elseif (check==0) then
        a_Mu=c_Mu
        b_Mu=c_Mu
    end if
    delt=abs(a_Mu-b_Mu)
    if (delt<E) exit
end do
res=c_Mu
end function

subroutine Fourier(X, Fo, Bi, E)
use functions_temp_find
implicit none
real*16 O, Bi, Mu, A, Fo, X, E, O_n, check,te,t0
integer :: i, n

n=0
O=0
O_n=0
check=1
do i=1,1000
    n=n+1
    Mu=bis(n, E, Bi)
    A=A_f(Mu)
    O_n=O_f(Mu, A, Fo, X)
    O=O+O_n
    check=abs(O)
    if (check<E) exit
end do
Te=130.0+273.15
T0=500.0+273.15
write(*,*)'T=', O

end subroutine

end module 

program main
use functions_temp_find
use methods_temp_find
implicit none
real*16 X, Fo, Bi, E

write(*,*)'Fo=1'
X=0
Fo=1.0
Bi=0.095
call Fourier(X, Fo, Bi, E)

Bi=1.09
call Fourier(X, Fo, Bi, E)

Bi=68.2
call Fourier(X, Fo, Bi, E)

write(*,*)'Fo=1.5'
X=0
Fo=1.5
Bi=0.095
call Fourier(X, Fo, Bi, E)

Bi=1.09
call Fourier(X, Fo, Bi, E)

Bi=68.2
call Fourier(X, Fo, Bi, E)

end program
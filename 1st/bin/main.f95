module functions
contains 
 
real*16 Function f(Mu, Bi) result(res)
real*16 :: Mu, Bi
 
res=(cos(Mu))/(sin(Mu))-(Mu/Bi) !function for bisection
 
end function
 
real*16 Function A_f(Mu) result(res)
real*16 :: Mu
 
res=(2*sin(Mu))/(Mu+sin(Mu)*cos(Mu)) !function for A
 
end function
 
real*16 Function O_f(Mu, A, Fo, X) result(res)
real*16 :: Mu, A, Fo, X
 
res=A*exp( (-1)*Fo*(Mu**2))*cos(Mu*X) !function for theta
 
end function
 
end module
 
module methods
contains
 
real*16 Function bis(n, E, Bi) result(res) !bisection
use functions
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
res=c_Mu  !saving result
end function
 
character(len=20) function str(k) !function for translating integer to string
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
 
end function
 
subroutine Fourier(check1, Te, T0, delta, alpha, lambda, cp, p, E)
use functions
implicit none
real*16 Te, T0, delta, lambda, cp, p, O, Bi, Mu, A, Fo, X, Fo_max, d, check, E, en, param, st, time_01
real*4 Bi_wr
integer :: i, j, n, check1, k, u, h, alpha, start, Fo_points
real*16, allocatable, dimension (:) :: O_n, T, tau, xx,tau_points
 
Bi=(alpha*delta)/(2*lambda)
write(*,'(F7.4)') Bi
if (Bi<0.5) then !Fo_max based on Bi
    Bi_wr=anint(Bi*10000)/10000
elseif (Bi>10) then
    Bi_wr=anint(Bi*10)/10
else
    Bi_wr=anint(Bi*1000)/1000
end if
if (Bi<1.25) then
    Fo_max=3.11*Bi**(-0.81)
elseif (Bi>20) then
    Fo_max=1.1
else
    Fo_max=2.76*Bi**(-0.31)
end if
if (check1==0) then
    d=Fo_max/100
    en=Fo_max
    write(13,*) alpha, Bi_wr
else
    d=1.0/99
    en=1.0
    write(14,*) alpha, Bi_wr
end if
 
if (check1==0) then
    allocate(tau_points(100))
    open(15, file='bin/res/time/fo_points_alpha'//trim(str(alpha))//'.txt') 
    do Fo_points=1,100
        if (Fo_points>Fo_max) exit
        tau_points(Fo_points)=(Fo_points*cp*p*(delta**2))/(4*lambda)
        write(15,*) Fo_points, tau_points(Fo_points), 400
    end do
    deallocate(tau_points)
end if
 
do k=1,3 !do for 3 points (x=0, x=d/4 and x=d/2 or Fo=0.1*Fo_max, Fo=0.5*Fo_max or Fo=0.9*Fo_max)  
    allocate(O_n(100))
    allocate(T(100))
    T(1)=(T0-Te)+Te
    if (check1==0) then !opening files for results
        allocate(tau(100))
        tau(1)=0
        open(10, file='bin/res/time/alpha'//trim(str(alpha))//'_point'//trim(str(k))//'_temp.txt')
        open(11, file='bin/res/time/alpha'//trim(str(alpha))//'_point'//trim(str(k))//'_time.txt')
        write(10,*) T(1)
        write(11,*) tau(1)
        start=2
    else
        allocate(xx(100))
        open(10, file='bin/res/coord/alpha'//trim(str(alpha))//'_point'//trim(str(k))//'_temp.txt')
        open(12, file='bin/res/coord/alpha'//trim(str(alpha))//'_point'//trim(str(k))//'_coord.txt')
        start=1
    end if
 
    if (k==1)then !param=static time or coordinate 
        if (check1==0) then
            param=(2*0)/delta
            write(*,'(a40,$)')'for point    (x=0)    iterations='
        else
            param=Fo_max*0.1
            write(*,'(a40,$)')'for point (0.1*Fo_max) iterations='
        end if
    elseif (k==2) then
        if (check1==0) then
            param=(2*(delta/4))/delta
            write(*,'(a40,$)')'for point (x=delta/4) iterations='
        else
            param=Fo_max*0.5
            write(*,'(a40,$)')'for point (0.5*Fo_max) iterations='
        end if
    elseif (k==3) then
        if (check1==0) then
            param=(2*(delta/2))/delta
            write(*,'(a40,$)')'for point (x=delta/2) iterations='
        else
            param=Fo_max*0.9
            write(*,'(a40,$)')'for point (0.9*Fo_max) iterations='
        end if
    end if
 
    st=0
    if (check1==1) then
        st=-d
    end if
    do j=start,100
        st=st+d !st=dynamic time or coordinate
        n=0
        O=0
        O_n(j)=0
        check=1
        do h=1,1000
            n=n+1
            Mu=bis(n, E, Bi)
            A=A_f(Mu)
            if (check1==0) then
                O=O_f(Mu, A, st, param)
            else
                O=O_f(Mu, A, param, st)
            end if
            O_n(j)=O_n(j)+O
            check=abs(O)
            if (check<E) then
                n=n-1
                exit
            end if
        end do
        T(j)=O_n(j)*(T0-Te)+Te
        if (check1==0) then
            tau(j)=(st*cp*p*(delta**2))/(4*lambda) !real time
            write(10,*) T(j)
            write(11,*) tau(j)
        else
            xx(j)=(st*delta)/2.0 !real coordinate 
            write(10,*) T(j)
            write(12,*) xx(j)
        end if
        if (st>en) exit
    end do
    close(10)
    if (check1==0) then
        close(11)
    else
        close(12)
    end if
    deallocate(O_n)
    deallocate(T)
    if (check1==0) then
        deallocate(tau)
    else
        deallocate(xx)
    end if
    write(*,'(i2)') n !number of iterations for |theta|<eps
end do
 
end subroutine
 
end module 
 
program main
use functions
use methods
implicit none
 
real*16 Te, T0, delta, lambda, cp, p, E
integer, allocatable, dimension (:) :: alpha_n
integer :: check1, n, i, alpha
 
Te=130.0+273.15
T0=500.0+273.15
delta=0.6
lambda=110.0 
cp=380.0
p=8600.0
E=0.000000000001
 
n=3
allocate(alpha_n(n))
    
alpha_n(1)=35
alpha_n(2)=400
alpha_n(3)=25000
 
open(13, file='bin/res/time/bi_alpha.txt')
open(14, file='bin/res/coord/bi_alpha.txt')
 
check1=1 !static Fo
write(*,*) 'T(x)'
write(*,*)
do i=1,n
    alpha=alpha_n(i)
    write(*,'(A7,$)') 'Bi='
    call Fourier(check1, Te, T0, delta, alpha, lambda, cp, p, E)
    write(*,*)
end do
 
check1=0 !static X
write(*,*) 'T(t)'
write(*,*)
alpha=alpha_n(1)
write(*,'(A7,$)') 'Bi='
call Fourier(check1, Te, T0, delta, alpha, lambda, cp, p, E)
write(*,*)
alpha=alpha_n(n)
write(*,'(A7,$)') 'Bi='
call Fourier(check1, Te, T0, delta, alpha, lambda, cp, p, E)
 
close(13)
close(14)
 
end program

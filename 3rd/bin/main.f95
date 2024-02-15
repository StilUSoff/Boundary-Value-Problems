module methods
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


subroutine finit_differences(p, Te, T0, delta, alpha, lambda, cp, po, check,t_res)
implicit none
real*16 Te, T0, delta, lambda, cp, po, Bi, Fo, Fo_max, st, delt_X, delt_Fo, one,delt_Fo_2,t_res
real*4 Bi_wr
integer :: i, j, n, alpha, p, p_Fo, Fo_points, popul, check, t_1
real*16, allocatable, dimension (:,:) ::tau_x
real*16, allocatable, dimension (:) :: points, AA,CC,BB,DD, O_n,T
integer, allocatable, dimension (:) :: points_dot

Bi=(alpha*delta)/(2*lambda)
write(*,*)Bi
if (Bi<=0.5) then !Fo_max based on Bi
    write(13,'(F5.3,$)') Bi
    write(*,'(F5.3)') Bi
elseif (Bi>=10) then
    write(13,'(F4.1,$)') Bi
    write(*,'(F4.1)') Bi
else
    write(13,'(F4.2,$)') Bi
    write(*,'(F4.2)') Bi
end if

if (Bi<1.25) then
    Fo_max=3.11*Bi**(-0.81)
elseif (Bi>20) then
    Fo_max=1.1
else
    Fo_max=2.76*Bi**(-0.31)
end if
write(*,'(A9,$)')'Fo_max= '
write(*,'(F12.6)')(Fo_max*cp*po*(delta**2))/(4*lambda)

allocate(points(6))
allocate(points_dot(6))
points(1)=0
points(2)=0.5 !delta/4
points(3)=1 !delta/2
points(4)=Fo_max*0.1 !(Fo_max*0.1*cp*po*(delta**2))/(4*lambda)
points(5)=Fo_max*0.5 !(Fo_max*0.5*cp*po*(delta**2))/(4*lambda)
points(6)=Fo_max*0.9 !(Fo_max*0.9*cp*po*(delta**2))/(4*lambda)

one=1
delt_X=one/(p-1)
delt_Fo_2=(delt_X**2)/2
delt_Fo=delt_Fo_2
if (check==1) delt_Fo=delt_Fo_2*1
if (check==2) delt_Fo=delt_Fo_2*5
if (check==3) delt_Fo=delt_Fo_2*20
if (check==4) delt_Fo=delt_Fo_2*100

p_Fo=int(Fo_max/delt_Fo)
write(13,*) p, p_Fo, alpha

allocate(O_n(p+1))
allocate(T(p+1))
allocate(tau_x(p_Fo+1,2))

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

st=0
popul=0
do i=1,p
    popul=popul+1
    O_n(i)=1
    T(i)=T0
    tau_x(i,2)=(st*delta)/2.0
    st=st+delt_X
    do j=1,3
        if ((points(j)>st-delt_x) .and. (points(j)<st+delt_x)) then
            points_dot(j)=popul
        end if
    end do
end do
tau_x(1,1)=0

if (check==0) then
    open(10, file='bin/res/Im='//trim(str(p))//'/alpha'//trim(str(alpha))//'_temp_t.txt')
    open(101, file='bin/res/Im='//trim(str(p))//'/alpha'//trim(str(alpha))//'_temp_x.txt')
    open(11, file='bin/res/Im='//trim(str(p))//'/alpha'//trim(str(alpha))//'_time.txt')
    open(21, file='bin/res/Im='//trim(str(p))//'/alpha'//trim(str(alpha))//'_coord.txt')
else
    open(10, file='bin/res/Im='//trim(str(p))//'_def/alpha'//trim(str(alpha))//'_dFo'//trim(str(check))//'_temp_t.txt')
    open(101, file='bin/res/Im='//trim(str(p))//'_def/alpha'//trim(str(alpha))//'_dFo'//trim(str(check))//'_temp_x.txt')
    open(11, file='bin/res/Im='//trim(str(p))//'_def/alpha'//trim(str(alpha))//'_dFo'//trim(str(check))//'_time.txt')
    open(21, file='bin/res/Im='//trim(str(p))//'_def/alpha'//trim(str(alpha))//'_dFo'//trim(str(check))//'_coord.txt')
end if

Fo=delt_Fo
n=0
tau_x(1,1)=0
do n=1,p_Fo+1 !time
    do i=4,6
        if ((points(i)>Fo-delt_fo) .and. (points(i)<Fo+delt_fo)) then
            points_dot(i)=n
        end if
    end do

    do i=1,p 
        DD(i)=O_n(i)/delt_Fo
    end do

    call sweep(p,AA,BB,CC,DD,O_n)

    do i=1,p
        T(i)=O_n(i)*(T0-Te)+Te
    end do
    tau_x(n+1,1)=(Fo*cp*po*(delta**2))/(4*lambda) !real time
    write(10,*) T(1), T(points_dot(2)), T(points_dot(3))
    write(10,*)

    do i=4,6
        if ((points(i)>Fo-delt_fo/2) .and. (points(i)<Fo+delt_fo/2)) then
            do j=1,p
                write(101,'(f18.10,$)') T(j)
            end do
            write(101,*)
        end if
    end do
    if ((abs(Fo-1.0)<=(delt_Fo/2.0)) .and. (check==0)) then
        t_1=n
        write(*,'(A10,$)')'For iter '
        write(*,'(i5,$)')t_1
        write(*,'(A11,$)')' and time '
        write(*,*) Fo, ':    T(0,Fo=1)=', O_n(1), (O_n(1)-t_res), ((O_n(1)-t_res)/t_res)*100, Bi
    elseif ((abs(Fo-1.5)<=(delt_Fo/2.0)) .and. (check/=0)) then
        t_1=n
        write(*,'(A10,$)')'For iter '
        write(*,'(i5,$)')t_1
        write(*,'(A11,$)')' and time '
        write(*,*) Fo, ':    T(0,Fo=1.5)=', O_n(1), (O_n(1)-0.345593785639), ((O_n(1)-0.345593785639)/0.345593785639)*100, Bi
    end if
    Fo=Fo+delt_Fo
    if (Fo>=Fo_max) exit
end do
do i=1,p_Fo
    write(11,*) tau_x(i,1)
end do
do i=1,p
    write(21,*) tau_x(i,2)
end do

write(20,*) 1, points_dot(2), points_dot(3), points_dot(4), points_dot(5), points_dot(6)
deallocate(points)
deallocate(points_dot)

close(10)
close(11)
close(101)
deallocate(O_n)
deallocate(T)
deallocate(tau_x)
deallocate(AA,BB,CC,DD)

end subroutine

end module 

program main
use methods
implicit none

real*16 Te, T0, delta, lambda, cp, po, E, t_res,T_res_n(3)
integer, allocatable, dimension (:) :: alpha_n, p_n, tau_points
integer :: n, i, alpha, p, j, check, Fo_points

Te=130.0+273.15
T0=500.0+273.15
delta=0.6
lambda=110.0 
cp=380.0
po=8600.0

allocate(p_n(5))! grid size
allocate(alpha_n(3))
allocate(tau_points(100))
alpha_n(1)=35
alpha_n(2)=400
alpha_n(3)=25000
p_n(1)=11
p_n(2)=21
p_n(3)=41
p_n(4)=81
p_n(5)=101

 T_res_n(1)= 0.926012747427353802229742345141718148      
 T_res_n(2)= 0.512375769314193418196748998307030519      
 T_res_n(3)= 0.115867936825182547473269331561302031  

write(*,*)'start program'

write(*,*) 'Main + Research 2'
check=0
do i=1,5
    p=p_n(i)
    write(*,*) '_________________________________________________________________________________________________________'
    write(*,'(a4,$)') ' Im='
    write(*,'(i3)') p
    write(*,*)
    open(13, file='bin/res/Im='//trim(str(p))//'/bi_iFORcoord_iFORtime_alpha.txt')
    open(20, file='bin/res/Im='//trim(str(p))//'/points.txt')
    do j=1,3
        alpha=alpha_n(j)
        t_res=T_res_n(j)
        write(*,'(A4,$)') ' Bi='
        call finit_differences(p, Te, T0, delta, alpha, lambda, cp, po,check,t_res)
        write(*,*)
    end do
    close(20)
    close(13)
    write(*,*) 'this Im done'
end do

write(*,*)
write(*,*) 'Research 1'
alpha=alpha_n(2)
do check=1,4
    p=21
    write(*,*) '_________________________________________________________________________________________________________'
        write(*,'(a4,$)') ' Im='
        write(*,'(i3)') p
        write(*,*)
        write(*,'(a12,$)') ' num of dFo:'
        write(*,'(i3)') check
        write(*,*)
    open(13, file='bin/res/Im='//trim(str(p))//'_def/bi_iFORcoord_iFORtime_alpha.txt')
    open(20, file='bin/res/Im='//trim(str(p))//'_def/points.txt')
    write(*,'(A4,$)') ' Bi='
    call finit_differences(p, Te, T0, delta, alpha, lambda, cp, po,check,t_res)
    write(*,*)
    close(20)
    close(13)
end do
write(*,*) 'this res done'

open(29, file='bin/res/fo_points.txt') 
do Fo_points=1,100
    tau_points(Fo_points)=(Fo_points*cp*po*(delta**2))/(4*lambda)
    write(29,*) Fo_points, tau_points(Fo_points), 400    
end do
close(29)

write(*,*) 'end program'

end program
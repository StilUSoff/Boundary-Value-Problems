module methods
contains

character(len=20) function str(k) !function for translating integer to string
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)

end function

subroutine finit_differences(p, Te, T0, delta, alpha, lambda, cp, po, check,T_res)
implicit none
real*16 Te, T0, delta, lambda, cp, po, Bi, Fo, Fo_max, st, delt_X, delt_Fo, one,T_res
real*4 Bi_wr
integer :: i, j, n, alpha, p, p_Fo, Fo_points, popul, check, t_1
real*16, allocatable, dimension (:,:) :: O_n,T,tau_x
real*16, allocatable, dimension (:) :: points, tau_points
integer, allocatable, dimension (:) :: points_dot

Bi=(alpha*delta)/(2*lambda)
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
delt_Fo=(delt_X**2)/2
if (check==1) delt_Fo=delt_Fo*1.1
p_Fo=int(Fo_max/delt_Fo)
write(13,*) p, p_Fo, alpha

allocate(O_n(p+1,2))
allocate(T(p+1,2))
allocate(tau_x(p_Fo+1,2))

st=0
popul=0
do i=1,p
    popul=popul+1
    O_n(i,2)=1
    T(i,2)=(T0-Te)+Te
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
    open(10, file='bin/res/Im='//trim(str(p))//'_10_perc/alpha'//trim(str(alpha))//'_temp_t.txt')
    open(101, file='bin/res/Im='//trim(str(p))//'_10_perc/alpha'//trim(str(alpha))//'_temp_x.txt')
    open(11, file='bin/res/Im='//trim(str(p))//'_10_perc/alpha'//trim(str(alpha))//'_time.txt')
    open(21, file='bin/res/Im='//trim(str(p))//'_10_perc/alpha'//trim(str(alpha))//'_coord.txt')
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
    O_n(:,1)=O_n(:,2)
    T(:,1)=T(:,2)
    do i=2,p-1
        O_n(i,2)=((delt_Fo)/(delt_X**2))*(O_n(i+1,1)+O_n(i-1,1))+O_n(i,1)*(1-(2*(delt_Fo)/(delt_X**2)))
        T(i,2)=O_n(i,2)*(T0-Te)+Te
    end do
    O_n(1,2)=O_n(2,2)
    O_n(p,2)=(O_n(p-1,2))/(1+(Bi*delt_X))
    T(1,2)=O_n(1,2)*(T0-Te)+Te
    T(p,2)=O_n(p,2)*(T0-Te)+Te
    tau_x(n+1,1)=(Fo*cp*po*(delta**2))/(4*lambda) !real time
    write(10,*) T(1,1), T(points_dot(2),1), T(points_dot(3),1)
    write(10,*)

    do i=4,6
        if ((points(i)>Fo-delt_fo/2) .and. (points(i)<Fo+delt_fo/2)) then
        do j=1,p
            write(101,'(f18.10,$)') T(j,1)
        end do
        write(101,*)
        end if
    end do
    if ((abs(Fo-1.0)<=(delt_Fo/2.0)) .and. (check==0)) then
        t_1=n
        write(*,'(A10,$)')'For iter '
        write(*,'(i5,$)')t_1
        write(*,'(A11,$)')' and time '
        write(*,*) Fo, ':    T(0,Fo=1)=', O_n(1,1), (O_n(1,1)-t_res), ((O_n(1,1)-t_res)/t_res)*100

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
allocate(tau_points(100))
if (alpha/=400) then 
    if (check==0) then
        open(19, file='bin/res/Im='//trim(str(p))//'/fo_points_alpha'//trim(str(alpha))//'.txt') 
    else
        open(19, file='bin/res/Im='//trim(str(p))//'_10_perc/fo_points_alpha'//trim(str(alpha))//'.txt') 
    end if
    do Fo_points=1,100
        tau_points(Fo_points)=(Fo_points*cp*po*(delta**2))/(4*lambda)
        write(19,*) Fo_points, tau_points(Fo_points)    
    end do
end if
deallocate(tau_points)
close(19)

write(20,*) 1, points_dot(2), points_dot(3), points_dot(4), points_dot(5), points_dot(6)
deallocate(points)
deallocate(points_dot)

close(10)
close(11)
close(101)
deallocate(O_n)
deallocate(T)
deallocate(tau_x)

end subroutine

end module 

program main
use methods
implicit none

real*16 Te, T0, delta, lambda, cp, po, E,T_res,T_res_n(3)
integer, allocatable, dimension (:) :: alpha_n, p_n
real*16, allocatable, dimension (:,:) :: points
integer :: n, i, alpha, p, j, check=0

Te=130.0+273.15
T0=500.0+273.15
delta=0.6
lambda=110.0 
cp=380.0
po=8600.0

allocate(p_n(5))! grid size
allocate(alpha_n(3))
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
do i=1,5
    p=p_n(i)
    write(*,*) '_______________________'
    write(*,'(a3,$)') 'Im='
    write(*,'(i3)') p
    write(*,*)
    open(13, file='bin/res/Im='//trim(str(p))//'/bi_iFORcoord_iFORtime_alpha.txt')
    open(20, file='bin/res/Im='//trim(str(p))//'/points.txt')
    do j=1,3
        alpha=alpha_n(j)
        T_res=T_res_n(j)
        write(*,'(A7,$)') 'Bi='
        call finit_differences(p, Te, T0, delta, alpha, lambda, cp, po,check,T_res)
        write(*,*)
    end do
    close(20)
    close(13)
    write(*,*) 'this Im done'
end do

p=21
write(*,*) '_______________________'
write(*,*) 'Fluct'
    write(*,'(a3,$)') 'Im='
    write(*,'(i3)') p
    write(*,*)
open(13, file='bin/res/Im='//trim(str(p))//'_10_perc/bi_iFORcoord_iFORtime_alpha.txt')
open(20, file='bin/res/Im='//trim(str(p))//'_10_perc/points.txt')
check=1
j=3
alpha=alpha_n(j)
write(*,'(A7,$)') 'Bi='
call finit_differences(p, Te, T0, delta, alpha, lambda, cp, po,check,T_res)
write(*,*)
close(20)
close(13)
write(*,*) 'this Im done'
write(*,*) 'end program'

end program


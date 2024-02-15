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


subroutine grid(L, Im, koeff_d, lambda_0, x, xx, delt_tau)
implicit none
real*16 :: L, koeff_d, lambda_0, alpha_1, alpha_2, T0, Te1, Te2, po, cp, delt_tau
real*16, dimension (Im) :: delta_x, x, xx, lambda
integer :: Im, i, n, Im_star

Im_star=(Im+1.0)/2.0
x(1)=0.0
x(Im)=L
delta_x(2)=(L*(koeff_d-1.0))/(2.0*((koeff_d**(Im_star-1.0))-1.0))
do i=2,Im_star
    if (i/=2) delta_x(i)=delta_x(i-1)*koeff_d
    x(i)=x(i-1)+delta_x(i)
end do
do i=Im_star+1,Im
    x(i)=x(i-1)+delta_x(Im+2-i)
end do
do i=2,Im
    xx(i)=(x(i)+x(i-1))/2.0
end do
do i=2,Im-1
    x(i)=(xx(i+1)+xx(i))/2.0
end do
delt_tau=(L**2)/(2*lambda_0*((Im-1)**2))

end subroutine


subroutine coefs(L, Im, koeff_d, Qv, lambda_0, koeff_t, alpha_1, alpha_2, T0, Te1, Te2, po, cp, eps)
implicit none
real*16 :: L, koeff_d, Qv, lambda_0, koeff_t, alpha_1, alpha_2, T0, Te1, Te2, po, cp
real*16 :: eps, delt_tau, balance
real*16, dimension (Im) :: x, xx, lambda, T, q,  AA, BB, CC, DD
integer :: Im, i, n, Im_star, check=0

call grid(L, Im, koeff_d, lambda_0, x, xx, delt_tau)

T=T0; AA=0; CC=0

do n=1,1000000
    
    lambda=lambda_0*(1.0+koeff_t*(T-T0))
    do i=2,Im-1
        AA(i)=(lambda(i+1)*lambda(i))/((xx(i+1)-xx(i))*(lambda(i)*(x(i+1)-xx(i+1))+lambda(i+1)*(xx(i+1)-x(i))))
        CC(i)=(lambda(i)*lambda(i-1))/((xx(i+1)-xx(i))*(lambda(i-1)*(x(i)-xx(i))+lambda(i)*(xx(i)-x(i-1))))
        BB(i)=-AA(i)-CC(i)-(1.0/delt_tau)
        DD(i)=(T(i)/delt_tau)+Qv
    end do
    AA(1)=(lambda(2)*lambda(1))/((xx(2)-x(1))*(lambda(1)*(x(2)-xx(2))+lambda(2)*(xx(2)-x(1))))
    BB(1)=-AA(1)-(1.0/delt_tau)-(alpha_1/(xx(2)-x(1)))
    DD(1)=(T(1)/delt_tau)+((alpha_1*Te1)/(xx(2)-x(1)))+Qv
    CC(Im)=(lambda(Im)*lambda(Im-1))/((x(Im)-xx(Im))*(lambda(Im-1)*(x(Im)-xx(Im))+lambda(Im)*(xx(Im)-x(Im-1))))
    BB(Im)=-CC(Im)-(1.0/delt_tau)-(alpha_2/(x(Im)-xx(Im)))
    DD(Im)=(T(Im)/delt_tau)+((alpha_2*Te2)/(x(Im)-xx(Im)))+Qv



    call sweep(Im,AA,BB,CC,DD,T)

    ! Pomeranc=((Qv*(L**2))/(minval(lambda)*abs(Te2-Te1)))
    ! if (Pomeranc>=2) then
    !     balance=abs(1.0-((alpha_2*(T(Im)-Te2)+alpha_1*(T(1)-Te1))/(Qv*L)))
    !     if (balance<=eps) then
    !         check=1
    !     end if
    ! else
    !     if (Te1<T(1)) then
    !         balance=abs(1.0-((Qv*L+alpha_1*(T(1)-Te1))/(alpha_2*(Te2-T(Im)))))
    !         if (balance<=eps) then
    !             check=1
    !         end if
    !     elseif (Te2<T(Im)) then
    !         balance=abs(1.0-((Qv*L+alpha_2*(T(Im)-Te2))/(alpha_1*(Te1-T(1)))))
    !         if (balance<=eps) then
    !             check=1
    !         end if
    !     end if
    ! end if

    ! if (check==1) exit

    balance=abs(1-(((alpha_2*(T(Im)-Te2))+alpha_1*(T(1)-Te1))/(Qv*L)))
    if (balance<=eps) exit

end do

write(16,*) T+273.15; write(17,*) lambda; write(14,*) x; write(15,*) xx
write(*,'(a12,$)') '  6) max T: ' ; write(*,'(f7.3,$)') maxval(T)+273.15; write(*,'(a19,$)') ' K ,  left wall T: ';
write(*,'(f7.3,$)') T(1)+273.15; write(*,'(a20,$)') ' K ,  right wall T: ' ; write(*,'(f7.3,$)') T(Im)+273.15; write(*,*) 'K'
write(*,'(a23,$)') '  8) heat dissipation: ' ; write(*,'(f8.3,$)') Qv*L ; write(*,*) 'W/m^2'
write(*,'(a16,$)') '  9) heat flow: ' ; write(*,'(f8.3,$)') (alpha_2*(T(Im)-Te2))+(alpha_1*(T(1)-Te1)); write(*,*) 'W/m^2'
write(*,'(a20,$)') ' 10) balance error: '; write(*,'(f14.12)') balance
write(*,'(a17,$)') ' left heat flow: '; write(*,'(f8.3,$)') alpha_1*(T(1)-Te1)
write(*,'(a20,$)') ' ,  right heat flow: '; write(*,'(f8.3,$)') alpha_2*(T(Im)-Te2)
write(*,'(a19,$)') ' ,  term balance: '; write(*,'(f14.12)') ((alpha_2*(T(Im)-Te2)+alpha_1*(T(1)-Te1))/(Qv*L))


end subroutine

end module 


program main
use methods
implicit none
real*16 :: L, koeff_d, Qv, lambda_0, koeff_t, alpha_1, alpha_2, T0, Te1, Te2, po, cp, eps
integer :: Im

open(201, file='bin/data.txt')
read(201,*) L; read(201,*) Im; read(201,*) koeff_d; read(201,*) Qv; read(201,*) lambda_0; read(201,*) koeff_t; read(201,*) alpha_1;
read(201,*) alpha_2; read(201,*) T0; read(201,*) Te1; read(201,*) Te2; read(201,*) po; read(201,*) cp; read(201,*) eps
close(201)

open(14, file='bin/results/x_i.txt'); open(15, file='bin/results/xx_i.txt')
open(16, file='bin/results/T.txt'); open(17, file='bin/results/lambda_i.txt')

write(*,*)'start program'; write(*,*); write(*,'(a4,$)') ' Im='; write(*,'(i3)') Im; write(*,*)

call coefs(L, Im, koeff_d, Qv, lambda_0, koeff_t, alpha_1, alpha_2, T0, Te1, Te2, po, cp, eps)

close(14); close(15); close(16); close(17)

write(*,*); write(*,*) 'end program'

end program
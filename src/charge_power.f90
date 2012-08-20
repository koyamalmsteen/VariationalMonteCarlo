!
! $Id: charge_power.f90,v 1.3 2003/11/18 08:05:41 k-yukino Exp $
!

#include "parameter.h"

subroutine charge_power(alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,tmp_charge0,&
     nume1,nume2,nume3_1,nume3_2,nume4,nume5,name,error)
  use global_variables
  implicit none

  include "mpif.h"
  
  real(8),intent(in) :: alpha_psi1,alpha_h_psi1,alpha_h_h_psi1
  real(8),dimension(max_number_xi),intent(in) :: tmp_charge0
  real(8),dimension(max_number_xi),intent(out) :: nume1,nume2,nume3_1,&
       nume3_2,nume4,nume5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
! init
  nume1=0 ; nume2=0 ; nume3_1=0 ; nume3_2=0 ; nume4=0 ; nume5=0
  name="charge_power" ; error=0
!
! $B%a%$%s%k!<%A%s(B
!
  do i_count=1,max_number_xi

!! <Charge> tmp_charge0$B$G7W;;$7$F$$$k(B
! $BCm0U(B($B$d$d$3$7$$(B)
     nume1(i_count)=nume1(i_count)+tmp_charge0(i_count) &
          *(alpha_psi1**2)

!! <H Charge>            (<Charge>$B$N7k2L$r;H$&(B)
     nume2(i_count)=nume2(i_count)+tmp_charge0(i_count) &
          *(alpha_h_psi1*alpha_psi1)

!! <H Charge H>           (<Charge>$B$N7k2L$r;H$&(B) 
     nume3_1(i_count)=nume3_1(i_count)+tmp_charge0(i_count) &
          *(alpha_h_psi1**2)
!
!
     if( POWER==2 ) then
!! <O H^2>   ---> <psi|O|alpha><alpha|H^2|psi>$B$H$9$k(B (<Charge>$B$N7k2L$r;H$&(B)
        nume3_2(i_count)=nume3_2(i_count)+tmp_charge0(i_count) &
             *(alpha_psi1*alpha_h_h_psi1)
          
!! <H O H^2> ---> <psi|HO|alpha><alpha|H^2|psi> (<Charge>$B$N7k2L$r;H$&(B)
        nume4(i_count)=nume4(i_count)+tmp_charge0(i_count) &
             *(alpha_h_psi1*alpha_h_h_psi1)

!! <H^2 O H^2> ---> <psi|H^2|alpha><alpha|OH^2|psi>
        nume5(i_count)=nume5(i_count)+tmp_charge0(i_count) &
             *(alpha_h_h_psi1**2)
     end if
  end do

end subroutine charge_power

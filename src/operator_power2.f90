!
! $Id: operator_power2.f90,v 1.1 2003/11/18 08:05:41 k-yukino Exp $ 
!
#include "parameter.h"

subroutine operator_power2(nume1,nume2,nume3_1,nume3_2,nume4,nume5,&
     pow2_c0,pow2_c1,pow2_c2,denomi_pow2,result,name,error)
  use global_variables
  implicit none

  real(8),dimension(max_number_xi),intent(in) :: nume1,nume2,nume3_1,nume3_2,nume4,nume5
  real(8),intent(in) :: pow2_c0,pow2_c1,pow2_c2
  real(8),intent(in) :: denomi_pow2
  real(8),dimension(max_number_xi),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count 
!
  result=0 ; name="operator_power2" ; error=0

  if( POWER/=2 ) then
     return
  else 
     if( denomi_pow2==0 ) then
        error=-1 ; error=0
     else
        do i_count=1,max_number_xi
           result(i_count)=&
                (nume5(i_count)*(pow2_c2**2) &
                + nume4(i_count)*(pow2_c1*pow2_c2) &
                + nume3_2(i_count)*(2*pow2_c0*pow2_c2) &
                + nume3_1(i_count)*(pow2_c1**2) &
                + nume2(i_count)*(2*pow2_c0*pow2_c1) &
                + nume1(i_count)*(pow2_c0**2) )/denomi_pow2
        end do
     end if
  end if

end subroutine operator_power2

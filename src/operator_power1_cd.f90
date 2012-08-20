!
! $Id: operator_power1_cd.f90,v 1.1 2003/11/18 08:05:41 k-yukino Exp $ 
!
#include "parameter.h"

subroutine operator_power1_cd(nume1,nume2,nume3_1,pow1_c0,pow1_c1,&
     denomi_pow1,result,name,error)
  use global_variables
  implicit none

  real(8),dimension(TOTAL_SITE_NUMBER),intent(in) :: nume1,nume2,nume3_1
  real(8),intent(in) :: pow1_c0,pow1_c1
  real(8),intent(in) :: denomi_pow1
  real(8),dimension(max_number_xi),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! 
  integer :: i_count
!
  result=0 ; name="operator_power1_cd" ; error=0

  if( POWER==0 ) then
     return
  else 
     if( denomi_pow1==0 ) then
        error=-1 ; return
     else        
        do i_count=1,TOTAL_SITE_NUMBER
           result(i_count)=(nume3_1(i_count)*(pow1_c1**2) &
                + nume2(i_count)*(2*pow1_c0*pow1_c1) &
                + nume1(i_count)*(pow1_c0**2) )/denomi_pow1
        end do
     end if
  end if
  
end subroutine operator_power1_cd

!
! $Id: operator_power1.f90,v 1.1 2003/11/18 08:05:41 k-yukino Exp $ 
!
#include "parameter.h"

subroutine operator_power1(nume1,nume2,nume3_1,pow1_c0,pow1_c1,&
     denomi_pow1,result,name,error)
  use global_variables
  implicit none

  real(8),dimension(max_number_xi),intent(in) :: nume1,nume2,nume3_1
  real(8),intent(in) :: pow1_c0,pow1_c1
  real(8),intent(in) :: denomi_pow1
  real(8),dimension(max_number_xi),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! 
  integer :: i_count
!
  result=0 ; name="operator_power1" ; error=0

  write(*,*) "nume1=",nume1
  write(*,*) "nume2=",nume2
  write(*,*) "nume3_1=",nume3_1
  write(*,*) "pow1_c0=",pow1_c0
  write(*,*) "pow1_c1=",pow1_c1

  stop

  if( POWER==0 ) then
     return
  else 
     if( denomi_pow1==0 ) then
        error=-1 ; return
     else        
        do i_count=1,max_number_xi
           result(i_count)=(nume3_1(i_count)*(pow1_c1**2) &
                /(dble(TOTAL_SITE_NUMBER)**2 ) &
                + nume2(i_count)*(2*pow1_c0*pow1_c1) &
                /dble(TOTAL_SITE_NUMBER) &
                + nume1(i_count)*(pow1_c0**2) &
                /(dble(TOTAL_SITE_NUMBER)**2) &
           )/denomi_pow1
        end do
     end if
  end if
  
end subroutine operator_power1

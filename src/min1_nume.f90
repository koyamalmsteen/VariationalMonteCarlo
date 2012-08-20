!
! $Id: min1_nume.f90,v 1.1 2004/02/15 10:52:35 k-yukino Exp $
!

#include "parameter.h"

subroutine min1_nume(h0,h1,h2,h3,pow1_c0,pow1_c1,result_pow1,name,error)
  use global_variables
  implicit none

  real(8),intent(in) :: h0,h1,h2,h3
  real(8),intent(out) :: pow1_c0,pow1_c1
  real(8),intent(out) :: result_pow1
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8),external :: evaluate_pow1
  integer :: i_count,j_count
  real(8) :: temp_result,min_pow1
!
  pow1_c0=0 ; pow1_c1=0 ; result_pow1=0
  name="min1_nume" ; error=0
!
  temp_result=0 ; min_pow1=99999

  do i_count=-100,100          ! C0
     do j_count=-100,100          ! C1
        temp_result=evaluate_pow1(h0,h1,h2,h3,dble(i_count*0.01),&
             dble(j_count*0.01))

        if( min_pow1>temp_result ) then
           min_pow1=temp_result
           pow1_c0=dble(i_count*0.01)
           pow1_c1=dble(j_count*0.01)
        end if
     end do
  end do

  result_pow1=min_pow1
  
end subroutine min1_nume


real(8) function evaluate_pow1(h0,h1,h2,h3,C0,C1)
  implicit none

  real(8),intent(in) :: h0,h1,h2,h3
  real(8),intent(in) :: C0,C1
  
  if( C0/=0 .and. C1/=0 ) then
     evaluate_pow1=( (C1**2)*h3+(2*C0*C1)*h2+(C0**2)*h1 ) &
          /( (C1**2)*h2+(2*C0*C1)*h1+(C0**2)*h0 )
  else
     evaluate_pow1=9999
  end if

end function evaluate_pow1


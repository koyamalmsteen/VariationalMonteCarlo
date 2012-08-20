!
! $Id: min1.f90,v 1.1 2004/02/15 10:52:35 k-yukino Exp $
!

#include "parameter.h"

subroutine min1(h0,h1,h2,h3,pow1_c0,pow1_c1,result_pow1,name,error)
  use global_variables
  implicit none

  real(8),intent(in) :: h0,h1,h2,h3
  real(8),intent(out) :: pow1_c0,pow1_c1
  real(8),intent(out) :: result_pow1 
  character(32),intent(out) :: name
  integer,intent(out) :: error
!
  pow1_c1=0
  name="min1" ; error=0
  result_pow1=0

  if( MINIMUM==0 ) then
     call min1_lag(h0,h1,h2,h3,pow1_c0,pow1_c1,result_pow1,name,error)
     call error_check(name,error)
  else
     call min1_nume(h0,h1,h2,h3,pow1_c0,pow1_c1,result_pow1,name,error)
     call error_check(name,error)
  end if

end subroutine min1

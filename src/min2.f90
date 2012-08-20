!
! $Id: min2.f90,v 1.1 2004/02/15 10:52:35 k-yukino Exp $
!

#include "parameter.h"

subroutine min2(h0,h1,h2,h3,h4,h5,pow2_c0,pow2_c1,pow2_c2,&
     result_pow2,name,error)
  use global_variables
  implicit none

  real(8),intent(in) :: h0,h1,h2,h3,h4,h5
  real(8),intent(out) :: pow2_c0,pow2_c1,pow2_c2
  real(8),intent(out) :: result_pow2
  character(32),intent(out) :: name
  integer,intent(out) :: error
!
  pow2_c0=0 ; pow2_c1=0 ; pow2_c2=0
  name="min2" ; error=0
  result_pow2=0

  if( MINIMUM==0 ) then
     call min2_lag(h0,h1,h2,h3,h4,h5,pow2_c0,pow2_c1,&
          pow2_c2,result_pow2,name,error)
     call error_check(name,error)
  else
     call min2_nume(h0,h1,h2,h3,h4,h5,pow2_c0,pow2_c1,&
          pow2_c2,result_pow2,name,error)
     call error_check(name,error)
  end if

end subroutine min2

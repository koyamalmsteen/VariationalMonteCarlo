!
! $Id: pow5_zero.f90,v 1.1 2003/11/05 09:21:14 k-yukino Exp $
!
#include "parameter.h"

subroutine pow5_zero(result,name,error)
  use global_variables
  implicit none

  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
!
  result=0 ; name="pow5_zero" ; error=0


end subroutine pow5_zero

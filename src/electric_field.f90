!
! $Id: electric_field.f90,v 1.1 2003/06/04 10:12:21 k-yukino Exp $
!
#include "parameter.h"

subroutine set_electric_field(result,name,error)
  use global_variables
  implicit none

  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  result=0
  name="electric_field" ; error=0

  call set_static1_field(name,error)
  call set_static2_field(name,error)

end subroutine set_electric_field

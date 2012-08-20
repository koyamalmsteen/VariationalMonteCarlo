!
! $Id: charge_density.f90,v 1.1 2003/08/14 09:27:37 k-yukino Exp $
!
#include "parameter.h"

subroutine charge_density(aaa,result,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,intent(in) :: aaa
  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  result=0 ; name="charge_density" ; error=0
! init local 

  if( global_site_table_up(aaa)==1 ) then
     result=1
  end if

  if( global_site_table_down(aaa)==1 ) then
     result=result+1
  end if

  result=result

end subroutine charge_density

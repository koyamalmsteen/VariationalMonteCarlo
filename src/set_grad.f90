!
! $Id: set_grad.f90,v 1.1 2003/11/05 09:22:17 k-yukino Exp $
!
#include "parameter.h"

subroutine set_grad(xxx,yyy,grad,name,error)
  implicit none

  integer,intent(in) :: xxx,yyy
  real(8),intent(out) :: grad
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  grad=0 ; name="set_grad" ; error=0

  if( xxx==0 ) then
     error=-1 ; return
  end if

  grad=dble(yyy)/dble(xxx)

end subroutine set_grad

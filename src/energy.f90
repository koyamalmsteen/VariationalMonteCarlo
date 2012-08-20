!
! $Id: energy.f90,v 1.3 2003/11/20 14:01:58 k-yukino Exp $
!
#include "parameter.h"

subroutine energy(vector_up,vector_down,alpha_psi0,lambda,alpha_h_psi1,&
     name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  real(8),intent(in) :: alpha_psi0,lambda
  real(8),intent(out) :: alpha_h_psi1
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  alpha_h_psi1=0
  name="energy" ; error=0

!
! 期待値計算ルーチン
!
  call calc_element(alpha_psi0,lambda,vector_up,vector_down,&
       alpha_h_psi1,name,error)
  call error_check(name,error)

end subroutine energy



#ifdef DEBUG_ENERGY
program main
  impclicit none

end program main
#endif /* DEBUG_ENERGY */

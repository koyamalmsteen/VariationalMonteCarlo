!
! $Id$
!
#include "parameter.h"

subroutine calc_alpha_p_p_psi1(vector_up,vector_down,alpha_psi0,&
     alpha_p_p_psi1,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  real(8),intent(in) :: alpha_psi0
  real(8),intent(out) :: alpha_p_p_psi1
  character(32),intent(out) :: name
  integer,intent(out) :: error
  real(8),external :: calc_total_pole,jastrow
! init
  name="calc_alpha_p_p_psi1" ; error=0 
  alpha_p_p_psi1=0
  
  alpha_p_p_psi1=calc_total_pole(vector_up,vector_down)**2*alpha_psi0 &
       *jastrow(vector_up,vector_down)
!
! koko deha lambda wo tsukae nai. nazenara, f(0)de naikara
!

end subroutine calc_alpha_p_p_psi1

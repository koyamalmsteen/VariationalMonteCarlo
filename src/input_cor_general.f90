!
! $Id: input_cor_general.f90,v 1.1 2004/02/15 10:52:35 k-yukino Exp $
!
#include "parameter.h"

subroutine input_cor_general(alpha_psi1,tmp_sxsy0,tmp_super0,tmp_sz0,&
     name,error)
  use global_variables
  implicit none

  real(8),intent(in) :: alpha_psi1 
  real(8),dimension(max_number_xi),intent(in) :: tmp_sxsy0,tmp_super0,tmp_sz0
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  name="input_cor_general" ; error=0

  result_sxsy0=result_sxsy0+tmp_sxsy0/alpha_psi1 &
       /dble(MAX_MONTECARLO_SAMPLE)/dble(TOTAL_SITE_NUMBER)
  result_s0=result_s0+(tmp_sz0+tmp_sxsy0/alpha_psi1) &
       /dble(MAX_MONTECARLO_SAMPLE)/dble(TOTAL_SITE_NUMBER)
  result_super0=result_super0+tmp_super0/alpha_psi1 &
       /dble(MAX_MONTECARLO_SAMPLE)/dble(TOTAL_SITE_NUMBER)

end subroutine input_cor_general

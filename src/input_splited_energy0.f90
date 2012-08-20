!
! $Id: input_splited_energy0.f90,v 1.2 2004/03/10 03:04:02 k-yukino Exp ju448 $
!
#include "parameter.h"

subroutine input_splited_energy0(alpha_psi1,alpha_h_psi1,sample_number,&
     name,error)
  use global_variables
  implicit none

  real(8),intent(in) :: alpha_psi1
  real(8),intent(in) :: alpha_h_psi1
  integer,intent(in) :: sample_number
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  name="input_splited_energy0" ; error=0

  if( mod(sample_number,5)==1 ) then
     splited_energy(1,sample_number/5+1)=alpha_h_psi1/alpha_psi1&
          /dble(MAX_MONTECARLO_SAMPLE/5)
  else if( mod(sample_number,5)==2 ) then
     splited_energy(2,sample_number/5+1)=&
          alpha_h_psi1/alpha_psi1/dble(MAX_MONTECARLO_SAMPLE/5)
  else if( mod(sample_number,5)==3 ) then
     splited_energy(3,sample_number/5+1)=&
          alpha_h_psi1/alpha_psi1/dble(MAX_MONTECARLO_SAMPLE/5)
  else if( mod(sample_number,5)==4 ) then
     splited_energy(4,sample_number/5+1)=&
          alpha_h_psi1/alpha_psi1/dble(MAX_MONTECARLO_SAMPLE/5)
  else if( mod(sample_number,5)==0 ) then
     splited_energy(5,sample_number/5+1)=&
          alpha_h_psi1/alpha_psi1/dble(MAX_MONTECARLO_SAMPLE/5)
  end if

end subroutine input_splited_energy0

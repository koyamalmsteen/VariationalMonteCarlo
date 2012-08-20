!
! $Id: energy_zero_power.f90,v 1.2 2003/11/20 16:10:27 k-yukino Exp $
!
#include "parameter.h"

subroutine energy_zero_power(name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8) :: h0_zero,h1_zero,h2_zero,h3_zero,h4_zero,h5_zero
  character(32) :: name
  integer :: error
! init
  name="energy_zero_power" ; error=0
  h0_zero=1 ; h1_zero=0 ; h2_zero=0 ; h3_zero=0 ; h4_zero=0 ; h5_zero=0

  call pow1_zero(h1_zero,name,error)
  call error_check(name,error)
  call pow2_zero(h2_zero,name,error)
  call error_check(name,error)
  call pow3_zero(h3_zero,name,error)
  call error_check(name,error)
  call pow4_zero(h4_zero,name,error)       ! 実は作ってない
  call error_check(name,error)     
  call pow5_zero(h5_zero,name,error)       ! 実は作ってない
  call error_check(name,error)     

  h0_zero=1

end subroutine energy_zero_power

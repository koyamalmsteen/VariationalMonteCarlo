!
! $Id: input_cor_number.f90,v 1.1 2004/02/15 10:52:35 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine input_cor_number(tmp_sz0,tmp_charge0,tmp_charge_density0,name,error)
  use global_variables
  implicit none

  real(8),dimension(max_number_xi),intent(in) :: tmp_sz0,tmp_charge0
  real(8),dimension(TOTAL_SITE_NUMBER),intent(in) :: tmp_charge_density0
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
! init
  name="input_cor_number" ; error=0

  result_sz0=result_sz0+tmp_sz0/dble(MAX_MONTECARLO_SAMPLE) &
       /dble(TOTAL_SITE_NUMBER)
  result_charge0=result_charge0 &
       +tmp_charge0/dble(MAX_MONTECARLO_SAMPLE) &
       /dble(TOTAL_SITE_NUMBER)
  result_charge_density0=result_charge_density0 &
       +tmp_charge_density0 &
       /dble(MAX_MONTECARLO_SAMPLE)
  do i_count=1,TOTAL_SITE_NUMBER
     result_pole0=result_pole0 &
          +xaxis(i_count)*(1-tmp_charge_density0(i_count)) &
          /dble(MAX_MONTECARLO_SAMPLE)
  end do

end subroutine input_cor_number

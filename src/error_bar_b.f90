!
! $Id: error_bar_b.f90,v 1.1 2004/03/10 03:04:02 k-yukino Exp $
!
#include "parameter.h"

subroutine error_bar_b(splited_charge_density,lowest_charge_density,&
     highest_charge_density,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),dimension(5,max_number_xi),intent(in) :: splited_charge_density
  real(8),dimension(max_number_xi),intent(out) :: lowest_charge_density,&
       highest_charge_density
  character(32),intent(out) :: name
  integer :: error
! local
  integer :: i_count,j_count
! init
  lowest_charge_density=0 ; highest_charge_density=0
  name="error_bar_b" ; error=0
  
  ! 仮に代入
  do i_count=1,TOTAL_SITE_NUMBER
     lowest_charge_density(i_count)=splited_charge_density(1,i_count)
     highest_charge_density(i_count)=splited_charge_density(1,i_count)
  end do

  do i_count=2,5
     do j_count=1,TOTAL_SITE_NUMBER
        if( lowest_charge_density(j_count)>splited_charge_density(i_count,j_count) ) then
           lowest_charge_density(j_count)=splited_charge_density(i_count,j_count)
        end if
        if( highest_charge_density(j_count)<splited_charge_density(i_count,j_count) ) then
           highest_charge_density(j_count)=splited_charge_density(i_count,j_count)
        end if
     end do
  end do

end subroutine error_bar_b

!
! $Id: error_bar_c.f90,v 1.1 2004/03/10 03:04:02 k-yukino Exp $
!
#include "parameter.h"

subroutine error_bar_c(splited_pole,lowest_pole,highest_pole,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),dimension(5),intent(in) :: splited_pole
  real(8),intent(out) :: lowest_pole,highest_pole
  character(32),intent(out) :: name
  integer :: error
! local
  integer :: i_count,j_count
! init
  lowest_pole=0 ; highest_pole=0
  name="error_bar_c" ; error=0
  
  ! 仮に代入
  lowest_pole=splited_pole(1)
  highest_pole=splited_pole(1)

  do i_count=2,5
     if( lowest_pole>splited_pole(i_count) ) then
        lowest_pole=splited_pole(i_count)
     end if
     if( highest_pole<splited_pole(i_count) ) then
        highest_pole=splited_pole(i_count)
     end if
  end do

end subroutine error_bar_c

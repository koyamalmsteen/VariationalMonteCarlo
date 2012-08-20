!
! $Id: error_bar_a.f90,v 1.1 2004/03/10 03:04:02 k-yukino Exp $
!
#include "parameter.h"

subroutine error_bar_a(splited_cor,lowest_cor,highest_cor,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),dimension(5,max_number_xi),intent(in) :: splited_cor
  real(8),dimension(max_number_xi),intent(out) :: lowest_cor,highest_cor
  character(32),intent(out) :: name
  integer :: error
! local
  integer :: i_count,j_count
! init
  lowest_cor=0 ; highest_cor=0
  name="error_bar_a" ; error=0
  
  ! 仮に代入
  do i_count=1,max_number_xi
     lowest_cor(i_count)=splited_cor(1,i_count)
     highest_cor(i_count)=splited_cor(1,i_count)
  end do

  do i_count=2,5
     do j_count=1,max_number_xi
        if( lowest_cor(j_count)>splited_cor(i_count,j_count) ) then
           lowest_cor(j_count)=splited_cor(i_count,j_count)
        end if
        if( highest_cor(j_count)<splited_cor(i_count,j_count) ) then
           highest_cor(j_count)=splited_cor(i_count,j_count)
        end if
     end do
  end do

end subroutine error_bar_a

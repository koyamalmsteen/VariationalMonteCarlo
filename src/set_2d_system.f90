!
! $Id: set_2d_system.f90,v 1.2 2003/11/05 09:21:45 k-yukino Exp $
!
#include "parameter.h"
 
subroutine set_2d_system(delta_x,delta_y,tmptmp_site_table,name,error)
  use global_variables
  implicit none

  type coordinate
     integer :: xaxis
     integer :: yaxis
  end type coordinate

  integer,intent(out) :: delta_x,delta_y
  type(coordinate),dimension(TOTAL_SITE_NUMBER),intent(out) :: tmptmp_site_table 
  character(32),intent(out) :: name
  integer,intent(out) :: error
! lcoal
  integer :: i_count,j_count
  integer :: site_number
! init
  tmptmp_site_table%xaxis=0 ; tmptmp_site_table%yaxis=0
  delta_x=0 ; delta_y=0
  name="set_2d_system" ; error=0
! init local
  site_number=0

  do i_count=1,9999                           ! 傾き
     do j_count=0,i_count
        call count_site_number(i_count,j_count,tmptmp_site_table,site_number,&
                               name,error)
        call error_check(name,error)

        if( TOTAL_SITE_NUMBER==site_number ) then
           delta_x=i_count ; delta_y=j_count
           return
        end if
     end do
  end do

  error=-1                  ! なんらかのエラー

end subroutine set_2d_system

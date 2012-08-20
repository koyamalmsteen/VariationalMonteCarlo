!
! $Id: charge.f90,v 1.4 2003/08/14 09:26:54 k-yukino Exp $ 
!
#include "parameter.h"

subroutine charge(aaa,bbb,result,name,error)
  use global_variables
  implicit none

  include "mpif.h"

! 電子状態はグローバルサイトテーブルにセットしておく
  integer,intent(in) :: aaa,bbb
  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  name="charge" ; error=0
  result=0

  if( global_site_table_up(aaa)==0 .and. global_site_table_down(aaa)==0 ) then
     if( global_site_table_up(bbb)==0 &
          .and. global_site_table_down(bbb)==0 ) then
        result=1
     else if( global_site_table_up(bbb)==1 &
          .and. global_site_table_down(bbb)==1 ) then
        result=-1
     end if
  else if( global_site_table_up(aaa)==1 &
       .and. global_site_table_down(aaa)==1 ) then
     if( global_site_table_up(bbb)==0 &
          .and. global_site_table_down(bbb)==0 ) then
        result=-1
     else if( global_site_table_up(bbb)==1 &
          .and. global_site_table_down(bbb)==1 ) then
        result=1
     end if
  else
     result=0
  end if

end subroutine charge



#ifdef DEBUG_CHARGE
program main
  implicit none

  character(32) :: name
  integer :: error
! MPI
  integer :: IERROR
! init
  name="main" ; error=0
  
  call MPI_INIT(IERROR)

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CHARGE */

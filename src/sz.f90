!
! $Id: sz.f90,v 1.1 2003/08/14 09:34:47 k-yukino Exp $
!
#include "parameter.h"

subroutine sz(aaa,bbb,result,name,error)
  use global_variables
  implicit none

  include "mpif.h"

! 電子状態はグローバルサイトテーブルにセットしておく
  integer,intent(in) :: aaa,bbb
  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8) :: result_a,result_b
! init
  name="sz" ; error=0
  result=0
  result_a=0 ; result_b=0

!!! A
  if( global_site_table_up(aaa)==1 ) then
     result_a=result_a+1
  end if
  if( global_site_table_down(aaa)==1 ) then
     result_a=result_a-1
  end if

!!! B
  if( global_site_table_up(bbb)==1 ) then
     result_b=result_b+1
  end if
  if( global_site_table_down(bbb)==1 ) then
     result_b=result_b-1
  end if

  result=0.25*result_a*result_b

end subroutine sz



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

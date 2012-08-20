!
! $Id: calc_pole0.f90,v 1.1 2004/03/10 03:04:02 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine calc_pole0(aaa,tmp_result,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,intent(in) :: aaa
  real(8),intent(out) :: tmp_result
  character(32),intent(out) :: name
  integer,intent(out) :: error  
! local
  integer :: qn 
  integer :: max_x
! init
  name="calc_pole0" ; error=0
  tmp_result=0
  max_x=0

! up電子があれば
  if( global_site_table_up(aaa)==1 .and. global_site_table_down(aaa)==1 ) then
     qn=-1
  else if( global_site_table_up(aaa)==1 .and. global_site_table_down(aaa)==0 ) then
     qn=0
  else if( global_site_table_up(aaa)==0 .and. global_site_table_down(aaa)==1 ) then
     qn=0
  else
     qn=1
  end if

  tmp_result=xaxis(aaa)*qn

end subroutine calc_pole0



#ifdef DEBUG_JASTROW
program main
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  real(8),external :: jastrow
  character(32) :: name
  integer :: error
! local
  integer :: i_count
  integer,external :: count_number_xi
  real(8) :: tmp
! MPI関連
  integer :: IERROR
! init
  name="main" ; error=0
  i_count=0 ; tmp=0

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

!  do i_count=1,TOTAL_UP_ELECTRON
!     write(*,*) "gamma_up >" ; read(*,*) gamma_up(i_count)
!  end do
!
!  do i_count=1,TOTAL_DOWN_ELECTRON
!     write(*,*) "gamma_down >" ; read(*,*) gamma_down(i_count)
!  end do
  
  gamma_up(1)=1 ; gamma_up(2)=2 ; gamma_up(3)=3 
  gamma_up(4)=4 ; gamma_up(5)=5 ; gamma_up(6)=6
  gamma_up(7)=7 ; gamma_up(8)=8 ; gamma_up(9)=9
  gamma_up(10)=10 ; gamma_up(11)=11 ; gamma_up(12)=12
  gamma_up(13)=13 ; gamma_up(14)=14 ; gamma_up(15)=15

  gamma_down(1)=16 ; gamma_down(2)=17 ; gamma_down(3)=18
  gamma_down(4)=19 ; gamma_down(5)=20 ; gamma_down(6)=21
  gamma_down(7)=22 ; gamma_down(8)=23 ; gamma_down(9)=24
  gamma_down(10)=25 ; gamma_down(11)=26 ; gamma_down(12)=27
  gamma_down(13)=28 ; gamma_down(14)=29 ; gamma_down(15)=1

  tmp=jastrow(gamma_up,gamma_down)
  write(*,*) "jastrow=",tmp

  call MPI_FINALIZE(IERROR)

end program main
#endif

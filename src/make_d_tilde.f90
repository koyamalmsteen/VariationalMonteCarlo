!
! $Id: make_d_tilde.f90,v 1.13 2003/11/20 14:05:22 k-yukino Exp $  
!
! up,downスピン共用
!

#include "parameter.h"

subroutine make_d_tilde(gamma_up,gamma_down,spin,algorithm,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: gamma_down
  integer,intent(in) :: spin     ! spin=1 up,spin=2 down
  integer,intent(in) :: algorithm ! algorithm=1 for montecarlo , algorithm=2
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count,j_count
! init
  name="make_d_tilde" ; error=0

! argument check

!
! このルーチンの入力部分にgamma_upとgamma_downの二つを要求しているが、
! 実は片方しか使わない。(spin変数を見て片方だけ使う)
!

! main
  if( algorithm==1 ) then
     if( spin==1 ) then
        do i_count=1,TOTAL_UP_ELECTRON
           do j_count=1,TOTAL_UP_ELECTRON
              d_tilde_up(i_count,j_count)= &
                   unitary_d_up(j_count,gamma_up(i_count))
           end do
        end do
     else if( spin==2 ) then
        do i_count=1,TOTAL_DOWN_ELECTRON
           do j_count=1,TOTAL_DOWN_ELECTRON
              d_tilde_down(i_count,j_count)= &
                   unitary_d_down(j_count,gamma_down(i_count))
           end do
        end do
     else
        error=-11 ; return
     end if
  else if( algorithm==2 ) then                 ! algorithm==2なら
     if( spin==1 ) then
        do i_count=1,TOTAL_UP_ELECTRON
           do j_count=1,TOTAL_UP_ELECTRON
              d_tilde_up_work(i_count,j_count)= &
                   unitary_d_up(j_count,gamma_up(i_count))
           end do
        end do
     else if( spin==2 ) then 
        do i_count=1,TOTAL_DOWN_ELECTRON
           do j_count=1,TOTAL_DOWN_ELECTRON
              d_tilde_down_work(i_count,j_count)= &
                   unitary_d_down(j_count,gamma_down(i_count))
           end do
        end do
     else
        error=-22 ; return
     end if
  else
     error=-44 ; return
  end if

end subroutine make_d_tilde



#ifdef DEBUG_MAKE_D_TILDE
program main
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON):: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON):: gamma_down
  character(32) :: name
  integer :: error
! local
  integer :: ier,i_count,j_count
! init
  name="main" ; error=0 ; ier=0 ; i_count=0 ; j_count=0

  call MPI_INIT(IERROR)

  write(*,*) "TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "TOTAL_UP_ELECTRON=",TOTAL_UP_ELECTRON

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        unitary_d_up(i_count,j_count)=i_count*10 + j_count
        unitary_d_down(i_count,j_count)=i_count*10 + j_count
        write(*,*) i_count,j_count,"=",unitary_d_up(i_count,j_count)
        write(*,*) i_count,j_count,"=",unitary_d_down(i_count,j_count)
     end do
  end do

  do i_count=1,TOTAL_UP_ELECTRON
     write(*,*) "gamma_up ",i_count,">"
     read(*,*) gamma_up(i_count)
  end do

  call make_d_tilde(gamma_up,gamma_down,1,1,name,error)
  call error_check(name,error)

  do i_count=1,TOTAL_UP_ELECTRON
     do j_count=1,TOTAL_UP_ELECTRON
        write(*,*) "d_tilde_up(",i_count,j_count,")=",&
                                  d_tilde_up(i_count,j_count)
     end do
  end do

  call MPI_FINALIZE(IERROR)

end program main

#endif /* DEBUG_MAKE_D_TILDE */ 

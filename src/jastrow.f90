!
! $Id: jastrow.f90,v 1.11 2003/11/05 15:31:58 k-yukino Exp k-yukino $
!
#include "parameter.h"

real(8) function jastrow(gamma_up,gamma_down)
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: gamma_down
! local
  integer :: i_count,j_count
  real(8) :: tmp_result
! init
  tmp_result=0 ; i_count=0 ; j_count=0

  if( NUMBER_XI/=0 ) then
! up$B$H(Bdown
     do i_count=1,TOTAL_UP_ELECTRON
        do j_count=1,TOTAL_DOWN_ELECTRON
           tmp_result=tmp_result &
                +xi_table(gamma_up(i_count),gamma_down(j_count))%xi
        end do
     end do
  end if



! up$B$H(Bup
  if( NUMBER_XI>1 ) then            ! $BJQJ,%Q%i%a!<%?$,(B2$B0J>e$N$H$-$O9T$J$&(B

     do i_count=1,TOTAL_UP_ELECTRON-1
        do j_count=i_count+1,TOTAL_UP_ELECTRON
           tmp_result=tmp_result &
                +xi_table(gamma_up(i_count),gamma_up(j_count))%xi
        end do
     end do

! down$B$H(Bdown
     do i_count=1,TOTAL_DOWN_ELECTRON-1
        do j_count=i_count+1,TOTAL_DOWN_ELECTRON
           tmp_result=tmp_result &
                +xi_table(gamma_down(i_count),gamma_down(j_count))%xi
        end do
     end do

  end if

   jastrow=exp(-0.5*tmp_result)

end function jastrow



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
! MPI$B4XO"(B
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

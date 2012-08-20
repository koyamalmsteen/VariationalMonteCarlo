!
! $Id: calc_new_rho.f90,v 1.5 2002/12/20 05:42:59 k-yukino Exp $
!
#include "parameter.h"

subroutine calc_new_rho(rho_up_new,rho_down_new,name,error)
  use global_variables
  implicit none

  real(8),dimension(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),intent(out) :: rho_up_new,rho_down_new
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count,j_count,k_count
! init 
  i_count=0 ; j_count=0 ; rho_up_new=0 ; rho_down_new=0 
  name="calc_rho_sigma" ; error=0

! argument check
  
!
! 今扱っているHFO係数の要素は実数であるので、D^*はDそのものである。
! よって以下でよい。
!

  rho_up_new=0 ; rho_down_new=0

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        do k_count=1,TOTAL_UP_ELECTRON
           rho_up_new(i_count,j_count)=rho_up_new(i_count,j_count) & 
                                       + unitary_d_up(i_count,k_count) &
                                       * unitary_d_up(j_count,k_count)
        end do

        do k_count=1,TOTAL_DOWN_ELECTRON
           rho_down_new(i_count,j_count)=rho_down_new(i_count,j_count) & 
                                         + unitary_d_down(i_count,k_count) &
                                         * unitary_d_down(j_count,k_count)
        end do
     end do
  end do

end subroutine calc_new_rho



#ifdef DEBUG_CALC_NEW_RHO
program main
  use global_variables
  implicit none

  real(8),dimension(TOTAL_SITE_NUMBER) :: rho_up_new,rho_down_new
  character(32) :: name
  integer :: error
! local
  integer :: ier,i_count,j_count
  real(8) :: ransuu
! MPI関連
  integer IERROR
! init
  rho_up_new=0
  name="main" ; error=0 ; ier=0 ; ransuu=0 ; i_count=0 ; j_count=0
  IERROR=0

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

  call random_init(name,error)
  call error_check(name,error)

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        unitary_d_up(i_count,j_count)=i_count*10+j_count
        write(*,*) unitary_d_up(i_count,j_count)
     end do
     do j_count=1,TOTAL_SITE_NUMBER
        unitary_d_up(i_count,j_count)=i_count*10+j_count
!        write(*,*) unitary_d_up(i_count,j_count)
     end do
  end do

  call calc_new_rho(rho_up_new,rho_down_new,name,error)
  call error_check(name,error)

  do i_count=1,TOTAL_SITE_NUMBER
     write(*,*) "rho_up_new(",i_count,")=",rho_up_new(i_count)
  end do

  write(*,*) "TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "TOTAL_UP_ELECTRON=",TOTAL_UP_ELECTRON

  call MPI_FINALIZE(IERROR) 
end program main
#endif /* DEBUG_CALC_NEW_RHO */

!
! $Id: sxsy_zero.f90,v 1.4 2003/11/20 14:49:12 k-yukino Exp $
!
#include "parameter.h"

subroutine sxsy_zero(aaa,bbb,result,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,intent(in) :: aaa,bbb
  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: k_count
  real(8) :: temp_a,temp_b,temp_c,temp_d
! init
  result=0 ; name="sxsy_zero" ; error=0
! init local
  temp_a=0 ; temp_b=0 ; temp_c=0 ; temp_d=0

! numerator

  do k_count=1,TOTAL_UP_ELECTRON
     temp_a=temp_a+unitary_d_up(k_count,aaa)*unitary_d_up(k_count,bbb)
  end do

  do k_count=1,TOTAL_DOWN_ELECTRON
     temp_b=temp_b+unitary_d_down(k_count,bbb)*unitary_d_down(k_count,aaa)
  end do

  do k_count=1,TOTAL_DOWN_ELECTRON
     temp_c=temp_c+unitary_d_down(k_count,aaa)*unitary_d_down(k_count,bbb)
  end do

  do k_count=1,TOTAL_UP_ELECTRON
     temp_d=temp_d+unitary_d_up(k_count,bbb)*unitary_d_up(k_count,aaa)
  end do

  result=0.5*(temp_a*(1-temp_b)+temp_c*(1-temp_d))

end subroutine sxsy_zero



#ifdef DEBUG_SXSY_ZERO
program main
  use global_variables
  implicit none

  character(32) :: name
  integer :: error
! local
  real(8) :: tmp
  integer,dimension(2) :: ic
  real(8) :: zero_approx_energy
  integer :: hf_iteration
  integer :: i_count,j_count
! name="main" ; error=0
  ic=1 ; zero_approx_energy=0
  tmp=0

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

! パラメーターの初期化
  call random_init(ic,name,error)
  call error_check(name,error)

! 平面波解
  call zero_approx(zero_approx_energy,hf_iteration,name,error)
  call error_check(name,error)

  write(*,*) "JIGEN=",JIGEN,"TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "INIT_WAVE_VECTOR=",INIT_WAVE_VECTOR
  write(*,*)

  result_sxsy_zero=0

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=i_count,TOTAL_SITE_NUMBER
        call sxsy_zero(i_count,j_count,tmp,name,error)

        result_sxsy_zero(distance_sequence(i_count,j_count))= &
             result_sxsy_zero(distance_sequence(i_count,j_count))+tmp
     end do
  end do

  write(*,*) "result_sxsy_zero=",result_sxsy_zero/dble(TOTAL_SITE_NUMBER)

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_SXSY_ZERO */

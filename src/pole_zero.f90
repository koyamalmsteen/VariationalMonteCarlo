!
! $Id: pole_zero.f90,v 1.1 2004/03/10 03:04:02 k-yukino Exp $
!
#include "parameter.h"

subroutine pole_zero(aaa,result,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,intent(in) :: aaa
  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local 
  real(8),dimension(2) :: temp
  integer :: k_count
! init
  result=0 ; name="pole_zero" ; error=0
  temp=0

  do k_count=1,TOTAL_UP_ELECTRON
     temp(1)=temp(1)+unitary_d_up(k_count,aaa)*unitary_d_up(k_count,aaa)
  end do

  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(2)=temp(2)+unitary_d_down(k_count,aaa)*unitary_d_down(k_count,aaa)
  end do

  result=xaxis(aaa)*(1-temp(1)-temp(2))

end subroutine pole_zero


#ifdef DEBUG_SZ_ZERO
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
  
  result_sz_zero=0

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=i_count,TOTAL_SITE_NUMBER
        call sz_zero(i_count,j_count,tmp,name,error)

        result_sz_zero(distance_sequence(i_count,j_count))= &
             result_sz_zero(distance_sequence(i_count,j_count))+tmp
     end do
  end do

  write(*,*) "result_sz_zero=",result_sz_zero/dble(TOTAL_SITE_NUMBER)

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_SZ_ZERO */

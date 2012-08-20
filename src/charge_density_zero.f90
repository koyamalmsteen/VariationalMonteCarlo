!
! $Id: charge_density_zero.f90,v 1.2 2003/11/20 14:56:35 k-yukino Exp $
!
#include "parameter.h"

subroutine charge_density_zero(aaa,result,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,intent(in) :: aaa
  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: k_count
  real(8) :: temp1,temp2
! init
  result=0 ; name="charge_density_zero" ; error=0
! init local 
  temp1=0 ; temp2=0

  do k_count=1,TOTAL_UP_ELECTRON
     temp1=temp1+unitary_d_up(k_count,aaa)*unitary_d_up(k_count,aaa)
  end do
     
  do k_count=1,TOTAL_DOWN_ELECTRON
     temp2=temp2+unitary_d_down(k_count,aaa)*unitary_d_down(k_count,aaa)
  end do

  result=temp1+temp2

end subroutine charge_density_zero



#ifdef DEBUG_CHARGE_DENSITY_ZERO
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
  integer :: i_count
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
  
  result_charge_density_zero=0

  do i_count=1,TOTAL_SITE_NUMBER
     call charge_density_zero(i_count,tmp,name,error)

     result_charge_density_zero(i_count)= &
          result_charge_density_zero(i_count)+tmp
  end do

  write(*,*) "result_charge_density_zero=",result_charge_density_zero

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CHARGE_DENSITY_ZERO */


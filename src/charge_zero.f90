!
! $Id: charge_zero.f90,v 1.5 2003/11/20 14:47:23 k-yukino Exp $
!
#include "parameter.h"

subroutine charge_zero(aaa,bbb,result,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,intent(in) :: aaa,bbb
  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: k_count
  real(8),dimension(16) :: temp
! init
  result=0 ; name="charge_zero" ; error=0
! init local 
  temp=0

! temp(1)
  do k_count=1,TOTAL_UP_ELECTRON
     temp(1)=temp(1)+unitary_d_up(k_count,aaa)*unitary_d_up(k_count,aaa)
  end do

! temp(2)
  do k_count=1,TOTAL_UP_ELECTRON
     temp(2)=temp(2)+unitary_d_up(k_count,bbb)*unitary_d_up(k_count,bbb)
  end do

! temp(3)
  do k_count=1,TOTAL_UP_ELECTRON
     temp(3)=temp(3)+unitary_d_up(k_count,aaa)*unitary_d_up(k_count,bbb)
  end do

! temp(4)
  do k_count=1,TOTAL_UP_ELECTRON
     temp(4)=temp(4)+unitary_d_up(k_count,bbb)*unitary_d_up(k_count,aaa)
  end do

!
!

! temp(5)
  do k_count=1,TOTAL_UP_ELECTRON
     temp(5)=temp(5)+unitary_d_up(k_count,aaa)*unitary_d_up(k_count,aaa)
  end do

! temp(6)
  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(6)=temp(6)+unitary_d_down(k_count,bbb)*unitary_d_down(k_count,bbb)
  end do

!
!

! temp(7)
  do k_count=1,TOTAL_UP_ELECTRON
     temp(7)=temp(7)+unitary_d_up(k_count,aaa)*unitary_d_up(k_count,aaa)
  end do
  
!
!

! temp(8)
  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(8)=temp(8)+unitary_d_down(k_count,aaa)*unitary_d_down(k_count,aaa)
  end do

! temp(9)
  do k_count=1,TOTAL_UP_ELECTRON
     temp(9)=temp(9)+unitary_d_up(k_count,bbb)*unitary_d_up(k_count,bbb)
  end do

!
!

! temp(10)
  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(10)=temp(10)+unitary_d_down(k_count,aaa)*unitary_d_down(k_count,aaa)
  end do

! temp(11)
  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(11)=temp(11)+unitary_d_down(k_count,bbb)*unitary_d_down(k_count,bbb)
  end do

! temp(12)
  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(12)=temp(12)+unitary_d_down(k_count,aaa)*unitary_d_down(k_count,bbb)
  end do

! temp(13)
  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(13)=temp(13)+unitary_d_down(k_count,bbb)*unitary_d_down(k_count,aaa)
  end do

!
!

! temp(14)
  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(14)=temp(14)+unitary_d_down(k_count,aaa)*unitary_d_down(k_count,aaa)
  end do

!
!

! temp(15)
  do k_count=1,TOTAL_UP_ELECTRON
     temp(15)=temp(15)+unitary_d_up(k_count,bbb)*unitary_d_up(k_count,bbb)
  end do

!
!

! temp(16)
  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(16)=temp(16)+unitary_d_down(k_count,bbb)*unitary_d_down(k_count,bbb)
  end do

  result=temp(1)*temp(2)+temp(3)*(1-temp(4))&
       +temp(5)*temp(6)&
       -temp(7) &
       +temp(8)*temp(9) &
       +temp(10)*temp(11)+temp(12)*(1-temp(13))&
       -temp(14)-temp(15)-temp(16)+1

end subroutine charge_zero



#ifdef DEBUG_CHARGE_ZERO
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

  result_charge_zero=0

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=i_count,TOTAL_SITE_NUMBER
        call charge_zero(i_count,j_count,tmp,name,error)

        result_charge_zero(distance_sequence(i_count,j_count))= &
             result_charge_zero(distance_sequence(i_count,j_count))+tmp
     end do
  end do

  write(*,*) "result_charge_zero=",result_charge_zero/dble(TOTAL_SITE_NUMBER)

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CHARGE_ZERO*/

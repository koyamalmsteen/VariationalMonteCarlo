!
! $Id: pow1_zero.f90,v 1.1 2003/08/14 09:24:18 k-yukino Exp $
!
#include "parameter.h"

subroutine pow1_zero(result,name,error)
  use global_variables
  implicit none

  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i1_count
  integer :: j1_count
  integer :: k_count
  real(8) :: t1
  real(8),dimension(2) :: temp
! init
  result=0
  name="pow1_zero"
  error=0
!
  temp=0
  

! 第1項
  do i1_count=1,TOTAL_SITE_NUMBER
     do j1_count=1,TOTAL_SITE_NUMBER
        if( neighbor_table(i1_count,j1_count)==1 ) then

! 第1項↑
           t1=0

           do k_count=1,TOTAL_UP_ELECTRON
              t1=t1+unitary_d_up(k_count,i1_count)&
                   *unitary_d_up(k_count,j1_count)
           end do

           temp(1)=t1

! 第1項↓
           t1=0

           do k_count=1,TOTAL_DOWN_ELECTRON
              t1=t1+unitary_d_down(k_count,i1_count)&
                   *unitary_d_down(k_count,j1_count)
           end do

           temp(2)=t1

           result=result+dble(TRANSFER)*( temp(1)+temp(2) )
           
        end if
     end do
  end do

end subroutine pow1_zero

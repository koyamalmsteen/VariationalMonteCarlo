!
! $Id: super_zero.f90,v 1.5 2003/11/20 14:51:21 k-yukino Exp $
!

!
! k_count=1,右方向
! k_count=2,下方向
! k_count=3,左方向
! k_count=4,上方向
!

#include "parameter.h"

subroutine super_zero(aaa,bbb,result,name,error)
  use global_variables
  implicit none

  include "mpif.h"
  
  integer,intent(in) :: aaa,bbb
  real(8),intent(out) :: result 
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: k_count,l_count,n_count
  integer :: fugou
  real(8) :: temp11,temp12,temp21,temp22,temp31,temp32,temp41,temp42
! init
  result=0 ; name="super_zero" ; error=0
  fugou=0
!
  temp11=0 ; temp12=0 ; temp21=0 ; temp22=0
  temp31=0 ; temp32=0 ; temp41=0 ; temp42=0

  do k_count=1,4
     do l_count=1,4
        temp11=0 ; temp12=0 ; temp21=0 ; temp22=0
        temp31=0 ; temp32=0 ; temp41=0 ; temp42=0

        fugou=(-1)**(k_count+l_count)

        if( neighbor_table2(aaa,k_count)/=0 .and. &
             neighbor_table2(bbb,l_count)/=0 ) then

! 第1項
           do n_count=1,TOTAL_UP_ELECTRON
              temp11=temp11+unitary_d_up(n_count,aaa)*unitary_d_up(n_count,bbb)
           end do

           do n_count=1,TOTAL_DOWN_ELECTRON
              temp12=temp12 &
                   +unitary_d_down(n_count,neighbor_table2(aaa,k_count))&
                   *unitary_d_down(n_count,neighbor_table2(bbb,l_count))
           end do

! 第2項
           do n_count=1,TOTAL_UP_ELECTRON
              temp21=temp21+unitary_d_up(n_count,aaa)&
                   *unitary_d_up(n_count,neighbor_table2(bbb,l_count))
           end do

           do n_count=1,TOTAL_DOWN_ELECTRON
              temp22=temp22&
                   +unitary_d_down(n_count,neighbor_table2(aaa,k_count))&
                   *unitary_d_down(n_count,bbb)
           end do

! 第3項
           do n_count=1,TOTAL_DOWN_ELECTRON
              temp31=temp31+unitary_d_down(n_count,aaa)&
                   *unitary_d_down(n_count,bbb)
           end do

           do n_count=1,TOTAL_UP_ELECTRON
              temp32=temp32+unitary_d_up(n_count,neighbor_table2(aaa,k_count))&
                   *unitary_d_up(n_count,neighbor_table2(bbb,l_count))
           end do

! 第4項
           do n_count=1,TOTAL_DOWN_ELECTRON
              temp41=temp41+unitary_d_down(n_count,aaa)&
                   *unitary_d_down(n_count,neighbor_table2(bbb,l_count))
           end do

           do n_count=1,TOTAL_UP_ELECTRON
              temp42=temp42+unitary_d_up(n_count,neighbor_table2(aaa,k_count))&
                   *unitary_d_up(n_count,bbb)
           end do

           result=result+fugou*(temp11*temp12+temp21*temp22+temp31*temp32&
                +temp41*temp42)*0.25*0.25

        end if

     end do
  end do

end subroutine super_zero



#ifdef DEBUG_SUPER_ZERO
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
  
  result_super_zero=0

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=i_count,TOTAL_SITE_NUMBER
        call super_zero(i_count,j_count,tmp,name,error)

        result_super_zero(distance_sequence(i_count,j_count))= &
             result_super_zero(distance_sequence(i_count,j_count))+tmp
     end do
  end do

  write(*,*) "result_super_zero=",result_super_zero/dble(TOTAL_SITE_NUMBER)

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_SUPER_ZERO */

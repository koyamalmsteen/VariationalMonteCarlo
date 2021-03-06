!
! $Id: sz_zero.f90,v 1.3 2003/11/20 14:45:31 k-yukino Exp $
!
#include "parameter.h"

subroutine sz_zero(aaa,bbb,result,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,intent(in) :: aaa,bbb
  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local 
  integer :: k_count
  real(8),dimension(12) :: temp
! init
  result=0 ; name="sz_zero" ; error=0
! init local
  temp=0

  do k_count=1,TOTAL_UP_ELECTRON
     temp(1)=temp(1)+unitary_d_up(k_count,aaa)*unitary_d_up(k_count,aaa)
  end do

  do k_count=1,TOTAL_UP_ELECTRON
     temp(2)=temp(2)+unitary_d_up(k_count,bbb)*unitary_d_up(k_count,bbb)
  end do

  do k_count=1,TOTAL_UP_ELECTRON
     temp(3)=temp(3)+unitary_d_up(k_count,aaa)*unitary_d_up(k_count,bbb)
  end do

  do k_count=1,TOTAL_UP_ELECTRON
     temp(4)=temp(4)+unitary_d_up(k_count,bbb)*unitary_d_up(k_count,aaa)
  end do

  do k_count=1,TOTAL_UP_ELECTRON
     temp(5)=temp(5)+unitary_d_up(k_count,aaa)*unitary_d_up(k_count,aaa)
  end do

  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(6)=temp(6)+unitary_d_down(k_count,bbb)*unitary_d_down(k_count,bbb)
  end do

  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(7)=temp(7)+unitary_d_down(k_count,aaa)*unitary_d_down(k_count,aaa)
  end do

  do k_count=1,TOTAL_UP_ELECTRON
     temp(8)=temp(8)+unitary_d_up(k_count,bbb)*unitary_d_up(k_count,bbb)
  end do

  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(9)=temp(9)+unitary_d_down(k_count,aaa)*unitary_d_down(k_count,aaa)
  end do

  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(10)=temp(10)+unitary_d_down(k_count,bbb)*unitary_d_down(k_count,bbb)
  end do

  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(11)=temp(11)+unitary_d_down(k_count,aaa)*unitary_d_down(k_count,bbb)
  end do

  do k_count=1,TOTAL_DOWN_ELECTRON
     temp(12)=temp(12)+unitary_d_down(k_count,bbb)*unitary_d_down(k_count,aaa)
  end do

  result=0.25*(temp(1)*temp(2)+temp(3)*(1-temp(4))&
       -temp(5)*temp(6)-temp(7)*temp(8)&
       +temp(9)*temp(10)+temp(11)*(1-temp(12)) )


end subroutine sz_zero


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

! $B%Q%i%a!<%?!<$N=i4|2=(B
  call random_init(ic,name,error)
  call error_check(name,error)

! $BJ?LLGH2r(B
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

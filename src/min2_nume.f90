!
! $Id: min2_nume.f90,v 1.1 2004/02/15 10:52:35 k-yukino Exp k-yukino $
!

#include "parameter.h"

subroutine min2_nume(h0,h1,h2,h3,h4,h5,pow2_c0,pow2_c1,pow2_c2,&
     result_pow2,name,error)
  use global_variables
  implicit none

  real(8),intent(in) :: h0,h1,h2,h3,h4,h5
  real(8),intent(out) :: pow2_c0,pow2_c1,pow2_c2
  real(8),intent(out) :: result_pow2
  character(32),intent(out) :: name
  integer,intent(out) :: error
!
  real(8),external :: evaluate_pow2
  integer :: i_count,j_count,k_count
  real(8) :: temp_result,min_pow2
! init
  pow2_c0=0 ; pow2_c1=0 ; pow2_c2=0
  name="min2_nume" ; error=0
  result_pow2=0

  min_pow2=99999

  if( POWER==2 ) then
     temp_result=0

     do i_count=-100,100           ! C0
        do j_count=-100,100        ! C1
           do k_count=-100,100     ! C2 
              temp_result=evaluate_pow2(h0,h1,h2,h3,&
                   h4,h5,dble(i_count*0.01),dble(j_count*0.01),&
                   dble(k_count*0.01))

              if( min_pow2>temp_result ) then
                 min_pow2=temp_result
                 pow2_c0=dble(i_count*0.01)
                 pow2_c1=dble(j_count*0.01)
                 pow2_c2=dble(k_count*0.01)
              end if
           end do
        end do
     end do
  end if

  result_pow2=min_pow2

end subroutine min2_nume


real(8) function evaluate_pow2(h0,h1,h2,h3,h4,h5,C0,C1,C2)
  implicit none

  real(8),intent(in) :: h0,h1,h2,h3,h4,h5
  real(8),intent(in) :: C0,C1,C2

  if( C0/=0 .and. C1/=0 .and. C2/=0 ) then
     evaluate_pow2=( (C2**2)*h5 &
          +2*C1*C2*h4 &
          +2*C0*C2*h3 &
          +(C1**2)*h3 &
          +2*C0*C1*h2 &
          +(C0**2)*h1 ) &
          /( (C2**2)*h4 &
          +2*C1*C2*h3 &
          +2*C0*C2*h2 & 
          +(C1**2)*h2 &
          +2*C0*C1*h1 &
          +(C0**2)*h0 )
  else
     evaluate_pow2=9999
  end if

end function evaluate_pow2

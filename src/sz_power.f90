!
! $Id: sz_power.f90,v 1.4 2004/02/15 10:52:35 k-yukino Exp $
!

#include "parameter.h"

subroutine sz_power(alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,tmp_sz0,&
     nume1,nume2,nume3_1,nume3_2,nume4,nume5,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(in) :: alpha_psi1,alpha_h_psi1,alpha_h_h_psi1
  real(8),dimension(max_number_xi),intent(in) :: tmp_sz0
  real(8),dimension(max_number_xi),intent(out) :: nume1,nume2,nume3_1,&
       nume3_2,nume4,nume5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
! init
  nume1=0 ; nume2=0 ; nume3_1=0 ; nume3_2=0 ; nume4=0 ; nume5=0
  name="sz_power" ; error=0
!
! メインルーチン
!

  do i_count=1,max_number_xi
!! <Sz> tmp_sz0で計算している
     nume1(i_count)=tmp_sz0(i_count)*(alpha_psi1**2)

!! <H Sz>            (<Sz>の結果を使う)
     nume2(i_count)=tmp_sz0(i_count)*(alpha_h_psi1*alpha_psi1)

!! <H Sz H>           (<Sz>の結果を使う) 
     nume3_1(i_count)=tmp_sz0(i_count)*(alpha_h_psi1**2)
!
!
! 
     if( POWER==2 ) then
!! <O H^2>   ---> <psi|O|alpha><alpha|H^2|psi>とする (<Sz>の結果を使う)
        nume3_2(i_count)=tmp_sz0(i_count)*(alpha_psi1*alpha_h_h_psi1)

!! <H O H^2> ---> <psi|HO|alpha><alpha|H^2|psi> (<Sz>の結果を使う)
        nume4(i_count)=tmp_sz0(i_count)*(alpha_h_psi1*alpha_h_h_psi1)

!! <H^2 O H^2> ---> <psi|H^2|alpha><alpha|OH^2|psi>
        nume5(i_count)=tmp_sz0(i_count)*(alpha_h_h_psi1**2)

     end if
  end do

end subroutine sz_power

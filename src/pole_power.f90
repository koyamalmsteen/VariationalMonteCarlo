!
! $Id: pole_power.f90,v 1.2 2004/03/10 03:04:02 k-yukino Exp k-yukino $
!

#include "parameter.h"

subroutine pole_power(alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,&
     tmp_pole0,nume1,nume2,nume3_1,nume3_2,nume4,nume5,name,error)
  use global_variables
  implicit none

  include "mpif.h"
  
  real(8),intent(in) :: alpha_psi1,alpha_h_psi1,alpha_h_h_psi1
  real(8),dimension(TOTAL_SITE_NUMBER),intent(in) :: tmp_pole0
  real(8),dimension(TOTAL_SITE_NUMBER),intent(out) :: nume1,nume2,nume3_1,&
       nume3_2,nume4,nume5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
! init
  nume1=0 ; nume2=0 ; nume3_1=0 ; nume3_2=0 ; nume4=0 ; nume5=0
  name="pole_power" ; error=0

!
! メインルーチン
!
  do i_count=1,TOTAL_SITE_NUMBER
!! <ChargeDensity> tmp_pole0で計算している 

! 注意(ややこしい)
     nume1(i_count)=tmp_pole0(i_count)*(alpha_psi1**2)

!! <H ChargeDensity>            (<ChargeDensity>の結果を使う)
     nume2(i_count)=tmp_pole0(i_count)*(alpha_h_psi1*alpha_psi1)

!! <H ChargeDensity H>           (<ChargeDensity>の結果を使う) 
     nume3_1(i_count)=tmp_pole0(i_count)*(alpha_h_psi1*alpha_h_psi1)

     if( POWER==2 ) then
!! <O H^2>   ---> <psi|O|alpha><alpha|H^2|psi>とする (<ChargeDensity>の結果を使う)
        nume3_2(i_count)=tmp_pole0(i_count)*alpha_psi1*alpha_h_h_psi1
          
!! <H O H^2> ---> <psi|HO|alpha><alpha|H^2|psi> (<ChargeDensity>の結果を使う)
        nume4(i_count)=tmp_pole0(i_count)*alpha_h_psi1*alpha_h_h_psi1

!! <H^2 O H^2> ---> <psi|H^2|alpha><alpha|OH^2|psi>
        nume5(i_count)=tmp_pole0(i_count)*alpha_h_h_psi1**2
     end if
  end do

end subroutine pole_power

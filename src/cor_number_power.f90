!
! $Id: cor_number_power.f90,v 1.4 2004/02/15 10:52:35 k-yukino Exp k-yukino $
!

#include "parameter.h"

subroutine cor_number_power(alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,&
     tmp_sz0,tmp_charge0,tmp_charge_density0,&
     sz_nume1,sz_nume2,sz_nume3_1,sz_nume3_2,sz_nume4,sz_nume5,&
     charge_nume1,charge_nume2,charge_nume3_1,&
     charge_nume3_2,charge_nume4,charge_nume5,&
     charge_density_nume1,charge_density_nume2,&
     charge_density_nume3_1,charge_density_nume3_2,&
     charge_density_nume4,charge_density_nume5,&
     name,error)
  use global_variables
  implicit none

  include "mpif.h"
  
  real(8),intent(in) :: alpha_psi1,alpha_h_psi1,alpha_h_h_psi1
  real(8),dimension(max_number_xi),intent(in) :: tmp_sz0
  real(8),dimension(max_number_xi),intent(in) :: tmp_charge0
  real(8),dimension(TOTAL_SITE_NUMBER),intent(in) :: tmp_charge_density0
  real(8),dimension(max_number_xi),intent(out) :: sz_nume1,sz_nume2,sz_nume3_1,sz_nume3_2,sz_nume4,sz_nume5
  real(8),intent(out) :: charge_nume1,charge_nume2,charge_nume3_1,charge_nume3_2,charge_nume4,charge_nume5
  real(8),intent(out) :: charge_density_nume1,charge_density_nume2,charge_density_nume3_1,charge_density_nume3_2,charge_density_nume4,charge_density_nume5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  sz_nume1=0 ; sz_nume2=0 ; sz_nume3_1=0 ; sz_nume3_2=0 
  sz_nume4=0 ; sz_nume5=0
! 
  charge_nume1=0 ; charge_nume2=0 ; charge_nume3_1=0 
  charge_nume3_2=0 ; charge_nume4=0 ; charge_nume5=0
!
  charge_density_nume1=0 ; charge_density_nume2=0 ; charge_density_nume3_1=0
  charge_density_nume3_2=0 ; charge_density_nume4=0 ; charge_density_nume5=0
!
  name="cor_number_power" ; error=0
! local init

!
! sz (szルーチンは状態ベクトルを送る必要がないようにできている)
!    (但し、要サンプル番号)
!
  call sz_power(alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,tmp_sz0,&
       sz_nume1,sz_nume2,sz_nume3_1,sz_nume3_2,sz_nume4,sz_nume5,name,error)
  call error_check(name,error)

!
! charge (szルーチンは状態ベクトルを送る必要がないようにできている)
!    (但し、要サンプル番号)
!

  call charge_power(alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,tmp_charge0,&
       charge_nume1,charge_nume2,charge_nume3_1,charge_nume3_2,charge_nume4,&
       charge_nume5,name,error)
  call error_check(name,error)

!
! charge_density (charge_densityルーチンは状態ベクトルを送る必要がないようにできている)
!    (但し、要サンプル番号)
!
  call charge_density_power(alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,&
       tmp_charge_density0,charge_density_nume1,charge_density_nume2,&
       charge_density_nume3_1,charge_density_nume3_2,charge_density_nume4,&
       charge_density_nume5,name,error)
  call error_check(name,error)

end subroutine cor_number_power

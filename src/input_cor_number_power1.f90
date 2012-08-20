!
! $Id: input_cor_number_power1.f90,v 1.2 2004/03/10 03:04:02 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine input_cor_number_power1(pow1_c0,pow1_c1,alpha_psi1,alpha_h_psi1,&
     sz_nume1,sz_nume2,sz_nume3_1,&
     sz_nume3_2,sz_nume4,sz_nume5,&
     charge_nume1,charge_nume2,charge_nume3_1,&
     charge_nume3_2,charge_nume4,charge_nume5,&
     charge_density_nume1,charge_density_nume2,&
     charge_density_nume3_1,charge_density_nume3_2,&
     charge_density_nume4,charge_density_nume5,&
     name,error)
  use global_variables
  implicit none

  real(8),intent(in) :: pow1_c0,pow1_c1
  real(8),intent(in) :: alpha_psi1,alpha_h_psi1
  real(8),dimension(max_number_xi),intent(in) :: sz_nume1,sz_nume2,sz_nume3_1,&
       sz_nume3_2,sz_nume4,sz_nume5
  real(8),dimension(max_number_xi),intent(in) :: charge_nume1,charge_nume2,&
       charge_nume3_1,charge_nume3_2,charge_nume4,charge_nume5
  real(8),dimension(TOTAL_SITE_NUMBER),intent(in) :: charge_density_nume1,&
       charge_density_nume2,charge_density_nume3_1,charge_density_nume3_2,&
       charge_density_nume4,charge_density_nume5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
! init
  name="input_cor_number_power1" ; error=0

! sz
  result_sz1=result_sz1+( sz_nume3_1*(pow1_c1**2)/(alpha_psi1**2) &
       +2*sz_nume2*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
       +sz_nume1*(pow1_c0**2)/(alpha_psi1**2) )
! charge
  result_charge1=result_charge1+( charge_nume3_1*(pow1_c1**2) &
       /(alpha_psi1**2) &
       +2*charge_nume2*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
       +charge_nume1*(pow1_c0**2)/(alpha_psi1**2) )
! charge_density
  result_charge_density1=result_charge_density1 &
       +( charge_density_nume3_1*(pow1_c1**2)/(alpha_psi1**2) &
       +2*charge_density_nume2*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
       +charge_density_nume1*(pow1_c0**2)/(alpha_psi1**2) )
! pole
  do i_count=1,TOTAL_SITE_NUMBER
     result_pole1=result_pole1 &
          +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count))&
          *(pow1_c1**2) /(alpha_psi1**2)  &
          +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count))&
          *(pow1_c0*pow1_c1)/(alpha_psi1**2) & 
          +xaxis(i_count)*(1*alpha_psi1**2-charge_density_nume1(i_count))&
          *(pow1_c0**2) /(alpha_psi1**2) 
  end do

end subroutine input_cor_number_power1

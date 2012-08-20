!
! $Id: input_cor_number_power2.f90,v 1.1 2004/02/15 10:52:35 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine input_cor_number_power2(pow2_c0,pow2_c1,pow2_c2,alpha_psi1,&
     alpha_h_psi1,alpha_h_h_psi1,&
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

  real(8),intent(in) :: pow2_c0,pow2_c1,pow2_c2
  real(8),intent(in) :: alpha_psi1,alpha_h_psi1,alpha_h_h_psi1
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
  name="input_cor_number_power2" ; error=0

! sz
  result_sz2=result_sz2+( sz_nume5*(pow2_c2**2)/(alpha_psi1**2) &
       +2*sz_nume4*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
       +2*sz_nume3_2*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
       +sz_nume3_1*(pow2_c1**2)/(alpha_psi1**2) &
       +2*sz_nume2*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
       +sz_nume1*(pow2_c0**2)/(alpha_psi1**2) )
! charge
  result_charge2=result_charge2+( charge_nume5*(pow2_c2**2) &
       /(alpha_psi1**2) &
       +2*charge_nume4*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
       +2*charge_nume3_2*(pow2_c0*pow2_c2)/(alpha_psi1**2 ) &
       +charge_nume3_1*(pow2_c1**2)/(alpha_psi1**2 ) &
       +2*charge_nume2*(pow2_c0*pow2_c1)/(alpha_psi1**2 ) &
       +charge_nume1*(pow2_c0**2)/(alpha_psi1**2) )
! charge_density
  result_charge_density2=result_charge_density2 &
       +( charge_density_nume5*(pow2_c2**2)/(alpha_psi1**2) &
       +2*charge_density_nume4*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
       +2*charge_density_nume3_2*(pow2_c0*pow2_c2)/(alpha_psi1**2)&
       +charge_density_nume3_1*(pow2_c1**2)/(alpha_psi1**2) &
       +2*charge_density_nume2*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
       +charge_density_nume1*(pow2_c0**2)/(alpha_psi1**2) )

! pole
  do i_count=1,TOTAL_SITE_NUMBER
     result_pole2=result_pole2 &
          +xaxis(i_count)*(1*alpha_h_h_psi1**2-charge_density_nume5(i_count)) &
          *(pow2_c2**2)/(alpha_psi1**2) &
          +2*xaxis(i_count)*(1*alpha_h_psi1*alpha_h_h_psi1-charge_density_nume4(i_count)) &
          *(pow2_c1*pow2_c2) /(alpha_psi1**2) &
          +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_h_psi1-charge_density_nume3_2(i_count)) &
          *(pow2_c0*pow2_c2) /(alpha_psi1**2)&
          +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count)) &
          *(pow2_c1**2) /(alpha_psi1**2)  &
          +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count)) &
          *(pow2_c0*pow2_c1) /(alpha_psi1**2) &
          +xaxis(i_count)*(1*alpha_psi1*alpha_psi1-charge_density_nume1(i_count)) &
          *(pow2_c0**2) /(alpha_psi1**2)
  end do

end subroutine input_cor_number_power2

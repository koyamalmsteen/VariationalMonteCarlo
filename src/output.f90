!
! $Id: output.f90,v 1.9 2004/03/10 03:04:02 k-yukino Exp $
!
#include "parameter.h"

subroutine output(zero_approx_energy,hf_iteration_result,result_pow0,&
     lowest_pow0,highest_pow0,result_pow1,lowest_pow1,highest_pow1,&
     result_pow2,lowest_pow2,highest_pow2,&
     lowest_sz0,highest_sz0,&
     lowest_sz1,highest_sz1,&
     lowest_sz2,highest_sz2,&
     lowest_charge0,highest_charge0,&
     lowest_charge1,highest_charge1,&
     lowest_charge2,highest_charge2,&
     lowest_sxsy0,highest_sxsy0,&
     lowest_sxsy1,highest_sxsy1,&
     lowest_sxsy2,highest_sxsy2,&
     lowest_s0,highest_s0,&
     lowest_s1,highest_s1,&
     lowest_s2,highest_s2,&
     lowest_super0,highest_super0,&
     lowest_super1,highest_super1,&
     lowest_super2,highest_super2,&
     lowest_charge_density0,highest_charge_density0,&
     lowest_charge_density1,highest_charge_density1,&
     lowest_charge_density2,highest_charge_density2,&
     lowest_pole0,highest_pole0,&
     lowest_pole1,highest_pole1,&
     lowest_pole2,highest_pole2,&
     pow1_c0,pow1_c1,pow2_c0,&
     pow2_c1,pow2_c2,name,error)
  use global_variables
  implicit none

  include "mpif.h"
  
  real(8),intent(in) :: zero_approx_energy
  integer,intent(in) :: hf_iteration_result
  real(8),intent(in) :: result_pow0,lowest_pow0,highest_pow0
  real(8),intent(in) :: result_pow1,lowest_pow1,highest_pow1
  real(8),intent(in) :: result_pow2,lowest_pow2,highest_pow2
  real(8),dimension(max_number_xi),intent(in) :: lowest_sz0,lowest_sz1,lowest_sz2
  real(8),dimension(max_number_xi),intent(in) :: highest_sz0,highest_sz1,highest_sz2

  real(8),dimension(max_number_xi),intent(in) :: lowest_charge0,lowest_charge1,lowest_charge2
  real(8),dimension(max_number_xi),intent(in) :: highest_charge0,highest_charge1,highest_charge2

  real(8),dimension(max_number_xi),intent(in) :: lowest_sxsy0,lowest_sxsy1,lowest_sxsy2
  real(8),dimension(max_number_xi),intent(in) :: highest_sxsy0,highest_sxsy1,highest_sxsy2

  real(8),dimension(max_number_xi),intent(in) :: lowest_s0,lowest_s1,lowest_s2
  real(8),dimension(max_number_xi),intent(in) :: highest_s0,highest_s1,highest_s2
  real(8),dimension(max_number_xi),intent(in) :: lowest_super0,lowest_super1,lowest_super2
  real(8),dimension(max_number_xi),intent(in) :: highest_super0,highest_super1,highest_super2
  real(8),dimension(TOTAL_SITE_NUMBER),intent(in) :: lowest_charge_density0,lowest_charge_density1,lowest_charge_density2
  real(8),dimension(TOTAL_SITE_NUMBER),intent(in) :: highest_charge_density0,highest_charge_density1,highest_charge_density2
  real(8),intent(in) :: lowest_pole0,lowest_pole1,lowest_pole2,&
       highest_pole0,highest_pole1,highest_pole2
  real(8),intent(in) :: pow1_c0,pow1_c1
  real(8),intent(in) :: pow2_c0,pow2_c1,pow2_c2
  character(32),intent(out) :: name
  integer,intent(out) :: error
  real(8),dimension(:),allocatable :: tmp_xi
!
  integer :: i_count
! init
  name="output" ; error=0

! HFO

  call output_hfo(values,unitary_d_up,unitary_d_down)

! phys_conditions

  call output_phys_conditions(values);

! etc
  call output_etc(values,hf_iteration_result);

! energy

  call output_energy(values,zero_approx_energy,result_pow0,&
       lowest_pow0,highest_pow0,result_pow1,lowest_pow1,highest_pow1,&
       result_pow2,lowest_pow2,highest_pow2)

! cor
  call output_cor(values,max_number_xi,&
       result_sz0,lowest_sz0,highest_sz0,&
       result_sz1,lowest_sz1,highest_sz1,&
       result_sz2,lowest_sz2,highest_sz2,&
       result_sxsy0,lowest_sxsy0,highest_sxsy0,&
       result_sxsy1,lowest_sxsy1,highest_sxsy1,&
       result_sxsy2,lowest_sxsy2,highest_sxsy2,&
       result_s0,lowest_s0,highest_s0,&
       result_s1,lowest_s1,highest_s1,&
       result_s2,lowest_s2,highest_s2,&
       result_charge0,lowest_charge0,highest_charge0,&
       result_charge1,lowest_charge1,highest_charge1,&
       result_charge2,lowest_charge2,highest_charge2,&
       result_super0,lowest_super0,highest_super0,&
       result_super1,lowest_super1,highest_super1,&
       result_super2,lowest_super2,highest_super2,&
       result_sz_zero,result_sxsy_zero,result_s_zero,&
       result_charge_zero,result_super_zero);

! output_pole
  call output_pole(values,&
       result_pole0,lowest_pole0,highest_pole0,&
       result_pole1,lowest_pole1,highest_pole1,&
       result_pole2,lowest_pole2,highest_pole2,&
       result_pole_zero)

! xi_table

  allocate(tmp_xi(NUMBER_XI))
  tmp_xi=0

  do i_count=1,NUMBER_XI
     tmp_xi(i_count)=seq_xi(i_count)
  end do

  call output_xi(values,tmp_xi)

! state_vector
  call output_vector(values,state_vector)

! charge_density
  call output_charge_density(values,&
       result_charge_density0,lowest_charge_density0,highest_charge_density0,&
       result_charge_density1,lowest_charge_density1,highest_charge_density1,&
       result_charge_density2,lowest_charge_density2,highest_charge_density2,&
       result_charge_density_zero)

! coef
  
  call output_coefficient(values,pow1_c0,pow1_c1,pow2_c0,pow2_c1,pow2_c2);

! calc_time
  call output_calc_time(values,time_start,time_end,NUMBER_PE)

end subroutine output

!
! $Id: input_cor_number_power.f90,v 1.1 2004/02/15 10:52:35 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine input_splited_cor_number_power(sample_number,&
     sz_nume1,sz_nume2,sz_nume3_1,sz_nume3_2,sz_nume4,sz_nume5,&
     charge_nume1,charge_nume2,charge_nume3_1,&
     charge_nume3_2,charge_nume4,charge_nume5,&
     charge_density_nume1,charge_density_nume2,&
     charge_density_nume3_1,charge_density_nume3_2,&
     charge_density_nume4,charge_density_nume5,&
     pole_nume1,pole_nume2,&
     pole_nume3_1,pole_nume3_2,&
     pole_nume4,pole_nume5,name,error)
  use global_variables
  implicit none

  integer,intent(in) :: sample_number
  real(8),dimension(max_number_xi),intent(in) :: sz_nume1,sz_nume2,sz_nume3_1,&
       sz_nume3_2,sz_nume4,sz_nume5
  real(8),dimension(max_number_xi),intent(in) :: charge_nume1,charge_nume2,&
       charge_nume3_1,charge_nume3_2,charge_nume4,charge_nume5
  real(8),dimension(TOTAL_SITE_NUMBER),intent(in) :: charge_density_nume1,&
       charge_density_nume2,charge_density_nume3_1,charge_density_nume3_2,&
       charge_density_nume4,charge_density_nume5
  real(8),dimension(TOTAL_SITE_NUMBER),intent(in) :: pole_nume1,&
       pole_nume2,pole_nume3_1,pole_nume3_2,pole_nume4,pole_nume5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
! init
  name="input_splited_cor_number" ; error=0

  if( mod(sample_number,5)==1 ) then
     do i_count=1,max_number_xi
        splited_sz1(1,i_count)=splited_sz1(1,i_count) &
             +sz_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz2(1,i_count)=splited_sz2(1,i_count) &
             +sz_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz3_1(1,i_count)=splited_sz3_1(1,i_count) &
             +sz_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz3_2(1,i_count)=splited_sz3_2(1,i_count) &
             +sz_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz4(1,i_count)=splited_sz4(1,i_count) &
             +sz_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz5(1,i_count)=splited_sz5(1,i_count) &
             +sz_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)

        splited_charge1(1,i_count)=splited_charge1(1,i_count) &
             +charge_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge2(1,i_count)=splited_charge2(1,i_count) &
             +charge_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge3_1(1,i_count)=splited_charge3_1(1,i_count) &
             +charge_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge3_2(1,i_count)=splited_charge3_2(1,i_count) &
             +charge_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge4(1,i_count)=splited_charge4(1,i_count) &
             +charge_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge5(1,i_count)=splited_charge5(1,i_count) &
             +charge_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density1(1,i_count)=splited_charge_density1(1,i_count) &
             +charge_density_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density2(1,i_count)=splited_charge_density2(1,i_count) &
             +charge_density_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density3_1(1,i_count)= &
             splited_charge_density3_1(1,i_count) &
             +charge_density_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density3_2(1,i_count)= &
             splited_charge_density3_2(1,i_count) &
             +charge_density_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density4(1,i_count)=splited_charge_density4(1,i_count) &
             +charge_density_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density5(1,i_count)=splited_charge_density5(1,i_count) &
             +charge_density_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)

        splited_pole1(1,i_count)=splited_pole1(1,i_count) &
             +pole_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole2(1,i_count)=splited_pole2(1,i_count) &
             +pole_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole3_1(1,i_count)=splited_pole3_1(1,i_count) &
             +pole_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole3_2(1,i_count)=splited_pole3_2(1,i_count) &
             +pole_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole4(1,i_count)=splited_pole4(1,i_count) &
             +pole_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole5(1,i_count)=splited_pole5(1,i_count) &
             +pole_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
     end do
  else if( mod(sample_number,5)==2 ) then
     do i_count=1,max_number_xi
        splited_sz1(2,i_count)=splited_sz1(2,i_count) &
             +sz_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz2(2,i_count)=splited_sz2(2,i_count) &
             +sz_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz3_1(2,i_count)=splited_sz3_1(2,i_count) &
             +sz_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz3_2(2,i_count)=splited_sz3_2(2,i_count) &
             +sz_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz4(2,i_count)=splited_sz4(2,i_count) &
             +sz_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz5(2,i_count)=splited_sz5(2,i_count) &
             +sz_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)

        splited_charge1(2,i_count)=splited_charge1(2,i_count) &
             +charge_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge2(2,i_count)=splited_charge2(2,i_count) &
             +charge_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge3_1(2,i_count)=splited_charge3_1(2,i_count) &
             +charge_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge3_2(2,i_count)=splited_charge3_2(2,i_count) &
             +charge_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge4(2,i_count)=splited_charge4(2,i_count) &
             +charge_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge5(2,i_count)=splited_charge5(2,i_count) &
             +charge_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density1(2,i_count)=splited_charge_density1(2,i_count) &
             +charge_density_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density2(2,i_count)=splited_charge_density2(2,i_count) &
             +charge_density_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density3_1(2,i_count)= &
             splited_charge_density3_1(2,i_count) &
             +charge_density_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density3_2(2,i_count)= &
             splited_charge_density3_2(2,i_count) &
             +charge_density_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density4(2,i_count)=splited_charge_density4(2,i_count) &
             +charge_density_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density5(2,i_count)=splited_charge_density5(2,i_count) &
             +charge_density_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)

        splited_pole1(2,i_count)=splited_pole1(2,i_count) &
             +pole_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole2(2,i_count)=splited_pole2(2,i_count) &
             +pole_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole3_1(2,i_count)=splited_pole3_1(2,i_count) &
             +pole_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole3_2(2,i_count)=splited_pole3_2(2,i_count) &
             +pole_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole4(2,i_count)=splited_pole4(2,i_count) &
             +pole_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole5(2,i_count)=splited_pole5(2,i_count) &
             +pole_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
     end do


  else if( mod(sample_number,5)==3 ) then
     do i_count=1,max_number_xi
        splited_sz1(3,i_count)=splited_sz1(3,i_count) &
             +sz_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz2(3,i_count)=splited_sz2(3,i_count) &
             +sz_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz3_1(3,i_count)=splited_sz3_1(3,i_count) &
             +sz_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz3_2(3,i_count)=splited_sz3_2(3,i_count) &
             +sz_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz4(3,i_count)=splited_sz4(3,i_count) &
             +sz_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz5(3,i_count)=splited_sz5(3,i_count) &
             +sz_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)

        splited_charge1(3,i_count)=splited_charge1(3,i_count) &
             +charge_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge2(3,i_count)=splited_charge2(3,i_count) &
             +charge_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge3_1(3,i_count)=splited_charge3_1(3,i_count) &
             +charge_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge3_2(3,i_count)=splited_charge3_2(3,i_count) &
             +charge_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge4(3,i_count)=splited_charge4(3,i_count) &
             +charge_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge5(3,i_count)=splited_charge5(3,i_count) &
             +charge_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density1(3,i_count)=splited_charge_density1(3,i_count) &
             +charge_density_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density2(3,i_count)=splited_charge_density2(3,i_count) &
             +charge_density_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density3_1(3,i_count)= &
             splited_charge_density3_1(3,i_count) &
             +charge_density_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density3_2(3,i_count)= &
             splited_charge_density3_2(3,i_count) &
             +charge_density_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density4(3,i_count)=splited_charge_density4(3,i_count) &
             +charge_density_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density5(3,i_count)=splited_charge_density5(3,i_count) &
             +charge_density_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)

        splited_pole1(3,i_count)=splited_pole1(3,i_count) &
             +pole_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole2(3,i_count)=splited_pole2(3,i_count) &
             +pole_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole3_1(3,i_count)=splited_pole3_1(3,i_count) &
             +pole_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole3_2(3,i_count)=splited_pole3_2(3,i_count) &
             +pole_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole4(3,i_count)=splited_pole4(3,i_count) &
             +pole_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole5(3,i_count)=splited_pole5(3,i_count) &
             +pole_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
     end do
  else if( mod(sample_number,5)==4 ) then
     do i_count=1,max_number_xi
        splited_sz1(4,i_count)=splited_sz1(4,i_count) &
             +sz_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz2(4,i_count)=splited_sz2(4,i_count) &
             +sz_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz3_1(4,i_count)=splited_sz3_1(4,i_count) &
             +sz_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz3_2(4,i_count)=splited_sz3_2(4,i_count) &
             +sz_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz4(4,i_count)=splited_sz4(4,i_count) &
             +sz_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz5(4,i_count)=splited_sz5(4,i_count) &
             +sz_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)

        splited_charge1(4,i_count)=splited_charge1(4,i_count) &
             +charge_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge2(4,i_count)=splited_charge2(4,i_count) &
             +charge_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge3_1(4,i_count)=splited_charge3_1(4,i_count) &
             +charge_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge3_2(4,i_count)=splited_charge3_2(4,i_count) &
             +charge_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge4(4,i_count)=splited_charge4(4,i_count) &
             +charge_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge5(4,i_count)=splited_charge5(4,i_count) &
             +charge_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density1(4,i_count)=splited_charge_density1(4,i_count) &
             +charge_density_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density2(4,i_count)=splited_charge_density2(4,i_count) &
             +charge_density_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density3_1(4,i_count)= &
             splited_charge_density3_1(4,i_count) &
             +charge_density_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density3_2(4,i_count)= &
             splited_charge_density3_2(4,i_count) &
             +charge_density_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density4(4,i_count)=splited_charge_density4(4,i_count) &
             +charge_density_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density5(4,i_count)=splited_charge_density5(4,i_count) &
             +charge_density_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)

        splited_pole1(4,i_count)=splited_pole1(4,i_count) &
             +pole_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole2(4,i_count)=splited_pole2(4,i_count) &
             +pole_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole3_1(4,i_count)=splited_pole3_1(4,i_count) &
             +pole_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole3_2(4,i_count)=splited_pole3_2(4,i_count) &
             +pole_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole4(4,i_count)=splited_pole4(4,i_count) &
             +pole_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole5(4,i_count)=splited_pole5(4,i_count) &
             +pole_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
     end do
  else if( mod(sample_number,5)==0 ) then
     do i_count=1,max_number_xi
        splited_sz1(5,i_count)=splited_sz1(5,i_count) &
             +sz_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz2(5,i_count)=splited_sz2(5,i_count) &
             +sz_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz3_1(5,i_count)=splited_sz3_1(5,i_count) &
             +sz_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz3_2(5,i_count)=splited_sz3_2(5,i_count) &
             +sz_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz4(5,i_count)=splited_sz4(5,i_count) &
             +sz_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sz5(5,i_count)=splited_sz5(5,i_count) &
             +sz_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)

        splited_charge1(5,i_count)=splited_charge1(5,i_count) &
             +charge_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge2(5,i_count)=splited_charge2(5,i_count) &
             +charge_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge3_1(5,i_count)=splited_charge3_1(5,i_count) &
             +charge_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge3_2(5,i_count)=splited_charge3_2(5,i_count) &
             +charge_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge4(5,i_count)=splited_charge4(5,i_count) &
             +charge_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge5(5,i_count)=splited_charge5(5,i_count) &
             +charge_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) &
             /dble(TOTAL_SITE_NUMBER)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density1(5,i_count)=splited_charge_density1(5,i_count) &
             +charge_density_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density2(5,i_count)=splited_charge_density2(5,i_count) &
             +charge_density_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density3_1(5,i_count)= &
             splited_charge_density3_1(5,i_count) &
             +charge_density_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density3_2(5,i_count)= &
             splited_charge_density3_2(5,i_count) &
             +charge_density_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density4(5,i_count)=splited_charge_density4(5,i_count) &
             +charge_density_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_charge_density5(5,i_count)=splited_charge_density5(5,i_count) &
             +charge_density_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)

        splited_pole1(5,i_count)=splited_pole1(5,i_count) &
             +pole_nume1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole2(5,i_count)=splited_pole2(5,i_count) &
             +pole_nume2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole3_1(5,i_count)=splited_pole3_1(5,i_count) &
             +pole_nume3_1(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole3_2(5,i_count)=splited_pole3_2(5,i_count) &
             +pole_nume3_2(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
        splited_pole4(5,i_count)=splited_pole4(5,i_count) &
             +pole_nume4(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
        splited_pole5(5,i_count)=splited_pole5(5,i_count) &
             +pole_nume5(i_count)/dble(MAX_MONTECARLO_SAMPLE/5) 
     end do
  end if

end subroutine input_splited_cor_number_power

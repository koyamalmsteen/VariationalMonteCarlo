!
! $Id: input_splited_cor_number.f90,v 1.2 2004/03/10 03:04:02 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine input_splited_cor_number(sample_number,tmp_sz0,tmp_charge0,&
     tmp_charge_density0,name,error)
  use global_variables
  implicit none

  integer,intent(in) :: sample_number
  real(8),dimension(max_number_xi),intent(in) :: tmp_sz0,tmp_charge0
  real(8),dimension(TOTAL_SITE_NUMBER),intent(in) :: tmp_charge_density0
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
! init
  name="input_splited_cor_number" ; error=0

  if( mod(sample_number,5)==1 ) then
     do i_count=1,max_number_xi
        splited_sz0(1,i_count)=splited_sz0(1,i_count) &
             +tmp_sz0(i_count)/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
        splited_charge0(1,i_count)=splited_charge0(1,i_count) &
             +tmp_charge0(i_count)/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density0(1,i_count)=splited_charge_density0(1,i_count) &
             +tmp_charge_density0(i_count)/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_pole0(1)=splited_pole0(1) &
             + xaxis(i_count)*(1-tmp_charge_density0(i_count)) &
             / dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
  else if( mod(sample_number,5)==2 ) then
     do i_count=1,max_number_xi
        splited_sz0(2,i_count)=splited_sz0(2,i_count) & 
             +tmp_sz0(i_count)/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
        splited_charge0(2,i_count)=splited_charge0(2,i_count) &
             +tmp_charge0(i_count)/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density0(2,i_count)=splited_charge_density0(2,i_count)&
             +tmp_charge_density0(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_pole0(2)=splited_pole0(2) &
             +xaxis(i_count)*(1-tmp_charge_density0(i_count)) &
             /dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
  else if( mod(sample_number,5)==3 ) then
     do i_count=1,max_number_xi
        splited_sz0(3,i_count)=splited_sz0(3,i_count) &
             +tmp_sz0(i_count)/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
        splited_charge0(3,i_count)=splited_charge0(3,i_count) &
             +tmp_charge0(i_count)/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density0(3,i_count)=splited_charge_density0(3,i_count)&
             +tmp_charge_density0(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_pole0(3)=splited_pole0(3) &
             +xaxis(i_count)*(1-tmp_charge_density0(i_count)) &
             /dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
  else if( mod(sample_number,5)==4 ) then
     do i_count=1,max_number_xi
        splited_sz0(4,i_count)=splited_sz0(4,i_count) &
             +tmp_sz0(i_count)/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
        splited_charge0(4,i_count)=splited_charge0(4,i_count) &
             +tmp_charge0(i_count)/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density0(4,i_count)=splited_charge_density0(4,i_count)&
             +tmp_charge_density0(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_pole0(4)=splited_pole0(4) &
             +xaxis(i_count)*(1-tmp_charge_density0(i_count)) &
             /dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
  else if( mod(sample_number,5)==0 ) then
     do i_count=1,max_number_xi
        splited_sz0(5,i_count)=splited_sz0(5,i_count) &
             +tmp_sz0(i_count)/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
        splited_charge0(5,i_count)=splited_charge0(5,i_count) &
             +tmp_charge0(i_count)/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density0(5,i_count)=splited_charge_density0(5,i_count)&
             +tmp_charge_density0(i_count)/dble(MAX_MONTECARLO_SAMPLE/5)
     end do
     do i_count=1,TOTAL_SITE_NUMBER
        splited_pole0(5)=splited_pole0(5) &
             +xaxis(i_count)*(1-tmp_charge_density0(i_count)) &
             /dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
  end if

end subroutine input_splited_cor_number

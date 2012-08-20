!
! $Id: input_splited_energy_power.f90,v 1.2 2004/03/10 03:04:02 k-yukino Exp $
!
#include "parameter.h"

subroutine input_splited_energy_power(sample_number,tmp_h0,tmp_h1,tmp_h2,&
     tmp_h3,tmp_h4,tmp_h5,name,error)
  use global_variables
  implicit none

  integer,intent(in) :: sample_number 
  real(8),intent(in) :: tmp_h0,tmp_h1,tmp_h2,tmp_h3,tmp_h4,tmp_h5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  name="input_splited_energy_power" ; error=0

  if( mod(sample_number,5)==1 ) then
     splited_h0(1)=splited_h0(1)+tmp_h0
     splited_h1(1)=splited_h1(1)+tmp_h1
     splited_h2(1)=splited_h2(1)+tmp_h2
     splited_h3(1)=splited_h3(1)+tmp_h3
     splited_h4(1)=splited_h4(1)+tmp_h4
     splited_h5(1)=splited_h5(1)+tmp_h5
  else if( mod(sample_number,5)==2 ) then
     splited_h0(2)=splited_h0(2)+tmp_h0
     splited_h1(2)=splited_h1(2)+tmp_h1
     splited_h2(2)=splited_h2(2)+tmp_h2
     splited_h3(2)=splited_h3(2)+tmp_h3
     splited_h4(2)=splited_h4(2)+tmp_h4
     splited_h5(2)=splited_h5(2)+tmp_h5
  else if( mod(sample_number,5)==3 ) then
     splited_h0(3)=splited_h0(3)+tmp_h0
     splited_h1(3)=splited_h1(3)+tmp_h1
     splited_h2(3)=splited_h2(3)+tmp_h2
     splited_h3(3)=splited_h3(3)+tmp_h3
     splited_h4(3)=splited_h4(3)+tmp_h4
     splited_h5(3)=splited_h5(3)+tmp_h5
  else if( mod(sample_number,5)==4 ) then
     splited_h0(4)=splited_h0(4)+tmp_h0
     splited_h1(4)=splited_h1(4)+tmp_h1
     splited_h2(4)=splited_h2(4)+tmp_h2
     splited_h3(4)=splited_h3(4)+tmp_h3
     splited_h4(4)=splited_h4(4)+tmp_h4
     splited_h5(4)=splited_h5(4)+tmp_h5
  else if( mod(sample_number,5)==0 ) then
     splited_h0(5)=splited_h0(5)+tmp_h0
     splited_h1(5)=splited_h1(5)+tmp_h1
     splited_h2(5)=splited_h2(5)+tmp_h2
     splited_h3(5)=splited_h3(5)+tmp_h3
     splited_h4(5)=splited_h4(5)+tmp_h4
     splited_h5(5)=splited_h5(5)+tmp_h5
  end if

end subroutine input_splited_energy_power

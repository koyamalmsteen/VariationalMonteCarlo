!
! $Id: cor_number.f90,v 1.4 2004/03/10 03:04:02 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine cor_number(tmp_sz0,tmp_charge0,tmp_charge_density0,name,error)
  use global_variables
  implicit none

  include "mpif.h"
  
  real(8),dimension(max_number_xi),intent(out) :: tmp_sz0,tmp_charge0
  real(8),dimension(TOTAL_SITE_NUMBER),intent(out) :: tmp_charge_density0
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8) :: tmp_result
  integer ::  i_count,j_count
! init
  tmp_sz0=0 ; tmp_charge0=0 ; tmp_charge_density0=0
  name="cor_number" ; error=0
  tmp_result=0

! $B$3$N%k!<%A%s$KF~$kA0$K(Bglobal_site_table$B$,%;%C%H$5$l$F$$$k(B

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=i_count,TOTAL_SITE_NUMBER ! $B%@%V%k%+%&%s%H$r$7$J$$$?$a(B
        tmp_result=0

        call sz(i_count,j_count,tmp_result,name,error)
        call error_check(name,error)

        tmp_sz0(distance_sequence(i_count,j_count))&
             &=tmp_sz0(distance_sequence(i_count,j_count)) &
             +tmp_result
     end do
  end do

! charge
  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=i_count,TOTAL_SITE_NUMBER ! $B%@%V%k%+%&%s%H$r$7$J$$$?$a(B
        tmp_result=0

        call charge(i_count,j_count,tmp_result,name,error)
        call error_check(name,error)

        tmp_charge0(distance_sequence(i_count,j_count))&
             &=tmp_charge0(distance_sequence(i_count,j_count))+tmp_result
     end do
  end do

!
! $BEE2YL)EY$N7W;;(B($BAj4X4X?t$G$O$J$$!#$h$C$F%5%$%H?t$NBg$-$5$r;}$D(B)
!

! charge_density
  tmp_charge_density0=0

  do i_count=1,TOTAL_SITE_NUMBER
     tmp_result=0
     call charge_density(i_count,tmp_result,name,error)
     call error_check(name,error)

     tmp_charge_density0(i_count)=tmp_charge_density0(i_count)+tmp_result
  end do

end subroutine cor_number

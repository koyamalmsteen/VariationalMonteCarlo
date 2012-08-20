!
! $Id: check_flat_site_number.f90,v 1.1 2003/06/04 10:11:35 k-yukino Exp $
!
#include "parameter.h"

subroutine check_flat_site_number(number_unit,name,error)
  use global_variables
  implicit none

  integer,intent(out) :: number_unit
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count,j_count,k_count
  integer :: flat_site_number
! init
  number_unit=0 ; name="check_flat_site_number" ; error=0
  flat_site_number=0

! フラットバンドに関しては、構造単位からサイト数を割り出す

  if( FLAT_BAND==1 )then
     do k_count=1,1000
        flat_site_number=0
        do i_count=0,2*k_count
           do j_count=0,2*k_count
              if( mod(i_count,2)==0 .or. mod(j_count,2)==0 ) then
                 flat_site_number=flat_site_number+1
              end if
           end do
        end do

        if( flat_site_number==TOTAL_SITE_NUMBER ) then
           number_unit=k_count ; error=0 ; return
        end if
     end do
  else if( FLAT_BAND==2 )then  
     do k_count=1,1000
        flat_site_number=0
        do i_count=0,2*k_count
           do j_count=0,2*k_count
              if( mod(i_count+j_count,2)/=0 ) then
                 flat_site_number=flat_site_number+1
              end if
           end do
        end do

        if( flat_site_number==TOTAL_SITE_NUMBER ) then
           number_unit=k_count ; error=0 ; return
        end if
     end do
  else if( FLAT_BAND==3 )then
     do k_count=1,1000
        flat_site_number=0
        do i_count=0,2*k_count
           do j_count=0,2*k_count
              if( mod(i_count+j_count,2)==0 ) then
                 flat_site_number=flat_site_number+1
              end if
           end do
        end do

        if( flat_site_number==TOTAL_SITE_NUMBER ) then
           number_unit=k_count ; error=0 ; return
        end if
     end do
  else
     error=-1 ; return
  end if

  error=-1

end subroutine check_flat_site_number

!
! $Id: init_global_variables_mielke.f90,v 1.1 2003/04/30 12:12:04 k-yukino Exp $
!
#include "parameter.h"

subroutine init_global_variables_mielke(number_unit,name,error)
  use global_variables
  implicit none

  type coordinate
     integer :: xaxis
     integer :: yaxis
  end type coordinate

  integer,intent(in) :: number_unit
  character(32),intent(out) :: name
  integer,intent(out) :: error
!
  type(coordinate),dimension(:),allocatable :: tmp_site_table
  integer :: i_count,j_count,k_count
  integer :: ier,temp1,temp2
! init
  name="init_global_variables_lieb" ; error=0
  ier=0 ; temp1=0 ; temp2=0

  allocate(tmp_site_table(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("tmp_site_table","init_global_variables_mielke",1,ier)
  tmp_site_table%xaxis=0 ; tmp_site_table%yaxis=0

!
! ラティス・ポイントの座標を書き出す
!
  k_count=1

  do i_count=0,2*number_unit
     do j_count=0,2*number_unit
        if( mod(i_count+j_count,2)/=0 ) then
           tmp_site_table(k_count)%xaxis=i_count
           tmp_site_table(k_count)%yaxis=j_count*(-1)
           k_count=k_count+1
        end if
     end do
  end do
!
! distance_tableを作る(やはりここでのdistance_tableもsqrtなしで定義)
!

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,i_count
        distance_table(i_count,j_count)=( (tmp_site_table(i_count)%xaxis &
             -tmp_site_table(j_count)%xaxis)**2&
             +(tmp_site_table(i_count)%yaxis-tmp_site_table(j_count)%yaxis)**2&
             )
        distance_table(j_count,i_count)=distance_table(i_count,j_count)
     end do
  end do

!
! 隣接テーブルを作る
!
  neighbor_table=0            ! 初期化

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        if( mod(tmp_site_table(i_count)%xaxis,2)==0 ) then  ! even number
           if( tmp_site_table(i_count)%xaxis==&
                tmp_site_table(j_count)%xaxis .and. &
                distance_table(i_count,j_count)==4 ) then
              neighbor_table(i_count,j_count)=1
              neighbor_table(j_count,i_count)=1
           end if

           if( tmp_site_table(i_count)%xaxis+1==tmp_site_table(j_count)%xaxis &
                .and. distance_table(i_count,j_count)==2 ) then
              neighbor_table(i_count,j_count)=1
              neighbor_table(j_count,i_count)=1
           end if

        else                                                ! odd number
           if( tmp_site_table(i_count)%yaxis==&
                tmp_site_table(j_count)%yaxis .and. &
                distance_table(i_count,j_count)==4 ) then
              neighbor_table(i_count,j_count)=1
              neighbor_table(j_count,i_count)=1
           end if

           if( tmp_site_table(i_count)%xaxis+1==tmp_site_table(j_count)%xaxis &
                .and. distance_table(i_count,j_count)==2 ) then
              neighbor_table(i_count,j_count)=1
              neighbor_table(j_count,i_count)=1
           end if

        end if
     end do
  end do

end subroutine init_global_variables_mielke

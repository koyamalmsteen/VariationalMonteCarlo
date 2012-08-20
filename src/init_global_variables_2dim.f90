!
! $Id: init_global_variables_2dim.f90,v 1.5 2004/03/10 03:04:02 k-yukino Exp k-yukino $
!

#include "parameter.h"

subroutine init_global_variables_2dim(zahyou,name,error)
  use global_variables
  implicit none

  type coordinate
     integer :: xaxis
     integer :: yaxis
  end type coordinate
  
  integer,dimension(TOTAL_SITE_NUMBER,2) :: zahyou
  character(32),intent(out) :: name
  integer,intent(out) :: error
!
  type(coordinate),dimension(:,:),allocatable :: tmp_site_table
  type(coordinate),dimension(TOTAL_SITE_NUMBER) :: tmptmp_site_table
  integer :: i_count,j_count,k_count
  integer :: ier,temp1,temp2
  integer :: delta_x,delta_y
! init
  name="init_global_variables_2dim" ; error=0
  ier=0 ; temp1=0 ; temp2=0
  delta_x=0 ; delta_y=0
!
  allocate(tmp_site_table(TOTAL_SITE_NUMBER,9),stat=ier)
  call stat_check("tmp_site_table","init_global_tables_for_2dim_normal",1,ier)
  tmp_site_table%xaxis=0 ; tmp_site_table%yaxis=0
  
  call set_2d_system(delta_x,delta_y,tmptmp_site_table,name,error)
  call error_check(name,error)

  do i_count=1,TOTAL_SITE_NUMBER
     tmp_site_table(i_count,5)%xaxis=tmptmp_site_table(i_count)%xaxis
     tmp_site_table(i_count,5)%yaxis=tmptmp_site_table(i_count)%yaxis
  end do

!! No.1
  do i_count=1,TOTAL_SITE_NUMBER
     tmp_site_table(i_count,1)%xaxis=tmp_site_table(i_count,5)%xaxis &
          +(delta_x-delta_y)
     tmp_site_table(i_count,1)%yaxis=tmp_site_table(i_count,5)%yaxis &
          +(delta_x+delta_y)
!! No.2
     tmp_site_table(i_count,2)%xaxis=tmp_site_table(i_count,5)%xaxis &
          -delta_y
     tmp_site_table(i_count,2)%yaxis=tmp_site_table(i_count,5)%yaxis &
          +delta_x

!! No.3
     tmp_site_table(i_count,3)%xaxis=tmp_site_table(i_count,5)%xaxis &
          -(delta_x+delta_y)
     tmp_site_table(i_count,3)%yaxis=tmp_site_table(i_count,5)%yaxis &
          +(delta_x-delta_y)

!! No.4
     tmp_site_table(i_count,4)%xaxis=tmp_site_table(i_count,5)%xaxis &
          +delta_x
     tmp_site_table(i_count,4)%yaxis=tmp_site_table(i_count,5)%yaxis &
          +delta_y

!! No.6
     tmp_site_table(i_count,6)%xaxis=tmp_site_table(i_count,5)%xaxis &
          -delta_x
     tmp_site_table(i_count,6)%yaxis=tmp_site_table(i_count,5)%yaxis &
          -delta_y

!! No.7
     tmp_site_table(i_count,7)%xaxis=tmp_site_table(i_count,5)%xaxis &
          +(delta_x+delta_y)
     tmp_site_table(i_count,7)%yaxis=tmp_site_table(i_count,5)%yaxis &
          -(delta_x-delta_y)

!! No.8
     tmp_site_table(i_count,8)%xaxis=tmp_site_table(i_count,5)%xaxis &
          +delta_y
     tmp_site_table(i_count,8)%yaxis=tmp_site_table(i_count,5)%yaxis &
          -delta_x

!! No.9
     tmp_site_table(i_count,9)%xaxis=tmp_site_table(i_count,5)%xaxis &
          -(delta_x-delta_y)
     tmp_site_table(i_count,9)%yaxis=tmp_site_table(i_count,5)%yaxis &
          -(delta_x+delta_y)
  end do

!
! distanceテーブルを作る
!

! 仮の再短距離を入力(この段階では周期的境界条件を仮定していない距離である)
! distance_tableはsqrtはしないという定義のもとプログラムを作成する。

  temp1=0

  do j_count=1,TOTAL_SITE_NUMBER
     do i_count=1,TOTAL_SITE_NUMBER
        temp1=(tmp_site_table(i_count,5)%xaxis&
             & -tmp_site_table(j_count,5)%xaxis)**2&
             & +(tmp_site_table(i_count,5)%yaxis&
             & -tmp_site_table(j_count,5)%yaxis)**2

        distance_table(i_count,j_count)=temp1
     end do
  end do

! ここからちゃんと比較する。
  do k_count=1,9
     temp1=0
     temp2=0
     do j_count=1,TOTAL_SITE_NUMBER
        do i_count=1,TOTAL_SITE_NUMBER
           temp1=(tmp_site_table(i_count,5)%xaxis&
                & -tmp_site_table(j_count,k_count)%xaxis)**2&
                & +(tmp_site_table(i_count,5)%yaxis&
                & -tmp_site_table(j_count,k_count)%yaxis)**2
           
           if( distance_table(i_count,j_count)>temp1 ) then
              distance_table(i_count,j_count)=temp1
           end if
        end do
     end do
  end do

!
! 座標の情報zahyou(TOTAL_SITE_NUMBER,2)にコピーしておく(電界で使う)
!     
  do i_count=1,TOTAL_SITE_NUMBER
     zahyou(i_count,1)=tmp_site_table(i_count,5)%xaxis
     zahyou(i_count,2)=tmp_site_table(i_count,5)%yaxis
     xaxis(i_count)=tmp_site_table(i_count,5)%xaxis+1     ! 原点を1とした
  end do

!
! 隣接テーブルを作る(隣あってたら1、違えば0)
!
  neighbor_table=0

  if( PERIODIC==0 ) then                        ! p-p 周期的境界条件を使う
     do k_count=1,9
        do j_count=1,TOTAL_SITE_NUMBER
           do i_count=1,TOTAL_SITE_NUMBER
              if( tmp_site_table(i_count,5)%xaxis== &
                   tmp_site_table(j_count,k_count)%xaxis & 
                   .and. abs(tmp_site_table(i_count,5)%yaxis &
                   -tmp_site_table(j_count,k_count)%yaxis)==1 .or. &
                   tmp_site_table(i_count,5)%yaxis== &
                   tmp_site_table(j_count,k_count)%yaxis & 
                   .and. abs(tmp_site_table(i_count,5)%xaxis &
                   -tmp_site_table(j_count,k_count)%xaxis)==1 ) then
                 neighbor_table(i_count,j_count)=1
              end if
           end do
        end do
     end do
  else if( PERIODIC==1 ) then                    ! a-a 使わない
     do j_count=1,TOTAL_SITE_NUMBER
        do i_count=1,TOTAL_SITE_NUMBER
           if( tmp_site_table(i_count,5)%xaxis== &
                tmp_site_table(j_count,5)%xaxis & 
                .and. abs(tmp_site_table(i_count,5)%yaxis &
                -tmp_site_table(j_count,5)%yaxis)==1 .or. &
                tmp_site_table(i_count,5)%yaxis== &
                tmp_site_table(j_count,5)%yaxis & 
                .and. abs(tmp_site_table(i_count,5)%xaxis &
                -tmp_site_table(j_count,5)%xaxis)==1 ) then
              neighbor_table(i_count,j_count)=1
           end if
        end do
     end do
  else if( PERIODIC==2 ) then                    ! p-a 一方向だけ(4,5,6)
     do k_count=4,6
        do j_count=1,TOTAL_SITE_NUMBER
           do i_count=1,TOTAL_SITE_NUMBER
              if( tmp_site_table(i_count,5)%xaxis== &
                   tmp_site_table(j_count,k_count)%xaxis & 
                   .and. abs(tmp_site_table(i_count,5)%yaxis &
                   -tmp_site_table(j_count,k_count)%yaxis)==1 .or. &
                   tmp_site_table(i_count,5)%yaxis== &
                   tmp_site_table(j_count,k_count)%yaxis & 
                   .and. abs(tmp_site_table(i_count,5)%xaxis &
                   -tmp_site_table(j_count,k_count)%xaxis)==1 ) then
                 neighbor_table(i_count,j_count)=1
              end if
           end do
        end do
     end do
  end if

!
! neighbor_table2を作る neighbor_table2(サイト番号,方向)
! 方向=1 : 右
! 方向=2 : 下
! 方向=3 : 左
! 方向=4 : 上
!
  neighbor_table2=0

  if( PERIODIC==0 ) then                        ! p-p 周期的境界条件を使う
     do k_count=1,9
        do j_count=1,TOTAL_SITE_NUMBER
           do i_count=1,TOTAL_SITE_NUMBER
! 方向=1 : 右
              if( tmp_site_table(i_count,5)%yaxis==&
                   tmp_site_table(j_count,k_count)%yaxis .and. &
                   tmp_site_table(i_count,5)%xaxis &
                   -tmp_site_table(j_count,k_count)%xaxis==-1 ) then
                 neighbor_table2(i_count,1)=j_count
              end if
! 方向=2 : 下
              if( tmp_site_table(i_count,5)%xaxis==&
                   tmp_site_table(j_count,k_count)%xaxis .and. &
                   tmp_site_table(i_count,5)%yaxis &
                   -tmp_site_table(j_count,k_count)%yaxis==1 ) then
                 neighbor_table2(i_count,2)=j_count
              end if
! 方向=3 : 左
              if( tmp_site_table(i_count,5)%yaxis==&
                   tmp_site_table(j_count,k_count)%yaxis .and. &
                   tmp_site_table(i_count,5)%xaxis &
                   -tmp_site_table(j_count,k_count)%xaxis==1 ) then
                 neighbor_table2(i_count,3)=j_count
              end if
! 方向=4 : 上
              if( tmp_site_table(i_count,5)%xaxis==&
                   tmp_site_table(j_count,k_count)%xaxis .and. &
                   tmp_site_table(i_count,5)%yaxis &
                   -tmp_site_table(j_count,k_count)%yaxis==-1 ) then
                 neighbor_table2(i_count,4)=j_count
              end if

           end do
        end do
     end do
  else if( PERIODIC==1 ) then                    ! a-a 使わない
     do i_count=1,TOTAL_SITE_NUMBER
        do j_count=1,TOTAL_SITE_NUMBER
! 方向=1 : 右
           if( tmp_site_table(i_count,5)%yaxis == &
                tmp_site_table(j_count,5)%yaxis .and. &
                tmp_site_table(i_count,5)%xaxis &
                -tmp_site_table(j_count,5)%xaxis==-1 ) then
              neighbor_table2(i_count,1)=j_count
           end if
! 方向=2 : 下
           if( tmp_site_table(i_count,5)%xaxis == &
                tmp_site_table(j_count,5)%xaxis .and. &
                tmp_site_table(i_count,5)%yaxis &
                -tmp_site_table(j_count,5)%yaxis==1 ) then
              neighbor_table2(i_count,2)=j_count
           end if
! 方向=3 : 左
           if( tmp_site_table(i_count,5)%yaxis == &
                tmp_site_table(j_count,5)%yaxis .and. &
                tmp_site_table(i_count,5)%xaxis &
                -tmp_site_table(j_count,5)%xaxis==1 ) then
              neighbor_table2(i_count,3)=j_count
           end if
! 方向=4 : 上
           if( tmp_site_table(i_count,5)%xaxis == &
                tmp_site_table(j_count,5)%xaxis .and. &
                tmp_site_table(i_count,5)%yaxis &
                -tmp_site_table(j_count,5)%yaxis==-1 ) then
              neighbor_table2(i_count,4)=j_count
           end if
        end do
     end do
  else if( PERIODIC==2 ) then                    ! p-a 一方向だけ(4,5,6)
!
! 出来ていない
!
!     do k_count=4,6
!        do j_count=1,TOTAL_SITE_NUMBER
!           do i_count=1,TOTAL_SITE_NUMBER
!              if( tmp_site_table(i_count,5)%xaxis== &
!                   tmp_site_table(j_count,k_count)%xaxis & 
!                   .and. abs(tmp_site_table(i_count,5)%yaxis &
!                   -tmp_site_table(j_count,k_count)%yaxis)==1 .or. &
!                   tmp_site_table(i_count,5)%yaxis== &
!                   tmp_site_table(j_count,k_count)%yaxis & 
!                   .and. abs(tmp_site_table(i_count,5)%xaxis &
!                   -tmp_site_table(j_count,k_count)%xaxis)==1 ) then
!                 neighbor_table(i_count,j_count)=1
!              end if
!           end do
!        end do
!     end do
  end if

  deallocate(tmp_site_table,stat=ier)
  call stat_check("tmp_site_table","init_global_tables_for_2dim_normal",2,ier)



end subroutine init_global_variables_2dim

!
! $Id: count_site_number.f90,v 1.1 2003/11/05 09:24:47 k-yukino Exp $
!
#include "parameter.h"

subroutine count_site_number(delta_x,delta_y,&
                             tmptmp_site_table,site_number,name,error)
  implicit none
  
  type coordinate
     integer :: xaxis
     integer :: yaxis
  end type coordinate

  integer,intent(in) :: delta_x,delta_y
  type(coordinate),dimension(TOTAL_SITE_NUMBER),intent(out) :: tmptmp_site_table
  integer,intent(out) :: site_number 
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count,j_count
  integer,dimension(4) :: xxx
  integer,dimension(4) :: yyy
  real(8),dimension(4) :: grad
! init
  tmptmp_site_table%xaxis=0 ; tmptmp_site_table%yaxis=0
  site_number=0
  name="count_site_number" ; error=0

! 座標の指定
  xxx(1)=0 ; yyy(1)=0
  xxx(2)=delta_x ; yyy(2)=delta_y
  xxx(4)=delta_y ; yyy(4)=-delta_x
  xxx(3)=xxx(2)+xxx(4) ; yyy(3)=yyy(2)+yyy(4)

  if( delta_y/=0 )then
     call set_grad(xxx(2),yyy(2),grad(1),name,error)
     call error_check(name,error)
     call set_grad(xxx(3)-xxx(2),yyy(3)-yyy(2),grad(2),name,error)
     call error_check(name,error)
     call set_grad(xxx(3)-xxx(4),yyy(3)-yyy(4),grad(3),name,error)
     call error_check(name,error)
     call set_grad(xxx(4),yyy(4),grad(4),name,error)
     call error_check(name,error)

     do i_count=0,delta_x+delta_y
        do j_count=delta_y,-delta_x,-1

           if( i_count==0 ) then
              if( j_count==0 ) then
                 site_number=1
                 tmptmp_site_table(site_number)%xaxis=i_count
                 tmptmp_site_table(site_number)%yaxis=j_count
              end if
           else                               ! i_countが0以外の場合
              if( delta_x==delta_y ) then      !! delta_x==delta_y
                 if( i_count<xxx(2) ) then
                    if( j_count<=grad(1)*i_count .and. &
                         j_count>=grad(4)*i_count ) then
                       site_number=site_number+1
                       tmptmp_site_table(site_number)%xaxis=i_count
                       tmptmp_site_table(site_number)%yaxis=j_count
                    end if
                 else if( xxx(2)==i_count ) then
                    if( j_count<grad(1)*i_count .and. &
                         j_count>grad(4)*i_count ) then
                       site_number=site_number+1
                       tmptmp_site_table(site_number)%xaxis=i_count
                       tmptmp_site_table(site_number)%yaxis=j_count
                    end if
                 else if( xxx(2)<i_count ) then
                    if( j_count<yyy(2)+grad(2)*(i_count-xxx(2)) .and. &
                         j_count>yyy(4)+grad(3)*(i_count-xxx(2)) ) then
                       site_number=site_number+1

                       tmptmp_site_table(site_number)%xaxis=i_count
                       tmptmp_site_table(site_number)%yaxis=j_count
                    end if
                 end if
              else                             !! delta_x/=delta_y
                 if( i_count<=xxx(4) ) then
                    if( j_count<=grad(1)*i_count .and. &
                         j_count>=grad(4)*i_count ) then
                       if( i_count==xxx(4) .and. j_count==yyy(4) ) then
! なにもしない
                       else
                          site_number=site_number+1
                          tmptmp_site_table(site_number)%xaxis=i_count
                          tmptmp_site_table(site_number)%yaxis=j_count
                       end if
                    end if
                 else if( i_count>xxx(4) .and. i_count<xxx(2) ) then
                    if( j_count<=grad(1)*i_count .and. &
                         j_count>yyy(4)+grad(3)*(i_count-xxx(4)) ) then

                       site_number=site_number+1
                       tmptmp_site_table(site_number)%xaxis=i_count
                       tmptmp_site_table(site_number)%yaxis=j_count
                    end if
                 else if( i_count>=xxx(2) ) then
                    if( j_count<yyy(2)+grad(2)*(i_count-xxx(2)) .and. &
                         j_count>yyy(4)+grad(3)*(i_count-xxx(4)) ) then
                       site_number=site_number+1
                       tmptmp_site_table(site_number)%xaxis=i_count
                       tmptmp_site_table(site_number)%yaxis=j_count
                    end if
                 end if
              end if
           end if

        end do
     end do

  else if( delta_y==0 ) then                    ! delta_y==0のときは
     site_number=0

     do i_count=0,delta_x-1
        do j_count=0,delta_x-1
           site_number=site_number+1
           tmptmp_site_table(site_number)%xaxis=i_count
           tmptmp_site_table(site_number)%yaxis=-j_count
        end do
     end do

  end if

end subroutine count_site_number




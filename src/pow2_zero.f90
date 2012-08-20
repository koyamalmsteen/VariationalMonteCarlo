!
! $Id: pow2_zero.f90,v 1.1 2003/08/14 09:24:24 k-yukino Exp $
!
#include "parameter.h"

subroutine pow2_zero(result,name,error)
  use global_variables
  implicit none

  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i1_count,i2_count
  integer :: j1_count,j2_count
  integer :: sigma1,sigma2
  integer :: k_count
  real(8),dimension(2) :: temp
  real(8) :: t1,t2,t3
! init
  result=0
  name="pow2_zero"
  error=0
!
  temp=0

  do i1_count=1,TOTAL_SITE_NUMBER
     do j1_count=1,TOTAL_SITE_NUMBER
        do i2_count=1,TOTAL_SITE_NUMBER
           do j2_count=1,TOTAL_SITE_NUMBER

              temp=0

              if( neighbor_table(i1_count,j1_count)==1 ) then
                 if( neighbor_table(i2_count,j2_count)==1 ) then
                    do sigma1=1,2
                       do sigma2=1,2

! 第1項
                          t1=0

                          if( sigma1==1 ) then
                             do k_count=1,TOTAL_UP_ELECTRON
                                t1=t1+unitary_d_up(k_count,i1_count)&
                                     *unitary_d_up(k_count,j1_count)
                             end do
                          else
                             do k_count=1,TOTAL_DOWN_ELECTRON
                                t1=t1+unitary_d_down(k_count,i1_count)&
                                     *unitary_d_down(k_count,j1_count)
                             end do
                          end if

                          t2=0
                          
                          if( sigma2==1 ) then
                             do k_count=1,TOTAL_UP_ELECTRON
                                t2=t2+unitary_d_up(k_count,i2_count)&
                                     *unitary_d_up(k_count,j2_count)
                             end do
                          else
                             do k_count=1,TOTAL_DOWN_ELECTRON
                                t2=t2+unitary_d_down(k_count,i2_count)&
                                     *unitary_d_down(k_count,j2_count)
                             end do
                          end if

                          temp(1)=t1*t2

! 第2項
                          t1=0
                          
                          if( sigma1==sigma2 ) then
                             if( sigma1==1 ) then
                                do k_count=1,TOTAL_UP_ELECTRON
                                   t1=t1+unitary_d_up(k_count,i1_count)&
                                        *unitary_d_up(k_count,j2_count)
                                end do
                             else
                                do k_count=1,TOTAL_DOWN_ELECTRON
                                   t1=t1+unitary_d_down(k_count,i1_count)&
                                        *unitary_d_down(k_count,j2_count)
                                end do
                             end if
                          else
                             t1=0
                          end if

                          t2=0

                          if( i2_count==j1_count .and. sigma1==sigma2 ) then
                             t2=1
                          else
                             t2=0
                          end if

                          t3=0 

                          if( sigma1==sigma2 ) then
                             if( sigma1==1 ) then
                                do k_count=1,TOTAL_UP_ELECTRON
                                   t3=t3+unitary_d_up(k_count,i2_count)&
                                        *unitary_d_up(k_count,j1_count)
                                end do
                             else
                                do k_count=1,TOTAL_DOWN_ELECTRON
                                   t3=t3+unitary_d_down(k_count,i2_count)&
                                        *unitary_d_down(k_count,j1_count)
                                end do
                             end if
                          end if

                          temp(2)=t1*(t2-t3)

                          result=result+dble(TRANSFER)*dble(TRANSFER)&
                               *(temp(1)+temp(2))
                       end do
                    end do
                 end if
              end if
           end do
        end do
     end do
  end do

end subroutine pow2_zero

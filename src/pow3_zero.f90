!
! $Id: pow3_zero.f90,v 1.1 2003/08/14 09:24:31 k-yukino Exp $
!
#include "parameter.h"

subroutine pow3_zero(result,name,error)
  use global_variables
  implicit none

  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i1_count,i2_count,i3_count
  integer :: j1_count,j2_count,j3_count
  integer :: sigma1,sigma2,sigma3
  integer :: k_count
  integer,external :: delta
  real(8),dimension(6) :: temp
  real(8) :: t1,t2,t3
! init
  result=0
  name="pow3_zero"
  error=0
  temp=0

! 第1項
  do i1_count=1,TOTAL_SITE_NUMBER
     do j1_count=1,TOTAL_SITE_NUMBER
        do i2_count=1,TOTAL_SITE_NUMBER
           do j2_count=1,TOTAL_SITE_NUMBER
              do i3_count=1,TOTAL_SITE_NUMBER
                 do j3_count=1,TOTAL_SITE_NUMBER

                    temp=0

                    if( neighbor_table(i1_count,j1_count)==1 ) then
                       if( neighbor_table(i2_count,j2_count)==1 ) then
                          if( neighbor_table(i3_count,j3_count)==1 ) then

                             do sigma1=1,2
                                do sigma2=1,2
                                   do sigma3=1,2

                                      temp=0
! 第1項
                                      t1=0 ; t2=0 ; t3=0

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

                                      if( sigma3==1 ) then
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t3=t3+unitary_d_up(k_count,i3_count)&
                                                 *unitary_d_up(k_count,j3_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t3=t3+unitary_d_down(k_count,i3_count)&
                                                 *unitary_d_down(k_count,j3_count)
                                         end do
                                      end if
                             
                                      temp(1)=temp(1)+t1*t2*t3

! 第2項
                                      t1=0 ; t2=0 ; t3=0

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

                                      if( sigma2==1 ) then
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t2=t2+unitary_d_up(k_count,i2_count)&
                                                 *unitary_d_up(k_count,j3_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t2=t2+unitary_d_down(k_count,i2_count)&
                                                 *unitary_d_down(k_count,j3_count)
                                         end do
                                      end if

                                      if( sigma3==1 ) then
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t3=t3+unitary_d_up(k_count,i3_count)&
                                                 *unitary_d_up(k_count,j2_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t3=t3+unitary_d_down(k_count,i3_count)&
                                                 *unitary_d_down(k_count,j2_count)
                                         end do
                                      end if
                             
                                      if( sigma2==sigma3 ) then
                                         if( i3_count==j2_count ) then
                                            temp(2)=temp(2)+t1*t2*(1-t3)
                                         else
                                            temp(2)=temp(2)+t1*t2*(-t3)
                                         end if
                                      end if

! 第3項
                                      t1=0 ; t2=0 ; t3=0
                                      
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

                                      if( sigma2==1 ) then 
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t2=t2+unitary_d_up(k_count,i2_count)&
                                                 *unitary_d_up(k_count,j1_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t2=t2+unitary_d_down(k_count,i2_count)&
                                                 *unitary_d_down(k_count,j1_count)
                                         end do
                                      end if

                                      if( sigma3==1 ) then
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t3=t3+unitary_d_up(k_count,i3_count)&
                                                 *unitary_d_up(k_count,j3_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t3=t3+unitary_d_down(k_count,i3_count)&
                                                 *unitary_d_down(k_count,j3_count)
                                         end do
                                      end if
                             
                                      if( sigma1==sigma2 ) then
                                         if( i2_count==j1_count )then
                                            temp(3)=temp(3)+t1*(1-t2)*t3
                                         else
                                            temp(3)=temp(3)+t1*(-t2)*t3
                                         end if
                                      end if

! 第4項
                                      t1=0 ; t2=0 ; t3=0
                                      
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

                                      if( sigma2==1 ) then
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t2=t2+unitary_d_up(k_count,i2_count)&
                                                 *unitary_d_up(k_count,j3_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t2=t2+unitary_d_down(k_count,i2_count)&
                                                 *unitary_d_down(k_count,j3_count)
                                         end do
                                      end if

                                      if( sigma3==1 ) then
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t3=t3+unitary_d_up(k_count,i3_count)&
                                                 *unitary_d_up(k_count,j1_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t3=t3+unitary_d_down(k_count,i3_count)&
                                                 *unitary_d_down(k_count,j1_count)
                                         end do
                                      end if
                             
                                      if( sigma1==sigma2 .and. &
                                           sigma2==sigma3 .and. &
                                           sigma1==sigma3 ) then
                                         if( i3_count==j1_count )then
                                            temp(4)=temp(4)-t1*t2*(1-t3)
                                         else
                                            temp(4)=temp(4)-t1*t2*(-t3)
                                         end if

                                      end if
! 第5項
                                      t1=0 ; t2=0 ; t3=0
                                      
                                      if( sigma1==1 ) then
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t1=t1+unitary_d_up(k_count,i1_count)&
                                                 *unitary_d_up(k_count,j3_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t1=t1+unitary_d_down(k_count,i1_count)&
                                                 *unitary_d_down(k_count,j3_count)
                                         end do
                                      end if

                                      if( sigma2==1 ) then
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t2=t2+unitary_d_up(k_count,i2_count)&
                                                 *unitary_d_up(k_count,j1_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t2=t2+unitary_d_down(k_count,i2_count)&
                                                 *unitary_d_down(k_count,j1_count)
                                         end do
                                      end if

                                      if( sigma3==1 ) then
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t3=t3+unitary_d_up(k_count,i3_count)&
                                                 *unitary_d_up(k_count,j2_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t3=t3+unitary_d_down(k_count,i3_count)&
                                                 *unitary_d_down(k_count,j2_count)
                                         end do
                                      end if

                                      if( sigma1==sigma3 .and. &
                                           sigma1==sigma2 .and. &
                                           sigma2==sigma3 ) then
                                         if( i2_count==j1_count ) then
                                            if( i3_count==j2_count ) then
                                               temp(5)=temp(5)+t1*(1-t2)&
                                                    *(1-t3)
                                            else
                                               temp(5)=temp(5)+t1*(1-t2)&
                                                    *(-t3)
                                            end if
                                         else
                                            if( i3_count==j2_count ) then
                                               temp(5)=temp(5)+t1*(-t2)&
                                                    *(1-t3)
                                            else
                                               temp(5)=temp(5)+t1*(-t2)&
                                                    *(-t3)
                                            end if
                                         end if
                                      end if

! 第6項
                                      t1=0 ; t2=0 ; t3=0

                                      if( sigma1==1 ) then
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t1=t1+unitary_d_up(k_count,i1_count)&
                                                 *unitary_d_up(k_count,j3_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t1=t1+unitary_d_down(k_count,i1_count)&
                                                 *unitary_d_down(k_count,j3_count)
                                         end do
                                      end if

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

                                      if( sigma3==1 ) then
                                         do k_count=1,TOTAL_UP_ELECTRON
                                            t3=t3+unitary_d_up(k_count,i3_count)&
                                                 *unitary_d_up(k_count,j1_count)
                                         end do
                                      else
                                         do k_count=1,TOTAL_DOWN_ELECTRON
                                            t3=t3+unitary_d_down(k_count,i3_count)&
                                                 *unitary_d_down(k_count,j1_count)
                                         end do
                                      end if
                             
                                      if( sigma1==sigma3 ) then
                                         if( i3_count==j1_count ) then
                                            temp(6)=temp(6)+t1*t2*(1-t3)
                                         else
                                            temp(6)=temp(6)+t1*t2*(-t3)
                                         end if
                                      end if

                                      result=result+(temp(1)+temp(2)+temp(3)&
                                           +temp(4)+temp(5)+temp(6))&
                                           *dble(TRANSFER)*dble(TRANSFER)&
                                           *dble(TRANSFER)
                                   end do
                                end do
                             end do
                          end if
                       end if
                    end if
                 end do
              end do
           end do
        end do
     end do
  end do

end subroutine pow3_zero

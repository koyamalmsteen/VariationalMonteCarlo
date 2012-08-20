!
! $Id: super_nume5.f90,v 1.4 2003/11/20 14:04:19 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine super_nume5(alpha_h_h_psi1,vector_up,vector_down,aaa,bbb,result234,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(in) :: alpha_h_h_psi1
  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  integer,intent(in) :: aaa,bbb
  real(8),intent(out) :: result234
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8),dimension(TOTAL_UP_ELECTRON,TOTAL_UP_ELECTRON) :: d_tilde_up_inverse_store
  real(8),dimension(TOTAL_DOWN_ELECTRON,TOTAL_DOWN_ELECTRON) :: d_tilde_down_inverse_store
  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  integer,dimension(TOTAL_SITE_NUMBER) :: global_site_table_up_store,global_site_table_down_store
  integer :: i_count,j_count,k_count
  real(8),external :: jastrow
  real(8) :: temp_a,temp_b
  integer,dimension(TOTAL_UP_ELECTRON) :: right_vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: right_vector_down
  integer :: number_ket_vector_up,number_ket_vector_down
  integer,dimension(2*JIGEN*TOTAL_UP_ELECTRON+1,TOTAL_UP_ELECTRON) :: available_ket_up_table
  integer,dimension(2*JIGEN*TOTAL_DOWN_ELECTRON+1,TOTAL_DOWN_ELECTRON) :: available_ket_down_table
  real(8) :: result4
  integer :: double_occupy 
  real(8) :: result2,result3
  integer :: tmp
  integer :: b2u_count,b2d_count
  real(8) :: lambda
  real(8) :: input_up,input_down,input
  integer :: sign                              ! $BFb@Q$N:]$N$b$N(B
! init
  name="super_nume5" ; error=0
! init local
  temp_a=0 ; temp_b=0
  tmp=0
  result234=0 ; result2=0 ; result3=0 ; result4=0

!
  global_site_table_up_store=global_site_table_up
  global_site_table_down_store=global_site_table_down
  d_tilde_up_inverse_store=d_tilde_up_inverse
  d_tilde_down_inverse_store=d_tilde_down_inverse
!
! <psi|H^2|alpha><alpha|O|beta1><beta1|H|beta2><beta2|H|psi>
!                 ^^^^^^^^^^^^^^^   ^^^^^^^^^^^^  ^^^^^^^^^^^^
!                  result2       result3       result4
!
  
!
! result234$B$r5a$a$k(B
!
! <alpha|$B$H7k$SIU$/(B|beta_1>$B$G7W;;$r9T$J$&!#(B
!
  do i_count=1,4
     do j_count=1,4

        if( neighbor_table2(aaa,i_count)/=0 .and. &
             neighbor_table2(bbb,j_count)/=0 ) then 

!!
!! $BBh(B1$B9`$N7W;;(B (aaa;i_count bbb;j_count)
!!

! $B%J%s%P!<%*%Z%l!<%?!<$N>l9g(B
           if( aaa==bbb .and. &
                neighbor_table2(aaa,i_count)== &
                neighbor_table2(bbb,j_count) ) then
              
              if( global_site_table_down_store(aaa)==1 .and. &
                   global_site_table_up_store&
                   &(neighbor_table2(aaa,i_count))==1 ) then
                 temp_a=-input

              end if
           else
              if( global_site_table_up_store(aaa)==0 .and. &
                   global_site_table_down_store(bbb)==1 .and. &
                   global_site_table_up_store&
                   &(neighbor_table2(aaa,i_count))==0 .and. &
                   global_site_table_up_store&
                   &(neighbor_table2(bbb,j_count))==1 &
                   ) then

                 gamma_up=vector_up
                 gamma_down=vector_down

 ! $B"-CV$-49$(:n6H(B
                 do k_count=1,TOTAL_DOWN_ELECTRON
                    if( gamma_down(k_count)==bbb ) then
                       gamma_down(k_count)=aaa
                       tmp=k_count
                       exit
                    end if
                 end do

! $B",CV$-49$(:n6H(B
                 do k_count=1,TOTAL_UP_ELECTRON
                    if( gamma_up(k_count)==&
                         neighbor_table2(bbb,j_count) ) then
                       gamma_up(k_count)=neighbor_table2(aaa,i_count)
                       tmp=k_count
                       exit
                    end if
                 end do

                 result2=1                         ! $B$3$l$,%_%=(B

! result3_4$B$r5a$a$k(B 
                 call choice_ket_vector(gamma_up,gamma_down,&
                      number_ket_vector_up,number_ket_vector_down,&
                      available_ket_up_table,available_ket_down_table,&
                      name,error)
                 call error_check(name,error)

                 do b2u_count=1,number_ket_vector_up
                    do k_count=1,TOTAL_UP_ELECTRON
                       right_vector_up(k_count)=available_ket_up_table(b2u_count,&
                            k_count)
                    end do

!
                    global_site_table_up=0
                    do k_count=1,TOTAL_UP_ELECTRON
                       global_site_table_up(right_vector_up(k_count))=1
                    end do

                    do b2d_count=1,number_ket_vector_down
                       if( b2u_count==number_ket_vector_up .or. &
                            b2d_count==number_ket_vector_down ) then
                          
                          do k_count=1,TOTAL_DOWN_ELECTRON
                             right_vector_down(k_count)=available_ket_down_table(&
                                  b2d_count,k_count)
                          end do
!
                          global_site_table_down=0
                          do k_count=1,TOTAL_DOWN_ELECTRON
                             global_site_table_down(right_vector_down(k_count))=1
                          end do

 ! $B%1%C%H$HA4$/F1$8$b$N$G64$s$@$H$-%/!<%m%s$,:nMQ(B
                          if( b2u_count==number_ket_vector_up .and. &  
                               b2d_count==number_ket_vector_down ) then
                                   
! $B%@%V%k%*%-%e%Q%$$N?t$rD4$Y$k!#(B
                             double_occupy=0

                             do k_count=1,TOTAL_SITE_NUMBER
                                if( global_site_table_up(k_count)==1 .and. &
                                     global_site_table_down(k_count)==1 ) then
                                   double_occupy=double_occupy+1
                                end if
                             end do

                             result3=double_occupy*COULOMB
                          else                                  ! $B$=$l0J30$G$O(Bt$B$,8z$/(B
                             result3=TRANSFER
                          end if

!
! right_vector_up$B$H(Bright_vector_down($B8GDj(B)$B$NFb@Q(B<beta_1 |phi>$B$r5a$a$k!#(B
!
                          call make_d_tilde(right_vector_up,right_vector_down,&
                               1,1,name,error)
                          call error_check(name,error)
                          call invert(1,1,name,error)
                          call error_check(name,error)
                          
                          call make_d_tilde(right_vector_up,right_vector_down,&
                               2,1,name,error)
                          call error_check(name,error)
                          call invert(2,1,name,error)
                          call error_check(name,error)
!
! $B4|BTCM7W;;$N(B<beta_1|phi>
!
                          input_up=1
                          sign=0
                          do k_count=1,TOTAL_UP_ELECTRON
                             input_up=input_up*after_lu_up(k_count,k_count)
                             if( ipiv_up(k_count)/=k_count ) then
                                sign=sign+1
                             end if
                          end do

                          input_up=input_up*((-1)**sign)

                          input_down=1
                          sign=0
                          do k_count=1,TOTAL_DOWN_ELECTRON
                             input_down=input_down*after_lu_down(k_count,k_count)
                             if( ipiv_down(k_count)/=k_count ) then
                                sign=sign+1
                             end if
                          end do

                          input_down=input_down*((-1)**sign)
        
                          lambda=jastrow(right_vector_up,right_vector_down)
                          input = input_up * input_down

                          call calc_element(input,lambda,right_vector_up,&
                               right_vector_down,result4,name,error)
                          call error_check(name,error)

! result234
                          result234=result234+result2*result3*result4

                       end if
                    end do
                 end do
              end if
           end if


!
! $BBh(B2$B9`$N7W;;(B (aaa;i_count bbb;j_count)
!
           if( aaa==neighbor_table2(bbb,j_count) .and. &
                bbb==neighbor_Table2(aaa,i_count) ) then

              if( global_site_table_up_store(aaa)==1 .and. &
                   global_site_table_down_store(bbb)==1 )then 
                 temp_b=-input
              end if
           else
              if( global_site_table_up_store(aaa)==0 .and. &
                   global_site_table_up_store&
                   &(neighbor_table2(bbb,j_count))==1 .and. &
                   global_site_table_down_store&
                   &(neighbor_table2(aaa,i_count))==0 .and. &
                   global_site_table_down_store(bbb)==1 )then 
! $B$b$70\F0@h$,6u$J$i0\F0$5$;$k(B

                 gamma_up=vector_up
                 gamma_down=vector_down
! $B",CV$-49$(:n6H(B
                 do k_count=1,TOTAL_UP_ELECTRON
                    if( gamma_up(k_count)==&
                         neighbor_table2(bbb,j_count) ) then
                       gamma_up(k_count)=aaa
                       tmp=k_count
                       exit
                    end if
                 end do

! $B"-CV$-49$(:n6H(B
                 do k_count=1,TOTAL_DOWN_ELECTRON
                    if( gamma_down(k_count)==bbb ) then
                       gamma_down(k_count)=neighbor_table2(aaa,i_count)
                       tmp=k_count
                       exit
                    end if
                 end do

                 result2=1                      ! $B$3$l$,%_%=(B

! result3_4$B$r5a$a$k(B        
                 call choice_ket_vector(gamma_up,gamma_down,&
                      number_ket_vector_up,number_ket_vector_down,&
                      available_ket_up_table,available_ket_down_table,&
                      name,error)
                 call error_check(name,error)

                 do b2u_count=1,number_ket_vector_up
                    do k_count=1,TOTAL_UP_ELECTRON
                       right_vector_up(k_count)=available_ket_up_table(b2u_count,k_count)
                    end do
!
                    global_site_table_up=0
                    do k_count=1,TOTAL_UP_ELECTRON
                       global_site_table_up(right_vector_up(k_count))=1
                    end do

                    do b2d_count=1,number_ket_vector_down
                       if( b2u_count==number_ket_vector_up .or. &
                            b2d_count==number_ket_vector_down ) then

                          do k_count=1,TOTAL_DOWN_ELECTRON
                             right_vector_down(k_count)=available_ket_down_table(b2d_count,k_count)
                          end do
!
                          global_site_table_down=0
                          do k_count=1,TOTAL_DOWN_ELECTRON
                             global_site_table_down(right_vector_down(k_count))=1
                          end do

 ! $B%1%C%H$HA4$/F1$8$b$N$G64$s$@$H$-%/!<%m%s$,:nMQ(B
                          if( b2u_count==number_ket_vector_up .and. &  
                               b2d_count==number_ket_vector_down ) then

! $B%@%V%k%*%-%e%Q%$$N?t$rD4$Y$k!#(B
                             double_occupy=0

                             do k_count=1,TOTAL_SITE_NUMBER
                                if( global_site_table_up(k_count)==1 .and. &
                                     global_site_table_down(k_count)==1 ) then
                                   double_occupy=double_occupy+1
                                end if
                             end do

                             result3=double_occupy*COULOMB
                          else                                  ! $B$=$l0J30$G$O(Bt$B$,8z$/(B
                             result3=TRANSFER
                          end if
!
! right_vector_up$B$H(Bright_vector_down($B8GDj(B)$B$NFb@Q(B<beta_1 |phi>$B$r5a$a$k!#(B
!
                          call make_d_tilde(right_vector_up,right_vector_down,&
                               1,1,name,error)
                          call error_check(name,error)
                          call invert(1,1,name,error)
                          call error_check(name,error)

                          call make_d_tilde(right_vector_up,right_vector_down,&
                               2,1,name,error)
                          call error_check(name,error)
                          call invert(2,1,name,error)
                          call error_check(name,error)
!
! $B4|BTCM7W;;$N(B<beta_1|phi>
!
                          input_up=1
                          sign=0
                          do k_count=1,TOTAL_UP_ELECTRON
                             input_up=input_up*after_lu_up(k_count,k_count)
                             if( ipiv_up(k_count)/=k_count ) then
                                sign=sign+1
                             end if
                          end do

                          input_up=input_up*((-1)**sign)
           
                          input_down=1
                          sign=0
                          do k_count=1,TOTAL_DOWN_ELECTRON
                             input_down=input_down*after_lu_down(k_count,k_count)
                             if( ipiv_down(k_count)/=k_count ) then
                                sign=sign+1
                             end if
                          end do

                          input_down=input_down*((-1)**sign)
                          
                          lambda=jastrow(right_vector_up,right_vector_down)
                          input = input_up * input_down

                          call calc_element(input,lambda,right_vector_up,&
                               right_vector_down,result4,name,error)
                          call error_check(name,error)
! result234
                          result234=result234+result2*result3*result4

                       end if
                    end do
                 end do
              end if
           end if

!!
!! $BBh(B3$B9`$N7W;;(B (aaa;i_count bbb;j_count)
!!

! $B%J%s%P!<%*%Z%l!<%?!<$N>l9g(B
           if( aaa==bbb .and. &
                neighbor_table2(aaa,i_count)== &
                neighbor_table2(bbb,j_count) ) then
              
              if( global_site_table_up_store(aaa)==1 .and. &
                   global_site_table_down_store&
                   &(neighbor_table2(aaa,i_count))==1 ) then
                 temp_a=-input

              end if
           else
              if( global_site_table_up_store(aaa)==0 .and. &
                   global_site_table_up_store(bbb)==1 .and. &
                   global_site_table_down_store&
                   &(neighbor_table2(aaa,i_count))==0 .and. &
                   global_site_table_down_store&
                   &(neighbor_table2(bbb,j_count))==1 &
                   ) then

                 gamma_up=vector_up
                 gamma_down=vector_down

 ! $B",CV$-49$(:n6H(B
                 do k_count=1,TOTAL_UP_ELECTRON
                    if( gamma_up(k_count)==bbb ) then
                       gamma_up(k_count)=aaa
                       tmp=k_count
                       exit
                    end if
                 end do

! $B"-CV$-49$(:n6H(B
                 do k_count=1,TOTAL_DOWN_ELECTRON
                    if( gamma_down(k_count)==&
                         neighbor_table2(bbb,j_count) ) then
                       gamma_down(k_count)=neighbor_table2(aaa,i_count)
                       tmp=k_count
                       exit
                    end if
                 end do

                 result2=1                         ! $B$3$l$,%_%=(B

! result3_4$B$r5a$a$k(B 
                 call choice_ket_vector(gamma_up,gamma_down,&
                      number_ket_vector_up,number_ket_vector_down,&
                      available_ket_up_table,available_ket_down_table,&
                      name,error)
                 call error_check(name,error)

                 do b2u_count=1,number_ket_vector_up
                    do k_count=1,TOTAL_UP_ELECTRON
                       right_vector_up(k_count)=available_ket_up_table(b2u_count,&
                            k_count)
                    end do
!
                    global_site_table_up=0
                    do k_count=1,TOTAL_UP_ELECTRON
                       global_site_table_up(right_vector_up(k_count))=1
                    end do

                    do b2d_count=1,number_ket_vector_down
                       if( b2u_count==number_ket_vector_up .or. &
                            b2d_count==number_ket_vector_down ) then
                          
                          do k_count=1,TOTAL_DOWN_ELECTRON
                             right_vector_down(k_count)=available_ket_down_table(&
                                  b2d_count,k_count)
                          end do
!
                          global_site_table_down=0
                          do k_count=1,TOTAL_DOWN_ELECTRON
                             global_site_table_down(right_vector_down(k_count))=1
                          end do

 ! $B%1%C%H$HA4$/F1$8$b$N$G64$s$@$H$-%/!<%m%s$,:nMQ(B
                          if( b2u_count==number_ket_vector_up .and. &  
                               b2d_count==number_ket_vector_down ) then
                                                       
! $B%@%V%k%*%-%e%Q%$$N?t$rD4$Y$k!#(B
                             double_occupy=0

                             do k_count=1,TOTAL_SITE_NUMBER
                                if( global_site_table_up(k_count)==1 .and. &
                                     global_site_table_down(k_count)==1 ) then
                                   double_occupy=double_occupy+1
                                end if
                             end do

                             result3=double_occupy*COULOMB
                          else                                  ! $B$=$l0J30$G$O(Bt$B$,8z$/(B
                             result3=TRANSFER
                          end if
!
! right_vector_up$B$H(Bright_vector_down($B8GDj(B)$B$NFb@Q(B<beta_1 |phi>$B$r5a$a$k!#(B
!
                          call make_d_tilde(right_vector_up,right_vector_down,&
                               1,1,name,error)
                          call error_check(name,error)
                          call invert(1,1,name,error)
                          call error_check(name,error)
                          
                          call make_d_tilde(right_vector_up,right_vector_down,&
                               2,1,name,error)
                          call error_check(name,error)
                          call invert(2,1,name,error)
                          call error_check(name,error)
!
! $B4|BTCM7W;;$N(B<beta_1|phi>
!
                          input_up=1
                          sign=0
                          do k_count=1,TOTAL_UP_ELECTRON
                             input_up=input_up*after_lu_up(k_count,k_count)
                             if( ipiv_up(k_count)/=k_count ) then
                                sign=sign+1
                             end if
                          end do

                          input_up=input_up*((-1)**sign)

                          input_down=1
                          sign=0
                          do k_count=1,TOTAL_DOWN_ELECTRON
                             input_down=input_down*after_lu_down(k_count,k_count)
                             if( ipiv_down(k_count)/=k_count ) then
                                sign=sign+1
                             end if
                          end do

                          input_down=input_down*((-1)**sign)
        
                          lambda=jastrow(right_vector_up,right_vector_down)
                          input = input_up * input_down

                          call calc_element(input,lambda,right_vector_up,&
                               right_vector_down,result4,name,error)
                          call error_check(name,error)

! result234
                          result234=result234+result2*result3*result4

                       end if
                    end do
                 end do
              end if
           end if

!
! $BBh(B4$B9`$N7W;;(B (aaa;i_count bbb;j_count)
!
           if( aaa==neighbor_table2(bbb,j_count) .and. &
                bbb==neighbor_Table2(aaa,i_count) ) then

              if( global_site_table_up_store(aaa)==1 .and. &
                   global_site_table_down_store(bbb)==1 )then 
                 temp_b=-input

              end if
           else
              if( global_site_table_up_store(aaa)==0 .and. &
                   global_site_table_up_store&
                   &(neighbor_table2(bbb,j_count))==1 .and. &
                   global_site_table_down_store&
                   &(neighbor_table2(aaa,i_count))==0 .and. &
                   global_site_table_down_store(bbb)==1 )then 
! $B$b$70\F0@h$,6u$J$i0\F0$5$;$k(B

                 gamma_up=vector_up
                 gamma_down=vector_down
! $B",CV$-49$(:n6H(B
                 do k_count=1,TOTAL_UP_ELECTRON
                    if( gamma_up(k_count)==&
                         neighbor_table2(bbb,j_count) ) then
                       gamma_up(k_count)=aaa
                       tmp=k_count
                       exit
                    end if
                 end do

                 call make_d_tilde(gamma_up,gamma_down,1,1,name,error)
                 call error_check(name,error)
! $B"-CV$-49$(:n6H(B
                 do k_count=1,TOTAL_DOWN_ELECTRON
                    if( gamma_down(k_count)==bbb ) then
                       gamma_down(k_count)=neighbor_table2(aaa,i_count)
                       tmp=k_count
                       exit
                    end if
                 end do

                 call make_d_tilde(gamma_up,gamma_down,2,1,name,error)
                 call error_check(name,error)

                 result2=1                      ! $B$3$l$,%_%=(B

! result3_4$B$r5a$a$k(B        
                 call choice_ket_vector(gamma_up,gamma_down,&
                      number_ket_vector_up,number_ket_vector_down,&
                      available_ket_up_table,available_ket_down_table,&
                      name,error)
                 call error_check(name,error)

                 do b2u_count=1,number_ket_vector_up
                    do k_count=1,TOTAL_UP_ELECTRON
                       right_vector_up(k_count)=available_ket_up_table(b2u_count,k_count)
                    end do
!
                    global_site_table_up=0
                    do k_count=1,TOTAL_UP_ELECTRON
                       global_site_table_up(right_vector_up(k_count))=1
                    end do

                    do b2d_count=1,number_ket_vector_down
                       if( b2u_count==number_ket_vector_up .or. &
                            b2d_count==number_ket_vector_down ) then

                          do k_count=1,TOTAL_DOWN_ELECTRON
                             right_vector_down(k_count)=available_ket_down_table(b2d_count,k_count)
                          end do
!
                          global_site_table_down=0
                          do k_count=1,TOTAL_DOWN_ELECTRON
                             global_site_table_down(right_vector_down(k_count))=1
                          end do

 ! $B%1%C%H$HA4$/F1$8$b$N$G64$s$@$H$-%/!<%m%s$,:nMQ(B
                          if( b2u_count==number_ket_vector_up .and. &  
                               b2d_count==number_ket_vector_down ) then

! $B%@%V%k%*%-%e%Q%$$N?t$rD4$Y$k!#(B
                             double_occupy=0

                             do k_count=1,TOTAL_SITE_NUMBER
                                if( global_site_table_up(k_count)==1 .and. &
                                     global_site_table_down(k_count)==1 ) then
                                   double_occupy=double_occupy+1
                                end if
                             end do

                             result3=double_occupy*COULOMB
                          else                                  ! $B$=$l0J30$G$O(Bt$B$,8z$/(B
                             result3=TRANSFER
                          end if
!
! right_vector_up$B$H(Bright_vector_down($B8GDj(B)$B$NFb@Q(B<beta_1 |phi>$B$r5a$a$k!#(B
!
                          call make_d_tilde(right_vector_up,right_vector_down,&
                               1,1,name,error)
                          call error_check(name,error)
                          call invert(1,1,name,error)
                          call error_check(name,error)

                          call make_d_tilde(right_vector_up,right_vector_down,&
                               2,1,name,error)
                          call error_check(name,error)
                          call invert(2,1,name,error)
                          call error_check(name,error)
!
! $B4|BTCM7W;;$N(B<beta_1|phi>
!
                          input_up=1
                          sign=0
                          do k_count=1,TOTAL_UP_ELECTRON
                             input_up=input_up*after_lu_up(k_count,k_count)
                             if( ipiv_up(k_count)/=k_count ) then
                                sign=sign+1
                             end if
                          end do
                          
                          input_up=input_up*((-1)**sign)
           
                          input_down=1
                          sign=0
                          do k_count=1,TOTAL_DOWN_ELECTRON
                             input_down=input_down*after_lu_down(k_count,k_count)
                             if( ipiv_down(k_count)/=k_count ) then
                                sign=sign+1
                             end if
                          end do

                          input_down=input_down*((-1)**sign)
        
                          lambda=jastrow(right_vector_up,right_vector_down)
                          input = input_up * input_down

                          call calc_element(input,lambda,right_vector_up,&
                               right_vector_down,result4,name,error)
                          call error_check(name,error)
! result234
                          result234=result234+result2*result3*result4

                       end if
                    end do
                 end do
              end if
           end if

        end if
     end do
  end do

  result234=0.5*result234

  global_site_table_up=global_site_table_up_store
  global_site_table_down=global_site_table_down_store
  d_tilde_up_inverse=d_tilde_up_inverse_store
  d_tilde_down_inverse=d_tilde_down_inverse_store

end subroutine super_nume5

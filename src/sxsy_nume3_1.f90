!
! $Id: sxsy_nume3_1.f90,v 1.4 2003/11/20 14:04:26 k-yukino Exp $
!
#include "parameter.h"

subroutine sxsy_nume3_1(alpha_h_psi1,vector_up,vector_down,aaa,bbb,result23,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(in) :: alpha_h_psi1
  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  integer,intent(in) :: aaa,bbb
  real(8),intent(out) :: result23
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer,dimension(TOTAL_SITE_NUMBER) :: global_site_table_up_in&
       ,global_site_table_down_in
  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  integer :: i_count
  real(8),external :: jastrow
  real(8) :: temp_a,temp_b
  real(8) :: result23_1,result23_2
  integer :: tmp
  real(8),dimension(TOTAL_UP_ELECTRON,TOTAL_UP_ELECTRON) :: d_tilde_up_inverse_store
  real(8),dimension(TOTAL_DOWN_ELECTRON,TOTAL_DOWN_ELECTRON) :: d_tilde_down_inverse_store
  real(8) :: lambda
  real(8) :: input_up,input_down,input
  integer :: sign
! init
  name="sxsy_nume3_1" ; error=0
! init local
  temp_a=0 ; temp_b=0
  tmp=0
  result23=0 ; result23_1=0 ; result23_2=0

  global_site_table_up_in=global_site_table_up
  global_site_table_down_in=global_site_table_down
  d_tilde_up_inverse_store=d_tilde_up_inverse
  d_tilde_down_inverse_store=d_tilde_down_inverse

!
! <psi|H|alpha><alpha|O|beta_1><beta_1|H|psi>
!               ^^^^^^^^^^^^^   ^^^^^^^^^^^^^^
!                 result2         result3     
!


!
! result23$B$r5a$a$k(B
!
! <alpha|$B$H7k$SIU$/(B|beta_1>$B$G7W;;$r9T$J$&!#(B
!
  if( aaa==bbb ) then                            ! $B%J%s%P!<%*%Z%l!<%?!<(B
! $BBh(B1$B9`(B
     if( global_site_table_up_in(aaa)==1 ) then
        temp_a=1
        if( global_site_table_down_in(aaa)==1 ) then
           temp_a=0
        end if
     end if

! $BBh(B2$B9`(B
     if( global_site_table_down_in(aaa)==1 ) then
        temp_b=1
        if( global_site_table_up_in(aaa)==1 ) then
           temp_b=0
        end if
     end if
     
     result23=(temp_a+temp_b)*alpha_h_psi1
  else                                           ! aaa/=bbb
! $BBh(B1$B9`$N7W;;(B (aaa;i_count bbb;j_count)
     if( global_site_table_up_in(aaa)==0 .and. &
          global_site_table_up_in(bbb)==1 .and. &
          global_site_table_down_in(bbb)==0 .and. &
          global_site_table_down_in(aaa)==1 ) then 
! $B$b$70\F0@h$,6u$J$i0\F0$5$;$k(B

! $B:n6HNN0h$K%3%T!<(B
        gamma_up=vector_up
        gamma_down=vector_down

! $B",CV$-49$(:n6H(B
        do i_count=1,TOTAL_UP_ELECTRON
           if( gamma_up(i_count)==bbb ) then
              gamma_up(i_count)=aaa
              tmp=i_count
              exit
           end if
        end do
        
! $B"-CV$-49$(:n6H(B

        do i_count=1,TOTAL_DOWN_ELECTRON
           if( gamma_down(i_count)==aaa ) then
              gamma_down(i_count)=bbb
              tmp=i_count
              exit
           end if
        end do

        call make_d_tilde(gamma_up,gamma_down,1,1,name,error) 
        call error_check(name,error)
        call invert(1,1,name,error)
        call error_check(name,error)

        call make_d_tilde(gamma_up,gamma_down,2,1,name,error)
        call error_check(name,error)
        call invert(2,1,name,error)
        call error_check(name,error)

!
! $BMWAG7W;;$HFb@Q(B
!
        input_up=1
        sign=0
        do i_count=1,TOTAL_UP_ELECTRON
           input_up=input_up*after_lu_up(i_count,i_count)
           if( ipiv_up(i_count)/=i_count ) then
              sign=sign+1
           end if
        end do

        input_up=input_up*((-1)**sign)

        input_down=1
        sign=0
        do i_count=1,TOTAL_DOWN_ELECTRON
           input_down=input_down*after_lu_down(i_count,i_count)
           if( ipiv_down(i_count)/=i_count ) then
              sign=sign+1
           end if
        end do

        input_down=input_down*((-1)**sign)

        input=input_up*input_down

        lambda=jastrow(gamma_up,gamma_down)

! (global_site_table_up,down$B$K%;%C%H$7$F$*$/(B)
        global_site_table_up=0
        do i_count=1,TOTAL_UP_ELECTRON
           global_site_table_up(gamma_up(i_count))=1
        end do

        global_site_table_down=0
        do i_count=1,TOTAL_DOWN_ELECTRON
           global_site_table_down(gamma_down(i_count))=1
        end do
        
        call calc_element(input,lambda,gamma_up,&
             gamma_down,result23_1,name,error)
        call error_check(name,error)
     else
        result23_1=0
     end if

!
! $BBh(B2$B9`$N7W;;(B (aaa;i_count bbb;j_count)
!
     if( global_site_table_up_in(bbb)==0 .and. &
          global_site_table_up_in(aaa)==1 .and. &
          global_site_table_down_in(aaa)==0 .and. &
          global_site_table_down_in(bbb)==1 )then   
! $B$b$70\F0@h$,6u$J$i0\F0$5$;$k(B

! $B:n6HNN0h$K%3%T!<(B        
        gamma_up=vector_up
        gamma_down=vector_down

! $B",CV$-49$(:n6H(B
        do i_count=1,TOTAL_UP_ELECTRON
           if( gamma_up(i_count)==aaa ) then
              gamma_up(i_count)=bbb
              tmp=i_count
              exit
           end if
        end do

! $B"-CV$-49$(:n6H(B
        do i_count=1,TOTAL_DOWN_ELECTRON
           if( gamma_down(i_count)==bbb ) then
              gamma_down(i_count)=aaa
              tmp=i_count
              exit
           end if
        end do

        call make_d_tilde(gamma_up,gamma_down,1,1,name,error)
        call error_check(name,error)
        call invert(1,1,name,error)
        call error_check(name,error)

        call make_d_tilde(gamma_up,gamma_down,2,1,name,error)
        call error_check(name,error)
        call invert(2,1,name,error)
        call error_check(name,error)
!
! $BMWAG7W;;$HFb@Q(B
!
        input_up=1
        sign=0
        do i_count=1,TOTAL_UP_ELECTRON
           input_up=input_up*after_lu_up(i_count,i_count)
           if( ipiv_up(i_count)/=i_count ) then
              sign=sign+1
           end if
        end do

        input_up=input_up*((-1)**sign)

        input_down=1
        sign=0
        do i_count=1,TOTAL_DOWN_ELECTRON
           input_down=input_down*after_lu_down(i_count,i_count)
           if( ipiv_down(i_count)/=i_count ) then
              sign=sign+1
           end if
        end do

        input_down=input_down*((-1)**sign)

        input=input_up*input_down
        
        lambda=jastrow(gamma_up,gamma_down)

! (global_site_table_up,down$B$K%;%C%H$7$F$*$/(B)
        global_site_table_up=0
        do i_count=1,TOTAL_UP_ELECTRON
           global_site_table_up(gamma_up(i_count))=1
        end do

        global_site_table_down=0
        do i_count=1,TOTAL_DOWN_ELECTRON
           global_site_table_down(gamma_down(i_count))=1
        end do

        call calc_element(input,lambda,gamma_up,&
             gamma_down,result23_2,name,error)
        call error_check(name,error)
     else
        result23_2=0
     end if

     result23=result23_1+result23_2
  end if

  result23=0.5*result23

  global_site_table_up=global_site_table_up_in
  global_site_table_down=global_site_table_down_in
  d_tilde_up_inverse=d_tilde_up_inverse_store
  d_tilde_down_inverse=d_tilde_down_inverse_store

end subroutine sxsy_nume3_1

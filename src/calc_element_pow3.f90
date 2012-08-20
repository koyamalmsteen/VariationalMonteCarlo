!
! $Id: calc_element_pow3.f90,v 1.6 2004/03/10 03:04:02 k-yukino Exp k-yukino $
!
!
! E=<phi_0 |H^3|phi_0 >
!
#include "parameter.h"

subroutine calc_element_pow3(vector_up,vector_down,alpha_h_h_psi1,&
     name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  real(8),intent(out) :: alpha_h_h_psi1
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: double_occupy
  real(8),external :: jastrow,lambda_field
  integer :: i_count,b3u_count,b3d_count
  integer,dimension(TOTAL_UP_ELECTRON) :: right_vector_up3
  integer,dimension(TOTAL_DOWN_ELECTRON) :: right_vector_down3
  real(8) :: result3,result4_5
  integer :: number_ket_vector_up,number_ket_vector_down,sign
  integer,dimension(2*JIGEN*TOTAL_UP_ELECTRON+1,TOTAL_UP_ELECTRON) :: available_ket_up_table
  integer,dimension(2*JIGEN*TOTAL_DOWN_ELECTRON+1,TOTAL_DOWN_ELECTRON) :: available_ket_down_table
  real(8) :: input,input_up,input_down,tmp_jastrow
  real(8),external :: calc_total_pole
  real(8) :: pole
! init
  name="calc_element_pow3" ; error=0
  available_ket_up_table=0 ; available_ket_down_table=0
! init [local]
  number_ket_vector_up=0 ; number_ket_vector_down=0
  result3=0
  alpha_h_h_psi1=0 ; result4_5=0
  input=0 ; input_up=0 ; input_down=0  ; sign=0 ; tmp_jastrow=0
  pole=0

!
! <phi_0|H^3|phi_0>/<phi_0|phi_0> -->  
! <phi_0 |beta_1><beta_1|H|alpha><alpha|H|beta_2><beta_2|H|beta3><beta3|phi_0>
!
!                                ^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^ ^^^^^^^^^^^ 
!                                    result3         result4         result5

! <alpha|$B$H7k$SIU$/(B|beta_2>$B$r5a$a$k!#(B

  call choice_ket_vector(vector_up,vector_down,&
       number_ket_vector_up,number_ket_vector_down,&
       available_ket_up_table,available_ket_down_table,&
       name,error)
  call error_check(name,error)

  do b3u_count=1,number_ket_vector_up         ! $B",$r8GDj(B
     do i_count=1,TOTAL_UP_ELECTRON
        right_vector_up3(i_count)=available_ket_up_table(b3u_count,i_count)
     end do

! calc_element$B$G;H$&$N$G(Bglobal_site_table_up$B$K=q$-JQ$((B
     global_site_table_up=0
     do i_count=1,TOTAL_UP_ELECTRON
        global_site_table_up(right_vector_up3(i_count))=1
     end do

! for (result4_5)
! <right_vector_up3|$B$H3N<B$K7k$SIU$/(B|right_vector_up3>$B$N5U9TNs$r7W;;$7$F$*$/(B
! ($B9bB.2=(B)
     call make_d_tilde(right_vector_up3,right_vector_down3,1,1,name,error)
     call error_check(name,error)
     call invert(1,1,name,error)
     call error_check(name,error)

! $BFb@Q(B(up$BJ,(B)
     input_up=1
     sign=0
     do i_count=1,TOTAL_UP_ELECTRON
        input_up=input_up*after_lu_up(i_count,i_count)
        if( ipiv_up(i_count)/=i_count ) then
           sign=sign+1
        end if
     end do

     input_up=input_up*((-1)**sign)

     do b3d_count=1,number_ket_vector_down           ! $B"-$r8GDj(B
        ! $B%O%_%k%H%K%"%s$N7A$h$j!"$I$A$i$+$N%9%T%s$,(Bfix$B$5$l$F$$$kI,MW$,$"$k$N$G(B
        if( b3u_count==number_ket_vector_up .or. &
             b3d_count==number_ket_vector_down ) then     

           do i_count=1,TOTAL_DOWN_ELECTRON
              right_vector_down3(i_count)=available_ket_down_table(&
                   b3d_count,i_count)
           end do

! calc_element$B$G;H$&$N$G(Bglobal_site_table_down$B$K=q$-JQ$((B
           global_site_table_down=0
           do i_count=1,TOTAL_DOWN_ELECTRON
              global_site_table_down(right_vector_down3(i_count))=1
           end do

! for (result4_5)
! <right_vector_up3|$B$H3N<B$K7k$SIU$/(B|right_vector_up3>$B$N5U9TNs$r7W;;$7$F$*$/(B
! ($B9bB.2=(B)
           call make_d_tilde(right_vector_up3,right_vector_down3,2,1,&
                name,error)
           call error_check(name,error)
           call invert(2,1,name,error)
           call error_check(name,error)

! $BFb@Q(B(down$BJ,(B)
           input_down=1
           sign=0
           do i_count=1,TOTAL_DOWN_ELECTRON
              input_down=input_down*after_lu_down(i_count,i_count)
              if( ipiv_down(i_count)/=i_count ) then
                 sign=sign+1
              end if
           end do

           input_down=input_down*((-1)**sign)
!
! result3$B$r5a$a$k(B
! 
           if( b3u_count==number_ket_vector_up .and. &   ! $B%1%C%H$HA4$/F1$8$b$N(B
                b3d_count==number_ket_vector_down ) then ! $B$G$O$5$s$@$H$-(B
                                                         ! $B%/!<%m%s$,:nMQ(B
! $B%@%V%k%*%-%e%Q%$$N?t$rD4$Y$k!#(B
              double_occupy=0

              do i_count=1,TOTAL_SITE_NUMBER
                 if( global_site_table_up(i_count)==1 .and. &
                      global_site_table_down(i_count)==1 ) then
                    double_occupy=double_occupy+1
                 end if
              end do

              result3=double_occupy*dble(COULOMB)

! electric_field
              if( electric_field/=0 .and. PERIODIC==1 .and. xi_field/=0 ) then
                 pole=calc_total_pole(right_vector_up3,right_vector_down3)
                 result3=result3+electric_field*pole
              end if
           else                                     ! $B$=$l0J30$G$O(Bt$B$,8z$/(B
              result3=dble(TRANSFER)
           end if
!
! $BEE3&$N9`(B
!
           tmp_jastrow=jastrow(right_vector_up3,right_vector_down3)
           if( electric_field/=0 .and. PERIODIC==1 .and. xi_field/=0 ) then
              tmp_jastrow=tmp_jastrow*lambda_field(right_vector_up3,right_vector_down3)
           end if

           input = input_up * input_down
!
! result4_5
!
           result4_5=0

           call calc_element(input,tmp_jastrow,right_vector_up3,&
                right_vector_down3,result4_5,name,error)
           call error_check(name,error)

! alpha_h_h_psi1$B$O(Bpow4$B$N7k2L$K;H$($k$+$i;H$&(B

           alpha_h_h_psi1=alpha_h_h_psi1+result3*result4_5
        end if
     end do
  end do

! available_ket_up_table$B$r5U=g$KF~$l$F$k$+$i!"(Bglobal_site_table$B$O(B
! $B$3$N%k!<%A%s$KF~$C$F$-$?$H$-$HF1$8$K$J$C$F$$$k!#(B

end subroutine calc_element_pow3



#ifdef DEBUG_CALC_ELEMENT_POW3

program main
  use global_variables
  implicit none

  integer,dimension(TOTAL_UP_ELECTRON) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: vector_down
  real(8) :: alpha_h_h_psi1,tmp_jastrow,input_up,input_down,alpha_psi0
  real(8),external :: jastrow
  character(32) :: name
  integer :: error
! local
  real(8) :: zero_approx_energy
  integer :: hf_iteration_result
  integer :: i_count
  integer,dimension(1) :: ic
! init
  alpha_h_h_psi1=0 ; tmp_jastrow=0 ; error=0 ; name="main"
  vector_up=0 ; vector_down=0
  input_up=0 ; input_down=0 ; alpha_psi0=0
  ic(1)=1
  zero_approx_energy=0 ; hf_iteration_result=0

  call MPI_INIT(IERROR)

! $B%Q%i%a!<%?!<$N=i4|2=(B
  call init_global_variables(name,error)
  call error_check(name,error)

! $BMp?t%F!<%V%k$N=i4|2=(B
  call random_init(ic,name,error)
  call error_check(name,error)

! $BBh(B0$B6a;w$N<B9T(B(INIT_WAVE_VECTOR = 0;Hartree-Fock , 1 ;Plane wave )
! (MYRANK==0$B$N$_$,<B9T$7!"7k2L$r(Bbcast$B$9$k(B)
  call zero_approx(zero_approx_energy,hf_iteration_result,name,error)
  call error_check(name,error)

  do i_count=1,TOTAL_UP_ELECTRON
     write(*,*) "input vector_up" ; read(*,*) vector_up(i_count)
  end do

  do i_count=1,TOTAL_DOWN_ELECTRON
     write(*,*) "input vector_down" ; read(*,*) vector_down(i_count)
  end do

  call make_d_tilde(vector_up,vector_down,1,1,name,error)
  call error_check(name,error)
  call make_d_tilde(vector_up,vector_down,2,1,name,error)
  call error_check(name,error)

  write(*,*) "d_tilde_up=",d_tilde_up
  write(*,*) "d_tilde_down=",d_tilde_down
  
! $B5U9TNs$r5a$a$k(B
  call invert(1,1,name,error)
  call invert(2,1,name,error)

  tmp_jastrow=jastrow(vector_up,vector_down)

  global_site_table_up=0
  do i_count=1,TOTAL_UP_ELECTRON
     global_site_table_up(vector_up(i_count))=1
  end do

  global_site_table_down=0
  do i_count=1,TOTAL_DOWN_ELECTRON
     global_site_table_down(vector_down(i_count))=1
  end do

! <alpha|psi0>$B$r5a$a$k(B
  input_up=1
  do i_count=1,TOTAL_UP_ELECTRON
     input_up=input_up*d_tilde_up(i_count,i_count)
  end do

  input_down=1
  do i_count=1,TOTAL_DOWN_ELECTRON
     input_down=input_down*d_tilde_down(i_count,i_count)
  end do

  alpha_psi0 = input_up * input_down

  write(*,*) "alpha_psi0=",alpha_psi0
  write(*,*) "tmp_jastrow=",tmp_jastrow
  write(*,*) "vector_up=",vector_up
  write(*,*) "vector_down=",vector_down
  write(*,*) 
  write(*,*) "COULOMB=",COULOMB,"TRANSFER=",TRANSFER
  write(*,*)

  call calc_element_pow3(vector_up,vector_down,alpha_h_h_psi1,name,error)

  write(*,*) "alpha_h_h_psi1=",alpha_h_h_psi1

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CALC_ELEMENT_POW3 */ 

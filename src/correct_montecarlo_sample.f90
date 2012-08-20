!
! $Id: correct_montecarlo_sample.f90,v 1.5 2004/02/15 10:52:35 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine correct_montecarlo_sample(sample_number,vector_up,vector_down,&
     alpha_psi0,alpha_psi1,lambda,name,error)
  use global_variables
  implicit none

  integer,intent(in) :: sample_number
  integer,dimension(TOTAL_UP_ELECTRON),intent(inout) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(inout) :: vector_down
  real(8),intent(out) :: alpha_psi0,alpha_psi1,lambda
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8) :: input_up,input_down
  integer :: i_count,j_count
  integer :: idou_spin,idou_saki,idou_moto
  integer :: sign
  real(8) :: ransuu
  real(8) :: denomi_lambda
  real(8),external :: jastrow,lambda_field
  real(8) :: acceptance_ratio
  integer,dimension(TOTAL_UP_ELECTRON) :: vector_up_work 
  integer,dimension(TOTAL_DOWN_ELECTRON) :: vector_down_work 
  real(8) :: tmp_pole
! init
  alpha_psi0=0 ; alpha_psi1=0 ; lambda=0
  name="correct_montecarlo_sample" ; error=0
! local init
  sign=0
  tmp_pole=0
  acceptance_ratio=0
  idou_spin=0 ; idou_saki=0 ; idou_moto=0
  ransuu=0
  denomi_lambda=0
  input_up=0 ; input_down=0

! $B0lHV=i$a$N(Bgamma$B$rMp?t$h$j7hDj$9$k(B(sort$B$bI,?\(B)
  if( sample_number==1 ) then
     vector_up=0 ; vector_down=0                ! $B:G=i$@$1=i4|2=(B

     if( PROJECTION==0 ) then
        call choice_gamma_sigma(TOTAL_UP_ELECTRON,vector_up,name,error)
        call error_check(name,error)
        call choice_gamma_sigma(TOTAL_DOWN_ELECTRON,vector_down,name,error)
        call error_check(name,error)
     else                  /* $B%W%m%8%'%/%7%g%s%*%Z%l!<%?!<IU$-(B */
        call choice_gamma_sigma_projection(TOTAL_UP_ELECTRON,vector_up,&
             name,error)
        call error_check(name,error)
        call choice_gamma_sigma_projection(TOTAL_DOWN_ELECTRON,&
             vector_down,name,error)
        call error_check(name,error)
     end if

! $B5U9TNs$r:n$k(B
     call make_d_tilde(vector_up,vector_down,1,1,name,error)
     call error_check(name,error)
     call make_d_tilde(vector_up,vector_down,2,1,name,error)
     call error_check(name,error)

     call invert(1,1,name,error)
     call error_check(name,error)
     call invert(2,1,name,error)
     call error_check(name,error)

! $B:G=i$N(BTRASH$B8D$N>uBV$r<N$F$k(B
     call trash(vector_up,vector_down,name,error)
     call error_check(name,error)

! $B=i4|$N;n9TEE;R>uBV(B
! jastrow$B4X?t$N7W;;$r9bB.2=$5$;$k$?$a!"0J2<$N$h$&$K%5%$%HI=<($K=q$-49$($k(B
! ($B=q$-49$(:n6H$r9M$($F$b9bB.2=$5$l$k(B)
     global_site_table_up=0
     do i_count=1,TOTAL_UP_ELECTRON
        global_site_table_up(vector_up(i_count))=1
     end do
     global_site_table_down=0
     do i_count=1,TOTAL_DOWN_ELECTRON
        global_site_table_down(vector_down(i_count))=1
     end do

!
! $B9bB.2=$N$?$a!":G=i$NEE;R>uBV$G$N%(%M%k%.!<4|BTCM7W;;$NJ,Jl$N(B
! <alpha|psi>=lambda_alpha<alpha|HF>$B$r$7$F$*$/!#(B
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
     stored_input_up=input_up                     ! $BJ]B8(B

     input_down=1
     sign=0
     do i_count=1,TOTAL_DOWN_ELECTRON
        input_down=input_down*after_lu_down(i_count,i_count)
        if( ipiv_down(i_count)/=i_count ) then
           sign=sign+1
        end if
     end do

     input_down=input_down*((-1)**sign)
     stored_input_down=input_down                 ! $BJ]B8(B

     lambda=jastrow(vector_up,vector_down)

     if( electric_field/=0 .and. PERIODIC==1 .and. xi_field/=0 ) then 
        lambda=lambda*lambda_field(vector_up,vector_down) 
     end if

     stored_lambda=lambda

     alpha_psi0 = input_up * input_down
     alpha_psi1 = input_up * input_down * lambda

     stored_alpha_psi0=alpha_psi0
     stored_alpha_psi1=alpha_psi1

/*
     state_vector(sample_number)%vector_up=vector_up
     state_vector(sample_number)%vector_down=vector_down
     state_vector(sample_number)%alpha_psi0=alpha_psi0
     state_vector(sample_number)%alpha_psi1=alpha_psi1
     state_vector(sample_number)%lambda=lambda
*/

     return
  end if

!
! $B%5%s%W%kHV9f(B2$BHVL\0J9_(B
!
  if( sample_number<=MAX_MONTECARLO_SAMPLE ) then
     do j_count=2,TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON

! $B?7$7$$EE;R>uBV(B($B$N8uJd(B)$B$N(Balpha'$B$r:n$k!#(B
! ($B8E$$EE;R>uBV$O(Bvector_up,vector_down$B$H$7$F;D$7$F$*$/(B($B%j%8%'%/%H;~$K;H$&$+$i(B))

! $B%o!<%/MQ$K0\$9(B
        vector_up_work=vector_up
        vector_down_work=vector_down

        if( PROJECTION==0 ) then 
           call choice_new_gamma(vector_up_work,vector_down_work,&
                idou_spin,idou_moto,idou_saki,name,error)
           call error_check(name,error)
        else
           call choice_new_gamma_projection(vector_up_work,vector_down_work,&
                idou_spin,idou_moto,idou_saki,name,error)
           call error_check(name,error)
        end if

! $B?7$7$$EE;R>uBV(B($B$N8uJd(B)$B$N(Balpha'$B$N9TNs$r:n@.(B
        if( idou_spin==1 ) then  ! $B",%9%T%s$,F0$$$?$i(B
           call make_d_tilde(vector_up_work,vector_down_work,1,1,name,error)
           call error_check(name,error)
        else if( idou_spin==2 ) then ! $B"-%9%T%s$,F0$$$?$i(B 
           call make_d_tilde(vector_up_work,vector_down_work,2,1,name,error)
           call error_check(name,error)
        end if
     
! $B%"%/%;%W%?%s%9!&%l%7%*$N7W;;(B
        call calc_acceptance_ratio(vector_up_work,vector_down_work,&
             stored_lambda,denomi_lambda,idou_spin,idou_moto,&
             acceptance_ratio,name,error)
        call error_check(name,error)
     
! $B%"%/%;%W%?%s%9!&%l%7%*$HHf3S$9$kMp?t$N@8@.(B
        call fortran_random(ransuu,name,error)
        call error_check(name,error)

        if( ransuu<=acceptance_ratio ) then
 
           if( idou_spin==1 )then
              vector_up=vector_up_work          ! $B@5<0:NMQ$N>l9g99?7$9$k(B

              global_site_table_up=0                   ! $B%5%$%HI=<($K(B
              do i_count=1,TOTAL_UP_ELECTRON
                 global_site_table_up(vector_up(i_count))=1
              end do

! $B5U9TNs:n@.(B
              call make_d_tilde(vector_up,vector_down,1,1,name,error)
              call error_check(name,error)
              call invert(1,1,name,error)
              call error_check(name,error)

              input_up=1
              sign=0
              do i_count=1,TOTAL_UP_ELECTRON
                 input_up=input_up*after_lu_up(i_count,i_count)
                 if( ipiv_up(i_count)/=i_count ) then
                    sign=sign+1
                 end if
              end do

              input_up=input_up*((-1)**sign)
              input_down=stored_input_down
              stored_input_up=input_up

              lambda=denomi_lambda
              alpha_psi0 = input_up * input_down
              alpha_psi1 = input_up * input_down * lambda

              stored_lambda=denomi_lambda
              stored_alpha_psi0=alpha_psi0
              stored_alpha_psi1=alpha_psi1
           else
              vector_down=vector_down_work

              global_site_table_down=0                 ! $B%5%$%HI=<($K(B
              do i_count=1,TOTAL_DOWN_ELECTRON
                 global_site_table_down(vector_down(i_count))=1
              end do

! $B5U9TNs:n@.(B
              call make_d_tilde(vector_up,vector_down,2,1,name,error)
              call error_check(name,error)
              call invert(2,1,name,error)
              call error_check(name,error)

              input_down=1
              sign=0
              do i_count=1,TOTAL_DOWN_ELECTRON
                 input_down=input_down*after_lu_down(i_count,i_count)
                 if( ipiv_down(i_count)/=i_count ) then
                    sign=sign+1
                 end if
              end do

              input_up=stored_input_up
              input_down=input_down*((-1)**sign)
              stored_input_down=input_down

              lambda=denomi_lambda
              alpha_psi0 = input_up * input_down
              alpha_psi1 = input_up * input_down * lambda

              stored_lambda=denomi_lambda
              stored_alpha_psi0=alpha_psi0
              stored_alpha_psi1=alpha_psi1

           end if
        else
!       $B$3$N>l9gEE;R>uBV$KJQ99$O$J$$!#A0$NEE;R>uBV$r$=$N$^$^:N$jF~$l$k(B 
!       $B$h$C$F%U%i%C%0$bJQ2=$J$7(B

           input_up=stored_input_up
           input_down=stored_input_down

           lambda=stored_lambda
           alpha_psi0=stored_alpha_psi0
           alpha_psi1=stored_alpha_psi1

        end if
     end do

     /*
     state_vector(sample_number)%vector_up=vector_up
     state_vector(sample_number)%vector_down=vector_down
     state_vector(sample_number)%alpha_psi0=alpha_psi0
     state_vector(sample_number)%alpha_psi1=alpha_psi1
     state_vector(sample_number)%lambda=lambda
     */
  end if

end subroutine correct_montecarlo_sample



#ifdef DEBUG_CORRECT_MONTECARLO_SAMPLE
program main
  use global_variables
  implicit none

  integer :: sample_number
  integer,dimension(TOTAL_UP_ELECTRON) :: vector_up 
  integer,dimension(TOTAL_DOWN_ELECTRON) :: vector_down
  real(8) :: alpha_psi0,alpha_psi1,lambda
  character(32) :: name
  integer :: error
! local
  integer,dimension(2) :: ic
  real(8) :: zero_approx_energy
  integer :: hf_iteration_result
! 
  sample_number=1 ; ic=1

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

  call random_init(ic,name,error)
  call error_check(name,error)

  call zero_approx(zero_approx_energy,hf_iteration_result,name,error)
  call error_check(name,error)
  
  do while( sample_number<=MAX_MONTECARLO_SAMPLE )
     call correct_montecarlo_sample(sample_number,vector_up,vector_down,&
          alpha_psi0,alpha_psi1,lambda,name,error)
     call error_check(name,error)
     write(*,*) "sample_number=",sample_number
     
     sample_number=sample_number+1
  end do

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CORRECT_MONTECARLO_SAMPLE */

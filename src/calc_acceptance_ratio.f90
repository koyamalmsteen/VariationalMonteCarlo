!
! $Id: calc_acceptance_ratio.f90,v 1.5 2004/02/15 10:52:35 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine calc_acceptance_ratio(gamma_up_work,gamma_down_work,tmp_lambda,&
                                 denomi_lambda,idou_spin,idou_moto,acceptance_ratio,name,error)
  use global_variables
  implicit none

  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: gamma_up_work
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: gamma_down_work
  real(8),intent(in) :: tmp_lambda
  real(8),intent(out) :: denomi_lambda
  integer,intent(in) :: idou_spin,idou_moto
  real(8),intent(out) :: acceptance_ratio
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
  real(8),external :: jastrow,lambda_field
  real(8) :: q
! init
  acceptance_ratio=0 ; name="calc_acceptance_ratio" ; error=0
! init [local]
  denomi_lambda=0
  q=0

  if( idou_spin==1 ) then         ! $B",$,F0$$$?$J$i(B
     do i_count=1,TOTAL_UP_ELECTRON
        q=q+d_tilde_up_inverse(i_count,idou_moto)&
             *d_tilde_up(idou_moto,i_count)
     end do

  else if( idou_spin==2 ) then     ! $B"-$,F0$$$?$J$i(B
     do i_count=1,TOTAL_DOWN_ELECTRON
        q=q+d_tilde_down_inverse(i_count,idou_moto) &
             *d_tilde_down(idou_moto,i_count)
     end do
  end if

! $BJQJ,%Q%i%a!<%?$N?t$,(B>0$B$J$i(B
  if( NUMBER_XI>0 ) then
     denomi_lambda=jastrow(gamma_up_work,gamma_down_work)
  end if

! $BEE>l$N8z2L$,$"$k$J$i(B
  if( electric_field/=0 .and. PERIODIC==1 .and. xi_field/=0 ) then
     denomi_lambda=denomi_lambda*lambda_field(gamma_up_work,gamma_down_work)
  end if

  q=q*denomi_lambda/tmp_lambda

  acceptance_ratio=(abs(q))**2

end subroutine calc_acceptance_ratio



#ifdef DEBUG_CALC_ACCEPTANCE_RATIO
program main
  use global_variables
  implicit none

  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up,gamma_up_work
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down,gamma_down_work
  character(32) :: name
  integer :: error
!
  real(8) :: ransuu
  integer :: number_accept,number_reject
! local
  real(8) :: denominator_up,denominator_down,denominator
  integer :: i_count,j_count,k_count
  real(8) :: acceptance_ratio
! local 
  real(8) :: temp
  real(8) :: zero_approx_energy
  integer :: hf_iteration_result
  real(8) :: bunbo,bunbo_up,bunbo_down
  integer :: idou_spin,idou_moto
  integer,dimension(1) :: ic
  character(32) :: filename
  integer :: kari
! NAG$B4XO"(B
  real(8),dimension(:),allocatable :: wk_space
  integer :: ifail 
! MPI$B4XO"(B
  integer :: IERROR
! init
  denominator_up=0 ; denominator_down=0
  ic(1)=1
!
  zero_approx_energy=0 ; hf_iteration_result=0

  number_accept=0 ; number_reject=0
  d_tilde_up_inverse=0 ; d_tilde_down_inverse=0

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

  call random_init(ic,name,error)
  call error_check(name,error)

  call zero_approx(zero_approx_energy,hf_iteration_result,name,error)
  call error_check(name,error)

! $B=i4|(Bgamma$B$rMp?t$h$j7hDj!#:G=i$@$1(Bsort$B$r$7$F$*$/(B($BJL$K$9$kI,MW$O$J$$$,(B)
  call choice_gamma_sigma(TOTAL_UP_ELECTRON,gamma_up,name,error)
  call error_check(name,error)
  call choice_gamma_sigma(TOTAL_DOWN_ELECTRON,gamma_down,name,error)
  call error_check(name,error)

!
! $B$5$F!"$3$N;~E@$G:G=i$NEE;R>uBV$O7hDj$5$l$?!#(B
! ($B$^$@%5%s%W%k?t$r%$%s%/%j%a%s%H$7$^$;$s!#8e$G$9$k$N$G!#(B)
!

! $B$5$F!"(B $B=i4|$NEE;R>uBV$rD4$Y$k!#(B

  k_count=0

  do i_count=1,TOTAL_UP_ELECTRON
     do j_count=1,TOTAL_DOWN_ELECTRON
        if( gamma_up(i_count)==gamma_down(j_count) ) then ! $B%@%V%k%*%-%e%Q%$(B
           k_count=k_count+1
        end if
     end do
  end do

! $B=i4|&C$N5U9TNs$r:n$k(B
  call make_d_tilde(gamma_up,gamma_down,1,1,name,error)
  call error_check(name,error)
  call make_d_tilde(gamma_up,gamma_down,2,1,name,error)
  call error_check(name,error)
  call invert(1,1,name,error) 
  call error_check(name,error)
  call invert(2,1,name,error) 
  call error_check(name,error)

  acceptance_ratio=1

! $B;n9TEE;R>uBV$r:n6HNN0h$K(B

  gamma_up_work=gamma_up ; gamma_down_work=gamma_down

! $B$3$3$+$i%a%$%s$N%k!<%W(B

  temp=0

!  open(21,file="aho")

  do while( number_accept<MAX_MONTECARLO_SAMPLE )
!
!     if( mod(number_accept,10000)==0 .and. temp/=number_accept ) then
!        write(21,*) number_accept,&
!                  &dble(100*number_accept)/dble(number_accept+number_reject)
!        temp=number_accept
!     end if

     call fortran_random(ransuu,name,error)
     call error_check(name,error)

     if( ransuu<=acceptance_ratio ) then

        number_accept=number_accept+1
        gamma_up=gamma_up_work ; gamma_down=gamma_down_work ! $B@5<0$K:NMQ$7$?(B
!
! $B$3$3$G!"%(%M%k%.!<4|BTCM7W;;$NJ,Jl$N(B<alpha|psi>$B$r$7$F$*$/!#(B
!        
        bunbo_up=1
        do i_count=1,TOTAL_UP_ELECTRON
           bunbo_up=bunbo_up*d_tilde_up(i_count,i_count)
        end do

        bunbo_down=1
        do i_count=1,TOTAL_DOWN_ELECTRON
           bunbo_down=bunbo_down*d_tilde_down(i_count,i_count)
        end do

        bunbo = bunbo_up * bunbo_down

! $BMWAG7W;;(B!!!

        if( idou_spin==1 )then        ! $B",%9%T%s$,F0$$$?$i(B
           call make_d_tilde(gamma_up,gamma_down,1,1,name,error)
           call error_check(name,error)

           call invert(1,1,name,error) 
           call error_check(name,error)
        end if
     else                                  ! $B%j%8%'%/%H$5$l$?$i(B 
        number_reject=number_reject+1
     end if

! $B?7$7$$EE;R>uBV&C(B'$B$r:n$k!#(B
! ($B8E$$EE;R>uBV$O(Bgamma_up,gamma_down$B$H$7$F;D$7$F$*$/(B($B%j%8%'%/%H;~$K;H$&$+$i(B))

     gamma_up_work=gamma_up ; gamma_down_work=gamma_down ! $B%o!<%/MQ$K0\$9(B

     call choice_new_gamma(gamma_up_work,gamma_down_work,&
                           idou_spin,&
                           idou_moto,name,error)
     call error_check(name,error)

! $B0J2<$G!"?7$7$$;n9TEE;R>uBV$NItJ,9TNs$r:n$k!#(B
! $BItJ,9TNs$r:n$k=hM}$O0\F0$7$?%9%T%s$N$_$G9T$J$($P$h$$!#(B
! $B0\F0$7$F$J$1$l$PJQ99$O$J$$$+$i!#(B

     if( idou_spin==1 ) then  ! $B",%9%T%s$,F0$$$?$i(B
        call make_d_tilde(gamma_up_work,gamma_down_work,1,1,name,error)
        call error_check(name,error)
! $B$3$l$O5U9TNs$r:n$kI,MW$,$J$$!#5U9TNs$H3]$19g$o$;$k9TNs$@$+$i(B
     else if( idou_spin==2 ) then ! $B"-%9%T%s$,F0$$$?$i(B 
        call make_d_tilde(gamma_up_work,gamma_down_work,2,1,name,error)
        call error_check(name,error)
! $B$3$l$O5U9TNs$r:n$kI,MW$,$J$$!#5U9TNs$H3]$19g$o$;$k9TNs$@$+$i(B
     end if

! $B%"%/%;%W%?%s%9!&%l%7%*$r5a$a$k!#(B

     call calc_acceptance_ratio(gamma_up,gamma_down,gamma_up_work,&
                                gamma_down_work,idou_spin,&
                                idou_moto,acceptance_ratio,name,error)
     call error_check(name,error)
  end do

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CALC_ACCEPTANCE_RATIO */

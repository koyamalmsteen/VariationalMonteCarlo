!
! $Id: trash.f90,v 1.5 2004/02/15 10:52:35 k-yukino Exp $
!
#include "parameter.h"

subroutine trash(vector_up,vector_down,name,error)
  use global_variables
  implicit none

  integer,dimension(TOTAL_UP_ELECTRON),intent(inout) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(inout) :: vector_down
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer,dimension(TOTAL_UP_ELECTRON) :: vector_up_work
  integer,dimension(TOTAL_DOWN_ELECTRON) :: vector_down_work
  integer :: number_accept,idou_spin,idou_saki,idou_moto
  real(8) :: acceptance_ratio,alpha_psi,alpha_psi_up,alpha_psi_down,&
             tmp_lambda,ransuu
  integer :: sign 
  real(8),external :: jastrow
  integer :: i_count
  real(8) :: denomi_lambda
! init
  name="trash" ; error=0
! 
  number_accept=0 ; idou_spin=0 ; idou_saki=0
  idou_moto=0 ; tmp_lambda=0 ; sign=0 ; ransuu=0
  alpha_psi=0 ; alpha_psi_up=0 ; alpha_psi_down=0
  denomi_lambda=0

  acceptance_ratio=1         ! $B:G=i$N>uBV$OI,$:%"%/%;%W%H$9$k$N$G(B

! TRASH$B%b%s%F%+%k%m%9%F%C%W$@$1EE;R>uBV$r<N$F$k!#(B

  vector_up_work=vector_up ; vector_down_work=vector_down
  number_accept=1

  do while( number_accept<=TRASH*(TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON) )
! ($B$3$3$G%b%s%F%+%k%m%5%s%W%k$N3NN(7W;;$r$9$k(B)
     call fortran_random(ransuu,name,error)
     call error_check(name,error)

     if( ransuu<=acceptance_ratio ) then
        number_accept=number_accept+1
        vector_up=vector_up_work          ! $B@5<0:NMQ$N>l9g99?7$9$k(B
        vector_down=vector_down_work      !
     else
        number_accept=number_accept+1
!!
!!      $B$3$N>l9gEE;R>uBV$KJQ99$O$J$$!#A0$NEE;R>uBV$r$=$N$^$^:N$jF~$l$k(B 
!!
     end if

! $B9bB.2=$N$?$a$K%5%$%HI=<($K$7$F$*$/(B
     if( idou_spin==1 ) then
        global_site_table_up=0
        do i_count=1,TOTAL_UP_ELECTRON
           global_site_table_up(vector_up(i_count))=1
        end do
     else
        global_site_table_down=0
        do i_count=1,TOTAL_DOWN_ELECTRON
           global_site_table_down(vector_down(i_count))=1
        end do
     end if

! $B5U9TNs$r:n$kJ}K!(B
     if( idou_spin==1 )then   ! $B",%9%T%s$,F0$$$?$i(B
        call make_d_tilde(vector_up,vector_down,1,1,name,error)
        call error_check(name,error)
        call invert(1,1,name,error) 
        call error_check(name,error)
     else                                   ! $B"-%9%T%s$,F0$$$?$i(B
        call make_d_tilde(vector_up,vector_down,2,1,name,error)
        call error_check(name,error)
        call invert(2,1,name,error)
        call error_check(name,error)
     end if
!
! $B$3$3$G!"%(%M%k%.!<4|BTCM7W;;$NJ,Jl$N(B<alpha|psi>=lambda_alpha<alpha|HF>
! $B$r$7$F$*$/!#(B
!
     alpha_psi_up=1 
     sign=0        
     do i_count=1,TOTAL_UP_ELECTRON
        alpha_psi_up=alpha_psi_up*after_lu_up(i_count,i_count)
        if( ipiv_up(i_count)/=i_count ) then
           sign=sign+1
        end if
     end do

     alpha_psi_up=alpha_psi_up*((-1)**sign)

     alpha_psi_down=1
     sign=0
     do i_count=1,TOTAL_DOWN_ELECTRON
        alpha_psi_down=alpha_psi_down*after_lu_down(i_count,i_count)
        if( ipiv_down(i_count)/=i_count ) then
           sign=sign+1
        end if
     end do

     alpha_psi_down=alpha_psi_down*((-1)**sign)

     tmp_lambda=jastrow(vector_up,vector_down)
     alpha_psi = alpha_psi_up * alpha_psi_down*tmp_lambda

! $B?7$7$$EE;R>uBV&C(B'$B$r:n$k!#(B
! ($B8E$$EE;R>uBV$O(Bvector_up,vector_down$B$H$7$F;D$7$F$*$/(B($B%j%8%'%/%H;~$K;H$&$+$i(B))

     vector_up_work=vector_up ; vector_down_work=vector_down   ! $B%o!<%/MQ$K0\$9(B

     if( PROJECTION==0 ) then 
        call choice_new_gamma(vector_up_work,vector_down_work,&
             idou_spin,idou_moto,idou_saki,name,error)
        call error_check(name,error)
     else
        call choice_new_gamma_projection(vector_up_work,vector_down_work,&
             idou_spin,idou_moto,idou_saki,name,error)
        call error_check(name,error)
     end if

! $B0J2<$G!"?7$7$$;n9TEE;R>uBV$NItJ,9TNs$r:n$k!#(B
! $BItJ,9TNs$r:n$k=hM}$O0\F0$7$?%9%T%s$N$_$G9T$J$($P$h$$!#(B
! $B0\F0$7$F$J$1$l$PJQ99$O$J$$$+$i!#(B

! $B5U9TNs$r:n$kJ}K!(B  
     if( idou_spin==1 ) then  ! $B",%9%T%s$,F0$$$?$i(B
        call make_d_tilde(vector_up_work,vector_down_work,1,1,name,error)
        call error_check(name,error)
! down$B$O0\F0$7$F$$$J$$$N$G!"2?$b$7$J$$(B
     else if( idou_spin==2 ) then ! $B"-%9%T%s$,F0$$$?$i(B 
        call make_d_tilde(vector_up_work,vector_down_work,2,1,name,error)
        call error_check(name,error)
! up$B$O0\F0$7$F$$$J$$$N$G!"2?$b$7$J$$(B
     end if

! $B%"%/%;%W%?%s%9!&%l%7%*$N7W;;(B
     call calc_acceptance_ratio(vector_up_work,vector_down_work,tmp_lambda,&
                                denomi_lambda,idou_spin,idou_moto,&
                                acceptance_ratio,name,error)
     call error_check(name,error)
  end do

end subroutine trash



#ifdef DEBUG_TRASH
program main
  implici nonte

end program main
#endif /* DEBUG_TRASH */

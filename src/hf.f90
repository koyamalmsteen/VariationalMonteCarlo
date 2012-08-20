!
! $Id: hf.f90,v 1.26 2004/02/15 10:52:35 k-yukino Exp $
!
#include "parameter.h"

subroutine hf(hf_energy,hf_iteration_result,name,error)
  use global_variables
  implicit none

  real(8),intent(out) :: hf_energy
  integer,intent(out) :: hf_iteration_result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8),dimension(:,:),allocatable :: rho_up,rho_down,rho_up_new,rho_down_new
  real(8) :: ransuu,ransuu_tmp,result_up,result_down
  integer :: ier
  integer :: i_count,j_count
  real(8) :: hf_energy_coulomb,hf_energy_transfer
! lapack
  integer :: lapack_nnn,lapack_lda,lapack_lwork,lapack_ifail
  real(8),dimension(:),allocatable :: lapack_www_up,lapack_www_down
  real(8),dimension(:),allocatable :: lapack_work
  character*1 :: lapack_job,lapack_uplo
! init
  hf_energy=0 ; hf_iteration_result=0
  name="hf" ; error=0
! int local
  i_count=0 ; j_count=0 ; result_up=0 ; result_down=0
  hf_energy_transfer=0 ; hf_energy_coulomb=0

! init
  lapack_job='v'                   ! eigen value and eigen vector are needed
  lapack_uplo='u'                  ! upper triangular part of A is stored
  lapack_nnn=TOTAL_SITE_NUMBER
  lapack_lda=TOTAL_SITE_NUMBER
  lapack_lwork=64*TOTAL_SITE_NUMBER
  lapack_ifail=0

  allocate(lapack_www_up(lapack_lda),stat=ier)
  call stat_check("lapack_www_up","hf",1,ier)
  allocate(lapack_www_down(lapack_lda),stat=ier)
  call stat_check("lapack_www_down","hf",1,ier)
  allocate(lapack_work(lapack_lwork),stat=ier)
  call stat_check("lapack_work","hf",1,ier)
! init
  lapack_www_up=0 ; lapack_www_down=0 ; lapack_work=0
!
  allocate(rho_up(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("rho_up","hf",1,ier)
  allocate(rho_down(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("rho_down","hf",1,ier)
  allocate(rho_up_new(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("rho_up_new","hf",1,ier)
  allocate(rho_down_new(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("rho_down_new","hf",1,ier)
! init
  unitary_d_up=0 ; unitary_d_down=0
  rho_up=0 ; rho_down=0
  rho_up_new=0 ; rho_down_new=0

!!! decision of primitive rho_up and rho_down using random number

!
! $B!y#1(B $B=i4|$N(Brho$B$r7h$a$k%k!<%A%s(B($BMp?t%F!<%V%k$N=i4|2=$O%a%$%s$G(B)
!
  if( HF_INITIAL_RHO==0 ) then                 ! 1,0,1,0...  
     call fortran_random(ransuu_tmp,name,error)! $B","-$+"-",$r$r7h$a$k$@$1(B
     call error_check(name,error)              ! $B0lHVDc%(%M%k%.!<$,F@$i$l$k(B
                                               ! $B$*$9$9$a!*(B
     do i_count=1,TOTAL_SITE_NUMBER
        if( mod(i_count,2)==1 ) then
           if( ransuu_tmp<0.5 ) then
              rho_up(i_count,i_count)=1 ; rho_down(i_count,i_count)=0 
           else
              rho_up(i_count,i_count)=0 ; rho_down(i_count,i_count)=1
           end if
        else
           if( ransuu_tmp<0.5 ) then
              rho_up(i_count,i_count)=0 ; rho_down(i_count,i_count)=1
           else
              rho_up(i_count,i_count)=1 ; rho_down(i_count,i_count)=0
           end if
        endif
     end do
  else if( HF_INITIAL_RHO==1 ) then     ! sdw$BE*$J2r$r2>Dj$9$k(B
                                        ! 0.5>,0.5<,0.5>,0.5<($B$b$7$/$O5U(B)$B$H(B
                                        ! $B$J$C$F$$$k(B
     call fortran_random(ransuu_tmp,name,error) ! $BJB$Y$k=gHV$r7h$a$k0Y(B
     call error_check(name,error)

     do i_count=1,TOTAL_SITE_NUMBER
        call fortran_random(ransuu,name,error)
        call error_check(name,error)
        if( ransuu_tmp>0.5 ) then
           rho_up(i_count,i_count)=0.5+ransuu/2*((-1)**i_count)
        else
           rho_up(i_count,i_count)=0.5-ransuu/2*((-1)**i_count)
        end if

        call fortran_random(ransuu,name,error)
        call error_check(name,error)
        if( ransuu_tmp<0.5 ) then
           rho_down(i_count,i_count)=0.5-ransuu/2*((-1)**(i_count-1))
        else
           rho_down(i_count,i_count)=0.5+ransuu/2*((-1)**(i_count-1))
        end if
     end do
  else if( HF_INITIAL_RHO==2 ) then            ! sdw$BE*$J2r$r2>Dj$7$J$$(B
     call fortran_random(ransuu,name,error)   ! $BEE;R$NB8:_3NN($O3F%5%$%H$G(B
     call error_check(name,error)         ! 1$B$K$J$k$h$&$K$J$C$F$$$k!#(B

     if( ransuu<0.5 ) then                ! $B","-$N=g$KJB$Y$k(B
        do i_count=1,TOTAL_SITE_NUMBER
           call fortran_random(ransuu,name,error)
           call error_check(name,error)

           rho_up(i_count,i_count)=ransuu
           rho_down(i_count,i_count)=1-ransuu
        end do
     else 
        do i_count=1,TOTAL_SITE_NUMBER    ! $B","-$N=g(B
           call fortran_random(ransuu,name,error)
           call error_check(name,error)

           rho_up(i_count,i_count)=1-ransuu
           rho_down(i_count,i_count)=ransuu
        end do
     end if
  else if( HF_INITIAL_RHO==3 ) then       ! $BA4$/%i%s%@%`$K:n$k(B($B<}B+$7$K$/$$(B)
     do i_count=1,TOTAL_SITE_NUMBER
        call fortran_random(ransuu,name,error)
        call error_check(name,error) 
        rho_up(i_count,i_count)=ransuu
        
        call fortran_random(ransuu,name,error)
        call error_check(name,error)
        rho_down(i_count,i_count)=ransuu
     end do
  end if

! $B$d$C$D$1(B
  if( JIGEN==2 .and. TOTAL_SITE_NUMBER==16 ) then
     rho_up(1,1)=1
     rho_up(2,2)=0
     rho_up(3,3)=1
     rho_up(4,4)=0

     rho_up(5,5)=0
     rho_up(6,6)=1
     rho_up(7,7)=0
     rho_up(8,8)=1

     rho_up(9,9)=1
     rho_up(10,10)=0
     rho_up(11,11)=1
     rho_up(12,12)=0

     rho_up(13,13)=0
     rho_up(14,14)=1
     rho_up(15,15)=0
     rho_up(16,16)=1
!
     rho_down(1,1)=0
     rho_down(2,2)=1
     rho_down(3,3)=0
     rho_down(4,4)=1

     rho_down(5,5)=1
     rho_down(6,6)=0
     rho_down(7,7)=1
     rho_down(8,8)=0

     rho_down(9,9)=0
     rho_down(10,10)=1
     rho_down(11,11)=0
     rho_down(12,12)=1

     rho_down(13,13)=1
     rho_down(14,14)=0
     rho_down(15,15)=1
     rho_down(16,16)=0
  end if

! $B:GBg7+$jJV$72s?t$^$G7W;;$9$k!#$=$l$^$G$K<}B+$7$?$i%k!<%W$+$iH4$1$k!#(B

  do i_count=1,MAX_HF_ITERATION
     call make_fock_matrix(rho_up,rho_down,name,error)
     call error_check(name,error)

! $BBP3Q2=(B($B<BBP>N9TNs(B)
     lapack_work=0
     call dsyev(lapack_job,lapack_uplo,lapack_nnn,unitary_d_up,lapack_lda,&
          lapack_www_up,lapack_work,lapack_lwork,lapack_ifail)
     lapack_work=0
     call dsyev(lapack_job,lapack_uplo,lapack_nnn,unitary_d_down,lapack_lda,&
          lapack_www_down,lapack_work,lapack_lwork,lapack_ifail)

! $B?7$7$$(Brho$B$r5a$a$k(B
     call calc_new_rho(rho_up_new,rho_down_new,name,error)
     call error_check(name,error)

! $B=i4|$N(Brho$B$H?7$7$$(Brho$B$H$NHf3S(B($B<}B+H=Dj(B)
     
     result_up=0 ; result_down=0             ! $B=i4|2=(B($BI,?\(B)

     do j_count=1,TOTAL_SITE_NUMBER
        result_up=result_up &
             +( rho_up_new(j_count,j_count)-rho_up(j_count,j_count) )**2
        result_down=result_down &
             +( rho_down_new(j_count,j_count)-rho_down(j_count,j_count) )**2
     end do

     result_up=(result_up/TOTAL_SITE_NUMBER)**0.5          ! $B<}B+H=Dj<0(B
     result_down=(result_down/TOTAL_SITE_NUMBER)**0.5      !

! $B0J2<$N>r7o$G<}B+$H$7$F!"(Bhf$B2r$,$b$H$^$C$?$3$H$H$9$k!#(B
     if( result_up<dble(HF_CONVERGENCE_CONDITION) .and. &
          &result_down<dble(HF_CONVERGENCE_CONDITION) ) then 
        exit
     end if
     
! $B<}B+$7$F$J$+$C$?$i!"?7$7$/:n$C$?(Brho_up_new,rho_down_new$B$r85$K!"$5$i$K(B
! $B<}B+$9$k$^$G7W;;$r7+$jJV$9!#(B
     rho_up=rho_up_new ; rho_down=rho_down_new
     rho_up_new=0 ; rho_down_new=0
  end do

  if( i_count>MAX_HF_ITERATION ) then
     hf_iteration_result=i_count-1  
     write(*,*) "i_count,MAX_HF_ITERATION=",i_count,MAX_HF_ITERATION
     write(*,*) "rho_up=",rho_up
     write(*,*) "rho_down=",rho_down
     name="hf" ; error=-2 ; return            ! $B<}B+$7$J$+$C$?>l9g%(%i!<$r=P$9(B
  else
     hf_iteration_result=i_count 
  end if

  deallocate(rho_up_new,stat=ier)
  call stat_check("rho_up_new","hf",2,ier)
  deallocate(rho_down_new,stat=ier)
  call stat_check("rho_down_new","hf",2,ier)

  deallocate(lapack_www_up,stat=ier)
  call stat_check("lapack_www_up","hf",2,ier)
  deallocate(lapack_www_down,stat=ier)
  call stat_check("lapack_www_down","hf",2,ier)
  deallocate(lapack_work,stat=ier)
  call stat_check("lapack_work","hf",2,ier)

!!! E_HF=<HF|H|HF>$B$r5a$a$k(B

!
! $B%H%i%s%9%U%!!<(B
!
  hf_energy_transfer=0
  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        if( neighbor_table(i_count,j_count)==1 ) then
           hf_energy_transfer=hf_energy_transfer+dble(TRANSFER)&
                           *(rho_up(i_count,j_count)+rho_down(i_count,j_count))
        end if
     end do
  end do
!
! $B%/!<%m%s(B
!
  hf_energy_coulomb=0
  do i_count=1,TOTAL_SITE_NUMBER
     hf_energy_coulomb=hf_energy_coulomb+dble(COULOMB)&
                          *rho_up(i_count,i_count)*rho_down(i_count,i_count)
  end do


  hf_energy=hf_energy_transfer+hf_energy_coulomb

end subroutine hf



#ifdef DEBUG_HF
program main
  use global_variables
  implicit none

  character(32) :: name
  integer :: error
! MPI$B4XO"(B
  integer :: IERROR
! local
  integer,dimension(1) :: ic
  real(8) :: hf_energy
  integer :: hf_iteration_result
! init
  name="main" ; error=0
  unitary_d_up=0 ; unitary_d_down=0 
  hf_energy=0 ; hf_iteration_result=0

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

  call random_init(ic,name,error)
  call error_check(name,error)

  call hf(hf_energy,hf_iteration_result,name,error)
  call error_check(name,error)

  write(*,*) "JIGEN=",JIGEN,"TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "TOTAL_UP_ELECTRON=",TOTAL_UP_ELECTRON
  write(*,*) "TOTAL_DOWN_ELECTRON=",TOTAL_DOWN_ELECTRON
  write(*,*) "COULOMB=",COULOMB,"TRANSFER=",TRANSFER
  write(*,*) "hf_energy=",hf_energy

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_HF */

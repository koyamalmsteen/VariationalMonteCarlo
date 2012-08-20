!
! $Id: choice_gamma_sigma_projection.f90,v 1.1 2004/02/15 10:52:35 k-yukino Exp $
!
#include "parameter.h"

subroutine choice_gamma_sigma_projection(total_sigma_electron,gamma_sigma,&
     name,error)
  use global_variables
  implicit none

  include "mpif.h"
  
  integer,intent(in) :: total_sigma_electron
  integer,dimension(total_sigma_electron),intent(inout) :: gamma_sigma
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer,dimension(total_sigma_electron) :: dummy     ! $B%@%_!<L$;HMQ(B
  integer :: i_count,j_count,k_count,l_count,ransuu_int,ier,&
             total_empty_site_number
  integer,dimension(TOTAL_SITE_NUMBER) :: site_table
  integer,dimension(:),allocatable :: empty_site_table
! init arguments
  name="choice_gamma_sigma" ; error=0
! local init
  i_count=0 ; j_count=0 ; k_count=0 ; ransuu_int=0 ; ier=0
  site_table=0 ; total_empty_site_number=TOTAL_SITE_NUMBER

! argument check
  if( total_sigma_electron<=0 ) then       ! $BEE;R?t$,(B0$B$N;~$bBP1~$7$F$$$J$$(B
     error=-92 ; return
  else if( TOTAL_SITE_NUMBER<total_sigma_electron ) then
     error=-912 ; return
  end if

! $B$3$l$+$i=i4|EE;R>uBV$rA*$V$N$G!"0l$D$bEE;R$OKd$a$i$l$F$$$J$$(B
! $B$@$+$i%5%$%H%F!<%V%k$b(Bempty$B$G$h$$(B

! $B6u%5%$%H$rJ];}$7$F$*$/NN0h3NJ](B

  do i_count=1,total_sigma_electron
! $B6u$N%5%$%H$rD4$Y$F!"$=$NCf$+$i0l$D$N%5%$%H$rEE;R$GKd$a$k(B
     allocate(empty_site_table(total_empty_site_number),stat=ier)
     call stat_check("empty_site_table","choice_gamma_sigma",1,ier)
     empty_site_table=0  

     total_empty_site_number=0
     do l_count=1,TOTAL_SITE_NUMBER
        if( site_table(l_count)==0 ) then
           total_empty_site_number=total_empty_site_number+1
        end if
     end do

!  $B0J2<$N$h$&$K$7$F6u%5%$%H$N%A%'%C%/$r9T$J$&!#(B

     k_count=1

     do j_count=1,TOTAL_SITE_NUMBER
        if( site_table(j_count)==0 ) then
           empty_site_table(k_count)=j_count
           k_count=k_count+1
        end if
     end do

! $B%A%'%C%/$7$?6u%5%$%H$+$iEE;R$rKd$a$k%5%$%H$r$R$H$D7h$a$k(B

     call fortran_random2(1,total_empty_site_number,ransuu_int,name,error)
     call error_check(name,error)

     gamma_sigma(i_count)=empty_site_table(ransuu_int)

!
! gamma_sigma$B$,$R$H$DEE;R$GKd$a$i$l$?$N$G!"?7$7$$(Bsite_table$B$r:n$kI,MW$,$G$-$?!#(B
! $B$D$^$j!"$"$?$i$7$$(Bsite_table$B$r:n$j!">e5-$N6u%5%$%H$rD4$Y$F!"?7$7$/EE;R$r(B
! $B5M$a$k:n6H$r$9$k$o$1$G$"$k!#(B
!  $B$3$3$G$O!"0J2<$K%5%$%H$X$NJQ49%W%m%0%i%`$r=q$/!#(B
!
     deallocate(empty_site_table,stat=ier)
     call stat_check("empty_site_table","choice_gamma_sigma",2,ier)

! $B%5%$%HI=<($X$NJQ49(B
     site_table=0

     do j_count=1,i_count
        site_table(gamma_sigma(j_count))=1
     end do

     total_empty_site_number=total_empty_site_number-1
  end do

! slatec library($B>:$Y$-(B(1)$B$KJB$Y$k!#HFMQ@-$N$?$a!"(Bcxml$B$+$i(Bslatec$B$K=q$-49$($?(B)
  call isort(gamma_sigma,dummy,total_sigma_electron,1)

end subroutine choice_gamma_sigma_projection



#ifdef DEBUG_CHOICE_GAMMA_SIGMA
program main
  implicit none

  include "mpif.h"

  integer :: total_sigma_electron
  integer,dimension(:),allocatable :: gamma_sigma
  character(32) :: name
  integer :: error
! local
  integer :: i_count,ier
  integer,external :: iargc
  character(32) :: string
! MPI$B4XO"(B
  integer :: IERROR
! init
  name="main" ; error=0
  i_count=0 ; ier=0

  call MPI_INIT(IERROR)

  if( iargc()==2 )then
     call getarg(2,string) ; read(string,*) total_sigma_electron
  else
     write(*,*) "total_sigma_electron >" ; read(*,*) total_sigma_electron
  end if

  allocate(gamma_sigma(total_sigma_electron),stat=ier)
  call stat_check("gamma_sigma","choice_gamma_sigma",1,ier)

! init
  gamma_sigma=0

! $BMp?t%F!<%V%k$N=i4|2=(B ($BMWCm0U2U=j(B)
  call random_init(name,error)
  call error_check(name,error)

  call choice_gamma_sigma(total_sigma_electron,gamma_sigma,name,error)
  call error_check(name,error)

  write(*,*) "gamma_sigma=",gamma_sigma

!
  deallocate(gamma_sigma,stat=ier)
  call stat_check("gamma_sigma","choice_gamma_sigma",2,ier)

  
  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CHOICE_GAMMA_SIGMA */ 

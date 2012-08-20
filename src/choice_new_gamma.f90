!
! $Id: choice_new_gamma.f90,v 1.18 2003/02/06 09:25:11 k-yukino Exp $
!
#include "parameter.h"

subroutine choice_new_gamma(gamma_up,gamma_down,idou_moto_electron_spin,&
                            idou_moto_electron,idou_saki_site,name,error)
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: gamma_down
  integer,intent(out) :: idou_moto_electron_spin
  integer,intent(out) :: idou_moto_electron
  integer,intent(out) :: idou_saki_site
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init 
  idou_moto_electron_spin=-1              ! $B$I$C$A$N%9%T%s$NEE;R$rF0$+$7$?$+(B
                                          !
  idou_moto_electron=0                    ! $B$3$l$r%A%'%C%/$7$F$*$/$H!"(B
                                          ! $B%"%/%;%W%?%s%9!&%l%7%*7W;;$r(B
  idou_saki_site=0                        ! $B9bB.2=$G$-$k$N$G!#(B
                                          ! $B$3$l$b=EMW$G$9!#$3$N9TNsMWAG$N$_(B
                                          ! $BJQ99$,$"$k$o$1$G$9!#(B
  name="choice_new_gamma" ; error=0

  if( GAMMA_ALGORITHM==0 ) then                ! $B$^$C$?$/%i%s%@%`(B
     do while( idou_moto_electron_spin==-1 )
        call choice_new_gamma_part1(gamma_up,gamma_down,&
                                    idou_moto_electron_spin,&
                                    idou_moto_electron,&
                                    idou_saki_site,&
                                    name,error)
        call error_check(name,error)
     end do
  else if( GAMMA_ALGORITHM==1 ) then           ! $B:F6a@\%5%$%H$J6u%5%$%H(B
     do while( idou_moto_electron_spin==-1 ) 
        call choice_new_gamma_part2(gamma_up,gamma_down,&
                                    idou_moto_electron_spin,&
                                    idou_moto_electron,&
                                    idou_saki_site,name,error)
        call error_check(name,error)
     end do
  end if

end subroutine choice_new_gamma



#ifdef DEBUG_CHOICE_NEW_GAMMA
program main
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  character(32) :: name
  integer :: error
! local
  integer :: ier,i_count
  integer :: idou_moto_electron_spin,idou_moto_electron,idou_saki_site
! MPI$B4XO"(B
  integer :: IERROR
! init
  gamma_up=0 ; gamma_down=0 ; name="main" ; error=0
! init [local]
  ier=0 ; i_count=0
  idou_moto_electron_spin=0 ; idou_moto_electron=0 ; idou_saki_site=0

  call MPI_INIT(IERROR)
  
  call init_global_variables(name,error)
  call error_check(name,error)

  call random_init(name,error)
  call error_check(name,error)

  do i_count=1,TOTAL_UP_ELECTRON
     write(*,*) "gamma_up",i_count,">" ; read(*,*) gamma_up(i_count)
  end do

  do i_count=1,TOTAL_DOWN_ELECTRON
     write(*,*) "gamma_down",i_count,">" ; read(*,*) gamma_down(i_count)
  end do

  do i_count=1,10 
     call choice_new_gamma(gamma_up,gamma_down,idou_moto_electron_spin,&
                           idou_moto_electron,idou_saki_site,name,error)
     call error_check(name,error)

     write(*,*) i_count,",",gamma_up,"-",gamma_down
  end do

  write(*,*) "JIGEN=",JIGEN
  write(*,*) "TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER


  call MPI_FINALIZE(IERROR)

end program main

#endif /* DEBUG_CHOICE_NEW_GAMMA */

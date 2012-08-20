! 
! $Id: choice_new_gamma_part2.f90,v 1.11 2002/12/20 06:15:13 k-yukino Exp k-yukino $
!

! $B:F6a@\%5%$%H$K!"EE;R$rF0$+$9!#:F6a@\$,A4$FKd$^$C$F$$$?$i!"(B
! idou_moto_electron_spin=-1$B$rJV$9!#8F$S=P$785$G!"$=$l$r%A%'%C%/$7$F(B
! $B:FEY$3$N%k!<%A%s$r8F$S=P$9!#(B

#include "parameter.h"

subroutine choice_new_gamma_part2(gamma_up,gamma_down,idou_moto_electron_spin,&
                                  idou_moto_electron,idou_saki_site,name,error)
  use global_variables
  implicit none

  include "mpif.h"

! $B$A$J$_$K(Bidou_moto_electron$B$O%5%$%HI=<($G$O$J$$$h!#%,%s%^I=<($M"v(B
! idou_saki_site$B$O%5%$%HI=<((B
  integer,intent(out) :: idou_moto_electron_spin,idou_moto_electron ,&
                         idou_saki_site 
  integer,dimension(TOTAL_UP_ELECTRON),intent(inout) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(inout) :: gamma_down
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: ier,ransuu_int
  integer,dimension(JIGEN*2) :: nearest_neighbor
  integer,dimension(TOTAL_SITE_NUMBER) :: site_table
  integer :: i_count,j_count,temp
  integer :: empty_site_number
  integer,external :: state_include
  integer,dimension(:),allocatable :: empty_site_table
! init
  idou_moto_electron_spin=0 ; idou_moto_electron=0
  idou_saki_site=0
  name="choice_new_gamma_part2" ; error=0
! local init
  i_count=0 ; site_table=0 ; temp=0

  do while (1)

! $B",$H"-$NEE;R$r$4$C$A$c$K$7$F!"$=$3$+$i0\F085EE;R$r0l$DA*$S=P$9(B
     call fortran_random2(1,TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON,&
                          idou_moto_electron,name,error)
     call error_check(name,error)

! 1$B!A(BTOTAL_UP_ELECTRON$B$^$G$J$i",$NEE;R$,A*$P$l$k!#(B
! TOTAL_UP_ELECTORN+1$B!A(BTOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON$B$^$G$J$i"-EE;R!#(B

     if( idou_moto_electron>=1 .and. &
         idou_moto_electron<=TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON ) then

        if( idou_moto_electron<=TOTAL_UP_ELECTRON ) then

           idou_moto_electron_spin=1               ! up$B%9%T%s(B 
           idou_moto_electron=idou_moto_electron   ! 1$B!A(BTOTAL_UP_ELECTRON$B$^$G(B
                                                   ! $B$J$N$G$=$N$^$^(B 
        else if( idou_moto_electron>=TOTAL_UP_ELECTRON+1 ) then
           idou_moto_electron_spin=2               ! down$B%9%T%s(B
           idou_moto_electron=idou_moto_electron  &! $B$3$3$G!":85-$N$h$&$K(B
                              -TOTAL_UP_ELECTRON   !  $BCV$-49$($k(B
        end if
     end if

!
! $B$3$NCJ3,$G0\F085$O7h$^$C$?!#0\F0@h$O$3$N0\F085$KNY@\$9$k(B(JIGEN*2)$B%5%$%H$+$i(B
! $B%i%s%@%`$KA*$V!#:F6a@\%5%$%H%F!<%V%k$r:n$k!#(B
!

     j_count=1                                   ! $B%+%&%s%?(B  

     if( idou_moto_electron_spin==1 ) then       ! up$B%9%T%s$,A*Br$5$l$?=hM}(B
        do i_count=1,TOTAL_SITE_NUMBER
           if( neighbor_table(gamma_up(idou_moto_electron),i_count)==1 ) then
              nearest_neighbor(j_count)=i_count
              j_count=j_count+1 
           end if
        end do
     else                                        ! down$B%9%T%s$,A*Br$5$l$?=hM}(B
        do i_count=1,TOTAL_SITE_NUMBER
           if( neighbor_table(gamma_down(idou_moto_electron),i_count)==1 ) then
              nearest_neighbor(j_count)=i_count
              j_count=j_count+1
           end if
        end do
     end if

!
! $B:F6a@\%5%$%H$rLVMe$7$?%F!<%V%k$,$G$-$?!#$=$NCf$+$i!"6u%5%$%H$N$_$N%F!<%V%k(B
! $B$r:n$j=P$9!#(B
!
     empty_site_number=0

     if( idou_moto_electron_spin==1 ) then         ! up$B%9%T%s$,A*Br$5$l$?=hM}(B
! $B$^$:!":F6a@\6u%5%$%H?t$rD4$Y$k(B
        do i_count=1,JIGEN*2
! moshi hashi de nakereba
           if( nearest_neighbor(i_count)/=0 ) then
              if( state_include(TOTAL_UP_ELECTRON,&
                   &nearest_neighbor(i_count),gamma_up)/=1 ) then
                 empty_site_number=empty_site_number+1
              end if
           end if
        end do
     else                                          ! down$B%9%T%s$,A*Br$5$l$?=hM}(B
! $B$^$:!":F6a@\6u%5%$%H?t$rD4$Y$k(B
        do i_count=1,JIGEN*2
! moshi hashi de nakereba
           if( nearest_neighbor(i_count)/=0 ) then
              if( state_include(TOTAL_DOWN_ELECTRON,&
                   &nearest_neighbor(i_count),gamma_down)/=1 ) then
                 empty_site_number=empty_site_number+1
              end if
           end if
        end do
     end if

     if( empty_site_number==0 ) then            ! $B$b$7:F6a@\$,A4$FKd$^$C$F(B
        idou_moto_electron_spin=-1 ; return     ! $B$$$?$i(B
     end if

! $B$b$7!"0\F085$N:F6a@\$,A4$F!"Kd$a$i$l$F$$$?$i!"0J2<$N=hM}$r9T$J$o$J$$$G!"(B
! $B0\F085$rA*$SD>$7!";O$a$+$i$d$j$J$*$9!#(B

     if( empty_site_number/=0 ) then               ! $B0\F0@h$,$"$l$P(B
  
        if( idou_moto_electron_spin==1 ) then      ! up$B%9%T%s$,A*Br$5$l$?=hM}(B
           allocate(empty_site_table(empty_site_number),stat=ier)
           call stat_check("empty_site_table","choice_new_gamma_part2",1,ier)
        
! $B<B:]$K:F6a@\6u%5%$%H$N%F!<%V%k$r:n$k(B
           
           j_count=1
        
           do i_count=1,JIGEN*2
! moshi hashi denai nara
              if( nearest_neighbor(i_count)/=0 ) then
                 if( state_include(TOTAL_UP_ELECTRON,&
                      &nearest_neighbor(i_count),gamma_up)/=1 ) then
                    empty_site_table(j_count)=nearest_neighbor(i_count)
                    j_count=j_count+1
                 end if
              end if
           end do
        else                                    ! down$B%9%T%s$,A*Br$5$l$?=hM}(B
           allocate(empty_site_table(empty_site_number),stat=ier)
           call stat_check("empty_site_table","choice_new_gamma_part2",1,ier)
        
! $B<B:]$K:F6a@\6u%5%$%H$N%F!<%V%k$r:n$k(B

           j_count=1
        
           do i_count=1,JIGEN*2
! moshi hashi de nainara
              if( nearest_neighbor(i_count)/=0 ) then
                 if( state_include(TOTAL_DOWN_ELECTRON,&
                      &nearest_neighbor(i_count),gamma_down)/=1 ) then
                    empty_site_table(j_count)=nearest_neighbor(i_count)
                    j_count=j_count+1
                 end if
              end if
           end do
        end if

!$B0J2<$NMp?t$GHt$P$9@h$r7h$a$k!#(B
        call fortran_random2(1,empty_site_number,ransuu_int,name,error)
        call error_check(name,error)

! $B$G$O<B:]$KCV$-49$($k(B($B$3$3$O%5%$%HI=<($NJ}$,JXMx$J$N$G!"I=<($rJQ99$7$?$$(B
! $B5$;}$A$O;3!9$@$,!"<B$OITET9g$,$"$k$N$G$7$J$$(B)$B"+%"%/%;%W%?%s%9%l%7%*$r(B
! $B5a$a$k9bB.2=%k!<%A%s$N;~:$$j$^$9!#(B
        if( idou_moto_electron_spin==1 ) then      ! up$B%9%T%s$,A*Br$5$l$?=hM}(B
           gamma_up(idou_moto_electron)=empty_site_table(ransuu_int)
           idou_saki_site=empty_site_table(ransuu_int)
           return
        else                                       ! down$B%9%T%s$,A*Br$5$l$?=hM}(B
           gamma_down(idou_moto_electron)=empty_site_table(ransuu_int)
           idou_saki_site=empty_site_table(ransuu_int)
           return
        end if
     end if
  end do

end subroutine choice_new_gamma_part2



#ifdef DEBUG_CHOICE_NEW_GAMMA_PART2
program main
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  character(32) :: name
  integer :: error
! local
  integer :: idou_moto_electron_spin,idou_moto_electron,idou_saki_site
  integer :: ier,i_count
! MPI$B4XO"(B
  integer :: IERROR
! init
  gamma_up=0 ; gamma_down=0 ; name="main" ; error=0
! init [local]
  idou_moto_electron_spin=0 ; idou_moto_electron=0 ; idou_saki_site=0
  ier=0 ; i_count=0
! init [MPI]
  IERROR=0

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

  call choice_new_gamma_part2(gamma_up,gamma_down,&
                              idou_moto_electron_spin,idou_moto_electron,&
                              idou_saki_site,name,error)
  call error_check(name,error)

  write(*,*) "JIGEN=",JIGEN
  write(*,*) "TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "gamma_up=",gamma_up
  write(*,*) "gamma_down=",gamma_down
  write(*,*) "idou_moto_electron_spin=",idou_moto_electron_spin
  write(*,*) "idou_moto_electron=",idou_moto_electron
  write(*,*) "idou_saki_site=",idou_saki_site

  call MPI_FINALIZE(IERROR)

end program main

#endif /* DEBUG_CHOICE_NEW_GAMMA_PART2 */

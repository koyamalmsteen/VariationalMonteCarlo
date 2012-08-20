!
! $Id: count_number_xi.f90,v 1.10 2003/11/20 14:04:02 k-yukino Exp $
!
#include "parameter.h"

integer function count_number_xi()
  use global_variables
  implicit none

  include "mpif.h"

! local
  integer :: i_count,j_count
  integer,dimension(:),allocatable :: table,new_table,new_temp_table
  integer :: number_accept,temp

  allocate(table(TOTAL_SITE_NUMBER*TOTAL_SITE_NUMBER))
  table=0

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        table((j_count-1)*TOTAL_SITE_NUMBER+i_count)=distance_table(j_count,i_count)
     end do
  end do

! $B$3$3$+$i=EJ#9`$r$^$S$/I,MW$,$"$k!#(B
! $B$$$A$P$s:G=i$N9`$OL5>r7o$GFHN)9`$G$"$k$N$G!"(B

  allocate(new_table(1))         ! $B$^$:$O(B1$B$@$1%a%b%j$r3NJ]$G$9$%!#(B
  new_table(1)=table(1)          ! $B0lHV;O$a$N9`$OFHN)9`$H$7$F:NMQ(B

  number_accept=1                ! $B$"$/$;$W$H$5$l$?$N$O(B1$B8D(B

  do i_count=1,TOTAL_SITE_NUMBER
     temp=0                        ! $BI,?\$G$9(B
     do j_count=1,number_accept
        if( table(i_count)==new_table(j_count) ) then  ! new_table$BCf$K(B
                                                       ! table(i_count)$B$,$"$l$P(B
           temp=temp+1                                 ! temp$B$,%;%C%H$5$l$k(B
        end if
     end do

     if( temp==0 ) then               ! $B$b$7FHN)$J$b$N$,$_$D$+$C$?$i2C$($k(B

! $B=t!9$N=i4|2=(B
        allocate(new_temp_table(number_accept+1))
        new_temp_table=0              ! init
        
        do j_count=1,number_accept
           new_temp_table(j_count)=new_table(j_count)
        end do
        
        new_temp_table(number_accept+1)=table(i_count) 

        deallocate(new_table)
        allocate(new_table(number_accept+1))
        new_table=new_temp_table             ! $B%5%$%:$,F1$8$J$N$G$3$NI=5-$G$h$$(B
        
        deallocate(new_temp_table)
        number_accept=number_accept+1
     end if
  end do

  deallocate(table)
  deallocate(new_table)

  count_number_xi=number_accept

end function count_number_xi



#ifdef DEBUG_COUNT_NUMBER_XI
program main
  use global_variables
  implicit none

  include "mpif.h"

  integer :: result
  integer,external :: count_number_xi 
  character(32) :: name
  integer :: error
! init
  name="main" ; error=0 

  call MPI_INIT(IERROR) 

  call init_global_variables(name,error)
  call error_check(name,error)

  result=count_number_xi()
  write(*,*) "result=",result
  write(*,*) "JIGEN=",JIGEN
  write(*,*) "TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER

  call MPI_FINALIZE(IERROR) 

end program main
#endif /* DEBUG_COUNT_NUMBER_XI */

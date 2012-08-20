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

! ここから重複項をまびく必要がある。
! いちばん最初の項は無条件で独立項であるので、

  allocate(new_table(1))         ! まずは1だけメモリを確保ですぅ。
  new_table(1)=table(1)          ! 一番始めの項は独立項として採用

  number_accept=1                ! あくせぷとされたのは1個

  do i_count=1,TOTAL_SITE_NUMBER
     temp=0                        ! 必須です
     do j_count=1,number_accept
        if( table(i_count)==new_table(j_count) ) then  ! new_table中に
                                                       ! table(i_count)があれば
           temp=temp+1                                 ! tempがセットされる
        end if
     end do

     if( temp==0 ) then               ! もし独立なものがみつかったら加える

! 諸々の初期化
        allocate(new_temp_table(number_accept+1))
        new_temp_table=0              ! init
        
        do j_count=1,number_accept
           new_temp_table(j_count)=new_table(j_count)
        end do
        
        new_temp_table(number_accept+1)=table(i_count) 

        deallocate(new_table)
        allocate(new_table(number_accept+1))
        new_table=new_temp_table             ! サイズが同じなのでこの表記でよい
        
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

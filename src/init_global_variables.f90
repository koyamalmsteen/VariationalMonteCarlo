!
! $Id: init_global_tables.f90,v 1.3 2002/05/31 00:31:10 k-yukino Exp
! k-yukino $
!

#include "parameter.h"

subroutine init_global_variables(name,error)
  use global_variables
  implicit none

  character(32),intent(out) :: name
  integer,intent(out) :: error
!
  integer :: i_count,j_count,k_count,l_count
  integer,external :: count_number_xi
  integer :: number_unit
  character(32) :: buff
  integer,dimension(:),allocatable :: tmp_distance,dummy
  integer,dimension(TOTAL_SITE_NUMBER,2) :: zahyou
! local
  integer :: ier,temp
! init
  name="init_global_variables" ; error=0
! init [local]
  temp=0 ; ier=0 
  d_tilde_up=0 ; d_tilde_down=0
  distance_table=0

! $B%U%i%C%H%P%s%I$N;~$O!"7W;;$G$-$k%5%$%H?t$,Ht$SHt$S$K$J$k!#(B
! $B$h$C$F%5%$%H?t$N;XDj$^$A$,$$$,$J$$$+$r0J2<$G3NG'$9$k!#(B

  if( FLAT_BAND/=0 ) then
     call check_flat_site_number(number_unit,name,error)
     call error_check(name,error)
  end if

  if( JIGEN==1 ) then
! distance$B%F!<%V%k$r:n$k(B
     if( PERIODIC==0 ) then                  ! p $B<~4|E*6-3&>r7o$r;H$&(B
        do j_count=1,TOTAL_SITE_NUMBER
           do i_count=1,TOTAL_SITE_NUMBER
              if( i_count==j_count ) then
                 temp=0
              else if( i_count>j_count ) then
                 temp=min(abs(i_count-j_count),(TOTAL_SITE_NUMBER&
                      & -i_count) +(j_count-1)+1)
              else if( i_count<j_count ) then
                 temp=min(abs(i_count-j_count),(TOTAL_SITE_NUMBER&
                      & -j_count) +(i_count-1)+1)
              end if
           
              distance_table(i_count,j_count)=temp
           end do
        end do
     else if( PERIODIC==1 ) then 
        do j_count=1,TOTAL_SITE_NUMBER
           do i_count=1,TOTAL_SITE_NUMBER
              distance_table(i_count,j_count)=abs(j_count-i_count)
           end do
        end do
     else
        error=-4 ; return
     end if
        
! $BNY@\%F!<%V%k$r:n$k!#(B
     neighbor_table=0

     if( PERIODIC==0 ) then                      ! p $B<~4|E*6-3&>r7o$r;H$&(B
        do i_count=1,TOTAL_SITE_NUMBER
           do j_count=1,TOTAL_SITE_NUMBER
              if( abs(i_count-j_count)==1 .or. &
                   ((i_count==1) .and. (j_count==TOTAL_SITE_NUMBER)) .or. &
                   ((i_count==TOTAL_SITE_NUMBER) .and. (j_count==1)) ) then
                 neighbor_table(i_count,j_count)=1
              end if
           end do
        end do
     else if( PERIODIC==1 ) then                 ! a $B<~4|E*6-3&>r7o$r;H$o$J$$(B
        do i_count=1,TOTAL_SITE_NUMBER
           do j_count=1,TOTAL_SITE_NUMBER
              if( abs(i_count-j_count)==1  ) then
                 neighbor_table(i_count,j_count)=1
              end if
           end do
        end do
     end if

! xaxis$B%"%l%$$r%;%C%H$9$k(B
     do i_count=1,TOTAL_SITE_NUMBER
        xaxis(i_count)=i_count
     end do
  end if

  if( JIGEN==2 ) then
     if( FLAT_BAND==0 ) then             ! $B%U%i%C%H%P%s%I$G$J$$$H$-(B
        call init_global_variables_2dim(zahyou,name,error)
        call error_check(name,error)
     else if( FLAT_BAND==1 ) then        ! lieb
        call init_global_variables_lieb(number_unit,name,error)
        call error_check(name,error)
     else if( FLAT_BAND==2 ) then        ! mielke
        call init_global_variables_mielke(number_unit,name,error)
        call error_check(name,error)
     else if( FLAT_BAND==3 ) then        ! tasaki
        call init_global_variables_tasaki(number_unit,name,error)
        call error_check(name,error)
     end if
  end if

! $BJQJ,%Q%i%a!<%?!<$N<B:]$NCM$rBeF~$9$kA0$K!"G0$N$?$a$K=i4|2=$9$k!#(B
  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        xi_table(i_count,j_count)%way=0.0
        xi_table(i_count,j_count)%xi=0.0
     end do
  end do

  max_number_xi=count_number_xi()

  if( max_number_xi<NUMBER_XI ) then
     write(*,*) "NUMBER_XI error"
     call MPI_FINALIZE(IERROR)
     stop
  end if

  allocate(seq_xi(max_number_xi),stat=ier)
  call stat_check("seq_xi","1",name,error)
  seq_xi=0

  call getarg(1,buff)                                       
  read(buff,'(f8.0)') electric_field

  call getarg(2,buff)
  read(buff,'(f8.0)') xi_field

  call getarg(3,buff)
  read(buff,'(f8.0)') coulomb

  do i_count=1,NUMBER_XI
     call getarg(i_count+3,buff)

     read(buff,'(f8.0)') seq_xi(i_count)
  end do

  allocate(tmp_distance(max_number_xi),stat=ier)     ! max_number_xi
  call stat_check("tmp_distance","1",name,error)
  tmp_distance=0

  allocate(dummy(max_number_xi),stat=ier)            ! max_number_xi
  call stat_check("dummy","1",name,error)
  dummy=0

! $BFHN)$J5wN%$N$_$rNs5s(B
  tmp_distance=-1
  k_count=0
  l_count=0

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,max_number_xi
        if( tmp_distance(j_count)/=distance_table(1,i_count) ) then
           k_count=k_count+1
        end if
     end do

     if( k_count==max_number_xi ) then
        tmp_distance(l_count+1)=distance_table(1,i_count)
        l_count=l_count+1
     end if

     k_count=0
  end do

! slatec 
  call isort(tmp_distance,dummy,max_number_xi,1)

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        xi_table(i_count,j_count)%way=distance_table(i_count,j_count)
     end do
  end do

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        do k_count=1,max_number_xi
           if( xi_table(i_count,j_count)%way==tmp_distance(k_count) ) then
              xi_table(i_count,j_count)%xi=seq_xi(k_count)
           end if
        end do
     end do
  end do

  allocate(independent_distance(max_number_xi),stat=ier)
  call stat_check("independent_distance","init_global_variables",1,ier)

  independent_distance=tmp_distance

! distance_sequence
  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER

        do k_count=1,max_number_xi
           if( independent_distance(k_count) == &
                &distance_table(i_count,j_count) ) then
              distance_sequence(i_count,j_count)=k_count
           end if
        end do

     end do
  end do

  deallocate(tmp_distance,stat=ier)
  call stat_check("tmp_distance","init_global_variables",2,ier)
  deallocate(dummy,stat=ier)
  call stat_check("dummy","init_global_variables",2,ier)

!
! $B0J2<$O(BMYRANK0$B$G$7$+7W;;$7$J$$$N$G!"$3$3$G3NJ](B
!
  allocate(result_s_zero(max_number_xi),stat=ier)
  call stat_check("result_s_zero","cor",1,ier)
  result_s_zero=0

  allocate(result_sxsy_zero(max_number_xi),stat=ier)
  call stat_check("result_sxsy_zero","cor",1,ier)
  result_sxsy_zero=0
  
  allocate(result_sz_zero(max_number_xi),stat=ier)
  call stat_check("result_sz_zero","cor",1,ier)
  result_sz_zero=0
  
  allocate(result_charge_zero(max_number_xi),stat=ier)
  call stat_check("result_charge_zero","cor",1,ier)
  result_charge_zero=0

  allocate(result_super_zero(max_number_xi),stat=ier)
  call stat_check("result_super_zero","cor",1,ier)
  result_super_zero=0
!
  result_pole_zero=0


  result_pole0=0 ; result_pole1=0 ; result_pole2=0 

  

!
! $B%(%i!<%P!<(B($B%(%M%k%.!<(B)
!
  splited_energy=0
  buff_splited=0
  splited_h0=0 ; splited_h1=0 ; splited_h2=0 
  splited_h3=0 ; splited_h4=0 ; splited_h5=0 
  buff5=0

end subroutine init_global_variables



#ifdef DEBUG_INIT_GLOBAL_VARIABLES
program main
  use global_variables
  implicit none
  
  character(32) :: name
  integer :: error
! local
  integer,dimension(:),allocatable :: aaa1,aaa2,aaa3,aaa4,aaa5
! local
  integer :: i_count,j_count
! init
  name="main" ; error=0

  call MPI_INIT(IERROR) 

  call init_global_variables(name,error)
  call error_check(name,error)

  write(*,*) "aaa1"
  allocate(aaa1(22))
  write(*,*) "aaa2"
  allocate(aaa2(22))
  write(*,*) "aaa3"
  allocate(aaa3(22))
  write(*,*) "aaa4"
  allocate(aaa4(22))
  write(*,*) "aaa5"
  allocate(aaa5(23))
  write(*,*) "seikou"

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_INIT_GLOBAL_VARIABLES */

!
! $Id: cor_zero.f90,v 1.5 2004/03/10 03:04:02 k-yukino Exp $
!
#include "parameter.h"

subroutine cor_zero(name,error)
  use global_variables
  implicit none

  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count,j_count
  integer :: ier
  real(8) :: tmp_result
! init
  name="cor_zero" ; error=0
  ier=0
  tmp_result=0

!        
! <$BBh(B0$B6a;w(B|O|$BBh(B0$B6a;w(B>$B$r7W;;(B
!
  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=i_count,TOTAL_SITE_NUMBER
        tmp_result=0

        call sz_zero(i_count,j_count,tmp_result,name,error)
        call error_check(name,error)
        result_sz_zero(distance_sequence(i_count,j_count))=&
             result_sz_zero(distance_sequence(i_count,j_count)) &
             +tmp_result/dble(TOTAL_SITE_NUMBER)

     end do
  end do

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=i_count,TOTAL_SITE_NUMBER
        tmp_result=0
        call charge_zero(i_count,j_count,tmp_result,name,error)
        call error_check(name,error)
        result_charge_zero(distance_sequence(i_count,j_count))=&
             result_charge_zero(distance_sequence(i_count,j_count)) &
             +tmp_result/dble(TOTAL_SITE_NUMBER)
     end do
  end do

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=i_count,TOTAL_SITE_NUMBER
        tmp_result=0
        call sxsy_zero(i_count,j_count,tmp_result,name,error)
        call error_check(name,error)
        result_sxsy_zero(distance_sequence(i_count,j_count))=&
             result_sxsy_zero(distance_sequence(i_count,j_count)) &
             +tmp_result/dble(TOTAL_SITE_NUMBER)
     end do
  end do

  if( JIGEN==2 ) then
     do i_count=1,TOTAL_SITE_NUMBER
        do j_count=i_count,TOTAL_SITE_NUMBER
           tmp_result=0
           call super_zero(i_count,j_count,tmp_result,name,error)
           call error_check(name,error)

           result_super_zero(distance_sequence(i_count,j_count))=&
                result_super_zero(distance_sequence(i_count,j_count)) &
                +tmp_result/dble(TOTAL_SITE_NUMBER)
        end do
     end do
  end if
  
! $BB-$79g$o$;$F(Bs$B$r:n$k(B
  do i_count=1,max_number_xi
     result_s_zero(i_count)=result_sxsy_zero(i_count)&
          +result_sz_zero(i_count)
  end do

!
! $BEE2YL)EY$N7W;;(B($BAj4X4X?t$G$O$J$$!#$h$C$F%5%$%H?t$NBg$-$5$r;}$D(B)
!
!
! charge_density
  do i_count=1,TOTAL_SITE_NUMBER
     tmp_result=0
     call charge_density_zero(i_count,tmp_result,name,error)
     call error_check(name,error)
     result_charge_density_zero(i_count)=&
          result_charge_density_zero(i_count) &
          +tmp_result
  end do

!
! pole$B$N7W;;(B
!
!
! pole
  do i_count=1,TOTAL_SITE_NUMBER
     tmp_result=0
     call pole_zero(i_count,tmp_result,name,error)
     call error_check(name,error)
     result_pole_zero=result_pole_zero+tmp_result
  end do

end subroutine cor_zero



#ifdef DEBUG_COR_ZERO
program main
  use global_variables
  implicit none

  include "mpif.h"

  character(32) :: name
  integer :: error
! local
  integer,dimension(2) :: ic
  real(8) :: zero_approx_energy
  integer :: hf_iteration_result
! init
  name="main" ; error=0
  unitary_d_up=0 ; unitary_d_down=0 
  zero_approx_energy=0 ; hf_iteration_result=0

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

  call random_init(ic,name,error)
  call error_check(name,error)

  call zero_approx(zero_approx_energy,hf_iteration_result,name,error)
  call error_check(name,error)

  call cor_zero(name,error)
  call error_check(name,error)

  write(*,*) "JIGEN=",JIGEN,"TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "INIT_WAVE_VECTOR=",INIT_WAVE_VECTOR
  write(*,*) 
  write(*,*) "result_sz_zero=",result_sz_zero
  write(*,*) "result_sxsy_zero=",result_sxsy_zero
  write(*,*) "result_s_zero=",result_s_zero
  write(*,*) "result_charge_zero=",result_charge_zero
  write(*,*) "result_super_zero=",result_super_zero
  write(*,*) "result_charge_density_zero=",result_charge_density_zero

  call MPI_FINALIZE(IERROR)


end program main
#endif /* DEBUG_COR_ZERO */

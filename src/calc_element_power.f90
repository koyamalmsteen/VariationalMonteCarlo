!
! $Id: calc_element_power.f90,v 1.6 2003/08/14 09:26:51 k-yukino Exp $
!

!
! すでに集められた<alpha|を使って、パワーをかける。
!

!
! E_pow1=<Psi_1 |H|Psi_1 >/<Psi_1 |Psi_1 >
!  =C1^2 <phi_0 |H^3 |phi_0 >+2C1 <phi_0 |H^2 |phi_0 >+<phi_0 |H|phi_0 >
!   --------------------------------------------------------------------
!   C1^2 <phi_0 |H^2 |phi_0 >+2C1 <phi_0 |H|phi_0 >+<phi_0 |phi_0 >  
!

!
! E_pow2=<Psi_2 |H|Psi_2 >/<Psi_2 |Psi_2 >
!

#include "parameter.h"

subroutine calc_element_power(result_pow0,result_pow1,result_pow2,&
                              pow1_c1,pow2_c1,pow2_c2,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(out) :: result_pow0,result_pow1,result_pow2,pow1_c1,pow2_c1,pow2_c2
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count,j_count
  real(8) :: temp0,temp1,temp2,temp3,temp4,temp5
  real(8) :: total_temp0,total_temp1,total_temp2,total_temp3,total_temp4,&
             total_temp5
  real(8) :: result3_4_5
  integer,dimension(TOTAL_UP_ELECTRON) :: left_vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: left_vector_down
  real(8),external :: evaluate_pow1,evaluate_pow2
  real(8) :: temp_result,temp_c1,temp_c2,min_pow1,min_pow2
  integer :: idou_moto_electron_spin
  real(8) :: pow_fn
! MPI関連
  integer :: COMM,NUMBER_PE,MYRANK
  integer :: IERROR 
! init
  result_pow1=0 ; result_pow2=0
  name="calc_element_power" ; error=0
  temp0=0 ; temp1=0 ; temp2=0 ; temp3=0 ; temp4=0 ; temp5=0
  total_temp0=0 ; total_temp1=0 ; total_temp2=0 ; total_temp3=0
  total_temp4=0 ; total_temp5=0
  pow1_c1=0 ; pow2_c1=0 ; pow2_c2=0
  temp_c1=0 ; temp_c2=0
  idou_moto_electron_spin=0
  min_pow1=0 ; min_pow2=0
  pow_fn=0

  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERROR)

  COMM=MPI_COMM_WORLD
  call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMBER_PE,IERROR)

!
! 初期の試行電子状態をセットする。
!
  do i_count=1,TOTAL_UP_ELECTRON
     left_vector_up(i_count)=local_vector1%vector_up(i_count)
  end do

  do i_count=1,TOTAL_DOWN_ELECTRON
     left_vector_down(i_count)=local_vector1%vector_down(i_count)
  end do

!
! メインルーチン
!

!!
!! 初期電子状態に関するエネルギー期待値だけ先に計算しておく
!!

  temp0=0
  temp1=0
  temp2=0
  temp3=0
  temp4=0
  temp5=0

  temp0=local_vector1%alpha_psi1**2
  temp1=local_vector1%alpha_h_psi1*local_vector1%alpha_psi1
  temp2=abs(local_vector1%alpha_h_psi1)**2
  call calc_element_pow3_dash(left_vector_up,left_vector_down,1,&
                              temp3,result3_4_5,name,error)
  call error_check(name,error)
! result3_4_5は<alpha|HH|psi0>なのでストアしておく(相関関数のときに使う)
  local_vector1%alpha_h_h_psi1=result3_4_5

  if( POWER==2 ) then
     temp4=abs(result3_4_5)**2

     call calc_element_pow5_dash(left_vector_up,left_vector_down,&
                                 result3_4_5,temp5,name,error)
     call error_check(name,error)
  end if

  total_temp0=0
  total_temp1=0
  total_temp2=0
  total_temp3=0
  total_temp4=0
  total_temp5=0
  
  total_temp0=total_temp0+temp0/temp0
  total_temp1=total_temp1+temp1/temp0
  total_temp2=total_temp2+temp2/temp0
  total_temp3=total_temp3+temp3/temp0
  total_temp4=total_temp4+temp4/temp0
  total_temp5=total_temp5+temp5/temp0

!!
!! (初期電子配置に関する計算終了)
!!


! 以下は必須です。 
  idou_moto_electron_spin=0

!
! 以下は2番目以降の電子状態での期待値計算
!
  do i_count=1,MAX_MONTECARLO_SAMPLE-1

     if( idou_moto_electron_spin/=3 ) then    ! 前と同じ電子状態ならさぼる

        temp0=0 ; temp1=0 ; temp2=0
        temp0=local_vector2(i_count)%alpha_psi1**2             ! pow0
        temp1=local_vector2(i_count)%alpha_h_psi1 &            ! pow1
              *local_vector2(i_count)%alpha_psi1
        temp2=abs(local_vector2(i_count)%alpha_h_psi1)**2      ! pow2

! pow3
        temp3=0

        if( HEEBRICE==0 ) then
           call calc_element_pow3_dash(left_vector_up,left_vector_down,&
                i_count+1,temp3,result3_4_5,name,error)
           call error_check(name,error)
! 
! result3_4_5は<alpha|HH|psi0>で相関関数のときに使うので、ストアしておく
!
           local_vector2(i_count)%alpha_h_h_psi1=result3_4_5
           
        else if( HEEBRICE==1 ) then
           call calc_element_pow3_dash_heebrice(left_vector_up,&
                left_vector_down,i_count+1,&
                temp3,result3_4_5,name,error)
           call error_check(name,error)
        end if

        if( POWER==2 ) then
! pow4
           temp4=0
           temp4=abs(result3_4_5)**2

! pow5
           temp5=0
           if( HEEBRICE==0 ) then
              call calc_element_pow5_dash(left_vector_up,left_vector_down,&
                   result3_4_5,temp5,name,error)
              call error_check(name,error)
           else if( HEEBRICE==1 ) then
              call calc_element_pow5_dash_heebrice(left_vector_up,&
                   left_vector_down,&
                   result3_4_5,temp5,name,error)
              call error_check(name,error)
           end if

        end if
     end if

     total_temp0=total_temp0+temp0/temp0
     total_temp1=total_temp1+temp1/temp0
     total_temp2=total_temp2+temp2/temp0
     total_temp3=total_temp3+temp3/temp0
     total_temp4=total_temp4+temp4/temp0
     total_temp5=total_temp5+temp5/temp0
/*
     powpow(i_count+1)%temp0=temp0/temp0
     powpow(i_count+1)%temp1=temp1/temp0
     powpow(i_count+1)%temp2=temp2/temp0
     powpow(i_count+1)%temp3=temp3/temp0
     powpow(i_count+1)%temp4=temp4/temp0
     powpow(i_count+1)%temp5=temp5/temp0
*/
!
! 次の電子状態を呼ぶ
! local_vector2(i_count)%idou_spin==3のときはなにしない
!
     if( local_vector2(i_count+1)%idou_spin==1 ) then           ! ↑が動いたら
        left_vector_up( local_vector2(i_count+1)%idou_moto ) = &
                                            local_vector2(i_count+1)%idou_saki
        idou_moto_electron_spin=local_vector2(i_count+1)%idou_spin
     else if( local_vector2(i_count+1)%idou_spin==2 ) then      ! ↓なら
        left_vector_down( local_vector2(i_count+1)%idou_moto ) = &
                                            local_vector2(i_count+1)%idou_saki
        idou_moto_electron_spin=local_vector2(i_count+1)%idou_spin
     end if

  end do

  call MPI_REDUCE(total_temp0,temp0,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,&
                  IERROR)
  call MPI_REDUCE(total_temp1,temp1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,&
                  IERROR)
  call MPI_REDUCE(total_temp2,temp2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,&
                  IERROR)
  call MPI_REDUCE(total_temp3,temp3,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,&
                  IERROR)
  call MPI_REDUCE(total_temp4,temp4,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,&
                  IERROR)
  call MPI_REDUCE(total_temp5,temp5,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,&
                  IERROR)

  if( MYRANK==0 ) then
     temp0=temp0/MAX_MONTECARLO_SAMPLE
     temp1=temp1/MAX_MONTECARLO_SAMPLE
     temp2=temp2/MAX_MONTECARLO_SAMPLE
     temp3=temp3/MAX_MONTECARLO_SAMPLE
     temp4=temp4/MAX_MONTECARLO_SAMPLE
     temp5=temp5/MAX_MONTECARLO_SAMPLE

     write(*,*) "temp0=",temp0
     write(*,*) "temp1=",temp1
     write(*,*) "temp2=",temp2
     write(*,*) "temp3=",temp3
     write(*,*) "temp4=",temp4
     write(*,*) "temp5=",temp5
     write(*,*)

     temp_result=0 ; min_pow1=99999

! 最小値化する(part1)
     do i_count=-100,100
        temp_c1=0 ; temp_c1=dble(i_count)*0.1
        temp_result=evaluate_pow1(temp0,temp1,temp2,temp3,temp_c1)

        if( min_pow1>temp_result ) then
           min_pow1=temp_result
           pow1_c1=temp_c1
        end if
     end do

     if( POWER==2 ) then
        temp_result=0 ; min_pow2=99999

! 最小値化する(part2)
        if( POWER==2 ) then
           do i_count=-100,100
              do j_count=-100,100
                 temp_c1=0 ; temp_c2=0
                 temp_c1=dble(i_count)*0.1
                 temp_c2=dble(j_count)*0.1
                 temp_result=evaluate_pow2(temp0,temp1,temp2,temp3,&
                                           temp4,temp5,temp_c1,temp_c2)
                 if( min_pow2>temp_result ) then
                    min_pow2=temp_result
                    pow2_c1=temp_c1
                    pow2_c2=temp_c2
                 end if
              end do
           end do
        end if
     end if

     result_pow0=temp1/temp0            ! 注意!! 
     result_pow1=min_pow1
     result_pow2=min_pow2
  end if

end subroutine calc_element_power



real(8) function evaluate_pow1(temp0,temp1,temp2,temp3,C1)
  implicit none

  real(8),intent(in) :: temp0,temp1,temp2,temp3,C1
  
  evaluate_pow1=( C1*C1*temp3+2*C1*temp2+temp1 )/&
                ( C1*C1*temp2+2*C1*temp1+temp0 )

end function evaluate_pow1



real(8) function evaluate_pow2(temp0,temp1,temp2,temp3,temp4,temp5,C1,C2)
  implicit none

  real(8),intent(in) :: temp0,temp1,temp2,temp3,temp4,temp5,C1,C2

  evaluate_pow2=( C2*C2*temp5+C1*C2*temp4+(C1*C1+2*C2+C1*C2)*temp3 &
                  +2*C1*temp2+temp1 )/&
                ( C2*C2*temp4+C1*C2*temp3+(C1*C1+2*C2+C1*C2)*temp2 &
                  +2*C1*temp1+temp0 )

end function evaluate_pow2





#ifdef DEBUG_CALC_ELEMENT_POWER
program main
  use global_variables 
  implicit none
  
  include "mpif.h"
  
  real(8) :: result
  character(32) :: name
  integer :: error
  integer :: i_count,j_count,tmp,ier,IERROR
  integer,external :: count_number_xi
!
  name="main" ; error=0
  
  call MPI_INIT(IERROR)

  call init_distance()

  tmp=count_number_xi()

  allocate(vtable(tmp),stat=ier)
  call stat_check("vtable","main",1,ier)
  do i_count=1,tmp
     vtable(i_count)%way=0
     vtable(i_count)%xi=0
  end do

! set variation parameter
  do i_count=1,tmp

  end do

  call calc_element_power(result,name,error)  
  call error_check(name,error)

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CALC_ELEMENT_POWER */

!
! $Id: global_variables.f90,v 1.20 2004/03/10 03:04:02 k-yukino Exp k-yukino $
!
#include "parameter.h"

module global_variables
  implicit none

  type state_vector_structure
     integer,dimension(TOTAL_UP_ELECTRON) :: vector_up
     integer,dimension(TOTAL_DOWN_ELECTRON) :: vector_down
     real(8) :: alpha_psi0
     real(8) :: alpha_psi1
     real(8) :: alpha_h_psi1
     real(8) :: alpha_h_h_psi1
     real(8) :: alpha_h_h_h_psi1
     real(8) :: lambda
  end type state_vector_structure

  type vparameter
     integer :: way
     real(8) :: xi
  end type vparameter

!  type(state_vector_structure),dimension(:),allocatable :: state_vector
  type(state_vector_structure),dimension(MAX_MONTECARLO_SAMPLE) :: state_vector
! BCAST$BI,MW(B(MYRANK==0$B$,=i4|2=$7$F(BBCAST$B$9$k(B)
  real(8),dimension(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER) :: unitary_d_up,&
       unitary_d_down
  integer,dimension(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER) :: distance_table
  integer,dimension(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER) :: neighbor_table
  integer,dimension(TOTAL_SITE_NUMBER,4) :: neighbor_table2
  integer,dimension(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER) :: distance_sequence
  type(vparameter),dimension(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER) :: xi_table
  integer :: max_number_xi
! BCAST$BITI,MW(B($B$?$@$N(Bwork$BNN0h(B)
  integer,dimension(TOTAL_SITE_NUMBER) :: global_site_table_up,&
       global_site_table_up2
  integer,dimension(TOTAL_SITE_NUMBER) :: global_site_table_down,&
       global_site_table_down2
  real(8),dimension(TOTAL_UP_ELECTRON,TOTAL_UP_ELECTRON) :: d_tilde_up,&
       d_tilde_up_inverse
  real(8),dimension(TOTAL_DOWN_ELECTRON,TOTAL_DOWN_ELECTRON) :: d_tilde_down,&
       d_tilde_down_inverse
  real(8),dimension(TOTAL_UP_ELECTRON,TOTAL_UP_ELECTRON) :: after_lu_up
  real(8),dimension(TOTAL_DOWN_ELECTRON,TOTAL_DOWN_ELECTRON) :: after_lu_down
  integer,dimension(TOTAL_UP_ELECTRON) :: ipiv_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: ipiv_down
  real(8),dimension(TOTAL_UP_ELECTRON,TOTAL_UP_ELECTRON) :: d_tilde_up_work
  real(8),dimension(TOTAL_DOWN_ELECTRON,TOTAL_DOWN_ELECTRON) :: d_tilde_down_work
! BCAST$BITI,MW(B(MYRANK==0$B$N$_;}$C$F$$$l$P$h$$$,!"JXMx$J$N$G%0%m!<%P%kJQ?t$K$7$?(B)
  real(8),dimension(:),allocatable :: seq_xi
  integer,dimension(:),allocatable :: independent_distance
! $BAj4X4X?t(B
  real(8),dimension(:),allocatable :: result_sz_zero,result_s_zero,result_charge_zero,result_sxsy_zero,result_super_zero
  real(8),dimension(:),allocatable :: result_sz0,result_sz1,result_sz2
  real(8),dimension(:),allocatable :: result_sxsy0,result_sxsy1,result_sxsy2
  real(8),dimension(:),allocatable :: result_s0,result_s1,result_s2
  real(8),dimension(:),allocatable :: result_super0,result_super1,result_super2
  real(8),dimension(:),allocatable :: result_charge0,result_charge1,result_charge2
! $B0J2<$OAj4X4X?t$G$J$$$N$G(BTOTAL_SITE_NUMBER$B$@$13NJ](B(allocate_global_variables
! $B$G=i4|2=$9$k!#(B)
  real(8),dimension(TOTAL_SITE_NUMBER) :: result_charge_density0,&
       result_charge_density1,result_charge_density2
  real(8) :: result_pole0,result_pole1,result_pole2
  real(8),dimension(TOTAL_SITE_NUMBER) ::  result_charge_density_zero
  real(8) :: result_pole_zero
! $B7W;;;~4V(B($B$3$3$O$J$<$+(Breal$B$GDj5A!#M}M3$OK:$l$?(B)
  real(8) :: time_start,time_end
! correct_montecarlo_sample$B$G;HMQ(B
  real(8) :: stored_alpha_psi0,stored_alpha_psi1,stored_lambda
  real(8) :: stored_input_up,stored_input_down
! exp_code$BMQ(B
  character(8) :: date
  character(10) :: time
  character(5) :: zone
  integer,dimension(8) :: values
  integer,dimension(TOTAL_SITE_NUMBER) :: xaxis
!
! MPI$B4XO"(B
! 
  integer :: NUMBER_PE,MYRANK,IERROR
!
! $B%(%i!<%P!<(B($B%(%M%k%.!<(B)
!
  real(8),dimension(5,MAX_MONTECARLO_SAMPLE/5) :: splited_energy,buff_splited
  real(8),dimension(5) :: splited_h0,splited_h1,splited_h2,splited_h3,&
       splited_h4,splited_h5,buff5
!
! $B%(%i!<%P!<(B
!
  real(8),dimension(:,:),allocatable :: buff_cor5 
! Sz
  real(8),dimension(:,:),allocatable :: splited_sz0,splited_sz1,splited_sz2
! charge
  real(8),dimension(:,:),allocatable :: splited_charge0,splited_charge1,splited_charge2
! charge_density
  real(8),dimension(:,:),allocatable :: splited_charge_density0,splited_charge_density1,splited_charge_density2
! pole
  real(8),dimension(5) :: splited_pole0,splited_pole1,splited_pole2
! sxsy
  real(8),dimension(:,:),allocatable :: splited_sxsy0,splited_sxsy1,splited_sxsy2
! s
  real(8),dimension(:,:),allocatable :: splited_s0,splited_s1,splited_s2
! super
  real(8),dimension(:,:),allocatable :: splited_super0,splited_super1,splited_super2
! getarg
  real(8) :: coulomb,electric_field,xi_field

end module global_variables



#ifdef DEBUG_GLOBAL_VARIABLES
program main
  use global_variables
  implicit none

end program main
#endif /* DEBUG_GLOBAL_VARIABLES */

!
! $Id: vmc.f90,v 1.7 2004/03/10 03:03:55 k-yukino Exp k-yukino $
!
#include "parameter.h"

program main
  use global_variables
  implicit none

  include "mpif.h"

! local
  real(8),dimension(5) :: error_pow0,error_pow1,error_pow2
  real(8) :: lowest_pow0,lowest_pow1,lowest_pow2
  real(8) :: highest_pow0,highest_pow1,highest_pow2
  integer :: sample_number
  real(8) :: zero_approx_energy
  integer :: hf_iteration_result
  character(32) :: name
  integer :: ier,error,i_count,j_count
  real(8) :: tmp,dummy0,dummy1,dummy2
! MPI関連
  integer :: POSITION,PACK_SIZE,TMP_PACK_SIZE_DOUBLE,&
       TMP_PACK_SIZE_INTEGER,TMP_PACK_SIZE_NEWTYPE
  real(8),dimension(20*(TOTAL_SITE_NUMBER**2)) :: BUFFER
  integer,dimension(2) :: BLOCKLENS,OLD_TYPES,INDICES
  integer :: NEWTYPE,COUNT
! local
  integer,dimension(TOTAL_UP_ELECTRON) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: vector_down
  real(8) :: result_pow0,result_pow1,result_pow2
  real(8) :: pow1_c0,pow1_c1
  real(8) :: pow2_c0,pow2_c1,pow2_c2
  real(8),dimension(MAX_MONTECARLO_SAMPLE) :: tmp_real,total_real
  integer,dimension(TOTAL_UP_ELECTRON) :: tmp_vector_up,total_vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: tmp_vector_down,total_vector_down
  integer,dimension(2) :: ic
! local
  real(8) :: denomi_pow1,denomi_pow2,buff
  real(8),dimension(5) :: denomi_splited_pow1,denomi_splited_pow2
!
! cor
!
  real(8),dimension(:),allocatable :: tmp_sz0,tmp_charge0,tmp_sxsy0,tmp_super0
  real(8),dimension(TOTAL_SITE_NUMBER) :: tmp_charge_density0
!
  integer :: need_memory 
!
  real(8),dimension(:),allocatable :: sz_nume1,sz_nume2,sz_nume3_1,sz_nume3_2,sz_nume4,sz_nume5
  real(8),dimension(:),allocatable :: charge_nume1,charge_nume2,charge_nume3_1,charge_nume3_2,charge_nume4,charge_nume5
  real(8),dimension(:),allocatable :: charge_density_nume1,charge_density_nume2,charge_density_nume3_1,charge_density_nume3_2
  real(8),dimension(:),allocatable :: charge_density_nume4,charge_density_nume5
  real(8),dimension(:),allocatable :: sxsy_nume1,sxsy_nume2,sxsy_nume3_1,sxsy_nume3_2,sxsy_nume4,sxsy_nume5
  real(8),dimension(:),allocatable :: super_nume1,super_nume2,super_nume3_1,super_nume3_2,super_nume4,super_nume5
  real(8),dimension(:),allocatable :: buff_cor,buff_cd
  real(8),dimension(:,:),allocatable :: buff_splited_cor,buff_splited_cd  
!
  real(8) :: alpha_psi0,alpha_psi1,alpha_h_psi1,alpha_h_h_psi1&
       ,alpha_h_h_h_psi1,lambda
  real(8) :: tmp_h0,tmp_h1,tmp_h2,tmp_h3,tmp_h4,tmp_h5
  real(8) :: h0,h1,h2,h3,h4,h5
! error_bar
  real(8),dimension(:),allocatable :: lowest_sz0,lowest_sz1,lowest_sz2
  real(8),dimension(:),allocatable :: lowest_charge0,lowest_charge1,lowest_charge2
  real(8),dimension(:),allocatable :: lowest_sxsy0,lowest_sxsy1,lowest_sxsy2
  real(8),dimension(:),allocatable :: lowest_s0,lowest_s1,lowest_s2
  real(8),dimension(:),allocatable :: lowest_super0,lowest_super1,lowest_super2
  real(8),dimension(:),allocatable :: lowest_charge_density0,lowest_charge_density1,lowest_charge_density2
  real(8) :: lowest_pole0,lowest_pole1,lowest_pole2 
  real(8),dimension(:),allocatable :: highest_sz0,highest_sz1,highest_sz2
  real(8),dimension(:),allocatable :: highest_charge0,highest_charge1,highest_charge2
  real(8),dimension(:),allocatable :: highest_sxsy0,highest_sxsy1,highest_sxsy2
  real(8),dimension(:),allocatable :: highest_s0,highest_s1,highest_s2
  real(8),dimension(:),allocatable :: highest_super0,highest_super1,highest_super2
  real(8),dimension(:),allocatable :: highest_charge_density0,highest_charge_density1,highest_charge_density2
  real(8) :: highest_pole0,highest_pole1,highest_pole2 
!
  real(8) :: alpha_p_psi1,alpha_p_p_psi1,alpha_h_p_psi1
  real(8) :: bunsi,bunbo,xi_field0 
! init [ MPI ]
  NUMBER_PE=0 ; MYRANK=0 ; IERROR=0
  PACK_SIZE=0 ; TMP_PACK_SIZE_DOUBLE=0 ; TMP_PACK_SIZE_INTEGER=0
  TMP_PACK_SIZE_NEWTYPE=0
! init [ pow0,pow ]
  result_pow0=0 ; result_pow1=0 ; result_pow2=0
! init [ local ]
  name="vmc" ; error=0
  BUFFER=0 ; POSITION=0
  pow1_c0=0 ; pow1_c1=0 ; pow2_c0=0 ; pow2_c1=0 ; pow2_c2=0
  zero_approx_energy=0 ; hf_iteration_result=0 
  ic=0
  vector_up=0 ; vector_down=0 
  alpha_psi0=0 ; alpha_psi1=0 ; alpha_h_psi1=0 ; alpha_h_h_psi1=0
  alpha_h_h_h_psi1=0 ; lambda=0
  tmp_h0=0 ; tmp_h1=0 ; tmp_h2=0 ; tmp_h3=0 ; tmp_h4=0 ; tmp_h5=0
  h0=0 ; h1=0 ; h2=0 ; h3=0 ; h4=0 ; h5=0

  denomi_pow1=0 ; denomi_pow2=0
  buff=0
  need_memory=0 ; sample_number=0
!
  denomi_splited_pow1=0 ; denomi_splited_pow2=0 ; buff5=0
!
  error_pow0=0
  lowest_pow0=0 ; highest_pow0=0

  lowest_pole0=0 ; lowest_pole1=0 ; lowest_pole2=0 
  highest_pole0=0 ; highest_pole1=0 ; highest_pole2=0 
!
  alpha_p_psi1=0 ; alpha_p_p_psi1=0 ; alpha_h_p_psi1=0
  bunsi=0 ; bunbo=0 ; xi_field0=0

! MPIの準備
  call MPI_INIT(IERROR) 
  call error_check("MPI_INIT",IERROR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERROR) 
  call error_check("MPI_COMM_RANK",IERROR)

  call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMBER_PE,IERROR)
  call error_check("MPI_COMM_SIZE",IERROR)

  if( MYRANK==0 ) then
! パラメーターの初期化
     call init_global_variables(name,error)
     call error_check(name,error)

     call parameter_check(name,error) ! 取得したパラメーターのチェックルーチン
     call error_check(name,error)

     time_start=MPI_WTIME()

! まずexp_codeを決定するために現在の日付などを得る。
! このvaluesを1〜8まで並べて主キーとする。
     call date_and_time(date,time,zone,values)

! これを乱数の種にする
     do i_count=1,4
        ic(1)=ic(1)+values(i_count)
     end do
     do i_count=5,8
        ic(2)=ic(2)+values(i_count)
     end do
! デバッグ用
!     ic=0

! メモリの使用量を見積もる(もし多かったら中断する)
     call estimate_use_memory(need_memory,name,error)
     call error_check(name,error)

! 第0近似の実行(INIT_WAVE_VECTOR = 0;Hartree-Fock , 1 ;Plane wave )
! (MYRANK==0のみが実行し、結果をbcastする)
     call zero_approx(zero_approx_energy,hf_iteration_result,name,error)
     call error_check(name,error)

! デバッグ用 - 相関関数の計算 (第0近似により<phi|O|phi>を計算している)
     call cor_zero(name,error)
     call error_check(name,error)

! デバッグ用 - エネルギー期待値の計算 (この計算は結構きつい。大きいサイズでは
! 無理)
#ifdef DEBUG     
     if( POWER==1 .or. POWER==2 ) then
        call energy_zero_power(name,error)
        call error_check(name,error)
     end if
#endif /* DEBUG */
  end if

! 乱数の種を全てのプロセッサに渡す。
  call MPI_BCAST(ic,2,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
  call error_check("MPI_BCAST",IERROR)

! 乱数の初期化
  call random_init(ic,name,error)
  call error_check(name,error)

! MPIで構造体を扱うための準備
  COUNT=2
  BLOCKLENS(1)=1
  BLOCKLENS(2)=1
  OLD_TYPES(1)=MPI_INTEGER
  OLD_TYPES(2)=MPI_DOUBLE_PRECISION
! INDICESの設定(アドレスの設定)
  call MPI_ADDRESS(xi_table%way,INDICES(1),IERROR)
  call MPI_ADDRESS(xi_table%xi,INDICES(2),IERROR)
! INDICES(1)=0としたアドレスにしておくと便利なので置き換える。
  INDICES(2)=INDICES(2)-INDICES(1)
  INDICES(1)=0

  call MPI_TYPE_STRUCT(COUNT,BLOCKLENS,INDICES,OLD_TYPES,NEWTYPE,IERROR)
  call MPI_TYPE_COMMIT(NEWTYPE,IERROR)

! 各プロセッサーに「ユニタリー行列」を配る
! 通信回数を減らすためにMPI_PACKを使う
  if( NUMBER_PE/=1 ) then
     POSITION=0                  ! こうしとかなループした時に問題ある

! unitary_d_up
     call MPI_PACK_SIZE(TOTAL_SITE_NUMBER**2,MPI_DOUBLE_PRECISION,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_DOUBLE,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_DOUBLE
! unitary_d_down
     call MPI_PACK_SIZE(TOTAL_SITE_NUMBER**2,MPI_DOUBLE_PRECISION,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_DOUBLE,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_DOUBLE
! distance_table
     call MPI_PACK_SIZE(TOTAL_SITE_NUMBER**2,MPI_INTEGER,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_INTEGER,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_INTEGER
! neighbor_table 
     call MPI_PACK_SIZE(TOTAL_SITE_NUMBER**2,MPI_INTEGER,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_INTEGER,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_INTEGER
! neighbor_table2
     call MPI_PACK_SIZE(4*TOTAL_SITE_NUMBER,MPI_INTEGER,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_INTEGER,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_INTEGER
! distance_sequence
     call MPI_PACK_SIZE(TOTAL_SITE_NUMBER**2,MPI_INTEGER,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_INTEGER,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_INTEGER
! max_number_xi
     call MPI_PACK_SIZE(1,MPI_INTEGER,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_INTEGER,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_INTEGER
! xi_table
     call MPI_PACK_SIZE(TOTAL_SITE_NUMBER**2,NEWTYPE,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_NEWTYPE,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_NEWTYPE
! electric_field
     call MPI_PACK_SIZE(1,MPI_DOUBLE_PRECISION,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_DOUBLE,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_DOUBLE
! xi_field
     call MPI_PACK_SIZE(1,MPI_DOUBLE_PRECISION,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_DOUBLE,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_DOUBLE
! coulomb
     call MPI_PACK_SIZE(1,MPI_DOUBLE_PRECISION,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_DOUBLE,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_DOUBLE
! xaxis
     call MPI_PACK_SIZE(TOTAL_SITE_NUMBER,MPI_INTEGER,&
          MPI_COMM_WORLD,TMP_PACK_SIZE_INTEGER,IERROR)
     PACK_SIZE=PACK_SIZE+TMP_PACK_SIZE_INTEGER


     if ( MYRANK==0 ) then
        call MPI_PACK(unitary_d_up,TOTAL_SITE_NUMBER**2,MPI_DOUBLE_PRECISION,&
             BUFFER,PACK_SIZE,POSITION,MPI_COMM_WORLD,IERROR)
        call MPI_PACK(unitary_d_down,TOTAL_SITE_NUMBER**2,&
             MPI_DOUBLE_PRECISION,BUFFER,PACK_SIZE,POSITION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_PACK(distance_table,TOTAL_SITE_NUMBER**2,&
             MPI_INTEGER,BUFFER,PACK_SIZE,POSITION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_PACK(neighbor_table,TOTAL_SITE_NUMBER**2,&
             MPI_INTEGER,BUFFER,PACK_SIZE,POSITION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_PACK(neighbor_table2,TOTAL_SITE_NUMBER*4,&
             MPI_INTEGER,BUFFER,PACK_SIZE,POSITION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_PACK(distance_sequence,TOTAL_SITE_NUMBER**2,&
             MPI_INTEGER,BUFFER,PACK_SIZE,POSITION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_PACK(max_number_xi,1,MPI_INTEGER,BUFFER,PACK_SIZE,POSITION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_PACK(electric_field,1,MPI_DOUBLE_PRECISION,BUFFER,PACK_SIZE,POSITION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_PACK(xi_field,1,MPI_DOUBLE_PRECISION,BUFFER,PACK_SIZE,POSITION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_PACK(coulomb,1,MPI_DOUBLE_PRECISION,BUFFER,PACK_SIZE,POSITION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_PACK(xaxis,TOTAL_SITE_NUMBER,MPI_INTEGER,BUFFER,PACK_SIZE,POSITION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_BCAST(BUFFER,PACK_SIZE,MPI_PACKED,0,MPI_COMM_WORLD,IERROR)
     else
        call MPI_BCAST(BUFFER,PACK_SIZE,MPI_PACKED,0,MPI_COMM_WORLD,IERROR)
        call MPI_UNPACK(BUFFER,PACK_SIZE,POSITION,unitary_d_up,&
             TOTAL_SITE_NUMBER**2,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
        call MPI_UNPACK(BUFFER,PACK_SIZE,POSITION,unitary_d_down,&
             TOTAL_SITE_NUMBER**2,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
        call MPI_UNPACK(BUFFER,PACK_SIZE,POSITION,distance_table,&
             TOTAL_SITE_NUMBER**2,MPI_INTEGER,MPI_COMM_WORLD,IERROR)
        call MPI_UNPACK(BUFFER,PACK_SIZE,POSITION,neighbor_table,&
             TOTAL_SITE_NUMBER**2,MPI_INTEGER,MPI_COMM_WORLD,IERROR)
        call MPI_UNPACK(BUFFER,PACK_SIZE,POSITION,neighbor_table2,&
             TOTAL_SITE_NUMBER*4,MPI_INTEGER,MPI_COMM_WORLD,IERROR)
        call MPI_UNPACK(BUFFER,PACK_SIZE,POSITION,distance_sequence,&
             TOTAL_SITE_NUMBER**2,MPI_INTEGER,MPI_COMM_WORLD,IERROR)
        call MPI_UNPACK(BUFFER,PACK_SIZE,POSITION,max_number_xi,1,MPI_INTEGER,&
             MPI_COMM_WORLD,IERROR)
        call MPI_UNPACK(BUFFER,PACK_SIZE,POSITION,electric_field,1,MPI_DOUBLE_PRECISION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_UNPACK(BUFFER,PACK_SIZE,POSITION,xi_field,1,MPI_DOUBLE_PRECISION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_UNPACK(BUFFER,PACK_SIZE,POSITION,coulomb,1,MPI_DOUBLE_PRECISION,&
             MPI_COMM_WORLD,IERROR)
        call MPI_UNPACK(BUFFER,PACK_SIZE,POSITION,xaxis,TOTAL_SITE_NUMBER,MPI_INTEGER,&
             MPI_COMM_WORLD,IERROR)
     end if
  end if

! max_number_xiを得てからじゃないとallocateできないもの
  call init_global_variables_after(name,error)

  if( MYRANK/=0 ) then
     allocate(seq_xi(max_number_xi),stat=ier)
     call stat_check("seq_xi","vmc",name,error)
     seq_xi=0    
  end if

  call MPI_BCAST(seq_xi,max_number_xi,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        xi_table(i_count,j_count)%way=distance_sequence(i_count,j_count)
        xi_table(i_count,j_count)%xi=seq_xi(distance_sequence(i_count,j_count))
     end do
  end do

! 相関関数用配列の確保
  call allocate_global_variables(name,error)
  call error_check(name,error)

!
! 相関関数power用
!

! sz
  allocate(sz_nume1(max_number_xi),stat=ier)
  call stat_check("sz_nume1","vmc",name,error) ; sz_nume1=0
  allocate(sz_nume2(max_number_xi),stat=ier)
  call stat_check("sz_nume2","vmc",name,error) ; sz_nume2=0
  allocate(sz_nume3_1(max_number_xi),stat=ier)
  call stat_check("sz_nume3_1","vmc",name,error) ; sz_nume3_1=0
  allocate(sz_nume3_2(max_number_xi),stat=ier)
  call stat_check("sz_nume3_2","vmc",name,error) ; sz_nume3_2=0
  allocate(sz_nume4(max_number_xi),stat=ier)
  call stat_check("sz_nume4","vmc",name,error) ; sz_nume4=0
  allocate(sz_nume5(max_number_xi),stat=ier)
  call stat_check("sz_nume5","vmc",name,error) ; sz_nume5=0
! charge
  allocate(charge_nume1(max_number_xi),stat=ier)
  call stat_check("charge_nume1","vmc",name,error) ; charge_nume1=0
  allocate(charge_nume2(max_number_xi),stat=ier)
  call stat_check("charge_nume2","vmc",name,error) ; charge_nume2=0
  allocate(charge_nume3_1(max_number_xi),stat=ier)
  call stat_check("charge_nume3_1","vmc",name,error) ; charge_nume3_1=0
  allocate(charge_nume3_2(max_number_xi),stat=ier)
  call stat_check("charge_nume3_2","vmc",name,error) ; charge_nume3_2=0
  allocate(charge_nume4(max_number_xi),stat=ier)
  call stat_check("charge_nume4","vmc",name,error) ; charge_nume4=0
  allocate(charge_nume5(max_number_xi),stat=ier)
  call stat_check("charge_nume5","vmc",name,error) ; charge_nume5=0
! charge_density
  allocate(charge_density_nume1(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("charge_density_nume1","vmc",name,error)
  charge_density_nume1=0
  allocate(charge_density_nume2(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("charge_density_nume2","vmc",name,error)
  charge_density_nume2=0
  allocate(charge_density_nume3_1(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("charge_density_nume3_1","vmc",name,error)
  charge_density_nume3_1=0
  allocate(charge_density_nume3_2(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("charge_density_nume3_2","vmc",name,error)
  charge_density_nume3_2=0
  allocate(charge_density_nume4(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("charge_density_nume4","vmc",name,error)
  charge_density_nume4=0
  allocate(charge_density_nume5(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("charge_density_nume5","vmc",name,error)
  charge_density_nume5=0
! sxsy
  allocate(sxsy_nume1(max_number_xi),stat=ier)
  call stat_check("sxsy_nume1","vmc",name,error) ; sxsy_nume1=0
  allocate(sxsy_nume2(max_number_xi),stat=ier)
  call stat_check("sxsy_nume2","vmc",name,error) ; sxsy_nume2=0
  allocate(sxsy_nume3_1(max_number_xi),stat=ier)
  call stat_check("sxsy_nume3_1","vmc",name,error) ; sxsy_nume3_1=0
  allocate(sxsy_nume3_2(max_number_xi),stat=ier)
  call stat_check("sxsy_nume3_2","vmc",name,error) ; sxsy_nume3_2=0
  allocate(sxsy_nume4(max_number_xi),stat=ier)
  call stat_check("sxsy_nume4","vmc",name,error) ; sxsy_nume4=0
  allocate(sxsy_nume5(max_number_xi),stat=ier)
  call stat_check("sxsy_nume5","vmc",name,error) ; sxsy_nume5=0
! super
  allocate(super_nume1(max_number_xi),stat=ier)
  call stat_check("super_nume1","vmc",name,error) ; super_nume1=0
  allocate(super_nume2(max_number_xi),stat=ier)
  call stat_check("super_nume2","vmc",name,error) ; super_nume2=0
  allocate(super_nume3_1(max_number_xi),stat=ier)
  call stat_check("super_nume3_1","vmc",name,error) ; super_nume3_1=0
  allocate(super_nume3_2(max_number_xi),stat=ier)
  call stat_check("super_nume3_2","vmc",name,error) ; super_nume3_2=0
  allocate(super_nume4(max_number_xi),stat=ier)
  call stat_check("super_nume4","vmc",name,error) ; super_nume4=0
  allocate(super_nume5(max_number_xi),stat=ier)
  call stat_check("super_nume5","vmc",name,error) ; super_nume5=0
!
! テンポラリ
!
  allocate(tmp_sz0(max_number_xi),stat=ier)
  call stat_check("tmp_sz0","vmc",name,error) ; tmp_sz0=0
  allocate(tmp_charge0(max_number_xi),stat=ier)
  call stat_check("tmp_charge0","vmc",name,error) ; tmp_charge0=0
  allocate(tmp_sxsy0(max_number_xi),stat=ier)
  call stat_check("tmp_sxsy0","vmc",name,error) ; tmp_sxsy0=0
  allocate(tmp_super0(max_number_xi),stat=ier)
  call stat_check("tmp_super0","vmc",name,error) ; tmp_super0=0
! MPI_REDUCEした際のバッファ
  allocate(buff_cor(max_number_xi),stat=ier)
  call stat_check("buff_cor","vmc",name,error) ; buff_cor=0
  allocate(buff_cd(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("buff_cd","vmc",name,error) ; buff_cd=0
!
! エラーバー用
!
  allocate(lowest_sz0(max_number_xi),stat=ier)
  call stat_check("lowest_sz0","vmc",name,error)
  lowest_sz0=0
  allocate(lowest_sz1(max_number_xi),stat=ier)
  call stat_check("lowest_sz1","vmc",name,error)
  lowest_sz1=0
  allocate(lowest_sz2(max_number_xi),stat=ier)
  call stat_check("lowest_sz2","vmc",name,error)
  lowest_sz2=0

  allocate(lowest_charge0(max_number_xi),stat=ier)
  call stat_check("lowest_charge0","vmc",name,error)
  lowest_charge0=0
  allocate(lowest_charge1(max_number_xi),stat=ier)
  call stat_check("lowest_charge1","vmc",name,error)
  lowest_charge1=0
  allocate(lowest_charge2(max_number_xi),stat=ier)
  call stat_check("lowest_charge2","vmc",name,error)
  lowest_charge2=0

  allocate(lowest_sxsy0(max_number_xi),stat=ier)
  call stat_check("lowest_sxsy0","vmc",name,error)
  lowest_sxsy0=0
  allocate(lowest_sxsy1(max_number_xi),stat=ier)
  call stat_check("lowest_sxsy1","vmc",name,error)
  lowest_sxsy1=0
  allocate(lowest_sxsy2(max_number_xi),stat=ier)
  call stat_check("lowest_sxsy2","vmc",name,error)
  lowest_sxsy2=0

  allocate(lowest_s0(max_number_xi),stat=ier)
  call stat_check("lowest_s0","vmc",name,error)
  lowest_s0=0
  allocate(lowest_s1(max_number_xi),stat=ier)
  call stat_check("lowest_s1","vmc",name,error)
  lowest_s1=0
  allocate(lowest_s2(max_number_xi),stat=ier)
  call stat_check("lowest_s2","vmc",name,error)
  lowest_s2=0

  allocate(lowest_super0(max_number_xi),stat=ier)
  call stat_check("lowest_super0","vmc",name,error)
  lowest_super0=0
  allocate(lowest_super1(max_number_xi),stat=ier)
  call stat_check("lowest_super1","vmc",name,error)
  lowest_super1=0
  allocate(lowest_super2(max_number_xi),stat=ier)
  call stat_check("lowest_super2","vmc",name,error)
  lowest_super2=0

  allocate(lowest_charge_density0(TOTAL_SITE_NUMBER),stat=ier) 
  call stat_check("lowest_charge_density0","vmc",name,error)
  lowest_charge_density0=0
  allocate(lowest_charge_density1(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("lowest_charge_density1","vmc",name,error)
  lowest_charge_density1=0
  allocate(lowest_charge_density2(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("lowest_charge_density2","vmc",name,error)
  lowest_charge_density2=0
!
  allocate(highest_sz0(max_number_xi),stat=ier)
  call stat_check("highest_sz0","vmc",name,error)
  highest_sz0=0
  allocate(highest_sz1(max_number_xi),stat=ier)
  call stat_check("highest_sz1","vmc",name,error)
  highest_sz1=0
  allocate(highest_sz2(max_number_xi),stat=ier)
  call stat_check("highest_sz2","vmc",name,error)
  highest_sz2=0

  allocate(highest_charge0(max_number_xi),stat=ier)
  call stat_check("highest_charge0","vmc",name,error)
  highest_charge0=0
  allocate(highest_charge1(max_number_xi),stat=ier)
  call stat_check("highest_charge1","vmc",name,error)
  highest_charge1=0
  allocate(highest_charge2(max_number_xi),stat=ier)
  call stat_check("highest_charge2","vmc",name,error)
  highest_charge2=0

  allocate(highest_sxsy0(max_number_xi),stat=ier)
  call stat_check("highest_sxsy0","vmc",name,error)
  highest_sxsy0=0
  allocate(highest_sxsy1(max_number_xi),stat=ier)
  call stat_check("highest_sxsy1","vmc",name,error)
  highest_sxsy1=0
  allocate(highest_sxsy2(max_number_xi),stat=ier)
  call stat_check("highest_sxsy2","vmc",name,error)
  highest_sxsy2=0

  allocate(highest_s0(max_number_xi),stat=ier)
  call stat_check("highest_s0","vmc",name,error)
  highest_s0=0
  allocate(highest_s1(max_number_xi),stat=ier)
  call stat_check("highest_s1","vmc",name,error)
  highest_s1=0
  allocate(highest_s2(max_number_xi),stat=ier)
  call stat_check("highest_s2","vmc",name,error)
  highest_s2=0

  allocate(highest_super0(max_number_xi),stat=ier)
  call stat_check("highest_super0","vmc",name,error)
  highest_super0=0
  allocate(highest_super1(max_number_xi),stat=ier)
  call stat_check("highest_super1","vmc",name,error)
  highest_super1=0
  allocate(highest_super2(max_number_xi),stat=ier)
  call stat_check("highest_super2","vmc",name,error)
  highest_super2=0

  allocate(highest_charge_density0(TOTAL_SITE_NUMBER),stat=ier) 
  call stat_check("highest_charge_density0","vmc",name,error)
  highest_charge_density0=0
  allocate(highest_charge_density1(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("highest_charge_density1","vmc",name,error)
  highest_charge_density1=0
  allocate(highest_charge_density2(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("highest_charge_density2","vmc",name,error)
  highest_charge_density2=0

!
! エネルギー期待値の計算
!
  sample_number=1

  do while( sample_number<=MAX_MONTECARLO_SAMPLE )
     ! モンテカルロ・サンプルの収集
     call correct_montecarlo_sample(sample_number,vector_up,vector_down,&
          alpha_psi0,alpha_psi1,lambda,name,error)
     call error_check(name,error)

     if( mod(sample_number+1,NUMBER_PE)==MYRANK ) then
        ! (1).エネルギー期待値計算

        call energy(vector_up,vector_down,alpha_psi0,lambda,&
             alpha_h_psi1,name,error)
        call error_check(name,error)

        result_pow0=result_pow0+alpha_h_psi1/alpha_psi1 &
             /dble(MAX_MONTECARLO_SAMPLE)

        if( PERTUR==1 ) then
           call calc_alpha_p_psi1(vector_up,vector_down,alpha_psi0,lambda,&
                alpha_p_psi1,name,error)
           call calc_alpha_p_p_psi1(vector_up,vector_down,alpha_psi0,&
                alpha_p_p_psi1,name,error)
           call calc_alpha_h_p_psi1(vector_up,vector_down,alpha_psi0,lambda,&
                alpha_h_p_psi1,name,error)
           result_pow0=result_pow0+2*electric_field**2*alpha_p_psi1**2 &
                +alpha_p_p_psi1*alpha_h_psi1+alpha_p_psi1*alpha_h_p_psi1
           bunsi=bunsi+alpha_p_psi1**2
           bunbo=bunbo+alpha_p_psi1*alpha_h_p_psi1 &
                +alpha_p_psi1*alpha_h_p_psi1
        end if

        call input_splited_energy0(alpha_psi1,alpha_h_psi1,sample_number,&
             name,error)
        call error_check(name,error) 

        ! (2).エネルギー期待値のpower(1次,2次)
        if( POWER==1 .or. POWER==2 )then
           
           tmp_h0=0 ; tmp_h1=0 ; tmp_h2=0
           tmp_h3=0 ; tmp_h4=0 ; tmp_h5=0

           call energy_power(vector_up,vector_down,alpha_psi1,&
                alpha_h_psi1,tmp_h0,tmp_h1,tmp_h2,tmp_h3,tmp_h4,tmp_h5,&
                alpha_h_h_psi1,alpha_h_h_h_psi1,name,error)
           call error_check(name,error)

           h0=h0+tmp_h0 ; h1=h1+tmp_h1 ; h2=h2+tmp_h2
           h3=h3+tmp_h3 ; h4=h4+tmp_h4 ; h5=h5+tmp_h5

           call input_splited_energy_power(sample_number,tmp_h0,tmp_h1,&
                tmp_h2,tmp_h3,tmp_h4,tmp_h5,name,error)
           call error_check(name,error)
        end if

        ! (3).相関関数計算1(ナンバーオペレーター)…対角性分のみ
        if( COR==1 .or. COR==2 ) then
           call cor_number(tmp_sz0,tmp_charge0,tmp_charge_density0,name,error)
           call error_check(name,error)

           call input_cor_number(tmp_sz0,tmp_charge0,tmp_charge_density0,name,error)
           call error_check(name,error)

           call input_splited_cor_number(sample_number,tmp_sz0,tmp_charge0,&
                tmp_charge_density0,name,error)
           call error_check(name,error)

        end if

        ! (5).相関関数計算2(ナンバーオペレーター以外)

        if( COR==2 ) then
           call cor_general(alpha_psi0,lambda,vector_up,vector_down,&
                tmp_sxsy0,tmp_super0,name,error)
           call error_check(name,error)

           call input_cor_general(alpha_psi1,&
                tmp_sxsy0,tmp_super0,tmp_sz0,name,error)
           call error_check(name,error)

           call input_splited_cor_general(sample_number,alpha_psi1,&
                tmp_sxsy0,tmp_super0,tmp_sz0,name,error)
           call error_check(name,error)
        end if
     end if

!
! もしPOWER==0ならvectorを出力(POW1,2の計算をするなら後で出力)
!
! exp_code,サンプル番号,<alpha|psi0>,<alpha|psi1>,<alpha|h|psi1>
! <alpha|hh|psi1>,<alpha|hhh|psi1>,lambdaを出力

     if( mod(sample_number+1,NUMBER_PE)==MYRANK ) then
        state_vector(sample_number)%vector_up=vector_up
        state_vector(sample_number)%vector_down=vector_down
        state_vector(sample_number)%alpha_psi0=alpha_psi0
        state_vector(sample_number)%alpha_psi1=alpha_psi1
        state_vector(sample_number)%alpha_h_psi1=alpha_h_psi1
        state_vector(sample_number)%alpha_h_h_psi1=alpha_h_h_psi1
        state_vector(sample_number)%alpha_h_h_h_psi1=alpha_h_h_h_psi1
        state_vector(sample_number)%lambda=lambda
     end if

     sample_number=sample_number+1
  end do

!
! REDUCE(POW1,POW2のための要素)
!
  buff=0
  call MPI_REDUCE(h0,buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  h0=buff ; buff=0
  call MPI_REDUCE(h1,buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  h1=buff ; buff=0
  call MPI_REDUCE(h2,buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  h2=buff ; buff=0
  call MPI_REDUCE(h3,buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  h3=buff ; buff=0
  call MPI_REDUCE(h4,buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  h4=buff ; buff=0
  call MPI_REDUCE(h5,buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  h5=buff ; buff=0
  
!
! REDUCE(POW0の相関関数)
!
  buff_cor=0
  call MPI_REDUCE(result_sz0,buff_cor,max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
  result_sz0=buff_cor

  buff_cor=0
  call MPI_REDUCE(result_sxsy0,buff_cor,max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
  result_sxsy0=buff_cor

  buff_cor=0
  call MPI_REDUCE(result_s0,buff_cor,max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
  result_s0=buff_cor

  buff_cor=0
  call MPI_REDUCE(result_charge0,buff_cor,max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
  result_charge0=buff_cor

  buff_cor=0
  call MPI_REDUCE(result_super0,buff_cor,max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
  result_super0=buff_cor

  buff_cd=0
  call MPI_REDUCE(result_charge_density0,buff_cd,&
       TOTAL_SITE_NUMBER,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
  result_charge_density0=buff_cd

  buff=0
  call MPI_REDUCE(result_pole0,buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  result_pole0=buff

!
! POWERが最小値になるように係数を決定する。
!
  if( MYRANK==0 ) then
     if( POWER==1 .or. POWER==2 ) then
! (pow1)エネルギーが最小値になるように系数C0,C1を選ぶ

        call min1(h0,h1,h2,h3,pow1_c0,pow1_c1,result_pow1,name,error)
        call error_check(name,error)

! (pow2)エネルギーが最小値になるように系数C0,C1,C2を選ぶ

        if( POWER==2 ) then
           call min2(h0,h1,h2,h3,h4,h5,pow2_c0,pow2_c1,pow2_c2,&
                result_pow2,name,error)
           call error_check(name,error)

           write(*,*) "h1/h0=",h1/h0
           write(*,*) "h2/h1=",h2/h1
           write(*,*) "h3/h2=",h3/h2
           write(*,*) "h4/h3=",h4/h3
           write(*,*) "h5/h4=",h5/h4
           write(*,*) 
           write(*,*) "h2/h0=",h2/h0
           write(*,*) "h3/h1=",h3/h1
           write(*,*) "h4/h2=",h4/h2
           write(*,*) "h5/h3=",h5/h3
           

        end if
     end if
  end if

!
! coefficientを全てのPCに配る
!
  call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

  call MPI_BCAST(pow1_c0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
  call error_check("MPI_BCAST",IERROR)
  call MPI_BCAST(pow1_c1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
  call error_check("MPI_BCAST",IERROR)
  call MPI_BCAST(pow2_c0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
  call error_check("MPI_BCAST",IERROR)
  call MPI_BCAST(pow2_c1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
  call error_check("MPI_BCAST",IERROR)
  call MPI_BCAST(pow2_c2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
  call error_check("MPI_BCAST",IERROR)

!
! REDUCE(POW0のエネルギー)
!

  buff=0
  call MPI_REDUCE(result_pow0,buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  result_pow0=buff
!
! REDUCE
!
  buff_splited=0
  call MPI_REDUCE(splited_energy,buff_splited,5*MAX_MONTECARLO_SAMPLE/5,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_energy=buff_splited                           ! 移す

!
! POW0のエラーバーを出す
!
  if( MYRANK==0 ) then
     error_pow0=0
  
     do i_count=1,5
        do j_count=1,MAX_MONTECARLO_SAMPLE/5
           error_pow0(i_count)=error_pow0(i_count) &
                +splited_energy(i_count,j_count)
        end do
     end do

     ! 仮に代入
     lowest_pow0=error_pow0(1)
     highest_pow0=error_pow0(1)

     do i_count=2,5
        if( lowest_pow0>error_pow0(i_count) ) then
           lowest_pow0=error_pow0(i_count)
        end if
        if( highest_pow0<error_pow0(i_count) ) then
           highest_pow0=error_pow0(i_count)
        end if
     end do

  end if

!
! REDUCE(splited_h0,h1,h2,h3,h4,h5)
!
  buff5=0
  call MPI_REDUCE(splited_h0,buff5,5,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_h0=buff5

  buff5=0
  call MPI_REDUCE(splited_h1,buff5,5,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_h1=buff5

  buff5=0
  call MPI_REDUCE(splited_h2,buff5,5,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_h2=buff5

  buff5=0
  call MPI_REDUCE(splited_h3,buff5,5,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_h3=buff5

  buff5=0
  call MPI_REDUCE(splited_h4,buff5,5,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_h4=buff5

  buff5=0
  call MPI_REDUCE(splited_h5,buff5,5,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_h5=buff5

!
! POW1のエラーバーを出す。
!
  if( MYRANK==0 ) then
     if( POWER==1 .or. POWER==2 ) then
        ! 1
        tmp_h0=splited_h0(1) ; tmp_h1=splited_h1(1)
        tmp_h2=splited_h2(1) ; tmp_h3=splited_h3(1)

        tmp=0

        call min1(tmp_h0,tmp_h1,tmp_h2,tmp_h3,&
             dummy0,dummy1,tmp,name,error)
        call error_check(name,error)

        error_pow1(1)=tmp
        ! 2
        tmp_h0=splited_h0(2) ; tmp_h1=splited_h1(2)
        tmp_h2=splited_h2(2) ; tmp_h3=splited_h3(2)

        tmp=0
        call min1(tmp_h0,tmp_h1,tmp_h2,tmp_h3,&
             dummy0,dummy1,tmp,name,error)
        call error_check(name,error)
        error_pow1(2)=tmp
        ! 3
        tmp_h0=splited_h0(3) ; tmp_h1=splited_h1(3)
        tmp_h2=splited_h2(3) ; tmp_h3=splited_h3(3)

        tmp=0
        call min1(tmp_h0,tmp_h1,tmp_h2,tmp_h3,&
             dummy0,dummy1,tmp,name,error)
        call error_check(name,error)
        error_pow1(3)=tmp
        ! 4
        tmp_h0=splited_h0(4) ; tmp_h1=splited_h1(4)
        tmp_h2=splited_h2(4) ; tmp_h3=splited_h3(4)

        tmp=0
        call min1(tmp_h0,tmp_h1,tmp_h2,tmp_h3,&
             dummy0,dummy1,tmp,name,error)
        call error_check(name,error)
        error_pow1(4)=tmp
        ! 5
        tmp_h0=splited_h0(5) ; tmp_h1=splited_h1(5)
        tmp_h2=splited_h2(5) ; tmp_h3=splited_h3(5)

        tmp=0
        call min1(tmp_h0,tmp_h1,tmp_h2,tmp_h3,&
             dummy0,dummy1,tmp,name,error)
        call error_check(name,error)
        error_pow1(5)=tmp

        ! 仮に代入
        lowest_pow1=error_pow1(1)
        highest_pow1=error_pow1(1)

        do i_count=2,5
           if( lowest_pow1>error_pow1(i_count) ) then
              lowest_pow1=error_pow1(i_count)
           end if
           if( highest_pow1<error_pow1(i_count) ) then
              highest_pow1=error_pow1(i_count)
           end if
        end do
     end if
     

     if( POWER==2 ) then
        ! 1
        tmp_h0=splited_h0(1) ; tmp_h1=splited_h1(1) ; tmp_h2=splited_h2(1) 
        tmp_h3=splited_h3(1) ; tmp_h4=splited_h4(1) ; tmp_h5=splited_h5(1)

        tmp=0

        call min2(tmp_h0,tmp_h1,tmp_h2,tmp_h3,tmp_h4,tmp_h5,&
             dummy0,dummy1,dummy2,tmp,name,error)
        call error_check(name,error)

        error_pow2(1)=tmp

        ! 2
        tmp_h0=splited_h0(2) ; tmp_h1=splited_h1(2) ; tmp_h2=splited_h2(2) 
        tmp_h3=splited_h3(2) ; tmp_h4=splited_h4(2) ; tmp_h5=splited_h5(2)
     
        tmp=0

        call min2(tmp_h0,tmp_h1,tmp_h2,tmp_h3,tmp_h4,tmp_h5,&
             dummy0,dummy1,dummy2,tmp,name,error)
        call error_check(name,error)

        error_pow2(2)=tmp

        ! 3
        tmp_h0=splited_h0(3) ; tmp_h1=splited_h1(3) ; tmp_h2=splited_h2(3) 
        tmp_h3=splited_h3(3) ; tmp_h4=splited_h4(3) ; tmp_h5=splited_h5(3)

        tmp=0

        call min2(tmp_h0,tmp_h1,tmp_h2,tmp_h3,tmp_h4,tmp_h5,&
             dummy0,dummy1,dummy2,tmp,name,error)
        call error_check(name,error)

        error_pow2(3)=tmp

        ! 4
        tmp_h0=splited_h0(4) ; tmp_h1=splited_h1(4) ; tmp_h2=splited_h2(4) 
        tmp_h3=splited_h3(4) ; tmp_h4=splited_h4(4) ; tmp_h5=splited_h5(4)
        
        tmp=0

        call min2(tmp_h0,tmp_h1,tmp_h2,tmp_h3,tmp_h4,tmp_h5,&
             dummy0,dummy1,dummy2,tmp,name,error)
        call error_check(name,error)

        error_pow2(4)=tmp

        ! 5
        tmp_h0=splited_h0(5) ; tmp_h1=splited_h1(5) ; tmp_h2=splited_h2(5) 
        tmp_h3=splited_h3(5) ; tmp_h4=splited_h4(5) ; tmp_h5=splited_h5(5)

        tmp=0

        call min2(tmp_h0,tmp_h1,tmp_h2,tmp_h3,tmp_h4,tmp_h5,&
             dummy0,dummy1,dummy2,tmp,name,error)
        call error_check(name,error)

        error_pow2(5)=tmp

        ! 仮に代入
        lowest_pow2=error_pow2(1)
        highest_pow2=error_pow2(1)

        do i_count=2,5
           if( lowest_pow2>error_pow2(i_count) ) then
              lowest_pow2=error_pow2(i_count)
           end if
           if( highest_pow2<error_pow2(i_count) ) then
              highest_pow2=error_pow2(i_count)
           end if
        end do
     end if
  end if

#ifdef DEBUG
  write(*,*) "delta_e=",delta_e/dble(MAX_MONTECARLO_SAMPLE)
#endif

!
! 一応エネルギー期待値の計算は終り。次にCOR==1or2 .and. POWER/=0の場合
!
  if( (POWER/=0 .and. COR==1) .or. (POWER/=0 .and. COR==2)  ) then
!!
!! 相関関数のPowerを求めるルーチン(もう一回同じ計算を繰り返す必要がある。)
!!

! 乱数の初期化(さっきと同じ結果を得るために必須、同じ種icで)
     call random_init(ic,name,error)
     call error_check(name,error)

     sample_number=1

     do while( sample_number<=MAX_MONTECARLO_SAMPLE )
        ! モンテカルロ・サンプルの収集
        call correct_montecarlo_sample(sample_number,vector_up,vector_down,&
             alpha_psi0,alpha_psi1,lambda,name,error)
        call error_check(name,error)

        if( mod(sample_number+1,NUMBER_PE)==MYRANK ) then
           if( mod(sample_number,10000)==0 ) then
!              write(*,*) "COR:sample_number=",sample_number
           end if
           ! (1).エネルギー期待値計算
           ! エネルギー期待値の計算はすでに完了している
           ! (必要なのは、alpha_h_psi1だけ)

!!!!          call energy(vector_up,vector_down,alpha_psi0,lambda,&
!!!!              alpha_h_psi1,name,error)
!!!!          call error_check(name,error)

!!!!  inserted below
           alpha_h_psi1=state_vector(sample_number)%alpha_h_psi1

           ! (2).エネルギー期待値のpower(1次,2次)
           tmp_h0=0 ; tmp_h1=0 ; tmp_h2=0
           tmp_h3=0 ; tmp_h4=0 ; tmp_h5=0

!!!!           call energy_power(vector_up,vector_down,alpha_psi1,&
!!!!                alpha_h_psi1,tmp_h0,tmp_h1,tmp_h2,tmp_h3,tmp_h4,&
!!!!                tmp_h5,alpha_h_h_psi1,&
!!!!                alpha_h_h_h_psi1,name,error)
!!!!           call error_check(name,error)
           alpha_psi1=state_vector(sample_number)%alpha_psi1
           alpha_h_psi1=state_vector(sample_number)%alpha_h_psi1
           alpha_h_h_psi1=state_vector(sample_number)%alpha_h_h_psi1
           alpha_h_h_h_psi1=state_vector(sample_number)%alpha_h_h_h_psi1
           tmp_h0=(alpha_psi1**2)/(alpha_psi1**2)
           tmp_h1=(alpha_psi1*alpha_h_psi1)/(alpha_psi1**2)
           tmp_h2=(alpha_h_psi1**2)/(alpha_psi1**2)
           tmp_h3=(alpha_h_psi1*alpha_h_h_psi1)/(alpha_psi1**2)
           if( POWER==2 ) then 
              tmp_h4=(alpha_h_h_psi1**2)/(alpha_psi1**2)
              tmp_h5=(alpha_h_h_psi1*alpha_h_h_h_psi1)/(alpha_psi1**2)
           end if
           alpha_h_h_psi1=state_vector(sample_number)%alpha_h_h_psi1
           alpha_h_h_h_psi1=state_vector(sample_number)%alpha_h_h_h_psi1

           h0=h0+tmp_h0 ; h1=h1+tmp_h1 ; h2=h2+tmp_h2
           h3=h3+tmp_h3 ; h4=h4+tmp_h4 ; h5=h5+tmp_h5

           denomi_pow1=denomi_pow1+(pow1_c1**2)*(alpha_h_psi1**2) &
                /(alpha_psi1**2) &
                +2*pow1_c0*pow1_c1*(alpha_psi1*alpha_h_psi1)/(alpha_psi1**2) &
                +(pow1_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)
           denomi_pow2=denomi_pow2+(pow2_c2**2)*(alpha_h_h_psi1**2) &
                /(alpha_psi1**2) &
                +2*pow2_c1*pow2_c2*(alpha_h_psi1*alpha_h_h_psi1) &
                /(alpha_psi1**2) &
                +2*pow2_c0*pow2_c2*(alpha_h_h_psi1*alpha_psi1)/(alpha_psi1**2)&
                +(pow2_c1**2)*(alpha_h_psi1**2)/(alpha_psi1**2) &
                +(2*pow2_c0*pow2_c1)*(alpha_psi1*alpha_h_psi1)/(alpha_psi1**2)&
                +(pow2_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)

           if( mod(sample_number,5)==1 ) then
              denomi_splited_pow1(1)=denomi_splited_pow1(1) &
                   +(pow1_c1**2)*(alpha_h_psi1**2)/(alpha_psi1**2) &
                   +2*pow1_c0*pow1_c1*(alpha_psi1*alpha_h_psi1) &
                   /(alpha_psi1**2) &
                   +(pow1_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)
              denomi_splited_pow2(1)=denomi_splited_pow2(1) &
                   +(pow2_c2**2)*(alpha_h_h_psi1**2)/(alpha_psi1**2) &
                   +2*pow2_c1*pow2_c2*(alpha_h_psi1*alpha_h_h_psi1) &
                   /(alpha_psi1**2) &
                   +2*pow2_c0*pow2_c2*(alpha_h_h_psi1*alpha_psi1) &
                   /(alpha_psi1**2)&
                   +(pow2_c1**2)*(alpha_h_psi1**2)/(alpha_psi1**2) &
                   +(2*pow2_c0*pow2_c1)*(alpha_psi1*alpha_h_psi1) &
                   /(alpha_psi1**2)&
                   +(pow2_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)
           else if( mod(sample_number,5)==2 ) then
              denomi_splited_pow1(2)=denomi_splited_pow1(2) &
                   +(pow1_c1**2)*(alpha_h_psi1**2)/(alpha_psi1**2) &
                   +2*pow1_c0*pow1_c1*(alpha_psi1*alpha_h_psi1) &
                   /(alpha_psi1**2) &
                   +(pow1_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)
              denomi_splited_pow2(2)=denomi_splited_pow2(2) &
                   +(pow2_c2**2)*(alpha_h_h_psi1**2)/(alpha_psi1**2) &
                   +2*pow2_c1*pow2_c2*(alpha_h_psi1*alpha_h_h_psi1) &
                   /(alpha_psi1**2) &
                   +2*pow2_c0*pow2_c2*(alpha_h_h_psi1*alpha_psi1) &
                   /(alpha_psi1**2)&
                   +(pow2_c1**2)*(alpha_h_psi1**2)/(alpha_psi1**2) &
                   +(2*pow2_c0*pow2_c1)*(alpha_psi1*alpha_h_psi1) &
                   /(alpha_psi1**2)&
                   +(pow2_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)
           else if( mod(sample_number,5)==3 ) then
              denomi_splited_pow1(3)=denomi_splited_pow1(3) &
                   +(pow1_c1**2)*(alpha_h_psi1**2)/(alpha_psi1**2) &
                   +2*pow1_c0*pow1_c1*(alpha_psi1*alpha_h_psi1) &
                   /(alpha_psi1**2) &
                   +(pow1_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)
              denomi_splited_pow2(3)=denomi_splited_pow2(3) &
                   +(pow2_c2**2)*(alpha_h_h_psi1**2)/(alpha_psi1**2) &
                   +2*pow2_c1*pow2_c2*(alpha_h_psi1*alpha_h_h_psi1) &
                   /(alpha_psi1**2) &
                   +2*pow2_c0*pow2_c2*(alpha_h_h_psi1*alpha_psi1) &
                   /(alpha_psi1**2)&
                   +(pow2_c1**2)*(alpha_h_psi1**2)/(alpha_psi1**2) &
                   +(2*pow2_c0*pow2_c1)*(alpha_psi1*alpha_h_psi1) &
                   /(alpha_psi1**2)&
                   +(pow2_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)
           else if( mod(sample_number,5)==4 ) then
              denomi_splited_pow1(4)=denomi_splited_pow1(4) &
                   +(pow1_c1**2)*(alpha_h_psi1**2)/(alpha_psi1**2) &
                   +2*pow1_c0*pow1_c1*(alpha_psi1*alpha_h_psi1) &
                   /(alpha_psi1**2) &
                   +(pow1_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)
              denomi_splited_pow2(4)=denomi_splited_pow2(4) &
                   +(pow2_c2**2)*(alpha_h_h_psi1**2)/(alpha_psi1**2) &
                   +2*pow2_c1*pow2_c2*(alpha_h_psi1*alpha_h_h_psi1) &
                   /(alpha_psi1**2) &
                   +2*pow2_c0*pow2_c2*(alpha_h_h_psi1*alpha_psi1) &
                   /(alpha_psi1**2)&
                   +(pow2_c1**2)*(alpha_h_psi1**2)/(alpha_psi1**2) &
                   +(2*pow2_c0*pow2_c1)*(alpha_psi1*alpha_h_psi1) &
                   /(alpha_psi1**2)&
                   +(pow2_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)
           else if( mod(sample_number,5)==0 ) then
              denomi_splited_pow1(5)=denomi_splited_pow1(5) &
                   +(pow1_c1**2)*(alpha_h_psi1**2)/(alpha_psi1**2) &
                   +2*pow1_c0*pow1_c1*(alpha_psi1*alpha_h_psi1) &
                   /(alpha_psi1**2) &
                   +(pow1_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)
              denomi_splited_pow2(5)=denomi_splited_pow2(5) &
                   +(pow2_c2**2)*(alpha_h_h_psi1**2)/(alpha_psi1**2) &
                   +2*pow2_c1*pow2_c2*(alpha_h_psi1*alpha_h_h_psi1) &
                   /(alpha_psi1**2) &
                   +2*pow2_c0*pow2_c2*(alpha_h_h_psi1*alpha_psi1) &
                   /(alpha_psi1**2)&
                   +(pow2_c1**2)*(alpha_h_psi1**2)/(alpha_psi1**2) &
                   +(2*pow2_c0*pow2_c1)*(alpha_psi1*alpha_h_psi1) &
                   /(alpha_psi1**2)&
                   +(pow2_c0**2)*(alpha_psi1**2)/(alpha_psi1**2)
           end if


           ! (3).相関関数計算1(ナンバーオペレーター)…対角性分のみ
           ! すでに先の計算で相関関数のpow0は得ているのでストアしない
           ! (必要なのはtmpの結果のみ)

           call cor_number(tmp_sz0,tmp_charge0,tmp_charge_density0,name,error)
           call error_check(name,error)

           ! エラーバー用のルーチンはすでに実行しているので、ここではしない

           ! (4).相関関数の計算1のpower

           call cor_number_power(alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,&
                tmp_sz0,tmp_charge0,tmp_charge_density0,sz_nume1,&
                sz_nume2,sz_nume3_1,sz_nume3_2,&
                sz_nume4,sz_nume5,&
                charge_nume1,charge_nume2,charge_nume3_1,&
                charge_nume3_2,charge_nume4,charge_nume5,&
                charge_density_nume1,charge_density_nume2,&
                charge_density_nume3_1,charge_density_nume3_2,&
                charge_density_nume4,charge_density_nume5,&
                name,error)
           call error_check(name,error)

           call input_cor_number_power1(pow1_c0,pow1_c1,alpha_psi1,alpha_h_psi1,&
                sz_nume1,sz_nume2,sz_nume3_1,&
                sz_nume3_2,sz_nume4,sz_nume5,&
                charge_nume1,charge_nume2,charge_nume3_1,&
                charge_nume3_2,charge_nume4,charge_nume5,&
                charge_density_nume1,charge_density_nume2,&
                charge_density_nume3_1,charge_density_nume3_2,&
                charge_density_nume4,charge_density_nume5,&
                name,error)
           call error_check(name,error)

           if( POWER==2 ) then
              call input_cor_number_power2(pow2_c0,pow2_c1,pow2_c2,alpha_psi1,&
                   alpha_h_psi1,alpha_h_h_psi1,&
                   sz_nume1,sz_nume2,sz_nume3_1,&
                   sz_nume3_2,sz_nume4,sz_nume5,&
                   charge_nume1,charge_nume2,charge_nume3_1,&
                   charge_nume3_2,charge_nume4,charge_nume5,&
                   charge_density_nume1,charge_density_nume2,&
                   charge_density_nume3_1,charge_density_nume3_2,&
                   charge_density_nume4,charge_density_nume5,&
                   name,error)
           end if

           call input_splited_cor_number_power(sample_number,&
                alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,&
                pow1_c0,pow1_c1,pow2_c0,pow2_c1,pow2_c2,&
                sz_nume1,sz_nume2,sz_nume3_1,sz_nume3_2,&
                sz_nume4,sz_nume5,&
                charge_nume1,charge_nume2,charge_nume3_1,&
                charge_nume3_2,charge_nume4,charge_nume5,&
                charge_density_nume1,charge_density_nume2,&
                charge_density_nume3_1,charge_density_nume3_2,&
                charge_density_nume4,charge_density_nume5,&
                name,error)
           call error_check(name,error)

           ! (5).相関関数計算2(ナンバーオペレーター以外)
           ! すでに先の計算で相関関数のpow0は得ているのでストアしない
           ! (必要なのはtmpの結果のみ)
           if( COR==2 ) then
              call cor_general(alpha_psi0,lambda,vector_up,vector_down,&
                   tmp_sxsy0,tmp_super0,name,error)
              call error_check(name,error)

              ! (6).相関関数計算2のpower
              call cor_general_power(alpha_psi0,alpha_psi1,alpha_h_psi1,&
                   alpha_h_h_psi1,lambda,vector_up,vector_down,&
                   sxsy_nume1,sxsy_nume2,sxsy_nume3_1,&
                   sxsy_nume3_2,sxsy_nume4,sxsy_nume5,&
                   super_nume1,super_nume2,super_nume3_1,&
                   super_nume3_2,super_nume4,super_nume5,&
                   name,error)
              call error_check(name,error)

              call input_cor_general_power1(pow1_c0,pow1_c1,alpha_psi0,&
                   alpha_psi1,alpha_h_psi1,&
                   alpha_h_h_psi1,lambda,vector_up,vector_down,&
                   sz_nume1,sz_nume2,sz_nume3_1,&
                   sz_nume3_2,sz_nume4,sz_nume5,&
                   sxsy_nume1,sxsy_nume2,sxsy_nume3_1,&
                   sxsy_nume3_2,sxsy_nume4,sxsy_nume5,&
                   super_nume1,super_nume2,super_nume3_1,&
                   super_nume3_2,super_nume4,super_nume5,&
                   name,error)

              if( POWER==2 ) then
                 call input_cor_general_power2(pow2_c0,pow2_c1,pow2_c2,&
                      alpha_psi0,alpha_psi1,alpha_h_psi1,&
                      alpha_h_h_psi1,lambda,vector_up,vector_down,&
                      sz_nume1,sz_nume2,sz_nume3_1,&
                      sz_nume3_2,sz_nume4,sz_nume5,&
                      sxsy_nume1,sxsy_nume2,sxsy_nume3_1,&
                      sxsy_nume3_2,sxsy_nume4,sxsy_nume5,&
                      super_nume1,super_nume2,super_nume3_1,&
                      super_nume3_2,super_nume4,super_nume5,&
                      name,error)
              end if

              call input_splited_cor_general_power(sample_number,pow1_c0,&
                   pow1_c1,pow2_c0,pow2_c1,pow2_c1,alpha_psi0,alpha_psi1,&
                   alpha_h_psi1,alpha_h_h_psi1,lambda,vector_up,vector_down,&
                   sz_nume1,sz_nume2,sz_nume3_1,sz_nume3_2,sz_nume4,sz_nume5,&
                   sxsy_nume1,sxsy_nume2,sxsy_nume3_1,sxsy_nume3_2,&
                   sxsy_nume4,sxsy_nume5,&
                   super_nume1,super_nume2,super_nume3_1,super_nume3_2,&
                   super_nume4,super_nume5,name,error)
           end if
        end if

! exp_code,サンプル番号,<alpha|psi0>,<alpha|psi1>,<alpha|h|psi1>
! <alpha|hh|psi1>,<alpha|hhh|psi1>,lambdaを出力

        sample_number=sample_number+1
     end do

!
! REDUCE(POW1の相関関数)
!
     buff_cor=0
     call MPI_REDUCE(result_sz1,buff_cor,max_number_xi,&
          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
     result_sz1=buff_cor ; buff_cor=0
     call MPI_REDUCE(result_sxsy1,buff_cor,max_number_xi,&
          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
     result_sxsy1=buff_cor ; buff_cor=0
     call MPI_REDUCE(result_s1,buff_cor,max_number_xi,&
          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
     result_s1=buff_cor ; buff_cor=0
     call MPI_REDUCE(result_charge1,buff_cor,max_number_xi,&
          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
     result_charge1=buff_cor ; buff_cor=0
     call MPI_REDUCE(result_super1,buff_cor,max_number_xi,&
          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
     result_super1=buff_cor ; buff_cd=0
     call MPI_REDUCE(result_charge_density1,buff_cd,&
          TOTAL_SITE_NUMBER,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,IERROR)
     result_charge_density1=buff_cd ; buff=0
     call MPI_REDUCE(result_pole1,buff,&
          1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,IERROR)
     result_pole1=buff     


!
! REDUCE(POW2の相関関数)
!
     buff_cor=0
     call MPI_REDUCE(result_sz2,buff_cor,max_number_xi,&
          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
     result_sz2=buff_cor ; buff_cor=0
     call MPI_REDUCE(result_sxsy2,buff_cor,max_number_xi,&
          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
     result_sxsy2=buff_cor ; buff_cor=0
     call MPI_REDUCE(result_s2,buff_cor,max_number_xi,&
          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
     result_s2=buff_cor ; buff_cor=0
     call MPI_REDUCE(result_charge2,buff_cor,max_number_xi,&
          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
     result_charge2=buff_cor ; buff_cor=0
     call MPI_REDUCE(result_super2,buff_cor,max_number_xi,&
          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR) 
     result_super2=buff_cor ; buff_cd=0
     call MPI_REDUCE(result_charge_density2,buff_cd,&
          TOTAL_SITE_NUMBER,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,IERROR)
     result_charge_density2=buff_cd ; buff=0
     call MPI_REDUCE(result_pole2,buff,&
          1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,IERROR)
     result_pole2=buff

!
! denomi_pow1,denomi_pow2
!
     call MPI_REDUCE(denomi_pow1,buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,IERROR)
     denomi_pow1=buff
     call MPI_REDUCE(denomi_pow2,buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,IERROR)
     denomi_pow2=buff

     if( MYRANK==0 ) then
        result_sz1=result_sz1/denomi_pow1/dble(TOTAL_SITE_NUMBER)
        result_charge1=result_charge1/denomi_pow1/dble(TOTAL_SITE_NUMBER)
        result_charge_density1=result_charge_density1/denomi_pow1
        result_sxsy1=result_sxsy1/denomi_pow1/dble(TOTAL_SITE_NUMBER)
        result_s1=result_s1/denomi_pow1/dble(TOTAL_SITE_NUMBER)
        result_super1=result_super1/denomi_pow1/dble(TOTAL_SITE_NUMBER)
        result_pole1=result_pole1/denomi_pow1
     end if

!
! denomi_splited_pow1,denomi_splited_pow2
!
     buff5=0
     call MPI_REDUCE(denomi_splited_pow1,buff5,5,MPI_DOUBLE_PRECISION,&
          MPI_SUM,0,MPI_COMM_WORLD,IERROR)
     denomi_splited_pow1=buff5

     buff5=0
     call MPI_REDUCE(denomi_splited_pow2,buff5,5,MPI_DOUBLE_PRECISION,&
          MPI_SUM,0,MPI_COMM_WORLD,IERROR)
     denomi_splited_pow2=buff5

     if( MYRANK==0 ) then
        if( POWER==2 ) then
           result_sz2=result_sz2/denomi_pow2/dble(TOTAL_SITE_NUMBER)
           result_charge2=result_charge2/denomi_pow2/dble(TOTAL_SITE_NUMBER)
           result_charge_density2=result_charge_density2/denomi_pow2
           result_sxsy2=result_sxsy2/denomi_pow2/dble(TOTAL_SITE_NUMBER)
           result_s2=result_s2/denomi_pow2/dble(TOTAL_SITE_NUMBER)
           result_super2=result_super2/denomi_pow2/dble(TOTAL_SITE_NUMBER)
           result_pole2=result_pole2/denomi_pow2
        end if
     end if

#ifdef DEBUG
     if( MYRANK==0 .and. (POWER==1 .or. POWER==2) ) then
        call display_h_pow(h0,h1,h2,h3,h4,h5,name,error)
     end if
#endif
  end if

! alpha_psi0をreduce
  tmp_real=0 ; total_real=0

  do i_count=1,MAX_MONTECARLO_SAMPLE
     tmp_real(i_count)=state_vector(i_count)%alpha_psi0
  end do
  call MPI_REDUCE(tmp_real,total_real,MAX_MONTECARLO_SAMPLE,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  do i_count=1,MAX_MONTECARLO_SAMPLE
     state_vector(i_count)%alpha_psi0=total_real(i_count)
  end do

! alpha_psi1をreduce
  tmp_real=0 ; total_real=0
  do i_count=1,MAX_MONTECARLO_SAMPLE
     tmp_real(i_count)=state_vector(i_count)%alpha_psi1
  end do
  call MPI_REDUCE(tmp_real,total_real,MAX_MONTECARLO_SAMPLE,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  do i_count=1,MAX_MONTECARLO_SAMPLE
     state_vector(i_count)%alpha_psi1=total_real(i_count)
  end do

! alpha_h_psi1をreduce
  tmp_real=0 ; total_real=0
  do i_count=1,MAX_MONTECARLO_SAMPLE
     tmp_real(i_count)=state_vector(i_count)%alpha_h_psi1
  end do
  call MPI_REDUCE(tmp_real,total_real,MAX_MONTECARLO_SAMPLE,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  do i_count=1,MAX_MONTECARLO_SAMPLE
     state_vector(i_count)%alpha_h_psi1=total_real(i_count)
  end do

! alpha_h_h_psi1をreduce
  tmp_real=0 ; total_real=0
  do i_count=1,MAX_MONTECARLO_SAMPLE
     tmp_real(i_count)=state_vector(i_count)%alpha_h_h_psi1
  end do
  call MPI_REDUCE(tmp_real,total_real,MAX_MONTECARLO_SAMPLE,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  do i_count=1,MAX_MONTECARLO_SAMPLE
     state_vector(i_count)%alpha_h_h_psi1=total_real(i_count)
  end do

! alpha_h_h_h_psi1をreduce
  tmp_real=0 ; total_real=0
  do i_count=1,MAX_MONTECARLO_SAMPLE
     tmp_real(i_count)=state_vector(i_count)%alpha_h_h_h_psi1
  end do
  call MPI_REDUCE(tmp_real,total_real,MAX_MONTECARLO_SAMPLE,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  do i_count=1,MAX_MONTECARLO_SAMPLE
     state_vector(i_count)%alpha_h_h_h_psi1=total_real(i_count)
  end do

! lambdaをreduce
  tmp_real=0 ; total_real=0
  do i_count=1,MAX_MONTECARLO_SAMPLE
     tmp_real(i_count)=state_vector(i_count)%lambda
  end do
  call MPI_REDUCE(tmp_real,total_real,MAX_MONTECARLO_SAMPLE,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  do i_count=1,MAX_MONTECARLO_SAMPLE
     state_vector(i_count)%lambda=total_real(i_count)
  end do

! up
/*
  do j_count=1,MAX_MONTECARLO_SAMPLE
     tmp_vector_up=0 ; total_vector_up=0

     tmp_vector_up=state_vector(j_count)%vector_up

     call MPI_REDUCE(tmp_vector_up,total_vector_up,TOTAL_UP_ELECTRON,&
          MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)

     state_vector(j_count)%vector_up=total_vector_up
  end do

! down
  do j_count=1,MAX_MONTECARLO_SAMPLE
     tmp_vector_down=0 ; total_vector_down=0

     tmp_vector_down=state_vector(j_count)%vector_down

     call MPI_REDUCE(tmp_vector_down,total_vector_down,TOTAL_DOWN_ELECTRON,&
          MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
     
     state_vector(j_count)%vector_down=total_vector_down
  end do

*/

!
! splited
!
  allocate(buff_splited_cor(5,max_number_xi),stat=ier)
  buff_splited_cor=0
  allocate(buff_splited_cd(5,TOTAL_SITE_NUMBER),stat=ier)
  buff_splited_cd=0

  buff_splited_cor=0
  call MPI_REDUCE(splited_sz0,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_sz0=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_charge0,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_charge0=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_sxsy0,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_sxsy0=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_s0,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_s0=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_super0,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_super0=buff_splited_cor

  buff_splited_cd=0
  call MPI_REDUCE(splited_charge_density0,buff_splited_cd,&
       5*TOTAL_SITE_NUMBER,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  splited_charge_density0=buff_splited_cd

  buff5=0
  call MPI_REDUCE(splited_pole0,buff5,5,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  splited_pole0=buff5



  buff_splited_cor=0
  call MPI_REDUCE(splited_sz1,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_sz1=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_charge1,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_charge1=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_sxsy1,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_sxsy1=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_s1,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_s1=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_super1,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_super1=buff_splited_cor

  buff_splited_cd=0
  call MPI_REDUCE(splited_charge_density1,buff_splited_cd,&
       5*TOTAL_SITE_NUMBER,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  splited_charge_density1=buff_splited_cd

  buff5=0
  call MPI_REDUCE(splited_pole1,buff5,5,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  splited_pole1=buff5


  buff_splited_cor=0
  call MPI_REDUCE(splited_sz2,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_sz2=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_charge2,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_charge2=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_sxsy2,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_sxsy2=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_s2,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_s2=buff_splited_cor

  buff_splited_cor=0
  call MPI_REDUCE(splited_super2,buff_splited_cor,5*max_number_xi,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  splited_super2=buff_splited_cor

  buff_splited_cd=0
  call MPI_REDUCE(splited_charge_density2,buff_splited_cd,&
       5*TOTAL_SITE_NUMBER,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  splited_charge_density2=buff_splited_cd

  buff5=0
  call MPI_REDUCE(splited_pole2,buff5,5,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,IERROR)
  splited_pole2=buff5

  
  if( MYRANK==0 ) then
     do i_count=1,max_number_xi
        splited_sz0(1,i_count)=splited_sz0(1,i_count)/dble(TOTAL_SITE_NUMBER)
        splited_sz0(2,i_count)=splited_sz0(2,i_count)/dble(TOTAL_SITE_NUMBER)
        splited_sz0(3,i_count)=splited_sz0(3,i_count)/dble(TOTAL_SITE_NUMBER)
        splited_sz0(4,i_count)=splited_sz0(4,i_count)/dble(TOTAL_SITE_NUMBER)
        splited_sz0(5,i_count)=splited_sz0(5,i_count)/dble(TOTAL_SITE_NUMBER)

        splited_charge0(1,i_count)=splited_charge0(1,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge0(2,i_count)=splited_charge0(2,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge0(3,i_count)=splited_charge0(3,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge0(4,i_count)=splited_charge0(4,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_charge0(5,i_count)=splited_charge0(5,i_count) &
             /dble(TOTAL_SITE_NUMBER)
! charge_density0は何もしなくてよい
        splited_sxsy0(1,i_count)=splited_sxsy0(1,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sxsy0(2,i_count)=splited_sxsy0(2,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sxsy0(3,i_count)=splited_sxsy0(3,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sxsy0(4,i_count)=splited_sxsy0(4,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_sxsy0(5,i_count)=splited_sxsy0(5,i_count) &
             /dble(TOTAL_SITE_NUMBER)

        splited_s0(1,i_count)=splited_sz0(1,i_count)+splited_sxsy0(1,i_count)
        splited_s0(2,i_count)=splited_sz0(2,i_count)+splited_sxsy0(2,i_count)
        splited_s0(3,i_count)=splited_sz0(3,i_count)+splited_sxsy0(3,i_count)
        splited_s0(4,i_count)=splited_sz0(4,i_count)+splited_sxsy0(4,i_count)
        splited_s0(5,i_count)=splited_sz0(5,i_count)+splited_sxsy0(5,i_count)

        splited_super0(1,i_count)=splited_super0(1,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_super0(2,i_count)=splited_super0(2,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_super0(3,i_count)=splited_super0(3,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_super0(4,i_count)=splited_super0(4,i_count) &
             /dble(TOTAL_SITE_NUMBER)
        splited_super0(5,i_count)=splited_super0(5,i_count) &
             /dble(TOTAL_SITE_NUMBER)

        if( POWER==1 .or. POWER==2 ) then
           splited_sz1(1,i_count)=splited_sz1(1,i_count) &
                /denomi_splited_pow1(1)/dble(TOTAL_SITE_NUMBER)
           splited_sz1(2,i_count)=splited_sz1(2,i_count) &
                /denomi_splited_pow1(2)/dble(TOTAL_SITE_NUMBER)
           splited_sz1(3,i_count)=splited_sz1(3,i_count) &
                /denomi_splited_pow1(3)/dble(TOTAL_SITE_NUMBER)
           splited_sz1(4,i_count)=splited_sz1(4,i_count) &
                /denomi_splited_pow1(4)/dble(TOTAL_SITE_NUMBER)
           splited_sz1(5,i_count)=splited_sz1(5,i_count) &
                /denomi_splited_pow1(5)/dble(TOTAL_SITE_NUMBER)

           splited_charge1(1,i_count)=splited_charge1(1,i_count) &
                /denomi_splited_pow1(1)/dble(TOTAL_SITE_NUMBER)
           splited_charge1(2,i_count)=splited_charge1(2,i_count) &
                /denomi_splited_pow1(2)/dble(TOTAL_SITE_NUMBER)
           splited_charge1(3,i_count)=splited_charge1(3,i_count) &
                /denomi_splited_pow1(3)/dble(TOTAL_SITE_NUMBER)
           splited_charge1(4,i_count)=splited_charge1(4,i_count) &
                /denomi_splited_pow1(4)/dble(TOTAL_SITE_NUMBER)
           splited_charge1(5,i_count)=splited_charge1(5,i_count) &
                /denomi_splited_pow1(5)/dble(TOTAL_SITE_NUMBER)

           splited_charge_density1(1,i_count)= &
                splited_charge_density1(1,i_count)/denomi_splited_pow1(1)
           splited_charge_density1(2,i_count)= &
                splited_charge_density1(2,i_count)/denomi_splited_pow1(2)
           splited_charge_density1(3,i_count)= &
                splited_charge_density1(3,i_count)/denomi_splited_pow1(3)
           splited_charge_density1(4,i_count)= &
                splited_charge_density1(4,i_count)/denomi_splited_pow1(4)
           splited_charge_density1(5,i_count)= &
                splited_charge_density1(5,i_count)/denomi_splited_pow1(5)

           splited_sxsy1(1,i_count)=splited_sxsy1(1,i_count) &
                /denomi_splited_pow1(1)/dble(TOTAL_SITE_NUMBER)
           splited_sxsy1(2,i_count)=splited_sxsy1(2,i_count) &
                /denomi_splited_pow1(2)/dble(TOTAL_SITE_NUMBER)
           splited_sxsy1(3,i_count)=splited_sxsy1(3,i_count) &
                /denomi_splited_pow1(3)/dble(TOTAL_SITE_NUMBER)
           splited_sxsy1(4,i_count)=splited_sxsy1(4,i_count) &
                /denomi_splited_pow1(4)/dble(TOTAL_SITE_NUMBER)
           splited_sxsy1(5,i_count)=splited_sxsy1(5,i_count) &
                /denomi_splited_pow1(5)/dble(TOTAL_SITE_NUMBER)

           splited_s1(1,i_count)=splited_s1(1,i_count) &
                /denomi_splited_pow1(1)/dble(TOTAL_SITE_NUMBER)
           splited_s1(2,i_count)=splited_s1(2,i_count) &
                /denomi_splited_pow1(2)/dble(TOTAL_SITE_NUMBER)
           splited_s1(3,i_count)=splited_s1(3,i_count) &
                /denomi_splited_pow1(3)/dble(TOTAL_SITE_NUMBER)
           splited_s1(4,i_count)=splited_s1(4,i_count) &
                /denomi_splited_pow1(4)/dble(TOTAL_SITE_NUMBER)
           splited_s1(5,i_count)=splited_s1(5,i_count) &
                /denomi_splited_pow1(5)/dble(TOTAL_SITE_NUMBER)

           splited_super1(1,i_count)=splited_super1(1,i_count) &
                /denomi_splited_pow1(1)/dble(TOTAL_SITE_NUMBER)
           splited_super1(2,i_count)=splited_super1(2,i_count) &
                /denomi_splited_pow1(2)/dble(TOTAL_SITE_NUMBER)
           splited_super1(3,i_count)=splited_super1(3,i_count) &
                /denomi_splited_pow1(3)/dble(TOTAL_SITE_NUMBER)
           splited_super1(4,i_count)=splited_super1(4,i_count) &
                /denomi_splited_pow1(4)/dble(TOTAL_SITE_NUMBER)
           splited_super1(5,i_count)=splited_super1(5,i_count) &
                /denomi_splited_pow1(5)/dble(TOTAL_SITE_NUMBER)
        end if

        if( POWER==2 ) then
           splited_sz2(1,i_count)=splited_sz2(1,i_count) &
                /denomi_splited_pow2(1)/dble(TOTAL_SITE_NUMBER)
           splited_sz2(2,i_count)=splited_sz2(2,i_count) &
                /denomi_splited_pow2(2)/dble(TOTAL_SITE_NUMBER)
           splited_sz2(3,i_count)=splited_sz2(3,i_count) &
                /denomi_splited_pow2(3)/dble(TOTAL_SITE_NUMBER)
           splited_sz2(4,i_count)=splited_sz2(4,i_count) &
                /denomi_splited_pow2(4)/dble(TOTAL_SITE_NUMBER)
           splited_sz2(5,i_count)=splited_sz2(5,i_count) &
                /denomi_splited_pow2(5)/dble(TOTAL_SITE_NUMBER)

           splited_charge2(1,i_count)=splited_charge2(1,i_count) &
                /denomi_splited_pow2(1)/dble(TOTAL_SITE_NUMBER)
           splited_charge2(2,i_count)=splited_charge2(2,i_count) &
                /denomi_splited_pow2(2)/dble(TOTAL_SITE_NUMBER)
           splited_charge2(3,i_count)=splited_charge2(3,i_count) &
                /denomi_splited_pow2(3)/dble(TOTAL_SITE_NUMBER)
           splited_charge2(4,i_count)=splited_charge2(4,i_count) &
                /denomi_splited_pow2(4)/dble(TOTAL_SITE_NUMBER)
           splited_charge2(5,i_count)=splited_charge2(5,i_count) &
                /denomi_splited_pow2(5)/dble(TOTAL_SITE_NUMBER)

           splited_charge_density2(1,i_count)= &
                splited_charge_density2(1,i_count)/denomi_splited_pow2(1)
           splited_charge_density2(2,i_count)= &
                splited_charge_density2(2,i_count)/denomi_splited_pow2(2)
           splited_charge_density2(3,i_count)= &
                splited_charge_density2(3,i_count)/denomi_splited_pow2(3)
           splited_charge_density2(4,i_count)= &
                splited_charge_density2(4,i_count)/denomi_splited_pow2(4)
           splited_charge_density2(5,i_count)= &
                splited_charge_density2(5,i_count)/denomi_splited_pow2(5)

           splited_sxsy2(1,i_count)=splited_sxsy2(1,i_count) &
                /denomi_splited_pow2(1)/dble(TOTAL_SITE_NUMBER)
           splited_sxsy2(2,i_count)=splited_sxsy2(2,i_count) &
                /denomi_splited_pow2(2)/dble(TOTAL_SITE_NUMBER)
           splited_sxsy2(3,i_count)=splited_sxsy2(3,i_count) &
                /denomi_splited_pow2(3)/dble(TOTAL_SITE_NUMBER)
           splited_sxsy2(4,i_count)=splited_sxsy2(4,i_count) &
                /denomi_splited_pow2(4)/dble(TOTAL_SITE_NUMBER)
           splited_sxsy2(5,i_count)=splited_sxsy2(5,i_count) &
                /denomi_splited_pow2(5)/dble(TOTAL_SITE_NUMBER)

           splited_s2(1,i_count)=splited_s2(1,i_count) &
                /denomi_splited_pow2(1)/dble(TOTAL_SITE_NUMBER)
           splited_s2(2,i_count)=splited_s2(2,i_count) &
                /denomi_splited_pow2(2)/dble(TOTAL_SITE_NUMBER)
           splited_s2(3,i_count)=splited_s2(3,i_count) &
                /denomi_splited_pow2(3)/dble(TOTAL_SITE_NUMBER)
           splited_s2(4,i_count)=splited_s2(4,i_count) &
                /denomi_splited_pow2(4)/dble(TOTAL_SITE_NUMBER)
           splited_s2(5,i_count)=splited_s2(5,i_count) &
                /denomi_splited_pow2(5)/dble(TOTAL_SITE_NUMBER)

           splited_super2(1,i_count)=splited_super2(1,i_count) &
                /denomi_splited_pow2(1)/dble(TOTAL_SITE_NUMBER)
           splited_super2(2,i_count)=splited_super2(2,i_count) &
                /denomi_splited_pow2(2)/dble(TOTAL_SITE_NUMBER)
           splited_super2(3,i_count)=splited_super2(3,i_count) &
                /denomi_splited_pow2(3)/dble(TOTAL_SITE_NUMBER)
           splited_super2(4,i_count)=splited_super2(4,i_count) &
                /denomi_splited_pow2(4)/dble(TOTAL_SITE_NUMBER)
           splited_super2(5,i_count)=splited_super2(5,i_count) &
                /denomi_splited_pow2(5)/dble(TOTAL_SITE_NUMBER)

        end if
     end do

     if( POWER==1 .or. POWER==2 ) then
        splited_pole1(1)=splited_pole1(1)/denomi_splited_pow1(1)
        splited_pole1(2)=splited_pole1(2)/denomi_splited_pow1(2)
        splited_pole1(3)=splited_pole1(3)/denomi_splited_pow1(3)
        splited_pole1(4)=splited_pole1(4)/denomi_splited_pow1(4)
        splited_pole1(5)=splited_pole1(5)/denomi_splited_pow1(5)
     end if

     if( POWER==2 ) then
        splited_pole2(1)=splited_pole2(1)/denomi_splited_pow2(1)
        splited_pole2(2)=splited_pole2(2)/denomi_splited_pow2(2)
        splited_pole2(3)=splited_pole2(3)/denomi_splited_pow2(3)
        splited_pole2(4)=splited_pole2(4)/denomi_splited_pow2(4)
        splited_pole2(5)=splited_pole2(5)/denomi_splited_pow2(5)
     end if

!
! エラーバー(MYRANK==0のみ)
!
     if( COR==1 .or. COR==2 ) then
        call error_bar_a(splited_sz0,lowest_sz0,highest_sz0,name,error)
        call error_bar_a(splited_charge0,lowest_charge0,highest_charge0,&
             name,error)
        call error_bar_b(splited_charge_density0,lowest_charge_density0,&
             highest_charge_density0,name,error)
        call error_bar_c(splited_pole0,lowest_pole0,highest_pole0,name,error)

        if( POWER==1 .or. POWER==2 ) then
           call error_bar_a(splited_sz1,lowest_sz1,highest_sz1,name,error)
           call error_bar_a(splited_charge1,lowest_charge1,highest_charge1,&
                name,error)
           call error_bar_b(splited_charge_density1,lowest_charge_density1,&
                highest_charge_density1,name,error)
           call error_bar_c(splited_pole1,lowest_pole1,highest_pole1,&
                name,error)
        end if
        if( POWER==2 ) then
           call error_bar_a(splited_sz2,lowest_sz2,highest_sz2,name,error)
           call error_bar_a(splited_charge2,lowest_charge2,highest_charge2,&
                name,error)
           call error_bar_b(splited_charge_density2,lowest_charge_density2,&
                highest_charge_density2,name,error)
           call error_bar_c(splited_pole2,lowest_pole2,highest_pole2,&
                name,error)
        end if
     end if

     if( COR==2 ) then
        call error_bar_a(splited_sxsy0,lowest_sxsy0,highest_sxsy0,name,error)
        call error_bar_a(splited_s0,lowest_s0,highest_s0,name,error)
        if( JIGEN==2 ) then
           call error_bar_a(splited_super0,lowest_super0,highest_super0,&
                name,error)
        end if
     
        if( POWER==1 .or. POWER==2 ) then
           call error_bar_a(splited_sxsy1,lowest_sxsy1,highest_sxsy1,&
                name,error)
           call error_bar_a(splited_s1,lowest_s1,highest_s1,name,error)
           if( JIGEN==2 ) then
              call error_bar_a(splited_super1,lowest_super1,highest_super1,&
                   name,error)
           end if
        end if
        if( POWER==2 ) then
           call error_bar_a(splited_sxsy2,lowest_sxsy2,highest_sxsy2,&
                name,error)
           call error_bar_a(splited_s2,lowest_s2,highest_s2,name,error)
           if( JIGEN==2 ) then
              call error_bar_a(splited_super2,lowest_super2,highest_super2,&
                   name,error)
           end if
        end if
     end if
  end if



!
!syojijyounotame
!
  if( MYRANK==0 ) then
     write(*,*) xi_field,result_pow0,result_pole0,electric_field

     if( bunbo/=0 ) then
        write(*,*) "bunsi=",bunsi
        write(*,*) "bunbo=",bunbo
        xi_field=bunsi/bunbo
        write(*,*) "xi_field0=",bunsi/bunbo
     end if

! 計算終了時刻
     time_end=MPI_WTIME()

     call output(zero_approx_energy,hf_iteration_result,&
          result_pow0,lowest_pow0,highest_pow0,result_pow1,&
          lowest_pow1,highest_pow1,result_pow2,lowest_pow2,highest_pow2,&
          lowest_sz0,highest_sz0,&
          lowest_sz1,highest_sz1,&
          lowest_sz2,highest_sz2,&
          lowest_charge0,highest_charge0,&
          lowest_charge1,highest_charge1,&
          lowest_charge2,highest_charge2,&
          lowest_sxsy0,highest_sxsy0,&
          lowest_sxsy1,highest_sxsy1,&
          lowest_sxsy2,highest_sxsy2,&
          lowest_s0,highest_s0,&
          lowest_s1,highest_s1,&
          lowest_s2,highest_s2,&
          lowest_super0,highest_super0,&
          lowest_super1,highest_super1,&
          lowest_super2,highest_super2,&
          lowest_charge_density0,highest_charge_density0,& 
          lowest_charge_density1,highest_charge_density1,&
          lowest_charge_density2,highest_charge_density2,&
          lowest_pole0,highest_pole0,&
          lowest_pole1,highest_pole1,&
          lowest_pole2,highest_pole2,&
          pow1_c0,pow1_c1,pow2_c0,pow2_c1,pow2_c2,name,error)
     call error_check(name,error)
  end if

  write(*,*) "independent_distance=",independent_distance**0.5

#ifdef DEBUG
  write(*,*) "delta_e=",delta_e/dble(MAX_MONTECARLO_SAMPLE)
#endif

  call MPI_FINALIZE(IERROR)

end program main
 

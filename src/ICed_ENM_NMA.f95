program run_ICed_ENM
implicit none
character(500) :: PDB, chain, arg, option_value, movie_filename, output_prefix, output_name
character(80), allocatable, dimension(:) :: refined_PDB,ref_PDB,ref_PDB_CA
integer :: i_natom,i_nang,cnt,i,j,k,i_ncore,i_nchains,i_nkpair,i_nCA,i_nmodes,i_full_dof
integer :: ncols,ios,i_nline,n_refined,argc,arg_index
integer :: i_nframe,i_movie_mode,i_movie_only,movie_mode_start,movie_mode_end,i_movie_nmodes
integer :: i_movie_actual
double precision, allocatable, dimension(:,:) :: arr_coordi_pre,arr_coordi,arr_coordi_IC
integer, allocatable, dimension(:,:) :: arr_atom_type_pre,arr_atom_type,ang_set
double precision, allocatable, dimension(:) :: arr_mass
real, allocatable, dimension(:) :: arr_atom_num
character(3), allocatable, dimension(:) :: arr_resid_type_pre,arr_resid_type
character(1), allocatable, dimension(:) :: arr_chain_label_pre,arr_chain_label
integer, allocatable, dimension(:) :: arr_chain_index_pre,arr_chain_index
double precision, dimension(3) :: center_mass,temp_movie_coord
double precision, allocatable, dimension(:,:) :: T_ic,T_rigid,T_cross,T_combined,T_final,T_backup
double precision, allocatable, dimension(:) :: MA,MC
double precision, allocatable, dimension(:,:,:) :: PA,PC,IA,IC
double precision :: d_cart_cutoff,C_seq,C_cart,P_seq,P_cart,sum_part,sum_total
double precision, allocatable, dimension(:,:) :: k_mat
double precision, allocatable, dimension(:) :: K_off,K_diag
double precision, allocatable, dimension(:,:) :: K_ic,R_ic,K_rigid,K_cross,K_combined,K_final,K_backup
double precision, allocatable, dimension(:) :: D,D_final
double precision, allocatable, dimension(:,:) :: Q,Q_cc_final,Q_cc_raw,Q_CA,Q_CA_orthonormal,Q_origin
double precision, allocatable, dimension(:,:) :: Q_ic,Q_cc_ic,Q_CA_raw,Q_cc_orthonormal,temp_CC_block
double precision, allocatable, dimension(:,:) :: Q_rigid,Q_cc_rigid
double precision, allocatable, dimension(:) :: var_percent_cumul,var_percent
double precision, allocatable, dimension(:,:) :: arr_chain_com,I_all
double precision, allocatable, dimension(:,:,:) :: arr_chain_inertia
double precision, allocatable, dimension(:) :: arr_theo_bfactor,arr_RMSF
double precision, allocatable, dimension(:,:,:,:) :: arr_output_deformed_coordi_CC
double precision, allocatable, dimension(:,:,:,:) :: arr_output_deformed_coordi_CA
character(len=100000) :: fmt
double precision :: d_weight
logical :: make_movie,movie_mode_on,movie_only_on
logical :: write_IC,write_CC,write_CA_raw,write_variance
integer(8) :: t_start, t_end, rate
real :: elapsed_time
external    dsygv
external    dsygvx
make_movie = .false.;movie_mode_on = .false.;movie_only_on = .false. 
write_IC = .false.;write_CC = .false.;write_CA_raw = .false.;write_variance = .false.
!******* input define ********************************************************
argc = command_argument_count()
if (argc == 1) then
    call get_command_argument(1, arg)
    if (trim(arg) == "-h" .or. trim(arg) == "--help") then
        call print_usage()
        stop 0, quiet=.true.
    endif
endif

if (argc < 2) then
    call print_usage()
    stop 1, quiet=.true.
endif

call get_command_argument(1, PDB)
call get_command_argument(2, chain)

i_ncore = 4
i_nmodes = 3
i_movie_mode = 1
d_cart_cutoff = 8.0d0
output_prefix = ""
arg_index = 3
do while (arg_index <= argc)
    call get_command_argument(arg_index, arg)
    select case (trim(arg))
    case ("--core")
        if (arg_index + 1 > argc) then
            call print_usage()
            stop 1, quiet=.true.
        endif
        call get_command_argument(arg_index + 1, option_value)
        read(option_value,*) i_ncore
        if (i_ncore < 1) then
            write(*,'(A)') "Error: --core must be a positive integer."
            stop 1, quiet=.true.
        endif
        arg_index = arg_index + 2
    case ("--mode")
        if (arg_index + 1 > argc) then
            call print_usage()
            stop 1, quiet=.true.
        endif
        call get_command_argument(arg_index + 1, option_value)
        read(option_value,*) i_nmodes
        if (i_nmodes < 0) then
            write(*,'(A)') "Error: --mode must be zero or a positive integer."
            stop 1, quiet=.true.
        endif
        arg_index = arg_index + 2
    case ("--cutoff")
        if (arg_index + 1 > argc) then
            call print_usage()
            stop 1, quiet=.true.
        endif
        call get_command_argument(arg_index + 1, option_value)
        read(option_value,*) d_cart_cutoff
        if (d_cart_cutoff < 0.0d0) then
            write(*,'(A)') "Error: --cutoff must be a positive double precision number."
            stop 1, quiet=.true.
        endif
        arg_index = arg_index + 2
    case ("--out-prefix")
        if (arg_index + 1 > argc) then
            call print_usage()
            stop 1, quiet=.true.
        endif
        call get_command_argument(arg_index + 1, option_value)
        output_prefix = trim(adjustl(option_value))
        if (len_trim(output_prefix) > 0) then
            if (output_prefix(len_trim(output_prefix):len_trim(output_prefix)) /= "_" .and. &
                output_prefix(len_trim(output_prefix):len_trim(output_prefix)) /= "-" .and. &
                output_prefix(len_trim(output_prefix):len_trim(output_prefix)) /= "/") then
                output_prefix = trim(output_prefix)//"_"
            endif
        endif
        arg_index = arg_index + 2
    case ("--write-IC")
        write_IC = .true.
        arg_index = arg_index + 1
    case ("--write-CC")
        write_CC = .true.
        arg_index = arg_index + 1
    case ("--write-CA-raw")
        write_CA_raw = .true.
        arg_index = arg_index + 1
    case ("--write-variance")
        write_variance = .true.
        arg_index = arg_index + 1
    case ("--write-all")
        write_IC = .true.
        write_CC = .true.
        write_CA_raw = .true.
        write_variance = .true.
        arg_index = arg_index + 1
    case ("--movie")
        if (arg_index + 2 > argc) then
            call print_usage()
            stop 1, quiet=.true.
        endif
        make_movie = .true.
        call get_command_argument(arg_index + 1, option_value)
        read(option_value,*) d_weight
        if (d_weight < 0.0d0) then
            write(*,'(A)') "Error: --movie weight must be a positive double precision number."
            stop 1, quiet=.true.
        endif
        call get_command_argument(arg_index + 2, option_value)
        read(option_value,*) i_nframe
        if (i_nframe < 3 .or. mod(i_nframe, 2) == 0) then
            write(*,'(A)') "Error: --movie frame number must be an odd integer >= 3."
            stop 1, quiet=.true.
        endif
        arg_index = arg_index + 3
    case ("--movie-mode")
        if (arg_index + 1 > argc) then
            call print_usage()
            stop 1, quiet=.true.
        endif
        movie_mode_on = .true.
        call get_command_argument(arg_index + 1, option_value)
        read(option_value,*) i_movie_mode
        if (i_movie_mode < 1) then
            write(*,'(A)') "Error: --movie-mode must be a positive integer."
            stop 1, quiet=.true.
        endif
        arg_index = arg_index + 2
    case ("--movie-only")
        if (arg_index + 1 > argc) then
            call print_usage()
            stop 1, quiet=.true.
        endif
        movie_only_on = .true.
        call get_command_argument(arg_index + 1, option_value)
        read(option_value,*) i_movie_only
        if (i_movie_only < 1) then
            write(*,'(A)') "Error: --movie-only must be a positive integer."
            stop 1, quiet=.true.
        endif
        arg_index = arg_index + 2
    case ("-h", "--help")
        call print_usage()
        stop 0, quiet=.true.
    case default
        call print_usage()
        stop 1, quiet=.true.
    end select
enddo

call omp_set_num_threads(i_ncore)
call openblas_set_num_threads(i_ncore)
i_nline = 0
OPEN (1, file = PDB)
DO
    READ (1,*, END=10)
    i_nline = i_nline + 1
END DO
10 CLOSE (1)
if (.not.allocated(refined_PDB)) allocate(refined_PDB(i_nline))
call refine_PDB(PDB,chain,refined_PDB,i_nline,n_refined)
if (.not.allocated(arr_coordi_pre)) allocate(arr_coordi_pre(n_refined,3))
if (.not.allocated(arr_atom_type_pre)) allocate(arr_atom_type_pre(n_refined,3))
if (.not.allocated(arr_resid_type_pre)) allocate(arr_resid_type_pre(n_refined))
if (.not.allocated(arr_chain_label_pre)) allocate(arr_chain_label_pre(n_refined))
if (.not.allocated(arr_chain_index_pre)) allocate(arr_chain_index_pre(n_refined))
if (.not.allocated(ref_PDB)) allocate(ref_PDB(n_refined+1))
if (.not.allocated(ref_PDB_CA)) allocate(ref_PDB_CA(n_refined+1))

call read_PDB(refined_PDB,n_refined,arr_coordi_pre,arr_atom_type_pre,arr_resid_type_pre,arr_chain_label_pre,&
arr_chain_index_pre,i_natom,i_nCA,ref_PDB,ref_PDB_CA)

if (.not.allocated(arr_coordi)) allocate(arr_coordi(i_natom,3))
if (.not.allocated(arr_atom_type)) allocate(arr_atom_type(i_natom,3))
if (.not.allocated(arr_resid_type)) allocate(arr_resid_type(i_natom))
if (.not.allocated(arr_chain_label)) allocate(arr_chain_label(i_natom))
if (.not.allocated(arr_chain_index)) allocate(arr_chain_index(i_natom))
arr_coordi=arr_coordi_pre(1:i_natom,:)
arr_atom_type=arr_atom_type_pre(1:i_natom,:)
arr_resid_type=arr_resid_type_pre(1:i_natom)
arr_chain_label=arr_chain_label_pre(1:i_natom)
arr_chain_index=arr_chain_index_pre(1:i_natom)
deallocate(arr_coordi_pre);deallocate(arr_atom_type_pre)
deallocate(arr_resid_type_pre);deallocate(arr_chain_label_pre);deallocate(arr_chain_index_pre)

write(*,'(A56,I3)') "The num of chains considering missing parts and HETATM: ", arr_chain_index(i_natom)
i_nchains=arr_chain_index(i_natom)
if (.not.allocated(arr_mass)) allocate( arr_mass(i_natom) )
call make_mass_array(arr_atom_type, arr_resid_type, arr_mass)
! re-centering the coordi with mass information *************************************
center_mass=0.0
do i=1,i_natom
    center_mass=center_mass+arr_mass(i)*arr_coordi(i,:)
enddo
center_mass=center_mass/sum(arr_mass(:))
do i=1,i_natom
    arr_coordi(i,:)=arr_coordi(i,:)-center_mass
enddo

cnt=0
do i=1,i_natom
    if (arr_atom_type(i,2)==2) then ! Check phi angle not in PROLINE
        if (i >= 3 .AND. i <= i_natom-1) then ! not in the first residue and not as the last atom
            if (arr_atom_type(i+1,2)==4 .AND. i <= i_natom-2) then ! this residue has a virtual bead
                if (arr_chain_label(i-2)==arr_chain_label(i+2) .AND. arr_atom_type(i-2,1)+1==arr_atom_type(i+2,1) .AND. &
                    (arr_atom_type(i-2,2)==3 .OR. arr_atom_type(i-2,2)==33) .AND. arr_atom_type(i+2,2)==3 ) then ! previous residue's C and this residue's C should be in the same chain and be connected 
                    cnt=cnt+1
                endif
            else if (arr_atom_type(i+1,2)==3) then ! this residue has no virtual bead
                if (arr_chain_label(i-2)==arr_chain_label(i+1) .AND. arr_atom_type(i-2,1)+1==arr_atom_type(i+1,1) .AND. &
                    (arr_atom_type(i-2,2)==3 .OR. arr_atom_type(i-2,2)==33)) then ! previous residue's C and this residue's C should be in the same chain and be connected 
                    cnt=cnt+1
                endif
            endif
        endif
    else if (arr_atom_type(i,2)==3 .OR. arr_atom_type(i,2)==33) then ! Check psi angle
        if (i >= 3 .AND. i <= i_natom-1) then 
            if ((arr_atom_type(i-1,2)==4 .OR. arr_atom_type(i-1,2)==44) .AND. i >= 4) then ! this residue has a virtual beada
                if (arr_chain_label(i-3)==arr_chain_label(i+1) .AND. arr_atom_type(i-3,1)+1==arr_atom_type(i+1,1) .AND. & 
                (arr_atom_type(i-3,2)==1 .OR. arr_atom_type(i-3,2)==11) .AND. (arr_atom_type(i+1,2)==1 &
                    .OR. arr_atom_type(i+1,2)==11)) then
                    cnt=cnt+1
                endif
            else if (arr_atom_type(i-1,2)==2 .OR. arr_atom_type(i-1,2)==22) then
                if (arr_chain_label(i-2)==arr_chain_label(i+1) .AND. arr_atom_type(i-2,1)+1==arr_atom_type(i+1,1) .AND. &
                (arr_atom_type(i-2,2)==1 .OR. arr_atom_type(i-2,2)==11) .AND. (arr_atom_type(i+1,2)==1 &
                    .OR. arr_atom_type(i+1,2)==11)) then
                    cnt=cnt+1
                endif
            endif
        endif
    endif
enddo
i_nang=cnt

if (.not.allocated(arr_coordi_IC)) allocate( arr_coordi_IC(i_nang,6))
if (.not.allocated(arr_atom_num)) allocate( arr_atom_num(i_nang))
call coordi_IC(arr_coordi, arr_atom_type, arr_chain_index, arr_coordi_IC, arr_atom_num) 
if(.not.allocated(ang_set)) allocate( ang_set(i_nang*(i_nang-1)/2,2))
call make_angle_set (arr_atom_num, ang_set)
if(.not.allocated(T_ic)) allocate( T_ic(i_nang,i_nang))
if(.not.allocated(PA)) then 
allocate( PA(3,3,i_nang)); allocate( PC(3,3,i_nang))
allocate( IA(3,3,i_nang)); allocate( IC(3,3,i_nang))
allocate( MA(i_nang)); allocate( MC(i_nang)); allocate( I_all(3,3) )
endif
!!!!!!!!!!!!!Mass matrix calculation!!!!!!!!!!!!!!!!!!!!!!
call system_clock(count_rate=rate)

call system_clock(t_start)
call final_T_ic(arr_mass,arr_coordi,arr_coordi_IC,arr_atom_num,ang_set,PA,PC,IA,IC,MA,MC,I_all,T_ic)
call system_clock(t_end)
elapsed_time = real(t_end - t_start) / real(rate)
print *, 'compute_T_ic: ', elapsed_time
!!!!!!!!!!!!!T_rigid calculation!!!!!!!!!!!!!!!!!!!!!!
if (i_nchains > 1) then ! for the multi-chain proteins
    allocate(arr_chain_com(i_nchains,3))
    allocate(arr_chain_inertia(3,3,i_nchains))
    call system_clock(t_start)
    call compute_chain_rigid(arr_mass,arr_coordi,arr_chain_index,arr_chain_com,arr_chain_inertia,i_nchains)
    call system_clock(t_end)
    elapsed_time = real(t_end - t_start) / real(rate)
    print *, 'compute_chain_rigid: ', elapsed_time
    allocate(T_rigid(6 * i_nchains, 6 * i_nchains))
    call system_clock(t_start)
    call compute_T_rigid(arr_mass, arr_chain_index, arr_chain_inertia, i_nchains, T_rigid)
    call system_clock(t_end)
    elapsed_time = real(t_end - t_start) / real(rate)
    print *, 'compute_T_rigid: ', elapsed_time
    allocate(T_cross(i_nang, 6 * i_nchains))
    call system_clock(t_start)
    call compute_T_cross(arr_mass,arr_coordi_IC,arr_coordi,arr_chain_index,arr_chain_com,&
    i_nchains,arr_atom_num,PA,PC,IA,IC,MA,MC,I_all,T_cross)
    call system_clock(t_end)
    elapsed_time = real(t_end - t_start) / real(rate)
    print *, 'compute_T_cross: ', elapsed_time
    !!!!!!!!!!!!!Combine T_ic and T_rigid!!!!!!!!!!!!!!!!!!!!!!
    allocate(T_combined(i_nang+6*i_nchains,i_nang+6*i_nchains))
    call combine_IC_and_rigid(T_ic,T_rigid,T_cross,T_combined,i_nang)
    allocate(T_final(i_nang+6*i_nchains,i_nang+6*i_nchains))
    T_final=T_combined
else
    allocate(T_final(i_nang,i_nang))
    T_final=T_ic
endif
!!!!!!!!!!!!!Hessian matrix calculation!!!!!!!!!!!!!!!!!!!!!!
write(*,'(A21,F6.1)') "Cutoff distance (A): ", d_cart_cutoff

C_seq=286.4;C_cart=147.6;P_seq=-1.92;P_cart=-1.66
if (.not.allocated(k_mat)) allocate( k_mat(i_natom,i_natom))
call system_clock(t_start)
call make_k_mat (arr_coordi,k_mat,i_nkpair,d_cart_cutoff,&
                C_seq,C_cart,P_seq,P_cart,arr_atom_type,arr_chain_index)
call system_clock(t_end)
elapsed_time = real(t_end - t_start) / real(rate)
print *, 'compute_k_mat: ', elapsed_time
if (.not.allocated(R_ic)) allocate( R_ic(6,6*(i_natom**2+5*i_natom+2)/2))
call system_clock(t_start)
call make_R_ic(arr_coordi,k_mat,R_ic)
call system_clock(t_end)
elapsed_time = real(t_end - t_start) / real(rate)
print *, 'compute_R_ic: ', elapsed_time
if (.not.allocated(K_off)) allocate( K_off(ubound(ang_set,1)))
call system_clock(t_start)
call make_K_off(i_natom,arr_coordi_IC,arr_atom_num,ang_set,R_ic,K_off)
call system_clock(t_end)
elapsed_time = real(t_end - t_start) / real(rate)
print *, 'compute_K_off: ', elapsed_time
if (.not.allocated(K_diag)) allocate( K_diag(ubound(arr_coordi_IC,1)))
call system_clock(t_start)
call make_K_diag(i_natom,arr_coordi_IC,arr_atom_num,R_ic,K_diag)
call system_clock(t_end)
elapsed_time = real(t_end - t_start) / real(rate)
print *, 'compute_K_diag: ', elapsed_time
if (.not.allocated(K_ic)) allocate( K_ic(ubound(arr_coordi_IC,1),ubound(arr_coordi_IC,1)))
do i=1,ubound(K_off,1)
    K_ic(ang_set(i,1),ang_set(i,2))=K_off(i)
    K_ic(ang_set(i,2),ang_set(i,1))=K_off(i)
enddo
do i=1,ubound(K_diag,1)
    K_ic(i,i)=K_diag(i)
enddo
deallocate(K_off); deallocate(K_diag); deallocate(R_ic)
!*****************************************************************************
!!!!!!!!!!!!!K_rigid calculation!!!!!!!!!!!!!!!!!!!!!!
if (i_nchains > 1) then ! for the multi-chain proteins
    allocate(K_rigid(6 * i_nchains, 6 * i_nchains))
    call system_clock(t_start)
    call make_K_rigid(arr_coordi,arr_chain_index,i_nchains,&
                            k_mat,arr_chain_com,K_rigid)
    call system_clock(t_end)
    elapsed_time = real(t_end - t_start) / real(rate)
    print *, 'compute_K_rigid: ', elapsed_time                        
    allocate(K_cross(i_nang, 6 * i_nchains))
    call system_clock(t_start)
    call make_K_cross(arr_coordi,arr_chain_index,i_nchains,arr_coordi_IC,arr_atom_num,&
                        k_mat,arr_chain_com,K_cross)
    call system_clock(t_end)
    elapsed_time = real(t_end - t_start) / real(rate)
    print *, 'compute_K_cross: ', elapsed_time                        
!!!!!!!!!!!!!Combine K_ic and K_rigid!!!!!!!!!!!!!!!!!!!!!!
    allocate(K_combined(i_nang+6*i_nchains,i_nang+6*i_nchains))
    call combine_IC_and_rigid(K_ic,K_rigid,K_cross,K_combined,i_nang)
    allocate(K_final(i_nang+6*i_nchains,i_nang+6*i_nchains))
    K_final=K_combined
else
    allocate(K_final(i_nang,i_nang))
    K_final=K_ic
endif

!compute eig-values and -vectors
if (i_nchains > 1) then
    i_full_dof = i_nang + 6 * i_nchains
else
    i_full_dof = i_nang
endif   

if (i_nmodes /= 0) then
    write(*,'(A44,I5)') "The number of normal modes to be calculated: ", i_nmodes
else
    i_nmodes = i_full_dof
    if (i_nchains > 1) then
        write(*,'(A44,I5,A17)') "The number of normal modes to be calculated: ", i_full_dof - 6, " (Full spectrum) "
    else
        write(*,'(A44,I5,A17)') "The number of normal modes to be calculated: ", i_full_dof, " (Full spectrum) "
    endif
endif

if (i_nchains > 1) then
    if (i_nmodes /= i_full_dof) then
        allocate(Q(i_full_dof,i_nmodes + 6)); allocate(D(i_nmodes + 6))
    else !if (i_nmodes == i_full_dof)
        allocate(Q(i_full_dof,i_nmodes)); allocate(D(i_nmodes))
    endif
else
    allocate(Q(i_full_dof,i_nmodes)); allocate(D(i_nmodes))
endif

allocate ( T_backup(ubound(T_final,1),ubound(T_final,2)))
allocate ( K_backup(ubound(K_final,1),ubound(K_final,2)))
K_backup = K_final
T_backup = T_final
call system_clock(t_start)
if (i_nmodes /= i_full_dof) then ! for the calculation of the part of normal modes
    if (i_nchains > 1 ) then ! for the multi-chain system
        call eig_func_part_dsygvx(K_final,T_final,i_full_dof,Q,D,i_nmodes + 6)
    else ! for the single-chain system
        call eig_func_part_dsygvx(K_final,T_final,i_full_dof,Q,D,i_nmodes)
    endif
else ! for the calculation of the full set of normal modes
    call eig_func(K_final,T_final,i_full_dof,Q,D) ! Q: EIGENVECTOR // D: EIGENVALUE
endif
call system_clock(t_end)
elapsed_time = real(t_end - t_start) / real(rate)
print *, 'compute_eig_problem: ', elapsed_time                   
!!!! Transformation of Q-IC dofs to Q-CC dofs
call system_clock(t_start)
if (i_nchains > 1) then ! for the multi-chain system
    if (i_nmodes /= i_full_dof) then ! for the calculation of the part of normal modes
        allocate( Q_origin(i_nang + i_nchains * 6, i_nmodes)) ! Remove redundant DOF
        allocate( Q_ic(i_nang,i_nmodes)); allocate( Q_cc_ic(3*i_natom,i_nmodes));
        allocate( Q_rigid(i_nchains * 6,i_nmodes)); allocate( Q_cc_rigid(3*i_natom,i_nmodes))
        allocate( Q_cc_final(3*i_natom,i_nmodes)); allocate( Q_cc_raw(3*i_natom,i_nmodes));
        allocate( D_final(i_nmodes))
        D_final(:) = D(7:)
        Q_origin(:,:) = Q(:,7:)
        Q_ic(:,:) = Q(:i_nang,7:)
        Q_rigid(:,:) = Q(i_nang+1:,7:)
        call ICtoCC(arr_mass,Q_ic,arr_coordi,arr_coordi_IC,arr_atom_num,PA &
            ,PC,IA,IC,MA,MC,Q_cc_ic,i_nmodes)    
        call RIGIDtoCC(Q_rigid,Q_cc_rigid,arr_chain_index,arr_chain_com,&
            arr_coordi)
    else ! for the calculation of the full set of normal modes (i_nmodes = i_full_dof)
        allocate( Q_origin(i_nang + i_nchains * 6, i_full_dof - 6)) ! Remove redundant DOF
        allocate( Q_ic(i_nang,i_full_dof - 6)); allocate( Q_cc_ic(3*i_natom,i_full_dof - 6))
        allocate( Q_rigid(i_nchains * 6,i_full_dof - 6)); allocate( Q_cc_rigid(3*i_natom,i_full_dof - 6))
        allocate( Q_cc_final(3*i_natom,i_full_dof - 6)); allocate( Q_cc_raw(3*i_natom,i_full_dof - 6));
        allocate( D_final(i_full_dof - 6))
        D_final(:) = D(7:)
        Q_origin(:,:) = Q(:,7:)
        Q_ic(:,:) = Q(:i_nang,7:)
        Q_rigid(:,:) = Q(i_nang+1:,7:)
        
        call ICtoCC(arr_mass,Q_ic,arr_coordi,arr_coordi_IC,arr_atom_num,PA &
            ,PC,IA,IC,MA,MC,Q_cc_ic,i_nmodes-6)    
        call RIGIDtoCC(Q_rigid,Q_cc_rigid,arr_chain_index,arr_chain_com,&
            arr_coordi)
    endif
    Q_cc_final = Q_cc_rigid + Q_cc_ic
    
else ! for the single-chain system
    allocate( Q_origin(i_nang,i_nmodes)); ! Remove redundant DOF
    allocate( Q_ic(i_nang,i_nmodes)); allocate( Q_cc_ic(3*i_natom,i_nmodes))
    allocate( Q_cc_final(3*i_natom,i_nmodes)); allocate( Q_cc_raw(3*i_natom,i_nmodes));
    allocate( D_final(i_nmodes))
    D_final = D
    Q_origin = Q
    Q_ic = Q
    call ICtoCC(arr_mass,Q_ic,arr_coordi,arr_coordi_IC,arr_atom_num,PA &
        ,PC,IA,IC,MA,MC,Q_cc_ic,i_nmodes)    
    Q_cc_final = Q_cc_ic
endif
allocate(Q_CA_raw(3*i_nCA,ubound(Q_cc_final,2)))
call CCtoCA(Q_cc_final,Q_CA_raw,arr_atom_type)

Q_cc_raw = Q_cc_final

allocate(Q_cc_orthonormal(ubound(Q_cc_raw,1),ubound(Q_cc_raw,2)))
do i=1,ubound(Q_cc_orthonormal,2)
    Q_cc_orthonormal(:,i)=Q_cc_raw(:,i)/norm2(Q_cc_raw(:,i))
    do j=1,i-1
        Q_cc_orthonormal(:,i)=Q_cc_orthonormal(:,i)-&
        dot_product(Q_cc_orthonormal(:,i),Q_cc_orthonormal(:,j))*Q_cc_orthonormal(:,j)
    enddo
    Q_cc_orthonormal(:,i)=Q_cc_orthonormal(:,i)/norm2(Q_cc_orthonormal(:,i))
enddo
call system_clock(t_end)
elapsed_time = real(t_end - t_start) / real(rate)
print *, 'transform from IC to CC: ', elapsed_time                        

allocate(Q_CA(3*i_nCA,ubound(Q_cc_final,2)))
Q_CA=Q_CA_raw
allocate(Q_CA_orthonormal(3*i_nCA,ubound(Q_cc_final,2)))
do i=1,ubound(Q_CA_orthonormal,2)
    Q_CA(:,i)=Q_CA(:,i)/norm2(Q_CA(:,i))
    Q_CA_orthonormal(:,i)=Q_CA(:,i)
    do j=1,i-1
        Q_CA_orthonormal(:,i)=Q_CA_orthonormal(:,i)-&
        dot_product(Q_CA_orthonormal(:,i),Q_CA_orthonormal(:,j))*Q_CA_orthonormal(:,j)
    enddo
    Q_CA_orthonormal(:,i)=Q_CA_orthonormal(:,i)/norm2(Q_CA_orthonormal(:,i))
enddo

allocate(arr_theo_bfactor(i_nCA))
allocate(arr_RMSF(i_nCA))
call cal_theo_bfactor(Q_CA_raw,D_final,arr_theo_bfactor,arr_RMSF)

allocate( var_percent_cumul(ubound(Q_cc_final,2)))
allocate( var_percent(ubound(Q_cc_final,2)))
sum_total=0.0
sum_part=0.0
do i=1,ubound(Q_cc_final,2)
    sum_total=sum_total+1/D_final(i)
enddo
do i=1,ubound(Q_cc_final,2)
    sum_part=sum_part+1/D_final(i)
    var_percent_cumul(i)=(sum_part/sum_total)*100
    var_percent(i)=(1/(D_final(i)*sum_total))*100
enddo
write(*,'(A26,I5)') "Internal coordinate DOFs: ", i_nang
if (i_nchains > 1) then
    write(*,'(A26,I5)') "  Rigid-body motion DOFs: ", 6 * i_nchains
    write(*,'(A26,I5,A27)') "              Total DOFs: ", i_full_dof - 6, &
                                    " (IC_DOFs + Rigid_DOFs - 6)"
endif

allocate(temp_CC_block(3,ubound(Q_cc_raw,2)))
i = 1
do while (i < i_natom)
    if (arr_chain_index(i)==arr_chain_index(i+1) .AND. arr_atom_type(i,1)==arr_atom_type(i+1,1) .AND. &
        (arr_atom_type(i,2)==4 .OR. arr_atom_type(i,2)==44) .AND. &
        (arr_atom_type(i+1,2)==3 .OR. arr_atom_type(i+1,2)==33)) then
        temp_CC_block = Q_cc_raw(3*i-2:3*i,:)
        Q_cc_raw(3*i-2:3*i,:) = Q_cc_raw(3*(i+1)-2:3*(i+1),:)
        Q_cc_raw(3*(i+1)-2:3*(i+1),:) = temp_CC_block

        temp_CC_block = Q_cc_orthonormal(3*i-2:3*i,:)
        Q_cc_orthonormal(3*i-2:3*i,:) = Q_cc_orthonormal(3*(i+1)-2:3*(i+1),:)
        Q_cc_orthonormal(3*(i+1)-2:3*(i+1),:) = temp_CC_block
        i = i + 2
    else
        i = i + 1
    endif
enddo
deallocate(temp_CC_block)

!*****Movie generation******************************************************
if (.not. make_movie .and. (movie_mode_on .or. movie_only_on)) then
    write(*,'(A)') "Error: --movie-mode/--movie-only requires --movie."
    stop 1, quiet=.true.
endif
if (movie_mode_on .and. movie_only_on) then
    write(*,'(A)') "Error: use either --movie-mode or --movie-only, not both."
    stop 1, quiet=.true.
endif

if (make_movie) then
    if (movie_only_on) then
        if (i_movie_only > ubound(Q_ic,2)) then
            write(*,'(A)') "Error: --movie-only N must be equal or lower than the number of calculated modes."
            stop 1, quiet=.true.
        endif
        movie_mode_start = i_movie_only
        movie_mode_end = i_movie_only
    else
        if (i_movie_mode > ubound(Q_ic,2)) then
            write(*,'(A)') "Error: --movie-mode N must be equal or lower than the number of calculated modes."
            stop 1, quiet=.true.
        endif
        movie_mode_start = 1
        movie_mode_end = i_movie_mode
    endif
    i_movie_nmodes = movie_mode_end - movie_mode_start + 1
    if (.not.allocated(arr_output_deformed_coordi_CA)) allocate(arr_output_deformed_coordi_CA(i_nCA,3,&
                                                                i_nframe,i_movie_nmodes))
    if (.not.allocated(arr_output_deformed_coordi_CC)) allocate(arr_output_deformed_coordi_CC(i_natom,3,&
                                                                i_nframe,i_movie_nmodes))
    if (i_nchains > 1 ) then
        call gen_multi_chain_path_frame(Q_ic,Q_rigid,arr_coordi,arr_mass,&
                                        arr_chain_index,arr_atom_num,arr_atom_type,&
                                        d_weight,movie_mode_start,movie_mode_end,arr_output_deformed_coordi_CC)
    else
        call gen_single_chain_path_frame(Q_ic,arr_coordi,arr_mass,&
                                        arr_atom_num,arr_atom_type,d_weight,&
                                        movie_mode_start,movie_mode_end,arr_output_deformed_coordi_CC)
    endif
    do k=1,i_movie_nmodes
        do j=1,i_nframe
            i = 1
            do while (i < i_natom)
                if (arr_chain_index(i)==arr_chain_index(i+1) .AND. &
                    arr_atom_type(i,1)==arr_atom_type(i+1,1) .AND. &
                    (arr_atom_type(i,2)==4 .OR. arr_atom_type(i,2)==44) .AND. &
                    (arr_atom_type(i+1,2)==3 .OR. arr_atom_type(i+1,2)==33)) then
                    temp_movie_coord = arr_output_deformed_coordi_CC(i,:,j,k)
                    arr_output_deformed_coordi_CC(i,:,j,k) = arr_output_deformed_coordi_CC(i+1,:,j,k)
                    arr_output_deformed_coordi_CC(i+1,:,j,k) = temp_movie_coord
                    i = i + 2
                else
                    i = i + 1
                endif
            enddo
            cnt = 0
            do i=1,i_natom
                if (arr_atom_type(i,2)==2 .OR. arr_atom_type(i,2)==22) then
                    cnt = cnt + 1
                    arr_output_deformed_coordi_CA(cnt,:,j,k) = arr_output_deformed_coordi_CC(i,:,j,k)
                endif
            enddo
        enddo
    enddo

    do k=1,i_movie_nmodes
        i_movie_actual = movie_mode_start + k - 1
        write(movie_filename,'(A,"movie_mode",I0,"_CC.pdb")') trim(output_prefix), i_movie_actual
        open(115, file=trim(movie_filename), status='unknown')
        do j=1,i_nframe
            write(115,'(A5,I9)') 'MODEL', j
            do i=1,i_natom
                write(115,'(A30,3F8.3,A26)') ref_PDB(i)(1:30), &
                    arr_output_deformed_coordi_CC(i,1:3,j,k), ref_PDB(i)(55:80)
            enddo
            write(115,'(A6)') 'ENDMDL'
        enddo
        write(115,'(A3)') 'END'
        close(115)

        write(movie_filename,'(A,"movie_mode",I0,"_CA.pdb")') trim(output_prefix), i_movie_actual
        open(116, file=trim(movie_filename), status='unknown')
        do j=1,i_nframe
            write(116,'(A5,I9)') 'MODEL', j
            do i=1,i_nCA
                write(116,'(A30,3F8.3,A26)') ref_PDB_CA(i)(1:30), &
                    arr_output_deformed_coordi_CA(i,1:3,j,k), ref_PDB_CA(i)(55:80)
            enddo
            write(116,'(A6)') 'ENDMDL'
        enddo
        write(116,'(A3)') 'END'
        close(116)
    enddo
endif

!****compute variance percent***********************************************
!*****************************************************************************
output_name = trim(output_prefix)//'eval.txt'
open(107, file=trim(output_name), status='unknown')
output_name = trim(output_prefix)//'RMSF_CA.txt'
open(110, file=trim(output_name), status='unknown')
output_name = trim(output_prefix)//'evec_CA.txt'
open(112, file=trim(output_name), status='unknown')
output_name = trim(output_prefix)//'ref_CA.pdb'
open(114, file=trim(output_name), status='unknown')

if (write_IC) then
    output_name = trim(output_prefix)//'evec_IC.txt'
    open(105, file=trim(output_name), status='unknown')
    ncols = ubound(Q_origin,2)
    write(fmt, '(A,I0,A)') '(', ncols, 'F14.10)'
    do i=1,ubound(Q_origin,1)
        write(105,trim(fmt)) Q_origin(i,:)
    enddo
    close(105)
endif

if (write_CC) then
    output_name = trim(output_prefix)//'evec_CC_raw.txt'
    open(104, file=trim(output_name), status='unknown')
    ncols = ubound(Q_cc_raw, 2)
    write(fmt, '(A,I0,A)') '(', ncols, 'F14.10)'
    do i=1,3*i_natom
        write(104,trim(fmt)) Q_cc_raw(i,:)
    enddo
    close(104)

    output_name = trim(output_prefix)//'evec_CC.txt'
    open(111, file=trim(output_name), status='unknown')
    ncols = ubound(Q_cc_orthonormal, 2)
    write(fmt, '(A,I0,A)') '(', ncols, 'F10.6)'
    do i=1,3*i_natom
        write(111,trim(fmt)) Q_cc_orthonormal(i,:)
    enddo
    close(111)

    output_name = trim(output_prefix)//'ref_CC.pdb'
    open(113, file=trim(output_name), status='unknown')
    do i=1,i_natom+1
        write(113,'(A80)') ref_pdb(i)
    enddo
    close(113)
endif

if (write_variance) then
    output_name = trim(output_prefix)//'variance_cumulative.txt'
    open(106, file=trim(output_name), status='unknown')
    output_name = trim(output_prefix)//'variance.txt'
    open(108, file=trim(output_name), status='unknown')
endif

do i=1,ubound(Q_cc_orthonormal,2)
    write(107,'(G20.12)') D_final(i)
    if (write_variance) then
        write(106,'(F7.2)') var_percent_cumul(i)
        write(108,'(F7.2)') var_percent(i)
    endif
enddo
if (write_variance) then
    close(106)
    close(108)
endif

if (write_CA_raw) then
    output_name = trim(output_prefix)//'evec_CA_raw.txt'
    open(109, file=trim(output_name), status='unknown')
    ncols = ubound(Q_CA_raw, 2)
    write(fmt, '(A,I0,A)') '(', ncols, 'F14.10)'
    do i=1,3*i_nCA
        write(109,trim(fmt)) Q_CA_raw(i,:)
    enddo
    close(109)
endif

ncols = ubound(Q_CA_orthonormal, 2)
write(fmt, '(A,I0,A)') '(', ncols, 'F10.6)'
do i=1,3*i_nCA
    write(112,fmt) Q_CA_orthonormal(i,:)
enddo

do i=1,i_nCA
    write(110,'(F10.6)') arr_RMSF(i)
enddo

do i=1,i_nCA+1
    write(114,'(A80)') ref_pdb_CA(i)
enddo
close(107);close(110);close(112);close(114)
!************************************************************************
! Deallocate all allocatable arrays
!************************************************************************
if (allocated(arr_coordi)) deallocate(arr_coordi)
if (allocated(arr_atom_type)) deallocate(arr_atom_type)
if (allocated(arr_resid_type)) deallocate(arr_resid_type)
if (allocated(arr_chain_label)) deallocate(arr_chain_label)
if (allocated(arr_chain_index)) deallocate(arr_chain_index)
if (allocated(arr_mass)) deallocate(arr_mass)
if (allocated(arr_coordi_IC)) deallocate(arr_coordi_IC)
if (allocated(arr_atom_num)) deallocate(arr_atom_num)
if (allocated(ang_set)) deallocate(ang_set)
if (allocated(PA)) deallocate(PA)
if (allocated(PC)) deallocate(PC)
if (allocated(IA)) deallocate(IA)
if (allocated(IC)) deallocate(IC)
if (allocated(MA)) deallocate(MA)
if (allocated(MC)) deallocate(MC)
if (allocated(I_all)) deallocate(I_all)
if (allocated(T_ic)) deallocate(T_ic)
if (allocated(T_rigid)) deallocate(T_rigid)
if (allocated(T_cross)) deallocate(T_cross)
if (allocated(T_combined)) deallocate(T_combined)
if (allocated(T_final)) deallocate(T_final)
if (allocated(k_mat)) deallocate(k_mat)
if (allocated(K_ic)) deallocate(K_ic)
if (allocated(K_rigid)) deallocate(K_rigid)
if (allocated(K_cross)) deallocate(K_cross)
if (allocated(K_combined)) deallocate(K_combined)
if (allocated(K_final)) deallocate(K_final)
if (allocated(K_off)) deallocate(K_off)
if (allocated(K_diag)) deallocate(K_diag)
if (allocated(R_ic)) deallocate(R_ic)
if (allocated(D)) deallocate(D)
if (allocated(D_final)) deallocate(D_final)
if (allocated(Q)) deallocate(Q)
if (allocated(Q_ic)) deallocate(Q_ic)
if (allocated(Q_rigid)) deallocate(Q_rigid)
if (allocated(Q_cc_ic)) deallocate(Q_cc_ic)
if (allocated(Q_cc_rigid)) deallocate(Q_cc_rigid)
if (allocated(Q_cc_final)) deallocate(Q_cc_final)
if (allocated(Q_CA)) deallocate(Q_CA)
if (allocated(Q_CA_orthonormal)) deallocate(Q_CA_orthonormal)
if (allocated(var_percent_cumul)) deallocate(var_percent_cumul)
if (allocated(var_percent)) deallocate(var_percent)
if (allocated(arr_chain_com)) deallocate(arr_chain_com)
if (allocated(arr_chain_inertia)) deallocate(arr_chain_inertia)
if (allocated(T_backup)) deallocate(T_backup)
if (allocated(K_backup)) deallocate(K_backup)
contains
!*****************************************************************************
subroutine refine_PDB(PDB, chain, refined_PDB, i_nline, n_refined)
    implicit none
    character(len=*), intent(in) :: PDB
    character(len=*), intent(in) :: chain
    integer, intent(in) :: i_nline
    integer, intent(out) :: n_refined
    character(len=80), intent(out) :: refined_PDB(i_nline)
    character(len=80) :: line
    character(len=500) :: chain_tmp
    character(len=2) :: elem
    integer :: i
    integer :: ios, OpenStatus
    refined_PDB = ''
    n_refined = 0
    open(11, file=trim(PDB), action="READ", status="old", iostat=OpenStatus)
    if (OpenStatus /= 0) stop "*** Cannot open the PDB file ***"
    chain_tmp = adjustl(chain)
    do i = 1, i_nline
        read(11,'(A)', iostat=ios) line
        if (ios < 0) exit
        if (ios > 0) stop "*** Error while reading the PDB file. ***"
        if (line(1:4) == 'ATOM') then
            if (trim(line(77:78)) == '') then
                stop "*** The atom type is not included in the PDB file. Please add the atom type information in columns 77-78 of the PDB file. ***"
            endif
            elem = adjustl(line(77:78))
            if (index(trim(chain_tmp), line(22:22)) > 0 .and. elem(1:1) /= 'H') then
                n_refined = n_refined + 1
                refined_PDB(n_refined) = line
            endif
        endif
    enddo
    close(11)
end subroutine refine_PDB

subroutine read_PDB (refined_PDB,n_refined,arr_coordi,arr_atom_type,arr_resid_type,arr_chain_label,&
    arr_chain_index,i_natom,i_nCA,ref_PDB,ref_PDB_CA)
    implicit none
    character(80), dimension(:), intent(in) :: refined_PDB
    character(80), dimension(:), intent(out) :: ref_PDB, ref_PDB_CA
    character(80) :: pdb_N, pdb_CA, pdb_C, pdb_VB
    character(8) :: x_coordi,y_coordi,z_coordi
    double precision :: coordi_N(3),coordi_CA(3),coordi_C(3),coordi_VB(3)
    integer :: atom_type_N(2),atom_type_CA(2),atom_type_C(2),atom_type_VB(2)
    character(3) :: resid_type_N,resid_type_CA,resid_type_C,resid_type_VB
    character(1) :: chain_label_N,chain_label_CA,chain_label_C,chain_label_VB
    logical :: has_N, has_CA, has_C, has_VB
    integer, intent(in) :: n_refined
    integer, intent(out) :: i_natom,i_nCA
    character(3), dimension(:), intent(out) :: arr_resid_type
    character(1), dimension(:), intent(out) :: arr_chain_label
    character(1) :: current_chain
    integer, dimension(:), intent(out) :: arr_chain_index ! if missing residues are in the chain, it regards it as the two chains
    integer :: i,OpenStatus,cnt,total_cnt,CA_cnt,onoff,resid_num,resid_num_present,cnt_dist,VB_exist,start_resid
    integer, dimension(:), allocatable :: count_chain_atom,count_chain_resid
    double precision, dimension(3) :: va_coordi,temp_coordi
    double precision, dimension(:,:), intent(out) :: arr_coordi
    double precision, dimension(ubound(arr_coordi,1),3) :: arr_coordi_pre
    integer, dimension(:,:), intent(out) :: arr_atom_type
    integer, dimension(ubound(arr_atom_type,1),2) :: arr_atom_type_pre
    double precision :: mass_N,mass_C,mass_O,mass_S,mass_total
    mass_N=14.0;mass_C=12.0;mass_O=16.0;mass_S=32.1
    cnt=0;total_cnt=0;CA_cnt=0;onoff=0
    va_coordi=0;mass_total=0.0
    read(refined_PDB(1)(23:26), '(I4)') resid_num
    current_chain=refined_PDB(1)(22:22)
    arr_coordi = 0.0d0
    arr_atom_type = 0
    arr_resid_type = ''
    arr_chain_label = ''
    arr_chain_index = 0
    arr_coordi_pre = 0.0d0
    arr_atom_type_pre = 0
    has_N = .false.; has_CA = .false.; has_C = .false.; has_VB = .false.
    coordi_N = 0.0d0; coordi_CA = 0.0d0; coordi_C = 0.0d0; coordi_VB = 0.0d0
    do i=1,n_refined
        read(refined_PDB(i)(23:26), '(I4)') resid_num_present
        if (resid_num_present .ne. resid_num .OR. refined_PDB(i)(22:22) .ne. current_chain) then
            if (has_VB) then
                va_coordi=va_coordi/mass_total
                write(x_coordi, '(F8.3)') va_coordi(1); x_coordi=adjustr(x_coordi)
                write(y_coordi, '(F8.3)') va_coordi(2); y_coordi=adjustr(y_coordi)
                write(z_coordi, '(F8.3)') va_coordi(3); z_coordi=adjustr(z_coordi)
                pdb_VB = refined_PDB(i-1)(1:12)//" VB "//refined_PDB(i-1)(17:30)//&
                                                x_coordi//y_coordi//z_coordi//refined_PDB(i-1)(55:80)
                pdb_VB(77:78) = " X"                                                
                coordi_VB=va_coordi
                atom_type_VB(1)=resid_num
                if (refined_PDB(i-1)(18:20)=="PRO") then
                    atom_type_VB(2)=44
                else
                    atom_type_VB(2)=4
                endif
                chain_label_VB=refined_PDB(i-1)(22:22)
                resid_type_VB=refined_PDB(i-1)(18:20)
                cnt=0;va_coordi=0;mass_total=0.0
            endif
            ! Always append in N / CA / C / VB order
            if (has_N) then
                total_cnt = total_cnt + 1
                arr_coordi_pre(total_cnt,:) = coordi_N
                arr_atom_type_pre(total_cnt,:) = atom_type_N
                arr_chain_label(total_cnt) = chain_label_N
                arr_resid_type(total_cnt) = resid_type_N
                call set_pdb_serial(pdb_N, total_cnt, ref_PDB(total_cnt))
            endif 
            if (has_CA) then
                total_cnt = total_cnt + 1
                CA_cnt = CA_cnt + 1
                arr_coordi_pre(total_cnt,:) = coordi_CA
                arr_atom_type_pre(total_cnt,:) = atom_type_CA
                arr_chain_label(total_cnt) = chain_label_CA
                arr_resid_type(total_cnt) = resid_type_CA
                call set_pdb_serial(pdb_CA, total_cnt, ref_PDB(total_cnt))
                call set_pdb_serial(pdb_CA, CA_cnt, ref_PDB_CA(CA_cnt))
            endif
            if (has_C) then
                total_cnt = total_cnt + 1
                arr_coordi_pre(total_cnt,:) = coordi_C
                arr_atom_type_pre(total_cnt,:) = atom_type_C
                arr_chain_label(total_cnt) = chain_label_C
                arr_resid_type(total_cnt) = resid_type_C
                call set_pdb_serial(pdb_C, total_cnt, ref_PDB(total_cnt))
            endif
            if (has_VB) then
                total_cnt = total_cnt + 1
                arr_coordi_pre(total_cnt,:) = coordi_VB
                arr_atom_type_pre(total_cnt,:) = atom_type_VB
                arr_chain_label(total_cnt) = chain_label_VB
                arr_resid_type(total_cnt) = resid_type_VB
                call set_pdb_serial(pdb_VB, total_cnt, ref_PDB(total_cnt))
            endif
            has_N = .false.; has_CA = .false.; has_C = .false.; has_VB = .false.
            coordi_N = 0.0d0; coordi_CA = 0.0d0; coordi_C = 0.0d0; coordi_VB = 0.0d0
            resid_num=resid_num_present
            current_chain=refined_PDB(i)(22:22)

            select case (refined_PDB(i)(13:16))
            case (" N  ")
                pdb_N=refined_PDB(i)
                has_N = .true.
                read(refined_PDB(i)(31:38), '(F8.3)') coordi_N(1)
                read(refined_PDB(i)(39:46), '(F8.3)') coordi_N(2)
                read(refined_PDB(i)(47:54), '(F8.3)') coordi_N(3)
                atom_type_N(1)=resid_num
                chain_label_N=refined_PDB(i)(22:22)
                resid_type_N=refined_PDB(i)(18:20)
                if (refined_PDB(i)(18:20)=="PRO") then
                    atom_type_N(2)=11
                else
                    atom_type_N(2)=1
                endif
            case (" CA ")
                pdb_CA=refined_PDB(i)
                has_CA = .true.
                read(refined_PDB(i)(31:38), '(F8.3)') coordi_CA(1)
                read(refined_PDB(i)(39:46), '(F8.3)') coordi_CA(2)
                read(refined_PDB(i)(47:54), '(F8.3)') coordi_CA(3)
                atom_type_CA(1)=resid_num
                chain_label_CA=refined_PDB(i)(22:22)
                resid_type_CA=refined_PDB(i)(18:20)
                if (refined_PDB(i)(18:20)=="PRO") then
                    atom_type_CA(2)=22
                else
                    atom_type_CA(2)=2
                endif
            case (" C  ")
                pdb_C=refined_PDB(i)
                has_C = .true.
                read(refined_PDB(i)(31:38), '(F8.3)') coordi_C(1)
                read(refined_PDB(i)(39:46), '(F8.3)') coordi_C(2)
                read(refined_PDB(i)(47:54), '(F8.3)') coordi_C(3)
                atom_type_C(1)=resid_num
                chain_label_C=refined_PDB(i)(22:22)
                resid_type_C=refined_PDB(i)(18:20)
                if (refined_PDB(i)(18:20)=="PRO") then
                    atom_type_C(2)=33
                else
                    atom_type_C(2)=3
                endif
            case (" O  ")
                ! skip O as a bead, but its mass is included in C later by make_mass_array subroutine
            case default
                read(refined_PDB(i)(31:38), '(F8.3)') temp_coordi(1)
                read(refined_PDB(i)(39:46), '(F8.3)') temp_coordi(2)
                read(refined_PDB(i)(47:54), '(F8.3)') temp_coordi(3)
                select case (adjustl(refined_PDB(i)(77:78)))
                case ("N")
                    va_coordi=va_coordi+temp_coordi*mass_N
                    mass_total=mass_total+mass_N
                    has_VB = .true.
                case ("C")
                    va_coordi=va_coordi+temp_coordi*mass_C
                    mass_total=mass_total+mass_C
                    has_VB = .true.
                case ("O")
                    va_coordi=va_coordi+temp_coordi*mass_O
                    mass_total=mass_total+mass_O
                    has_VB = .true.
                case ("S")
                    va_coordi=va_coordi+temp_coordi*mass_S
                    mass_total=mass_total+mass_S
                    has_VB = .true.
                case default
                    ! skip H or other undefined atoms
                end select
            end select
        else
            select case (refined_PDB(i)(13:16))
            case (" N  ")
                pdb_N=refined_PDB(i)
                has_N = .true.
                read(refined_PDB(i)(31:38), '(F8.3)') coordi_N(1)
                read(refined_PDB(i)(39:46), '(F8.3)') coordi_N(2)
                read(refined_PDB(i)(47:54), '(F8.3)') coordi_N(3)
                atom_type_N(1)=resid_num
                chain_label_N=refined_PDB(i)(22:22)
                resid_type_N=refined_PDB(i)(18:20)
                if (refined_PDB(i)(18:20)=="PRO") then
                    atom_type_N(2)=11
                else
                    atom_type_N(2)=1
                endif
            case (" CA ")
                pdb_CA=refined_PDB(i)
                has_CA = .true.
                read(refined_PDB(i)(31:38), '(F8.3)') coordi_CA(1)
                read(refined_PDB(i)(39:46), '(F8.3)') coordi_CA(2)
                read(refined_PDB(i)(47:54), '(F8.3)') coordi_CA(3)
                atom_type_CA(1)=resid_num
                chain_label_CA=refined_PDB(i)(22:22)
                resid_type_CA=refined_PDB(i)(18:20)
                if (refined_PDB(i)(18:20)=="PRO") then
                    atom_type_CA(2)=22
                else
                    atom_type_CA(2)=2
                endif
            case (" C  ")
                pdb_C=refined_PDB(i)
                has_C = .true.
                read(refined_PDB(i)(31:38), '(F8.3)') coordi_C(1)
                read(refined_PDB(i)(39:46), '(F8.3)') coordi_C(2)
                read(refined_PDB(i)(47:54), '(F8.3)') coordi_C(3)
                atom_type_C(1)=resid_num
                chain_label_C=refined_PDB(i)(22:22)
                resid_type_C=refined_PDB(i)(18:20)
                if (refined_PDB(i)(18:20)=="PRO") then
                    atom_type_C(2)=33
                else
                    atom_type_C(2)=3
                endif
            case (" O  ")
                ! skip O as a bead, but its mass is included in C later by make_mass_array subroutine
            case default
                read(refined_PDB(i)(31:38), '(F8.3)') temp_coordi(1)
                read(refined_PDB(i)(39:46), '(F8.3)') temp_coordi(2)
                read(refined_PDB(i)(47:54), '(F8.3)') temp_coordi(3)
                select case (adjustl(refined_PDB(i)(77:78)))
                case ("N")
                    va_coordi=va_coordi+temp_coordi*mass_N
                    mass_total=mass_total+mass_N
                    has_VB = .true.
                case ("C")
                    va_coordi=va_coordi+temp_coordi*mass_C
                    mass_total=mass_total+mass_C
                    has_VB = .true.
                case ("O")
                    va_coordi=va_coordi+temp_coordi*mass_O
                    mass_total=mass_total+mass_O
                    has_VB = .true.
                case ("S")
                    va_coordi=va_coordi+temp_coordi*mass_S
                    mass_total=mass_total+mass_S
                    has_VB = .true.
                case default
                    ! skip H or other undefined atoms
                end select
            end select
        endif
    enddo
    if (has_VB) then
        va_coordi=va_coordi/mass_total
        write(x_coordi, '(F8.3)') va_coordi(1); x_coordi=adjustr(x_coordi)
        write(y_coordi, '(F8.3)') va_coordi(2); y_coordi=adjustr(y_coordi)
        write(z_coordi, '(F8.3)') va_coordi(3); z_coordi=adjustr(z_coordi)
        pdb_VB = refined_PDB(n_refined)(1:12)//" VB "//refined_PDB(n_refined)(17:30)//&
                                        x_coordi//y_coordi//z_coordi//refined_PDB(n_refined)(55:80)
        pdb_VB(77:78) = " X"
        coordi_VB=va_coordi
        atom_type_VB(1)=resid_num
        if (refined_PDB(n_refined)(18:20)=="PRO") then
            atom_type_VB(2)=44
        else
            atom_type_VB(2)=4
        endif
        chain_label_VB=refined_PDB(n_refined)(22:22)
        resid_type_VB=refined_PDB(n_refined)(18:20)
    endif
    if (has_N) then
        total_cnt = total_cnt + 1
        arr_coordi_pre(total_cnt,:) = coordi_N
        arr_atom_type_pre(total_cnt,:) = atom_type_N
        arr_chain_label(total_cnt) = chain_label_N
        arr_resid_type(total_cnt) = resid_type_N
        call set_pdb_serial(pdb_N, total_cnt, ref_PDB(total_cnt))
    endif

    if (has_CA) then
        total_cnt = total_cnt + 1
        CA_cnt = CA_cnt + 1
        arr_coordi_pre(total_cnt,:) = coordi_CA
        arr_atom_type_pre(total_cnt,:) = atom_type_CA
        arr_chain_label(total_cnt) = chain_label_CA
        arr_resid_type(total_cnt) = resid_type_CA
        call set_pdb_serial(pdb_CA, total_cnt, ref_PDB(total_cnt))
        call set_pdb_serial(pdb_CA, CA_cnt, ref_PDB_CA(CA_cnt))
    endif

    if (has_C) then
        total_cnt = total_cnt + 1
        arr_coordi_pre(total_cnt,:) = coordi_C
        arr_atom_type_pre(total_cnt,:) = atom_type_C
        arr_chain_label(total_cnt) = chain_label_C
        arr_resid_type(total_cnt) = resid_type_C
        call set_pdb_serial(pdb_C, total_cnt, ref_PDB(total_cnt))
    endif

    if (has_VB) then
        total_cnt = total_cnt + 1
        arr_coordi_pre(total_cnt,:) = coordi_VB
        arr_atom_type_pre(total_cnt,:) = atom_type_VB
        arr_chain_label(total_cnt) = chain_label_VB
        arr_resid_type(total_cnt) = resid_type_VB
        call set_pdb_serial(pdb_VB, total_cnt, ref_PDB(total_cnt))
    endif
    ref_PDB(total_cnt+1) = "END"
    ref_PDB_CA(CA_cnt+1) = "END"

    do i=1,total_cnt
        if (i/=total_cnt) then
            if (arr_atom_type_pre(i,1)==arr_atom_type_pre(i+1,1) .AND. &
            ((arr_atom_type_pre(i,2)==3 .OR. arr_atom_type_pre(i,2)==33) .AND. (arr_atom_type_pre(i+1,2)==4&
            .OR. arr_atom_type_pre(i+1,2)==44))) then
                arr_coordi(i,:)=arr_coordi_pre(i+1,:)
                arr_coordi(i+1,:)=arr_coordi_pre(i,:)
                arr_atom_type(i,1:2)=arr_atom_type_pre(i+1,:)
                arr_atom_type(i+1,1:2)=arr_atom_type_pre(i,:)
                onoff=1
            else
                if (onoff.eq.0) then
                    arr_coordi(i,:)=arr_coordi_pre(i,:)
                    arr_atom_type(i,1:2)=arr_atom_type_pre(i,:)
                elseif (onoff.eq.1) then
                    onoff=0
                endif
            endif
        else
            if (onoff.eq.0) then
                arr_coordi(i,:)=arr_coordi_pre(i,:)
                arr_atom_type(i,1:2)=arr_atom_type_pre(i,:)
            endif
        endif
    enddo
    !!! Check there are missing parts in a protein
    cnt=1 ! first_chain starts as '1'
    do i=1,total_cnt
        if (i/=1) then
            if (arr_chain_label(i-1)/=arr_chain_label(i)) then ! different chain
                cnt=cnt+1
            else ! in the same chain
                if (arr_atom_type(i,2)==1 .OR. arr_atom_type(i,2)==11) then
                    if (.NOT.((arr_atom_type(i-1,2)==3.OR.arr_atom_type(i-1,2)==33) .AND. &
                    arr_atom_type(i-1,1)+1==arr_atom_type(i,1))) then ! NOT when i atom is N and (i-1) atom is C 
                    cnt=cnt+1
                    write(*,'(A16,A1,I5,A5,A1,I5,A9,I5)') "Missing between ", arr_chain_label(i-1), arr_atom_type(i-1,1), " and ", &
                                arr_chain_label(i-1), arr_atom_type(i,1), " ** Gap: ", arr_atom_type(i,1)-arr_atom_type(i-1,1)-1 
                    endif
                else if (arr_atom_type(i,2)==2 .OR. arr_atom_type(i,2)==22) then
                    if (.NOT.((arr_atom_type(i-1,2)==1.OR.arr_atom_type(i-1,2)==11) .AND. &
                    arr_atom_type(i-1,1)==arr_atom_type(i,1))) then ! NOT when i atom is CA and (i-1) atom is N 
                    cnt=cnt+1
                    write(*,'(A16,A1,I5,A5,A1,I5,A9,I5)') "Missing between ", arr_chain_label(i-1), arr_atom_type(i-1,1), " and ", &
                                arr_chain_label(i-1), arr_atom_type(i,1), " ** Gap: ", arr_atom_type(i,1)-arr_atom_type(i-1,1)-1
                    endif
                else if (arr_atom_type(i,2)==3 .OR. arr_atom_type(i,2)==33) then
                    if (.NOT.((arr_atom_type(i-1,2)==2 .OR. arr_atom_type(i-1,2)==22 .OR. &
                    arr_atom_type(i-1,2)==4 .OR. arr_atom_type(i-1,2)==44) .AND. &
                    arr_atom_type(i-1,1)==arr_atom_type(i,1))) then ! NOT when i atom is C and (i-1) atom is CA or VB 
                    cnt=cnt+1
                    write(*,'(A16,A1,I5,A5,A1,I5,A9,I5)') "Missing between ", arr_chain_label(i-1), arr_atom_type(i-1,1), " and ", &
                                arr_chain_label(i-1), arr_atom_type(i,1), " ** Gap: ", arr_atom_type(i,1)-arr_atom_type(i-1,1)-1
                    endif
                else if (arr_atom_type(i,2)==4 .OR. arr_atom_type(i,2)==44) then
                    if (.NOT.((arr_atom_type(i-1,2)==2 .OR. arr_atom_type(i-1,2)==22) .AND. &
                    arr_atom_type(i-1,1)==arr_atom_type(i,1))) then ! NOT when i atom is VB and (i-1) atom is CA 
                    cnt=cnt+1
                    write(*,'(A16,A1,I5,A5,A1,I5,A9,I5)') "Missing between ", arr_chain_label(i-1), arr_atom_type(i-1,1), " and ", &
                                arr_chain_label(i-1), arr_atom_type(i,1), " ** Gap: ", arr_atom_type(i,1)-arr_atom_type(i-1,1)-1
                    endif
                else if (arr_atom_type(i,2)==5) then
                    cnt=cnt+1
                    write(*,*) "Hetero atom"
                endif
            endif
        endif
        arr_chain_index(i)=cnt
    enddo
    allocate(count_chain_atom(arr_chain_index(total_cnt)))
    allocate(count_chain_resid(arr_chain_index(total_cnt)))
    cnt=0
    start_resid=0
    do i = 1, total_cnt
        if (i==1) then
            cnt=cnt+1
            start_resid = arr_atom_type(i,1)
        else if (i==total_cnt) then
            cnt=cnt+1
            count_chain_atom(arr_chain_index(i))=cnt
            count_chain_resid(arr_chain_index(i))=arr_atom_type(i,1) - start_resid + 1
        else
            if (arr_chain_index(i-1)==arr_chain_index(i)) then
                cnt=cnt+1
            else
                count_chain_atom(arr_chain_index(i-1))=cnt
                count_chain_resid(arr_chain_index(i-1))=arr_atom_type(i-1,1) - start_resid + 1
                cnt=0
                cnt=cnt+1
                start_resid=arr_atom_type(i,1)
            endif
        endif
    enddo

    !if (arr_chain_index(total_cnt)>1) then
        do i = 1,ubound(count_chain_atom,1)
            write(*,'(A9,I3,A5,I5,A9,I5,A5)') "Fragment ", i, " has ", count_chain_resid(i), " resids, ",&
                                                 count_chain_atom(i), " atoms"
        enddo
        cnt=0
        do i = 1,total_cnt
            if (arr_atom_type(i,2)==2 .OR. arr_atom_type(i,2)==22) then
                cnt = cnt + 1
            endif
        enddo
    !endif
    
    cnt_dist=0
    VB_exist=0
    do i=1,total_cnt
        if (arr_atom_type(i,2).eq.1 .OR. arr_atom_type(i,2).eq.11 .OR. arr_atom_type(i,2).eq.2 &
        .OR. arr_atom_type(i,2).eq.22 ) then
            cnt_dist=cnt_dist+1
        else if (arr_atom_type(i,2).eq.4 .OR. arr_atom_type(i,2).eq.44) then
            cnt_dist=cnt_dist+1
            VB_exist=1
        else if (arr_atom_type(i,2).eq.3 .or. arr_atom_type(i,2).eq.33) then
            if (VB_exist == 0) then 
                cnt_dist=cnt_dist+1
            else 
                VB_exist=0
            endif
        endif
        arr_atom_type(i,3)=cnt_dist ! Distance index in the chain
    enddo
    i_natom=total_cnt
    i_nCA=CA_cnt
end subroutine
!***************************************************************************************************
subroutine set_pdb_serial(line_in, serial, line_out)
    implicit none
    character(80), intent(in) :: line_in
    integer, intent(in) :: serial
    character(80), intent(out) :: line_out
    line_out = line_in
    write(line_out(7:11), '(I5)') serial
end subroutine set_pdb_serial
!***************************************************************************************************
subroutine compute_chain_rigid(arr_mass,arr_coordi,arr_chain_index,arr_chain_com,arr_chain_inertia,i_nchains)
    double precision, dimension(:), intent(in) :: arr_mass
    double precision, dimension(:,:), intent(in) :: arr_coordi
    integer, dimension(:), intent(in) :: arr_chain_index
    double precision, dimension(:,:), intent(out) :: arr_chain_com
    double precision, dimension(:,:,:), intent(out) :: arr_chain_inertia
    integer :: i_nchains,i,j,i_natom
    double precision :: chain_mass
    double precision, dimension(3) :: r
    double precision, dimension(3,3) :: EYE3
    EYE3(:,:)=0.0;EYE3(1,1)=1.0;EYE3(2,2)=1.0;EYE3(3,3)=1.0
    i_natom = ubound(arr_coordi,1)
    do i = 1,i_nchains
        arr_chain_com(i,:)=0.0
        chain_mass=0.0
        do j=1, i_natom
            if (arr_chain_index(j) == i) then
                arr_chain_com(i,:) = arr_chain_com(i,:)+arr_mass(j)*arr_coordi(j,:)
                chain_mass=chain_mass+arr_mass(j)
            endif
        enddo
        arr_chain_com(i,:) = arr_chain_com(i,:) / chain_mass
        arr_chain_inertia(:,:,i) = 0.0
        do j=1, i_natom
            if (arr_chain_index(j) == i) then
                r = arr_coordi(j,:) - arr_chain_com(i,:)
                arr_chain_inertia(:,:,i) = arr_chain_inertia(:,:,i)+arr_mass(j)*&
                (dot_product(r,r)*EYE3 - matmul(reshape(r,[3,1]),reshape(r,[1,3])))
            endif
        enddo
    enddo
end subroutine
!***************************************************************************************************
subroutine coordi_IC (arr_coordi, arr_atom_type, arr_chain_index, arr_coordi_IC, arr_atom_num)
    integer :: i_natom, i_nang
    integer :: i, cnt
    double precision, dimension(:,:), intent(in) :: arr_coordi
    integer, dimension(:,:), intent(in) :: arr_atom_type
    integer, dimension(:), intent(in) :: arr_chain_index
    double precision, dimension(3) :: diff_mat, cross_mat, temp_coordi
    double precision, dimension(:,:), intent(out) :: arr_coordi_IC
    real, dimension(:), intent(out) :: arr_atom_num
    i_natom=ubound(arr_coordi,1)
    i_nang=ubound(arr_coordi_IC,1)
    cnt=0
    do i=1,i_natom
        if (arr_atom_type(i,2)==2) then ! Check phi angle not in PROLINE
            if (i >= 3 .AND. i <= i_natom-1) then ! not in the first residue and not as the last atom
                if (arr_atom_type(i+1,2)==4 .AND. i <= i_natom-2) then ! this residue has a virtual bead
                    if (arr_chain_index(i-2)==arr_chain_index(i+2) .AND. arr_atom_type(i-2,1)+1==arr_atom_type(i+2,1) .AND. &
                        (arr_atom_type(i-2,2)==3 .OR. arr_atom_type(i-2,2)==33) .AND. arr_atom_type(i+2,2)==3 ) then ! previous residue's C and this residue's C should be in the same chain and be connected 
                        cnt=cnt+1
                        arr_atom_num(cnt)=i-0.5
                        diff_mat=arr_coordi(i,:)-arr_coordi(i-1,:)
                        diff_mat=diff_mat/norm2(diff_mat)
                        temp_coordi(:)=0d0
                        temp_coordi=arr_coordi(i,:)
                        cross_mat=cross1D(diff_mat,temp_coordi)
                        arr_coordi_IC(cnt,1:3)=diff_mat
                        arr_coordi_IC(cnt,4:6)=cross_mat
                    endif
                else if (arr_atom_type(i+1,2)==3) then ! this residue has no virtual bead
                    if (arr_chain_index(i-2)==arr_chain_index(i+1) .AND. arr_atom_type(i-2,1)+1==arr_atom_type(i+1,1) .AND. &
                        (arr_atom_type(i-2,2)==3 .OR. arr_atom_type(i-2,2)==33)) then ! previous residue's C and this residue's C should be in the same chain and be connected 
                        cnt=cnt+1
                        arr_atom_num(cnt)=i-0.5
                        diff_mat=arr_coordi(i,:)-arr_coordi(i-1,:)
                        diff_mat=diff_mat/norm2(diff_mat)
                        temp_coordi(:)=0d0
                        temp_coordi=arr_coordi(i,:)
                        cross_mat=cross1D(diff_mat,temp_coordi)
                        arr_coordi_IC(cnt,1:3)=diff_mat
                        arr_coordi_IC(cnt,4:6)=cross_mat
                    endif
                endif
            endif
        else if (arr_atom_type(i,2)==3 .OR. arr_atom_type(i,2)==33) then ! Check psi angle
            if (i >= 3 .AND. i <= i_natom-1) then 
                if ((arr_atom_type(i-1,2)==4 .OR. arr_atom_type(i-1,2)==44) .AND. i >= 4) then ! this residue has a virtual beada
                    if (arr_chain_index(i-3)==arr_chain_index(i+1) .AND. arr_atom_type(i-3,1)+1==arr_atom_type(i+1,1) .AND. & 
                    (arr_atom_type(i-3,2)==1 .OR. arr_atom_type(i-3,2)==11) .AND. (arr_atom_type(i+1,2)==1 &
                        .OR. arr_atom_type(i+1,2)==11)) then
                        cnt=cnt+1
                        arr_atom_num(cnt)=i-0.5
                        diff_mat=arr_coordi(i,:)-arr_coordi(i-2,:)
                        diff_mat=diff_mat/norm2(diff_mat)
                        temp_coordi(:)=0d0
                        temp_coordi=arr_coordi(i,:)
                        cross_mat=cross1D(diff_mat,temp_coordi)
                        arr_coordi_IC(cnt,1:3)=diff_mat
                        arr_coordi_IC(cnt,4:6)=cross_mat
                    endif
                else if (arr_atom_type(i-1,2)==2 .OR. arr_atom_type(i-1,2)==22) then
                    if (arr_chain_index(i-2)==arr_chain_index(i+1) .AND. arr_atom_type(i-2,1)+1==arr_atom_type(i+1,1) .AND. &
                    (arr_atom_type(i-2,2)==1 .OR. arr_atom_type(i-2,2)==11) .AND. (arr_atom_type(i+1,2)==1 &
                        .OR. arr_atom_type(i+1,2)==11)) then
                        cnt=cnt+1
                        arr_atom_num(cnt)=i-0.5
                        diff_mat=arr_coordi(i,:)-arr_coordi(i-1,:)
                        diff_mat=diff_mat/norm2(diff_mat)
                        temp_coordi(:)=0d0
                        temp_coordi=arr_coordi(i,:)
                        cross_mat=cross1D(diff_mat,temp_coordi)
                        arr_coordi_IC(cnt,1:3)=diff_mat
                        arr_coordi_IC(cnt,4:6)=cross_mat
                    endif
                endif
            endif
        endif
    enddo
end subroutine coordi_IC
!***************************************************************************************************
subroutine make_angle_set(arr_atom_num,ang_set)
    integer :: i,j, cnt, cnt1, i_nang
    real, dimension(:), intent(in) :: arr_atom_num
    integer, dimension(:,:), intent(out) :: ang_set
    i_nang=ubound(arr_atom_num,1)
    cnt=0
    cnt1=0
    do i=1,i_nang-1
        do j=i+1,i_nang
            cnt=cnt+1
            ang_set(cnt,1)=i
            ang_set(cnt,2)=j
        enddo
    enddo
end subroutine make_angle_set
!***************************************************************************************************
subroutine final_T_ic(arr_mass,arr_coordi,arr_coordi_IC,arr_atom_num,ang_set,PA,PC,IA,IC,MA,MC,I_all,T_ic)
    double precision, dimension(:), intent(in) :: arr_mass
    double precision, dimension(:,:), intent(in) :: arr_coordi
    double precision, dimension(:,:), intent(in) :: arr_coordi_IC
    real, dimension(:), intent(in) :: arr_atom_num
    integer, dimension(:,:), intent(in) :: ang_set
    double precision, dimension(:), intent(out) :: MA, MC
    double precision, dimension(:,:,:), intent(out) :: PA,PC,IA,IC
    double precision, dimension(:,:), intent(out) :: I_all,T_ic
    double precision, dimension(ubound(arr_coordi_IC,1),ubound(arr_coordi_IC,1)) :: T_off,T_diag
    double precision, dimension(3,3,ubound(arr_coordi,1)) :: PA_all,PC_all,IA_all,IC_all
    call make_sub_mat(arr_mass,arr_coordi,PA_all,PC_all,IA_all,IC_all)
    call make_T_off(arr_mass,arr_coordi_IC,ang_set,arr_atom_num,PA_all,PC_all,IA_all,IC_all,T_off)
    call make_T_diag(arr_mass,arr_coordi_IC,arr_atom_num,PA_all,PC_all,IA_all,IC_all,&
    PA,PC,IA,IC,MA,MC,I_all,T_diag)
    T_ic=T_off+T_diag
end subroutine final_T_ic
!***************************************************************************************************
subroutine combine_IC_and_rigid(mat_ic,mat_rigid,mat_cross,mat_combined,i_nang)
    double precision, dimension(:,:), intent(in) :: mat_ic,mat_rigid,mat_cross
    double precision, dimension(:,:), intent(out) :: mat_combined
    integer :: i_nang
    mat_combined = 0.0
    mat_combined(:i_nang,:i_nang) = mat_ic
    mat_combined(i_nang+1:,i_nang+1:) = mat_rigid
    mat_combined(:i_nang,i_nang+1:) = mat_cross
    mat_combined(i_nang+1:,:i_nang) = transpose(mat_cross)
end subroutine
!***************************************************************************************************
subroutine make_sub_mat(arr_mass,arr_coordi,PA_all,PC_all,IA_all,IC_all)
    double precision, dimension(:), intent(in) :: arr_mass
    double precision, dimension(:,:), intent(in) :: arr_coordi
    double precision, dimension(:,:,:), intent(out) :: PA_all,PC_all,IA_all,IC_all
    double precision, dimension(3,3,ubound(arr_coordi,1)) :: P_all,I_all
    double precision :: norm_val, temp_coordi(3)
    double precision, dimension(3,3) :: EYE3
    integer i, i_natom
    i_natom = ubound(arr_coordi,1)
    !$OMP parallel do schedule(static) private(i, norm_val, EYE3, temp_coordi)
    do i=1, i_natom
        temp_coordi(:)=0d0
        temp_coordi=arr_coordi(i,:)
        P_all(:,:,i)=arr_mass(i)*P_sub(temp_coordi)
        EYE3(:,:)=0d0
        norm_val=arr_coordi(i,1)**2+arr_coordi(i,2)**2+arr_coordi(i,3)**2
        EYE3(1,1)=arr_mass(i)*norm_val
        EYE3(2,2)=arr_mass(i)*norm_val
        EYE3(3,3)=arr_mass(i)*norm_val
        I_all(:,:,i)=EYE3-arr_mass(i)*index_mul(arr_coordi(i,:))
    enddo
    PA_all(:,:,1) = P_all(:,:,1)
    IA_all(:,:,1) = I_all(:,:,1)
    PC_all(:,:,i_natom) = P_all(:,:,i_natom)
    IC_all(:,:,i_natom) = I_all(:,:,i_natom)
    do i = 2, i_natom
        PA_all(:,:,i) = PA_all(:,:,i-1) + P_all(:,:,i)
        IA_all(:,:,i) = IA_all(:,:,i-1) + I_all(:,:,i)
        PC_all(:,:,i_natom-i+1) = PC_all(:,:,i_natom-i+2) + P_all(:,:,i_natom-i+1)
        IC_all(:,:,i_natom-i+1) = IC_all(:,:,i_natom-i+2) + I_all(:,:,i_natom-i+1)
    end do
end subroutine make_sub_mat
!***************************************************************************************************
subroutine make_T_off(arr_mass,arr_coordi_IC,ang_set,arr_atom_num,PA_all,PC_all,IA_all,IC_all,T_off)
    double precision, dimension(:), intent(in) :: arr_mass
    double precision, dimension(:,:), intent(in) :: arr_coordi_IC
    integer, dimension(:,:), intent(in) :: ang_set
    real, dimension(:), intent(in) :: arr_atom_num
    double precision, dimension(:,:,:), intent(in) :: PA_all,PC_all,IA_all,IC_all
    double precision, dimension(:,:), intent(out) :: T_off
    integer, dimension(ubound(ang_set,1)) :: MA_atom_num,MC_atom_num
    double precision, dimension(ubound(ang_set,1)) :: MA,MC
    double precision :: M
    integer :: i,ll,rr,i_nang_set,i_natom
    double precision, dimension(3,3,ubound(ang_set,1)) :: PA,PC,IA,IC
    double precision, dimension(3,3) :: I_mat,inv_I,temp1,temp2,temp3,temp4,EYE3
    double precision ,dimension(6,6) :: temp
    double precision, dimension(1,6) :: left_side
    double precision, dimension(6,1) :: right_side
    double precision, dimension(1,1) :: value1
    double precision, dimension(ubound(arr_mass,1)) :: prefix
    EYE3(:,:)=0d0;EYE3(1,1)=1d0;EYE3(2,2)=1d0;EYE3(3,3)=1d0
    i_nang_set = ubound(ang_set,1)
    i_natom = ubound(arr_mass,1)
    M=sum(arr_mass)
    MA_atom_num=floor(arr_atom_num(ang_set(:,1)))
    MC_atom_num=ceiling(arr_atom_num(ang_set(:,2)))
    do i=1,i_natom
        if (i==1) then
            prefix(i)=arr_mass(1)
        else
            prefix(i) = prefix(i-1) + arr_mass(i)
        endif
    enddo
    !$OMP parallel do schedule(static) private(i)
    do i=1,i_nang_set
        MA(i)=prefix(MA_atom_num(i))
        MC(i)=prefix(i_natom)-prefix(MC_atom_num(i)-1)
    enddo
    PA=PA_all(:,:,floor(arr_atom_num(ang_set(:,1))))
    PC=PC_all(:,:,ceiling(arr_atom_num(ang_set(:,2))))
    IA=IA_all(:,:,floor(arr_atom_num(ang_set(:,1))))
    IC=IC_all(:,:,ceiling(arr_atom_num(ang_set(:,2))))
    I_mat=IC_all(:,:,1)
    inv_I=inverse_3d(I_mat)
    !$OMP parallel do schedule(static) private(i, ll, rr, left_side, right_side, temp1, temp2, temp3, temp4, temp, value1)
    do i=1,i_nang_set
        ll=ang_set(i,1)
        rr=ang_set(i,2)
        left_side(1,:)=arr_coordi_IC(ll,:)
        right_side(:,1)=arr_coordi_IC(rr,:)
        temp1=matmul(matmul(IA(:,:,i),inv_I),IC(:,:,i))+matmul(transpose(PA(:,:,i)),PC(:,:,i))/dble(M)
        temp2=matmul(matmul(IA(:,:,i),inv_I),transpose(PC(:,:,i)))+transpose(PA(:,:,i))*dble(MC(i))/dble(M)
        temp3=matmul(matmul(PA(:,:,i),inv_I),IC(:,:,i))+PC(:,:,i)*dble(MA(i))/dble(M)
        temp4=matmul(matmul(PA(:,:,i),inv_I),transpose(PC(:,:,i)))+EYE3*dble(MA(i))*dble(MC(i))/dble(M)
        temp(1:3,1:3)=temp1; temp(1:3,4:6)=temp2; temp(4:6,1:3)=temp3; temp(4:6,4:6)=temp4
        value1=matmul(matmul(left_side,temp),right_side)
        T_off(ll,rr)=value1(1,1)
        T_off(rr,ll)=value1(1,1)
    enddo
end subroutine make_T_off
!***************************************************************************************************
subroutine make_T_diag(arr_mass,arr_coordi_IC,arr_atom_num,PA_all,PC_all,IA_all,IC_all,&
    PA,PC,IA,IC,MA,MC,I_all,T_diag)
    double precision, dimension(:), intent(in) :: arr_mass
    double precision, dimension(:,:), intent(in) :: arr_coordi_IC
    real, dimension(:), intent(in) :: arr_atom_num
    double precision, dimension(:,:,:), intent(in) :: PA_all,PC_all,IA_all,IC_all
    double precision, dimension(:,:), intent(out) :: T_diag
    integer, dimension(ubound(arr_coordi_IC,1)) :: MA_atom_num,MC_atom_num
    double precision, dimension(ubound(arr_coordi_IC,1)) :: MA,MC
    integer, dimension(ubound(arr_coordi_IC,1)) :: ang_list
    double precision :: M
    integer :: i, i_nang, i_nang_set, i_natom
    double precision, dimension(3,3,ubound(arr_coordi_IC,1)) :: PA,PC,IA,IC
    double precision, dimension(:,:), intent(out) :: I_all
    double precision, dimension(3,3) :: inv_I_all,temp1,temp2,temp3,temp4,EYE3
    double precision ,dimension(6,6) :: temp
    double precision, dimension(1,6) :: left_side
    double precision, dimension(6,1) :: right_side
    double precision, dimension(1,1) :: value1
    EYE3(:,:)=0d0;EYE3(1,1)=1d0;EYE3(2,2)=1d0;EYE3(3,3)=1d0
    i_nang = ubound(arr_coordi_IC, 1)
    i_nang_set = ubound(ang_list, 1)
    i_natom = ubound(arr_mass, 1)
    do i=1, i_nang
        ang_list(i)=i
    enddo

    M=sum(arr_mass)
    MA_atom_num=floor(arr_atom_num(ang_list(:)))
    MC_atom_num=ceiling(arr_atom_num(ang_list(:)))
    do i=1, i_nang_set
        MA(i)=sum(arr_mass(1:MA_atom_num(i)))
        MC(i)=sum(arr_mass(MC_atom_num(i):i_natom))
    enddo

    PA=PA_all(:,:,floor(arr_atom_num(ang_list(:))))
    PC=PC_all(:,:,ceiling(arr_atom_num(ang_list(:))))
    IA=IA_all(:,:,floor(arr_atom_num(ang_list(:))))
    IC=IC_all(:,:,ceiling(arr_atom_num(ang_list(:))))
    I_all=IC_all(:,:,1)
    inv_I_all=inverse_3d(I_all)

    !$OMP parallel do schedule(static) private(i, left_side, right_side, temp1, temp2, temp3, temp4, temp, value1)
    DO i=1,i_nang
        left_side(1,:)=arr_coordi_IC(i,:)
        right_side(:,1)=arr_coordi_IC(i,:)
        temp1=matmul(matmul(IA(:,:,i),inv_I_all),IC(:,:,i))+matmul(transpose(PA(:,:,i)),PC(:,:,i))/dble(M)
        temp2=matmul(matmul(IA(:,:,i),inv_I_all),transpose(PC(:,:,i)))+transpose(PA(:,:,i))*dble(MC(i))/dble(M)
        temp3=matmul(matmul(PA(:,:,i),inv_I_all),IC(:,:,i))+PC(:,:,i)*dble(MA(i))/dble(M)
        temp4=matmul(matmul(PA(:,:,i),inv_I_all),transpose(PC(:,:,i)))+EYE3*dble(MA(i))*dble(MC(i))/dble(M)
        temp(1:3,1:3)=temp1; temp(1:3,4:6)=temp2; temp(4:6,1:3)=temp3; temp(4:6,4:6)=temp4
        value1=matmul(matmul(left_side,temp),right_side)
        T_diag(i,i)=value1(1,1)
    enddo
end subroutine make_T_diag
!***************************************************************************************************
subroutine compute_T_rigid(arr_mass,arr_chain_index,arr_chain_inertia,i_nchains,T_rigid)
    double precision, dimension(:), intent(in) :: arr_mass
    integer, dimension(:), intent(in) :: arr_chain_index
    double precision, dimension(:,:,:), intent(in) :: arr_chain_inertia
    integer :: i_nchains
    double precision, dimension(:,:), intent(out) :: T_rigid
    double precision :: chain_mass
    integer :: i, j, i_natom
    double precision, dimension(3,3) :: EYE3
    EYE3(:,:)=0d0;EYE3(1,1)=1d0;EYE3(2,2)=1d0;EYE3(3,3)=1d0
    i_natom = ubound(arr_mass, 1)
    T_rigid(:,:) = 0d0
    do i=1,i_nchains
        chain_mass=0d0
        do j=1, i_natom
            if (arr_chain_index(j) == i) then
                chain_mass = chain_mass + arr_mass(j)
            endif   
        enddo
        T_rigid(6*(i-1)+1:6*(i-1)+3,6*(i-1)+1:6*(i-1)+3) = chain_mass*EYE3
        T_rigid(6*(i-1)+4:6*i, 6*(i-1)+4:6*i) = arr_chain_inertia(:,:,i)
    enddo
end subroutine
!***************************************************************************************************
!T_cross = off-diagonal term with IC and rigid
subroutine compute_T_cross(arr_mass,arr_coordi_IC,arr_coordi,arr_chain_index,arr_chain_com,&
    i_nchains,arr_atom_num,PA,PC,IA,IC,MA,MC,I_all,T_cross)
    double precision, dimension(:), intent(in) :: arr_mass
    integer, dimension(:), intent(in) :: arr_chain_index
    real, dimension(:), intent(in) :: arr_atom_num
    double precision, dimension(:,:), intent(in) :: arr_coordi_IC,arr_coordi,arr_chain_com
    double precision, dimension(:,:,:), intent(in) :: PA,PC,IA,IC
    double precision, dimension(:), intent(in) :: MA,MC
    double precision, dimension(3,3), intent(in) :: I_all
    double precision, dimension(3,3) :: inv_I_all
    double precision, dimension(3,3) :: temp1_A,temp2_A,temp3_A,temp4_A
    double precision, dimension(3,3) :: temp1_C,temp2_C,temp3_C,temp4_C
    double precision :: left_side(1,6),temp_A(6,6),temp_C(6,6),temp_coordi(3)
    double precision :: term_A(1,6),term_C(1,6),right_side(6,3),rot_change(3,3),r(3)
    double precision :: trans_A(6,3),trans_C(6,3),rot_A(6,3),rot_C(6,3)
    double precision, dimension(:,:), intent(out) :: T_cross
    double precision :: M
    integer :: i, j, k, i_nchains, i_nang, i_natom
    double precision, dimension(3,3) :: EYE3
    EYE3(:,:) = 0d0; EYE3(1,1) = 1d0; EYE3(2,2) = 1d0; EYE3(3,3) = 1d0
    i_nang = ubound(arr_coordi_IC, 1); i_natom = ubound(arr_mass, 1)
    M=sum(arr_mass)
    inv_I_all=inverse_3d(I_all)
    !$OMP parallel do schedule(static) private(i, j, k, left_side, right_side, r, rot_change, &
    !$OMP& temp1_A, temp2_A, temp3_A, temp4_A, temp_A, term_A, &
    !$OMP& temp1_C, temp2_C, temp3_C, temp4_C, temp_C, term_C, &
    !$OMP& trans_A, trans_C, rot_A, rot_C, temp_coordi)
    do i = 1, i_nang ! 1 to i_nang
        temp1_A(:,:)=0d0;temp2_A(:,:)=0d0;temp3_A(:,:)=0d0;temp4_A(:,:)=0d0;temp_A(:,:)=0d0;term_A(:,:)=0d0
        temp1_C(:,:)=0d0;temp2_C(:,:)=0d0;temp3_C(:,:)=0d0;temp4_C(:,:)=0d0;temp_C(:,:)=0d0;term_C(:,:)=0d0
        left_side(1,:)=arr_coordi_IC(i,:)
        temp1_A=matmul(I_all-IA(:,:,i),inv_I_all)
        temp2_A=-transpose(PA(:,:,i))/dble(M)
        temp3_A=-matmul(PA(:,:,i),inv_I_all)
        temp4_A=(dble(M)-dble(MA(i)))/dble(M)*EYE3
        temp_A(1:3,1:3)=temp1_A;temp_A(1:3,4:6)=temp2_A;temp_A(4:6,1:3)=temp3_A;temp_A(4:6,4:6)=temp4_A
        term_A=matmul(left_side,temp_A)
        temp1_C=matmul(I_all-IC(:,:,i),inv_I_all)
        temp2_C=-transpose(PC(:,:,i))/dble(M)
        temp3_C=-matmul(PC(:,:,i),inv_I_all)
        temp4_C=(dble(M)-dble(MC(i)))/dble(M)*EYE3
        temp_C(1:3,1:3)=temp1_C;temp_C(1:3,4:6)=temp2_C;temp_C(4:6,1:3)=temp3_C;temp_C(4:6,4:6)=temp4_C
        term_C=-1d0 * matmul(left_side,temp_C)
        do j = 1, i_nchains
            trans_A(:,:) = 0d0; trans_C(:,:) = 0d0; rot_A(:,:) = 0d0; rot_C(:,:) = 0d0
            do k = 1, i_natom ! 1 to i_natom
                if (arr_chain_index(k) == j) then
                    temp_coordi(:)=0d0
                    temp_coordi=arr_coordi(k,:)
                    right_side(1:3,1:3) = -1d0 * (P_sub(temp_coordi))
                    right_side(4:6,1:3) = EYE3
                    r = arr_coordi(k,:) - arr_chain_com(j,:)
                    rot_change = -1d0 * (P_sub(r))
                    if (k <= floor(arr_atom_num(i))) then ! when the atom is included in block A and chain j                        
                        trans_A = trans_A + arr_mass(k) * right_side
                        rot_A = rot_A + arr_mass(k) * matmul(right_side,rot_change)
                    elseif (k >= ceiling(arr_atom_num(i))) then ! when the atom is included in block C and chain j
                        trans_C = trans_C + arr_mass(k) * right_side
                        rot_C = rot_C + arr_mass(k) * matmul(right_side,rot_change)
                    endif
                endif
            enddo
            ! T_cross = size(i_nang, 6 * i_nchains)
            T_cross(i,6*(j-1)+1:6*(j-1)+3) = reshape(matmul(term_A,trans_A) + matmul(term_C,trans_C),[3])
            T_cross(i,6*(j-1)+4:6*(j-1)+6) = reshape(matmul(term_A,rot_A) + matmul(term_C,rot_C), [3])
        enddo
    enddo
end subroutine
!***************************************************************************************************
function cross(a,b)
    integer i
    double precision, dimension(:,:), intent(in) :: a, b
    double precision, dimension(1:ubound(a,1),1:ubound(a,2)) :: cross
    do i=1,ubound(cross,1)
        cross(i,1)=a(i,2)*b(i,3)-a(i,3)*b(i,2)
        cross(i,2)=a(i,3)*b(i,1)-a(i,1)*b(i,3)
        cross(i,3)=a(i,1)*b(i,2)-a(i,2)*b(i,1)
    enddo
end function
!***************************************************************************************************
function skew_sym_mat(a)
    integer i
    double precision, dimension(3), intent(in) :: a
    double precision, dimension(3,3) :: skew_sym_mat
    skew_sym_mat=0.0
    do i=1,3
        skew_sym_mat(1,2)=-a(3);skew_sym_mat(1,3)=a(2)
        skew_sym_mat(2,1)=a(3);skew_sym_mat(2,3)=-a(1)
        skew_sym_mat(3,1)=-a(2);skew_sym_mat(3,2)=a(1)
    enddo
end function
!***************************************************************************************************
function cross1D(a,b)
    double precision, dimension(3), intent(in) :: a, b
    double precision, dimension(3) :: cross1D
    cross1D(1)=a(2)*b(3)-a(3)*b(2)
    cross1D(2)=a(3)*b(1)-a(1)*b(3)
    cross1D(3)=a(1)*b(2)-a(2)*b(1)
end function
!***************************************************************************************************
function index_mul(a)
    double precision, dimension(:), intent(in) :: a
    double precision, dimension(1:ubound(a,1),1:ubound(a,1)) :: index_mul
    integer :: i,j
    do i=1,ubound(a,1)
        do j=1,ubound(a,1)
            index_mul(i,j)=a(i)*a(j)
        enddo
    enddo
end function
!***************************************************************************************************
function P_sub(a)
    double precision, dimension(3), intent(in) :: a
    double precision, dimension(3,3) :: P_sub
    P_sub(1,1)=0d0
    P_sub(1,2)=-a(3)
    P_sub(1,3)=a(2)
    P_sub(2,1)=a(3)
    P_sub(2,2)=0d0
    P_sub(2,3)=-a(1)
    P_sub(3,1)=-a(2)
    P_sub(3,2)=a(1)
    P_sub(3,3)=0d0
end function
!***************************************************************************************************
function inverse_3d(a)
    double precision, dimension(3,3), intent(in) :: a
    double precision, dimension(3,3):: inverse_3d
    double precision, dimension(3,3) :: adj_mat
    double precision, dimension(6,6) :: temp_mat
    double precision :: det_val
    integer :: i,j
    det_val=0
    temp_mat(1:3,1:3)=a
    temp_mat(1:3,4:6)=a
    temp_mat(4:6,1:3)=a
    temp_mat(4:6,4:6)=a
    do i=1,3
        det_val=det_val+temp_mat(1,i)*(temp_mat(2,i+1)*temp_mat(3,i+2)-temp_mat(2,i+2)*temp_mat(3,i+1))
    enddo
    do i=1,3
        do j=1,3
        adj_mat(i,j)=(temp_mat(i+1,j+1)*temp_mat(i+2,j+2)-temp_mat(i+1,j+2)*temp_mat(i+2,j+1))
        enddo
    enddo
    adj_mat=transpose(adj_mat)
    inverse_3d=adj_mat/det_val
end function
!***************************************************************************************************
subroutine make_k_mat (arr_coordi,k_mat,i_nkpair,d_cart_cutoff,&
    C_seq,C_cart,P_seq,P_cart,arr_atom_type,arr_chain_index)
integer :: i_natom, i, j
double precision, intent(in) :: d_cart_cutoff,C_seq,C_cart,P_seq,P_cart
double precision :: chain_index
integer, intent(out) :: i_nkpair
integer, dimension(:,:), intent(in) :: arr_atom_type
integer, dimension(:), intent(in) :: arr_chain_index
double precision :: d_seq_dist
!double precision :: proportion
double precision, dimension(:,:), intent(in) :: arr_coordi
double precision, dimension(:,:), intent(out) :: k_mat
i_natom=ubound(arr_coordi,1)
i_nkpair=0
k_mat=0.0      
do i=1,i_natom-1
    do j=i+1,i_natom
        if( norm2(arr_coordi(i,1:3)-arr_coordi(j,1:3))<=d_cart_cutoff) then
            if (arr_chain_index(i)/=arr_chain_index(j)) then
                chain_index=0.0 ! two atoms are in different chains
                d_seq_dist=1000.0 !
            else
                chain_index=1.0 ! two atoms are within a chain
                if (arr_atom_type(i,2)==4 .OR. arr_atom_type(i,2)==44) then ! if first atom is VB
                    if (arr_atom_type(i,3)==arr_atom_type(j,3)) then ! it means one is VB and the other is C in the same amino acid
                        d_seq_dist=2.0
                    else
                        d_seq_dist=dble(arr_atom_type(j,3)-arr_atom_type(i,3)+2)
                    endif
                else ! if first atom is not VB
                    d_seq_dist=dble(arr_atom_type(j,3)-arr_atom_type(i,3))
                endif
            endif
            i_nkpair=i_nkpair+1
            k_mat(i,j)=chain_index*C_seq*(d_seq_dist**P_seq)+&
                    C_cart*(norm2(arr_coordi(i,1:3)-arr_coordi(j,1:3))**P_cart)
            k_mat(j,i)=k_mat(i,j)
        endif
    enddo
enddo
!write(*,*) "the number of all k pairs: ", i_nkpair
end subroutine make_k_mat
!***************************************************************************************************
subroutine make_R_ic(arr_coordi,k_mat,R_ic)
    integer :: i_natom, i, j, block_index, block_index2
    double precision, dimension(:,:), intent(in) :: arr_coordi
    double precision, dimension(:,:), intent(in) :: k_mat
    double precision, dimension(:,:), intent(out) :: R_ic
    i_natom = ubound(arr_coordi,1)
    R_ic(:,:)=0d0
    do i = 1,i_natom
        do j = i_natom,i,-1
            block_index=i*(2*(i_natom+1)+1-i)/2+j
            block_index2=(i-1)*(2*(i_natom+1)+2-i)/2+j
            R_ic(:,6*(block_index-1)+1:6*(block_index))= &
            R_ic(:,6*(block_index)+1:6*(block_index+1))+ &
            R_ic(:,6*(block_index2-1)+1:6*(block_index2))- &
            R_ic(:,6*(block_index2)+1:6*(block_index2+1))+ & 
            make_S(arr_coordi,k_mat,i,j)
        enddo
    enddo
end subroutine make_R_ic
!***************************************************************************************************
function make_S(arr_coordi,k_mat,i,j) result(S_mat)
    double precision, dimension(:,:), intent(in) :: arr_coordi
    double precision, dimension(:,:), intent(in) :: k_mat
    integer :: i,j
    double precision, dimension(6,6) :: S_mat
    double precision, dimension(3) :: rij,cross_ij,temp_coordi1,temp_coordi2
    double precision, dimension(6) :: temp
    if (i.LT.1 .OR. j.GT.ubound(arr_coordi,1)) then 
        S_mat(:,:)=0d0
    else
        if (k_mat(i,j) == 0d0) then
            S_mat(:,:)=0d0
        else
            rij = arr_coordi(i,:)-arr_coordi(j,:)
            temp_coordi1(:)=0d0; temp_coordi2(:)=0d0
            temp_coordi1=arr_coordi(i,:)
            temp_coordi2=arr_coordi(j,:)
            cross_ij = cross1D(temp_coordi1,temp_coordi2)
            temp(1:3) = cross_ij; temp(4:6)=rij
            S_mat=(k_mat(i,j)/norm2(rij)**2)*index_mul(temp)
        endif
    endif
end function make_S
!***************************************************************************************************
subroutine make_K_off(i_natom,arr_coordi_IC,arr_atom_num,ang_set,R_ic,K_off)
    double precision, dimension(:,:), intent(in) :: arr_coordi_IC ! i_nang X 6
    real, dimension(:), intent(in) :: arr_atom_num ! i_nang
    integer, dimension(:,:), intent(in) :: ang_set 
    integer :: i_natom, i_nang_set, block_index
    double precision, dimension(1,6) :: left_side
    double precision, dimension(6,1) :: right_side
    double precision :: ls(6),rs(6),M(6,6),tmp(6)
    double precision, dimension(:,:), intent(in) :: R_ic ! 6*(i_natom+1) x 6*(i_natom+1)
    double precision, dimension(:), intent(out) :: K_off ! (i_nang_set_reduced)
    double precision, dimension(1,1) :: value1
    integer i,lower,upper
    i_nang_set = ubound(ang_set, 1)
    K_off(:) = 0d0
    !$OMP parallel do schedule(static) private(i, ls, rs, M, lower, upper, block_index)
    do i = 1, i_nang_set
        ls = arr_coordi_IC(ang_set(i,1),:)
        rs = arr_coordi_IC(ang_set(i,2),:)
        lower=floor(arr_atom_num(ang_set(i,1)))
        upper=ceiling(arr_atom_num(ang_set(i,2)))
        block_index=lower*(2*(i_natom+1)+1-lower)/2+upper
        M=R_ic(:,6*(block_index-1)+1:6*block_index)
        K_off(i)= sum(ls(:) * matmul(M,rs))
    enddo
end subroutine
!***************************************************************************************************
subroutine make_K_diag(i_natom,arr_coordi_IC,arr_atom_num,R_ic,K_diag)
    double precision, dimension(:,:), intent(in) :: arr_coordi_IC ! i_nang X 6
    real, dimension(:), intent(in) :: arr_atom_num ! i_nang
    double precision, dimension(1,6) :: left_side
    double precision, dimension(6,1) :: right_side
    double precision, dimension(:,:), intent(in) :: R_ic ! 6*(i_natom+1) x 6*(i_natom+1)
    double precision, dimension(:), intent(out) :: K_diag ! (i_nang)
    double precision, dimension(1,1) :: value1       
    integer :: i,lower,upper,i_natom,i_nang,block_index
    i_nang = ubound(arr_coordi_IC, 1)
    K_diag(:) = 0d0
    !$OMP PARALLEL DO private(i,left_side,right_side,lower,upper,value1,block_index)
    do i = 1, i_nang
        left_side(1,:) = arr_coordi_IC(i,:)
        right_side(:,1) = arr_coordi_IC(i,:)
        lower=floor(arr_atom_num(i))
        upper=ceiling(arr_atom_num(i))
        block_index=lower*(2*(i_natom+1)+1-lower)/2+upper
        value1=matmul(matmul(left_side,R_ic(:,6*(block_index-1)+1:6*block_index))&
        ,right_side)
        K_diag(i)=value1(1,1)
    enddo
end subroutine
!***************************************************************************************************
subroutine make_K_rigid(arr_coordi,arr_chain_index,i_nchains,&
                            k_mat,arr_chain_com,K_rigid)
    double precision, dimension(:,:), intent(in) :: arr_coordi
    double precision, dimension(:,:), intent(in) :: k_mat
    double precision, dimension(:,:), intent(in) :: arr_chain_com
    integer, dimension(:), intent(in) :: arr_chain_index
    double precision, dimension(:,:), intent(out) :: K_rigid
    integer :: i, j, ci, cj, i_natom, i_nchains, idx, i_npairs, pair
    integer, allocatable :: chain_pairlist(:,:)
    double precision :: dist
    double precision, dimension(3,3) :: Y,Ji,Jj
    double precision, dimension(6,6) :: block_A,block_B,block_C
    double precision, dimension(6,6) :: temp_A,temp_B,temp_C
    double precision, dimension(3) :: r
    i_natom = ubound(arr_coordi, 1)
    K_rigid(:,:) = 0d0
    i_npairs = i_nchains * (i_nchains -1) / 2
    allocate( chain_pairlist(i_npairs, 2))
    idx = 0
    do ci = 1, i_nchains
        do cj = ci + 1, i_nchains
            idx = idx + 1
            chain_pairlist(idx,1) = ci
            chain_pairlist(idx,2) = cj
        end do
    end do

    do pair = 1, i_npairs
        ci = chain_pairlist(pair, 1)
        cj = chain_pairlist(pair, 2)
        block_A = 0.0;block_B = 0.0;block_C = 0.0
        do i = 1, i_natom
            if (arr_chain_index(i) == ci) then
                do j = 1, i_natom
                    if (arr_chain_index(j) == cj) then
                        if (k_mat(i,j) > 0.0) then
                            r = arr_coordi(i,:) - arr_coordi(j,:)
                            dist = norm2(r)
                            Y = matmul(reshape(r,[3,1]),reshape(r,[1,3]))/dist**2
                            Ji = skew_sym_mat(arr_coordi(i,:)-arr_chain_com(ci,:))
                            Jj = skew_sym_mat(arr_coordi(j,:)-arr_chain_com(cj,:))
                            temp_A=0.0;temp_B=0.0;temp_C=0.0
                            temp_A(1:3,1:3)=Y
                            temp_A(1:3,4:6)= transpose(matmul(Ji,Y))
                            temp_A(4:6,1:3)= matmul(Ji,Y)
                            temp_A(4:6,4:6)=-matmul(matmul(Ji,Y),Ji)
                            temp_B(1:3,1:3)=-transpose(Y)
                            temp_B(1:3,4:6)=-transpose(matmul(Jj,Y))
                            temp_B(4:6,1:3)= transpose(matmul(Y,Ji))
                            temp_B(4:6,4:6)= transpose(matmul(matmul(Jj,Y),Ji))
                            temp_C(1:3,1:3)=Y
                            temp_C(1:3,4:6)= transpose(matmul(Jj,Y))
                            temp_C(4:6,1:3)= matmul(Jj,Y)
                            temp_C(4:6,4:6)=-matmul(matmul(Jj,Y),Jj)
                            block_A = block_A + k_mat(i,j)*temp_A 
                            block_B = block_B + k_mat(i,j)*temp_B 
                            block_C = block_C + k_mat(i,j)*temp_C   
                        endif
                    endif
                enddo
            endif
        enddo
        K_rigid(6*(ci-1)+1:6*(ci-1)+6,6*(ci-1)+1:6*(ci-1)+6)=&
        K_rigid(6*(ci-1)+1:6*(ci-1)+6,6*(ci-1)+1:6*(ci-1)+6)+block_A
        K_rigid(6*(cj-1)+1:6*(cj-1)+6,6*(cj-1)+1:6*(cj-1)+6)=&
        K_rigid(6*(cj-1)+1:6*(cj-1)+6,6*(cj-1)+1:6*(cj-1)+6)+block_C
        K_rigid(6*(ci-1)+1:6*(ci-1)+6,6*(cj-1)+1:6*(cj-1)+6)=block_B
        K_rigid(6*(cj-1)+1:6*(cj-1)+6,6*(ci-1)+1:6*(ci-1)+6)=transpose(block_B)
    enddo
end subroutine
!***************************************************************************************************
subroutine make_K_cross(arr_coordi,arr_chain_index,i_nchains,arr_coordi_IC,arr_atom_num,&
    k_mat,arr_chain_com,K_cross)
implicit none
double precision, dimension(:,:), intent(in) :: arr_coordi
double precision, dimension(:,:), intent(in) :: k_mat
double precision, dimension(:,:), intent(in) :: arr_coordi_IC ! i_nang X 6
real, dimension(:), intent(in) :: arr_atom_num ! i_nang
double precision, dimension(:,:), intent(in) :: arr_chain_com
integer, dimension(:), intent(in) :: arr_chain_index
integer :: i_nchains
double precision, dimension(:,:), intent(out) :: K_cross
double precision, allocatable, dimension(:,:) :: hess_sum_list,hess_rot_sum_list
integer :: i,j,ci,ang,upper,lower,i_nang,i_natom,i_chain_idx
double precision :: dist,left(1,6),kij,alpha,c11,c12,c13,c21,c22,c23,c31,c32,c33
double precision :: H11,H12,H13,H21,H22,H23,H31,H32,H33
double precision :: H_cart(3 * ubound(arr_coordi, 1), 3 * ubound(arr_coordi, 1))
double precision, dimension(3,3) :: EYE3,Hii,temp_hess,temp_hess_rot
double precision, dimension(6,3) :: angle_change,temp_trans,temp_rot,S,SR,blkT,blkR
double precision, dimension(3,3) :: rigid_change,Hij
double precision, dimension(1,3) :: trans_term, rot_term
double precision, dimension(3) :: r,r_chain,temp_coordi
double precision, allocatable :: angle_change_table(:,:,:),pref_sum(:,:),pref_rot(:,:)
integer, allocatable :: upper_of(:),lower_of(:)
double precision :: l1,l2,l3,l4,l5,l6,t1,t2,t3,r1,r2,r3
EYE3(:,:) = 0d0; EYE3(1,1) = 1d0; EYE3(2,2) = 1d0; EYE3(3,3) = 1d0
i_nang = ubound(arr_coordi_IC, 1);i_natom = ubound(arr_coordi, 1)
!!! Cartesian coordinate HESSIAN
H_cart(:,:) = 0d0
!$omp parallel do schedule(guided,16) default(none)                                &
!$omp& private(i,j,kij,r,dist,alpha,c11,c12,c13,c21,c22,c23,c31,c32,c33)           &
!$omp& shared(H_cart,arr_coordi,k_mat,i_natom)
do i = 1, i_natom
    do j = i + 1, i_natom
        if (k_mat(i,j) > 1.0d-8) then
            kij = k_mat(i,j)
            r = arr_coordi(i,:) - arr_coordi(j,:)
            dist = norm2(r)
            alpha = -kij/(dist*dist)
            c11 = alpha*r(1)*r(1); c12 = alpha*r(1)*r(2); c13 = alpha*r(1)*r(3)
            c21 = alpha*r(2)*r(1); c22 = alpha*r(2)*r(2); c23 = alpha*r(2)*r(3)
            c31 = alpha*r(3)*r(1); c32 = alpha*r(3)*r(2); c33 = alpha*r(3)*r(3)
            H_cart(3*i-2,3*j-2) = c11; H_cart(3*i-2,3*j-1) = c12; H_cart(3*i-2,3*j) = c13
            H_cart(3*i-1,3*j-2) = c21; H_cart(3*i-1,3*j-1) = c22; H_cart(3*i-1,3*j) = c23
            H_cart(3*i  ,3*j-2) = c31; H_cart(3*i  ,3*j-1) = c32; H_cart(3*i  ,3*j) = c33
        endif
    enddo
enddo
H_cart = H_cart + transpose(H_cart)

!$omp parallel do schedule(guided,16) default(none)                               &
!$omp& private(i,j,H11,H12,H13,H21,H22,H23,H31,H32,H33)                           &
!$omp& shared(H_cart,i_natom)
do i = 1, i_natom
    H11=0d0; H12=0d0; H13=0d0
    H21=0d0; H22=0d0; H23=0d0
    H31=0d0; H32=0d0; H33=0d0
    do j = 1, i_natom
        if (j == i) cycle
        H11 = H11 - H_cart(3*i-2,3*j-2);  H12 = H12 - H_cart(3*i-2,3*j-1);  H13 = H13 - H_cart(3*i-2,3*j)
        H21 = H21 - H_cart(3*i-1,3*j-2);  H22 = H22 - H_cart(3*i-1,3*j-1);  H23 = H23 - H_cart(3*i-1,3*j)
        H31 = H31 - H_cart(3*i  ,3*j-2);  H32 = H32 - H_cart(3*i  ,3*j-1);  H33 = H33 - H_cart(3*i  ,3*j)
    end do
    H_cart(3*i-2,3*i-2) = H11;  H_cart(3*i-2,3*i-1) = H12;  H_cart(3*i-2,3*i) = H13
    H_cart(3*i-1,3*i-2) = H21;  H_cart(3*i-1,3*i-1) = H22;  H_cart(3*i-1,3*i) = H23
    H_cart(3*i  ,3*i-2) = H31;  H_cart(3*i  ,3*i-1) = H32;  H_cart(3*i  ,3*i) = H33
end do

allocate(angle_change_table(6,3,i_natom))
!$omp parallel do default(none) private(i) schedule(static) &
!$omp& shared(angle_change_table,EYE3,arr_coordi,i_natom)
do i = 1, i_natom
    angle_change_table(1:3,1:3,i) = P_sub(arr_coordi(i,:))  ! 3x3 skew
    angle_change_table(4:6,1:3,i) = -EYE3                   ! -I
end do

allocate( hess_sum_list(6 * i_natom, 3 * i_nchains)); allocate( hess_rot_sum_list(6 * i_natom, 3 * i_nchains))
hess_sum_list(:,:) = 0d0
hess_rot_sum_list(:,:) = 0d0
!$omp parallel do schedule(guided,16) default(none)                                      &
!$omp& private(i, j, i_chain_idx, r_chain, rigid_change, Hij, S, SR)                     &
!$omp& shared(hess_sum_list, hess_rot_sum_list, H_cart, arr_coordi, arr_chain_com,       &
!$omp&        arr_chain_index, i_natom, i_nchains, angle_change_table, k_mat)
do i = 1, i_natom
    do j = 1, i_natom
        if (k_mat(i,j) > 1.0d-8 .OR. i==j) then
            i_chain_idx = arr_chain_index(j)
            r_chain = arr_coordi(j,:) - arr_chain_com(i_chain_idx,:)
            Hij(1,1) = H_cart(3*i-2,3*j-2); Hij(1,2) = H_cart(3*i-2,3*j-1); Hij(1,3) = H_cart(3*i-2,3*j)
            Hij(2,1) = H_cart(3*i-1,3*j-2); Hij(2,2) = H_cart(3*i-1,3*j-1); Hij(2,3) = H_cart(3*i-1,3*j)
            Hij(3,1) = H_cart(3*i  ,3*j-2); Hij(3,2) = H_cart(3*i  ,3*j-1); Hij(3,3) = H_cart(3*i  ,3*j)
            S = matmul( angle_change_table(:,:,i), Hij )
            rigid_change = -1.0 * (P_sub(r_chain))

            hess_sum_list(6*(i-1)+1:6*i, 3*(i_chain_idx-1)+1:3*i_chain_idx) =  &
                hess_sum_list(6*(i-1)+1:6*i, 3*(i_chain_idx-1)+1:3*i_chain_idx) + S
            
            SR = matmul( S, rigid_change )
            hess_rot_sum_list(6*(i-1)+1:6*i, 3*(i_chain_idx-1)+1:3*i_chain_idx) =  &
                 hess_rot_sum_list(6*(i-1)+1:6*i, 3*(i_chain_idx-1)+1:3*i_chain_idx) + SR
        endif
    enddo
enddo

allocate(pref_sum(6*i_natom, 3*i_nchains))
allocate(pref_rot(6*i_natom, 3*i_nchains))
pref_sum(1:6, :) = hess_sum_list(1:6, :)
pref_rot(1:6, :) = hess_rot_sum_list(1:6, :)
do i = 2, i_natom
    pref_sum(6*(i-1)+1:6*i, :) = pref_sum(6*(i-2)+1:6*(i-1), :) + hess_sum_list(6*(i-1)+1:6*i, :)
    pref_rot(6*(i-1)+1:6*i, :) = pref_rot(6*(i-2)+1:6*(i-1), :) + hess_rot_sum_list(6*(i-1)+1:6*i, :)
end do
allocate(upper_of(i_nang), lower_of(i_nang))
do ang = 1, i_nang
  upper_of(ang) = ceiling(arr_atom_num(ang))
  lower_of(ang) = floor(arr_atom_num(ang))
end do
K_cross(:,:) = 0d0
!$omp parallel do collapse(2) schedule(guided,16) default(none)                           &
!$omp& private(ci, ang, upper, lower, blkT, blkR, l1,l2,l3,l4,l5,l6, t1,t2,t3, r1,r2,r3) &
!$omp& shared(K_cross, pref_sum, pref_rot, arr_coordi_IC, upper_of, lower_of, i_natom, i_nang, i_nchains)
do ci = 1, i_nchains
    do ang = 1, i_nang
        upper = upper_of(ang)
        lower = lower_of(ang)
        
        if (upper <= lower) then ! choose domain C
            if (upper > 1) then
                blkT(:,:) = pref_sum(6*(i_natom-1)+1:6*i_natom, 3*(ci-1)+1:3*ci) - &
                    pref_sum(6*(upper-2)+1:6*(upper-1), 3*(ci-1)+1:3*ci)
                blkR(:,:) = pref_rot(6*(i_natom-1)+1:6*i_natom, 3*(ci-1)+1:3*ci) - &
                    pref_rot(6*(upper-2)+1:6*(upper-1), 3*(ci-1)+1:3*ci)
            else
                blkT(:,:) = pref_sum(6*(i_natom-1)+1:6*i_natom, 3*(ci-1)+1:3*ci)
                blkR(:,:) = pref_rot(6*(i_natom-1)+1:6*i_natom, 3*(ci-1)+1:3*ci)
            endif
        else
            if (lower >= 1) then
                blkT(:,:) =  pref_sum(6*(lower-1)+1:6*lower,  3*(ci-1)+1:3*ci)
                blkR(:,:) =  pref_rot(6*(lower-1)+1:6*lower,  3*(ci-1)+1:3*ci)
            else
                blkT(:,:) = 0d0
                blkR(:,:) = 0d0
            end if
            blkT(:,:) = -blkT(:,:)
            blkR(:,:) = -blkR(:,:)
        endif
        ! left = 1x6 : arr_coordi_IC(ang,1:6)
        l1 = arr_coordi_IC(ang,1); l2 = arr_coordi_IC(ang,2); l3 = arr_coordi_IC(ang,3)
        l4 = arr_coordi_IC(ang,4); l5 = arr_coordi_IC(ang,5); l6 = arr_coordi_IC(ang,6)
        ! trans_term = left * blkT  (1x6 · 6x3)  → (1x3)
        t1 = l1*blkT(1,1) + l2*blkT(2,1) + l3*blkT(3,1) + l4*blkT(4,1) + l5*blkT(5,1) + l6*blkT(6,1)
        t2 = l1*blkT(1,2) + l2*blkT(2,2) + l3*blkT(3,2) + l4*blkT(4,2) + l5*blkT(5,2) + l6*blkT(6,2)
        t3 = l1*blkT(1,3) + l2*blkT(2,3) + l3*blkT(3,3) + l4*blkT(4,3) + l5*blkT(5,3) + l6*blkT(6,3)
        ! rot_term = left * blkR  (1x6 · 6x3)  → (1x3)
        r1 = l1*blkR(1,1) + l2*blkR(2,1) + l3*blkR(3,1) + l4*blkR(4,1) + l5*blkR(5,1) + l6*blkR(6,1)
        r2 = l1*blkR(1,2) + l2*blkR(2,2) + l3*blkR(3,2) + l4*blkR(4,2) + l5*blkR(5,2) + l6*blkR(6,2)
        r3 = l1*blkR(1,3) + l2*blkR(2,3) + l3*blkR(3,3) + l4*blkR(4,3) + l5*blkR(5,3) + l6*blkR(6,3)
        ! Write K_cross block.
        K_cross(ang,6*(ci-1)+1) = t1
        K_cross(ang,6*(ci-1)+2) = t2
        K_cross(ang,6*(ci-1)+3) = t3
        K_cross(ang,6*(ci-1)+4) = r1
        K_cross(ang,6*(ci-1)+5) = r2
        K_cross(ang,6*(ci-1)+6) = r3
    enddo
enddo
end subroutine
!***************************************************************************************************
Subroutine eig_func(K_ic_mat,T_ic_mat,i_dof,Q_order,D_order)
    character JOBZ, UPLO
    integer :: ITYPE,LDA,LDB,N,LWORK,INFO,i_dof
    double precision, dimension(:,:), intent(in) :: T_ic_mat
    double precision, dimension(:,:), intent(out) :: K_ic_mat,Q_order
    double precision, dimension(i_dof,i_dof) :: Q
    double precision, dimension(:) :: D_order
    double precision, dimension(i_dof) :: W
    double precision, dimension(i_dof*10) :: WORK
    ITYPE=1
    N=i_dof
    LDA=N; LDB=N
    LWORK=i_dof*10
    JOBZ = 'V'
    UPLO = 'U'
    !example of DSYGV: https://github.com/numericalalgorithmsgroup/LAPACK_Examples/blob/master/examples/source/dsygv_example.f90
    call DSYGV(ITYPE,JOBZ,UPLO,N,K_ic_mat,LDA,T_ic_mat,LDB,W,WORK,LWORK,INFO)
    !do i=1,i_nang
    !    if (W(i) .lt. 0) STOP "*** Minus Eigenvalue error !!! ***"     
    !enddo
    Q=K_ic_mat
    Q_order=Q
    D_order=W
    !call sort_ascending(D,Q,i_nang,D_order,Q_order)
end subroutine
!***************************************************************************************************
Subroutine eig_func_part_dsygvx(K_ic_mat,T_ic_mat,i_dof,Q,D,lowest_num)
    character JOBZ, UPLO, range
    integer :: ITYPE,LDA,LDB,LDZ,N,INFO,i_dof,il,iu,lowest_num,i,j,LWORK
    double precision, dimension(:,:) :: T_ic_mat
    double precision, dimension(:,:) :: K_ic_mat
    double precision, dimension(i_dof) :: D
    double precision, dimension(:,:) :: Q
    double precision, dimension(:), allocatable :: WORK
    integer, dimension(i_dof*5) :: IWORK
    integer, dimension(i_dof) :: IFAIL
    double precision :: vl, vu, abstol, WORKQ(1)
    double precision, external :: dlamch
    INFO=0; IFAIL(:)=0; IWORK(:)=0; ITYPE=1; range='I'
    N=i_dof; LDA=N; LDB=N
    LDZ=i_dof; JOBZ = 'V' ;UPLO = 'U'
    vl=0.0; vu=0.0; abstol = 2.0d0 * dlamch('S'); il=1; iu=lowest_num
    !example of DSYGV: https://github.com/numericalalgorithmsgroup/LAPACK_Examples/blob/master/examples/source/dsygv_example.f90
    LWORK = -1
    call dsygvx(ITYPE,JOBZ,range,UPLO,N,K_ic_mat,LDA,T_ic_mat,LDB,vl,vu,il,iu,&
                abstol,lowest_num,D,Q,LDZ,WORKQ,LWORK,IWORK,IFAIL,INFO)
    if (INFO .ne. 0) return

    LWORK = int(WORKQ(1))
    allocate(WORK(LWORK))
    call dsygvx(ITYPE,JOBZ,range,UPLO,N,K_ic_mat,LDA,T_ic_mat,LDB,vl,vu,il,iu,&
                abstol,lowest_num,D,Q,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
end subroutine
!***************************************************************************************************
subroutine ICtoCC(arr_mass,Q,arr_coordi,arr_coordi_IC,arr_atom_num,PA_mat &
    ,PC_mat,IA_mat,IC_mat,MA_mat,MC_mat,Q_cc,i_nmodes)
    double precision, dimension(:), intent(in) :: arr_mass
    double precision, dimension(:,:) :: arr_coordi, arr_coordi_IC
    real, dimension(:) :: arr_atom_num
    double precision, dimension(:,:,:) :: PA_mat,PC_mat,IA_mat,IC_mat
    double precision, dimension(:) :: MA_mat,MC_mat
    double precision, dimension(:,:) :: Q
    double precision, dimension(:,:), allocatable :: Q_local_A,Q_local_C
    double precision, dimension(:,:), intent(out) :: Q_cc
    !double precision, dimension(:,:), intent(out) :: Q_cc_weighted
    double precision, dimension(3,6) :: left_mat=0
    double precision, dimension(3,3) :: I_mat,I_inv,temp11,temp12,temp13,temp14,EYE3
    double precision, dimension(3,3) :: temp21,temp22,temp23,temp24
    double precision, dimension(6,1) :: data_ic_A,data_ic_C
    double precision, dimension(6,6) :: temp1,temp2 ! temp = [temp1, temp2]
    integer :: i_natom,i_nang,i,s,i_nmodes          !        [temp3, temp4]
    double precision :: M
    double precision, dimension(6,i_nmodes) :: local_sum
    double precision, dimension(:,:), allocatable :: local_C,local_A
    double precision, dimension(:,:,:), allocatable :: mat_domain_A, mat_domain_C
    integer, dimension(1) :: min_pos
    double precision, dimension(1) :: min_val
    i_natom=ubound(arr_coordi,1)
    i_nang=ubound(arr_coordi_IC,1)
    M=sum(arr_mass)
    EYE3=0
    EYE3(1,1)=1; EYE3(2,2)=1; EYE3(3,3)=1
    allocate( Q_local_A(1,i_nmodes));allocate( Q_local_C(1,i_nmodes))
    I_mat=IA_mat(:,:,1)+IC_mat(:,:,1)
    I_inv=inverse_3d(I_mat)
    allocate ( local_C(6,i_nmodes)); allocate ( local_A(6,i_nmodes))
    allocate ( mat_domain_A(6,i_nmodes,i_nang))
    allocate ( mat_domain_C(6,i_nmodes,i_nang))
    do i=1,i_nang
        local_C = 0
        local_A = 0
        !! When an atom in domain C
        Q_local_C(1,:)=Q(i,:)
        data_ic_C(:,1)=arr_coordi_IC(i,:)
        temp11=matmul(I_inv,IA_mat(:,:,i))
        temp12=matmul(I_inv,transpose(PA_mat(:,:,i)))
        temp13=PA_mat(:,:,i)/M
        temp14=EYE3*MA_mat(i)/M
        temp1(1:3,1:3)=temp11
        temp1(1:3,4:6)=temp12
        temp1(4:6,1:3)=temp13
        temp1(4:6,4:6)=temp14
        local_C = matmul(matmul(temp1,data_ic_C),Q_local_C)
        if(i.NE.1) then
            mat_domain_C(:,:,i)=-local_C(:,:)+mat_domain_C(:,:,i-1)
        else
            mat_domain_C(:,:,i)=-local_C(:,:)
        endif
        !! When an atom in domain A
        s=i_nang-i+1
        Q_local_A(1,:)=Q(s,:)
        data_ic_A(:,1)=arr_coordi_IC(s,:)
        temp21=matmul(I_inv,IC_mat(:,:,s))
        temp22=matmul(I_inv,transpose(PC_mat(:,:,s)))
        temp23=PC_mat(:,:,s)/M
        temp24=EYE3*MC_mat(s)/M
        temp2(1:3,1:3)=temp21
        temp2(1:3,4:6)=temp22
        temp2(4:6,1:3)=temp23
        temp2(4:6,4:6)=temp24
        local_A = matmul(matmul(temp2,data_ic_A),Q_local_A)
        if(i.NE.1) then
            mat_domain_A(:,:,s)=local_A(:,:)+mat_domain_A(:,:,s+1)
        else
            mat_domain_A(:,:,s)=local_A(:,:)
        endif
    enddo
    
    do i=1,i_natom
        local_sum(:,:) = 0d0
        ! make [P E] matrix
        left_mat(1,2)=-arr_coordi(i,3); left_mat(1,3)=arr_coordi(i,2)
        left_mat(2,1)=arr_coordi(i,3); left_mat(2,3)=-arr_coordi(i,1)
        left_mat(3,1)=-arr_coordi(i,2); left_mat(3,2)=arr_coordi(i,1)
        left_mat(1,4)=1; left_mat(2,5)=1; left_mat(3,6)=1
        !******************************************************************
        ! find the nearest dihedral angle to atom i 
        min_pos=minloc(abs(arr_atom_num-i))
        min_val=arr_atom_num(min_pos)-i
        if (min_val(1).GT.0) then ! the first ~ (min_pos-1)th angles make i be in domain C // (min_pos)th ~ nang th angles make i be in domain A
            if (min_pos(1).NE.1) then
                !Q_cc(3*i-2:3*i,:)=matmul(left_mat,mat_domain_A(:,:,min_pos)+mat_domain_C(:,:,min_pos-1))
                local_sum = reshape(mat_domain_A(:,:,min_pos) + mat_domain_C(:,:,min_pos-1), [6,i_nmodes])
            else
                !Q_cc(3*i-2:3*i,:)=matmul(left_mat,mat_domain_A(:,:,min_pos))
                local_sum = reshape(mat_domain_A(:,:,min_pos), [6,i_nmodes])
            endif
        else ! the first ~ (min_pos)th angles make i be in domain C // (min_pos+1)th ~ nang th angles make i be in domain A
            if (min_pos(1).NE.i_nang) then
                local_sum = reshape(mat_domain_A(:,:,min_pos+1) + mat_domain_C(:,:,min_pos), [6,i_nmodes])
            else
                local_sum = reshape(mat_domain_C(:,:,min_pos), [6,i_nmodes])
            endif
        endif
        Q_cc(3*i-2:3*i,:)=matmul(left_mat,local_sum(:,:))
    enddo
    
    !!$OMP parallel do schedule(static) private(i)
    !do i=1,i_natom
    !    Q_cc_weighted(3*i-2:3*i,:) = sqrt(arr_mass(i)) * Q_cc(3*i-2:3*i,:)
    !enddo
end subroutine
!***************************************************************************************************
subroutine RIGIDtoCC(Q_rigid,Q_cc_rigid,arr_chain_index,arr_chain_com,arr_coordi)
    double precision, dimension(:,:), intent(in) :: Q_rigid
    double precision, dimension(:,:), intent(out) :: Q_cc_rigid
    integer, dimension(:), intent(in) :: arr_chain_index
    double precision, dimension(:,:), intent(in) :: arr_chain_com,arr_coordi
    integer :: i,nmode,c
    double precision :: r(3), t(3), rot(3), d(3), J(3,3)
    Q_cc_rigid = 0.0d0
    !$OMP parallel do collapse(2) schedule(static) private(i, nmode, c, r, t, rot, J, d)
    do nmode = 1, size(Q_rigid,2)
        do i = 1, size(arr_coordi,1)
            c = arr_chain_index(i)
            r = arr_coordi(i,:) - arr_chain_com(c,:)    ! The relative position from the chain's COM position
            t = Q_rigid(6*(c-1)+1:6*(c-1)+3, nmode)     ! The translation vector
            rot = Q_rigid(6*(c-1)+4:6*(c-1)+6, nmode)     ! The translation vector

            ! the displacement vector d = t + (w x r) >> J = skew-symmetric matrix for expressing "w x { }"  
            J(1,1) = 0.0d0; J(1,2) = -rot(3); J(1,3) = rot(2)
            J(2,1) = rot(3); J(2,2) = 0.0d0; J(2,3) = -rot(1)
            J(3,1) = -rot(2); J(3,2) = rot(1); J(3,3) = 0.0d0
            d = t + matmul(J,r)
            
            Q_cc_rigid(3*i-2:3*i, nmode) = d
        enddo
    enddo
    !!$OMP parallel do schedule(static) private(i)
    !do i = 1, size(arr_coordi,1)
    !    Q_cc_rigid_weighted(3*i-2:3*i, :) = sqrt(arr_mass(i)) * Q_cc_rigid(3*i-2:3*i, :)
    !enddo
end subroutine
!***************************************************************************************************
subroutine CCtoCA(Q_cc,Q_CA,arr_atom_type)
    double precision, dimension(:,:), intent(in) :: Q_cc
    double precision, dimension(:,:), intent(out) :: Q_CA
    integer, dimension(:,:), intent(in) :: arr_atom_type
    integer :: i, cnt, i_nline
    cnt=0; i_nline = ubound(arr_atom_type, 1)
    do i = 1, i_nline
        if (arr_atom_type(i,2).eq.2 .OR. arr_atom_type(i,2).eq.22) then
            cnt=cnt+1
            Q_CA(cnt*3-2:cnt*3,:)=Q_cc(i*3-2:i*3,:)
        endif
    enddo
end subroutine
!***************************************************************************************************
!!!! Paper Vicky Choi (2006) J.Chem.Inf.Model introduced the simple rotation method
subroutine gen_single_chain_path_frame(arr_eig_vec_IC,arr_ref_coordi,arr_mass,&
                                        arr_atom_num,arr_atom_type,d_weight,&
                                        movie_mode_start,movie_mode_end,arr_output_deformed_coordi_CC)
    double precision, dimension(:,:), intent(in) :: arr_eig_vec_IC,arr_ref_coordi
    double precision, dimension(:), intent(in) :: arr_mass
    real, dimension(:), intent(in) :: arr_atom_num
    integer, dimension(:,:), intent(in) :: arr_atom_type
    integer :: i_nframe,movie_mode_start,movie_mode_end
    double precision :: d_weight,norm_minus,norm_plus
    double precision, dimension(:,:,:,:), intent(out) :: arr_output_deformed_coordi_CC
    double precision, dimension(ubound(arr_ref_coordi,1),3) :: arr_deformed_coordi_minus
    double precision, dimension(ubound(arr_ref_coordi,1),3) :: arr_deformed_coordi_plus
    double precision, dimension(ubound(arr_ref_coordi,1),3) :: arr_coordi_diff_plus
    double precision, dimension(ubound(arr_ref_coordi,1),3) :: arr_coordi_diff_minus
    integer :: i,j,k,s,mode_num,out_mode,cnt,i_nang,i_natom,mid
    integer, dimension(ubound(arr_eig_vec_IC,1)) :: lower,upper,axis_lower,axis_upper
    integer, dimension(ubound(arr_eig_vec_IC,1),2) :: fragment_range
    double precision, dimension(:,:), allocatable :: unit_theta
    double precision, dimension(3) :: Q_ref,D,D_rotm,D_rotp,n_ref
    double precision, dimension(3) :: qvm,qvp,qvm_acc_temp
    double precision, dimension(3) :: tm,tp,qvp_acc_temp,tm_prev,tp_prev
    double precision, dimension(3) :: Qm_world_rot,Qp_world_rot,Qm_rot,Qp_rot
    double precision :: thetap,thetam,M
    double precision, dimension(3) :: center,Qm_world,Qp_world
    double precision, dimension(3) :: qvp_prev,qvm_prev,nm_world,np_world,tm_prev_rot,tp_prev_rot
    double precision, dimension(3,3) :: rot_mat
    double precision, dimension(ubound(arr_eig_vec_IC,1),3) :: tm_acc,tp_acc
    double precision, dimension(ubound(arr_eig_vec_IC,1),3) :: qvm_acc,qvp_acc
    double precision, dimension(ubound(arr_eig_vec_IC,1)) :: q0m_acc,q0p_acc
    double precision :: q0m,q0p,q0m_acc_temp,q0p_acc_temp,q0m_prev,q0p_prev
    allocate(unit_theta(ubound(arr_eig_vec_IC,1),ubound(arr_eig_vec_IC,2)))
    out_mode = 0
    i_nframe=ubound(arr_output_deformed_coordi_CC,3)
    unit_theta=d_weight*arr_eig_vec_IC
    i_nang=ubound(arr_eig_vec_IC,1)
    i_natom=ubound(arr_ref_coordi,1)    
    mid=(i_nframe-1)/2+1
    M=sum(arr_mass(:))
    do i=1,i_nang
        upper(i)=ceiling(arr_atom_num(i))
        lower(i)=floor(arr_atom_num(i))
        axis_upper(i)=upper(i)
        axis_lower(i)=lower(i)
        if ((arr_atom_type(upper(i),2)==3 .OR. arr_atom_type(upper(i),2)==33) .AND. &
            (arr_atom_type(upper(i)-1,2)==4 .OR. arr_atom_type(upper(i)-1,2)==44)) then
            axis_lower(i)=upper(i)-2
        endif
    enddo
    do i=1,i_nang
        if (i/=i_nang) then
            fragment_range(i,1)=ceiling(arr_atom_num(i))+1
            fragment_range(i,2)=ceiling(arr_atom_num(i+1))
        else
            fragment_range(i,1)=ceiling(arr_atom_num(i))+1
            fragment_range(i,2)=i_natom
        endif
    enddo
    
    do mode_num=movie_mode_start,movie_mode_end
        out_mode = out_mode + 1
        write(*,*) "mode num: ", mode_num
        arr_output_deformed_coordi_CC(:,:,mid,out_mode)=arr_ref_coordi(:,:)
        do i=1,(i_nframe-1)/2
            arr_deformed_coordi_minus(:,:)=arr_ref_coordi(:,:)
            arr_deformed_coordi_plus(:,:)=arr_ref_coordi(:,:)
            ! accumulated q0, qv, and t
            q0m_acc=0.0d0;q0p_acc=0.0d0;
            qvm_acc=0.0d0;qvp_acc=0.0d0;
            tm_acc=0.0d0;tp_acc=0.0d0
            do j=1,i_nang
                arr_coordi_diff_plus=0.0d0
                arr_coordi_diff_minus=0.0d0
                n_ref(:)=arr_ref_coordi(axis_upper(j),:)-arr_ref_coordi(axis_lower(j),:)
                n_ref=n_ref/norm2(n_ref)
                Q_ref=arr_ref_coordi(axis_upper(j),:)
                thetam=unit_theta(j,mode_num)*(-i)
                thetap=unit_theta(j,mode_num)*i

                 if (j == 1) then
                    q0m_prev = 1.0d0; qvm_prev = 0.0d0
                    q0p_prev = 1.0d0; qvp_prev = 0.0d0
                    tm_prev = 0.0d0;  tp_prev = 0.0d0
                else
                    q0m_prev = q0m_acc(j-1); qvm_prev = qvm_acc(j-1,:)
                    q0p_prev = q0p_acc(j-1); qvp_prev = qvp_acc(j-1,:)
                    tm_prev = tm_acc(j-1,:); tp_prev = tp_acc(j-1,:)
                end if

                call cal_qt(q0m_prev, qvm_prev, n_ref, nm_world)
                call cal_qt(q0p_prev, qvp_prev, n_ref, np_world)
                call cal_qt(q0m_prev, qvm_prev, Q_ref, Qm_rot)
                call cal_qt(q0p_prev, qvp_prev, Q_ref, Qp_rot)

                Qm_world = Qm_rot + tm_prev
                Qp_world = Qp_rot + tp_prev

                q0m=cos(0.5d0*thetam)
                qvm=sin(0.5d0*thetam)*nm_world   
                q0p=cos(0.5d0*thetap)
                qvp=sin(0.5d0*thetap)*np_world
                
                call quat_mul(q0m, qvm, q0m_prev, qvm_prev, q0m_acc_temp, qvm_acc_temp)
                call quat_mul(q0p, qvp, q0p_prev, qvp_prev, q0p_acc_temp, qvp_acc_temp)
                
                norm_minus = sqrt(q0m_acc_temp**2 + dot_product(qvm_acc_temp,qvm_acc_temp))
                norm_plus = sqrt(q0p_acc_temp**2 + dot_product(qvp_acc_temp,qvp_acc_temp))
                if (abs(norm_minus - 1.0d0) > 1.0d-12) then
                    q0m_acc_temp = q0m_acc_temp / norm_minus
                    qvm_acc_temp = qvm_acc_temp / norm_minus
                endif
                if (abs(norm_plus - 1.0d0) > 1.0d-12) then
                    q0p_acc_temp = q0p_acc_temp / norm_plus
                    qvp_acc_temp = qvp_acc_temp / norm_plus
                endif
                
                call cal_qt(q0m,qvm,Qm_world,Qm_world_rot)
                call cal_qt(q0p,qvp,Qp_world,Qp_world_rot)
                tm = Qm_world-Qm_world_rot
                tp = Qp_world-Qp_world_rot

                call cal_qt(q0m, qvm, tm_prev, tm_prev_rot)
                call cal_qt(q0p, qvp, tp_prev, tp_prev_rot)
                tm_acc(j,:) = tm + tm_prev_rot
                tp_acc(j,:) = tp + tp_prev_rot                    

                q0m_acc(j)=q0m_acc_temp                    
                q0p_acc(j)=q0p_acc_temp 
                qvm_acc(j,:)=qvm_acc_temp
                qvp_acc(j,:)=qvp_acc_temp
            enddo

            do j=1,i_nang
                if (fragment_range(j,2) < fragment_range(j,1)) cycle
                
                do k=fragment_range(j,1),fragment_range(j,2)
                    D=arr_ref_coordi(k,:)
                    call cal_qt(q0m_acc(j),qvm_acc(j,:),D,D_rotm)
                    call cal_qt(q0p_acc(j),qvp_acc(j,:),D,D_rotp)
                    arr_deformed_coordi_minus(k,1:3)=D_rotm+tm_acc(j,:)
                    arr_deformed_coordi_plus(k,1:3)=D_rotp+tp_acc(j,:)
                enddo
            enddo

            !!!!centering and rotating deformed coordi_BB
            !!!!minus direction
            center=0.0
            do s=1,ubound(arr_deformed_coordi_minus,1)
                center=center+arr_mass(s)*arr_deformed_coordi_minus(s,:)
            enddo
            center=center/M
            do s=1,ubound(arr_deformed_coordi_minus,1)
                arr_deformed_coordi_minus(s,:)=arr_deformed_coordi_minus(s,:)-center
            enddo
            rot_mat=0.0!;temp1=0.0;temp2=0.0
            call make_a_hat_kabsch(arr_ref_coordi,arr_deformed_coordi_minus,rot_mat) ! Kabasch method is used ***************
            arr_deformed_coordi_minus=transpose(matmul&
                    (rot_mat,transpose(arr_deformed_coordi_minus)))
            !!!!plus direction
            center=0.0
            do s=1,ubound(arr_deformed_coordi_plus,1)
                center=center+arr_mass(s)*arr_deformed_coordi_plus(s,:)
            enddo
            center=center/M
            do s=1,ubound(arr_deformed_coordi_plus,1)
                arr_deformed_coordi_plus(s,:)=arr_deformed_coordi_plus(s,:)-center
            enddo
            rot_mat=0.0!;temp1=0.0;temp2=0.0
            call make_a_hat_kabsch(arr_ref_coordi,arr_deformed_coordi_plus,rot_mat) ! Kabasch method is used ***************
            arr_deformed_coordi_plus=transpose(matmul&
                    (rot_mat,transpose(arr_deformed_coordi_plus)))
            arr_output_deformed_coordi_CC(:,:,mid - i,out_mode)=arr_deformed_coordi_minus(:,:)
            arr_output_deformed_coordi_CC(:,:,mid + i,out_mode)=arr_deformed_coordi_plus(:,:)
        enddo
    enddo
end subroutine

subroutine cal_qt(q0,qv,Q,Q_rot) ! [q0, qv, t]
    double precision, dimension(3), intent(in) :: Q,qv
    double precision, dimension(3), intent(out) :: Q_rot
    double precision, intent(in) :: q0
    double precision, dimension(3) :: a,b
    a = cross1D(qv,Q) ! a = (qv ​× Q)
    b = cross1D(qv,a) ! b = qv × (qv x Q)
    ! Q' = Q + 2q0 ​(qv ​× Q) + 2qv​ × (qv ​× Q)
    Q_rot = Q + 2*q0*a + 2 * b
    ! t = Q - Q_rot
end subroutine

subroutine quat_mul(q0, qv, q0p, qvp, q0acc, qvacc) ! q0,qv = current q0,qv / q0p,qvp = parent's q0,qv, / q0acc, qvacc = accumuated ones
    double precision, intent(in)  :: q0, qv(3), q0p, qvp(3)
    double precision, intent(out) :: q0acc, qvacc(3)
    ! (qa * qb)
    q0acc = q0*q0p - dot_product(qv, qvp)
    qvacc = q0*qvp + q0p*qv + cross1D(qv, qvp)
end subroutine
!***************************************************************************************************
!!!! Paper Vicky Choi (2006) J.Chem.Inf.Model introduced the simple rotation method
subroutine gen_multi_chain_path_frame(arr_eig_vec_IC,arr_eig_vec_rigid,arr_ref_coordi,arr_mass,&
                                        arr_chain_index,arr_atom_num,arr_atom_type,&
                                        d_weight,movie_mode_start,movie_mode_end,arr_output_deformed_coordi_CC)
    double precision, dimension(:,:), intent(in) :: arr_eig_vec_IC,arr_ref_coordi
    double precision, dimension(:,:), intent(in) :: arr_eig_vec_rigid
    double precision, dimension(:), intent(in) :: arr_mass
    integer, dimension(:), intent(in) :: arr_chain_index
    real, dimension(:), intent(in) :: arr_atom_num
    integer, dimension(:,:), intent(in) :: arr_atom_type
    integer, dimension(:), allocatable :: target_chain_start,target_chain_end
    integer :: i_nframe
    double precision :: d_weight,norm_minus,norm_plus
    double precision, dimension(:,:,:,:), intent(out) :: arr_output_deformed_coordi_CC
    double precision, dimension(ubound(arr_ref_coordi,1),3) :: arr_deformed_coordi_minus
    double precision, dimension(ubound(arr_ref_coordi,1),3) :: arr_deformed_coordi_plus
    double precision, dimension(ubound(arr_ref_coordi,1),3) :: arr_coordi_diff_plus
    double precision, dimension(ubound(arr_ref_coordi,1),3) :: arr_coordi_diff_minus
    integer :: i,j,k,s,mode_num,cnt,i_nang,i_natom,mid
    integer :: c_start,c_end,c_leng,movie_mode_start,movie_mode_end,out_mode
    integer, dimension(ubound(arr_eig_vec_IC,1)) :: lower,upper,axis_lower,axis_upper
    integer, dimension(ubound(arr_eig_vec_IC,1),2) :: fragment_range
    double precision, dimension(3,3) :: rodp_rigid,rodm_rigid
    double precision, dimension(3,3) :: jacop_rigid,jacom_rigid
    double precision, dimension(:,:), allocatable :: unit_theta
    double precision, dimension(:,:), allocatable :: arr_chain_comp,arr_chain_comm
    double precision, dimension(:,:), allocatable :: unit_rigid,rigid_trans_vec
    double precision, dimension(3) :: Q_ref,D,D_rotm,D_rotp,n_vec,n_ref
    double precision, dimension(3) :: qvm,qvp,qvm_acc_temp
    double precision, dimension(3) :: tm,tp,qvp_acc_temp,tm_prev,tp_prev
    double precision, dimension(3) :: Qm_world_rot,Qp_world_rot,Qm_rot,Qp_rot
    double precision :: thetap,thetam,M
    double precision, dimension(3) :: center,Qm_world,Qp_world,u
    double precision, dimension(3) :: tm_rigid,tp_rigid,n_rigid
    double precision, dimension(3) :: qvp_prev,qvm_prev,nm_world,np_world,tm_prev_rot,tp_prev_rot
    double precision, dimension(3,3) :: rot_mat
    double precision, dimension(ubound(arr_eig_vec_IC,1),3) :: tm_acc,tp_acc
    double precision, dimension(ubound(arr_eig_vec_IC,1),3) :: qvm_acc,qvp_acc
    double precision, dimension(ubound(arr_eig_vec_IC,1)) :: q0m_acc,q0p_acc
    double precision :: q0m,q0p,q0m_acc_temp,q0p_acc_temp,q0m_prev,q0p_prev
    double precision :: theta_rigid_minus,theta_rigid_plus
    allocate(unit_theta(ubound(arr_eig_vec_IC,1),ubound(arr_eig_vec_IC,2)))
    allocate(unit_rigid(ubound(arr_eig_vec_rigid,1),ubound(arr_eig_vec_rigid,2)))
    out_mode = 0
    i_nframe = ubound(arr_output_deformed_coordi_CC,3)
    unit_theta=d_weight*arr_eig_vec_IC
    unit_rigid=d_weight*arr_eig_vec_rigid
    i_nang=ubound(arr_eig_vec_IC,1)
    i_natom=ubound(arr_ref_coordi,1)    
    mid=(i_nframe-1)/2+1
    M=sum(arr_mass(:))
    do i=1,i_nang
        upper(i)=ceiling(arr_atom_num(i))
        lower(i)=floor(arr_atom_num(i))
        axis_upper(i)=upper(i)
        axis_lower(i)=lower(i)
        if ((arr_atom_type(upper(i),2)==3 .OR. arr_atom_type(upper(i),2)==33) .AND. &
            (arr_atom_type(upper(i)-1,2)==4 .OR. arr_atom_type(upper(i)-1,2)==44)) then
            axis_lower(i)=upper(i)-2
        endif
    enddo
    do i=1,i_nang
        if (i/=i_nang) then
            fragment_range(i,1)=ceiling(arr_atom_num(i))+1
            fragment_range(i,2)=ceiling(arr_atom_num(i+1))
        else
            fragment_range(i,1)=ceiling(arr_atom_num(i))+1
            fragment_range(i,2)=i_natom
        endif
    enddo
    i_nchains=arr_chain_index(i_natom)
    allocate(arr_chain_comp(i_nchains,3))
    allocate(arr_chain_comm(i_nchains,3))
    allocate(target_chain_start(i_nchains));allocate(target_chain_end(i_nchains))
    allocate(rigid_trans_vec(i_nchains,3))
    do i=1,i_nchains
        do s=1,i_natom
            if (arr_chain_index(s)==i) then
                target_chain_end(i)=s
            endif
            if (arr_chain_index(i_natom-s+1)==i) then
                target_chain_start(i)=i_natom-s+1
            endif
        enddo
    enddo
    do mode_num=movie_mode_start,movie_mode_end
        out_mode = out_mode + 1
        write(*,*) "mode num: ", mode_num
        arr_output_deformed_coordi_CC(:,:,mid,out_mode)=arr_ref_coordi(:,:)
        do i=1,(i_nframe-1)/2
            arr_deformed_coordi_minus(:,:)=arr_ref_coordi(:,:)
            arr_deformed_coordi_plus(:,:)=arr_ref_coordi(:,:)
            ! accumulated q0, qv, and t
            q0m_acc=0.0d0;q0p_acc=0.0d0;
            qvm_acc=0.0d0;qvp_acc=0.0d0;
            tm_acc=0.0d0;tp_acc=0.0d0
            do j=1,i_nang
                arr_coordi_diff_plus=0.0d0
                arr_coordi_diff_minus=0.0d0
                n_ref(:)=arr_ref_coordi(axis_upper(j),:)-arr_ref_coordi(axis_lower(j),:)
                n_ref=n_ref/norm2(n_ref)
                n_vec(:)=arr_ref_coordi(axis_upper(j),:)-arr_ref_coordi(axis_lower(j),:)
                n_vec=n_vec/norm2(n_vec)
                Q_ref=arr_ref_coordi(axis_upper(j),:)
                thetam=unit_theta(j,mode_num)*(-i)
                thetap=unit_theta(j,mode_num)*i
                 if (j == 1) then
                    q0m_prev = 1.0d0; qvm_prev = 0.0d0
                    q0p_prev = 1.0d0; qvp_prev = 0.0d0
                    tm_prev = 0.0d0; tp_prev = 0.0d0
                else
                    q0m_prev = q0m_acc(j-1); qvm_prev = qvm_acc(j-1,:)
                    q0p_prev = q0p_acc(j-1);  qvp_prev = qvp_acc(j-1,:)
                    tm_prev = tm_acc(j-1,:); tp_prev = tp_acc(j-1,:)
                end if

                call cal_qt(q0m_prev, qvm_prev, n_ref, nm_world)
                call cal_qt(q0p_prev, qvp_prev, n_ref, np_world)
                call cal_qt(q0m_prev, qvm_prev, Q_ref, Qm_rot)
                call cal_qt(q0p_prev, qvp_prev, Q_ref, Qp_rot)

                Qm_world = Qm_rot + tm_prev
                Qp_world = Qp_rot + tp_prev

                q0m=cos(0.5d0*thetam)
                qvm=sin(0.5d0*thetam)*nm_world   
                q0p=cos(0.5d0*thetap)
                qvp=sin(0.5d0*thetap)*np_world
                
                call quat_mul(q0m, qvm, q0m_prev, qvm_prev, q0m_acc_temp, qvm_acc_temp)
                call quat_mul(q0p, qvp, q0p_prev, qvp_prev, q0p_acc_temp, qvp_acc_temp)

                norm_minus = sqrt(q0m_acc_temp**2 + dot_product(qvm_acc_temp,qvm_acc_temp))
                norm_plus = sqrt(q0p_acc_temp**2 + dot_product(qvp_acc_temp,qvp_acc_temp))
                if (abs(norm_minus - 1.0d0) > 1.0d-12) then
                    q0m_acc_temp = q0m_acc_temp / norm_minus
                    qvm_acc_temp = qvm_acc_temp / norm_minus
                endif
                if (abs(norm_plus - 1.0d0) > 1.0d-12) then
                    q0p_acc_temp = q0p_acc_temp / norm_plus
                    qvp_acc_temp = qvp_acc_temp / norm_plus
                endif

                call cal_qt(q0m,qvm,Qm_world,Qm_world_rot)
                call cal_qt(q0p,qvp,Qp_world,Qp_world_rot)
                tm = Qm_world-Qm_world_rot
                tp = Qp_world-Qp_world_rot

                call cal_qt(q0m, qvm, tm_prev, tm_prev_rot)
                call cal_qt(q0p, qvp, tp_prev, tp_prev_rot)
                tm_acc(j,:) = tm + tm_prev_rot
                tp_acc(j,:) = tp + tp_prev_rot                    

                q0m_acc(j)=q0m_acc_temp                    
                q0p_acc(j)=q0p_acc_temp 
                qvm_acc(j,:)=qvm_acc_temp
                qvp_acc(j,:)=qvp_acc_temp
            enddo
            do j=1,i_nang
                if (fragment_range(j,2) < fragment_range(j,1)) cycle
                
                do k=fragment_range(j,1),fragment_range(j,2)
                    D=arr_ref_coordi(k,:)
                    call cal_qt(q0m_acc(j),qvm_acc(j,:),D,D_rotm)
                    call cal_qt(q0p_acc(j),qvp_acc(j,:),D,D_rotp)
                    arr_deformed_coordi_minus(k,1:3)=D_rotm+tm_acc(j,:)
                    arr_deformed_coordi_plus(k,1:3)=D_rotp+tp_acc(j,:)
                enddo
            enddo

            call compute_chain_com(arr_mass,arr_deformed_coordi_minus,arr_chain_index,arr_chain_comm,i_nchains)
            call compute_chain_com(arr_mass,arr_deformed_coordi_plus,arr_chain_index,arr_chain_comp,i_nchains)
            do j=1,i_nchains
                c_start= target_chain_start(j)
                c_end  = target_chain_end(j)
                c_leng = c_end - c_start + 1
                tm_rigid=unit_rigid(6*(j-1)+1:6*(j-1)+3,mode_num)*(-dble(i))
                tp_rigid=unit_rigid(6*(j-1)+1:6*(j-1)+3,mode_num)*(dble(i))
                n_rigid=unit_rigid(6*(j-1)+4:6*j,mode_num)
                if (norm2(n_rigid) > 1.0d-12) then
                    u = n_rigid/norm2(n_rigid)
                    theta_rigid_minus=norm2(n_rigid)*(-dble(i))
                    theta_rigid_plus=norm2(n_rigid)*(dble(i))
                    rodm_rigid=rodrigues_rot(u,theta_rigid_minus)
                    rodp_rigid=rodrigues_rot(u,theta_rigid_plus)
                    jacom_rigid=left_jacobian(u,theta_rigid_minus)/theta_rigid_minus
                    jacop_rigid=left_jacobian(u,theta_rigid_plus)/theta_rigid_plus
                else
                    rodm_rigid = 0.0d0
                    rodm_rigid(1,1)=1.0d0; rodm_rigid(2,2)=1.0d0; rodm_rigid(3,3)=1.0d0
                    rodp_rigid = 0.0d0
                    rodp_rigid(1,1)=1.0d0; rodp_rigid(2,2)=1.0d0; rodp_rigid(3,3)=1.0d0
                    jacom_rigid = 0.0d0
                    jacom_rigid(1,1)=1.0d0; jacom_rigid(2,2)=1.0d0; jacom_rigid(3,3)=1.0d0
                    jacop_rigid = 0.0d0
                    jacop_rigid(1,1)=1.0d0; jacop_rigid(2,2)=1.0d0; jacop_rigid(3,3)=1.0d0
                    theta_rigid_minus = 0.0d0; theta_rigid_plus = 0.0d0
                endif             
                arr_deformed_coordi_minus(c_start:c_end,:)=transpose(matmul(&
                rodm_rigid,transpose(arr_deformed_coordi_minus(c_start:c_end,:)-&
                spread(arr_chain_comm(j,:), dim=1, ncopies=c_leng))))+&
                spread(matmul(jacom_rigid,tm_rigid),dim=1,ncopies=c_leng)+&
                spread(arr_chain_comm(j,:), dim=1, ncopies=c_leng)
                arr_deformed_coordi_plus(c_start:c_end,:)=transpose(matmul(&
                rodp_rigid,transpose(arr_deformed_coordi_plus(c_start:c_end,:)-&
                spread(arr_chain_comp(j,:), dim=1, ncopies=c_leng))))+&
                spread(matmul(jacop_rigid,tp_rigid),dim=1,ncopies=c_leng)+&
                spread(arr_chain_comp(j,:), dim=1, ncopies=c_leng)
            enddo
            
            !!!!centering and rotating deformed coordi_BB
            !!!!minus direction
            center=0.0
            do s=1,ubound(arr_deformed_coordi_minus,1)
                center=center+arr_mass(s)*arr_deformed_coordi_minus(s,:)
            enddo
            center=center/M
            do s=1,ubound(arr_deformed_coordi_minus,1)
                arr_deformed_coordi_minus(s,:)=arr_deformed_coordi_minus(s,:)-center
            enddo
            rot_mat=0.0!;temp1=0.0;temp2=0.0
            call make_a_hat_kabsch(arr_ref_coordi,arr_deformed_coordi_minus,rot_mat) ! Kabasch method is used ***************
            arr_deformed_coordi_minus=transpose(matmul&
                    (rot_mat,transpose(arr_deformed_coordi_minus)))
            !!!!plus direction
            center=0.0
            do s=1,ubound(arr_deformed_coordi_plus,1)
                center=center+arr_mass(s)*arr_deformed_coordi_plus(s,:)
            enddo
            center=center/M
            do s=1,ubound(arr_deformed_coordi_plus,1)
                arr_deformed_coordi_plus(s,:)=arr_deformed_coordi_plus(s,:)-center
            enddo
            rot_mat=0.0!;temp1=0.0;temp2=0.0
            call make_a_hat_kabsch(arr_ref_coordi,arr_deformed_coordi_plus,rot_mat) ! Kabasch method is used ***************
            arr_deformed_coordi_plus=transpose(matmul&
                    (rot_mat,transpose(arr_deformed_coordi_plus)))
            arr_output_deformed_coordi_CC(:,:,mid - i,out_mode)=arr_deformed_coordi_minus(:,:)
            arr_output_deformed_coordi_CC(:,:,mid + i,out_mode)=arr_deformed_coordi_plus(:,:)
        enddo
    enddo
end subroutine
!***************************************************************************************************
function rodrigues_rot(n_vec,theta)
    double precision, dimension(3), intent(in) :: n_vec
    double precision, intent(in) :: theta
    double precision, dimension(3) :: normed_vec
    double precision, dimension(3,3) :: rodrigues_rot
    double precision, dimension(3,3) :: EYE3,skew_sym
    double precision, dimension(3,1) :: vec1
    EYE3=0.0d0; EYE3(1,1)=1.0d0; EYE3(2,2)=1.0d0; EYE3(3,3)=1.0d0
    skew_sym=0.0d0;
    normed_vec = n_vec/norm2(n_vec)
    skew_sym(1,2)=-normed_vec(3); skew_sym(1,3)=normed_vec(2)
    skew_sym(2,1)=normed_vec(3); skew_sym(2,3)=-normed_vec(1)
    skew_sym(3,1)=-normed_vec(2); skew_sym(3,2)=normed_vec(1)
    vec1(:,1)=normed_vec
    rodrigues_rot = EYE3*cos(theta)+(1-cos(theta))*matmul(vec1,transpose(vec1))+skew_sym*sin(theta)
end function
!***************************************************************************************************
function left_jacobian(n_vec,theta) ! actually get G (= theta * left_Jacovian)
    double precision, dimension(3), intent(in) :: n_vec
    double precision, intent(in) :: theta
    double precision, dimension(3) :: normed_vec
    double precision, dimension(3,3) :: left_jacobian
    double precision, dimension(3,3) :: EYE3,skew_sym
    double precision, dimension(3,1) :: vec1
    EYE3=0.0d0; EYE3(1,1)=1.0d0; EYE3(2,2)=1.0d0; EYE3(3,3)=1.0d0
    skew_sym=0.0d0;
    normed_vec = n_vec/norm2(n_vec)
    skew_sym(1,2)=-normed_vec(3); skew_sym(1,3)=normed_vec(2)
    skew_sym(2,1)=normed_vec(3); skew_sym(2,3)=-normed_vec(1)
    skew_sym(3,1)=-normed_vec(2); skew_sym(3,2)=normed_vec(1)
    vec1(:,1)=normed_vec
    if (theta < 1d-6) then
        left_jacobian = (theta-(theta**2)/6d0)*EYE3+(theta**2)/6d0*matmul(vec1,transpose(vec1))+&
                        0.5d0*theta*skew_sym
    else
        left_jacobian = sin(theta)*EYE3+(theta-sin(theta))*matmul(vec1,transpose(vec1))+&
                        (1d0-cos(theta))*skew_sym
    endif
end function
!***************************************************************************************************
subroutine compute_chain_com(arr_mass,arr_coordi,arr_chain_index,arr_chain_com,i_nchains)
    double precision, dimension(:), intent(in) :: arr_mass
    double precision, dimension(:,:), intent(in) :: arr_coordi
    integer, dimension(:), intent(in) :: arr_chain_index
    double precision, dimension(:,:), intent(out) :: arr_chain_com
    integer :: i_nchains,i,j,i_natom
    double precision :: chain_mass
    double precision, dimension(3,3) :: EYE3
    EYE3(:,:)=0.0;EYE3(1,1)=1.0;EYE3(2,2)=1.0;EYE3(3,3)=1.0
    i_natom = ubound(arr_coordi,1)
    do i = 1,i_nchains
        arr_chain_com(i,:)=0.0
        chain_mass=0.0
        do j=1, i_natom
            if (arr_chain_index(j) == i) then
                arr_chain_com(i,:) = arr_chain_com(i,:)+arr_mass(j)*arr_coordi(j,:)
                chain_mass=chain_mass+arr_mass(j)
            endif
        enddo
        arr_chain_com(i,:) = arr_chain_com(i,:) / chain_mass
    enddo
end subroutine
!***************************************************************************************************
subroutine make_a_hat_kabsch(arr_ref, arr_target, rot_mat)
    implicit none
    double precision, intent(in) :: arr_ref(:,:), arr_target(:,:)
    double precision, intent(out) :: rot_mat(3,3)
    double precision :: C(3,3), U(3,3), VT(3,3), S(3)
    double precision :: work(100)
    integer :: i, j, k, info, lwork

    ! Step 1: Compute covariance matrix C = B^T * A
    C = 0.0
    do i = 1, size(arr_ref,1)
        do j = 1, 3
            do k = 1, 3
                C(j,k) = C(j,k) + arr_target(i,j) * arr_ref(i,k)
            enddo
        enddo
    enddo

    ! Step 2: Compute SVD: C = U * diag(S) * VT
    lwork = 100
    call dgesvd('A','A',3,3,C,3,S,U,3,VT,3,work,lwork,info)
    if (info /= 0) then
        print *, 'Error in DGESVD, info = ', info
        stop
    endif

    ! Step 3: Check for reflection, correct if necessary
    if (det3x3(matmul(VT,transpose(U))) < 0.0d0) then
        VT(3,:) = -VT(3,:)  ! flip sign of last row of V^T
    endif

    ! Step 4: Compute rotation matrix: R = V * U^T = (VT^T) * U^T
    rot_mat = matmul(transpose(VT), transpose(U))
end subroutine
!***************************************************************************************************
double precision function det3x3(mat)
double precision, intent(in) :: mat(3,3)
det3x3 = mat(1,1)*(mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) &
       - mat(1,2)*(mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) &
       + mat(1,3)*(mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))
end function
!***************************************************************************************************
subroutine cal_theo_bfactor(evec,eval,arr_theo_bfactor,arr_RMSF)
    double precision, dimension(:,:), intent(in) :: evec
    double precision, dimension(:), intent(in) :: eval
    double precision, dimension(:), intent(out) :: arr_theo_bfactor
    double precision, dimension(:), intent(out) :: arr_RMSF
    double precision :: T=303.15,pi=4.0*atan(1.0),kb,eigval_scale_factor    !kb=1.380649*(10.0**(-3)) ! with unit: kg*A^2/(s^2*T)
    double precision :: temp,CA_Da_mass,amino_mass
    integer :: i,j,i_nmodes
    eigval_scale_factor=4.184*(10.0**26) ! Kcal/(mol*A^2*Da) to 1/s^2
    kb=8.314*(10.0**23) ! with unit: Da*A^2/(s^2*T)
    CA_Da_mass=12.0
    amino_mass=105.0
    i_nmodes=ubound(eval,1)
    do i=1,ubound(evec,1)/3 ! the number of atoms
        temp=0
        do j=1,i_nmodes ! the number of dihedral angles having DOF
            temp=temp+(eval(j)**(-1))*( evec(3*i-2,j)**2+evec(3*i-1,j)**2+evec(3*i,j)**2)
        enddo
        arr_theo_bfactor(i)=8.0/3.0*(pi**2)*kb*T*temp/(CA_Da_mass*eigval_scale_factor)
    enddo
    
    do i=1,ubound(arr_theo_bfactor,1)
        arr_RMSF(i)=sqrt(3.0/(8.0*(pi**2))*arr_theo_bfactor(i))
    enddo
end subroutine
!***************************************************************************************************
subroutine make_mass_array(arr_atom_type,arr_resid_type,arr_mass)
    integer, dimension(:,:), intent(in) :: arr_atom_type
    character(3), dimension(:), intent(in) :: arr_resid_type
    double precision, dimension(:), intent(out) :: arr_mass
    double precision :: mass_N,mass_C,mass_O,mass_ZN
    integer :: i, i_nline
    ! atomic weights is expressed atomic mass units (Daltons). C including Oxygen mass, CA including side chain's mass
    i_nline = ubound(arr_resid_type, 1)
    mass_N=14.0;mass_C=12.0;mass_O=16.0;mass_ZN=65.4
    do i = 1, i_nline
        if (arr_atom_type(i,2).eq.1 .or. arr_atom_type(i,2).eq.11) then
            arr_mass(i) = mass_N
        else if (arr_atom_type(i,2).eq.2 .or. arr_atom_type(i,2).eq.22) then
            arr_mass(i) = mass_C
        else if (arr_atom_type(i,2).eq.3 .or. arr_atom_type(i,2).eq.33) then
            arr_mass(i) = mass_C+mass_O
        else if (arr_atom_type(i,2).eq.5) then
            if (arr_resid_type(i) == 'ZN') then
                arr_mass(i) = mass_ZN
            endif
        else
            if (arr_resid_type(i) == 'ALA') then
                arr_mass(i) = 82.0-mass_N-2*mass_C-2*mass_O ! whole weight - N CA COOH
            else if (arr_resid_type(i) == 'ARG') then
                arr_mass(i) = 160.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'ASN') then
                arr_mass(i) = 124.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'ASP') then
                arr_mass(i) = 126.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'CYS') then
                arr_mass(i) = 114.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'GLU') then
                arr_mass(i) = 138.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'GLN') then
                arr_mass(i) = 136.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'GLY') then
                arr_mass(i) = 70.0-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'HIS') then
                arr_mass(i) = 146.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'HYP') then
                arr_mass(i) = 122.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'ILE') then
                arr_mass(i) = 118.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'LEU') then
                arr_mass(i) = 118.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'LYS') then
                arr_mass(i) = 132.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'MET') then
                arr_mass(i) = 138.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'PHE') then
                arr_mass(i) = 154.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'PRO') then
                arr_mass(i) = 106.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'SER') then
                arr_mass(i) = 98.0-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'THR') then
                arr_mass(i) = 110.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'TRP') then
                arr_mass(i) = 192.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'TYR') then
                arr_mass(i) = 170.1-mass_N-2*mass_C-2*mass_O
            else if (arr_resid_type(i) == 'VAL') then
                arr_mass(i) = 106.1-mass_N-2*mass_C-2*mass_O
            else
                arr_mass(i) = 102.0-mass_N-2*mass_C-2*mass_O ! average weight of amino acid is 110 DA > Heavy atom only: 110 * 0.93 (average proportion of heavy atom's weight to whole weight)
            end if
        endif
    enddo   
end subroutine
!***************************************************************************************************
subroutine print_usage()
    implicit none
    write(*,'(A)') "Usage:"
    write(*,'(A)') "  ICed_ENM_NMA <PDB> <chain> [options]"
    write(*,'(A)') ""
    write(*,'(A)') "Required arguments:"
    write(*,'(A)') "  <PDB>       Input PDB file path."
    write(*,'(A)') "  <chain>     Chain ID(s) to analyze, e.g. A or ABC."
    write(*,'(A)') ""
    write(*,'(A)') "Options:"
    write(*,'(A)') "  --core N         Number of CPU cores to use. Default: 4."
    write(*,'(A)') "  --mode N         Number of lowest normal modes to calculate. Default: 3."
    write(*,'(A)') "                   Use --mode 0 to calculate the full set of normal modes."
    write(*,'(A)') "  --cutoff X       Cartesian cutoff distance in Angstrom. Default: 8.0."
    write(*,'(A)') "  --out-prefix S   Prefix output file names with S. Default: no prefix."
    write(*,'(A)') "                   If S does not end with _, -, or /, '_' is added automatically."
    write(*,'(A)') "  --write-IC       Also write evec_IC.txt."
    write(*,'(A)') "  --write-CC       Also write evec_CC_raw.txt, evec_CC.txt, and ref_CC.pdb."
    write(*,'(A)') "  --write-CA-raw   Also write evec_CA_raw.txt."
    write(*,'(A)') "  --write-variance Also write variance_cumulative.txt and variance.txt."
    write(*,'(A)') "  --write-all      Write all optional output files."
    write(*,'(A)') "  --movie W F      Generate movie coordinates with weight W and frame count F."
    write(*,'(A)') "                   F must be an odd integer >= 3."
    write(*,'(A)') "  --movie-mode N   With --movie, generate movies from mode 1 to mode N."
    write(*,'(A)') "                   Default: 1. N must be <= the number of calculated modes."
    write(*,'(A)') "  --movie-only N   With --movie, generate a movie only for mode N."
    write(*,'(A)') "                   N must be <= the number of calculated modes."
    write(*,'(A)') "                   Do not use together with --movie-mode."
    write(*,'(A)') "  -h, --help       Show this help message."
    write(*,'(A)') ""
    write(*,'(A)') "Examples:"
    write(*,'(A)') "  ICed_ENM_NMA 1CRN.pdb A"
    write(*,'(A)') "  ICed_ENM_NMA 1CRN.pdb A --core 8 --mode 20 --cutoff 10.0"
    write(*,'(A)') "  ICed_ENM_NMA mutant.pdb A --mode 50 --out-prefix A_35"
    write(*,'(A)') "  ICed_ENM_NMA 1CRN.pdb A --write-CC --write-variance"
    write(*,'(A)') "  ICed_ENM_NMA 1CRN.pdb A --mode 0"
    write(*,'(A)') "  ICed_ENM_NMA 1CRN.pdb A --movie 1.0 21"
    write(*,'(A)') "  ICed_ENM_NMA 1CRN.pdb A --mode 10 --movie 1.0 21 --movie-only 7"
end subroutine
end program

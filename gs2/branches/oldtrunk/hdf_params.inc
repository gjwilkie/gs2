!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                        Grid parameters                         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! file name
  character (500)           :: fname_grids
  character (6),  parameter :: grids_filepref = '_grids'
  ! group name
  character (6),  parameter :: gname_grids = '/Grids'
  character (11), parameter :: gname_params = '/Parameters'
  ! datasets name
  character (2),  parameter :: dname_nx = 'nx'
  character (2),  parameter :: dname_ny = 'ny'
  character (2),  parameter :: dname_nz = 'nz'
  character (6),  parameter :: dname_negrid = 'negrid'
  character (7),  parameter :: dname_nlambda = 'nlambda'
  character (5),  parameter :: dname_nspec = 'nspec'
  character (7),  parameter :: dname_ntheta0 = 'ntheta0'
  character (4),  parameter :: dname_naky = 'naky'
  character (6),  parameter :: dname_ntgrid = 'ntgrid'
  character (2),  parameter :: dname_x0 = 'x0'
  character (2),  parameter :: dname_y0 = 'y0'
  character (2),  parameter :: dname_z0 = 'z0'
  character (2),  parameter :: dname_kx = 'kx'
  character (2),  parameter :: dname_ky = 'ky'
  character (1),  parameter :: dname_energy = 'e'
  character (2),  parameter :: dname_al = 'al'
  character (1),  parameter :: dname_w = 'w'
  character (2),  parameter :: dname_wl = 'wl'
  character (6),  parameter :: dname_layout = 'layout'
  character (11), parameter :: dname_accel = 'accelerated'

  ! file identifier
  integer (HID_T) :: grids_file_id
  ! group identifier
  integer (HID_T) :: grp_id


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                 Distribution function parameters               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! file name
  character (500)          :: fname_dist
  character (6), parameter :: dist_filepref = '_dist_'
  character (500)          :: fname_field
  character (7), parameter :: field_filepref = '_field_'
  ! datasets name
  character (10)           :: dname_pdf
  character (3), parameter :: dname_gk = 'g_k'
  character (3), parameter :: dname_gx = 'g_x'
  character (3), parameter :: dname_hk = 'h_k'
  character (4), parameter :: dname_time = 'time'
  character (5), parameter :: dname_phik = 'phi_k'
  character (6), parameter :: dname_apark = 'Apar_k'
  character (6), parameter :: dname_bpark = 'Bpar_k'

  integer (HID_T) :: dist_file_id    ! File identifier
  integer (HID_T) :: field_file_id    ! File identifier
  integer (HID_T) :: popen_id        ! property list id
  integer (HID_T) :: pwrite_id       ! property list id
  integer (HID_T) :: dsp_pdf         ! data space id
  integer (HID_T) :: dst_pdf         ! datasets id
  integer (HID_T) :: mem_pdf         ! memory id

  ! rank and dimensions
! TT>
!  integer, parameter :: rank_pdf= 7    ! Dataset rank
!  integer, parameter :: rank_mem= 2    ! curent_datum rank
  integer            :: rank_pdf      ! Dataset rank
  integer            :: rank_mem      ! memory data rank
! <TT
! TT>
!  integer (HSSIZE_T), dimension (rank_pdf) :: offset
!  integer (HSIZE_T), dimension (rank_pdf) :: count
!  integer (HSIZE_T), dimension (rank_pdf) :: dim_pdf ! Dataset dimensions
!  integer (HSSIZE_T), dimension (rank_mem) :: offset_mem, count_mem
!  integer (HSIZE_T), dimension (rank_mem) :: offset_mem, count_mem
  integer (HSIZE_T), dimension (:), allocatable :: offset, count
  integer (HSIZE_T), dimension (:), allocatable :: stride, block
  integer (HSIZE_T), dimension (:), allocatable :: offset_mem, count_mem
  integer (HSIZE_T), dimension (:), allocatable :: dim_pdf, dim_g
! <TT

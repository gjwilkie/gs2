# include "define.inc"

module gs2_dist_io

!
! This module ported from AstroGK by Tomo Tatsuno
!

! To do: flush right before serial file close
!        add TEST for test message head
!        support non-square array for real space output

# ifdef HDF
# if FCOMPILER == _GFORTRAN_
  use hdf5
# else
  use hdf5, only: HID_T, HSIZE_T!, HSSIZE_T
# endif
# endif
  implicit none

  public :: hdf_write_dist, hdf_finish_dist

  private

!!$ Defines names, params for the HDF5 datasets, files etc.
!!$ ==George==
# ifdef HDF
  include 'hdf_params.inc'
# endif

  logical, save :: accelerated
  logical, parameter :: export_complex = .true.
  logical, parameter :: export_h = .true.
  logical, save :: square ! non-square output works only if export_complex = T.
  logical :: test = .false.

contains

  subroutine hdf_write_dist (g0)

    use theta_grid, only: ntgrid
    use gs2_layouts, only: yxf_lo, accelx_lo, g_lo
    use file_utils, only: error_unit, run_name
    use mp, only: proc0
# ifdef HDF
    use mp, only: nproc, barrier, iproc
    use convert, only: c2r
    use gs2_time, only: user_time
    use fields_arrays, only: phinew, aparnew, bparnew
    use gs2_transforms, only: transform2
    use run_parameters, only: fphi, fapar, fbpar
    use collisions, only: g_adjust
! TT> for test
    use gs2_layouts, only: is_idx, ie_idx, il_idx, ik_idx, it_idx
    use le_grids, only: e, al
! <TT
# if FCOMPILER != _GFORTRAN_
    use hdf5, only: H5F_ACC_TRUNC_F, H5F_ACC_RDWR_F, H5S_SELECT_SET_F
    use hdf5, only: H5F_SCOPE_LOCAL_F, h5fflush_f
    use hdf5, only: h5fcreate_f, h5fopen_f, h5fclose_f
    use hdf5, only: h5screate_simple_f, h5sselect_hyperslab_f, h5sclose_f
    use hdf5, only: h5dcreate_f, h5dwrite_f, h5dclose_f
# endif
    use hdf_wrapper, only: hdf_write, hdf_error, hdf_file_real, hdf_mem_real
!    use mpi,  only: mpi_wtime
    use job_manage, only: timer => timer_local
# endif
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g0
    integer, save :: error=0
# ifdef HDF
    integer, save :: cnt=0
    integer :: ia, ierr
    real :: t_start, t_finish, zero
    real,    dimension (:,:),     allocatable :: gx
    real,    dimension (:,:,:),   allocatable :: agx
    real,    dimension (:,:,:,:), allocatable :: gkr
    complex, dimension (:,:,:),   allocatable :: g1
    character (3) :: suffix
! TT> test
    integer :: iac, ik, it, il, ie, iglo, iyxf
! <TT
# endif

    call hdf_init_dist (error)
    if (error < 0) return
    if (proc0 .and. test) print *, 'hdf_init_dist done'

# ifdef HDF
    allocate (g1 (-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))
    g1 = g0

    if (export_h) call g_adjust (g1, phinew, bparnew, fphi, fbpar)

! For now, assume that we have a box domain and a square-form layout array.

    ! determine file name
    cnt = cnt + 1
    if (cnt > 999) then
       write (error_unit(),*) &
            'WARNING: Too many hdf_write_dist calls.  Skip writes.'
       error = -1
       return
    end if
    write (suffix, '(i3.3)') cnt
    write (fname_dist, '(a)') trim(run_name) // trim(dist_filepref) &
         & // trim(suffix) // '.hdf'
!!$    write (fname_field, '(a)') trim(run_name) // trim(field_filepref) &
!!$         & // trim(suffix) // '.hdf'

!!$    ! create file, write time and field and close it
!!$    ! TT: may not need to be closed?
!!$    if (proc0) then
!!$       call h5fcreate_f (fname_field, H5F_ACC_TRUNC_F, field_file_id, error)
!!$       if (error >= 0) call hdf_write (field_file_id, dname_time, user_time, error)
!!$       if (export_complex) then
!!$          if (use_Phi .and. error>=0) &
!!$               call hdf_write (field_file_id, dname_phik, phinew, error)
!!$               ! intentionally create error
!!$!               call hdf_write (grp_id, dname_phik, phinew, error)
!!$          if (use_Apar .and. error>=0) &
!!$               call hdf_write (field_file_id, dname_apark, aparnew, error)
!!$          if (use_Bpar .and. error>=0) &
!!$               call hdf_write (field_file_id, dname_bpark, bparnew, error)
!!$       end if
!!$!       if (error >= 0) call h5fflush_f (dist_file_id, H5F_SCOPE_LOCAL_F, error)
!!$!       if (error < 0) call hdf_error (message='hdf serial file flush error')
!!$       call h5fclose_f (field_file_id, error)
!!$       if (test) print *, 'hdf_write_dist: serial write done'
!!$    end if

    if (export_complex) then
       allocate (gkr(2, -ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))

       call c2r (g1, gkr)
       if (test .and. proc0) print *, 'hdf_write_dist: c2r done'

    else

       ! transform data to real space
       if (accelerated) then
          allocate (agx(-ntgrid:ntgrid, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
          agx = 0.
          call transform2 (g1, agx, ia)
       else
          allocate (gx(yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
          gx = 0.
          call transform2 (g1, gx)
       end if

    end if

    deallocate (g1)

! TT> write artificial data for test: agreed.
!!$    if (test) then
!!$       if (export_complex) then
!!$          gkr = 0.0
!!$          do iglo = g_lo%llim_proc, g_lo%ulim_proc
!!$             ik = ik_idx (g_lo, iglo)
!!$             it = it_idx (g_lo, iglo)
!!$             il = il_idx (g_lo, iglo)
!!$             ie = ie_idx (g_lo, iglo)
!!$             if ((ik /= 1) .or. (it /= 1)) cycle
!!$             gkr(1,:,:,iglo) = cos(al(il)) * sin(e(ie,1))
!!$             gkr(2,:,:,iglo) = sin(al(il)) * sin(e(ie,1))
!!$          end do
!!$       else
!!$          if (accelerated) then
!!$             do iac = accelx_lo%llim_proc, accelx_lo%ulim_proc
!!$!                ik = ik_idx (accelx_lo, iac)
!!$!                it = it_idx (accelx_lo, iac)
!!$                il = il_idx (accelx_lo, iac)
!!$                ie = ie_idx (accelx_lo, iac)
!!$                agx(:,:,iac) = cos(al(il)) * sin(e(ie,1))
!!$             end do
!!$          else
!!$             do iyxf = yxf_lo%llim_proc, yxf_lo%ulim_proc
!!$!                it = it_idx (yxf_lo, iyxf)
!!$                il = il_idx (yxf_lo, iyxf)
!!$                ie = ie_idx (yxf_lo, iyxf)
!!$                gx(:,iyxf) = cos(al(il)) * sin(e(ie,1))
!!$             end do
!!$          end if
!!$       end if
!!$    end if
! <TT

!    call barrier

    ! open file for parallel writing
    call h5fcreate_f (fname_dist, H5F_ACC_TRUNC_F, dist_file_id, error, &
         access_prp = popen_id)
!!$    call h5fopen_f (fname_dist, H5F_ACC_RDWR_F, dist_file_id, error, &
!!$         access_prp = popen_id)
    if (error < 0) call hdf_error &
         (message='ERROR: in parallel file open', iproc=iproc)
    if (test .and. proc0) print *, 'hdf_write_dist: parallel file open done'

    ! create datasets for dist func
    if (error >= 0) then
       call h5dcreate_f (dist_file_id, dname_pdf, hdf_file_real, dsp_pdf, &
            dst_pdf, error)
       if (error < 0) call hdf_error &
            (message='ERROR: in dist func dataset creation', iproc=iproc)
       if (test .and. proc0) print *, 'hdf_write_dist: dist dataset creation done'
    end if

!    if (test .and. proc0) t_start = mpi_wtime()
    if (test .and. proc0) t_start = timer()

    if (export_complex) then

       if (error >= 0) then
          call h5dwrite_f (dst_pdf, hdf_mem_real, gkr, dim_g, error, &
               file_space_id=dsp_pdf, mem_space_id=mem_pdf, xfer_prp=pwrite_id)
          if (error < 0) call hdf_error &
               (message='ERROR: in parallel write', iproc=iproc)
          if (test .and. proc0) print *, 'hdf_write_dist: parallel write done'
       end if

! minor bug fix here, associated with memory diagnostics.  8/12/09  BD
       deallocate (gkr)

    else

       if (accelerated) then

          if (error >= 0) then
             call h5dwrite_f (dst_pdf, hdf_mem_real, agx, dim_g, error, &
                  file_space_id=dsp_pdf, mem_space_id=mem_pdf, &
                  xfer_prp=pwrite_id)
             if (error < 0) call hdf_error &
                  (message='ERROR: in hdf parallel write', iproc=iproc)
          end if

          deallocate (agx)

       else

          if (error >= 0) then
             call h5dwrite_f (dst_pdf, hdf_mem_real, gx, dim_g, error, &
                  file_space_id=dsp_pdf, mem_space_id=mem_pdf, &
                  xfer_prp=pwrite_id)
             if (error < 0) call hdf_error &
                  (message='ERROR: in hdf parallel write', iproc=iproc)
          end if

          deallocate (gx)

       end if

    end if

    if (test .and. proc0) then
!       t_finish = mpi_wtime()
       t_finish = timer()
       print *, 'HDF total I/O time = ', t_finish - t_start, ' seconds'
    end if

!    call barrier

    call h5dclose_f (dst_pdf, ierr)
    call h5fclose_f (dist_file_id, ierr)

    ! open file, write time and field and close it
    ! TT: may not need to be closed?
    if (proc0) then
       zero = epsilon(0.0)
       call h5fopen_f (fname_dist, H5F_ACC_RDWR_F, dist_file_id, error)
       if (error >= 0) call hdf_write (dist_file_id, dname_time, user_time, error)
       if (export_complex) then
          if (fphi>zero .and. error>=0) &
               call hdf_write (dist_file_id, dname_phik, phinew, error)
               ! intentionally create error
!               call hdf_write (grp_id, dname_phik, phinew, error)
          if (fapar>zero .and. error>=0) &
               call hdf_write (dist_file_id, dname_apark, aparnew, error)
          if (fbpar>zero .and. error>=0) &
               call hdf_write (dist_file_id, dname_bpark, bparnew, error)
       end if
!       if (error >= 0) call h5fflush_f (dist_file_id, H5F_SCOPE_LOCAL_F, error)
!       if (error < 0) call hdf_error (message='hdf serial file flush error')
       call h5fclose_f (dist_file_id, error)
       if (test) print *, 'hdf_write_dist: serial write done'
    end if

# endif
  end subroutine hdf_write_dist

  subroutine hdf_init_dist (error)

    use file_utils, only: error_unit
    use mp, only: proc0
# ifdef HDF
    use mpi, only: MPI_INFO_NULL
    use mp, only: communicator, nproc, iproc
    use kt_grids, only: nx, ny, ntheta0, naky, box
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use gs2_layouts, only: layout, yxf_lo, g_lo, accelx_lo
    use gs2_layouts, only: is_idx, ie_idx, il_idx, ik_idx, it_idx, isign_idx, ig_idx
    use gs2_transforms, only: init_transforms
# if FCOMPILER != _GFORTRAN_
    use hdf5, only: H5P_FILE_ACCESS_F, H5S_SELECT_SET_F
    use hdf5, only: H5P_DATASET_XFER_F, H5FD_MPIO_INDEPENDENT_F
    use hdf5, only: h5pcreate_f, h5pset_fapl_mpio_f, h5pset_dxpl_mpio_f
    use hdf5, only: h5screate_simple_f, h5dcreate_f, h5sselect_hyperslab_f
# endif
# endif
    implicit none
    integer, intent (inout) :: error
# ifdef HDF
    logical, save :: initialized=.false.
    integer :: ik_min, ik_max, it_min, it_max, ig_min, ig_max
    integer :: il_min, il_max, ie_min, ie_max
    integer :: isign_min, isign_max, ispec_min, ispec_max, g_size

    if (initialized .or. error < 0) return

    if (.not.export_complex .and. .not.box) then
       if (proc0) write (error_unit(),*) &
            'ERROR: We support real space dist_fn output for box domain only.'
       error = -1
       return
    end if

    if (export_complex) then

       ! set dataname
       dname_pdf = dname_gk
       if (export_h) dname_pdf = dname_hk

       ! memory dist func dimensions
       rank_mem = 4
       allocate (dim_g(rank_mem), offset_mem(rank_mem), count_mem(rank_mem))
       g_size = max (g_lo%ulim_proc - g_lo%llim_proc + 1, 0)
       dim_g = (/ 2, ntgrid*2+1, 2, g_size /)
       offset_mem = 0
       count_mem = dim_g

       ! file dist func dimensions
       ! first check if array structure is square type
       select case (layout)
       case ('yxels')
          if (nproc <= nspec) then
             square = mod(nspec,nproc)==0
          else if (nproc <= nspec*nlambda) then
               square = mod(nspec*nlambda, nproc)==0 .and. mod(nproc, nspec)==0
          else if (nproc <= nspec*nlambda*negrid) then
               square = mod(nspec*nlambda*negrid, nproc)==0 &
                    .and. mod(nproc, nspec*nlambda)==0
          else if (nproc <= nspec*nlambda*negrid*ntheta0) then
               square = mod(nspec*nlambda*negrid*ntheta0, nproc)==0 &
                    .and. mod(nproc, nspec*nlambda*negrid)==0
          else if (nproc <= nspec*nlambda*negrid*ntheta0*naky) then
               square = mod(nspec*nlambda*negrid*ntheta0*naky, nproc)==0 &
                    .and. mod(nproc, nspec*nlambda*negrid*ntheta0)==0
          end if
       case ('yxles')
          if (nproc <= nspec) then
             square = mod(nspec,nproc)==0
          else if (nproc <= nspec*negrid) then
               square = mod(nspec*negrid, nproc)==0 .and. mod(nproc, nspec)==0
          else if (nproc <= nspec*negrid*nlambda) then
               square = mod(nspec*negrid*nlambda, nproc)==0 &
                    .and. mod(nproc, nspec*negrid) == 0
          else if (nproc <= nspec*negrid*nlambda*ntheta0) then
               square = mod(nspec*negrid*nlambda*ntheta0, nproc)==0 &
                    .and. mod(nproc, nspec*negrid*nlambda)==0
          else if (nproc <= nspec*negrid*nlambda*ntheta0*naky) then
               square = mod(nspec*negrid*nlambda*ntheta0*naky, nproc)==0 &
                    .and. mod(nproc, nspec*negrid*nlambda*ntheta0)==0
          end if
       case ('lexys')
          if (nproc <= nspec) then
             square = mod(nspec,nproc)==0
          else if (nproc <= nspec*naky) then
               square = mod(nspec*naky, nproc)==0 .and. mod(nproc, nspec)==0
          else if (nproc <= nspec*naky*ntheta0) then
               square = mod(nspec*naky*ntheta0, nproc)==0 &
                    .and. mod(nproc, nspec*naky)==0
          else if (nproc <= nspec*naky*ntheta0*negrid) then
               square = mod(nspec*naky*ntheta0*negrid, nproc)==0 &
                    .and. mod(nproc, nspec*naky*ntheta0)==0
          else if (nproc <= nspec*naky*ntheta0*negrid*nlambda) then
               square = mod(nspec*naky*ntheta0*negrid*nlambda, nproc)==0 &
                    .and. mod(nproc, nspec*naky*ntheta0*negrid)==0
          end if
       case ('lxyes')
          if (nproc <= nspec) then
             square = mod(nspec,nproc)==0
          else if (nproc <= nspec*negrid) then
               square = mod(nspec*negrid, nproc)==0 .and. mod(nproc, nspec)==0
          else if (nproc <= nspec*negrid*naky) then
               square = mod(nspec*negrid*naky, nproc)==0 &
                    .and. mod(nproc, nspec*negrid)==0
          else if (nproc <= nspec*negrid*naky*ntheta0) then
               square = mod(nspec*negrid*naky*ntheta0, nproc)==0 &
                    .and. mod(nproc, nspec*negrid*naky)==0
          else if (nproc <= nspec*negrid*naky*ntheta0*nlambda) then
               square = mod(nspec*negrid*naky*ntheta0*nlambda, nproc)==0 &
                    .and. mod(nproc, nspec*negrid*naky*ntheta0)==0
          end if
       case ('lyxes')
          if (nproc <= nspec) then
             square = mod(nspec,nproc)==0
          else if (nproc <= nspec*negrid) then
               square = mod(nspec*negrid, nproc)==0 .and. mod(nproc, nspec)==0
          else if (nproc <= nspec*negrid*ntheta0) then
               square = mod(nspec*negrid*ntheta0, nproc)==0 &
                    .and. mod(nproc, nspec*negrid)==0
          else if (nproc <= nspec*negrid*ntheta0*naky) then
               square = mod(nspec*negrid*ntheta0*naky, nproc)==0 &
                    .and. mod(nproc, nspec*negrid*ntheta0)==0
          else if (nproc <= nspec*negrid*ntheta0*naky*nlambda) then
               square = mod(nspec*negrid*ntheta0*naky*nlambda, nproc)==0 &
                    .and. mod(nproc, nspec*negrid*ntheta0*naky)==0
          end if
       end select

       if (test .and. proc0) print *, 'square = ', square

       if (square) then
          rank_pdf = 8
          allocate (dim_pdf(rank_pdf), offset(rank_pdf), count(rank_pdf))

          ! TT: following should work if array structure is square
          ik_min = ik_idx (g_lo, g_lo%llim_proc)
          ik_max = ik_idx (g_lo, g_lo%ulim_proc)
          it_min = it_idx (g_lo, g_lo%llim_proc)
          it_max = it_idx (g_lo, g_lo%ulim_proc)
          il_min = il_idx (g_lo, g_lo%llim_proc)
          il_max = il_idx (g_lo, g_lo%ulim_proc)
          ie_min = ie_idx (g_lo, g_lo%llim_proc)
          ie_max = ie_idx (g_lo, g_lo%ulim_proc)
          ispec_min = is_idx (g_lo, g_lo%llim_proc)
          ispec_max = is_idx (g_lo, g_lo%ulim_proc)

          dim_pdf(1:3) = (/ 2, ntgrid*2+1, 2 /)
          offset(1:3) = 0
          count (1:3) = dim_pdf(1:3)
          select case (layout)
          case ('yxels')
             dim_pdf(4:7) = (/ naky, ntheta0, negrid, nlambda /)
             offset(4:7) = (/ ik_min-1, it_min-1, ie_min-1, il_min-1 /)
             count (4:7) = (/ ik_max, it_max, ie_max, il_max /) - offset(4:7)
          case ('yxles')
             dim_pdf(4:7) = (/ naky, ntheta0, nlambda, negrid /)
             offset(4:7) = (/ ik_min-1, it_min-1, il_min-1, ie_min-1 /)
             count(4:7) = (/ ik_max, it_max, il_max, ie_max /) - offset(4:7)
          case ('lexys')
             dim_pdf(4:7) = (/ nlambda, negrid, ntheta0, naky /)
             offset(4:7) = (/ il_min-1, ie_min-1, it_min-1, ik_min-1 /)
             count(4:7) = (/ il_max, ie_max, it_max, ik_max /) - offset(4:7)
          case ('lxyes')
             dim_pdf(4:7) = (/ nlambda, ntheta0, naky, negrid /)
             offset(4:7) = (/ il_min-1, it_min-1, ik_min-1, ie_min-1 /)
             count(4:7) = (/ il_max, it_max, ik_max, ie_max /) - offset(4:7)
          case ('lyxes')
             dim_pdf(4:7) = (/ nlambda, naky, ntheta0, negrid /)
             offset(4:7) = (/ il_min-1, ik_min-1, it_min-1, ie_min-1 /)
             count(4:7) = (/ il_max, ik_max, it_max, ie_max /) - offset(4:7)
          end select
          dim_pdf(8) = nspec
          offset(8) = ispec_min - 1
          count (8) = ispec_max - offset(8)

          if (test) print *, 'iproc=', iproc, ' offset=', offset, &
               ' count=', count

       else

          ! if not square, write g_lo as it is
          rank_pdf = rank_mem
          allocate (dim_pdf(rank_pdf), offset(rank_pdf), count(rank_pdf))
          dim_pdf(1:3) = dim_g(1:3)
          dim_pdf (4) = g_lo%ulim_world - g_lo%llim_world + 1
          offset(1:3) = 0
          offset(4) = g_lo%llim_proc
          count = count_mem

          if (test) print *, 'iproc=', iproc, ' offset=', offset, &
               ' count=', count

       end if

    else ! if not export_complex

       ! set dataname
       dname_pdf = dname_gx

! TT> moved to write_grids
       ! initialize FFT and accel flag
       call init_transforms &
            (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)
! <TT

       if (accelerated) then

          ! memory dist func dimensions
          rank_mem = 3
          allocate (dim_g(rank_mem), offset_mem(rank_mem), count_mem(rank_mem))
          g_size = max (accelx_lo%ulim_proc - accelx_lo%llim_proc + 1, 0)
          dim_g = (/ ntgrid*2+1, 2, g_size /)
          offset_mem = 0
          count_mem = dim_g

          ! file dist func dimensions
          rank_pdf = 7
          allocate (dim_pdf(rank_pdf), offset(rank_pdf), count(rank_pdf))

          ! TT: following should work if array structure is square
          ik_min = ik_idx (accelx_lo, accelx_lo%llim_proc)
          ik_max = ik_idx (accelx_lo, accelx_lo%ulim_proc)
          it_min = it_idx (accelx_lo, accelx_lo%llim_proc)
          it_max = it_idx (accelx_lo, accelx_lo%ulim_proc)
          il_min = il_idx (accelx_lo, accelx_lo%llim_proc)
          il_max = il_idx (accelx_lo, accelx_lo%ulim_proc)
          ie_min = ie_idx (accelx_lo, accelx_lo%llim_proc)
          ie_max = ie_idx (accelx_lo, accelx_lo%ulim_proc)
          ispec_min = is_idx (accelx_lo, accelx_lo%llim_proc)
          ispec_max = is_idx (accelx_lo, accelx_lo%ulim_proc)

          dim_pdf(1:4) = (/ ntgrid*2+1, 2, ny, nx /)
          offset(1:4) = (/ 0, 0, ik_min-1, it_min-1 /)
          count (1:4) = (/ ntgrid*2+1, 2, ik_max-ik_min+1, it_max-it_min+1 /)
          select case (layout)
          case ('yxels')
             dim_pdf(5:6) = (/ negrid, nlambda /)
             offset(5:6) = (/ ie_min-1, il_min-1 /)
             count (5:6) = (/ ie_max, il_max /) - offset(5:6)
          case ('yxles')
             dim_pdf(5:6) = (/ nlambda, negrid /)
             offset(5:6) = (/ il_min-1, ie_min-1 /)
             count (5:6) = (/ il_max, ie_max /) - offset(5:6)
          case default
             ! There should be no other accelerated layout
             print *, 'accelerated layout? ', trim(layout)
          end select
          dim_pdf(7) = nspec
          offset(7) = ispec_min - 1
          count (7) = ispec_max - offset(7)

          if ( g_size /= product(count(3:7)) ) print*, 'what to do?'

       else ! if not accelerated

          ! memory dist func dimensions
          rank_mem = 2
          allocate (dim_g(rank_mem), offset_mem(rank_mem), count_mem(rank_mem))
          g_size = max (yxf_lo%ulim_proc - yxf_lo%llim_proc + 1, 0)
          dim_g = (/ ny, g_size /)
          offset_mem = 0
          count_mem = dim_g

          ! file dist func dimensions
          rank_pdf = 7
          allocate (dim_pdf(rank_pdf), offset(rank_pdf), count(rank_pdf))

          ! TT: following should work if array structure is square
          it_min = it_idx (yxf_lo, yxf_lo%llim_proc)
          it_max = it_idx (yxf_lo, yxf_lo%ulim_proc)
          ig_min = ig_idx (yxf_lo, yxf_lo%llim_proc)
          ig_max = ig_idx (yxf_lo, yxf_lo%ulim_proc)
          isign_min = isign_idx(yxf_lo, yxf_lo%llim_proc)
          isign_max = isign_idx(yxf_lo, yxf_lo%ulim_proc)
          il_min = il_idx (yxf_lo, yxf_lo%llim_proc)
          il_max = il_idx (yxf_lo, yxf_lo%ulim_proc)
          ie_min = ie_idx (yxf_lo, yxf_lo%llim_proc)
          ie_max = ie_idx (yxf_lo, yxf_lo%ulim_proc)
          ispec_min = is_idx (yxf_lo, yxf_lo%llim_proc)
          ispec_max = is_idx (yxf_lo, yxf_lo%ulim_proc)

          dim_pdf(1:4) = (/ ny, nx, ntgrid*2+1, 2/)
          offset(1:4) = (/ 0, it_min-1, ig_min+ntgrid, isign_min-1 /)
          count (1:4) = (/ ny, it_max-it_min+1, ig_max-ig_min+1, &
               isign_max-isign_min+1 /)
          select case (layout)
          case ('yxels')
             dim_pdf(5:6) = (/ negrid, nlambda /)
             offset(5:6) = (/ ie_min-1, il_min-1 /)
             count (5:6) = (/ ie_max, il_max /) - offset(5:6)
          case default
             dim_pdf(5:6) = (/ nlambda, negrid /)
             offset(5:6) = (/ il_min-1, ie_min-1 /)
             count (5:6) = (/ il_max, ie_max /) - offset(5:6)
          end select
          dim_pdf(7) = nspec
          offset(7) = ispec_min - 1
          count (7) = ispec_max - offset(7)

          if ( g_size /= product(count(2:7)) ) print*, 'what to do?'

       end if

    end if

    ! define hyperslab selection for memory
    call h5screate_simple_f (rank_mem, dim_g, mem_pdf, error)
    call h5sselect_hyperslab_f (mem_pdf, H5S_SELECT_SET_F, &
         offset_mem, count_mem, error)
    if (error /= 0) print*, 'error in select memory hyperslab for pdf'

    ! define hyperslab selection for file
    call h5screate_simple_f (rank_pdf, dim_pdf, dsp_pdf, error)
    call h5sselect_hyperslab_f (dsp_pdf, H5S_SELECT_SET_F, & 
         offset, count, error)
    if (error /= 0) print*, 'error in select file hyperslab for pdf'

! TT> dim_g is refered in h5dwrite_f later
!    deallocate (dim_pdf, dim_g, offset, offset_mem, count, count_mem)
    deallocate (dim_pdf, offset, offset_mem, count, count_mem)
! <TT

    ! Set up parallel open-file property
    call h5pcreate_f (H5P_FILE_ACCESS_F, popen_id, error)
    if (proc0 .and. error<0) print *, 'error creating popen id'
    call h5pset_fapl_mpio_f (popen_id, communicator, MPI_INFO_NULL, error)
    if (proc0 .and. error<0) print *, 'error creating popen'

    ! Set up parallel write-to-file property
    call h5pcreate_f (H5P_DATASET_XFER_F, pwrite_id, error)
    if (proc0 .and. error<0) print *, 'error creating pwrite id'
    call h5pset_dxpl_mpio_f (pwrite_id, H5FD_MPIO_INDEPENDENT_F, error)
    if (proc0 .and. error<0) print *, 'error creating pwrite'

    ! proc0 writes grids
    if (proc0) call hdf_write_grids (error)

    initialized = .true.
# else
    if (error >= 0) then
       if (proc0) write (error_unit(),*) &
            'WARNING: hdf_init_dist is called without hdf5 compilation'
       error = -1
    end if
# endif
  end subroutine hdf_init_dist

  subroutine hdf_finish_dist !(error)
# ifdef HDF
# if FCOMPILER != _GFORTRAN_
    use hdf5, only: h5pclose_f, h5sclose_f
# endif
# endif
    implicit none
    integer :: error
# ifdef HDF
    call h5pclose_f (popen_id, error)
    call h5pclose_f (pwrite_id, error)
    call h5sclose_f (dsp_pdf, error)
    call h5sclose_f (mem_pdf, error)
# endif
  end subroutine hdf_finish_dist

  !
  ! write grids
  !
  subroutine hdf_write_grids (error)
# ifdef HDF
    use file_utils, only: error_unit
    use mp, only: proc0
    use file_utils, only: run_name
    use kt_grids, only: nx, ny, ntheta0, naky, box, akx, aky
    use kt_grids_box, only: x0, y0
    use theta_grid, only: ntgrid!, z0
    use le_grids, only: nlambda, negrid, e, al, w, wl
    use species, only: nspec
    use gs2_layouts, only: layout
!    use gs2_transforms, only: init_transforms
# if FCOMPILER != _GFORTRAN_
    use hdf5, only: h5fcreate_f, h5fclose_f, h5gcreate_f, h5gclose_f
    use hdf5, only: H5F_ACC_TRUNC_F
# endif
    use hdf_wrapper, only: hdf_write, hdf_error
# endif
    implicit none
    integer, intent (out) :: error 
# ifdef HDF

    ! create file
    write (fname_grids, '(a)') trim(run_name) // trim(grids_filepref) // '.hdf'
    call h5fcreate_f (fname_grids, H5F_ACC_TRUNC_F, grids_file_id, error)
    if (error < 0) call hdf_error (file=fname_grids)

    ! write layout and accel flag
    if (error >= 0) call hdf_write (grids_file_id, dname_layout, layout, error)
    if (error >= 0 .and. .not.export_complex) then
!!$       ! initialize FFT and accel flag
!!$       call init_transforms &
!!$            (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)
       call hdf_write (grids_file_id, dname_accel, accelerated, error)
    end if

    ! create parameter group
    if (error >= 0) then
       call h5gcreate_f (grids_file_id, gname_params, grp_id, error)
       if (error < 0) call hdf_error (grp=gname_params)
    end if

    ! write number of grid points
    if (error >= 0) call hdf_write (grp_id, dname_nx, nx, error)
    if (error >= 0) call hdf_write (grp_id, dname_ny, ny, error)
    if (error >= 0) call hdf_write (grp_id, dname_nz, 2*ntgrid+1, error)
    if (error >= 0) call hdf_write (grp_id, dname_nlambda, nlambda, error)
    if (error >= 0) call hdf_write (grp_id, dname_negrid, negrid, error)
    if (error >= 0) call hdf_write (grp_id, dname_nspec, nspec, error)
    if (error >= 0) call hdf_write (grp_id, dname_ntheta0, ntheta0, error)
    if (error >= 0) call hdf_write (grp_id, dname_naky, naky, error)
    if (error >= 0) call hdf_write (grp_id, dname_ntgrid, ntgrid, error)

    ! write spatial length if box
    if (box) then
       if (error >= 0) call hdf_write (grp_id, dname_x0, x0, error)
       if (error >= 0) call hdf_write (grp_id, dname_y0, y0, error)
    end if
!    if (error >= 0) call hdf_write (grp_id, dname_z0, z0, error)

    ! close parameter group
    if (error >= 0) then
       call h5gclose_f (grp_id, error)
       if (error < 0) &
            call hdf_error (message='error closing group in write_grids')
    end if

    ! create grids group
    if (error >= 0) then
       call h5gcreate_f (grids_file_id, gname_grids, grp_id, error)
       if (error < 0) call hdf_error (grp=gname_grids)
    end if

    ! write velocity-space grid points and wave number if not box
    if (error >= 0) call hdf_write (grp_id, dname_energy, e, error)
    if (error >= 0) call hdf_write (grp_id, dname_al, al, error)
    if (error >= 0) call hdf_write (grp_id, dname_w, w, error)
    if (error >= 0) call hdf_write (grp_id, dname_wl, wl, error)
    if (.not.box) then
       if (error >= 0) call hdf_write (grp_id, dname_kx, akx(1), error)
       if (error >= 0) call hdf_write (grp_id, dname_ky, aky(1), error)
    end if

    ! close group and file
    call h5gclose_f (grp_id, error)
    if (error < 0) &
         call hdf_error (message='error closing group in write_grids')

    call h5fclose_f (grids_file_id, error)
    if (error < 0) &
         call hdf_error (message='error finishing grids_hdf')
# else
    error = -1
# endif
  end subroutine hdf_write_grids

end module gs2_dist_io

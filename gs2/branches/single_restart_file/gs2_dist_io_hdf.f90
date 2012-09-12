module hdfparams

!!$ Defines names, params for the HDF5 datasets, files etc.
!!$ ==George==

  use hdf5
  Implicit None

  CHARACTER(LEN=20), PARAMETER :: dsetname_pdf = "Density Function" 
  CHARACTER(LEN=20), PARAMETER :: dsetname_energy = "Energy grid" 
  CHARACTER(LEN=20), PARAMETER :: filepref = "_density_gs2"  
  CHARACTER(LEN=32) :: filename
  INTEGER(HID_T) :: file_id    !! File identifier
  INTEGER(HID_T) :: dsp_pdf_id  !! Data space id
  INTEGER(HID_T) :: pdf_id, mem_pdf_id  !! Dataset id; memory id for current_datum
  INTEGER(HID_T) :: dsp_energy_id  !! Data space id
  INTEGER(HID_T) :: energy_id
  INTEGER, parameter ::   rank_pdf= 7      ! Dataset rank
  INTEGER, parameter ::   rank_mem= 1      ! curent_datum rank
  INTEGER(HSIZE_T), DIMENSION(rank_pdf)  :: dim_pdf  ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(rank_mem) :: count_mem= (/1/)
  INTEGER(HID_T) :: plist_id    !! property list id
  INTEGER(HSSIZE_T) :: offset(rank_pdf)

  integer, parameter :: rank1 = 1
  integer (hsize_t), dimension(rank1) :: dim1

End Module hdfparams

module gs2_dist_io

  implicit none

  public :: write_dist

contains
  
  subroutine write_dist (g0)

    use mp, only: proc0
    use gs2_transforms, only: init_transforms, transform2
    use gs2_layouts, only: yxf_lo, accelx_lo, g_lo
    use gs2_layouts, only: is_idx, ie_idx, il_idx, ik_idx, it_idx, isign_idx, ig_idx
    use kt_grids, only: naky, ntheta0, nx, ny
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use hdf5
    use hdfparams

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0   

    real, dimension (:,:,:), allocatable :: agx
    real, dimension (:,:),   allocatable ::  gx

    integer :: iglo, ik, it, is, ig, il, ie, ia, isgn, i, j, k
    logical :: alloc = .true., accelerated
    integer :: error

    integer(HSIZE_T), parameter, dimension(7) :: count = (/1, 1, 1, 1, 1, 1, 1/) 
    integer(HSIZE_T), parameter, dimension(1) :: dim_datum = (/1/)
    real, dimension(1) ::  current_datum

! For now, assume that we have a box domain and have initialized the nonlinear 
! (FFT) arrays.  

! Need to initialize output with HDF; suggest call from here
! GEORGE: NEED TO FILL IN THIS PROCEDURE
!!$ Done; needs testing

    call init_hdf(error)
    ! call write_grids
    if (error /= 0) then
       print *, "error with initializing hdf5, exiting..."
       return
    end if

! Initialize accelerated variable
! DONE

    call init_transforms (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)

! allocate space for transformed data
! DONE
    
    if (alloc) then
       if (accelerated) then
          allocate (agx(-ntgrid:ntgrid, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
          agx = 0.
       else
          allocate ( gx(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
          gx = 0.
       end if
       alloc = .false.
    end if
   
! transform data to real space
! ia is only used to get the compiler to recognize need for accelerated transforms
! DONE

    if (accelerated) then
       call transform2 (g0, agx, ia)
    else
       call transform2 (g0, gx)
    end if

! Format of data at this point depends on whether the accelerated transforms were used 

! write data to file    
! GEORGE: NEED TO FILL IN THESE ELEMENTS
!!$ Done; not tested ==George==

!!$ Set up parallel write-to-file property
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)

    CALL h5screate_simple_f(rank_mem, dim_datum, mem_pdf_id, error)
    CALL h5sselect_hyperslab_f(mem_pdf_id, H5S_SELECT_SET_F, &
         offset, count_mem, error) 

    if (accelerated) then
!!$       do k = accelx+lo%llim_proc, accelx_lo%ulim_proc
       !!$ probably typo here; ==George==
       do k = accelx_lo%llim_proc, accelx_lo%ulim_proc
          do isgn = 1, 2     ! sigma index
             offset(4) = isgn - 1
!!$ ig takes negative values : is that OK with HDF? need to look at this
             do ig = -ntgrid, ntgrid ! z (a.k.a. theta) index

! Index functions for accelerated layouts are not written!  (wow)
! Format will be 

                is = is_idx(accelx_lo, k)   ! species index
                ie = ie_idx(accelx_lo, k)   ! energy index
                il = il_idx(accelx_lo, k)   ! lambda index
                ik = ik_idx(accelx_lo, k)   ! y index 
                it = it_idx(accelx_lo, k)   ! x index 

                offset(1) = it - 1 
                offset(2) = ik - 1 
                offset(3) = ig + ntgrid
                offset(5) = ie - 1
                offset(6) = il - 1
                offset(7) = is - 1 

              ! write agx (ig, isgn, k)  [HDF FUNCTION NEEDED HERE]
                current_datum(1) = agx (ig, isgn, k)

                CALL h5sselect_hyperslab_f(dsp_pdf_id, H5S_SELECT_SET_F, & 
                     offset, count, error)

                CALL H5dwrite_f(pdf_id, H5T_NATIVE_DOUBLE, current_datum, &
                     dim_datum, error, file_space_id = dsp_pdf_id, mem_space_id = mem_pdf_id, & 
                     xfer_prp=plist_id)
             end do
          end do
       end do
    else
       do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
          do i = 1, yxf_lo%ny

! Index functions for real space layouts are not written!  (wow)
! Format will be 

                is = is_idx(yxf_lo, j)   ! species index
                ie = ie_idx(yxf_lo, j)   ! energy index
                il = il_idx(yxf_lo, j)   ! lambda index
                isgn = isign_idx(yxf_lo, j) ! sigma index
                ig = ig_idx(yxf_lo, j)   ! theta (z) index
                it = it_idx(yxf_lo, j)   ! x index 

             ! write gx (i, j)   [HDF FUNCTION NEEDED HERE]

                offset(1) = it - 1
                offset(2) = i - 1        !!$ presumably y index
                offset(3) = ig + ntgrid
                offset(4) = isgn - 1
                offset(5) = ie - 1
                offset(6) = il - 1
                offset(7) = is - 1

                current_datum(1) = gx (i, j)
                CALL h5sselect_hyperslab_f(dsp_pdf_id, H5S_SELECT_SET_F, & 
                     offset, count, error)
                CALL H5dwrite_f(pdf_id, H5T_NATIVE_DOUBLE, current_datum, &
                     dim_datum, error, file_space_id = dsp_pdf_id, mem_space_id = mem_pdf_id, & 
                     xfer_prp=plist_id)

          end do
       end do
    end if

    call finish_hdf (error)

  end subroutine write_dist

  subroutine write_grids (error)

    use mp, only: proc0
    use file_utils, only: run_name
    use kt_grids, only: nx, ny
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid, e
    use species, only: nspec
    use hdf5
    use hdfparams

    Implicit None

    logical :: initialized == .false.

    if (initialized) return
    initialized = .true.

    if (proc0) then
       dim1(1) = negrid
       call h5dwrite_f(energy_id, H5T_NATIVE_DOUBLE, e(:,1), dim1, error)
    

       

    end if


  end subroutine write_grids

  subroutine init_hdf(error)

    use mp
    use file_utils, only: run_name
    use kt_grids, only: nx, ny
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use hdf5
    use hdfparams

    Implicit None

    INTEGER, intent(out) :: error

    INTEGER :: info, comm

! initialize output with HDF library from here.  
! GEORGE: NEED TO FILL IN THIS PROCEDURE
!!$ Done; still not tested ==George==
    
    info = MPIINFONULL
    comm = communicator
    dim_pdf(1) = nx
    dim_pdf(2) = ny
    dim_pdf(3) = 2*ntgrid + 1
    dim_pdf(4) = 2 !!$ velocity orientations
    dim_pdf(5) = negrid
    dim_pdf(6) = nlambda
    dim_pdf(7) = nspec
    
    write(filename,fmt="(a)") trim(run_name)//trim(filepref)//'.hdf'
 
    CALL h5open_f(error)
    
    if (error /= 0) then 
       return
    end if
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, info,  error)
!!$ Open file for parallel writing    
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, & 
         access_prp = plist_id)
    CALL h5pclose_f(plist_id, error)
!!$ Create the PDF dataspace
    CALL h5screate_simple_f(rank_pdf, dim_pdf, dsp_pdf_id, error)
!!$ Create the energy grid dataspace
    dim1(1) = negrid
    CALL h5screate_simple_f(rank1, dim1, dsp_energy_id, error)
!!$ Create the PDF dataset
    CALL h5dcreate_f(file_id, dsetname_pdf , H5T_NATIVE_DOUBLE, dsp_pdf_id, &
         pdf_id, error)
!!$ Create the energy grid dataset
    CALL h5dcreate_f(file_id, dsetname_energy , H5T_NATIVE_DOUBLE, dsp_energy_id, &
         energy_id, error)

  end subroutine init_hdf

  subroutine finish_hdf(error)
    use hdf5
    use hdfparams

    Implicit None

    INTEGER, intent(out) :: error
    CALL h5sclose_f(dsp_pdf_id, error)
    CALL h5dclose_f(pdf_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

  end subroutine finish_hdf
  
end module gs2_dist_io

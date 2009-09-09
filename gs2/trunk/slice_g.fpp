# include "define.inc"

program slice_g

  !
  ! ASCII file generator for the slice of dist func hdf file
  !   written by Tomo Tatsuno (Apr 28 2008)
  ! USAGE: slice_g <hdf5 dist filename> [it] [ik] [ig] [is]
  !   [it], [ik], [ig], [is] are optional integers for slice taken.
  !   If not specified, default (1,1,0,1) are used.
  !   Slice is written out to standard output.
  !
  ! To do: support multi-species case --- done (6/2/08)
  !        support h output with a proper g -> h conversion
  !                 --- done with hk read (7/23/08)
  !        support non-square array for real space output
  !

# if FCOMPILER == _GFORTRAN_
  use hdf5
# else
  use hdf5, only: HID_T, HSIZE_T
#endif
  use hdf_wrapper, only: hdf_finish

  implicit none
  include 'hdf_params.inc'

  logical :: test=.false.           ! test flag
  logical :: write_vplane = .true.  ! write out vpar-vperp plane

! TT: below may be used from other agk f90 files?
  character (5) :: layout
  integer :: nx, ny, nz, ntheta0, naky, ntgrid, negrid, nlambda, nspec
! TT: above may be used from other agk f90 files?
  logical :: import_complex, import_h, accelerated
  integer :: ik_out=1, it_out=1, ig_out=0, is_out=1
  integer :: ie, il
  real :: time, vpar, vperp
  real, dimension (:), allocatable :: al
  real, dimension (:,:), allocatable :: e
  real, dimension (:,:,:), allocatable :: gx
  real, dimension (:,:,:,:), allocatable :: gkr
  complex, dimension (:,:,:), allocatable :: gk

  ! check command line args etc
  call check_command
  if (test) print *, 'check command done'

  ! check dist and grids file status
  call check_file
  if (test) print *, 'file check done'

  ! read grids data
  call read_grid_data

  if (import_complex) then
     write (0,*) 'ntheta0 = ', ntheta0, 'naky = ', naky, 'ntgrid = ', ntgrid
  else
     write (0,*) 'nx = ', nx, 'ny = ', ny, 'ntgrid = ', ntgrid
  end if

  ! read dist func data in slice
  call read_dist_data

  ! write to standard output
  if (import_complex) then
     write (*,'(2(a,i5,3x))') '# it = ', it_out, 'ik = ', ik_out
     if (write_vplane) then
        if (import_h) then
           write (*,'(a)') '#  vpar    vperp    Re(hk)   Im(hk)'
        else
           write (*,'(a)') '#  vpar    vperp    Re(gk)   Im(gk)'
        end if
     else
        if (import_h) then
           write (*,'(a)') '#  al    e    Re(hk(1))   Im(hk(1)) &
                & Re(hk(2))   Im(hk(2))'
        else
           write (*,'(a)') '#  al    e    Re(gk(1))   Im(gk(1)) &
                & Re(gk(2))   Im(gk(2))'
        end if
     end if
  else
     write (*,'(2(a,i5))') '# ik = ', ik_out, ' it = ', it_out
     if (write_vplane) then
        write (*,'(a)') '#  vpar    vperp    gx'
     else
        write (*,'(a)') '#  al    e    gx(1)    gx(2)'
     end if
  end if
  if (write_vplane) then
     do ie=1, negrid
        do il=1, nlambda
           vpar = sqrt( (1.0 - al(il)) * e(ie,1) )
           vperp = sqrt( al(il) * e(ie,1) )
           if (import_complex) then
              if (layout == 'yxels') then
                 write (*,'(4es15.5)') vpar, vperp, gk(1,ie,il)
              else
                 write (*,'(4es15.5)') vpar, vperp, gk(1,il,ie)
              end if
           else
              if (layout == 'yxels') then
                 write (*,'(3es15.5)') vpar, vperp, gx(1,ie,il)
              else
                 write (*,'(3es15.5)') vpar, vperp, gx(1,il,ie)
              end if
           end if
        end do
        do il=nlambda, 1, -1
           vpar = -sqrt( (1.0 - al(il)) * e(ie,1) )
           vperp = sqrt( al(il) * e(ie,1) )
           if (import_complex) then
              if (layout == 'yxels') then
                 write (*,'(4es15.5)') vpar, vperp, gk(2,ie,il)
              else
                 write (*,'(4es15.5)') vpar, vperp, gk(2,il,ie)
              end if
           else
              if (layout == 'yxels') then
                 write (*,'(3es15.5)') vpar, vperp, gx(2,ie,il)
              else
                 write (*,'(3es15.5)') vpar, vperp, gx(2,il,ie)
              end if
           end if
        end do
        print *
     end do
  else
     do il=1, nlambda
        do ie=1, negrid
           if (import_complex) then
              if (layout == 'yxels') then
                 write(*,'(6es15.5)') al(il), e(ie,1), gk(1,ie,il), gk(2,ie,il)
              else
                 write(*,'(6es15.5)') al(il), e(ie,1), gk(1,il,ie), gk(2,il,ie)
              end if
           else
              if (layout == 'yxels') then
                 write(*,'(4es15.5)') al(il), e(ie,1), gx(1,ie,il), gx(2,ie,il)
              else
                 write(*,'(4es15.5)') al(il), e(ie,1), gx(1,il,ie), gx(2,il,ie)
              end if
           end if
        end do
        print*
     end do
  end if

  deallocate (offset, count, stride, block)
  deallocate (offset_mem, count_mem)
  deallocate (dim_g)
  if (allocated(gx)) deallocate (gx)
  if (allocated(gk)) deallocate (gk)
  deallocate (al, e)

  call hdf_finish

contains

  subroutine check_command

# if FCOMPILER == _GFORTRAN_
    use command_line, only: cl_getarg
# else
    use command_line, only: iargc, cl_getarg
# endif
    logical :: flag
    character (100) :: arg
    integer :: no_arg, ier

    no_arg = iargc()
    if ( (no_arg < 1) .or. (no_arg > 5) ) then
       write (0,*) 'USAGE: slice_g <hdf5 dist filename> [it] [ik] [ig] [is]'
       stop
    end if
    if (test) print *, 'no_arg is ', no_arg

    ! check dist file exist
    call cl_getarg (1, fname_dist, 0, ier)
    inquire (file=fname_dist, exist=flag)
    if (.not.flag) then
       write (0,*) 'ERROR: File not found: ', trim(fname_dist)
       stop
    end if
    if (test) print *, trim(fname_dist), flag

    ! set up i{t,k,g}_out
    if (no_arg > 1) then
       call cl_getarg (2, arg, 0, ier)
       read (arg,*) it_out
       if (test) print *, 'it_out = ', it_out
    end if

    if (no_arg > 2) then
       call cl_getarg (3, arg, 0, ier)
       read (arg,*) ik_out
       if (test) print *, 'ik_out = ', ik_out
    end if

    if (no_arg > 3) then
       call cl_getarg (4, arg, 0, ier)
       read (arg,*) ig_out
    end if

    if (no_arg > 4) then
       call cl_getarg (5, arg, 0, ier)
       read (arg,*) is_out
    end if

  end subroutine check_command

  subroutine check_file

# if FCOMPILER != _GFORTRAN_
    use hdf5, only: h5fis_hdf5_f
# endif
    use hdf_wrapper, only: hdf_init
    use command_line, only: cl_getarg
    use file_utils, only: init_file_utils, error_unit, run_name
    logical :: flag
    integer :: ier

    ! setup run_name and error_unit
    call init_file_utils (flag, input=.false., error=.false., &
         name=trim(fname_dist(1:index(fname_dist, dist_filepref)-1)))
    if (test) print *, 'run name is ', trim(run_name)
    if (test) print *, 'error_unit is ', error_unit()

    ! check grid file exist
    write (fname_grids, '(a)') trim(run_name) // trim(grids_filepref) // '.hdf'
    inquire (file=fname_grids, exist=flag)
    if (.not.flag) then
       write (0,*) 'ERROR: File not found: ', trim(fname_grids)
       stop
    end if

    call hdf_init (stop=.true.)
    if (test) print *, 'hdf_init done'

    call h5fis_hdf5_f (fname_dist, flag, ier)
    if (.not.flag) then
       write (0,*) 'ERROR: File not hdf5 format: ', trim(fname_dist)
       stop
    end if

    call h5fis_hdf5_f (fname_grids, flag, ier)
    if (.not.flag) then
       write (0,*) 'ERROR: File not hdf5 format: ', trim(fname_grids)
       stop
    end if

  end subroutine check_file

  subroutine read_grid_data

# if FCOMPILER != _GFORTRAN_
    use hdf5, only: H5F_ACC_RDONLY_F, h5fopen_f, h5fclose_f
    use hdf5, only: h5gn_members_f, h5gget_obj_info_idx_f, h5gopen_f, h5gclose_f
# endif
    use hdf_wrapper, only: hdf_read, hdf_error
    integer (HID_T) :: otp
    character (20) :: name
    logical :: found=.false., box
    integer :: nmem, imem, ier

    ! open grid file
    call h5fopen_f (fname_grids, H5F_ACC_RDONLY_F, grids_file_id, ier)
    if (ier /= 0) call hdf_error (file=fname_grids)
    if (test) print *, 'opened grid file'

    ! read layout
    call hdf_read (grids_file_id, dname_layout, layout)

    ! inquire if accel flag exists
    call h5gn_members_f (grids_file_id, '/', nmem, ier)
    do imem=1, nmem
       call h5gget_obj_info_idx_f (grids_file_id, '/', imem-1, name, otp, ier)
       if (name == dname_accel) then
          found = .true.
          if (test) print *, 'Yes, it is! ', trim(name), otp
          exit
       end if
       if (test) print *, 'No, it isn''t. ', imem, trim(name), otp
    end do
    ! existence of accel means real space output
    if (found) then
       call hdf_read (grids_file_id, dname_accel, accelerated)
    else
       if (test) print *, 'No accel flag found.'
    end if
    import_complex = .not.found
    if (test) print *, 'accel = ', accelerated
    if (test) print *, 'import_complex = ', import_complex

    ! open parameter group and inquire if x0 exists
    call h5gopen_f (grids_file_id, gname_params, grp_id, ier)
    if (ier /= 0) call hdf_error (grp=gname_params)
    call h5gn_members_f (grp_id, gname_params, nmem, ier)
    do imem=1, nmem
       call h5gget_obj_info_idx_f &
            (grids_file_id, gname_params, imem-1, name, otp, ier)
       if (name == dname_x0) then
          found = .true.
          if (test) print *, 'Yes, it is! ', trim(name), otp
          exit
       end if
       if (test) print *, 'No, it isn''t. ', imem, trim(name), otp
    end do
    box = found
    if (test) print *, 'box = ', box

    ! read parameter data
    call hdf_read (grp_id, dname_nx, nx)
    call hdf_read (grp_id, dname_ny, ny)
    call hdf_read (grp_id, dname_nz, nz)
    call hdf_read (grp_id, dname_ntheta0, ntheta0)
    call hdf_read (grp_id, dname_naky, naky)
    call hdf_read (grp_id, dname_ntgrid, ntgrid)
    call hdf_read (grp_id, dname_negrid, negrid)
    call hdf_read (grp_id, dname_nlambda, nlambda)
    call hdf_read (grp_id, dname_nspec, nspec)
    call h5gclose_f (grp_id, ier)

    ! open grid group and read data
    call h5gopen_f (grids_file_id, gname_grids, grp_id, ier)
    if (ier /= 0) call hdf_error (grp=gname_grids)
    allocate (al(nlambda), e(negrid,nspec))
    call hdf_read (grp_id, dname_al, al)
    call hdf_read (grp_id, dname_energy, e)
    call h5gclose_f (grp_id, ier)

    ! close grid file
    call h5fclose_f (grids_file_id, ier)

  end subroutine read_grid_data

  subroutine read_dist_data

# if FCOMPILER != _GFORTRAN_
    use hdf5, only: H5F_ACC_RDONLY_F, h5fopen_f, h5fclose_f
    use hdf5, only: h5dopen_f, h5dget_space_f, h5dread_f, h5dclose_f
    use hdf5, only: H5S_SELECT_SET_F, h5screate_simple_f
    use hdf5, only: h5sselect_hyperslab_f, h5sclose_f
    use hdf5, only: h5sget_simple_extent_ndims_f
    use hdf5, only: h5gn_members_f, h5gget_obj_info_idx_f
# endif
    use hdf_wrapper, only: hdf_mem_real, hdf_read, hdf_error
    use convert, only: r2c
    integer :: ier
    logical :: square
    integer (HID_T) :: otp
    character (20) :: name
    logical :: found=.false.
    integer :: nmem, imem

    ! open distribution function and read time
    call h5fopen_f (fname_dist, H5F_ACC_RDONLY_F, dist_file_id, ier)
    if (ier < 0) call hdf_error (file=fname_dist)
    call hdf_read (dist_file_id, dname_time, time)
    write (*,'(a,es15.5)') '# time = ', time

    ! define file dist func selection
    if (import_complex) then
       ! inquire g or h output
       call h5gn_members_f (dist_file_id, '/', nmem, ier)
       do imem=1, nmem
          call h5gget_obj_info_idx_f (dist_file_id, '/', imem-1, name, otp, ier)
          if (name == dname_hk) then
             found = .true.
             if (test) print *, 'Yes, it is! ', trim(name), otp
             exit
          end if
          if (test) print *, 'No, it isn''t. ', imem, trim(name), otp
       end do
       import_h = found
       dname_pdf = dname_gk
       if (import_h) dname_pdf = dname_hk
    else
       dname_pdf = dname_gx
    end if
    ! open dist func dataset
    call h5dopen_f (dist_file_id, dname_pdf, dst_pdf, ier)
    if (ier < 0) call hdf_error (dset=dname_pdf)
    call h5dget_space_f (dst_pdf, dsp_pdf, ier)
    if (import_complex) then
       call h5sget_simple_extent_ndims_f (dsp_pdf, rank_pdf, ier)
       if (test) print *, 'rank_pdf: ', rank_pdf
       if (ier < 0) call hdf_error (message='failed obtaining rank_pdf')
       square = rank_pdf==8
    else
       rank_pdf = 7
    end if

    allocate (offset(rank_pdf), count(rank_pdf))
    allocate (stride(rank_pdf), block(rank_pdf))
    offset = 0
    count = 1
    stride = 1
    block = 1

    if (import_complex) then
       ! (/ ri, theta, sign, g_lo /)
       offset(2) = ig_out+ntgrid
       count(1) = 2
       count(3) = 2
       if (square) then
          select case (layout)
          case ('yxels')
             offset(4:5) = (/ ik_out-1, it_out-1 /)
             count(6:7) = (/ negrid, nlambda /)
          case ('yxles')
             offset(4:5) = (/ ik_out-1, it_out-1 /)
             count(6:7) = (/ nlambda, negrid /)
          case ('lexys')
             offset(6:7) = (/ it_out-1, ik_out-1 /)
             count(4:5) = (/ nlambda, negrid /)
          case ('lxyes')
             offset(5:6) = (/ it_out-1, ik_out-1 /)
             count(4) = nlambda
             count(7) = negrid
          case ('lyxes')
             offset(5:6) = (/ ik_out-1, it_out-1 /)
             count(4) = nlambda
             count(7) = negrid
          end select
          offset(8) = is_out-1
       else ! if not square
          select case (layout)
          case ('yxels', 'yxles')
             offset (4) = ik_out-1 + naky * (it_out-1)
             stride (4) = naky*ntheta0
             count (4) = nlambda*negrid
          case ('lexys')
             offset (4) = nlambda * negrid * (it_out-1 + ntheta0 * (ik_out-1))
             stride (4) = 1
             count (4) = nlambda*negrid
          case ('lxyes')
             offset (4) = nlambda * (it_out-1 + ntheta0 * (ik_out-1))
             stride (4) = nlambda * ntheta0 * naky
             block (4) = nlambda
             count (4) = negrid
          case ('lyxes')
             offset (4) = nlambda * (ik_out-1 + naky * (it_out-1))
             stride (4) = nlambda * naky * ntheta0
             block (4) = nlambda
             count (4) = negrid
          end select
          offset(4) = offset(4) + naky*ntheta0*nlambda*negrid*(is_out-1)
          if (test) then
             print *, 'layout: ', trim(layout)
             print *, 'offset: ', offset
             print *, 'stride: ', stride
             print *, 'count: ', count
             print *, 'block: ', block
          end if
       end if
    else ! if not import_complex
       if (accelerated) then
          ! accel_lo case follows
          offset(1) = ig_out + ntgrid
          offset(3) = ik_out - 1
          offset(4) = it_out - 1
       else
          ! The following is for non-accel layout
          offset(1:3) = (/ ik_out-1, it_out-1, ig_out+ntgrid /)
       end if
       count(2) = 2
       select case (layout)
       case ('yxels')
          count(5:6) = (/ negrid, nlambda /)
       case default
          count(5:6) = (/ nlambda, negrid /)
       end select
       offset(7) = is_out - 1
    end if

    ! select file hyperslab
    call h5sselect_hyperslab_f &
         (dsp_pdf, H5S_SELECT_SET_F, offset, count, ier, stride, block)
    if (ier < 0) call hdf_error (message='error in select file hyperslab')

    ! define memory hyperslab
    ! This is needed since dimension is different from file
    if (import_complex) then
       select case (layout)
       case ('yxels')
          allocate (gk(2, negrid, nlambda), gkr(2, 2, negrid, nlambda))
       case default
          allocate (gk(2, nlambda, negrid), gkr(2, 2, nlambda, negrid))
       end select
       rank_mem = size(shape(gkr))
       allocate (dim_g(rank_mem), offset_mem(rank_mem), count_mem(rank_mem))
       dim_g = shape(gkr)
    else
       select case (layout)
       case ('yxels')
          allocate (gx(2, negrid, nlambda))
       case default
          allocate (gx(2, nlambda, negrid))
       end select
       rank_mem = size(shape(gx))
       allocate (dim_g(rank_mem), offset_mem(rank_mem), count_mem(rank_mem))
       dim_g = shape(gx)
    end if
    if (test) print *, 'dim_g: ', dim_g
    call h5screate_simple_f (rank_mem, dim_g, mem_pdf, ier)
    offset_mem = 0
    count_mem = dim_g
    call h5sselect_hyperslab_f (mem_pdf, H5S_SELECT_SET_F, offset_mem, &
         count_mem, ier)

    ! check if number of elements agrees
    if (product(count*block) /= product(dim_g)) then
       write (0,*) 'ERROR: array mismatch.'
       write (0,*) 'dim_g: ', dim_g
       write (0,*) 'count: ', count
       write (0,*) 'block: ', block
       goto 1
    end if

    ! read !
    if (import_complex) then
       call h5dread_f (dst_pdf, hdf_mem_real, gkr, dim_g, ier, &
            mem_space_id=mem_pdf, file_space_id=dsp_pdf)
       call r2c (gk, gkr)
       deallocate (gkr)
    else
       call h5dread_f (dst_pdf, hdf_mem_real, gx, dim_g, ier, &
            mem_space_id=mem_pdf, file_space_id=dsp_pdf)
    end if
    if (ier < 0) call hdf_error (dset=dname_pdf)

1   call h5sclose_f (mem_pdf, ier)
    call h5sclose_f (dsp_pdf, ier)
    call h5dclose_f (dst_pdf, ier)
    call h5fclose_f (dist_file_id, ier)

  end subroutine read_dist_data

end program slice_g

!> This module contains functions for writing quantities
!! that have the dimensions of the fields, including
!! volume averaged quantities, etc.
module diagnostics_fields
   use simpledataio
   use simpledataio_write
   use theta_grid, only: ntgrid
   use kt_grids, only: naky, ntheta0
   use diagnostics_config, only: diagnostics_type
   implicit none
   !integer, parameter :: gnostics%rtype = SDATIO_DOUBLE

contains

  subroutine write_fields(gnostics)
    use fields_arrays, only: phinew, aparnew, bparnew
    use run_parameters, only: fphi, fapar, fbpar
    type(diagnostics_type), intent(inout) :: gnostics

    
    if (fphi >epsilon(0.0)) call write_standard_field_properties(gnostics, &
      'phi',  'The electrostatic potential', 'T_r/e', phinew, gnostics%distributed)
    if (fapar>epsilon(0.0)) call write_standard_field_properties(gnostics, &
      'apar', 'The parallel magnetic potential', '...', aparnew, gnostics%distributed)
    if (fbpar>epsilon(0.0)) call write_standard_field_properties(gnostics, &
      'bpar', 'The parallel magnetic potential', '...', bparnew, gnostics%distributed)
  end subroutine write_fields

  subroutine write_standard_field_properties(gnostics, field_name, field_description, &
    field_units, field_value, distributed)
    use diagnostics_create_and_write, only: create_and_write_variable
    use volume_averages
    type(diagnostics_type), intent(inout) :: gnostics
    character(*), intent(in) :: field_name, field_description, field_units
    complex, dimension(:,:,:), intent(in) :: field_value
    logical, intent(in) :: distributed
    real, dimension(ntheta0, naky) :: field2_by_mode
    real, dimension(naky) :: field2_by_ky
    real, dimension(ntheta0) :: field2_by_kx

    call average_theta(field_value, field_value, field2_by_mode, distributed)
    !call create_and_write_field(gnostics%sfile, field_name, field_description, field_units, field_value)
    call create_and_write_field_by_mode(gnostics, field_name, field_description, field_units, &
      field_value, field2_by_mode, distributed)

    call average_kx(field2_by_mode, field2_by_ky, distributed)
    call create_and_write_variable(gnostics, gnostics%rtype, field_name//"2_by_ky", "Yt", &
      field_description//" squared and averaged over theta and kx, as a function of time", &
      "("//field_units//")^2", field2_by_ky)

    call average_ky(field2_by_mode, field2_by_kx, distributed)
    call create_and_write_variable(gnostics, gnostics%rtype, field_name//"2_by_kx", "Xt", &
      field_description//" squared and averaged over theta and ky, as a function of time", &
      "("//field_units//")^2", field2_by_kx)

    call create_and_write_variable(gnostics, gnostics%rtype, field_name//"2", "t", &
      field_description//" squared and averaged over theta, kx and ky, as a function of time", &
      "("//field_units//")^2", sum(field2_by_kx))

    if (field_name .eq. 'phi') gnostics%current_results%phi2 = sum(field2_by_kx)
    if (field_name .eq. 'apar') gnostics%current_results%apar2 = sum(field2_by_kx)
    if (field_name .eq. 'bpar') gnostics%current_results%bpar2 = sum(field2_by_kx)

  end subroutine write_standard_field_properties

  !> Dump the fields to a time-labelled ascii file:
  !! suggestion... use the netcdf file with 
  !! write_phi_over_time etc on instead :-)
  subroutine dump_fields_periodically(gnostics)
    use file_utils, only: get_unused_unit
    use kt_grids, only: naky, ntheta0, aky, theta0, akx
    use theta_grid, only: theta
    use fields_arrays, only: phinew, aparnew, bparnew
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real :: t
    character(200) :: filename
    integer :: ik, it, ig, unit
    integer :: ntg_out 

    ! EGH I'm pretty sure that there is no difference
    ! between ntg_out and ntgrid in the old diagnostics 
    ! module... anyone disagree?
    ntg_out = ntgrid

    t = gnostics%user_time

    call get_unused_unit (unit)
    write (filename, "('dump.fields.t=',e13.6)") t
    open (unit=unit, file=filename, status="unknown")
    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntg_out, ntg_out
             write (unit=unit, fmt="(20(1x,e12.5))") &
                  theta(ig), aky(ik), theta0(it,ik), &
                  phinew(ig,it,ik), aparnew(ig,it,ik), &
                  bparnew(ig,it,ik), t, akx(it)
          end do
          write (unit, "()")
       end do
    end do
    close (unit=unit)
  end subroutine dump_fields_periodically

  subroutine create_and_write_field_by_mode(gnostics, field_name, field_description, field_units, &
      val, field2_by_mode, distributed)
   use fields_parallelization, only: field_k_local
   type(diagnostics_type), intent(in) :: gnostics
   character(*), intent(in) :: field_name
   character(*), intent(in) :: field_description
   character(*), intent(in) :: field_units
   character(len=len(field_name)+9) :: field2_by_mode_name
   character(len=len(field_name)+2) :: field_t_name
   complex, intent(in), dimension(-ntgrid:ntgrid, ntheta0, naky) :: val
   real, dimension(ntheta0, naky), intent(in) :: field2_by_mode
   logical, intent(in) :: distributed
   complex, dimension(-ntgrid:ntgrid, 1, 1) :: dummyc
   real, dimension(1,1) :: dummyr
   logical :: write_field_by_time
   integer :: it,ik

   field2_by_mode_name =  field_name//"2_by_mode" 
   field_t_name = field_name//"_t"


   write_field_by_time = .false.
   if (field_name .eq. 'phi'  .and. gnostics%write_phi_over_time ) write_field_by_time = .true.
   if (field_name .eq. 'apar' .and. gnostics%write_apar_over_time) write_field_by_time = .true.
   if (field_name .eq. 'bpar' .and. gnostics%write_bpar_over_time) write_field_by_time = .true.
   
   if (gnostics%create) then 
     call create_variable(gnostics%sfile, gnostics%rtype, field_name, "rzXY", field_description, field_units)
     call create_variable(gnostics%sfile, gnostics%rtype, field2_by_mode_name, "XYt", &
       field_description//" squared and averaged over theta, as a function of kx and ky" , "("//field_units//")^2")
     if (write_field_by_time) then 
       call create_variable(gnostics%sfile, gnostics%rtype, field_t_name, "rzXYt", &
         field_description//" as a function of time", field_units)
     end if


   end if
   
   if (gnostics%create .or. .not. gnostics%wryte) return
   
   if (.not. distributed) then

     call write_variable(gnostics%sfile, field_name, val)
     call write_variable(gnostics%sfile, field2_by_mode_name, field2_by_mode)
     if (write_field_by_time) call write_variable(gnostics%sfile, field_t_name, val)

   else

     call set_count(gnostics%sfile, field_name, "X", 1)
     call set_count(gnostics%sfile, field_name, "Y", 1)
     call set_count(gnostics%sfile, field2_by_mode_name, "X", 1)
     call set_count(gnostics%sfile, field2_by_mode_name, "Y", 1)
     if (write_field_by_time) then
       call set_count(gnostics%sfile, field_t_name, "X", 1)
       call set_count(gnostics%sfile, field_t_name, "Y", 1)
     end if

     ! For some reason every process has to make at least
     ! one write to a variable with an infinite dimension.
     ! Here we make some dummy writes to satisfy that
     if (write_field_by_time) call write_variable_with_offset(gnostics%sfile, field_t_name, dummyc)
     call write_variable_with_offset(gnostics%sfile, field2_by_mode_name, dummyr)

     do ik = 1,naky
       do it = 1,ntheta0
         if (field_k_local(it,ik)) then

           call set_start(gnostics%sfile, field_name, "X", it)
           call set_start(gnostics%sfile, field_name, "Y", ik)
           call write_variable_with_offset(gnostics%sfile, field_name, val)

           call set_start(gnostics%sfile, field2_by_mode_name, "X", it)
           call set_start(gnostics%sfile, field2_by_mode_name, "Y", ik)
           call write_variable_with_offset(gnostics%sfile, field2_by_mode_name, field2_by_mode)

           if (write_field_by_time) then
             call set_start(gnostics%sfile, field_t_name, "X", it)
             call set_start(gnostics%sfile, field_t_name, "Y", ik)
             call write_variable_with_offset(gnostics%sfile, field_t_name, val)
           end if
         end if
       end do
     end do
   end if
  end subroutine create_and_write_field_by_mode

  subroutine write_movie(gnostics)
    use gs2_transforms, only: init_transforms
    use gs2_layouts, only: yxf_lo
    use theta_grid, only: ntgrid
    use run_parameters, only: fphi, fapar, fbpar
    use gs2_transforms, only: transform2
    use gs2_io, only: nc_loop_movie
    use fields_arrays, only: phinew, aparnew, bparnew
    use mp, only: proc0, mp_abort
    use diagnostics_create_and_write, only: create_and_write_variable
    use le_grids, only: nlambda, negrid
    use species, only: nspec
    use kt_grids, only: nx, ny
    use file_utils, only: error_unit
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real :: t
    integer :: nnx, nny
    real, dimension(:,:,:), allocatable :: yxphi, yxapar, yxbpar
    logical :: accelerated ! dummy variable here
    logical, save :: warning = .false.

    call init_transforms (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)

   

    t = gnostics%user_time

    ! EAB 09/17/03 -- modify dump_fields_periodically to print out inverse fft of fields in x,y
    nnx = yxf_lo%nx
    nny = yxf_lo%ny

    !write (*,*) 'Grid sizes', nnx, nny

    if (proc0 .and. .not. warning) then
      write (error_unit(),*) "WARNING: write_movie writes arrays of size nx*ny which contain a lot of &
        & redundancy, leading to unnecessarily large output files... consider doing &
        & the Fourier transforms of the fields in post-processing instead"
      warning = .true.
    endif


    !<DD>Commented as removed writing of ntot in favour of apar for consistency
    !call getmoms (phinew, bparnew, ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)

    if (fphi > epsilon(0.0)) then
       allocate (yxphi(nnx,nny,-ntgrid:ntgrid))
       call transform2 (phinew, yxphi, nny, nnx)
       call create_and_write_variable(gnostics, gnostics%rtype, "phi_by_xmode", "xyzt", &
          "The electric potential in real space", "Tr/e", yxphi)
    end if

    if (fapar > epsilon(0.0)) then
       allocate (yxapar(nnx,nny,-ntgrid:ntgrid))
       call transform2 (aparnew, yxapar, nny, nnx)
       call create_and_write_variable(gnostics, gnostics%rtype, "apar_by_xmode", "xyzt", &
          "The parallel vector potential in real space", "TBC", yxapar)
    end if

    if (fbpar > epsilon(0.0)) then 
       allocate (yxbpar(nnx,nny,-ntgrid:ntgrid))
       call transform2 (bparnew, yxbpar, nny, nnx)
       call create_and_write_variable(gnostics, gnostics%rtype, "bpar_by_xmode", "xyzt", &
          "The parallel magnetic field in real space", "TBC", yxapar)
    end if

    !if (proc0) then
       !call nc_loop_movie(nout_movie, t, yxphi, yxapar, yxbpar)
    !end if

    if (fphi > epsilon(0.0)) deallocate (yxphi)
    if (fapar > epsilon(0.0)) deallocate (yxapar)
    if (fbpar > epsilon(0.0)) deallocate (yxbpar)
    !nout_movie = nout_movie + 1
  end subroutine write_movie


end module diagnostics_fields


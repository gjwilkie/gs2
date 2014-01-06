!> This module contains functions for writing quantities
!! that have the dimensions of the fields, including
!! volume averaged quantities, etc.
module diagnostics_write_fields
	 use simpledataio
	 use simpledataio_write
   use theta_grid, only: ntgrid
   use kt_grids, only: naky, ntheta0
   use diagnostics_config, only: diagnostics_type
   implicit none
   logical :: fields_local
   integer, parameter :: REAL_TYPE = SDATIO_DOUBLE

contains
  !function field_k_local(ik,it)
    !use mp, only: iproc, nproc
    !integer, intent(in) :: ik, it
    !logical :: field_k_local

    !! This is temporary while the fields are being parallelised
    !!write (*,*) 'it', it, ' nproc ', nproc, 'iproc ', iproc, 'mod', mod(it,nproc)
    !field_k_local = (mod(iproc,ntheta0) == it-1)
    !!field_k_local = (mod(it,nproc) == iproc)
  !end function field_k_local

  subroutine write_standard_field_properties(gnostics, field_name, field_description, &
    field_units, field_value)
    use kt_grids, only: ntheta0, naky
    use diagnostics_create_and_write, only: create_and_write_variable
    use volume_averages
    type(diagnostics_type), intent(in) :: gnostics
    character(*), intent(in) :: field_name, field_description, field_units
    complex, dimension(:,:,:), intent(in) :: field_value
    real, dimension(ntheta0, naky) :: field2_by_mode
    real, dimension(naky) :: field2_by_ky
    real, dimension(ntheta0) :: field2_by_kx

    call average_theta(field_value, field_value, field2_by_mode, .true.)
    !call create_and_write_field(gnostics%sfile, field_name, field_description, field_units, field_value)
    call create_and_write_field_by_mode(gnostics, field_name, field_description, field_units, &
      field_value, field2_by_mode, .true.)

    call average_kx(field2_by_mode, field2_by_ky, .true.)
    call create_and_write_variable(gnostics%sfile, REAL_TYPE, field_name//"2_by_ky", "yt", &
      field_description//" squared and averaged over theta and kx, as a function of time", &
      "("//field_units//")^2", field2_by_ky)

    call average_ky(field2_by_mode, field2_by_kx, .true.)
    call create_and_write_variable(gnostics%sfile, REAL_TYPE, field_name//"2_by_kx", "xt", &
      field_description//" squared and averaged over theta and ky, as a function of time", &
      "("//field_units//")^2", field2_by_kx)

    call create_and_write_variable(gnostics%sfile, REAL_TYPE, field_name//"2", "t", &
      field_description//" squared and averaged over theta, kx and ky, as a function of time", &
      "("//field_units//")^2", sum(field2_by_kx))

  end subroutine write_standard_field_properties

  subroutine create_and_write_field_by_mode(gnostics, field_name, field_description, field_units, &
      val, field2_by_mode, write_field_by_time)
   use fields_parallelization, only: field_k_local
   type(diagnostics_type), intent(in) :: gnostics
   character(*), intent(in) :: field_name
   character(*), intent(in) :: field_description
   character(*), intent(in) :: field_units
   character(len=len(field_name)+9) :: field2_by_mode_name
   character(len=len(field_name)+2) :: field_t_name
   complex, intent(in), dimension(-ntgrid:ntgrid, ntheta0, naky) :: val
   real, dimension(ntheta0, naky), intent(in) :: field2_by_mode
   complex, dimension(-ntgrid:ntgrid, 1, 1) :: dummyc
   real, dimension(1,1) :: dummyr
   logical, intent(in) :: write_field_by_time
   integer :: it,ik

   field2_by_mode_name =  field_name//"2_by_mode" 
   field_t_name = field_name//"_t"


	 
	 if (.not. variable_exists(gnostics%sfile, field_name)) then 
	   call create_variable(gnostics%sfile, SDATIO_DOUBLE, field_name, "rzxy", field_description, field_units)
	 end if
	 if (.not. variable_exists(gnostics%sfile, field2_by_mode_name)) then 
	   call create_variable(gnostics%sfile, SDATIO_DOUBLE, field2_by_mode_name, "xyt", &
       field_description//" squared and averaged over theta, as a function of kx and ky" , "("//field_units//")^2")
	 end if
	 if (write_field_by_time .and. .not. variable_exists(gnostics%sfile, field_t_name)) then 
	   call create_variable(gnostics%sfile, SDATIO_DOUBLE, field_t_name, "rzxyt", &
       field_description//" as a function of time", field_units)
	 end if
	 
   if (fields_local) then

     call write_variable(gnostics%sfile, field_name, val)
     call write_variable(gnostics%sfile, field2_by_mode_name, field2_by_mode)
     if (write_field_by_time) call write_variable(gnostics%sfile, field_t_name, val)

   else

     call set_count(gnostics%sfile, field_name, "x", 1)
     call set_count(gnostics%sfile, field_name, "y", 1)
     call set_count(gnostics%sfile, field2_by_mode_name, "x", 1)
     call set_count(gnostics%sfile, field2_by_mode_name, "y", 1)
     if (write_field_by_time) then
       call set_count(gnostics%sfile, field_t_name, "x", 1)
       call set_count(gnostics%sfile, field_t_name, "y", 1)
     end if

     ! For some reason every process has to make at least
     ! one write to a variable with an infinite dimension.
     ! Here we make some dummy writes to satisfy that
     call write_variable(gnostics%sfile, field_t_name, dummyc)
     call write_variable(gnostics%sfile, field2_by_mode_name, dummyr)

     do ik = 1,naky
       do it = 1,ntheta0
         if (field_k_local(it,ik)) then

           call set_start(gnostics%sfile, field_name, "x", it)
           call set_start(gnostics%sfile, field_name, "y", ik)
           call write_variable(gnostics%sfile, field_name, val)

           call set_start(gnostics%sfile, field2_by_mode_name, "x", it)
           call set_start(gnostics%sfile, field2_by_mode_name, "y", ik)
           call write_variable(gnostics%sfile, field2_by_mode_name, field2_by_mode)

           if (write_field_by_time) then
             call set_start(gnostics%sfile, field_t_name, "x", it)
             call set_start(gnostics%sfile, field_t_name, "y", ik)
             call write_variable(gnostics%sfile, field_t_name, val)
           end if
         end if
       end do
     end do
   end if
  end subroutine create_and_write_field_by_mode

end module diagnostics_write_fields


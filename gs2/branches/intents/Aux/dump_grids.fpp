program dump_grids
  use job_manage, only: job_fork, time_message
  use mp, only: init_mp, proc0, nproc, broadcast
  use file_utils, only: init_file_utils, run_name
  use kt_grids, only: init_kt_grids
  use theta_grid, only: init_theta_grid
  use le_grids, only: init_le_grids

  implicit none
  logical :: list
  logical :: accel_x, accel_y
  logical, parameter :: quiet=.false.
  character (500), target :: cbuff
  character (len=500) :: filename
  real :: time_measure(2)

  ! <doc> Initialize message passing </doc>
  call init_mp
  
  ! <doc> Report # of processors being used </doc>
  if (proc0) then
     if(.not.quiet)then
        if (nproc == 1) then
           write(*,'("Running on ",I0," processor")') nproc
        else
           write(*,'("Running on ",I0," processors")') nproc
        end if
     endif

     ! <doc> Call init_file_utils, ie. initialize the inputs and outputs, checking 
     !  whether we are doing a run or a list of runs. </doc>
     call init_file_utils (list, name="gs")
  end if

  call broadcast (list)
  
  ! <doc> If given a list of jobs, fork </doc>
  if (list) call job_fork
  
  if (proc0) cbuff = trim(run_name)
  call broadcast (cbuff)
  if (.not. proc0) run_name => cbuff

  !Initialise
  time_measure=0.
  !/Theta
  if(proc0.and.(.not.quiet)) then 
     write(*,'("Init theta_grids")',advance='no') 
     call time_message(.false.,time_measure,'init-theta')
  endif
  call init_theta_grid
  if(proc0.and.(.not.quiet)) then
     call time_message(.false.,time_measure,'init-theta')
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.
  !/KT
  if(proc0.and.(.not.quiet)) then 
     write(*,'("Init kt_grids   ")',advance='no') 
     call time_message(.false.,time_measure,'init-kt')
  endif
  call init_kt_grids
  if(proc0.and.(.not.quiet)) then
     call time_message(.false.,time_measure,'init-kt')
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.
  !/LE
  if(proc0.and.(.not.quiet)) then 
     write(*,'("Init le_grids   ")',advance='no') 
     call time_message(.false.,time_measure,'init-le')
  endif
  call init_le_grids(accel_x,accel_y)
  if(proc0.and.(.not.quiet)) then
     call time_message(.false.,time_measure,'init-le')
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.

  !Now write to file
  if(proc0.and.(.not.quiet)) then 
     write(*,'("Write to file   ")',advance='no') 
     call time_message(.false.,time_measure,'write-file')
  endif
  filename=trim(adjustl(run_name))
  if(proc0) call write_grids_to_file(filename)
  if(proc0.and.(.not.quiet)) then
     call time_message(.false.,time_measure,'write-file')
     write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
  endif
  time_measure=0.

contains

#ifdef NETCDF
!Prefer to write to netcdf file if possible
  subroutine write_grids_to_file(fname)
    use netcdf, only: NF90_CLOBBER
    use netcdf, only: nf90_create, nf90_close
    use netcdf, only: nf90_def_dim, nf90_def_var, nf90_enddef
    use netcdf, only: nf90_put_var
    
    use netcdf_utils, only: get_netcdf_code_precision
    use netcdf_utils, only: netcdf_real

    use kt_grids, only: aky, akx, theta0
    use theta_grid, only: theta
    use le_grids, only: al, energy

    implicit none
    character(len=*), intent(in) :: fname
    integer :: ncid, ierr
    integer :: ky_dimid, kx_dimid, theta_dimid, al_dimid, energy_dimid
    integer :: ky_varid, kx_varid, theta_varid, al_varid, energy_varid, theta0_varid

    !First get precisions
    if(netcdf_real==0) netcdf_real=get_netcdf_code_precision()

    !First create a file
    ierr=nf90_create(trim(adjustl(fname))//".grids.nc",NF90_CLOBBER,ncid)

    !Define dimensions
    ierr=nf90_def_dim(ncid,'aky',size(aky),ky_dimid)
    ierr=nf90_def_dim(ncid,'akx',size(akx),kx_dimid)
    ierr=nf90_def_dim(ncid,'theta',size(theta),theta_dimid)
    ierr=nf90_def_dim(ncid,'lambda',size(al),al_dimid)
    ierr=nf90_def_dim(ncid,'energy',size(energy),energy_dimid)

    !Define variables
    ierr=nf90_def_var(ncid,'ky',netcdf_real,(/ky_dimid/),ky_varid)
    ierr=nf90_def_var(ncid,'kx',netcdf_real,(/kx_dimid/),kx_varid)
    ierr=nf90_def_var(ncid,'theta0',netcdf_real,(/kx_dimid,ky_dimid/),theta0_varid)
    ierr=nf90_def_var(ncid,'theta',netcdf_real,(/theta_dimid/),theta_varid)
    ierr=nf90_def_var(ncid,'lambda',netcdf_real,(/al_dimid/),al_varid)
    ierr=nf90_def_var(ncid,'energy',netcdf_real,(/energy_dimid/),energy_varid)

    !End definitions
    ierr=nf90_enddef(ncid)

    !Now fill in variable data
    ierr=nf90_put_var(ncid,ky_varid,aky)
    ierr=nf90_put_var(ncid,kx_varid,akx)
    ierr=nf90_put_var(ncid,theta0_varid,theta0)
    ierr=nf90_put_var(ncid,theta_varid,theta)
    ierr=nf90_put_var(ncid,al_varid,al)
    ierr=nf90_put_var(ncid,energy_varid,energy)

    !Now close file
    ierr=nf90_close(ncid)
  end subroutine write_grids_to_file
#else
!Fall back to binary output when netcdf unavailable
  subroutine write_grids_to_file(fname)
    use kt_grids, only: aky, akx, theta0
    use theta_grid, only: theta
    use le_grids, only: al, energy
    use file_utils, only: get_unused_unit

    implicit none
    character(len=*), intent(in) :: fname
    integer :: unit

    call get_unused_unit(unit)

    open(unit=unit,file=trim(adjustl(fname))//".ky",form="unformatted")
    write(unit) aky
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".kx",form="unformatted")
    write(unit) akx
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".theta0",form="unformatted")
    write(unit) theta0
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".theta",form="unformatted")
    write(unit) theta
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".lambda",form="unformatted")
    write(unit) al
    close(unit)

    open(unit=unit,file=trim(adjustl(fname))//".energy",form="unformatted")
    write(unit) energy
    close(unit)
  end subroutine write_grids_to_file
#endif
end program dump_grids

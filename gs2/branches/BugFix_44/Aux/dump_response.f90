!> A program to generate response matrix data for a given input file 
!! and dump the resulting arrays to file to be used by subsequent runs
!!  Written by : David Dickinson (ddickinson@users.sourceforge.net)
program dump_response
  use job_manage, only: job_fork
  use mp, only: init_mp, proc0, nproc, broadcast
  use file_utils, only: init_file_utils, run_name, input_unit_exist
  implicit none
  logical :: list, exist, jnk_exit
  logical, parameter :: quiet=.false.
  character (500), target :: cbuff
  integer :: istep, n_time_steps, in_file
  namelist /dump_response_knobs/n_time_steps

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
  
  if(proc0) then
     n_time_steps=1
     in_file = input_unit_exist ("dump_response_knobs", exist)
     if(exist) read(unit=in_file,nml=dump_response_knobs)
  endif
  call broadcast(n_time_steps)

  !This routine initialises the response matrix and dumps it to file.
  call init_and_dump(quiet)
  do istep=2,n_time_steps
     !Reduce the time step
     call next_time_step

     !Dump response
     call init_and_dump(quiet)
  enddo

  !Clean up
  call finish_dump_response

contains
  subroutine next_time_step
    use gs2_time, only: code_dt, save_dt
    use gs2_reinit, only: reduce_time_step
    use fields, only: f_reset => reset_init
    use collisions, only: c_reset => reset_init
    use dist_fn, only: d_reset => reset_init
    implicit none

    !First reduce the time step
    call reduce_time_step

    !Save the time step
    call save_dt(code_dt)

    !Now reset the modules
    call d_reset
    call c_reset
    call f_reset
  end subroutine next_time_step

  subroutine init_and_dump(quiet)
    use fields, only: fields_pre_init, fields_init_response, set_dump_and_read_response, dump_response_to_file
    use job_manage, only: time_message
    use gs2_time, only: code_dt
    use mp, only: proc0
    implicit none
    logical, intent(in) :: quiet
    real :: time_measure(2)
    character (len=64) :: suffix, str_code_dt

    if(.not.quiet) write(*,'("")')

    !Prepare to initialise the fields
    if(proc0.and.(.not.quiet)) then 
       write(*,'("Pre-init")',advance='no') 
       call time_message(.false.,time_measure,'Pre-init')
    endif
    call fields_pre_init
    if(proc0.and.(.not.quiet)) then
       call time_message(.false.,time_measure,'Pre-init')
       write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
    endif
    time_measure=0.

    !Now need to disable dump_response and read_response before
    !we try to initialise the response matrix
    call set_dump_and_read_response(.false.,.false.)

    !Now call the initialisation routine
    if(proc0.and.(.not.quiet)) then
       write(*,'("Init")', advance='no') 
       call time_message(.false.,time_measure,'Init')
    endif
    call fields_init_response
    if(proc0.and.(.not.quiet)) then
       call time_message(.false.,time_measure,'Init')
       write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
    endif
    time_measure=0.

    !Now we make the file suffix
    write(str_code_dt,'(F12.10)') code_dt
    write(suffix,'(".dt_",A,".response")') trim(adjustl(str_code_dt))

    !Now write response to files
    if(proc0.and.(.not.quiet)) then
       write(*,'("Dumping to file with suffix ",A,A,A)',advance='no') "'",trim(suffix),"'"
       call time_message(.false.,time_measure,'Dump')
    endif
    call dump_response_to_file(suffix)
    if(proc0.and.(.not.quiet)) then
       call time_message(.false.,time_measure,'Dump')
       write(*,'("  --> Done in : ",F12.6," seconds")') time_measure(1)
    endif
    time_measure=0.
  end subroutine init_and_dump

  subroutine finish_dump_response
    use antenna, only: finish_antenna
    use collisions, only: finish_collisions
    use dist_fn, only: finish_dist_fn
    use fields, only: finish_fields
    use file_utils, only: finish_file_utils
    use hyper, only: finish_hyper
    use init_g, only: finish_init_g
    use kt_grids, only: finish_kt_grids
    use le_grids, only: finish_le_grids
    use mp, only: proc0, finish_mp
    use run_parameters, only: finish_run_parameters
    use species, only: finish_species

    implicit none

    call finish_antenna
    call finish_collisions
    call finish_dist_fn
    call finish_fields
    call finish_hyper
    call finish_init_g
    call finish_kt_grids
    call finish_le_grids
    call finish_run_parameters
    call finish_species
    if (proc0) call finish_file_utils
    call finish_mp

  end subroutine finish_dump_response

end program dump_response


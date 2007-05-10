program gs2
  use mp, only: init_mp, finish_mp, proc0, nproc, broadcast
  use file_utils, only: init_file_utils, finish_file_utils, run_name
  use fields, only: init_fields
  use gs2_diagnostics, only: init_gs2_diagnostics, finish_gs2_diagnostics

  use run_parameters, only: nstep, delt, tnorm
  use fields, only: advance
  use gs2_diagnostics, only: loop_diagnostics
  use gs2_reinit, only: reset_time_step, check_time_step
  use gs2_time, only: init_time, write_dt
  use init_g, only: tstart
  use check, only: checkstop

  implicit none
  real :: dt_cfl
  integer :: istep, istep_end, unit
  logical :: exit, reset

  call init_mp

  if (nproc > 1 .and. proc0) write(*,*) 'Running on ',nproc,' processors'
  if (nproc == 1) write(*,*) 'Running on ',nproc,' processor'

  if (proc0) call init_file_utils (name="gs")
  call broadcast (run_name)

  call init_fields
  call init_gs2_diagnostics
  call init_time (tstart*tnorm, delt)

  istep_end = nstep
  do istep = 1, nstep

     call advance (istep, dt_cfl)
     call loop_diagnostics (istep, exit)
     call check_time_step (istep, delt, dt_cfl, reset, exit)
     if (reset) call reset_time_step (istep, dt_cfl, exit)

     if (mod(istep,5) == 0) call checkstop(exit)

     if (exit) then
        istep_end = istep
        exit
     end if
  end do
  if (proc0) call write_dt(tnorm)

  call finish_gs2_diagnostics (istep_end)
  if (proc0) call finish_file_utils
  call finish_mp

end program gs2





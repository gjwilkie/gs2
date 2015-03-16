!> This module sits on top of gs2_main and provides
!! a toolkit for measuring and optimising performance
module gs2_optimisation
  use gs2_main, only: gs2_program_state_type


  contains

    subroutine initialize_gs2_optimisation(state)
      use optimisation_config, only: init_optimisation_config
      use gs2_main, only: initialize_gs2, finalize_gs2
      use mp, only: init_mp, mp_comm
      type(gs2_program_state_type), intent(inout) :: state
      call init_mp
      state%mp_comm_external = .true.
      state%mp_comm = mp_comm
      allocate(state%optim%sorted_optimisations(0))
      allocate(state%optim%sorted_results(0))
      write (*,*) 'size', size(state%optim%sorted_results)
      ! We have to initialize_gs2 so that we 
      ! can read the optimisation_config namelist
      call initialize_gs2(state)
      call init_optimisation_config(state%optim)
      call finalize_gs2(state)
    end subroutine initialize_gs2_optimisation

    subroutine finalize_gs2_optimisation(state)
      use mp, only: finish_mp
      use optimisation_config, only: finish_optimisation_config
      type(gs2_program_state_type), intent(inout) :: state
      call finish_optimisation_config(state%optim)
      call finish_mp
    end subroutine finalize_gs2_optimisation

    subroutine optimise_gs2(state)
      type(gs2_program_state_type), intent(inout) :: state

      state%print_times = .false.
      state%is_external_job = .true.
      ! Initialize optimisation results
      !call measure_timestep(state)
      state%optim%results%optimal_time = -1.0
      state%optim%results%optimal = .false.

      if (state%optim%max_unused_procs .eq. 0 &
          .or. state%optim%max_imbalance .eq. 0.0) then
        call optimise_layout(state)
      else
        call optimise_nprocs(state)
      end if

      call output_results(state)

    end subroutine optimise_gs2

    subroutine output_results(state)
      use mp, only: proc0
      implicit none
      type(gs2_program_state_type), intent(inout) :: state
      integer :: i,n

      if (proc0) then
        write (state%optim%outunit, '(A16," ",A6)') 'timestep-time', 'layout' 
        n = size(state%optim%sorted_results)
        do i = 1,n
            write(state%optim%outunit, &
              '(E16.9," ",A6)') &
            state%optim%sorted_results(i)%time, &
            state%optim%sorted_optimisations(i)%layout
        end do
      end if

    end subroutine output_results

    subroutine optimise_nprocs(state)
      type(gs2_program_state_type), intent(inout) :: state
    end subroutine optimise_nprocs

    subroutine optimise_layout(state)
      use gs2_main, only: prepare_optimisations_overrides
      use mp, only: proc0
      type(gs2_program_state_type), intent(inout) :: state
      call prepare_optimisations_overrides(state)
      state%init%opt_ov%override_layout = .true.
      state%init%opt_ov%layout = 'lxyes'
      call optimise_flags(state)
      state%init%opt_ov%layout = 'lexys'
      call optimise_flags(state)
      state%init%opt_ov%layout = 'xyles'
      call optimise_flags(state)
      state%init%opt_ov%layout = 'yxles'
      call optimise_flags(state)
    end subroutine optimise_layout

    subroutine optimise_flags(state)
      type(gs2_program_state_type), intent(inout) :: state
      call measure_timestep(state)
    end subroutine optimise_flags

    subroutine measure_timestep(state)
      use gs2_main, only: gs2_program_state_type
      use gs2_main, only: initialize_gs2
      use gs2_main, only: initialize_equations
      use gs2_main, only: initialize_diagnostics
      use gs2_main, only: evolve_equations
      use gs2_main, only: run_eigensolver
      use gs2_main, only: finalize_diagnostics
      use gs2_main, only: finalize_equations
      use gs2_main, only: finalize_gs2
      use overrides, only: optimisations_overrides_type
      use optimisation_config, only: optimisation_results_type

      implicit none
      type(gs2_program_state_type), intent(inout) :: state
      type(optimisations_overrides_type), &
        dimension(:), allocatable :: sorted_opts_temp
      type(optimisation_results_type), &
        dimension(:), allocatable :: sorted_res_temp
      integer :: i,n
      real :: t

      call initialize_gs2(state)
      call initialize_equations(state)
      call initialize_diagnostics(state)
      call evolve_equations(state, state%optim%nstep_measure)
      call finalize_diagnostics(state)
      call finalize_equations(state)
      call finalize_gs2(state)

      if (state%optim%measure_all) then 
        t = state%timers%advance(1)/real(state%optim%nstep_measure)
      else
        t = state%timers%timestep(1)/real(state%optim%nstep_measure)
        !t = state%timers%timestep(1)
      endif

      state%optim%results%time = t
      if (t .lt.  state%optim%results%optimal_time .or. &
        state%optim%results%optimal_time .lt. 0.0) then
        state%optim%results%optimal_time = t
        state%optim%results%optimal = .true.
      end if

      n = size(state%optim%sorted_results)
      write (*,*) 'size2', size(state%optim%sorted_results)
      allocate(sorted_opts_temp(n), sorted_res_temp(n))
      if (n>0) then
        do i = 1,n
          sorted_opts_temp(i) = state%optim%sorted_optimisations(i)
          sorted_res_temp(i) = state%optim%sorted_results(i)
          end do
      end if

      deallocate(state%optim%sorted_optimisations)
      deallocate(state%optim%sorted_results)
      allocate(state%optim%sorted_optimisations(n+1))
      allocate(state%optim%sorted_results(n+1))

      i=1
      do 
        if (i>n) exit
        if (sorted_res_temp(i)%time > t) exit
        state%optim%sorted_optimisations(i) = sorted_opts_temp(i)
        state%optim%sorted_results(i) = sorted_res_temp(i)
        i = i+1
      end do
      state%optim%sorted_optimisations(i) = state%init%opt_ov
      state%optim%sorted_results(i) = state%optim%results
      i = i + 1
      do while (i < n+2)
        state%optim%sorted_optimisations(i) = sorted_opts_temp(i-1)
        state%optim%sorted_results(i) = sorted_res_temp(i-1)
        i = i + 1
      end do
        
      deallocate(sorted_opts_temp, sorted_res_temp)
      
    end subroutine measure_timestep

end module gs2_optimisation

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

      ! Initialize optimisation results
      call measure_timestep(state)
      state%optim%results%optimal_time = state%optim%results%last_time
      state%optim%results%optimal = .true.

      if (state%optim%max_unused_procs .eq. 0 &
          .or. state%optim%max_imbalance .eq. 0.0) then
        call optimise_layout(state)
      else
        call optimise_nprocs(state)
      end if

    end subroutine optimise_gs2

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
      if (proc0) write (*,*) 'layout', state%init%opt_ov%layout, &
        'time', state%optim%results%last_time
      state%init%opt_ov%layout = 'lexys'
      call optimise_flags(state)
      if (proc0) write (*,*) 'layout', state%init%opt_ov%layout, &
        'time', state%optim%results%last_time
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

      implicit none
      type(gs2_program_state_type), intent(inout) :: state

      call initialize_gs2(state)
      call initialize_equations(state)
      call initialize_diagnostics(state)
      call evolve_equations(state, 20)
      call finalize_diagnostics(state)
      call finalize_equations(state)
      call finalize_gs2(state)

      state%optim%results%last_time = state%timers%advance(1)/10.0
      if (state%optim%results%last_time .le. &
          state%optim%results%optimal_time) then
        state%optim%results%optimal_time = &
        state%optim%results%last_time
        state%optim%results%optimal = .true.
      end if
      
    end subroutine measure_timestep

end module gs2_optimisation

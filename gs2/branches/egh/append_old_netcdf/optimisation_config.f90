
! DO NOT EDIT THIS FILE
! This file has been automatically generated using generate_optimisation_config.rb

!> A module for handling the configuration of the optimisation
!! module via the namelist optimisation_config.
module optimisation_config
  use overrides, only: optimisations_overrides_type
  implicit none

  private

  public :: init_optimisation_config
  public :: finish_optimisation_config
  public :: optimisation_type
  public :: optimisation_results_type

  type optimisation_results_type
    ! Configuration

    ! Results
    real :: time
    real :: optimal_time
    real :: cost
    real :: optimal_cost
    real :: efficiency
    integer :: nproc
    logical :: optimal = .true.

  end type optimisation_results_type


  !> A type for storing the optimisation configuration,
  !! the results
  type optimisation_type
    integer :: nproc_max
    type(optimisation_results_type) :: results
    type(optimisations_overrides_type), &
      dimension(:), pointer :: sorted_optimisations
    type(optimisation_results_type), dimension(:), pointer :: sorted_results
    integer :: outunit
     logical :: on
     logical :: auto
     logical :: measure_all
     logical :: warm_up
     integer :: nstep_measure
     real :: max_imbalance
     integer :: max_unused_procs
     real :: min_efficiency
  end type optimisation_type

contains
  subroutine init_optimisation_config(optim)
    use file_utils, only: open_output_file
    use mp, only: nproc
    implicit none
    type(optimisation_type), intent(inout) :: optim
    call read_parameters(optim)
    call open_output_file(optim%outunit, '.optim')
    optim%nproc_max = nproc
  end subroutine init_optimisation_config

  subroutine finish_optimisation_config(optim)
    implicit none
    type(optimisation_type), intent(inout) :: optim
  end subroutine finish_optimisation_config


  subroutine read_parameters(optim)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type(optimisation_type), intent(inout) :: optim
    logical :: on
    logical :: auto
    logical :: measure_all
    logical :: warm_up
    integer :: nstep_measure
    real :: max_imbalance
    integer :: max_unused_procs
    real :: min_efficiency
    namelist /optimisation_config/ &
         on, &
         auto, &
         measure_all, &
         warm_up, &
         nstep_measure, &
         max_imbalance, &
         max_unused_procs, &
         min_efficiency

    integer :: in_file
    logical :: exist

    if (proc0) then
       on = .false.
       auto = .true.
       measure_all = .false.
       warm_up = .false.
       nstep_measure = 5
       max_imbalance = -1
       max_unused_procs = 0
       min_efficiency = -1.0

       in_file = input_unit_exist ("optimisation_config", exist)
       if (exist) read (unit=in_file, nml=optimisation_config)

       optim%on = on
       optim%auto = auto
       optim%measure_all = measure_all
       optim%warm_up = warm_up
       optim%nstep_measure = nstep_measure
       optim%max_imbalance = max_imbalance
       optim%max_unused_procs = max_unused_procs
       optim%min_efficiency = min_efficiency

    end if

    call broadcast (optim%on)
    call broadcast (optim%auto)
    call broadcast (optim%measure_all)
    call broadcast (optim%warm_up)
    call broadcast (optim%nstep_measure)
    call broadcast (optim%max_imbalance)
    call broadcast (optim%max_unused_procs)
    call broadcast (optim%min_efficiency)
    
  end subroutine read_parameters
end module optimisation_config

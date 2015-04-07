# This Ruby script automatically generates the module diagnostics_config. It is compatible with any version of Ruby, which means that it should work out of the box on any system, since even the most decrepit systems usually have Ruby 1.8.7. It is automatically invoked by the Makefile; thus, after editing it a simple make will trigger the generation of the file diagnostics_config.f90.
#

# This is a list of all the input parameters to the namelist
# diagnostics_config, with their default values.
# Leave third element of array empty to use the default default
#
# This list is used in generating diagnostics_config.f90
input_variables_for_diagnostics_config = [
  ['integer', 'nwrite', '10'],
  ['integer', 'nwrite_large', '100'],
  ['logical', 'write_any', '.true.'],

  # If built with parallel IO capability 
  # enable parallel IO. Currently disabled by 
  # default because of problems on Helios:
  # parallel IO seems to work on 1 node but not
  # on >1. This problem does not affect parallel
  # IO for the single restart file because it only 
  # affects output with infinite dimensions
  ['logical', 'enable_parallel', '.false.'],

  # If true, produce netcdf-4 files when running
  # with serial IO. Otherwise produce classic file
  # format. Parallel IO always produces netcdf-4 files.
  ['logical', 'serial_netcdf4', '.false.'],

  # Controls the theta location of omega calculation,
  # also of any quantities that are written out for
  # a given value of theta
  ['integer', 'igomega', '0'],

  # Write instantaneous quantities to screen
  ['logical', 'print_line', '.false.'],
  ['logical', 'print_flux_line', '.false.'],

  # Write instantaneous quantites to .new.out
  ['logical', 'write_line', '.true.'],
  ['logical', 'write_flux_line', '.true.'],

  # Parameters for writing out fields
  ['logical', 'write_fields', '.true.'],
  ['logical', 'write_phi_over_time'],
  ['logical', 'write_apar_over_time'],
  ['logical', 'write_bpar_over_time'],
  ['logical', 'write_movie', '.false.'],
  ['logical', 'dump_fields_periodically', '.false.'],

  # Parameters for writing out moments such as density etc
  ['logical', 'write_moments', '.true.'],
  ['logical', 'write_full_moments_notgc', '.false.'],
  # Write out 4-D moments as a function of time ... gives
  # LARGE data files!
  ['logical', 'write_ntot_over_time', '.false.'],
  ['logical', 'write_density_over_time', '.false.'],
  ['logical', 'write_upar_over_time', '.false.'],
  ['logical', 'write_tperp_over_time', '.false.'],


  # Parameters for writing out fluxes
  ['logical', 'write_fluxes', '.true.'],
  ['logical', 'write_fluxes_by_mode', '.false.'],
  ['logical', 'write_symmetry', '.false.'],
  ['logical', 'write_parity', '.false.'],

  # Parameters for writing out growth rates and frequencies
  ['logical', 'write_omega', '.true.'],
  ['integer', 'navg', '10'],
  ['real', 'omegatinst', '1.0e6'],
  ['real', 'omegatol', '-0.001'],
  ['logical', 'exit_when_converged', '.true.'],

  # Parameters for writing out velocity space diagnostics
  ['logical', 'write_verr', '.true.'],
  ['logical', 'write_max_verr', '.false.'],
  ['integer', 'ncheck', '10'],

  # Parameters controlling heating diagnositics
  ['logical', 'write_heating', '.false.'],


  # If true, write out old-style text files
  ['logical', 'write_ascii', '.true.'],

  # If true, write the dist_fn at a range of points in space
  # as a function of velocity to an output file. Do not enable
  # in the new and old diagnostics modules at the same time
  ['logical', 'write_gyx', '.false.'],
  ['logical', 'write_g', '.false.'],
  ['logical', 'write_lpoly', '.false.'],


  # Switches on write collision_error
  ['logical', 'write_cerr', '.false.'],

  # Parameters controlling Trinity convergence tests
  ['integer', 'conv_nstep_av', '4000'],
  ['real', 'conv_test_multiplier', '4e-1'],
  ['integer', 'conv_min_step', '4000'],
  ['integer', 'conv_max_step', '80000'],
  ['integer', 'conv_nsteps_converged', '10000'],
  ['logical', 'use_nonlin_convergence', '.false.'],

  # Parameters determining what turbulence characteristics are calculated
  ['logical', 'write_cross_phase', '.false.'],
  ['logical', 'write_correlation', '.true.'],
  ['logical', 'write_correlation_extend', '.false.'],

  # Parameters controlling diagnostics for the antennna
  ['logical', 'write_jext', '.false.'],
  ['logical', 'write_lorentzian', '.false.'],

  # Parameters controlling some old routines which dump
  # stuff to text files at the end of the simulation
  ['logical', 'write_eigenfunc', '.false.'],
  ['logical', 'write_final_fields', '.false.'],
  ['logical', 'write_kpar', '.false.'],
  ['logical', 'write_final_epar', '.false.'],
  ['logical', 'write_final_db', '.false.'],
  ['logical', 'write_final_moments', '.false.'],
  ['logical', 'write_final_antot', '.false.'],
  ['logical', 'write_gs', '.false.'],

  # Save the current state of the simulation so 
  # that it can be restarted
  # At the moment the default is false because
  # we don't want it to conflict with the old
  # module, but eventually I think it should
  # default to true. EGH
  ['integer', 'nsave', '1000'],
  ['logical', 'save_for_restart', '.false.'],
  ['logical', 'file_safety_check', '.true.'],

  # Save the distribution function
  ['logical', 'save_distfn', '.false.'],

]


class Generator
  attr_accessor :name, :type
  def initialize(type, name, default)
    @type = type
    @name = name
    @default = default
  end
  def declaration
    "#@type :: #@name"
  end
  def set_default
    "#@name = #{default}"
  end
  def default
    @default || case type
                when /logical/
                  '.false.'
                when /integer/
                  '1'
                end
  end
  def parameters_type_value
    "gnostics%#@name"
  end
  def set_parameters_type_value
    "#{parameters_type_value} = #@name"
  end
  def broadcast
    "call broadcast (#{parameters_type_value})"
  end
end

begin
  4.times.map{|i|}
rescue
  puts "You appear to be running ruby 1.8.6 or lower... suggest you upgrade your ruby version!"
  class Integer
    def times(&block)
      if block
        (0...self).to_a.each{|i| yield(i)}
      else
        return  (0...self).to_a
      end
    end
  end
end
generators = input_variables_for_diagnostics_config.map{|type, name, default| Generator.new(type,name,default)}

string = <<EOF

! DO NOT EDIT THIS FILE
! This file has been automatically generated using generate_diagnostics_config.rb

!> A module for handling the configuration of the diagnostics
!! module via the namelist diagnostics_config.
module diagnostics_config
  use simpledataio, only: sdatio_file
  use diagnostics_ascii, only: diagnostics_ascii_type
  use diagnostics_dimensions, only: diagnostics_dimension_list_type
  implicit none

  private

  public :: init_diagnostics_config
  public :: finish_diagnostics_config
  public :: diagnostics_type
  public :: results_summary_type
  public :: override_screen_printout_options

  !> A type for storing the current results of the simulation
  type results_summary_type
     real :: phi2
     real :: apar2
     real :: bpar2
     real :: total_heat_flux
     real :: total_momentum_flux
     real :: total_particle_flux
     real :: max_growth_rate

     ! Individual heat fluxes
     real, dimension(:), allocatable :: species_es_heat_flux
     real, dimension(:), allocatable :: species_apar_heat_flux
     real, dimension(:), allocatable :: species_bpar_heat_flux

     ! Total fluxes
     real, dimension(:), allocatable :: species_heat_flux
     real, dimension(:), allocatable :: species_momentum_flux
     real, dimension(:), allocatable :: species_particle_flux
     real, dimension(:), allocatable :: species_energy_exchange

     ! Average total fluxes
     real, dimension(:), allocatable :: species_heat_flux_avg
     real, dimension(:), allocatable :: species_momentum_flux_avg
     real, dimension(:), allocatable :: species_particle_flux_avg

     ! Heating
     real, dimension(:), allocatable :: species_heating
     real, dimension(:), allocatable :: species_heating_avg
  end type results_summary_type

  !> A type for storing the diagnostics configuration,
  !! a reference to the output file, and current 
  !! results of the simulation
  type diagnostics_type
     type(sdatio_file) :: sfile
     type(diagnostics_ascii_type) :: ascii_files
     type(results_summary_type) :: current_results
     type(diagnostics_dimension_list_type) :: dims
     !> Integer below gives the sdatio type 
     !! which corresponds to a gs2 real
     integer :: rtype
     integer :: itype
     integer :: istep
     logical :: create
     logical :: wryte
     logical :: distributed
     logical :: parallel
     logical :: exit
     logical :: vary_vnew_only
     logical :: calculate_fluxes
     logical :: is_trinity_run
     real :: user_time
     real :: user_time_old
     real, dimension(:), allocatable :: fluxfac
     #{generators.map{|g| g.declaration}.join("\n     ") }
  end type diagnostics_type

  !> Used for testing... causes screen printout to be 
  !! generated regardless of the values of print_line 
  !! and print_flux_line if set to true
  logical :: override_screen_printout_options = .false.

contains
  subroutine init_diagnostics_config(gnostics)
    implicit none
    type(diagnostics_type), intent(out) :: gnostics
    call read_parameters(gnostics)
    call allocate_current_results(gnostics)
  end subroutine init_diagnostics_config

  subroutine finish_diagnostics_config(gnostics)
    implicit none
    type(diagnostics_type), intent(out) :: gnostics
    call deallocate_current_results(gnostics)
  end subroutine finish_diagnostics_config

  subroutine allocate_current_results(gnostics)
    use species, only: nspec
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics

    allocate(gnostics%current_results%species_es_heat_flux(nspec))
    allocate(gnostics%current_results%species_apar_heat_flux(nspec))
    allocate(gnostics%current_results%species_bpar_heat_flux(nspec))

    allocate(gnostics%current_results%species_heat_flux(nspec))
    allocate(gnostics%current_results%species_momentum_flux(nspec))
    allocate(gnostics%current_results%species_particle_flux(nspec))
    allocate(gnostics%current_results%species_energy_exchange(nspec))
    allocate(gnostics%current_results%species_heat_flux_avg(nspec))
    allocate(gnostics%current_results%species_momentum_flux_avg(nspec))
    allocate(gnostics%current_results%species_particle_flux_avg(nspec))
    allocate(gnostics%current_results%species_heating(nspec))
    allocate(gnostics%current_results%species_heating_avg(nspec))
  end subroutine allocate_current_results

  subroutine deallocate_current_results(gnostics)
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    ! This routine needs to be fixed: I don't know 
    ! how to correctly deallocate the derived type
    !<DD>What's written below should work fine, just need to
    !deallocate anything explicitly allocated.
    return
    
    ! One call deallocates gnostics and all allocatable arrays 
    ! within it
    !deallocate(gnostics%current_results)
    write (*,*) "DEALLOCATING"
    deallocate(gnostics%current_results%species_heat_flux)
    deallocate(gnostics%current_results%species_momentum_flux)
    deallocate(gnostics%current_results%species_particle_flux)
    deallocate(gnostics%current_results%species_heat_flux_avg)
    deallocate(gnostics%current_results%species_momentum_flux_avg)
    deallocate(gnostics%current_results%species_particle_flux_avg)
    deallocate(gnostics%current_results%species_heating)
    deallocate(gnostics%current_results%species_heating_avg)
  end subroutine deallocate_current_results


  subroutine read_parameters(gnostics)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type(diagnostics_type), intent(out) :: gnostics
    #{generators.map{|g| g.declaration}.join("\n    ") }
    namelist /diagnostics_config/ &
         #{generators.map{|g| g.name}.join(", &\n         ")}

    integer :: in_file
    logical :: exist

    if (proc0) then
       #{generators.map{|g| g.set_default}.join("\n       ")}

       in_file = input_unit_exist ("diagnostics_config", exist)
       if (exist) read (unit=in_file, nml=diagnostics_config)

       #{generators.map{|g| g.set_parameters_type_value}.join("\n       ")}

    end if

    #{generators.map{|g| g.broadcast}.join("\n    ")}
    
    if (override_screen_printout_options) then 
       gnostics%print_line = .true.
       gnostics%print_flux_line = .true.
    end if
  end subroutine read_parameters
end module diagnostics_config

EOF

File.open(ARGV[-1], 'w'){|f| f.puts string}

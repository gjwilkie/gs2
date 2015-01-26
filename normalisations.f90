module normalisations
  implicit none
  
  private

  public :: norms
  public :: init_normalisations, finish_normalisations
  public :: check_normalisations, wnml_normalisations

  !>Define a type to hold the set of normalisations, we
  !!could imagine extending this in the future to provide
  !!a normalise/denormalise routine which can be used to
  !!automatically convert quantities
  type norms_type
     private
     !Internal parameters
     real :: def_val = -999.0 !Default value, used to check if parameter has been set.

     !Note def_val is really a parameter, but not allowed to define params in derived type
     logical :: some_norms_set=.false.
     logical :: all_norms_set=.false.
     logical :: initialised = .false.

     !The following is used to help automate loops over normalisations
     integer :: nnorm
     character(len=6), dimension(:), allocatable, public :: names

     !The normalisations
     real :: aref !Reference mass in atomic mass units
     real :: zref !Reference charge in units of the elementary charge
     real :: nref !Reference density in m^3
     real :: tref !Reference temperature in eV
     real :: lref !Reference length in m
     real :: vref !Reference (thermal) velocity in m/s
     real :: bref !Reference magnetic field in Tesla
     real :: rhoref !Reference Larmor radius in m
   contains
     private
     procedure, public :: init => norms_init
     procedure, public :: finish => norms_finish
     procedure, public :: reset => norms_reset
     procedure :: read_parameters => norms_read_parameters
     procedure :: set_value => norms_set_value
     procedure, public :: get_value => norms_get_value
     procedure :: check_if_set => norms_check_if_set
     procedure :: set_logicals => norms_set_logicals
  end type norms_type

  type(norms_type) :: norms !The normalisation object, there should only be one of these
  logical :: exist !Has this namelist been specified?

contains
  !/////////////////////////////
  !// TYPE BOUND PROCEDURES
  !/////////////////////////////
  subroutine norms_set_value(self, val_name, val)
    use runtime_tests, only: verbosity
    use mp, only: proc0
    implicit none
    class(norms_type), intent(in out) :: self
    character(len=*), intent(in) :: val_name
    real, intent(in) :: val
    
    !Here we can have a message in very verbose runs
    if((verbosity().gt.5).and.proc0)then
       if(val.eq.self%def_val)then
          write(6,'("The ",A," normalisation has not been set. This may prevent conversion of some quantities.")') trim(val_name)
       else
          write(6,'("Setting the ",A," normalisation to ",F12.5)') trim(val_name), val
       endif
    endif

    !Should probably convert to lower case here, but until we
    !add some string utils to do this sort of stuff we'll rely
    !on developer doing the right thing.
    select case (trim(val_name))
    case("aref")
       self%aref=val
    case("zref")
       self%zref=val
    case("nref")
       self%nref=val
    case("tref")
       self%tref=val
    case("lref")
       self%lref=val
    case("vref")
       self%vref=val
    case("bref")
       self%bref=val
    case("rhoref")
       self%rhoref=val
    case default
       if(proc0) write(6,'("Warning : Attempt to set unknown normalisation ",A," --> Ignoring")')
    end select
  end subroutine norms_set_value

  elemental function norms_get_value(self, val_name)
    !use mp, only: mp_abort
    implicit none
    class(norms_type), intent(in) :: self
    character(len=*), intent(in) :: val_name
    real :: norms_get_value

    !Should probably convert to lower case here, but until we
    !add some string utils to do this sort of stuff we'll rely
    !on developer doing the right thing.
    select case (trim(val_name))
    case("aref")
       norms_get_value=self%aref
    case("zref")
       norms_get_value=self%zref
    case("nref")
       norms_get_value=self%nref
    case("tref")
       norms_get_value=self%tref
    case("lref")
       norms_get_value=self%lref
    case("vref")
       norms_get_value=self%vref
    case("bref")
       norms_get_value=self%bref
    case("rhoref")
       norms_get_value=self%rhoref
    !case default
    !   call mp_abort("Invalid normalisation requested")
    end select
  end function norms_get_value

  subroutine norms_read_parameters(self)
    use file_utils, only: input_unit_exist
    use mp, only: proc0, broadcast
    implicit none
    class(norms_type), intent(in out) :: self
    real :: aref,zref,nref,tref,lref,vref,bref,rhoref
    integer :: in_file
    namelist/normalisations/aref,zref,nref,tref,lref,vref,bref,rhoref

    !Set defaults, check if namelist present and read
    if(proc0) then
       aref = self%def_val
       zref = self%def_val
       nref = self%def_val
       tref = self%def_val
       lref = self%def_val
       vref = self%def_val
       bref = self%def_val
       rhoref = self%def_val

       in_file = input_unit_exist("normalisations", exist)
       if(exist) read(unit=in_file, nml=normalisations)
    endif

    !Now broadcast values
    call broadcast(aref)
    call broadcast(zref)
    call broadcast(nref)
    call broadcast(tref)
    call broadcast(lref)
    call broadcast(vref)
    call broadcast(bref)
    call broadcast(rhoref)

    !Now copy parameters into holder type
    call self%set_value("aref",aref)
    call self%set_value("zref",zref)
    call self%set_value("nref",nref)
    call self%set_value("tref",tref)
    call self%set_value("lref",lref)
    call self%set_value("vref",vref)
    call self%set_value("bref",bref)
    call self%set_value("rhoref",rhoref)
  end subroutine norms_read_parameters

  !>Initialise the norms object
  subroutine norms_init(self)
    implicit none
    class(norms_type), intent(in out) :: self

    !First setup allowed normalisations
    self%nnorm=8
    if(.not.allocated(self%names)) allocate(self%names(self%nnorm))
    self%names(1)="aref"
    self%names(2)="zref"
    self%names(3)="nref"
    self%names(4)="tref"
    self%names(5)="lref"
    self%names(6)="vref"
    self%names(7)="bref"
    self%names(8)="rhoref"

    if(self%initialised) return
    self%initialised = .true.

    call self%read_parameters
  end subroutine norms_init

  !>Reset the properties
  subroutine norms_reset(self)
    implicit none
    class(norms_type), intent(in out) :: self
    integer :: i
    
    !Loop over parameters and set them to def_val
    do i=1,len(self%names)
       call self%set_value(self%names(i),self%def_val)
    enddo

    !Set the logical vars
    self%initialised=.false.
    call self%set_logicals
  end subroutine norms_reset

  !>Reset and free memory
  subroutine norms_finish(self)
    implicit none
    class(norms_type), intent(in out) :: self
    call self%reset
    if(allocated(self%names)) deallocate(self%names)
  end subroutine norms_finish

  !Decide if a given normalisation has been set
  elemental function norms_check_if_set(self,val_name)
    implicit none
    class(norms_type), intent(in) :: self
    character(len=*), intent(in) :: val_name
    logical :: norms_check_if_set
    norms_check_if_set = self%get_value(val_name).ne.self%def_val
  end function norms_check_if_set

  !Decide if all/some of the normalisations have been set
  subroutine norms_set_logicals(self)
    implicit none
    class(norms_type), intent(in out) :: self
    integer :: i
    logical :: some_set, all_set

    !Init internals
    some_set=.false.
    all_set=.true.

    !Loop over parameters and set them to def_val
    do i=1,len(self%names)
       all_set=all_set.and.self%check_if_set(self%names(i))
       some_set=some_set.or.self%check_if_set(self%names(i))
    enddo

    !Update object parameters
    self%some_norms_set=some_set
    self%all_norms_set=all_set
  end subroutine norms_set_logicals

  !/////////////////////////////
  !// MODULE LEVEL PROCDURES
  !/////////////////////////////
  subroutine check_normalisations(report_unit)
    implicit none
    integer, intent(in) :: report_unit
    character(len=7) :: msg
    if(norms%all_norms_set)then
       msg = "All of"
    else
       if(norms%some_norms_set)then
          msg = "Some of"
       else
          msg = "None of"
       endif
    endif

    write(report_unit,fmt='(A," the normalisation parameters have been provided.")') trim(msg)
  end subroutine check_normalisations

  subroutine wnml_normalisations(unit)
    implicit none
    integer, intent(in) :: unit
    integer :: i
    if(.not.exist) return

    write(unit,*)
    do i=1,len(norms%names)
       write(unit,fmt='(A," = ", F12.5)') trim(norms%names(i)),norms%get_value(norms%names(i))
    enddo
    write(unit,'("/")')
  end subroutine wnml_normalisations

  !>Read input file and populate the norms object
  subroutine init_normalisations
    implicit none
    call norms%init
  end subroutine init_normalisations

  !>Free memory etc. associated with normalisations
  subroutine finish_normalisations
    implicit none
    call norms%finish
  end subroutine finish_normalisations

  !>Reset normalisations object
  subroutine reset_normalisations
    implicit none
    call norms%reset
  end subroutine reset_normalisations
end module normalisations

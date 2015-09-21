!> A module for finding eigenvalues and eigenmodes
!! Requires SLEPC and PETSC
#ifdef WITH_EIG
module eigval
  use petscvec
  use petscmat
  use slepceps

  implicit none

!Allows us to user SLEPC_VERSION_MAJOR etc. to automatically
!disable unsupportted features
#include <slepcversion.h>
#include <finclude/petscvecdef.h>
#include <finclude/petscmatdef.h>
#include <finclude/slepcepsdef.h>

  private

  public :: eigval_functional
  public :: init_eigval, finish_eigval, BasicSolve
  public :: read_parameters, wnml_eigval, check_eigval
  public :: time_eigval

  !//////////////////
  !Slepc related configuration parameters
  !//////////////////
  !Solver type
  integer, parameter :: &        !Note not all of these may be available
       SolverTypePower=1,&
       SolverTypeSubspace=2,&
       SolverTypeArnoldi=3,&
       SolverTypeLanczos=4,&
       SolverTypeKrylovSchur=5,&
       SolverTypeGD=6,&
       SolverTypeJD=7,&
       SolverTypeRQCG=8,&
       SolverTypeCISS=9,&
       SolverTypeLapack=10,&
       SolverTypeArpack=11,&
       SolverTypeBlzpack=12,&
       SolverTypeTrlan=13,&
       SolverTypeBlopex=14,&
       SolverTypePrimme=15,&
       SolverTypeFeast=16,&
       SolverTypeNotSpecified=17 !=>Use slepc default
  integer :: solver_option_switch

  !Extraction type
  integer, parameter :: & !Note not all extractions are supported by every solver
       ExtractionRitz=1,&
       ExtractionHarmonic=2,&
       ExtractionHarmonicRelative=3,&
       ExtractionHarmonicRight=4,&
       ExtractionHarmonicLargest=5,&
       ExtractionRefined=6,&
       ExtractionRefinedHarmonic=7,&
       ExtractionNotSpecified=8 !=>Use slepc default
  integer :: extraction_option_switch

  !Which eigenpairs do we search for
  integer, parameter :: &
       WhichLargestMagnitude=1,&
       WhichSmallestMagnitude=2,&
       WhichLargestReal=3,&
       WhichSmallestReal=4,&
       WhichLargestImaginary=5,&
       WhichSmallestImaginary=6,&
       WhichTargetMagnitude=7,&
       WhichTargetReal=8,&
       WhichTargetImaginary=9,& !Complex builds only
       WhichAll=10,&   !Requires an interval to be set
       WhichUser=11,&  !Requires a used comparison routine to be defined
       WhichNotSpecified=12 !=>Use slepc default
  integer :: which_option_switch

  !What sort of spectral transform do we use
  integer, parameter :: &
       TransformShell=1,&
       TransformShift=2,&
       TransformInvert=3,&
       TransformCayley=4,&
       TransformFold=5,&
       TransformPrecond=6,&
       TransformNotSpecified=7
  integer :: transform_option_switch

  !Number of eigenvalues to seek
  integer :: n_eig

  !Maximum number of iterations (calls to advance) allowed
  integer :: max_iter

  !Tolerance to converge to
  real :: tolerance

  !Real and imaginary components of target (roughly where we expect eigvals to be)
  real :: targ_re, targ_im

  !Do we specify the initial condition (through ginit)?
  logical :: use_ginit

  !Number of GS2 advance steps per SLEPc call to advance_eigen
  !May be useful when looking at marginal modes etc.
  integer :: nadv

  !If true then save restart files for each eigenmode
  logical :: save_restarts

  !A custom type to make it easy to encapsulate specific settings
  type EpsSettings
     logical :: use_ginit
     integer :: n_eig
     integer :: max_iter
     integer :: solver_option_switch
     integer :: extraction_option_switch
     integer :: which_option_switch
     integer :: transform_option_switch
     PetscInt :: local_size, global_size
     real :: tolerance
     real :: targ_re, targ_im
     complex :: target
     PetscScalar :: target_slepc
  end type EpsSettings

  !//////////////////////
  !Operations parameters
  !//////////////////////
  logical, parameter :: allow_command_line_settings=.true.    
  character(len=12),parameter :: nml_name="eigval_knobs"
  real, save :: time_eigval(2)=0.

contains

  !>Returns true if GS2 was compiled with WITH_EIG defined
  function eigval_functional()
    logical :: eigval_functional
    eigval_functional=.true.
  end function eigval_functional
  
  !> Initialise eigenvalue problem
  subroutine init_eigval
    use mp, only: mp_comm, proc0
    use file_utils, only: error_unit
    implicit none
    PetscErrorCode :: ierr
#ifndef PETSC_USE_COMPLEX
    integer :: unit, ii
    integer, dimension(2) :: units
#endif

    !If we don't have eigenvalue support then ignore any settings in input file
    if(.not.eigval_functional())then
       return
    endif

    !If petsc has not been compiled with complex support then display
    !a warning that the results are most likely incorrect.
#ifndef PETSC_USE_COMPLEX
    units(1)=error_unit() !Error file
    units(2)=6            !Screen
    if(proc0)then
       do ii=1,2
          unit=units(ii)
          write(unit,'(66("#"))')
          write(unit,'("# ","WARNING : The PETSC library used does not have complex support"," #")')
          write(unit,'("# ","          --> The results are most likely incorrect.          "," #")')
          write(unit,'(66("#"))')
       enddo
    endif
#endif

    !Get the input parameters
    call read_parameters
    
    !Set PETSC_COMM_WORLD to be mp_comm so we can use whatever sub-comm
    !needed for list mode to work
    PETSC_COMM_WORLD=mp_comm

    !Initialise slepc
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  end subroutine init_eigval

  !> Clean up eigenvalue problem
  subroutine finish_eigval
    implicit none
    PetscErrorCode :: ierr

    !Finish up slepc
    call SlepcFinalize(ierr)
  end subroutine finish_eigval

  !> Read the eigval_knobs namelist
  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: broadcast, proc0
    implicit none

    !NOTE: text_option is case sensitive (currently)
    type(text_option), dimension(18), parameter :: solver_type_opts = &
         (/&
         text_option('slepc_default',SolverTypeNotSpecified),&
         text_option('default',SolverTypeKrylovSchur),&
         text_option('power',SolverTypePower),&
         text_option('subspace',SolverTypeSubspace),&
         text_option('arnoldi',SolverTypeArnoldi),&
         text_option('lanczos',SolverTypeLanczos),&
         text_option('krylov',SolverTypeKrylovSchur),&
         text_option('GD',SolverTypeGD),&
         text_option('JD',SolverTypeJD),&
         text_option('RQCG',SolverTypeRQCG),&
         text_option('CISS',SolverTypeCISS),&
         text_option('lapack',SolverTypeLapack),&
         text_option('arpack',SolverTypeArpack),&
         text_option('blzpack',SolverTypeBlzpack),&
         text_option('trlan',SolverTypeTrlan),&
         text_option('blopex',SolverTypeBlopex),&
         text_option('primme',SolverTypePrimme),&
         text_option('feast',SolverTypeFeast)&
         /)

    type(text_option), dimension(9), parameter :: extraction_type_opts = &
         (/&
         text_option('slepc_default',ExtractionNotSpecified),&
         text_option('default',ExtractionNotSpecified),&
         text_option('ritz',ExtractionRitz),&
         text_option('harmonic',ExtractionHarmonic),&
         text_option('harmonic_relative',ExtractionHarmonicRelative),&
         text_option('harmonic_right',ExtractionHarmonicRight),&
         text_option('harmonic_largest',ExtractionHarmonicLargest),&
         text_option('refined',ExtractionRefined),&
         text_option('refined_harmonic',ExtractionRefinedHarmonic)&
         /)

    type(text_option), dimension(13), parameter :: which_type_opts = &
         (/&
         text_option('slepc_default',WhichNotSpecified),&
         text_option('default',WhichTargetMagnitude),&
         text_option('largest_magnitude',WhichLargestMagnitude),&
         text_option('smallest_magnitude',WhichSmallestMagnitude),&
         text_option('largest_real',WhichLargestReal),&
         text_option('smallest_real',WhichSmallestReal),&
         text_option('largest_imaginary',WhichLargestImaginary),&
         text_option('smallest_imaginary',WhichSmallestImaginary),&
         text_option('target_magnitude',WhichTargetMagnitude),&
         text_option('target_real',WhichTargetReal),&
         text_option('target_imaginary',WhichTargetImaginary),&
         text_option('all',WhichAll),&
         text_option('user',WhichUser)&
         /)

    type(text_option), dimension(8), parameter :: transform_type_opts = &
         (/&
         text_option('slepc_default',TransformNotSpecified),&
         text_option('default',TransformNotSpecified),&
         text_option('shell',TransformShell),&
         text_option('shift',TransformShift),&
         text_option('invert',TransformInvert),&
         text_option('cayley',TransformCayley),&
         text_option('fold',TransformFold),&
         text_option('precond',TransformPrecond)&
         /)

    character(len=20) :: solver_option, extraction_option, which_option, transform_option

    namelist /eigval_knobs/ n_eig, max_iter, tolerance, &
         solver_option, extraction_option, which_option, transform_option,&
         targ_re, targ_im, use_ginit, nadv, save_restarts
    integer :: ierr, in_file
    logical :: exist

    if(proc0)then
       !Defaults
       n_eig=1
       max_iter=PETSC_DECIDE
       tolerance=1.0d-6
       targ_re=0.5
       targ_im=0.5
       solver_option='default'
       extraction_option='default'
       which_option='default'
       transform_option='default'
       use_ginit=.false.
       nadv=1
       save_restarts=.false.

       !Check if name list exists and if so read it
       in_file=input_unit_exist(nml_name,exist)
       if(exist) read(unit=in_file,nml=eigval_knobs)
       
       !Get error unit for output
       ierr = error_unit()
       
       !Convert string options to integer flags
       call get_option_value &
            (solver_option,solver_type_opts,solver_option_switch,&
            ierr,"solver_option in "//nml_name,.true.)
       
       call get_option_value &
            (extraction_option,extraction_type_opts,extraction_option_switch,&
            ierr,"extraction_option in "//nml_name,.true.)

       call get_option_value &
            (which_option,which_type_opts,which_option_switch,&
            ierr,"which_option in "//nml_name,.true.)

       call get_option_value &
            (transform_option,transform_type_opts,transform_option_switch,&
            ierr,"transform_option in "//nml_name,.true.)
    endif

    !Broadcasts
    call broadcast(n_eig)
    call broadcast(max_iter)
    call broadcast(tolerance)
    call broadcast(solver_option_switch)
    call broadcast(extraction_option_switch)
    call broadcast(which_option_switch)
    call broadcast(transform_option_switch)
    call broadcast(targ_re)
    call broadcast(targ_im)
    call broadcast(use_ginit)
    call broadcast(nadv)
    call broadcast(save_restarts)
  end subroutine read_parameters

  !> Write the eigval_knobs namelist
  subroutine wnml_eigval(unit)
    use mp, only: mp_abort
    implicit none
    integer, intent(in) :: unit
    character(len=4) :: inden="    "
    character(len=30) :: choice
    write(unit,*)
    write(unit,'(" &",A)') nml_name
    !Basic pars
    write(unit,'(A,A,"=",L1)') inden,'use_ginit',use_ginit
    write(unit,'(A,A,"=",I0)') inden,'n_eig',n_eig
    write(unit,'(A,A,"=",I0)') inden,'max_iter',max_iter
    write(unit,'(A,A,"=",I0)') inden,'nadv',nadv
    write(unit,'(A,A,"=",F12.5)') inden,'toleranace',tolerance
    write(unit,'(A,A,"=",F12.5)') inden,'targ_re',targ_re
    write(unit,'(A,A,"=",F12.5)') inden,'targ_im',targ_im
    !string pars
    !Solver
    select case(solver_option_switch)
    case(SolverTypePower)
       choice='power'
    case(SolverTypeSubspace)
       choice='subspace'
    case(SolverTypeArnoldi)
       choice='arnoldi'
    case(SolverTypeLanczos)
       choice='lanczos'
    case(SolverTypeKrylovSchur)
       choice='krylov'
    case(SolverTypeGD)
       choice='GD'
    case(SolverTypeJD)
       choice='JD'
    case(SolverTypeRQCG)
       choice='RQCG'
    case(SolverTypeCISS)
       choice='CISS'
    case(SolverTypeLapack)
       choice='lapack'
    case(SolverTypeArpack)
       choice='arpack'
    case(SolverTypeBlzpack)
       choice='blzpack'
    case(SolverTypeTrlan)
       choice='trlan'
    case(SolverTypeBlopex)
       choice='blopex'
    case(SolverTypePrimme)
       choice='primme'
    case(SolverTypeFeast)
       choice='feast'
    case(SolverTypeNotSpecified)
       choice='slepc_default'
    case default
       !Should never get here
       call mp_abort("Unknown value of solver_option_switch")
    end select
    write(unit,'(A,A,"=",A,A,A)') inden,'solver_option','"',choice,'"'

    !Extraction
    select case(extraction_option_switch)
    case(ExtractionRitz)
       choice='ritz'
    case(ExtractionHarmonic)
       choice='harmonic'
    case(ExtractionHarmonicRelative)
       choice='harmonic_relative'
    case(ExtractionHarmonicRight)
       choice='harmonic_right'
    case(ExtractionHarmonicLargest)
       choice='harmonic_largest'
    case(ExtractionRefined)
       choice='refined'
    case(ExtractionRefinedHarmonic)
       choice='refined_harmonic'
    case(ExtractionNotSpecified)
       choice='slepc_default'
    case default
       !Should never get here
       call mp_abort("Unknown value of extraction_option_switch")
    end select
    write(unit,'(A,A,"=",A,A,A)') inden,'extraction_option','"',choice,'"'
    
    !Which type
    select case(which_option_switch)
    case(WhichLargestMagnitude)
       choice='largest_magnitude'
    case(WhichSmallestMagnitude)
       choice='smallest_magnitude'
    case(WhichLargestReal)
       choice='largest_real'
    case(WhichSmallestReal)
       choice='smallest_real'
    case(WhichLargestImaginary)
       choice='largest_imaginary'
    case(WhichSmallestImaginary)
       choice='smallest_imaginary'
    case(WhichTargetMagnitude)
       choice='target_magnitude'
    case(WhichTargetReal)
       choice='target_real'
    case(WhichTargetImaginary)
       choice='target_imaginary'
    case(WhichAll)
       choice='all'
    case(WhichUser)
       choice='user'
    case(WhichNotSpecified)
       choice='slepc_default'
    case default
       !Should never get here
       call mp_abort("Unknown value of which_option_switch")
    end select
    write(unit,'(A,A,"=",A,A,A)') inden,'which_option','"',choice,'"'

    !Transform type
    select case(transform_option_switch)
    case(TransformShell)
       choice='shell'
    case(TransformShift)
       choice='shift'
    case(TransformInvert)
       choice='invert'
    case(TransformCayley)
       choice='cayley'
    case(TransformFold)
       choice='fold'
    case(TransformPrecond)
       choice='precond'
    case(TransformNotSpecified)
       choice='slepc_default'
    case default
       !Should never get here
       call mp_abort("Unknown value of transform_option_switch")
    end select
    write(unit,'(A,A,"=",A,A,A)') inden,'transform_option','"',choice,'"'

    !Done
    write(unit,'(" /")')

  end subroutine wnml_eigval

  !> Check the eigval settings
  subroutine check_eigval
  end subroutine check_eigval

  !> Create a default settings object
  subroutine InitSettings(eps_settings)
    use gs2_time, only: code_dt
    implicit none
    type(EpsSettings), intent(inout) :: eps_settings
    eps_settings%use_ginit=use_ginit
    eps_settings%n_eig=n_eig
    eps_settings%max_iter=max_iter
    eps_settings%solver_option_switch=solver_option_switch
    eps_settings%extraction_option_switch=extraction_option_switch
    eps_settings%which_option_switch=which_option_switch
    eps_settings%transform_option_switch=transform_option_switch
    eps_settings%tolerance=tolerance
    eps_settings%targ_re=targ_re
    eps_settings%targ_im=targ_im
    eps_settings%target=cmplx(targ_re,targ_im)
    !Convert GS2 eigenvalue to slepc eigenvalue
    eps_settings%target_slepc=Eig_Gs2ToSlepc(eps_settings%target)
    !Initialise the dimension attributes to -1 to indicate unset
    eps_settings%local_size=-1
    eps_settings%global_size=-1
  end subroutine InitSettings

  !> Setup the matrix operator
  subroutine InitMatrixOperator(mat_operator, Adv_routine, eps_settings)
    implicit none
    Mat, intent(inout) :: mat_operator
    external :: Adv_routine !This returns the result of mat_operator.X (where X is a vector)
    type(EpsSettings), intent(in) :: eps_settings
    PetscErrorCode :: ierr
    !Note, whilst we're defining a matrix (i.e. 2D array) we only specify
    !the local and global length of one side. This is because it is a shell
    !matrix where we only ever deal with vectors which can interact with the
    !matrix and never the matrix itself (hence why Sum(n_loc*n_loc)/=n_glob*n_glob
    !is ok!).

    !Make a shell matrix operator (AdvMat)
    call MatCreateShell(PETSC_COMM_WORLD,eps_settings%local_size,eps_settings%local_size,&
         eps_settings%global_size,eps_settings%global_size,PETSC_NULL_INTEGER,mat_operator,ierr)

    !Set the shell MATOP_MULT operation, i.e. which routine returns the result
    !of a AdvMat.x
    call MatShellSetOperation(mat_operator,MATOP_MULT,Adv_routine,ierr)
  end subroutine InitMatrixOperator

  !> Create a solver and set parameters based on eps_settings
  subroutine InitEPS(eps_solver,mat_operator,eps_settings)
    use init_g, only: ginit
    implicit none
    EPS, intent(inout) :: eps_solver !The solver object to initialise
    Mat, intent(in) :: mat_operator !The matrix to find eigenpairs of
    type(EpsSettings), intent(in) :: eps_settings
    Vec :: init_vec
    PetscErrorCode :: ierr
    logical :: restarted

    !Now create the eps solver
    call EPSCreate(PETSC_COMM_WORLD,eps_solver,ierr)

    !Set the matrix operator we want to find eigenvalues of
    call EPSSetOperators(eps_solver,mat_operator,PETSC_NULL_OBJECT,ierr)

    !Set the type of eigenproblem -->Always non-hermittian for us
    call EPSSetProblemType(eps_solver,EPS_NHEP,ierr)

    !Now setup from options
    call SetupEPS(eps_solver,eps_settings)

    !Now initialise if we requested this
    if(eps_settings%use_ginit)then
       !Initialise g
       call ginit(restarted)

       !Setup vector
       call MatGetVecs(mat_operator,PETSC_NULL_OBJECT,init_vec,ierr)
       call GnewToVec(init_vec)

       !Set the initial vector
       call EPSSetInitialSpace(eps_solver,1,init_vec,ierr)

       !Now destroy the vector
       call VecDestroy(init_vec,ierr)
    endif
  end subroutine InitEPS

  !> Setup the passed solver
  subroutine SetupEPS(eps_solver,eps_settings)
    use mp, only: mp_abort
    implicit none
    EPS, intent(inout) :: eps_solver
    type(EpsSettings), intent(in) :: eps_settings
    ST :: st
    PetscErrorCode :: ierr
    EPSType :: SolverType, TransformType
    PetscInt :: ExtractionType, WhichType

    !NOTE: There is no consistency checking here (yet) so it is up to the user
    !to make sure they request a valid combination of parameters

    !First set the number of eigenpairs to find
    call EPSSetDimensions(eps_solver,eps_settings%n_eig,&
         PETSC_DECIDE,PETSC_DECIDE,ierr)

    !Now set tolerance and iteration limits
    call EPSSetTolerances(eps_solver,eps_settings%tolerance,&
         eps_settings%max_iter,ierr)

    !Now set what type of eigenpairs we're looking for if we don't ask for slepc_default
    if(eps_settings%which_option_switch.ne.WhichNotSpecified)then
       select case(eps_settings%which_option_switch)
       case(WhichLargestMagnitude)
          WhichType=EPS_LARGEST_MAGNITUDE
       case(WhichSmallestMagnitude)
          WhichType=EPS_SMALLEST_MAGNITUDE
       case(WhichLargestReal)
          WhichType=EPS_LARGEST_REAL
       case(WhichSmallestReal)
          WhichType=EPS_SMALLEST_REAL
       case(WhichLargestImaginary)
          WhichType=EPS_LARGEST_IMAGINARY
       case(WhichSmallestImaginary)
          WhichType=EPS_SMALLEST_IMAGINARY
       case(WhichTargetMagnitude)
          WhichType=EPS_TARGET_MAGNITUDE
       case(WhichTargetReal)
          WhichType=EPS_TARGET_REAL
       case(WhichTargetImaginary)
          WhichType=EPS_TARGET_IMAGINARY
       case(WhichAll)
          WhichType=EPS_ALL
       case(WhichUser)
          WhichType=EPS_WHICH_USER
       case default
          !Should never get here
          call mp_abort("Unknown value of which_option_switch")
       end select

       call EpsSetWhichEigenpairs(eps_solver,WhichType,ierr)
    endif
    
    !Set the target
    call EpsSetTarget(eps_solver,eps_settings%target_slepc,ierr)

    !Now set the solver type if we don't ask for slepc_default
    if(eps_settings%solver_option_switch.ne.SolverTypeNotSpecified)then
       select case(eps_settings%solver_option_switch)
       case(SolverTypePower)
          SolverType=EPSPOWER
       case(SolverTypeSubspace)
          SolverType=EPSSUBSPACE
       case(SolverTypeArnoldi)
          SolverType=EPSARNOLDI
       case(SolverTypeLanczos)
          SolverType=EPSLANCZOS
       case(SolverTypeKrylovSchur)
          SolverType=EPSKRYLOVSCHUR
#if SLEPC_VERSION_GE(3,1,0)
       case(SolverTypeGD)
          SolverType=EPSGD
       case(SolverTypeJD)
          SolverType=EPSJD
#endif
#if SLEPC_VERSION_GE(3,3,0)
       case(SolverTypeRQCG)
          SolverType=EPSRQCG
#endif
#if SLEPC_VERSION_GE(3,4,0)
       case(SolverTypeCISS)
          SolverType=EPSCISS
#endif
       case(SolverTypeLapack)
          SolverType=EPSLAPACK
       case(SolverTypeArpack)
          SolverType=EPSARPACK
       case(SolverTypeBlzpack)
          SolverType=EPSBLZPACK
       case(SolverTypeTrlan)
          SolverType=EPSTRLAN
       case(SolverTypeBlopex)
          SolverType=EPSBLOPEX
       case(SolverTypePrimme)
          SolverType=EPSPRIMME
#if SLEPC_VERSION_GE(3,4,0)
       case(SolverTypeFeast)
          SolverType=EPSFEAST
#endif
       case default
          !Should never get here
          call mp_abort("Unknown value of solver_option_switch")
       end select

       call EPSSetType(eps_solver,SolverType,ierr)
    endif

    !Now set the extraction method
    if(eps_settings%extraction_option_switch.ne.ExtractionNotSpecified)then
       select case(eps_settings%extraction_option_switch)
       case(ExtractionRitz)
          ExtractionType=EPS_RITZ
       case(ExtractionHarmonic)
          ExtractionType=EPS_HARMONIC
       case(ExtractionHarmonicRelative)
          ExtractionType=EPS_HARMONIC_RELATIVE
       case(ExtractionHarmonicRight)
          ExtractionType=EPS_HARMONIC_RIGHT
       case(ExtractionHarmonicLargest)
          ExtractionType=EPS_HARMONIC_LARGEST
       case(ExtractionRefined)
          ExtractionType=EPS_REFINED
       case(ExtractionRefinedHarmonic)
          ExtractionType=EPS_REFINED_HARMONIC
       case default
          !Should never get here
          call mp_abort("Unknown value of extraction_option_switch")
       end select

       call EPSSetExtraction(eps_solver,ExtractionType,ierr)
    endif

    !Now set the spectral transform
    if(eps_settings%transform_option_switch.ne.TransformNotSpecified)then        
       !Get spectral transform object
       call EPSGetSt(eps_solver,st,ierr)

       select case(eps_settings%transform_option_switch)
       case(TransformShell)
          TransformType=STSHELL
       case(TransformShift)
          TransformType=STSHIFT
       case(TransformInvert)
          TransformType=STSINVERT
       case(TransformCayley)
          TransformType=STCAYLEY
#if SLEPC_VERSION_LE(3,4,4)
       case(TransformFold)
          TransformType=STFOLD
#endif
       case(TransformPrecond)
          TransformType=STPRECOND
       case default
          !Should never get here
          call mp_abort("Unknown value of transform_option_switch")
       end select

       !Set the shift
       call STSetShift(st,eps_settings%target_slepc,ierr)

       !Set the type
       call STSetType(st,TransformType,ierr)

       !Set the matrix mode to shell -- This is needed to use a lot of the
       !spectral transforms but is commented out for the moment because a
       !few other settings would need to be set including the KSP and PC
       !types as without implementing other matrix methods (like MATOP_GET_DIAGONAL)
       !the default KSP/PC won't work. To enable this at run time try something like
       ! -st_type sinvert -st_ksp_type gmres -st_pc_type none -st_matmode shell
       !call STSetMatMode(st,ST_MATMODE_SHELL,ierr)
    endif

    !Now we can allow users to override settings from command line if we like
    if(allow_command_line_settings)then
       !Do the solver object (also updates st object)
       call EPSSetFromOptions(eps_solver,ierr)
    endif
  end subroutine SetupEPS

  subroutine ReportSolverSettings(eps_solver,fext)
    use mp, only: proc0
    use file_utils, only: open_output_file, close_output_file
    implicit none
    EPS, intent(in) :: eps_solver
    character(len=*), intent(in), optional :: fext !Optional text to add to file name
    character(len=80) :: extension
    integer :: ounit
    EPSType :: TmpType !String
    EPSExtraction :: Extract !Integer
    PetscInt :: TmpInt
    PetscScalar :: TmpScal !Complex
    PetscReal :: TmpReal
    ST :: st
    PetscErrorCode :: ierr

    !Only proc0 does the writing
    if(.not.proc0) return

    !Set the default extension
    extension=".eig.solver_settings"

    !Add optional text is passed
    if(present(fext)) then
       extension=trim(fext)//trim(extension)
    endif

    !Open an output file for writing
    call open_output_file(ounit,extension)
    write(ounit,'("Slepc eigensolver settings:")')

    !Now we fetch settings and write them to file
#ifdef PETSC_USE_COMPLEX
    write(ounit,'("   ",A30)') "Compiled with COMPLEX support"
#else
    write(ounit,'("   ",A30)') "Compiled with REAL support"
#endif

    !1. Solver type
    call EPSGetType(eps_solver,TmpType,ierr)
    write(ounit,'("   ",A22,2X,":",2X,A22)') "Type", adjustl(TmpType)

    !2. Number of eigenvalues to look for
    call EPSGetDimensions(eps_solver,TmpInt,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)
    write(ounit,'("   ",A22,2X,":",2X,I0)') "Number of eigenvalues", TmpInt

    !3. Tolerances
    call EPSGetTolerances(eps_solver,TmpReal,TmpInt,ierr)
    write(ounit,'("   ",A22,2X,":",2X,I0)') "Max iterations", TmpInt
    write(ounit,'("   ",A22,2X,":",1X,E11.4)') "Tolerance", TmpReal

    !4. ExtractionType
    call EPSGetExtraction(eps_solver,Extract,ierr)
    !Could have a select case (in function) to convert integer to string name
    write(ounit,'("   ",A22,2X,":",2X,I0)') "Extraction type", Extract !Note integer, not string

    !5. WhichType
    call EPSGetWhichEigenPairs(eps_solver,TmpInt,ierr)
    !Could have a select case (in function) to convert integer to string name
    write(ounit,'("   ",A22,2X,":",2X,I0)') "Target type", TmpInt !Note integer, not string
    !5b. Target value
    call EPSGetTarget(eps_solver,TmpScal,ierr)
#ifdef PETSC_USE_COMPLEX
    write(ounit,'("   ",A22,2X,":",2X,F12.5)') "Real target (slepc)", PetscRealPart(TmpScal)
    write(ounit,'("   ",A22,2X,":",2X,F12.5)') "Imag target (slepc)", PetscImaginaryPart(TmpScal)
#else
    write(ounit,'("   ",A22,2X,":",2X,F12.5)') "Mag. target (slepc)", PetscRealPart(TmpScal)
#endif

    !6. Transform type
    call EPSGetST(eps_solver,st,ierr)
    call STGetType(st,TmpType,ierr)
    write(ounit,'("   ",A22,2X,":",2X,A22)') "Transform type", adjustl(TmpType)

    !Close output file
    call close_output_file(ounit)
  end subroutine ReportSolverSettings

  subroutine BasicSolve
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    implicit none
    EPS :: my_solver
    Mat :: my_operator
    PetscInt :: d1_size, d2_size, d3_size_global, d3_size_local
    PetscInt :: n_loc, n_glob
    type(EpsSettings) :: my_settings
    PetscErrorCode :: ierr

    !Define dimensions for problem
    d1_size=2*ntgrid+1 !Size of theta grid
    d2_size=2          !Size of sigma grid
    d3_size_local=(1+g_lo%ulim_alloc-g_lo%llim_proc) !Size of local distributed domain
    if(g_lo%ulim_proc.lt.g_lo%llim_proc) d3_size_local=0
    d3_size_global=(1+g_lo%ulim_world-g_lo%llim_world) !Size of global distributed domain

    !How big is the "advance operator matrix"
    n_loc  = d1_size*d2_size*d3_size_local
    n_glob = d1_size*d2_size*d3_size_global

    !Pack the settings type
    call InitSettings(my_settings)

    !Set the sizes
    my_settings%local_size=n_loc
    my_settings%global_size=n_glob

    !Initialise the matrix operator
    call InitMatrixOperator(my_operator, advance_eigen, my_settings)

    !Create the solver
    call InitEPS(my_solver,my_operator,my_settings)

    !Attach a monitor
    !call EPSMonitorSet(my_solver,SimpleMonitor,PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr)

    !Write settings to file
    call ReportSolverSettings(my_solver)

    !Now we're ready to solve
    call EPSSolve(my_solver,ierr)

    !Now do the reporting and diagnostic output
    call ReportResults(my_solver,my_operator,my_settings)

    !Now destroy objects
    call EPSDestroy(my_solver,ierr)
    call MatDestroy(my_operator,ierr)
  end subroutine BasicSolve

  !> Routine to write results to screen and netcdf
  subroutine ReportResults(my_solver,my_operator,my_settings)
    use mp, only: proc0
    use gs2_time, only: code_dt, user_time, user_dt
    use collisions, only: vnmult
    use run_parameters, only: fphi, fapar, fbpar
    use fields, only: set_init_fields
    use fields_arrays, only: phinew, bparnew
    use gs2_save, only: EigNetcdfID, init_eigenfunc_file, finish_eigenfunc_file
    use gs2_save, only: add_eigenpair_to_file, gs2_save_for_restart, finish_gs2_save
    use file_utils, only: run_name
    use dist_fn_arrays, only: gnew, g_adjust
    use gs2_diagnostics, only: save_distfn
    implicit none
    EPS, intent(in) :: my_solver
    Mat, intent(in) :: my_operator
    type(EpsSettings), intent(in) :: my_settings
    PetscErrorCode :: ierr
    PetscInt :: iteration_count, n_converged
    Vec :: eig_vec_r, eig_vec_i
    PetscScalar :: eig_val_r, eig_val_i
    complex :: EigVal
    logical :: all_found
    PetscInt :: ieig
    type(EigNetcdfID) :: io_ids
    character(len=20) :: restart_unique_string
    integer :: istatus

    !Find out how many iterations were performed
    call EPSGetIterationNumber(my_solver,iteration_count,ierr)

    !Find out how many converged eigenpairs where found
    call EPSGetConverged(my_solver,n_converged,ierr)
    all_found=(n_converged.ge.my_settings%n_eig)
    if(proc0) write(6,'("Found ",I0," eigenvalues in ",I0," iterations.")') n_converged,iteration_count

    if(n_converged.gt.0)then
       !Initialise the eigenvector arrays
       call MatGetVecs(my_operator,PETSC_NULL_OBJECT,eig_vec_r,ierr)
       call MatGetVecs(my_operator,PETSC_NULL_OBJECT,eig_vec_i,ierr)

       if(proc0)call init_eigenfunc_file(trim(run_name)//"_eig.out.nc",fphi,fapar,fbpar,io_ids)

       !Now loop over converged values
       do ieig=0,n_converged-1
          !Get eigenpair
          call EPSGetEigenpair(my_solver,ieig,eig_val_r,&
               eig_val_i,eig_vec_r,eig_vec_i,ierr)

          !NOTE: If petsc has been compiled with complex support then _i variables
          !are empty, else contains imaginary component
#ifdef PETSC_USE_COMPLEX
          EigVal=cmplx(PetscRealPart(eig_val_r),PetscImaginaryPart(eig_val_r))
#else
          EigVal=cmplx(PetscRealPart(eig_val_r),PetscRealPart(eig_val_i))
#endif
          !Convert to GS2 eigenvalue
          Eigval=Eig_SlepcToGs2(EigVal)

          !Get field eigenmodes
          !/First set g
#ifdef PETSC_USE_COMPLEX
          call VecToGnew(eig_vec_r)
#else
          call VecToGnew2(eig_vec_r,eig_vec_i)
#endif
          !/Then calculate fields
          call set_init_fields

          !Add to file
          if(proc0)call add_eigenpair_to_file(EigVal,fphi,fapar,fbpar,io_ids)

          !Save restart file
          if(save_restarts) then
             !First make the unique bit of the name
             write(restart_unique_string,'(".eig_",I0)') ieig

             !Save restart files
             call gs2_save_for_restart(gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar,.true.,fileopt=restart_unique_string)

             !Save distfn files
             if(save_distfn)then
                !Convert h to distribution function
                call g_adjust(gnew,phinew,bparnew,fphi,fbpar)
                
                !Save dfn, fields and velocity grids to file
                call gs2_save_for_restart(gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar,.true.,distfn=.true.,fileopt=restart_unique_string)

                !Convert distribution function back to h
                call g_adjust(gnew,phinew,bparnew,-fphi,-fbpar)
             endif

             !Reset the status of save_for_restart
             call finish_gs2_save
          endif
       enddo

       !Close netcdf file
       if(proc0)call finish_eigenfunc_file(io_ids)
 
       call VecDestroy(eig_vec_r,ierr)
       call VecDestroy(eig_vec_i,ierr)
    else
       if(proc0) write(6,'("No converged eigenvalues found.")')
    endif
  end subroutine ReportResults

  elemental complex function Eig_SlepcToGs2(eig)
    use gs2_time, only: code_dt
    implicit none
    complex, intent(in) :: eig
    Eig_SlepcToGs2=log(eig)*cmplx(0.0,1.0)/(code_dt*nadv)
  end function Eig_SlepcToGs2

  elemental PetscScalar function Eig_Gs2ToSlepc(eig)
    use gs2_time, only: code_dt
    implicit none
    complex, intent(in) :: eig
    !If PETSC doesn't have a complex type then return the
    !magnitude of the eigenvalue.
#ifdef PETSC_USE_COMPLEX
    Eig_Gs2ToSlepc=exp(-cmplx(0.0,1.0)*code_dt*eig*nadv)
#else
    Eig_Gs2ToSlepc=abs(exp(-cmplx(0.0,1.0)*code_dt*eig*nadv))
#endif
  end function Eig_Gs2ToSlepc

  subroutine advance_eigen(MatOperator,VecIn,Res,ierr)
    use fields, only: set_init_fields, advance
    implicit none
    Mat, intent(in) :: MatOperator
    Vec, intent(inout) :: VecIn, Res
    PetscErrorCode, intent(inout) :: ierr
    integer, parameter :: istep=2
    integer :: is

    !First unpack input vector into gnew
    call VecToGnew(VecIn)

    !Now set fields to be consistent with gnew
    call set_init_fields

    !Now do a number of time steps
    do is=1,nadv
       !Note by using a fixed istep we
       !force the explicit terms to be excluded
       !except for the very first call.
       call advance(istep)
    enddo

    !Now pack gnew into output
    call GnewToVec(Res)
  end subroutine advance_eigen

  subroutine VecToGnew(VecIn)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use dist_fn_arrays, only: gnew
    implicit none
    Vec, intent(inout) :: VecIn
    PetscScalar, pointer :: VecInArr(:)
    PetscInt :: d1_size, d2_size, d3_size_global, d3_size_local
    PetscErrorCode :: ierr
    integer :: ig, isgn, iglo, local_index

    !No local data
    if(g_lo%ulim_proc.lt.g_lo%llim_proc) return

    !Define dimensions
    d1_size=2*ntgrid+1 !Size of theta grid
    d2_size=2          !Size of sigma grid
    d3_size_local=(1+g_lo%ulim_alloc-g_lo%llim_proc) !Size of local distributed domain
    d3_size_global=(1+g_lo%ulim_world-g_lo%llim_world) !Size of global distributed domain

    !Get a pointer to the data
    call VecGetArrayF90(VecIn,VecInArr,ierr)

    !Extract
    do iglo=g_lo%llim_proc,g_lo%ulim_alloc
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             !Form local index (note we could just having a running counter which we
             !increment on each loop
             local_index=1+((iglo-g_lo%llim_proc)*d2_size*d1_size)+(isgn-1)*d1_size+(ig+ntgrid)
             gnew(ig,isgn,iglo)=VecInArr(local_index)
          enddo
       enddo
    enddo

    !Get rid of pointer
    call VecRestoreArrayF90(VecIn,VecInArr,ierr)
  end subroutine VecToGnew

#ifndef PETSC_USE_COMPLEX
  subroutine VecToGnew2(VecIn,VecInI)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use dist_fn_arrays, only: gnew
    implicit none
    Vec, intent(inout) :: VecIn,VecInI
    PetscScalar, pointer :: VecInArr(:),VecInIArr(:)
    PetscInt :: d1_size, d2_size, d3_size_global, d3_size_local
    PetscErrorCode :: ierr
    integer :: ig, isgn, iglo, local_index

    !No local data
    if(g_lo%ulim_proc.lt.g_lo%llim_proc) return

    !Define dimensions
    d1_size=2*ntgrid+1 !Size of theta grid
    d2_size=2          !Size of sigma grid
    d3_size_local=(1+g_lo%ulim_alloc-g_lo%llim_proc) !Size of local distributed domain
    d3_size_global=(1+g_lo%ulim_world-g_lo%llim_world) !Size of global distributed domain

    !Get a pointer to the data
    call VecGetArrayF90(VecIn,VecInArr,ierr)
    call VecGetArrayF90(VecInI,VecInIArr,ierr)

    !Extract
    do iglo=g_lo%llim_proc,g_lo%ulim_alloc
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             !Form local index (note we could just having a running counter which we
             !increment on each loop
             local_index=1+((iglo-g_lo%llim_proc)*d2_size*d1_size)+(isgn-1)*d1_size+(ig+ntgrid)
             gnew(ig,isgn,iglo)=cmplx(PetscRealPart(VecInArr(local_index)),PetscRealPart(VecInIArr(local_index)))
          enddo
       enddo
    enddo

    !Get rid of pointer
    call VecRestoreArrayF90(VecIn,VecInArr,ierr)
    call VecRestoreArrayF90(VecInI,VecInIArr,ierr)
  end subroutine VecToGnew2
#endif
 
  subroutine GnewToVec(VecIn)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use dist_fn_arrays, only: gnew
    implicit none
    Vec, intent(inout) :: VecIn
    PetscScalar, pointer :: VecInArr(:)
    PetscInt :: d1_size, d2_size, d3_size_global, d3_size_local
    PetscErrorCode :: ierr
    integer :: ig, isgn, iglo, local_index

    !No local data
    if(g_lo%ulim_proc.lt.g_lo%llim_proc) return

    !Define dimensions
    d1_size=2*ntgrid+1 !Size of theta grid
    d2_size=2          !Size of sigma grid
    d3_size_local=(1+g_lo%ulim_alloc-g_lo%llim_proc) !Size of local distributed domain
    d3_size_global=(1+g_lo%ulim_world-g_lo%llim_world) !Size of global distributed domain

    !Get a pointer to the data
    call VecGetArrayF90(VecIn,VecInArr,ierr)

    !Extract
    do iglo=g_lo%llim_proc,g_lo%ulim_alloc
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             !Form local index (note we could just having a running counter which we
             !increment on each loop
             local_index=1+((iglo-g_lo%llim_proc)*d2_size*d1_size)+(isgn-1)*d1_size+(ig+ntgrid)
             VecInArr(local_index)=gnew(ig,isgn,iglo)
          enddo
       enddo
    enddo

    !Get rid of pointer
    call VecRestoreArrayF90(VecIn,VecInArr,ierr)
  end subroutine GnewToVec

  PetscErrorCode function SimpleMonitor(solver,iters,nconv,eig_r,eig_i,err_est,n_est,mm)
    use mp, only: proc0
    use gs2_time, only: code_dt
    implicit none
    EPS, intent(inout) :: solver
    PetscInt, intent(inout) :: iters, nconv, n_est
    PetscScalar,dimension(:), intent(inout) :: eig_r, eig_i
    PetscReal, dimension(:), intent(inout) :: err_est
    integer, pointer, intent(inout), optional :: mm
    integer, save :: nconv_prev=0, iter_prev=0
    integer :: new_conv, new_iter

    !Set return value
    SimpleMonitor=0

    !If no change in nconv then return
    if(nconv_prev.eq.nconv)then
       return
    endif

    !Calculate how many extra modes found and how many iterations it took
    new_conv=nconv-nconv_prev
    new_iter=iters-iter_prev
    if(proc0)write(6,'("Found ",I0," eigenvalues in ",I0," iterations:")') new_conv, new_iter

    !Update previous value
    nconv_prev=nconv
    iter_prev=iters
  end function SimpleMonitor
end module eigval

#else

!###############################################################
!THE SECTION BELOW PROVIDES THE EIGENVAL MODULE IN THE ABSENCE OF
!SLEPC/PETSC (i.e. WITH_EIG not defined).
module eigval
  implicit none
  private
  public :: eigval_functional
contains
  !>Returns true if GS2 was compiled with WITH_EIG defined
  function eigval_functional()
    logical :: eigval_functional
    eigval_functional=.false.
  end function eigval_functional
end module eigval

#endif

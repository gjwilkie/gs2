# include "define.inc"

MODULE gs2_save

!==============================================================================================
!MODULE IMPORTS
!==============================================================================================
# ifdef NETCDF
!   use netcdf, only: NF90_FLOAT, NF90_DOUBLE
    USE netcdf, ONLY: NF90_NOWRITE, NF90_CLOBBER, NF90_NOERR
    USE netcdf, ONLY: nf90_create, nf90_open, nf90_sync, nf90_close
    USE netcdf, ONLY: nf90_def_dim, nf90_def_var, nf90_enddef
    USE netcdf, ONLY: nf90_put_var, nf90_get_var, nf90_strerror
    USE netcdf, ONLY: nf90_inq_dimid, nf90_inquire_dimension
    USE netcdf, ONLY: nf90_inq_varid, nf90_inquire_variable
    USE netcdf_utils, ONLY: get_netcdf_code_precision
    USE netcdf_utils, ONLY: check_netcdf_file_precision
    USE netcdf_utils, ONLY: netcdf_error
    USE netcdf_utils, ONLY: netcdf_real, kind_nf
# endif
!----------------------------------------------------------------------------------------------

!==============================================================================================
!VARIABLE DECLARATIONS
!==============================================================================================
    !No implicits
    IMPLICIT NONE

    !Public access definitions
    PUBLIC :: gs2_restore, gs2_save_for_restart
    PUBLIC :: init_save, init_dt, init_tstart, init_ant_amp
    PUBLIC :: init_vnm
    !<DD 27-08-2010> Make save_distfn public
    PUBLIC :: gs2_save_distfn
    !</DD>

    !Interface to restore routines
    INTERFACE gs2_restore
        MODULE PROCEDURE gs2_restore_many, gs2_restore_one
    END INTERFACE

    !Internal variables
    PRIVATE

    CHARACTER (300), SAVE :: restart_file

# ifdef NETCDF
    !Arrays used for temporaries
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: tmpr, tmpi, ftmpr, ftmpi
    REAL, ALLOCATABLE, DIMENSION(:) :: stmp  ! MR: tmp var for kx_shift
    REAL, ALLOCATABLE, DIMENSION(:) :: atmp

    !Netcdf identifiers
    INTEGER (KIND_NF) :: ncid, thetaid, signid, gloid, kyid, kxid, nk_stir_dim
    INTEGER (KIND_NF) :: phir_id, phii_id, aparr_id, apari_id, bparr_id, bpari_id
    INTEGER (KIND_NF) :: kx_shift_id   ! MR: added to save kx_shift variable
    INTEGER (KIND_NF) :: t0id, gr_id, gi_id, vnm1id, vnm2id, delt0id
    INTEGER (KIND_NF) :: a_antr_id, b_antr_id, a_anti_id, b_anti_id
    !<DD 29-08-2010> Netcdf id defintions for energy and lambda dimensions/variables
    INTEGER (KIND_NF) :: egridid,lgridid
    INTEGER (KIND_NF) :: energy_id, lambda_id
    !</DD>
    !<DD 01-09-2010> Netcdf id definitions for nspecies dimensions/variables
    INTEGER (KIND_NF) :: nspecid, spec_id
    !</DD>
    !<DD 02-09-2010> Netcdf id definitions for vparallel and vperpendicular variables
    INTEGER (KIND_NF) :: vpa_id, vperp2_id
    !</DD>

    !Switches
    LOGICAL :: initialized = .false.
    !<DD 04-09-2010> Added initialized_dfn switch to determine if distfn file
    !has been initialized
    LOGICAL :: initialized_dfn = .false.
    !</DD>
    LOGICAL :: test = .false.
# endif
!----------------------------------------------------------------------------------------------

 CONTAINS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE gs2_save_for_restart &
    (g, t0, delt0, vnm, istatus, fphi, fapar, fbpar, exit_in)
!GS2_SAVE_FOR_RESTART(gs2_save.fpp): Subroutine to save internally used
!distribution function(h) and EM fields (as well as a couple of other
!variable) to a restart file (1 for each processor).

!==============================================================================================
!MODULE IMPORTS
!==============================================================================================
# ifdef NETCDF
    USE constants, ONLY: kind_rs, kind_rd, pi
    USE fields_arrays, ONLY: phinew, aparnew, bparnew
    USE dist_fn_arrays, ONLY: kx_shift  !MR
    USE kt_grids, ONLY: naky, ntheta0
    USE antenna_data, ONLY: nk_stir, a_ant, b_ant, ant_on
# endif
    USE mp, ONLY: nproc, iproc, proc0
    USE theta_grid, ONLY: ntgrid
! Must include g_layout_type here to avoid obscure bomb while compiling
! gs2_diagnostics.f90 (which uses this module) with the Compaq F90 compiler:
! TT>
    USE gs2_layouts, ONLY: g_lo
    USE layouts_type, ONLY: g_layout_type
! <TT
    USE file_utils, ONLY: error_unit
!----------------------------------------------------------------------------------------------

!==============================================================================================
!VARIABLE DECLARATIONS
!==============================================================================================
    !No implicits
    IMPLICIT NONE

    !Interface variables
    COMPLEX, DIMENSION (-ntgrid:,:,g_lo%llim_proc:), INTENT (IN) :: g
    REAL, INTENT (IN) :: t0, delt0
    REAL, DIMENSION (2), INTENT (IN) :: vnm
    REAL, INTENT (IN) :: fphi, fapar, fbpar
    INTEGER, INTENT (OUT) :: istatus
    LOGICAL, INTENT (IN), OPTIONAL :: exit_in

# ifdef NETCDF
    !Local variables
    CHARACTER (306) :: file_proc
    CHARACTER (10) :: suffix
    INTEGER :: i, n_elements, ierr
    LOGICAL :: exit
!----------------------------------------------------------------------------------------------

!==============================================================================================
!START OF WORK SECTION
!==============================================================================================
    !Decide what to do with file at end of routine
    IF (PRESENT(exit_in)) THEN
       exit = exit_in
    ELSE
       exit = .false.
    END IF

    !Get the number of elements for this processor
    n_elements = g_lo%ulim_proc-g_lo%llim_proc+1

    !If no elements then nothing to write so escape routine
    IF (n_elements <= 0) RETURN

    !#####################
    !#File initialisation#
    !#####################
    !If file not initialised then create file, define dimensions and variables
    IF (.not.initialized) THEN
        !Don't perform initialisation twice
        initialized = .true.

        !Define file prefix
        file_proc = TRIM(restart_file)

        !Define file suffix (.procnum)
        WRITE (suffix,'(a1,i0)') '.', iproc

        !Define full file name
        file_proc = TRIM(TRIM(file_proc)//ADJUSTL(suffix))

        !Create file
        istatus = nf90_create (file_proc, NF90_CLOBBER, ncid)

        !Check file was created successfully
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_create error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !###################
        !#Define dimensions#
        !###################
        !<note> this if statement shouldn't be needed due to earlier check which
        !results in subroutine exit if n_elements <=0 </note>
        IF (n_elements > 0) THEN
            !Define theta dimension
            istatus = nf90_def_dim (ncid, "theta", 2*ntgrid+1, thetaid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim theta error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define sign dimension (sign of v_par)
            istatus = nf90_def_dim (ncid, "sign", 2, signid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim sign error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define glo dimension (corresponds to layout grid i.e. kx,ky,lambda,energy,species mapped to 1d array)
            istatus = nf90_def_dim (ncid, "glo", n_elements, gloid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim glo error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define aky dimension (number of ky modes)
            istatus = nf90_def_dim (ncid, "aky", naky, kyid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim aky error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define akx dimension (number of kx modes)
            istatus = nf90_def_dim (ncid, "akx", ntheta0, kxid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim akx error: ", nf90_strerror(istatus)
                GOTO 1
            END IF
        END IF !End of IF (n_elements > 0) THEN

        !###################
        !#Define variables #
        !###################

        !Get real variable data type specifier for netcdf routines
        IF (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()

        !Define t0 variable (user time)
        istatus = nf90_def_var (ncid, "t0", netcdf_real, t0id)

        !Check variable created successfully
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_def_var t0 error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !Define delt0 variable (time step)
        istatus = nf90_def_var (ncid, "delt0", netcdf_real, delt0id)

        !Check variable created successfully
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_def_var delt0 error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !Define vnm1 variable (??)
        istatus = nf90_def_var (ncid, "vnm1", netcdf_real, vnm1id)

        !Check variable created successfully
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_def_var vnm(1) error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !Define vnm2 variable (??)
        istatus = nf90_def_var (ncid, "vnm2", netcdf_real, vnm2id)

        !Check variable created successfully
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_def_var vnm(2) error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !Create antenna related variables if ant_on is true
        IF (ant_on) THEN
            istatus = nf90_def_dim (ncid, "nk_stir", nk_stir, nk_stir_dim)
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim nk_stir error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            istatus = nf90_def_var (ncid, "a_ant_r", netcdf_real, nk_stir_dim, a_antr_id)
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var a_ant_r error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            istatus = nf90_def_var (ncid, "a_ant_i", netcdf_real, nk_stir_dim, a_anti_id)
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var a_ant_i error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            istatus = nf90_def_var (ncid, "b_ant_r", netcdf_real, nk_stir_dim, b_antr_id)
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var b_ant_r error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            istatus = nf90_def_var (ncid, "b_ant_i", netcdf_real, nk_stir_dim, b_anti_id)
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var b_ant_i error: ", nf90_strerror(istatus)
                GOTO 1
            END IF
        END IF

        !<note> this if statement shouldn't be needed due to earlier check which
        !results in subroutine exit if n_elements <=0 </note>
        IF (n_elements > 0) THEN
            !Define gr variable (real component of g)
            istatus = nf90_def_var (ncid, "gr", netcdf_real, &
                (/ thetaid, signid, gloid /), gr_id)

            !Check variable created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var g error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define gi variable (imaginary component of gi)
            istatus = nf90_def_var (ncid, "gi", netcdf_real, &
                (/ thetaid, signid, gloid /), gi_id)

            !Check variable created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var g error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !If there are finite electrostatic potential perturbations
            !then write them to file
            IF (fphi > epsilon(0.)) THEN
                !Define phi_r variable (real component of \phi)
                istatus = nf90_def_var (ncid, "phi_r", netcdf_real, &
                    (/ thetaid, kxid, kyid /), phir_id)

                !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var phi error: ", nf90_strerror(istatus)
                    GOTO 1
                END IF

                !Define phi_i variable (imaginary component of \phi)
                istatus = nf90_def_var (ncid, "phi_i", netcdf_real, &
                    (/ thetaid, kxid, kyid /), phii_id)

                 !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var phi error: ", nf90_strerror(istatus)
                    GOTO 1
                END IF
            END IF

            !If there are finite parallel magnetic potential perturbations
            !then write them to file
            IF (fapar > epsilon(0.)) THEN
                !Define apar_r variable (real component A_\parallel)
                istatus = nf90_def_var (ncid, "apar_r", netcdf_real, &
                    (/ thetaid, kxid, kyid /), aparr_id)

                !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var apar error: ", nf90_strerror(istatus)
                    GOTO 1
                END IF

                !Define apar_i variable (imaginary component A_\parallel)
                istatus = nf90_def_var (ncid, "apar_i", netcdf_real, &
                    (/ thetaid, kxid, kyid /), apari_id)

                !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var apar error: ", nf90_strerror(istatus)
                    GOTO 1
                END IF
            END IF

            !If there are finite perpendicular magnetic potential perturbations
            !then write them to file
            IF (fbpar > epsilon(0.)) THEN
                !Define bpar_r variable (real component B_\parallel)
                istatus = nf90_def_var (ncid, "bpar_r", netcdf_real, &
                    (/ thetaid, kxid, kyid /), bparr_id)

                !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var bparr error: ", nf90_strerror(istatus)
                    GOTO 1
                END IF

                !Define bpar_i variable (imaginary component B_\parallel)
                istatus = nf90_def_var (ncid, "bpar_i", netcdf_real, &
                    (/ thetaid, kxid, kyid /), bpari_id)

                !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var bpari error: ", nf90_strerror(istatus)
                    GOTO 1
                END IF
            END IF
        END IF !End of IF (n_elements > 0) THEN

        !Define kx_shift variable (gives shift in kx due to flow shear)
        istatus = nf90_def_var (ncid, "kx_shift", netcdf_real, &
            (/ kyid /), kx_shift_id)

        !Check variable created correctly
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_def_var kx_shift error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !######################
        !#End definition stage#
        !######################

        !End definition stage of file
        istatus = nf90_enddef (ncid)

        !Check end was successful
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE (ierr,*) "nf90_enddef error: ", nf90_strerror(istatus)
            GOTO 1
        END IF
    END IF !End of IF (.not. initialized) THEN

    !#################
    !#File population#
    !#################

    !#################
    !#Store variables#
    !#################

    !Store delt0
    istatus = nf90_put_var (ncid, delt0id, delt0)

    !Check store was successful
    IF (istatus /= NF90_NOERR) THEN
        ierr = error_unit()
        WRITE (ierr,*) "nf90_put_var delt0 error: ", nf90_strerror(istatus)
        GOTO 1
    END IF

    !Store t0
    istatus = nf90_put_var (ncid, t0id, t0)

    !Check store was successful
    IF (istatus /= NF90_NOERR) THEN
        ierr = error_unit()
        WRITE (ierr,*) "nf90_put_var t0 error: ", nf90_strerror(istatus)
        GOTO 1
    END IF

    !Store vnm(1)
    istatus = nf90_put_var (ncid, vnm1id, vnm(1))

    !Check store was successful
    IF (istatus /= NF90_NOERR) THEN
        ierr = error_unit()
        WRITE (ierr,*) "nf90_put_var vnm(1) error: ", nf90_strerror(istatus)
        GOTO 1
    END IF

    !Store vnm(2)
    istatus = nf90_put_var (ncid, vnm2id, vnm(2))

    !Check store was successful
    IF (istatus /= NF90_NOERR) THEN
        ierr = error_unit()
        WRITE (ierr,*) "nf90_put_var vnm(2) error: ", nf90_strerror(istatus)
        GOTO 1
    END IF

    !Escape point in case of netcdf error at any point
1   CONTINUE

    !If an error has been encountered then close the file and exit routine
    IF (istatus /= NF90_NOERR) THEN
       i = nf90_close (ncid)
       RETURN
    END IF

    !<note> this if statement shouldn't be needed due to earlier check which
    !results in subroutine exit if n_elements <=0 </note>
    IF (n_elements > 0) THEN

        !Fill antenna variables if ant_on set
        IF (ant_on) THEN
            IF (.not. ALLOCATED(atmp)) ALLOCATE (atmp(nk_stir))
            atmp = REAL(a_ant)
            istatus = nf90_put_var (ncid, a_antr_id, atmp)

            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE (ierr,*) "nf90_put_var a_antr error: ", &
                    nf90_strerror(istatus), ' ', iproc
            END IF

            atmp = AIMAG(a_ant)
            istatus = nf90_put_var (ncid, a_anti_id, atmp)

            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE (ierr,*) "nf90_put_var a_anti error: ", &
                    nf90_strerror(istatus), ' ', iproc
            END IF

            atmp = REAL(b_ant)
            istatus = nf90_put_var (ncid, b_antr_id, atmp)
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE (ierr,*) "nf90_put_var b_antr error: ", &
                    nf90_strerror(istatus), ' ', iproc
            END IF

            atmp = AIMAG(b_ant)
            istatus = nf90_put_var (ncid, b_anti_id, atmp)
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE (ierr,*) "nf90_put_var b_anti error: ", &
                    nf90_strerror(istatus), ' ', iproc
            END IF

            DEALLOCATE (atmp)
        END IF

        !Allocate array to hold components of g individually
        IF (.not. ALLOCATED(tmpr)) &
            ALLOCATE (tmpr(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))

        !Make temporary copy of the real component of g
        tmpr = REAL(g)

        !Store variable real(g)
        istatus = nf90_put_var (ncid, gr_id, tmpr)

        !Check store was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, gr_id)

        !Make temporary copy of the imaginary component of g
        tmpr = AIMAG(g)

        !Store variable aimag(g)
        istatus = nf90_put_var (ncid, gi_id, tmpr)

        !Check store was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, gi_id)

        !Allocate arrays to hold components of fields individually
        IF (.not. ALLOCATED(ftmpr)) ALLOCATE (ftmpr(2*ntgrid+1,ntheta0,naky))
        IF (.not. ALLOCATED(ftmpi)) ALLOCATE (ftmpi(2*ntgrid+1,ntheta0,naky))

        !If finite electrostatic potential perturbation then write \phi
        IF (fphi > epsilon(0.)) THEN
            !Nake temporary copy of real(\phi)
            ftmpr = REAL(phinew)

            !Store variable real(\phi)
            istatus = nf90_put_var (ncid, phir_id, ftmpr)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, phir_id)

            !Make temporary copy of imaginary(\phi)
            ftmpi = AIMAG(phinew)

            !Store variable aimag(\phi)
            istatus = nf90_put_var (ncid, phii_id, ftmpi)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, phii_id)
        END IF

        !If finite parallel magnetic potential perturbations then write A_\parallel
        IF (fapar > epsilon(0.)) THEN
            !Make temporary copy of real(A_\parallel)
            ftmpr = REAL(aparnew)

            !Store variable real(A_\parallel)
            istatus = nf90_put_var (ncid, aparr_id, ftmpr)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, aparr_id)

            !Make temporary copy of imaginary(A_\parallel)
            ftmpi = AIMAG(aparnew)

            !Store variable imaginary(A_\parallel)
            istatus = nf90_put_var (ncid, apari_id, ftmpi)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, apari_id)
        END IF

        !If finite perpendicular magnetic potential perturbations then write B_\parallel
        IF (fbpar > epsilon(0.)) THEN
            !Make temporary copy of real(B_\parallel)
            ftmpr = REAL(bparnew)

            !Store variable real(B_\parallel)
            istatus = nf90_put_var (ncid, bparr_id, ftmpr)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, bparr_id)

            !Make temporary copy of imaginary(B_\parallel)
            ftmpi = AIMAG(bparnew)

            !Store variable imaginary(B_\parallel)
            istatus = nf90_put_var (ncid, bpari_id, ftmpi)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, bpari_id)
        END IF
    END IF !End of IF (n_elements > 0) THEN

    !If kx_shift is allocated then output kx otherwise output zeros
    IF (allocated(kx_shift)) THEN ! MR begin
        !Make temporary array to hold kx_shift if allocated, zeros if not
        IF (.not. ALLOCATED(stmp)) ALLOCATE (stmp(naky))

        !Fill in output array
        stmp = kx_shift

        !Store variable
        istatus = nf90_put_var (ncid, kx_shift_id, stmp)

        !Check store was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, kx_shift_id)
    ELSE
        !Make temporary array to hold kx_shift if allocated, zeros if not
        IF (.not. ALLOCATED(stmp)) ALLOCATE (stmp(naky))

        !Fill in output array
        stmp = 0.

        !Store variable
        istatus = nf90_put_var (ncid, kx_shift_id, stmp)

        !Check store was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, kx_shift_id)
    END IF ! MR end

    !If exit_in flag set then close the file, if not then sync the file
    IF (exit) THEN
        !Close file
        i = nf90_close (ncid)
    ELSE
        !Sync file
        i = nf90_sync (ncid)

        !Check sync was successful
        IF (i /= NF90_NOERR) &
            CALL netcdf_error (istatus, message='nf90_sync error')
    END IF

# else
    !Statement to display if try to save distfn without netcdf available
    IF (proc0) WRITE (error_unit(),*) &
        'WARNING: gs2_save_for_restart is called without netcdf library'

# endif
!----------------------------------------------------------------------------------------------
END SUBROUTINE gs2_save_for_restart
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE gs2_restore_many (g, scale, istatus, fphi, fapar, fbpar, many)
!GS2_RESTORE_MANY(gs2_save.fpp): Subroutine to restore distribution
!function, EM fields and a couple of other variables from a restart
!file on each processor.
!==============================================================================================
!MODULE IMPORTS
!==============================================================================================
!MR, 2007: restore kx_shift array if already allocated
# ifdef NETCDF
    USE mp, ONLY: proc0, iproc, nproc
    USE fields_arrays, ONLY: phinew, aparnew, bparnew
    USE fields_arrays, ONLY: phi, apar, bpar
    USE dist_fn_arrays, ONLY: kx_shift   ! MR
    USE kt_grids, ONLY: naky, ntheta0
# endif
    USE theta_grid, ONLY: ntgrid
    USE gs2_layouts, ONLY: g_lo
    USE file_utils, ONLY: error_unit
!----------------------------------------------------------------------------------------------

!==============================================================================================
!VARIABLE DECLARATIONS
!==============================================================================================
    !No implicits
    IMPLICIT NONE

    !Interface variables
    COMPLEX, DIMENSION (-ntgrid:,:,g_lo%llim_proc:), INTENT (OUT) :: g
    REAL, INTENT (IN) :: scale
    INTEGER, INTENT (OUT) :: istatus
    REAL, INTENT (IN) :: fphi, fapar, fbpar
    LOGICAL, INTENT (IN) :: many
# ifdef NETCDF

    !Internal variables
    CHARACTER (306) :: file_proc
    CHARACTER (10) :: suffix
    INTEGER :: i, n_elements, ierr
    REAL :: fac
!----------------------------------------------------------------------------------------------

!==============================================================================================
!WORK SECTION
!==============================================================================================

    !Get number of elements on this processor
    n_elements = g_lo%ulim_proc-g_lo%llim_proc+1

    !If no elements then nothing to read so escape routine
    IF (n_elements <= 0) RETURN

    !Check that the file has not been initialised already
    IF (.not.initialized) THEN
        !######################
        !#File Initialisations#
        !######################

        !Set initialised state to True so only do section once
!       initialized = .true.

        !Get file prefix
        file_proc = TRIM(restart_file)

        !Define file suffix (*.iproc)
        WRITE (suffix,'(a1,i0)') '.', iproc

        !Define full file name
        file_proc = TRIM(TRIM(file_proc)//ADJUSTL(suffix))

        !Open specified file for reading
        istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)

        !Check open was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, file=file_proc)

        !Load netcdf precision variable if not done already
        IF (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()

        !Call check_netcdf_file_precision
        !   Opens file and tries to read a variable 't0' (typically representing time?)
        !   Then checks preicision of data in file and compares this to precision currently
        !   being used in the code.
        !   If these don't match then will print a warning to screen
        CALL check_netcdf_file_precision (ncid)

        !######################
        !#Obtain dimension ids#
        !######################

        !Obtain id for dimension theta
        istatus = nf90_inq_dimid (ncid, "theta", thetaid)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, dim='theta')

        !Obtain id for dimension sign
        istatus = nf90_inq_dimid (ncid, "sign", signid)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, dim='sign')

        !Obtain id for dimension glo
        istatus = nf90_inq_dimid (ncid, "glo", gloid)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, dim='glo')

        !Obtain id for dimension aky
        istatus = nf90_inq_dimid (ncid, "aky", kyid)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, dim='aky')

        !Obtain id for dimension akx
        istatus = nf90_inq_dimid (ncid, "akx", kxid)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, dim='akx')

        !#######################
        !#Obtain dimension info#
        !#######################

        !Get information about dimension
        istatus = nf90_inquire_dimension (ncid, thetaid, len=i)

        !Check information retrieved successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, dimid=thetaid)

        !Perform consistency check
        IF (i /= 2*ntgrid + 1) WRITE(*,*) 'Restart error: ntgrid=? ',i,' : ',ntgrid,' : ',iproc

        !Get information about dimension
        istatus = nf90_inquire_dimension (ncid, signid, len=i)

        !Check information retrieved successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, dimid=signid)

        !Perform consistency check
        IF (i /= 2) WRITE(*,*) 'Restart error: sign=? ',i,' : ',iproc

        !Get information about dimension
        istatus = nf90_inquire_dimension (ncid, gloid, len=i)

        !Check information retrieved successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, dimid=gloid)

        !Perform consistency check
        IF (i /= g_lo%ulim_proc-g_lo%llim_proc+1) WRITE(*,*) 'Restart error: glo=? ',i,' : ',iproc

        !Get information about dimension
        istatus = nf90_inquire_dimension (ncid, kyid, len=i)

        !Check information retrieved successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, dimid=kyid)

        !Perform consistency check
        IF (i /= naky) WRITE(*,*) 'Restart error: naky=? ',i,' : ',naky,' : ',iproc

        !Get information about dimension
        istatus = nf90_inquire_dimension (ncid, kxid, len=i)

        !Check information retrieved successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, dimid=kxid)

        !Perform consistency check
        IF (i /= ntheta0) WRITE(*,*) 'Restart error: ntheta0=? ',i,' : ',ntheta0,' : ',iproc

        !#####################
        !#Obtain variable ids#
        !#####################

        !If finite electrostatic potential read in phi
        IF (fphi > epsilon(0.)) THEN
            !Get variable id
            istatus = nf90_inq_varid (ncid, "phi_r", phir_id)

            !Check id obtained successfully
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='phi_r')

            !Get variable id
            istatus = nf90_inq_varid (ncid, "phi_i", phii_id)

            !Check id obtained successfully
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='phi_i')
        END IF

        !If finite parallel magnetic potential read in apar
        IF (fapar > epsilon(0.)) THEN
            !Get variable id
            istatus = nf90_inq_varid (ncid, "apar_r", aparr_id)

            !Check id obtained successfully
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='apar_r')

            !Get variable id
            istatus = nf90_inq_varid (ncid, "apar_i", apari_id)

            !Check id obtained successfully
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='apar_i')
        END IF

        !If finite perpendicular magnetic potential read in bpar
        IF (fbpar > epsilon(0.)) THEN
            !Get variable id
            istatus = nf90_inq_varid (ncid, "bpar_r", bparr_id)

            !Check id obtained successfully
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='bpar_r')

            !Get variable id
            istatus = nf90_inq_varid (ncid, "bpar_i", bpari_id)

            !Check id obtained successfully
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='bpar_i')
        END IF

        !Get kx_shift id if kx_shift allocated (i.e. if this is a flow shear run)
        IF (allocated(kx_shift)) THEN   ! MR begin
            !Get variable id
            istatus = nf90_inq_varid (ncid, "kx_shift", kx_shift_id)

            !Check id obtained successfully
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='kx_shift')
        END IF   ! MR end

        !Get variable id
        istatus = nf90_inq_varid (ncid, "gr", gr_id)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='gr')

        !Get variable id
        istatus = nf90_inq_varid (ncid, "gi", gi_id)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='gi')
    END IF !End of IF (.not.initialized) THEN

    !###################
    !#Read in variables#
    !###################

    !Allocate temporary arrays to hold dist fn data
    IF (.not. ALLOCATED(tmpr)) &
         ALLOCATE (tmpr(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))
    IF (.not. ALLOCATED(tmpi)) &
         ALLOCATE (tmpi(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))

    !Intitialise temporary arrays to 0
    tmpr = 0.; tmpi = 0.

    !Get real component of distfn
    istatus = nf90_get_var (ncid, gr_id, tmpr)

    !Check read was successful
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, gr_id)

    !Get imaginary component of distfn
    istatus = nf90_get_var (ncid, gi_id, tmpi)

    !Check read was successful
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, gi_id)

    !Recombine components to form complex distribution function
    g = CMPLX(tmpr, tmpi)

    !Read in kx_shift if this is allocated
    IF (allocated(kx_shift)) THEN   ! MR begin
        !Allocate temporary array if needed
        IF (.not. ALLOCATED(stmp)) ALLOCATE (stmp(naky))   ! MR

        !Read in kx_shift
        istatus = nf90_get_var (ncid, kx_shift_id, stmp)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, kx_shift_id)

        !Update kx_shift variable with data from file
        kx_shift = stmp
    END IF   ! MR end

    !Allocate temporary arrays to hold field data
    IF (.not. ALLOCATED(ftmpr)) ALLOCATE (ftmpr(2*ntgrid+1,ntheta0,naky))
    IF (.not. ALLOCATED(ftmpi)) ALLOCATE (ftmpi(2*ntgrid+1,ntheta0,naky))

    !If finite electrostatic potential read in phi
    IF (fphi > epsilon(0.)) THEN
        !Read in real component of phi
        istatus = nf90_get_var (ncid, phir_id, ftmpr)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, phir_id)

        !Read in imaginary component of phi
        istatus = nf90_get_var (ncid, phii_id, ftmpi)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, phii_id)

        !Initialise phi
        phi = 0.

        !Recombine components to create complex phi
        phinew = CMPLX(ftmpr, ftmpi)
    END IF

    !If finite parallel magnetic potential read in apar
    IF (fapar > epsilon(0.)) THEN
        !Read in real component of apar
        istatus = nf90_get_var (ncid, aparr_id, ftmpr)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, aparr_id)

        !Read in imaginary component of apar
        istatus = nf90_get_var (ncid, apari_id, ftmpi)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, apari_id)

        !Initialise apar
        apar = 0.

        !Recombine components to create complex apar
        aparnew = CMPLX(ftmpr, ftmpi)
    END IF

    !If finite perpendicular magnetic potential read in bpar
    IF (fbpar > epsilon(0.)) THEN
        !Read in real component of bpar
        istatus = nf90_get_var (ncid, bparr_id, ftmpr)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, bparr_id)

        !Read in imaginary component of bpar
        istatus = nf90_get_var (ncid, bpari_id, ftmpi)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, bpari_id)

        !Initialise bpar
        bpar = 0.

        !Recombine components to form complex bpar
        bparnew = CMPLX(ftmpr, ftmpi)
    END IF

    !Scale fields if required
    IF (scale > 0.) THEN
        g = g*scale
        phinew = phinew*scale
        aparnew = aparnew*scale
        bparnew = bparnew*scale
    ELSE
        fac = - scale/(MAXVAL(ABS(phinew)))
        g = g*fac
        phinew = phinew*fac
        aparnew = aparnew*fac
        bparnew = bparnew*fac
    END IF

    !######################
    !#Final file functions#
    !######################
    ! RN 2008/05/23: this was commented out. why?
    !Close file
!    istatus = nf90_close (ncid)

    !Check close was successful
    IF (istatus /= NF90_NOERR) THEN
        ierr = error_unit()
        WRITE(ierr,*) "nf90_close error: ", nf90_strerror(istatus),' ',iproc
    END IF

# else
    !Warning message to display if try to call subroutine without netcdf libraries
    WRITE (error_unit(),*) &
        'ERROR: gs2_restore_many is called without netcdf'
# endif
!----------------------------------------------------------------------------------------------

END SUBROUTINE gs2_restore_many
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


! RN 2008/05/23:
!  This can be removed. restore_many seems to work for single proc.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE gs2_restore_one (g, scale, istatus, fphi, fapar, fbpar)
!GS2_RESTORE_ONE(gs2_save.fpp): Subroutine to restore distribution
!function, EM fields and a couple of other variables on one proc.
!Almost identical to GS2_RESTORE_MANY and can probably be replaced
!by this routine.
!==============================================================================================
!MODULE IMPORTS
!==============================================================================================
!MR, 2007: restore kx_shift array if allocated
# ifdef NETCDF
    USE mp, ONLY: proc0, iproc, nproc
    USE fields_arrays, ONLY: phinew, aparnew, bparnew
    USE fields_arrays, ONLY: phi, apar, bpar
    USE dist_fn_arrays, ONLY: kx_shift   ! MR
    USE kt_grids, ONLY: naky, ntheta0
# endif
    USE theta_grid, ONLY: ntgrid
    USE gs2_layouts, ONLY: g_lo
    USE file_utils, ONLY: error_unit
!----------------------------------------------------------------------------------------------

!==============================================================================================
!VARIABLE DECLARATIONS
!==============================================================================================
    !No implicits
    IMPLICIT NONE

    !Interface variables
    COMPLEX, DIMENSION (-ntgrid:,:,g_lo%llim_proc:), INTENT (OUT) :: g
    REAL, INTENT (IN) :: scale
    INTEGER, INTENT (OUT) :: istatus
    REAL, INTENT (IN) :: fphi, fapar, fbpar

# ifdef NETCDF
    !Internal variables
    INTEGER :: n_elements, ierr, i
!----------------------------------------------------------------------------------------------

!==============================================================================================
!WORK SECTION
!==============================================================================================

    !######################
    !#File Initialisations#
    !######################

    !Open specified file for reading
    istatus = nf90_open (TRIM(restart_file), NF90_NOWRITE, ncid)

    !Check open was successful
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, file=restart_file)

    !######################
    !#Obtain dimension ids#
    !######################

    !Obtain id for dimension theta
    istatus = nf90_inq_dimid (ncid, "theta", thetaid)

    !Check id obtained successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, dim='theta')

    !Obtain id for dimension sign
    istatus = nf90_inq_dimid (ncid, "sign", signid)

    !Check id obtained successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, dim='sign')

    !Obtain id for dimension glo
    istatus = nf90_inq_dimid (ncid, "glo", gloid)

    !Check id obtained successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, dim='glo')

    !Obtain id for dimension aky
    istatus = nf90_inq_dimid (ncid, "aky", kyid)

    !Check id obtained successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, dim='aky')

    !Obtain id for dimension akx
    istatus = nf90_inq_dimid (ncid, "akx", kxid)

    !Check id obtained successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, dim='akx')

    !#######################
    !#Obtain dimension info#
    !#######################

    !Get information about dimension
    istatus = nf90_inquire_dimension (ncid, thetaid, len=i)

    !Check information retrieved successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, dimid=thetaid)

    !Perform consistency check
    IF (i /= 2*ntgrid + 1) WRITE(*,*) 'Restart error: ntgrid=? ',i,' : ',ntgrid,' : ',iproc

    !Get information about dimension
    istatus = nf90_inquire_dimension (ncid, signid, len=i)

    !Check information retrieved successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, signid)

    !Perform consistency check
    IF (i /= 2) WRITE(*,*) 'Restart error: sign=? ',i,' : ',iproc

    !Get information about dimension
    istatus = nf90_inquire_dimension (ncid, gloid, len=i)

    !Check information retrieved successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, gloid)

    !Perform consistency check
    IF (i /= g_lo%ulim_world-g_lo%llim_world+1) WRITE(*,*) 'Restart error: glo=? ',i,' : ',iproc

    !Get information about dimension
    istatus = nf90_inquire_dimension (ncid, kyid, len=i)

    !Check information retrieved successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, kyid)

    !Perform consistency check
    IF (i /= naky) WRITE(*,*) 'Restart error: naky=? ',i,' : ',naky,' : ',iproc

    !Get information about dimension
    istatus = nf90_inquire_dimension (ncid, kxid, len=i)

    !Check information retrieved successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, kxid)

    !Perform consistency check
    IF (i /= ntheta0) WRITE(*,*) 'Restart error: ntheta0=? ',i,' : ',ntheta0,' : ',iproc

    !#####################
    !#Obtain variable ids#
    !#####################

    !If finite electrostatic potential read in phi
    IF (fphi > epsilon(0.)) THEN
        !Get variable id
        istatus = nf90_inq_varid (ncid, "phi_r", phir_id)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='phi_r')

        !Get variable id
        istatus = nf90_inq_varid (ncid, "phi_i", phii_id)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='phi_i')
    END IF

    !If finite parallel magnetic potential read in apar
    IF (fapar > epsilon(0.)) THEN
        !Get variable id
        istatus = nf90_inq_varid (ncid, "apar_r", aparr_id)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='apar_r')

        !Get variable id
        istatus = nf90_inq_varid (ncid, "apar_i", apari_id)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='apar_i')
    END IF

    !If finite perpendicular magnetic potential read in bpar
    IF (fbpar > epsilon(0.)) THEN
        !Get variable id
        istatus = nf90_inq_varid (ncid, "bpar_r", bparr_id)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='bpar_r')

        !Get variable id
        istatus = nf90_inq_varid (ncid, "bpar_i", bpari_id)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='bpar_i')
    END IF

    !Get kx_shift id if kx_shift allocated (i.e. if this is a flow shear run)
    IF (allocated(kx_shift)) THEN   ! MR begin
        !Get variable id
        istatus = nf90_inq_varid (ncid, "kx_shift", kx_shift_id)

        !Check id obtained successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='kx_shift')
    END IF   ! MR end

    !Get variable id
    istatus = nf90_inq_varid (ncid, "gr", gr_id)

    !Check id obtained successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='gr')

    !Get variable id
    istatus = nf90_inq_varid (ncid, "gi", gi_id)

    !Check id obtained successfully
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='gi')

    !Define number of elements
    n_elements = g_lo%ulim_proc - g_lo%llim_proc + 1

    !If number of elements not > 0 then processor has nothing to read
    !so exit subroutine
    IF (n_elements <= 0) THEN
        phinew = 0.
        aparnew = 0.
        bparnew = 0.
        RETURN
    END IF

    !###################
    !#Read in variables#
    !###################

    !Allocate temporary arrays to hold dist fn data
    IF (.not. ALLOCATED(tmpr)) ALLOCATE (tmpr(2*ntgrid+1, 2, n_elements))
    IF (.not. ALLOCATED(tmpi)) ALLOCATE (tmpi(2*ntgrid+1, 2, n_elements))

    !Intitialise temporary arrays to 0
    tmpr = 0.;  tmpi = 0.

    !Get real component of distfn
    istatus = nf90_get_var (ncid, gr_id, tmpr, &
        start=(/ 1, 1, g_lo%llim_proc /), &
        count=(/ 2*ntgrid + 1, 2, n_elements /))

    !Check read was successful
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, gr_id)

    !Get imaginary component of distfn
    istatus = nf90_get_var (ncid, gi_id, tmpi, &
        start=(/ 1, 1, g_lo%llim_proc /), &
        count=(/ 2*ntgrid + 1, 2, n_elements /))

    !Check read was successful
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, gi_id)

    !Recombine components to form complex distribution function
    g = CMPLX(tmpr, tmpi)*scale

    !Read in kx_shift if this is allocated
    IF (allocated(kx_shift)) THEN   ! MR begin
        !Allocate temporary array if needed
        IF (.not. ALLOCATED(stmp)) ALLOCATE (stmp(naky))   ! MR

        !Read in kx_shift
        istatus = nf90_get_var (ncid, kx_shift_id, stmp)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, kx_shift_id)

        !Update kx_shift variable with data from file
        kx_shift = stmp
    END IF   ! MR end

    !Allocate temporary arrays to hold field data
    IF (.not. ALLOCATED(ftmpr)) ALLOCATE (ftmpr(2*ntgrid+1,ntheta0,naky))
    IF (.not. ALLOCATED(ftmpi)) ALLOCATE (ftmpi(2*ntgrid+1,ntheta0,naky))

    !If finite electrostatic potential read in phi
    IF (fphi > epsilon(0.)) THEN
        !Read in real component of phi
        istatus = nf90_get_var (ncid, phir_id, ftmpr)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, phir_id)

        !Read in imaginary component of phi
        istatus = nf90_get_var (ncid, phii_id, ftmpi)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, phii_id)

        !Initialise phi
        phi = 0.

        !Recombine components to create complex phi
        phinew = CMPLX(ftmpr, ftmpi)*scale
    END IF

    !If finite parallel magnetic potential read in apar
    IF (fapar > epsilon(0.)) THEN
        !Read in real component of bpar
        istatus = nf90_get_var (ncid, aparr_id, ftmpr)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, aparr_id)

        !Read in imaginary component of apar
        istatus = nf90_get_var (ncid, apari_id, ftmpi)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, apari_id)

        !Initialise apar
        apar = 0.

        !Recombine components to create complex apar
        aparnew = CMPLX(ftmpr, ftmpi)*scale
    END IF

    !If finite perpendicular magnetic potential read in bpar
    IF (fbpar > epsilon(0.)) THEN
        !Read in real component of bpar
        istatus = nf90_get_var (ncid, bparr_id, ftmpr)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, bparr_id)

        !Read in imaginary component of bpar
        istatus = nf90_get_var (ncid, bpari_id, ftmpi)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, bpari_id)

        !Initialise bpar
        bpar = 0.

        !Recombine components to form complex bpar
        bparnew = CMPLX(ftmpr, ftmpi)*scale
    END IF

    !######################
    !#Final file functions#
    !######################

    ! RN 2008/05/23: this was commented out. why?
    !Close file
!    istatus = nf90_close (ncid)

    !Check close was successful
    IF (istatus /= NF90_NOERR) THEN
        ierr = error_unit()
        WRITE(ierr,*) "nf90_close error: ", nf90_strerror(istatus),' ',iproc
    END IF

# else
    !Warning message to display if try to call subroutine without netcdf libraries
    WRITE (error_unit(),*) &
        'ERROR: gs2_restore_one is called without netcdf'
# endif
!----------------------------------------------------------------------------------------------

END SUBROUTINE gs2_restore_one
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE init_save (file)
!INIT_SAVE(gs2_save): Subroutine to set the restart_file variable used
!to define the output files prefix

!==============================================================================================
!VARIABLE DECLARATIONS
!==============================================================================================
    !Interface variables
    CHARACTER(300), INTENT (IN) :: file
!----------------------------------------------------------------------------------------------

!==============================================================================================
!WORK SECTION
!==============================================================================================

    !Update internal variable with passed arguement
    restart_file = file
!----------------------------------------------------------------------------------------------

END SUBROUTINE init_save
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE init_dt (delt0, istatus)
!INIT_DT(gs2_save.fpp): Subroutine to read in the time step
!from a restart file. Runs on proc0 only and then broadcasts
!result to all other procs
!==============================================================================================
!MODULE IMPORTS
!==============================================================================================
# ifdef NETCDF
    USE mp, ONLY: nproc, proc0, broadcast
    USE file_utils, ONLY: error_unit
# endif
!----------------------------------------------------------------------------------------------

!==============================================================================================
!VARIABLE DECLARATIONS
!==============================================================================================
    !No implicits
    IMPLICIT NONE

    !Interface variables
    REAL, INTENT (IN OUT) :: delt0
    INTEGER, INTENT (OUT) :: istatus

# ifdef NETCDF
    !Internal variables
    CHARACTER (306) :: file_proc
!----------------------------------------------------------------------------------------------

!==============================================================================================
!WORK SECTION
!==============================================================================================
    !If on master processor perform below else don't do anything
    IF (proc0) THEN
        !If file not initialised
        IF (.not. initialized) THEN

            !Define restart file name (for proc0)
            file_proc=TRIM(TRIM(restart_file)//'.0')

            !Open restart file
            istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)

            !Check open was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus,file=file_proc)

            !Get id for variable 'delt0'
            istatus = nf90_inq_varid (ncid, "delt0", delt0id)

            !Check id retrieved successfully
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='delt0')
        END IF

        !Read in variable delt0
        istatus = nf90_get_var (ncid, delt0id, delt0)

        !Check delt0 was successfully read
        IF (istatus /= NF90_NOERR) THEN
            CALL netcdf_error (istatus, ncid, delt0id, message=' in init_dt')
            delt0 = -1.
        END IF

        !If file not initialised (i.e. if first section was performed) close file (assume there is no error so far)
        IF (.not.initialized) istatus = nf90_close (ncid)
    END IF

    !Send error status to all procs
    CALL broadcast (istatus)

    !Send delt0 to all procs
    CALL broadcast (delt0)

# endif
!----------------------------------------------------------------------------------------------

END SUBROUTINE init_dt
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE init_vnm (vnm, istatus)
!INIT_VNM(gs2_save.fpp): Subroutine to read in vnm[1|2]
!from a restart file. Runs on proc0 only and then broadcasts
!result to all other procs
!==============================================================================================
!MODULE IMPORTS
!==============================================================================================
# ifdef NETCDF
    USE mp, ONLY: nproc, proc0, broadcast
    USE file_utils, ONLY: error_unit
# endif
!----------------------------------------------------------------------------------------------

!==============================================================================================
!VARIABLE DECLARATIONS
!==============================================================================================
    !No implicits
    IMPLICIT NONE

    !Interface variables
    REAL, DIMENSION(2), INTENT (IN OUT) :: vnm
    INTEGER, INTENT (OUT) :: istatus

# ifdef NETCDF
    !Internal variables
    CHARACTER (306) :: file_proc
!----------------------------------------------------------------------------------------------

!==============================================================================================
!WORK SECTION
!==============================================================================================
    !If on master processor perform below else don't do anything
    IF (proc0) THEN
        !If file not initialised
        IF (.not.initialized) THEN
            !Define restart file name (for proc0)
            file_proc=TRIM(TRIM(restart_file)//'.0')

            !Open restart file
            istatus = nf90_open (file_proc, 0, ncid)

            !Check open was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, file=file_proc)

            !Get id for variable 'vnm1'
            istatus = nf90_inq_varid (ncid, "vnm1", vnm1id)

            !Check id retrieved successfully
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='vnm1')

            !Get id for variable 'vnm2'
            istatus = nf90_inq_varid (ncid, "vnm2", vnm2id)

            !Check id retrieved successfully
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='vnm2')
        END IF

        !Read in variable vnm1
        istatus = nf90_get_var (ncid, vnm1id, vnm(1))

        !Check vnm1 was successfully read
        IF (istatus /= NF90_NOERR) THEN
            CALL netcdf_error (istatus, ncid, vnm1id, message=' in init_vnm')
            vnm(1) = 0.
        END IF

        !Read in variable vnm2
        istatus = nf90_get_var (ncid, vnm2id, vnm(2))

        !Check vnm2 was successfully read
        IF (istatus /= NF90_NOERR) THEN
            CALL netcdf_error (istatus, ncid, vnm2id, message=' in init_vnm')
            vnm(2) = 0.
        END IF

        !If file not initialised (i.e. if first section was performed) close file (assume there is no error so far)
        IF (.not. initialized) istatus = nf90_close (ncid)
    END IF

    !Send error status to all procs
    CALL broadcast (istatus)

    !Send delt0 to all procs
    CALL broadcast (vnm)

# endif
!----------------------------------------------------------------------------------------------

END SUBROUTINE init_vnm
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! This routine gets a_ant and b_ant for proc 0 only!!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE init_ant_amp (a_ant, b_ant, nk_stir, istatus)
!INIT_ANT_AMP(gs2_save.fpp): Subroutine to read in the antenna
!setting from a restart file. Runs on proc0 only and then
!broadcasts result to all other procs
!==============================================================================================
!MODULE IMPORTS
!==============================================================================================
# ifdef NETCDF
    USE file_utils, ONLY: error_unit
    USE constants, ONLY: zi
# endif
    USE mp, ONLY: proc0
!----------------------------------------------------------------------------------------------

!==============================================================================================
!VARIABLE DECLARATIONS
!==============================================================================================
    !No implicits
    IMPLICIT NONE

    !Interface variables
    COMPLEX, DIMENSION(:), INTENT (IN OUT) :: a_ant, b_ant
    INTEGER, INTENT (IN) :: nk_stir
    INTEGER, INTENT (OUT) :: istatus

# ifdef NETCDF
    !Internal variables
    CHARACTER (306) :: file_proc
    INTEGER :: ierr, i
!----------------------------------------------------------------------------------------------

!==============================================================================================
!WORK SECTION
!==============================================================================================
    !If on master processor perform below else don't do anything
    IF (proc0) THEN
        !Initialise a_and and b_ant
        a_ant = 0. ; b_ant = 0.

        !Define restart file name (for proc0)
        file_proc=TRIM(TRIM(restart_file)//'.0')

        !Open restart file
        istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)

        !Check open was successful
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_open in init_ant_amp error: ", nf90_strerror(istatus)
            WRITE(ierr,*) "If you did not intend for this to be a restarted run with an external antenna,"
            WRITE(ierr,*) "you may ignore the error message above."
            RETURN
        END IF

        !Get id for dimension 'nk_stir'
        istatus = nf90_inq_dimid (ncid, "nk_stir", nk_stir_dim)

        !Check id retrieved successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, dim='nk_stir')

        !Get information about dimension
        istatus = nf90_inquire_dimension (ncid, nk_stir_dim, len=i)

        !Check information retrieved successfully
        IF (istatus /= NF90_NOERR) &
            CALL netcdf_error (istatus, ncid, dimid=nk_stir_dim)

        !Perform a consistency check
        IF (i /= nk_stir) WRITE(*,*) 'Restart error: nk_stir=? ',i,' : ',nk_stir

        !Get id for variable 'a_ant_r'
        istatus = nf90_inq_varid (ncid, "a_ant_r", a_antr_id)

        !Check id retrieved successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='a_ant_r')

        !Get id for variable 'a_ant_i'
        istatus = nf90_inq_varid (ncid, "a_ant_i", a_anti_id)

        !Check id retrieved successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='a_ant_i')

        !Get id for variable 'b_ant_r'
        istatus = nf90_inq_varid (ncid, "b_ant_r", b_antr_id)

        !Check id retrieved successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='b_ant_r')

        !Get id for variable 'b_ant_i'
        istatus = nf90_inq_varid (ncid, "b_ant_i", b_anti_id)

        !Check id retrieved successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='b_ant_i')

        !If temporary array not allocated allocate it
        IF (.not. ALLOCATED(atmp)) ALLOCATE (atmp(nk_stir))

        !Initialise temporary array
        atmp = 0.

        !Read variable into atmp
        istatus = nf90_get_var (ncid, a_antr_id, atmp)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, a_antr_id)

        !Put real component into array
        a_ant = atmp

        !Read variable into atmp
        istatus = nf90_get_var (ncid, a_anti_id, atmp)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, a_anti_id)

        !Update array with complex component
        a_ant = a_ant + zi * atmp

        !Read variable into atmp
        istatus = nf90_get_var (ncid, b_antr_id, atmp)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, b_antr_id)

        !Put real component into array
        b_ant = atmp

        !Read variable into atmp
        istatus = nf90_get_var (ncid, b_anti_id, atmp)

        !Check read was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, b_anti_id)

        !Update array with complex component
        b_ant = b_ant + zi * atmp

        !Deallocate temp array
        DEALLOCATE (atmp)

        !Close restart file (assume there are no errors)
        istatus = nf90_close (ncid)

    END IF

# else
    !If no netcdf set error status to 2 (on proc0 only)
    IF (proc0) istatus = 2

# endif
!----------------------------------------------------------------------------------------------

END SUBROUTINE init_ant_amp
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE init_tstart (tstart, istatus)
!INIT_TSTART(gs2_save.fpp): Subroutine to read in t0, the current
!simulation time, from a restart file. Runs on proc0 only and
!then broadcasts result to all other procs
!==============================================================================================
!MODULE IMPORTS
!==============================================================================================
# ifdef NETCDF
    USE mp, ONLY: nproc, proc0, broadcast
    USE file_utils, ONLY: error_unit
# endif
!----------------------------------------------------------------------------------------------

!==============================================================================================
!VARIABLE DEFINITIONS
!==============================================================================================
    !No implicits
    IMPLICIT NONE

    !Interface variables
    REAL, INTENT (IN OUT) :: tstart
    INTEGER, INTENT (OUT) :: istatus

# ifdef NETCDF
    !Internal variables
    CHARACTER (306) :: file_proc
    INTEGER :: ierr
!----------------------------------------------------------------------------------------------

!==============================================================================================
!WORK SECTION
!==============================================================================================
    !If on master processor perform below else don't do anything
    IF (proc0) THEN
        !If file not initialised
        IF (.not.initialized) THEN
            !Define restart file name (for proc0)
            file_proc=TRIM(TRIM(restart_file)//'.0')

            !Open restart file
            istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)

            !Check open was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, file=file_proc)
        END IF

        !Get id for variable 't0'
        istatus = nf90_inq_varid (ncid, "t0", t0id)

        !Check id retrieved successfully
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, var='t0')

        !Read in variable t0
        istatus = nf90_get_var (ncid, t0id, tstart)

        !Check variable successfully read
        IF (istatus /= NF90_NOERR) THEN
            CALL netcdf_error (istatus, ncid, t0id, message=' in init_tstart')
            tstart = -1.
        END IF

        !If file not initialised (i.e. if first section was performed) close file (assume there is no error so far)
        IF (.not.initialized) istatus = nf90_close (ncid)
    END IF

    !Broadcast istatus to all procs
    CALL broadcast (istatus)

    !Broadcast tstart to all procs
    CALL broadcast (tstart)

# endif
!----------------------------------------------------------------------------------------------

END SUBROUTINE init_tstart
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!!$
!!$# ifdef NETCDF
!!$  function get_netcdf_code_precision () result (code_real)
!!$
!!$    use file_utils, only: error_unit
!!$    use constants, only: pi, kind_rs, kind_rd
!!$    integer :: code_real
!!$    integer :: ierr
!!$
!!$    ! second condition for Cray
!!$    if ( (kind(pi)==kind_rs) .or. (kind_rs==kind_rd) ) then
!!$       code_real = NF90_FLOAT
!!$    else if (kind(pi)==kind_rd) then
!!$       code_real = NF90_DOUBLE
!!$    else
!!$       ierr = error_unit()
!!$       write (ierr,*) 'ERROR: precision mismatch in get_netcdf_code_precision'
!!$    end if
!!$
!!$  end function get_netcdf_code_precision
!!$
!!$  subroutine check_netcdf_file_precision (ncid, filename)
!!$
!!$    use file_utils, only: error_unit
!!$
!!$    integer (kind_nf), intent (in), optional :: ncid
!!$    character (*), intent (in), optional :: filename
!!$    integer (kind_nf) :: file_real
!!$    integer (kind_nf) :: ist, ncid_private, tid
!!$    integer :: ierr
!!$
!!$    ist = NF90_NOERR
!!$    file_real = -1
!!$
!!$    if (present(ncid)) then
!!$       if (present(filename)) then
!!$          ierr = error_unit()
!!$          write (ierr,*) 'WARNING: in calling check_netcdf_file_precision'
!!$          write (ierr,*) &
!!$               'WARNING: both filename and ncid given -- filename ignored'
!!$       end if
!!$       ncid_private = ncid
!!$    else
!!$       if (present(filename)) then
!!$          ist = nf90_open (filename, NF90_NOWRITE, ncid_private)
!!$          if (test) write (error_unit(),*) &
!!$               'opened netcdf file ', trim(filename), ' with ncid: ', &
!!$               ncid_private, ' in check_netcdf_file_precision'
!!$          if (ist /= NF90_NOERR) then
!!$             call netcdf_error (ist, file=filename)
!!$             return
!!$          end if
!!$       else
!!$          ierr = error_unit()
!!$          write (ierr,*) 'ERROR: in calling check_netcdf_file_precision'
!!$          write (ierr,*) 'ERROR: either filename or ncid should be given'
!!$          return
!!$       end if
!!$    end if
!!$
!!$    ist = nf90_inq_varid (ncid_private, 't0', tid)
!!$    if (ist /= NF90_NOERR) call netcdf_error (ist, var='t0')
!!$
!!$    ! get file_real
!!$    if (ist == NF90_NOERR) then
!!$       ist = nf90_inquire_variable (ncid_private, tid, xtype=file_real)
!!$       if (ist /= NF90_NOERR) call netcdf_error (ist, ncid_private, tid)
!!$    end if
!!$
!!$    if (.not.present(ncid)) then
!!$       ist = nf90_close (ncid_private)
!!$       if (ist /= NF90_NOERR) call netcdf_error (ist, file=filename)
!!$    end if
!!$
!!$    ! check if file_real == code_real
!!$    if (file_real /= netcdf_real) then
!!$       ierr = error_unit()
!!$       write (ierr,*) 'WARNING: precision mismatch in input netcdf file and running code'
!!$       if (file_real == NF90_FLOAT) then
!!$          write (ierr,*) 'WARNING: file_real = NF90_FLOAT'
!!$       else if (file_real == NF90_DOUBLE) then
!!$          write (ierr,*) 'WARNING: file_real = NF90_DOUBLE'
!!$       else
!!$          write (ierr,*) 'WARNING: unknown file_real', file_real
!!$       end if
!!$       if (netcdf_real == NF90_FLOAT) then
!!$          write (ierr,*) 'WARNING: code_real = NF90_FLOAT'
!!$       else if (netcdf_real == NF90_DOUBLE) then
!!$          write (ierr,*) 'WARNING: code_real = NF90_DOUBLE'
!!$       else
!!$          write (ierr,*) 'WARNING: unknown code_real'
!!$       end if
!!$    end if
!!$
!!$  end subroutine check_netcdf_file_precision
!!$
!!$  subroutine netcdf_error &
!!$       (istatus, ncid, varid, dimid, file, dim, var, att, message)
!!$
!!$    use mp, only: proc0, iproc
!!$    use file_utils, only: error_unit
!!$    use netcdf, only: NF90_GLOBAL
!!$
!!$    integer (kind_nf), intent (in) :: istatus
!!$    integer (kind_nf), intent (in), optional :: ncid
!!$    integer (kind_nf), intent (in), optional :: varid
!!$    integer (kind_nf), intent (in), optional :: dimid
!!$    character (*), intent (in), optional :: file
!!$    character (*), intent (in), optional :: dim
!!$    character (*), intent (in), optional :: var
!!$    character (*), intent (in), optional :: att
!!$    character (*), intent (in), optional :: message
!!$    integer (kind_nf) :: ist
!!$    integer :: ierr
!!$    character (20) :: varname, dimname
!!$
!!$    ierr = error_unit()
!!$
!!$    write (ierr, '(2a,$)') 'ERROR: ', trim (nf90_strerror (istatus))
!!$    ! TT: If $ control fails, there is an alternative advance='no' specifier
!!$
!!$    if (present(file)) &
!!$         write (ierr, '(2a,$)') ' in file: ', trim (file)
!!$
!!$    if (present(dim)) &
!!$         write (ierr, '(2a,$)') ' in dimension: ', trim (dim)
!!$
!!$    if (present(var)) &
!!$         write (ierr, '(2a,$)') ' in variable: ', trim (var)
!!$
!!$    if (present(varid)) then
!!$       if (present(ncid)) then
!!$          if ( (varid == NF90_GLOBAL) .and. present(att) ) then
!!$             write (ierr, '(2a)') ' in global attribute: ', trim(att)
!!$             return
!!$          else
!!$             ist = nf90_inquire_variable (ncid, varid, varname)
!!$             if (ist == NF90_NOERR) then
!!$                write (ierr, '(a,i8,2a,$)') ' in varid: ', varid, &
!!$                     & ' variable name: ', trim (varname)
!!$             else
!!$                write (ierr, *) ''
!!$                write (ierr, '(3a,i8,a,i8,$)') 'ERROR in netcdf_error: ', &
!!$                     trim (nf90_strerror(ist)), ' in varid: ', varid, &
!!$                     ', ncid: ', ncid
!!$             end if
!!$          end if
!!$          if (present(att)) &
!!$               write (ierr, '(2a)') ' with the attribute: ', trim(att)
!!$       else
!!$          write (ierr, *) ''
!!$          write (ierr, '(2a,$)') 'ERROR in netcdf_error: ', &
!!$               & 'ncid missing while varid present in the argument'
!!$       end if
!!$    end if
!!$
!!$    if (present(dimid)) then
!!$       if (present(ncid)) then
!!$          ist = nf90_inquire_dimension (ncid, dimid, dimname)
!!$          if (ist == NF90_NOERR) then
!!$             write (ierr, '(a,i8,2a,$)') ' in dimid: ', dimid, &
!!$                  & ' dimension name: ', trim (dimname)
!!$          else
!!$             write (ierr, *) ''
!!$             write (ierr, '(3a,i8,a,i8,$)') 'ERROR in netcdf_error: ', &
!!$                  trim (nf90_strerror(ist)), ' in dimid: ', dimid, &
!!$                  ', ncid: ', ncid
!!$          end if
!!$       else
!!$          write (ierr, *) ''
!!$          write (ierr, '(2a,$)') 'ERROR in netcdf_error: ', &
!!$               & 'ncid missing while dimid present in the argument'
!!$       end if
!!$    end if
!!$
!!$    if (present(message)) write (ierr, '(a,$)') trim(message)
!!$
!!$    write (ierr, '(a,i8)') ' on iproc: ', iproc
!!$
!!$  end subroutine netcdf_error
!!$# endif
!!$! <TT
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE gs2_save_distfn &
     (g, t0, delt0, vnm, istatus, fphi, fapar, fbpar, exit_in)
!<DD 27-08-2010>
!GS2_SAVE_DISTFN: Saves a copy of the distribution function adjusted to convert
!from internal GS2 distribution function (h) to actual component of interest (g).
!Produces a file with name prefix".dfn."proc_no.
!Based on copy of gs2_save_for_restart hence probably saves more than desired
!</DD>

!==============================================================================================
!MODULE IMPORTS
!==============================================================================================
    USE constants, ONLY: kind_rs, kind_rd, pi           !Standard constants
    USE fields_arrays, ONLY: phinew, aparnew, bparnew   !Fields
    USE dist_fn_arrays, ONLY: kx_shift                  !Important for runs with flow shear
    USE kt_grids, ONLY: naky, ntheta0                   !Number of ky and kx used
    USE mp, ONLY: nproc, iproc, proc0                   !MPI related constants
    USE theta_grid, ONLY: ntgrid                        !Number of theta grid points

    ! Must include g_layout_type here to avoid obscure bomb while compiling
    ! gs2_diagnostics.f90 (which uses this module) with the Compaq F90 compiler:
    ! TT>
    USE gs2_layouts, ONLY: g_lo                         !Distfn index
    USE layouts_type, ONLY: g_layout_type               !Layout
    ! <TT

    USE file_utils, ONLY: error_unit                    !Error file unit

    !<DD 29-08-2010> Added use of energy and lambda grids so they can be written
    USE le_grids, ONLY: e,al                            !Energy and lambda grids
    USE le_grids, ONLY: negrid, nlambda                 !Number of energy and lambda points
    !</DD>

    !<DD 01-09-2010> Added number of species
    USE species, ONLY: nspec                            !Number of species
    !</DD>

    !<DD 02-09-2010> Added use of parallel and perpendicular velocities
    USE dist_fn_arrays, ONLY: vpa, vperp2               !Parallel and perpendicular velocities
    !</DD>


!----------------------------------------------------------------------------------------------

!==============================================================================================
!VARIABLE DECLARATIONS
!==============================================================================================
    !No implicits
    IMPLICIT NONE

    !Interface variables
    REAL, INTENT (IN) :: t0, delt0                  !User time and time step
    REAL, DIMENSION (2), INTENT (IN) :: vnm         !
    REAL, INTENT (IN) :: fphi, fapar, fbpar         !Fields
    INTEGER, INTENT (OUT) :: istatus                !Write status variable
    LOGICAL, INTENT (IN), OPTIONAL :: exit_in       !True=>Close file at end, False=>Sync file at end
    COMPLEX, DIMENSION (-NTGRID:,:,G_LO%LLIM_PROC:), INTENT (IN) :: g !Adjusted Distribution function

# ifdef NETCDF
    !Internal variables
    CHARACTER (306) :: file_proc
    CHARACTER (10) :: suffix
    INTEGER :: i, n_elements, ierr
    LOGICAL :: exit
!----------------------------------------------------------------------------------------------

!==============================================================================================
!WORK SECTION
!==============================================================================================
    !Decide what to do with file at end of routine
    IF (PRESENT(exit_in)) THEN
       exit = exit_in
    ELSE
       exit = .false.
    END IF

    !Get the number of elements for this processor
    n_elements = g_lo%ulim_proc-g_lo%llim_proc+1

    !If no elements then nothing to write so escape routine
    IF (n_elements <= 0) RETURN

    !<DD 27-08-2010> Translate h back to g
    !CALL g_adjust (g, phinew, bparnew, fphi, fbpar)
    !</DD>

    !#####################
    !#File initialisation#
    !#####################
    !If file not initialised then create file, define dimensions and variables
    IF (.not. initialized_dfn) THEN
        !Don't perform initialisation twice
        initialized_dfn = .true.

        !Define file prefix
        file_proc = TRIM(restart_file)

        !Define file suffix (.dfn.procnum)
        WRITE (suffix,'(a5,i0)') '.dfn.', iproc

        !Define full file name
        file_proc = TRIM(TRIM(file_proc)//ADJUSTL(suffix))

        !Create file
        istatus = nf90_create (file_proc, NF90_CLOBBER, ncid)

        !Check file was created successfully
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_create error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !###################
        !#Define dimensions#
        !###################
        !<note> this if statement shouldn't be needed due to earlier check which
        !results in subroutine exit if n_elements <=0 </note>
        IF (n_elements > 0) THEN
            !Define theta dimension
            istatus = nf90_def_dim (ncid, "theta", 2*ntgrid+1, thetaid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim theta error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define sign dimension (sign of v_par)
            istatus = nf90_def_dim (ncid, "sign", 2, signid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim sign error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define glo dimension (corresponds to layout grid i.e. kx,ky,lambda,energy,species mapped to 1d array)
            istatus = nf90_def_dim (ncid, "glo", n_elements, gloid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim glo error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define aky dimension (number of ky modes)
            istatus = nf90_def_dim (ncid, "aky", naky, kyid)
            IF (istatus /= NF90_NOERR) THEN
                 ierr = error_unit()
                 WRITE(ierr,*) "nf90_def_dim aky error: ", nf90_strerror(istatus)
                 GOTO 1
            END IF

            !Define akx dimension (number of kx modes)
            istatus = nf90_def_dim (ncid, "akx", ntheta0, kxid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                 ierr = error_unit()
                 WRITE(ierr,*) "nf90_def_dim akx error: ", nf90_strerror(istatus)
                 GOTO 1
            END IF

            !<DD 29-08-2010> Define energy and lambda dimensions

            !Define negrid dimension (number of energy grid points)
            istatus = nf90_def_dim (ncid, "negrid", negrid, egridid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim negrid error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define nlambda dimension (number of pitch angles)
            istatus = nf90_def_dim (ncid, "nlambda", nlambda, lgridid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim nlambda error: ", nf90_strerror(istatus)
                GOTO 1
            END IF
            !</DD>

            !<DD 01-09-2010> Define species dimension

            !Define nspec dimension (number of species)
            istatus = nf90_def_dim (ncid, "nspec", nspec, nspecid)

            !Check dimension created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_dim nspec error: ", nf90_strerror(istatus)
                GOTO 1
            END IF
            !</DD>

        END IF !End of IF (n_elements > 0) THEN


        !###################
        !#Define variables #
        !###################

        !Get real variable data type specifier for netcdf routines
        if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()

        !Define t0 variable (user time)
        istatus = nf90_def_var (ncid, "t0", netcdf_real, t0id)

        !Check variable created successfully
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_def_var t0 error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !Define delt0 variable (time step)
        istatus = nf90_def_var (ncid, "delt0", netcdf_real, delt0id)

        !Check variable created successfully
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_def_var delt0 error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !Define vnm1 variable (??)
        istatus = nf90_def_var (ncid, "vnm1", netcdf_real, vnm1id)

        !Check variable created successfully
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_def_var vnm(1) error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !Define vnm2 variable (??)
        istatus = nf90_def_var (ncid, "vnm2", netcdf_real, vnm2id)

        !Check variable created successfully
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_def_var vnm(2) error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !Define kx_shift variable (gives shift in kx due to flow shear
        istatus = nf90_def_var (ncid, "kx_shift", netcdf_real, &
            (/ kyid /), kx_shift_id)

        !Check variable created successfully
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE(ierr,*) "nf90_def_var kx_shift error: ", nf90_strerror(istatus)
            GOTO 1
        END IF

        !<note> this if statement shouldn't be needed due to earlier check which
        !results in subroutine exit if n_elements <=0 </note>
        IF (n_elements > 0) THEN

            !Define gr variable (real component of g)
            istatus = nf90_def_var (ncid, "gr", netcdf_real, &
                (/ thetaid, signid, gloid /), gr_id)

            !Check variable created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var g error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define gi variable (imaginary component of gi)
            istatus = nf90_def_var (ncid, "gi", netcdf_real, &
                (/ thetaid, signid, gloid /), gi_id)

            !Check variable created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var g error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !If there are finite electrostatic potential perturbations
            !then write them to file
            IF (fphi > epsilon(0.)) THEN
                !Define phi_r variable (real component of \phi)
                istatus = nf90_def_var (ncid, "phi_r", netcdf_real, &
                    (/ thetaid, kxid, kyid /), phir_id)

                !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var phi error: ", nf90_strerror(istatus)
                    GOTO 1
                END IF

                !Define phi_i variable (imaginary component of \phi)
                istatus = nf90_def_var (ncid, "phi_i", netcdf_real, &
                    (/ thetaid, kxid, kyid /), phii_id)

                !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var phi error: ", nf90_strerror(istatus)
                    GOTO 1
                END IF
            END IF

            !If there are finite parallel magnetic potential perturbations
            !then write them to file
            IF (fapar > epsilon(0.)) THEN
                !Define apar_r variable (real component A_\parallel)
                istatus = nf90_def_var (ncid, "apar_r", netcdf_real, &
                    (/ thetaid, kxid, kyid /), aparr_id)

                !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var apar error: ", nf90_strerror(istatus)
                    GOTO 1
                END if

                !Define apar_i variable (imaginary component A_\parallel)
                istatus = nf90_def_var (ncid, "apar_i", netcdf_real, &
                    (/ thetaid, kxid, kyid /), apari_id)

                !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var apar error: ", nf90_strerror(istatus)
                    GOTO 1
                END if
            END IF

            !If there are finite perpendicular magnetic potential perturbations
            !then write them to file
            IF (fbpar > epsilon(0.)) THEN
                !Define bpar_r variable (real component B_\parallel)
                istatus = nf90_def_var (ncid, "bpar_r", netcdf_real, &
                    (/ thetaid, kxid, kyid /), bparr_id)

                !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var bparr error: ", nf90_strerror(istatus)
                    GOTO 1
                END IF

                !Define bpar_i variable (imaginary component B_\parallel)
                istatus = nf90_def_var (ncid, "bpar_i", netcdf_real, &
                    (/ thetaid, kxid, kyid /), bpari_id)

                !Check variable created successfully
                IF (istatus /= NF90_NOERR) THEN
                    ierr = error_unit()
                    WRITE(ierr,*) "nf90_def_var bpari error: ", nf90_strerror(istatus)
                    GOTO 1
                END IF
            END IF

            !<DD 29-08-2010> Define energy and lambda variables
            !Define energy variable (energy for each species)
            istatus = nf90_def_var (ncid, "energy", netcdf_real, &
                (/ egridid, nspecid /), energy_id)

            !Check variable created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var energy error: ", nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define lambda variable (pitch angles)
            istatus = nf90_def_var (ncid, "lambda", netcdf_real, &
                (/ lgridid /), lambda_id)

            !Check variable created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var lambda error: ", nf90_strerror(istatus)
                GOTO 1
            END IF
            !</DD>

            !<DD 01-09-2010> define species variable
            !Define species variable
            istatus = nf90_def_var (ncid, "species", netcdf_real, &
                (/ nspecid /), spec_id)
            !Check variable created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var species error: ",nf90_strerror(istatus)
                GOTO 1
            END IF
            !</DD>

            !<DD 02-09-2010> Define velocity variables
            !Define vpa variable
            istatus = nf90_def_var (ncid, "vpa", netcdf_real, &
                (/ thetaid, signid, gloid /), vpa_id)

            !Check variable created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var vpa error: ",nf90_strerror(istatus)
                GOTO 1
            END IF

            !Define vperp2 variable
            istatus = nf90_def_var (ncid, "vperp2", netcdf_real, &
                (/ thetaid, gloid /), vperp2_id)

            !Check variable created successfully
            IF (istatus /= NF90_NOERR) THEN
                ierr = error_unit()
                WRITE(ierr,*) "nf90_def_var vperp2 error: ",nf90_strerror(istatus)
                GOTO 1
            END IF
            !</DD>
        END IF !End of IF (n_elements > 0) THEN

        !######################
        !#End definition stage#
        !######################

        !End definition stage of file
        istatus = nf90_enddef (ncid)

        !Check end was successful
        IF (istatus /= NF90_NOERR) THEN
            ierr = error_unit()
            WRITE (ierr,*) "nf90_enddef error: ", nf90_strerror(istatus)
            GOTO 1
        END IF
    END IF !End of IF (.not. initialized) THEN


    !#################
    !#File population#
    !#################

    !#################
    !#Store variables#
    !#################

    !Store delt0
    istatus = nf90_put_var (ncid, delt0id, delt0)

    !Check store was successful
    IF (istatus /= NF90_NOERR) THEN
        ierr = error_unit()
        WRITE (ierr,*) "nf90_put_var delt0 error: ", nf90_strerror(istatus)
        GOTO 1
    END IF

    !Store t0
    istatus = nf90_put_var (ncid, t0id, t0)

    !Check store was successful
    IF (istatus /= NF90_NOERR) THEN
        ierr = error_unit()
        WRITE (ierr,*) "nf90_put_var t0 error: ", nf90_strerror(istatus)
        GOTO 1
    END IF

    !Store vnm(1)
    istatus = nf90_put_var (ncid, vnm1id, vnm(1))

    !Check store was successful
    IF (istatus /= NF90_NOERR) THEN
        ierr = error_unit()
        WRITE (ierr,*) "nf90_put_var vnm(1) error: ", nf90_strerror(istatus)
        GOTO 1
    END IF

    !Store vnm(2)
    istatus = nf90_put_var (ncid, vnm2id, vnm(2))

    !Check store was successful
    IF (istatus /= NF90_NOERR) THEN
        ierr = error_unit()
        WRITE (ierr,*) "nf90_put_var vnm(2) error: ", nf90_strerror(istatus)
        GOTO 1
    END IF


    !Make temporary array to hold kx_shift if allocated, zeros if not
    IF (.not. ALLOCATED(stmp)) ALLOCATE (stmp(naky))

    !If kx_shift is allocated then output that else output zeros
    IF (ALLOCATED(kx_shift)) THEN
       stmp = kx_shift
    ELSE
       stmp = 0.
    END IF

    !Store variable (kx_shift/zeros)
    istatus = nf90_put_var (ncid, kx_shift_id, stmp)

    !Check store was successful
    IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, kx_shift_id)


    !<note> this if statement shouldn't be needed due to earlier check which
    !results in subroutine exit if n_elements <=0 </note>
    IF (n_elements > 0) THEN
        !Allocate array to hold components of g individually
        IF (.not. ALLOCATED(tmpr)) &
            ALLOCATE (tmpr(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))

        !Make temporary copy of the real component of g
        tmpr = REAL(g)

        !Store variable real(g)
        istatus = nf90_put_var (ncid, gr_id, tmpr)

        !Check store was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, gr_id)

        !Make temporary copy of the imaginary component of g
        tmpr = AIMAG(g)

        !Store variable aimag(g)
        istatus = nf90_put_var (ncid, gi_id, tmpr)

        !Check store was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, gi_id)

        !<DD 29-08-2010> Fill energy and lambda information

        !Store variable energy
        istatus = nf90_put_var (ncid, energy_id, e(:,1:nspec))

        !Check store was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, energy_id)

        !Store variable lambda
        istatus = nf90_put_var (ncid, lambda_id, al)

        !Check store was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, lambda_id)
        !</DD>

        !<DD 02-09-2010> Fill velocity variables

        !Store variable vpa
        istatus = nf90_put_var (ncid, vpa_id, vpa)

        !Check store was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, vpa_id)

        !Store variable vperp2
        istatus = nf90_put_var (ncid, vperp2_id, vperp2)

        !Check store was successful
        IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, vperp2_id)
        !</DD>

        !Allocate array to hold components of fields individually
        IF (.not. ALLOCATED(ftmpr)) ALLOCATE (ftmpr(2*ntgrid+1,ntheta0,naky))

        !If finite electrostatic potential perturbation then write \phi
        IF (fphi > epsilon(0.)) THEN
            !Nake temporary copy of real(\phi)
            ftmpr = REAL(phinew)

            !Store variable real(\phi)
            istatus = nf90_put_var (ncid, phir_id, ftmpr)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, phir_id)

            !Make temporary copy of imaginary(\phi)
            ftmpr = AIMAG(phinew)

            !Store variable aimag(\phi)
            istatus = nf90_put_var (ncid, phii_id, ftmpr)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, phii_id)
        END IF

        !If finite parallel magnetic potential perturbations then write A_\parallel
        IF (fapar > epsilon(0.)) THEN
            !Make temporary copy of real(A_\parallel)
            ftmpr = REAL(aparnew)

            !Store variable real(A_\parallel)
            istatus = nf90_put_var (ncid, aparr_id, ftmpr)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, aparr_id)

            !Make temporary copy of imaginary(A_\parallel)
            ftmpr = AIMAG(aparnew)

            !Store variable imaginary(A_\parallel)
            istatus = nf90_put_var (ncid, apari_id, ftmpr)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, apari_id)
        END IF

        !If finite perpendicular magnetic potential perturbations then write B_\parallel
        IF (fbpar > epsilon(0.)) THEN
            !Make temporary copy of real(B_\parallel)
            ftmpr = REAL(bparnew)

            !Store variable real(B_\parallel)
            istatus = nf90_put_var (ncid, bparr_id, ftmpr)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, bparr_id)

            !Make temporary copy of imaginary(B_\parallel)
            ftmpr = AIMAG(bparnew)

            !Store variable imaginary(B_\parallel)
            istatus = nf90_put_var (ncid, bpari_id, ftmpr)

            !Check store was successful
            IF (istatus /= NF90_NOERR) CALL netcdf_error (istatus, ncid, bpari_id)
        END IF
    END IF !End of IF (n_elements > 0) THEN


    !Escape point in case of netcdf error at any point
1   CONTINUE

    !If an error has been encountered then close the file and exit routine
    IF (istatus /= NF90_NOERR) THEN
        i = nf90_close (ncid)
        RETURN
    END IF

    !If exit_in flag set then close the file, if not then sync the file
    IF (exit) THEN
        !Close file
        i = nf90_close (ncid)
    ELSE
        !Sync file
        i = nf90_sync (ncid)

        !Check sync was successful
        IF (i /= NF90_NOERR) &
            CALL netcdf_error (istatus, message='nf90_sync error')
    END IF
!----------------------------------------------------------------------------------------------

# else
    !Statement to display if try to save distfn without netcdf available
    IF (proc0) WRITE (error_unit(),*) &
         'WARNING: gs2_save_distfn is called without netcdf library'
# endif

END SUBROUTINE gs2_save_distfn
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

END MODULE gs2_save

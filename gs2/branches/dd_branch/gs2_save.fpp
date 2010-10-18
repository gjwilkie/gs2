# include "define.inc"

module gs2_save

# ifdef NETCDF
!  use netcdf, only: NF90_FLOAT, NF90_DOUBLE
  use netcdf, only: NF90_NOWRITE, NF90_CLOBBER, NF90_NOERR
  use netcdf, only: nf90_create, nf90_open, nf90_sync, nf90_close
  use netcdf, only: nf90_def_dim, nf90_def_var, nf90_enddef
  use netcdf, only: nf90_put_var, nf90_get_var, nf90_strerror
  use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid, nf90_inquire_variable
  
  use netcdf_utils, only: get_netcdf_code_precision
  use netcdf_utils, only: check_netcdf_file_precision
  use netcdf_utils, only: netcdf_error
  use netcdf_utils, only: netcdf_real, kind_nf
# endif

  implicit none

  public :: gs2_restore, gs2_save_for_restart
  public :: init_save, init_dt, init_tstart, init_ant_amp
  public :: init_vnm
!# ifdef NETCDF
!  public :: netcdf_real, kind_nf, get_netcdf_code_precision, netcdf_error
!# endif

  interface gs2_restore
     module procedure gs2_restore_many, gs2_restore_one
  end interface
  
  private
  character (300), save :: restart_file

# ifdef NETCDF
  real, allocatable, dimension(:,:,:) :: tmpr, tmpi, ftmpr, ftmpi
  real, allocatable, dimension(:) :: stmp  ! MR: tmp var for kx_shift
  real, allocatable, dimension(:) :: atmp
!  integer, parameter :: kind_nf = kind (NF90_NOERR)
  integer (kind_nf) :: ncid, thetaid, signid, gloid, kyid, kxid, nk_stir_dim
  integer (kind_nf) :: phir_id, phii_id, aparr_id, apari_id, bparr_id, bpari_id
  integer (kind_nf) :: kx_shift_id   ! MR: added to save kx_shift variable
  integer (kind_nf) :: t0id, gr_id, gi_id, vnm1id, vnm2id, delt0id
  integer (kind_nf) :: a_antr_id, b_antr_id, a_anti_id, b_anti_id
!<DD> Added for saving distribution function
  INTEGER (KIND_NF) :: egridid,lgridid
  INTEGER (KIND_NF) :: energy_id, lambda_id
  INTEGER (KIND_NF) :: nspecid, spec_id
  LOGICAL :: initialized_dfn= .false.
!</DD> Added for saving distribution function
!  integer (kind_nf) :: netcdf_real=0

  logical :: initialized = .false.
  logical :: test = .false.
# endif

contains

  subroutine gs2_save_for_restart &
       (g, t0, delt0, vnm, istatus, fphi, fapar, fbpar, exit_in,distfn)
!<DD 18-10-2010> Added flag distfn to gs2_save_for_restart
!If present then will save to "rootname.nc.dfn.proc" and will
!include extra information including velocity and energy grids
!</DD>
 
!MR, 2007: save kx_shift array in restart file if allocated    
# ifdef NETCDF
    use constants, only: kind_rs, kind_rd, pi
    use fields_arrays, only: phinew, aparnew, bparnew
    use dist_fn_arrays, only: kx_shift  !MR
    use kt_grids, only: naky, ntheta0
    use antenna_data, only: nk_stir, a_ant, b_ant, ant_on
# endif    
    use mp, only: nproc, iproc, proc0
    use theta_grid, only: ntgrid
! Must include g_layout_type here to avoid obscure bomb while compiling
! gs2_diagnostics.f90 (which uses this module) with the Compaq F90 compiler:
! TT>
!    use gs2_layouts, only: g_lo, g_layout_type
    use gs2_layouts, only: g_lo
    use layouts_type, only: g_layout_type
! <TT
    use file_utils, only: error_unit
    !<DD> Added for saving distribution function
    USE LE_GRIDS, ONLY: e,al,negrid,nlambda
    USE SPECIES, ONLY: nspec
    USE DIST_FN_ARRAYS, ONLY: vpa,vperp2
    !</DD> Added for saving distribution function
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, intent (in) :: t0, delt0
    real, dimension (2), intent (in) :: vnm
    real, intent (in) :: fphi, fapar, fbpar
    integer, intent (out) :: istatus
    logical, intent (in), optional :: exit_in
    LOGICAL, INTENT (in), optional :: distfn !<DD> Added for saving distribution function
# ifdef NETCDF
    character (306) :: file_proc
    character (10) :: suffix
    integer :: i, n_elements, ierr
    logical :: exit
    logical :: local_init !<DD> Added for saving distribution function
    if (present(exit_in)) then
       exit = exit_in
    else
       exit = .false.
    end if

    n_elements = g_lo%ulim_proc-g_lo%llim_proc+1
    if (n_elements <= 0) return

    !<DD> Added for saving distribution function
    IF (PRESENT(distfn)) THEN
    	local_init=initialized_dfn
    ELSE
    	local_init=initialized
    END IF
    !</DD> Added for saving distribution function
    
    if (.not.local_init) then

       !<DD> Added for saving distribution function
       IF (PRESENT(distfn)) THEN
        	initialized_dfn=.true.
       ELSE
    		initialized = .true.
       END IF
       !</DD> Added for saving distribution function
       
       file_proc = trim(restart_file)
       
       !<DD> Added for saving distribution function
       IF (PRESENT(distfn)) THEN
        	WRITE (suffix,'(a5,i0)') '.dfn.', iproc	
       ELSE
            WRITE (suffix,'(a1,i0)') '.', iproc
       END IF
       !</DD> Added for saving distribution function
       

       write (suffix,'(a1,i0)') '.', iproc
       file_proc = trim(trim(file_proc)//adjustl(suffix))
       
       istatus = nf90_create (file_proc, NF90_CLOBBER, ncid)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_create error: ", nf90_strerror(istatus)
          goto 1
       end if
       
       if (n_elements > 0) then
          istatus = nf90_def_dim (ncid, "theta", 2*ntgrid+1, thetaid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim theta error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          istatus = nf90_def_dim (ncid, "sign", 2, signid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim sign error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          istatus = nf90_def_dim (ncid, "glo", n_elements, gloid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim glo error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          istatus = nf90_def_dim (ncid, "aky", naky, kyid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim aky error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          istatus = nf90_def_dim (ncid, "akx", ntheta0, kxid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim akx error: ", nf90_strerror(istatus)
             goto 1
          end if
       end if
       
       !<DD> Added for saving distribution function
       	IF (PRESENT(distfn)) THEN
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
	        END IF 
       !</DD> Added for saving distribution function

       if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()

       istatus = nf90_def_var (ncid, "t0", netcdf_real, t0id)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_def_var t0 error: ", nf90_strerror(istatus)
          goto 1
       end if
       
       istatus = nf90_def_var (ncid, "delt0", netcdf_real, delt0id)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_def_var delt0 error: ", nf90_strerror(istatus)
          goto 1
       end if
       
       istatus = nf90_def_var (ncid, "vnm1", netcdf_real, vnm1id)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_def_var vnm(1) error: ", nf90_strerror(istatus)
          goto 1
       end if
       
       istatus = nf90_def_var (ncid, "vnm2", netcdf_real, vnm2id)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_def_var vnm(2) error: ", nf90_strerror(istatus)
          goto 1
       end if
       
       if (ant_on) then
          istatus = nf90_def_dim (ncid, "nk_stir", nk_stir, nk_stir_dim)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim nk_stir error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          istatus = nf90_def_var (ncid, "a_ant_r", netcdf_real, nk_stir_dim, a_antr_id)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_var a_ant_r error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          istatus = nf90_def_var (ncid, "a_ant_i", netcdf_real, nk_stir_dim, a_anti_id)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_var a_ant_i error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          istatus = nf90_def_var (ncid, "b_ant_r", netcdf_real, nk_stir_dim, b_antr_id)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_var b_ant_r error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          istatus = nf90_def_var (ncid, "b_ant_i", netcdf_real, nk_stir_dim, b_anti_id)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_var b_ant_i error: ", nf90_strerror(istatus)
             goto 1
          end if
       end if
       
       if (n_elements > 0) then
          istatus = nf90_def_var (ncid, "gr", netcdf_real, &
               (/ thetaid, signid, gloid /), gr_id)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_var g error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          istatus = nf90_def_var (ncid, "gi", netcdf_real, &
               (/ thetaid, signid, gloid /), gi_id)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_var g error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          if (fphi > epsilon(0.)) then
             istatus = nf90_def_var (ncid, "phi_r", netcdf_real, &
                  (/ thetaid, kxid, kyid /), phir_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var phi error: ", nf90_strerror(istatus)
                goto 1
             end if
             
             istatus = nf90_def_var (ncid, "phi_i", netcdf_real, &
                  (/ thetaid, kxid, kyid /), phii_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var phi error: ", nf90_strerror(istatus)
                goto 1
             end if
          end if

          if (fapar > epsilon(0.)) then
             istatus = nf90_def_var (ncid, "apar_r", netcdf_real, &
                  (/ thetaid, kxid, kyid /), aparr_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var apar error: ", nf90_strerror(istatus)
                goto 1
             end if
             
             istatus = nf90_def_var (ncid, "apar_i", netcdf_real, &
                  (/ thetaid, kxid, kyid /), apari_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var apar error: ", nf90_strerror(istatus)
                goto 1
             end if
          end if

          if (fbpar > epsilon(0.)) then
             istatus = nf90_def_var (ncid, "bpar_r", netcdf_real, &
                  (/ thetaid, kxid, kyid /), bparr_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var bparr error: ", nf90_strerror(istatus)
                goto 1
             end if

             istatus = nf90_def_var (ncid, "bpar_i", netcdf_real, &
                  (/ thetaid, kxid, kyid /), bpari_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var bpari error: ", nf90_strerror(istatus)
                goto 1
             end if
          end if

          !<DD> Added for saving distribution function
           IF (PRESENT(distfn)) THEN
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
            END IF
            !</DD> Added for saving distribution function
       end if

! remove allocated conditional because we want to be able to restart
! using exb shear from a case which does not have exb shear (i.e.
! we need kx_shift variable defined in netcdf file even if no exb
! shear present in simulation) -- MAB + CMR
!       if (allocated(kx_shift)) then   ! MR begin
       istatus = nf90_def_var (ncid, "kx_shift", netcdf_real, &
            (/ kyid /), kx_shift_id)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_def_var kx_shift error: ", nf90_strerror(istatus)
          goto 1
       endif
          
!       endif   ! MR end 
       
       istatus = nf90_enddef (ncid)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write (ierr,*) "nf90_enddef error: ", nf90_strerror(istatus)
          goto 1
       end if
    end if

    istatus = nf90_put_var (ncid, delt0id, delt0)
    if (istatus /= NF90_NOERR) then
       ierr = error_unit()
       write (ierr,*) "nf90_put_var delt0 error: ", nf90_strerror(istatus)
       goto 1
    end if
 
    istatus = nf90_put_var (ncid, t0id, t0)
    if (istatus /= NF90_NOERR) then
       ierr = error_unit()
       write (ierr,*) "nf90_put_var t0 error: ", nf90_strerror(istatus)
       goto 1
    end if
 
    istatus = nf90_put_var (ncid, vnm1id, vnm(1))
    if (istatus /= NF90_NOERR) then
       ierr = error_unit()
       write (ierr,*) "nf90_put_var vnm(1) error: ", nf90_strerror(istatus)
       goto 1
    end if
 
    istatus = nf90_put_var (ncid, vnm2id, vnm(2))
    if (istatus /= NF90_NOERR) then
       ierr = error_unit()
       write (ierr,*) "nf90_put_var vnm(2) error: ", nf90_strerror(istatus)
       goto 1
    end if

1   continue
    if (istatus /= NF90_NOERR) then
       i = nf90_close (ncid)
       return
    end if

    if (n_elements > 0) then

       if (ant_on) then

          if (.not. allocated(atmp)) allocate (atmp(nk_stir))
          atmp = real(a_ant)
          istatus = nf90_put_var (ncid, a_antr_id, atmp)

          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write (ierr,*) "nf90_put_var a_antr error: ", &
                  nf90_strerror(istatus), ' ', iproc
          end if

          atmp = aimag(a_ant)
          istatus = nf90_put_var (ncid, a_anti_id, atmp)

          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write (ierr,*) "nf90_put_var a_anti error: ", &
                  nf90_strerror(istatus), ' ', iproc
          end if

          atmp = real(b_ant)
          istatus = nf90_put_var (ncid, b_antr_id, atmp)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write (ierr,*) "nf90_put_var b_antr error: ", &
                  nf90_strerror(istatus), ' ', iproc
          end if

          atmp = aimag(b_ant)
          istatus = nf90_put_var (ncid, b_anti_id, atmp)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write (ierr,*) "nf90_put_var b_anti error: ", &
                  nf90_strerror(istatus), ' ', iproc
          end if
          deallocate (atmp)
       end if

       if (.not. allocated(tmpr)) &
            allocate (tmpr(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))

       tmpr = real(g)
       istatus = nf90_put_var (ncid, gr_id, tmpr)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gr_id)
       
       tmpr = aimag(g)
       istatus = nf90_put_var (ncid, gi_id, tmpr)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gi_id)
       
       !<DD> Added for saving distribution function
       IF (PRESENT(distfn)) THEN
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
       END IF
       !</DD> Added for saving distribution function

       if (.not. allocated(ftmpr)) allocate (ftmpr(2*ntgrid+1,ntheta0,naky))
       if (.not. allocated(ftmpi)) allocate (ftmpi(2*ntgrid+1,ntheta0,naky))

       if (fphi > epsilon(0.)) then
          ftmpr = real(phinew)
          istatus = nf90_put_var (ncid, phir_id, ftmpr)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phir_id)
          
          ftmpi = aimag(phinew)
          istatus = nf90_put_var (ncid, phii_id, ftmpi)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phii_id)
       end if

       if (fapar > epsilon(0.)) then
          ftmpr = real(aparnew)
          istatus = nf90_put_var (ncid, aparr_id, ftmpr)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, aparr_id)
          
          ftmpi = aimag(aparnew)
          istatus = nf90_put_var (ncid, apari_id, ftmpi)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, apari_id)
       end if

       if (fbpar > epsilon(0.)) then
          ftmpr = real(bparnew)
          istatus = nf90_put_var (ncid, bparr_id, ftmpr)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bparr_id)
          
          ftmpi = aimag(bparnew)
          istatus = nf90_put_var (ncid, bpari_id, ftmpi)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bpari_id)
       end if
    end if

    if (allocated(kx_shift)) then ! MR begin
       if (.not. allocated(stmp)) allocate (stmp(naky))   
       stmp = kx_shift
       istatus = nf90_put_var (ncid, kx_shift_id, stmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, kx_shift_id)
    else
       if (.not. allocated(stmp)) allocate (stmp(naky))
       stmp = 0.
       istatus = nf90_put_var (ncid, kx_shift_id, stmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, kx_shift_id)
    endif ! MR end

    if (exit) then
       i = nf90_close (ncid)
    else
       i = nf90_sync (ncid)
       if (i /= NF90_NOERR) &
            call netcdf_error (istatus, message='nf90_sync error')
    end if

# else

    if (proc0) write (error_unit(),*) &
         'WARNING: gs2_save_for_restart is called without netcdf library'

# endif

  end subroutine gs2_save_for_restart

  subroutine gs2_restore_many (g, scale, istatus, fphi, fapar, fbpar, many)
!MR, 2007: restore kx_shift array if already allocated
# ifdef NETCDF
    use mp, only: proc0, iproc, nproc
    use fields_arrays, only: phinew, aparnew, bparnew
    use fields_arrays, only: phi, apar, bpar
    use dist_fn_arrays, only: kx_shift   ! MR
    use kt_grids, only: naky, ntheta0
# endif
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use file_utils, only: error_unit
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    real, intent (in) :: scale
    integer, intent (out) :: istatus
    real, intent (in) :: fphi, fapar, fbpar
    logical, intent (in) :: many
# ifdef NETCDF
    character (306) :: file_proc
    character (10) :: suffix
    integer :: i, n_elements, ierr
    real :: fac

    n_elements = g_lo%ulim_proc-g_lo%llim_proc+1
    if (n_elements <= 0) return
    
    if (.not.initialized) then
!       initialized = .true.
       file_proc = trim(restart_file)

       write (suffix,'(a1,i0)') '.', iproc
       file_proc = trim(trim(file_proc)//adjustl(suffix))
       
       istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, file=file_proc)
       
       ! check precision
       if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()
       call check_netcdf_file_precision (ncid)

       istatus = nf90_inq_dimid (ncid, "theta", thetaid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='theta')
       
       istatus = nf90_inq_dimid (ncid, "sign", signid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='sign')
       
       istatus = nf90_inq_dimid (ncid, "glo", gloid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='glo')
       
       istatus = nf90_inq_dimid (ncid, "aky", kyid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='aky')
       
       istatus = nf90_inq_dimid (ncid, "akx", kxid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='akx')
       
       istatus = nf90_inquire_dimension (ncid, thetaid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=thetaid)
       if (i /= 2*ntgrid + 1) write(*,*) 'Restart error: ntgrid=? ',i,' : ',ntgrid,' : ',iproc
       
       istatus = nf90_inquire_dimension (ncid, signid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=signid)
       if (i /= 2) write(*,*) 'Restart error: sign=? ',i,' : ',iproc
       
       istatus = nf90_inquire_dimension (ncid, gloid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=gloid)
       if (i /= g_lo%ulim_proc-g_lo%llim_proc+1) write(*,*) 'Restart error: glo=? ',i,' : ',iproc
       
       istatus = nf90_inquire_dimension (ncid, kyid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=kyid)
       if (i /= naky) write(*,*) 'Restart error: naky=? ',i,' : ',naky,' : ',iproc
       
       istatus = nf90_inquire_dimension (ncid, kxid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=kxid)
       if (i /= ntheta0) write(*,*) 'Restart error: ntheta0=? ',i,' : ',ntheta0,' : ',iproc
       
       if (fphi > epsilon(0.)) then
          istatus = nf90_inq_varid (ncid, "phi_r", phir_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='phi_r')
          
          istatus = nf90_inq_varid (ncid, "phi_i", phii_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='phi_i')
       end if

       if (fapar > epsilon(0.)) then
          istatus = nf90_inq_varid (ncid, "apar_r", aparr_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='apar_r')
          
          istatus = nf90_inq_varid (ncid, "apar_i", apari_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='apar_i')
       end if

       if (fbpar > epsilon(0.)) then
          istatus = nf90_inq_varid (ncid, "bpar_r", bparr_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='bpar_r')

          istatus = nf90_inq_varid (ncid, "bpar_i", bpari_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='bpar_i')
       end if

       if (allocated(kx_shift)) then   ! MR begin
          istatus = nf90_inq_varid (ncid, "kx_shift", kx_shift_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='kx_shift')
       endif   ! MR end

       istatus = nf90_inq_varid (ncid, "gr", gr_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='gr')
       
       istatus = nf90_inq_varid (ncid, "gi", gi_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='gi')
    end if
    
    if (.not. allocated(tmpr)) &
         allocate (tmpr(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))
    if (.not. allocated(tmpi)) &
         allocate (tmpi(2*ntgrid+1,2,g_lo%llim_proc:g_lo%ulim_alloc))

    tmpr = 0.; tmpi = 0.
    istatus = nf90_get_var (ncid, gr_id, tmpr)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gr_id)

    istatus = nf90_get_var (ncid, gi_id, tmpi)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gi_id)

    g = cmplx(tmpr, tmpi)

    if (.not. allocated(ftmpr)) allocate (ftmpr(2*ntgrid+1,ntheta0,naky))
    if (.not. allocated(ftmpi)) allocate (ftmpi(2*ntgrid+1,ntheta0,naky))

    if (allocated(kx_shift)) then   ! MR begin
       if (.not. allocated(stmp)) allocate (stmp(naky))   ! MR 
       istatus = nf90_get_var (ncid, kx_shift_id, stmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, kx_shift_id)
       kx_shift = stmp
    endif   ! MR end

    if (fphi > epsilon(0.)) then
       istatus = nf90_get_var (ncid, phir_id, ftmpr)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phir_id)
       
       istatus = nf90_get_var (ncid, phii_id, ftmpi)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phii_id)
       
       phi = 0.
       phinew = cmplx(ftmpr, ftmpi)
    end if

    if (fapar > epsilon(0.)) then
       istatus = nf90_get_var (ncid, aparr_id, ftmpr)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, aparr_id)
       
       istatus = nf90_get_var (ncid, apari_id, ftmpi)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, apari_id)
       
       apar = 0.
       aparnew = cmplx(ftmpr, ftmpi)
    end if

    if (fbpar > epsilon(0.)) then
       istatus = nf90_get_var (ncid, bparr_id, ftmpr)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bparr_id)
       
       istatus = nf90_get_var (ncid, bpari_id, ftmpi)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bpari_id)
       
       bpar = 0.
       bparnew = cmplx(ftmpr, ftmpi)
    end if

    if (scale > 0.) then
       g = g*scale
       phinew = phinew*scale
       aparnew = aparnew*scale
       bparnew = bparnew*scale
    else
       fac = - scale/(maxval(abs(phinew)))
       g = g*fac
       phinew = phinew*fac
       aparnew = aparnew*fac
       bparnew = bparnew*fac
    end if

    ! RN 2008/05/23: this was commented out. why?
!    istatus = nf90_close (ncid)
    if (istatus /= NF90_NOERR) then
       ierr = error_unit()
       write(ierr,*) "nf90_close error: ", nf90_strerror(istatus),' ',iproc
    end if

# else
    
    write (error_unit(),*) &
         'ERROR: gs2_restore_many is called without netcdf'

# endif

  end subroutine gs2_restore_many

! RN 2008/05/23:
!  This can be removed. restore_many seems to work for single proc.
  subroutine gs2_restore_one (g, scale, istatus, fphi, fapar, fbpar)
!MR, 2007: restore kx_shift array if allocated
# ifdef NETCDF
    use mp, only: proc0, iproc, nproc
    use fields_arrays, only: phinew, aparnew, bparnew
    use fields_arrays, only: phi, apar, bpar
    use dist_fn_arrays, only: kx_shift   ! MR
    use kt_grids, only: naky, ntheta0
# endif
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use file_utils, only: error_unit
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    real, intent (in) :: scale
    integer, intent (out) :: istatus
    real, intent (in) :: fphi, fapar, fbpar
# ifdef NETCDF
    integer :: n_elements, ierr
    integer :: i

    istatus = nf90_open (trim(restart_file), NF90_NOWRITE, ncid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, file=restart_file)

    istatus = nf90_inq_dimid (ncid, "theta", thetaid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='theta')

    istatus = nf90_inq_dimid (ncid, "sign", signid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='sign')

    istatus = nf90_inq_dimid (ncid, "glo", gloid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='glo')

    istatus = nf90_inq_dimid (ncid, "aky", kyid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='aky')

    istatus = nf90_inq_dimid (ncid, "akx", kxid)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='akx')

    istatus = nf90_inquire_dimension (ncid, thetaid, len=i)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=thetaid)
    if (i /= 2*ntgrid + 1) write(*,*) 'Restart error: ntgrid=? ',i,' : ',ntgrid,' : ',iproc

    istatus = nf90_inquire_dimension (ncid, signid, len=i)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, signid)
    if (i /= 2) write(*,*) 'Restart error: sign=? ',i,' : ',iproc

    istatus = nf90_inquire_dimension (ncid, gloid, len=i)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gloid)
    if (i /= g_lo%ulim_world-g_lo%llim_world+1) write(*,*) 'Restart error: glo=? ',i,' : ',iproc

    istatus = nf90_inquire_dimension (ncid, kyid, len=i)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, kyid)
    if (i /= naky) write(*,*) 'Restart error: naky=? ',i,' : ',naky,' : ',iproc

    istatus = nf90_inquire_dimension (ncid, kxid, len=i)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, kxid)
    if (i /= ntheta0) write(*,*) 'Restart error: ntheta0=? ',i,' : ',ntheta0,' : ',iproc

    if (fphi > epsilon(0.)) then
       istatus = nf90_inq_varid (ncid, "phi_r", phir_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='phi_r')
       
       istatus = nf90_inq_varid (ncid, "phi_i", phii_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='phi_i')
    end if

    if (fapar > epsilon(0.)) then
       istatus = nf90_inq_varid (ncid, "apar_r", aparr_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='apar_r')
       
       istatus = nf90_inq_varid (ncid, "apar_i", apari_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='apar_i')
    end if

    if (fbpar > epsilon(0.)) then
       istatus = nf90_inq_varid (ncid, "bpar_r", bparr_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='bpar_r')

       istatus = nf90_inq_varid (ncid, "bpar_i", bpari_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='bpar_i')
    end if

    if (allocated(kx_shift)) then   ! MR begin
       istatus = nf90_inq_varid (ncid, "kx_shift", kx_shift_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='kx_shift')
    endif   ! MR end

    istatus = nf90_inq_varid (ncid, "gr", gr_id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='gr')

    istatus = nf90_inq_varid (ncid, "gi", gi_id)
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='gi')

    n_elements = g_lo%ulim_proc - g_lo%llim_proc + 1

    if (n_elements <= 0) then
       phinew = 0.
       aparnew = 0.
       bparnew = 0.
       return
    endif

    if (.not. allocated(tmpr)) allocate (tmpr(2*ntgrid+1, 2, n_elements))    
    if (.not. allocated(tmpi)) allocate (tmpi(2*ntgrid+1, 2, n_elements))

    tmpr = 0.;  tmpi = 0.
    istatus = nf90_get_var (ncid, gr_id, tmpr, &
         start=(/ 1, 1, g_lo%llim_proc /), &
         count=(/ 2*ntgrid + 1, 2, n_elements /))
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gr_id)
    
    istatus = nf90_get_var (ncid, gi_id, tmpi, &
         start=(/ 1, 1, g_lo%llim_proc /), &
         count=(/ 2*ntgrid + 1, 2, n_elements /))
    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gi_id)

    g = cmplx(tmpr, tmpi)*scale

    if (.not. allocated(ftmpr)) allocate (ftmpr(2*ntgrid+1,ntheta0,naky))
    if (.not. allocated(ftmpi)) allocate (ftmpi(2*ntgrid+1,ntheta0,naky))

    if (allocated(kx_shift)) then   ! MR begin
       if (.not. allocated(stmp)) allocate (stmp(naky))   ! MR
       istatus = nf90_get_var (ncid, kx_shift_id, stmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, kx_shift_id)
       kx_shift = stmp
    endif   ! MR end


    if (fphi > epsilon(0.)) then
       istatus = nf90_get_var (ncid, phir_id, ftmpr)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phir_id)
       
       istatus = nf90_get_var (ncid, phii_id, ftmpi)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phii_id)
       
       phi = 0.
       phinew = cmplx(ftmpr, ftmpi)*scale
    end if

    if (fapar > epsilon(0.)) then
       istatus = nf90_get_var (ncid, aparr_id, ftmpr)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, aparr_id)
       
       istatus = nf90_get_var (ncid, apari_id, ftmpi)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, apari_id)
       
       apar = 0.
       aparnew = cmplx(ftmpr, ftmpi)*scale
    end if

    if (fbpar > epsilon(0.)) then
       istatus = nf90_get_var (ncid, bparr_id, ftmpr)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bparr_id)
       
       istatus = nf90_get_var (ncid, bpari_id, ftmpi)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, bpari_id)
       
       bpar = 0.
       bparnew = cmplx(ftmpr, ftmpi)*scale
    end if

    ! RN 2008/05/23: this was commented out. why?
!    istatus = nf90_close (ncid)       
    if (istatus /= NF90_NOERR) then
       ierr = error_unit()
       write(ierr,*) "nf90_close error: ", nf90_strerror(istatus),' ',iproc
    end if

# else

    write (error_unit(),*) &
         'ERROR: gs2_restore_one is called without netcdf'

# endif

  end subroutine gs2_restore_one

  subroutine init_save (file)
    character(300), intent (in) :: file
    
    restart_file = file

  end subroutine init_save

  subroutine init_dt (delt0, istatus)

# ifdef NETCDF
    use mp, only: nproc, proc0, broadcast
    use file_utils, only: error_unit
# endif
    implicit none
    real, intent (in out) :: delt0
    integer, intent (out) :: istatus
# ifdef NETCDF
    character (306) :: file_proc

    if (proc0) then

       if (.not. initialized) then

          file_proc=trim(trim(restart_file)//'.0')

          istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus,file=file_proc)

          istatus = nf90_inq_varid (ncid, "delt0", delt0id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='delt0')
       end if

       istatus = nf90_get_var (ncid, delt0id, delt0)

       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, delt0id, message=' in init_dt')
          delt0 = -1.
       endif           

       if (.not.initialized) istatus = nf90_close (ncid)
    endif

    call broadcast (istatus)
    call broadcast (delt0)

# endif

  end subroutine init_dt

  subroutine init_vnm (vnm, istatus)

# ifdef NETCDF
    use mp, only: nproc, proc0, broadcast
    use file_utils, only: error_unit
# endif
    implicit none
    real, dimension(2), intent (in out) :: vnm
    integer, intent (out) :: istatus
# ifdef NETCDF
    character (306) :: file_proc

    if (proc0) then
       if (.not.initialized) then

          file_proc=trim(trim(restart_file)//'.0')

          istatus = nf90_open (file_proc, 0, ncid)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, file=file_proc)

          istatus = nf90_inq_varid (ncid, "vnm1", vnm1id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='vnm1')

          istatus = nf90_inq_varid (ncid, "vnm2", vnm2id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='vnm2')
       end if

       istatus = nf90_get_var (ncid, vnm1id, vnm(1))

       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, vnm1id, message=' in init_vnm')
          vnm(1) = 0.
       endif           

       istatus = nf90_get_var (ncid, vnm2id, vnm(2))

       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, vnm2id, message=' in init_vnm')
          vnm(2) = 0.
       endif           

       if (.not. initialized) istatus = nf90_close (ncid)
    endif

    call broadcast (istatus)
    call broadcast (vnm)

# endif

  end subroutine init_vnm

! This routine gets a_ant and b_ant for proc 0 only!!
  subroutine init_ant_amp (a_ant, b_ant, nk_stir, istatus)

# ifdef NETCDF
    use file_utils, only: error_unit
    use constants, only: zi
# endif
    use mp, only: proc0
    implicit none
    complex, dimension(:), intent (in out) :: a_ant, b_ant
    integer, intent (in) :: nk_stir
    integer, intent (out) :: istatus
# ifdef NETCDF
    character (306) :: file_proc
    integer :: ierr, i

    if (proc0) then
       a_ant = 0. ; b_ant = 0.

       file_proc=trim(trim(restart_file)//'.0')
       
       istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_open in init_ant_amp error: ", nf90_strerror(istatus) 
          write(ierr,*) "If you did not intend for this to be a restarted run with an external antenna,"
          write(ierr,*) "you may ignore the error message above."
          return
       endif

       istatus = nf90_inq_dimid (ncid, "nk_stir", nk_stir_dim)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='nk_stir')

       istatus = nf90_inquire_dimension (ncid, nk_stir_dim, len=i)
       if (istatus /= NF90_NOERR) &
            call netcdf_error (istatus, ncid, dimid=nk_stir_dim)
       if (i /= nk_stir) write(*,*) 'Restart error: nk_stir=? ',i,' : ',nk_stir

       istatus = nf90_inq_varid (ncid, "a_ant_r", a_antr_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='a_ant_r')

       istatus = nf90_inq_varid (ncid, "a_ant_i", a_anti_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='a_ant_i')

       istatus = nf90_inq_varid (ncid, "b_ant_r", b_antr_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='b_ant_r')

       istatus = nf90_inq_varid (ncid, "b_ant_i", b_anti_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='b_ant_i')

       if (.not. allocated(atmp)) allocate (atmp(nk_stir))
       atmp = 0.

       istatus = nf90_get_var (ncid, a_antr_id, atmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, a_antr_id)
       a_ant = atmp

       istatus = nf90_get_var (ncid, a_anti_id, atmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, a_anti_id)
       a_ant = a_ant + zi * atmp

       istatus = nf90_get_var (ncid, b_antr_id, atmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, b_antr_id)
       b_ant = atmp

       istatus = nf90_get_var (ncid, b_anti_id, atmp)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, b_anti_id)
       b_ant = b_ant + zi * atmp

       deallocate (atmp)
       istatus = nf90_close (ncid)       
    endif

# else

    if (proc0) istatus = 2

# endif

  end subroutine init_ant_amp

  subroutine init_tstart (tstart, istatus)

# ifdef NETCDF
    use mp, only: nproc, proc0, broadcast
    use file_utils, only: error_unit
# endif
    implicit none
    real, intent (in out) :: tstart
    integer, intent (out) :: istatus
# ifdef NETCDF
    character (306) :: file_proc
    integer :: ierr

    if (proc0) then
       
       if (.not.initialized) then

          file_proc=trim(trim(restart_file)//'.0')

          istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, file=file_proc)
       end if

       istatus = nf90_inq_varid (ncid, "t0", t0id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='t0')

       istatus = nf90_get_var (ncid, t0id, tstart)
       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, t0id, message=' in init_tstart')
          tstart = -1.
       end if           
       if (.not.initialized) istatus = nf90_close (ncid)
    endif

    call broadcast (istatus)
    call broadcast (tstart)

# endif

  end subroutine init_tstart
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

end module gs2_save

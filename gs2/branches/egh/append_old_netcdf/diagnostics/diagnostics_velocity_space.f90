!> A module which writes out quantities which writes out 
!! diagnostic quantities which assess whether the velocity 
!! space resolution is sufficient.
module diagnostics_velocity_space
  implicit none
  private

  public :: init_diagnostics_velocity_space, write_velocity_space_checks
  public :: write_collision_error
  public :: finish_diagnostics_velocity_space

contains
  
  subroutine init_diagnostics_velocity_space(gnostics)
    use le_grids, only: init_weights
    use mp, only: proc0
    use diagnostics_config, only: diagnostics_type
    use collisions, only: collision_model_switch
    use collisions, only: init_lorentz_error
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics 
    
    if (.not. gnostics%write_verr) return
    
    ! initialize weights for less accurate integrals used
    ! to provide an error estimate for v-space integrals (energy and untrapped)
    if (proc0) call init_weights
    
    if (gnostics%write_cerr) then
       if (collision_model_switch == 1 .or. collision_model_switch == 5) then
          call init_lorentz_error
       else
          gnostics%write_cerr = .false.
       end if
    end if
  end subroutine init_diagnostics_velocity_space

  subroutine finish_diagnostics_velocity_space(gnostics)
    use le_grids, only: finish_weights
    use mp, only: proc0
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics 
    if (proc0) call finish_weights
  end subroutine finish_diagnostics_velocity_space

  
  
  subroutine write_velocity_space_checks(gnostics)
    use dist_fn, only: get_verr, get_gtran
    use mp, only: proc0
    use le_grids, only: nlambda, ng2
    use fields_arrays, only: phinew, bparnew
    use gs2_time, only: user_time
    use collisions, only: vnmult
    use species, only: spec
    use diagnostics_config, only: diagnostics_type
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_create_and_write, only: create_and_write_variable_noread
    use diagnostics_dimensions, only: dim_string
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real, dimension (:,:), allocatable :: errest
    integer, dimension (:,:), allocatable :: erridx
    real :: geavg, glavg, gtavg
    
    allocate(errest(5,2), erridx(5,3))
    errest = 0.0; erridx = 0
    
    if (.not. gnostics%replay) then
      ! error estimate obtained by comparing standard integral with less-accurate integral
      call get_verr (errest, erridx, phinew, bparnew)
      
      ! error estimate based on monitoring amplitudes of legendre polynomial coefficients
      call get_gtran (geavg, glavg, gtavg, phinew, bparnew)
    end if

    
    !errest(5,:) = -4
    
    if (.not. gnostics%vary_vnew_only) then
       call create_and_write_variable_noread(gnostics, gnostics%rtype, "vspace_lpcfrac", &
          dim_string([gnostics%dims%generic_3,gnostics%dims%time]), &
          "Fraction of free energy contained in the high order coefficients of &
          & the Legendre polynomial transform of (1) energy space, (2) untrapped &
          & pitch angles and (3) trapped pitch angles (each should ideally be < 0.1).  &
          & Note that there are no trapped pitch angles for certain geometries", &
          "1", (/geavg, glavg, gtavg/))
       call create_and_write_variable(gnostics, gnostics%rtype, "vspace_err", &
          dim_string([gnostics%dims%generic_5,gnostics%dims%generic_2,gnostics%dims%time]), &
          "Estimate of the (1) absolute and (2) relative errors resulting from &
          & velocity space integrals in the calculation of the following quantities &
          & in the given dimensions: (1) k phi, energy (2) k phi, untrapped pitch angles &
          & (3) k phi, trapped pitch angles, (4) k apar, energy, (5) k apar, untrapped &
          & angles. Relative errors should be < 0.1. ", &
          "absolute error measures have units T_r/(e rho_r)", errest)
       call create_and_write_variable_noread(gnostics, gnostics%rtype, "vspace_vnewk", &
          dim_string([gnostics%dims%generic_2,gnostics%dims%time]), &
          "If the simulation is set to vary the collisionality in order to keep &
          & error in velocity integrals to acceptable levels, contains species 1 &
          & collisionality in (1) pitch angle and (2) energy  ", &
          "v_thr/a", (/vnmult(1)*spec(1)%vnewk, vnmult(2)*spec(1)%vnewk/))
       ! This next statement causes annoying printout because of an error in netcdf 4.1
       ! The error is fixed in 4.2
       ! See
       ! http://www.unidata.ucar.edu/software/netcdf/docs/known_problems.html#f90-debug-segfault
      if (gnostics%write_max_verr) &
           call create_and_write_variable(gnostics, gnostics%itype, "vspace_err_maxindex", &
            dim_string([gnostics%dims%generic_5,gnostics%dims%generic_3,gnostics%dims%time]), &
           "Gives the (1) theta index, (2) ky index and (3) kx index of the maximum &
           & error resulting from the &
           & velocity space integrals in the calculation of the following quantities &
           & in the given dimensions: (1) k phi, energy (2) k phi, untrapped pitch angles &
           & (3) k phi, trapped pitch angles, (4) k apar, energy, (5) k apar, untrapped &
           & angles. Relative errors should be < 0.1. ", &
           "1", erridx)
   end if
   
   if (proc0 .and. gnostics%write_ascii) call write_ascii
   deallocate(errest,erridx)
   
 contains
   subroutine write_ascii
     if (nlambda - ng2 > 1) then
        write(gnostics%ascii_files%lpc,"(4(1x,e13.6))") user_time, geavg, glavg, gtavg
     else
        write(gnostics%ascii_files%lpc,"(3(1x,e13.6))") user_time, geavg, glavg
     end if
     write(gnostics%ascii_files%vres,"(8(1x,e13.6))") user_time, errest(1,2), errest(2,2), errest(3,2), &
          errest(4,2), errest(5,2), vnmult(1)*spec(1)%vnewk, vnmult(2)*spec(1)%vnewk
     if (gnostics%write_max_verr) then
        write(gnostics%ascii_files%vres2,"(3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6),3(i8),(1x,e13.6))") &
             erridx(1,1), erridx(1,2), erridx(1,3), errest(1,1), &
             erridx(2,1), erridx(2,2), erridx(2,3), errest(2,1), &
             erridx(3,1), erridx(3,2), erridx(3,3), errest(3,1), &
             erridx(4,1), erridx(4,2), erridx(4,3), errest(4,1), &
             erridx(5,1), erridx(5,2), erridx(5,3), errest(5,1)
     end if
   end subroutine write_ascii
 end subroutine write_velocity_space_checks

 subroutine write_collision_error (gnostics)
    
    use mp, only: proc0, send, receive, barrier
    use le_grids, only: ng2, jend, nlambda, lambda_map
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gnew, g_adjust
    use fields_arrays, only: phinew, bparnew
    use run_parameters, only: fphi, fbpar
    use gs2_layouts, only: lz_lo, ig_idx, idx_local, proc_id
    use gs2_layouts, only: ik_idx, ie_idx, is_idx, it_idx, il_idx, g_lo
    use collisions, only: dtot, fdf, fdb
    use redistribute, only: gather, scatter
    use file_utils, only: open_output_file, close_output_file
    use gs2_time, only: user_time
    use diagnostics_config, only: diagnostics_type
    implicit none

    type(diagnostics_type), intent(in) :: gnostics


    integer :: je, te, ig, il, ip, ilz, ie, is, ik, it
    integer :: igmax, ikmax, itmax, iemax, ilmax, ismax
    integer, save :: unit
    complex, dimension (:), allocatable :: ltmp, ftmp
    complex, dimension (:,:), allocatable :: lcoll, fdcoll, glze
    complex, dimension (:,:,:), allocatable :: g0
!    logical :: first = .true.
    real :: etmp, emax, etot, eavg, edenom, ltmax
    real :: time

    allocate (ltmp(2*nlambda), ftmp(2*nlambda))
    allocate (lcoll(2*nlambda,lz_lo%llim_proc:lz_lo%ulim_alloc))
    allocate (fdcoll(2*nlambda,lz_lo%llim_proc:lz_lo%ulim_alloc))
    allocate (glze(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
    allocate (g0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))

    lcoll = 0.0; fdcoll = 0.0; glze = 0.0; ltmp = 0.0; ftmp = 0.0
    etmp = 0.0; emax = 0.0; etot = 0.0; eavg = 0.0; edenom = 1.0; ltmax = 1.0

!    if (first .and. proc0) then
    !if (.not.cerrinit .and. proc0) then
!!       first = .false.
       !cerrinit = .true.
    !end if

! convert gnew from g to h
    call g_adjust(gnew,phinew,bparnew,fphi,fbpar)

    g0 = gnew

! convert gnew from h back to g
    call g_adjust(gnew,phinew,bparnew,-fphi,-fbpar)

! map from g0(ig,isgn,iglo) to glze(il,ilz)
    call gather (lambda_map, g0, glze)

! loop over ig, isign, ik, it, ie, is
    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
       ig = ig_idx(lz_lo,ilz)
       
       if (jend(ig) == 0) then      ! no trapped particles
! je = number of + vpa grid pts
! te = number of + and - vpa grid pts
          je = ng2
          te = 2*ng2
          
! find d/d(xi) ((1+xi**2)( d g(xi)/ d(xi) )) at each xi
! using lagrange (lcoll) and finite difference (fdcoll)
          il = 1
          do ip = il, il+2
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip)*glze(ip,ilz)
          end do

          il = 2
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+2)*glze(ip,ilz)
          end do

          do il=3,ng2
             do ip=il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+3)*glze(ip,ilz)
             end do
          end do

          do il=ng2+1, 2*ng2-2
             do ip = il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2*ng2-il+1,il-ip+3)*glze(ip,ilz)
             end do
          end do

          il = 2*ng2-1
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2,il-ip+2)*glze(ip,ilz)
          end do

          il = 2*ng2
          do ip = il-2, il
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,1,il-ip+1)*glze(ip,ilz)
          end do

! deal with xi from 1-eps -> eps
          do il=2,ng2
             fdcoll(il,ilz) = fdf(ig,il)*(glze(il+1,ilz) - glze(il,ilz)) - fdb(ig,il)*(glze(il,ilz) - glze(il-1,ilz))
          end do
! deal with xi from -eps -> -1+eps
          do il=ng2+1, 2*ng2-1
             fdcoll(il,ilz) = fdb(ig,2*ng2-il+1)*(glze(il+1,ilz) - glze(il,ilz)) - fdf(ig,2*ng2-il+1)*(glze(il,ilz) - glze(il-1,ilz))
          end do

          fdcoll(1,ilz) = fdf(ig,1)*(glze(2,ilz) - glze(1,ilz))
          fdcoll(2*ng2,ilz) = -fdf(ig,1)*(glze(2*ng2,ilz) - glze(2*ng2-1,ilz))

       else       ! trapped particle runs          
          je = jend(ig)
          te = 2*je - 1

          il = 1
          do ip = il, il+2
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip)*glze(ip,ilz)
          end do

          il = 2
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+2)*glze(ip,ilz)
          end do

          do il=3,je
             do ip=il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+3)*glze(ip,ilz)
             end do
          end do

          do il=je+1, te-2
             do ip = il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,te-il+1,il-ip+3)*glze(ip,ilz)
             end do
          end do

          il = te-1
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2,il-ip+2)*glze(ip,ilz)
          end do

          il = te
          do ip = il-2, il
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,1,il-ip+1)*glze(ip,ilz)
          end do

! is il=je handled correctly here?
          do il=2,je
             fdcoll(il,ilz) = fdf(ig,il)*(glze(il+1,ilz) - glze(il,ilz)) - fdb(ig,il)*(glze(il,ilz) - glze(il-1,ilz))
          end do
          do il=je+1,te-1
             fdcoll(il,ilz) = fdb(ig,te-il+1)*(glze(il+1,ilz) - glze(il,ilz)) - fdf(ig,te-il+1)*(glze(il,ilz) - glze(il-1,ilz))
          end do
          
          fdcoll(1,ilz) = fdf(ig,1)*(glze(2,ilz) - glze(1,ilz))
          fdcoll(te,ilz) = -fdf(ig,1)*(glze(te,ilz) - glze(te-1,ilz))
          
       end if
    end do

    time = user_time

    do ilz=lz_lo%llim_world, lz_lo%ulim_world
       ig = ig_idx(lz_lo, ilz)
       ik = ik_idx(lz_lo, ilz)
       it = it_idx(lz_lo, ilz)
       ie = ie_idx(lz_lo, ilz)
       is = is_idx(lz_lo, ilz)
       je = jend(ig)

       if (je == 0) then
          te = 2*ng2
       else
          te = 2*je-1
       end if

       if (idx_local (lz_lo, ilz)) then
          if (proc0) then 
             ltmp = lcoll(:,ilz)
             ftmp = fdcoll(:,ilz)
          else
             call send (lcoll(:,ilz), 0)
             call send (fdcoll(:,ilz), 0)
          endif
       else if (proc0) then
          call receive (ltmp, proc_id(lz_lo, ilz))
          call receive (ftmp, proc_id(lz_lo, ilz))
       endif
       call barrier

       do il=1,te
          etmp = cabs(ltmp(il) - ftmp(il))
          
          if (etmp > emax) then
             emax = etmp
             ltmax = cabs(ltmp(il))
             ikmax = ik
             itmax = it
             iemax = ie
             ismax = is
             ilmax = il
             igmax = ig
          end if
          
          etot = etot + etmp
          edenom = edenom + cabs(ltmp(il))
       end do
    end do

    eavg = etot/edenom
    emax = emax/ltmax

    if (proc0 .and. gnostics%write_ascii) then
       write(unit,"((1x,e13.6),6(i8),2(1x,e13.6))") time, &
            igmax, ikmax, itmax, iemax, ilmax, ismax, emax, eavg
    end if

    deallocate (lcoll, fdcoll, glze, ltmp, ftmp, g0)
    
  end subroutine write_collision_error
end module diagnostics_velocity_space

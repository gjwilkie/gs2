module lowflow
  
  implicit none

  public :: get_lowflow_terms, dphidth
  public :: phineo
  public :: mach_lab    

  real, dimension (:), allocatable :: dphidth
  real, dimension (:,:,:,:,:), allocatable :: coefs
  real, dimension (:,:), allocatable :: phineo, dcoefsdr
  real, dimension (:,:,:), allocatable :: dcoefsdth
  real, dimension (:), allocatable :: rad_neo, theta_neo
  real, dimension (:), allocatable :: mach_lab

contains
  
  subroutine get_lowflow_terms (theta, al, energy, bmag, dHdEc, dHdxic, vpadHdEc, dHdrc, &
       dHdthc, hneoc, dphidrc, dphidthc, phi_neo, lf_default, lf_decompose)
    
    use mp, only: proc0, broadcast
    use le_grids, only: w, wl
    use theta_grid, only: delthet, jacob, ntgrid, Rplot
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use geometry, only: rmaj
    use species, only: spec

    implicit none
    
    real, dimension (:), intent (in) :: theta, bmag
    real, dimension (:), intent (in) ::  al
    real, dimension (:), intent (in) :: energy
    real, dimension (:,:,:,:,:), intent (out) :: dHdec, dHdxic, dHdrc, dHdthc, vpadHdEc, hneoc
    real, dimension (:), intent (out) :: dphidthc, dphidrc, phi_neo
    logical, intent (in) :: lf_default, lf_decompose

    real :: rgeo2, drlt, drln

    real, dimension (:,:,:,:,:), allocatable :: hneo
    real, dimension (:,:,:,:), allocatable :: dHdxi, dHdE, vpadHdE, dHdr, dHdth
    real, dimension (:,:,:), allocatable :: legp
    real, dimension (:,:), allocatable :: xi, emax, dphidr
    real, dimension (:,:,:,:), allocatable :: chebyp1, chebyp2
    real, dimension (:), allocatable :: dl_over_b, drad_neo
    real, dimension (:,:), allocatable :: transport 
    real, dimension (:), allocatable :: phineo2
    real, dimension (:,:), allocatable :: upar2, qpar2
    real, dimension (:,:,:), allocatable :: upar, qpar, afac, bfac
    real, dimension (:,:), allocatable :: pflx, qflx, vflx, qparB, upar1, uparB, upar_over_B
    real, dimension (:), allocatable ::  radius, jboot, phi2, vtor, upar0
    real, dimension (size(energy)) :: w_r 

    integer :: il, ie, is, ns, nc, nl, nr, ig, ixi, ir, ir_loc
    integer :: ntheta, nlambda, nenergy, nxi
    integer, save :: neo_unit, neot_unit

!    logical, save :: initialized = .false.
    logical, dimension (:,:), allocatable :: forbid
    
!    if (initialized) return
!    initialized = .true.

    ntheta = size(theta)
    nlambda = size(al)
    nenergy = size(energy)
    nxi = 2*nlambda-1

    if (.not. allocated(xi)) allocate (xi(ntheta,nxi))
    if (.not. allocated(forbid)) allocate(forbid(ntheta,nxi))
    forbid = .false. ; xi = 0.0

    ! get xi=vpar/v from lambda=(mu/E)*Bnorm
    do il = 1, nlambda-1
       do ig = 1, ntheta
          if (1.0-al(il)*bmag(ig) < 0.0) then
             forbid(ig,il) = .true.
             forbid(ig,2*nlambda-il) = .true.
          else
             xi(ig,il) = -sqrt(1.0-al(il)*bmag(ig))
             xi(ig,2*nlambda-il) = -xi(ig,il)
          end if
       end do
    end do

    ! ns is nspec, nc is nenergy, nl is nlambda,
    ! nr is number of radii
    ! taken from neo
    call read_neocoefs (theta, ns, nc, nl, nr, ir_loc)
    
    if (.not. allocated(chebyp1)) allocate (chebyp1(nr,nenergy,0:nc-1,ns), chebyp2(nr,nenergy,0:nc-1,ns))
    if (.not. allocated(legp)) allocate (legp(ntheta,nxi,0:nl+1))
    if (.not. allocated(hneo)) allocate (hneo(nr,ntheta,nxi,nenergy,ns))
    if (.not. allocated(dHdr)) allocate (   dHdr(ntheta,nxi,nenergy,ns))
    if (.not. allocated(dHdth)) allocate (  dHdth(ntheta,nxi,nenergy,ns))
    if (.not. allocated(dHdxi)) allocate (  dHdxi(ntheta,nxi,nenergy,ns))
    if (.not. allocated(dHdE)) allocate (   dHdE(ntheta,nxi,nenergy,ns))
    if (.not. allocated(vpadHdE)) allocate (vpadHdE(ntheta,nxi,nenergy,ns))
    if (.not. allocated (dphidr)) allocate (dphidr(ntheta,ns))
    if (.not. allocated (dl_over_b)) allocate (dl_over_b(ntheta))
    if (.not. allocated (emax)) allocate (emax(nr,ns))

    ! better to be taken from neo
    if (proc0) write (*,*) '# make sure ENERGY_MAX=16.0 in NEO INPUT file'
    ! emax is ENERGY_MAX (from NEO input file) times v_{ts}^2

    ! two ways to deal with radial variation of temperature in 
    ! NEO energy variable: 1st (lf_default) is to take it into
    ! account when constructing distribution function, i.e. H(EGS2,r);
    ! 2nd is to construct H(ENEO) and deal with it later by
    ! including a temp gradient term in wstarfac in dist_fn.fpp
    if (lf_default) then
       ! v_{ts} is a function of radius, so we need to convert emax
       ! to its equivalent value using v_{ts} from the center radius.
       ! this is necessary because we will be taking radial derivatives of 
       ! the distribution function with v fixed, not v/v_t(r) fixed.
       do is = 1, ns
          emax(2,is) = 16.0 ! this is EMAX for center grid point
          emax(1,is) = emax(2,is)*(1.0-spec(is)%tprim*(rad_neo(1)-rad_neo(2)))
          emax(3,is) = emax(2,is)*(1.0-spec(is)%tprim*(rad_neo(3)-rad_neo(2)))
       end do
    else
       emax = 16.0
    end if

    allocate (drad_neo(nr)) ; drad_neo = 0.
    drad_neo = rad_neo - rad_neo(2)

    ! get legendre polynomials on gs2 pitch-angle grid
    legp = 0.0
    do ixi = 1, nxi
       do ig = 1, ntheta
          if (.not. forbid(ig,ixi)) call legendre (xi(ig,ixi), legp(ig,ixi,:))
       end do
    end do

    ! get chebyshev polynomials of the first and second kinds on the gs2 energy grid
    ! note the argument of the chebyshev polynomials is z = 2 * sqrt(energy/emax) - 1
    do is = 1, ns
       do ie = 1, nenergy
          do ir = 1, nr
             call chebyshev (zfnc(energy(ie),emax(ir,is)), chebyp1(ir,ie,:,is), 1)
             call chebyshev (zfnc(energy(ie),emax(ir,is)), chebyp2(ir,ie,:,is), 2)
          end do
       end do
    end do

    ! first, get hneo = F_1/F_0 = H_1/F_0 - Z_s * e * Phi / T_s
    ! H_1/F_0 is constructed from Legendre and Chebyshev polynomials
    do is = 1, ns
       do ig = 1, ntheta
          do ir = 1, nr
             call get_H (coefs(ir,ig,:,:,is), legp(ig,:,:), chebyp1(ir,:,:,is), &
                  hneo(ir,ig,:,:,is), phineo(ir,ig), ig, is)
          end do
       end do
    end do

    ! calculate and write to file some useful quantities like Upar and Qpar
    dl_over_b = delthet*jacob
    dl_over_b = dl_over_b/sum(dl_over_b)

    if (.not. allocated(pflx)) allocate (pflx(nr,ns), qflx(nr,ns), vflx(nr,ns), upar1(nr,ns))  
    if (.not. allocated(upar)) then
       allocate (phineo2(nr))
       allocate (upar2(nr,ns), qpar2(nr,ns))
       allocate (uparB(nr,ns), qparB(nr,ns))
       allocate (upar(nr,ns,ntheta), qpar(nr,ns,ntheta))
       allocate (upar_over_b(nr,ns))
       upar = 0. ; qpar = 0.
       upar2 = 0. ; qpar2 = 0. ; phineo2 = 0.
       uparB = 0. ; qparB = 0.
       upar_over_B = 0.
    end if
    if (.not. allocated(transport)) allocate (transport(nr,(5+ns*8)))  !JPL, NEO output for all radius
    if (.not. allocated(radius)) allocate (radius(nr), phi2(nr), jboot(nr), vtor(nr), upar0(nr))
    if (.not. allocated(mach_lab)) allocate (mach_lab(ns))  

    if (proc0) then
       call get_unused_unit (neo_unit)
       ! read the diagnosis in NEO from "neo_transport.out"
       open (unit=neo_unit, file='neo_transport.out', status='old', action='read')
       do ir=1, nr
          read (neo_unit,*) transport(ir,:)
          radius(ir) = transport(ir,1)
          phi2(ir) = transport(ir,2)
          jboot(ir) = transport(ir,3)
          vtor(ir) = transport(ir,4)
          upar0(ir) = transport(ir,5)
          do is = 1, ns
             pflx(ir,is) = transport(ir,5+(is-1)*8+1)/sqrt(2.0)  ! NEO uses vt=sqrt(T/m) while GS2 uses vt=sqrt(2T/m)
             qflx(ir,is) = transport(ir,5+(is-1)*8+2)/(sqrt(2.0)**3)
             vflx(ir,is) = transport(ir,5+(is-1)*8+3)/(sqrt(2.0)**2)
             upar1(ir,is) = transport(ir,5+(is-1)*8+4)/sqrt(2.0)
          end do
       end do
       close (neo_unit)
       
       ! evaluate the neoclassical parallel heat flux and parallel velcotiy in GS2 geometry
       do ir=1, nr 
          do is=1, ns
             drlt = drad_neo(ir)*spec(is)%tprim ; drln = drad_neo(ir)*spec(is)%fprim
             do ie = 1, nenergy
                w_r(ie) = w(ie) * exp(energy(ie)*(1.0-1.0/(1.0-drlt))) / (1.0-drlt)**1.5
                do ixi = 1, nxi
                   il = min(ixi, nxi+1-ixi)
                   do ig = 1, ntheta
                      upar(ir,is,ig) = upar(ir,is,ig) + hneo(ir,ig,ixi,ie,is) &
                           * sqrt(energy(ie))*xi(ig,ixi) * w_r(ie) * wl(-ntgrid+ig-1,il)
                      qpar(ir,is,ig) = qpar(ir,is,ig) + hneo(ir,ig,ixi,ie,is) &
                           * sqrt(energy(ie))*xi(ig,ixi)*(energy(ie)-2.5*(1.0-drlt)) &
                           * (1.0-drln) * w_r(ie) * wl(-ntgrid+ig-1,il)
                   end do
                end do
             end do
          end do
       end do

       do is = 1, ns
          do ir = 1, nr
             uparB(ir,is) = sum(upar(ir,is,:)*dl_over_b*bmag)
             qparB(ir,is) = sum(qpar(ir,is,:)*dl_over_b*bmag)
             upar_over_B(ir,is) = sum(upar(ir,is,:)*dl_over_b/bmag)
             upar2(ir,is) = sum(upar(ir,is,:)**2*dl_over_b)
             qpar2(ir,is) = sum(qpar(ir,is,:)**2*dl_over_b)
             phineo2(ir) = sum(phineo(ir,:)**2*dl_over_b)
          end do
       end do

       rgeo2 = 0.
       do ig = 1, ntheta
          rgeo2 = rgeo2 + dl_over_b(ig)*Rplot(-ntgrid+ig-1)**2
       end do

       ! note that for total toroidal angular momentum (including diamagnetic and ExB) to be zero,
       ! upar_over_b/rgeo2 = -(omega_phi * a / v_{t,ref}) * (a/R0)

       call open_output_file (neot_unit,".neotransp")
      ! print out in ".neotransp" file
       write (neot_unit,fmt='(a14,a8,12a14)') "# 1) rad", "2) spec", "3) pflx", "4) qflx", "5) vflx", &
            "6) qparflx", "7) uparB(GS2)", "8) upar1(NEO)", "9) <phi**2>", "10) bootstrap", "11) upar_o_B", &
            "12) upar2", "13) qpar2", "14) phi2"
       do ir=1, nr
          do is = 1, ns
             write (neot_unit,fmt='(e14.5,i8,12e14.5)') radius(ir), is, pflx(ir,is), qflx(ir,is), vflx(ir,is), &
                  qparB(ir,is),uparB(ir,is), upar1(ir,is)*sqrt(1.0-drlt), &
                  phi2(ir), jboot(ir), rmaj*upar_over_B(ir,is)/rgeo2, upar2(ir,is), qpar2(ir,is), phineo2(ir)
             if (ir==2)  mach_lab(is)=rmaj*upar_over_B(ir,is)/rgeo2 
          end do
       end do
       call close_output_file (neot_unit)
    end if

    call broadcast(mach_lab)

    do ir = 1, nr
       do is = 1, ns
          call broadcast (upar(ir,is,:))
          call broadcast (qpar(ir,is,:))
       end do
    end do

    allocate (afac(nr,ns,ntheta), bfac(nr,ns,ntheta))

    ! Here only to set up F_1/F_0 = xi*(v/vth)*(a + b*(v/vth)^2)
    if (lf_decompose) then

       ! for testing
!       upar = 0.
!       qpar(1,:,:) = qpar(2,:,:) ; qpar(3,:,:) = qpar(2,:,:)
!       qpar = 0.
       ! added after doing the qpar=0 test (to decouple upar and dupar/dr effects)
!       upar(1,:,:) = upar(2,:,:) ; upar(3,:,:) = upar(2,:,:)
       ! added after doing the qpar=0,du/dr=0 test...am eliminating the du/dtheta effect
!       do is = 1, ns
!          do ir = 1, nr
!             upar(ir,is,:) = sum(upar(ir,is,:)*dl_over_b)
!          end do
!       end do

       do is = 1, ns
          do ir = 1, nr
             drlt = drad_neo(ir)*spec(is)%tprim ; drln = drad_neo(ir)*spec(is)%fprim
             afac(ir,is,:) = 2.*(upar(ir,is,:)-qpar(ir,is,:)/(1.0-drln) - drlt*upar(ir,is,:))/(1.0-drlt)**2
             bfac(ir,is,:) = 0.8*qpar(ir,is,:) / ((1.0-drln) * (1.0-drlt)**3)
!             afac(ir,is,:) = 2.*(upar(ir,is,:)-qpar(ir,is,:)/(1.0-drln)/(1.0-drlt))/sqrt(1.0-drlt)
!             bfac(ir,is,:) = 0.8*qpar(ir,is,:) / ((1.0-drln) * (1.0-drlt)**1.5)
          end do
       end do

       coefs = 0.0
       do ig = 1, ntheta
          coefs(:,ig,1,0,:) = sqrt(emax)*(0.5*afac(:,:,ig) + 0.3125*bfac(:,:,ig)*emax)
          coefs(:,ig,1,1,:) = 2.*sqrt(emax)*(0.25*afac(:,:,ig) + 0.234375*bfac(:,:,ig)*emax)
          coefs(:,ig,1,2,:) = 2.*0.09375*emax**1.5*bfac(:,:,ig)
          coefs(:,ig,1,3,:) = 2.*0.015625*emax**1.5*bfac(:,:,ig)
       end do

       ! first, get hneo = F_1/F_0 = H_1/F_0 - Z_s * e * Phi / T_s
       ! H_1/F_0 is constructed from Legendre and Chebyshev polynomials
       do is = 1, ns
          do ig = 1, ntheta
             do ir = 1, nr
                call get_H (coefs(ir,ig,:,:,is), legp(ig,:,:), chebyp1(ir,:,:,is), &
                     hneo(ir,ig,:,:,is), phineo(ir,ig), ig, is)
             end do
          end do
       end do
    
       upar = 0. ; qpar = 0.
       do ir = 1, nr
          do is = 1, ns
             drlt = drad_neo(ir)*spec(is)%tprim ; drln = drad_neo(ir)*spec(is)%fprim
             do ie = 1, nenergy
                w_r(ie) = w(ie) * exp(energy(ie)*(1.0-1.0/(1.0-drlt))) / (1.0-drlt)**1.5
                do ixi = 1, nxi
                   il = min(ixi, nxi+1-ixi)
                   do ig = 1, ntheta
                      upar(ir,is,ig) = upar(ir,is,ig) + hneo(ir,ig,ixi,ie,is) &
                           * sqrt(energy(ie))*xi(ig,ixi) * w_r(ie) * wl(-ntgrid+ig-1,il)
                      qpar(ir,is,ig) = qpar(ir,is,ig) + hneo(ir,ig,ixi,ie,is) &
                           * sqrt(energy(ie))*xi(ig,ixi)*(energy(ie)-2.5*(1.0-drlt)) &
                           * (1.0-drln) * w_r(ie) * wl(-ntgrid+ig-1,il)
                   end do
                end do
             end do
          end do
       end do
       
       do is = 1, ns
          do ir = 1, nr
             upar2(ir,is) = sum(upar(ir,is,:)**2*dl_over_b)
             qpar2(ir,is) = sum(qpar(ir,is,:)**2*dl_over_b)
          end do
       end do

    end if

    ! get dphi/dr at GS2 radius (center of 3 NEO radii)
    ! note that radgrad takes derivative of e*phi/Tref,
    ! not the derivative of Z_s*e*phi/T_s
    do ig = 1, ntheta
       ! Note that we are neglecting the term prop. to dTref/dr, which is fine
       ! as long as Tref in NEO is the same for all radii.
       ! This will require a changing TEMP with radius if the species temperature
       ! changes with radius.
       call get_radgrad (phineo(:,ig), rad_neo, ir_loc, dphidr(ig,1))
    end do

    ! get radial derivative of F_1/F_0 at fixed (xi,E)
    ! note that F_1/F_0 = H_1/F_0 - Z_s * e * Phi / T_s
!     do is = 1, ns
!        do ig = 1, ntheta
!           do ie = 0, nc-1
!              do ixi = 0, nl
!                 ! get radial derivative of spectral coefficients of H_1/H_0
!                 call get_radgrad (coefs(:,ig,ixi,ie,is), rad_neo, ir_loc, dcoefsdr(ixi,ie))
!              end do
!           end do
!           call get_gradH (dcoefsdr, dphidr(ig,is), legp(ig,:,:), chebyp1(ir_loc,:,:,is), dHdr(ig,:,:,is))
!        end do
!     end do

    ! get theta derivative of F_1/F_0 at fixed (xi,E)
    do is = 1, ns
       do ie = 0, nc-1
          do ixi = 0, nl
             ! get theta derivative of spectral coefficients of F_1/F_0
             call get_thgrad (coefs(ir_loc,:,ixi,ie,is), theta, dcoefsdth(:,ixi,ie))
          end do
       end do
       ! note that dphidth is calculated in read_neocoefs
       do ig = 1, ntheta
          call get_gradH (dcoefsdth(ig,:,:), dphidth(ig)*spec(is)%zt, legp(ig,:,:), chebyp1(ir_loc,:,:,is), dHdth(ig,:,:,is))
       end do
    end do

    do is = 1, ns
       do ig = 1, ntheta
          ! get dH/dxi at fixed E and dH/dE at fixed xi
          call get_dHdxi (coefs(ir_loc,ig,:,:,is), legp(ig,:,:), chebyp1(ir_loc,:,:,is), xi(ig,:), dHdxi(ig,:,:,is))
          call get_dHdE (coefs(ir_loc,ig,:,:,is), legp(ig,:,:), chebyp1(ir_loc,:,:,is), chebyp2(ir_loc,:,:,is), &
               energy(:), emax(ir_loc,is), dHdE(ig,:,:,is))
       end do
    end do

    ! get radial derivative of F_1/F_0 at fixed (xi,E)
    ! note that F_1/F_0 = H_1/F_0 - Z_s * e * Phi / T_s
    do is = 1, ns
       do ig = 1, ntheta
!          do ie = 0, nc-1
!             do ixi = 0, nl
!                ! get radial derivative of spectral coefficients of H_1/H_0
!                call get_radgrad (coefs(:,ig,ixi,ie,is), rad_neo, ir_loc, dcoefsdr(ixi,ie))
          do ie = 1, nenergy
             do ixi = 1, nxi
                call get_radgrad (hneo(:,ig,ixi,ie,is), rad_neo, ir_loc, dHdr(ig,ixi,ie,is))
             end do
          end do
!          call get_gradH (dcoefsdr, dphidr(ig,is), legp(ig,:,:), chebyp1(ir_loc,:,:,is), dHdr(ig,:,:,is))
       end do
    end do

    ! get dH/dtheta and dH/dr
!     do is = 1, ns
!        do ie = 1, nenergy
!           do ixi = 1, nxi
!              do ig = 1, ntheta
!                 call get_dHdr (hneo(:,ig,ixi,ie,is), rad_neo, ir_loc, dHdr(ig,ixi,ie,is))
!              end do
!              call get_dHdth (hneo(ir_loc,:,ixi,ie,is), theta, dHdth(:,ixi,ie,is))
!           end do
!        end do
!     end do

! deprecated ->    ! this is actually dphi/dr - phi*d(ln T)/dr, where phi = e*Phi/Tref
!    dphidrc(1:ntheta-1) = 0.5*(dphidr(1:ntheta-1,1) + dphidr(2:ntheta,1))/spec(1)%zt
    dphidrc(1:ntheta-1) = 0.5*(dphidr(1:ntheta-1,1) + dphidr(2:ntheta,1))
    dphidrc(ntheta) = dphidrc(1)
    dphidthc(1:ntheta-1) = 0.5*(dphidth(1:ntheta-1) + dphidth(2:ntheta))
    dphidthc(ntheta) = dphidth(1)

    phi_neo = phineo(ir_loc,:)

    ! vpadHdE is the derivative of F1/F0 with respect to energy at fixed mu (not xi)
    do ie = 1, nenergy
       do ixi = 1, nxi
          do ig = 1, ntheta
             vpadHdE(ig,ixi,ie,:) = sqrt(energy(ie))*xi(ig,ixi)*dHdE(ig,ixi,ie,:) &
                  + (1.-xi(ig,ixi)**2)*dHdxi(ig,ixi,ie,:)/(2.*sqrt(energy(ie)))
          end do
       end do
    end do

    do is = 1, ns
       do ie = 1, nenergy
          do il = 1, nlambda
             ! should be weighted with backdif for consistency -- MAB
             dHdrc(1:ntheta-1,il,ie,1,is) = 0.5*(dHdr(1:ntheta-1,2*nlambda-il,ie,is) &
                  + dHdr(2:ntheta,2*nlambda-il,ie,is))
             dHdrc(1:ntheta-1,il,ie,2,is) = 0.5*(dHdr(1:ntheta-1,il,ie,is) &
                  + dHdr(2:ntheta,il,ie,is))
             dHdthc(1:ntheta-1,il,ie,1,is) = 0.5*(dHdth(1:ntheta-1,2*nlambda-il,ie,is) &
                  + dHdth(2:ntheta,2*nlambda-il,ie,is))
             dHdthc(1:ntheta-1,il,ie,2,is) = 0.5*(dHdth(1:ntheta-1,il,ie,is) &
                  + dHdth(2:ntheta,il,ie,is))
             dHdEc(1:ntheta-1,il,ie,1,is) = 0.5*(dHdE(1:ntheta-1,2*nlambda-il,ie,is) &
                  + dHdE(2:ntheta,2*nlambda-il,ie,is))
             dHdEc(1:ntheta-1,il,ie,2,is) = 0.5*(dHdE(1:ntheta-1,il,ie,is) &
                  + dHdE(2:ntheta,il,ie,is))
             dHdxic(1:ntheta-1,il,ie,1,is) = 0.5*(dHdxi(1:ntheta-1,2*nlambda-il,ie,is) &
                  + dHdxi(2:ntheta,2*nlambda-il,ie,is))
             dHdxic(1:ntheta-1,il,ie,2,is) = 0.5*(dHdxi(1:ntheta-1,il,ie,is) &
                  + dHdxi(2:ntheta,il,ie,is))
             vpadHdEc(1:ntheta-1,il,ie,1,is) = 0.5*(vpadHdE(1:ntheta-1,2*nlambda-il,ie,is) &
                  + vpadHdE(2:ntheta,2*nlambda-il,ie,is))
             vpadHdEc(1:ntheta-1,il,ie,2,is) = 0.5*(vpadHdE(1:ntheta-1,il,ie,is) &
                  + vpadHdE(2:ntheta,il,ie,is))
             hneoc(1:ntheta-1,il,ie,1,is) = 0.5*(hneo(ir_loc,1:ntheta-1,2*nlambda-il,ie,is) &
                  + hneo(ir_loc,2:ntheta,2*nlambda-il,ie,is))
             hneoc(1:ntheta-1,il,ie,2,is) = 0.5*(hneo(ir_loc,1:ntheta-1,il,ie,is) &
                  + hneo(ir_loc,2:ntheta,il,ie,is))
          end do
       end do
    end do
    dHdrc(ntheta,:,:,:,:) = 0.0 ; dHdthc(ntheta,:,:,:,:) = 0.0 ; vpadHdEc(ntheta,:,:,:,:) = 0.0
    dHdEc(ntheta,:,:,:,:) = 0.0 ; dHdxic(ntheta,:,:,:,:) = 0.0 ; hneoc(ntheta,:,:,:,:) = 0.0
    dphidrc(ntheta) = 0.0 ; dphidthc(ntheta) = 0.0

    deallocate (xi, emax, chebyp1, chebyp2, legp, dHdr, dHdth, dHdxi, dHdE, vpadHdE, hneo, forbid)
    deallocate (coefs, drad_neo)
    deallocate (dphidr, dl_over_b, transport)
    deallocate (pflx, qflx, vflx, upar1, qparB)
    deallocate (upar, qpar, uparB)
    deallocate (upar2, qpar2)
    deallocate (afac, bfac)

  end subroutine get_lowflow_terms
  
  function zfnc (enrgy, enrgymax)
    
    implicit none
    
    real :: enrgy, enrgymax
    real :: zfnc
    
    zfnc = 2.*sqrt(enrgy/enrgymax)-1.
    
    return
    
  end function zfnc
  
  ! knd = 1 (2) for cheb polys of first (second) kind
  subroutine chebyshev (x, chebyp, knd)
    
    implicit none
    
    integer, intent (in) :: knd
    real, intent (in) :: x
    real, dimension (0:), intent (out) :: chebyp
    
    integer :: n, idx
    
    n = size(chebyp)-1
    
    chebyp(0) = 1.0
    
    if (knd == 2) then
       chebyp(1) = 2.*x
    else
       chebyp(1) = x
    end if
    do idx = 2, n
       chebyp(idx) = 2.*x*chebyp(idx-1) - chebyp(idx-2)
    end do
    
  end subroutine chebyshev
  
  subroutine legendre (x, legp)
    
    implicit none
    
    real, intent (in) :: x
    real, dimension (0:), intent (out) :: legp
    
    integer :: n, idx
    
    n = size(legp)-1
    
    legp(0) = 1.0
    legp(1) = x
    
    do idx = 2, n
       legp(idx) = ((2.*idx-1.)*x*legp(idx-1) + (1.-idx)*legp(idx-2))/idx
    end do
    
  end subroutine legendre
  
  subroutine get_H (gjk, legdre, chebyshv, h, phi, ig, is)

    use species, only: spec

    implicit none

    real, dimension (0:,0:), intent (in) :: gjk
    real, dimension (:,0:), intent (in) :: legdre, chebyshv
    real, intent (in) :: phi
    real, dimension (:,:), intent (out) :: h
    integer, intent (in) :: ig, is

    integer :: ix, ij, ik

    h = 0.0

    do ix = 1, size(h,1)
       do ik = 0, size(gjk,2)-1
          do ij = 0, size(gjk,1)-1
             h(ix,:) = h(ix,:) + chebyshv(:,ik)*legdre(ix,ij)*gjk(ij,ik)
          end do
       end do
    end do

    ! F_1s = H_1s - Z_s*e*Phi_1/T_s
    ! want Z_s*e*Phi/T_s, but phi from neo is e*Phi/Tnorm
    ! at center (GS2) radius, this is a simple multiply by Z_s * Tref/Ts
    ! note that this assumes the Tref used in NEO is the same as that used in GS2
    h = h - spec(is)%zt*phi

  end subroutine get_H

  subroutine get_dHdxi (gjk, legdre, chebyshv, x, dH)
    
    implicit none
    
    real, dimension (0:,0:), intent (in) :: gjk
    real, dimension (:,0:), intent (in) :: legdre
    real, dimension (:), intent (in) :: x
    real, dimension (:,0:), intent (in) :: chebyshv
    real, dimension (:,:), intent (out) :: dH
    
    integer :: ij, ik, ix
    
    dH = 0.0
    
    ! calculate dH/dxi = sum_{j,k} g_{j,k} * T_{k}(z) * d(P_{j}(xi))/d(xi), where z=2*sqrt(E/EMAX)-1
    ! note that dP_{j}/dxi = (j+1)*(P_{j+1} - xi*P_{j})/(xi^2-1)
    do ix = 1, size(x)
       do ik = 0, size(gjk,2)-1
          do ij = 0, size(gjk,1)-1
             dH(ix,:) = dH(ix,:) + (chebyshv(:,ik)*(ij+1)*(legdre(ix,ij+1)-x(ix)*legdre(ix,ij))*gjk(ij,ik))/(x(ix)**2-1.)
          end do
       end do
    end do

  end subroutine get_dHdxi
  
  subroutine get_dHdE (gjk, legdre, chebyshv1, chebyshv2, x, xmax, dH)
    
    implicit none
    
    real, dimension (0:,0:), intent (in) :: gjk
    real, dimension (:,0:), intent (in) :: legdre
    real, dimension (:), intent (in) :: x
    real, intent (in) :: xmax
    real, dimension (:,0:), intent (in) :: chebyshv1, chebyshv2
    real, dimension (:,:), intent (out) :: dH
    
    integer :: ij, ik, ix
    real :: cfac
    
    dH = 0.0
    
    ! dH/dE at fixed xi is sum_{j,k} c_{j,k}*P_j(xi)*d(T_{n}(z))/dz / sqrt(E*EMAX)
    ! z=2*sqrt(E/EMAX)-1
    ! note that d(T_{k})/dz = k*U_{k-1}, where U_k is the Chebyshev polynomial of the 2nd kind
    do ix = 1, size(x)
       do ik = 0, size(gjk,2)-1
          if (ik==0) then
! cfac nonzero here if hneo = F_1 instead of F_1/F_0
!             cfac = -chebyshv1(ix,ik)
             cfac = 0.0
          else
! extra chebyshv1 term not needed since hneo = F_1/F_0 instead of F_1
!             cfac = ik/sqrt(xmax*x(ix))*chebyshv2(ix,ik-1)-chebyshv1(ix,ik)
             cfac = ik/sqrt(xmax*x(ix))*chebyshv2(ix,ik-1)
          end if
          do ij = 0, size(gjk,1)-1
             dH(:,ix) = dH(:,ix) + (legdre(:,ij)*cfac*gjk(ij,ik))
          end do
       end do
    end do
    
  end subroutine get_dHdE

  ! returns radial derivative at center of 3 points
  subroutine get_radgrad (h, rad, ir, dh)

    implicit none
    
    real, dimension (:), intent (in) :: h, rad
    integer, intent (in) :: ir
    real, intent (out) :: dh
    
    dh = (h(ir+1)-h(ir-1))/(rad(ir+1)-rad(ir-1))
    
  end subroutine get_radgrad

  subroutine get_thgrad (h, th, dh)

    implicit none

    real, dimension (:), intent (in) :: h, th
    real, dimension (:), intent (out) :: dh

    integer :: ig, nth

    nth = size(th)

    do ig = 2, nth-1
       dh(ig) = (h(ig+1)-h(ig-1))/(th(ig+1)-th(ig-1))
    end do
    ! note that H_neo is periodic in theta
    dh(1) = (h(2)-h(nth))/(2.*(th(2)-th(1)))
    dh(nth) = (h(1)-h(nth-1))/(2.*(th(nth)-th(nth-1)))

  end subroutine get_thgrad

  subroutine get_gradH (dgjk, dphi, legdre, chebyshv, dH)
    
    implicit none
    
    real, dimension (0:,0:), intent (in) :: dgjk
    real, intent (in) :: dphi
    real, dimension (:,0:), intent (in) :: legdre
    real, dimension (:,0:), intent (in) :: chebyshv
    real, dimension (:,:), intent (out) :: dH
    
    integer :: ij, ik, ix
    
    dH = 0.0
    
    ! calculate dH/drho or dH/dtheta at fixed (xi,E) as sum_{j,k} d(g_{j,k})/drho * T_{k}(z) * P_{j}(xi)
    ! similarly for dH/dtheta
    ! where z=2*sqrt(E/EMAX)-1
    do ix = 1, size(legdre,1)
       do ik = 0, size(dgjk,2)-1
          do ij = 0, size(dgjk,1)-1
             dH(ix,:) = dH(ix,:) + chebyshv(:,ik)*legdre(ix,ij)*dgjk(ij,ik)
          end do
       end do
    end do

    dH = dH - dphi

  end subroutine get_gradH

  subroutine read_neocoefs (theta, nspec_neo, nenergy_neo, nxi_neo, nrad_neo, ir_neo)

    use mp, only: proc0, broadcast
    use splines, only: lf_spline
    use file_utils, only: get_unused_unit
    use constants, only: pi,twopi

    implicit none

    real, dimension (:), intent (in) :: theta
    integer, intent (out) :: nspec_neo, nenergy_neo, nxi_neo, nrad_neo, ir_neo

    integer :: is, ik, ij, ig, ir, idx, ntheta, ntheta_neo
    integer :: ptheta_neo, ip 
    integer, save :: neo_unit, neof_unit, neophi_unit
   
    real, dimension (:), allocatable :: tmp, neo_coefs, dum, neo_phi, dneo_phi
    real, dimension (:), allocatable :: theta_neo_ext,neo_coefs_ext, neo_phi_ext, dneo_phi_ext !JPL
    ntheta = size(theta)

    if (proc0) then
       call get_unused_unit (neo_unit)

       ! read in number of grid points from neo's grid.out file
       open (unit=neo_unit, file='neo_grid.out', status="old", action="read")
       read (neo_unit,*) nspec_neo
       read (neo_unit,*) nenergy_neo
       read (neo_unit,*) nxi_neo
       read (neo_unit,*) ntheta_neo
       if (.not. allocated(theta_neo)) allocate (theta_neo(ntheta_neo))
       do ig = 1, ntheta_neo
          read (neo_unit,*) theta_neo (ig)
       end do
       read (neo_unit,*) nrad_neo
       if (.not. allocated(rad_neo)) allocate (rad_neo(nrad_neo))
       do ir = 1, nrad_neo
          read (neo_unit,*) rad_neo(ir)
       end do
       close (neo_unit)
    end if

    call broadcast (nspec_neo)
    call broadcast (nenergy_neo)
    call broadcast (nxi_neo)
    call broadcast (ntheta_neo)
    call broadcast (nrad_neo)
    if (.not. allocated(theta_neo)) allocate (theta_neo(ntheta_neo))
    if (.not. allocated(rad_neo)) allocate (rad_neo(nrad_neo))
    call broadcast (theta_neo)
    call broadcast (rad_neo)

    ! for now, set ir_neo by hand, but best to derive it from neo output in future
    ir_neo = 2

    if (.not. allocated(tmp)) allocate (tmp(ntheta_neo*(nxi_neo+1)*nenergy_neo*nspec_neo*nrad_neo))
    if (.not. allocated(neo_coefs)) allocate (neo_coefs(ntheta_neo), neo_phi(ntheta_neo), dneo_phi(ntheta_neo))
    if (.not. allocated(dum)) allocate (dum(ntheta))
    if (.not. allocated(coefs)) allocate (coefs(nrad_neo,ntheta,0:nxi_neo,0:nenergy_neo-1,nspec_neo))
    if (.not. allocated(phineo)) allocate (phineo(nrad_neo,ntheta))
    if (.not. allocated(dphidth)) allocate (dphidth(ntheta))
    if (.not. allocated(dcoefsdr)) then
       allocate (dcoefsdr(0:nxi_neo,0:nenergy_neo-1)) ; dcoefsdr = 0.
       allocate (dcoefsdth(ntheta,0:nxi_neo,0:nenergy_neo-1)) ; dcoefsdth = 0.
    end if

    if (proc0) then
       
       !JPL: if the range of theta grid in GS2 is beyond [-pi:pi](e.g. "nperiod > 1"),
       !     extend the theta grid in NEO ([-pi:pi]) to the periodic one in the theta range of GS2.
       ptheta_neo=ceiling((maxval(theta)+pi)/twopi)  !the period of 2pi in theta range
       if (ptheta_neo .gt. 1) then 
          if (.not. allocated(theta_neo_ext)) allocate (theta_neo_ext(ntheta_neo*(2*ptheta_neo-1)))
          if (.not. allocated(neo_coefs_ext)) allocate (neo_coefs_ext(ntheta_neo*(2*ptheta_neo-1)))
          if (.not. allocated(neo_phi_ext)) allocate (neo_phi_ext(ntheta_neo*(2*ptheta_neo-1)))
          if (.not. allocated(dneo_phi_ext)) allocate (dneo_phi_ext(ntheta_neo*(2*ptheta_neo-1)))
       end if

       ! read in H1^{nc} (adiabatic piece of F1^{nc}) from neo's f.out file
       call get_unused_unit (neof_unit)
       open (unit=neof_unit, file='neo_f.out', status="old", action="read")

       read (neof_unit,*) tmp

       idx = 1
       do ir = 1, nrad_neo
          do is = 1, nspec_neo
             do ik = 0, nenergy_neo-1
                do ij = 0, nxi_neo
                   do ig = 1, ntheta_neo
                      neo_coefs(ig) = tmp(idx)
                      idx = idx+1
                   end do

                   if ( ptheta_neo .gt. 1) then
                      do ip = 1, ntheta_neo*(2*ptheta_neo-1)
                         theta_neo_ext(ip)=theta_neo(mod((ip-1),ntheta_neo)+1)+twopi*(int((ip-1)/ntheta_neo)-(ptheta_neo-1))
                         neo_coefs_ext(ip)=neo_coefs(mod((ip-1),ntheta_neo)+1)
                      end do
                      call lf_spline (theta_neo_ext, neo_coefs_ext, theta, coefs(ir,:,ij,ik,is), dum)  !JPL
                   else
                   ! need to interpolate coefficients from neo's theta grid to gs2's
                      call lf_spline (theta_neo, neo_coefs, theta, coefs(ir,:,ij,ik,is), dum)
                   end if
                end do
             end do
          end do
       end do

       close (neof_unit)
    end if

    deallocate (tmp)

    if (.not. allocated(tmp)) allocate (tmp(ntheta_neo*nrad_neo))

    if (proc0) then
       ! read in phi1^{nc} from neo's phi.out file
       call get_unused_unit (neophi_unit)
       open (unit=neophi_unit, file='neo_phi.out', status="old", action="read")

       read (neophi_unit,*) tmp

       idx = 1
       do ir = 1, nrad_neo
          do ig = 1, ntheta_neo
             neo_phi(ig) = tmp(idx)
             idx = idx+1
          end do

          ! TMP FOR TESTING -- MAB
!          neo_phi = 0.

          ! need to interpolate coefficients from neo's theta grid to gs2's
          if ( ptheta_neo .gt. 1) then
             do ip = 1, ntheta_neo*(2*ptheta_neo-1)
                neo_phi_ext(ip)=neo_phi(mod((ip-1),ntheta_neo)+1)
             end do
             call lf_spline (theta_neo_ext, neo_phi_ext, theta, phineo(ir,:), dum)  !JPL
          else          
             call lf_spline (theta_neo, neo_phi, theta, phineo(ir,:), dum)
          end if

          ! at central radius, calculate dphi/dth and interpolate onto gs2 grid
          if (ir == 2) then
             if ( ptheta_neo .gt. 1) then  
                do ip = 1, ntheta_neo*(2*ptheta_neo-1)
                   dneo_phi_ext(ip)=dneo_phi(mod((ip-1),ntheta_neo)+1)
                end do
                call get_thgrad (neo_phi_ext, theta_neo_ext, dneo_phi_ext)
                call lf_spline (theta_neo_ext, dneo_phi_ext, theta, dphidth, dum) 
             else
                call get_thgrad (neo_phi, theta_neo, dneo_phi)
                call lf_spline (theta_neo, dneo_phi, theta, dphidth, dum)
             end if
          end if
       end do
       
       close (neophi_unit)

    end if

    call broadcast (dphidth)
    do ir = 1, nrad_neo
       call broadcast (phineo(ir,:))
    end do
    do ir = 1, nrad_neo
       do is = 1, nspec_neo
          do ik = 0, nenergy_neo-1
             do ij = 0, nxi_neo
                call broadcast (coefs(ir,:,ij,ik,is))
             end do
          end do
       end do
    end do

    deallocate (tmp, neo_coefs, theta_neo, neo_phi, dneo_phi, dum)
    if (proc0 .and. (ptheta_neo .gt. 1)) then 
       deallocate (neo_coefs_ext, theta_neo_ext, neo_phi_ext, dneo_phi_ext) 
    endif
  end subroutine read_neocoefs

end module lowflow

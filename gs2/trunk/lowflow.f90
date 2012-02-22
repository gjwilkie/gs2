module lowflow
  
  implicit none

  public :: get_lowflow_terms, dphidth

  real, dimension (:), allocatable :: dphidth
  real, dimension (:,:,:,:,:), allocatable :: coefs
  real, dimension (:,:), allocatable :: phineo
  real, dimension (:), allocatable :: rad_neo, theta_neo

contains
  
  subroutine get_lowflow_terms (theta, al, energy, bmag, dHdEc, dHdxic, vpadHdEc, dHdrc, &
       dHdthc, hneoc, dphidrc, dphidthc, phi_neo)
    
    use mp, only: proc0
    use le_grids, only: w, wl
    use theta_grid, only: delthet, jacob, ntgrid
    use file_utils, only: open_output_file, close_output_file, get_unused_unit

    implicit none
    
    real, dimension (:), intent (in) :: theta, bmag
    real, dimension (:), intent (in) ::  al
    real, dimension (:), intent (in) :: energy
    real, dimension (:,:,:,:,:), intent (out) :: dHdec, dHdxic, dHdrc, dHdthc, vpadHdEc, hneoc
    real, dimension (:), intent (out) :: dphidthc, dphidrc, phi_neo

    real, dimension (:,:,:,:,:), allocatable :: hneo
    real, dimension (:,:,:,:), allocatable :: dHdxi, dHdE, vpadHdE, dHdr, dHdth
    real, dimension (:,:,:), allocatable :: legp
    real, dimension (:,:), allocatable :: xi, chebyp1, chebyp2
    real, dimension (:), allocatable :: dphidr, dl_over_b, transport
    real, dimension (:), allocatable :: pflx, qflx, vflx, qpar, upar1
    real :: emax, radius, jboot, phi2, vtor, upar0

    integer :: il, ie, is, ns, nc, nl, nr, ig, ixi, ir, ir_loc
    integer :: ntheta, nlambda, nenergy, nxi
    integer, save :: neo_unit, neot_unit

    logical, save :: initialized = .false.
    logical, dimension (:,:), allocatable :: forbid
    
    if (initialized) return
    initialized = .true.

    ntheta = size(theta)
    nlambda = size(al)
    nenergy = size(energy)
    nxi = 2*nlambda-1

    allocate (xi(ntheta,nxi), forbid(ntheta,nxi))
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

    ! ns is nspec, nc is nenergy, nl is nlambda
    ! taken from neo
    call read_neocoefs (theta, ns, nc, nl, nr, ir_loc)

    ! TMP FOR TESTING -- MAB
!    coefs = 0.0
!    coefs(:,:,1,0,:) = 0.6 ; coefs(:,:,3,0,:) = 0.4
!    coefs(:,0,1,:) = 1.0
    
    ! better to be taken from neo
    if (proc0) write (*,*) '# make sure ENERGY_MAX=16.0 in NEO INPUT file'
    emax = 16.0
    
    allocate (chebyp1(nenergy,0:nc-1), chebyp2(nenergy,0:nc-1))
    allocate (legp(ntheta,nxi,0:nl+1))
    allocate (hneo(nr,ntheta,nxi,nenergy,ns))
    allocate (   dHdr(ntheta,nxi,nenergy,ns))
    allocate (  dHdth(ntheta,nxi,nenergy,ns))
    allocate (  dHdxi(ntheta,nxi,nenergy,ns))
    allocate (   dHdE(ntheta,nxi,nenergy,ns))
    allocate (vpadHdE(ntheta,nxi,nenergy,ns))
    allocate (dphidr(ntheta))
    allocate (dl_over_b(ntheta))

    legp = 0.0
    do ixi = 1, nxi
       do ig = 1, ntheta
          if (.not. forbid(ig,ixi)) call legendre (xi(ig,ixi), legp(ig,ixi,:))
       end do
    end do
    
    do ie = 1, nenergy
       call chebyshev (zfnc(energy(ie),emax), chebyp1(ie,:), 1)
       call chebyshev (zfnc(energy(ie),emax), chebyp2(ie,:), 2)
    end do
    
! BD:  Switched order of first two loops for efficiency.  MAB should double-check for correctness
    do is = 1, ns
       do ir = 1, nr
          do ig = 1, ntheta
             ! get_H returns hneo = F_1 / F_0
             call get_H (coefs(ir,ig,:,:,is), legp(ig,:,:), chebyp1, hneo(ir,ig,:,:,is), phineo(ir,ig))
             call get_dHdxi (coefs(ir_loc,ig,:,:,is), legp(ig,:,:), chebyp1, xi(ig,:), dHdxi(ig,:,:,is))
             call get_dHdE (coefs(ir_loc,ig,:,:,is), legp(ig,:,:), chebyp1, chebyp2, energy(:), emax, dHdE(ig,:,:,is))
!             call get_dHdE (hneo(ig,ig,:,:,is),energy(:),dHdE(ig,:,:,is))
          end do
       end do
    end do

    ! get dH/dtheta and dH/dr
    do is = 1, ns
       do ie = 1, nenergy
          do ixi = 1, nxi
             do ig = 1, ntheta
                call get_dHdr (hneo(:,ig,ixi,ie,is), rad_neo, ir_loc, dHdr(ig,ixi,ie,is))
             end do
             call get_dHdth (hneo(ir_loc,:,ixi,ie,is), theta, dHdth(:,ixi,ie,is))
          end do
       end do
    end do

    ! get dphi/dr
    do ig = 1, ntheta
       call get_dHdr (phineo(:,ig), rad_neo, ir_loc, dphidr(ig))
    end do

    dphidrc(1:ntheta-1) = 0.5*(dphidr(1:ntheta-1) + dphidr(2:ntheta))
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

    dl_over_b = delthet*jacob
    dl_over_b = dl_over_b/sum(dl_over_b)

    ! get parallel heat flux (int d3v vpa * v^2 * F1^{nc})                                                                                                                  
    allocate (pflx(ns), qflx(ns), vflx(ns), upar1(ns), qpar(ns))
    qpar = 0.
    do ie = 1, nenergy
       do ixi = 1, nxi
          il = min(ixi, nxi+1-ixi)
          do ig = 1, ntheta
             qpar = qpar + hneo(ir_loc,ig,ixi,ie,:)*energy(ie) &
                  *(xi(ig,ixi)*sqrt(energy(ie))-upar1)*w(ie)*wl(-ntgrid+ig-1,il)*dl_over_b(ig)
          end do
       end do
    end do

    allocate (transport(5+ns*8))

    if (proc0) then
       call get_unused_unit (neo_unit)
       open (unit=neo_unit, file='neo_transport.out', status='old', action='read')
       read (neo_unit,*) transport
       radius = transport(1)
       phi2 = transport(2)
       jboot = transport(3)
       vtor = transport(4)
       upar0 = transport(5)
       do is = 1, ns
          pflx(is) = transport(5+(is-1)*8+1)
          qflx(is) = transport(5+(is-1)*8+2)
          vflx(is) = transport(5+(is-1)*8+3)
          upar1(is) = transport(5+(is-1)*8+4)
       end do
       close (neo_unit)

       call open_output_file (neot_unit,".neotransp")
       write (neot_unit,fmt='(a110)') "# 1) rad,    2) spec, 3) pflx,    4) qflx,    5) vflx   , 6) qparflx, 7) upar1   , 8) <phi**2>, 9) bootstrap"
       do is = 1, ns
          write (neot_unit,fmt='(e14.5,i4,7e14.5)') radius, is, pflx(is), qflx(is), vflx(is), qpar(is), upar1(is), &
               phi2, jboot
       end do
       call close_output_file (neot_unit)
    end if

    deallocate (xi, chebyp1, chebyp2, legp, coefs, dHdr, dHdth, dHdxi, dHdE, vpadHdE, hneo, forbid)
    deallocate (dphidr, dl_over_b, transport)
    deallocate (pflx, qflx, vflx, upar1, qpar)

  end subroutine get_lowflow_terms
  
  function zfnc (enrgy, enrgymax)
    
    implicit none
    
    real :: enrgy, enrgymax
    real :: zfnc
    
    zfnc = 2.*sqrt(enrgy/enrgymax)-1.
    
    return
    
  end function zfnc
  
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
  
  subroutine get_H (gjk, legdre, chebyshv, h, phi)

    implicit none

    real, dimension (0:,0:), intent (in) :: gjk
    real, dimension (:,0:), intent (in) :: legdre, chebyshv
    real, intent (in) :: phi
    real, dimension (:,:), intent (out) :: h

    integer :: ix, ij, ik

    h = 0.0

    do ix = 1, size(h,1)
       do ik = 0, size(gjk,2)-1
          do ij = 0, size(gjk,1)-1
             h(ix,:) = h(ix,:) + chebyshv(:,ik)*legdre(ix,ij)*gjk(ij,ik)
          end do
       end do
    end do

    ! F_1 = H_1 - Phi_1
    h = h - phi

!     do ik = 0, size(gjk,2)-1
!        do ij = 0, size(gjk,1)-1
!           write (*,*) 'chebyshv', ij, ik, gjk(ij,ik)
!        end do
!     end do

!    do ij = 1, size(h,2)
!       do ix = 1, size(h,1)
!          write (*,*) 'hneo, ix, ie', ix, ij, h(ix,ij)
!       end do
!    end do

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

  subroutine get_dHdr (h, rad, ir, dh)

    implicit none

    real, dimension (:), intent (in) :: h, rad
    integer, intent (in) :: ir
    real, intent (out) :: dh

    dh = (h(ir+1)-h(ir-1))/(rad(ir+1)-rad(ir-1))

  end subroutine get_dHdr

  subroutine get_dHdth (h, th, dh)

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

  end subroutine get_dHdth

  subroutine read_neocoefs (theta, nspec_neo, nenergy_neo, nxi_neo, nrad_neo, ir_neo)

    use mp, only: proc0, broadcast
    use splines, only: lf_spline
    use file_utils, only: get_unused_unit

    implicit none

    real, dimension (:), intent (in) :: theta
    integer, intent (out) :: nspec_neo, nenergy_neo, nxi_neo, nrad_neo, ir_neo

    integer :: is, ik, ij, ig, ir, idx, ntheta, ntheta_neo
    integer, save :: neo_unit, neof_unit, neophi_unit
   
    real, dimension (:), allocatable :: tmp, neo_coefs, dum, neo_phi, dneo_phi

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

    allocate (tmp(ntheta_neo*(nxi_neo+1)*nenergy_neo*nspec_neo*nrad_neo))
    allocate (neo_coefs(ntheta_neo), neo_phi(ntheta_neo), dneo_phi(ntheta_neo))
    allocate (dum(ntheta))
    if (.not. allocated(coefs)) allocate (coefs(nrad_neo,ntheta,0:nxi_neo,0:nenergy_neo-1,nspec_neo))
    if (.not. allocated(phineo)) allocate (phineo(nrad_neo,ntheta))
    if (.not. allocated(dphidth)) allocate (dphidth(ntheta))

    if (proc0) then
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
                   ! need to interpolate coefficients from neo's theta grid to gs2's
                   call lf_spline (theta_neo, neo_coefs, theta, coefs(ir,:,ij,ik,is), dum)
                end do
             end do
          end do
       end do

       close (neof_unit)
    end if

    deallocate (tmp)

    allocate (tmp(ntheta_neo*nrad_neo))

    if (proc0) then
       ! read in phi1^{nc} from neo's phi.out file
       open (unit=neophi_unit, file='neo_phi.out', status="old", action="read")

       read (neophi_unit,*) tmp

       idx = 1
       do ir = 1, nrad_neo
          do ig = 1, ntheta_neo
             neo_phi(ig) = tmp(idx)
             idx = idx+1
          end do

          ! need to interpolate coefficients from neo's theta grid to gs2's
          call lf_spline (theta_neo, neo_phi, theta, phineo(ir,:), dum)

          ! at central radius, calculate dphi/dth and interpolate onto gs2 grid
          if (ir == 2) then
             call get_dHdth (neo_phi, theta_neo, dneo_phi)
             call lf_spline (theta_neo, dneo_phi, theta, dphidth, dum)
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

  end subroutine read_neocoefs

end module lowflow

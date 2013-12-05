module leq

  implicit none

  integer :: nr, nt
  private

  integer :: ntg

  real :: dpsidrho, dIdrho, bi, dqdr, dI2dr2
  real, dimension (:), allocatable :: d2Rdrdth, d2Zdrdth, grho, bmag, dBdrho, d2Bdrdth
  real, dimension (:), allocatable :: dgradpardrho, gradpar, dgradparBdrho, dBdth, gradparb
  real, dimension (:), allocatable :: cvdrift0, dcvdrift0drho, dgbdrift0drho, theta
  real, dimension (:), allocatable :: varthet, dvarthdr, gradrho_gradthet, cross, d2varthdr2
  real, dimension (:), allocatable :: gradthet2, gradalph_gradthet, gradrho_gradalph, gradalph2
  real, dimension (:), allocatable :: cvdrift, gbdrift, d2Bdr2, d2Rdr2, d2Zdr2, drz, drzdth
  real, dimension (:), allocatable :: d2Rdr2dth, d2Zdr2dth, d2gpsidr2, d2Idr2, dcrossdr
  real, dimension (:), allocatable :: dcvdriftdrho, dgbdriftdrho, dgds2dr, dgds21dr, dgds22dr
  real, dimension (:,:), allocatable :: Rr, Zr
  real, allocatable, dimension (:)     :: eqpsi, fp, beta, pressure
  real, allocatable, dimension (:,:)   :: R_psi, Z_psi
  real, allocatable, dimension (:,:,:) :: drm, dzm, dbtm, dpm, dtm
  real, allocatable, dimension (:,:,:) :: dpcart, dtcart, dbtcart
  real, allocatable, dimension (:,:,:) :: dpbish, dtbish, dbtbish

  real, dimension (:), allocatable :: jacrho, delthet, djacdrho, d2jacdr2, dRdrho, dZdrho, dRdth, dZdth

  real :: beta_0
  
  type :: flux_surface
     real :: R_center, R_geo, k, kp, d, dp, r, dr, delp, q, shat, pp, a, ap, &
          betaprim, betadbprim, d2qdr2
     integer :: nt
  end type flux_surface

  type (flux_surface) :: surf

  public :: leq_init, leqin, gradient, eqitem, bgradient, leqcoefs

  public :: invR, Rpos, Zpos, diameter, btori, dbtori,  qfun, pfun, &
       dpfun, betafun, psi, rcenter, dpdrhofun

contains

  subroutine leqin(R0, Ra, k, kp, d, dp, r, dr, s, qq, qs, a, ap, bp, bpp, d2q, nt_used)
        
    real :: R0, Ra, k, kp, d, dp, r, dr, s, qq, qs, a, ap, bp, bpp, d2q
    integer :: nt_used

    surf%R_center = R0
    surf%R_geo = Ra
    surf%delp = s
    surf%k = k
    surf%kp = kp
    surf%q = qq
    surf%shat = qs
    surf%d = d
    surf%dp = dp
    surf%a = a
    surf%ap = ap
    surf%r = r
    surf%dr = dr
    surf%pp = 0.
    surf%betaprim = bp
    surf%betadbprim = bpp
    surf%d2qdr2 = d2q

    beta_0 = 1.

    nr = 3
    nt = nt_used
    ntg = nt-1
    if(.not.allocated(beta)) call alloc_arrays(3, nt)
    surf%nt = nt
    dqdr = surf%shat*surf%q/surf%r
    call leq_init

  end subroutine leqin

  subroutine alloc_arrays(nr, nt)

    integer :: nr, nt

    allocate(eqpsi(nr), fp(nr), beta(nr), pressure(nr))
    allocate(R_psi(nr, nt), Z_psi(nr, nt))
    allocate(drm(nr, nt, 2), dzm(nr, nt, 2), dbtm(nr, nt, 2), &
         dpm(nr, nt, 2), dtm(nr, nt, 2))
    allocate(dpcart(nr, nt, 2), dtcart(nr, nt, 2), dbtcart(nr, nt, 2))
    allocate(dpbish(nr, nt, 2), dtbish(nr, nt, 2), dbtbish(nr, nt, 2))
    allocate (Rr(3,-ntg:ntg), Zr(3,-ntg:ntg))
    allocate (jacrho(-ntg:ntg), djacdrho(-ntg:ntg), d2jacdr2(-ntg:ntg))
    allocate (d2Rdrdth(-ntg:ntg), d2Zdrdth(-ntg:ntg))
    allocate (grho(-ntg:ntg), bmag(-ntg:ntg), dBdrho(-ntg:ntg), dgradpardrho(-ntg:ntg), gradpar(-ntg:ntg))
    allocate (d2Bdrdth(-ntg:ntg), dgradparBdrho(-ntg:ntg), dBdth(-ntg:ntg), gradparb(-ntg:ntg))
    allocate (cvdrift0(-ntg:ntg), dcvdrift0drho(-ntg:ntg), dgbdrift0drho(-ntg:ntg), theta(-ntg:ntg))
    allocate (dRdrho(-ntg:ntg), dZdrho(-ntg:ntg), dRdth(-ntg:ntg), dZdth(-ntg:ntg))
    allocate (varthet(-ntg:ntg), dvarthdr(-ntg:ntg), gradrho_gradthet(-ntg:ntg), cross(-ntg:ntg))
    allocate (gradthet2(-ntg:ntg), gradalph2(-ntg:ntg), gradalph_gradthet(-ntg:ntg))
    allocate (gradrho_gradalph(-ntg:ntg))
    allocate (cvdrift(-ntg:ntg), gbdrift(-ntg:ntg), d2Rdr2(-ntg:ntg), d2Zdr2(-ntg:ntg), d2Bdr2(-ntg:ntg))
    allocate (drz(-ntg:ntg), drzdth(-ntg:ntg), d2Rdr2dth(-ntg:ntg), d2Zdr2dth(-ntg:ntg))
    allocate (d2gpsidr2(-ntg:ntg), d2Idr2(-ntg:ntg), d2varthdr2(-ntg:ntg), dcrossdr(-ntg:ntg))
    allocate (dcvdriftdrho(-ntg:ntg), dgbdriftdrho(-ntg:ntg))
    allocate (dgds21dr(-ntg:ntg), dgds2dr(-ntg:ntg), dgds22dr(-ntg:ntg))

  end subroutine alloc_arrays

  subroutine leq_init

    implicit none
    real, dimension(nr, nt) :: eqpsi1, eqth, eqbtor
    
    real, dimension (-ntg:ntg) :: d2Rdth2, d2Zdth2

    real dr(3)
    real pi, t, r
    integer i, j

    pi=2*acos(0.)
    dr(1) = -surf%dr
    dr(2) = 0.
    dr(3) = surf%dr
    
    do j=-ntg,ntg
       theta(j) = j*pi/real(ntg)
       do i=1,3
          r = surf%r + dr(i)
          Rr(i,j) = Rpos(r,theta(j))
          Zr(i,j) = Zpos(r,theta(j))
       end do
    end do

    do j=1,nt
       do i=1,nr
          r = surf%r + dr(i)
          t = (j-1)*pi/real(nt-1)
          R_psi(i,j) = Rpos(r, t) 
          Z_psi(i,j) = Zpos(r, t)
          eqth(i,j) = t
          eqpsi1(i,j) = 1 + dr(i)
          eqbtor(i,j) = surf%r_geo/R_psi(i,j)
       enddo
    enddo

    do i=1,nr
       pressure(i) = -dr(i)
    enddo

    eqpsi(:) = eqpsi1(:,1)

    call derm(eqth,   dtm,  'T')
    call derm(R_psi,  drm,  'E')
    call derm(Z_psi,  dzm,  'O')
    call derm(eqbtor, dbtm, 'E')
    call derm(eqpsi1, dpm,  'E')

    allocate (delthet(-ntg:ntg-1))
    ! get delta theta as a function of theta
    delthet = theta(-ntg+1:)-theta(:ntg-1)

    ! get dR/drho and dZ/drho
    call get_drho (Rr, dRdrho)
    call get_drho (Zr, dZdrho)

    ! get dR/dtheta and dZ/dtheta
    call get_dthet(Rr(2,:), dRdth)
    call get_dthet(Zr(2,:), dZdth)

    ! I=Btor*R is a flux function
    ! bi = I/(Btor(psi,theta of Rgeo)*a) = Rgeo/a
    bi = surf%R_geo

    ! get second derivatives of R and Z with respect to theta
    call get_d2dthet2 (Rr(2,:), d2Rdth2)
    call get_d2dthet2 (Zr(2,:), d2Zdth2)
    ! get mixed theta and rho derivatives of R and Z
    call get_dthet (dRdrho, d2Rdrdth)
    call get_dthet (dZdrho, d2Zdrdth)

    ! get the Jacobian of the transformation from (rho,theta,zeta) to (R,Z,zeta)
    call get_jacrho

    ! get dpsinorm/drho
    call get_dpsidrho

    ! get grad rho
    call get_gradrho

    ! quantity needed in calculation of dI/drho and djacrho/drho
    drz = (dRdrho*dRdth + dZdrho*dZdth)/jacrho
    call get_dthet (drz, drzdth)

    ! get dI/drho
    call get_dIdrho

    ! get djacrho/drho
    call get_djacdrho

    ! get d2R/drho2 and d2Z/drho2
    call get_d2RZdr2

    ! get theta derivative of d2R/drho2 and d2Z/drho2
    call get_dthet (d2Rdr2, d2Rdr2dth)
    call get_dthet (d2Zdr2, d2Zdr2dth)

    ! calculate the magnitude of B (normalized by B(psi,theta corresponding to Rgeo))
    ! B/B0 = sqrt(I**2 + |grad psi|**2)/R
    bmag = sqrt(bi**2 + (grho*dpsidrho)**2)/Rr(2,:)

    ! get dB/dtheta
    call get_dthet (bmag, dbdth)
    
    ! calculate b . grad theta
    gradpar = dpsidrho/(bmag*jacrho)
    ! b . grad B
    gradparb = gradpar*dBdth

    ! get dB/drho and d2B/drho2
    call get_dBdrho

    ! d (b . grad theta) / drho
    dgradpardrho = -gradpar*(dBdrho/bmag + djacdrho/jacrho)

    ! get d/dtheta (dB/drho)
    call get_dthet (dBdrho, d2Bdrdth)
    
    dgradparBdrho = dgradpardrho*dBdrho + gradpar*d2Bdrdth

    ! this is 2 *(bhat/B x grad B / B) . (grad q) * dpsiN/drho / (bhat . grad B)
    ! same as usual GS2 definition
    cvdrift0 = -2.*bi*surf%shat*surf%q/(bmag**2*surf%r)

    ! this is the rho derivative of cvdrift0
    dcvdrift0drho = cvdrift0*(gradparb*dIdrho/bi + dgradparbdrho - 2.*gradparb*dBdrho/bmag) &
         - 2.*bi*gradparb**2*surf%d2qdr2/bmag**2
!    dcvdrift0drho = cvdrift0*(dIdrho/bi + dgradparbdrho/gradparb - 2.*dBdrho/bmag) &
!         - 2.*bi*gradparb*surf%d2qdr2/bmag**2

    ! this is the rho derivative of gbdrift0 = (bhat x gradB/B) . (grad q)
    dgbdrift0drho = cvdrift0*bmag*(gradparb*dIdrho/bi + dgradparbdrho - gradparb*dBdrho/bmag) &
         - 2.*bi*gradparb**2*surf%d2qdr2/bmag
!    dgbdrift0drho = cvdrift0*bmag*(dIdrho/bi + dgradparbdrho/gradparb - dBdrho/bmag) &
!         - 2.*bi*gradparb*surf%d2qdr2/bmag

    cvdrift0 = cvdrift0*gradparb

    ! obtain varthet = (I/(q*(dpsi/dr)) * int_0^theta dtheta' jacrho/R^2
    call get_varthet

    ! obtain dvarthet/drho
    call get_dvarthdr

    ! grad theta . grad theta
    gradthet2 = (Rr(2,:)/jacrho)**2*(dRdrho**2 + dZdrho**2)
    ! grad rho . grad theta
    gradrho_gradthet = -(Rr(2,:)/jacrho)**2*(dRdrho*dRdth+dZdrho*dZdth)
    ! grad alpha . grad theta
    gradalph_gradthet = -(varthet*dqdr + surf%q*dvarthdr)*gradrho_gradthet &
         - bi*jacrho/(dpsidrho*Rr(2,:)**2)*gradthet2
    ! grad rho . grad alpha
    gradrho_gradalph = -(varthet*dqdr + surf%q*dvarthdr)*grho**2 &
         - bi*jacrho/(dpsidrho*Rr(2,:)**2)*gradrho_gradthet

    ! grad alpha . grad alpha
    gradalph2 = (1./Rr(2,:)**2) + ((varthet*dqdr+surf%q*dvarthdr)*grho)**2 &
         + 2.*bi*jacrho*(varthet*dqdr+surf%q*dvarthdr)*gradrho_gradthet/(dpsidrho*Rr(2,:)**2) &
         + (bi*jacrho/(dpsidrho*Rr(2,:)**2))**2*gradthet2

    ! this is (grad alpha x B) . grad theta
    cross = dpsidrho*(gradrho_gradalph*gradalph_gradthet - gradalph2*gradrho_gradthet)

    ! this is bhat/B x (grad B) . grad alpha * 2 * dpsiN/drho    
    gbdrift = 2.0*(-dBdrho + cross*dBdth*dpsidrho/bmag**2)
    ! this is bhat/B x (bhat . grad bhat) . grad alpha * 2 * dpsiN/drho
    ! this is assuming betaprim = 4*pi*ptot/B0^2 * (-d ln ptot / drho)
    cvdrift = (gbdrift + 2.0*surf%betaprim/bmag)/bmag

    ! this is d/drho (bhat/B x (grad B) . grad alpha) * 2 * dpsiN/drho
    dgbdriftdrho = 2.0*(-d2Bdr2 + dpsidrho*(dcrossdr*dBdth+cross(d2Bdrdth-2.*dBdth))/bmag**2)
    dcvdriftdrho = (dgbdriftdrho - dBdrho/bmag)/bmag + 2.0*surf%betadbprim/bmag**2 &
         - 4.0*surf%betaprim/bmag**3

    ! get d^2I/drho^2 and d^2 Jac / dr^2
    call get_d2Idr2_d2jacdr2

    ! get d/dr [(grad alpha x B) . grad theta]
    call get_dcrossdr

!    do j=-ntg,ntg
!       write (*,'(a4,7e12.4)') 'leq', theta(j), gbdrift(j)/bmag(j), cvdrift(j), gradalph2(j)*dpsidrho**2, &
!            gradrho_gradalph(j)*dqdr*dpsidrho**2, (grho(j)*dqdr*dpsidrho)**2, dvarthdr(j)
!    end do
!    write (*,*)

! below is actually grad(rho) instead of grad(psi),
! and 'cartesian' refers to (R,Z) coordinates -- MAB
! grad(psi) in cartesian form 
    call eqdcart(dpm, dpcart)
! grad(psi) in Bishop form 
    call eqdbish(dpcart, dpbish)

! grad(BT) in cartesian form
    call eqdcart(dbtm, dbtcart)
! grad(BT) in Bishop form
    call eqdbish(dbtcart, dbtbish)

! grad(theta) in cartesian form
    call eqdcart(dtm, dtcart)
! grad(theta) in Bishop form
    call eqdbish(dtcart, dtbish)

  end subroutine leq_init

  ! takes in f(r), with r given at three radial locations
  ! and returns df = df/dr at the middle radius
  subroutine get_drho (f, df)

    implicit none

    real, dimension (:,-ntg:), intent (in) :: f
    real, dimension (-ntg:), intent (out) :: df

    df = 0.5*(f(3,:)-f(1,:))/surf%dr

  end subroutine get_drho

  ! given function f(theta), calculate second derivative
  ! of f with respect to theta
  ! second order accurate, with equal grid spacing assumed
  subroutine get_d2dthet2 (f, d2f)

    implicit none

    real, dimension (-ntg:), intent (in) :: f
    real, dimension (-ntg:), intent (out) :: d2f

    ! assuming equal grid spacing in theta here
    d2f(-ntg+1:ntg-1) = (f(:ntg-2)-2.*f(-ntg+1:ntg-1)+f(-ntg+2:))/delthet(-ntg+1:ntg-1)**2

    ! use periodicity at boundary
    d2f(-ntg) = (f(ntg-1)-2.*f(-ntg)+f(-ntg+1))/delthet(-ntg+1)**2
    d2f(ntg) = d2f(-ntg)

  end subroutine get_d2dthet2

  ! given function f(theta:-pi->pi), calculate theta derivative
  ! second order accurate, with equal grid spacing assumed
  ! assumes periodic in theta -- may need to change this in future
  subroutine get_dthet (f, df)

    implicit none

    real, dimension (-ntg:), intent (in) :: f
    real, dimension (-ntg:), intent (out) :: df

    ! assuming equal grid spacing in theta here
    df(-ntg+1:ntg-1) = (f(-ntg+2:)-f(:ntg-2))/(delthet(:ntg-2)+delthet(-ntg+1:))

    ! use periodicity at boundary
    df(-ntg) = (f(-ntg+1)-f(ntg-1))/(delthet(-ntg)+delthet(ntg-1))
    df(ntg) = df(-ntg)

  end subroutine get_dthet

  subroutine get_jacrho

    implicit none

    ! jacrho = R*(dR/drho * dZ/dtheta - dR/dtheta * dZ/drho)
    jacrho = Rr(2,:)*(dRdrho*dZdth - dRdth*dZdrho)

  end subroutine get_jacrho

  ! get dpsinorm/drho = (I/2*pi*q)*int_0^{2*pi} dthet jacrho/R**2
  subroutine get_dpsidrho

    use constants, only: pi

    implicit none
    
    ! theta_integrate returns integral from 0 -> 2*pi
    call theta_integrate (jacrho/Rr(2,:)**2, dpsidrho)

    ! integration done using trapezoidal rule
    dpsidrho = dpsidrho*bi/(2.*pi*surf%q)

  end subroutine get_dpsidrho

  subroutine get_gradrho

    implicit none

    grho = Rr(2,:)*sqrt(dRdth**2 + dZdth**2)/jacrho

  end subroutine get_gradrho

  subroutine get_dIdrho

    implicit none

    real :: num1, num2, denom
    real, dimension (:), allocatable :: dum

    allocate (dum(-ntg:ntg)) ; dum = 0.

    ! betaprim below is (4*pi*ptot/B0^2)*(-d ln ptot / drho)
    dum = (-2.*(dRdth*d2Rdrdth + dZdth*d2Zdrdth)/jacrho &
         - drzdth + surf%betaprim*jacrho/dpsidrho**2 ) / grho**2

    call theta_integrate (jacrho*(2.*dRdrho/Rr(2,:) + surf%shat/surf%r)/Rr(2,:)**2, num1)
    call theta_integrate (jacrho/Rr(2,:)**2*(1. + (bi/(grho*dpsidrho))**2), denom)
    call theta_integrate (dum, num2)

    dIdrho = bi*(num2 + num1)/denom

    deallocate (dum)

  end subroutine get_dIdrho

  subroutine get_djacdrho

    implicit none

    ! betaprim below is (4*pi*ptot/B0^2)*(-d ln ptot / drho)
    djacdrho = (Rr(2,:)/grho)**2*(2.*(dRdth*d2Rdrdth+dZdth*d2Zdrdth)/jacrho &
         + drzdth + jacrho*(bi*dIdrho/Rr(2,:)**2 - surf%betaprim)/dpsidrho**2)

  end subroutine get_djacdrho

  subroutine get_d2RZdr2

    implicit none

    ! get factor common to both d2R/drho2 and d2Z/drho2
    d2Rdr2 = ((djacdrho-dRdrho/Rr(2,:))/Rr(2,:) &
         - dRdrho*d2Zdrdth + dZdrho*d2Rdrdth)/(dRdth**2+dZdth**2)

    d2Zdr2 = d2Rdr2*dRdth
    d2Rdr2 = d2Rdr2*dZdth

  end subroutine get_d2RZdr2

  subroutine get_d2Idr2_d2jacdr2

    implicit none

    real :: denom, num1, num2, num3, num4, num5
    real, dimension (-ntg:ntg) :: tmp, fac1, fac2, fac3, mod, jacth, d2Rdr2dth, d2Zdr2dth

    tmp = jacrho/Rr(2,:)**2 + bi**2*jacrho/(Rr(2,:)*dpsidrho*grho)**2
    call theta_integrate (tmp, denom)
    denom = denom/bi

    call get_dthet (djacdrho/jacrho, jacth)

    mod = dRdth**2+dZdth**2
    fac1 = dRdth*d2Rdrdth+dZdth*d2Zdrdth
    fac2 = (d2Rdr2*dRdth+dRdrho*d2Rdrdth+d2Zdr2*dZdth+dZdrho*d2Zdrdth)/jacrho
    call get_dthet (fac2, fac3)

    tmp = -( (djacdrho-2.*jacrho*fac1/mod)/mod * (2.*fac1+jacrho*drzdth) &
         + jacrho/mod*(2.*(d2Rdrdth**2+dRdth*d2Rdr2dth+d2Zdrdth**2+dZdth*d2Zdr2dth) &
         - jacrho*drz*jacth + jacrho*fac3) ) &
         / Rr(2,:)**2
    call theta_integrate (tmp, num1)
    
    ! this is the first term in the 2nd radial derivative of the Jacobian
    d2jacdr2 = -tmp*Rr(2,:)**2

    ! betadbprim is (4*pi*ptot/B0^2)*(-d^2ptot/drho^2 / ptot)
    tmp = surf%betadbprim*jacrho**3/(mod*dpsidrho**2*Rr(2,:)**2)
    call theta_integrate (tmp, num2)

    d2jacdr2 = d2jacdr2 - tmp*Rr(2,:)**2

    tmp = surf%betaprim*jacrho**3/(mod*dpsidrho**2*Rr(2,:)**2) &
         * (3.*djacdrho/jacrho - 2.*fac1/mod)
    call theta_integrate (tmp, num3)

    d2jacdr2 = d2jacdr2 - tmp*Rr(2,:)**2

    tmp = -dIdrho*jacrho**3*bi*(2.*fac1/mod + dIdrho/bi + 3.*djacdrho/jacrho - 2.*dRdrho/Rr(2,:)) &
         / (mod*dpsidrho**2*Rr(2,:)**4)
    call theta_integrate (tmp, num4)

    d2jacdr2 = d2jacdr2 - tmp*Rr(2,:)**2

    tmp = 2.*jacrho/Rr(2,:)**3*(djacdrho/jacrho - 3.*dRdrho/Rr(2,:) + d2Rdr2) &
         + dqdr*jacrho/(surf%q*Rr(2,:)**2)*(djacdrho/jacrho - dqdr/surf%q - 2.*dRdrho/Rr(2,:)) &
         + jacrho*surf%d2qdr2/(surf%q*Rr(2,:)**2) + 2.*dRdrho*djacdrho/Rr(2,:)**2 &
         - dIdrho*jacrho/(Rr(2,:)**2*bi)*(djacdrho/jacrho - dIdrho/bi - 2.*dRdrho/Rr(2,:))
    call theta_integrate (tmp, num5)

    d2Idr2 = (num1+num2+num3+num4+num5)/denom
    d2jacdr2 = d2jacdr2 + d2Idr2*bi*jacrho**3/(mod*Rr(2,:)**2*dpsidrho**2)

  end subroutine get_d2Idr2_d2jacdr2

  subroutine get_dBdrho

    implicit none

    d2gpsidr2 = 2.*(dpsidrho*Rr(2,:)/jacrho)**2 &
         * (2.*dRdrho/Rr(2,:)-2.*djacdrho/jacrho) &
         * (dRdrho/Rr(2,:)-djacdrho/jacrho+dRdth*d2Rdrdth+dZdth*d2Zdrdth) &
         - dRdrho**2/Rr(2,:)**2 + d2Rdr2/Rr(2,:) + djacdrho**2/jacrho**2 &
         - d2jacdr2/jacrho + d2Rdrdth**2 + dRdth*d2Rdr2dth &
         + d2Zdrdth**2 + dZdth*d2Zdr2dth

    ! dB/drho
    dBdrho = ( (bi*dIdrho + dpsidrho**2*(grho**2*(dRdrho/Rr(2,:)-djacdrho/jacrho) &
         + (Rr(2,:)/jacrho)**2*(dRdth*d2Rdrdth+dZdth*d2Zdrdth))) &
         / (bmag*Rr(2,:)) - bmag*dRdrho) / Rr(2,:)

    ! get d/drho (dB/drho)
    d2Bdr2 = -dBdrho*dRdrho/Rr(2,:) + bmag*(dRdrho/Rr(2,:))**2 &
         - bmag*d2Rdr2/Rr(2,:) + 0.5*(2.*dIdrho**2 + 2.*bi*d2Idr2 &
         + d2gpsidr2)/(bmag*Rr(2,:)**2) &
         - (dBdrho + bmag*dRdrho/Rr(2,:))*(2.*dRdrho/Rr(2,:)+dBdrho/bmag)

  end subroutine get_dBdrho

  subroutine get_dcrossdr

    implicit none

    integer :: i
    real, dimension (-ntg:ntg) :: dgr2, dgrgt, dgt2, dgagr, dgagt, dga2

    dgr2 = 2.*(Rr(2,:)/jacrho)**2*(dRdrho/Rr(2,:)-djacdrho/jacrho &
         + dRdth*d2Rdrdth + dZdth*d2Zdrdth)
    dgrgt = -(Rr(2,:)/jacrho)**2*((2.*dRdrho/Rr(2,:)-2.*djacdrho/jacrho) &
         + d2Rdr2*dRdth+dRdrho*d2Rdrdth+d2Zdr2*dZdth+dZdrho*d2Zdrdth)
    dgt2 = (Rr(2,:)/jacrho)**2*(2.*(dRdrho/Rr(2,:)-djacdrho/jacrho) &
         + 2.*(dRdrho*d2Rdrdth + dZdrho*d2Zdrdth))
    ! this is d/drho (|grad alph|^2)
    ! will later multiply it by 0.5*dpsidrho**2
    dga2 = -2*dRdrho/Rr(2,:)**3 + dgr2*(varthet*dqdr+surf%q*dvarthdr) &
         + (2.0*grho**2*(varthet*dqdr+surf%q*dvarthdr) &
         + 2.*bi*jacrho*gradrho_gradthet/(dpsidrho*Rr(2,:)**2)) &
         *(surf%d2qdr2*varthet+2.*dqdr*dvarthdr+surf%q*d2varthdr2) &
         + 2.*(varthet*dqdr+surf%q*dvarthdr)*bi*jacrho/(dpsidrho*Rr(2,:)**2) &
         * (dgrgt + gradrho_gradthet*(dIdrho/bi + djacdrho/jacrho - 2.*dRdrho/Rr(2,:))) &
         + (bi*jacrho/(dpsidrho*Rr(2,:)**2))**2*(dgt2 + gradthet2*(2.*dIdrho/bi + 2.*djacdrho/jacrho &
         - 4.*dRdrho/Rr(2,:)))
    dgagr = -grho**2*(dvarthdr*dqdr+2.*varthet*dqdr+surf%q*d2varthdr2) &
         - dgr2*(varthet*dqdr+surf%q*dvarthdr) - bi*jacrho/(dpsidrho*Rr(2,:)**2) &
         * (dgrgt + gradrho_gradthet*(dIdrho/bi + djacdrho/jacrho - 2.*dRdrho/Rr(2,:)))

!    do i = -ntg, ntg
!       write (*,'(a3,6e12.4)') 'dg', dgr2(i), varthet(i)*dqdr, dvarthdr(i), d2varthdr2(i), dgrgt(i), bi*jacrho(i)/(dpsidrho*Rr(2,i)**2)
!    end do
    dgagt = -gradrho_gradthet*(2.*dvarthdr*dqdr+varthet*surf%d2qdr2+surf%q*d2varthdr2) &
         - dgrgt*(varthet*dqdr+surf%q*dvarthdr) - bi*jacrho/(dpsidrho*Rr(2,:)**2) &
         * (dgt2 + gradthet2*(dIdrho/bi + djacdrho/jacrho - 2.*dRdrho/Rr(2,:)))

    dcrossdr = dpsidrho*(dgagr*gradalph_gradthet+gradrho_gradalph*dgagt &
         - dga2*gradrho_gradthet - gradalph2*dgrgt)

    ! this is (dpsi/drho)^2*d|grad alpha|^2/dr
    dgds2dr = dga2*dpsidrho**2
    ! this is (dpsi/drho)^2*d(grad alpha . grad q)/dr
    dgds21dr = -dgagr*dqdr*dpsidrho**2
    ! this is (dpsi/drho)^2*d(|grad q|^2)/dr
    dgds22dr = (dqdr**2*dgr2 + 2.*grho**2*dqdr*surf%d2qdr2)*dpsidrho**2

    ! note that dkperp2/dr = (dpsi/drho*n0/a)^2*(dgds2dr+theta0*dgds21dr+theta0^2*dgds22dr)

  end subroutine get_dcrossdr

  subroutine get_varthet

    use constants, only: pi

    implicit none

    integer :: i
    real, dimension (-ntg:ntg) :: tmp

    call theta_integrate_indef(jacrho/Rr(2,:)**2, varthet)
    varthet = bi*varthet/(dpsidrho*surf%q)

  end subroutine get_varthet

  subroutine get_dvarthdr

    implicit none

    real, dimension (-ntg:ntg) :: dum

    dum = jacrho*(dIdrho/bi - dqdr/surf%q + djacdrho/jacrho &
         - 2.*dRdrho/Rr(2,:))/Rr(2,:)**2
    call theta_integrate_indef(dum, dvarthdr)
    dvarthdr = bi*dvarthdr/(dpsidrho*surf%q)

    dum = bi*jacrho/(surf%q*dpsidrho*Rr(2,:)**2)*((dIdrho/bi - dqdr/surf%q &
         + djacdrho/jacrho - 2.*dRdrho/Rr(2,:))**2 &
         + d2Idr2/bi - (dIdrho/bi)**2 - surf%d2qdr2/surf%q &
         + (dqdr/surf%q)**2 + d2jacdr2/jacrho - (djacdrho/jacrho)**2 &
         - 2.*d2Rdr2/Rr(2,:) + 2.*(dRdrho/Rr(2,:))**2)
    call theta_integrate_indef(dum, d2varthdr2)

  end subroutine get_dvarthdr

  subroutine theta_integrate (integrand, integral)

    implicit none

    real, dimension (-ntg:), intent (in) :: integrand
    real, intent (out) :: integral

    ! use trapezoidal rule to integrate in theta
    integral = 0.5*sum(delthet(:ntg-1)*(integrand(:ntg-1) + integrand(-ntg+1:ntg)))

  end subroutine theta_integrate

  ! get indefinite integral of integrand
  subroutine theta_integrate_indef (integrand, integral)

    implicit none

    real, dimension (-ntg:), intent (in) :: integrand
    real, dimension (-ntg:), intent (out) :: integral

    integer :: i

    ! use trapezoidal rule to integrate in theta
    integral(0) = 0.0
    do i = 1, ntg
       integral(i) = integral(i-1)+0.5*delthet(i-1)*(integrand(i-1)+integrand(i))
    end do
    do i = -1, -ntg, -1
       integral(i) = integral(i+1)-0.5*delthet(i)*(integrand(i+1)+integrand(i))
    end do

  end subroutine theta_integrate_indef

  subroutine derm(f, dfm, char)

    implicit none
    integer i, j
    character*1 :: char
    real f(:,:), dfm(:,:,:), pi

    pi = 2*acos(0.)
    
    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))         
    
    i=nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))
   
! assume up-down symmetry for now:
 
    select case (char)
    case ('E') 
       j=1
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j+1))
       
       j=nt
       dfm(:,j,2) = -0.5*(f(:,j-1)-f(:,j-1))
    case ('O')
       j=1
       dfm(:,j,2) = 0.5*(f(:,j+1)+f(:,j+1))
       
       j=nt      
       dfm(:,j,2) = -0.5*(f(:,j-1)+f(:,j-1))
    case ('T')
       j=1
       dfm(:,j,2) = f(:,j+1)
       
       j=nt      
       dfm(:,j,2) = pi - f(:,j-1)
    end select

    do i=2,nr-1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))
    enddo
    
    do j=2,nt-1
       dfm(:,j,2)=0.5*(f(:,j+1)-f(:,j-1))
    enddo
    
  end subroutine derm

  subroutine gradient(rgrid, theta, grad, char, rp, nth_used, ntm)

    use splines, only: inter_d_cspl
    implicit none
    
    integer nth_used, ntm
    character*1 char
    real rgrid(-ntm:), theta(-ntm:), grad(-ntm:,:)
    real tmp(2), aa(1), daa(1), rp, rpt(1)
    real, dimension(nr,nt,2) :: dcart
    integer i
    
    select case(char)
    case('D')  ! diagnostic 
       dcart = dbtcart
    case('P') 
       dcart = dpcart
    case('R') 
       dcart = dpcart  ! dpcart is correct for 'R'
    case('T')
       dcart = dtcart
    end select
    
    do i=-nth_used,-1
       call eqitem(rgrid(i), theta(i), dcart(:,:,1), tmp(1), 'R')
       call eqitem(rgrid(i), theta(i), dcart(:,:,2), tmp(2), 'Z')
       if(char == 'T') then
          grad(i,1)=-tmp(1)
          grad(i,2)=-tmp(2)
       else
          grad(i,1)=tmp(1)
          grad(i,2)=tmp(2)
       endif
    enddo

    do i=0,nth_used
       call eqitem(rgrid(i), theta(i), dcart(:,:,1), tmp(1), 'R')
       call eqitem(rgrid(i), theta(i), dcart(:,:,2), tmp(2), 'Z')
       grad(i,1)=tmp(1)
       grad(i,2)=tmp(2)
    enddo

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1)*0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1)*0.5*beta_0
       enddo
    endif

  end subroutine gradient

  subroutine bgradient(rgrid, theta, grad, char, rp, nth_used, ntm)

    use splines, only: inter_d_cspl
    implicit none
    
    integer nth_used, ntm
    character*1 char
    real rgrid(-ntm:), theta(-ntm:), grad(-ntm:,:)
    real :: aa(1), daa(1), rp, rpt(1)
    real, dimension(nr,nt,2) ::  dbish
    integer i

    select case(char)
    case('D')  ! diagnostic
       dbish = dbtbish
    case('P') 
       dbish = dpbish
    case('R') 
       dbish = dpbish  ! dpcart is correct for 'R'
    case('T')
       dbish = dtbish
    end select

    do i=-nth_used,nth_used
       call eqitem(rgrid(i), theta(i), dbish(:,:,1), grad(i,1), 'R')
       call eqitem(rgrid(i), theta(i), dbish(:,:,2), grad(i,2), 'Z')
    enddo

    if (char == 'T') then
       where (theta(-nth_used:nth_used) < 0.0)
          grad(-nth_used:nth_used,1) = -grad(-nth_used:nth_used,1)
          grad(-nth_used:nth_used,2) = -grad(-nth_used:nth_used,2)
       end where
    end if

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1) * 0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1) * 0.5*beta_0
       enddo
    endif

  end subroutine bgradient

  subroutine eqitem(r, theta_in, f, fstar, char)
      
    integer :: j, istar, jstar
    character*1 :: char
    real :: r, thet, fstar, tp, tps, theta_in
    real :: st, dt, pi
    real, dimension(:,:) :: f
    real, dimension(size(f,2)) :: mtheta
    
    pi = 2.*acos(0.)
! find r on psi mesh
    
    istar = 2

! Now do theta direction

    thet = mod2pi(theta_in)

! assume up-down symmetry

    tp=abs(thet)
    tps=1.
    if(char == 'Z' .and. thet /= 0.) tps=thet/abs(thet)
        
! get thet on theta mesh

    mtheta = (/ ( real(j-1)*pi/real(nt-1), j=1,nt) /)
  
! note that theta(1)=0 for local_eq theta 

    jstar=-1    
    do j=1,nt
       if(tp < mtheta(j)) then
          dt = tp - mtheta(j-1)
          st = mtheta(j) - tp
          jstar=j-1
          exit
       endif
       if(jstar /= -1) write(*,*) 'exit error j'
    enddo
      
! treat theta = pi separately
  
    if(jstar == -1) then
       jstar=nt-1
       dt=mtheta(jstar+1)-mtheta(jstar)
       st=0.
    endif

! use opposite area stencil to interpolate

    fstar=f(istar    , jstar    )  * st &
         +f(istar    , jstar + 1)  * dt 
    fstar=fstar*tps &
         /(mtheta(jstar+1)-mtheta(jstar))
!     write(*,*) i, dr, dt, sr, st
!     write(*,*) f(istar,jstar+1),f(istar+1,jstar+1)
!     write(*,*) f(istar,jstar),f(istar+1,jstar)
!     write(*,*) eqpsi(istar),eqpsi(istar+1)
!     write(*,*) mtheta(jstar),mtheta(jstar+1)
      
  end subroutine eqitem

  subroutine eqdcart(dfm, dfcart)
      
    implicit none

    real, dimension (:,:,:), intent(in)  :: dfm
    real, dimension (:,:,:), intent(out) :: dfcart
    real, dimension (size(dfm,1),size(dfm,2)) :: denom
    integer :: i, j
      
    denom(:,:) = drm(:,:,1)*dzm(:,:,2) - drm(:,:,2)*dzm(:,:,1)

    dfcart = 0.
    
    dfcart(:,:,1) =   dfm(:,:,1)*dzm(:,:,2) - dzm(:,:,1)*dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1)*drm(:,:,2) + drm(:,:,1)*dfm(:,:,2)
    
    do j=1,nt
       do i=2,nr
          dfcart(i,j,:)=dfcart(i,j,:)/denom(i,j)
       enddo
    enddo
    
  end subroutine eqdcart

  subroutine eqdbish(dcart, dbish)
    
    implicit none
    real, dimension(:, :, :), intent (in) :: dcart
    real, dimension(:, :, :), intent(out) :: dbish
    real, dimension(size(dcart,1),size(dcart,2)) :: denom
    integer :: i, j

    denom(:,:) = sqrt(dpcart(:,:,1)**2 + dpcart(:,:,2)**2)

    dbish(:,:,1) = dcart(:,:,1)*dpcart(:,:,1) + dcart(:,:,2)*dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1)*dpcart(:,:,2) + dcart(:,:,2)*dpcart(:,:,1)
    
    do j=1,nt
       do i=2,nr
          dbish(i,j,:) = dbish(i,j,:)/denom(i,j)
      enddo
    enddo

  end subroutine eqdbish

  subroutine leqcoefs (dgpdr, dgpbdr, dcvd0dr, dgbd0dr, dcvdr, dgbdr, dg2dr, dg21dr, dg22dr, dbdr)

    implicit none

    real, dimension (-ntg:), intent (out) :: dgpdr, dgpbdr, dcvd0dr, dgbd0dr, dcvdr, dgbdr
    real, dimension (-ntg:), intent (out) :: dg2dr, dg21dr, dg22dr, dbdr

    integer :: i

    dgpdr = dgradpardrho
    dgpbdr = dgradparbdrho
    dcvd0dr = dcvdrift0drho
    dgbd0dr = dgbdrift0drho
    dcvdr = dcvdriftdrho
    dgbdr = dgbdriftdrho
    dg2dr = dgds2dr
    dg21dr = dgds21dr
    dg22dr = dgds22dr
    dbdr = dBdrho

!    do i = -ntg, ntg
!       write (*,'(a9,17e12.4)') 'leqcoefs', theta(i), dgpdr(i), dgpbdr(i), dcvd0dr(i), dgbd0dr(i), dcvdr(i), dgbdr(i), &
!            dg2dr(i), dg21dr(i), dg22dr(i), dbdr(i), cvdrift0(i), dIdrho, dgradparbdrho(i), gradparb(i), dqdr, dpsidrho**2
!    end do

  end subroutine leqcoefs

  function invR (r, theta)
   
    real, intent (in) :: r, theta
    real :: invR
    
    invR=1./Rpos(r, theta)
    
  end function invR

  function rcenter(rp)

    real, intent(in) :: rp
    real :: rcenter

    rcenter = surf%R_center
    
  end function rcenter

  function Rpos (r, theta)
   
    use constants, only: pi

    real, intent (in) :: r, theta
    real :: Rpos
    real :: g, gp, dr
    
    dr = r - surf%r

! For Y Xiao: 
!    g = surf%delp/surf%r + surf%d * sin(theta)**2
!    Rpos = surf%R_center*(1.+r*(cos(theta)-g)-g*dr)

    g = cos(theta + surf%d * sin(theta-0.5*pi*surf%a))
    gp = -sin(theta + surf%d * sin(theta-0.5*pi*surf%a))*(surf%dp*sin(theta-0.5*surf%a)-surf%d*0.5*pi*surf%ap*cos(theta-0.5*pi*surf%a))
    
    Rpos = surf%R_center + surf%delp*dr + g*surf%r + (g+surf%r*gp)*dr
    
  end function Rpos

  function Zpos (r, theta)
   
    real, intent (in) :: r, theta
    real :: Zpos, dr
    
    dr = r - surf%r
    Zpos = surf%k*sin(theta)*surf%r + (surf%r*surf%kp + surf%k)*sin(theta)*dr
    
  end function Zpos

  function psi (r, theta)
   
    real, intent (in) :: r, theta
    real :: psi

    psi = r - surf%r
    
  end function psi

  function mod2pi (theta)
    
    real, intent(in) :: theta
    real :: pi, th, mod2pi
    real, parameter :: theta_tol = 1.e-6
    logical :: out
    
    pi=2.*acos(0.)
    
    if(theta <= pi .and. theta >= -pi) then
       mod2pi = theta
       return
    endif
    
    if(theta - theta_tol <= pi .and. theta >= -pi) then
       mod2pi = pi
       return
    endif

    if(theta <= pi .and. theta + theta_tol >= -pi) then
       mod2pi = -pi
       return
    endif

    th=theta
    out=.true.
    do while(out)
       if(th > pi) th = th - 2.*pi
       if(th <-pi) th = th + 2.*pi
       if(th <= pi .and. th >= -pi) out=.false.
    enddo
    mod2pi=th
    
  end function mod2pi
   
  function diameter (rp)
  
    real :: rp, diameter

    diameter = 2.*rp

  end function diameter

  function dbtori (pbar)
    real :: pbar, dbtori
    dbtori = 1.
  end function dbtori

  function btori (pbar)
    real :: pbar, btori
    btori = surf%r_geo
  end function btori

  function qfun (pbar)
    real :: pbar, qfun
    qfun = surf%q
  end function qfun

  function pfun (pbar)
    real :: pbar, pfun
    pfun = 0.5*beta_0
  end function pfun
  
  function dpfun (pbar)  
    real :: pbar, dpfun    

       dpfun = -1.

  end function dpfun

  function dpdrhofun(rho)

    real :: rho, dpdrhofun

    dpdrhofun = surf%pp

  end function dpdrhofun
  
  function betafun (pbar)  
    real :: pbar, betafun
    betafun = beta_0
  end function betafun

end module leq

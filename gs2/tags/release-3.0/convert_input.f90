program convert_input
  implicit none
  character(500) :: name

  integer :: loop1
  character(50) :: string
  integer :: ntheta, nperiod, ngauss, negrid, islow, nthetag
  real :: fcv, tpfac
  integer :: icv, kperiod
  real :: apfac
  real :: dif(5)
  integer :: nscreen, nout, nstep, ngamstep, isvar, nwrite
  integer :: iphi
  real :: width0
  real :: delt, ominst, phydif, power, tol
  real :: gamma, absom, fv, gridfac
  complex :: fexp1, fexp3
  real :: test1, test2, test3
  real :: bakdif(5)
  complex :: alpha, alpha1
  complex :: ala, ala1
  integer :: il00, ie00
  integer :: ncspec(5), ietg
  integer :: naky
  real :: akyarr(16)
  complex :: omegacom
  real :: ecut, we
  integer :: iwe
  real :: epsq, phiinit
  integer :: nbd
  real :: bdtheta0inp
  integer :: ith0mat
  integer :: nam, lnqdam
  integer :: ngk, lnqdgk
  real :: coefland
  integer :: itmax
  real :: err
  real :: rmin, rmax
  real :: delrho
  integer :: irho
  integer :: itor, iflux, ifulls, ivmom
  integer :: iq, ieik
  real :: alpmhdfac
  integer :: igridgen, npadd
  integer :: ivnew
  real :: vncoef
  integer ::iadd
  real :: filter, filtperp
  real :: driftknob

  real :: eps, epsl, pk
  real :: qinp, shat, shift, beta
  real :: rmaj, rhoc, airat, akappri, akappa, tri, tripri
  real :: ushift1, uprim(5)
  real :: feprim, zeff, zeffp
  real :: dbeam, fprim4, an5, fprim5
  real :: tprim(5)
  real :: vnewk(5)
  real :: amass(5)
  real :: temp(5)
  real :: z(5)
  real :: shiftpri, akappri2, tripri2, r_geo

  complex :: fexp(5)
  real :: fprim(5), an(5)
  real :: alpmhd, dpdrho
  real :: theta0
  real :: pi

  open (unit=1, file="name", status="old")
  read (unit=1, fmt="(a)") name

  open (unit=2, file=name, status="old")
  call read_gy0inp
  close (unit=2)

  read (unit=1, fmt="(a)") name
  open (unit=2, file=name, status="old")
  call read_gyinp
  close (unit=2)
  close (unit=1)

  call fiddle_with_parameters

  call print_conversion

contains

  subroutine read_gy0inp
    implicit none
    real :: fexp1r, fexp1i, fexp3r, fexp3i

    read (2,*) loop1
    read (2,"(a)") string
    read (2,*) ntheta, nperiod, ngauss, negrid, islow, nthetag
    read (2,*) fcv, tpfac, icv, kperiod, apfac
    read (2,*) dif(1), dif(2), dif(3), dif(4), dif(5)
    read (2,*) nscreen, nout, nstep, ngamstep, isvar, nwrite
    read (2,*) iphi, width0
    read (2,*) delt, ominst, phydif, power, tol
    read (2,*) gamma, absom, fv, gridfac
    read (2,*) fexp1r, fexp1i, fexp3r, fexp3i, test1, test2, test3
    fexp1 = cmplx(fexp1r,fexp1i)
    fexp3 = cmplx(fexp3r,fexp3i)
    read (2,*) bakdif(1), bakdif(2), bakdif(3), bakdif(4), bakdif(5)
    read (2,*) alpha, alpha1
    read (2,*) ala, ala1, il00, ie00
    read (2,*) ncspec(1), ncspec(2), ncspec(3), ncspec(4), ncspec(5), ietg
    read (2,*) naky
    read (2,*) akyarr(1),  akyarr(2),  akyarr(3),  akyarr(4)
    read (2,*) akyarr(5),  akyarr(6),  akyarr(7),  akyarr(8)
    read (2,*) akyarr(9),  akyarr(10), akyarr(11), akyarr(12)
    read (2,*) akyarr(13), akyarr(14), akyarr(15), akyarr(16)
    read (2,*) omegacom
    read (2,*) ecut, we, iwe
    read (2,*) epsq, phiinit
    read (2,*) nbd, bdtheta0inp
    read (2,*) nam, lnqdam
    read (2,*) ngk, lnqdgk
    read (2,*) coefland
    read (2,*) itmax, err
    read (2,*) rmin, rmax
    read (2,*) delrho, irho
    read (2,*) itor, iflux, ifulls, ivmom
    read (2,*) iq, ieik, alpmhdfac, igridgen, npadd
    read (2,*) ivnew, vncoef, iadd
    read (2,*) filter, filtperp
    read (2,*) driftknob
  end subroutine read_gy0inp

  subroutine read_gyinp
    implicit none
    read (2,*) eps, epsl, pk
    read (2,*) qinp, shat, shift, beta
    read (2,*) rmaj, rhoc, airat, akappri, akappa, tri, tripri
    read (2,*) ushift1, uprim(1), uprim(2), uprim(3), uprim(4), uprim(5)
    read (2,*) feprim, zeff, zeffp
    read (2,*) dbeam, fprim4, an5, fprim5
    read (2,*) tprim(1), tprim(2), tprim(3), tprim(4), tprim(5)
    read (2,*) vnewk(1), vnewk(2), vnewk(3), vnewk(4), vnewk(5)
    read (2,*) amass(1), amass(2), amass(3), amass(4), amass(5)
    read (2,*)  temp(1),  temp(2),  temp(3),  temp(4),  temp(5)
    read (2,*)     z(1),     z(2),     z(3),     z(4),     z(5)
    read (2,*) shiftpri, akappri2, tripri2, r_geo
  end subroutine read_gyinp

  subroutine fiddle_with_parameters
    implicit none
    real :: nelect, nduet, nimp, nduetp, nimpp

    if (ncspec(3) /= 0) then
       z(3) = -1.0
    else
       z(3) = 1.0
    end if

    if (ifulls == 0) then
       igridgen = 1
       npadd = 0
    end if

    eps = rhoc/rmaj
    epsl = fcv*2.0/rmaj
    pk = 2.0/(qinp*rmaj)

    fexp(1) = fexp1; fexp(2) = fexp1; fexp(5) = fexp1
    fexp(3) = fexp3; fexp(4) = fexp3

    nelect = 1.0
    nduetp = nelect*(zeffp + (zeff-z(2))*feprim)/(1.0-z(2))
    nduetp = nduetp - nelect*dbeam*fprim4
    nimpp = nelect*(zeffp + (zeff-1.0)*feprim)/(z(2)*(z(2)-1.0))
    nduet = nelect*(z(2)-zeff)/(z(2)-1.0) - nelect*dbeam
    nimp = nelect*(zeff-1.0)/(z(2)*(z(2)-1.0))

    an(1) = nduet
    fprim(1) = nduetp/nduet
    tprim(1) = tprim(1)*tpfac

    an(2) = nimp
    if (zeff > 1.0 + 2.0*epsilon(0.0)) then
       fprim(2) = nimpp/nimp
    else
       fprim(2) = 0.0
    end if
    tprim(2) = tprim(2)*tprim(1)

    an(3) = nelect
    fprim(3) = feprim
    tprim(3) = tprim(3)*tpfac
    an(4) = dbeam
    if (islow == 1) an(4) = 0.5*dbeam
    fprim(4) = fprim4
    tprim(4) = tprim(4)*tpfac
    an(5) = an5
    fprim(5) = fprim5
    tprim(5) = tprim(5)*tprim(1)

    alpmhd = sum((fprim(1:4) + tprim(1:4))*an(1:4)*temp(1:4))

    dpdrho = alpmhd
    if (alpmhdfac > 2.0*epsilon(0.0)) then
       alpmhd = alpmhd*beta*qinp**2*rmaj
       shift = -alpmhd*alpmhdfac
    end if

    pi = 2.0*acos(0.0)
    theta0 = bdtheta0inp*pi

  end subroutine fiddle_with_parameters

  subroutine print_conversion
    implicit none
    integer :: ik, is, ispec

    print "('&collisions_knobs')"
    if (vnewk(3) < -epsilon(0.0)) then
       print *, "collision_model='krook'"
    else
       print *, "collision_model='lorentz'"
    end if
    vnewk = abs(vnewk)
    print *, "vncoef=", vncoef
    print *, "absom=", absom
    print *, "ivnew=", ivnew
    print *, "conserve_number=", iadd /= 0
    print *, "conserve_momentum=.true."
    print "('/')"

    print "('&dist_fn_knobs')"
    print *, "kperiod=", kperiod
    print *, "gridfac=", gridfac
    print *, "apfac=", apfac
    print *, "driftknob=", driftknob
    print "('/')"

    print "('&source_knobs')"
    print *, "source_option='default'"

    print "('&init_g_knobs')"
    print *, "ginit_option='default'"
    print *, "width0=", width0
    print *, "phiinit=", phiinit
    print "('/')"

    print "('&fields_knobs')"
    print *, "field_option='implicit'"
    print "('/')"

    print "('&gs2_diagnostics_knobs')"
    if (nwrite > 0) then
       string = ".true."
    else
       string = ".false."
    end if
    print *, "print_line=", string
    print *, "write_line=", string
    print *, "write_phi=", string
    print *, "write_apar=", string
    print *, "write_aperp=", string
    print *, "write_omega=", string
    print *, "write_omavg=", string
    print *, "write_qheat=", string
    print *, "write_pflux=", string
    print *, "write_vflux=", string
    print *, "write_qmheat=", string
    print *, "write_pmflux=", string
    print *, "write_vmflux=", string
    print *, "write_dmix=", string
    print *, "write_kperpnorm=", string
    print *, "write_phitot=", string
    print *, "nwrite=", nscreen
    print *, "navg=", max(nscreen,10)
    print *, "omegatol=", tol
    print *, "print_old_units=", .true.
    print *, "write_eigenfunc=", .true.
    print *, "write_final_fields=", .false.
    print *, "igomega=", nout
    print "('/')"

    print "('&le_grids_knobs')"
    print *, "ngauss=", ngauss
    print *, "negrid=", negrid
    print *, "ecut=", ecut
    print "('/')"

    print "('&parameters')"
    print *, "beta=", beta
    print *, "zeff=", zeff
    print "('/')"

    print "('&kt_grids_knobs')"
    print *, "grid_option='specified'"
    print "('/')"

    print "('&kt_grids_specified_parameters')"
    print *, "naky=", naky
    print *, "ntheta0=", 1
    print "('/')"

    do ik = 1, naky
       write (string, *) ik
       string = trim(adjustl(string))
       print "('&kt_grids_specified_element_', a)", string
       print *, "aky=", akyarr(ik)
       print *, "theta0=", theta0
       print "('/')"
    end do

    print "('&knobs')"
    print *, "fphi=", test1
    print *, "fapar=", test2
    print *, "faperp=", test3
    print *, "delt=", delt
    print *, "nstep=", nstep
    print "('/')"

    is = 0
    do ispec = 1, 5
       if (ncspec(ispec) == 0) cycle
       is = is + 1
       write (string, *) is
       string = trim(adjustl(string))
       print "('&species_parameters_', a)", string
       print *, "z=", z(ispec)
       print *, "mass=", amass(ispec)
       print *, "dens=", an(ispec)
       print *, "temp=", temp(ispec)
       print *, "tprim=", tprim(ispec)
       print *, "fprim=", fprim(ispec)
       print *, "uprim=", uprim(ispec)
       print *, "vnewk=", vnewk(ispec)
       if (ispec == 3) then
          print *, "type='electron'"
       else if (ispec == 4 .and. islow /= 0) then
          print *, "type='slowing-down'"
       else
          print *, "type='ion'"
       end if
       print "('/')"
       print "('&dist_fn_species_knobs_', a)", string
       print *, "fexpr=", real(fexp(ispec))
       print *, "fexpi=", aimag(fexp(ispec))
       print *, "bakdif=", bakdif(ispec)
       print "('/')"
    end do

    print "('&species_knobs')"
    print *, "nspec=", is
    print "('/')"

    print "('&theta_grid_parameters')"
    print *, "rhoc=", rhoc
    print *, "rmaj=", rmaj
    print *, "r_geo=", r_geo
    print *, "eps=", abs(eps)
    print *, "epsl=", epsl
    print *, "qinp=", qinp
    print *, "shat=", shat
    print *, "pk=", pk
    print *, "shift=", shift
    print *, "akappa=", akappa
    print *, "akappri=", akappri
    print *, "tri=", tri
    print *, "tripri=", tripri
    print *, "ntheta=", nthetag
    print *, "nperiod=", nperiod
    print "('/')"

    print "('&theta_grid_knobs')"
    if (ifulls == 0) then
       print *, "equilibrium_option='s-alpha'"
    else if (igridgen == 0) then
       print *, "equilibrium_option='grid.out'"
    else
       print *, "equilibrium_option='eik'"
    end if
    print "('/')"

    print "('&theta_grid_eik_knobs')"
    print *, "itor=", itor
    print *, "iflux=", iflux
    print *, "irho=", irho
    print *, "equal_arc=.true."
    print *, "ppl_eq=.false."
    print *, "gen_eq=.false."
    print *, "vmom_eq=.false."
    print *, "efit_eq=.false."
    print *, "bishop=5"
    print *, "s_hat_input=", shat
    print *, "alpha_input=", alpmhd
    print *, "!invLp_input="
    print *, "!beta_prime_input="
    print *, "dp_mult=1.0"
    print *, "delrho=", delrho
    print *, "rmin=", rmin
    print *, "rmax=", rmax
    print *, "ismooth=", 0
    print *, "!ak0"
    print *, "!k1"
    print *, "!k2"
    print *, "isym=0"
    print *, "writelots=.false."
    print "('/')"

    print "('&theta_grid_file_knobs')"
    print *, "gridout_file='grid.out'"
    print "('/')"

    print "('&theta_grid_gridgen_knobs')"
    print *, "npadd=", npadd
    print *, "alknob=0.0"
    print *, "epsknob=1e-5"
    print *, "extrknob=0.0"
    print *, "tension=1.0"
    print *, "thetamax=0.0"
    print *, "deltaw=0.0"
    print *, "widthw=1.0"
    print "('/')"

    print "('&theta_grid_salpha_knobs')"
    if (icv == 0) then
       print *, "model_option='no-curvature'"
    else if (eps < 0.0) then
       print *, "model_option='alpha1'"
    else
       print *, "model_option='s-alpha'"
    end if
    print *, "alpmhdfac=", alpmhdfac
    print *, "alpha1=", real(alpha1)
    print "('/')"

  end subroutine print_conversion

end program convert_input

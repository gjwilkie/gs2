program test
  use simpledataio
  use simpledataio_write
  implicit none
  type (sdatio_file) :: sdatfile
!  type (sdatio_variable), pointer :: svar
!  type (sdatio_variable) :: svar2
  real :: parameter1 = 0.5
  double precision, dimension(2) :: yvar = (/0.1d0, 0.3d0/)
  real, dimension(2) :: floatvar = (/0.1,0.3/)
  integer, dimension(2) :: iy = (/1,2/)
  !double precision, dimension(3,2) ::  phivar = (/(/0.1,0.3/), (/2.0, 4.0/), (/-1.0, 3.6/)/);
  double precision, dimension(3,2) ::  phivar = reshape((/0.1d0,2.0d0,-1.0d0, 0.3d0,4.0d0, 3.6d0/), (/3,2/))
  complex*16 , dimension(3,2) ::  phicomp
  complex :: parameter_comp 
  double precision, dimension(2) :: phi_tvar
  !write (*, *) 'phivar', phivar(2,1)
  integer :: i
  double precision :: t
  integer, dimension(2) ::  idxs
  integer, dimension(2) :: idxs2 = (/2,1/)
  double precision ::  val = 32.9
  !complex :: ZI = 

  phicomp = cmplx(phivar, phivar*2.0d0) ! for some reason this results in a loss
    ! of precision
  phicomp =  phivar + phivar*cmplx(0.0,2.0)
  parameter_comp =  cmplx(4.0, 5.0)
  write (*,*) 'phicomp(1,1) is ', phicomp(1,1)

  call createfile(sdatfile, "test.cdf")
  call add_dimension(sdatfile, "x", 3, "The x coordinate", "m")
  call add_dimension(sdatfile, "y", 2, "The y coordinate", "m")
  call add_dimension(sdatfile, "r", 2, "Real and imaginary parts", "")
  call add_dimension(sdatfile, "t", SDATIO_UNLIMITED, "The time coordinate", "s")
  call print_dimensions(sdatfile)

  call create_variable(sdatfile, SDATIO_DOUBLE, "phi", "xy", "Some potential", "Vm")
  call create_variable(sdatfile, SDATIO_DOUBLE, "phicomp", "rxy", "Some complex potential", "Vm")
  call create_variable(sdatfile, SDATIO_DOUBLE, "parameter", "", "A scalar parameter", "(none)")
  call create_variable(sdatfile, SDATIO_DOUBLE, "parameter_comp", "r", "A complex parameter", "(none)")
  call create_variable(sdatfile, SDATIO_FLOAT, "floatvar", "y", "A single precision variable.", "Vm")
  call create_variable(sdatfile, SDATIO_DOUBLE, "phi_t", "yt", "Some potential as a function of y and time", "Vm")
  call create_variable(sdatfile, SDATIO_DOUBLE, "phi_txy", "rxyt", "Some complex potential as a function of y and time", "Vm")
  call create_variable(sdatfile, SDATIO_DOUBLE, "y", "y", "Values of the y coordinate", "m")
  call create_variable(sdatfile, SDATIO_DOUBLE, "t", "t", "Values of the time coordinate", "m")

  call write_variable(sdatfile, "parameter", parameter1)
  call write_variable(sdatfile, "parameter_comp", parameter_comp)

  call print_variables(sdatfile)
  call create_variable(sdatfile, SDATIO_INT, "iky", "y", "y index values", "(none)")

  call write_variable(sdatfile, "y", yvar)
  call write_variable(sdatfile, "iky", iy)
  call write_variable(sdatfile, "phi", phivar)
  call write_variable(sdatfile, "phicomp", phicomp)
  call write_variable(sdatfile, "floatvar", floatvar)
  
  call set_count(sdatfile, "phi_txy", "x", 1)
  call set_count(sdatfile, "phi_txy", "y", 1)
  call write_variable(sdatfile, "phi_txy", phicomp)


  do i = 1,6
    t = 0.3d0 + real(i);
    phi_tvar(1) = 4 + real(i)/2.0;
    phi_tvar(2) = 6 + real(i)*3.0; 
    call write_variable(sdatfile, "t", t);
    call write_variable(sdatfile, "phi_t", phi_tvar);
    
    call set_offset(sdatfile, "phi_txy", "x", 1)

    call set_start(sdatfile, "phi_txy", "x", 2)
    call write_variable_with_offset(sdatfile, "phi_txy", phicomp)

    call set_start(sdatfile, "phi_txy", "y", 2)
    call write_variable(sdatfile, "phi_txy", parameter_comp)
    call set_start(sdatfile, "phi_txy", "y", 1)

    call set_start(sdatfile, "phi_txy", "x", 3)
    call write_variable_with_offset(sdatfile, "phi_txy", phicomp)


    call increment_start(sdatfile, "t");
    ! if (i>2) stop
  end do


  write (*,*) 'y', yvar
  write(*,*) "This should be t: ", variable_exists(sdatfile, "t")
  write(*,*) "This should be t: ", variable_exists(sdatfile, "phi_t")
  write(*,*) "This should be f: ", variable_exists(sdatfile, "tbbb")

  !find_variable(sdatfile, "parameter")
  !svar2 = svar(1)
  !write (*,*) 'id', svar%type_size
  !call write_variable(sdatfile, "parameter", c_address(parameter1));
  call closefile(sdatfile)

end program test

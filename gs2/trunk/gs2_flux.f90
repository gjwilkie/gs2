module gs2_flux
  
  use regression, only: reg_type
  implicit none
  public :: init_gs2_flux, reset_init
  
  real :: delta, eps
  integer :: nflux
  integer :: flux_option_switch
  integer, parameter :: flux_option_none = 1
  integer, parameter :: flux_option_adjust = 2

  real, dimension (:), allocatable :: qflux, t_hist, qflux_last, tprim_last, qin
  real, dimension (:), allocatable :: qhi, tphi, qlo, tplo
  real, dimension (:,:), allocatable :: qflux_hist

  real, dimension (:,:), allocatable :: qold, tpold
  
  integer, parameter :: nfluxsave = 20
  integer :: out_unit, out2_unit
  logical, private :: initialized = .false.
  logical :: gs2_flux_adjust = .false.
  
  type (reg_type), dimension(:), allocatable :: reg

contains
  
  subroutine init_gs2_flux

    use species, only: nspec, spec, init_species
    use file_utils, only: open_output_file
    use mp
    
    if (initialized) return
    initialized = .true.
    
    call init_species
    allocate (qin(nspec)) ; qin = 0.

    call read_parameters

    select case (flux_option_switch)
    case (flux_option_none)
       return
    case (flux_option_adjust)
       if (proc0) then
          call open_output_file (out_unit, ".cf")
          call open_output_file (out2_unit, ".cf2")
       end if

       gs2_flux_adjust = .true.
       allocate (qflux(nspec)) ; qflux = 0.
       allocate (qflux_last(nspec)) ; qflux_last = 0.
       allocate (tprim_last(nspec)) ; tprim_last = spec%tprim
       allocate (t_hist(0:nflux-1)) ; t_hist = 0.
       allocate (qflux_hist(0:nflux-1,nspec)) ; qflux_hist = 0.
       allocate (qold(nfluxsave, nspec)) ; qold = 0.
       allocate (tpold(nfluxsave, nspec)) ; tpold = 0.
       allocate (reg(nspec))
       allocate (qhi(nspec), tphi(nspec)) ; qhi = 0. ; tphi = 0.
       allocate (qlo(nspec), tplo(nspec)) ; qlo = 0. ; tplo = 0.
    end select
    
  end subroutine init_gs2_flux
  
  subroutine read_parameters 
    use mp
    use text_options, only: text_option, get_option_value
    use file_utils, only: input_unit, input_unit_exist, error_unit
    use file_utils, only: get_indexed_namelist_unit
    use species, only: nspec
    integer :: ierr, in_file, is, unit
    logical :: exist
    real :: heat
    type (text_option), dimension (2), parameter :: fluxopts = &
         (/ text_option('default', flux_option_none), &
            text_option('adjust', flux_option_adjust) /)
    character (20) :: flux_option

    namelist /gs2_flux_knobs/ nflux, flux_option, delta, eps
    namelist /flux_target/ heat
    
    if (proc0) then
       eps = 0.5
       delta = 0.1
       nflux = 50
       flux_option = 'default'
       
       in_file=input_unit_exist("gs2_flux_knobs", exist)
       if (exist) read (unit=input_unit("gs2_flux_knobs"), nml=gs2_flux_knobs)

       ierr = error_unit()
       call get_option_value &
            (flux_option, fluxopts, flux_option_switch, &
            ierr, "flux_option in gs2_flux_knobs")
    end if
       
    call broadcast (eps)
    call broadcast (delta)
    call broadcast (nflux)
    call broadcast (flux_option_switch)
    
    if (flux_option_switch == flux_option_none) return

    if (proc0) then
       do is=1,nspec
          call get_indexed_namelist_unit (unit, "flux_target", is)
          heat = 50.
          read (unit=unit, nml=flux_target)
          close (unit=unit)
          qin(is) = heat
       end do
    end if

    call broadcast (qin)

  end subroutine read_parameters
  
  subroutine check_flux (i, t, heat_fluxes)

    use mp, only: proc0, broadcast
    use dist_fn, only: reset_physics
    use species, only: nspec, spec
    use file_utils, only: flush_output_file
    use regression
    integer, intent (in) :: i
    real, intent (in) :: t
    real, dimension(:), intent (in) :: heat_fluxes
    real, dimension (nspec) :: tp
    real :: t_interval
    integer :: j
    integer :: isave = 0, iuse
    logical :: first = .true.
    logical :: second = .false., third = .false.
    logical :: bi

    if (flux_option_switch == flux_option_none) return
    
    if (proc0) then
       qflux_hist(mod(i-1,nflux),:) = heat_fluxes
       t_hist(mod(i-1,nflux)) = t
    end if

    if (mod(i, nflux) /= 0) return
       
    if (proc0) then

       qflux = 0.
       t_interval = 0.
       do j=nflux/2,nflux-1
          qflux = qflux + qflux_hist(j,:)*(t_hist(j)-t_hist(j-1))
          t_interval = t_interval + t_hist(j)-t_hist(j-1)
       end do
       qflux = qflux / t_interval

       write (out_unit, "('t= ',e10.4,' heat flux= ',1(1x,e10.4),&
            &' tprim= ',1(1x,e10.4))") t, qflux(1), spec(1)%tprim
       call flush_output_file(out_unit, ".cf")

       if (first) then
          tp = (1.+delta) * spec%tprim
          first = .false. ; second = .true. ; third = .false.
       elseif (second) then
          isave = isave + 1
          iuse = mod(isave-1,nfluxsave)+1
          qold(iuse,:) = (qflux+qflux_last)/2.
          tpold(iuse,:) = (spec%tprim+tprim_last)/2.

! estimate upper and lower bounds
          do j=1,nspec
             if (qold(iuse,j) > qin(j)) then
                if (qhi(j) == 0.) then
                   qhi(j) = qold(iuse,j)
                   tphi(j) = tpold(iuse,j)
                elseif (qold(iuse,j) < qhi(j)) then
                   qhi(j) = 0.5*(qhi(j)+qold(iuse,j))
                   tphi(j) = 0.5*(tphi(j)+tpold(iuse,j))
                end if
             elseif (qold(iuse,j) < qin(j)) then
                if (qlo(j) == 0.) then
                   qlo(j) = qold(iuse,j)
                   tplo(j) = tpold(iuse,j)
                elseif (qold(iuse,j) > qlo(j)) then
                   qlo(j) = 0.5*(qlo(j)+qold(iuse,j))
                   tplo(j) = 0.5*(tplo(j)+tpold(iuse,j))
                end if
             end if
          end do

!          write (*,*) 'inserting tp= ',tpold(iuse,1),' and q = ',qold(iuse,1),' into regression.'
          call regress (reg, tpold(iuse,:), qold(iuse,:))

! Newton estimate averaged against present value
          tp = spec%tprim*(1.-delta*(qold(iuse,:)-qin)/(qflux-qflux_last))
          tp = eps*tp+(1.-eps)*spec%tprim
          delta = -delta*0.95

          first = .true.; second = .false. ; third = .false.

          if (iuse > 1) then
             do j=1,nspec
                bi = .true.
                if (qhi(j) == 0. .or. qlo(j) == 0.) bi = .false.
                if (abs(qold(iuse,j)) < abs(qin(j)/5.)) then
                   tp(j)=1.5*spec(j)%tprim
                   write (*,*) 'Attempting to blast out'
                   cycle
                end if

! try regression estimate occasionally.
                if (mod(iuse,5) == 0) bi = .false.

                if (bi) then
                   tp(j) = 0.5*(tplo(j)+tphi(j))
                   write (*,*) 'Used bisection estimate: ',tplo(j), tphi(j), qlo(j), qhi(j)
                else
                   tp(j) = xp(reg(j), qin(j))
                   write (*,*) 'Used regression estimate'
                end if
                tp(j) = eps*tp(j)+(1.-eps)*spec(j)%tprim
             end do
             first = .true.; second = .false. ; third = .false.
             write (out2_unit, "('t= ',e10.4,' a= ',1(1x,e10.4),&
                  &' b= ',1(1x,e10.4),' Tpred= ',1(1x,e10.4))") t, reg%a, reg%b, xp(reg, qin)
             call flush_output_file(out2_unit, ".cf2")
          end if
       elseif (third) then
          ! nothing
       end if
       
       write (*,*) qflux, spec%tprim, tp
       qflux_last = qflux
    end if

    call broadcast (tp)
    tprim_last = spec%tprim
    spec%tprim=tp

    call reset_physics


  end subroutine check_flux
  
  subroutine reset_init

    implicit none

    deallocate (qin)
    initialized = .false.

  end subroutine reset_init

end module gs2_flux

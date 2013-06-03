module dist_fn

  implicit none

  public :: init_dist_fn, finish_dist_fn
  public :: M_class, N_class, i_class

  ! source term in implicit solve
  complex, dimension (:,:,:), allocatable :: source

  ! response matrix in g
  real, dimension (:,:,:), allocatable :: gresponse

  ! matrix needed for reponse matrix approach
  real, dimension (:,:), allocatable :: m_mat

  ! coefficients multiplying g_{i,j}^{n+1}, g_{i+1,j}^{n+1}, g_{i-1,j}^{n+1}, etc.
  real, dimension (:,:), allocatable :: dvpafac
  real, dimension (:,:,:), allocatable :: dthetfac

  ! needed for twist and shift BC during field solve
  integer :: i_class
  integer, dimension(:), allocatable :: M_class, N_class

contains

  subroutine init_dist_fn

    use gs3_input, only: nvgrid, ntgrid
    use dist_fn_arrays, only: gold, gnew, g
    use gs3_layouts, only: g_lo, imu_idx
    use constants, only: zi
    use vpa_grid, only: vpa
    use mu_grid, only: mu
    use theta_grid, only: theta, bmag
    use kt_grids, only: naky, ntheta0

    implicit none

    integer :: ig, iv, imu, iglo

    if (.not. allocated(gold)) then
       allocate (gold(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; gold = 0.0
       allocate (gnew(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; gnew = 0.0
       allocate (g(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; g = 0.0
       allocate (source(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; source = 0.0
       allocate (gresponse(nvgrid,nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; gresponse = 0.0
       allocate (m_mat(nvgrid,nvgrid)) ; m_mat = 0.0
    end if

    call init_implicit_solve

    i_class = 1
    if (.not. allocated(M_class)) then
       allocate (M_class(i_class))
       allocate (N_class(i_class))
    end if
    M_class = naky*ntheta0 ; N_class = 1

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       imu = imu_idx(g_lo,iglo)
       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             gold(ig,iv,iglo) = exp(-vpa(iv)**2-mu(imu)*bmag(ig))*exp(-theta(ig)**2)
          end do
       end do
    end do

    gnew = gold

  end subroutine init_dist_fn

  subroutine advance_dist_fn

    use dist_fn_arrays, only: gold, gnew

    implicit none

    ! solve dg/dt + vpa dg/dtheta = 0 implicitly with bidiagonal scheme
    call implicit_solve

    ! update gold
    gold = gnew

  end subroutine advance_dist_fn

  subroutine init_implicit_solve

    use gs3_input, only: t_impfac, thet_impfac, vpa_impfac, dt, ntgrid, nvgrid
    use gs3_layouts, only: g_lo, imu_idx
    use theta_grid, only: delthet, theta
    use vpa_grid, only: dvpa, vpa
    use mu_grid, only: mu

    implicit none

    integer :: iglo, imu

    if (.not. allocated(dthetfac)) then
       allocate (dthetfac(-ntgrid:ntgrid,-nvgrid:nvgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; dthetfac = 0.0
       allocate (dvpafac(-ntgrid:ntgrid,-nvgrid:nvgrid)) ; dvpafac = 0.0
    end if

    ! need to add in multiplication by b . grad theta
    dvpafac = dt*spread(vpa,1,2*ntgrid+1)/spread(delthet,2,2*nvgrid+1)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       ! get im corresponding to iglo
       imu = imu_idx(g_lo,iglo)

       ! 0.1*sin(theta) is currently a placeholder for (b . grad B)
       dthetfac(:,:,iglo) = 0.5*dt*mu(imu)/spread(dvpa,1,2*ntgrid+1)*(0.1*spread(sin(theta),2,2*nvgrid+1))

    end do

    call get_gresponse_matrix

  end subroutine init_implicit_solve

  subroutine get_gresponse_matrix

    use gs3_input, only: ntgrid, nvgrid
    use gs3_layouts, only: g_lo
    use luinverse, only: invert_matrix
    use dist_fn_arrays, only: gnew

    implicit none

    integer :: ig, iv, iglo, i
    real, dimension (:,:), allocatable :: p_mat

    allocate (p_mat(nvgrid,nvgrid)) ; p_mat = 0.0

    ! ensure that gnew and source are initialized to zero
    gnew = 0.0 ; source = 0.0

    m_mat = 0.0

    ! get response of g(least negative theta, vpa > 0) to unit impulse in g(theta=0, vpa > 0)
    ! this will be used to calculate g(theta=0, vpa > 0) each time step, which will
    ! start the implicit sweep in theta and vpa;
    ! i.e. g(theta=-delthet,vpa>=0)^{n+1} = gresponse * g(theta=0,vpa>0)^{n+1} + r(theta=-delthet)^{n}
    ! the r(theta=-delthet)^{n} term will be obtained at each time step with a sweep setting g(theta=0,vpa>0)^{n+1}=0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do i = 1, nvgrid
          ! give a unit impulse to g at this vpa
          gnew(0,i,iglo) = 1.0

          ! sweep over (theta,vpa) plane and get values of gnew at theta=-delthet, vpa>0
          call implicit_sweep (iglo)

          ! fill in ith column of initial (nvgrid) x (nvgrid) x (nmu) response matrix
          gresponse(:,i,iglo) = real(gnew(-1,1:nvgrid,iglo))

          ! reset g to zero everywhere
          gnew(:,:,iglo) = 0.0
       end do

       ! get the matrices multiplying g(theta=0,vpa>0) and g(theta=-delthet,vpa>0) 
       ! in the implicit solve. these are necessary to obtain the final response matrix

       do iv = 1, nvgrid
          m_mat(iv,iv) = -dvpafac(-1,iv)
          p_mat(iv,iv) = 1.0 + dvpafac(-1,iv)
       end do

       gresponse(:,:,iglo) = p_mat + matmul(m_mat,gresponse(:,:,iglo))

       ! take the inverse to get the final response matrix
       call invert_matrix (gresponse(:,:,iglo))

       p_mat = 0.0
    end do

    deallocate (p_mat)

  end subroutine get_gresponse_matrix

  ! solve dg/dt + vpa dg/delthet = 0 implicitly with bidiagonal scheme; i.e.,
  ! g_{i+1}^{n+1}*(t_impfac + thet_impfac*(vpa*dt/dz)) + g_{i}^{n+1}*(1-t_impfac - thet_impfac*(vpa*dt/dz))
  ! = g_{i+1}^{n}*(t_impfac - (1-thet_impfac)*(vpa*dt/dz)) + g_{i}^{n}*(1-t_impfac + (1-thet_impfac)*(vpa*dt/dz))
  subroutine implicit_solve

    use gs3_input, only: ntgrid, nvgrid
    use gs3_layouts, only: g_lo
    use dist_fn_arrays, only: gnew

    implicit none

    integer :: iglo

    real, dimension (:), allocatable :: dgdgn, source_mod

    allocate (dgdgn(nvgrid)) ; dgdgn = 0.0
    allocate (source_mod(nvgrid)) ; source_mod = 0.0

    call get_source

    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       ! initially set g^{n+1}(theta=0,vpa>=0) to zero
       gnew(0,1:nvgrid,iglo) = 0.0
       
       ! sweep through vpa and theta once to obtain the response of g(theta=-delthet,vpa>0)^{n+1}
       ! to g^{n} when g^{n+1}(theta=0,vpa>0) = 0
       call implicit_sweep (iglo, dgdgn)
       
       ! first obtain g(theta=0,vpa>0) using pre-computed response matrix
       source_mod = source(0,1:nvgrid,iglo) - dgdgn
       
       ! get g^{n+1}(theta=0,vpa>0) using response matrix
       gnew(0,1:nvgrid,iglo) = matmul(gresponse(:,:,iglo), source_mod)
       
       ! with g^{n+1}(theta=0,vpa>0) specified, sweep once more to get g^{n+1} everywhere else
       call implicit_sweep (iglo)
       
       ! reset dgdgn to zero
       dgdgn = 0.0

    end do

    deallocate (dgdgn, source_mod)

  end subroutine implicit_solve

  ! implicit_sweep starts with a boundary condition for g^{n+1} along the line (theta=0,vpa>=0)
  ! as well as zero BCs at (theta=-theta_max,vpa>0), (theta=theta_max,vpa<0), (vpa=-vpa_max,theta<0),
  ! and (vpa=vpa_max,theta>0).  note that dvpa/dt<0 for theta>0 and dvpa/dt>0 for theta<0.
  ! and solves for g^{n+1}
  subroutine implicit_sweep (iglo, dgdgn)

    use gs3_input, only: ntgrid, nvgrid
    use dist_fn_arrays, only: gnew

    implicit none

    integer, intent (in) :: iglo
    real, dimension (:), intent (in out), optional :: dgdgn

    integer :: ig, iv

    ! initialize gnew to zero for all (theta,vpa) except (theta=0,vpa>0)
    ! which was set earlier as initial condition
    gnew(-ntgrid:-1,:,iglo) = 0.0 ; gnew(1:,:,iglo) = 0.0 ; gnew(0,-nvgrid:0,iglo) = 0.0

    ! first calculate gnew at theta=vpa=0, which is decoupled from other 
    ! theta-vpa points
    gnew(0,0,iglo) = source(0,0,iglo)

    ! boundary condition for vpa > 0 particles is g(theta=-infinity) = 0
    ! boundary condition for theta > 0 (corresponds to dvpa/dt < 0) is 
    ! g(vpa=infinity) = 0
    do iv = nvgrid-1, 0, -1
       do ig = 1, ntgrid
          gnew(ig,iv,iglo) = ( source(ig,iv,iglo) + gnew(ig-1,iv,iglo)*dvpafac(ig-1,iv) &
               + gnew(ig,iv+1,iglo)*dthetfac(ig,iv,iglo) ) / (1.0 + dvpafac(ig-1,iv) + dthetfac(ig,iv,iglo))
       end do
    end do
    
    ! boundary condition for vpa < 0 particles is g(theta=infinity) = 0
    ! boundary condition for theta > 0 (dvpa/dt < 0) is g(vpa=infinity)=0
    do iv = -1, -nvgrid, -1
       do ig = ntgrid-1, 0, -1
          gnew(ig,iv,iglo) = ( source(ig,iv,iglo) - gnew(ig+1,iv,iglo)*dvpafac(ig,iv) &
               + gnew(ig,iv+1,iglo)*dthetfac(ig,iv,iglo) ) / (1.0 - dvpafac(ig,iv) + dthetfac(ig,iv,iglo))
       end do
    end do
    
    ! boundary condition for vpa < 0 is g(theta=infinity) = 0
    ! boundary condition for theta < 0 (dvpa/dt > 0) is g(vpa=-infinity)=0
    do iv = -nvgrid+1, 0
       do ig = -1, -ntgrid, -1
          gnew(ig,iv,iglo) = ( source(ig,iv,iglo) - gnew(ig+1,iv,iglo)*dvpafac(ig,iv) &
               - gnew(ig,iv-1,iglo)*dthetfac(ig,iv-1,iglo) ) / (1.0 - dvpafac(ig,iv) - dthetfac(ig,iv-1,iglo))
       end do
    end do
    
    ! boundary condition for vpa > 0 is g(theta=-infinity) = 0
    ! boundary condition for theta < 0 (dvpa/dt > 0) is g(vpa=-infinity)=0
    do iv = 1, nvgrid
       do ig = -ntgrid+1, -1
          gnew(ig,iv,iglo) = ( source(ig,iv,iglo) + gnew(ig-1,iv,iglo)*dvpafac(ig-1,iv) &
               - gnew(ig,iv-1,iglo)*dthetfac(ig,iv-1,iglo) ) / (1.0 + dvpafac(ig-1,iv) - dthetfac(ig,iv-1,iglo))
       end do
    end do

    ! will need to remove 'real' operator below later -- MAB
    ! this will require tackling the inversion of a complex matrix gresponse as well
    if (present(dgdgn)) dgdgn = matmul(m_mat,real(gnew(-1,1:nvgrid,iglo)))
    
  end subroutine implicit_sweep
  
  subroutine get_source

    use gs3_input, only: ntgrid, nvgrid, dt
    use gs3_layouts, only: g_lo
    use dist_fn_arrays, only: gold

    implicit none

    integer :: ig, iv, iglo

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do iv = -nvgrid, nvgrid
          do ig = -ntgrid, ntgrid
             source(ig,iv,iglo) = gold(ig,iv,iglo)
          end do
       end do
    end do

  end subroutine get_source

  subroutine finish_dist_fn

    use dist_fn_arrays, only: gold, gnew, g

    implicit none

    if (allocated(gold)) then
       deallocate (gold)
       deallocate (gnew)
       deallocate (g)
       deallocate (source)
       deallocate (gresponse)
       deallocate (m_mat)
    end if

    if (allocated(dvpafac)) then
       deallocate (dvpafac)
       deallocate (dthetfac)
    end if

    if (allocated(M_class)) deallocate (M_class, N_class)

  end subroutine finish_dist_fn

end module dist_fn

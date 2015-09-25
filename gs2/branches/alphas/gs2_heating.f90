module gs2_heating

  implicit none
  private

  public :: heating_diagnostics
  public :: init_htype, zero_htype, del_htype
  public :: avg_h, avg_hk, avg_he
  public :: hk_repack, he_repack
  !GGH Density-velocity perturbations
  public :: dens_vel_diagnostics
  public :: init_dvtype, zero_dvtype, del_dvtype
  public :: avg_dv, avg_dvk
  public :: init_hetype, zero_hetype, del_hetype
   
  interface init_htype
     module procedure init_htype_0, init_htype_1, init_htype_2, init_htype_3
  end interface

  interface zero_htype
     module procedure zero_htype_0, zero_htype_1, zero_htype_2, zero_htype_3
  end interface

  interface del_htype
     module procedure del_htype_0, del_htype_1, del_htype_2, del_htype_3
  end interface

  interface init_hetype
     module procedure init_hetype_1, init_hetype_2
  end interface

  interface zero_hetype
     module procedure zero_hetype_1, zero_hetype_2
  end interface

  interface del_hetype
     module procedure del_hetype_1, del_hetype_2
  end interface

  !GGH
  interface init_dvtype
     module procedure init_dvtype_0, init_dvtype_1, init_dvtype_2, init_dvtype_3
  end interface

  interface zero_dvtype
     module procedure zero_dvtype_0, zero_dvtype_1, zero_dvtype_2, zero_dvtype_3
  end interface

  interface del_dvtype
     module procedure del_dvtype_0, del_dvtype_1, del_dvtype_2, del_dvtype_3
  end interface

  type :: heating_diagnostics
! total quantities:
     real :: energy
     real :: energy_dot
     real :: antenna
     real :: eapar                                       !int k_perp^2 A_par^2/8 pi
     real :: ebpar                                       !int B_par^2/8 pi
! species by species:
     real, dimension(:), pointer :: delfs2     => null() !int T/F0 dfs^2/2
     real, dimension(:), pointer :: hs2        => null() !int T/F0 hs^2/2
     real, dimension(:), pointer :: phis2      => null() !int q^2 n/T phi^2/2
     real, dimension(:), pointer :: hypervisc  => null() 
     real, dimension(:), pointer :: hyperres   => null()
     real, dimension(:), pointer :: hypercoll  => null()
     real, dimension(:), pointer :: collisions => null()
     real, dimension(:), pointer :: imp_colls => null()
     real, dimension(:), pointer :: gradients  => null()
!     real, dimension(:), pointer :: curvature  => null()
     real, dimension(:), pointer :: heating    => null()
     real, dimension(:), pointer :: dh2dt    => null()
     real, dimension(:), pointer :: parstream    => null()
  end type heating_diagnostics

!GGH>
  type :: dens_vel_diagnostics
     !GGH NOTE: Dimension here is for species
     real, dimension(:), pointer :: dvpar  => null()
     real, dimension(:), pointer :: dvperp   => null()
     real, dimension(:), pointer :: dn => null()
  end type dens_vel_diagnostics
!<GGH
contains

  subroutine init_htype_0 (h, nspec)

    type (heating_diagnostics) :: h
    integer, intent (in) :: nspec
    
    allocate (h % delfs2(nspec))      
    allocate (h % hs2(nspec))      
    allocate (h % phis2(nspec))      
    allocate (h % hypervisc(nspec))   
    allocate (h % hyperres(nspec))    
    allocate (h % hypercoll(nspec))  
    allocate (h % collisions(nspec))  
    allocate (h % imp_colls(nspec))  
    allocate (h % gradients(nspec))   
!    allocate (h % curvature(nspec))   
    allocate (h % heating(nspec))     

    call zero_htype (h)

  end subroutine init_htype_0

  subroutine init_htype_1 (h, nspec)

    type (heating_diagnostics), dimension(:) :: h
    integer, intent (in) :: nspec
    integer :: n, nmax

    nmax = size(h)
    
    do n=1,nmax
       allocate (h(n) % delfs2(nspec))      
       allocate (h(n) % hs2(nspec))      
       allocate (h(n) % phis2(nspec))      
       allocate (h(n) % hypervisc(nspec))   
       allocate (h(n) % hyperres(nspec))    
       allocate (h(n) % hypercoll(nspec))  
       allocate (h(n) % collisions(nspec))  
       allocate (h(n) % imp_colls(nspec))  
       allocate (h(n) % gradients(nspec))   
       !    allocate (h(n) % curvature(nspec))   
       allocate (h(n) % heating(nspec))     
    end do

    call zero_htype (h)

  end subroutine init_htype_1

  subroutine init_htype_2 (h, nspec)

    type (heating_diagnostics), dimension(:,:) :: h
    integer, intent (in) :: nspec
    integer :: m, n, mmax, nmax

    mmax = size(h, 1)
    nmax = size(h, 2)

    do n = 1, nmax
       do m = 1, mmax
          allocate (h(m,n) % delfs2(nspec))      
          allocate (h(m,n) % hs2(nspec))      
          allocate (h(m,n) % phis2(nspec))      
          allocate (h(m,n) % hypervisc(nspec))   
          allocate (h(m,n) % hyperres(nspec))    
          allocate (h(m,n) % hypercoll(nspec))  
          allocate (h(m,n) % collisions(nspec))  
          allocate (h(m,n) % imp_colls(nspec))  
          allocate (h(m,n) % gradients(nspec))   
          !    allocate (h(m,n) % curvature(nspec))   
          allocate (h(m,n) % heating(nspec))     
       end do
    end do

    call zero_htype (h)

  end subroutine init_htype_2
    
  subroutine init_htype_3 (h, nspec)

    type (heating_diagnostics), dimension(:,:,:) :: h
    integer, intent (in) :: nspec
    integer :: l, m, n, lmax, mmax, nmax

    lmax = size(h, 1)
    mmax = size(h, 2)
    nmax = size(h, 3)

    do n = 1, nmax
       do m = 1, mmax
          do l = 1, lmax    
             allocate (h(l,m,n) % delfs2(nspec))      
             allocate (h(l,m,n) % hs2(nspec))      
             allocate (h(l,m,n) % phis2(nspec))      
             allocate (h(l,m,n) % hypervisc(nspec))   
             allocate (h(l,m,n) % hyperres(nspec))    
             allocate (h(l,m,n) % hypercoll(nspec))  
             allocate (h(l,m,n) % collisions(nspec))  
             allocate (h(l,m,n) % imp_colls(nspec))  
             allocate (h(l,m,n) % gradients(nspec))   
             !    allocate (h(l,m,n) % curvature(nspec))   
             allocate (h(l,m,n) % heating(nspec))     
          end do
       end do 
    end do

    call zero_htype (h)

  end subroutine init_htype_3
 
  subroutine init_hetype_1 (h, nspec)

    type (heating_diagnostics), dimension(:) :: h
    integer, intent (in) :: nspec
    integer :: n, nmax

    nmax = size(h)
    
    do n=1,nmax
       allocate (h(n) % dh2dt(nspec))      
       allocate (h(n) % hs2(nspec))      
       allocate (h(n) % gradients(nspec))   
       allocate (h(n) % heating(nspec))     
       allocate (h(n) % parstream(nspec))     
    end do

    call zero_hetype (h)

  end subroutine init_hetype_1

  subroutine init_hetype_2 (h, nspec)

    type (heating_diagnostics), dimension(:,:) :: h
    integer, intent (in) :: nspec
    integer :: n, nmax, m, mmax

    mmax = size(h,1)
    nmax = size(h,2)
    
    do n=1,nmax
    do m=1,mmax
       allocate (h(m,n) % dh2dt(nspec))      
       allocate (h(m,n) % hs2(nspec))      
       allocate (h(m,n) % gradients(nspec))   
       allocate (h(m,n) % heating(nspec))     
       allocate (h(m,n) % parstream(nspec))     
    end do
    end do

    call zero_hetype (h)

  end subroutine init_hetype_2

  subroutine zero_htype_0 (h)

    type (heating_diagnostics) :: h
    
    h % energy = 0.
    h % energy_dot = 0.
    h % antenna = 0.
    h % eapar = 0.
    h % ebpar = 0.
    h % delfs2 = 0. 
    h % hs2 = 0. 
    h % phis2 = 0. 
    h % hypervisc = 0. 
    h % hyperres = 0.  
    h % hypercoll = 0.  
    h % collisions = 0.
    h % imp_colls = 0.
    h % gradients = 0. 
!    h % curvature = 0. 
    h % heating = 0.   

  end subroutine zero_htype_0

  subroutine zero_htype_1 (h)

    type (heating_diagnostics), dimension(:) :: h
    integer :: n, nmax

    nmax = size(h)
    
    h % energy = 0.
    h % energy_dot = 0.
    h % antenna = 0.
    h % eapar = 0.
    h % ebpar = 0.
    do n=1,nmax
       h(n) % delfs2 = 0. 
       h(n) % hs2 = 0. 
       h(n) % phis2 = 0. 
       h(n) % hypervisc = 0. 
       h(n) % hyperres = 0.  
       h(n) % hypercoll = 0.  
       h(n) % collisions = 0.
       h(n) % imp_colls = 0.
       h(n) % gradients = 0. 
       !    h(n) % curvature = 0. 
       h(n) % heating = 0.   
    end do

  end subroutine zero_htype_1

  subroutine zero_htype_2 (h)

    type (heating_diagnostics), dimension(:,:) :: h
    integer :: n, nmax, m, mmax

    mmax = size(h, 1)
    nmax = size(h, 2)
    
    h % energy = 0.
    h % energy_dot = 0.
    h % antenna = 0.
    h % eapar = 0.
    h % ebpar = 0.
    do n=1,nmax
       do m=1,mmax
          h(m,n) % delfs2 = 0. 
          h(m,n) % hs2 = 0. 
          h(m,n) % phis2 = 0. 
          h(m,n) % hypervisc = 0. 
          h(m,n) % hyperres = 0.  
          h(m,n) % hypercoll = 0.  
          h(m,n) % collisions = 0.
          h(m,n) % imp_colls = 0.
          h(m,n) % gradients = 0. 
          !    h(m,n) % curvature = 0. 
          h(m,n) % heating = 0.   
       end do
    end do

  end subroutine zero_htype_2

  subroutine zero_htype_3 (h)

    type (heating_diagnostics), dimension(:,:,:) :: h
    integer :: n, nmax, m, mmax, l, lmax

    lmax = size(h, 1)
    mmax = size(h, 2)
    nmax = size(h, 3)
        
    h % energy = 0.
    h % energy_dot = 0.
    h % antenna = 0.
    h % eapar = 0.
    h % ebpar = 0.
    do n=1,nmax
       do m=1,mmax
          do l=1,lmax
             h(l,m,n) % delfs2 = 0. 
             h(l,m,n) % hs2 = 0. 
             h(l,m,n) % phis2 = 0. 
             h(l,m,n) % hypervisc = 0. 
             h(l,m,n) % hyperres = 0.  
             h(l,m,n) % hypercoll = 0.  
             h(l,m,n) % collisions = 0.
             h(l,m,n) % imp_colls = 0.
             h(l,m,n) % gradients = 0. 
             !    h(l,m,n) % curvature = 0. 
             h(l,m,n) % heating = 0.   
          end do
       end do
    end do

  end subroutine zero_htype_3

  subroutine zero_hetype_1 (h)

    type (heating_diagnostics), dimension(:) :: h
    integer :: n, nmax

    nmax = size(h)
    
    do n=1,nmax
       h(n) % dh2dt = 0. 
       h(n) % hs2 = 0. 
       h(n) % gradients = 0. 
       h(n) % heating = 0.   
       h(n) % parstream = 0.   
    end do

  end subroutine zero_hetype_1


  subroutine zero_hetype_2 (h)

    type (heating_diagnostics), dimension(:,:) :: h
    integer :: n, nmax, m, mmax

    mmax = size(h,1)
    nmax = size(h,2)
    
    do n=1,nmax
    do m=1,mmax
       h(m,n) % dh2dt = 0. 
       h(m,n) % hs2 = 0. 
       h(m,n) % gradients = 0. 
       h(m,n) % heating = 0.   
       h(m,n) % parstream = 0.   
    end do
    end do

  end subroutine zero_hetype_2

  subroutine del_htype_0 (h)

    type (heating_diagnostics) :: h
    
    deallocate (h % delfs2)
    deallocate (h % hs2)
    deallocate (h % phis2)
    deallocate (h % hypervisc)
    deallocate (h % hyperres)
    deallocate (h % hypercoll)
    deallocate (h % collisions)
    deallocate (h % imp_colls)
    deallocate (h % gradients)
!    deallocate (h % curvature)
    deallocate (h % heating)

  end subroutine del_htype_0

  subroutine del_htype_1 (h)

    type (heating_diagnostics), dimension(:) :: h
    integer :: m, mmax
    
    mmax = size (h, 1)
    
    do m=1,mmax
       deallocate (h(m) % delfs2)
       deallocate (h(m) % hs2)
       deallocate (h(m) % phis2)
       deallocate (h(m) % hypervisc)
       deallocate (h(m) % hyperres)
       deallocate (h(m) % hypercoll)
       deallocate (h(m) % collisions)
       deallocate (h(m) % imp_colls)
       deallocate (h(m) % gradients)
       !       deallocate (h(m,n) % curvature)
       deallocate (h(m) % heating)
    end do

  end subroutine del_htype_1

  subroutine del_htype_2 (h)

    type (heating_diagnostics), dimension(:,:) :: h
    integer :: m, mmax, n, nmax
    
    mmax = size (h, 1)
    nmax = size (h, 2)
    
    do n=1,nmax
       do m=1,mmax
          deallocate (h(m,n) % delfs2)
          deallocate (h(m,n) % hs2)
          deallocate (h(m,n) % phis2)
          deallocate (h(m,n) % hypervisc)
          deallocate (h(m,n) % hyperres)
          deallocate (h(m,n) % hypercoll)
          deallocate (h(m,n) % collisions)
          deallocate (h(m,n) % imp_colls)
          deallocate (h(m,n) % gradients)
          !    deallocate (h(m,n) % curvature)
          deallocate (h(m,n) % heating)
       end do
    end do

  end subroutine del_htype_2

  subroutine del_htype_3 (h)

    type (heating_diagnostics), dimension(:,:,:) :: h
    integer :: l, m, n, lmax, mmax, nmax

    lmax = size(h, 1)
    mmax = size(h, 2)
    nmax = size(h, 3)

    do n = 1, nmax
       do m = 1, mmax
          do l = 1, lmax    
             deallocate (h(l,m,n) % delfs2)      
             deallocate (h(l,m,n) % hs2)      
             deallocate (h(l,m,n) % phis2)      
             deallocate (h(l,m,n) % hypervisc)   
             deallocate (h(l,m,n) % hyperres)    
             deallocate (h(l,m,n) % hypercoll)  
             deallocate (h(l,m,n) % collisions)  
             deallocate (h(l,m,n) % imp_colls)  
             deallocate (h(l,m,n) % gradients)   
             !    deallocate (h(l,m,n) % curvature)   
             deallocate (h(l,m,n) % heating)     
          end do
       end do 
    end do

  end subroutine del_htype_3

  subroutine del_hetype_1 (h)

    type (heating_diagnostics), dimension(:) :: h
    integer :: m, mmax
    
    mmax = size (h, 1)
    
    do m=1,mmax
       deallocate (h(m) % dh2dt)
       deallocate (h(m) % hs2)
       deallocate (h(m) % gradients)
       deallocate (h(m) % heating)
       deallocate (h(m) % parstream)
    end do

  end subroutine del_hetype_1

  subroutine del_hetype_2 (h)

    type (heating_diagnostics), dimension(:,:) :: h
    integer :: m, mmax, n, nmax
    
    mmax = size (h, 1)
    nmax = size (h, 2)
    
    do n=1,nmax
    do m=1,mmax
       deallocate (h(m,n) % dh2dt)
       deallocate (h(m,n) % hs2)
       deallocate (h(m,n) % gradients)
       deallocate (h(m,n) % heating)
       deallocate (h(m,n) % parstream)
    end do
    end do

  end subroutine del_hetype_2

  subroutine avg_h (h, h_hist, istep, navg)

    use mp, only: proc0
    use species, only: nspec
    type (heating_diagnostics) :: h
    type (heating_diagnostics), dimension (0:) :: h_hist
    integer, intent (in) :: istep, navg
    integer :: is, i

    if (proc0) then
       if (navg > 1) then
          if (istep > 1) then
             h_hist(mod(istep,navg)) % energy     = h % energy
             h_hist(mod(istep,navg)) % energy_dot = h % energy_dot
             h_hist(mod(istep,navg)) % antenna    = h % antenna
             h_hist(mod(istep,navg)) % eapar      = h % eapar
             h_hist(mod(istep,navg)) % ebpar      = h % ebpar
             h_hist(mod(istep,navg)) % delfs2     = h % delfs2
             h_hist(mod(istep,navg)) % hs2        = h % hs2
             h_hist(mod(istep,navg)) % phis2      = h % phis2
             h_hist(mod(istep,navg)) % hypervisc  = h % hypervisc
             h_hist(mod(istep,navg)) % hyperres   = h % hyperres
             h_hist(mod(istep,navg)) % hypercoll  = h % hypercoll
             h_hist(mod(istep,navg)) % collisions = h % collisions
             h_hist(mod(istep,navg)) % imp_colls  = h % imp_colls
             h_hist(mod(istep,navg)) % gradients  = h % gradients
!             h_hist(mod(istep,navg)) % curvature  = h % curvature
             h_hist(mod(istep,navg)) % heating    = h % heating
          end if
          
          if (istep >= navg) then
             call zero_htype(h)
             do i=0,navg-1
                h % energy     = h%energy     + h_hist(i) % energy / real(navg)
                h % energy_dot = h%energy_dot + h_hist(i) % energy_dot / real(navg)
                h % antenna    = h%antenna    + h_hist(i) % antenna / real(navg)
                h % eapar      = h%eapar      + h_hist(i) % eapar    / real(navg)
                h % ebpar      = h%ebpar      + h_hist(i) % ebpar    / real(navg)
                do is = 1,nspec
                   h % delfs2(is)     = h%delfs2(is)     + h_hist(i) % delfs2(is) / real(navg)
                   h % hs2(is)        = h%hs2(is)        + h_hist(i) % hs2(is) / real(navg)
                   h % phis2(is)      = h%phis2(is)      + h_hist(i) % phis2(is) / real(navg)
                   h % hypervisc(is)  = h%hypervisc(is)  + h_hist(i) % hypervisc(is) / real(navg)
                   h % hyperres(is)   = h%hyperres(is)   + h_hist(i) % hyperres(is) / real(navg)
                   h % hypercoll(is)  = h%hypercoll(is)  + h_hist(i) % hypercoll(is) / real(navg)
                   h % collisions(is) = h%collisions(is) + h_hist(i) % collisions(is) / real(navg)
                   h % imp_colls(is) = h%imp_colls(is) + h_hist(i) % imp_colls(is) / real(navg)
                   h % gradients(is)  = h%gradients(is)  + h_hist(i) % gradients(is) / real(navg)
!                  h % curvature(is)  = h%curvature(is)  + h_hist(i) % curvature(is) / real(navg)
                   h % heating(is)    = h%heating(is)    + h_hist(i) % heating(is) / real(navg)
                end do
             end do
          end if
       end if
    end if

  end subroutine avg_h

  subroutine avg_hk (hk, hk_hist, istep, navg)

    use mp, only: proc0
    use species, only: nspec
    type (heating_diagnostics), dimension(:,:) :: hk
    type (heating_diagnostics), dimension (:,:,0:) :: hk_hist
    integer, intent (in) :: istep, navg
    integer :: is, i, m, n, mmax, nmax

    mmax = size(hk, 1)
    nmax = size(hk, 2)

    if (proc0) then
       if (navg > 1) then
          if (istep > 1) then
             hk_hist(:,:,mod(istep,navg)) % energy     = hk % energy
             hk_hist(:,:,mod(istep,navg)) % energy_dot = hk % energy_dot
             hk_hist(:,:,mod(istep,navg)) % antenna    = hk % antenna
             hk_hist(:,:,mod(istep,navg)) % eapar      = hk % eapar
             hk_hist(:,:,mod(istep,navg)) % ebpar      = hk % ebpar
            do n=1,nmax
                do m=1,mmax
                   do is = 1,nspec
                      hk_hist(m,n,mod(istep,navg)) % delfs2(is)     = hk(m,n) % delfs2(is)
                      hk_hist(m,n,mod(istep,navg)) % hs2(is)        = hk(m,n) % hs2(is)
                      hk_hist(m,n,mod(istep,navg)) % phis2(is)      = hk(m,n) % phis2(is)
                      hk_hist(m,n,mod(istep,navg)) % hypervisc(is)  = hk(m,n) % hypervisc(is)
                      hk_hist(m,n,mod(istep,navg)) % hyperres(is)   = hk(m,n) % hyperres(is)
                      hk_hist(m,n,mod(istep,navg)) % hypercoll(is)  = hk(m,n) % hypercoll(is)
                      hk_hist(m,n,mod(istep,navg)) % collisions(is) = hk(m,n) % collisions(is)
                      hk_hist(m,n,mod(istep,navg)) % imp_colls(is)  = hk(m,n) % imp_colls(is)
                      hk_hist(m,n,mod(istep,navg)) % gradients(is)  = hk(m,n) % gradients(is)
!                     hk_hist(m,n,mod(istep,navg)) % curvature(is)  = hk(m,n) % curvature(is)
                      hk_hist(m,n,mod(istep,navg)) % heating(is)    = hk(m,n) % heating(is)
                   end do
                end do
             end do
          end if
          
          if (istep >= navg) then
             call zero_htype (hk)
             do n=1,nmax
                do m=1,mmax
                   do i=0,navg-1
                      hk(m,n) % energy     = hk(m,n)%energy     + hk_hist(m,n,i) % energy / real(navg)
                      hk(m,n) % energy_dot = hk(m,n)%energy_dot + hk_hist(m,n,i) % energy_dot / real(navg)
                      hk(m,n) % antenna    = hk(m,n)%antenna    + hk_hist(m,n,i) % antenna / real(navg)
                      hk(m,n) % eapar      = hk(m,n)%eapar      + hk_hist(m,n,i) % eapar / real(navg)
                      hk(m,n) % ebpar      = hk(m,n)%ebpar      + hk_hist(m,n,i) % ebpar / real(navg)
                      do is = 1,nspec
                         hk(m,n)%delfs2(is)     = hk(m,n)%delfs2(is)    + hk_hist(m,n,i) % delfs2(is) / real(navg)
                         hk(m,n)%hs2(is)        = hk(m,n)%hs2(is)       + hk_hist(m,n,i) % hs2(is) / real(navg)
                         hk(m,n)%phis2(is)      = hk(m,n)%phis2(is)     + hk_hist(m,n,i) % phis2(is) / real(navg)
                         hk(m,n)%hypervisc(is)  = hk(m,n)%hypervisc(is) + hk_hist(m,n,i) % hypervisc(is) / real(navg)
                         hk(m,n)%hyperres(is)   = hk(m,n)%hyperres(is)  + hk_hist(m,n,i) % hyperres(is) / real(navg)
                         hk(m,n)%hypercoll(is)  = hk(m,n)%hypercoll(is) + hk_hist(m,n,i) % hypercoll(is) / real(navg)
                         hk(m,n)%collisions(is) = hk(m,n)%collisions(is)+ hk_hist(m,n,i) % collisions(is) / real(navg)
                         hk(m,n)%imp_colls(is)  = hk(m,n)%imp_colls(is) + hk_hist(m,n,i) % imp_colls(is) / real(navg)
                         hk(m,n)%gradients(is)  = hk(m,n)%gradients(is) + hk_hist(m,n,i) % gradients(is) / real(navg)
!                        hk(m,n)%curvature(is)  = hk(m,n)%curvature(is) + hk_hist(m,n,i) % curvature(is) / real(navg)
                         hk(m,n)%heating(is)    = hk(m,n)%heating(is)   + hk_hist(m,n,i) % heating(is) / real(navg)
                      end do
                   end do
                end do
             end do
          end if
       end if
    end if

  end subroutine avg_hk

  subroutine avg_he (he, he_hist, istep, navg)

    use mp, only: proc0
    use species, only: nspec
    type (heating_diagnostics), dimension(:) :: he
    type (heating_diagnostics), dimension (:,0:) :: he_hist
    integer, intent (in) :: istep, navg
    integer :: is, i, m, n, mmax, nmax

    nmax = size(he)

    if (proc0) then
       if (navg > 1) then
          if (istep > 1) then
            do n=1,nmax
                do is = 1,nspec
                   he_hist(n,mod(istep,navg)) % hs2(is)        = he(n) % hs2(is)
                   he_hist(n,mod(istep,navg)) % gradients(is)  = he(n) % gradients(is)
                   he_hist(n,mod(istep,navg)) % heating(is)    = he(n) % heating(is)
                   he_hist(n,mod(istep,navg)) % dh2dt(is)    = he(n) % dh2dt(is)
                   he_hist(n,mod(istep,navg)) % parstream(is)    = he(n) % parstream(is)
                end do
             end do
          end if
          
          if (istep >= navg) then
             call zero_hetype (he)
             do n=1,nmax
                do i=0,navg-1
                   do is = 1,nspec
                      he(n)%hs2(is)        = he(n)%hs2(is)       + he_hist(n,i) % hs2(is) / real(navg)
                      he(n)%gradients(is)  = he(n)%gradients(is) + he_hist(n,i) % gradients(is) / real(navg)
                      he(n)%heating(is)    = he(n)%heating(is)   + he_hist(n,i) % heating(is) / real(navg)
                      he(n)%dh2dt(is)    = he(n)%dh2dt(is)   + he_hist(n,i) % dh2dt(is) / real(navg)
                      he(n)%parstream(is)    = he(n)%parstream(is)   + he_hist(n,i) % parstream(is) / real(navg)
                   end do
                end do
             end do
          end if
       end if
    end if

  end subroutine avg_he

  subroutine hk_repack (hk, i, tmp)

    use species, only: nspec
    type (heating_diagnostics), dimension (:,:), intent(in) :: hk
    real, dimension(:,:,:), intent (out) :: tmp
    integer, intent (in) :: i
    integer :: is, m, n, mmax, nmax

    mmax = size(tmp, 1)
    nmax = size(tmp, 2)

    select case (i)
    case (1) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%hypervisc(is)
             end do
          end do
       end do
    case (2) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%hyperres(is)
             end do
          end do
       end do
    case (3) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%collisions(is)
             end do
          end do
       end do
    case (4) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%gradients(is)
             end do
          end do
       end do
    case (5) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%heating(is)
             end do
          end do
       end do
    case (6) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%hypercoll(is)
             end do
          end do
       end do
    case (7) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%delfs2(is)
             end do
          end do
       end do
    case (8) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%hs2(is)
             end do
          end do
       end do
    case (9) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%phis2(is)
             end do
          end do
       end do
    case (10) 
       do is=1,nspec
          do n=1,nmax
             do m=1,mmax
                tmp(m,n,is) = hk(m,n)%imp_colls(is)
             end do
          end do
       end do
    end select
    
  end subroutine hk_repack


  subroutine he_repack (he, i, tmp)

    use species, only: nspec
    type (heating_diagnostics), dimension (:), intent(in) :: he
    real, dimension(:,:), intent (out) :: tmp
    integer, intent (in) :: i
    integer :: is, m, n, mmax, nmax

    nmax = size(tmp, 1)

    select case (i)
    case (1) 
       do is=1,nspec
          do n=1,nmax
             tmp(n,is) = he(n)%gradients(is)
          end do
       end do
    case (2) 
       do is=1,nspec
          do n=1,nmax
             tmp(n,is) = he(n)%heating(is)
          end do
       end do
    case (3) 
       do is=1,nspec
          do n=1,nmax
             tmp(n,is) = he(n)%hs2(is)
          end do
       end do
    case (4) 
       do is=1,nspec
          do n=1,nmax
             tmp(n,is) = he(n)%dh2dt(is)
          end do
       end do
    case (5) 
       do is=1,nspec
          do n=1,nmax
             tmp(n,is) = he(n)%parstream(is)
          end do
       end do

    end select
    
    
  end subroutine he_repack

!=============================================================================
!<GGH
! Density velocity perturbation routines
!=============================================================================
    subroutine init_dvtype_0 (dv, nspec)

    type (dens_vel_diagnostics) :: dv
    integer, intent (in) :: nspec

       allocate (dv % dvpar(nspec))
       allocate (dv % dvperp(nspec))    
       allocate (dv % dn(nspec))  

    call zero_dvtype (dv)

  end subroutine init_dvtype_0
!=============================================================================
  subroutine init_dvtype_1 (dv, nspec)
    type (dens_vel_diagnostics), dimension(:) :: dv
    integer, intent (in) :: nspec
    integer :: n, nmax

    nmax = size(dv)
    
    do n=1,nmax
       allocate (dv(n) % dvpar(nspec))
       allocate (dv(n) % dvperp(nspec))    
       allocate (dv(n) % dn(nspec))  
    end do

    call zero_dvtype (dv)

  end subroutine init_dvtype_1
!=============================================================================
  subroutine init_dvtype_2 (dv, nspec)
    type (dens_vel_diagnostics), dimension(:,:) :: dv
    integer, intent (in) :: nspec
    integer :: m, n, mmax, nmax

    mmax = size(dv, 1)
    nmax = size(dv, 2)

    do n = 1, nmax
       do m = 1, mmax
          allocate (dv(m,n) % dvpar(nspec))
          allocate (dv(m,n) % dvperp(nspec))    
          allocate (dv(m,n) % dn(nspec))  
       end do
    end do

    call zero_dvtype (dv)

  end subroutine init_dvtype_2
!=============================================================================
  subroutine init_dvtype_3 (dv, nspec)
    type (dens_vel_diagnostics), dimension(:,:,:) :: dv
    integer, intent (in) :: nspec
    integer :: l, m, n, lmax, mmax, nmax

    lmax = size(dv, 1)
    mmax = size(dv, 2)
    nmax = size(dv, 3)

    do n = 1, nmax
       do m = 1, mmax
          do l = 1, lmax    
             allocate (dv(l,m,n) % dvpar(nspec))
             allocate (dv(l,m,n) % dvperp(nspec))    
             allocate (dv(l,m,n) % dn(nspec))  
          end do
       end do 
    end do

    call zero_dvtype (dv)

  end subroutine init_dvtype_3
!=============================================================================
  subroutine zero_dvtype_0 (dv)

    type (dens_vel_diagnostics) :: dv
    integer :: n, nmax

    dv % dvpar = 0. 
    dv % dvperp = 0.  
    dv % dn = 0.

  end subroutine zero_dvtype_0
!=============================================================================
  subroutine zero_dvtype_1 (dv)

    type (dens_vel_diagnostics), dimension(:) :: dv
    integer :: n, nmax

    nmax = size(dv)
    
    do n=1,nmax
       dv(n) % dvpar = 0. 
       dv(n) % dvperp = 0.  
       dv(n) % dn = 0.
    end do

  end subroutine zero_dvtype_1
!=============================================================================
 subroutine zero_dvtype_2 (dv)

    type (dens_vel_diagnostics), dimension(:,:) :: dv
    integer :: n, nmax, m, mmax

    mmax = size(dv, 1)
    nmax = size(dv, 2)
        
    do n=1,nmax
       do m=1,mmax
          dv(m,n) % dvpar = 0. 
          dv(m,n) % dvperp = 0.  
          dv(m,n) % dn = 0.
       end do
    end do

  end subroutine zero_dvtype_2
!=============================================================================
 subroutine zero_dvtype_3 (dv)

    type (dens_vel_diagnostics), dimension(:,:,:) :: dv
    integer :: n, nmax, m, mmax, l, lmax

    lmax = size(dv, 1)
    mmax = size(dv, 2)
    nmax = size(dv, 3)
        
    do n=1,nmax
       do m=1,mmax
          do l=1,lmax
             dv(l,m,n) % dvpar = 0. 
             dv(l,m,n) % dvperp = 0.  
             dv(l,m,n) % dn = 0.
          end do
       end do
    end do

  end subroutine zero_dvtype_3
!=============================================================================
  subroutine del_dvtype_0 (dv)
    type (dens_vel_diagnostics) :: dv
    
    deallocate (dv % dvpar)
    deallocate (dv % dvperp)
    deallocate (dv % dn)

  end subroutine del_dvtype_0
!=============================================================================
  subroutine del_dvtype_1 (dv)
    type (dens_vel_diagnostics), dimension(:) :: dv
    integer :: m, mmax

    mmax = size (dv, 1)

    do m=1,mmax
       deallocate (dv(m) % dvpar)
       deallocate (dv(m) % dvperp)
       deallocate (dv(m) % dn)
    end do

  end subroutine del_dvtype_1
!=============================================================================
 subroutine del_dvtype_2 (dv)

    type (dens_vel_diagnostics), dimension(:,:) :: dv
    integer :: m, mmax, n, nmax
    
    mmax = size (dv, 1)
    nmax = size (dv, 2)
    
    do n=1,nmax
       do m=1,mmax
          deallocate (dv(m,n) % dvpar)
          deallocate (dv(m,n) % dvperp)
          deallocate (dv(m,n) % dn)
       end do
    end do

  end subroutine del_dvtype_2
!=============================================================================
 subroutine del_dvtype_3 (dv)

    type (dens_vel_diagnostics), dimension(:,:,:) :: dv
    integer :: l, lmax, m, mmax, n, nmax
    
    lmax = size (dv, 1)
    mmax = size (dv, 2)
    nmax = size (dv, 3)
    
    do l=1,lmax
       do n=1,nmax
          do m=1,mmax
             deallocate (dv(l,m,n) % dvpar)
             deallocate (dv(l,m,n) % dvperp)
             deallocate (dv(l,m,n) % dn)
          end do
       end do
    end do

  end subroutine del_dvtype_3
!=============================================================================
  subroutine avg_dv (dv, dv_hist, istep, navg)
    use mp, only: proc0
    use species, only: nspec
    type (dens_vel_diagnostics) :: dv
    type (dens_vel_diagnostics), dimension (0:) :: dv_hist
    integer, intent (in) :: istep, navg
    integer :: is, i

    if (proc0) then
       if (navg > 1) then
          if (istep > 1) then
             dv_hist(mod(istep,navg)) % dvpar(:)     = dv % dvpar(:)
             dv_hist(mod(istep,navg)) % dvperp(:)    = dv % dvperp(:)
             dv_hist(mod(istep,navg)) % dn(:)        = dv % dn(:)
          end if
          
          if (istep >= navg) then
             call zero_dvtype(dv)
             do i=0,navg-1
                dv % dvpar(:)  = dv % dvpar(:)  + dv_hist(i) % dvpar(:)  / real(navg)
                dv % dvperp(:) = dv % dvperp(:) + dv_hist(i) % dvperp(:) / real(navg)
                dv % dn(:)     = dv % dn(:)     + dv_hist(i) % dn (:)    / real(navg)
             end do
          end if
       end if
    end if

  end subroutine avg_dv
!=============================================================================
  subroutine avg_dvk (dvk, dvk_hist, istep, navg)
    use mp, only: proc0
    use species, only: nspec
    type (dens_vel_diagnostics), dimension(:,:) :: dvk
    type (dens_vel_diagnostics), dimension (:,:,0:) :: dvk_hist
    integer, intent (in) :: istep, navg
    integer :: is, i, m, n, mmax, nmax

    mmax = size(dvk, 1)
    nmax = size(dvk, 2)

    if (proc0) then
       if (navg > 1) then
          if (istep > 1) then
             do n=1,nmax
                do m=1,mmax
                   dvk_hist(m,n,mod(istep,navg)) % dvpar(:)     = dvk(m,n) % dvpar(:)
                   dvk_hist(m,n,mod(istep,navg)) % dvperp(:)    = dvk(m,n) % dvperp(:)
                   dvk_hist(m,n,mod(istep,navg)) % dn(:)        = dvk(m,n) % dn(:)
             enddo
          enddo
          end if
          
          if (istep >= navg) then
             call zero_dvtype (dvk)
             do n=1,nmax
                do m=1,mmax
                   do i=0,navg-1
                      dvk(m,n) % dvpar(:)  = dvk(m,n) % dvpar(:)  + dvk_hist(m,n,i) % dvpar(:)  / real(navg)
                      dvk(m,n) % dvperp(:) = dvk(m,n) % dvperp(:) + dvk_hist(m,n,i) % dvperp(:) / real(navg)
                      dvk(m,n) % dn(:)     = dvk(m,n) % dn(:)     + dvk_hist(m,n,i) % dn (:)    / real(navg)
                   end do
                enddo
             enddo
          end if
       end if
    end if

  end subroutine avg_dvk
!>GGH
!=============================================================================

end module gs2_heating

module gs2_heating

  implicit none
  private

  public :: heating_diagnostics
  public :: init_htype, zero_htype, del_htype
  public :: avg_h, avg_hk
  public :: hk_repack
  
  interface init_htype
     module procedure init_htype_0, init_htype_1, init_htype_2, init_htype_3
  end interface

  interface zero_htype
     module procedure zero_htype_0, zero_htype_1, zero_htype_2, zero_htype_3
  end interface

  interface del_htype
     module procedure del_htype_0, del_htype_2
  end interface

  type :: heating_diagnostics
! total quantities:
     real :: energy
     real :: energy_dot
     real :: antenna
! species by species:
     real, dimension(:), pointer :: hypervisc  => null()
     real, dimension(:), pointer :: hyperres   => null()
     real, dimension(:), pointer :: collisions => null()
     real, dimension(:), pointer :: gradients  => null()
!     real, dimension(:), pointer :: curvature  => null()
     real, dimension(:), pointer :: heating    => null()
  end type heating_diagnostics

contains

  subroutine init_htype_0 (h, nspec)

    type (heating_diagnostics) :: h
    integer, intent (in) :: nspec
    
    allocate (h % hypervisc(nspec))   
    allocate (h % hyperres(nspec))    
    allocate (h % collisions(nspec))  
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
       allocate (h(n) % hypervisc(nspec))   
       allocate (h(n) % hyperres(nspec))    
       allocate (h(n) % collisions(nspec))  
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
          allocate (h(m,n) % hypervisc(nspec))   
          allocate (h(m,n) % hyperres(nspec))    
          allocate (h(m,n) % collisions(nspec))  
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
             allocate (h(l,m,n) % hypervisc(nspec))   
             allocate (h(l,m,n) % hyperres(nspec))    
             allocate (h(l,m,n) % collisions(nspec))  
             allocate (h(l,m,n) % gradients(nspec))   
             !    allocate (h(l,m,n) % curvature(nspec))   
             allocate (h(l,m,n) % heating(nspec))     
          end do
       end do 
    end do

    call zero_htype (h)

  end subroutine init_htype_3
    
  subroutine zero_htype_0 (h)

    type (heating_diagnostics) :: h
    
    h % energy = 0.
    h % energy_dot = 0.
    h % antenna = 0.
    h % hypervisc = 0. 
    h % hyperres = 0.  
    h % collisions = 0.
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
    do n=1,nmax
       h(n) % hypervisc = 0. 
       h(n) % hyperres = 0.  
       h(n) % collisions = 0.
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
    do n=1,nmax
       do m=1,mmax
          h(m,n) % hypervisc = 0. 
          h(m,n) % hyperres = 0.  
          h(m,n) % collisions = 0.
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
    do n=1,nmax
       do m=1,mmax
          do l=1,lmax
             h(l,m,n) % hypervisc = 0. 
             h(l,m,n) % hyperres = 0.  
             h(l,m,n) % collisions = 0.
             h(l,m,n) % gradients = 0. 
             !    h(l,m,n) % curvature = 0. 
             h(l,m,n) % heating = 0.   
          end do
       end do
    end do

  end subroutine zero_htype_3

  subroutine del_htype_0 (h)

    type (heating_diagnostics) :: h
    
    deallocate (h % hypervisc)
    deallocate (h % hyperres)
    deallocate (h % collisions)
    deallocate (h % gradients)
!    deallocate (h % curvature)
    deallocate (h % heating)

  end subroutine del_htype_0

  subroutine del_htype_2 (h)

    type (heating_diagnostics), dimension(:,:) :: h
    integer :: m, mmax, n, nmax
    
    m = size (h, 1)
    n = size (h, 2)
    
    do n=1,nmax
       do m=1,mmax
          deallocate (h(m,n) % hypervisc)
          deallocate (h(m,n) % hyperres)
          deallocate (h(m,n) % collisions)
          deallocate (h(m,n) % gradients)
          !    deallocate (h(m,n) % curvature)
          deallocate (h(m,n) % heating)
       end do
    end do

  end subroutine del_htype_2

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
             h_hist(mod(istep,navg)) % hypervisc  = h % hypervisc
             h_hist(mod(istep,navg)) % hyperres   = h % hyperres
             h_hist(mod(istep,navg)) % collisions = h % collisions
             h_hist(mod(istep,navg)) % gradients  = h % gradients
!             h_hist(mod(istep,navg)) % curvature  = h % curvature
             h_hist(mod(istep,navg)) % heating    = h % heating
          end if
          
          if (istep >= navg) then
             call zero_htype(h)
             do i=0,navg-1
                h % energy     = h%energy     + h_hist(i) % energy / real(navg)
                h % energy_dot = h%energy_dot + h_hist(i) % energy_dot / real(navg)
                h % antenna     = h%antenna    + h_hist(i) % antenna / real(navg)
                do is = 1,nspec
                   h % hypervisc(is)  = h%hypervisc(is)  + h_hist(i) % hypervisc(is) / real(navg)
                   h % hyperres(is)   = h%hyperres(is)   + h_hist(i) % hyperres(is) / real(navg)
                   h % collisions(is) = h%collisions(is) + h_hist(i) % collisions(is) / real(navg)
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
             do n=1,nmax
                do m=1,mmax
                   do is = 1,nspec
                      hk_hist(m,n,mod(istep,navg)) % hypervisc(is)  = hk(m,n) % hypervisc(is)
                      hk_hist(m,n,mod(istep,navg)) % hyperres(is)   = hk(m,n) % hyperres(is)
                      hk_hist(m,n,mod(istep,navg)) % collisions(is) = hk(m,n) % collisions(is)
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
                      do is = 1,nspec
                         hk(m,n)%hypervisc(is)  = hk(m,n)%hypervisc(is) + hk_hist(m,n,i) % hypervisc(is) / real(navg)
                         hk(m,n)%hyperres(is)   = hk(m,n)%hyperres(is)  + hk_hist(m,n,i) % hyperres(is) / real(navg)
                         hk(m,n)%collisions(is) = hk(m,n)%collisions(is)+ hk_hist(m,n,i) % collisions(is) / real(navg)
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
    end select
    
  end subroutine hk_repack

end module gs2_heating

!> \\file grid.f03
module grid
  
  use fft

  implicit none
  
  double precision :: pi = 4*atan(1.0)

  contains

  ! generate the 1D array that defines the 3D real space
  function space_grid(N,dr) result(r)
    
    integer, intent(in) :: N
    double precision, intent(in) :: dr
    double precision, allocatable :: r(:)
   
    ! local variables
    integer :: i
    
    allocate(r(N))

    ! spatial array
    do i = 1,N
      r(i) = (-N/2.0 + i - 1.0)*dr
    end do

  end function space_grid
  
  ! generate the 1D array that defines the 3D momentum space
  function mom_grid(Nr,dr) result(kr)
    
    integer, intent(in) :: Nr

    double precision, intent(in) :: dr
    
    double precision, allocatable :: kr(:)

    ! local variables
    integer :: i

    double precision :: Lr
    double precision :: dkr

    double precision, allocatable :: kr_low(:), kr_high(:)

    ! box size
    Lr = Nr*dr

    ! momentum space step
    dkr = 2.0*pi/Lr

    allocate(kr_low(Nr/2))
    allocate(kr_high(Nr/2))
    allocate(kr(Nr))

    ! spectral array defined as [0,...,(Nr-1)*(dkr/2),-Nr*(dkr/2),...,-dkr2]
    do i = 1,Nr/2
      kr_low(i) = (i-1.0)*dkr
      kr_high(i) = (-Nr/2.0 + i - 1.0)*dkr
    end do
    kr(1:Nr/2) = kr_low
    kr(Nr/2+1:Nr) = kr_high

    deallocate(kr_low)
    deallocate(kr_high)

  end function mom_grid

  ! generate the 3D array for the exponential form of the laplacian operator
  function exp_lap(kx,ky,kz,dt)

    implicit none

    integer :: Nx, Ny, Nz
    integer :: i, j, k

    double complex :: dt

    double precision :: kx(:), ky(:), kz(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: exp_lap(:,:,:)

    Nx = size(kx)
    Ny = size(ky)
    Nz = size(kz)
    
    allocate(exp_lap(Nx,Ny,Nz))

    ! laplacian operator defined in the ssfm form
    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          exp_lap(i,j,k) = exp(-0.5*dt*(kx(i)**2 + ky(j)**2 + kz(k)**2))
        end do
      end do
    end do

  end function exp_lap
end module grid

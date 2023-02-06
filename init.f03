!> \\file init.f03
module init
  
  use HDF5
  use fft
  use grid

  implicit none
  
  contains

  ! generate the initial form of the wavefunction 
  function init_wav(x,y,z,init_type,gauss_sig,x_shift)

    double precision, intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: init_type    
    double precision, intent(in) :: gauss_sig
    double precision, intent(in) :: x_shift


    complex(C_DOUBLE_COMPLEX), allocatable :: init_wav(:,:,:)
   
    ! Local variables 
    integer :: Nx, Ny, Nz
    integer :: i, j, k
   
    Nx = size(x)
    Ny = size(y)
    Nz = size(z)

    allocate(init_wav(Nx,Ny,Nz))
    
    ! define initial wavefunction profile as a Gaussian
    if (init_type == 1) then
      write(*,*) "initial input to imaginary time: Gaussian"
      do k = 1, Nz
        do j = 1, Ny
          do i = 1, Nx
            init_wav(i,j,k) = exp(-((x(i)-x_shift)**2.0/(2.0*gauss_sig**2.0) &
                                  + y(j)**2.0/(2.0*gauss_sig**2.0) &
                                  + z(k)**2.0/(2.0*gauss_sig**2.0)))**0.5
          end do
        end do
      end do
    ! define initial wavefunction profile as a Super-Gaussian
    elseif (init_type == 2) then
      write(*,*) "initial input to imaginary time: Super-Gaussian"
      do k = 1, Nz
        do j = 1, Ny
          do i = 1, Nx
            init_wav(i,j,k) = exp(-((x(i)-x_shift)**2.0/(2.0*gauss_sig**2.0) &
                                  + y(j)**2.0/(2.0*gauss_sig**2.0) &
                                  + z(k)**2.0/(2.0*gauss_sig**2.0))**3.0)**0.5
          end do
        end do
      end do
    elseif (init_type == 3) then
      write(*,*) "initial input to imaginary time: Super-Gaussian + const."
      do k = 1, Nz
        do j = 1, Ny
          do i = 1, Nx
            init_wav(i,j,k) = exp(-(x(i)**2.0/(2.0*gauss_sig**2.0) &
                                  + y(j)**2.0/(2.0*gauss_sig**2.0) &
                                  + z(k)**2.0/(2.0*gauss_sig**2.0))**3.0)**0.5 &
                              + 0.001
          end do
        end do
      end do
    end if

  end function init_wav

  function init_pot(x, y, z, omgx, omgy, omgz, x_shift)

    double precision, intent(in) :: x(:), y(:), z(:)
    double precision, intent(in) :: omgx, omgy, omgz
    double precision, intent(in) :: x_shift

    double precision, allocatable :: init_pot(:, :, :)
   
    ! Local variables 
    integer :: Nx, Ny, Nz
    integer :: i, j, k
   
    Nx = size(x)
    Ny = size(y)
    Nz = size(z)

    write(*, *) "trapping potential with shift:", x_shift

    allocate(init_pot(Nx, Ny, Nz))
    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          init_pot(i,j,k) = 0.5*omgx**2.0*(x(i)-x_shift)**2.0 &
                          + 0.5*omgy**2.0*y(j)**2.0 &
                          + 0.5*omgz**2.0*z(k)**2.0 
        end do
      end do
    end do
  end function init_pot


  function readin_wav(x,y,z,comp)
    
    double precision, intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: comp
    complex(C_DOUBLE_COMPLEX), allocatable :: readin_wav(:,:,:)
    
    ! Local variables 
    type(C_PTR) :: f_ptr
    double precision, allocatable, target :: psi_imag(:,:,:), psi_real(:,:,:)

    integer :: Nx, Ny, Nz

    integer :: errors
    integer(HID_T) :: file_id, dset_id
    integer(HSIZE_T), dimension(3) :: dims

    dims(1) = size(x)
    dims(2) = size(y)
    dims(3) = size(z)

    allocate(psi_imag(dims(1),dims(2),dims(3)))
    allocate(psi_real(dims(1),dims(2),dims(3)))
    allocate(readin_wav(dims(1),dims(2),dims(3)))

    write(*,*) "initial wavefunction read from file"
    ! open hdf5 API
    call h5open_f(errors)

    ! open wavefunction file
    if (comp == 1) then
      call h5fopen_f('psi1_init.h5', H5F_ACC_RDONLY_F, file_id, errors)
    elseif(comp == 2) then
      call h5fopen_f('psi2_init.h5', H5F_ACC_RDONLY_F, file_id, errors)
    end if
    ! open imaginary component dataset
    call h5dopen_f(file_id, 'psi_imag', dset_id, errors)
    f_ptr = C_LOC(psi_imag(1,1,1))
    ! read in imaginary component
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, errors)
    call h5dclose_f(dset_id, errors)

    ! open real component dataset
    call h5dopen_f(file_id, 'psi_real', dset_id, errors)
    f_ptr = C_LOC(psi_real(1,1,1))
    ! read in real component
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, errors)
    call h5dclose_f(dset_id, errors)

    call h5fclose_f(file_id, errors)

    ! close hdf5 API
    call h5close_f(errors)

    ! construct wavefunction from real and imaginary components
    readin_wav = psi_real + cmplx(0.0,1.0)*psi_imag

    deallocate(psi_imag)
    deallocate(psi_real)
  end function readin_wav

  function readin_chempot(comp)
    
    integer, intent(in) :: comp
    double precision :: readin_chempot
    ! Local variables 

    integer :: errors
    integer(HID_T) :: file_id, dset_id
    integer(HSIZE_T), dimension(1) :: scal_dim=1

    write(*,*) "initial chemical potential read from file"
    ! open hdf5 API
    call h5open_f(errors)

    ! open wavefunction file
    if (comp == 1) then
      call h5fopen_f('psi1_init.h5', H5F_ACC_RDONLY_F, file_id, errors)
    elseif(comp == 2) then
      call h5fopen_f('psi2_init.h5', H5F_ACC_RDONLY_F, file_id, errors)
    end if

    ! open chemical potential dataset
    call h5dopen_f(file_id, 'mu', dset_id, errors)
    ! read in imaginary component
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, readin_chempot, scal_dim, errors)
    call h5dclose_f(dset_id, errors)

    call h5fclose_f(file_id, errors)

    ! close hdf5 API
    call h5close_f(errors)
    
    write(*,*) "using mu:"
    write(*,*) readin_chempot

  end function readin_chempot
end module init

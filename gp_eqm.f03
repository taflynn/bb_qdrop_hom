!> \\file gp_lck.f03
program gp_lck
  use OMP_LIB
  use json_module
  use HDF5
  use fft
  use grid
  use init
  use time

  implicit none

  complex(C_DOUBLE_COMPLEX), allocatable :: psi1(:,:,:), psi2(:,:,:)

  integer :: Nx, Ny, Nz
  double precision :: dx, dy, dz 
  double precision :: im_dt_coef, re_dt_coef
  
  integer :: im_t_steps, re_t_steps
  integer :: im_t_save, re_t_save
  
  integer :: init_type1, init_type2
  double precision :: gauss_sig1, gauss_sig2

  double precision :: N1, N2
  double precision :: alpha, beta, eta

  integer :: pot_type1=1, pot_type2=2
  double precision :: omgx1, omgy1, omgz1
  double precision :: omgx2, omgy2, omgz2

  integer :: im_real
  double complex :: dt
  double precision :: mu1, mu2

  double precision, allocatable :: x(:), y(:), z(:)
  double precision, allocatable :: kx(:), ky(:), kz(:)

  complex(C_DOUBLE_COMPLEX), allocatable :: dk2(:,:,:)
  
  double precision, allocatable :: V1(:, :, :), V2(:, :, :)

  ! json-fortran parameters
  type(json_file) :: json
  type(json_core) ::jCore
  logical :: is_found
  type(json_value), pointer :: parent_ptr, child_ptr, c_child_ptr

  ! HDF5 parameters
  integer :: error
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dspace_id
  integer(HSIZE_T), dimension(1) :: dims_r
  character(len=7) :: filename_grid = 'grid.h5'

  ! initialising the json_file object
  call json%initialize()

  ! loading in the input file
  call json%load(filename='config.json') 

  call json%print()

  call json%get_core(jCore)
  
  ! Grid dimensions
  call json%get('grid_size', parent_ptr, is_found)
  call jCore%get_child(parent_ptr, 'Nx', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, Nx); end if

  call jCore%get_child(parent_ptr, 'Ny', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, Ny); end if

  call jCore%get_child(parent_ptr, 'Nz', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, Nz); end if
  ! Grid spacing
  call json%get('grid_space', parent_ptr, is_found)
  call jCore%get_child(parent_ptr, 'dx', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, dx); end if

  call jCore%get_child(parent_ptr, 'dy', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, dy); end if

  call jCore%get_child(parent_ptr, 'dz', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, dz); end if

  ! Time data
  call json%get('time', parent_ptr, is_found)
  call jCore%get_child(parent_ptr, 'im_dt_coef', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, im_dt_coef); end if

  call jCore%get_child(parent_ptr, 're_dt_coef', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, re_dt_coef); end if

  call jCore%get_child(parent_ptr, 'im_t_steps', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, im_t_steps); end if

  call jCore%get_child(parent_ptr, 're_t_steps', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, re_t_steps); end if

  call jCore%get_child(parent_ptr, 'im_t_save', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, im_t_save); end if

  call jCore%get_child(parent_ptr, 're_t_save', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, re_t_save); end if

  ! Initial wavefunction
  call json%get('wav_init', parent_ptr, is_found)
  call jCore%get_child(parent_ptr, 'init_type1', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, init_type1); end if

  call jCore%get_child(parent_ptr, 'init_type2', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, init_type2); end if

  call jCore%get_child(parent_ptr, 'gauss_sig1', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, gauss_sig1); end if

  call jCore%get_child(parent_ptr, 'gauss_sig2', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, gauss_sig2); end if

  ! Population numbers
  call json%get('pop_nums', parent_ptr, is_found)
  call jCore%get_child(parent_ptr, 'N1', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, N1); end if

  call jCore%get_child(parent_ptr, 'N2', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, N2); end if
 
  ! Parameters
  call json%get('params', parent_ptr, is_found)
  call jCore%get_child(parent_ptr, 'alpha', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, alpha); end if

  call jCore%get_child(parent_ptr, 'beta', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, beta); end if

  call jCore%get_child(parent_ptr, 'eta', child_ptr, is_found)
  if (is_found) then; call jCore%get(child_ptr, eta); end if
    ! Potentials
  call json%get('potentials', parent_ptr, is_found)
  call jCore%get_child(parent_ptr, 'pot1', child_ptr, is_found)
  call jCore%get_child(child_ptr, 'omgx1', c_child_ptr, is_found)
  if (is_found) then; call jCore%get(c_child_ptr, omgx1); end if
  call jCore%get_child(child_ptr, 'omgy1', c_child_ptr, is_found)
  if (is_found) then; call jCore%get(c_child_ptr, omgy1); end if
  call jCore%get_child(child_ptr, 'omgz1', c_child_ptr, is_found)
  if (is_found) then; call jCore%get(c_child_ptr, omgz1); end if

  call jCore%get_child(parent_ptr, 'pot2', child_ptr, is_found)
  call jCore%get_child(child_ptr, 'omgx2', c_child_ptr, is_found)
  if (is_found) then; call jCore%get(c_child_ptr, omgx2); end if
  call jCore%get_child(child_ptr, 'omgy2', c_child_ptr, is_found)
  if (is_found) then; call jCore%get(c_child_ptr, omgy2); end if
  call jCore%get_child(child_ptr, 'omgz2', c_child_ptr, is_found)
  if (is_found) then; call jCore%get(c_child_ptr, omgz2); end if


  call json%destroy()

  ! set up 3D spatial grid
  x = space_grid(Nx,dx)
  y = space_grid(Ny,dy)
  z = space_grid(Nz,dz)

  ! set up 3D momentum space grid
  kx = mom_grid(Nx,dx)
  ky = mom_grid(Ny,dy)
  kz = mom_grid(Nz,dz)

  ! initiate the hdf5 environment
  call h5open_f(error)
  ! create the grid.h5 file
  call h5fcreate_f(filename_grid, H5F_ACC_TRUNC_F, file_id, error)
  ! saving x-array
  dims_r = Nx
  call h5screate_simple_f(1, dims_r, dspace_id, error); if (Nx .ne. Ny .and. Nx .ne. Nz) stop 
  call h5dcreate_f(file_id, 'x', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, dims_r, error)
  call h5dclose_f(dset_id, error)
  ! saving y-array
  dims_r = Ny
  call h5dcreate_f(file_id, 'y', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, y, dims_r, error)
  call h5dclose_f(dset_id, error)
  ! saving z-array
  dims_r = Nz
  call h5dcreate_f(file_id, 'z', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, z, dims_r, error)
  call h5dclose_f(dset_id, error)

  ! saving kx-array
  dims_r = Nx
  call h5dcreate_f(file_id, 'kx', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, kx, dims_r, error)
  call h5dclose_f(dset_id, error)
  ! saving ky-array
  dims_r = Ny
  call h5dcreate_f(file_id, 'ky', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ky, dims_r, error)
  call h5dclose_f(dset_id, error)
  ! saving kz-array
  dims_r = Nz
  call h5dcreate_f(file_id, 'kz', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, kz, dims_r, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  ! close the grid.h5 file
  call h5fclose_f(file_id, error)
  
  ! compute initial profile of wavefunction 1
  if (init_type1 == 0) then
    ! read in the initial wavefunction from a .h5 file
    psi1 = readin_wav(x,y,z,1)
  elseif (init_type1 /= 0) then
    ! calculate the initial wavefunction
    psi1 = init_wav(x,y,z,init_type1,gauss_sig1)
    call renorm(psi1,dx,dy,dz,N1)
  end if
  
  ! compute initial profile of wavefunction 2
  if (init_type2 == 0) then
    ! read in the initial wavefunction from a .h5 file
    psi2 = readin_wav(x,y,z,2)
  elseif (init_type2 /= 0) then
    ! calculate the initial wavefunction
    psi2 = init_wav(x,y,z,init_type2,gauss_sig2)
    call renorm(psi2,dx,dy,dz,N2)
  end if

  V1 = init_pot(x, y, z, omgx1, omgy1, omgz1, pot_type1)
  V2 = init_pot(x, y, z, omgx2, omgy2, omgz2, pot_type2)

  ! begin time-stepping
  if (im_t_steps > 0) then
    write(*,*) "beginning imaginary time"
  
    ! imaginary time step
    dt = im_dt_coef*min(dx,dy,dz)**2
    
    ! state that the time-stepping should expect imaginary time
    im_real = 0
    
    ! initialise ssfm laplacian term and chemical potentials
    dk2 = exp_lap(kx,ky,kz,dt)
    if (init_type1 == 0) then
      mu1 = readin_chempot(1)
    else
      mu1 = -1.0
    end if
    if (init_type2 == 0) then
      mu2 = readin_chempot(2)
    else
      mu2 = -1.0
    end if

    ! imaginary time function
    call ssfm(psi1,psi2,dk2,im_t_steps,im_t_save,dt,dx,dy,dz,V1,V2,N1,N2,alpha,beta,eta,mu1,mu2,im_real)
  end if
  if (re_t_steps > 0) then
    write(*,*) "beginning real time"
    
    ! real time step
    dt = complex(0.0,1.0)*re_dt_coef*min(dx,dy,dz)**2
    
    ! state that the time-stepping should expect real time
    im_real = 1

    ! initilise ssfm laplacian term and chemical potential
    dk2 = exp_lap(kx,ky,kz,dt)
    if (init_type1 == 0) then
      mu1 = readin_chempot(1)
    end if
    if (init_type2 == 0) then
      mu2 = readin_chempot(2)
    end if
    
    ! real time function
    call ssfm(psi1,psi2,dk2,im_t_steps,im_t_save,dt,dx,dy,dz,V1,V2,N1,N2,alpha,beta,eta,mu1,mu2,im_real)
  end if
  if (im_t_steps == 0 .and. re_t_steps == 0) then
    ! if there are no time-steps for imaginary and real time then stop program 
    stop
  end if
end program

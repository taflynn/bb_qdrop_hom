!> \\file time.f03
module time
  use OMP_LIB
  use HDF5
  use fft
  use rhs

  implicit none

  contains

  ! main ssfm time-stepping function
  subroutine ssfm(psi1,psi2,dk2,t_steps,t_save,dt,dx,dy,dz,V1,V2,N1,N2,alpha,beta,eta,mu1,mu2,im_real)

    integer, intent(in) :: t_steps, t_save, im_real
    double precision, intent(in) :: dx, dy, dz
    double complex, intent(in) :: dt
    double precision, intent(in) :: N1, N2
    double precision, intent(in) :: alpha, beta, eta
    complex(C_DOUBLE_COMPLEX), intent(in) :: dk2(:,:,:)
    double precision, intent(in) :: V1(:,:,:), V2(:,:,:)
    
    complex(C_DOUBLE_COMPLEX), intent(inout) :: psi1(:,:,:), psi2(:,:,:)
    double precision, intent(inout) :: mu1, mu2

    ! local variables 
    integer :: l
    integer :: Nx, Ny, Nz, Nxyz

    double complex :: t = 0.0

    type(C_PTR) :: plan_forw, plan_back
    
    complex(C_DOUBLE_COMPLEX), allocatable :: psi1_k(:,:,:), psi2_k(:,:,:)
    complex(C_DOUBLE_COMPLEX), allocatable :: psi1_prev(:,:,:), psi2_prev(:,:,:)
    
    integer(C_INT) :: nthreads
    integer(C_INT) :: error
  
    character(len=4) :: t_current 
    integer :: h5_error 
    integer(HID_T) :: file_id
    integer(HID_T) :: dset_id
    integer(HID_T) :: dspace_id
    integer(HSIZE_T), dimension(1) :: scal_dim=(/0/)
    integer(HSIZE_T), dimension(3) :: dims
    character(len=19) :: filename_wav

    call h5open_f(h5_error)

    dims = shape(psi1)
    Nx = dims(1)
    Ny = dims(2)
    Nz = dims(3)

    Nxyz = Nx*Ny*Nz

    allocate(psi1_k(Nx,Ny,Nz))
    allocate(psi1_prev(Nx,Ny,Nz))
    allocate(psi2_k(Nx,Ny,Nz))
    
    ! initialising FFTW with threads
    error = fftw_init_threads()
    nthreads = omp_get_max_threads()
    call fftw_plan_with_nthreads(int(nthreads,C_INT))
  
    if (im_real == 0) then 
      write(*,*) "number of threads:" 
      write(*,*) nthreads
    end if 

    ! constructing FFTW plans (just using component 1)
    plan_forw = fftw_plan_dft_3d(Nz,Ny,Nx,psi1,psi1_k,FFTW_FORWARD,FFTW_ESTIMATE)
    plan_back = fftw_plan_dft_3d(Nz,Ny,Nx,psi1_k,psi1,FFTW_BACKWARD,FFTW_ESTIMATE)

    do l = 1, t_steps

      psi1_prev = psi1

      ! first component 1 half-step (non-linear terms)
      call V_rhs1(psi1,psi2,mu1,V1,alpha,beta,eta,dt,Nx,Ny,Nz)
      
      ! first component 2 half-step (non-linear terms)
      call V_rhs2(psi1_prev,psi2,mu2,V2,alpha,beta,eta,dt,Nx,Ny,Nz)

      ! FFT wavefunction 1 to momentum space
      call fftw_execute_dft(plan_forw,psi1,psi1_k)
      
      ! FFT wavefunction 2 to momentum space
      call fftw_execute_dft(plan_forw,psi2,psi2_k)

      ! kinetic energy step of component 1
      call T_rhs(psi1_k,dk2,Nx,Ny,Nz)

      ! kinetic energy step of component 2
      call T_rhs(psi2_k,dk2,Nx,Ny,Nz)

      ! IFFT wavefunction 1 to real space
      call fftw_execute_dft(plan_back,psi1_k,psi1)
      !$omp parallel workshare
      psi1(:,:,:) = psi1(:,:,:)/dble(Nxyz)
      !$omp end parallel workshare

      ! IFFT wavefunction 2 to real space
      call fftw_execute_dft(plan_back,psi2_k,psi2)
      !$omp parallel workshare
      psi2(:,:,:) = psi2(:,:,:)/dble(Nxyz)
      !$omp end parallel workshare

      psi1_prev = psi1

      ! last component 1 half-step (non-linear terms)
      call V_rhs1(psi1,psi2,mu1,V1,alpha,beta,eta,dt,Nx,Ny,Nz)

      ! last component 2 half-step (non-linear terms)
      call V_rhs2(psi1_prev,psi2,mu2,V2,alpha,beta,eta,dt,Nx,Ny,Nz)

      ! in imaginary time: renormalise and compute chemical potentials
      if (im_real == 0) then
        ! renormalise wavefunction 1
        call renorm(psi1,dx,dy,dz,N1)
        ! renomalise wavefunction 2
        call renorm(psi2,dx,dy,dz,N2)
        ! FFT wavefunction 1 to real space
        call fftw_execute_dft(plan_forw,psi1,psi1_k)
        ! FFT wavefunction 2 to real space
        call fftw_execute_dft(plan_forw,psi2,psi2_k)
        ! chemical potential 1
        mu1 = chem_pot1(psi1,psi2,psi1_k,dk2,plan_back,V1,Nx,Ny,Nz,dt,alpha,beta,eta)
        ! chemical potential 2
        mu2 = chem_pot2(psi1,psi2,psi2_k,dk2,plan_back,V2,Nx,Ny,Nz,dt,alpha,beta,eta)
      end if

      t = t + dt 

      ! data outputting
      if (mod(l,t_save) == 0) then
        ! writing to screen
        write(*,*) "Percentage Completed"
        write(*,*) 100*(dble(l)/dble(T_STEPS))
        write(*,*) "Central Density 1"
        write(*,*) abs(psi1(1+Nx/2,1+Ny/2,1+Nz/2))**2
        write(*,*) "Central Density 2"
        write(*,*) abs(psi2(1+Nx/2,1+Ny/2,1+Nz/2))**2
        write(*,*) "Chemical Potential 1"
        write(*,*) mu1
        write(*,*) "Chemical Potential 2"
        write(*,*) mu2
        ! outputting
        write(t_current,'(I4.4)') 100*l/T_STEPS
        if (im_real == 0) then
          filename_wav = 'psi_im_t_'//trim(t_current)//'.h5'
        elseif (im_real == 1) then
          filename_wav = 'psi_re_t_'//trim(t_current)//'.h5'
        end if
        ! create wavefunction output file
        call h5fcreate_f(filename_wav, H5F_ACC_TRUNC_F, file_id, h5_error)
        ! save chemical potential 1
        call h5screate_simple_f(0, scal_dim, dspace_id, h5_error)
        call h5dcreate_f(file_id, 'mu1', H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5_error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mu1, scal_dim, h5_error)
        call h5dclose_f(dset_id, h5_error)
        ! save chemical potential 2
        call h5dcreate_f(file_id, 'mu2', H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5_error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mu2, scal_dim, h5_error)
        call h5dclose_f(dset_id, h5_error)
        ! save current time step
        call h5dcreate_f(file_id, 't', H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5_error)
        if (im_real == 0) then
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, real(t), scal_dim, h5_error)
        elseif (im_real == 1) then
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, aimag(t), scal_dim, h5_error)
        end if
        call h5dclose_f(dset_id, h5_error)
        call h5sclose_f(dspace_id, h5_error)
        ! save real component of wavefunction 1
        call h5screate_simple_f(3, dims, dspace_id, h5_error)
        call h5dcreate_f(file_id, 'psi1_real', H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5_error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, real(psi1), dims, h5_error)
        call h5dclose_f(dset_id, h5_error)
        ! save imaginary component of wavefunction 1
        call h5dcreate_f(file_id, 'psi1_imag', H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5_error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, aimag(psi1), dims, h5_error)
        call h5dclose_f(dset_id, h5_error)
        ! save real component of wavefunction 2
        call h5dcreate_f(file_id, 'psi2_real', H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5_error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, real(psi2), dims, h5_error)
        call h5dclose_f(dset_id, h5_error)
        ! save imaginary component of wavefunction 2
        call h5dcreate_f(file_id, 'psi2_imag', H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5_error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, aimag(psi2), dims, h5_error)
        call h5dclose_f(dset_id, h5_error)
        call h5sclose_f(dspace_id, h5_error)
        call h5fclose_f(file_id, h5_error)
      end if
    end do
    call h5close_f(h5_error)
  end subroutine ssfm

  subroutine renorm(psi,dx,dy,dz,N)
    
    double precision, intent(in) :: dx, dy, dz
    double precision, intent(in) :: N
   
    complex(C_DOUBLE_COMPLEX), intent(inout) :: psi(:,:,:)

    ! local variables
    double precision :: norm
    
    !$omp parallel workshare
    norm = sum(abs(psi(:,:,:))**2.0)*dx*dy*dz
    !$omp end parallel workshare

    !$omp parallel workshare
    psi(:,:,:) = psi(:,:,:)*sqrt(N/norm)
    !$omp end parallel workshare

  end subroutine renorm

  function chem_pot1(psi1,psi2,psi1_k,dk2,plan_back,V1,Nx,Ny,Nz,dt,alpha,beta,eta)
    
    type(C_PTR), intent(in) :: plan_back
    
    integer, intent(in) :: Nx, Ny, Nz
    
    double precision, intent(in) :: alpha, beta, eta

    double complex, intent(in) :: dk2(:,:,:)

    complex(C_DOUBLE_COMPLEX), intent(in) :: psi1(:,:,:), psi2(:,:,:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: psi1_k(:,:,:)
    double precision, intent(in) :: V1(:,:,:)

    double precision :: chem_pot1

    ! local variables
    integer :: i, j, k
   
    double complex :: dt

    complex(C_DOUBLE_COMPLEX), allocatable :: lap_psi1(:,:,:), lap_psi1k(:,:,:)

    allocate(lap_psi1(Nx,Ny,Nz))
    allocate(lap_psi1k(Nx,Ny,Nz))

    ! compute second derivative of wavefunction in momentum space
    !$omp parallel do collapse(3)   
    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          lap_psi1k(i,j,k) = (2.0/dt)*log(dk2(i,j,k))*psi1_k(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do
    
    ! transform kinetic energy term back into real space
    call fftw_execute_dft(plan_back,lap_psi1k,lap_psi1)
    !$omp parallel workshare
    lap_psi1 = lap_psi1/dble(Nx*Ny*Nz)
    !$omp end parallel workshare
    
    deallocate(lap_psi1k)

    ! compute chemical potential
    !$omp parallel workshare
    chem_pot1 = sum(-0.5*conjg(psi1(:,:,:))*lap_psi1(:,:,:) &
                    +abs(psi1(:,:,:))**4.0 + eta*abs(psi1(:,:,:))**2.0*abs(psi2(:,:,:))**2.0 &
                    +alpha*(abs(psi1(:,:,:))**2.0 &
                    +beta*abs(psi2(:,:,:))**2.0)**1.5*abs(psi1(:,:,:))**2.0)/sum(abs(psi1(:,:,:))**2.0)
    !$omp end parallel workshare

  end function chem_pot1

  function chem_pot2(psi1,psi2,psi2_k,dk2,plan_back,V2,Nx,Ny,Nz,dt,alpha,beta,eta)
    
    type(C_PTR), intent(in) :: plan_back
    
    integer, intent(in) :: Nx, Ny, Nz
    
    double precision, intent(in) :: alpha, beta, eta

    double complex, intent(in) :: dk2(:,:,:)

    complex(C_DOUBLE_COMPLEX), intent(in) :: psi1(:,:,:), psi2(:,:,:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: psi2_k(:,:,:)
    double precision, intent(in) :: V2(:,:,:)

    double precision :: chem_pot2

    ! local variables
    integer :: i, j, k
   
    double complex :: dt

    complex(C_DOUBLE_COMPLEX), allocatable :: lap_psi2(:,:,:), lap_psi2k(:,:,:)

    allocate(lap_psi2(Nx,Ny,Nz))
    allocate(lap_psi2k(Nx,Ny,Nz))

    ! compute second derivative of wavefunction in momentum space
    !$omp parallel do collapse(3)   
    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          lap_psi2k(i,j,k) = (2.0/dt)*log(dk2(i,j,k))*psi2_k(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do
    
    ! transform kinetic energy term back into real space
    call fftw_execute_dft(plan_back,lap_psi2k,lap_psi2)
    !$omp parallel workshare
    lap_psi2(:,:,:) = lap_psi2(:,:,:)/dble(Nx*Ny*Nz)
    !$omp end parallel workshare
    
    deallocate(lap_psi2k)

    ! compute chemical potential
    !$omp parallel workshare
    chem_pot2 = sum(-0.5*conjg(psi2(:,:,:))*lap_psi2(:,:,:) &
                    +beta*abs(psi2(:,:,:))**4.0 + eta*beta*abs(psi1(:,:,:))**2.0*abs(psi2(:,:,:))**2.0 &
                    +alpha*beta**2.0*(abs(psi1(:,:,:))**2.0 &
                    +beta*abs(psi2(:,:,:))**2.0)**1.5*abs(psi2(:,:,:))**2.0)/sum(abs(psi2(:,:,:))**2.0)
    !$omp end parallel workshare

  end function chem_pot2

end module time

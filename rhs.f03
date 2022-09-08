!> \\file rhs.f03
module rhs
 
  use fft

  implicit none

  contains

  ! generate the initial form of the wavefunction 
  subroutine V_rhs1(psi1,psi2,mu1,alpha,beta,eta,dt,Nx,Ny,Nz)
    implicit none
    
    integer, intent(in) :: Nx, Ny, Nz
    double precision, intent(in) :: mu1
    double complex, intent(in) :: dt
    double precision, intent(in) :: alpha, beta, eta

    complex(C_DOUBLE_COMPLEX), intent(in) :: psi2(:,:,:)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: psi1(:,:,:)

    ! local variables
    integer :: i, j, k
    
    ! iterate across x, y and z of the two-body and Lee-Huang-Yang terms
    !$omp parallel do collapse(3)
    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          psi1(i,j,k) = psi1(i,j,k)*exp(-0.5*dt*(abs(psi1(i,j,k))**2.0 + eta*abs(psi2(i,j,k))**2.0 &
                                       + alpha*(abs(psi1(i,j,k))**2.0 + beta*abs(psi2(i,j,k))**2.0)**1.5 - mu1))
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine V_rhs1
  
  subroutine V_rhs2(psi1,psi2,mu2,alpha,beta,eta,dt,Nx,Ny,Nz)
    implicit none
    
    integer, intent(in) :: Nx, Ny, Nz
    double precision, intent(in) :: mu2
    double complex, intent(in) :: dt
    double precision, intent(in) :: alpha, beta, eta

    complex(C_DOUBLE_COMPLEX), intent(in) :: psi1(:,:,:)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: psi2(:,:,:)

    ! local variables
    integer :: i, j, k
    
    ! iterate across x, y and z of the two-body and Lee-Huang-Yang terms
    !$omp parallel do collapse(3)
    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          psi2(i,j,k) = psi2(i,j,k)*exp(-0.5*dt*(beta*abs(psi2(i,j,k))**2.0 + eta*beta*abs(psi1(i,j,k))**2.0 &
                                       + alpha*beta**2.0*(abs(psi1(i,j,k))**2.0 + beta*abs(psi2(i,j,k))**2.0)**1.5 - mu2))
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine V_rhs2

  subroutine T_rhs(psi_k,dk2,Nx,Ny,Nz)
    implicit none

    integer, intent(in) :: Nx, Ny, Nz
    complex(C_DOUBLE_COMPLEX), intent(in) :: dk2(:,:,:)

    complex(C_DOUBLE_COMPLEX), intent(inout) :: psi_k(:,:,:)

    integer :: i, j, k
    
    ! iterate across kx, ky and kz of the kinetic energy term
    !$omp parallel do collapse(3)
    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          psi_k(i,j,k) = dk2(i,j,k)*psi_k(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine T_rhs
end module rhs

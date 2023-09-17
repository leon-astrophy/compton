!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Setting up everything to compute the xray polarization
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SETUP
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose between GRMHD or analytic models !
IF(analytic) THEN
  
  ! Analytic model !
  CALL setup_analytic

ELSEIF(grmhd) THEN
 
  ! load data file from GRMHD snapshots !
  call setup_grmhd

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Setup the computational domain for analytic model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SETUP_ANALYTIC
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, assign nx, ny, nz !

nx = nr
ny = nth
nz = nphi

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Setup the computational domain for grmhd simulation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SETUP_GRMHD
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, assign nx, ny, nz !

CALL HDF5_OPEN

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Populate the arrays for doing the integration
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BUILD_ARRAYS
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate arrays !

! Viewing angle !
ALLOCATE(theta_0(1:n_angle))
ALLOCATE(pdegree(1:n_angle))
ALLOCATE(stokes_total(1:3,1:n_angle))

! Grid !
ALLOCATE(r_grid(1:nx,1:ny,1:nz))
ALLOCATE(th_grid(1:nx,1:ny,1:nz))
ALLOCATE(phi_grid(1:nx,1:ny,1:nz))
ALLOCATE(vol(1:nx,1:ny,1:nz))
ALLOCATE(gdet(1:nx,1:ny,1:nz))

! Fluid !
ALLOCATE(rho(1:nx,1:ny,1:nz))
ALLOCATE(beta(1:nx,1:ny,1:nz))
ALLOCATE(gamma(1:nx,1:ny,1:nz))
ALLOCATE(gammam1(1:nx,1:ny,1:nz))
ALLOCATE(beta_vel(1:nx,1:ny,1:nz))
ALLOCATE(nebar(1:nx,1:ny,1:nz))

! Electron !
ALLOCATE(nhat_x(1:nx,1:ny,1:nz))
ALLOCATE(nhat_y(1:nx,1:ny,1:nz))
ALLOCATE(nhat_z(1:nx,1:ny,1:nz))

! Angular factor !
ALLOCATE(ang_fac_isc(1:nx,1:ny,1:nz))
ALLOCATE(ang_fac_qsc(1:nx,1:ny,1:nz))
ALLOCATE(ang_fac_usc(1:nx,1:ny,1:nz))

! Angular factor !
ALLOCATE(ang_fac_out(1:2,1:3,1:nx,1:ny,1:nz))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign viewing angle !

DO i = 1, n_angle
  theta_0(i) = (DBLE(i) - 0.5d0)*pi/DBLE(n_angle)
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Setting up everything to compute the xray polarization
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SETUP
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONVERT ALL ANGLES INTO RADIANS !

theta_j = theta_j*pi/180.0d0
theta_observe = theta_observe*pi/180.0d0
phi_j = phi_j*pi/180.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONVERT FREQUENCY INTO HZ, CONVERT TEMPERATURE INTO KELVIN
nu_light = 10.0d0**(nu_light)
bbody_temp = 10.0d0**(bbody_temp)

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
! Number of angles 

! Set the number of viewing angles to 1 for discrete !
n_angle = 1

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

! Assign according to parameter file !
nx = nr
ny = nth
nz = nphi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

! Assign nx, ny, nz from GRMHD snapshots !
CALL HDF5_OPEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Populate the arrays for doing the integration
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ARRAY
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate arrays !

! Viewing angle !
ALLOCATE(theta_0(1:n_angle))

! Grid !
ALLOCATE(r_grid(1:nx,1:ny,1:nz))
ALLOCATE(th_grid(1:nx,1:ny,1:nz))
ALLOCATE(phi_grid(1:nx,1:ny,1:nz))

! Fluid !
ALLOCATE(nrho(1:nx,1:ny,1:nz))
ALLOCATE(gamma(1:nx,1:ny,1:nz))
ALLOCATE(beta_vel(1:nx,1:ny,1:nz))

! Electron !
ALLOCATE(nhat_r(1:nx,1:ny,1:nz))
ALLOCATE(nhat_th(1:nx,1:ny,1:nz))
ALLOCATE(nhat_phi(1:nx,1:ny,1:nz))

! Photon !
ALLOCATE(khat_r(1:nx,1:ny,1:nz))
ALLOCATE(khat_th(1:nx,1:ny,1:nz))
ALLOCATE(khat_phi(1:nx,1:ny,1:nz))

! Angular factor !
ALLOCATE(fsc(1:nx,1:ny,1:nz))
ALLOCATE(gsc(1:nx,1:ny,1:nz))
ALLOCATE(hsc(1:nx,1:ny,1:nz))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For double integral !

ALLOCATE(rhop_dir(0:n_integral))
ALLOCATE(alphap_dir(0:n_integral))
ALLOCATE(cosalpha_p(0:n_integral))
ALLOCATE(sinalpha_p(0:n_integral))
ALLOCATE(cosrho_p(0:n_integral))
ALLOCATE(sinrho_p(0:n_integral))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign viewing angle !

! Choose between full-angle calculation or single angle !
DO i = 1, n_angle
  theta_0(i) = theta_observe
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
drhop = (x_high - x_low)/DBLE(n_integral)
dalphap = (y_high - y_low)/DBLE(n_integral)

! Get x-coordinate !
DO i = 0, n_integral
  rhop_dir(i) = x_low + i*drhop
  cosrho_p(i) = DCOS(rhop_dir(i))
  sinrho_p(i) = DSIN(rhop_dir(i))
END DO

! Get y-coordinate !
DO i = 0, n_integral
  alphap_dir(i) = y_low + i*dalphap
  cosalpha_p(i) = DCOS(alphap_dir(i))
  sinalpha_p(i) = DSIN(alphap_dir(i))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

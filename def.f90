!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Definition file containing global variables and arrays
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DEFINITION
IMPLICIT NONE
INCLUDE "param.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The dimension of the system !
INTEGER :: nx, ny, nz

! Current viewing angle !
REAL*8 :: theta_view 

! Step in the double integral !
REAL*8 :: drhop, dalphap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for black body !
REAL*8 :: hkt
REAL*8, PARAMETER :: k_boltz = 1.3807D-16 !cm2 g s-2 K-1
REAL*8, PARAMETER :: hplanck = 6.6261D-27 !cm2 g s-1
REAL*8, PARAMETER :: clight = 2.99792458D10 !cms-1
REAL*8, PARAMETER :: bbody_factor = 2.0d0*hplanck/clight**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for fake disc model 

! fixed parameters !
REAL*8, PARAMETER :: r_sh = 2.0d0

! free parameters !
REAL*8 :: gamma_0 = 0.5d0
REAL*8 :: s_gamma = 5.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define the integration range !
REAL*8, PARAMETER :: x_low = 0 
REAL*8, PARAMETER :: x_high = PI
REAL*8, PARAMETER :: y_low = 0 
REAL*8, PARAMETER :: y_high = 2*PI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid, in general should be 3D !

! r !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: r_grid

! theta !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: th_grid

! phi !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: phi_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fluid variables !

! Electron density !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: nrho

! Relativistic Gamma !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: gamma

! Relativistic beta !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: beta_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Electron directional vector !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: nhat_r
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: nhat_th
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: nhat_phi

! Photon directional vector !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: khat_r
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: khat_th
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: khat_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Angular factor !

! Viewing angle !
REAL*8, ALLOCATABLE, DIMENSION(:) :: theta_0

! Scattering factor for stokes I !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: fsc

! Scattering factor for stokes Q !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: gsc

! Scattering factor for stokes U !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: hsc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For doing the double integral !

! Arrays !
REAL*8, ALLOCATABLE, DIMENSION(:) :: rhop_dir
REAL*8, ALLOCATABLE, DIMENSION(:) :: alphap_dir
REAL*8, ALLOCATABLE, DIMENSION(:) :: cosalpha_p
REAL*8, ALLOCATABLE, DIMENSION(:) :: sinalpha_p
REAL*8, ALLOCATABLE, DIMENSION(:) :: cosrho_p
REAL*8, ALLOCATABLE, DIMENSION(:) :: sinrho_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output !

REAL*8 :: isctaui0
REAL*8 :: qsctaui0
REAL*8 :: usctaui0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function !
contains

	! dot product !
  REAL*8 function rdotx(theta_in, phi_in)
  implicit none
  REAL*8 :: theta_in, phi_in
  rdotx = DSIN(theta_in)*DCOS(phi_in)
  end function

	! dot product !
  REAL*8 function rdoty(theta_in, phi_in)
  implicit none
  REAL*8 :: theta_in, phi_in
  rdoty = DSIN(theta_in)*DSIN(phi_in)
  end function

	! dot product !
  REAL*8 function rdotz(theta_in, phi_in)
  implicit none
  REAL*8 :: theta_in, phi_in
  rdotz = DCOS(theta_in)
  end function

	! dot product !
  REAL*8 function thdotx(theta_in, phi_in)
  implicit none
  REAL*8 :: theta_in, phi_in
  thdotx = DCOS(theta_in)*DCOS(phi_in)
  end function

	! dot product !
  REAL*8 function thdoty(theta_in, phi_in)
  implicit none
  REAL*8 :: theta_in, phi_in
  thdoty = DCOS(theta_in)*DSIN(phi_in)
  end function

	! dot product !
  REAL*8 function thdotz(theta_in, phi_in)
  implicit none
  REAL*8 :: theta_in, phi_in
  thdotz = -DSIN(theta_in)
  end function

	! dot product !
  REAL*8 function phidotx(theta_in, phi_in)
  implicit none
  REAL*8 :: theta_in, phi_in
  phidotx = -DSIN(phi_in)
  end function

	! dot product !
  REAL*8 function phidoty(theta_in, phi_in)
  implicit none
  REAL*8 :: theta_in, phi_in
  phidoty = DCOS(phi_in)
  end function

	! dot product !
  REAL*8 function phidotz(theta_in, phi_in)
  implicit none
  REAL*8 :: theta_in, phi_in
  phidotz = 0.0D0
  end function

END MODULE

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

! Total stokes vector !
REAL*8 :: isc_total, qsc_total, usc_total

! Current viewing angle !
REAL*8 :: theta_view 

! Step in the double integral !
REAL*8 :: drhop, dalphap

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

! Volume !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fluid variables !

! Relativistic Gamma !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: gamma

! One over gamma !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: gammam1

! Gamma factor !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: gam_fac

! Relativistic beta !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: beta_vel

! Dimensionless electron density !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: nebar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Electron unit direction vector !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: nhat_x
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: nhat_y
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: nhat_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Angular factor !

! Viewing angle !
REAL*8, ALLOCATABLE, DIMENSION(:) :: theta_0

! Polarization degree
REAL*8, ALLOCATABLE, DIMENSION(:) :: pdegree

! For stokes I !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ang_fac_isc

! For stokes Q !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ang_fac_qsc

! For stokes U !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ang_fac_usc

! Output angular map !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: ang_fac_out

! Total integrated Stokes 
REAL*8, ALLOCATABLE, DIMENSION(:,:) :: stokes_total

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

END MODULE
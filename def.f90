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

! determinant of the black hole metric !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: gdet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fluid variables !

! Mass density !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: rho

! Plasma beta !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: beta

! Relativistic Gamma !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: gamma

! One over gamma !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: gammam1

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

END MODULE
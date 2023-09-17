!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get the initial conditions according to the selected model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MODEL
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose between GRMHD or analytic models !
IF(analytic) THEN
  
  ! Analytic model !
  CALL model_analytic

ELSEIF(grmhd) THEN
 
  ! load data file from GRMHD snapshots !
  CALL model_grmhd

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Setup the initial conditions for the analytic model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MODEL_ANALYTIC
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k
INTEGER :: j_th

! Real !
REAL*8 :: dr, dtheta, dphi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, setup grid !

! r-direction !
dr = 1.0d0/DBLE(nx)
DO i = 1, nx
  r_grid(i,:,:) = (DBLE(i) - 0.5d0)*dr
END DO

! theta-direction !
dtheta = pi/DBLE(ny)
DO j = 1, ny
  th_grid(:,j,:) = (DBLE(j) - 0.5d0)*dtheta
END DO

! phi-direction !
dphi = 2.0d0*pi/DBLE(nz)
DO k = 1, nz
  phi_grid(:,:,k) = (DBLE(k) - 0.5d0)*dphi
END DO

! Volume !
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nz
      vol(i,j,k) = r_grid(i,j,k)**2*DSIN(th_grid(i,j,k))*dr*dtheta*dphi
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Then, setup electron number density !

! Choose according to model, dirac-delta model !
IF(delta) THEN

  ! Electron is located at a particular angle !
  DO j = 1, ny
    IF((th_grid(1,j,1) > theta_j)) THEN
      j_th = j
      EXIT
    END IF
  END DO

  ! Assign dirac-delta distribution !
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz

        ! Electron number density and 
        IF(j == j_th) THEN
          nebar(i,j,k) = dphi/(2.0d0*pi)/DBLE(nx)/vol(i,j,k)
          gamma(i,j,k) = gamma_j
          gammam1(i,j,k) = 1.0d0/gamma_j
          beta_vel(i,j,k) = DSQRT(1.0d0 - 1.0d0/gamma_j**2)
        ELSE
          nebar(i,j,k) = 0.0d0
          gamma(i,j,k) = 1.0d0
          gammam1(i,j,k) = 1.0d0
          beta_vel(i,j,k) = 0.0D0
        END IF

        ! Electron directional vector, orthonomral basis !
        nhat_x(i,j,k) = DSIN(th_grid(i,j,k))*DCOS(phi_grid(i,j,k))
        nhat_y(i,j,k) = DSIN(th_grid(i,j,k))*DSIN(phi_grid(i,j,k))
        nhat_z(i,j,k) = DCOS(th_grid(i,j,k))

      END DO
    END DO
  END DO

ELSE IF(gaussian) THEN

  ! Not yet implemented 
  STOP 'gaussian beam not implemented'

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Setup the initial conditions for the grmhd model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MODEL_GRMHD
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, load from HDF5 !
CALL HDF5_LOAD

! Then, assign !
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nz

      ! one over gamma and beta !
      gammam1(i,j,k) = 1.0d0/gamma(i,j,k)
      beta_vel(i,j,k) = DSQRT(1.0d0 - gammam1(i,j,k)**2)

      ! Electron directional vector, orthonomral basis !
      nhat_x(i,j,k) = DSIN(th_grid(i,j,k))*DCOS(phi_grid(i,j,k))
      nhat_y(i,j,k) = DSIN(th_grid(i,j,k))*DSIN(phi_grid(i,j,k))
      nhat_z(i,j,k) = DCOS(th_grid(i,j,k))

    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

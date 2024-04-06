!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get the initial conditions according to the selected model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INITIAL
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose between GRMHD or analytic models !
IF(analytic) THEN
  
  ! Analytic model !
  CALL initial_analytic

ELSEIF(grmhd) THEN
 
  ! load data file from GRMHD snapshots !
  CALL initial_grmhd

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set black body constant !
IF(bbody) hkt = hplanck/k_boltz/bbody_temp

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Setup the initial conditions for the analytic model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INITIAL_ANALYTIC
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k
INTEGER :: j_th

! Real !
REAL*8 :: dr, dtheta, dphi, s_eq, volume

! Real !
REAL*8 :: kdotr, kdotth, kdotphi
REAL*8 :: xdotr, xdotth, xdotphi
REAL*8 :: zdotr, zdotth, zdotphi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, setup grid !

! r-direction !
IF(disc) THEN
  dr = (40.0d0 - r_sh)/DBLE(nx)
  DO i = 1, nx
    r_grid(i,:,:) = r_sh + (DBLE(i) - 0.5d0)*dr
  END DO
ELSE
  dr = 1.0d0/DBLE(nx)
  DO i = 1, nx
    r_grid(i,:,:) = (DBLE(i) - 0.5d0)*dr
  END DO
END IF

! theta-direction !
dtheta = pi/DBLE(ny)
DO j = 1, ny
  th_grid(:,j,:) = theta_j !(DBLE(j) - 0.5d0)*dtheta
END DO

! phi-direction !
dphi = 2.0d0*pi/DBLE(nz)
DO k = 1, nz
  phi_grid(:,:,k) = (DBLE(k) - 0.5d0)*dphi
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose by disc or beamed jet model !

IF(disc) THEN

  ! Assign dirac-delta distribution !
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz

        ! polar radius !
        s_eq = r_grid(i,j,k)*DSIN(th_grid(i,j,k))

        ! Inititalize !
        gamma(i,j,k) = gamma_0*DEXP(-(s_eq - r_sh)**2/s_gamma**2) + 1.0d0
        beta_vel(i,j,k) = DSQRT(1.0d0 - 1.0d0/gamma(i,j,k)**2)

        ! Electron directional vector, orthonomral basis !
        nhat_r(i,j,k) = 1.0d0
        nhat_th(i,j,k) = 0.0d0
        nhat_phi(i,j,k) = 0.0d0

        ! Photon directional vector, orthonomral basis !
        xdotr = rdotx(th_grid(i,j,k), phi_grid(i,j,k))
        xdotth = thdotx(th_grid(i,j,k), phi_grid(i,j,k))
        xdotphi = phidotx(th_grid(i,j,k), phi_grid(i,j,k))
        zdotr = rdotz(th_grid(i,j,k), phi_grid(i,j,k))
        zdotth = thdotz(th_grid(i,j,k), phi_grid(i,j,k))
        zdotphi = phidotz(th_grid(i,j,k), phi_grid(i,j,k))

        ! Components !
        khat_r(i,j,k) = DSIN(theta_observe)*xdotr + DCOS(theta_observe)*zdotr
        khat_th(i,j,k) = DSIN(theta_observe)*xdotth + DCOS(theta_observe)*zdotth
        khat_phi(i,j,k) = DSIN(theta_observe)*xdotphi + DCOS(theta_observe)*zdotphi

      END DO
    END DO 
  END DO

ELSE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Then, setup electron number density !

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

        ! Electron number density, the absolute scale is not important 
        !IF(i == 1 .AND. j == j_th) THEN
          nrho(i,j,k) = 1.0D0
          gamma(i,j,k) = gamma_j
          beta_vel(i,j,k) = DSQRT(1.0d0 - 1.0d0/gamma(i,j,k)**2)
        !ELSE
        !  nrho(i,j,k) = 0.0D0
        !  gamma(i,j,k) = 1.0d0
        !  beta_vel(i,j,k) = 0.0D0
        !END IF

        ! Electron directional vector, orthonomral basis !
        nhat_r(i,j,k) = 0.0d0
        nhat_th(i,j,k) = 0.0d0
        nhat_phi(i,j,k) = 1.0d0

        ! Photon directional vector, orthonomral basis !
        xdotr = rdotx(th_grid(i,j,k), phi_grid(i,j,k))
        xdotth = thdotx(th_grid(i,j,k), phi_grid(i,j,k))
        xdotphi = phidotx(th_grid(i,j,k), phi_grid(i,j,k))
        zdotr = rdotz(th_grid(i,j,k), phi_grid(i,j,k))
        zdotth = thdotz(th_grid(i,j,k), phi_grid(i,j,k))
        zdotphi = phidotz(th_grid(i,j,k), phi_grid(i,j,k))

        ! Components !
        khat_r(i,j,k) = DSIN(theta_observe)*xdotr + DCOS(theta_observe)*zdotr
        khat_th(i,j,k) = DSIN(theta_observe)*xdotth + DCOS(theta_observe)*zdotth
        khat_phi(i,j,k) = DSIN(theta_observe)*xdotphi + DCOS(theta_observe)*zdotphi

      END DO
    END DO 
  END DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Setup the initial conditions for the grmhd model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INITIAL_GRMHD
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Real !
REAL*8 :: kdotr, kdotth, kdotphi
REAL*8 :: xdotr, xdotth, xdotphi
REAL*8 :: zdotr, zdotth, zdotphi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, load from HDF5 !
CALL HDF5_LOAD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Then, assign !
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nz

      ! one over gamma and beta !
      beta_vel(i,j,k) = DSQRT(1.0d0 - 1.0d0/gamma(i,j,k)**2)

      ! Photon directional vector, orthonomral basis !
      xdotr = rdotx(th_grid(i,j,k), phi_grid(i,j,k))
      xdotth = thdotx(th_grid(i,j,k), phi_grid(i,j,k))
      xdotphi = phidotx(th_grid(i,j,k), phi_grid(i,j,k))
      zdotr = rdotz(th_grid(i,j,k), phi_grid(i,j,k))
      zdotth = thdotz(th_grid(i,j,k), phi_grid(i,j,k))
      zdotphi = phidotz(th_grid(i,j,k), phi_grid(i,j,k))

      ! Components !
      khat_r(i,j,k) = DSIN(theta_observe)*xdotr + DCOS(theta_observe)*zdotr
      khat_th(i,j,k) = DSIN(theta_observe)*xdotth + DCOS(theta_observe)*zdotth
      khat_phi(i,j,k) = DSIN(theta_observe)*xdotphi + DCOS(theta_observe)*zdotphi

    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

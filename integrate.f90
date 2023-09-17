!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Integrate the whole space to obtain the scattering stokes parameter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTEGRATE
USE OMP_LIB
USE DEFINITION
IMPLICIT NONE

! Integer ! 
INTEGER :: n, i, j, k

! Angles !
REAL*8 :: sin_0, cos_0

! Angles !
REAL*8 :: cosz_0, sinz_0

! Angles !
REAL*8 :: cosrho_0, sinrho_0

! Angles !
REAL*8 :: cosrho0_p, sinrho0_p

! Angles !
REAL*8 :: cos2zeta, sin2zeta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Loop over the viewing angle !
DO n = 1, n_angle

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Print to debug !
  WRITE (*,*) 'step', n

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize !
  theta_view = theta_0(n)
  isc_total = 0.0d0
  qsc_total = 0.0d0
  usc_total = 0.0d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign !
  sin_0 = DSIN(theta_0(n))
  cos_0 = DCOS(theta_0(n))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We loop over the computational domain !
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) FIRSTPRIVATE(cos_0, sin_0, theta_view) &
  !$OMP PRIVATE(cosz_0, sinz_0, cosrho_0, sinrho_0, cosrho0_p, sinrho0_p, cos2zeta, sin2zeta) &
  !$OMP REDUCTION(+:isc_total, qsc_total, usc_total)
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angles !
        CALL RHO_0(nhat_x(i,j,k), nhat_z(i,j,k), cos_0, sin_0, cosrho_0, sinrho_0)
        CALL Z_0(nhat_z(i,j,k), cosz_0, sinz_0)
        CALL TWOZETA(phi_grid(i,j,k), cosz_0, cos_0, cosrho_0, sin_0, sinrho_0, cos2zeta, sin2zeta)
        CALL RHO0_P(cosrho_0, beta_vel(i,j,k), cosrho0_p, sinrho0_p)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angular contribution !
        IF(powerlaw) THEN
          CALL ANG_FAC_POWERLAW(cosrho0_p, cosrho_0, cos2zeta, sin2zeta, beta_vel(i,j,k), ang_fac_isc(i,j,k), ang_fac_qsc(i,j,k), ang_fac_usc(i,j,k))
        ELSEIF(custom) THEN
          ! Not yet implemented 
        END IF

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Sum up the contribution !
        isc_total = isc_total + ang_fac_isc(i,j,k)*nebar(i,j,k)*vol(i,j,k)*gammam1(i,j,k)**2
        qsc_total = qsc_total + ang_fac_qsc(i,j,k)*nebar(i,j,k)*vol(i,j,k)*gammam1(i,j,k)**2
        usc_total = usc_total + ang_fac_usc(i,j,k)*nebar(i,j,k)*vol(i,j,k)*gammam1(i,j,k)**2
      
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign total integrated stokes parameter
  stokes_total(1,n) = isc_total
  stokes_total(2,n) = qsc_total
  stokes_total(3,n) = usc_total
  pdegree(n) = (qsc_total/isc_total)*DSQRT(1.0d0 + (usc_total/qsc_total)**2)

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Loop over the viewing angle !
DO n = 1, n_angle

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize !
  theta_view = theta_0(n)
  isc_total = 0.0d0
  qsc_total = 0.0d0
  usc_total = 0.0d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign !
  sin_0 = DSIN(theta_0(n))
  cos_0 = DCOS(theta_0(n))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get map for minimum PD !  
  IF(pdegree(n) == MINVAL(pdegree)) THEN

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We loop over the computational domain !
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) FIRSTPRIVATE(cos_0, sin_0, theta_view) &
    !$OMP PRIVATE(cosz_0, sinz_0, cosrho_0, sinrho_0, cosrho0_p, sinrho0_p, cos2zeta, sin2zeta) &
    !$OMP REDUCTION(+:isc_total, qsc_total, usc_total)
    DO i = 1, nx
      DO j = 1, ny
        DO k = 1, nz

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Get angles !
          CALL RHO_0(nhat_x(i,j,k), nhat_z(i,j,k), cos_0, sin_0, cosrho_0, sinrho_0)
          CALL Z_0(nhat_z(i,j,k), cosz_0, sinz_0)
          CALL TWOZETA(phi_grid(i,j,k), cosz_0, cos_0, cosrho_0, sin_0, sinrho_0, cos2zeta, sin2zeta)
          CALL RHO0_P(cosrho_0, beta_vel(i,j,k), cosrho0_p, sinrho0_p)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Get angular contribution !
          IF(powerlaw) THEN
            CALL ANG_FAC_POWERLAW(cosrho0_p, cosrho_0, cos2zeta, sin2zeta, beta_vel(i,j,k), ang_fac_out(1,1,i,j,k), ang_fac_out(1,2,i,j,k), ang_fac_out(1,3,i,j,k))
          ELSEIF(custom) THEN
            ! Not yet implemented 
          END IF

        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get map for maximum PD !  
  ELSEIF(pdegree(n) == MAXVAL(pdegree)) THEN

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We loop over the computational domain !
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) FIRSTPRIVATE(cos_0, sin_0, theta_view) &
    !$OMP PRIVATE(cosz_0, sinz_0, cosrho_0, sinrho_0, cosrho0_p, sinrho0_p, cos2zeta, sin2zeta) &
    !$OMP REDUCTION(+:isc_total, qsc_total, usc_total)
    DO i = 1, nx
      DO j = 1, ny
        DO k = 1, nz

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Get angles !
          CALL RHO_0(nhat_x(i,j,k), nhat_z(i,j,k), cos_0, sin_0, cosrho_0, sinrho_0)
          CALL Z_0(nhat_z(i,j,k), cosz_0, sinz_0)
          CALL TWOZETA(phi_grid(i,j,k), cosz_0, cos_0, cosrho_0, sin_0, sinrho_0, cos2zeta, sin2zeta)
          CALL RHO0_P(cosrho_0, beta_vel(i,j,k), cosrho0_p, sinrho0_p)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Get angular contribution !
          IF(powerlaw) THEN
            CALL ANG_FAC_POWERLAW(cosrho0_p, cosrho_0, cos2zeta, sin2zeta, beta_vel(i,j,k), ang_fac_out(2,1,i,j,k), ang_fac_out(2,2,i,j,k), ang_fac_out(2,3,i,j,k))
          ELSEIF(custom) THEN
            ! Not yet implemented 
          END IF

        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO

  END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! HDF5 output !
CALL HDF5_OUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE
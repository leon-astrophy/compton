!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Assuming isoptropic and power law radiation, choose between full or approx calculation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTEGRATE_POWERLAW
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose !

! Choose between exact or head-on approx !
IF(exact) THEN
  CALL EXACT_POWERLAW
ELSEIF(headon) THEN
  CALL HEADON_POWERLAW
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Integrate the whole space to obtain the scattering stokes parameter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HEADON_POWERLAW
USE OMP_LIB
USE DEFINITION
IMPLICIT NONE

! Integer ! 
INTEGER :: n, i, j, k

! Integer !
INTEGER :: n_min, n_max

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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Loop over the viewing angle !
DO n = 1, n_angle

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Print to debug !
  WRITE (*,*) 'step', n

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize !
  theta_view = theta_0(n)
  isc_total = 0.0d0
  qsc_total = 0.0d0
  usc_total = 0.0d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign !
  sin_0 = DSIN(theta_0(n))
  cos_0 = DCOS(theta_0(n))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We loop over the computational domain !
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) FIRSTPRIVATE(cos_0, sin_0, theta_view) &
  !$OMP PRIVATE(cosz_0, sinz_0, cosrho_0, sinrho_0, cosrho0_p, sinrho0_p, cos2zeta, sin2zeta) &
  !$OMP REDUCTION(+:isc_total, qsc_total, usc_total)
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angles !
        CALL RHO_0(nhat_x(i,j,k), nhat_z(i,j,k), cos_0, sin_0, cosrho_0, sinrho_0)
        CALL Z_0(nhat_z(i,j,k), cosz_0, sinz_0)
        CALL TWOZETA(phi_grid(i,j,k), cosz_0, cos_0, cosrho_0, sin_0, sinrho_0, cos2zeta, sin2zeta)
        CALL RHO0_P(cosrho_0, beta_vel(i,j,k), cosrho0_p, sinrho0_p)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angular contribution !
        CALL ANGFAC_POWERLAW_HEADON(cosrho0_p, cosrho_0, cos2zeta, sin2zeta, beta_vel(i,j,k), ang_fac_isc(i,j,k), ang_fac_qsc(i,j,k), ang_fac_usc(i,j,k))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Sum up the contribution !
        gam_fac(i,j,k) = gammam1(i,j,k)**2
        isc_total = isc_total + ang_fac_isc(i,j,k)*nebar(i,j,k)*vol(i,j,k)*gam_fac(i,j,k)
        qsc_total = qsc_total + ang_fac_qsc(i,j,k)*nebar(i,j,k)*vol(i,j,k)*gam_fac(i,j,k)
        usc_total = usc_total + ang_fac_usc(i,j,k)*nebar(i,j,k)*vol(i,j,k)*gam_fac(i,j,k)
      
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Scale up !
  isc_total = isc_total*1.5d0*(2.0d0**(s_power)/(2.0d0+s_power))
  qsc_total = qsc_total*1.5d0*(2.0d0**(s_power)/(2.0d0+s_power))
  usc_total = usc_total*1.5d0*(2.0d0**(s_power)/(2.0d0+s_power))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign total integrated stokes parameter
  stokes_total(1,n) = isc_total
  stokes_total(2,n) = qsc_total
  stokes_total(3,n) = usc_total
  pdegree(n) = (qsc_total/isc_total)

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Look for minimum and maximum PD position angle !
IF(full) THEN

  ! Locate minimum and maximum index !
  n_min = MINLOC(pdegree, DIM=1)
  n_max = MAXLOC(pdegree, DIM=1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Look for angular factor at minimum !

  ! Initialize !
  n = n_min
  theta_view = theta_0(n)
  isc_total = 0.0d0
  qsc_total = 0.0d0
  usc_total = 0.0d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign !
  sin_0 = DSIN(theta_0(n))
  cos_0 = DCOS(theta_0(n))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We loop over the computational domain !
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) FIRSTPRIVATE(cos_0, sin_0, theta_view) &
  !$OMP PRIVATE(cosz_0, sinz_0, cosrho_0, sinrho_0, cosrho0_p, sinrho0_p, cos2zeta, sin2zeta) &
  !$OMP REDUCTION(+:isc_total, qsc_total, usc_total)
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angles !
        CALL RHO_0(nhat_x(i,j,k), nhat_z(i,j,k), cos_0, sin_0, cosrho_0, sinrho_0)
        CALL Z_0(nhat_z(i,j,k), cosz_0, sinz_0)
        CALL TWOZETA(phi_grid(i,j,k), cosz_0, cos_0, cosrho_0, sin_0, sinrho_0, cos2zeta, sin2zeta)
        CALL RHO0_P(cosrho_0, beta_vel(i,j,k), cosrho0_p, sinrho0_p)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angular contribution !
        CALL ANGFAC_POWERLAW_HEADON(cosrho0_p, cosrho_0, cos2zeta, sin2zeta, beta_vel(i,j,k), ang_fac_out(1,1,i,j,k), ang_fac_out(1,2,i,j,k), ang_fac_out(1,3,i,j,k))
        
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Look for angular factor at maximum !

  ! Initialize !
  n = n_max
  theta_view = theta_0(n)
  isc_total = 0.0d0
  qsc_total = 0.0d0
  usc_total = 0.0d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign !
  sin_0 = DSIN(theta_0(n))
  cos_0 = DCOS(theta_0(n))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We loop over the computational domain !
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) FIRSTPRIVATE(cos_0, sin_0, theta_view) &
  !$OMP PRIVATE(cosz_0, sinz_0, cosrho_0, sinrho_0, cosrho0_p, sinrho0_p, cos2zeta, sin2zeta) &
  !$OMP REDUCTION(+:isc_total, qsc_total, usc_total)
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angles !
        CALL RHO_0(nhat_x(i,j,k), nhat_z(i,j,k), cos_0, sin_0, cosrho_0, sinrho_0)
        CALL Z_0(nhat_z(i,j,k), cosz_0, sinz_0)
        CALL TWOZETA(phi_grid(i,j,k), cosz_0, cos_0, cosrho_0, sin_0, sinrho_0, cos2zeta, sin2zeta)
        CALL RHO0_P(cosrho_0, beta_vel(i,j,k), cosrho0_p, sinrho0_p)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angular contribution !
        CALL ANGFAC_POWERLAW_HEADON(cosrho0_p, cosrho_0, cos2zeta, sin2zeta, beta_vel(i,j,k), ang_fac_out(2,1,i,j,k), ang_fac_out(2,2,i,j,k), ang_fac_out(2,3,i,j,k))
        
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Integrate the whole space to obtain the scattering stokes parameter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EXACT_POWERLAW
USE OMP_LIB
USE DEFINITION
IMPLICIT NONE

! Integer ! 
INTEGER :: n, i, j, k

! Integer !
INTEGER :: n_min, n_max

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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Loop over the viewing angle !
DO n = 1, n_angle

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Print to debug !
  WRITE (*,*) 'step', n

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize !
  theta_view = theta_0(n)
  isc_total = 0.0d0
  qsc_total = 0.0d0
  usc_total = 0.0d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign !
  sin_0 = DSIN(theta_0(n))
  cos_0 = DCOS(theta_0(n))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We loop over the computational domain !
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) FIRSTPRIVATE(cos_0, sin_0, theta_view) &
  !$OMP PRIVATE(cosz_0, sinz_0, cosrho_0, sinrho_0, cosrho0_p, sinrho0_p, cos2zeta, sin2zeta) &
  !$OMP REDUCTION(+:isc_total, qsc_total, usc_total)
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angles !
        CALL RHO_0(nhat_x(i,j,k), nhat_z(i,j,k), cos_0, sin_0, cosrho_0, sinrho_0)
        CALL Z_0(nhat_z(i,j,k), cosz_0, sinz_0)
        CALL TWOZETA(phi_grid(i,j,k), cosz_0, cos_0, cosrho_0, sin_0, sinrho_0, cos2zeta, sin2zeta)
        CALL RHO0_P(cosrho_0, beta_vel(i,j,k), cosrho0_p, sinrho0_p)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angular contribution !
        CALL ANGFAC_POWERLAW_FULL(cosrho0_p, sinrho0_p, cosrho_0, cos2zeta, sin2zeta, beta_vel(i,j,k), ang_fac_isc(i,j,k), ang_fac_qsc(i,j,k), ang_fac_usc(i,j,k))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Sum up the contribution !
        gam_fac(i,j,k) = gammam1(i,j,k)**(6.0d0 + 2.0d0*s_power)
        isc_total = isc_total + ang_fac_isc(i,j,k)*nebar(i,j,k)*vol(i,j,k)*gam_fac(i,j,k)
        qsc_total = qsc_total + ang_fac_qsc(i,j,k)*nebar(i,j,k)*vol(i,j,k)*gam_fac(i,j,k)
        usc_total = usc_total + ang_fac_usc(i,j,k)*nebar(i,j,k)*vol(i,j,k)*gam_fac(i,j,k)
      
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Scale up !
  isc_total = isc_total*3.0d0/(16.0d0*pi)
  qsc_total = qsc_total*3.0d0/(16.0d0*pi)
  usc_total = usc_total*3.0d0/(16.0d0*pi)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign total integrated stokes parameter
  stokes_total(1,n) = isc_total
  stokes_total(2,n) = qsc_total
  stokes_total(3,n) = usc_total
  pdegree(n) = (qsc_total/isc_total)

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Look for minimum and maximum PD position angle !
IF(full) THEN

  ! Locate minimum and maximum index !
  n_min = MINLOC(pdegree, DIM=1)
  n_max = MAXLOC(pdegree, DIM=1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Look for angular factor at minimum pd

  ! Initialize
  n = n_min
  theta_view = theta_0(n)
  isc_total = 0.0d0
  qsc_total = 0.0d0
  usc_total = 0.0d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign !
  sin_0 = DSIN(theta_0(n))
  cos_0 = DCOS(theta_0(n))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We loop over the computational domain !
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) FIRSTPRIVATE(cos_0, sin_0, theta_view) &
  !$OMP PRIVATE(cosz_0, sinz_0, cosrho_0, sinrho_0, cosrho0_p, sinrho0_p, cos2zeta, sin2zeta) &
  !$OMP REDUCTION(+:isc_total, qsc_total, usc_total)
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angles !
        CALL RHO_0(nhat_x(i,j,k), nhat_z(i,j,k), cos_0, sin_0, cosrho_0, sinrho_0)
        CALL Z_0(nhat_z(i,j,k), cosz_0, sinz_0)
        CALL TWOZETA(phi_grid(i,j,k), cosz_0, cos_0, cosrho_0, sin_0, sinrho_0, cos2zeta, sin2zeta)
        CALL RHO0_P(cosrho_0, beta_vel(i,j,k), cosrho0_p, sinrho0_p)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angular contribution !
        CALL ANGFAC_POWERLAW_FULL(cosrho0_p, sinrho0_p, cosrho_0, cos2zeta, sin2zeta, beta_vel(i,j,k), ang_fac_out(1,1,i,j,k), ang_fac_out(1,2,i,j,k), ang_fac_out(1,3,i,j,k))
      
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Look for angular factor at minimum pd

  ! Initialize !
  n = n_max
  theta_view = theta_0(n)
  isc_total = 0.0d0
  qsc_total = 0.0d0
  usc_total = 0.0d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign !
  sin_0 = DSIN(theta_0(n))
  cos_0 = DCOS(theta_0(n))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We loop over the computational domain !
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) FIRSTPRIVATE(cos_0, sin_0, theta_view) &
  !$OMP PRIVATE(cosz_0, sinz_0, cosrho_0, sinrho_0, cosrho0_p, sinrho0_p, cos2zeta, sin2zeta) &
  !$OMP REDUCTION(+:isc_total, qsc_total, usc_total)
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angles !
        CALL RHO_0(nhat_x(i,j,k), nhat_z(i,j,k), cos_0, sin_0, cosrho_0, sinrho_0)
        CALL Z_0(nhat_z(i,j,k), cosz_0, sinz_0)
        CALL TWOZETA(phi_grid(i,j,k), cosz_0, cos_0, cosrho_0, sin_0, sinrho_0, cos2zeta, sin2zeta)
        CALL RHO0_P(cosrho_0, beta_vel(i,j,k), cosrho0_p, sinrho0_p)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angular contribution !
        CALL ANGFAC_POWERLAW_FULL(cosrho0_p, sinrho0_p, cosrho_0, cos2zeta, sin2zeta, beta_vel(i,j,k), ang_fac_out(2,1,i,j,k), ang_fac_out(2,2,i,j,k), ang_fac_out(2,3,i,j,k))
      
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate angular factor in the polarization integral, for the headon approx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ANGFAC_POWERLAW_HEADON(cosrho0p_in, cosrho0_in, cos2zeta_in, sin2zeta_in, betavel_in, isc_out, qsc_out, usc_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: cosrho0p_in, cosrho0_in, cos2zeta_in, sin2zeta_in, betavel_in

! Output !
REAL*8, INTENT(OUT) :: isc_out, qsc_out, usc_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate !

isc_out = (1.0d0 + cosrho0p_in**2)/(1.0d0 - betavel_in*cosrho0_in)**(2.0d0 + s_power)
qsc_out = (1.0d0 - cosrho0p_in**2)*cos2zeta_in/(1.0d0 - betavel_in*cosrho0_in)**(2.0d0 + s_power)
usc_out = (1.0d0 - cosrho0p_in**2)*sin2zeta_in/(1.0d0 - betavel_in*cosrho0_in)**(2.0d0 + s_power)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate angular factor in the polarization integral, for full integral
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ANGFAC_POWERLAW_FULL(cosrho0p_in, sinrho0p_in, cosrho0_in, cos2zeta_in, sin2zeta_in, betavel_in, isc_out, qsc_out, usc_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: cosrho0p_in, sinrho0p_in, cosrho0_in, cos2zeta_in, sin2zeta_in, betavel_in

! Output !
REAL*8, INTENT(OUT) :: isc_out, qsc_out, usc_out

! Real !
REAL*8 :: iprime_out, qprime_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate !

CALL DBINTEGRAL_POWERLAW(cosrho0p_in, sinrho0p_in, cosrho0_in, betavel_in, iprime_out, qprime_out)
isc_out = iprime_out/(1.0d0 - betavel_in*cosrho0_in)**(2.0d0 + s_power)
qsc_out = qprime_out*cos2zeta_in/(1.0d0 - betavel_in*cosrho0_in)**(2.0d0 + s_power)
usc_out = qprime_out*sin2zeta_in/(1.0d0 - betavel_in*cosrho0_in)**(2.0d0 + s_power)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Fortran subroutine for performing double integral
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DBINTEGRAL_POWERLAW(cosrho0p_in, sinrho0p_in, cosrho0_in, betavel_in, iprime_out, qprime_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: cosrho0p_in, sinrho0p_in, cosrho0_in, betavel_in
REAL*8, INTENT(OUT) :: iprime_out, qprime_out

! INTEGER !
INTEGER :: i, j, k

! REAL !
REAL*8 :: uqsc_out
REAL*8 :: fin_m1, fin_c, fin_p1
REAL*8 :: fout_m1, fout_c, fout_p1
REAL*8 :: gin_m1, gin_c, gin_p1
REAL*8 :: gout_m1, gout_c, gout_p1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! evaluate the integral !

! Initialize !
iprime_out = 0
qprime_out = 0

! Loop over the domain !
DO i = 1, n_integral, 2

  ! Initialize !
  fout_m1 = 0
  fout_c = 0
  fout_p1 = 0
  gout_m1 = 0
  gout_c = 0
  gout_p1 = 0

  ! Loop over the domain !
  DO j = 1, n_integral, 2

    ! For the scattering intensity !
    CALL INTEGRAND_IP_POWERLAW(cosalpha_p(j-1), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, betavel_in, fin_m1)
    CALL INTEGRAND_IP_POWERLAW(cosalpha_p(j), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, betavel_in, fin_c)
    CALL INTEGRAND_IP_POWERLAW(cosalpha_p(j+1), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, betavel_in, fin_p1)
    fout_m1 = fout_m1 + dalphap*(fin_m1 + 4.0d0*fin_c + fin_p1)/3.0d0

    CALL INTEGRAND_IP_POWERLAW(cosalpha_p(j-1), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, betavel_in, fin_m1)
    CALL INTEGRAND_IP_POWERLAW(cosalpha_p(j), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, betavel_in, fin_c)
    CALL INTEGRAND_IP_POWERLAW(cosalpha_p(j+1), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, betavel_in, fin_p1)
    fout_c = fout_c + dalphap*(fin_m1 + 4.0d0*fin_c + fin_p1)/3.0d0

    CALL INTEGRAND_IP_POWERLAW(cosalpha_p(j-1), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, betavel_in, fin_m1)
    CALL INTEGRAND_IP_POWERLAW(cosalpha_p(j), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, betavel_in, fin_c)
    CALL INTEGRAND_IP_POWERLAW(cosalpha_p(j+1), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, betavel_in, fin_p1)
    fout_p1 = fout_p1 + dalphap*(fin_m1 + 4.0d0*fin_c + fin_p1)/3.0d0

    ! For the stokes U and Q !
    CALL INTEGRAND_QP_POWERLAW(alphap_dir(j-1), cosalpha_p(j-1), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, betavel_in, gin_m1)
    CALL INTEGRAND_QP_POWERLAW(alphap_dir(j), cosalpha_p(j), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, betavel_in, gin_c)
    CALL INTEGRAND_QP_POWERLAW(alphap_dir(j+1), cosalpha_p(j+1), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, betavel_in, gin_p1)
    gout_m1 = gout_m1 + dalphap*(gin_m1 + 4.0d0*gin_c + gin_p1)/3.0d0

    CALL INTEGRAND_QP_POWERLAW(alphap_dir(j-1), cosalpha_p(j-1), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, betavel_in, gin_m1)
    CALL INTEGRAND_QP_POWERLAW(alphap_dir(j), cosalpha_p(j), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, betavel_in, gin_c)
    CALL INTEGRAND_QP_POWERLAW(alphap_dir(j+1), cosalpha_p(j+1), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, betavel_in, gin_p1)
    gout_c = gout_c + dalphap*(gin_m1 + 4.0d0*gin_c + gin_p1)/3.0d0    

    CALL INTEGRAND_QP_POWERLAW(alphap_dir(j-1), cosalpha_p(j-1), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, betavel_in, gin_m1)
    CALL INTEGRAND_QP_POWERLAW(alphap_dir(j), cosalpha_p(j), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, betavel_in, gin_c)
    CALL INTEGRAND_QP_POWERLAW(alphap_dir(j+1), cosalpha_p(j+1), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, betavel_in, gin_p1)
    gout_p1 = gout_p1 + dalphap*(gin_m1 + 4.0d0*gin_c + gin_p1)/3.0d0  

  END DO
  
  ! Do the outer integration !
  iprime_out = iprime_out + drhop*(fout_m1 + 4.0d0*fout_c + fout_p1)/3.0d0
  qprime_out = qprime_out + drhop*(gout_m1 + 4.0d0*gout_c + gout_p1)/3.0d0

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computing the inner integrand of the scattering intensity function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTEGRAND_IP_POWERLAW(cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, betavel_in, int_out)
USE DEFINITION
IMPLICIT NONE

! REAL !
REAL*8, INTENT(IN) :: cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, betavel_in
REAL*8, INTENT(OUT) :: int_out

! REAL !
REAL*8 :: cosw_p, sinw_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
CALL W_P (cosalphap_in, cosrho0p_in, sinrho0p_in, cosrhop_in, sinrhop_in, cosw_p, sinw_p)

! Assign !
int_out = sinrhop_in*(1.0d0+cosw_p**2)/(1.0d0+betavel_in*cosrhop_in)**(3+s_power)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computing the inner integrand of the Stokes Q and U parameter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTEGRAND_QP_POWERLAW(alphap_in, cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, betavel_in, int_out)
USE DEFINITION
IMPLICIT NONE

! REAL !
REAL*8, INTENT(IN) :: alphap_in, cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, betavel_in
REAL*8, INTENT(OUT) :: int_out

! REAL !
REAL*8 :: cosw_p, sinw_p
REAL*8 :: cos2eta, sin2eta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
CALL W_P (cosalphap_in, cosrho0p_in, sinrho0p_in, cosrhop_in, sinrhop_in, cosw_p, sinw_p)
CALL TWOETA(alphap_in, cosrhop_in, cosrho0p_in, sinrho0p_in, cosw_p, sinw_p, cos2eta, sin2eta)

! Assign !
int_out = sinrhop_in*(1.0d0-cosw_p**2)*cos2eta/(1.0d0+betavel_in*cosrhop_in)**(3+s_power)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE

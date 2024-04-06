!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Assuming isoptropic and power law radiation, choose between full or approx calculation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTEGRATE
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose !

! Choose between exact or head-on approx !
IF(exact) THEN
  CALL EXACT_ANGFAC
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
USE DEFINITION
IMPLICIT NONE

! Integer ! 
INTEGER :: n, i, j, k

! Integer !
INTEGER :: n_min, n_max

! Angles !
REAL*8 :: sin_0, cos_0

! Angles !
REAL*8 :: cosrho_0, sinrho_0

! Angles !
REAL*8 :: cosrho0_p, sinrho0_p

! Angles !
REAL*8 :: cos2zeta, sin2zeta

! Directional angle for electron !
REAL*8 :: cos_thetae, edoty, edotz
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Loop over the viewing angle !
DO n = 1, n_angle

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Print to debug !
  WRITE (*,*) 'step', n

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize !
  theta_view = theta_0(n)
  isctaui0 = 0.0d0
  qsctaui0 = 0.0d0
  usctaui0 = 0.0d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign !
  sin_0 = DSIN(theta_0(n))
  cos_0 = DCOS(theta_0(n))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We loop over the computational domain !
  !!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) FIRSTPRIVATE(cos_0, sin_0, theta_view) &
  !!$OMP PRIVATE(cos_thetae, cosrho_0, sinrho_0, cosrho0_p, sinrho0_p, cos2zeta, sin2zeta) &
  !!$OMP REDUCTION(+:isc_total, qsc_total, usc_total)
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angles !
        CALL THETA_E(th_grid(i,j,k), phi_grid(i,j,k), nhat_r(i,j,k), nhat_th(i,j,k), nhat_phi(i,j,k), cos_thetae, edoty, edotz)
        CALL RHO_0(khat_r(i,j,k), khat_th(i,j,k), khat_phi(i,j,k), nhat_r(i,j,k), nhat_th(i,j,k), nhat_phi(i,j,k), cosrho_0, sinrho_0)
        CALL TWOZETA(edoty, edotz, cos_0, cosrho_0, sin_0, sinrho_0, cos2zeta, sin2zeta)
        CALL RHO0_P(cosrho_0, beta_vel(i,j,k), cosrho0_p, sinrho0_p)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angular contribution !
        CALL ANGFAC_POWERLAW_HEADON(cosrho0_p, cosrho_0, cos2zeta, sin2zeta, beta_vel(i,j,k), fsc(i,j,k), gsc(i,j,k), hsc(i,j,k))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get the scattering factor !
        fsc(i,j,k) = fsc(i,j,k)*1.5d0*(2.0d0**(s_power)/(2.0d0+s_power))/gamma(i,j,k)**2
        gsc(i,j,k) = gsc(i,j,k)*1.5d0*(2.0d0**(s_power)/(2.0d0+s_power))/gamma(i,j,k)**2
        hsc(i,j,k) = hsc(i,j,k)*1.5d0*(2.0d0**(s_power)/(2.0d0+s_power))/gamma(i,j,k)**2

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Sum it up !
        isctaui0 = isctaui0 + nrho(i,j,k)*fsc(i,j,k)/DBLE(nphi)
        qsctaui0 = qsctaui0 + nrho(i,j,k)*gsc(i,j,k)/DBLE(nphi)
        usctaui0 = usctaui0 + nrho(i,j,k)*hsc(i,j,k)/DBLE(nphi)

      END DO
    END DO
  END DO
  !!$OMP END PARALLEL DO

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Integrate the whole space to obtain the scattering stokes parameter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EXACT_ANGFAC
USE DEFINITION
IMPLICIT NONE

! Integer ! 
INTEGER :: n, i, j, k

! Integer !
INTEGER :: n_min, n_max

! Angles !
REAL*8 :: sin_0, cos_0

! Angles !
REAL*8 :: cosrho_0, sinrho_0

! Angles !
REAL*8 :: cosrho0_p, sinrho0_p

! Angles !
REAL*8 :: cos2zeta, sin2zeta

! Directional angle for electron !
REAL*8 :: cos_thetae, edoty, edotz
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Loop over the viewing angle !
DO n = 1, n_angle

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Print to debug !
  WRITE (*,*) 'step', n

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize !
  theta_view = theta_0(n)
  isctaui0 = 0.0d0
  qsctaui0 = 0.0d0
  usctaui0 = 0.0d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Assign !
  sin_0 = DSIN(theta_0(n))
  cos_0 = DCOS(theta_0(n))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We loop over the computational domain !
  !!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) FIRSTPRIVATE(cos_0, sin_0, theta_view) &
  !!$OMP PRIVATE(cos_thetae, cosrho_0, sinrho_0, cosrho0_p, sinrho0_p, cos2zeta, sin2zeta) &
  !!$OMP REDUCTION(+:isc_total, qsc_total, usc_total)
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angles !
        CALL THETA_E(th_grid(i,j,k), phi_grid(i,j,k), nhat_r(i,j,k), nhat_th(i,j,k), nhat_phi(i,j,k), cos_thetae, edoty, edotz)
        CALL RHO_0(khat_r(i,j,k), khat_th(i,j,k), khat_phi(i,j,k), nhat_r(i,j,k), nhat_th(i,j,k), nhat_phi(i,j,k), cosrho_0, sinrho_0)
        CALL TWOZETA(edoty, edotz, cos_0, cosrho_0, sin_0, sinrho_0, cos2zeta, sin2zeta)
        CALL RHO0_P(cosrho_0, beta_vel(i,j,k), cosrho0_p, sinrho0_p)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get angular contribution !
        CALL ANGFAC_FULL(cosrho0_p, sinrho0_p, cosrho_0, cos2zeta, sin2zeta, gamma(i,j,k), beta_vel(i,j,k), fsc(i,j,k), gsc(i,j,k), hsc(i,j,k))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Get the scattering factor !
        fsc(i,j,k) = fsc(i,j,k)*3.0d0/(16.0d0*pi)
        gsc(i,j,k) = gsc(i,j,k)*3.0d0/(16.0d0*pi)
        hsc(i,j,k) = hsc(i,j,k)*3.0d0/(16.0d0*pi)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Sum it up !
        IF(powerlaw) THEN
          isctaui0 = isctaui0 + nrho(i,j,k)*fsc(i,j,k)/DBLE(nphi)
          qsctaui0 = qsctaui0 + nrho(i,j,k)*gsc(i,j,k)/DBLE(nphi)
          usctaui0 = usctaui0 + nrho(i,j,k)*hsc(i,j,k)/DBLE(nphi)
        END IF

      END DO
    END DO
  END DO
  !!$OMP END PARALLEL DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Scale it down by Cv-s
  IF(powerlaw) THEN
    isctaui0 = isctaui0/(c_power*nu_light**(-s_power))
    qsctaui0 = qsctaui0/(c_power*nu_light**(-s_power))
    usctaui0 = usctaui0/(c_power*nu_light**(-s_power))
  END IF
  
END DO

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

SUBROUTINE ANGFAC_FULL(cosrho0p_in, sinrho0p_in, cosrho0_in, cos2zeta_in, sin2zeta_in, gamma_in, betavel_in, isc_out, qsc_out, usc_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: cosrho0p_in, sinrho0p_in, cosrho0_in, cos2zeta_in, sin2zeta_in, gamma_in, betavel_in

! Output !
REAL*8, INTENT(OUT) :: isc_out, qsc_out, usc_out

! Real !
REAL*8 :: iprime_out, qprime_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate !

CALL DBINTEGRAL(cosrho0p_in, sinrho0p_in, cosrho0_in, gamma_in, betavel_in, iprime_out, qprime_out)
isc_out = iprime_out/gamma_in**(3.0d0)/(1.0d0 - betavel_in*cosrho0_in)**(2.0d0)
qsc_out = qprime_out*cos2zeta_in/gamma_in**(3.0d0)/(1.0d0 - betavel_in*cosrho0_in)**(2.0d0)
usc_out = qprime_out*sin2zeta_in/gamma_in**(3.0d0)/(1.0d0 - betavel_in*cosrho0_in)**(2.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Fortran subroutine for performing double integral
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DBINTEGRAL(cosrho0p_in, sinrho0p_in, cosrho0_in, gamma_in, betavel_in, iprime_out, qprime_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: cosrho0p_in, sinrho0p_in, cosrho0_in, gamma_in, betavel_in
REAL*8, INTENT(OUT) :: iprime_out, qprime_out

! INTEGER !
INTEGER :: i, j, k

! REAL !
REAL*8 :: uqsc_out
REAL*8 :: fin_m1, fin_c, fin_p1
REAL*8 :: fout_m1, fout_c, fout_p1
REAL*8 :: gin_m1, gin_c, gin_p1
REAL*8 :: gout_m1, gout_c, gout_p1

! REAL ! 
REAL*8 :: d1_m1, d1_c, d1_p1
REAL*8 :: nu_m1, nu_c, nu_p1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

  ! Compute doppler shifted frequency and the d1 factor !
  CALL doppler_d1(gamma_in, betavel_in, cosrho_p(i-1), d1_m1)
  CALL doppler_d1(gamma_in, betavel_in, cosrho_p(i), d1_c)
  CALL doppler_d1(gamma_in, betavel_in, cosrho_p(i+1), d1_p1)
  CALL frequency_shift(gamma_in, betavel_in, cosrho0_in, d1_m1, nu_m1)
  CALL frequency_shift(gamma_in, betavel_in, cosrho0_in, d1_c, nu_c)
  CALL frequency_shift(gamma_in, betavel_in, cosrho0_in, d1_p1, nu_p1)

  ! Loop over the domain !
  DO j = 1, n_integral, 2

    ! For the scattering intensity !
    CALL INTEGRAND_IP(cosalpha_p(j-1), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, fin_m1)
    CALL INTEGRAND_IP(cosalpha_p(j), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, fin_c)
    CALL INTEGRAND_IP(cosalpha_p(j+1), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, fin_p1)
    fout_m1 = fout_m1 + dalphap*(fin_m1 + 4.0d0*fin_c + fin_p1)/3.0d0

    CALL INTEGRAND_IP(cosalpha_p(j-1), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, d1_c, nu_c, fin_m1)
    CALL INTEGRAND_IP(cosalpha_p(j), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, d1_c, nu_c, fin_c)
    CALL INTEGRAND_IP(cosalpha_p(j+1), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, d1_c, nu_c, fin_p1)
    fout_c = fout_c + dalphap*(fin_m1 + 4.0d0*fin_c + fin_p1)/3.0d0

    CALL INTEGRAND_IP(cosalpha_p(j-1), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, fin_m1)
    CALL INTEGRAND_IP(cosalpha_p(j), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, fin_c)
    CALL INTEGRAND_IP(cosalpha_p(j+1), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, fin_p1)
    fout_p1 = fout_p1 + dalphap*(fin_m1 + 4.0d0*fin_c + fin_p1)/3.0d0

    ! For the stokes U and Q !
    CALL INTEGRAND_QP(alphap_dir(j-1), cosalpha_p(j-1), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, gin_m1)
    CALL INTEGRAND_QP(alphap_dir(j), cosalpha_p(j), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, gin_c)
    CALL INTEGRAND_QP(alphap_dir(j+1), cosalpha_p(j+1), cosrho_p(i-1), sinrho_p(i-1), cosrho0p_in, sinrho0p_in, d1_m1, nu_m1, gin_p1)
    gout_m1 = gout_m1 + dalphap*(gin_m1 + 4.0d0*gin_c + gin_p1)/3.0d0

    CALL INTEGRAND_QP(alphap_dir(j-1), cosalpha_p(j-1), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, d1_c, nu_c, gin_m1)
    CALL INTEGRAND_QP(alphap_dir(j), cosalpha_p(j), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, d1_c, nu_c, gin_c)
    CALL INTEGRAND_QP(alphap_dir(j+1), cosalpha_p(j+1), cosrho_p(i), sinrho_p(i), cosrho0p_in, sinrho0p_in, d1_c, nu_c, gin_p1)
    gout_c = gout_c + dalphap*(gin_m1 + 4.0d0*gin_c + gin_p1)/3.0d0    

    CALL INTEGRAND_QP(alphap_dir(j-1), cosalpha_p(j-1), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, gin_m1)
    CALL INTEGRAND_QP(alphap_dir(j), cosalpha_p(j), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, gin_c)
    CALL INTEGRAND_QP(alphap_dir(j+1), cosalpha_p(j+1), cosrho_p(i+1), sinrho_p(i+1), cosrho0p_in, sinrho0p_in, d1_p1, nu_p1, gin_p1)
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

SUBROUTINE INTEGRAND_IP(cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, d1_in, nu_in, int_out)
USE DEFINITION
IMPLICIT NONE

! REAL !
REAL*8, INTENT(IN) :: cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, d1_in, nu_in
REAL*8, INTENT(OUT) :: int_out

! REAL !
REAL*8 :: cosw_p, sinw_p, inu_shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
CALL W_P (cosalphap_in, cosrho0p_in, sinrho0p_in, cosrhop_in, sinrhop_in, cosw_p, sinw_p)

! Get intensity !
CALL RADIATION_FIELD(nu_in, inu_shift)

! Assign !
int_out = sinrhop_in*(1.0d0+cosw_p**2)*inu_shift*d1_in**(3.0d0) !/(1.0d0+betavel_in*cosrhop_in)**(3+s_power)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computing the inner integrand of the Stokes Q and U parameter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTEGRAND_QP(alphap_in, cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, d1_in, nu_in, int_out)
USE DEFINITION
IMPLICIT NONE

! REAL !
REAL*8, INTENT(IN) :: alphap_in, cosalphap_in, cosrhop_in, sinrhop_in, cosrho0p_in, sinrho0p_in, d1_in, nu_in
REAL*8, INTENT(OUT) :: int_out

! REAL !
REAL*8 :: inu_shift
REAL*8 :: cosw_p, sinw_p
REAL*8 :: cos2eta, sin2eta
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
CALL W_P (cosalphap_in, cosrho0p_in, sinrho0p_in, cosrhop_in, sinrhop_in, cosw_p, sinw_p)
CALL TWOETA(alphap_in, cosrhop_in, cosrho0p_in, sinrho0p_in, cosw_p, sinw_p, cos2eta, sin2eta)

! Get intensity !
CALL RADIATION_FIELD(nu_in, inu_shift)

! Assign !
int_out = sinrhop_in*(1.0d0-cosw_p**2)*cos2eta*inu_shift*d1_in**(3.0d0)  !/(1.0d0+betavel_in*cosrhop_in)**(3+s_power)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computing doppler shifted photon frequency 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE frequency_shift(gamma_in, betavel_in, cosrho0_in, d1_in, nu_shift)
USE DEFINITION
IMPLICIT NONE

! REAL !
REAL*8, INTENT(IN) :: gamma_in, betavel_in, cosrho0_in, d1_in
REAL*8, INTENT(OUT) :: nu_shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
nu_shift = nu_light*gamma_in*(1.0d0 - betavel_in*cosrho0_in)/d1_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computing the doppler factor d1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE doppler_d1(gamma_in, betavel_in, cosrhop_in, d1_out)
USE DEFINITION
IMPLICIT NONE

! REAL !
REAL*8, INTENT(IN) :: gamma_in, betavel_in, cosrhop_in
REAL*8, INTENT(OUT) :: d1_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
d1_out = 1.0d0/(gamma_in*(1.0d0 + betavel_in*cosrhop_in))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE
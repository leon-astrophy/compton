!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate cosrho_0 and sinrho_0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RHO_0(nhatx_in, nhatz_in, cos0_in, sin0_in, cosrho0_out, sinrho0_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: nhatx_in, nhatz_in, cos0_in, sin0_in

! Output !
REAL*8, INTENT(OUT) :: cosrho0_out, sinrho0_out

! Angles !
REAL*8 :: angle

! Calculate 
cosrho0_out = nhatx_in*sin0_in + nhatz_in*cos0_in
cosrho0_out = MIN(cosrho0_out, 1.0d0)
cosrho0_out = MAX(cosrho0_out, -1.0d0)
sinrho0_out = DSQRT(1.0d0 - cosrho0_out**2)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate cosz_0 and sinz_0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Z_0(nhatz_in, cosz0_out, sinz0_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: nhatz_in

! Output !
REAL*8, INTENT(OUT) :: cosz0_out, sinz0_out

! Calculate 
cosz0_out = nhatz_in
cosz0_out = MIN(cosz0_out, 1.0d0)
cosz0_out = MAX(cosz0_out, -1.0d0)
sinz0_out = DSQRT(1.0d0 - cosz0_out**2)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate cos2zeta and sin2zeta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TWOZETA(phie_in, cosz0_in, cos0_in, cosrho0_in, sin0_in, sinrho0_in, cos2zeta_out, sin2zeta_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: phie_in, cosz0_in, cos0_in, cosrho0_in, sin0_in, sinrho0_in

! Output !
REAL*8, INTENT(OUT) :: cos2zeta_out, sin2zeta_out

! Local !
REAL*8 :: coszeta, sinzeta, zeta, num, den

! Calculate 
num = (cosz0_in - cos0_in*cosrho0_in)
den = sin0_in*sinrho0_in
IF(den == 0) THEN
  coszeta = 0.0d0
ELSE
  coszeta = num/den
END IF
coszeta = MIN(coszeta, 1.0d0)
coszeta = MAX(coszeta, -1.0d0)
zeta = DACOS(coszeta)

! Angle !
IF(phie_in >= pi) THEN
  zeta = 2.0D0*PI - zeta
END IF
sinzeta = DSIN(zeta)

! Get angles !
cos2zeta_out = 2.0d0*coszeta**2 - 1.0d0
sin2zeta_out = 2.0d0*sinzeta*coszeta

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate cosrho0_p and sinrho0_p
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RHO0_P(cosrho0_in, betavel_in, cosrho0p_out, sinrho0p_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: cosrho0_in, betavel_in

! Output !
REAL*8, INTENT(OUT) :: cosrho0p_out, sinrho0p_out

! Local !
REAL*8 :: coszeta, sinzeta, num, den

! Calculate 
num = cosrho0_in - betavel_in
den = 1.0d0 - betavel_in*cosrho0_in
IF(den == 0) THEN
  cosrho0p_out = 0.0d0
ELSE
  cosrho0p_out = num/den
END IF        
cosrho0p_out = MIN(cosrho0p_out, 1.0d0)
cosrho0p_out = MAX(cosrho0p_out, -1.0d0)
sinrho0p_out = DSQRT(1.0d0 - cosrho0p_out**2)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate cow_p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE W_P(cosalphap_in, cosrho0p_in, sinrho0p_in, cosrhop_in, sinrhop_in, coswp_out, sinwp_out)
USE DEFINITION
IMPLICIT NONE
  
! Input !
REAL*8, INTENT(IN) :: cosalphap_in, cosrho0p_in, sinrho0p_in, cosrhop_in, sinrhop_in

! Output !
REAL*8, INTENT(OUT) :: coswp_out, sinwp_out
  
! Calculate 
coswp_out = cosalphap_in*sinrhop_in*sinrho0p_in + cosrhop_in*cosrho0p_in
coswp_out = MIN(coswp_out, 1.0d0)
coswp_out = MAX(coswp_out, -1.0d0)
sinwp_out = DSQRT(1.0d0 - coswp_out**2)
  
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate cos2eta and sin2eta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TWOETA(alphap_in, cosrhop_in, cosrho0p_in, sinrho0p_in, coswp_in, sinwp_in, cos2eta_out, sin2eta_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: alphap_in, cosrhop_in, cosrho0p_in, sinrho0p_in, coswp_in, sinwp_in

! Output !
REAL*8, INTENT(OUT) :: cos2eta_out, sin2eta_out

! Local !
REAL*8 :: coseta, sineta, eta, num, den

! Calculate 
num = (cosrhop_in - cosrho0p_in*coswp_in)
den = sinrho0p_in*sinwp_in
IF(den == 0) THEN
  coseta = 0.0d0
ELSE
  coseta = num/den
END IF
coseta = MIN(coseta, 1.0d0)
coseta = MAX(coseta, -1.0d0)
eta = DACOS(coseta)

! Angle !
IF(alphap_in >= pi) THEN
  eta = 2.0D0*PI - eta
END IF
sineta = DSIN(eta)

! Get angles !
cos2eta_out = 2.0d0*coseta**2 - 1.0d0
sin2eta_out = 2.0d0*sineta*coseta

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
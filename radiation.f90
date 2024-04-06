!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Choose radiation law according to user input 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RADIATION_FIELD(nu_in, inu_out)
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: nu_in

! Output !
REAL*8, INTENT(OUT) :: inu_out

! Choose !
IF(powerlaw) THEN
  inu_out = c_power*nu_in**(-s_power)
ELSEIF(bbody) THEN
  inu_out = bbody_factor*nu_in**(3.0d0)/(DEXP(nu_in*hkt) - 1.0d0)
END IF

END SUBROUTINE
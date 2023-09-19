!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Choose the model that used to integrate to get the scattering stokes parameters
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTEGRATE
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose according to parameter.h !

! Isotropic plus power-law radiation or custom spectrum !
IF(powerlaw) THEN
  CALL INTEGRATE_POWERLAW
ELSEIF(custom) THEN

  ! Not yet implemetned !
  STOP 'not yet implemented'
  !CALL INTEGRATE_CUSTOM

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ouput after integration !

! HDF5 output !
CALL HDF5_OUT

END SUBROUTINE
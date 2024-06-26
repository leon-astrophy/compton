!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Computing polarization of an isotropic background radiation field
! using the formalism by Begelman et al. 1987 and Dexter et al. 2023 
! Written by H.S. Leon Chan on fall 2023
! JILA and APS, University of Colorado
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM XRAYPOL
#ifdef DCPU
USE OMP_LIB
#endif
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set openmp threads !

#ifdef DCPU
CALL OMP_SET_NUM_THREADS(64)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get terminal input (if any)
CALL TERMINAL

! Setup the computational domain !
CALL SETUP

! Build arrays !
CALL ARRAY

! Then, setup initial conditions !
CALL INITIAL

! Do the integraiton !
CALL INTEGRATE

! Ouput after integration !
CALL HDF5_OUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM

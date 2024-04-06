!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Parameter file containing various parameters for xraypol
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! pi constant 
REAL*8, PARAMETER :: pi = 4.D0*DATAN(1.D0)

! define small number to avoid coordinate singularity !
REAL*8, PARAMETER :: small_num = TINY(1.0D0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Integration density !
INTEGER :: n_integral = 32

! Number of grid along three axis !
INTEGER :: nr = 1
INTEGER :: nth = 1
INTEGER :: nphi = 256

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Electron distribution gamma for dirac-delta function !
REAL*8 :: gamma_j = 2.0d0

! Electron distribution angle for dirac-delta function, in DEGREE !
REAL*8 :: theta_j = 90.0D0

! Electron distribution angle for blob model, in DEGREE !
REAL*8 :: phi_j = 270.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Index for power law radiation field !
REAL*8, PARAMETER :: s_power = 1.0d0

! Proportionality constant for power law radiation field !
REAL*8, PARAMETER :: c_power = 1.0d0

! Target photon frequency (Log10 Hz) !
REAL*8 :: nu_light = DLOG10(1.0D18)

! Black body temperature (Log10 Kelvin) !
REAL*8 :: bbody_temp = DLOG10(1.0D9)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Target viewing angle IN DEGREE !
REAL*8 :: theta_observe = 90.0d0

! Number of angles !
INTEGER :: n_angle = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Output map? !
LOGICAL, PARAMETER :: map = .true.

! Toy disc model? !
LOGICAL, PARAMETER :: disc = .false.

! Electron model, choose from analytic or GRMHD !
LOGICAL, PARAMETER :: analytic = .true.
LOGICAL, PARAMETER :: grmhd = .false.

! Integration mode, head-one approximation or full calculation !
LOGICAL, PARAMETER :: headon = .false.
LOGICAL, PARAMETER :: exact = .true.

! The radiation field, power law or black body !
! Black body radiation is available only for exact calculation !
LOGICAL, PARAMETER :: powerlaw = .true.
LOGICAL, PARAMETER :: bbody = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! GRMHD snapshot name !
CHARACTER (len = 999) :: grmhd_file = 'pp-torus.h5'

! Ouput file name !
CHARACTER (len = 999) :: output_file = 'data'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

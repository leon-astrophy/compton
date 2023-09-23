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
INTEGER :: n_integral = 100

! Number of grid along three axis !
INTEGER :: nr = 288
INTEGER :: nth = 128
INTEGER :: nphi = 128

! Viewing angle density !
INTEGER :: n_angle = 100

! Electron distribution angle for dirac-delta function !
REAL*8, PARAMETER :: theta_j = 0.6d0

! Electron distribution gamma for dirac-delta function !
REAL*8, PARAMETER :: gamma_j = 2.0d0

! Index for power law radiation field !
REAL*8, PARAMETER :: s_power = 1.0d0

! Target photon frequency (Hz) !
REAL*8 :: nu_light = 1.0D17

! Target viewing angle !
REAL*8 :: theta_in = pi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculating full polar viewing angle or discrete angles !
LOGICAL, PARAMETER :: discrete = .true.
LOGICAL, PARAMETER :: full = .false.

! Electron model, choose from analytic or GRMHD !
LOGICAL, PARAMETER :: analytic = .true.
LOGICAL, PARAMETER :: grmhd = .false.

! Integration mode, head-one approximation or full calculation !
LOGICAL, PARAMETER :: headon = .false.
LOGICAL, PARAMETER :: exact = .true.

! The radiation field, power law or customized !
LOGICAL, PARAMETER :: powerlaw = .true.
LOGICAL, PARAMETER :: custom = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! GRMHD snapshot name !
CHARACTER (len = 999) :: grmhd_file = 'pp-torus.h5'

! Ouput file name !
CHARACTER (len = 999) :: output_file = 'data'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

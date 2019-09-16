module parameters
  implicit none
  ! Precision
  integer, parameter                    :: dp = kind(0.d0)
  integer, parameter                    :: sp = kind(0.0) 
  integer, parameter                    :: wp = dp ! Choose presiscion, dp or sp
  ! integer, parameter                    :: long_double = selected_real_kind(16)

  ! Natural constants as.
  real(wp), parameter                   :: PI     = acos(-1._wp)
  real(wp), parameter                   :: planck = 6.62607004e-34_wp
  real(wp), parameter                   :: redpla = planck/(2._wp*PI)
  real(wp), parameter                   :: me     = 9.10938356e-9_wp
  real(wp), parameter                   :: c      = 3e8_wp
  real(wp), parameter                   :: eps    = 1.0e-14_wp
  real(wp), parameter                   :: eps15  = 1.0e-15_wp
  complex(wp), parameter                :: i      = (0._wp, 1._wp)
  integer(wp), parameter                :: prog   = 10 ! Progress report frequency
  
  ! Parameters
  real(wp), parameter                   :: width   = 5.0_wp ! Width in reduced units
  real(wp), parameter                   :: height  = 8.0_wp  ! Height in reduced units
  real(wp), parameter                   :: lambda  = 500e-3_wp ! e-3 is eqv to e-9 meter using red. units
  real(wp), parameter                   :: problim = 1e-6_wp  ! Ray termination limmit
  real(wp), parameter                   :: nrlim   = 1e-13_wp ! Newton-Raphson limit
  integer(wp), parameter                :: PSOSlim = 0
  integer(wp), parameter                :: bouncelim = 900000000
  integer(wp), parameter                :: boxcntlim = 5000
  real(wp), parameter                   :: n0      = 1.0_wp
  real(wp)                              :: n1      = 1.9_wp  ! ellipse
  real(wp), parameter                   :: n2      = 2.0_wp  ! film
  character(len=4), parameter           :: potentialtype = "1p1c" ! n points, m lines = npml
  real(wp), parameter                   :: ni      = 0.0054_wp

  ! Initial parameters
  real(wp), parameter                   :: xi      = 0.0_wp+eps
  real(wp), parameter                   :: xf      = width/2._wp
  !                                                  0.1234567891ABCDEF | A ruler for number of decimals
  real(wp), parameter                   :: iang    = -89.99_wp
  real(wp), parameter                   :: fang    =  89.99_wp
  ! real(wp), parameter                   :: initn   = n0
  ! character(len=4), parameter           :: initloc = "C   "
  real(wp), parameter                   :: initn   = n0
  character(len=4), parameter           :: initloc = "C   "

  real(wp), parameter                   :: iskew   = 1.01_wp
  real(wp), parameter                   :: fskew   = 2.5_wp
  
  real(wp), parameter                   :: inang   = 0.0_wp ! In degrees
  real(wp), parameter                   :: xloc    = 0.5_wp
  integer(wp), parameter                :: nrays   = 10000
  integer(wp), parameter                :: vars    = 41 ! Shape variations
  integer(wp), parameter                :: randomlines = 14
  real(wp), parameter                   :: roughamp = 0.08_wp ! Amplitude of random structure
  real(wp), parameter                   :: normalizedthickness = 0.4_wp ! Double check width and height
  integer(wp), parameter                :: skiprandom = 24
  integer(wp), parameter                :: inst = 300 ! Instances of one type
  
  integer(wp), parameter                :: specshap = 40 ! Specify a specific "random" shape
  integer(wp), parameter                :: nshapes  = vars
  
  ! BCM
  integer(wp), parameter                :: grid = 500

  ! Switches
  logical                               :: recording = .false.
  logical, parameter                    :: optical   = .true. ! if not optical it's geometric
  logical, parameter                    :: weighted  = .true. ! probability weigths
  logical, parameter                    :: raysplitting = .true.
  logical, parameter                    :: storephasspc = .false.
  logical, parameter                    :: periodicbcs  = .false.
  logical, parameter                    :: trackrays    = .true.
  logical, parameter                    :: absorption   = .false.
  logical, parameter                    :: savefractal  = .true.
  logical, parameter                    :: addredundancy= .false.
  logical, parameter                    :: startagrid   = .false.
  integer(wp), parameter                :: specray = 1508218 ! Store a specific ray
  integer(wp), parameter                :: savepaths = 100+100 ! number of actual saved paths is minus 100

  ! Globals
  integer(wp)                           :: erays = 0
  integer(wp)                           :: deserter = 0
  integer(wp)                           :: PSOSrays = 0
  integer(wp)                           :: globaltime = 0
  real(wp)                              :: probloss, truncray
  real(wp)                              :: truncloss = 0._wp
  
  ! Define Poincare surface of section
  real(wp), dimension(2), parameter     :: PSOSnorm = [0,1] ! Normal pointing into the bucket
  character(len=4), parameter           :: PSOS = "B   "

  ! addBounce limits
  integer(wp), parameter    :: clim1 = 3
  integer(wp), parameter    :: clim2 = 10
  integer(wp), parameter    :: clim3 = 20
  integer(wp), parameter    :: clim4 = 40

  ! Do not change, these are initial values
  integer(wp)               :: splittings = 100 ! Do not change
  integer(wp)               :: nextid  = 100 ! Do not change
  integer(wp), parameter    :: pid     = 99 ! Do not change
  integer(wp), parameter    :: oplid   = 98 ! Do not change
  integer(wp), parameter    :: effid   = 97 ! Do not change
  integer(wp), parameter    :: cplid   = 96 ! Do not change
  integer(wp), parameter    :: varid   = 95 ! Do not change

end module parameters

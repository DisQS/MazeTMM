module RNG_RLFSR113
  USE SystemParameters
  implicit none
  
  ! accessibility
  private
  public :: rlfsr113
  public :: lfsrinit
  
  ! variables
  integer(kind=ikind) :: z1, z2, z3, z4
 
  ! parameter
  real(kind=rkind),    parameter :: AM = 4.656612873077d-10
  integer(kind=ikind), parameter :: IA = 16807
  integer(kind=ikind), parameter :: IM = 2147483647
  integer(kind=ikind), parameter :: IQ = 127773
  integer(kind=ikind), parameter :: IR = 2836
  
contains

  ! -----------------------------------------------------------------
  ! rlfsr113()  returns a random number of interval (0, 1)
  ! -----------------------------------------------------------------
  function rlfsr113() result(dRet)
    real(kind=rkind)        :: dRet

    integer(kind=ikind) :: b

    b  = ishft(ieor(ishft(z1,6),z1),-13)
    z1 = ieor(ishft(iand(z1,-2),18),b)

    b  = ishft(ieor(ishft(z2,2),z2),-27)
    z2 = ieor(ishft(iand(z2,-8),2),b)

    b  = ishft(ieor(ishft(z3,13),z3),-21)
    z3 = ieor(ishft(iand(z3,-16),7),b)

    b  = ishft(ieor(ishft(z4,3),z4),-12)
    z4 = ieor(ishft(iand(z4,-128),13),b)

    dRet=ishft(ieor(ieor(ieor(z1,z2),z3),z4),-1)*AM

  end function rlfsr113

  ! -----------------------------------------------------------------
  ! lfsrinit()  initialize rlfsr113 (z1,z2,z3,z4) 
  ! -----------------------------------------------------------------
  subroutine lfsrinit(idum)
    integer(kind=ikind), intent(inout) :: idum

    integer(kind=ikind) :: k,c1,c2,c3,c4
    real(kind=rkind)    :: dummy
    

    ! Check whether the FORTRAN integers can be used as unsigned long !

    ! data c1 /B'11111111111111111111111111111110'/
    !  data c1 /X'FFFFFFFE'/
    c1 = Z"FFFFFFFE"
    ! data c2 /B'11111111111111111111111111111000'/
    !  data c2 /X'FFFFFFF8'/
    c2 = Z"FFFFFFF8"
    ! data c3 /B'11111111111111111111111111110000'/
    !  data c3 /X'FFFFFFF0'/
    c3 = Z"FFFFFFF0"
    ! data c4 /B'11111111111111111111111110000000'/
    !  data c4 /X'FFFFFF80'/
    c4 = Z"FFFFFF80"
    
    if ((c1.ne.-2).or.(c2.ne.-8).or.(c3.ne.-16).or.(c4.ne.-128)) then
      print *,"c1,c2,c3,c4", c1,c2,c3,c4
      print *,'Nonstandard integer representation. Stoped.'
      stop
    endif

    ! Initialize z1,z2,z3,z4

    if (idum.le.0) idum=1
    k=(idum)/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum = idum + IM
    if (idum.lt.2) then
      z1=idum+2 
    else 
      z1=idum
    endif
    k=(idum)/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum = idum + IM
    if (idum.lt.8) then 
      z2=idum+8 
    else 
      z2=idum
    endif
    k=(idum)/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum = idum + IM
    if (idum.lt.16) then
      z3=idum+16 
    else 
      z3=idum
    endif
    k=(idum)/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum = idum + IM
    if (idum.lt.128) then
      z4=idum+128 
    else 
      z4=idum
    endif

    ! Make a single call to rlfsr113() to achieve a valid state
    dummy=rlfsr113()
      
  end subroutine lfsrinit

end module RNG_RLFSR113


module RNG
  use RNG_RLFSR113
  use SystemParameters

  implicit none

  ! accessibility
  private
  public :: SRANDOM
  public :: DRANDOM 
  public :: GRANDOM
  
  ! parameter

!  ! kind parameter for double precision
!  integer, parameter :: PRECISION = 8
  ! number of random number inside [0,1] for routine gauss
  integer, parameter :: GAUSS_N   = 20
contains

  ! ------------------------------------------------------------------
  ! SRANDOM() 
  !
  ! Random number generator SEED interface for use with any old RND
  ! ------------------------------------------------------------------
  subroutine SRANDOM( ISeed )
    integer, intent(in) :: ISeed

    integer idum

    ! change following lines to incorporate different RND generators
    idum = ISeed
    call lfsrinit(idum)
  
  end subroutine SRANDOM

  ! ------------------------------------------------------------------
  ! DRANDOM() 
  !
  ! Random number generator interface for use with any old RND
  ! ------------------------------------------------------------------
  function DRANDOM( ISeed ) result(dRet)
    integer, intent(in) :: ISeed
    real(kind=rkind):: dRet

    ! change following lines to incorporate different RND generators
    dRet = rlfsr113()  ! NOTE that ISeed is never used

  end function DRANDOM

  ! ------------------------------------------------------------------
  ! GRANDOM() 
  !
  ! Gaussian random number generator interface
  ! ------------------------------------------------------------------
  function GRANDOM( ISeed, avg, sigma) result (dRet)
    integer,          intent(in) :: ISeed
    real(kind=rkind), intent(in) :: avg
    real(kind=rkind), intent(in) :: sigma
    real(kind=rkind)             :: dRet
    
  
    ! change following lines to incorporate different RND generators
    call gauss(dRet, sigma, avg) ! NOTE that ISeed is never used
  
  end function GRANDOM


! TODO
  ! ------------------------------------------------------------------
  ! GAUSS()
  !
  ! THE ROUTINE GAUSS GENERATES A RANDOM NUMBER
  ! IN A GAUSSIAN DISTRIBUTION
  !
  ! VARIABLES:
  !
  !   X     - THE OUTPUT GAUSSIAN VARIABLE
  !   sigma - standard deviation
  !   mu    - average
  ! 
  ! NOTE: The random number generator rlfsr113 should be 
  !   initialised by calling
  !   the subroutine lfsrinit
  ! ------------------------------------------------------------------
  subroutine gauss(X,sigma,mu)
    real(kind=rkind), intent(out):: X
    real(kind=rkind), intent(in) :: sigma
    real(kind=rkind), intent(in) :: mu

    real(kind=rkind) :: Y, SUM
    integer i

    SUM=0.0
    do i=1,GAUSS_N
      Y=rlfsr113()
      Y=2.d0*(Y-0.5d0)
      SUM=SUM+Y
    end do

    X=mu+sigma*SUM* DSQRT( 3.0d0 /DBLE(GAUSS_N) )
   
  end subroutine gauss

end module RNG


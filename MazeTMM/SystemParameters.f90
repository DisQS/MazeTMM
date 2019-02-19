MODULE SystemParameters

  IMPLICIT NONE

  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,307)
  REAL(KIND=RKIND), PARAMETER  :: ZERO=0.0, ONE=1.0
  REAL(KIND=RKIND), PARAMETER  :: TINY=1.0D-15

  INTEGER(KIND=IKIND) :: Width ! The width of cross section for 2D systems, for 3D systerms the cross sction is in widthXwidth
  INTEGER(KIND=IKIND) :: Dim ! The dimension of the system
  INTEGER(KIND=IKIND) :: WidthSquared
  INTEGER(KIND=IKIND) :: NofOrtho !the number of TMM steps before doing orthogonalization
  INTEGER(KIND=IKIND) :: SubSteps !the numbner of substeps in a TMM step
  INTEGER(KIND=IKIND) :: MaxOrthoNo=10000000 !If not convegent, the TMM runs MaxOrthoNo times orthogonalization at most
  INTEGER(KIND=IKIND) :: MazeType
  
  
  CHARACTER :: BoundaryCon ! Boundary condition in non-transformed direction, "O" for open, "P" for periodic, "A" for anti-periodic
  CHARACTER :: LoopDirection ! Determine the program loop type, "E" means running energy loop, "D" means running disorder loop
  
  
  REAL(KIND=RKIND)  :: ConvCriterion  !TMM convergence criterion
  REAL(KIND=RKIND)  :: DiagDis, Energy0, Energy1, dEnergy ! for the energy loop
  REAL(KIND=RKIND)  :: Energy, DiagDis0, DiagDis1, dDiagDis ! for the disorder loop
  REAL(KIND=RKIND)  :: MazePotential
  INTEGER(KIND=IKIND) :: Conv !output parameter to indicate the convegence of a TMM
  
    
    
END MODULE SystemParameters
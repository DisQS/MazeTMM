PROGRAM TMMxD
    USE SystemParameters
    USE util
    USE RNG
    USE io_module

    IMPLICIT NONE
    
    REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE :: PsiA, PsiB, Maze
    REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: Gammas, GammasSquared, AveGammas,Variance
    INTEGER(KIND=IKIND):: PrintSwitch=1 ! controlling the print
    INTEGER(KIND=IKIND):: Index, IndexIn, LayerNo,Seed=71277
    REAL(KIND=RKIND):: Flux0,Flux1,dFlux,Flux, tick, tock
    
    !Read input parameters
    NAMELIST /SystemSettings/ Width, Dim, BoundaryCon, LoopDirection,MazeType,MazePotential
    NAMELIST /RunControl/ ConvCriterion, NofOrtho, SubSteps
    NAMELIST /EnergyLoop/ Energy0, Energy1, dEnergy, DiagDis
    NAMELIST /DisorderLoop/ DiagDis0, DiagDis1, dDiagDis, Energy
    CALL openUnit('Parameters.nml',138,'O')
    READ(138,SystemSettings)
    READ(138,RunControl)
    READ(138,EnergyLoop)
    READ(138,DisorderLoop)
    CLOSE(138)
    
    !only 3D open boundary condition and odd length case is implemented!
    IF ((3/=Dim).OR.('O'/=BoundaryCon).OR.(1/=mod(Width,2))) THEN
        PRINT *,'Only 3D open boundary condition and odd length case is implemented!'
        STOP
    END IF
    
    !Begin timing
    CALL CPU_TIME(tick)
    
    !Load Maze
    CALL AllocateMat(Maze,Width)
    CALL LoadMaze(Maze,Width,MazeType,MazePotential)

    ! Check and Print parameters
    IF ((Dim .NE. 2).AND.(Dim .NE. 3)) THEN
        PRINT *, 'Dimension ',Dim,' is not allowed!'
        STOP
    END IF
    IF ((BoundaryCon .NE. 'O') .AND. (BoundaryCon .NE. 'P') .AND. (BoundaryCon .NE. 'A')) THEN
        PRINT *, 'Boundary condition ',BoundaryCon,' is not allowed!'
        STOP
    END IF
    IF ((LoopDirection .NE. 'E') .AND. (LoopDirection .NE. 'D')) THEN
        PRINT *, 'Loop Direction ',LoopDirection,' is not allowed!'
        STOP
    END IF
    IF (PrintSwitch) THEN
        PRINT *, 'Beginning A Transform Matrix Case study!'
        PRINT '(A11,I1,A8,I2,A20,A1,A10,I2,A15,F6.1)', ' Dimension=',Dim,', Width=',Width, ' Boundary Condition=',BoundaryCon,' MazeType=',MazeType,' MazePotential=',MazePotential
        PRINT '(A15,E7.1,A12,I2,A30,I3)', ' ConvCriterion=',ConvCriterion,' SubSteps=',SubSteps,' Steps in each ortho=',NofOrtho
        IF ('E'==LoopDirection) THEN
            PRINT '(A26,F6.1,A4,F6.1,A11,F6.1,A22,F6.1)', ' TMM changing energy from ',Energy0,' to ',Energy1,' with step ',dEnergy,' and disorder strength ',DiagDis
        ELSE
            PRINT '(A28,F6.1,A4,F6.1,A11,F6.1,A12,F6.1)', ' TMM changing disorder from ',DiagDis0,' to ',DiagDis1,' with step ',dDiagDis,' and energy ',Energy
        END IF
        PRINT *,' '
        CALL ShowMaze(Maze,Width)
        PRINT *,' '
    END IF
    
    !Allocate memeory
    WidthSquared=Width*Width ! only needed in 3D case
    SELECT CASE (Dim)
    CASE(2)
        CALL AllocateMat(PsiA,Width)
        CALL AllocateMat(PsiB,Width)
        CALL AllocateVec(Gammas,Width)
        CALL AllocateVec(GammasSquared,Width)
        CALL AllocateVec(AveGammas,Width)
        CALL AllocateVec(Variance,Width)
    CASE(3)
        CALL AllocateMat(PsiA,WidthSquared)
        CALL AllocateMat(PsiB,WidthSquared)
        CALL AllocateVec(Gammas,WidthSquared)
        CALL AllocateVec(GammasSquared,WidthSquared)
        CALL AllocateVec(AveGammas,WidthSquared)
        CALL AllocateVec(Variance,WidthSquared) 
    CASE DEFAULT
        STOP
    END SELECT
    
    !Set up the loop
    SELECT CASE(LoopDirection)
    CASE('D')
        Flux0=DiagDis0
        Flux1=DiagDis1
        dFlux=dDiagDis
    CASE('E')
        Flux0=Energy0
        Flux1=Energy1
        dFlux=dEnergy
    CASE DEFAULT
        STOP
    ENDSELECT
    
    !set up the random number generator
    CALL SRANDOM(Seed)
    
    !Flux loop
    DO Flux=Flux0,Flux1,dFlux
        
        IF (LoopDirection=='D') THEN
            DiagDis=Flux
        ELSE
            Energy=Flux
        END IF
        IF (PrintSwitch) THEN
            PRINT '(A17, F6.1,A9,F6.1,A12)', ' Start Disorder ',DiagDis,' Energy ',Energy,' TMM study!'
        END IF
        
        !reset vectors and matrices
        Conv=" "
        PsiA=ZERO
        PsiB=ZERO
        Gammas=ZERO
        GammasSquared=ZERO
        AveGammas=ZERO
        Variance=ZERO
        IF (Dim==2) THEN
            DO Index=1,Width
                PsiA(Index,Index)= ONE
            END DO
        ELSE    
            DO Index=1,WidthSquared
                PsiA(Index,Index)= ONE
            END DO 
        END IF
        
        !TMM loops
        DO Index=1,MaxOrthoNo
            Do IndexIn=1, (NofOrtho*SubSteps), SubSteps
                
                LayerNo=(Index-1)*NofOrtho*SubSteps+IndexIn
                IF (Dim==2) THEN
                    PRINT *, '2D TMM is not implemented!'
                    STOP
                    !TMMMaze2D(PsiA,PsiB,LayerNo,Energy,DiagDis,Width)
                    !TMMMaze2D(PsiB,PsiA,LayerNo,Energy,DiagDis,Width)
                ELSE
                    CALL TMMMaze3D(PsiA,PsiB,Maze,LayerNo,Energy,DiagDis)
                    CALL TMMMaze3D(PsiB,PsiA,Maze,LayerNo+1,Energy,DiagDis)
                END IF
            END DO !IndexIn loop
            
            ! renormalize via Gram-Schmidt
            CALL ReNorm(PsiA,PsiB,Gammas,GammasSquared,WidthSquared)
            !sort the eigenvalues by LARGEST first AND also resort the eigenvectors accordingly
            CALL ReSort(PsiA,PsiB,Gammas,GammasSquared,WidthSquared)
            ! do the gamma computations
            CALL DoGamma(Gammas,GammasSquared,AveGammas,Variance,WidthSquared,Index,NofOrtho,SubSteps)
            !Wrtite AveGammas and Variance to file for later check of the convegence process
            !CALL WriteAveGammasAndVarToFile(Index,AveGammas(1),Variance(1))
            
            IF (PrintSwitch) THEN
                PRINT '(A20, F6.2,A10,F6.2,A14,I8,A23,F8.4,A8,F8.4)', ' TMM With Disorder ', DiagDis,' Energy ',Energy,' finish the ',Index,' ortho with avegamma ',AveGammas(1),' var ',Variance(1)
            END IF
            
            ! check accuracy 
            IF((Index>5) .AND. (Variance(1) .LE. ConvCriterion) .AND. (Variance(1).GE.TINY)) THEN
                 Conv=1
                 EXIT
            ENDIF
            
        END DO !TMM loop
        !If not convergent
        IF (Index>MaxOrthoNo) THEN
            PRINT '(A20, F6.2,A10,F6.2,A20)', ' TMM for disorder ',DiagDis,' Energy ',Energy,' is not convegent!'
            Conv=0
        END IF
        
        !Write results to file
        CALL WriteOutput(DiagDis,Energy,AveGammas,Variance,Index,Conv,PsiA,WidthSquared)
        
        
        IF (PrintSwitch) THEN
            PRINT '(A20, F6.2,A10,F6.2,A12)', ' TMM for disorder ',DiagDis,' Energy ',Energy,' is done!'
        END IF
         
    END DO !Flux loop
    
    DEALLOCATE(PsiA,PsiB,Gammas,GammasSquared,AveGammas,Variance,Maze)
    
    IF (PrintSwitch) THEN
        PRINT *, ' TMM is done!'
    END IF
    CALL CPU_TIME(tock)
    IF (PrintSwitch) THEN
        PRINT *, 'Time taken=',(tock-tick)/3600.0,' hours!'
    END IF
    
    PAUSE
END PROGRAM TMMxD
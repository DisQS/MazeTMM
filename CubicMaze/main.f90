!Use each forward and backword block Gamma to judge convergence
!Store every sample, but not average samples with negetive LE or error larger than 2 
!Update 22 May 2019
PROGRAM CubicMazeTMM
    USE SystemParameters
    USE RNG
    USE util
    USE io_module
    
    IMPLICIT NONE
    
    REAL(KIND=RKIND):: temp,tick,tock, Flux0,Flux1,dFlux,Flux
    INTEGER(KIND=IKIND):: Length, High, Width, NoOfOrtho, Realization, tempIndex, Seed=71277
    REAL(KIND=RKIND), DIMENSION(:,:,:), ALLOCATABLE :: Disorder, Maze
    REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE :: PsiA, PsiB
    REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: GammasPerSweep, OldGammas, accVariance
    REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: AveGammas,AveGammasSquared,StdErr,ReletiveErr, accVarianceSum
    INTEGER(KIND=IKIND):: PrintSwitch=1 ! controlling the print
    INTEGER(KIND=IKIND)::MazeWallNo!Store the wall lattice number in a maze
    REAL(KIND=RKIND)::MeanMazeWallRatio, CondG_PerSample, OldCondG_PerSample, accCondG_PerSample, accCondG_Sum, AveCondG, AveCondGSquare, StdErrCondG,ReletiveErrCondG 
    INTEGER(KIND=IKIND) :: Index, IndexFB, IndexTMM,OrthoedNo,EffSamples
    INTEGER(KIND=IKIND) :: Conv !output parameter to indicate the convegence of a TMM
    INTEGER(KIND=IKIND) :: MaxSweepNo=100 !If not convegent, the TMM runs MaxOrthoNo times orthogonalization at most
    
    INTEGER(KIND=IKIND) :: NoOfCell !The 3D systerms is in NoOfCell^3 form
    INTEGER(KIND=IKIND) :: OrthoSteps !the number of TMM steps before doing orthogonalization 
    INTEGER(KIND=IKIND) :: MazeType
    INTEGER(KIND=IKIND) :: NofDisConfig !the number of the disorder realizations
    CHARACTER :: LoopDirection ! Determine the program loop type, "E" means running energy loop, "D" means running disorder loop
    REAL(KIND=RKIND)  :: ConvCriterion  !TMM convergence criterion
    REAL(KIND=RKIND)  :: Ratio !the ratio of disorder strengths between the sites in path and wall
    REAL(KIND=RKIND)  :: DiagDis, Energy0, Energy1, dEnergy ! for the energy loop
    REAL(KIND=RKIND)  :: Energy, DiagDis0, DiagDis1, dDiagDis ! for the disorder loop
    
    !Read input parameters
    NAMELIST /SystemSettings/ NoOfCell, LoopDirection, MazeType, Ratio
    NAMELIST /RunControl/ ConvCriterion, OrthoSteps, NofDisConfig
    NAMELIST /EnergyLoop/ Energy0, Energy1, dEnergy, DiagDis
    NAMELIST /DisorderLoop/ DiagDis0, DiagDis1, dDiagDis, Energy
    CALL openUnit('parameters.nml',138,'O')
    READ(138,SystemSettings)
    READ(138,RunControl)
    READ(138,EnergyLoop)
    READ(138,DisorderLoop)
    CLOSE(138)
    
    
    SELECT CASE(MazeType)
    CASE(1)
        Width=NoOfCell*2-1
        High=NoOfCell*2-1
        Length=NoOfCell*2-1
    CASE DEFAULT
        PRINT *, 'Type ',MazeType,' has not implemented!'
        STOP
    END SELECT
    
    IF (mod(Length,OrthoSteps)>(OrthoSteps/2)) THEN
        NoOfOrtho=Length/OrthoSteps+1
    ELSE
        NoOfOrtho=Length/OrthoSteps
    END IF
        
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
    
    CALL AllocateR3Tensor(Maze,Width,High,Length)
    CALL AllocateR3Tensor(Disorder,Width,High,Length)
    CALL AllocateMat(PsiA,Width*High,Width*High)
    CALL AllocateMat(PsiB,Width*High,Width*High)
    CALL AllocateVec(GammasPerSweep,Width*High)
    CALL AllocateVec(OldGammas,Width*High)
    CALL AllocateVec(accVariance,Width*High)
    CALL AllocateVec(AveGammas,Width*High)
    CALL AllocateVec(AveGammasSquared,Width*High)
    CALL AllocateVec(StdErr,Width*High)
    CALL AllocateVec(ReletiveErr,Width*High)
    CALL AllocateVec(accVarianceSum,Width*High)

    
    
    !Flux loop
    DO Flux=Flux0,Flux1,dFlux
        
        IF (LoopDirection=='D') THEN
            DiagDis=Flux
        ELSE
            Energy=Flux
        END IF
        IF (1==PrintSwitch) THEN
            PRINT '(A17, F6.1,A9,F6.1,A12)', ' Start Disorder ',DiagDis,' Energy ',Energy,' TMM study!'
        END IF
        
        CALL CPU_TIME(tick)
        EffSamples=0;
        AveGammas=ZERO
        AveGammasSquared=ZERO
        accVarianceSum=ZERO
        MeanMazeWallRatio=ZERO
        accCondG_Sum=ZERO
        AveCondG=ZERO
        AveCondGSquare=ZERO
        !Different Disorder realization loop
        DO Realization=1,NofDisConfig
            
            !Reset parameters
            Conv=0
            MazeWallNo=0
            PsiA=ZERO
            PsiB=ZERO
            OldGammas=ZERO
            OldCondG_PerSample=ZERO
            DO tempIndex=1,Width*High
                PsiA(tempIndex,tempIndex)=ONE
            END DO
            !set up the random number generator
            CALL SRANDOM(Seed+(Realization-1)*100)
            
            CALL Generate3DMaze(Maze,NoOfCell,MazeType, Width, High, Length,MazeWallNo)
            CALL GenerateDisorder(Disorder,DiagDis,Maze,Width,High,Length,Ratio)
            
            !Convegence loop
            DO Index=1,MaxSweepNo
                GammasPerSweep=ZERO
                ! Forward and backward loop
                DO IndexFB=1,2 
                    !TMM loop
                    OrthoedNo=0
                    DO IndexTMM=1,(Length-1),2
                        CALL TMM3D(PsiA,PsiB,IndexTMM,Energy,Disorder,IndexFB,Width,High,Length)
                        CALL TMM3D(PsiB,PsiA,IndexTMM+1,Energy,Disorder,IndexFB,Width,High,Length)
                        
                        IF ((IndexTMM>(OrthoSteps-1)).AND.((0==mod(IndexTMM,OrthoSteps)).OR.(1==mod(IndexTMM,OrthoSteps))).AND.((OrthoedNo+1)<NoOfOrtho) ) THEN
                            CALL ReNorm(PsiA,PsiB,GammasPerSweep,Width*High)
                            OrthoedNo=OrthoedNo+1
                        END IF
                    END DO
                    IF (1==mod(Length,2)) THEN
                        CALL TMM3D(PsiA,PsiB,Length,Energy,Disorder,IndexFB,Width,High,Length)
                        CALL Swap(PsiA,PsiB,Width*High)
                    END IF
                    CALL WriteTest1(Width,Energy,Disorder)
                    CALL ReNorm(PsiA,PsiB,GammasPerSweep,Width*High)
                    CALL WriteTest2(Width*High,GammasPerSweep)
                    STOP
                END DO !Forward and backward
                !CALL ReSort(PsiA,PsiB,GammasPerSweep,Width*High)
                DO tempIndex=1,Width*High
                    accVariance(Width*High+1-tempIndex)=ABS(GammasPerSweep(tempIndex)/2/Length-OldGammas(tempIndex))/(GammasPerSweep(tempIndex)/2/Length)
                    OldGammas(tempIndex)=GammasPerSweep(tempIndex)/2/Length
                END DO
                CondG_PerSample=ZERO
                DO tempIndex=1,Width*High
                    CondG_PerSample=CondG_PerSample+2.0D0/COSH(OldGammas(tempIndex)*Length)/COSH(OldGammas(tempIndex)*Length)
                END DO
                accCondG_PerSample=ABS(CondG_PerSample-OldCondG_PerSample)/CondG_PerSample
                OldCondG_PerSample=CondG_PerSample
                
                !PRINT *, Index, OldGammas(Width*High)

                ! check accuracy 
                IF( (abs(accVariance(1)) .LE. ConvCriterion) .AND. (abs(accVariance(1)).GE.TINY) .AND. (Index>10)) THEN
                     Conv=1
                     EXIT
                END IF
            END DO !Convegence loop, repeating sweep. 
            !After this we get the convegent 'OldGammas, accVariance, CondG_PerSample, 
            !accCondG_PerSample, MazeWallNo' for each sample
            CALL WriteSampleData(Realization,Index,MazeType,DiagDis,Energy,MazeWallNo,Width, High, Length,&
                                        OldGammas,accVariance,CondG_PerSample,accCondG_PerSample,Ratio)
            CALL PrintSampleData(PrintSwitch,Realization,DiagDis,Energy,MazeWallNo,&
                                        OldGammas(Width*High),accVariance(1))
            
            !Count the effective samples
            IF ((accVariance(1).GE. 0.0).AND.(accVariance(1) .LE. 2.0)) THEN
                EffSamples=EffSamples+1 
                MeanMazeWallRatio=MeanMazeWallRatio+REAL(MazeWallNo)/REAL(Width*High*Length)
                AveCondG=AveCondG+CondG_PerSample
                AveCondGSquare= AveCondGSquare+CondG_PerSample**2
                accCondG_Sum=accCondG_Sum+accCondG_PerSample
                DO tempIndex=1,Width*High
                    AveGammas(tempIndex)=AveGammas(tempIndex)+OldGammas(Width*High+1-tempIndex)
                    AveGammasSquared(tempIndex)=AveGammasSquared(tempIndex)+OldGammas(Width*High+1-tempIndex)**2
                    accVarianceSum(tempIndex)=accVarianceSum(tempIndex)+accVariance(tempIndex)
                END DO
            END IF
        END DO !different sample loop
        IF (0==EffSamples) THEN
            PRINT *, 'No sample convegent for Dis',DiagDis,'En',Energy
            STOP
        END IF
        MeanMazeWallRatio=MeanMazeWallRatio/REAL(EffSamples)
        accCondG_Sum=accCondG_Sum/REAL(EffSamples)
        StdErr=ZERO
        ReletiveErr=ZERO
        DO tempIndex=1,Width*High
            StdErr(tempIndex)=SQRT(ABS(AveGammasSquared(tempIndex)/REAL(EffSamples)-(AveGammas(tempIndex)/REAL(EffSamples))**2)/REAL(MAX(EffSamples-1,1)))
            AveGammas(tempIndex)=AveGammas(tempIndex)/REAL(EffSamples)
            ReletiveErr(tempIndex)=StdErr(tempIndex)/AveGammas(tempIndex)
            accVarianceSum(tempIndex)=accVarianceSum(tempIndex)/REAL(EffSamples)
        END DO
        StdErrCondG=SQRT(ABS(AveCondGSquare/REAL(EffSamples)-(AveCondG/REAL(EffSamples))**2)/REAL(MAX(EffSamples-1,1)))
        AveCondG=AveCondG/REAL(EffSamples)
        ReletiveErrCondG=StdErrCondG/AveCondG
        CALL WriteAveData(LoopDirection,EffSamples,MazeType,DiagDis,Energy,MeanMazeWallRatio,Width, High, Length,&
                          AveGammas,StdErr,ReletiveErr,accVarianceSum,AveCondG,StdErrCondG,ReletiveErrCondG,accCondG_Sum,Ratio)
        CALL CPU_TIME(tock)
        CALL PrintAveData(PrintSwitch,EffSamples,DiagDis,Energy,tock-tick,MeanMazeWallRatio,AveGammas(1),StdErr(1),ReletiveErr(1))
    END DO !flux loop
    
                 
    DEALLOCATE(Maze,Disorder,PsiA,PsiB,GammasPerSweep,OldGammas,accVariance,&
               AveGammas,AveGammasSquared,StdErr,ReletiveErr, accVarianceSum)

    
END PROGRAM CubicMazeTMM
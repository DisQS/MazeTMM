MODULE io_module
    USE SystemParameters
    IMPLICIT NONE
    
    INTERFACE appendBaseName
        MODULE PROCEDURE appendBaseName_r, appendBaseName_i, appendBaseName_c
    END INTERFACE appendBaseName
    
    CONTAINS
    
    SUBROUTINE openUnit(fileName,myUnit,openKind)
        !Purpose: Open a file 'filename' as UNIT=myUnit
        !An error message is returned if the specified file cannot be opened
        IMPLICIT NONE
        INTEGER, INTENT(IN) ::  myUnit
        INTEGER :: fileStatus
        CHARACTER(len=*), INTENT(IN) :: fileName
        CHARACTER(132) ::  stopname
        CHARACTER, INTENT(IN), OPTIONAL :: openKind

	    stopname='*** Cannot open file named '//filename//'***'
        IF(PRESENT(openKind)) THEN
            IF(openKind=='N') THEN
	            OPEN(UNIT=myUnit, FILE=filename, STATUS='NEW', ACTION='WRITE',IOSTAT=fileStatus)
		            IF(fileStatus>0) THEN
			            PRINT *,stopname
			            STOP
		            END IF
            ELSE IF(openKind=='O') THEN
	            OPEN(UNIT=myUnit, FILE=filename, STATUS='OLD', ACTION='READWRITE',IOSTAT=fileStatus)
		            IF(fileStatus>0) THEN
			            PRINT *,stopname
			            STOP
		            END IF
            ELSE IF(openKind=='A') THEN
	            OPEN(UNIT=myUnit, FILE=filename, ACTION='WRITE',POSITION='APPEND',IOSTAT=fileStatus)
		            IF(fileStatus>0) THEN
			            PRINT *,stopname
			            STOP
                    END IF
            ELSE IF(openKind=='R') THEN
	            OPEN(UNIT=myUnit, FILE=filename, STATUS='REPLACE',ACTION='WRITE',IOSTAT=fileStatus)
		            IF(fileStatus>0) THEN
			            PRINT *,stopname
			            STOP
                    END IF
            ELSE
                STOP "Unknown option in openUnit!"
            END IF
        ELSE
	        OPEN(UNIT=myUnit, FILE=filename, STATUS='UNKNOWN', ACTION='READWRITE',IOSTAT=fileStatus)
		        IF(fileStatus>0) THEN
			        PRINT *,stopname
			        STOP
		        END IF
        END IF
    END SUBROUTINE openUnit
    
    SUBROUTINE createFileName(basename,diRectory)
        !Purpose: Begin a file name in the directory diRectory
        IMPLICIT NONE
        CHARACTER(len=*), INTENT(INOUT) :: baseName
        CHARACTER(len=*), INTENT(IN) :: diRectory

        baseName=diRectory
    END SUBROUTINE createFileName
    
    SUBROUTINE appendBaseName_r(basename,partDigs,partValue)
        IMPLICIT NONE
        CHARACTER(len=*), INTENT(INOUT) :: baseName
        INTEGER, INTENT(IN) :: partDigs
        REAL(KIND=rKind), INTENT(IN) :: partValue
        CHARACTER(16) :: iWstring, specString

        !Write the number of decimal places wanted into a string
	        WRITE(iWstring,'(I4)') partDigs
        !Generate the format type for a float string with partDigs decimal places	
	        specString="(1f16."//TRIM(ADJUSTL(iWstring))//")"
        !Write the value partValue to a string using partDigs decimal places
	        WRITE(iWstring,specString) partValue
        !append partname and the value to the given basename
        basename=TRIM(basename)//TRIM(ADJUSTL(iWstring))

    END SUBROUTINE appendBaseName_r
    
    SUBROUTINE appendBaseName_i(basename,partValue)
        !Purpose: Append to a file name the character string partName followed by 
        !the value partValue to partDigs digits for integers
        IMPLICIT NONE
        CHARACTER(len=*), INTENT(INOUT) :: baseName
        INTEGER, INTENT(IN) :: partValue
        CHARACTER(16) :: iWstring, specString

        !Write the number of digits wanted into a string
	        WRITE(iWstring,'(I4)') 16
        !Generate the format type for an integer string with partDigs digits	
	        specString="(I"//TRIM(ADJUSTL(iWstring))//")"
        !Write the value partValue to a string using partDigs digits
	        WRITE(iWstring,specString) partValue
        !append partname and the value to the given basename
        basename=TRIM(basename)//TRIM(ADJUSTL(iWstring))

    END SUBROUTINE appendBaseName_i

    SUBROUTINE appendBaseName_c(basename,partName)
        !Purpose: Append to a file name the character string partName
        IMPLICIT NONE
        CHARACTER(len=*), INTENT(INOUT) :: baseName
        CHARACTER(len=*), INTENT(IN) :: partname

        basename=TRIM(basename)//TRIM(partName)

    END SUBROUTINE appendBaseName_c
    
    LOGICAL FUNCTION CheckName(baseName)
        !Purpose: Returns TRUE if file name baseName exists and FALSE otherwise
        IMPLICIT NONE
        CHARACTER(len=*), INTENT(IN) :: baseName

        INQUIRE(FILE=baseName, EXIST=CheckName)

    END FUNCTION CheckName
    
    
    SUBROUTINE WriteSampleData(EffSamples,Ind,MzType,Dis,En,MzWallNo,Width, High, Length,&
        G,Var,CondG,accCondG,Ratio)
        IMPLICIT NONE
        
        INTEGER(KIND=IKIND),INTENT(IN):: EffSamples,MzType,MzWallNo,Width, High, Length,Ind
        REAL(KIND=RKIND),INTENT(IN)::Dis,En,G(Width*High),Var(Width*High),CondG,accCondG,Ratio
        CHARACTER(len=132) :: Filename
        INTEGER:: tempIndex
        
        CALL createFileName(Filename,'Output/')
        CALL appendBaseName(Filename,'CubeMazeNew_')
        CALL appendBaseName(Filename,Width)
        CALL appendBaseName(Filename,'X')
        CALL appendBaseName(Filename,High)
        CALL appendBaseName(Filename,'X')
        CALL appendBaseName(Filename,Length)
        CALL appendBaseName(Filename,'Mz')
        CALL appendBaseName(Filename,MzType)
        CALL appendBaseName(Filename,'En')
        CALL appendBaseName(Filename,2,En)
        CALL appendBaseName(Filename,'Dis')
        CALL appendBaseName(Filename,2,Dis)
        CALL appendBaseName(Filename,'Ratio')
        CALL appendBaseName(Filename,2,Ratio)
        CALL appendBaseName(Filename,'_Samples.Dat')
        
        IF (1_IKIND==EffSamples) THEN
            CALL openUnit(Filename,110,'R')
        ELSE
            CALL openUnit(Filename,110,'A')
        END IF
        WRITE(110,*) EffSamples,Ind,MzWallNo,CondG,accCondG,&
                     (G(Width*High+1-tempIndex),Var(tempIndex),tempIndex=1,5)
        CLOSE(110)
    
    END SUBROUTINE WriteSampleData
    
    
    SUBROUTINE PrintSampleData(PrintSwitch,EffSamples,Dis,En,MzWallNo,G,Var)
        IMPLICIT NONE
        
        INTEGER(KIND=IKIND),INTENT(IN):: EffSamples,MzWallNo,PrintSwitch
        REAL(KIND=RKIND),INTENT(IN)::Dis,En,G,Var
        
        
        IF(1==PrintSwitch) THEN
            PRINT '(A20,F8.2,A7,F8.2,A15,I6,A13)', ' TMM Cubic with en ',En,' dis ',&
                   Dis, ' finish the ',EffSamples,'-th sample!'
            PRINT '(A12,I6,A10,F10.6,A7,F11.8)',' MzWallNo:',MzWallNo, ' Gamma:',G,' Var:',Var
        END IF
        
        
    END SUBROUTINE PrintSampleData
    
    SUBROUTINE WriteAveData( LoopD,EffSamples,MzType,Dis,En,MeanMzWallRatio,Width, High, Length,&
        AveG,SErr,ReErr,accVariSum,ACondG,SECondG,RECondG,accCondGSum,Ratio)
        IMPLICIT NONE
        
        INTEGER(KIND=IKIND),INTENT(IN):: EffSamples,Width, High, Length, MzType
        CHARACTER ::LoopD
        REAL(KIND=RKIND),INTENT(IN)::Dis,En,AveG(Width*High),SErr(Width*High),ReErr(Width*High),accVariSum(Width*High),&
                                     MeanMzWallRatio,ACondG,SECondG,RECondG,accCondGSum,Ratio
        CHARACTER(len=132) :: Filename
        INTEGER:: tempIndex
        LOGICAL ::Existence
        
        SELECT CASE(LoopD)
        CASE('D')
            CALL createFileName(Filename,'Output/')
            CALL appendBaseName(Filename,'CubeMazeNew_')
            CALL appendBaseName(Filename,Width)
            CALL appendBaseName(Filename,'X')
            CALL appendBaseName(Filename,High)
            CALL appendBaseName(Filename,'X')
            CALL appendBaseName(Filename,Length)
            CALL appendBaseName(Filename,'Mz')
            CALL appendBaseName(Filename,MzType)
            CALL appendBaseName(Filename,'En')
            CALL appendBaseName(Filename,2,En)
            CALL appendBaseName(Filename,'Ratio')
            CALL appendBaseName(Filename,2,Ratio)
            CALL appendBaseName(Filename,'.Dat')
            
            Existence=CheckName(Filename)
            IF (.FALSE.==Existence) THEN
                CALL openUnit(Filename,110,'N')
            ELSE
                CALL openUnit(Filename,110,'A')
            END IF
            
            WRITE(110,*) Dis, EffSamples, MeanMzWallRatio,ACondG,SECondG,RECondG,accCondGSum, &
                         (AveG(tempIndex),SErr(tempIndex),ReErr(tempIndex),accVariSum(tempIndex),tempIndex=1,5)
            CLOSE(110)
        CASE('E')
            CALL createFileName(Filename,'Output/')
            CALL appendBaseName(Filename,'CubeMazeNew_')
            CALL appendBaseName(Filename,Width)
            CALL appendBaseName(Filename,'X')
            CALL appendBaseName(Filename,High)
            CALL appendBaseName(Filename,'X')
            CALL appendBaseName(Filename,Length)
            CALL appendBaseName(Filename,'Mz')
            CALL appendBaseName(Filename,MzType)
            CALL appendBaseName(Filename,'Dis')
            CALL appendBaseName(Filename,2,Dis)
            CALL appendBaseName(Filename,'Ratio')
            CALL appendBaseName(Filename,2,Ratio)
            CALL appendBaseName(Filename,'.Dat')
            
            Existence=CheckName(Filename)
            IF (.FALSE.==Existence) THEN
                CALL openUnit(Filename,110,'N')
            ELSE
                CALL openUnit(Filename,110,'A')
            END IF
            
            WRITE(110,*) En, EffSamples, MeanMzWallRatio,ACondG,SECondG,RECondG,accCondGSum, &
                         (AveG(tempIndex),SErr(tempIndex),ReErr(tempIndex),accVariSum(tempIndex),tempIndex=1,5)
            CLOSE(110)
        CASE DEFAULT
            STOP
        END SELECT
    
    END SUBROUTINE WriteAveData
    
    SUBROUTINE PrintAveData(PSwitch,EffSamples,Dis,En,timeUsed,MeanMzWallRatio,MG,SErr,ReErr)
        IMPLICIT NONE
        
        INTEGER(KIND=IKIND),INTENT(IN):: PSwitch,EffSamples
        REAL(KIND=RKIND),INTENT(IN)::Dis,En,MeanMzWallRatio,MG,SErr,ReErr
        REAL(KIND=RKIND):: timeUsed
        
        IF (1==PSwitch) THEN
            PRINT '(A15,F8.2,A5,F8.2,A25)',' TMM with En ',En,' Dis ',Dis,' finish all samples!'
            PRINT '(A13,I6,A12,F8.2,A16,F8.3)',' Sample No:',EffSamples, 'Time used:',timeUsed/3600.0,&
                  'MeanWallRatio:',MeanMzWallRatio
            PRINT '(A12,F10.6,A10,F10.6,A16,F10.6)',' AveGamma:',MG,'Std Err:',SErr,'Reletive Err:',ReErr
        END IF
    
    END SUBROUTINE PrintAveData
    
    SUBROUTINE WriteTest1(W,En,Dis)
        IMPLICIT NONE
        
        CHARACTER(len=132) :: Filename
        INTEGER(KIND=IKIND),INTENT(IN)::W
        INTEGER(KIND=IKIND)::Index1,Index2,Index3,Index4,iSite,jState,LayerNo
        REAL(KIND=RKIND)::En,Dis(W,W,W)
        REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE ::T,Ts
        INTEGER(KIND=IKIND):: statInt
        
        ALLOCATE(T(2*W*W,2*W*W), STAT=statInt)
		IF(statInt.ne.0) THEN
		    PRINT *, 'Failed to allocate Matrix!'
        END IF 
        ALLOCATE(Ts(2*W*W,2*W*W), STAT=statInt)
		IF(statInt.ne.0) THEN
		    PRINT *, 'Failed to allocate Matrix!'
        END IF 
        
        T=ZERO
        DO Index1=1,2*W*W
            T(Index1,Index1)=ONE
        END DO
        DO LayerNo=1,W
            Ts=ZERO
            DO Index1=1,W
                DO Index2=1,W
                    iSite =(Index1-1)*W+Index2
                    Ts(iSite,iSite)=-En+Dis(Index1,Index2,LayerNo)
                    IF (Index2>1) THEN
                        Ts(iSite,iSite-1)=-ONE
                    END IF
                    IF (Index2<W) THEN
                        Ts(iSite,iSite+1)=-ONE
                    END IF
                    IF(Index1>1) THEN
                        Ts(iSite,iSite-W)=-ONE
                    END IF
                    IF (Index1<W) THEN
                        Ts(iSite,iSite+W)=-ONE
                    END IF
                    Ts(iSite,W*W+iSite)=-ONE
                    Ts(W*W+iSite,iSite)=ONE
                END DO
            END DO
            T=MATMUL(Ts,T)
        END DO
        
        CALL createFileName(Filename,'Output/TMM.DAT')
        
        CALL openUnit(Filename,110,'R')
        DO Index1=1,2*W*W
                WRITE(110,*) (T(Index1,Index2), Index2=1,2*W*W)
        END DO
        CLOSE(110)
        
        DEALLOCATE(T,Ts)
        
    END SUBROUTINE WriteTest1
    
    SUBROUTINE WriteTest2(M,G)
        IMPLICIT NONE
        
        CHARACTER(len=132) :: Filename
        INTEGER(KIND=IKIND),INTENT(IN)::M
        INTEGER(KIND=IKIND)::Index1,Index2
        REAL(KIND=RKIND),INTENT(IN)::G(M)
    
        CALL createFileName(Filename,'Output/Gamma.DAT')
        
        CALL openUnit(Filename,110,'R')
        DO Index1=1,M
                WRITE(110,*) G(Index1)
        END DO
        CLOSE(110)
    
    END SUBROUTINE WriteTest2
    
END MODULE io_module
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
    
    SUBROUTINE LoadMaze(mz,d,Index1,potential)
        IMPLICIT NONE
        
        REAL(KIND=RKIND),DIMENSION(:,:) :: mz
        INTEGER(KIND=IKIND) :: d,Index1,temp1,temp2,tempVal
        REAL(KIND=RKIND):: potential
        CHARACTER(len=132) :: Filename
        
        CALL createFileName(Filename,'Mazes/')
        CALL appendBaseName(Filename,d)
        CALL appendBaseName(Filename,'X')
        CALL appendBaseName(Filename,d)
        CALL appendBaseName(Filename,'maze')
        CALL appendBaseName(Filename,Index1)
        CALL appendBaseName(Filename,'.Dat')
    
        mz=ZERO
        CALL openUnit(Filename,110,'O')
        DO temp1=1,d 
           DO temp2=1,d
               READ(110,*) tempVal
               IF (1==tempVal) THEN
                   mz(temp2,temp1)=potential
               END IF
           END DO
        END DO
        CLOSE(110)
        
    END SUBROUTINE LoadMaze
    
    
    SUBROUTINE ShowMaze(mz,d)
        IMPLICIT NONE
        REAL(KIND=RKIND),DIMENSION(:,:) :: mz
        INTEGER(KIND=IKIND) :: d,index,temp1,temp2
        REAL(KIND=RKIND):: MzWithWall(d+2,d+2)
        
        
        MzWithWall=1.0
        MzWithWall(2, 1) = 0.0
        MzWithWall(d+1, d+2) = 0.0
        DO temp1=1,d
            DO temp2=1,d
                MzWithWall(temp2+1,temp1+1)=mz(temp2,temp1)
            END DO
        END DO
        
        PRINT *, 'Maze used in this TMM:'
        DO temp1=1,(d+2)
            DO temp2=1,(d+2)
                IF (abs(MzWithWall(temp2,temp1))<0.1_RKIND) THEN
                    WRITE(*,'(A2)',ADVANCE='NO') '  '
                ELSE
                    WRITE(*,'(A2)',ADVANCE='NO') '[]'
                END IF
            END DO
            WRITE(*,*)' '
        END DO
    END SUBROUTINE ShowMaze
    
    
    SUBROUTINE WriteAveGammasAndVarToFile(out1,out2,out3)
        IMPLICIT NONE
        CHARACTER(len=132) :: Filename
        INTEGER(KIND=IKIND):: out1
        REAL(KIND=RKIND)::out2,out3
        
        CALL createFileName(Filename,'Output/')
        CALL appendBaseName(Filename,'TMM_')
        CALL appendBaseName(Filename,Width)
        CALL appendBaseName(Filename,'X')
        CALL appendBaseName(Filename,Width)
        CALL appendBaseName(Filename,'Maze')
        CALL appendBaseName(Filename,MazeType)
        CALL appendBaseName(Filename,'_En')
        CALL appendBaseName(Filename,2,Energy)
        CALL appendBaseName(Filename,'_Dis')
        CALL appendBaseName(Filename,2,DiagDis)
        CALL appendBaseName(Filename,'_convProc.Dat')
        
        IF (1_IKIND==out1) THEN
            CALL openUnit(Filename,110,'R')
        ELSE
            CALL openUnit(Filename,110,'A')
        END IF
        WRITE (110,*) out1,out2,out3
        CLOSE(110)
              
    END SUBROUTINE WriteAveGammasAndVarToFile
    
    
    LOGICAL FUNCTION CheckName(baseName)
    !Purpose: Returns TRUE if file name baseName exists and FALSE otherwise
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: baseName

    INQUIRE(FILE=baseName, EXIST=CheckName)

    END FUNCTION CheckName
    
    
    SUBROUTINE WriteOutput(Diag_Dis,En,AveG,Vari,tempIndex,Conv1,Psi_A,M)
        IMPLICIT NONE
        
        CHARACTER(len=132) :: Filename
        LOGICAL ::Existence
        REAL(KIND=RKIND)::Diag_Dis,En,AveG(M),Vari(M),Psi_A(M,M)
        INTEGER(KIND=IKIND)::tempIndex,Conv1,M,Index1,Index2, PsiA_Flag=1
        
        CALL createFileName(Filename,'Output/')
        CALL appendBaseName(Filename,'TMM_')
        CALL appendBaseName(Filename,Width)
        CALL appendBaseName(Filename,'X')
        CALL appendBaseName(Filename,Width)
        CALL appendBaseName(Filename,'Maze')
        CALL appendBaseName(Filename,MazeType)
        CALL appendBaseName(Filename,'.DAT')
        
        Existence=CheckName(Filename)
        IF (Existence) THEN
            CALL openUnit(Filename,110,'N')
        ELSE
            CALL openUnit(Filename,110,'A')
        END IF
        WRITE(110,*) Diag_Dis,En,Conv1,tempIndex,AveG(1),Vari(1),AveG(2),Vari(2),AveG(3),Vari(3),AveG(4),Vari(4),AveG(5),Vari(5)
        CLOSE(110)
        
        !write last collumn of PsiA to file
        IF (1==PsiA_Flag) THEN
            CALL createFileName(Filename,'Output/')
            CALL appendBaseName(Filename,'TMM_')
            CALL appendBaseName(Filename,Width)
            CALL appendBaseName(Filename,'X')
            CALL appendBaseName(Filename,Width)
            CALL appendBaseName(Filename,'Maze')
            CALL appendBaseName(Filename,MazeType)
            CALL appendBaseName(Filename,'_En')
            CALL appendBaseName(Filename,2,En)
            CALL appendBaseName(Filename,'_Dis')
            CALL appendBaseName(Filename,2,Diag_Dis)
            CALL appendBaseName(Filename,'WaveFunc.DAT')
        
            CALL openUnit(Filename,110,'R')
            DO Index1=1,Width
                    WRITE(110,*) (Psi_A((Index1-1)*Width+Index2,M)**2, Index2=1,Width)
            END DO
            CLOSE(110)
        END IF
        
    END SUBROUTINE WriteOutput
    
    
END MODULE io_module
    

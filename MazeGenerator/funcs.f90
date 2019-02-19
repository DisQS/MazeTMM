MODULE funcs
    IMPLICIT NONE
    
    CONTAINS
    
    SUBROUTINE AllocateMat(Mat,Size)
        IMPLICIT NONE
        
        INTEGER, DIMENSION(:,:), ALLOCATABLE ::Mat
        INTEGER::Size, statInt
        
        ALLOCATE(Mat(Size,Size), STAT=statInt)
		IF(statInt.ne.0) THEN
		    PRINT *, 'Failed to allocate Vertor!'
		END IF 
    END SUBROUTINE AllocateMat
    
    SUBROUTINE GenerateMaze(width, mz)

       IMPLICIT NONE
       INTEGER, INTENT(IN)  :: width
       INTEGER, INTENT(INOUT)  :: mz(width, width)
       INTEGER ::MzWithWall(width+2,width+2)
       INTEGER :: x, y

       MzWithWall=1
       MzWithWall(2, 2) = 0
       do y = 2, (width+2), 2
          do x = 2, (width+2), 2
             call carve_maze(width+2, width+2, MzWithWall, x, y)
          end do
       end do
       DO y=1,width
           DO x=1,width
               mz(x,y)=MzWithWall(x+1,y+1)
           END DO
       END DO
       
    END SUBROUTINE GenerateMaze
    
    ! Carve the maze at the specified coordinates.
    subroutine carve_maze(width, height, mz, x, y)

       implicit none
       integer, intent(in)     :: width, height
       integer, intent(inout)  :: mz(width, height)
       integer, intent(in)     :: x, y
       real     :: rand
       integer  :: dir, cnt
       integer  :: dx, dy, localx, localy
       integer  :: x1, y1, x2, y2

       call random_number(rand)  
       cnt = 0
       dir = rand * 4
       localx = x
       localy = y
       do
          dx = 0
          dy = 0
          select case (dir)
             case (0)
                dx = 1
             case (1)
                dy = 1
             case (2)
                dx = -1
             case default
                dy = -1
          end select
          x1 = localx + dx
          y1 = localy + dy
          x2 = x1 + dx
          y2 = y1 + dy
          if ( x2 .gt. 1 .and. x2 .lt. width .and. y2 .gt. 1 .and. y2 .lt. height) THEN
              IF (abs(mz(x1, y1)-1) .lt. 0.01 .and. abs(mz(x2, y2)-1) .lt. 0.01) then
                 mz(x1, y1) = 0
                 mz(x2, y2) = 0
                 localx = x2
                 localy = y2
                 call random_number(rand)
                 dir = rand * 4
                 cnt = 0
              ELSE
                  cnt = cnt + 1
                 if (cnt .gt. 3) then
                    exit
                 end if
                 dir = mod(dir + 1, 4)
              END IF
          else
             cnt = cnt + 1
             if (cnt .gt. 3) then
                exit
             end if
             dir = mod(dir + 1, 4)
          end if
       end do

    end subroutine carve_maze
    
    
    SUBROUTINE createFileName(basename,diRectory)
    !
    !Purpose: Begin a file name in the directory diRectory
    !
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(INOUT) :: baseName
    CHARACTER(len=*), INTENT(IN) :: diRectory

    baseName=diRectory

    END SUBROUTINE createFileName
    
    SUBROUTINE appendBaseName_i(basename,partValue)
    !
    !Purpose: Append to a file name the character string partName followed by 
    !the value partValue to partDigs digits for integers
    !
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
    !
    !Purpose: Append to a file name the character string partName
    !
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(INOUT) :: baseName
    CHARACTER(len=*), INTENT(IN) :: partname
    CHARACTER(16) :: iWstring, specString

    basename=TRIM(basename)//TRIM(partName)

    END SUBROUTINE appendBaseName_c
    
    SUBROUTINE openUnit(fileName,myUnit,openKind)
    !
    !Purpose: Open a file 'filename' as UNIT=myUnit
    !An error message is returned if the specified file cannot be opened
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) ::  myUnit
    INTEGER:: fileStatus
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
    
    SUBROUTINE WriteMazeToFile(Width,mz,Index)
    
    IMPLICIT NONE
    
    INTEGER::Width,mz(Width,Width),Index,temp1,temp2,MzWithWall(Width+2,Width+2)
    CHARACTER(len=132) :: Filename
    
    
    CALL createFileName(Filename,'GeneratedMazes/')
    CALL appendBaseName_i(Filename,Width)
    CALL appendBaseName_c(Filename,'X')
    CALL appendBaseName_i(Filename,Width)
    CALL appendBaseName_c(Filename,'maze')
    CALL appendBaseName_i(Filename,Index)
    CALL appendBaseName_c(Filename,'.Dat')
    
    CALL openUnit(Filename,110)
    DO temp1=1,Width
        DO temp2=1,Width
            WRITE(110,'(I6)') mz(temp2,temp1)
        END DO
    END DO
    WRITE(110,*) ' '
    WRITE(110,*) ' '
    
    
    MzWithWall=1
    MzWithWall(2, 1) = 0
    MzWithWall(Width+1, Width+2) = 0
    DO temp1=1,Width
        DO temp2=1,Width
            MzWithWall(temp2+1,temp1+1)=mz(temp2,temp1)
        END DO
    END DO
    
    DO temp1=1,(Width+2)
        DO temp2=1,(Width+2)
            IF (0==MzWithWall(temp2,temp1)) THEN
                WRITE(110,'(A2)',ADVANCE='NO') '  '
            ELSE
                WRITE(110,'(A2)',ADVANCE='NO') '[]'
            END IF
        END DO
        WRITE(110,*)' '
    END DO
    
    CLOSE(110)
    
    END SUBROUTINE WriteMazeToFile
    
    
END MODULE funcs
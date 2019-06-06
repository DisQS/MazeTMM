MODULE util
    USE SystemParameters
    USE RNG
    
    IMPLICIT NONE
    
    CONTAINS
    
    !Allocate real vector
    SUBROUTINE AllocateVec(Vec,Size)
        IMPLICIT NONE
        
        REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE ::Vec
        INTEGER(KIND=IKIND)::Size, statInt
        
        ALLOCATE(Vec(Size), STAT=statInt)
		IF(statInt.ne.0) THEN
		    PRINT *, 'Failed to allocate Vertor!'
		END IF 
    END SUBROUTINE AllocateVec
    
    !Allocate real matrix
    SUBROUTINE AllocateMat(Mat,Size1,Size2)
        IMPLICIT NONE
        
        REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE ::Mat
        INTEGER(KIND=IKIND)::Size1, Size2, statInt
        
        ALLOCATE(Mat(Size1,Size2), STAT=statInt)
		IF(statInt.ne.0) THEN
		    PRINT *, 'Failed to allocate Matrix!'
		END IF 
    END SUBROUTINE AllocateMat
    
    !Allocate real rank 3 tensor
    SUBROUTINE AllocateR3Tensor(Ten,Size1,Size2,Size3)
        IMPLICIT NONE
        
        REAL(KIND=RKIND), DIMENSION(:,:,:), ALLOCATABLE ::Ten
        INTEGER(KIND=IKIND)::Size1, Size2, Size3, statInt
        
        ALLOCATE(Ten(Size1,Size2,Size3), STAT=statInt)
		IF(statInt.ne.0) THEN
		    PRINT *, 'Failed to allocate Tensor!'
		END IF 
    END SUBROUTINE AllocateR3Tensor
    
    SUBROUTINE Generate3DMaze(mz, CellNo, MzType, Width, High, Length, MazeWallNo)

       IMPLICIT NONE
       INTEGER, INTENT(IN)  :: CellNo,MzType,Width, High, Length
       INTEGER, INTENT(INOUT)  ::MazeWallNo
       REAL(KIND=RKIND), INTENT(INOUT)  :: mz(Width, High, Length)
       INTEGER ::MzWithWall(width+2,High+2,Length+2)
       INTEGER :: x, y, z
       CHARACTER(len=132) ::TempStr
       
       WRITE(TempStr,'(I4)') MzType
       SELECT CASE (MzType)
       CASE(1)
           MzWithWall=1
           MzWithWall(2, 2, 2) = 0
           do z = 2, (Length+2), 2
               do y = 2, (High+2), 2
                  do x = 2, (Width+2), 2
                     call carve_maze(Width+2, High+2, Length+2, MzWithWall, x, y, z)
                  end do
               end do
           end do
           MazeWallNo=0 
           DO z=1,Length
               DO y=1,High
                   DO x=1,width
                       mz(x,y,z)=REAL(MzWithWall(x+1,y+1,z+1))
                       IF (ABS(mz(x,y,z)-1.0)<0.01) THEN
                           MazeWallNo=MazeWallNo+1
                       END IF
                   END DO
               END DO
           END DO
      CASE DEFAULT
           PRINT *, 'Type ',TRIM(TempStr),' has not implemented!'
           STOP
      END SELECT
       
    END SUBROUTINE Generate3DMaze
    
    ! Carve the maze at the specified coordinates.
    subroutine carve_maze(width, height, length, mz, x, y, z)

       implicit none
       integer, intent(in)     :: width, height, length
       integer, intent(inout)  :: mz(width, height, length)
       integer, intent(in)     :: x, y, z
       real(KIND=RKIND)     :: rand
       integer  :: dir, cnt
       integer  :: dx, dy, dz, localx, localy, localz
       integer  :: x1, y1, z1, x2, y2, z2

       
       rand=DRANDOM(1) 
       cnt = 0
       dir = rand * 6
       localx = x
       localy = y
       localz = z
       do
          dx = 0
          dy = 0
          dz = 0
          select case (dir)
             case (0)
                dx = 1
             case (1)
                dy = 1
             case (2)
                dz = 1
             case (3)
                dx = -1
             case (4)
                dy = -1
             case default
                dz = -1
          end select
          x1 = localx + dx
          y1 = localy + dy
          z1 = localz + dz
          x2 = x1 + dx
          y2 = y1 + dy
          z2 = z1 + dz
          if ( x2 .gt. 1 .and. x2 .lt. width .and. y2 .gt. 1 .and. y2 .lt. height .and. z2 .gt. 1 .and. z2 .lt. length) THEN
              IF (abs(REAL(mz(x1, y1, z1)-1)) .lt. 0.01 .and. abs(REAL(mz(x2, y2, z2)-1)) .lt. 0.01 ) then
                 mz(x1, y1, z1) = 0
                 mz(x2, y2, z2) = 0
                 localx = x2
                 localy = y2
                 localz = z2
                 rand=DRANDOM(1)
                 dir = rand * 6
                 cnt = 0
              ELSE
                  cnt = cnt + 1
                 if (cnt .gt. 5) then
                    exit
                 end if
                 dir = mod(dir + 1, 6)
              END IF
          else
             cnt = cnt + 1
             if (cnt .gt. 5) then
                exit
             end if
             dir = mod(dir + 1, 6)
          end if
       end do

    end subroutine carve_maze
    
    
    SUBROUTINE GenerateDisorder(Disorder,DisStr,Maze,Width,High,Length,R)
        IMPLICIT NONE
        
        INTEGER(KIND=IKIND),INTENT(IN)::Width,High,Length
        REAL(KIND=RKIND),INTENT(IN)::DisStr,R
        REAL(KIND=RKIND):: rand
        REAL(KIND=RKIND), INTENT(INOUT)::Disorder(Width,High,Length)
        REAL(KIND=RKIND), INTENT(IN)::Maze(Width,High,Length)
        INTEGER(KIND=IKIND)::x,y,z
        
        DO x=1,Width
            DO y=1,High
                DO z=1,Length
                    rand=DRANDOM(1)
                    IF (ABS(Maze(x,y,z)-1.0)<0.1) THEN
                        Disorder(x,y,z)=DisStr*(rand-0.5D0)
                    ELSE
                        Disorder(x,y,z)=R*DisStr*(rand-0.5D0)
                    END IF
                END DO
            END DO
        END DO
        
    END SUBROUTINE GenerateDisorder
    
     !TMM for the 3D system under OBC
    SUBROUTINE TMM3D(PsiA,PsiB,LayerNo,Energy,Disorder,IndexFB,Width,High,Length)
        IMPLICIT NONE
        
        INTEGER(KIND=IKIND)::LayerNo,IndexFB,Width,High,Length,TempLayerNo,index1,index2,index3,index4,jState,iSite
        REAL(KIND=RKIND):: PsiA(Width*High,Width*High),PsiB(Width*High,Width*High),Energy,Disorder(Width,High,Length),temp
        REAL(KIND=RKIND):: OnsitePot,PsiLeft,PsiRight,PsiUp,PsiDown
        
        IF (1==IndexFB) THEN
            TempLayerNo=LayerNo
        ELSE
            TempLayerNo=Length+1-LayerNo
        END IF
        
        DO index1=1,Width
            DO index2=1,High
                OnsitePot=-Energy+Disorder(index1,index2,TempLayerNo)
                
                DO index3=1,Width
                    DO index4=1,High
                        jState=(index3-1)*High+index4
                        iSite =(index1-1)*High+index2
                        IF (1==index2) THEN
                            PsiLeft=ZERO
                        ELSE
                            PsiLeft=PsiA(iSite-1,jState)
                        END IF
                        IF (High==index2) THEN
                            PsiRight=ZERO
                        ELSE
                            PsiRight=PsiA(iSite+1,jState)
                        END IF
                        IF (1==index1) THEN
                            PsiDown=ZERO
                        ELSE
                            PsiDown=PsiA(iSite-High,jState)
                        END IF
                        IF (Width==index1) THEN
                            PsiUp=ZERO
                        ELSE
                            PsiUp=PsiA(iSite+High,jState)
                        END IF
                        
                        temp=OnsitePot*PsiA(iSite,jState)-(Psileft+PsiRight+PsiUp+PsiDown)+PsiB(iSite,jState)*REAL(2*IndexFB-3)
                        PsiB(iSite,jState)=temp
                    END DO
                END DO
            END DO
        END DO
        
    END SUBROUTINE TMM3D
    
    
    SUBROUTINE Swap( PSI_A, PSI_B, M)
        IMPLICIT NONE
  
        INTEGER(KIND=IKIND) M,jState, tempindex
        REAL(KIND=RKIND) PSI_A(M,M), PSI_B(M,M),dummy

        DO jState=1,M
            DO tempindex=1,M     
                dummy= PSI_B(tempindex,jState)
                PSI_B(tempindex,jState)= PSI_A(tempindex,jState)
                PSI_A(tempindex,jState)= dummy
            ENDDO
        ENDDO

    END SUBROUTINE Swap
    
    
    SUBROUTINE ReNorm(PSI_A,PSI_B,G,M)
          IMPLICIT NONE
  
          INTEGER(KIND=IKIND):: M,IVec,JVec,KIndex
          REAL(KIND=RKIND) PSI_A(M,M), PSI_B(M,M)
          REAL(KIND=RKIND) G(M)
          REAL(KIND=RKIND) sum
          REAL(KIND=RKIND) dummy,norm
          EQUIVALENCE (dummy,norm)
  
          DO IVec=1,M
            DO JVec=1,(IVec-1)
                sum= ZERO
                DO KIndex=1,M
                    sum= sum + (PSI_A(KIndex,JVec))*PSI_A(KIndex,IVec)+ (PSI_B(KIndex,JVec))*PSI_B(KIndex,IVec)
                ENDDO !KIndex loop
                DO KIndex=1,M
                    PSI_A(KIndex,IVec)= PSI_A(KIndex,IVec) - sum * PSI_A(KIndex,JVec)
                    PSI_B(KIndex,IVec)= PSI_B(KIndex,IVec) - sum * PSI_B(KIndex,JVec)
                ENDDO !KIndex loop
            ENDDO !JVec loop
     
             ! calculation of norm
             norm= REAL(0.D0,RKIND)
             DO KIndex=1,M                      
                norm= norm + (PSI_A(KIndex,IVec)) * PSI_A(KIndex,IVec)+ (PSI_B(KIndex,IVec)) * PSI_B(KIndex,IVec)
             ENDDO !KIndex loop
             dummy= 1.D0/SQRT(norm)
             DO KIndex=1,M
                PSI_A(KIndex,IVec)= dummy * PSI_A(KIndex,IVec)
                PSI_B(KIndex,IVec)= dummy * PSI_B(KIndex,IVec)
             ENDDO !KIndex loop
             
             dummy = LOG(dummy)
             G(IVec) = G(IVec) - dummy
          ENDDO !IVec loop
  
    END SUBROUTINE ReNorm
    
    !Find the position with smallest value
    SUBROUTINE FindSmall(tempArray,N,number,position)
        IMPLICIT NONE
        
        INTEGER(KIND=IKIND):: N,number,position,Index1
        REAL(KIND=RKIND):: tempArray(N),temp
        
        position=1
        temp=tempArray(1)
        DO Index1=2,number
            IF (tempArray(Index1)<temp) THEN
                temp=tempArray(Index1)
                position=Index1
            END IF
        END DO
    
    END SUBROUTINE FindSmall
          
    !Bubble sort
    SUBROUTINE ReSort( PSI_A, PSI_B, array0, N)
          IMPLICIT NONE

          INTEGER(KIND=IKIND) N
          REAL(KIND=RKIND) PSI_A(N,N),PSI_B(N,N)
          REAL(KIND=RKIND) array0(N)
          INTEGER Index1,Index2,PS,tempIndex
          INTEGER::Sorted=1
          REAL(KIND=RKIND) temp
          
          DO tempIndex=1,(N-1)
              IF (array0(tempIndex)<array0(tempIndex+1)) THEN
                  Sorted=0
                  EXIT
              END IF
          END DO
          IF (1==Sorted) THEN
              RETURN
          END IF
          
          DO tempIndex=1,(N-1)
              Index1=N+1-tempIndex
              CALL FindSmall(array0,N,Index1,PS)
              IF (PS/=Index1) THEN
                  temp=array0(PS)
                  array0(PS)=array0(Index1)
                  array0(Index1)=temp
                  DO Index2=1,N
                      temp=PSI_A(Index2,PS);
                      PSI_A(Index2,PS)=PSI_A(Index2,Index1)
                      PSI_A(Index2,Index1)=temp
                      temp=PSI_B(Index2,PS);
                      PSI_B(Index2,PS)=PSI_B(Index2,Index1)
                      PSI_B(Index2,Index1)=temp
                  END DO
              END IF
          END DO

    END SUBROUTINE ReSort
    
END MODULE util
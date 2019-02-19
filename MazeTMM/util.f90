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
    SUBROUTINE AllocateMat(Mat,Size)
        IMPLICIT NONE
        
        REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE ::Mat
        INTEGER(KIND=IKIND)::Size, statInt
        
        ALLOCATE(Mat(Size,Size), STAT=statInt)
		IF(statInt.ne.0) THEN
		    PRINT *, 'Failed to allocate Matrix!'
		END IF 
    END SUBROUTINE AllocateMat
    
    !TMM for the 3D system with a 2D maze under OBC
    SUBROUTINE TMMMaze3D(Psi_A,Psi_B,Mz,Layer_No,En,Diag_Dis)
        IMPLICIT NONE
        
        INTEGER(KIND=IKIND)::Layer_No,index1,index2,index3,index4,jState,iSite
        REAL(KIND=RKIND):: En,Diag_Dis,temp,Psi_A(Width*Width,Width*Width),Psi_B(Width*Width,Width*Width),Mz(Width,Width)
        REAL(KIND=RKIND):: OnsitePot,PsiLeft,PsiRight,PsiUp,PsiDown
        
        
        DO index1=1,Width
            DO index2=1,Width
                OnsitePot=-En+Mz(index1,index2)+Diag_Dis*(DRANDOM(1)-0.5D0)
                
                DO index3=1,Width
                    DO index4=1,Width
                        jState=(index3-1)*Width+index4
                        iSite =(index1-1)*Width+index2
                        IF (1==index2) THEN
                            PsiLeft=ZERO
                        ELSE
                            PsiLeft=Psi_A(iSite-1,jState)
                        END IF
                        IF (Width==index2) THEN
                            PsiRight=ZERO
                        ELSE
                            PsiRight=Psi_A(iSite+1,jState)
                        END IF
                        IF (1==index1) THEN
                            PsiDown=ZERO
                        ELSE
                            PsiDown=Psi_A(iSite-Width,jState)
                        END IF
                        IF (Width==index1) THEN
                            PsiUp=ZERO
                        ELSE
                            PsiUp=Psi_A(iSite+Width,jState)
                        END IF
                        
                        temp=(OnsitePot*Psi_A(iSite,jState)-(Psileft+PsiRight+PsiUp+PsiDown)-Psi_B(iSite,jState))
                        Psi_B(iSite,jState)=temp
                    END DO
                END DO
            END DO
        END DO
        
    END SUBROUTINE TMMMaze3D
    
    
    SUBROUTINE ReNorm(PSI_A,PSI_B,G,G2,M)
          IMPLICIT NONE
  
          INTEGER(KIND=IKIND):: M,IVec,JVec,KIndex
          REAL(KIND=RKIND) PSI_A(M,M), PSI_B(M,M)
          REAL(KIND=RKIND) G(M), G2(M)
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
             G2(IVec)= G2(IVec) + dummy*dummy
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
    SUBROUTINE ReSort( PSI_A, PSI_B, array0, array1, N)
          IMPLICIT NONE

          INTEGER(KIND=IKIND) N
          REAL(KIND=RKIND) PSI_A(N,N),PSI_B(N,N)
          REAL(KIND=RKIND) array0(N), array1(N)
          INTEGER Index1,Index2,PS,tempIndex
          INTEGER::Sorted=1
          REAL(KIND=RKIND) temp
          
          DO tempIndex=1,(N-1)
              IF (array0(tempIndex)<array0(tempIndex+1)) THEN
                  Sorted=0
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
                  temp=array1(PS)
                  array1(PS)=array1(Index1)
                  array1(Index1)=temp
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
          
          
    SUBROUTINE DoGamma(G,G2,Ave_G,Vari,M,No,N_ofOrtho,Sub_Steps)
          IMPLICIT NONE
  
          INTEGER(KIND=IKIND) M, Index1,No,N_ofOrtho,Sub_Steps
          REAL(KIND=RKIND) G(M), G2(M),Ave_G(M),Vari(M)
          
          DO Index1=1,M
              Ave_G(M+1-Index1)=G(Index1)/REAL(No*N_ofOrtho*Sub_Steps)
              Vari(M+1-Index1)= SQRT( ABS((G2(Index1)/REAL(No) -(G(Index1)/REAL(No))**2 ) &
                      / REAL( MAX(No-1,1) ) )) / ABS( G(Index1)/REAL(No) )
          END DO
          
    END SUBROUTINE DoGamma
    
END MODULE util
PROGRAM MAZEGENERATER
    USE funcs
    IMPLICIT NONE
    
    
    INTEGER :: Width=39 !Maze is in Width X Width
    INTEGER, DIMENSION(:,:),ALLOCATABLE:: mz
    INTEGER:: Index !To index different mazes
    INTEGER :: seed(33)

    CALL itime(seed)
    CALL random_seed(put = seed)
    
    !Check width, only accept odd number
    IF (0==mod(Width,2)) THEN
        PRINT *, 'Only accept odd number in Width!'
        PAUSE
        STOP
    END IF
    
    CALL AllocateMat(mz,Width)
    
    DO Index=1,20
        CALL GenerateMaze(Width,mz)
        CALL WriteMazeToFile(Width,mz,Index)
    END DO
    DEALLOCATE(mz)
    
END PROGRAM MAZEGENERATER
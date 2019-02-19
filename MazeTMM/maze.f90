MODULE mazes
    USE SystemParameters
    
    IMPLICIT NONE
    
    CONTAINS
    
    SUBROUTINE GenerateMaze(width, mz)

       IMPLICIT NONE
       INTEGER (KIND=IKIND), INTENT(IN)  :: width
       REAL(KIND=RKIND), INTENT(INOUT)  :: mz(width, width)
       INTEGER (KIND=IKIND) :: x, y
       INTEGER :: seed(33)

       CALL itime(seed)
       CALL random_seed(put = seed)
       mz=ONE
       mz(2, 2) = ZERO
       do y = 2, width, 2
          do x = 2, width, 2
             call carve_maze(width, width, mz, x, y)
          end do
       end do
       mz(2, 1) = ZERO
       mz(width - 1, width) = ZERO

    END SUBROUTINE GenerateMaze
    
    ! Carve the maze at the specified coordinates.
    subroutine carve_maze(width, height, mz, x, y)

       implicit none
       integer(KIND=IKIND), intent(in)     :: width, height
       REAL(KIND=RKIND), intent(inout)  :: mz(width, height)
       integer(KIND=IKIND), intent(in)     :: x, y
       real     :: rand
       integer (KIND=IKIND) :: dir, cnt
       integer (KIND=IKIND) :: dx, dy, localx, localy
       integer (KIND=IKIND) :: x1, y1, x2, y2

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
              IF (abs(mz(x1, y1)-ONE) .lt. 0.01 .and. abs(mz(x2, y2)-ONE) .lt. 0.01) then
                 mz(x1, y1) = ZERO
                 mz(x2, y2) = ZERO
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
    
    subroutine ShowMaze(width, mz)

       implicit none
       integer(KIND=IKIND), intent(in) :: width
       integer(KIND=IKIND) :: height
       REAL(KIND=RKIND), intent(in) :: mz(width, width)
       integer(KIND=IKIND) :: x, y

       height=width
       do y = 1, height
          do x = 1, width
             if (abs(mz(x, y)) .lt. 0.01) then
                write (*,'(A2)',advance='no') '  '
             else
                write (*,'(A2)',advance='no') '[]'
             end if
          end do
          print *, ''
       end do

    end subroutine ShowMaze
    
END MODULE mazes
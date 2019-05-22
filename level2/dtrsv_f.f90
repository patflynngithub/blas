
! Adapted by Patrick Flynn : 5/16/19

! Example of using BLAS Level 2 DTRSV to get solution x to Ax=b,
!    where A is a triangular square matrix (upper in this case)
!  
!   ./dtrsv N          size of Ax=b: A is NxN, b is Nx1!   
!     
!     Outputs timing to file
! 
! =============================================================================================================
! 

! Example using the BLAS (Level 2) DTRSV procedure
program dtrsv_f

    implicit none
     
    ! local variable declarations
    integer               :: i,j          ! for looping
    integer               :: num_args     ! # of command line arguments
    character(len=32)     :: arg          ! command line argument
    integer               :: n            ! matrix/vector dimension(s)

    double precision, allocatable, dimension(:,:) :: A       ! the matrix of interest
    double precision, allocatable, dimension(:)   :: b_x     ! is both vectors of interest
                                                             ! initially is givne b of Ax=b
                                                             ! after DTRSV is solution x of Ax=b
    ! timing variables
    integer :: start_count, end_count, count_rate
    
    ! BLAS DTRSV variables
    character        :: uplo
    character        :: trans
    character        :: diag
    integer          :: lda
    integer          :: incx
    
    ! check command line argument count
    num_args = command_argument_count()
    if (num_args /= 1) then
        write(*,*) "Incorrect format; should be   ./dtrsv_f N     (note: A=NxN, b=Nx1)"
        stop 1
    end if
    
    ! Get matrix/vector dimension(s) from command line    
    call get_command_argument(1, arg)
    if (len_trim(arg) == 0) stop 1
    read(arg , *) n  ! convert string to integer

    allocate(A(n,n), b_x(n))

    ! initialize matrices
    ! not doing in column major order b/c want
    ! to have same matrix example as dgemm C code
    do i = 1,n
        do j = 1,n
            if (j < i) then
                A(i,j) = 0.0
            else            
                A(i,j) = (i-1)*n + j;
            end if
        end do
    end do
    
    do i = 1,n
        b_x(i) = i + 9
    end do

    ! output given matrix and vector if their dimensions are small enough
    if (n .le. 10) then
        write(*,*) "A ="
        call print_matrix(A,n,n)
        write(*,*) "b = "
        call print_vector(b_x,n)
    end if
    
    ! DTRSV setup
    uplo  = "u"
    trans = 'n'
    diag  = 'n'
    lda   = n
    incx  = 1

    call SYSTEM_CLOCK(start_count, count_rate)
    
    ! ***************************************************************************

    call DTRSV(uplo, trans, diag, n, A, lda, b_x, incx)
    
    ! ***************************************************************************
 
    call SYSTEM_CLOCK(end_count, count_rate)

    ! output solution x if vector size is small enough
    if (n <= 10) then
        write(*,*)
        write(*,*) "x = "
        call print_vector(b_x,n)
    end if

    write(*,'(A,I0)') "matrix/vector dimension(s) = ", n

    ! output timing data into a file 
    open(1, file = 'timingFortranDTRSV.txt', position = 'append')
    write (1,'(A,I0,A,F0.2,A)') "matrix dimension(s) = ", n, "  ", (end_count - start_count) * 1.0 / count_rate, " seconds"
           
    deallocate(A, b_x)
    
end program dtrsv_f

! -------------------------------------------------------------------

! Output vector
subroutine print_vector(v,n)

    implicit none
    
    integer, intent(in)          :: n
    double precision, intent(in) :: v(n)
    
    ! local variable
    integer            :: i
    
    do i = 1,n
        write(*,"(F6.2,A)", advance="no") v(i), " "
    end do
    write(*,*)
    
    return
    
end subroutine print_vector

! -------------------------------------------------------------------

! Output matrix
subroutine print_matrix(Mat,m,n)

    implicit none
    
    double precision, dimension(m,n), intent(in) :: Mat
    integer, intent(in)          :: m,n

    ! local variable
    integer            :: i,j
    
    do i = 1,m
        write(*,*) Mat(i,:)
    end do
    
    return
    
end subroutine print_matrix

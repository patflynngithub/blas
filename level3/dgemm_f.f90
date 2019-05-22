
! Adapted by Patrick Flynn : 5/16/19

! Example of using BLAS Level 3 DGEMM to get product (AB) of two matrices
!   
!     ./dgemm_f m k n               (note: A=mxk, B=kxn)
!     
!     Outputs timing to file
! 
! =============================================================================================================
! 

! Example using the BLAS (Level 3) DGEMM procedure
program dgemm_f

    implicit none
     
    ! local variable declarations
    integer               :: i,j          ! for looping
    integer               :: num_args     ! # of command line arguments
    character(len=32)     :: arg          ! command line argument
    integer               :: m,k,n        ! matrix dimensions

    double precision, allocatable, dimension(:,:)  :: A,B,C       ! the matrices of interest

    ! timing variables
    integer :: start_count, end_count, count_rate
    
    ! BLAS DGEMM variables
    character        :: transA, transB
    double precision :: alpha, beta;
    
    ! check command line argument count
    num_args = command_argument_count()
    if (num_args /= 3) then
        write(*,*) "Incorrect format; should be   ./dgemm_f m k n     (note: A=mxk, B=kxn)"
        stop 1
    end if
    
    ! Get matrix dimensions from command line    
    call get_command_argument(1, arg)
    if (len_trim(arg) == 0) stop 1
    read(arg , *) m  ! convert string to integer
    call get_command_argument(2, arg)
    if (len_trim(arg) == 0) stop 1
    read(arg , *) k  ! convert string to integer
    call get_command_argument(3, arg)
    if (len_trim(arg) == 0) stop 1
    read(arg , *) n  ! convert string to integer

    allocate(A(m,k), B(k,n), C(m,n))

    ! initialize matrices
    ! not doing in column major order b/c want
    ! to have same matrix example as dgemm C code
    do i = 1,m
        do j = 1,k
            A(i,j) = (i-1)*k + j
        end do
    end do
    
    do i = 1,k
        do j = 1,n
            B(i,j) = (i-1)*n + j
        end do
    end do
    
    do j = 1,n
        do i = 1,m
            C(i,j) = 0.0
        end do
    end do

    ! output A,B matrices if small enough
    if (m .le. 10 .AND. k .le. 10 .AND. n .le. 10) then
        write(*,*) "A ="
        call print_matrix(A,m,k)
        write(*,*) "B = "
        call print_matrix(B,k,n)
    end if
    
    ! DGEMM setup
    transA = 'n'
    transB = 'n'
    alpha  = 1.0;
    beta   = 0.0;

    call SYSTEM_CLOCK(start_count, count_rate)
    
    ! ***************************************************************************

    call dgemm(transA, transB, m, n, k, alpha, A, m, B, k, beta, C, m)
    
    ! ***************************************************************************

    call SYSTEM_CLOCK(end_count, count_rate)

    ! output C=AB matrix if small enough
    if (m .le. 10 .AND. n .le. 10) then
        write(*,*) "C = "
        call print_matrix(C,m,n)
    end if

    write(*,'(A,I0,A,I0,A,I0,A,I0)') "matrix dimensions (mxkxn) = ", m, "x",k,"x",n, "  size of C = ", m*n
    write(*,*) "alpha = ", alpha, "    beta = ", beta

    ! output timing data into a file 
    open(1, file = 'timingFortranDGEMM.txt', position = 'append')
    write (1,'(A,I0,A,I0,A,I0,A,I0,A,F0.2,A)') "matrix dimensions (mxkxn)= ", m, "x",k,"x",n, "  size of C = ", m*n, &
                "  ", (end_count - start_count) * 1.0 / count_rate, " seconds"
           
    deallocate(A,B,C)
    
end program dgemm_f

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
    write(*,*)
    
    return
    
end subroutine print_matrix

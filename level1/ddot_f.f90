! Adapted by Patrick Flynn : 5/16/19

! Example of using BLAS Level 1 DDOT to get dot product of two double vectors
!   
!     ./ddot N          N = size of the two vectors
!     
!     Outputs timing to file
! 
! =============================================================================================================
! 

! Example using the BLAS (Level 1) DDOT procedure
program ddot_f

    implicit none

    ! BLAS function
    external DDOT
    double precision :: DDOT
     
    ! variable declarations
    integer               :: i            ! for looping
    integer               :: num_args     ! # of command line arguments
    character(len=32)     :: arg          ! command line argument
    integer               :: n            ! vector size
    double precision      :: dotprod      ! dot product

    double precision, allocatable, dimension(:)      :: va, vb       ! the vectors of interest

    ! timing variables
    integer :: start_count, end_count, count_rate
    
    ! BLAS DDOT variables
    integer :: inca = 1, incb = 1 
    
    ! check command line argument count
    num_args = command_argument_count()
    if (num_args /= 1) then
        write(*,*) "Incorrect format; should be   ./ddot_f vector_size"
        stop 1
    end if
    
    ! Get vector size from command line    
    call get_command_argument(1, arg)
    if (len_trim(arg) == 0) stop 1
    read(arg , *) n  ! convert string to integer
    write (*,*) "vector size = ", n

    allocate(va(n), vb(n))

    ! initialize vectors
    do i = 1,n
        va(i) = i
    end do
    do i = 1,n
        vb(i) = i+2
    end do
    
    ! output vectors if small enough
    if (n .le. 10) then
        write(*,*) "vector a:"
        call print_vector(va,n)
        write(*,*) "vector b:"
        call print_vector(vb,n)
    end if
    
    call SYSTEM_CLOCK(start_count, count_rate)

    dotprod = DDOT(n, va, inca, vb, incb)

    call SYSTEM_CLOCK(end_count, count_rate)

    write (*,*) "dot product = ", dotprod

    ! output timing data into a file 
    open(1, file = 'timingFortranDDOT.txt', position = 'append')  
    write (1,*) "vector size n = ", n, "  ", (end_count - start_count) * 1.0 / count_rate, " seconds"
           
    deallocate(va,vb)
    
end program ddot_f

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

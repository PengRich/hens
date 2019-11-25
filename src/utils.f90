module random_numbers
    implicit none
    real(kind=4), parameter :: pi=4.0*ATAN(1.0)
    ! real(kind=8) :: rn=0.1d0
    contains
        function rand(rn)
            implicit none
			real(kind=8) :: ax, am, ac
			real(kind=8) :: rand, rn

            if(abs(rn-0.0d0).gt.1.d-6) then
                ac = dble(16807)
                rn = ac * rn
                rn = rn - dble(idint(rn))
            else
                ax = dble(8388607)
                am = dble(2147483647)
                rn = ax / am
            endif

            rand = rn

            return
        end function rand

        subroutine generate_normal_rand(n, rn, mean, std, array)
            implicit none
            integer(kind=4) :: n
            real(kind=8) :: rn, mean, std, array(n)
            integer(kind=4) :: i
            real(kind=8) :: temp

            do i=1, n
                array(i) = rand(rn)
            enddo
 
            DO i = 1, n-1, 2
              temp = std * SQRT(-2.0*LOG(array(i))) * COS(2*pi*array(i+1)) + mean
              array(i+1) = std * SQRT(-2.0*LOG(array(i))) * SIN(2*pi*array(i+1)) + mean
              array(i) = temp
            END DO
 
        end subroutine generate_normal_rand

end module random_numbers

module utility
    implicit none
    character(len=:), public, allocatable :: log_filename 
    contains
        subroutine get_log_filename(perfix)
            character(len=8)  :: date
            character(len=10) :: time
            character(len=4) :: folder 
            character(len=4) :: suffix
            character(len=12) :: f1
            character(len=16) :: f2
            character(len=*), optional :: perfix 
            ! using keyword arguments
            call date_and_time(DATE=date, TIME=time)
            folder = "log/"
            suffix = ".log"
            if(present(perfix))then
                allocate(character(len=len(perfix)+21) :: log_filename)
                f1 = trim(date) // trim(time)
                f2 = trim(f1) // suffix
                log_filename = trim(folder)//trim(perfix)//"_"//trim(f2)
            else
                allocate(character(len=20) :: log_filename)
                f1 = trim(date) // trim(time)
                f2 = trim(f1) // suffix
                log_filename = trim(folder) // trim(f2)
            endif
        end subroutine get_log_filename

        subroutine sort(n, a)
            implicit none
            integer(kind=4) :: n, i, j
            real(kind=8) :: a(n), x
         
            do i = 2, n
                x = a(i)
                j = i - 1
                do while (j >= 1)
                    if (a(j) <= x) exit
                    a(j + 1) = a(j)
                    j = j - 1
                end do
                a(j + 1) = x
            end do
        end subroutine

        function median(n, a)
            implicit none
            integer(kind=4) :: n
            real(kind=8) :: a(n)
            real(kind=8) :: median

            if(mod(n, 2) == 0) then
                median = (a(n/2+1) + a(n/2))/2.d0
            else
                median = a(n/2+1)
            end if
 
        end function median

end module utility

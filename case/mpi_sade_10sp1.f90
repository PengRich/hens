program sade_9sp1
    use mpi_sade
    implicit none
    call mpi_run_sade("10sp1", 4, 5*5*4*10, 100, 100.d0, 100, 100, 5)
end program sade_9sp1

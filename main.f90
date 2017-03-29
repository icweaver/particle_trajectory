program main
  use sho 
  implicit none
  integer :: status = 0
  call init
  print*, "Tmax (in seconds)                        : ", sho_tmax
  print*, "Tmax (in binary periods)                 : ", sho_tmax/sho_period
  do while (sho_t.le.sho_tmax)
    call sho_interrupt
    call sho_rk4
    !call sho_collision_test
    if (sho_test == 1) then ! Ends progam if test mass gets too close to accretor 
      print*, "Binary Periods completed before collision: ", sho_t/(2.0*sho_pi/omega)
      call exit(status)
    end if
    call sho_update
  end do
    call sho_exit
end program main



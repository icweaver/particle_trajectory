module sho
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Variables of Module
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer            :: sho_nstep          ! step number
    integer, parameter :: sho_unit   = 10    ! file identifier
    integer, parameter :: sho_unitt  = 11    ! file identifier
    integer, parameter :: sho_circle = 13    ! file identifier
    integer, parameter :: sho_line   = 14    ! file identifier
    integer, parameter :: sho_scratch= 15    ! file identifier
    integer, parameter :: sho_impact= 16    ! file identifier
    integer, parameter :: sho_acceleration   = 17    ! file identifier
    integer, parameter :: sho_region   = 18    ! file identifier
    integer, parameter :: sho_velimpact   = 19    ! file identifier
    real               :: sho_COM            ! M1r1 + M2r2 = 0, definition check
    real               :: sho_U              ! potential energy
    real               :: sho_Cj,sho_Cj0     ! jacobi constant
    real, dimension(6) :: sho_uold           ! solution vector for satellite (old)
    real, dimension(6) :: sho_unew           ! solution vector for satellite (new)
    real :: sho_x0,sho_y0,sho_vx0,sho_vy0
    real               :: sho_u0        
    real, dimension(6) :: sho_uscale         ! scalings used for dt
    real, dimension(6) :: sho_scaled         ! actual scaled vector uscale/udot
    real               :: d10,d20,davg,r0,v0 
    real, dimension(3) :: sho_r0, sho_v0 
    integer            :: sho_test = 0
    real               :: v_nozz, sho_period
    real               :: sho_rc             ! circularization radius
    real               :: sho_K,sho_K0       ! kinetic energy, initial kinetic energy
    real               :: sho_dt             ! time step
    real               :: sho_t              ! time
    real               :: sho_tmax           ! max time
    real, dimension(3) :: sho_L, sho_Lk,sho_L0,sho_L_sample
    real :: Lk,L, L_sample
    real, dimension(6) :: datas,abc,def
    real :: vv, v_min,b1,alpha,beta,L0
    real :: A, sho_theta_s
    real :: b1_norm
    real :: sho_original
    real :: impact_angle
    real :: v_r, sho_alpha, sho_beta, v_theta, vr_angle, v_frac

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Constants
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real, parameter :: sho_pi          = 4*atan(1.0)    ! pi
    real, parameter :: sho_G           = 6.67259e-8     ! gravitational constant
    real, parameter :: sho_kbol        = 1.3806488e-16  ! boltzman constant
    real, parameter :: sho_temperature = 2525.0         ! wasp-12b surface temperature 
    real, parameter :: sho_proton      = 1.67262178e-24 ! mass of proton
    real, parameter :: sho_threshold   = 1.0E-9         ! maximum allowed acceleration at L1
    real, parameter :: sho_terminate   = 1.092E11    ! terminates code if test mass gets too close
    real, parameter, dimension(3) :: iunit = (/1.0,0.0,0.0/) ! i unit vector

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! System Properties
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real, parameter :: sho_a      = 0.0229*1.5e13     ! binary separation (multiplier*au)
    real, parameter :: sho_R      = sho_a !6.955E10          ! solar radius. used for plotting normalization.
    real, parameter :: sho_M1     = 1.35*1.9891e33    ! M1 >= M2 (convention: M1 is the accretor)  
    real, parameter :: sho_M2     = 1.41*1.8986e30    ! M2 <= M1 
    real, parameter :: sho_mu     = sho_M1/(sho_M1+sho_M2)
    real, parameter :: sho_m      = 0.01*5.97219e27   ! multiplier*mass of Earth in g
    real, parameter, dimension(3) :: sho_M1r0 = (/ (sho_M2/(sho_M1+sho_M2))*sho_a, 0.0, 0.0/) 
    real, parameter, dimension(3) :: sho_M2r0 = (/ (-sho_M1/(sho_M1+sho_M2))*sho_a, 0.0, 0.0/) 
    real, parameter, dimension(3) :: sho_Omega  = (/ 0.0,0.0,(sho_G*(sho_M1+sho_M2))**0.5/(sho_a)**1.5 /)
    real, parameter :: omega      = sqrt(dot_product(sho_Omega, sho_Omega))
    integer :: gate = 0
    character (len=*), parameter :: sho_system = "WASP_12/b"
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! The Good Stuff
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        ! Subroutines
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine init ! Sets initial conditions 
            implicit none
            open(unit=sho_scratch,file="./data/scratch.dat",status="unknown")
            open(unit=sho_unit,file="./data/test.dat",status="unknown") ! main data file opened
            open(unit=sho_acceleration,file="./data/acceleration.dat",status="unknown") ! main data file opened
            open(unit=sho_impact,file="./data/impact.dat",status="unknown") ! impact data file opened
            open(unit=sho_region,file="./data/region.dat",status="unknown") ! impact data file opened
            open(unit=sho_velimpact,file="./data/velimpact.dat", Access = "append", status="old") ! impact data file opened
            write(sho_region,*) "binary_time angle"
            sho_period = 2.0*sho_pi/omega 
            v_nozz      = sqrt((3.0*sho_kbol*sho_temperature)/sho_proton)
            print*, 'v_nozz: ', v_nozz
            sho_uscale = sho_a*(1.0E-3)*(/ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 /) 
            b1         = sho_L1(sho_M1r0(1),sho_M2r0(1)) 
            b1_norm    = b1/sho_a
            print*, "sho_mu: ", sho_mu
            print*, "check: ", b1_norm+sho_mu/(b1_norm-1+sho_mu)**2-(1-sho_mu)/(b1_norm+sho_mu)**2
            A      = sho_mu/abs(b1_norm-1+sho_mu)**3+(1-sho_mu)/abs(b1_norm+sho_mu)**3
            print*, "A: ", A
            sho_theta_s = sho_angle(A)
            print*, "theta_s:", sho_theta_s
            alpha = 1.1
            beta = 1.1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        ! Start your engines 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        sho_theta_s = 0.0
            v_frac = 1
            !sho_uold = (/ b1 ,-0.1, 0.0, v_nozz*cos(sho_theta_s), -v_nozz*sin(sho_theta_s), 0.0 /)
            sho_uold = (/ b1 ,-0.1, 0.0, v_nozz*cos(sho_theta_s), -v_nozz*sin(sho_theta_s), 0.0 /)
            !sho_uold   = (/ b1,  -0.01 , 0.0, v_frac*v_nozz*cos(sho_theta_s)&
            !  &, -v_frac*(v_nozz)*sin(sho_theta_s), 0.0 /) ! Don't let y=0, see sho_rk4 
            print*, 'uold: ', sho_uold
            print*, 'acc: ', sho_func(0.0, sho_uold)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            sho_original = atan2(sho_uold(5),sho_uold(4))
            print*, 'v_x: ', sho_uold(4)
            print*, 'v_y: ', sho_uold(5)
            L0     = b1**2*omega !sho_magnitude(sho_cross(sho_uold(1:3),sho_uold(4:6))) 
            print*, "L0: ", L0
            !alpha  = atan2(sho_uold(2),sho_uold(1)) ! angle between x-axis and r
            !beta   = atan2(sho_uold(5),sho_uold(4)) ! angle between x-axis and v
            print*, "Angular momentum at b1: ", sho_angl(sho_uold)   
            sho_unew   = sho_uold ! so sho_interrupt includes initial values
            sho_r0     = sho_uold(1:3) 
            sho_v0     = sho_uold(4:6) 
            v0         = sqrt(dot_product(sho_v0,sho_v0))   ! |vx0, vy0, vz0| 
            d10        = sqrt(dot_product(sho_r0-sho_M1r0, sho_r0-sho_M1r0)) ! distance between tm and m1
            d20        = sqrt(dot_product(sho_r0-sho_M2r0, sho_r0-sho_M2r0)) ! distance between tm and m2
            davg       = (d10+d20)/2.0
            print*, 'v0                : ', sho_v0
            print*, 'M1 (x,y,z)        : ', sho_M1r0
            print*, 'M2 (x,y,z)        : ', sho_M2r0
            print*, 'b1(x,y,z)         : ', b1 
            print*, 'test mass (x,y,z) : ', sho_r0
            sho_COM    = (-sho_M1*sho_M2*sho_a)/(sho_M1+sho_M2) + (sho_M2*sho_M1*sho_a)/(sho_M1+sho_M2)
            print*, 'M1r1 + M2r2       : ', sho_COM
            print*, 'Omega             : ', omega
            print*, 'binary period     : ', (2.0*sho_pi)/omega 
            sho_tmax   = sho_period*4.0 ! period*multiplier 
            sho_U0     = sho_energy(sho_uold)
            sho_K0     = 0.5*dot_product(sho_uold(4:6),sho_uold(4:6))
            sho_Cj0    = -2.0*(sho_U0)
            open(unit=sho_circle,file="./data/circle.dat",status="unknown")
            write(sho_circle,*) sho_M1r0(1)/sho_R, sho_M1r0(2)/sho_R, 0.7 
            close(sho_circle)
            open(unit=sho_unitt,file="./data/masses.dat",status="unknown")
            write(sho_unitt,*) "x1 y1 z1 x2 y2 z2 b1 -- -- -- -- --" 
            write(sho_unitt,*) sho_M1r0/sho_R,sho_M2r0/sho_R,b1/sho_R,0.0,0.0,0.0,0.0,0.0 ! com at origin  
            write(sho_scratch,*) "x y z time(binary_period)" ! for artificial energy loss
            close(sho_unitt)
            open(unit=sho_line,file="./data/line.dat",status="unknown")
            sho_nstep  = 1.0            ! iteration count
            sho_t      = 0.0            ! total elapsed time
            v_min       = sqrt((sho_G*sho_M1)/(sho_a*1.0E-5))
            sho_L_sample = sqrt(sho_G*sho_M1*(sho_a/3.0))
            L_sample = sho_magnitude(sho_L_sample)
            write(sho_unit,*) "step time x y z vx vy vz U K Cj L Lk L0 distance dt"
            print *, 'Distance from planet to L1: ', abs(sho_M2r0(1) - b1)/(1.204E10), 'planet radii'
          end subroutine init

          subroutine sho_interrupt ! Writes calculated info to data file
            ! writes data, computes normalized conserved quantities
            implicit none
            real :: r,vec(3),v
            sho_U  = sho_energy(sho_unew)
            sho_K  = (0.5*dot_product(sho_unew(4:6),sho_unew(4:6)))
            sho_Cj = (sho_U + sho_K)/(sho_U0+sho_K0) ! -2 times the total normalized energy
            sho_L  = sho_angl(sho_unew)
            L      = sho_L(3)
            vec    = sho_unew(1:3) - sho_M1r0
            r      = sqrt(dot_product(vec,vec))
            Lk     = sqrt(sho_G*sho_M1*r)
            write(sho_unit,*) sho_nstep,sho_t/sho_period,sho_unew(1:3)/sho_R,sho_unew(4:6) &
                              &,sho_U,sho_K,sho_Cj,L,Lk,L0,sho_distance(sho_unew(1:3)), sho_dt
          end subroutine sho_interrupt

         subroutine sho_rk4  ! Takes a step 
            implicit none
            real, dimension(6) :: sho_K1,sho_K2, sho_K3, sho_K4  ! holds slopes used in rk4
            real, dimension(3) :: sho_r
            real :: kappa,phi,betaPrime,vdubs,r,v,ratio

            r        = sho_magnitude(sho_unew(1:3) - sho_M1r0)
            v        = sho_magnitude(sho_unew(4:6))

            ! calculates ks and takes a step in runge-kutta 4th order
            !sho_scaled = abs(sho_uscale / sho_func(sho_t, sho_uold))
            !sho_dt = 300 !minval(sho_scaled(1:2))
            sho_dt = 1.0E-3*sho_a / v

            sho_K1 = sho_func(sho_t, sho_uold)
            sho_K2 = sho_func(sho_t + sho_dt/2, sho_uold + sho_K1*sho_dt/2)
            sho_K3 = sho_func(sho_t + sho_dt/2, sho_uold + sho_K2*sho_dt/2)
            sho_K4 = sho_func(sho_t + sho_dt, sho_uold + sho_K3*sho_dt)
            sho_unew = sho_uold + (sho_dt/6)*(sho_K1 + 2*(sho_K2+sho_K3) + sho_K4)
            alpha  = atan2(sho_unew(2),sho_unew(1)) ! angle between x-axis and r
            beta   = atan2(sho_unew(5),sho_unew(4)) ! angle between x-axis and v


            if ((sho_unew(1) .lt. sho_M1r0(1)) .and. (sho_uold(2)/sho_unew(2)) .le. 0.0) then ! lap requirement 
              write(sho_scratch,*) sho_unew(1:3)/sho_R, sho_t/sho_period ! Writes adjustment coordinate
              close(sho_scratch)
              gate = 1
            end if

            if (gate .eq. 1) then
            !if ((sho_uold(1)/sho_unew(1)) .le. 0.0) then
              write(sho_line,*) sho_unew(1:3) 
              phi = beta-alpha ! angle between extended r and v

                if (phi .lt. 0.0) then ! makes phi continous 
                  phi = phi + 2.0*sho_pi
                end if

              kappa = (1.0+sin(phi))/2.0  
              !kappa = (1+L0/(r*v))/2.0
              vdubs = kappa*v
              betaPrime = alpha + asin(sin(phi)/kappa)
              !print*, L0/(r*v)
              !print*, L0/(r*kappa*v)
              !betaPrime = alpha + asin(L0/(r*kappa*v))
              sho_unew(4) = vdubs*cos(betaPrime)
              sho_unew(5) = vdubs*sin(betaPrime)
              sho_unew(6) = 0.0
              alpha  = atan2(sho_unew(2),sho_unew(1)) ! angle between x-axis and r
              beta   = atan2(sho_unew(5),sho_unew(4)) ! angle between x-axis and v
              datas = sho_func(sho_t,sho_unew)
              write(sho_acceleration,*) sho_t,sho_magnitude(datas(4:6))
            !end if
            end if

            if (sho_magnitude(sho_unew(1:3)-sho_M1r0)/sho_a .le. 0.7) then ! within c_R
              write(sho_impact,*) "binary_time x y z vx vy vz vr v_total v_kepler K_impact kT_L1 mass_ratio system angle"
              v_r =  sho_magnitude(sho_unew(4:6))*cos(alpha+beta) ! not sure if plus or minus
              v_theta =  sho_magnitude(sho_unew(4:6))*sin(alpha+beta)
              write(sho_impact,*) sho_t/sho_period, sho_unew(1:3)/sho_R, sho_unew(4:6)/sho_magnitude(sho_unew(4:6)), v_r, &
              & sho_magnitude(sho_unew(4:6)), sqrt(sho_G*sho_M1/sho_magnitude(sho_unew(1:3))), &
              & (0.5*sho_proton*dot_product(sho_unew(4:6),sho_unew(4:6))), sho_kbol*sho_temperature, sho_M1/sho_M2,&
              & sho_system, sho_pi/2. - sho_vrangle(-1.0*(sho_unew(1:3)-sho_M1r0),sho_unew(4:6)) !
              !makes tails of vectors meet
              close(sho_impact)

              write(sho_velimpact,*) v_frac ,v_frac*v_nozz/sqrt(sho_G*sho_M1/sho_a), &
              & (sho_pi/2. - sho_vrangle(-1.0*(sho_unew(1:3)-sho_M1r0),sho_unew(4:6)))&
              & * 360.0 / (2.0 *sho_pi)   
              close(sho_velimpact)
              !sho_test = 1 ! hops out of code
            end if  
            if (sho_t/sho_period .ge. 0.315490-0.315490*0.4 &
             & .and. sho_t/sho_period .le. 0.315490+0.315490*0.4) then ! within region 
                write(sho_region,*) sho_t/sho_period, sho_vrangle(-1.0*(sho_unew(1:3)),sho_unew(4:6)) !angles here
                ! negative sign flips r, sho_unew(1:3), so that the r and v
                ! vectors are tail to tail for the dot product definition
            end if


          end subroutine sho_rk4

          subroutine sho_update ! New value of last iteration becomes old one for the next. Time updated.
            ! updates vars
            implicit none
            sho_uold   = sho_unew
            sho_nstep  = sho_nstep + 1
            sho_t      = sho_t + sho_dt
          end subroutine sho_update

          subroutine sho_collision_test
            sho_test =  sho_collision(sho_unew(1:3))
          end subroutine sho_collision_test 
        
          subroutine sho_exit ! Ends program, clears memory.
            implicit none
            close(sho_line)
            close(sho_unit)
            close(sho_acceleration)
            close(sho_impact)
            close(sho_region)
          end subroutine sho_exit 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        ! Functions
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          function sho_angle(A)
            implicit none
            real :: sho_angle,A,dummy

            dummy     = -4/(3*A) + sqrt(1-8/(9*A))
            !dummy     = -sqrt(8/(9*A))*sqrt(1-2/A+3*sqrt(1-8/(9*A)))
            sho_angle = 0.5*asin(dummy) 
          end function

          function sho_energy(uin) ! Calculates energy each time step
            implicit none
            real, dimension(6), intent(in) :: uin
            real, dimension(3) :: pos
            real :: sho_energy
            pos         = uin(1:3)
            sho_energy  = sho_pot(pos) ! jacobi constant
          end function sho_energy
          
          function sho_pot(pos) ! Used by sho_energy in calculating potential energy
            implicit none
            real, dimension(3), intent(in) :: pos 
            real :: dist0,dist1,dist2,dummy
            real :: sho_pot
            dist0   = sqrt(dot_product(pos,pos)) ! distance of test mass from com (origin)
            dist1   = sqrt(dot_product(pos-sho_M1r0,pos-sho_M1r0)) ! distance from m1
            dist2   = sqrt(dot_product(pos-sho_M2r0,pos-sho_M2r0)) ! distance from m2
            sho_pot = -sho_G*sho_M1/dist1 - sho_G*sho_M2/dist2 - 0.5*(omega**2)*(dist0**2)
          end function 

          function sho_angl(uin) ! Calculates angular momentum
            implicit none
            real, dimension(3) :: sho_angl
            real, dimension(6), intent(in)  :: uin
            real, dimension(3) :: pos, vel

            pos = uin(1:3) - sho_M1r0
            vel = uin(4:6)

            sho_angl = sho_cross(pos,vel)
          end function

          function sho_collision(pos) ! Checks if particle gets too close to accretor 
            implicit none
            integer :: sho_collision
            integer :: test
            real    :: dist0,dist1,dist2
            real, dimension(3), intent(in) :: pos ! sho_unew(1:3)

            dist1   = sqrt(dot_product(pos-sho_M1r0,pos-sho_M1r0)) ! distance from m1
            !dist2   = sqrt(dot_product(pos-sho_M2r0,pos-sho_M2r0)) ! distance from m2

            if (dist1 .le. sho_terminate) then
              sho_collision = 1
            else
              sho_collision = 0
            end if
          end function sho_collision
   
          function sho_func(t,uin) ! Finds time derivatives
            ! function that determines the derivatives of input array uin
            ! for the satellite: 
            ! dx/dt =   vx
            ! dy/dt =   vy
            ! dz/dt =   vz
            ! dvx/dt = - (GM1/r^3)*x
            ! dvy/dt = - (GM1/r^3)*y
            ! dvx/dt = - (GM1/r^3)*z
            implicit none
            ! declaring input and output vars
            real, dimension(6), intent(in)  :: uin
            real, dimension(6) :: sho_func
            real, dimension(3)              :: g_centrifugal,g1,g2,g_coriolis,r,v
            real                            :: t, d1, d2
            
            r = uin(1:3)
            v = uin(4:6)
            ! local vars
            ! translating
            d1 = sqrt(dot_product(sho_M1r0-r,sho_M1r0-r))
            d2 = sqrt(dot_product(sho_M2r0-r,sho_M2r0-r))

            ! getting distances to masses
            ! getting accelerations
            g_centrifugal = (omega**2)*r
            g1 = (sho_G*sho_M1/(d1**3))*(sho_M1r0-r)
            g2 = (sho_G*sho_M2/(d2**3))*(sho_M2r0-r)
            g_coriolis = -2.0*sho_cross(sho_omega, v)
            ! computing derivatives
            sho_func(1:3) = v
            sho_func(4:6) = g1 + g2 + g_centrifugal + g_coriolis
          end function sho_func

          function sho_L1(M1x0, M2x0)
            implicit none
            real :: sho_L1
            real, intent(in)   :: M1x0, M2x0 ! primary and secondary stay on x-axis by convention
            integer            :: i
            real :: dummy(6),udum(6),up,down,acc
            udum(1)   = (M1x0 + M2x0)/2.0 ! location of first test cut on x-axis
            udum(2:6) = 0.0  
            dummy     = sho_func(0.0,udum) ! return acceleration at test cut
            acc       = dummy(4)
            up        = M1x0
            down      = M2x0
            do while (abs(acc) .gt. sho_threshold)
              if (acc < 0.0) then
                down    = udum(1)
                udum(1) = (up+udum(1))/2.0  ! moves cut over half the distance over
              else
                up      = udum(1)
                udum(1) = (down+udum(1))/2.0
              end if 
              dummy = sho_func(0.0,udum) ! updated accelerations
              acc   = dummy(4)  ! acceleration in x direction
            enddo
            sho_L1 =  udum(1) 
            print*, "acceleration at L1: ", acc
          end function sho_L1

          function sho_cross(u,v)
            implicit none
            real, dimension(3) :: sho_cross
            real, dimension(3), intent(in) :: u,v
            sho_cross(1) = u(2)*v(3) - u(3)*v(2)
            sho_cross(2) = u(3)*v(1) - u(1)*v(3)
            sho_cross(3) = u(1)*v(2) - u(2)*v(1)
          end function

          function sho_magnitude(vec)
            implicit none
            real :: sho_magnitude
            real, dimension(3) :: vec
            sho_magnitude = sqrt(dot_product(vec,vec))
          end function sho_magnitude

          function sho_distance(pos) ! Calculates distance from M1 in Rsol
            implicit none
            real :: sho_distance
            real, dimension(3) :: pos
            
            sho_distance = sho_magnitude(pos-sho_M1r0)/sho_R
         end function sho_distance

         function sho_vrangle(u,v) ! computes angle between two vectors
            implicit none
            real :: sho_vrangle
            real :: sho_dot
            real, dimension(3) :: u,v

            sho_vrangle = acos((u(1)*v(1)+u(2)*v(2)+u(3)*v(3))/(sho_magnitude(u)*sho_magnitude(v)))
        end function sho_vrangle

end module sho








module solvers
  use precision
  use functions
  implicit none
  private

  public :: rk2 
  public :: rk4
  
  public :: dopri54
  public :: dopri87

  public :: solver

  interface 
    subroutine  solver(ff, t, dt, u0, u, err)
       use precision
       interface 
         function ff(t,u) result(up)   
           use precision    
           real(dp), intent(in) :: t    
           real(dp), intent(in) :: u(:)    
           real(dp), allocatable :: up(:)    
         end function
       end interface
       real(dp), intent(in) :: t
       real(dp), intent(in) :: dt
       real(dp), intent(in) :: u0(:)
       real(dp), intent(inout) :: u(:)
       real(dp), intent(inout) :: err
    end subroutine solver
  end interface 

  
  contains

  subroutine rk2(f, t, dt, u0, u)
    procedure(func) :: f    
    real(dp), intent(in) :: t
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: u0(:)
    real(dp), intent(inout) :: u(:)

    real(dp), allocatable :: k1(:), k2(:)
 
    k1 = f(t,u0)
    k2 = f(t+dt, u0 + dt*k1)
    u = u0 + (k1+k2)*dt*0.5_dp  

  end subroutine rk2

  subroutine rk4(f, t, dt, u0, u)
    procedure(func) :: f    
    real(dp), intent(in) :: t
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: u0(:)
    real(dp), intent(inout) :: u(:)

    real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:)
 
    k1 = f(t          , u0             )
    k2 = f(t+dt*0.5_dp, u0+dt*k1*0.5_dp)
    k3 = f(t+dt*0.5_dp, u0+dt*k2*0.5_dp)
    k4 = f(t+dt       , u0+dt*k3       )
    u = u0 + (0.16666666666666666666_dp*k1+0.33333333333333333333_dp*k2+ &
          & 0.33333333333333333333_dp*k3+0.16666666666666666666_dp*k4)*dt

  end subroutine rk4
 
  ! --------------------------------------------------
  subroutine dopri54(f, t, h, u0, u, err)
    procedure(func) :: f    
    real(dp), intent(in) :: t, h
    real(dp), dimension(:), intent(in) :: u0 
    real(dp), dimension(:), intent(inout) :: u 
    real(dp), intent(inout) :: err
 
    real(dp), dimension(:,:), allocatable :: k 
    
    allocate(k(size(u0) ,7))
    
    k(:,1) = f(t, u0)*h
    k(:,2) = f(t+h/5.0_dp, u0+k(:,1)/5.0_dp)*h
    k(:,3) = f(t+3.0_dp*h/10.0_dp, u0+3.0_dp/40.0_dp*k(:,1)+9.0_dp/40.0_dp*k(:,2))*h
    k(:,4) = f(t+h*4.0_dp/5.0_dp, u0+44.0_dp/45.0_dp*k(:,1)+160.0_dp/45.0_dp*k(:,3)-168.0_dp/45.0_dp*k(:,2))*h
    k(:,5) = f(t+h*8.d0/9.d0, u0+38744.d0/13122.d0*k(:,1)+128896.d0/13122.d0*k(:,3)-152160.d0/13122.d0*k(:,2) &
            -3816.d0/13122.d0*k(:,4))*h
    k(:,6) = f(t+h, u0+9017.d0/3168.d0*k(:,1)+46732.d0/5247.d0*k(:,3)+5194.d0/18656.d0*k(:,4) &
            -34080.d0/3168.d0*k(:,2)-5103.d0/18656.d0*k(:,5))*h

    u = u0+(35.d0/384.d0*k(:,1)+500.d0/1113.d0*k(:,3)+125.d0/192.d0*k(:,4)+&
              11.d0/84.d0*k(:,6)-2187.d0/6784.d0*k(:,5))

    k(:,7) = f(t+h, u)*h

    ! just estimate the error as  ( z - y ):
    ! err = max(|u(5)-u(4)|)
    ! maxval(u(:)) maxloc(u(:)) indice del valore massimo 
    err = maxval(abs( 71.d0/57600.d0*k(:,1)+22.d0/525.d0*k(:,6) &
        & +71.d0/1920.d0*k(:,4) &
        & -71.d0/16695.d0*k(:,3)-17253.d0/339200.d0*k(:,5)-1.d0/40.d0*k(:,7) ))
 
    deallocate(k)
 
  end subroutine dopri54
  ! --------------------------------------------------


  subroutine dopri87(f, t, dt, u0, u, err)
    procedure(func) :: f    
    real(dp), intent(in) :: t
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: u0(:)
    real(dp), intent(inout) :: u(:)
    real(dp), intent(inout) :: err    

    real(dp), parameter :: c1 = 0                      
    real(dp), parameter :: c2 = 1.0_dp/18.0_dp                   
    real(dp), parameter :: c3 = 1.0_dp/12.0_dp                   
    real(dp), parameter :: c4 = 1.0_dp/8.0_dp                    
    real(dp), parameter :: c5 = 5.0_dp/16.0_dp                   
    real(dp), parameter :: c6 = 3.0_dp/8.0_dp                    
    real(dp), parameter :: c7 = 59.0_dp/400.0_dp                 
    real(dp), parameter :: c8 = 93.0_dp/200.0_dp                 
    real(dp), parameter :: c9 = 5490023248.0_dp/9719169821.0_dp  
    real(dp), parameter :: c10 = 13.0_dp/20                 
    real(dp), parameter :: c11 = 1201146811.0_dp/1299019798.0_dp 
    real(dp), parameter :: c12 = 1.0_dp
    real(dp), parameter :: c13 = 1.0_dp

    real(dp), parameter :: a21 = 1.0_dp /18.0_dp  
    real(dp), parameter :: a31 = 1.0_dp /48.0_dp  
    real(dp), parameter :: a32 = 1.0_dp /16.0_dp  
    real(dp), parameter :: a41 = 1.0_dp /32.0_dp   
    real(dp), parameter :: a43 = 3.0_dp /32.0_dp   
    real(dp), parameter :: a51 = 5.0_dp /16.0_dp   
    real(dp), parameter :: a53 =-75.0_dp /64.0_dp
    real(dp), parameter :: a54 = 75.0_dp /64.0_dp  
    real(dp), parameter :: a61 = 3.0_dp /80.0_dp
    real(dp), parameter :: a64 = 3.0_dp /16.0_dp  
    real(dp), parameter :: a65 = 3.0_dp /20.0_dp
    real(dp), parameter :: a71 = 29443841.0_dp /  614563906.0_dp  
    real(dp), parameter :: a74 = 77736538.0_dp / 692538347.0_dp  
    real(dp), parameter :: a75 =-28693883.0_dp / 1125000000.0_dp    
    real(dp), parameter :: a76 = 23124283.0_dp /  1800000000.0_dp  
   
    real(dp), parameter :: a81 = 16016141.0_dp / 946692911.0_dp  
    real(dp), parameter :: a84 =  61564180.0_dp / 158732637.0_dp
    real(dp), parameter :: a85 =  22789713.0_dp / 633445777.0_dp
    real(dp), parameter :: a86 =  545815736.0_dp / 2771057229.0_dp  
    real(dp), parameter :: a87 = -180193667.0_dp / 1043307555.0_dp 

    real(dp), parameter :: a91 = 39632708.0_dp  / 573591083.0_dp
    real(dp), parameter :: a94 = -433636366.0_dp / 683701615.0_dp
    real(dp), parameter :: a95 = -421739975.0_dp /  2616292301.0_dp
    real(dp), parameter :: a96 =  100302831.0_dp / 723423059.0_dp  
    real(dp), parameter :: a97 =  790204164.0_dp /  839813087.0_dp   
    real(dp), parameter :: a98 =  800635310.0_dp / 3783071287.0_dp

    real(dp), parameter :: a101 =  246121993.0_dp  / 1340847787.0_dp
    real(dp), parameter :: a104 = -37695042795.0_dp / 15268766246.0_dp  
    real(dp), parameter :: a105 = -309121744.0_dp / 1061227803.0_dp                            
    real(dp), parameter :: a106 = -12992083.0_dp / 490766935.0_dp  
    real(dp), parameter :: a107 = 6005943493.0_dp /  2108947869.0_dp
    real(dp), parameter :: a108 = 393006217.0_dp  /  1396673457.0_dp
    real(dp), parameter :: a109 = 123872331.0_dp  / 1001029789.0_dp 

    real(dp), parameter :: a111 = -1028468189.0_dp / 846180014.0_dp             
    real(dp), parameter :: a114 = 8478235783.0_dp / 508512852.0_dp  
    real(dp), parameter :: a115 = 1311729495.0_dp  /1432422823.0_dp  
    real(dp), parameter :: a116 = -10304129995.0_dp / 1701304382.0_dp   
    real(dp), parameter :: a117 = -48777925059.0_dp /  3047939560.0_dp  
    real(dp), parameter :: a118 =  15336726248.0_dp /  1032824649.0_dp
    real(dp), parameter :: a119 =  -45442868181.0_dp / 3398467696.0_dp  
    real(dp), parameter :: a1110 = 3065993473.0_dp / 597172653.0_dp

    real(dp), parameter :: a121 = 185892177.0_dp / 718116043.0_dp  
    real(dp), parameter :: a124 =  -3185094517.0_dp / 667107341.0_dp
    real(dp), parameter :: a125 =  -477755414.0_dp /  1098053517.0_dp  
    real(dp), parameter :: a126 =  -703635378.0_dp /  230739211.0_dp
    real(dp), parameter :: a127 = 5731566787.0_dp  /1027545527.0_dp
    real(dp), parameter :: a128 =  5232866602.0_dp / 850066563.0_dp
    real(dp), parameter :: a129 =  -4093664535.0_dp /  808688257.0_dp  
    real(dp), parameter :: a1210 = 3962137247.0_dp /  1805957418.0_dp
    real(dp), parameter :: a1211 = 65686358.0_dp / 487910083.0_dp 

    real(dp), parameter :: a131 = 403863854.0_dp / 491063109.0_dp  
    real(dp), parameter :: a134 =  -5068492393.0_dp / 434740067.0_dp  
    real(dp), parameter :: a135 =  -411421997.0_dp / 543043805.0_dp
    real(dp), parameter :: a136 =  652783627.0_dp / 914296604.0_dp
    real(dp), parameter :: a137 = 11173962825.0_dp / 925320556.0_dp  
    real(dp), parameter :: a138 = -13158990841.0_dp / 6184727034.0_dp  
    real(dp), parameter :: a139 = 3936647629.0_dp /  1978049680.0_dp  
    real(dp), parameter :: a1310 = -160528059.0_dp / 685178525.0_dp
    real(dp), parameter :: a1311 = 248638103.0_dp / 1413531060.0_dp  

    real(dp), parameter :: b81 = 14005451.0_dp/335480064.0_dp
    real(dp), parameter :: b86 = -59238493.0_dp/1068277825.0_dp
    real(dp), parameter :: b87 =  181606767.0_dp/758867731.0_dp 
    real(dp), parameter :: b88 = 561292985.0_dp/797845732.0_dp
    real(dp), parameter :: b89 = -1041891430.0_dp/1371343529.0_dp
    real(dp), parameter :: b810 = 760417239.0_dp/1151165299.0_dp
    real(dp), parameter :: b811=  118820643.0_dp/751138087.0_dp
    real(dp), parameter :: b812= -528747749.0_dp/2220607170.0_dp
    real(dp), parameter :: b813=  1.0_dp/4.0_dp

    real(dp), parameter :: b71 = 13451932.0_dp/ 455176623.0_dp
    real(dp), parameter :: b76 = -808719846.0_dp/976000145.0_dp
    real(dp), parameter :: b77 = 1757004468.0_dp/5645159321.0_dp
    real(dp), parameter :: b78 = 656045339.0_dp/265891186.0_dp
    real(dp), parameter :: b79 = -3867574721.0_dp/1518517206.0_dp
    real(dp), parameter :: b710 = 465885868.0_dp/322736535.0_dp
    real(dp), parameter :: b711= 53011238.0_dp/667516719.0_dp
    real(dp), parameter :: b712= 2.0_dp/45.0_dp

    real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), k5(:), k6(:)
    real(dp), allocatable :: k7(:), k8(:), k9(:), k10(:), k11(:), k12(:), k13(:)
    real(dp), allocatable :: y(:)
    
    k1 = dt*f(t       , u0 )
    k2 = dt*f(t+dt*c2 , u0+ a21*k1)
    k3 = dt*f(t+dt*c3 , u0+ a31*k1+ a32*k2)
    k4 = dt*f(t+dt*c4 , u0+ a41*k1+ a43*k3)
    k5 = dt*f(t+dt*c5 , u0+ a51*k1+ a53*k3+ a54*k4)
    k6 = dt*f(t+dt*c6 , u0+ a61*k1+ a64*k4+ a65*k5)
    k7 = dt*f(t+dt*c7 , u0+ a71*k1+ a74*k4+ a75*k5+ a76*k6)
    k8 = dt*f(t+dt*c8 , u0+ a81*k1+ a84*k4+ a85*k5+ a86*k6+ a87*k7)
    k9 = dt*f(t+dt*c9 , u0+ a91*k1+ a94*k4+ a95*k5+ a96*k6+ a97*k7+ a98*k8)
    k10= dt*f(t+dt*c10, u0+a101*k1+a104*k4+a105*k5+a106*k6+a107*k7+a108*k8+a109*k9)
    k11= dt*f(t+dt*c11, u0+a111*k1+a114*k4+a115*k5+a116*k6+a117*k7+a118*k8+a119*k9+a1110*k10)
    k12= dt*f(t+dt*c12, u0+a121*k1+a124*k4+a125*k5+a126*k6+a127*k7+a128*k8+a129*k9+a1210*k10+a1211*k11)
    k13= dt*f(t+dt*c13, u0+a131*k1+a134*k4+a135*k5+a136*k6+a137*k7+a138*k8+a139*k9+a1310*k10+a1311*k11)

    u = u0 + b81*k1+b86*k6+b87*k7+b88*k8+b89*k9+b810*k10+b811*k11+b812*k12+b813*k13

    y = u0 + b71*k1+b76*k6+b77*k7+b78*k8+b79*k9+b710*k10+b711*k11+b712*k12

    !err = sqrt(dot_product(u-y, u-y)) 
    err = maxval(abs(u-y)) 

  end subroutine dopri87

end module solvers

   

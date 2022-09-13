module functions
  use precision
  use interpolators, only : poly, poly1
  implicit none
  private

  public :: func
  public :: duffing
  public :: linear_duffing
  public :: variant
  public :: linear_variant
  public :: lorenz
  public :: lorenz_linear
  public :: sys0
  public :: sys1
  public :: sys2
  public :: sol0_duffing
  public :: sol0_variant
  public :: sol0_lorenz

  real(dp), public :: eps, ss, rr, bb
  integer, parameter :: neqs = 3
  integer :: i
  real(dp), public :: u1, u2, u3, x_g, y_g, z_g 
  real(dp), public :: w0
  real(dp), public :: dw
  !real(dp), parameter, public :: cc = 1.0_dp/12.0_dp
  real(dp), parameter, public :: cc = 0.07_dp  

  interface 
    function func(t,u) result(up)   
      use precision    
      real(dp), intent(in) :: t    
      real(dp), intent(in) :: u(:)
      real(dp), allocatable :: up(:)
      
    end function func
    
    function sol0(t) result(up)   
      use precision    
      real(dp), intent(in) :: t
      real(dp), allocatable :: up(:)

    end function sol0

  !  function sol0_lorenz(t) result(up)
   !   use precision
    !  real(dp), intent(in) :: t
     ! real(dp), allocatable :: up(:,:)
   ! end function sol0_lorenz
  end interface
    
  procedure(func), pointer, public :: ff
  procedure(func), pointer, public :: linear_ff
  procedure(sol0), pointer, public :: p_sol0
  
  contains
      
  function duffing(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))

    up(1) = u(2)/w0 
    up(2) = (-1.0_dp - eps*u(1)**2) * u(1)/w0 

  end function duffing

  
  function linear_duffing(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1) = u(2)/w0
    up(2) = (-1.0_dp - 3.0_dp*eps*u1*u1) * u(1)/w0

  end function linear_duffing

  ! dx/dt = y - y^2 - x * (x^2 - y^2 + 2/3 y^3 + c )
  ! dy/dt = x + (y - y^2)*(x^2 - y^2 + 2/3 y^3 + c )
  !
  ! c = 1/12  x=1/2; y=1 sta sull' orbita
  !
  ! c = 0.07 y=1, => x^2 = 1-2/3-0.07 
  !                  x = 0.5131601439446884141  
  function variant(t,u) result(up)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    real(dp) :: p, qy  
    allocate(up(size(u)))
    
    p = u(1)*u(1) - u(2)*u(2) + 2.0_dp/3.0_dp*u(2)**3 + cc 
    qy = u(2) - u(2)*u(2)

    up(1) = (qy - u(1)*p)/w0
    up(2) = (u(1) + qy * p)/w0

  end function variant

  ! p(x,y) = (x^2 - y^2 + 2/3 y^3 + c )
  ! dp/dx = 2x;  dp/dy = -2y + 2y^2 = -2(y-y^2)
  !
  ! fx = qy - x*p
  ! fy = x + qy*p
  ! 
  ! dfx/dx = - p(x,y) - 2x^2;   dfx/dy = 1 - 2y + 2*x*(y - y^2)
  ! dfy/dx = 1 + 2*x*(y - y^2); dfy/dy = (1 - 2y)*p - 2*(y - y^2)^2
  function linear_variant(t,u) result(up)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
   
    real(dp) :: p, qy  
    allocate(up(size(u)))
   
    p = u1*u1 - u2*u2 + 2.0_dp/3.0_dp*u2**3 + cc 
    qy = u2 - u2*u2

    up(1) = (-2.0_dp*u1*u1 - p) * u(1)/w0 + (1.0_dp - 2.0_dp*u2 + 2*u1*qy) * u(2)/w0
    up(2) = (1.0_dp + 2.0_dp*u1*qy) * u(1)/w0 + ((1.0_dp - 2.0_dp*u2)*p - 2.0_dp*qy*qy) * u(2)/w0
    
  end function linear_variant

  ! dx/dt = ss(y - x)
  ! dy/dt = x(rr - z) - y
  ! dz/dt = xy - bb*z

  function lorenz(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = ss*(u(2) - u(1))/w0
    up(2) = (u(1)*( rr - u(3) ) - u(2))/w0
    up(3) = (u(1)*u(2) - bb*u(3))/w0

  end function lorenz

  function lorenz_linear(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = (ss*(u(2) - u(1)))/w0
    up(2) = ((rr-u3)*u(1) -u(2) - u1*u(3))/w0
    up(3) = (u(1)*u2 +u1*u(2)- bb*u(3))/w0

  end function lorenz_linear

  ! dM/dt = A/w0 M
  function sys0(t,u) result(up)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    real(dp), allocatable :: u0(:)
    allocate(up(size(u)))
    allocate(u0(size(u)))

    u0 = poly(t)
    u1 = u0(1) 
    u2 = u0(2)
    if (neqs == 3) then
       u3 = u0(3)
    end if 
    up(:) = linear_ff(t,u)
    !up(:) = linear_variant(t,u)

  end function


  !  dy1/dt = A/w0 y1 + r(x0(t))/w0   
  ! r(t)/w0 = f(x0)/w0 - d/dt x0 
  function sys1(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
   
    real(dp), allocatable :: u0(:) 
    allocate(up(size(u)))
    allocate(u0(size(u)))
    ! polynomial interpolation of u0 
    u0 = poly(t)
    u1 = u0(1)
    u2 = u0(2)
    if (neqs == 3) then
       u3 = u0(3)
    end if

    
    up(:) = linear_ff(t,u) + ff(t,u0) - poly1(t)
    !up(:) = linear_variant(t,u) + variant(t,u0) - poly1(t) 
    !up(:) = linear_duffing(t,u) + duffing(t,u0) - poly1(t) 

  end function sys1

  !  d/dt y2 = A/w0 y2 - 1/w0 d/dt x0
  !
  !  In this way    y = y1 + dw * y2 is the solution to
  !  dy/dt = A/w0 y + r/w0 - dw/w0 d/dt x0  
  function sys2(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    real(dp), allocatable :: u0(:) 
    allocate(up(size(u)))
    allocate(u0(size(u)))
    ! interpolation of u0 
    u0 = poly(t)
    u1 = u0(1)
    u2 = u0(2)
    if (neqs == 3) then
       u3 = u0(3)
    end if


    up(:) = linear_ff(t,u) + poly1(t)/w0
    !up = linear_variant(t,u) + poly1(t)/w0
    !up = linear_duffing(t,u) + poly1(t)/w0
    
  end function sys2

  
  function sol0_variant(t) result(u0)
    real(dp), intent(in) :: t    
    real(dp), allocatable :: u0(:)
    
    !real(dp), parameter :: rrr = 0.5131601439446884141_dp  
    real(dp), parameter :: rrr = 0.5_dp  
    allocate(u0(2))
    
    u0(1) = rrr*cos(t)
    u0(2) = 1.0_dp+rrr*sin(t)

  end function sol0_variant

  function sol0_duffing(t) result(u0)
    real(dp), intent(in) :: t    
    real(dp), allocatable :: u0(:)
    
    allocate(u0(2))

    u0(1) =  cos(t) 
    u0(2) = -sin(t)

  end function sol0_duffing

  function sol0_lorenz(t) result(u0)
    real(dp), intent(in) :: t    
    real(dp), allocatable :: u0(:,:)
    integer :: N = 777

    allocate(u0(3,N+1))

    open (222, file = 'x0_guess.dat', status = 'old')
    open (223, file = 'y0_guess.dat', status = 'old')
    open (224, file = 'z0_guess.dat', status = 'old')
    do i = 1, N+1  
       read(222,*) u0(1, i) !, u0(2, i), u0(3, i)
       read(223,*) u0(2, i)
       read(224,*) u0(3, i)
  
       !print*, u0(2, i) !, u0(2,i), u0(3,i)
    end do
    
    close(222)
    

  end function sol0_lorenz


end module functions

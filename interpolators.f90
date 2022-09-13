module interpolators
  use precision
  implicit none
  private
  public :: set_points
  public :: clean_points 
  public :: poly
  public :: poly1 

  real(dp), allocatable :: x(:)
  real(dp), allocatable :: c(:,:)
  character(1) :: polytype = 'N'

  contains

  ! Performs cubic interpolation
  !        t
  ! |---|--x|---|
  ! t0  t1  t2  t3
  ! 
  subroutine set_points(tt, yy)
    real(dp), intent(in) :: tt(:)
    real(dp), intent(in) :: yy(:,:)

    integer :: nn, nrhs

    if (size(tt) /= size(yy,2)) then
       stop "ERROR: set_points size mismatch"
    end if
    nrhs = size(yy,1)
    nn = size(tt)
    call clean_points()
    allocate(c(nn,nrhs))
    allocate(x(nn))
    c = transpose(yy)
    x = tt
    polytype = 'L'
  end subroutine set_points

  subroutine solve_vandermonde()
    
    real(dp), allocatable :: A(:,:)
    integer :: ii, jj, nrhs, nn, err
    integer, allocatable :: ipv(:)

    nrhs = size(c,2)
    nn = size(x)
    
    allocate(A(nn,nn))   !matrice 5x5
    allocate(ipv(nn))
   
    ! Setup linear system:  A c = y
    ! c1 * x1^0 + c2 * x1^1 + c3 * x1^2 + ... = y1
    ! c1 * x2^0 + c2 * x2^1 + c3 * x2^2 + ... = y2
    ! c1 * x3^0 + c2 * x3^1 + c3 * x3^2 + ... = y3
    ! ...

    do ii = 1, nn
      A(:,ii) = 1.0_dp
      do jj = 2, ii
        A(:,ii) = A(:,ii)*x(:) 
      end do  
    end do
        
    ! Solve A c = y using LAPACK
    call dgesv(nn,nrhs,A,nn,ipv,c,nn,err)

    if (err /= 0) then
       print*, err
       stop "ERROR in dgesv"
    end if

    deallocate(A)
    deallocate(ipv)
    polytype = 'V' 

  end subroutine solve_vandermonde

  
  subroutine clean_points()
    if (allocated(c)) deallocate(c)    
    if (allocated(x)) deallocate(x)
  end subroutine clean_points


  ! Perform actual polynomial interpolation
  ! f = c1 + c2*t + c3*t*t + c4*t*t*t + ...
  function poly(t) result(f)
    real(dp), intent(in) :: t
    real(dp), allocatable :: f(:)

    integer :: ii, jj, nn
    real(dp) :: P

    nn = size(c,1)
    
    allocate(f(size(c,2)))
    f(:) = 0.0_dp
    
    ! Vendermonde or Lagrange interpolation
    select case(polytype)
    case('V')
      do ii = nn, 2, -1
        f(:) = (f(:) + c(ii,:))*t
      end do
      f(:) = f(:) + c(1,:)
    case('L')
      do ii = 1, nn
        P = 1.0_dp
        do jj = 1, nn
          if (jj == ii) cycle
          P = P*(t-x(jj))/(x(ii)-x(jj)) 
        end do  
        f(:) = f(:) + P*c(ii,:)
      end do 
    case default
      stop 'ERROR: polytype not selected'
    end select

  end function poly
  
  ! Perform polynomial interpolation of derivative
  ! f = c2 + 2*c3*t + 3*c4*t*t + ...
  function poly1(t) result(f)
    real(dp), intent(in) :: t
    real(dp), allocatable :: f(:)

    integer :: ii, jj, kk, nn
    real(dp) :: PP, SS
    
    nn = size(c,1)
    allocate(f(size(c,2)))
    f(:) = 0.0_dp
    
    ! Vendermonde or Lagrange interpolation
    select case(polytype)
    case('V')
      do ii = nn, 3, -1
        f(:) = (f(:) + (ii-1)*c(ii,:))*t
      end do
      f(:) = f(:) + c(2,:)
    case('L')
      do ii = 1, nn
        SS = 0.0_dp
        do jj = 1, nn
          if (jj == ii) cycle
          PP = 1.0_dp
          do kk = 1, nn
            if (kk == ii .or. kk == jj) cycle
            PP = PP*(t-x(kk))/(x(ii)-x(kk)) 
          end do
          SS = SS + PP/(x(ii)-x(jj))
        end do
        f(:) = f(:) + SS*c(ii,:)
      end do 
    case default
      stop 'ERROR: polytype not selected'
    end select

  end function poly1



end module interpolators

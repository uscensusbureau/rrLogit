!#####################################################################
! modified on 12/8/22 to use kind=long_int integers for dimensioning
! and indexing
!#######################################################################
module matrix_methods
!##  3/16: included new procedures for QR factorization via Householder
!##     reflections, with pivoting to handle rank deficiency
   use program_constants
   use error_handler
   implicit none
   private ! by default
   public :: cholesky_in_place, invert_lower, &
        premult_lower_by_transpose, sweep_forward, sweep_reverse, &
        householder_qr, householder_qr2, householder_qr_pivot, &
        householder_qr2_pivot, householder_ols, kronecker_symm, &
        kronecker_lower_tri, kronecker
   character(len=*), parameter :: modname = "matrix_methods"
   !####################################################################
contains 
   !####################################################################
   integer(kind=our_int) function cholesky_in_place(a, err, &
        suppress_msg ) result(answer)
      !### Overwrites lower triangle of a symmetric, pos.-def.
      !### matrix a with its cholesky factor.
      !### If suppress_msg = .true., does not store a message in err
      !### if a is not pos def
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(inout) :: a(:,:)
      type(error_type), intent(inout) :: err
      logical, intent(in), optional :: suppress_msg
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "cholesky_saxpy"
      integer(kind=long_int) :: p, j, k
      real(kind=our_dble) :: den
      logical :: suppress_msg_local
      ! begin
      answer = RETURN_FAIL
      if( present( suppress_msg ) ) then
         suppress_msg_local = suppress_msg
      else
         suppress_msg_local = .false.
      end if
      p = size(a,1,kind=long_int)
      if( p /= size(a,2,kind=long_int) ) goto 700
      do j = 1, p
         do k = 1, j-1
            a(j:p,j) = a(j:p,j) - a(j:p,k) * a(j,k)
         end do
         if( a(j,j) <= 0.D0 ) goto 710
         den = sqrt( a(j,j) )
         a(j:p,j) = a(j:p,j) / den
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
700   call err_handle(err, 1, &
           comment = "Non-square matrix encountered; square expected" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
710   if( .not. suppress_msg_local ) then
         call err_handle(err, 1, &
              comment = "Matrix not positive definite" )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      return
   end function cholesky_in_place
   !####################################################################
   integer(kind=our_int) function invert_lower(a, err, suppress_msg ) &
        result(answer)
      !### Overwrites a lower-triangular matrix a with its inverse
      !### by forward substitution. The upper triangle is untouched.
      !### If suppress_msg = .true., does not store a message in err
      !### if a is apparently singular
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(inout) :: a(:,:)
      type(error_type), intent(inout) :: err
      logical, intent(in), optional :: suppress_msg
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "invert_lower"
      integer(kind=long_int) :: p, i, j, k
      real(kind=our_dble) :: sum
      logical :: suppress_msg_local
      ! begin
      answer = RETURN_FAIL
      if( present( suppress_msg ) ) then
         suppress_msg_local = suppress_msg
      else
         suppress_msg_local = .false.
      end if
      p = size(a,1,kind=long_int)
      if( p /= size(a,2,kind=long_int) ) goto 700
      if( a(1,1) == 0.D0 ) goto 710
      a(1,1) = 1.D0 / a(1,1)
      do i = 2, p
         if( a(i,i) == 0.D0 ) goto 710
         a(i,i) = 1.D0 / a(i,i)
         do j = 1, i-1
            sum = 0.D0
            do k = j, i-1
               sum = sum + a(k,j) * a(i,k)
            end do
            a(i,j) = - sum * a(i,i)
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
700   call err_handle(err, 1, &
           comment = "Non-square matrix encountered; square expected" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
710   if( .not. suppress_msg_local ) then
         call err_handle(err, 1, &
           comment = "Matrix apparently singular" )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      return
   end function invert_lower
   !####################################################################
   integer(kind=our_int) function premult_lower_by_transpose(a, b, &
        err) result(answer)
      !### Premultiplies a lower-triangular matrix a by its upper-
      !### triangular transpose to produce a symmetric matrix b.
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(in) :: a(:,:)
      real(kind=our_dble), intent(out) :: b(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = &
           "premult_lower_by_transpose"
      integer(kind=long_int) :: p, i, j, k
      ! begin
      answer = RETURN_FAIL
      p = size(a,1,kind=long_int)
      if( p /= size(a,2,kind=long_int) ) goto 700
      if( ( p /= size(b,1,kind=long_int) ) .or. &
           ( p /= size(b,2,kind=long_int) ) ) goto 710
      do i = 1, p
         do j = 1, i
            b(i,j) = 0.D0
            do k = max(i,j), p   ! skip zero elements
               b(i,j) = b(i,j) + a(k,i) * a(k,j)
            end do
            b(j,i) = b(i,j)      ! copy upper triangle from lower
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps 
700   call err_handle(err, 1, &
           comment = "Non-square matrix encountered; square expected" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
710   call err_handle(err, 1, &
           comment = "Dimensions of matrix arguments not conformable" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
  end function premult_lower_by_transpose
   !####################################################################
   integer(kind=our_int) function sweep_forward(mat, p, err) &
        result(answer)
      !### Sweeps a symmetric matrix mat on position p. Only the 
      !### lower triangle is used; the upper triangle is ignored
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(inout) :: mat(:,:)
      integer(kind=long_int), intent(in) :: p
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "sweep_forward"
      integer(kind=long_int) :: n, i, j
      ! begin
      answer = RETURN_FAIL
      n = size(mat,1,kind=long_int)
      if( size(mat,2,kind=long_int) /= n ) goto 600
      if( (p < 0) .or. (p > n) ) goto 650
      if( abs(mat(p,p)) <= tiny(mat) ) goto 700
      mat(p,p) = -1.D0/mat(p,p)
      do i = 1, p-1
         mat(p,i)=-mat(p,i)*mat(p,p)
      end do
      do j = p+1,n
         mat(j,p)=-mat(j,p)*mat(p,p)
      end do
      if( abs(mat(p,p)) <= tiny(mat) ) goto 700
      do i=1,p-1
         do j=i,p-1
            mat(j,i)=mat(j,i)+mat(p,i)*mat(p,j)/mat(p,p)
         end do
         do j=p+1,n
            mat(j,i)=mat(j,i)+mat(p,i)*mat(j,p)/mat(p,p)
         end do
      end do
      do i=p+1,n
         do j=i,n
            mat(j,i)=mat(j,i)+mat(i,p)*mat(j,p)/mat(p,p)
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
600   call err_handle(err, 1, &
           comment = "Non-square matrix encountered; square expected" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
650   call err_handle(err, 1, &
           comment = "Pivot out of bounds")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
700   call err_handle(err, 1, &
           comment = "Matrix apparently singular" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function sweep_forward
   !####################################################################
   integer(kind=our_int) function sweep_reverse(mat, p, err) &
        result(answer)
      !### Reverse-sweeps a symmetric matrix mat on position p. Only the 
      !### lower triangle is used; the upper triangle is ignored
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(inout) :: mat(:,:)
      integer(kind=long_int), intent(in) :: p
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "sweep_reverse"
      integer(kind=long_int) :: n, i, j
      ! begin
      answer = RETURN_FAIL
      n = size(mat,1,kind=long_int)
      if( size(mat,2,kind=long_int) /= n ) goto 600
      if( (p < 0) .or. (p > n) ) goto 650
      if( abs(mat(p,p)) <= tiny(mat) ) goto 700
      mat(p,p) = -1.D0/mat(p,p)
      do i = 1, p-1
         mat(p,i)=mat(p,i)*mat(p,p)
      end do
      do j = p+1,n
         mat(j,p)=mat(j,p)*mat(p,p)
      end do
      if( abs(mat(p,p)) <= tiny(mat) ) goto 700
      do i=1,p-1
         do j=i,p-1
            mat(j,i)=mat(j,i)+mat(p,i)*mat(p,j)/mat(p,p)
         end do
         do j=p+1,n
            mat(j,i)=mat(j,i)+mat(p,i)*mat(j,p)/mat(p,p)
         end do
      end do
      do i=p+1,n
         do j=i,n
            mat(j,i)=mat(j,i)+mat(i,p)*mat(j,p)/mat(p,p)
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
600   call err_handle(err, 1, &
           comment = "Non-square matrix encountered; square expected" )
       call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
650   call err_handle(err, 1, &
           comment = "Pivot out of bounds")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
700   call err_handle(err, 1, &
           comment = "Matrix apparently singular" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function sweep_reverse
   !####################################################################
   integer(kind=our_int) function house( x, v, beta ) &
        result(answer)
      !### Computes the Householder vector, based on Algorithm 5.1.1
      !### in Golub and Van Loan (1996)
      !### ***Does not check dimensions of arguments***
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(in) :: x(:)
      real(kind=our_dble), intent(out) :: v(:), beta
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "house"
      real(kind=our_dble) :: mu, sigma, tmp
      integer(kind=long_int) :: n, i
      ! begin
      answer = RETURN_FAIL
      n = size(x,kind=long_int)
      if( n > 0 ) then
         sigma = 0.D0
         do i = 2, n
            sigma = sigma + x(i)**2
         end do
         v(1) = 1.D0
         do i = 2, n
            v(i) = x(i)
         end do
         if( sigma == 0.D0 ) then
            beta = 0.D0
         else
            mu = sqrt( x(1)**2 + sigma )
            if( x(1) <= 0.D0 ) then
               v(1) = x(1) - mu
            else
               v(1)= - sigma / ( x(1) + mu )
            end if
            beta = ( 2.D0 * v(1)**2 ) / ( sigma + v(1)**2 )
            tmp = v(1)
            v(:) = v(:) / tmp
            v(1) = 1.D0
         end if
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
    end function house
   !####################################################################
   integer(kind=our_int) function row_house( a, v, beta, w, err) &
        result(answer)
      !### Householder pre-multiplication, based on Algorithm 5.1.2
      !### in Golub and Van Loan (1989)
      !###   v has length nrow(a)
      !###   beta = 2.D0 / t(v)%*%v is assumed
      !###   w is a workspace of length ncol(a)
      !### ***Does not check dimensions of arguments***
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(inout) :: a(:,:)
      real(kind=our_dble), intent(in) :: v(:), beta
      real(kind=our_dble), intent(inout) :: w(:) ! workspace
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "row_house"
      real(kind=our_dble) :: sum
      integer(kind=long_int) :: n, m, i, j
      ! begin
      answer = RETURN_FAIL
      m = size(a,1,kind=long_int)
      n = size(a,2,kind=long_int)
      if( v(1) /= 1.D0 ) goto 600
      if( (m > 0) .and. (n > 0) ) then
         do j = 1, n
            sum = 0.D0
            do i = 1, m
               sum = sum + a(i,j) * v(i)
            end do
            w(j) = - ( beta * sum )
         end do
         do i = 1, m
            do j = 1, n
               a(i,j) = a(i,j) + v(i) * w(j)
            end do
         end do
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
600   call err_handle(err, 1, &
           comment = "First element of v is not 1.D0")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function row_house
   !####################################################################
   integer(kind=our_int) function row_house_vec( a, v, beta, err) &
        result(answer)
      !### Like row_house, except that it assumes a has a single column
      !###   beta = 2.D0 / t(v)%*%v is assumed
      !### ***Does not check dimensions of arguments***
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(inout) :: a(:)
      real(kind=our_dble), intent(in) :: v(:), beta
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "row_house_vec"
      real(kind=our_dble) :: w, sum
      integer(kind=long_int) :: m, i
      ! begin
      answer = RETURN_FAIL
      m = size(a,kind=long_int)
      if( v(1) /= 1.D0 ) goto 600
      if( m > 0 ) then
         sum = 0.D0
         do i = 1, m
            sum = sum + a(i) * v(i)
         end do
         w = - ( beta * sum )
         do i = 1, m
            a(i) = a(i) + v(i) * w
         end do
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
600   call err_handle(err, 1, &
           comment = "First element of v is not 1.D0")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function row_house_vec
   !####################################################################
   integer(kind=our_int) function householder_qr( x, wkn, wkp, err) &
        result(answer)
      !### Householder QR decomposition, based on Algorithm 5.2.1
      !### in Golub and Van Loan (1996).
      !### Does NOT handle rank deficiency. If rank deficiency is
      !### possible, use the pivoting function instead.
      !###   x(n,p) = matrix with n >= p
      !###   wkn is a workspace of length n
      !###   wkp is a workspace of length p
      !### Upon completion, x is overwritten with its QR decomposition
      !###         X = Q %*% R
      !### in factored-form representation. The upper-triangular
      !### matrix R is stored in the upper-right corner of x, and
      !### the remainder of x contains the Householder vectors
      !### normalized to have their first elements equal to 1.D0
      !### (the ones are not stored)
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(inout) :: x(:,:)
      real(kind=our_dble), intent(inout) :: wkn(:), wkp(:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "householder_qr"
      real(kind=our_dble) :: beta
      integer(kind=long_int) :: n, p, j
      ! begin
      answer = RETURN_FAIL
      n = size(x,1,kind=long_int)
      p = size(x,2,kind=long_int)
      if( p > n ) goto 500
      if( size(wkn,kind=long_int) /= n ) goto 600
      if( size(wkp,kind=long_int) /= p ) goto 610
      do j = 1, p
         if( house( x(j:n,j) , wkn(j:n), beta ) == RETURN_FAIL ) &
              goto 800
         if( row_house( x(j:n,j:p) , wkn(j:n), beta, wkp(j:p), &
              err ) == RETURN_FAIL ) goto 800
         if( j < n ) then
            x( (j+1):n, j ) = wkn( (j+1):n )
         end if
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
500   call err_handle(err, 1, &
           comment = "Argument x has more columns than rows")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
600   call err_handle(err, 1, &
           comment = "Workspace wkn has incorrect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
610   call err_handle(err, 1, &
           comment = "Workspace wkp has incorrect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
800   call err_handle(err, 1, &
           comment = "Attempted division by zero" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function householder_qr
   !####################################################################
   integer(kind=our_int) function householder_qr2( x, q, r, &
        wkn, wkp, wknp, err) result(answer)
      !### Like householder_qr, except that it leaves the input
      !### matrix unchanged and returns the R and Q matrices.
      !### Householder QR decomposition: 
      !###         X = Q %*% R where Q orthonormal, R upper-tri
      !###   x(n,p) = X-matrix (input)
      !###   q(n,p) = Q-matrix (output)
      !###   r(p,p) = R-matrix (output)
      !###   wkn is a workspace of length n
      !###   wkp is a workspace of length p
      !###   wknp is a workspace of size (n,p)
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(in) :: x(:,:)
      real(kind=our_dble), intent(out) :: q(:,:), r(:,:)
      real(kind=our_dble), intent(inout) :: wkn(:), wkp(:), wknp(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "householder_qr2"
      integer(kind=long_int) :: n, p, i, j
      real(kind=our_dble) :: beta
      ! begin
      answer = RETURN_FAIL
      n = size(x,1,kind=long_int)
      p = size(x,2,kind=long_int)
      if( ( size(q,1,kind=long_int) /= n ) .or. &
           ( size(q,2,kind=long_int) /= p ) ) goto 500
      if( ( size(r,1,kind=long_int) /= p ) .or. &
           ( size(r,2,kind=long_int) /= p ) ) goto 510
      if( size(wkn,kind=long_int) /= n ) goto 520
      if( size(wkp,kind=long_int) /= p ) goto 530
      if( ( size(wknp,1,kind=long_int) /= n ) .or. &
           ( size(wknp,2,kind=long_int) /= p ) ) goto 540
      wknp(:,:) = x(:,:)
      if( householder_qr( wknp, wkn, wkp, err) == RETURN_FAIL ) goto 800
      ! extract r
      r(:,:) = 0.D0
      do i = 1, p
         do j = i, p
            r(i,j) = wknp(i,j)
         end do
      end do
      ! fill upper part of q with identity
      q(:,:) = 0.D0
      do i = 1, p
         q(i,i) = 1.D0
      end do
      ! premultiply q by sequence of householder matrices, beginning
      ! the last one
      do j = p, 1, -1
         if( j == n ) cycle
         ! extract householder vector
         wkn(j) = 1.D0
         wkn((j+1):n) = wknp((j+1):n,j)
         ! compute beta
         beta = 2.D0 / sum( wkn(j:n)**2 )
         ! premultiply
         if( row_house( q(j:n,:), wkn(j:n), beta, wkp, err) &
             == RETURN_FAIL ) goto 800
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
500   call err_handle(err, 1, &
           comment = "Arguments x and q have different dimensions")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
510   call err_handle(err, 1, &
           comment = "Argument r has incorect dimensions")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
520   call err_handle(err, 1, &
           comment = "Argument wkn has incorect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
530   call err_handle(err, 1, &
           comment = "Argument wkp has incorect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
540   call err_handle(err, 1, &
           comment = "Argument wknp has incorect dimensions")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
800   call err_handle(err, 1, &
           comment = "Attempted division by zero" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function householder_qr2
   !####################################################################
   integer(kind=our_int) function householder_qr_pivot( x, &
        piv, rank, wkn, wkp1, wkp2, err) result(answer)
      !### Householder QR decomposition with column pivoting
      !### (rank-revealing) based on Algorithm 5.4.1
      !### in Golub and Van Loan (4th ed)
      !###   x(n,p) = matrix with n >= p
      !###   piv(p) = integer vector to hold column permutation
      !###   rank = apparent rank of x
      !###   wkn is a workspace of length n
      !###   wkp is a workspace of length p
      !### Upon completion, x is overwritten with its QR decomposition
      !###         X = Q %*% R
      !### in factored-form representation. The upper-triangular
      !### matrix R is stored in the upper-right corner of x, and
      !### the remainder of x contains the Householder vectors
      !### normalized to have their first elements equal to 1.D0
      !### (the ones are not stored)
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(inout) :: x(:,:)
      integer(kind=long_int), intent(out) :: piv(:), rank
      real(kind=our_dble), intent(inout) :: wkn(:), wkp1(:), wkp2(:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "householder_qr_pivot"
      real(kind=our_dble) :: sum, tau, beta, tmp, tolerance
      integer(kind=long_int) :: n, p, j, i, k, itmp
      logical :: stop_me
      ! begin
      answer = RETURN_FAIL
      n = size(x,1,kind=long_int)
      p = size(x,2,kind=long_int)
      if( size(wkn,kind=long_int) /= n ) goto 600
      if( size(wkp1,kind=long_int) /= p ) goto 610
      if( size(wkp2,kind=long_int) /= p ) goto 610
      if( size(piv,kind=long_int) /= p ) goto 620
      do j = 1, p
         piv(j) = j
         sum = 0.D0
         do i = 1, n
            sum = sum + x(i,j)**2
         end do
         wkp1(j) = sum
      end do
      rank = 0
      tau = maxval(wkp1)
      tolerance = 1.D-13 * maxval( abs(x) ) * sqrt( real(p,our_dble) )
      stop_me = .false.
      do
         if( ( tau <= 0.D0 ) .or. ( rank >= p ) .or. &
              ( rank >= n ) .or. stop_me ) exit
         rank = rank + 1
         do k = rank, p
            if( wkp1(k) == tau ) exit
         end do
         itmp = piv(rank)
         piv(rank) = piv(k)
         piv(k) = itmp
         do i = 1, n
            tmp = x(i,rank)
            x(i,rank) = x(i,k)
            x(i,k) = tmp
         end do
         tmp = wkp1(rank)
         wkp1(rank) = wkp1(k)
         wkp1(k) = tmp
         if( house( x(rank:n,rank), wkn(rank:n), beta ) &
              == RETURN_FAIL ) goto 800
         if( row_house( x(rank:n,rank:p), wkn(rank:n), beta, &
              wkp2(rank:p), err ) == RETURN_FAIL ) goto 800
         x( (rank+1):n, rank ) = wkn( (rank+1):n )
         do i = (rank+1), p
            wkp1(i) = wkp1(i) - x(rank,i)**2
         end do
         tau = maxval( wkp1( (rank+1):p ) )
         if( rank < p ) then
            ! check tolerance of remaining part of R
            if( maxval(abs(x((rank+1):p,(rank+1):p))) * &
                 sqrt(real(p-rank,our_dble)) <= tolerance ) then
               stop_me = .true.
            end if
         end if
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
600   call err_handle(err, 1, &
           comment = "Workspace wkn has incorrect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
610   call err_handle(err, 1, &
           comment = "Workspace wkp1 or wkp2 has incorrect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
620   call err_handle(err, 1, &
           comment = "Argument piv has incorrect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
800   call err_handle(err, 1, &
           comment = "Attempted division by zero" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function householder_qr_pivot
  !####################################################################
   integer(kind=our_int) function householder_qr2_pivot( x, &
        piv, rank, q, r, wkn, wkp1, wkp2, wknp, err) result(answer)
      !### Like householder_qr_pivot, except that it leaves the input
      !### matrix unchanged and returns the R and Q matrices
      !###    X %*% P = Q %*% R where Q orthonormal, R upper-tri, and
      !###        P is a permutation
      !###   x(n,p) = matrix (input)
      !###   piv(p) = integer vector for column permutation (output)
      !###   rank = apparent rank of x (output)
      !###   q(n,p) = Q-matrix (output)
      !###   r(p,p) = R-matrix (output)
      !###   wkn is a workspace of length n
      !###   wkp1 is a workspace of length p
      !###   wkp2 is a workspace of length p
      !###   wknp is a workspace of size (n,p)
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(in) :: x(:,:)
      integer(kind=long_int), intent(out) :: piv(:), rank
      real(kind=our_dble), intent(out) :: q(:,:), r(:,:)
      real(kind=our_dble), intent(inout) :: wkn(:), wkp1(:), &
           wkp2(:), wknp(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "householder_qr2_pivot"
      real(kind=our_dble) :: beta
      integer(kind=long_int) :: n, p, j, i
      ! begin
      answer = RETURN_FAIL
      n = size(x,1,kind=long_int)
      p = size(x,2,kind=long_int)
      if( size(piv,kind=long_int) /= p ) goto 300
      if( ( size(q,1,kind=long_int) /= n ) .or. &
           ( size(q,2,kind=long_int) /= p ) ) goto 500
      if( ( size(r,1,kind=long_int) /= p ) .or. &
           ( size(r,2,kind=long_int) /= p ) ) goto 510
      if( size(wkn,kind=long_int) /= n ) goto 520
      if( size(wkp1,kind=long_int) /= p ) goto 530
      if( size(wkp2,kind=long_int) /= p ) goto 530
      if( ( size(wknp,1,kind=long_int) /= n ) .or. &
           ( size(wknp,2,kind=long_int) /= p ) ) goto 540
      wknp(:,:) = x(:,:)
      if( householder_qr_pivot( wknp, piv, rank, wkn, wkp1, &
           wkp2, err) == RETURN_FAIL ) goto 800
      ! extract r
      r(:,:) = 0.D0
      do i = 1, p
         do j = i, p
            r(i,j) = wknp(i,j)
         end do
      end do
      ! fill upper part of q with identity
      q(:,:) = 0.D0
      do i = 1, p
         q(i,i) = 1.D0
      end do
      ! premultiply q by sequence of householder matrices, beginning
      ! the last one
      do j = p, 1, -1
         ! extract householder vector
         wkn(j) = 1.D0
         wkn((j+1):n) = wknp((j+1):n,j)
         ! compute beta
         beta = 2.D0 / sum( wkn(j:n)**2 )
         ! premultiply
         if( row_house( q(j:n,:), wkn(j:n), beta, wkp1, err) &
             == RETURN_FAIL ) goto 800
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
300   call err_handle(err, 1, &
           comment = "Argument piv has incorrect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
500   call err_handle(err, 1, &
           comment = "Arguments x and q have different dimensions")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
510   call err_handle(err, 1, &
           comment = "Argument r has incorect dimensions")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
520   call err_handle(err, 1, &
           comment = "Argument wkn has incorect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
530   call err_handle(err, 1, &
           comment = "Argument wkp1 or wkp2 has incorect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
540   call err_handle(err, 1, &
           comment = "Argument wknp has incorect dimensions")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
800   call err_handle(err, 1, &
           comment = "Attempted division by zero" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function householder_qr2_pivot
  !####################################################################
   integer(kind=our_int) function householder_ols( x, y, &
        rank, omit, coef, wkn1, wkn2, wkp1, wkp2, wknp, wkpp, &
        iwkp, err) result(answer)
      !### Ordinary least squares, possibly of deficient rank,
      !### using Householder QR
      !###   x(n,p) = matrix (input)
      !###   y(n) = response vector (input)
      !###   rank = apparent column rank of x (output)
      !###   omit(p) = logical array (output) (T in position j if
      !###        column j of X was omitted due to rank deficiency)
      !###   coef(p) = least squares coefficients (output)
      !###   wkn1 = real workspace of length n
      !###   wkn2 = real workspace of length n
      !###   wkp1 = real workspace of length p
      !###   wkp2 = real workspace of length p
      !###   wknp = real workspace of size (n,p)
      !###   wkpp = real workspace of size (p,p)
      !###   iwkp = integer workspace of length p
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(in) :: x(:,:), y(:)
      integer(kind=long_int), intent(out) :: rank
      logical, intent(out) :: omit(:)
      real(kind=our_dble), intent(out) :: coef(:)
      real(kind=our_dble), intent(inout) :: wkn1(:), wkn2(:), &
           wkp1(:), wkp2(:), wknp(:,:), wkpp(:,:)
      integer(kind=long_int), intent(inout) :: iwkp(:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      integer(kind=long_int) :: n, p, i, j
      real(kind=our_dble) :: beta
      character(len=*), parameter :: subname = "householder_ols"
      ! begin
      answer = RETURN_FAIL
      n = size(x,1,kind=long_int)
      p = size(x,2,kind=long_int)
      if( size(y,kind=long_int) /= n ) goto 500
      if( size(omit,kind=long_int) /= p ) goto 510
      if( size(coef,kind=long_int) /= p ) goto 520
      if( size(wkn1,kind=long_int) /= n ) goto 530
      if( size(wkn2,kind=long_int) /= n ) goto 530
      if( size(wkp1,kind=long_int) /= p ) goto 540
      if( size(wkp2,kind=long_int) /= p ) goto 540
      if( ( size(wknp,1,kind=long_int) /= n ) &
           .or. ( size(wknp,2,kind=long_int) /= p ) ) goto 550
      if( ( size(wkpp,1,kind=long_int) /= p ) .or. &
           ( size(wknp,2,kind=long_int) /= p ) ) goto 560
      if( size(iwkp,kind=long_int) /= p ) goto 570
      wknp(:,:) = x(:,:)
      if( householder_qr_pivot( wknp, iwkp, rank, wkn1, wkp1, &
           wkp2, err) == RETURN_FAIL ) goto 800
      ! premultiply y by sequence of householder matrices, beginning
      ! with the first one, storing result in wkn2
      wkn2(:) = y(:)
      do j = 1, p
         if( j == n ) exit
         wkn1(j) = 1.D0
         wkn1((j+1):n) = wknp((j+1):n,j)
         ! compute beta
         beta = 2.D0 / sum( wkn1(j:n)**2 )
         if( row_house_vec( wkn2(j:n), wkn1(j:n), beta, err) &
             == RETURN_FAIL ) goto 800
      end do
      ! extract transpose of r, put into wkpp, then invert
      do i = 1, rank
         do j = i, rank
            wkpp(j,i) = wknp(i,j)
         end do
      end do
      if( invert_lower( wkpp(1:rank,1:rank), err ) &
           == RETURN_FAIL ) goto 800
      ! store least squares coefficients in wkp1
      do j = 1, rank
         wkp1(j) = sum( wkpp(j:rank,j) * wkn2(j:rank) )
      end do
      ! restore the coefficients to original order
      do j = 1, rank
         coef( iwkp(j) ) = wkp1(j)
         omit( iwkp(j) ) = .false.
      end do
      do j = (rank+1), p
         coef( iwkp(j) ) = 0.D0
         omit( iwkp(j) ) = .true.
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
500   call err_handle(err, 1, &
           comment = "Arguments x and y have incompatible dimensions")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
510   call err_handle(err, 1, &
           comment = "Argument omit has incorect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
520   call err_handle(err, 1, &
           comment = "Argument beta has incorect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
530   call err_handle(err, 1, &
           comment = "Argument wkn1 or wkn2 has incorect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
540   call err_handle(err, 1, &
           comment = "Argument wkp1 or wkp2 has incorect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
550   call err_handle(err, 1, &
           comment = "Argument wknp has incorect dimensions")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
560   call err_handle(err, 1, &
           comment = "Argument wkpp has incorect dimensions")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
570   call err_handle(err, 1, &
           comment = "Argument iwkp has incorect length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
800   call err_handle(err, 1, &
           comment = "Attempted division by zero" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function householder_ols
  !####################################################################
   integer(kind=our_int) function kronecker_symm(a, b, ab, &
        err) result(answer)
      !### Kronecker product of two symmetric matrices
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(in) :: a(:,:), b(:,:)
      real(kind=our_dble), intent(out) :: ab(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "kronecker_symm"
      integer(kind=long_int) :: na, nb, ia, ja, ib, jb, i, j, ii, jj
      ! begin
      answer = RETURN_FAIL
      na = size( a, 1, kind=long_int )
      if( size( a, 2, kind=long_int ) /= na ) goto 200
      nb = size( b, 1, kind=long_int )
      if( size( b, 2, kind=long_int ) /= nb ) goto 200
      if( size( ab, 1, kind=long_int ) /= &
           size( ab, 2, kind=long_int ) ) goto 200
      if( size( ab, 1, kind=long_int ) /= na * nb ) goto 300
      ! compute lower triangle of ab
      do ia = 1, na
         ! off-diagonal blocks
         do ja = 1, ( ia - 1 )
            i = nb * ( ia - 1 )
            jj = nb * ( ja - 1 )
            do ib = 1, nb
               i = i + 1
               jj = jj + 1
               j = nb * ( ja - 1 )
               ii = nb * ( ia - 1 )
               do jb = 1, ib
                  j = j + 1
                  ii = ii + 1
                  ab( i, j ) = a( ia, ja ) * b( ib, jb )
                  if( ib /= jb ) then
                     ab( ii, jj ) = ab( i, j )
                  end if
               end do
            end do
         end do
         ! diagonal block
         ja = ia
         i = nb * ( ia - 1 )
         do ib = 1, nb
            i = i + 1
            j = nb * ( ja - 1 )
            do jb = 1, ib
               j = j + 1
               ab( i, j ) = a( ia, ja ) * b( ib, jb )
            end do
         end do
      end do
      ! fill in upper triangle
      do i = 1, na * nb
         do j = 1, i - 1
            ab( j, i ) = ab( i, j )
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps 
200   call err_handle(err, 1, &
           comment = "Non-square matrix encountered; square expected" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
300   call err_handle(err, 1, &
           comment = "Dimensions of matrix arguments not conformable" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function kronecker_symm
  !####################################################################
   integer(kind=our_int) function kronecker_lower_tri(a, b, ab, &
        err) result(answer)
      !### Kronecker product of two lower triangular matrices
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(in) :: a(:,:), b(:,:)
      real(kind=our_dble), intent(out) :: ab(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "kronecker_lower_tri"
      integer(kind=long_int) :: na, nb, ia, ja, ib, jb, i, j
      ! begin
      answer = RETURN_FAIL
      na = size( a, 1, kind=long_int )
      if( size( a, 2, kind=long_int ) /= na ) goto 200
      nb = size( b, 1, kind=long_int )
      if( size( b, 2, kind=long_int ) /= nb ) goto 200
      if( size( ab, 1, kind=long_int ) /= &
           size( ab, 2, kind=long_int ) ) goto 200
      if( size( ab, 1, kind=long_int ) /= na * nb ) goto 300
      ab(:,:) = 0.D0
      do ia = 1, na
         do ja = 1, ia
            i = ( ia - 1 ) * nb
            do ib = 1, nb
               i = i + 1
               j = ( ja - 1 ) * nb
               do jb = 1, ib
                  j = j + 1
                  ab( i, j ) = a( ia, ja ) * b( ib, jb )
               end do
            end do
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps 
200   call err_handle(err, 1, &
           comment = "Non-square matrix encountered; square expected" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
300   call err_handle(err, 1, &
           comment = "Dimensions of matrix arguments not conformable" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function kronecker_lower_tri
  !####################################################################
   integer(kind=our_int) function kronecker(a, b, ab, err) &
        result(answer)
      !### Kronecker product of two arbitrary matrices
      implicit none
      ! declare arguments
      real(kind=our_dble), intent(in) :: a(:,:), b(:,:)
      real(kind=our_dble), intent(out) :: ab(:,:)
      type(error_type), intent(inout) :: err
      ! declare local variables and parameters
      character(len=*), parameter :: subname = "kronecker"
      integer(kind=long_int) :: ma, na, mb, nb, &
           ia, ja, ib, jb, i, j
      ! begin
      answer = RETURN_FAIL
      ma = size( a, 1, kind=long_int )
      na = size( a, 2, kind=long_int )
      mb = size( b, 1, kind=long_int )
      nb = size( b, 2, kind=long_int )
      if( size( ab, 1, kind=long_int ) /= ma * mb ) goto 300
      if( size( ab, 2, kind=long_int ) /= na * nb ) goto 300
      do ia = 1, ma
         do ja = 1, na
            i = ( ia - 1 ) * mb
            do ib = 1, mb
               i = i + 1
               j = ( ja - 1 ) * nb
               do jb = 1, nb
                  j = j + 1
                  ab( i, j ) = a( ia, ja ) * b( ib, jb )
               end do
            end do
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps 
300   call err_handle(err, 1, &
           comment = "Dimensions of matrix arguments not conformable" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function kronecker
  !####################################################################
end module matrix_methods
!#######################################################################

!integrate: a function to return the integral of x^(exps[0])*y^(exps[1])*z^(exps[2]) over a given sample shape
real(8) function integrate(exps,shape,dims)
  implicit none
  integer, parameter :: dp=kind(0.d0)
  integer, intent(in), dimension(0:2) :: exps
  integer, intent(in) :: shape
  real(dp), intent(in), dimension(0:2) :: dims
  
  real(dp) :: Lx, Ly, Lz
  integer, dimension(0:2) :: exps_1

  !variables needed for qawo (numerical integration; part of QUADPACK)
  !IMPORTANT: quadpack only works with (kind = 4), NOT double precision
  integer ::  ier, neval
  real(kind=4) :: abserr, result, x

  !made f01 contained within function, rather than calling it as an external function
  !real(8), external :: f01


  !rectangular parallelepiped
  if (shape == 1) then
     Lx = dims(0)
     Ly = dims(1)
     Lz = dims(2)
    if (mod(exps(0),2) == 1 .OR. mod(exps(1),2) == 1 .OR. mod(exps(2),2) == 1) then
       integrate = 0     !odd powers will give an integral of 0 across the RPR
    else
       exps_1 = exps+1
       integrate = 8.0*Lx**(exps_1(0))*Ly**(exps_1(1))*Lz**(exps_1(2))/(exps_1(0)*exps_1(1)*exps_1(2))
    end if
  !diamond-faced prism
  else if (shape == 2) then
      Lx = dims(0)
      Ly = dims(1)
      Lz = dims(2)
      if (mod(exps(0),2) == 1 .OR. mod(exps(1),2) == 1 .OR. mod(exps(2),2) == 1) then
        integrate = 0     !odd powers will give an integral of 0 across the diamond
      else
        exps_1 = exps+1
        integrate = 8.0*Lx**(exps_1(0))*Ly**(exps_1(1))*Lz**(exps_1(2))&
             *Gamma(exps(0)+1.0)*Gamma(exps(1)+1.0)/((exps(2)+1.0)*Gamma(3.0+exps(0)+exps(1)))
      end if
  !cylinder
  else if (shape == 3) then
     Lx = dims(0) !radius
     Lz = dims(2) !half-height
     if (mod(exps(0),2) == 1 .OR. mod(exps(1),2) == 1 .OR. mod(exps(2),2) == 1) then
        integrate = 0     !odd powers will give an integral of 0 across the cylinder
     else
        exps_1 = exps+1
        !call qawo(function, lowerbound, upperbound, omega, integer, epsabs, epsrel, result, abserr, neval, ier)
        !call qawo (f01, 0.0E+00, 2*3.141592653589793E+00, 1.0E+00, 2, 0.0E+00, 0.001E+00, result, abserr, neval, ier )

        !important: upperbound, lowerbound, etc. must all be of the form n.nE+00
        call qags(f01, 0.0E+00, 2*3.141592653589793E+00, 0.0E+00, 0.001E+00, result, abserr, neval, ier)
        !^integral of cos(theta)^l*sin(theta)^m from 0 to 2pi

        integrate = 2*Lz**(exps_1(2))*Lx**(exps_1(0)+exps_1(1))/((exps_1(0)+exps_1(1))*exps_1(2))*result

     end if
  end if

  contains
    real(kind=4) function f01(x)
      real(kind=4) :: x
      f01 = cos(x)**exps(0)*sin(x)**exps(1)
      !previously had the following, which I believe was just a mistake:
      !f01 = cos(x)**exps(0)*sin(x)**(exps(1)-1)
      
      !f01 = cos(x)
      !if ( x <= 0.0E+00 ) then
      !   f01 = 0.0E+00
      !else
      !   f01 = log ( x ) / sqrt ( x )
      !end if

      return
    end function f01

end function integrate



!tried to make creation of E and Gamma as efficient as possible:
!put E into same do-loop as Gamma;
!moved multiplication by rho OUTSIDE of loops;
!initialized E to 0 so didn't need to set any elements individually to 0;
!made it so we only calculate a term for Gamma if C not equal to 0 for that element;
!taking advantage of the fact that E and Gamma are symmetric: build lower-left triangles of each, then copy over.

!build the matrices E and Gamma used in the eigenvalue problem, then solve for eigenvalues and eigenvectors
subroutine forwardsolver(basis,mass,C,shape,dims,eigenvals,eigenvects,lenbasis)
  implicit none
  integer, parameter :: dp=kind(0.d0)     !used for defining "float" variables
  integer, intent(in) :: lenbasis         !the length of the basis function
  integer, intent(in), dimension(0:lenbasis-1,0:2) :: basis
  real(dp), intent(in) :: mass
  real(dp), intent(in), dimension(0:2,0:2,0:2,0:2) :: C
  integer, intent(in) :: shape
  real(dp), intent(in), dimension(0:2) :: dims

  real(dp) :: integrate

  integer, dimension(0:2) :: ooo    !the matrix [0,0,0]; used to find volume
  real(dp) :: volume, rho
  integer :: ii, kk, i, k, j, l     !indices for "do" loops
  integer, dimension(0:3*lenbasis-1,0:2) :: multbasis    !the array "basis" repeated x3
  integer, dimension(0:3*lenbasis-1) :: component        !helps us to construct E, Gamma
  integer, dimension(0:2) :: exps   ![l,m,n]
  real(dp) :: integral, integral_sum
  real(dp), dimension(0:3*lenbasis-1,0:3*lenbasis-1) :: E, Gamma, ET, GammaT

  real(dp), dimension(102*lenbasis) :: WORK     !things used by the eigensolver; dsygv tells me that this is the optimum size
  integer :: INFO

  real(dp), intent(out), dimension(0:3*lenbasis-1) :: eigenvals
  real(dp), intent(out), dimension(0:3*lenbasis-1,0:3*lenbasis-1) :: eigenvects
  

  !by calculating rho in buildmatrices, I can make "integrate" only called from here.
  ooo(:) = 0
  volume = integrate(ooo,shape,dims)  !volume of sample in cm^3
  rho = mass/volume              !density of sample in g/cm^3

  multbasis(0:lenbasis-1,:) = basis
  multbasis(lenbasis:2*lenbasis-1,:) = basis
  multbasis(2*lenbasis:3*lenbasis-1,:) = basis

  component(0:lenbasis-1) = 0
  component(lenbasis:2*lenbasis-1) = 1
  component(2*lenbasis:3*lenbasis-1) = 2
   
  !create E and Gamma (only lower-left triangles; will use fact that they're both symmetric)
  E(:,:) = 0
  
  do kk = 0,3*lenbasis-1
     do ii = kk,3*lenbasis-1
      i = component(ii)  !0, 1, or 2
      k = component(kk)  !0, 1, or 2

      !create E
      if (i == k) then
        exps = multbasis(ii,:)+multbasis(kk,:)
        E(ii,kk) = integrate(exps,shape,dims)
      end if

      !create Gamma
      integral_sum = 0
      do j = 0,2
        do l = 0,2
          if (multbasis(ii,j) /= 0 .AND. multbasis(kk,l) /= 0 .AND. C(i,j,k,l) /= 0) then
            exps = multbasis(ii,:)+multbasis(kk,:)
            exps(j) = exps(j) - 1 !derivative on multbasis[ii]
            exps(l) = exps(l) - 1 !derivative on multbasis[kk]
            integral =  C(i,j,k,l)*multbasis(ii,j)*multbasis(kk,l)*integrate(exps,shape,dims)
            integral_sum = integral_sum + integral
          end if
        end do
      end do
      Gamma(ii,kk) = integral_sum
    end do
  end do

  !creating upper-right portions of matrices; using transposes to access elements most efficiently
  GammaT = transpose(Gamma)
  ET = transpose(E)
  do kk = 0,3*lenbasis-1
     do ii = 0,kk-1
        Gamma(ii,kk) = GammaT(ii,kk)
        E(ii,kk) = ET(ii,kk)
     end do
  end do
  
  !overall factor for the matrix; multiplying at the end for efficiency
  E = rho*E 

  !call rsg(3*lenbasis,3*lenbasis,Gamma,E,eigenvals,1,eigenvects,fv1,fv2,ierr) (obsolete)
  call dsygv(1,'V','L',3*lenbasis,Gamma,3*lenbasis,E,3*lenbasis,eigenvals,WORK,102*lenbasis,INFO)
  
  eigenvects = Gamma

end subroutine


!build only the matrices Gamma and E (need this subroutine for derivatives of matrices with respect to parameters)
subroutine matrixderivs(basis,mass,C,shape,dims,switch,derivnum,Gamma,E,lenbasis)
  implicit none
  integer, parameter :: dp=kind(0.d0)     !used for defining "float" variables
  integer, intent(in) :: lenbasis         !the length of the basis function
  integer, intent(in), dimension(0:lenbasis-1,0:2) :: basis
  real(dp), intent(in) :: mass
  real(dp), intent(in), dimension(0:2,0:2,0:2,0:2) :: C
  integer, intent(in) :: shape
  real(dp), intent(in), dimension(0:2) :: dims   !the dimensions of the sample
  integer, intent(in) :: switch    !0 for derivatives w.r.t. elastic moduli or alignment; 1 for derivatives of dimensions
  integer, intent(in) :: derivnum  !which dimension we're taking the derivative of (0,1,or 2)

  real(dp) :: integrate

  integer, dimension(0:2) :: ooo    !the matrix [0,0,0]; used to find volume
  real(dp) :: volume, rho
  integer :: ii, kk, i, k, j, l     !indices for "do" loops
  integer, dimension(0:3*lenbasis-1,0:2) :: multbasis    !the array "basis" repeated x3
  integer, dimension(0:3*lenbasis-1) :: component        !helps us to construct E, Gamma
  integer, dimension(0:2) :: exps   ![l,m,n]
  real(dp) :: integral, integral_sum
  real(dp), dimension(0:3*lenbasis-1,0:3*lenbasis-1) :: ET, GammaT

  real(dp), intent(out), dimension(0:3*lenbasis-1,0:3*lenbasis-1) :: Gamma, E

  !by calculating rho in buildmatrices, I can make "integrate" only called from here.
  ooo(:) = 0
  volume = integrate(ooo,shape,dims)  !volume of sample in cm^3
  rho = mass/volume              !density of sample in g/cm^3

  multbasis(0:lenbasis-1,:) = basis
  multbasis(lenbasis:2*lenbasis-1,:) = basis
  multbasis(2*lenbasis:3*lenbasis-1,:) = basis

  component(0:lenbasis-1) = 0
  component(lenbasis:2*lenbasis-1) = 1
  component(2*lenbasis:3*lenbasis-1) = 2

  !create E and Gamma (only lower-left triangles; will use fact that they're both symmetric)
  E(:,:) = 0

  do kk = 0,3*lenbasis-1
     do ii = kk,3*lenbasis-1
      i = component(ii)  !0, 1, or 2
      k = component(kk)  !0, 1, or 2

      !create E derivative (all 0's for derivatives w.r.t. elastic moduli or alignment)
      if (switch == 1) then
         if (i == k) then
            exps = multbasis(ii,:)+multbasis(kk,:)
            E(ii,kk) = integrate(exps,shape,dims)*(exps(derivnum)+1)/dims(derivnum)
         end if
      end if

      !create Gamma derivative
      integral_sum = 0
      do j = 0,2
        do l = 0,2
          if (multbasis(ii,j) /= 0 .AND. multbasis(kk,l) /= 0 .AND. C(i,j,k,l) /= 0) then
            exps = multbasis(ii,:)+multbasis(kk,:)
            exps(j) = exps(j) - 1 !derivative on multbasis[ii]
            exps(l) = exps(l) - 1 !derivative on multbasis[kk]
            if (switch == 0) then  !derivative of elastic moduli or alignment; done through input C
              integral =  C(i,j,k,l)*multbasis(ii,j)*multbasis(kk,l)*integrate(exps,shape,dims)
            else   !derivate of dimension; done "manually" here
              integral =  C(i,j,k,l)*multbasis(ii,j)*multbasis(kk,l)*integrate(exps,shape,dims)*(exps(derivnum)+1)/dims(derivnum)
            end if
            integral_sum = integral_sum + integral
          end if
        end do
      end do
      Gamma(ii,kk) = integral_sum
    end do
  end do

  !creating upper-right portion of Gamma; using transpose to access elements most efficiently
  GammaT = transpose(Gamma)
  ET = transpose(E)
  do kk = 0,3*lenbasis-1
     do ii = 0,kk-1
        Gamma(ii,kk) = GammaT(ii,kk)
        E(ii,kk) = ET(ii,kk)
     end do
  end do

  !overall factor for the matrix; multiplying at the end for efficiency
  E = rho*E 
  
end subroutine

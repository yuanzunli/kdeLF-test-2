!
!  f2py3 --fcompiler=intelem --f90flags=-openmp -liomp5 -c -m kde kdeLF.f90 --quiet
!
!  version: 2021_08_31_16_57 
!

!***********************************************************************************************************

! module dqag and dqags for f2py

!!!! version 2019_11_16_16_01

module dqag_and_dqags
implicit none
contains


subroutine dqage ( f, a, b, epsabs, epsrel, key, limit, result, abserr, &
  neval, ier, alist, blist, rlist, elist, iord, last )

!*****************************************************************************80
!
!! DQAGE estimates a definite integral.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral   i = integral of f over (a,b),
!      hopefully satisfying following claim for accuracy
!      abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!      key    - integer ( kind = 4 )
!               key for choice of local integration rule
!               a gauss-kronrod pair is used with
!                    7 - 15 points if key.lt.2,
!                   10 - 21 points if key = 2,
!                   15 - 31 points if key = 3,
!                   20 - 41 points if key = 4,
!                   25 - 51 points if key = 5,
!                   30 - 61 points if key.gt.5.
!
!      limit  - integer ( kind = 4 )
!               gives an upperbound on the number of subintervals
!               in the partition of (a,b), limit.ge.1.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for result and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more
!                       subdivisions by increasing the value
!                       of limit.
!                       however, if this yields no improvement it
!                       is rather advised to analyze the integrand
!                       in order to determine the integration
!                       difficulties. if the position of a local
!                       difficulty can be determined(e.g.
!                       singularity, discontinuity within the
!                       interval) one will probably gain from
!                       splitting up the interval at this point
!                       and calling the integrator on the
!                       subranges. if possible, an appropriate
!                       special-purpose integrator should be used
!                       which is designed for handling the type of
!                       difficulty involved.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 6 the input is invalid, because
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                       result, abserr, neval, last, rlist(1) ,
!                       elist(1) and iord(1) are set to zero.
!                       alist(1) and blist(1) are set to a and b
!                       respectively.
!
!      alist   - real ( kind = 8 )
!                vector of dimension at least limit, the first
!                 last  elements of which are the left
!                end points of the subintervals in the partition
!                of the given integration range (a,b)
!
!      blist   - real ( kind = 8 )
!                vector of dimension at least limit, the first
!                 last  elements of which are the right
!                end points of the subintervals in the partition
!                of the given integration range (a,b)
!
!      rlist   - real ( kind = 8 )
!                vector of dimension at least limit, the first
!                 last  elements of which are the
!                integral approximations on the subintervals
!
!      elist   - real ( kind = 8 )
!                vector of dimension at least limit, the first
!                 last  elements of which are the moduli of the
!                absolute error estimates on the subintervals
!
!      iord    - integer ( kind = 4 )
!                vector of dimension at least limit, the first k
!                elements of which are pointers to the
!                error estimates over the subintervals,
!                such that elist(iord(1)), ...,
!                elist(iord(k)) form a decreasing sequence,
!                with k = last if last.le.(limit/2+2), and
!                k = limit+1-last otherwise
!
!      last    - integer ( kind = 4 )
!                number of subintervals actually produced in the
!                subdivision process
!
!  Local Parameters:
!
!     alist     - list of left end points of all subintervals
!                 considered up to now
!     blist     - list of right end points of all subintervals
!                 considered up to now
!     rlist(i)  - approximation to the integral over
!                (alist(i),blist(i))
!     elist(i)  - error estimate applying to rlist(i)
!     maxerr    - pointer to the interval with largest
!                 error estimate
!     errmax    - elist(maxerr)
!     area      - sum of the integrals over the subintervals
!     errsum    - sum of the errors over the subintervals
!     errbnd    - requested accuracy max(epsabs,epsrel*
!                 abs(result))
!     *****1    - variable for the left subinterval
!     *****2    - variable for the right subinterval
!     last      - index for subdivision
!
!
!     machine dependent constants
!
!     epmach  is the largest relative spacing.
!     uflow  is the smallest positive magnitude.
!
  implicit none
  real(8) temp  !by Zunli for using f2py

  real ( kind = 8 ) a,abserr,alist,area,area1,area12,area2,a1,a2,b, &
    blist,b1,b2,defabs,defab1,defab2,elist,epmach, &
    epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f, &
    resabs,result,rlist,uflow
  integer ( kind = 4 ) ier,iord,iroff1,iroff2,k,key,keyf,last,limit, &
    maxerr, nrmax, neval

  dimension alist(limit),blist(limit),elist(limit),iord(limit), &
    rlist(limit)

  external f

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
!
!  test on validity of parameters
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0D+00
  elist(1) = 0.0D+00
  iord(1) = 0

  if(epsabs.le.0.0D+00.and. &
    epsrel.lt. max ( 0.5D+02*epmach,0.5d-28)) then
    ier = 6
    return
  end if
!
!  first approximation to the integral
!
  keyf = key
  if(key.le.0) keyf = 1
  if(key.ge.7) keyf = 6
  neval = 0
  if(keyf.eq.1) call dqk15(f,a,b,result,abserr,defabs,resabs)
  if(keyf.eq.2) call dqk21(f,a,b,result,abserr,defabs,resabs)
  if(keyf.eq.3) call dqk31(f,a,b,result,abserr,defabs,resabs)
  if(keyf.eq.4) call dqk41(f,a,b,result,abserr,defabs,resabs)
  if(keyf.eq.5) call dqk51(f,a,b,result,abserr,defabs,resabs)
  if(keyf.eq.6) call dqk61(f,a,b,result,abserr,defabs,resabs)
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
!
!  test on accuracy.
!
  errbnd =  max ( epsabs, epsrel* abs ( result ) )

  if(abserr.le.0.5D+02* epmach * defabs .and. &
    abserr.gt.errbnd) then
    ier = 2
  end if

  if(limit.eq.1) then
    ier = 1
  end if

  if ( ier .ne. 0 .or. &
    (abserr .le. errbnd .and. abserr .ne. resabs ) .or. &
    abserr .eq. 0.0D+00 ) then

    if(keyf.ne.1) then
      neval = (10*keyf+1)*(2*neval+1)
    else
      neval = 30*neval+15
    end if

    return

  end if
!
!  initialization
!
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  iroff1 = 0
  iroff2 = 0
!
!  main do-loop
!
  do last = 2, limit
!
!  bisect the subinterval with the largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5D+00*(alist(maxerr)+blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)

    if(keyf.eq.1) call dqk15(f,a1,b1,area1,error1,resabs,defab1)
    if(keyf.eq.2) call dqk21(f,a1,b1,area1,error1,resabs,defab1)
    if(keyf.eq.3) call dqk31(f,a1,b1,area1,error1,resabs,defab1)
    if(keyf.eq.4) call dqk41(f,a1,b1,area1,error1,resabs,defab1)
    if(keyf.eq.5) call dqk51(f,a1,b1,area1,error1,resabs,defab1)
    if(keyf.eq.6) call dqk61(f,a1,b1,area1,error1,resabs,defab1)

    if(keyf.eq.1) call dqk15(f,a2,b2,area2,error2,resabs,defab2)
    if(keyf.eq.2) call dqk21(f,a2,b2,area2,error2,resabs,defab2)
    if(keyf.eq.3) call dqk31(f,a2,b2,area2,error2,resabs,defab2)
    if(keyf.eq.4) call dqk41(f,a2,b2,area2,error2,resabs,defab2)
    if(keyf.eq.5) call dqk51(f,a2,b2,area2,error2,resabs,defab2)
    if(keyf.eq.6) call dqk61(f,a2,b2,area2,error2,resabs,defab2)
!
!  improve previous approximations to integral
!  and error and test for accuracy.
!
    neval = neval+1
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)

    if ( defab1 .ne. error1 .and. defab2 .ne. error2 ) then

      if( abs ( rlist(maxerr)-area12).le.0.1D-04* abs ( area12) &
        .and. erro12.ge.0.99D+00*errmax) then
        iroff1 = iroff1+1
      end if

      if(last.gt.10.and.erro12.gt.errmax) then
        iroff2 = iroff2+1
      end if

    end if

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd =  max ( epsabs,epsrel* abs ( area))

    if ( errbnd .lt. errsum ) then
!
!  test for roundoff error and eventually set error flag.
!
      if(iroff1.ge.6.or.iroff2.ge.20) then
        ier = 2
      end if
!
!  set error flag in the case that the number of subintervals
!  equals limit.
!
      if(last.eq.limit) then
        ier = 1
      end if
!
!  set error flag in the case of bad integrand behaviour
!  at a point of the integration range.
!
      if( max (  abs ( a1), abs ( b2)).le.(0.1D+01+0.1D+03* &
        epmach)*( abs ( a2)+0.1D+04*uflow)) then
        ier = 3
      end if

    end if
!
!  append the newly-created intervals to the list.
!
    if(error2.le.error1) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  call dqpsrt to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with the largest error estimate (to be bisected next).
!
    call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)

    if(ier.ne.0.or.errsum.le.errbnd) then
      exit
    end if

  end do
!
!  compute final result.
!
  result = 0.0D+00
  do k=1,last
    result = result+rlist(k)
  end do
  abserr = errsum

  if(keyf.ne.1) then
    neval = (10*keyf+1)*(2*neval+1)
  else
    neval = 30*neval+15
  end if

  if(temp==1.0) temp=f(a)  !by Zunli for using f2py
  return
end subroutine dqage

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine dqag ( f, a, b, epsabs, epsrel, key, result, abserr, neval, ier, &
  limit, lenw, last, iwork, work )

!*****************************************************************************80
!
!! DQAG approximates an integral over a finite interval.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral i = integral of f over (a,b),
!      hopefully satisfying following claim for accuracy
!      abs(i-result)le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!      f      - real ( kind = 8 )
!               function subprogam defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      epsabs - real ( kind = 8 )
!               absolute accoracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!      key    - integer ( kind = 4 )
!               key for choice of local integration rule
!               a gauss-kronrod pair is used with
!                 7 - 15 points if key.lt.2,
!                10 - 21 points if key = 2,
!                15 - 31 points if key = 3,
!                20 - 41 points if key = 4,
!                25 - 51 points if key = 5,
!                30 - 61 points if key.gt.5.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for result and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!                error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more
!                       subdivisions by increasing the value of
!                       limit (and taking the according dimension
!                       adjustments into account). however, if
!                       this yield no improvement it is advised
!                       to analyze the integrand in order to
!                       determine the integration difficulaties.
!                       if the position of a local difficulty can
!                       be determined (i.e.singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling the
!                       integrator on the subranges. if possible,
!                       an appropriate special-purpose integrator
!                       should be used which is designed for
!                       handling the type of difficulty involved.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 6 the input is invalid, because
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                       or limit.lt.1 or lenw.lt.limit*4.
!                       result, abserr, neval, last are set
!                       to zero.
!                       except when lenw is invalid, iwork(1),
!                       work(limit*2+1) and work(limit*3+1) are
!                       set to zero, work(1) is set to a and
!                       work(limit+1) to b.
!
!   dimensioning parameters
!      limit - integer ( kind = 4 )
!              dimensioning parameter for iwork
!              limit determines the maximum number of subintervals
!              in the partition of the given integration interval
!              (a,b), limit.ge.1.
!              if limit.lt.1, the routine will end with ier = 6.
!
!      lenw  - integer ( kind = 4 )
!              dimensioning parameter for work
!              lenw must be at least limit*4.
!              if lenw.lt.limit*4, the routine will end with
!              ier = 6.
!
!      last  - integer ( kind = 4 )
!              on return, last equals the number of subintervals
!              produced in the subdivision process, which
!              determines the number of significant elements
!              actually in the work arrays.
!
!   work arrays
!      iwork - integer ( kind = 4 )
!              vector of dimension at least limit, the first k
!              elements of which contain pointers to the error
!              estimates over the subintervals, such that
!              work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
!              form a decreasing sequence with k = last if
!              last.le.(limit/2+2), and k = limit+1-last otherwise
!
!      work  - real ( kind = 8 )
!              vector of dimension at least lenw
!              on return
!              work(1), ..., work(last) contain the left end
!              points of the subintervals in the partition of
!               (a,b),
!              work(limit+1), ..., work(limit+last) contain the
!               right end points,
!              work(limit*2+1), ..., work(limit*2+last) contain
!               the integral approximations over the subintervals,
!              work(limit*3+1), ..., work(limit*3+last) contain
!               the error estimates.
!
  implicit none
  real(8) temp  !by Zunli for using f2py


  integer ( kind = 4 ) lenw
  integer ( kind = 4 ) limit

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) f
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iwork(limit)
  integer ( kind = 4 ) key
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lvl
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) l3
  integer ( kind = 4 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ) work(lenw)
  external f
!
!  check validity of lenw.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  if(limit.lt.1.or.lenw.lt.limit*4) go to 10
!
!  prepare call for dqage.
!
  l1 = limit+1
  l2 = limit+l1
  l3 = limit+l2

  call dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval, &
    ier,work(1),work(l1),work(l2),work(l3),iwork,last)
!
!  call error handler if necessary.
!
  lvl = 0
10    continue

  if(ier.eq.6) lvl = 1
  if(ier.ne.0) call xerror('abnormal return from dqag ',26,ier,lvl)

  if(temp==1.0) temp=f(a)  !by Zunli for using f2py
  return
end subroutine dqag

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dqagse(f,a,b,epsabs,epsrel,limit,result,abserr,neval, &
     ier,alist,blist,rlist,elist,iord,last)

!*****************************************************************************80
!
!! DQAGSE estimates the integral of a function.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral i = integral of f over (a,b),
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!      limit  - integer ( kind = 4 )
!               gives an upperbound on the number of subintervals
!               in the partition of (a,b)
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for integral and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!                   = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more sub-
!                       divisions by increasing the value of limit
!                       (and taking the according dimension
!                       adjustments into account). however, if
!                       this yields no improvement it is advised
!                       to analyze the integrand in order to
!                       determine the integration difficulties. if
!                       the position of a local difficulty can be
!                       determined (e.g. singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling the
!                       integrator on the subranges. if possible,
!                       an appropriate special-purpose integrator
!                       should be used, which is designed for
!                       handling the type of difficulty involved.
!                   = 2 the occurrence of roundoff error is detec-
!                       ted, which prevents the requested
!                       tolerance from being achieved.
!                       the error may be under-estimated.
!                   = 3 extremely bad integrand behaviour
!                       occurs at some points of the integration
!                       interval.
!                   = 4 the algorithm does not converge.
!                       roundoff error is detected in the
!                       extrapolation table.
!                       it is presumed that the requested
!                       tolerance cannot be achieved, and that the
!                       returned result is the best which can be
!                       obtained.
!                   = 5 the integral is probably divergent, or
!                       slowly convergent. it must be noted that
!                       divergence can occur with any other value
!                       of ier.
!                   = 6 the input is invalid, because
!                       epsabs.le.0 and
!                       epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
!                       result, abserr, neval, last, rlist(1),
!                       iord(1) and elist(1) are set to zero.
!                       alist(1) and blist(1) are set to a and b
!                       respectively.
!
!      alist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the left end points
!               of the subintervals in the partition of the
!               given integration range (a,b)
!
!      blist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the right end points
!               of the subintervals in the partition of the given
!               integration range (a,b)
!
!      rlist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the integral
!               approximations on the subintervals
!
!      elist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the moduli of the
!               absolute error estimates on the subintervals
!
!      iord   - integer ( kind = 4 )
!               vector of dimension at least limit, the first k
!               elements of which are pointers to the
!               error estimates over the subintervals,
!               such that elist(iord(1)), ..., elist(iord(k))
!               form a decreasing sequence, with k = last
!               if last.le.(limit/2+2), and k = limit+1-last
!               otherwise
!
!      last   - integer ( kind = 4 )
!               number of subintervals actually produced in the
!               subdivision process
!
!  Local parameters:
!
!      the dimension of rlist2 is determined by the value of
!      limexp in routine dqelg (rlist2 should be of dimension
!      (limexp+2) at least).
!
!      list of major variables
!
!     alist     - list of left end points of all subintervals
!                 considered up to now
!     blist     - list of right end points of all subintervals
!                 considered up to now
!     rlist(i)  - approximation to the integral over
!                 (alist(i),blist(i))
!     rlist2    - array of dimension at least limexp+2 containing
!                 the part of the epsilon table which is still
!                 needed for further computations
!     elist(i)  - error estimate applying to rlist(i)
!     maxerr    - pointer to the interval with largest error
!                 estimate
!     errmax    - elist(maxerr)
!     erlast    - error on the interval currently subdivided
!                 (before that subdivision has taken place)
!     area      - sum of the integrals over the subintervals
!     errsum    - sum of the errors over the subintervals
!     errbnd    - requested accuracy max(epsabs,epsrel*
!                 abs(result))
!     *****1    - variable for the left interval
!     *****2    - variable for the right interval
!     last      - index for subdivision
!     nres      - number of calls to the extrapolation routine
!     numrl2    - number of elements currently in rlist2. if an
!                 appropriate approximation to the compounded
!                 integral has been obtained it is put in
!                 rlist2(numrl2) after numrl2 has been increased
!                 by one.
!     small     - length of the smallest interval considered up
!                 to now, multiplied by 1.5
!     erlarg    - sum of the errors over the intervals larger
!                 than the smallest interval considered up to now
!     extrap    - logical variable denoting that the routine is
!                 attempting to perform extrapolation i.e. before
!                 subdividing the smallest interval we try to
!                 decrease the value of erlarg.
!     noext     - logical variable denoting that extrapolation
!                 is no longer allowed (true value)
!
!      machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!     oflow is the largest positive magnitude.
!
  implicit none
  real(8) temp  !by Zunli for using f2py

  real ( kind = 8 ) a,abseps,abserr,alist,area,area1,area12,area2,a1, &
    a2,b,blist,b1,b2,correc,defabs,defab1,defab2, &
    dres,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd,errmax, &
    error1,error2,erro12,errsum,ertest,f,oflow,resabs,reseps,result, &
    res3la,rlist,rlist2,small,uflow
  integer ( kind = 4 ) id,ier,ierro,iord,iroff1,iroff2,iroff3,jupbnd, &
    k,ksgn,ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
  logical extrap,noext
  dimension alist(limit),blist(limit),elist(limit),iord(limit), &
   res3la(3),rlist(limit),rlist2(52)

  external f

  epmach = epsilon ( epmach )
!
!  test on validity of parameters
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0D+00
  elist(1) = 0.0D+00

  if(epsabs.le.0.0D+00.and.epsrel.lt. max ( 0.5D+02*epmach,0.5d-28)) then
    ier = 6
    return
  end if
!
!  first approximation to the integral
!
  uflow = tiny ( uflow )
  oflow = huge ( oflow )
  ierro = 0
  call dqk21(f,a,b,result,abserr,defabs,resabs)
!
!  test on accuracy.
!
  dres =  abs ( result)
  errbnd =  max ( epsabs,epsrel*dres)
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  if(abserr.le.1.0D+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
  if(limit.eq.1) ier = 1
  if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or. &
    abserr.eq.0.0D+00) go to 140
!
!  initialization
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = oflow
  nrmax = 1
  nres = 0
  numrl2 = 2
  ktmin = 0
  extrap = .false.
  noext = .false.
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ksgn = -1
  if(dres.ge.(0.1D+01-0.5D+02*epmach)*defabs) ksgn = 1
!
!  main do-loop
!
  do 90 last = 2,limit
!
!  bisect the subinterval with the nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5D+00*(alist(maxerr)+blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call dqk21(f,a1,b1,area1,error1,resabs,defab1)
    call dqk21(f,a2,b2,area2,error2,resabs,defab2)
!
!  improve previous approximations to integral
!  and error and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)
    if(defab1.eq.error1.or.defab2.eq.error2) go to 15
    if( abs ( rlist(maxerr)-area12).gt.0.1D-04* abs ( area12) &
    .or.erro12.lt.0.99D+00*errmax) go to 10
    if(extrap) iroff2 = iroff2+1
    if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
    rlist(last) = area2
    errbnd =  max ( epsabs,epsrel* abs ( area))
!
!  test for roundoff error and eventually set error flag.
!
    if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
    if(iroff2.ge.5) ierro = 3
!
!  set error flag in the case that the number of subintervals
!  equals limit.
!
    if(last.eq.limit) ier = 1
!
!  set error flag in the case of bad integrand behaviour
!  at a point of the integration range.
!
    if( max (  abs ( a1), abs ( b2)).le.(0.1D+01+0.1D+03*epmach)* &
    ( abs ( a2)+0.1D+04*uflow)) ier = 4
!
!  append the newly-created intervals to the list.
!
    if(error2.gt.error1) go to 20
    alist(last) = a2
    blist(maxerr) = b1
    blist(last) = b2
    elist(maxerr) = error1
    elist(last) = error2
    go to 30
   20   alist(maxerr) = a2
    alist(last) = a1
    blist(last) = b1
    rlist(maxerr) = area2
    rlist(last) = area1
    elist(maxerr) = error2
    elist(last) = error1
!
!  call dqpsrt to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
   30   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
    if(errsum.le.errbnd) go to 115
    if(ier.ne.0) go to 100
    if(last.eq.2) go to 80
    if(noext) go to 90
    erlarg = erlarg-erlast
    if( abs ( b1-a1).gt.small) erlarg = erlarg+erro12
    if(extrap) go to 40
!
!  test whether the interval to be bisected next is the
!  smallest interval.
!
    if( abs ( blist(maxerr)-alist(maxerr)).gt.small) go to 90
    extrap = .true.
    nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
!
!  the smallest interval has the largest error.
!  before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    id = nrmax
    jupbnd = last
    if(last.gt.(2+limit/2)) jupbnd = limit+3-last
    do k = id,jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if( abs ( blist(maxerr)-alist(maxerr)).gt.small) go to 90
      nrmax = nrmax+1
    end do
!
!  perform extrapolation.
!
   60   numrl2 = numrl2+1
    rlist2(numrl2) = area
    call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
    ktmin = ktmin+1
    if(ktmin.gt.5.and.abserr.lt.0.1D-02*errsum) ier = 5
    if(abseps.ge.abserr) go to 70
    ktmin = 0
    abserr = abseps
    result = reseps
    correc = erlarg
    ertest =  max ( epsabs,epsrel* abs ( reseps))
    if(abserr.le.ertest) go to 100
!
!  prepare bisection of the smallest interval.
!
   70   if(numrl2.eq.1) noext = .true.
    if(ier.eq.5) go to 100
    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small*0.5D+00
    erlarg = errsum
    go to 90
   80   small =  abs ( b-a)*0.375D+00
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area
   90 continue
!
!  set final result and error estimate.
!
  100 if(abserr.eq.oflow) go to 115
  if(ier+ierro.eq.0) go to 110
  if(ierro.eq.3) abserr = abserr+correc
  if(ier.eq.0) ier = 3
  if(result.ne.0.0D+00.and.area.ne.0.0D+00) go to 105
  if(abserr.gt.errsum) go to 115
  if(area.eq.0.0D+00) go to 130
  go to 110
  105 if(abserr/ abs ( result).gt.errsum/ abs ( area)) go to 115
!
!  test on divergence.
!
  110 if(ksgn.eq.(-1).and. max (  abs ( result), abs ( area)).le. &
   defabs*0.1D-01) go to 130
  if(0.1D-01.gt.(result/area).or.(result/area).gt.0.1D+03 &
   .or.errsum.gt. abs ( area)) ier = 6
  go to 130
!
!  compute global integral sum.
!
  115 result = 0.0D+00
  do k = 1,last
     result = result+rlist(k)
  end do
  abserr = errsum
  130 if(ier.gt.2) ier = ier-1
  140 neval = 42*last-21

  if(temp==1.0) temp=f(a)  !by Zunli for using f2py
  return
end subroutine dqagse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine dqags ( f, a, b, epsabs, epsrel, result, abserr, neval, ier, &
  limit, lenw, last, iwork, work )

!*****************************************************************************80
!
!! DQAGS estimates the integral of a function.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral  i = integral of f over (a,b),
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for integral and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more sub-
!                       divisions by increasing the value of limit
!                       (and taking the according dimension
!                       adjustments into account. however, if
!                       this yields no improvement it is advised
!                       to analyze the integrand in order to
!                       determine the integration difficulties. if
!                       the position of a local difficulty can be
!                       determined (e.g. singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling the
!                       integrator on the subranges. if possible,
!                       an appropriate special-purpose integrator
!                       should be used, which is designed for
!                       handling the type of difficulty involved.
!                   = 2 the occurrence of roundoff error is detec-
!                       ted, which prevents the requested
!                       tolerance from being achieved.
!                       the error may be under-estimated.
!                   = 3 extremely bad integrand behaviour
!                       occurs at some points of the integration
!                       interval.
!                   = 4 the algorithm does not converge.
!                       roundoff error is detected in the
!                       extrapolation table. it is presumed that
!                       the requested tolerance cannot be
!                       achieved, and that the returned result is
!                       the best which can be obtained.
!                   = 5 the integral is probably divergent, or
!                       slowly convergent. it must be noted that
!                       divergence can occur with any other value
!                       of ier.
!                   = 6 the input is invalid, because
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28)
!                       or limit.lt.1 or lenw.lt.limit*4.
!                       result, abserr, neval, last are set to
!                       zero.except when limit or lenw is invalid,
!                       iwork(1), work(limit*2+1) and
!                       work(limit*3+1) are set to zero, work(1)
!                       is set to a and work(limit+1) to b.
!
!   dimensioning parameters
!      limit - integer ( kind = 4 )
!              dimensioning parameter for iwork
!              limit determines the maximum number of subintervals
!              in the partition of the given integration interval
!              (a,b), limit.ge.1.
!              if limit.lt.1, the routine will end with ier = 6.
!
!      lenw  - integer ( kind = 4 )
!              dimensioning parameter for work
!              lenw must be at least limit*4.
!              if lenw.lt.limit*4, the routine will end
!              with ier = 6.
!
!      last  - integer ( kind = 4 )
!              on return, last equals the number of subintervals
!              produced in the subdivision process, detemines the
!              number of significant elements actually in the work
!              arrays.
!
!   work arrays
!      iwork - integer ( kind = 4 )
!              vector of dimension at least limit, the first k
!              elements of which contain pointers
!              to the error estimates over the subintervals
!              such that work(limit*3+iwork(1)),... ,
!              work(limit*3+iwork(k)) form a decreasing
!              sequence, with k = last if last.le.(limit/2+2),
!              and k = limit+1-last otherwise
!
!      work  - real ( kind = 8 )
!              vector of dimension at least lenw
!              on return
!              work(1), ..., work(last) contain the left
!               end-points of the subintervals in the
!               partition of (a,b),
!              work(limit+1), ..., work(limit+last) contain
!               the right end-points,
!              work(limit*2+1), ..., work(limit*2+last) contain
!               the integral approximations over the subintervals,
!              work(limit*3+1), ..., work(limit*3+last)
!               contain the error estimates.
!
  implicit none
  real(8) temp  !by Zunli for using f2py


  real ( kind = 8 ) a,abserr,b,epsabs,epsrel,f,result,work
  integer ( kind = 4 ) ier,iwork,last,lenw,limit,lvl,l1,l2,l3,neval
  dimension iwork(limit),work(lenw)

  external f
!
!  check validity of limit and lenw.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  if(limit.lt.1.or.lenw.lt.limit*4) go to 10
!
!  prepare call for dqagse.
!
  l1 = limit+1
  l2 = limit+l1
  l3 = limit+l2

  call dqagse(f,a,b,epsabs,epsrel,limit,result,abserr,neval, &
    ier,work(1),work(l1),work(l2),work(l3),iwork,last)
!
!  call error handler if necessary.
!
  lvl = 0
10    if(ier.eq.6) lvl = 1
  if(ier.ne.0) call xerror('abnormal return from dqags',26,ier,lvl)


  if(temp==1.0) temp=f(a)  !by Zunli for using f2py
  return
end subroutine dqags

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine dqelg ( n, epstab, result, abserr, res3la, nres )

!*****************************************************************************80
!
!! DQELG carries out the Epsilon extrapolation algorithm.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine determines the limit of a given sequence of
!      approximations, by means of the epsilon algorithm of
!      p.wynn. an estimate of the absolute error is also given.
!      the condensed epsilon table is computed. only those
!      elements needed for the computation of the next diagonal
!      are preserved.
!
!  Parameters:
!
!        n      - integer ( kind = 4 )
!                 epstab(n) contains the new element in the
!                 first column of the epsilon table.
!
!        epstab - real ( kind = 8 )
!                 vector of dimension 52 containing the elements
!                 of the two lower diagonals of the triangular
!                 epsilon table. the elements are numbered
!                 starting at the right-hand corner of the
!                 triangle.
!
!        result - real ( kind = 8 )
!                 resulting approximation to the integral
!
!        abserr - real ( kind = 8 )
!                 estimate of the absolute error computed from
!                 result and the 3 previous results
!
!        res3la - real ( kind = 8 )
!                 vector of dimension 3 containing the last 3
!                 results
!
!        nres   - integer ( kind = 4 )
!                 number of calls to the routine
!                 (should be zero at first call)
!
!  Local Parameters:
!
!     e0     - the 4 elements on which the computation of a new
!     e1       element in the epsilon table is based
!     e2
!     e3                 e0
!                  e3    e1    new
!                        e2
!     newelm - number of elements to be computed in the new
!              diagonal
!     error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
!     result - the element in the new diagonal with least value
!              of error
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     oflow is the largest positive magnitude.
!     limexp is the maximum number of elements the epsilon
!     table can contain. if this number is reached, the upper
!     diagonal of the epsilon table is deleted.
!
  implicit none

  real ( kind = 8 ) abserr,delta1,delta2,delta3, &
    epmach,epsinf,epstab,error,err1,err2,err3,e0,e1,e1abs,e2,e3, &
    oflow,res,result,res3la,ss,tol1,tol2,tol3
  integer ( kind = 4 ) i,ib,ib2,ie,indx,k1,k2,k3,limexp,n,newelm
  integer ( kind = 4 ) nres
  integer ( kind = 4 ) num
  dimension epstab(52),res3la(3)

  epmach = epsilon ( epmach )
  oflow = huge ( oflow )
  nres = nres+1
  abserr = oflow
  result = epstab(n)
  if(n.lt.3) go to 100
  limexp = 50
  epstab(n+2) = epstab(n)
  newelm = (n-1)/2
  epstab(n) = oflow
  num = n
  k1 = n

  do 40 i = 1,newelm

    k2 = k1-1
    k3 = k1-2
    res = epstab(k1+2)
    e0 = epstab(k3)
    e1 = epstab(k2)
    e2 = res
    e1abs =  abs ( e1)
    delta2 = e2-e1
    err2 =  abs ( delta2)
    tol2 =  max (  abs ( e2),e1abs)*epmach
    delta3 = e1 - e0
    err3 =  abs ( delta3)
    tol3 =  max ( e1abs, abs ( e0))*epmach
    if(err2.gt.tol2.or.err3.gt.tol3) go to 10
!
!  if e0, e1 and e2 are equal to machine accuracy, convergence is assumed.
!
    result = res
    abserr = err2+err3
    go to 100
   10   e3 = epstab(k1)
    epstab(k1) = e1
    delta1 = e1-e3
    err1 =  abs ( delta1)
    tol1 =  max ( e1abs, abs ( e3))*epmach
!
!  if two elements are very close to each other, omit
!  a part of the table by adjusting the value of n
!
    if(err1.le.tol1.or.err2.le.tol2.or.err3.le.tol3) go to 20
    ss = 0.1D+01/delta1+0.1D+01/delta2-0.1D+01/delta3
    epsinf =  abs ( ss*e1)
!
!  test to detect irregular behaviour in the table, and
!  eventually omit a part of the table adjusting the value
!  of n.
!
    if(epsinf.gt.0.1D-03) go to 30
   20   n = i+i-1
    go to 50
!
!  compute a new element and eventually adjust
!  the value of result.
!
   30   res = e1+0.1D+01/ss
    epstab(k1) = res
    k1 = k1-2
    error = err2 + abs ( res-e2 ) + err3

    if ( error .le. abserr ) then
      abserr = error
      result = res
    end if

   40 continue
!
!  shift the table.
!
   50 if(n.eq.limexp) n = 2*(limexp/2)-1
  ib = 1
  if((num/2)*2.eq.num) ib = 2
  ie = newelm+1
  do i=1,ie
    ib2 = ib+2
    epstab(ib) = epstab(ib2)
    ib = ib2
  end do
  if(num.eq.n) go to 80
  indx = num-n+1
  do i = 1,n
    epstab(i)= epstab(indx)
    indx = indx+1
  end do
   80 if(nres.ge.4) go to 90
  res3la(nres) = result
  abserr = oflow
  go to 100
!
!  compute error estimate
!
   90 abserr =  abs ( result-res3la(3))+ abs ( result-res3la(2)) &
    + abs ( result-res3la(1))
  res3la(1) = res3la(2)
  res3la(2) = res3la(3)
  res3la(3) = result
  100 continue

  abserr =  max ( abserr, 0.5D+01*epmach* abs ( result))

  return
end subroutine dqelg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dqk15(f,a,b,result,abserr,resabs,resasc)

!*****************************************************************************80
!
!! DQK15 carries out a 15 point Gauss-Kronrod quadrature rule.
!
!     the abscissae and weights are given for the interval (-1,1).
!     because of symmetry only the positive abscissae and their
!     corresponding weights are given.
!
!     xgk    - abscissae of the 15-point kronrod rule
!              xgk(2), xgk(4), ...  abscissae of the 7-point
!              gauss rule
!              xgk(1), xgk(3), ...  abscissae which are optimally
!              added to the 7-point gauss rule
!
!     wgk    - weights of the 15-point kronrod rule
!
!     wg     - weights of the 7-point gauss rule
!
!
!   gauss quadrature weights and kronron quadrature abscissae and weights
!   as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
!   bell labs, nov. 1981.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b), with error
!                     estimate
!                 j = integral of abs(f) over (a,b)
!  Parameters:
!
!      on entry
!        f      - real ( kind = 8 )
!                 function subprogram defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the calling program.
!
!        a      - real ( kind = 8 )
!                 lower limit of integration
!
!        b      - real ( kind = 8 )
!                 upper limit of integration
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 15-point
!                 kronrod rule (resk) obtained by optimal addition
!                 of abscissae to the7-point gauss rule(resg).
!
!        abserr - real ( kind = 8 )
!                 estimate of the modulus of the absolute error,
!                 which should not exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral j
!
!        resasc - real ( kind = 8 )
!                 approximation to the integral of abs(f-i/(b-a))
!                 over (a,b)
!
!  Local Parameters:
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc   - abscissa
!     fval*  - function value
!     resg   - result of the 7-point gauss formula
!     resk   - result of the 15-point kronrod formula
!     reskh  - approximation to the mean value of f over (a,b),
!              i.e. to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none
  real(8) temp  !by Zunli for using f2py

  
  real ( kind = 8 ) a,absc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f
  dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)

  data wg  (  1) / 0.129484966168869693270611432679082d0 /
  data wg  (  2) / 0.279705391489276667901467771423780d0 /
  data wg  (  3) / 0.381830050505118944950369775488975d0 /
  data wg  (  4) / 0.417959183673469387755102040816327d0 /

  data xgk (  1) / 0.991455371120812639206854697526329d0 /
  data xgk (  2) / 0.949107912342758524526189684047851d0 /
  data xgk (  3) / 0.864864423359769072789712788640926d0 /
  data xgk (  4) / 0.741531185599394439863864773280788d0 /
  data xgk (  5) / 0.586087235467691130294144838258730d0 /
  data xgk (  6) / 0.405845151377397166906606412076961d0 /
  data xgk (  7) / 0.207784955007898467600689403773245d0 /
  data xgk (  8) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.022935322010529224963732008058970d0 /
  data wgk (  2) / 0.063092092629978553290700663189204d0 /
  data wgk (  3) / 0.104790010322250183839876322541518d0 /
  data wgk (  4) / 0.140653259715525918745189590510238d0 /
  data wgk (  5) / 0.169004726639267902826583426598550d0 /
  data wgk (  6) / 0.190350578064785409913256402421014d0 /
  data wgk (  7) / 0.204432940075298892414161999234649d0 /
  data wgk (  8) / 0.209482141084727828012999174891714d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 15-point kronrod approximation to
!  the integral, and estimate the absolute error.
!
  fc = f(centr)
  resg = fc*wg(4)
  resk = fc*wgk(8)
  resabs =  abs ( resk)

  do j=1,3
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j = 1,4
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(8)* abs ( fc-reskh)
  do j=1,7
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) &
    abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)


  if(temp==1.0) temp=f(a)  !by Zunli for using f2py
  return
end subroutine dqk15

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dqk21(f,a,b,result,abserr,resabs,resasc)

!*****************************************************************************80
!
!! DQK21 carries out a 21 point Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b), with error
!                     estimate
!                 j = integral of abs(f) over (a,b)
!
!  Parameters:
!
!      on entry
!        f      - real ( kind = 8 )
!                 function subprogram defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the driver program.
!
!        a      - real ( kind = 8 )
!                 lower limit of integration
!
!        b      - real ( kind = 8 )
!                 upper limit of integration
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 21-point
!                 kronrod rule (resk) obtained by optimal addition
!                 of abscissae to the 10-point gauss rule (resg).
!
!        abserr - real ( kind = 8 )
!                 estimate of the modulus of the absolute error,
!                 which should not exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral j
!
!        resasc - real ( kind = 8 )
!                 approximation to the integral of abs(f-i/(b-a))
!                 over (a,b)
!
!  Local Parameters:
!
!
!     the abscissae and weights are given for the interval (-1,1).
!     because of symmetry only the positive abscissae and their
!     corresponding weights are given.
!
!     xgk    - abscissae of the 21-point kronrod rule
!              xgk(2), xgk(4), ...  abscissae of the 10-point
!              gauss rule
!              xgk(1), xgk(3), ...  abscissae which are optimally
!              added to the 10-point gauss rule
!
!     wgk    - weights of the 21-point kronrod rule
!
!     wg     - weights of the 10-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc   - abscissa
!     fval*  - function value
!     resg   - result of the 10-point gauss formula
!     resk   - result of the 21-point kronrod formula
!     reskh  - approximation to the mean value of f over (a,b),
!              i.e. to i/(b-a)
!
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none
  real(8) temp  !by Zunli for using f2py


  real ( kind = 8 ) a,absc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f
  dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)

  data wg  (  1) / 0.066671344308688137593568809893332d0 /
  data wg  (  2) / 0.149451349150580593145776339657697d0 /
  data wg  (  3) / 0.219086362515982043995534934228163d0 /
  data wg  (  4) / 0.269266719309996355091226921569469d0 /
  data wg  (  5) / 0.295524224714752870173892994651338d0 /

  data xgk (  1) / 0.995657163025808080735527280689003d0 /
  data xgk (  2) / 0.973906528517171720077964012084452d0 /
  data xgk (  3) / 0.930157491355708226001207180059508d0 /
  data xgk (  4) / 0.865063366688984510732096688423493d0 /
  data xgk (  5) / 0.780817726586416897063717578345042d0 /
  data xgk (  6) / 0.679409568299024406234327365114874d0 /
  data xgk (  7) / 0.562757134668604683339000099272694d0 /
  data xgk (  8) / 0.433395394129247190799265943165784d0 /
  data xgk (  9) / 0.294392862701460198131126603103866d0 /
  data xgk ( 10) / 0.148874338981631210884826001129720d0 /
  data xgk ( 11) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.011694638867371874278064396062192d0 /
  data wgk (  2) / 0.032558162307964727478818972459390d0 /
  data wgk (  3) / 0.054755896574351996031381300244580d0 /
  data wgk (  4) / 0.075039674810919952767043140916190d0 /
  data wgk (  5) / 0.093125454583697605535065465083366d0 /
  data wgk (  6) / 0.109387158802297641899210590325805d0 /
  data wgk (  7) / 0.123491976262065851077958109831074d0 /
  data wgk (  8) / 0.134709217311473325928054001771707d0 /
  data wgk (  9) / 0.142775938577060080797094273138717d0 /
  data wgk ( 10) / 0.147739104901338491374841515972068d0 /
  data wgk ( 11) / 0.149445554002916905664936468389821d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 21-point kronrod approximation to
!  the integral, and estimate the absolute error.
!
  resg = 0.0D+00
  fc = f(centr)
  resk = wgk(11)*fc
  resabs =  abs ( resk)
  do j=1,5
    jtw = 2*j
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j = 1,5
    jtwm1 = 2*j-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(11)* abs ( fc-reskh)

  do j=1,10
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) &
    abserr = resasc*min(0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)

  if(temp==1.0) temp=f(a)  !by Zunli for using f2py
  return
end subroutine dqk21

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dqk31(f,a,b,result,abserr,resabs,resasc)

!*****************************************************************************80
!
!! DQK31 carries out a 31 point Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b) with error
!                     estimate
!                 j = integral of abs(f) over (a,b)
!
!  Parameters:
!
!      on entry
!        f      - real ( kind = 8 )
!                 function subprogram defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the calling program.
!
!        a      - real ( kind = 8 )
!                 lower limit of integration
!
!        b      - real ( kind = 8 )
!                 upper limit of integration
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 31-point
!                 gauss-kronrod rule (resk), obtained by optimal
!                 addition of abscissae to the 15-point gauss
!                 rule (resg).
!
!        abserr - double precison
!                 estimate of the modulus of the modulus,
!                 which should not exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral j
!
!        resasc - real ( kind = 8 )
!                 approximation to the integral of abs(f-i/(b-a))
!                 over (a,b)
!
!  Local Parameters:
!
!
!     the abscissae and weights are given for the interval (-1,1).
!     because of symmetry only the positive abscissae and their
!     corresponding weights are given.
!
!     xgk    - abscissae of the 31-point kronrod rule
!              xgk(2), xgk(4), ...  abscissae of the 15-point
!              gauss rule
!              xgk(1), xgk(3), ...  abscissae which are optimally
!              added to the 15-point gauss rule
!
!     wgk    - weights of the 31-point kronrod rule
!
!     wg     - weights of the 15-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc   - abscissa
!     fval*  - function value
!     resg   - result of the 15-point gauss formula
!     resk   - result of the 31-point kronrod formula
!     reskh  - approximation to the mean value of f over (a,b),
!              i.e. to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none
  real(8) temp  !by Zunli for using f2py

  real ( kind = 8 ) a,absc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f

  dimension fv1(15),fv2(15),xgk(16),wgk(16),wg(8)

  data wg  (  1) / 0.030753241996117268354628393577204d0 /
  data wg  (  2) / 0.070366047488108124709267416450667d0 /
  data wg  (  3) / 0.107159220467171935011869546685869d0 /
  data wg  (  4) / 0.139570677926154314447804794511028d0 /
  data wg  (  5) / 0.166269205816993933553200860481209d0 /
  data wg  (  6) / 0.186161000015562211026800561866423d0 /
  data wg  (  7) / 0.198431485327111576456118326443839d0 /
  data wg  (  8) / 0.202578241925561272880620199967519d0 /

  data xgk (  1) / 0.998002298693397060285172840152271d0 /
  data xgk (  2) / 0.987992518020485428489565718586613d0 /
  data xgk (  3) / 0.967739075679139134257347978784337d0 /
  data xgk (  4) / 0.937273392400705904307758947710209d0 /
  data xgk (  5) / 0.897264532344081900882509656454496d0 /
  data xgk (  6) / 0.848206583410427216200648320774217d0 /
  data xgk (  7) / 0.790418501442465932967649294817947d0 /
  data xgk (  8) / 0.724417731360170047416186054613938d0 /
  data xgk (  9) / 0.650996741297416970533735895313275d0 /
  data xgk ( 10) / 0.570972172608538847537226737253911d0 /
  data xgk ( 11) / 0.485081863640239680693655740232351d0 /
  data xgk ( 12) / 0.394151347077563369897207370981045d0 /
  data xgk ( 13) / 0.299180007153168812166780024266389d0 /
  data xgk ( 14) / 0.201194093997434522300628303394596d0 /
  data xgk ( 15) / 0.101142066918717499027074231447392d0 /
  data xgk ( 16) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.005377479872923348987792051430128d0 /
  data wgk (  2) / 0.015007947329316122538374763075807d0 /
  data wgk (  3) / 0.025460847326715320186874001019653d0 /
  data wgk (  4) / 0.035346360791375846222037948478360d0 /
  data wgk (  5) / 0.044589751324764876608227299373280d0 /
  data wgk (  6) / 0.053481524690928087265343147239430d0 /
  data wgk (  7) / 0.062009567800670640285139230960803d0 /
  data wgk (  8) / 0.069854121318728258709520077099147d0 /
  data wgk (  9) / 0.076849680757720378894432777482659d0 /
  data wgk ( 10) / 0.083080502823133021038289247286104d0 /
  data wgk ( 11) / 0.088564443056211770647275443693774d0 /
  data wgk ( 12) / 0.093126598170825321225486872747346d0 /
  data wgk ( 13) / 0.096642726983623678505179907627589d0 /
  data wgk ( 14) / 0.099173598721791959332393173484603d0 /
  data wgk ( 15) / 0.100769845523875595044946662617570d0 /
  data wgk ( 16) / 0.101330007014791549017374792767493d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 31-point kronrod approximation to
!  the integral, and estimate the absolute error.
!
  fc = f(centr)
  resg = wg(8)*fc
  resk = wgk(16)*fc
  resabs =  abs ( resk)

  do j=1,7
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j = 1,8
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(16)* abs ( fc-reskh)

  do j=1,15
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) &
    abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)


  if(temp==1.0) temp=f(a)  !by Zunli for using f2py
  return
end subroutine dqk31

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dqk41 ( f, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! DQK41 carries out a 41 point Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b), with error
!                     estimate
!                 j = integral of abs(f) over (a,b)
!
!  Parameters:
!
!      on entry
!        f      - real ( kind = 8 )
!                 function subprogram defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the calling program.
!
!        a      - real ( kind = 8 )
!                 lower limit of integration
!
!        b      - real ( kind = 8 )
!                 upper limit of integration
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 41-point
!                 gauss-kronrod rule (resk) obtained by optimal
!                 addition of abscissae to the 20-point gauss
!                 rule (resg).
!
!        abserr - real ( kind = 8 )
!                 estimate of the modulus of the absolute error,
!                 which should not exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral j
!
!        resasc - real ( kind = 8 )
!                 approximation to the integal of abs(f-i/(b-a))
!                 over (a,b)
!
!  Local Parameters:
!
!
!     the abscissae and weights are given for the interval (-1,1).
!     because of symmetry only the positive abscissae and their
!     corresponding weights are given.
!
!     xgk    - abscissae of the 41-point gauss-kronrod rule
!              xgk(2), xgk(4), ...  abscissae of the 20-point
!              gauss rule
!              xgk(1), xgk(3), ...  abscissae which are optimally
!              added to the 20-point gauss rule
!
!     wgk    - weights of the 41-point gauss-kronrod rule
!
!     wg     - weights of the 20-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc   - abscissa
!     fval*  - function value
!     resg   - result of the 20-point gauss formula
!     resk   - result of the 41-point kronrod formula
!     reskh  - approximation to mean value of f over (a,b), i.e.
!              to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none
  real(8) temp  !by Zunli for using f2py

  
  real ( kind = 8 ) a,absc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f

  dimension fv1(20),fv2(20),xgk(21),wgk(21),wg(10)

  data wg  (  1) / 0.017614007139152118311861962351853d0 /
  data wg  (  2) / 0.040601429800386941331039952274932d0 /
  data wg  (  3) / 0.062672048334109063569506535187042d0 /
  data wg  (  4) / 0.083276741576704748724758143222046d0 /
  data wg  (  5) / 0.101930119817240435036750135480350d0 /
  data wg  (  6) / 0.118194531961518417312377377711382d0 /
  data wg  (  7) / 0.131688638449176626898494499748163d0 /
  data wg  (  8) / 0.142096109318382051329298325067165d0 /
  data wg  (  9) / 0.149172986472603746787828737001969d0 /
  data wg  ( 10) / 0.152753387130725850698084331955098d0 /

  data xgk (  1) / 0.998859031588277663838315576545863d0 /
  data xgk (  2) / 0.993128599185094924786122388471320d0 /
  data xgk (  3) / 0.981507877450250259193342994720217d0 /
  data xgk (  4) / 0.963971927277913791267666131197277d0 /
  data xgk (  5) / 0.940822633831754753519982722212443d0 /
  data xgk (  6) / 0.912234428251325905867752441203298d0 /
  data xgk (  7) / 0.878276811252281976077442995113078d0 /
  data xgk (  8) / 0.839116971822218823394529061701521d0 /
  data xgk (  9) / 0.795041428837551198350638833272788d0 /
  data xgk ( 10) / 0.746331906460150792614305070355642d0 /
  data xgk ( 11) / 0.693237656334751384805490711845932d0 /
  data xgk ( 12) / 0.636053680726515025452836696226286d0 /
  data xgk ( 13) / 0.575140446819710315342946036586425d0 /
  data xgk ( 14) / 0.510867001950827098004364050955251d0 /
  data xgk ( 15) / 0.443593175238725103199992213492640d0 /
  data xgk ( 16) / 0.373706088715419560672548177024927d0 /
  data xgk ( 17) / 0.301627868114913004320555356858592d0 /
  data xgk ( 18) / 0.227785851141645078080496195368575d0 /
  data xgk ( 19) / 0.152605465240922675505220241022678d0 /
  data xgk ( 20) / 0.076526521133497333754640409398838d0 /
  data xgk ( 21) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.003073583718520531501218293246031d0 /
  data wgk (  2) / 0.008600269855642942198661787950102d0 /
  data wgk (  3) / 0.014626169256971252983787960308868d0 /
  data wgk (  4) / 0.020388373461266523598010231432755d0 /
  data wgk (  5) / 0.025882133604951158834505067096153d0 /
  data wgk (  6) / 0.031287306777032798958543119323801d0 /
  data wgk (  7) / 0.036600169758200798030557240707211d0 /
  data wgk (  8) / 0.041668873327973686263788305936895d0 /
  data wgk (  9) / 0.046434821867497674720231880926108d0 /
  data wgk ( 10) / 0.050944573923728691932707670050345d0 /
  data wgk ( 11) / 0.055195105348285994744832372419777d0 /
  data wgk ( 12) / 0.059111400880639572374967220648594d0 /
  data wgk ( 13) / 0.062653237554781168025870122174255d0 /
  data wgk ( 14) / 0.065834597133618422111563556969398d0 /
  data wgk ( 15) / 0.068648672928521619345623411885368d0 /
  data wgk ( 16) / 0.071054423553444068305790361723210d0 /
  data wgk ( 17) / 0.073030690332786667495189417658913d0 /
  data wgk ( 18) / 0.074582875400499188986581418362488d0 /
  data wgk ( 19) / 0.075704497684556674659542775376617d0 /
  data wgk ( 20) / 0.076377867672080736705502835038061d0 /
  data wgk ( 21) / 0.076600711917999656445049901530102d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 41-point gauss-kronrod approximation to
!  the integral, and estimate the absolute error.
!
  resg = 0.0D+00
  fc = f(centr)
  resk = wgk(21)*fc
  resabs =  abs ( resk)

  do j=1,10
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j = 1,10
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(21)* abs ( fc-reskh)

  do j=1,20
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.D+00) &
    abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)


  if(temp==1.0) temp=f(a)  !by Zunli for using f2py
  return
end subroutine dqk41

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dqk51(f,a,b,result,abserr,resabs,resasc)

!*****************************************************************************80
!
!! DQK51 carries out a 51 point Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b) with error
!                     estimate
!                 j = integral of abs(f) over (a,b)
!
!  Parameters:
!
!      on entry
!        f      - real ( kind = 8 )
!                 function defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the calling program.
!
!        a      - real ( kind = 8 )
!                 lower limit of integration
!
!        b      - real ( kind = 8 )
!                 upper limit of integration
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 51-point
!                 kronrod rule (resk) obtained by optimal addition
!                 of abscissae to the 25-point gauss rule (resg).
!
!        abserr - real ( kind = 8 )
!                 estimate of the modulus of the absolute error,
!                 which should not exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral j
!
!        resasc - real ( kind = 8 )
!                 approximation to the integral of abs(f-i/(b-a))
!                 over (a,b)
!
!  Local Parameters:
!
!     the abscissae and weights are given for the interval (-1,1).
!     because of symmetry only the positive abscissae and their
!     corresponding weights are given.
!
!     xgk    - abscissae of the 51-point kronrod rule
!              xgk(2), xgk(4), ...  abscissae of the 25-point
!              gauss rule
!              xgk(1), xgk(3), ...  abscissae which are optimally
!              added to the 25-point gauss rule
!
!     wgk    - weights of the 51-point kronrod rule
!
!     wg     - weights of the 25-point gauss rule
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc   - abscissa
!     fval*  - function value
!     resg   - result of the 25-point gauss formula
!     resk   - result of the 51-point kronrod formula
!     reskh  - approximation to the mean value of f over (a,b),
!              i.e. to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none
  real(8) temp  !by Zunli for using f2py

  real ( kind = 8 ) a,absc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f

  dimension fv1(25),fv2(25),xgk(26),wgk(26),wg(13)

  data wg  (  1) / 0.011393798501026287947902964113235d0 /
  data wg  (  2) / 0.026354986615032137261901815295299d0 /
  data wg  (  3) / 0.040939156701306312655623487711646d0 /
  data wg  (  4) / 0.054904695975835191925936891540473d0 /
  data wg  (  5) / 0.068038333812356917207187185656708d0 /
  data wg  (  6) / 0.080140700335001018013234959669111d0 /
  data wg  (  7) / 0.091028261982963649811497220702892d0 /
  data wg  (  8) / 0.100535949067050644202206890392686d0 /
  data wg  (  9) / 0.108519624474263653116093957050117d0 /
  data wg  ( 10) / 0.114858259145711648339325545869556d0 /
  data wg  ( 11) / 0.119455763535784772228178126512901d0 /
  data wg  ( 12) / 0.122242442990310041688959518945852d0 /
  data wg  ( 13) / 0.123176053726715451203902873079050d0 /

  data xgk (  1) / 0.999262104992609834193457486540341d0 /
  data xgk (  2) / 0.995556969790498097908784946893902d0 /
  data xgk (  3) / 0.988035794534077247637331014577406d0 /
  data xgk (  4) / 0.976663921459517511498315386479594d0 /
  data xgk (  5) / 0.961614986425842512418130033660167d0 /
  data xgk (  6) / 0.942974571228974339414011169658471d0 /
  data xgk (  7) / 0.920747115281701561746346084546331d0 /
  data xgk (  8) / 0.894991997878275368851042006782805d0 /
  data xgk (  9) / 0.865847065293275595448996969588340d0 /
  data xgk ( 10) / 0.833442628760834001421021108693570d0 /
  data xgk ( 11) / 0.797873797998500059410410904994307d0 /
  data xgk ( 12) / 0.759259263037357630577282865204361d0 /
  data xgk ( 13) / 0.717766406813084388186654079773298d0 /
  data xgk ( 14) / 0.673566368473468364485120633247622d0 /
  data xgk ( 15) / 0.626810099010317412788122681624518d0 /
  data xgk ( 16) / 0.577662930241222967723689841612654d0 /
  data xgk ( 17) / 0.526325284334719182599623778158010d0 /
  data xgk ( 18) / 0.473002731445714960522182115009192d0 /
  data xgk ( 19) / 0.417885382193037748851814394594572d0 /
  data xgk ( 20) / 0.361172305809387837735821730127641d0 /
  data xgk ( 21) / 0.303089538931107830167478909980339d0 /
  data xgk ( 22) / 0.243866883720988432045190362797452d0 /
  data xgk ( 23) / 0.183718939421048892015969888759528d0 /
  data xgk ( 24) / 0.122864692610710396387359818808037d0 /
  data xgk ( 25) / 0.061544483005685078886546392366797d0 /
  data xgk ( 26) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.001987383892330315926507851882843d0 /
  data wgk (  2) / 0.005561932135356713758040236901066d0 /
  data wgk (  3) / 0.009473973386174151607207710523655d0 /
  data wgk (  4) / 0.013236229195571674813656405846976d0 /
  data wgk (  5) / 0.016847817709128298231516667536336d0 /
  data wgk (  6) / 0.020435371145882835456568292235939d0 /
  data wgk (  7) / 0.024009945606953216220092489164881d0 /
  data wgk (  8) / 0.027475317587851737802948455517811d0 /
  data wgk (  9) / 0.030792300167387488891109020215229d0 /
  data wgk ( 10) / 0.034002130274329337836748795229551d0 /
  data wgk ( 11) / 0.037116271483415543560330625367620d0 /
  data wgk ( 12) / 0.040083825504032382074839284467076d0 /
  data wgk ( 13) / 0.042872845020170049476895792439495d0 /
  data wgk ( 14) / 0.045502913049921788909870584752660d0 /
  data wgk ( 15) / 0.047982537138836713906392255756915d0 /
  data wgk ( 16) / 0.050277679080715671963325259433440d0 /
  data wgk ( 17) / 0.052362885806407475864366712137873d0 /
  data wgk ( 18) / 0.054251129888545490144543370459876d0 /
  data wgk ( 19) / 0.055950811220412317308240686382747d0 /
  data wgk ( 20) / 0.057437116361567832853582693939506d0 /
  data wgk ( 21) / 0.058689680022394207961974175856788d0 /
  data wgk ( 22) / 0.059720340324174059979099291932562d0 /
  data wgk ( 23) / 0.060539455376045862945360267517565d0 /
  data wgk ( 24) / 0.061128509717053048305859030416293d0 /
  data wgk ( 25) / 0.061471189871425316661544131965264d0 /
  data wgk ( 26) / 0.061580818067832935078759824240066d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 51-point kronrod approximation to
!  the integral, and estimate the absolute error.
!
  fc = f(centr)
  resg = wg(13)*fc
  resk = wgk(26)*fc
  resabs =  abs ( resk)

  do j=1,12
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j = 1,13
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(26)* abs ( fc-reskh)

  do j=1,25
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) &
    abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)


  if(temp==1.0) temp=f(a)  !by Zunli for using f2py
  return
end subroutine dqk51

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine dqk61(f,a,b,result,abserr,resabs,resasc)

!*****************************************************************************80
!
!! DQK61 carries out a 61 point Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b) with error
!                     estimate
!                 j = integral of  abs ( f) over (a,b)
!
!  Parameters:
!
!   on entry
!     f      - real ( kind = 8 )
!              function subprogram defining the integrand
!              function f(x). the actual name for f needs to be
!              declared e x t e r n a l in the calling program.
!
!     a      - real ( kind = 8 )
!              lower limit of integration
!
!     b      - real ( kind = 8 )
!              upper limit of integration
!
!   on return
!     result - real ( kind = 8 )
!              approximation to the integral i
!              result is computed by applying the 61-point
!              kronrod rule (resk) obtained by optimal addition of
!              abscissae to the 30-point gauss rule (resg).
!
!     abserr - real ( kind = 8 )
!              estimate of the modulus of the absolute error,
!              which should equal or exceed  abs ( i-result)
!
!     resabs - real ( kind = 8 )
!              approximation to the integral j
!
!     resasc - real ( kind = 8 )
!              approximation to the integral of  abs ( f-i/(b-a))
!
!  Local Parameters:
!
!     the abscissae and weights are given for the
!     interval (-1,1). because of symmetry only the positive
!     abscissae and their corresponding weights are given.
!
!     xgk   - abscissae of the 61-point kronrod rule
!             xgk(2), xgk(4)  ... abscissae of the 30-point
!             gauss rule
!             xgk(1), xgk(3)  ... optimally added abscissae
!             to the 30-point gauss rule
!
!     wgk   - weights of the 61-point kronrod rule
!
!     wg    - weigths of the 30-point gauss rule
!
!
!   gauss quadrature weights and kronron quadrature abscissae and weights
!   as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
!   bell labs, nov. 1981.
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     dabsc  - abscissa
!     fval*  - function value
!     resg   - result of the 30-point gauss rule
!     resk   - result of the 61-point kronrod rule
!     reskh  - approximation to the mean value of f
!              over (a,b), i.e. to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none
  real(8) temp  !by Zunli for using f2py

  
  real ( kind = 8 ) a,dabsc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f

  dimension fv1(30),fv2(30),xgk(31),wgk(31),wg(15)

  data wg  (  1) / 0.007968192496166605615465883474674d0 /
  data wg  (  2) / 0.018466468311090959142302131912047d0 /
  data wg  (  3) / 0.028784707883323369349719179611292d0 /
  data wg  (  4) / 0.038799192569627049596801936446348d0 /
  data wg  (  5) / 0.048402672830594052902938140422808d0 /
  data wg  (  6) / 0.057493156217619066481721689402056d0 /
  data wg  (  7) / 0.065974229882180495128128515115962d0 /
  data wg  (  8) / 0.073755974737705206268243850022191d0 /
  data wg  (  9) / 0.080755895229420215354694938460530d0 /
  data wg  ( 10) / 0.086899787201082979802387530715126d0 /
  data wg  ( 11) / 0.092122522237786128717632707087619d0 /
  data wg  ( 12) / 0.096368737174644259639468626351810d0 /
  data wg  ( 13) / 0.099593420586795267062780282103569d0 /
  data wg  ( 14) / 0.101762389748405504596428952168554d0 /
  data wg  ( 15) / 0.102852652893558840341285636705415d0 /

  data xgk (  1) / 0.999484410050490637571325895705811d0 /
  data xgk (  2) / 0.996893484074649540271630050918695d0 /
  data xgk (  3) / 0.991630996870404594858628366109486d0 /
  data xgk (  4) / 0.983668123279747209970032581605663d0 /
  data xgk (  5) / 0.973116322501126268374693868423707d0 /
  data xgk (  6) / 0.960021864968307512216871025581798d0 /
  data xgk (  7) / 0.944374444748559979415831324037439d0 /
  data xgk (  8) / 0.926200047429274325879324277080474d0 /
  data xgk (  9) / 0.905573307699907798546522558925958d0 /
  data xgk ( 10) / 0.882560535792052681543116462530226d0 /
  data xgk ( 11) / 0.857205233546061098958658510658944d0 /
  data xgk ( 12) / 0.829565762382768397442898119732502d0 /
  data xgk ( 13) / 0.799727835821839083013668942322683d0 /
  data xgk ( 14) / 0.767777432104826194917977340974503d0 /
  data xgk ( 15) / 0.733790062453226804726171131369528d0 /
  data xgk ( 16) / 0.697850494793315796932292388026640d0 /
  data xgk ( 17) / 0.660061064126626961370053668149271d0 /
  data xgk ( 18) / 0.620526182989242861140477556431189d0 /
  data xgk ( 19) / 0.579345235826361691756024932172540d0 /
  data xgk ( 20) / 0.536624148142019899264169793311073d0 /
  data xgk ( 21) / 0.492480467861778574993693061207709d0 /
  data xgk ( 22) / 0.447033769538089176780609900322854d0 /
  data xgk ( 23) / 0.400401254830394392535476211542661d0 /
  data xgk ( 24) / 0.352704725530878113471037207089374d0 /
  data xgk ( 25) / 0.304073202273625077372677107199257d0 /
  data xgk ( 26) / 0.254636926167889846439805129817805d0 /
  data xgk ( 27) / 0.204525116682309891438957671002025d0 /
  data xgk ( 28) / 0.153869913608583546963794672743256d0 /
  data xgk ( 29) / 0.102806937966737030147096751318001d0 /
  data xgk ( 30) / 0.051471842555317695833025213166723d0 /
  data xgk ( 31) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.001389013698677007624551591226760d0 /
  data wgk (  2) / 0.003890461127099884051267201844516d0 /
  data wgk (  3) / 0.006630703915931292173319826369750d0 /
  data wgk (  4) / 0.009273279659517763428441146892024d0 /
  data wgk (  5) / 0.011823015253496341742232898853251d0 /
  data wgk (  6) / 0.014369729507045804812451432443580d0 /
  data wgk (  7) / 0.016920889189053272627572289420322d0 /
  data wgk (  8) / 0.019414141193942381173408951050128d0 /
  data wgk (  9) / 0.021828035821609192297167485738339d0 /
  data wgk ( 10) / 0.024191162078080601365686370725232d0 /
  data wgk ( 11) / 0.026509954882333101610601709335075d0 /
  data wgk ( 12) / 0.028754048765041292843978785354334d0 /
  data wgk ( 13) / 0.030907257562387762472884252943092d0 /
  data wgk ( 14) / 0.032981447057483726031814191016854d0 /
  data wgk ( 15) / 0.034979338028060024137499670731468d0 /
  data wgk ( 16) / 0.036882364651821229223911065617136d0 /
  data wgk ( 17) / 0.038678945624727592950348651532281d0 /
  data wgk ( 18) / 0.040374538951535959111995279752468d0 /
  data wgk ( 19) / 0.041969810215164246147147541285970d0 /
  data wgk ( 20) / 0.043452539701356069316831728117073d0 /
  data wgk ( 21) / 0.044814800133162663192355551616723d0 /
  data wgk ( 22) / 0.046059238271006988116271735559374d0 /
  data wgk ( 23) / 0.047185546569299153945261478181099d0 /
  data wgk ( 24) / 0.048185861757087129140779492298305d0 /
  data wgk ( 25) / 0.049055434555029778887528165367238d0 /
  data wgk ( 26) / 0.049795683427074206357811569379942d0 /
  data wgk ( 27) / 0.050405921402782346840893085653585d0 /
  data wgk ( 28) / 0.050881795898749606492297473049805d0 /
  data wgk ( 29) / 0.051221547849258772170656282604944d0 /
  data wgk ( 30) / 0.051426128537459025933862879215781d0 /
  data wgk ( 31) / 0.051494729429451567558340433647099d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(b+a)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 61-point kronrod approximation to the
!  integral, and estimate the absolute error.
!
  resg = 0.0D+00
  fc = f(centr)
  resk = wgk(31)*fc
  resabs =  abs ( resk)

  do j=1,15
    jtw = j*2
    dabsc = hlgth*xgk(jtw)
    fval1 = f(centr-dabsc)
    fval2 = f(centr+dabsc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j=1,15
    jtwm1 = j*2-1
    dabsc = hlgth*xgk(jtwm1)
    fval1 = f(centr-dabsc)
    fval2 = f(centr+dabsc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(31)* abs ( fc-reskh)

  do j=1,30
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) &
    abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)


  if(temp==1.0) temp=f(a)  !by Zunli for using f2py
  return
end subroutine dqk61

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dqpsrt ( limit, last, maxerr, ermax, elist, iord, nrmax )

!*****************************************************************************80
!
!! DQPSRT maintains the order of a list of local error estimates.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  this routine maintains the descending ordering in the
!      list of the local error estimated resulting from the
!      interval subdivision process. at each call two error
!      estimates are inserted using the sequential search
!      method, top-down for the largest error estimate and
!      bottom-up for the smallest error estimate.
!
!  Parameters:
!
!        limit  - integer ( kind = 4 )
!                 maximum number of error estimates the list
!                 can contain
!
!        last   - integer ( kind = 4 )
!                 number of error estimates currently in the list
!
!        maxerr - integer ( kind = 4 )
!                 maxerr points to the nrmax-th largest error
!                 estimate currently in the list
!
!        ermax  - real ( kind = 8 )
!                 nrmax-th largest error estimate
!                 ermax = elist(maxerr)
!
!        elist  - real ( kind = 8 )
!                 vector of dimension last containing
!                 the error estimates
!
!        iord   - integer ( kind = 4 )
!                 vector of dimension last, the first k elements
!                 of which contain pointers to the error
!                 estimates, such that
!                 elist(iord(1)),...,  elist(iord(k))
!                 form a decreasing sequence, with
!                 k = last if last.le.(limit/2+2), and
!                 k = limit+1-last otherwise
!
!        nrmax  - integer ( kind = 4 )
!                 maxerr = iord(nrmax)
!
  implicit none

  real ( kind = 8 ) elist,ermax,errmax,errmin
  integer ( kind = 4 ) i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last, &
    lim
  integer ( kind = 4 ) limit
  integer ( kind = 4 ) maxerr
  integer ( kind = 4 ) nrmax
  dimension elist(last),iord(last)
!
!  check whether the list contains more than
!  two error estimates.
!
  if(last.gt.2) go to 10
  iord(1) = 1
  iord(2) = 2
  go to 90
!
!  this part of the routine is only executed if, due to a
!  difficult integrand, subdivision increased the error
!  estimate. in the normal case the insert procedure should
!  start after the nrmax-th largest error estimate.
!
   10 errmax = elist(maxerr)

  ido = nrmax-1
  do i = 1,ido
    isucc = iord(nrmax-1)
    if(errmax.le.elist(isucc)) go to 30
    iord(nrmax) = isucc
    nrmax = nrmax-1
  end do
!
!  compute the number of elements in the list to be maintained
!  in descending order. this number depends on the number of
!  subdivisions still allowed.
!
   30 jupbn = last
  if(last.gt.(limit/2+2)) jupbn = limit+3-last
  errmin = elist(last)
!
!  insert errmax by traversing the list top-down,
!  starting comparison from the element elist(iord(nrmax+1)).
!
  jbnd = jupbn-1
  ibeg = nrmax+1

  do i=ibeg,jbnd
    isucc = iord(i)
    if(errmax.ge.elist(isucc)) go to 60
    iord(i-1) = isucc
  end do

  iord(jbnd) = maxerr
  iord(jupbn) = last
  go to 90
!
!  insert errmin by traversing the list bottom-up.
!
   60 iord(i-1) = maxerr
  k = jbnd

  do j=i,jbnd
    isucc = iord(k)
    if(errmin.lt.elist(isucc)) go to 80
    iord(k+1) = isucc
    k = k-1
  end do

  iord(i) = last
  go to 90
   80 iord(k+1) = last
!
!     set maxerr and ermax.
!
   90 maxerr = iord(nrmax)
  ermax = elist(maxerr)

  return
end subroutine dqpsrt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine xerror ( xmess, nmess, nerr, level )

!*****************************************************************************80
!
!! XERROR replaces the SLATEC XERROR routine.
!
!  Modified:
!
!    12 September 2015
!
  implicit none

  integer ( kind = 4 ) level
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) nmess
  character ( len = * ) xmess

  if ( 1 <= LEVEL ) then
    WRITE ( *,'(1X,A)') XMESS(1:NMESS)
    WRITE ( *,'('' ERROR NUMBER = '',I5,'', MESSAGE LEVEL = '',I5)') &
        NERR,LEVEL
  end if

  return
end subroutine xerror

end module dqag_and_dqags


!%%%%%%%%%%%%%%%%%%%%%%%% module myquad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!************************               ******************************
! module myquad, estimates a double definite integral.
! Author: Zunli Yuan
! !!!! version 2019_11_16_16_01

!   main subroutine : myquad2d
!
!   Parameters:
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      g1     - real ( kind = 8 )
!               function subprogram defining the lower limit 
!               function g1(x). the actual name for g1 needs to be
!               declared 'external' in the driver program.
!      g2     - real ( kind = 8 )
!               function subprogram defining the upper limit
!               function g2(x). the actual name for g2 needs to be
!               declared 'external' in the driver program.

!      f2d    - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f2d(x,y). the actual name for f2d needs to be
!               declared 'external' in the driver program.
module myquad  
  use dqag_and_dqags
  implicit none
  procedure(f1),pointer :: p1
  procedure(f2),pointer :: p2
  procedure(f3),pointer :: p3  
  
  interface    
    function f1(x)
      implicit none
      real(8),intent(in) :: x
      real(8) f1
    end function
    
    function f2(x)
      implicit none
      real(8),intent(in) :: x
      real(8) f2
    end function

    function f3(x,y)
      implicit none
      real(8),intent(in) :: x,y
      real(8) f3
    end function
  end interface  
 
  contains
  
  subroutine myquad2d(a,b,g1,g2,f2d,ans)
    implicit none
    !real (8), external :: g1,g2,f2d  !Declaring variables in this way is not supported here by f2py
    real(8) temp
    real(8) g1,g2,f2d
    external g1,g2,f2d  !should be this for running f2py    
    real(8),intent(in) :: a,b
    real(8),intent(out) :: ans
    integer ( kind = 4 ), parameter :: limit = 500
    real ( kind = 8 ) abserr
    real ( kind = 8 ), parameter :: epsabs = 0.0D+00
    real ( kind = 8 ), parameter :: epsrel = 0.001D+00
    integer ( kind = 4 ) ier
    integer ( kind = 4 ) last
    integer ( kind = 4 ) neval
    !%%%%%%%%%%%%%%  same  
    integer ( kind = 4 ), parameter :: lenw = 4 * limit
    integer ( kind = 4 ) iwork(limit)   
    real ( kind = 8 ) work(lenw)
    integer ( kind = 4 ), parameter :: key = 6
    !%%%%%%%%%%%%%%  different
    
    if(temp==1.0) then  ! for using f2py
      temp=g1(a)        ! redundant code, but it seems to be necessary for using f2py 
      temp=g2(a)
      temp=f2d(a,b)
    end if              ! for using f2py
    
    p1 => g1
    p2 => g2
    p3 => f2d 
    
    call dqag ( H, a, b, epsabs, epsrel, key, ans, abserr, neval, ier, &
      limit, lenw, last, iwork, work )  
  
  end subroutine myquad2d 
  
  real(8) function H(x)
    implicit none
    !real ( kind = 8 ), external :: f
    real(8) a,b,ans
    integer ( kind = 4 ), parameter :: limit = 500
    real ( kind = 8 ) abserr
    real ( kind = 8 ), parameter :: epsabs = 0.0D+00
    real ( kind = 8 ), parameter :: epsrel = 0.001D+00
    integer ( kind = 4 ) ier
    integer ( kind = 4 ) last
    integer ( kind = 4 ) neval
    !%%%%%%%%%%%%%%  same  
    integer ( kind = 4 ), parameter :: lenw = 4 * limit
    integer ( kind = 4 ) iwork(limit)   
    real ( kind = 8 ) work(lenw)
    integer ( kind = 4 ), parameter :: key = 6
    !%%%%%%%%%%%%%%  different
    real(8) x,xx
    common /groupx/ xx
    
    xx=x
    a=p1(x)
    b=p2(x)
    call dqag ( f, a, b, epsabs, epsrel, key, ans, abserr, neval, ier, &
      limit, lenw, last, iwork, work )
    H=ans
  end function
  
  real(8) function f(y)
    implicit none
    real(8),external :: p
    real(8) x,y
    common /groupx/ x 
    
    f=p3(x,y)  
  end function

end module myquad
!***********************                 ***********************************
!                       end module myquad
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  
  
  
  
  
  
  

























!include "User_defined_func.f90"
!***********************************************************************************************************


module params
 implicit none
  real(8),parameter :: pi=3.141592653589793
  real(8), allocatable, dimension(:) :: xz,red,ylum,hi,f_pilot,weight,ni 
  real(8),save :: H0,matter,lambda
  real(8),save :: z1,z2,Lmin,Lmax
  real(8),save :: L1,L2
  integer,save :: ndata,process
  real(8),save :: h1,h2
  real(8),save :: h10,h20,beta 
  logical,save :: set_beta_fixed,small_sample_approximation,adaptive,weighting,absolute_magnitude,unbounded_lum
  logical, allocatable :: xdif(:,:),ydif(:,:)
  real(8), allocatable :: limx(:),limy(:)
  real(8),save :: nw
  !!$omp THREADPRIVATE(L1,L2)
end module



MODULE prob_functions
  implicit none
  CONTAINS  
!*****************************************
real(8) function p2d(z,L)
  use params
  implicit none
  real(8) x,y,z,L,bandwidth
  real(8),external :: f_lim,f_ref

  x=log( (z-z1)/(z2-z) )
  y=L-f_lim(z)
  p2d=f_ref(x,y) * ( 1/(z-z1) + 1/(z2-z) )
  return 
end function
!*****************************************
real(8) function p2da(z,L)
  use params
  implicit none
  real(8) x,y,z,L
  real(8),external :: f_lim,f_ada

  x=log( (z-z1)/(z2-z) )
  y=L-f_lim(z)
  p2da=f_ada(x,y) * ( 1/(z-z1) + 1/(z2-z) )
  return
end function
!*****************************************
real(8) function p1d(z,L)
  use params
  implicit none
  real(8) y,z,L
  real(8),external :: f_lim,f1d

  y=L-f_lim(z)
  p1d=f1d(y)/(z2-z1)
end function
!*****************************************
real(8) function p1da(z,L)
  use params
  implicit none
  real(8) y,z,L
  real(8),external :: f_lim,f1da
  
  y=L-f_lim(z)
  p1da=f1da(y)/(z2-z1)
  return   
end function
!*****************************************
real(8) function pw2d(z,L)
  use params
  implicit none
  real(8) x,y,z,L
  real(8),external :: f_lim,fw_ref

  x=log( (z-z1)/(z2-z))
  y=L-f_lim(z)
  pw2d=fw_ref(x,y) * ( 1/(z-z1) + 1/(z2-z) )
  return 
end function
!*****************************************
real(8) function pw2da(z,L)
  use params
  implicit none
  real(8) x,y,z,L
  real(8),external :: f_lim,fw_ada

  x=log( (z-z1)/(z2-z) )
  y=L-f_lim(z)
  pw2da=fw_ada(x,y) * ( 1/(z-z1) + 1/(z2-z) )
  if ( (pw2da>10000.0) .or. (pw2da<0.0) ) then
    write(110,"(f20.16,2X,3f9.4,2X,f12.3,2X,2f8.4)") z,L,x,y,pw2da,h10,h20    
  end if  
  return
end function

END MODULE prob_functions

!################################################################################################

MODULE fhi_functions           ! module for leave-more-out functions  
  implicit none
  CONTAINS
 !***************************************** 
real(8) function fhi(i,xi,yi)
  use params
  implicit none
  real(8) xi,yi
  real(8) x_m,y_m,y_p,xj,yj,temp
  real(8),external :: K
  integer i,j
  
  fhi=0.0
  do j=1,ndata
    xj=xz(j)
    yj=ylum(j)            	    
    x_m=(xi-xj)/h1            
    y_m=(yi-yj)/h2
    y_p=(yi+yj)/h2

    if ( xdif(j,i) ) then
      if ( ydif(j,i) ) then  
        temp=( K(x_m,y_m) + K(x_m,y_p) ) 
      else
        temp=K(x_m,y_p)
      end if        
      fhi=fhi + temp
    end if 
  end do
  fhi=fhi/(h1*h2) 
end function
!*****************************************

real(8) function fhi_ada(i,xi,yi)
  use params
  implicit none
  real(8) xi,yi,ha1,ha2
  real(8) x_m,y_m,y_p,xj,yj,temp
  real(8),external :: K
  integer i,j  
  
  fhi_ada=0.0
  do j=1,ndata
    xj=xz(j)
    yj=ylum(j)            	    
    ha1=h10*hi(j)
    ha2=h20*hi(j)
    x_m=(xi-xj)/ha1            
    y_m=(yi-yj)/ha2
    y_p=(yi+yj)/ha2

    if ( xdif(j,i) ) then
      if ( ydif(j,i) ) then  
        temp=( K(x_m,y_m) + K(x_m,y_p) ) / (ha1*ha2)  
      else
        temp=K(x_m,y_p) / (ha1*ha2)
      end if        
      fhi_ada=fhi_ada + temp
    end if
    
  end do   
  return
end function
!*****************************************
real(8) function fwhi(i,xi,yi)
  use params
  implicit none
  real(8) xi,yi
  real(8) x_m,y_m,y_p,xj,yj,temp
  real(8),external :: K
  integer i,j
  
  fwhi=0.0
  do j=1,ndata
    xj=xz(j)
    yj=ylum(j)            	    
    x_m=(xi-xj)/h1            
    y_m=(yi-yj)/h2
    y_p=(yi+yj)/h2
    if ( xdif(j,i) ) then
      if ( ydif(j,i) ) then  
        temp=( K(x_m,y_m) + K(x_m,y_p) ) 
      else
        temp=K(x_m,y_p)
      end if        
      fwhi=fwhi + temp*weight(j)
    end if          
  end do
  fwhi=fwhi/(h1*h2) 
end function
!*****************************************

real(8) function fwhi_ada(i,xi,yi)
  use params
  implicit none
  real(8) xi,yi,ha1,ha2
  real(8) x_m,y_m,y_p,xj,yj,temp
  real(8),external :: K
  integer i,j  
  
  fwhi_ada=0.0
  do j=1,ndata
    xj=xz(j)
    yj=ylum(j)            	    
    ha1=h10*hi(j)
    ha2=h20*hi(j)
    x_m=(xi-xj)/ha1            
    y_m=(yi-yj)/ha2
    y_p=(yi+yj)/ha2

    if ( xdif(j,i) ) then
      if ( ydif(j,i) ) then  
        temp=( K(x_m,y_m) + K(x_m,y_p) ) / (ha1*ha2)  
      else
        temp=K(x_m,y_p) / (ha1*ha2)
      end if        
      fwhi_ada=fwhi_ada + temp*weight(j)
    end if
  end do   
  return
end function

END MODULE fhi_functions

!################################################################################################
!               
!************************************************************************************************


subroutine f2py_value(nd,input_python,md,process_python)
  use params
  implicit none
  integer :: nd,md,process_python
  real(8) :: input_python(md)
!f2py intent(in) :: nd,input_python,process_python 
!integer intent(hide),depend(input_python) :: md=shape(input_python)   
  ndata=nd 
  process=process_python
  H0=input_python(1)
  matter=input_python(2)
  lambda=1-matter
  z1=input_python(3)
  z2=input_python(4)
  Lmin=input_python(5)
  Lmax=input_python(6)     
  return
end subroutine f2py_value

!*****************************************

subroutine initialize(epson,ans)
  use params
  implicit none
  integer i,j,num
  real(8) epson,ans,dif
  !f2py intent(in) :: epson
  !f2py intent(out) :: ans 
  !real(8),parameter :: epson=1e-9
  
  if (allocated(ni)) deallocate(ni)   
  if (allocated(xdif)) deallocate(xdif)
  if (allocated(ydif)) deallocate(ydif)  
  
  allocate(xdif(ndata,ndata))
  allocate(ydif(ndata,ndata))  
  open(unit=101,file="temp1.txt")
  open(unit=102,file="temp2.txt")  
  open(unit=103,file="temp3.txt")
  num=0
  do i=1,ndata
    do j=1,ndata
      if ( abs(xz(i)-xz(j))<epson ) then
        xdif(j,i)=.false.
        if ( (i /= j) .and.  ( abs(xz(i)-xz(j))>1e-9 ) ) then
          write(101,*) 'x,j,i',j,i, xz(i),xz(j)
        end if        
      else
        xdif(j,i)=.true.
      end if  
      if ( abs(ylum(i)-ylum(j))<epson ) then
        ydif(j,i)=.false.
        if ( (i /= j) .and. ( abs(ylum(i)-ylum(j))>1e-9 ) )then
          write(102,*) 'y,j,i',j,i, ylum(i),ylum(j)
        end if        
      else
        ydif(j,i)=.true.
      end if                
      
      dif=abs(ylum(i)-ylum(j))
      if ( dif > 1e-50 .and. dif < epson ) then
        num=num+1
      end if      
      !write(101,*) 'j,i',j,i,xdif(j,i),ydif(j,i)
    end do 
  end do 
   
  ans=num*1.0/ndata/ndata 
  !print*,num,ans 
   
  allocate(ni(ndata))
  if (weighting) then    
    do i=1,ndata
      ni(i)=0.0
      do j=1,ndata
        if ( xdif(j,i)==.false. ) then
          ni(i)=ni(i) + 2*weight(j)
        end if      
        if ( xdif(j,i) .and. ydif(j,i)==.false. )  then
          ni(i)= ni(i) + weight(j)     
        end if  
      end do
      write(103,*) 'i',i,ni(i)
    end do 
  end if 
  
  if (.not. weighting) then    
    do i=1,ndata
      ni(i)=0.0
      do j=1,ndata
        if ( xdif(j,i)==.false. ) then
          ni(i)=ni(i) + 2.0
        end if      
        if ( xdif(j,i) .and. ydif(j,i)==.false. )  then
          ni(i)= ni(i) + 1.0     
        end if  
      end do
      write(103,*) 'i',i,ni(i)
    end do 
  end if 
  
  
  !print*,ni(1)
  !print*,'nw',nw

  return    
end subroutine initialize


subroutine check()
  use params

  if (allocated(xz)) deallocate(xz)  
  if (allocated(red)) deallocate(red)
  if (allocated(ylum)) deallocate(ylum)  
  if (allocated(hi)) deallocate(hi)
  if (allocated(f_pilot)) deallocate(f_pilot)
  if (allocated(weight)) deallocate(weight) 
  if (allocated(ni)) deallocate(ni)   
  if (allocated(xdif)) deallocate(xdif)
  if (allocated(ydif)) deallocate(ydif)
  if (allocated(limx)) deallocate(limx)  
  if (allocated(limy)) deallocate(limy)    
  return    
end subroutine check

!################################################################################################
!************************************************************************************************
subroutine lnlike(h1_py,h2_py,ans)
  use params
  use myquad
  use prob_functions
  use fhi_functions
  implicit none
  real(8) :: h1_py,h2_py,ans
!f2py intent(in) :: h1_py,h2_py
!f2py intent(out) :: ans 
  real(8) :: f_hi(ndata)
  real(8) temp,x_m,y_m,y_p
  real(8) xi,yi,ftemp,xj,yj
  real(8) int2d,log_lik
  integer i,j
  !real(8),external :: fhi,fwhi
  real(8),external :: X1,X2   !p,pw    
  procedure(p2d), pointer :: p
  procedure(fhi), pointer :: fi
    
  h1=h1_py
  h2=h2_py 
  ftemp=0.0
  log_lik=0.0

  if (weighting) then
    p => pw2d
    fi => fwhi
  else
    p => p2d
    fi => fhi
    nw=ndata*1.0
  end if

  !!$omp parallel NUM_THREADS(process) private(i,xi,yi,ftemp) shared(log_lik)  	
  !!$omp do  
  do i=1,ndata
    f_hi(i)=0.0
    xi=xz(i)           
    yi=ylum(i)           	
    f_hi(i)=fi(i,xi,yi)    
    f_hi(i)=f_hi(i)*2/(2*nw-ni(i)) * ( 1/(red(i)-z1) + 1/(z2-red(i)) ) 
    ftemp = ftemp + log( f_hi(i) )
    !print*,'here',i,f_hi(i)    	
  end do    
  !!$omp end do
  !!$omp ATOMIC
  log_lik=log_lik + ftemp
  !!$omp end parallel  

  !***********************************************************************************************

  if(unbounded_lum) then
    ans = log_lik
  else
    if (absolute_magnitude) then
      L2=Lmin
      call myquad2d(z1,z2,X2,X1,p,int2d)
    else
      L2=Lmax
      call myquad2d(z1,z2,X1,X2,p,int2d)    
    end if 
    !print*,'here3',int2d
    ans = log_lik - nw*int2d
  end if    
  return
end subroutine lnlike

!################################################################################################
!                Calculate expected cumulative distribution functions (CDFs)  
!************************************************************************************************
!subroutine czf(z,nd,ans)
!  use params
!  use myquad
!  implicit none
!  integer :: nd,i
!  real(8) :: ans(nd),z(nd),area
!!f2py intent(in) :: z 
!!integer intent(hide),depend(z) :: nd=shape(z) 
!!f2py intent(out) :: ans
!  real(8),external :: p,p_ada,p1da,X1,X2,pp

!  L1=Lmin
!  L2=Lmax  
!  do i=1,ndata
!    if (small_sample_approximation) then
!      call myquad2d(z1,z(i),X1,X2,p1da,area)
!    else
!      if (adaptive) then        
!        call myquad2d(z1,z(i),X1,X2,p_ada,area)
!      else             
!        call myquad2d(z1,z(i),X1,X2,p,area)
!      end if
!    end if  
!    ans(i)=area
!    !write(*,'(5f10.6)'),z,area
!  end do
!  return

!end subroutine czf


!subroutine clf1d(lum,nd,ans)
!  use params
!  use dqag_and_dqags
!  implicit none
!  integer :: nd,i
!  real(8) :: lum(nd),ans(nd),area
!!f2py intent(in) :: lum,z 
!!integer intent(hide),depend(lum) :: nd=shape(lum) 
!!f2py intent(out) :: ans
!  real(8),external :: f1d,f1da
!  integer ( kind = 4 ), parameter :: limit = 500
!  integer ( kind = 4 ), parameter :: lenw = limit * 4
!  real ( kind = 8 ) a
!  real ( kind = 8 ) abserr
!  real ( kind = 8 ) b
!  real ( kind = 8 ), parameter :: epsabs = 0.0D+00
!  real ( kind = 8 ), parameter :: epsrel = 0.001D+00
!  integer ( kind = 4 ) ier
!  integer ( kind = 4 ) iwork(limit)
!  integer ( kind = 4 ) last
!  integer ( kind = 4 ) neval
!  real ( kind = 8 ) work(lenw)
!  
!  do i=1,ndata
!    a=0
!    b=lum(i)
!    if (adaptive) then     
!      call dqags ( f1da, a, b, epsabs, epsrel, area, abserr, neval, ier, limit, lenw, last, iwork, work )        
!    else
!      call dqags ( f1d, a, b, epsabs, epsrel, area, abserr, neval, ier, limit, lenw, last, iwork, work )
!    end if
!    ans(i)=area
!    !write(*,'(5f10.6)'),L2,area
!  end do
!  return
!end subroutine clf1d


!subroutine clf(lum,nd,z2s,ans)
!  use params
!  use myquad
!  implicit none
!  integer :: nd,i
!  real(8) :: lum(nd),ans(nd),z2s(nd),area
!!f2py intent(in) :: lum,z 
!!integer intent(hide),depend(lum) :: nd=shape(lum) 
!!f2py intent(out) :: ans
!  real(8),external :: p,p_ada,p1d,p1da,X1,X2,pp,pw

!  do i=1,ndata
!    L1=Lmin
!    L2=lum(i)
!    if (small_sample_approximation) then
!      if (adaptive) then     
!        call myquad2d(z1,z2s(i),X1,X2,p1da,area)        
!      else
!        call myquad2d(z1,z2s(i),X1,X2,p1d,area)
!      end if    
!    else
!      if (adaptive) then        
!        call myquad2d(z1,z2s(i),X1,X2,p_ada,area)
!      else        
!        if (absolute_magnitude) then
!          call myquad2d(z1,z2s(i),X2,X1,pw,area)      
!        else
!          call myquad2d(z1,z2s(i),X1,X2,p,area)
!        end if
!      end if
!    end if  
!    ans(i)=area
!    write(133,'(I4,2f10.6)'),i,L2,area
!  end do
!  return

!end subroutine clf



subroutine clf_new(lum,nd,z2s,ans)
  use params
  use myquad
  use prob_functions
  implicit none
  integer :: nd
  real(8) :: lum(nd),ans(nd),z2s(nd),ai
!f2py intent(in) :: lum,z 
!integer intent(hide),depend(lum) :: nd=shape(lum) 
!f2py intent(out) :: ans
  !real(8),external :: p,pw
  real(8),external :: X1,X2,g_lim,f_lim
  real(8), allocatable :: st(:)
  real(8),parameter :: dl=0.0002
  real(8),parameter :: dz=0.0002
  real(8) dlj,lc,z,zc,zt,gj,gjj,dli,f_lim_z2
  integer i,j,n
  procedure(p2d), pointer :: p

  if (weighting) then
    p => pw2d
  else
    p => p2d
  end if

  f_lim_z2=f_lim(z2)
  DO i=1,nd    
    IF (i==1) then
      L2=lum(i)
      if (absolute_magnitude) then
        call myquad2d(z1,z2s(i),X2,X1,p,ai)
      else
        call myquad2d(z1,z2s(i),X1,X2,p,ai)
      end if  
      ans(i)=ai      
    ELSE
      dli=abs(lum(i)-lum(i-1))
      if ( dli>dl ) then
        n= max( floor(dli/dl),2 )
      else
        n=1
      end if        
      
      write(121,*) i,dli,n
      
      allocate(st(n+1))
      dlj=dli/n
      if (absolute_magnitude) then   !*******absolute_magnitude is .true.
          st(1)=lum(i-1)
          do j=2,n+1
            st(j)=st(j-1)-dlj
          end do      
          gj=0.0
          do j=1,n
            lc=st(j)-dlj/2      
            if (lc>f_lim_z2) then       
              zt=g_lim(lc)
            else
              zt=z2
            end if    
            z=z1+dz
            gjj=0.0
            do while(z<=zt)
              zc=z-dz/2
              gjj=gjj+p(zc,lc)*dz*dlj         
              z=z+dz
            end do
            gj=gj+gjj
          end do                       
      
      else                         !*******absolute_magnitude is .false.             
          st(1)=lum(i-1)
          do j=2,n+1
            st(j)=st(j-1)+dlj
          end do       
          gj=0.0
          do j=1,n
            lc=st(j)+dlj/2      
            if (lc<=f_lim_z2) then       
              zt=g_lim(lc)
            else
              zt=z2
            end if    
            z=z1+dz
            gjj=0.0
            do while(z<=zt)
              zc=z-dz/2
              gjj=gjj+p(zc,lc)*dz*dlj    
              z=z+dz
            end do
            gj=gj+gjj
          end do
      end if                      !*******absolute_magnitude judgement end if   
      deallocate(st)
      ans(i)=ans(i-1)+gj    
    END IF
    write(113,'(I4,2f10.6)') i,lum(i),ans(i)
  END DO  
  return
end subroutine clf_new


!subroutine clf1(lum,nd,z2s,ans)
!  use params
!  use myquad
!  implicit none
!  integer :: nd,i,n
!  real(8) :: lum(nd),ans(nd),z2s(nd),area
!!f2py intent(in) :: lum,z 
!!integer intent(hide),depend(lum) :: nd=shape(lum) 
!!f2py intent(out) :: ans
!  real(8),external :: p,p_ada,p1da,X1,X2,pp

!  process=10
!  n=3
!  !!$omp parallel NUM_THREADS(n) private(i,area)  	
!  !!$omp do
!  do i=1,ndata
!    L1=Lmin
!    L2=lum(i)
!    if (small_sample_approximation) then
!      call myquad2d(z1,z2s(i),X1,X2,p1da,area)
!    else
!      if (adaptive) then        
!        call myquad2d(z1,z2s(i),X1,X2,p_ada,area)
!      else        
!        !call myquad2d(z1,z2,X1,X2,pp,area)      
!        call myquad2d(z1,z2s(i),X1,X2,p,area)
!      end if
!    end if  
!    ans(i)=area
!    !write(*,'(I5,5f10.6)'),i,L2,area
!  end do
!  !!$omp end do
!  !!$omp end parallel
!  return

!end subroutine clf1

!subroutine clf2(lum,nd,z2s,ans)
!  use params
!  use myquad
!  implicit none
!  integer :: nd,i,n
!  real(8) :: lum(nd),ans(nd),z2s(nd),area
!!f2py intent(in) :: lum,z 
!!integer intent(hide),depend(lum) :: nd=shape(lum) 
!!f2py intent(out) :: ans
!  real(8),external :: p,p_ada,p1da,X1,X2,pp
!  
!  n=10
!  do i=1,ndata
!    if((i-n)<=0) then
!      L1=Lmin
!    else
!      L1=lum(i-n)
!    end if
!    L2=lum(i)     

!    if (small_sample_approximation) then
!      call myquad2d(z1,z2s(i),X1,X2,p1da,area)
!    else
!      if (adaptive) then        
!        call myquad2d(z1,z2s(i),X1,X2,p_ada,area)
!      else
!        !call myquad2d(z1,z2,X1,X2,pp,area)      
!        call myquad2d(z1,z2s(i),X1,X2,p,area)
!      end if
!    end if    
!    
!    if((i-n)<=0) then
!      ans(i)=area
!    else
!      ans(i)=ans(i-n)+area
!    end if    
!    
!    !write(*,'(5f10.6)'),L2,area
!  end do
!  return

!end subroutine clf2



!subroutine clf2d(lum,nd,z,ans)
!  use params
!  use myquad
!  implicit none
!  integer :: nd,i
!  real(8) :: lum(nd),ans(nd),z(nd),area
!!f2py intent(in) :: lum,z 
!!integer intent(hide),depend(lum) :: nd=shape(lum) 
!!f2py intent(out) :: ans
!  real(8),external :: p,p_ada,p1da,X1,X2,pp
!    
!  !l2lim=f_lim(z2)
!  do i=1,ndata
!    !if(lum(i)<l2lim) then
!    !  L1=lum(i)
!    !end if     
!    L1=Lmin
!    L2=lum(i)    
!    if (small_sample_approximation) then
!      call myquad2d(z1,z(i),X1,X2,p1da,area)
!    else
!      if (adaptive) then        
!        call myquad2d(z1,z(i),X1,X2,p_ada,area)
!      else        
!        !call myquad2d(z1,z2,X1,X2,pp,area)      
!        call myquad2d(z1,z(i),X1,X2,p,area)
!      end if
!    end if  
!    ans(i)=area
!    !write(*,'(3f10.6)'),lum(i),z(i),area
!  end do
!  return

!end subroutine clf2d


!************************************************************************************************
!################################################################################################
real(8) function K(x,y)
  implicit none
  real(8),parameter :: pi=3.141592653589793  
  real(8) x,y
  
  K=1/(2*pi)*exp(-0.5*(x*x+y*y))
end function

!***************************************
subroutine prob(z,L,ans)
  use params
  use prob_functions
  implicit none
  real(8) z,L,ans
!f2py intent(in) :: z,L
!f2py intent(out) :: ans
  procedure(p2d), pointer :: p
  
  if (adaptive) then
    hi=f_pilot**(-beta)    
    !hi=(1+f_pilot)**(-beta)
    !hi=exp(-beta*f_pilot)  
    if (weighting) then
      p => pw2da
    else
      p => p2da
    end if     
  
  else  
    if (weighting) then
      p => pw2d
    else
      p => p2d
    end if  
  end if  

  ans=p(z,L)
  return
end subroutine prob



subroutine f_ref_f2py(x,y,ans)
  use params
  implicit none
  real(8) x,y,ans,x_minus,y_minus,y_plus,temp,ftemp
  real(8),external :: K
  integer i,num  
  real(8) xi,yi
!f2py intent(in) :: x,y
!f2py intent(out) :: ans

  num=ndata
  ans=0.0
  ftemp=0.0
  !!$omp parallel NUM_THREADS(process) private(i,xi,yi,x_minus,y_minus,y_plus,temp,ftemp) shared(ans)  	
  !!$omp do
  do i=1,num
    xi=xz(i)                 
    yi=ylum(i)                 
    x_minus=(x-xi)/h1
    y_minus=(y-yi)/h2
    y_plus=(y+yi)/h2
    temp= K(x_minus,y_minus)+K(x_minus,y_plus)
    ftemp=ftemp + temp
  end do
  !!$omp end do
  !!$omp ATOMIC
  ans=ans + ftemp
  !!$omp end parallel
  ans=ans/(h1*h2*num)
  return 
end subroutine f_ref_f2py

!***************************************
real(8) function f_ref(x,y)
  use params
  implicit none
  real(8) x,y,x_minus,y_minus,y_plus,temp,ftemp
  real(8),external :: K
  integer i,num  
  real(8) xi,yi

  num=ndata
  f_ref=0.0
  ftemp=0.0
  !!$omp parallel NUM_THREADS(process) private(i,xi,yi,x_minus,y_minus,y_plus,temp,ftemp) shared(f_ref)  	
  !!$omp do
  do i=1,num
    xi=xz(i)                 
    yi=ylum(i)                
    x_minus=(x-xi)/h1
    y_minus=(y-yi)/h2
    y_plus=(y+yi)/h2
    temp= K(x_minus,y_minus)+K(x_minus,y_plus)
    ftemp=ftemp + temp
  end do
  !!$omp end do
  !!$omp ATOMIC
  f_ref=f_ref + ftemp
  !!$omp end parallel
  f_ref=f_ref/(h1*h2*num)
  return  
end function

!***************************************
real(8) function pp(z,L)
  use params
  implicit none
  real(8) x,y,z,L,bandwidth
  real(8),external :: f_lim,f_ref

  x=log( (z-z1)/(z2-z) )
  y=L-f_lim(z)
  if (z>z1 .and. z<z2 .and. L>f_lim(z)) then 
    pp=f_ref(x,y) * ( 1/(z-z1) + 1/(z2-z) )
  else
    pp=0.0
  end if    
end function


!#############################################################################

real(kind=8) function X1(z)
  use params
  implicit none
  real(kind=8) z
  real(8),external :: f_lim 

  X1=f_lim(z)
  return
end function

real(kind=8) function X2(z)
  use params
  implicit none
  real(kind=8) z
  
  X2 = L2
  return
end function 
!************************************************************************************************
!################################################################################################
real(8) function f_ada(x,y)
  use params
  implicit none
  real(8) x,y,ha1,ha2
  real(8) xi,yi,x_m,y_m,y_p,ftemp
  real(8),external :: K
  integer i

  f_ada=0.0
  ftemp=0.0
  !!$omp parallel NUM_THREADS(process) private(i,xi,yi,ha1,ha2,x_m,y_m,y_p,ftemp) shared(f_ada)  	
  !!$omp do
  do i=1,ndata
    xi=xz(i)                   
    yi=ylum(i)            
    ha1=h10*hi(i)
    ha2=h20*hi(i)
    x_m=(x-xi)/ha1
    y_m=(y-yi)/ha2
    y_p=(y+yi)/ha2
    ftemp=ftemp + ( K(x_m,y_m)+K(x_m,y_p) )/(ha1*ha2)
  end do  	
  !!$omp end do
  !!$omp ATOMIC
  f_ada=f_ada + ftemp
  !!$omp end parallel
  f_ada=f_ada/ndata
end function

!################################################################################################
!************************************************
subroutine lnlike_ada(h1_py,h2_py,beta_py,ans)
  use params
  use myquad
  use prob_functions
  use fhi_functions
  implicit none
  real(8) :: h1_py,h2_py,beta_py,ans
!f2py intent(in) :: h1_py,h2_py,beta_py
!f2py intent(out) :: ans 
  real(8) :: f_hi(ndata)
  real(8) temp,x_m,y_m,y_p
  real(8) log_lik,int2d
  real(8) xi,yi,ftemp,xj,yj
  integer i,j
  !real(8),external :: fhi_ada,fwhi_ada
  real(8),external :: X1,X2   !,p_ada,pw_ada
  procedure(p2d), pointer :: p
  procedure(fhi), pointer :: fi

  if (weighting) then
    p => pw2da
    fi => fwhi_ada 
  else
    p => p2da
    fi => fhi_ada
    nw=ndata*1.0
  end if
    
  h10=h1_py
  h20=h2_py  
  beta=beta_py
  if (.not. set_beta_fixed) then  
    hi=f_pilot**(-beta)
    !hi=(1+f_pilot)**(-beta)
    !hi=exp(-beta*f_pilot) 
  end if  
  
  ftemp=0.0
  log_lik=0.0

  !******************************************************************************
  !!$omp parallel NUM_THREADS(process) private(i,xi,yi,ftemp) shared(log_lik)  	
  !!$omp do  
  do i=1,ndata
    f_hi(i)=0.0
    xi=xz(i)           
    yi=ylum(i)           	
    f_hi(i)=fi(i,xi,yi) 
    f_hi(i)=f_hi(i)*2/(2*nw-ni(i)) * ( 1/(red(i)-z1) + 1/(z2-red(i)) )
    ftemp = ftemp + log( f_hi(i) )   	
  end do    
  !!$omp end do
  !!$omp ATOMIC
  log_lik=log_lik + ftemp
  !!$omp end parallel  
  
  if(unbounded_lum) then
    ans = log_lik
  else
    if (absolute_magnitude) then
      L2=Lmin
      call myquad2d(z1,z2,X2,X1,p,int2d)
    else
      L2=Lmax
      call myquad2d(z1,z2,X1,X2,p,int2d)    
    end if 
    !print*,'here3',int2d
    ans = log_lik - nw*int2d
  end if 

  return
end subroutine lnlike_ada


!#######################################################################################################################

subroutine lnlike_1d(h2_py,ans)
  use params
  use myquad
  implicit none
  real(8) :: h2_py,ans
!f2py intent(in) :: h2_py
!f2py intent(out) :: ans 
  real(8) :: f_hi(ndata)
  real(8) temp
  real(8) yi,ftemp,yj
  real(8) int2d,log_lik
  integer i,j
  real(8),external :: fhi_1d
  real(8),external :: p1d,X1,X2    

  h2=h2_py  
  ftemp=0.0
  log_lik=0.0
  !!$omp parallel NUM_THREADS(process) private(i,yi,ftemp) shared(log_lik)  	
  !!$omp do
  do i=1,ndata
    f_hi(i)=0.0          
    yi=ylum(i)           	
    f_hi(i)=fhi_1d(i,yi)    
    f_hi(i)=f_hi(i)*2/(2*ndata-1)/(z2-z1)
    ftemp = ftemp + log( f_hi(i) )    	
  end do    
  !!$omp end do
  !!$omp ATOMIC
  log_lik=log_lik + ftemp
  !!$omp end parallel
  
  L1=Lmin
  L2=Lmax
  call myquad2d(z1,z2,X1,X2,p1d,int2d)  
  ans = log_lik - ndata*int2d  
  return
end subroutine lnlike_1d


real(8) function fhi_1d(i,yi)
  use params
  implicit none
  real(8) yi,y_m,y_p,yj,temp
  integer i,j
  
  fhi_1d=0.0
  do j=1,ndata
    yj=ylum(j)           
    y_m=(yi-yj)/h2
    y_p=(yi+yj)/h2
    if(j /= i) then
    	temp= ( exp(-y_m**2/2) + exp(-y_p**2/2) )    
    else
    	temp=exp(-y_p**2/2)	  
    end if
    fhi_1d=fhi_1d + temp
  end do
  fhi_1d=fhi_1d/(h2*sqrt(2*pi))
end function


!***************************************
real(8) function f1d(y)
  use params
  implicit none
  real(8) y,y_minus,y_plus,temp,ftemp,f_ref
  integer i,num  
  real(8) xi,yi

  num=ndata
  f_ref=0.0
  ftemp=0.0
  !!$omp parallel NUM_THREADS(process) private(i,yi,y_minus,y_plus,temp,ftemp) shared(f_ref)  	
  !!$omp do
  do i=1,num                     
    yi=ylum(i)    
    y_minus=(y-yi)/h2
    y_plus=(y+yi)/h2
    temp=( exp(-y_minus**2/2) + exp(-y_plus**2/2) )/sqrt(2*pi)    
    ftemp=ftemp + temp
  end do
  !!$omp end do
  !!$omp ATOMIC
  f_ref=f_ref + ftemp
  !!$omp end parallel
  f1d=f_ref/(h2*num)
  return  
end function

!***************************************
real(8) function p1d(z,L)
  use params
  implicit none
  real(8) y,z,L
  real(8),external :: f_lim,f1d

  y=L-f_lim(z)
  p1d=f1d(y)/(z2-z1)
end function

!#############################################################################
subroutine lnlike_1da(h2_py,beta_py,ans)
  use params
  use myquad
  implicit none
  real(8) :: h2_py,beta_py,ans
!f2py intent(in) :: h2_py,beta_py
!f2py intent(out) :: ans 
  real(8) :: f_hi(ndata)
  real(8) temp
  real(8) yi,ftemp,yj
  real(8) int2d,log_lik
  integer i,j
  real(8),external :: fhi_1da
  real(8),external :: p1da,X1,X2    

  h20=h2_py
  beta=beta_py  
  hi=f_pilot**(-beta)
  
  ftemp=0.0
  log_lik=0.0
  !!$omp parallel NUM_THREADS(process) private(i,yi,ftemp) shared(log_lik)  	
  !!$omp do
  do i=1,ndata              
    yi=ylum(i)           	
    f_hi(i)=fhi_1da(i,yi)    
    f_hi(i)=f_hi(i)*2/(2*ndata-1)/(z2-z1)
    ftemp = ftemp + log( f_hi(i) )    	
  end do    
  !!$omp end do
  !!$omp ATOMIC
  log_lik=log_lik + ftemp
  !!$omp end parallel
  
  L1=Lmin
  L2=Lmax
  call myquad2d(z1,z2,X1,X2,p1da,int2d)  
  ans = log_lik - ndata*int2d  
  return
end subroutine lnlike_1da


!*****************************************************
real(8) function f1da(y)
  use params
  implicit none
  real(8) y,h
  real(8) yi,y_m,y_p,ftemp
  !real(8),external :: K
  integer i
  
  f1da=0.0
  ftemp=0.0
  !!$omp parallel NUM_THREADS(process) private(i,yi,h,y_m,y_p,ftemp) shared(f1da)  	
  !!$omp do
  do i=1,ndata
    yi=ylum(i) 	            
    h=h20*hi(i)
    y_m=(y-yi)/h
    y_p=(y+yi)/h
    ftemp=ftemp + ( exp(-y_m**2/2) + exp(-y_p**2/2) )/sqrt(2*pi) /h       
  end do  	
  !!$omp end do
  !!$omp ATOMIC
  f1da=f1da + ftemp
  !!$omp end parallel
  f1da=f1da/ndata
end function

real(8) function fhi_1da(i,yi)
  use params
  implicit none
  real(8) yi,yj,h
  real(8) y_m,y_p,temp
  integer i,j
  
  fhi_1da=0.0
  do j=1,ndata
    yj=ylum(j)           	    
    h=h20*hi(j)
    y_m=(yi-yj)/h
    y_p=(yi+yj)/h
    if(j /= i) then
    	temp=( exp(-y_m**2/2) + exp(-y_p**2/2) )/sqrt(2*pi) /h        !( K(y_m) + K(y_p) ) / h      
    else
    	temp= exp(-y_p**2/2) /sqrt(2*pi) /h  	  
    end if
    fhi_1da=fhi_1da + temp
  end do
end function

real(8) function p1da(z,L)
  use params
  implicit none
  real(8) y,z,L
  real(8),external :: f_lim,f1da
  
  y=L-f_lim(z)
  p1da=f1da(y)/(z2-z1)
  return   
end function

!################################################################################################
!             The KDE Method Considering the Weighting Due to the Selection Function  
!************************************************************************************************
real(8) function fw_ref(x,y)
  use params
  implicit none
  real(8) x,y,x_minus,y_minus,y_plus,temp,ftemp
  real(8),external :: K
  integer i,num  
  real(8) xi,yi

  fw_ref=0.0
  ftemp=0.0
  !!$omp parallel NUM_THREADS(process) private(i,xi,yi,x_minus,y_minus,y_plus,temp,ftemp) shared(fw_ref)  	
  !!$omp do
  do i=1,ndata
    xi=xz(i)                 
    yi=ylum(i)                
    x_minus=(x-xi)/h1
    y_minus=(y-yi)/h2
    y_plus=(y+yi)/h2
    temp= K(x_minus,y_minus)+K(x_minus,y_plus)
    ftemp=ftemp + temp*weight(i)
  end do
  !!$omp end do
  !!$omp ATOMIC
  fw_ref=fw_ref + ftemp
  !!$omp end parallel
  fw_ref=fw_ref/(h1*h2*nw)
  return  
end function

!#############################################################################

real(8) function fw_ada(x,y)
  use params
  implicit none
  real(8) x,y,ha1,ha2
  real(8) xi,yi,x_m,y_m,y_p,ftemp
  real(8),external :: K
  integer i,num

  fw_ada=0.0
  ftemp=0.0
  !!$omp parallel NUM_THREADS(process) private(i,xi,yi,ha1,ha2,x_m,y_m,y_p,ftemp) shared(fw_ada)  	
  !!$omp do
  do i=1,ndata
    xi=xz(i)                   
    yi=ylum(i)            
    ha1=h10*hi(i)
    ha2=h20*hi(i)
    x_m=(x-xi)/ha1
    y_m=(y-yi)/ha2
    y_p=(y+yi)/ha2
    ftemp=ftemp + ( K(x_m,y_m)+K(x_m,y_p) )*weight(i)/(ha1*ha2)
  end do  	
  !!$omp end do
  !!$omp ATOMIC
  fw_ada=fw_ada + ftemp
  !!$omp end parallel
  fw_ada=fw_ada/nw
end function

!################################################################################################
!         !!! lscv method, not used in current version  
!************************************************************************************************

subroutine lscv(h,ans)
  use params
  !use myquad
  implicit none
  real(8) :: h,ans
!f2py intent(in) :: h
!f2py intent(out) :: ans 
  real(8) xi,xj,x_m,x_p
  real(8) K2hi_m,K2hi_p,Khi_m,Khi_p
  real(8) lscv1,lscv2
  real(8) f_hi(ndata),f_2hi(ndata)   
  real(8) xm2,xp2
  integer i,j 

  do i=1,ndata
    f_hi(i)=0.0
    f_2hi(i)=0.0
    xi=ylum(i)
    do j=1,ndata               
      xj=ylum(j)  	    
      x_m=(xi-xj)/h
      x_p=(xi+xj)/h
      xm2=-(x_m**2)/2
      xp2=-(x_p**2)/2      
      Khi_m=exp(xm2)               
      Khi_p=exp(xp2)            
      K2hi_m=exp(xm2/2)         
      K2hi_p=exp(xp2/2)         
      f_2hi(i) = f_2hi(i) + (K2hi_m + K2hi_p) 
      f_hi(i) =f_hi(i) + (Khi_m + Khi_p)
    end do  
    f_hi(i)=f_hi(i)-1.0
  end do 
  lscv1=sum(f_2hi)/(2*sqrt(pi))/(2*ndata*ndata*h)
  lscv2=sum(f_hi)/(sqrt(2*pi)) /(ndata*(2*ndata-1)*h) 
  !print*,'lscv1,lscv2',lscv1,lscv2
  ans=lscv1 - 2*lscv2
  return
end subroutine lscv

subroutine lscv2d(h1_py,h2_py,ans)
  use params
  !use myquad
  implicit none
  real(8) :: h1_py,h2_py,ans
!f2py intent(in) :: h1_py,h2_py
!f2py intent(out) :: ans 
  real(8) xi,yi,xj,yj,x_m,y_m,y_p
  real(8) K2hi_m,K2hi_p,Khi_m,Khi_p
  real(8) lscv1,lscv2
  real(8) f_hi(ndata),f_2hi(ndata)   
  real(8) temp1,temp2,xmym,xmyp
  integer i,j  
    
  h1=h1_py
  h2=h2_py

  do i=1,ndata
    f_hi(i)=0.0
    f_2hi(i)=0.0
    xi=xz(i)                 
    yi=ylum(i)
    do j=1,ndata
      xj=xz(j)                 
      yj=ylum(j)  	    
      x_m=(xi-xj)/h1  
      y_m=(yi-yj)/h2
      y_p=(yi+yj)/h2
      xmym=-(x_m**2 + y_m**2)/2
      xmyp=-(x_m**2 + y_p**2)/2        
      Khi_m=exp(xmym)                 !Khi_m=exp(-(x_m**2 + y_m**2)/2)
      Khi_p=exp(xmyp)                 !Khi_p=exp(-(x_m**2 + y_p**2)/2)
      K2hi_m=exp(xmym/2)              !K2hi_m=exp(-(x_m**2 + y_m**2)/4)
      K2hi_p=exp(xmyp/2)              !K2hi_p=exp(-(x_m**2 + y_p**2)/4)
      f_2hi(i) = f_2hi(i) + (K2hi_m + K2hi_p) 
      f_hi(i) =f_hi(i) + (Khi_m + Khi_p)
    end do  
    f_hi(i)=f_hi(i)-1.0
  end do 
  lscv1=sum(f_2hi)/(4*h1*h2*pi)/(2*ndata*ndata)
  lscv2=sum(f_hi)/(2*h1*h2*pi) *2/(ndata*(2*ndata-1)) 
  !print*,'3:lscv1,lscv2',lscv1,lscv2
  ans=lscv1 - lscv2
  return
end subroutine lscv2d  


!################################################################################################
!         !!! interpolation to calculate the f_lim function  
!************************************************************************************************
real(8) function f_lim(z)
  use params
  implicit none
  integer,external :: locate
  real(8) z
  integer n,nlim
  
  nlim=size(limx)
  n=locate(limx,z,nlim)
  f_lim=(z-limx(n))*(limy(n+1)-limy(n)) / (limx(n+1)-limx(n)) + limy(n)    

  !print*,z,n
end function


real(8) function g_lim(L)   ! g_lim is the inverse function of f_lim
  use params
  implicit none
  integer,external :: locate
  real(8) L
  integer n,nlim
  
  nlim=size(limy)
  n=locate(limy,L,nlim)
  g_lim=(L-limy(n))*(limx(n+1)-limx(n)) / (limy(n+1)-limy(n)) + limx(n)    

  !print*,z,n
end function

	FUNCTION locate(xx,x,nn)	
	IMPLICIT NONE
	integer nn
	real(8) xx(nn)
	!REAL(8), DIMENSION(:), INTENT(IN) :: xx
	
	REAL(8), INTENT(IN) :: x
	INTEGER :: locate
	INTEGER :: n,jl,jm,ju
	LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	jl=0
	ju=n+1
	do
		if (ju-jl <= 1) exit
		jm=(ju+jl)/2
		if (ascnd .eqv. (x >= xx(jm))) then
			jl=jm
		else
			ju=jm
		end if
	end do
	if (x == xx(1)) then
		locate=1
	else if (x == xx(n)) then
		locate=n-1
	else
		locate=jl
	end if
	END FUNCTION locate



!subroutine time_lim()
!  use params
!  real(8),external :: f_lim1,f_lim2,f_lim
!  real(8) t1,t2,t3,t4,y
!  integer i,j
!  
!  call cpu_time(t1)
!  do j=1,100
!    do i=1,ndata
!      y=f_lim1(red(i))
!    end do
!  end do  
!  call cpu_time(t2)
!  
!  do j=1,100
!    do i=1,ndata
!      y=f_lim2(red(i))
!    end do
!  end do      
!  call cpu_time(t3)
!  do j=1,100
!    do i=1,ndata
!      y=f_lim(red(i))
!    end do    
!  end do
!  call cpu_time(t4)  
!  print*,'f_lim1,','f_lim2,','f_lim'
!  print*,t2-t1,t3-t2,t4-t3
!  return    
!end subroutine time_lim


real(8) function f_lim1(z)
  use params
  implicit none
  integer :: locate
  integer :: n,jl,jm,ju
  LOGICAL :: ascnd
  integer :: nlim
  real(8) z
  
  nlim=size(limx)
  ascnd = (limx(nlim) >= limx(1))
  jl=0
  ju=n+1
  do while(.true.)
      if (ju-jl <= 1) exit
      jm=(ju+jl)/2
      if (ascnd .eqv. (z >= limx(jm))) then
        jl=jm
      else
        ju=jm
      end if
  end do
  if (z == limx(1)) then
	  locate=1
  else if (z == limx(nlim)) then
	  locate=n-1
  else
	  locate=jl
  end if

  n=locate
  f_lim1=(z-limx(n))*(limy(n+1)-limy(n)) / (limx(n+1)-limx(n)) + limy(n) 
  !print*,z,n
end function

subroutine trap2d(pun,ans)
  use params
  implicit none
  real(8) dz,ans,a,z0,zc,lc,Ai,Le,h,s,vj
  integer i,j,nlim
  real(8),external :: pun 
  
  nlim=size(limx)
  z0=limx(1)
  dz=(z2-z0)/(nlim-1)   
  ans=0.0
  do i=1,nlim-1
    if (i==1) then
      a=dz+z0-z1
      zc=z0+dz/2
    else
      a=dz
      zc=limx(i)+dz/2
    end if
    j=1
    Ai=0.0
    Le=0.0
    
    IF (absolute_magnitude) then    
      do while (Le>=L2)      
        if (j==1) then
          h=limy(i)-limy(i+1)
          s=a*h/2       
          lc=limy(i)-0.7*h
          Le=limy(i+1)
          j=2
        else
          h=0.005
          s=a*h
          lc=Le-h/2   
          Le=Le-h
        end if  
        vj=s*pun(zc,lc)
        Ai=Ai+vj
      end do   
  
    ELSE  
      do while (Le<=L2)      
        if (j==1) then
          h=limy(i+1)-limy(i)
          s=a*h/2       
          lc=limy(i)+0.7*h
          Le=limy(i+1)
          j=2
        else
          h=0.005
          s=a*h
          lc=Le+h/2   
          Le=Le+h
        end if  
        vj=s*pun(zc,lc)
        Ai=Ai+vj
      end do
    END IF  
      ans=ans+Ai 
  end do  
  return
end subroutine trap2d  
    
    































      





  

module ucb_heaviside

  use shared_input_parameters, only: SOURCE_T1, SOURCE_T2, SOURCE_T3, &
                               SOURCE_T4, TAU


  implicit none

  private

  ! expose only initialization and
  public :: ucb_stf, init_ucb_heaviside

  ! ============ Modified by <FM> January, 2022
  ! hard-wired source passband, origin time, and (optional) amplitude are set
  ! here:
!   double precision, parameter :: &
!        f1h = 1.d0 / 400.d0, &
!        f2h = 1.d0 / 250.d0, &
!        f3h = 1.d0 /  53.d0, &
!        f4h = 1.d0 /  40.d0, &
!        t0  = 400.d0, &
!        amp = 1.d7
  ! =======================================================
  ! source amplitude (optional)
  double precision, parameter ::  amp = 1.d0
  double precision :: f1h, f2h, f3h, f4h, t0
  ! =================== end modification ===================

  ! the STF type (displacement vs. velocity vs. acceleration) is set below in init_ucb_heaviside

  logical :: initialized = .false.

  integer :: nstep, itsource_max

  double precision :: dt

  double precision, dimension(:), allocatable :: g

contains

  function ucb_stf(t)
    implicit none
    double precision :: ucb_stf
    double precision, intent(in) :: t
    ! --
    integer :: it
    ! --
    if (.not. initialized) stop "Error [ucb_heaviside]: call to ucb_stf _before_ STF init"
    ! infer the time-step index, always starting from zero time
    it = 1 + NINT(max(0.d0,t) / dt)
    ucb_stf = g(min(it,nstep))
  end function ucb_stf

  subroutine init_ucb_heaviside(nstep_,dt_)
    implicit none
    ! --
    doubleprecision, parameter :: zero = 0.d0
    ! --
    integer, intent(in) :: nstep_
    doubleprecision, intent(in) :: dt_
    ! --
    integer :: it0,j,nstep2,i,ic,istat
    doubleprecision :: wt,freq,t1,t2,tmax,seuil,pi
    complex*16 :: dphi
    complex*16, dimension(:), allocatable :: spectre
    !
    ! Added <FM>, January 2022
    f1h = 1.d0/SOURCE_T1
    f2h = 1.d0/SOURCE_T2
    f3h = 1.d0/SOURCE_T3
    f4h = 1.d0/SOURCE_T4
    t0  = TAU
    write(*,*) "=============================="
    WRITE(*,*) "Heaviside Source T1:", SOURCE_T1
    WRITE(*,*) "Heaviside Source T2:", SOURCE_T2
    WRITE(*,*) "Heaviside Source T3:", SOURCE_T3
    WRITE(*,*) "Heaviside Source T4:", SOURCE_T4
    WRITE(*,*) "Source time-shift:", t0
    write(*,*) "=============================="
    ! ============================== 

    !
    nstep = nstep_
    dt = dt_
    allocate(g(nstep))
    !
    pi = 4.d0 * datan(1.d0)
    g(:) = zero
    !
    nstep2=int(2.d0**(int(log(dble(nstep))/log(2.d0))+1))
    !
    allocate(spectre(nstep2))
    spectre(:)=cmplx(0.d0,0.d0)
    do j=1,nstep2
       if (j<=nstep2/2) then
          freq=(j-1)/(dt*nstep2)
       else if (j==nstep2/2+1) then
          freq=1/(2.d0*dt)
       else
          freq=-(nstep2-j+1)/(dt*nstep2)
       endif
       dphi=exp(-2.d0*pi*freq*t0*cmplx(0.d0,1.d0))
       call wtcoef(abs(freq),f1h,f2h,f3h,f4h,wt)
       ! >>> STF TYPE IS SET HERE <<<
       ! displacement
       ! if (j/=1) spectre(j)=wt*dphi / cmplx(0.d0, 2.d0 * pi * freq)
       ! velocity
       if (j/=1) spectre(j)=wt*dphi
       ! acceleration
       !if (j/=1) spectre(j)=wt*dphi * cmplx(0.d0, 2.d0 * pi * freq)
    enddo
    call dfour1(spectre,nstep2,1)
    g(:)=amp*real(spectre(1:nstep))/nstep2/dt
    !on met les premier pas de temps a zero:
    tmax=nstep2*dt
    t1=0.d0
    t2=t0/5.d0
    do i=1,nstep
       call wtcoef((i-1)*dt,t1,t2,tmax,tmax,wt)
       g(i)=g(i)*wt
    enddo
    !
    deallocate(spectre,stat=istat)
    if (istat/=0) stop 'time_function deallocate error'
    !
    ! calcul du temps de fin de la source
    !
    seuil=1.d-3*maxval(abs(g(:)))
    itsource_max=-1
    ic=50 !50 pas de temps
    it0=t0/dt+1
    do i=it0,nstep-ic
       if (itsource_max<0) then
          if (maxval(abs(g(i:i+ic)))<=seuil) itsource_max=i
       endif
    enddo
    if (itsource_max<=0) itsource_max=nstep-ic
    !
    initialized = .true.
    !
  contains
    !----------------------------------------------------------------------
    subroutine wtcoef(f,f1,f2,f3,f4,wt_out)
      !----------------------------------------------------------------------
      implicit none
      !
      doubleprecision, intent(in) ::  f,f1,f2,f3,f4
      doubleprecision, intent(out)::  wt_out
      !
      if (f3.gt.f4) stop 'wtcoef: f3>f4 '
      if (f1.gt.f2) stop 'wtcoef: f1>f2 '
      if (f.le.f3.and.f.ge.f2) then
         wt_out=1.d0
      else if (f.gt.f4.or.f.lt.f1 ) then
         wt_out=0.d0
      else if (f.gt.f3.and.f.le.f4) then
         wt_out=0.5d0*(1.d0+cos(pi*(f-f3)/(f4-f3)))
      else if (f.ge.f1.and.f.lt.f2) then
         wt_out=0.5d0*(1.d0+cos(pi*(f-f2)/(f2-f1)))
      endif
      !----------------------------------------------------------------------
    end subroutine wtcoef
    !----------------------------------------------------------------------
  end subroutine init_ucb_heaviside

end module ucb_heaviside

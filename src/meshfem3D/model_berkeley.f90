!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
! Berkeley 3D model data
!--------------------------------------------------------------------------------------------------
  !
  module model_berkeley_par

  implicit none
  !
  double precision, dimension(:), allocatable :: mdl,kntrad,aknot,oknot,aknot2,oknot2
  integer, dimension(:), allocatable :: level,level2
  integer, parameter :: MAXPAR = 4  
  integer, dimension(MAXPAR) :: nknotA1,nknotA2
  character(1), dimension(:), allocatable :: parblock
  integer :: npar,ndisc,surface,NBPARAM
  logical :: hknots2_exist=.false.,unconformal=.false.
  !

  doubleprecision :: mohodepth = -1.0
    
   contains
   !! subroutine to read earth model
   subroutine det_moho_radius(model_file, mohodp)
      implicit none

      !! Declare vars
      integer :: FID = 10
      character (len=100), intent(in) :: model_file
      character (len=100) :: CTMP
      doubleprecision, allocatable :: radius(:),density(:),vpv(:),vsv(:),qkappa(:),qmu(:),vph(:),vsh(:),eta(:)
      doubleprecision, allocatable :: derivdensity(:), derivVp(:), derivVs(:)
      doubleprecision :: tol = 2.0d0**(-5), maxderivdensity = 0.
      integer :: i = 0, IERR = 0
      integer :: totalheadersval = 3
      integer :: num_lines, totalDiscont
      doubleprecision, intent(inout) :: mohodp

      ! model_file = berk_model1D
      if (mohodp < 0) then
         num_lines = 0 !initialize num of nodes in the file
         open(FID, file=trim(model_file), status="old",iostat=IERR)
         
      
         ! Get number of lines
         do i = 1, totalheadersval
               read( FID, * )   !! skip the header
         end do
         
         do while (IERR == 0)
               num_lines = num_lines + 1
               read(FID,*,iostat=IERR) CTMP
         end do
         num_lines = num_lines - 1
      
      
         ! Allocate array of strings
         allocate(radius(num_lines), density(num_lines), vpv(num_lines), vsv(num_lines), qkappa(num_lines))
         allocate(qmu(num_lines), vph(num_lines),vsh(num_lines),eta(num_lines))
         allocate(derivdensity(num_lines), derivVp(num_lines), derivVs(num_lines))
         
         
      
         ! Read the file content
         rewind(FID)
         do i = 1, totalheadersval
               read( FID, * )   !! skip the header
         end do
         do i = 1, num_lines
               read(FID,*) radius(i), density(i), vpv(i), vsv(i), qkappa(i),qmu(i), vph(i),vsh(i),eta(i)  !Read the data
         end do
      
      
         ! find the discontinuities
         totalDiscont = 0
         do i = 1, num_lines-1
               if (abs(radius(i+1)-radius(i)) < tol)  then
                  derivdensity(i) = density(i+1)-density(i)
                  derivVp(i) = vpv(i+1)-vpv(i)
                  derivVs(i) = vsv(i+1)-vsv(i)
                  
                  totalDiscont = totalDiscont + 1
               else
                  derivdensity(i) = (density(i+1)-density(i))/(radius(i+1)-radius(i))
                  derivVp(i) = (vpv(i+1)-vpv(i))/(radius(i+1)-radius(i))
                  derivVs(i) = (vsv(i+1)-vsv(i))/(radius(i+1)-radius(i))
               end if
         end do
      
         ! Determine the Mohorovicic discontinuity node
         ! Conditions to select the moho discontinuity: 
         ! 1. Radius don't change, hence discontinuity
         ! 2. Vsv(i) and Vsv(i+1) > 0
         ! 3. delta Vsv > 0
         ! 4. max depth of 90km
         ! 5. max density change within 90 km from surface
      
         do i = 1, num_lines-1
               if (abs(radius(i+1)-radius(i)) < tol) then
                  if ((abs(vsv(i)) > tol) .and. (abs(vsv(i+1)) > tol) .and. (abs(derivVs(i)) > tol) &
                  .and. (abs(radius(i)-6371000.) < 90000.)) then
                     
                     if (abs(derivdensity(i)) > maxderivdensity) then
                           maxderivdensity = abs(derivdensity(i))
                           mohodp = radius(i+1)
                     end if
                  end if
               end if
         end do
         
         close(FID)
         deallocate(radius,density,vpv,vsv,qkappa,qmu,vph,vsh,eta)
         deallocate(derivdensity, derivVp, derivVs)

         mohodp = mohodp/1000.d0 !convert to km
         mohodp = mohodp-0.1d0 !adjust output
      
      end if


   end subroutine det_moho_radius
  end module model_berkeley_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_berkeley_broadcast()

! standard routine to setup model

  use model_berkeley_par
  use constants

  implicit none

  integer :: unit1=51,unit2=52,dum2,i,j,k,n,mdim
  character :: trash
  character(len=100), parameter :: A3d_dat            = trim(A3d_folder)//'A3d.dat'
  character(len=100), parameter :: hknots_dat         = trim(A3d_folder)//'hknots.dat'
  character(len=100), parameter :: hknots2_dat        = trim(A3d_folder)//'hknots2.dat'
  !
  ! master reads
  !

  if(myrank==0)then

   print*,'[model_berkeley_broadcast] inside',myrank

  open(unit1,file=trim(hknots_dat),status='old')
  read(unit1,*) nknotA2(1)
  if(allocated(oknot).or.allocated(aknot).or.allocated(level)) then
     print*,'A3d already initiated'
     return
  endif
  allocate(oknot(nknotA2(1)),aknot(nknotA2(1)),level(nknotA2(1)))
  do i=1,nknotA2(1)
     read(unit1,*) oknot(i),aknot(i),level(i)
  enddo
  close(unit1)
  inquire(file=trim(hknots2_dat),exist=hknots2_exist)
  if(hknots2_exist) then
     open(unit1,file=trim(hknots2_dat),status='old')
     read(unit1,*) nknotA2(2)
     if(allocated(oknot2).or.allocated(aknot2).or.allocated(level2)) then
        print*,'A3d already initiated'
        return
     endif
     allocate(oknot2(nknotA2(2)),aknot2(nknotA2(2)),level2(nknotA2(2)))
     do i=1,nknotA2(2)
        read(unit1,*) oknot2(i),aknot2(i),level2(i)
     enddo
     close(unit1)
    else
     nknotA2(2)=nknotA2(1)
  endif

  !open model file and read in model
  open(unit2,file=A3d_dat,status='old')
  read(unit2,*) npar
  NBPARAM=npar
  if(npar>MAXPAR)   stop 'npar greater than MAXPAR'
  allocate(parblock(npar))
  do i=1,npar
     read(unit2,*) dum2,nknotA1(i),parblock(i)
     if(i>1.and.nknotA1(i)/=nknotA1(1)) then
       stop 'Inconsistent A1 splines between parameters'
     endif
     if(i>2) then
        nknotA2(i) = dum2
        if (nknotA2(i)/=nknotA2(2))   stop 'Param 3 and 4 need the same A2 splines than param 2'
     elseif (dum2/=nknotA2(i)) then
        stop 'Inconsistent hknots.dat and A3d.dat'
     endif
     if(i==2.and.nknotA2(i)/=nknotA2(1)) then
        unconformal=.true.
        if(.not.hknots2_exist)   stop 'unconformal grid requires hknots2.dat'
     endif
  enddo
  read(unit2,*) ndisc
  surface=0
  if (ndisc>0) then
     surface=ndisc
     if(unconformal)   print*,'discontinuities assumed same grid as first par'
     do i=1,surface
        read(unit2,*) dum2, trash
     enddo
  endif
  allocate(kntrad(nknotA1(1)))
  read(unit2,*) (kntrad(i),i=1,nknotA1(1))
  mdim=0
  do i=1,npar
     mdim=mdim+nknotA1(i)*nknotA2(i)
  enddo
  mdim=mdim+ndisc*nknotA2(1)
  allocate(mdl(mdim))
  n=0
  do i=1,npar
     do j=1,nknotA1(i)
        read(unit2,*) (mdl(k+n),k=1,nknotA2(i))
        n=n+nknotA2(i)
     enddo
  enddo
  do i=1,ndisc
     read(unit2,*) (mdl(k+n),k=1,nknotA2(1))
     n=n+nknotA2(1)
  enddo

  if(n/=mdim) stop 'init_A3d dimension error'
  close(unit2)

  endif

  !
  ! broadcast
  ! 
  !
  call BCAST_ALL_SINGLEI(nknotA2(1))
  !
  if(.not.allocated(oknot))allocate(oknot(nknotA2(1)))
  if(.not.allocated(aknot))allocate(aknot(nknotA2(1)))
  if(.not.allocated(level))allocate(level(nknotA2(1)))
  !
  call BCAST_ALL_I(level,nknotA2(1))
  call BCAST_ALL_DP(oknot,nknotA2(1))
  call BCAST_ALL_DP(aknot,nknotA2(1))
  !
  call BCAST_ALL_SINGLEL(hknots2_exist)
  !
  if(hknots2_exist) then
    !
    call BCAST_ALL_SINGLEI(nknotA2(2))
    !
    if(.not.allocated(oknot2))allocate(oknot2(nknotA2(2)))
    if(.not.allocated(aknot2))allocate(aknot2(nknotA2(2)))
    if(.not.allocated(level2))allocate(level2(nknotA2(2)))
    !
    call BCAST_ALL_I(level2,nknotA2(2))
    call BCAST_ALL_DP(oknot2,nknotA2(2))
    call BCAST_ALL_DP(aknot2,nknotA2(2))
    !
  else
    !
    nknotA2(2)=nknotA2(1)
    !
  endif
  !
  call BCAST_ALL_SINGLEI(npar)
  NBPARAM=npar

  if(.not.allocated(parblock))allocate(parblock(npar))
  call BCAST_ALL_CH_ARRAY(parblock,npar,1)
  call BCAST_ALL_I(nknotA1,MAXPAR)
  call BCAST_ALL_I(nknotA2,MAXPAR)
  call BCAST_ALL_SINGLEL(unconformal)
  call BCAST_ALL_SINGLEI(ndisc)

  if(.not.allocated(kntrad))allocate(kntrad(nknotA1(1)))
  call BCAST_ALL_DP(kntrad,nknotA1(1))
  
  call BCAST_ALL_SINGLEI(mdim)

  if(.not.allocated(mdl))allocate(mdl(mdim))
  call BCAST_ALL_DP(mdl,mdim)

  end subroutine model_berkeley_broadcast

!---------------------------------------------------
subroutine model_berkeley_shsv(r,theta,phi,dvsh,dvsv,dvph,dvpv,drho,eta_aniso,iregion_code,CRUSTAL)
!---------------------------------------------------
!returns isotropic vs, vp, and rho assuming scaling dlnVs/dlnVp=2 dlnVs/dlnrho=3
!also returns anisotropic parameters xi,fi,eta,Gc,Gs,Hc,Hs,Bc,Bs if ifanis=1

    use model_berkeley_par
    use constants

    implicit none

    double precision :: x,rho1d,vpv1d,vph1d,vsv1d,vsh1d,eta1d,Qmu1d,Qkappa1d
    integer :: iregion_code
    logical :: CRUSTAL

    double precision, intent(in) :: r,theta,phi
    double precision :: vs,vp,rho,Qs
    double precision :: xi,fi,eta,Gc,Gs,eta_aniso

    integer :: jump,effnknot,i,j,k
    integer, dimension(:), allocatable :: kindex
    double precision :: lat,lon,del,dr,dv,AA,CC,FF,LL,NN,eta1,adel1,r_
    double precision, dimension(8) :: adel
    double precision, dimension(:), allocatable :: dh
    double precision :: rad2deg,getdel, spbsp
    double precision :: vsv,vsh,vpv,vph, drho,scaleval
    real(kind=4) :: dvsv,dvsh,dvpv,dvph 
    double precision :: aa1, bb1
    character (len=256) :: earth_model_file = trim(A3d_folder)//'model1D.dat'

    rad2deg = 180.d0/PI
   !  adel=(/63.4, 31.7, 15.8, 7.9, 3.95, 1.98, 0.99/)
    ! Include level 8 (FM, April 2021 - Mod. courtesy of Dan Frost)
    adel=(/63.4, 31.7, 15.8, 7.9, 3.95, 1.98, 0.99, 0.495/)
    xi = 1.d0
    fi = 1.d0

    !if(.not.model1D_exist) &
    !    stop 'no 1D model file'

    !if(ifanis_berk==1) then
    !    if(.not.(present(xi).and.present(fi).and.present(eta))) &
    !        stop 'A3d_full: ifanis_berk inconsistent'
    !endif

    ! if(r>6340.9d0) then
    !     r_=6340.9d0
    ! else
    !     r_=r
    ! endif
    
   call det_moho_radius(earth_model_file, mohodepth)
   if(r>mohodepth) then
      r_=mohodepth
   else
      r_=r
   endif
   

    x = r_/EARTH_R_KM

    call model_1dberkeley(x,rho1d,vpv1d,vph1d,vsv1d,vsh1d,eta1d,Qkappa1d,Qmu1d,iregion_code,CRUSTAL)
    if (rho1d<1200.d0)call model_1dberkeley(r_/6367.999d0,rho1d,vpv1d,vph1d,vsv1d,vsh1d,eta1d,&
                                            Qkappa1d,Qmu1d,iregion_code,CRUSTAL) ! No water in RegSEM please

    scaleval=dsqrt(PI*GRAV*EARTH_RHOAV)
    rho1d=rho1d*EARTH_RHOAV
    vpv1d=vpv1d*(EARTH_R*scaleval)
    vph1d=vph1d*(EARTH_R*scaleval)
    vsv1d=vsv1d*(EARTH_R*scaleval)
    vsh1d=vsh1d*(EARTH_R*scaleval)

    rho = rho1d
    AA  = vph1d**2*rho
    CC  = vpv1d**2*rho
    LL  = vsv1d**2*rho
    NN  = vsh1d**2*rho
    FF  = eta1d*(AA-2.d0*LL)
    qs  = qmu1d

    !call get_1Dmodel_TI(AA,CC,FF,LL,NN,rho,Qs,r_*1000.d0)
    !if (rho<1200.d0)   call get_1Dmodel_TI(AA,CC,FF,LL,NN,rho,Qs,6367999.d0)   ! No water in RegSEM please

    eta1=FF/(AA-2.*LL)
    !Voigt average
    vp=sqrt((3.*CC+(8.+4.*eta1)*AA+8.*(1.-eta1)*LL)/(15.*rho))
    vs=sqrt((CC+(1.-2.*eta1)*AA+(6.+4.*eta1)*LL+5.*NN)/(15.*rho))

    eta=eta1
    xi=NN/LL
    fi=CC/AA

    ! Vs perturbation
    if (r_>kntrad(nknotA1(1)).or.r_<kntrad(1)) then
       dv=0.0
    else
       jump=0
       allocate(dh(nknotA2(1)),kindex(nknotA2(1)))
       lat=90-rad2deg*theta
       lon=rad2deg*phi
       effnknot=0
       do i=1,nknotA2(1)
          del=getdel(lat,lon,aknot(i),oknot(i))

          if(del<=adel(level(i))*2.0) then
             effnknot=effnknot+1
             kindex(effnknot)=i
             dh(effnknot)=spbsp(del,adel(level(i)))
          endif
       enddo

       dv=0.0
       do i=1,nknotA1(1)
          call fspl(i,nknotA1(1),kntrad,r_,dr)
          do j=1,effnknot
             dv=dv+dr*dh(j)*mdl(jump+kindex(j)+nknotA2(1)*(i-1))
          enddo
       enddo
       deallocate(dh,kindex)
    endif

    ! Scaling
    vs=vs+dv*vs
    vp=vp+0.5d0*dv*vp
    rho=rho+0.33333d0*dv*rho

    
    ! Perturbation of other parameters
    if(npar<2) then   ! no other parameters
          dv=0.
       else
          do k=2,npar
             if (r_>kntrad(nknotA1(k)).or.r_<kntrad(1)) then
                dv=0.0
             else
                jump = jump + nknotA1(k-1)*nknotA2(k-1)
                allocate(dh(nknotA2(k)),kindex(nknotA2(k)))
                lat=90-rad2deg*theta
                lon=rad2deg*phi
                effnknot=0
                do i=1,nknotA2(k)
                   if(unconformal) then
                      del=getdel(lat,lon,aknot2(i),oknot2(i))
                      adel1=adel(level2(i))
                   else
                      del=getdel(lat,lon,aknot(i),oknot(i))
                      adel1=adel(level(i))
                   endif
                   if(del<=adel1*2.0) then
                      effnknot=effnknot+1
                      kindex(effnknot)=i
                      dh(effnknot)=spbsp(del,adel1)
                   endif
                enddo

                dv=0.0
                do i=1,nknotA1(k)
                   call fspl(i,nknotA1(k),kntrad,r_,dr)
                   do j=1,effnknot
                      dv=dv+dr*dh(j)*mdl(jump+kindex(j)+nknotA2(k)*(i-1))
                   enddo
                enddo
                deallocate(dh,kindex)
             endif

             if (k==2) then
                xi=xi+dv*xi
                fi=fi-1.5*dv*fi
                eta=eta-2.5*dv*eta
             elseif (k==3) then
                Gc = dv
             elseif (k==4) then
                Gs = dv
             endif
             ! Here we can add a scaling to get Hc, Hs, Bc and Bs
          enddo
       endif
       
       ! ============= commented by <FM> on Feb 3, 2020 =====
       !vsv = sqrt(3.d0/(xi+2.d0))*vs
       !vsh = sqrt(xi)*vsv
       !vph = sqrt(5.d0/(fi+4.d0))*vp
       !vpv = sqrt(fi)*vph
       !rho = rho
       ! ====================================================
       ! New conversion relationships <FM> - Feb 3, 2020
       ! Auxiliar values
       aa1 = 3.d0 + ( 8.d0 + 4.d0 * eta ) / fi
       bb1 = 1.d0 + ( 1.d0 - 2.d0 * eta ) / fi
       vsv = sqrt( 15.d0 * (vp**2 * bb1 - vs ** 2 * aa1) /&
            (8.d0 * (1.d0-eta) * bb1 - (6.d0 + 4.d0 * eta + 5.d0 * xi) * aa1 ) )
       vsh = sqrt( xi ) * vsv
       vpv = sqrt( (15.d0 * vp ** 2 - 8.d0 * (1.d0-eta) * vsv ** 2 ) /&
                   (3.d0 + ( 8.d0 + 4.d0 * eta ) / fi ) )
       vph = vpv * sqrt(1.d0/fi)
       rho = rho
       ! ===================================================

   !
   dvsv = vsv/vsv1d-1.d0
   if(vsv1d==0)dvsv = 0.d0
   dvsh = vsh/vsh1d-1.d0
   if(vsh1d==0)dvsh = 0.d0
   dvpv = vpv/vpv1d-1.d0
   dvph = vph/vph1d-1.d0
   drho  = rho/rho1d-1.d0
   eta_aniso = eta
   !
!---------------------------------------------------
end subroutine model_berkeley_shsv
!---------------------------------------------------

!---------------------------------------------------
doubleprecision function getdel(a0,o0,a,o)
!---------------------------------------------------

    doubleprecision :: a0,o0,a,o
    doubleprecision :: q0,sq0,cq0,q,sq,cq,ff,sff,cff,arg    
    double precision :: deg2rad,rad2deg

    deg2rad = (4.d0*datan(1.d0))/180.d0
    rad2deg = 180.d0/(4.d0*datan(1.d0))

    q0=(90.0-a0)*deg2rad;
    sq0=sin(q0);
    cq0=cos(q0);

    q=(90.0-a)*deg2rad;
    sq=sin(q);
    cq=cos(q);

    ff=(o-o0)*deg2rad;
    sff=sin(ff);
    cff=cos(ff);

    arg=cq*cq0+sq*sq0*cff;
    if(arg > 1.) arg= 1.;
    if(arg < -1.) arg=-1.;
    getdel=rad2deg*acos(arg);

!---------------------------------------------------
end function getdel
!---------------------------------------------------


!---------------------------------------------------
doubleprecision function spbsp(hdel,ahdel)
!---------------------------------------------------

    doubleprecision :: hdel,ahdel

    if(hdel < ahdel) then
       spbsp=0.75*(hdel/ahdel)*(hdel/ahdel)*(hdel/ahdel)-1.5*(hdel/ahdel)*(hdel/ahdel)+1.
    else if(hdel <= ahdel*2.) then
       spbsp=0.25*(2.-hdel/ahdel)*(2.-hdel/ahdel)*(2.-hdel/ahdel)
       if(spbsp < 0) spbsp=0.0
    else
       spbsp=0.0
    endif

!---------------------------------------------------
end function spbsp
!---------------------------------------------------






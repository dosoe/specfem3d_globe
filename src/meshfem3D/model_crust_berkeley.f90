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

!---------------------------------------!
! berkeley crust implemented as in csem !
!---------------------------------------!
!
!==============================
module model_crust_berkeley_par
!==============================
  !
  implicit none
  !
  integer :: NBT,NBP
  real :: drin
  real :: drfiltre = 2.e0
  doubleprecision, parameter :: dr_ = 2.d0
  real, dimension(:,:), allocatable :: crust_array
  real, dimension(:,:)  , allocatable :: moho_start
  doubleprecision :: mohodepth = -1.0d0
  !
  !! Utpal Kumar, 2022-03-17
  contains
  !! subroutine to read earth model
  subroutine det_moho_depth(model_file, mohodp)
   implicit none
 
   !! Declare vars
   integer :: FID = 10
   character (len=100), intent(in) :: model_file
   character (len=100) :: CTMP
   doubleprecision, allocatable :: radius(:),density(:),vpv(:),vsv(:),qkappa(:),qmu(:),vph(:),vsh(:),eta(:)
   doubleprecision, allocatable :: derivdensity(:), derivVp(:), derivVs(:)
   doubleprecision :: tol = 2.0d0**(-5), maxderivdensity = 0., earthradius = 6371.
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
       mohodp = earthradius - mohodp !get moho depth
   
   end if
 end subroutine det_moho_depth
!==================================
end module model_crust_berkeley_par
!==================================
!
!================================================
subroutine model_berkeley_crust_broadcast(myrank)
!================================================
  !
  ! standard routine to setup model
  !
  use model_crust_berkeley_par
  use constants, only : A3d_folder
  !
  implicit none
  !
  integer :: myrank
  character(len=100), parameter :: file_crust = trim(A3d_folder)//'crust2cru2av_2x2.dat'
  character(len=100), parameter :: file_moho = trim(A3d_folder)//'crust2moho_2x2.dat'
  !
  print*,'I am initializing SMOOTH crust!',myrank
  open(52,file=file_moho,status='old',action='read')
  call read_crustmoho_filtre(52)
  drfiltre=dr_
  close(52) 
  open(52,file=file_crust,status='old',action='read')
  call read_crust_smooth_variable(52)
  close(52)   
  !
!============================================
end subroutine model_berkeley_crust_broadcast
!============================================
!
!========================================================
subroutine model_berkeley_crust(lat,lon,x,vp,vs,rho,moho)
!========================================================
  !
  use model_crust_berkeley_par
  use constants
  !
  implicit none
  !
  double precision,intent(in) :: lat,lon,x
  double precision,intent(out) :: vp,vs,rho,moho
  !
  double precision :: vsv,vsh
  !
  doubleprecision :: depth, theta, phi
  double precision :: moho_depth, scaleval
  real, parameter :: deg2rad=pi/180.,rad2deg=180./pi
  !
  theta = (90-lat)*deg2rad ! assumed lat range: [-90,90]
  phi = lon*deg2rad ! assumed lon range: [-180,180]
  !
  depth = (1-x)*EARTH_R_KM
  !
  call get_crust_val_csem(theta,phi,depth,rho,vp,vsv,vsh,moho_depth)
  !
  ! get equivalent isotropic vs
  !
  vs = dsqrt((2.d0*vsv**2.d0+vsh**2.d0)/3.d0)
  !
  ! scale values for specfem
  !
  scaleval=dsqrt(PI*GRAV*EARTH_RHOAV)
  vp=vp/(EARTH_R*scaleval)
  vs=vs/(EARTH_R*scaleval)
  rho=rho/EARTH_RHOAV
  moho = moho_depth/EARTH_R_KM
  !
!==================================
end subroutine model_berkeley_crust
!==================================
!
!==================================================================================
subroutine model_berkeley_crust_aniso(lat,lon,x,vpv,vph,vsv,vsh,eta_aniso,rho,moho)
!==================================================================================
  !
  use model_crust_berkeley_par
  use constants
  !
  implicit none
  !
  double precision,intent(in) :: lat,lon,x
  double precision,intent(out) :: vpv,vph,vsv,vsh,rho,moho
  !
  double precision :: vp,eta_aniso
  !
  doubleprecision :: depth, theta,phi
  double precision :: moho_depth, scaleval
  real, parameter :: deg2rad=pi/180.,rad2deg=180./pi
  !
  theta = (90-lat)*deg2rad ! assumed lat range: [-90,90]
  phi = lon*deg2rad ! assumed lon range: [-180,180]
  !
  depth = (1-x)*EARTH_R_KM
  !
  call get_crust_val_csem(theta,phi,depth,rho,vp,vsv,vsh,moho_depth)
  !
  ! scale values for specfem
  !
  scaleval=dsqrt(PI*GRAV*EARTH_RHOAV)
  vph=vp/(EARTH_R*scaleval)
  vpv=vp/(EARTH_R*scaleval)
  vsh=vsh/(EARTH_R*scaleval)
  vsv=vsv/(EARTH_R*scaleval)
  rho=rho/EARTH_RHOAV
  moho = moho_depth/EARTH_R_KM
  eta_aniso = 1.d0
  !
!========================================
end subroutine model_berkeley_crust_aniso
!========================================
!
!===================================================================
subroutine get_crust_val_csem(theta,phi,z,rho,vp,vsv,vsh,moho_depth)
!===================================================================
!
use model_crust_berkeley_par
use constants, only : A3d_folder
implicit none
!
double precision :: xi(5),wi(5),lagrange
double precision :: theta,phi,z
double precision :: rho,vp,vsv,vsh,moho_depth
double precision :: rho_cr(5),vp_cr(5),vsv_cr(5),vsh_cr(5)
double precision :: x
real :: moho_filtre
integer :: i
! double precision, parameter :: depth_moho_berkeley = 30.d0
! double precision, parameter :: depth_moho_berkeley = 24.4d0
double precision :: depth_moho_berkeley
character (len=256) :: earth_model_file = trim(A3d_folder)//'model1D.dat'
!
! Utpal Kumar, 2022-03-17
call det_moho_depth(earth_model_file, mohodepth)
depth_moho_berkeley = mohodepth
! write(*,*) "earth_model_file, depth_moho_berkeley", trim(earth_model_file), depth_moho_berkeley
! get moho depth
!
! moho_depth0 = 30.d0-dble(moho_filtre(theta,phi))
moho_depth = depth_moho_berkeley-dble(moho_filtre(theta,phi))
! write(*,*) "Moho depth: ", moho_depth, moho_depth0
!
! horizontal interpolation for all registered depths
! 
call crust_bilinear_variable(theta,phi,1,rho_cr(1),vp_cr(1),vsv_cr(1),vsh_cr(1))
call crust_bilinear_variable(theta,phi,2,rho_cr(2),vp_cr(2),vsv_cr(2),vsh_cr(2))
call crust_bilinear_variable(theta,phi,3,rho_cr(3),vp_cr(3),vsv_cr(3),vsh_cr(3))
call crust_bilinear_variable(theta,phi,4,rho_cr(4),vp_cr(4),vsv_cr(4),vsh_cr(4))
call crust_bilinear_variable(theta,phi,5,rho_cr(5),vp_cr(5),vsv_cr(5),vsh_cr(5))
!
! get gll nodes position
!
call LGL_NODES(4,1.d-6,XI,WI)
xi = moho_depth-(xi+1.d0)/2.d0*moho_depth
x = z
if(x>maxval(xi))x = maxval(xi)
if(x<minval(xi))x = minval(xi)
!
! init values
!
vp  = 0.d0
vsv = 0.d0
vsh = 0.d0
rho = 0.d0
!
! depth interpolation
!
do i = 1,5
   rho = rho+lagrange(i,xi,5,x)*rho_cr(i)
   vp  = vp +lagrange(i,xi,5,x)*vp_cr(i)
   vsv = vsv+lagrange(i,xi,5,x)*vsv_cr(i)
   vsh = vsh+lagrange(i,xi,5,x)*vsh_cr(i)
enddo
!
rho = rho * 1000.d0
vp  = vp  * 1000.d0
vsv = vsv * 1000.d0
vsh = vsh * 1000.d0
!
return
!
!================================
end subroutine get_crust_val_csem
!================================
!
!===============================================
double precision function lagrange(ind,xi,nxi,x)
!===============================================
!
implicit none
!
integer :: ind,nxi,i
double precision :: xi(nxi),x
!      
lagrange=1.0
!
do i = 1,nxi
   if(i.ne.ind)lagrange=lagrange*((x-xi(i))/(xi(ind)-xi(i)))
enddo
!
return
!
!====================
end function lagrange
!====================
!
!================================
SUBROUTINE LGL_NODES(N,EPS,XI,WI)
!================================
!
IMPLICIT NONE
!
INTEGER :: N ! POLYNOMIAL ORDER
DOUBLE PRECISION   :: EPS ! DESIRED ERROR
!
DOUBLE PRECISION   :: XI(N+1) ! GAUSS-LOBATTO-LEGENDRE INTERPOLATION POINTS
DOUBLE PRECISION   :: WI(N+1) ! GAUSS-LOBATTO-LEGENDRE QUADRATURE WEIGHT
!      
DOUBLE PRECISION   :: P(N+1,N+1),PI,XOLD(N+1)
DOUBLE PRECISION   :: ERMAX
INTEGER :: I,J,N1
!
PI=4.D0*DATAN(1.D0)
N1=N+1
!
DO I = 1,N1
   XI(I)=-COS(PI*real(I-1)/real(N1-1))
ENDDO
!
ERMAX=2.*EPS
!
DO WHILE(ERMAX.GT.EPS)
   DO I = 1,N1
      XOLD(I)=XI(I)
      P(I,1)=1.0
      P(I,2)=XI(I)
   ENDDO
   DO J = 2,N
      DO I = 1,N1
         P(I,J+1)=((2*J-1)*XI(I)*P(I,J)-(J-1)*P(I,J-1))/REAL(J)
      ENDDO
   ENDDO
   DO I = 1,N1
      XI(I)=XOLD(I)-(XI(I)*P(I,N1)-P(I,N))/(N1*P(I,N1))
   ENDDO
   ERMAX=0.0
   DO I = 1,N1
      ERMAX=MAX(ERMAX,ABS(XI(I)-XOLD(I)))
   ENDDO
ENDDO
!     
DO I = 1,N1
   WI(I)=2.0/(N*N1*P(I,N1)**2.0)
ENDDO
!
RETURN
!
!=======================
END SUBROUTINE LGL_NODES
!=======================

!-----------------------------------------------------------------
  real function moho_filtre(theta,phi)
!theta phi en radians
!dr en degre
!reponse en metre
!-----------------------------------------------------------------
    use model_crust_berkeley_par
    implicit none
    doubleprecision :: theta,phi
    real::t,p
    real, parameter :: pi=3.141592653589793,deg2rad=pi/180.,rad2deg=180./pi
    t=theta/deg2rad
    p=phi  /deg2rad
    moho_filtre=gauss_filtre1(moho_start,t,p,drfiltre)
contains

!-----------------------------------------------------------------
  function gauss_filtre1(tin,theta,phi,dr)
!-----------------------------------------------------------------
    real :: gauss_filtre1
    real, dimension(:,:)  :: tin
    real :: theta,phi,dr,thetar,phir,tmpnorm,inte
    real :: tmp
    integer :: i,ii,j,jj,LARG
    real, parameter :: pi=3.141592653589793,deg2rad=pi/180.,rad2deg=180./pi
    !
    tmp=0.
    tmpnorm=0. 
    LARG=10
    do i=1,LARG+1
       do j=1,LARG+1
          call get_indexloc(phi,theta,i,j,dr,LARG,ii,jj,phir,thetar)
          inte=cos_cylindre(theta,phi,dr,thetar,phir)*(dr/real(LARG/2)*deg2rad)**2*sin(thetar*deg2rad)
          tmp=tmp+tin(ii,jj)*inte
          tmpnorm=tmpnorm+inte
       enddo
    enddo
    gauss_filtre1=tmp/tmpnorm

!-----------------------------------------------------------------
  end function gauss_filtre1
!-----------------------------------------------------------------


!----------------------------------------------------------------------
  real function cos_cylindre(t0_,p0_,d0_,theta_,phi_)
!----------------------------------------------------------------------
    implicit none
    real :: t0,p0,d0,theta,phi, d_ang
    real :: t0_,p0_,d0_,theta_,phi_,cosa
!
    t0=t0_*deg2rad
    p0=p0_*deg2rad
    theta=theta_*deg2rad
    phi=phi_*deg2rad
    d0=d0_*deg2rad
!distance angulaire au centre du cylindre:
    cosa=cos(theta)*cos(t0)+sin(theta)*sin(t0)*cos(phi-p0)
    if (cosa >= 1. ) then
       d_ang=0.
    else if (cosa <= -1. ) then
       d_ang=4.*atan(1.)
    else
       d_ang=acos(cos(theta)*cos(t0)+sin(theta)*sin(t0)*cos(phi-p0))
    endif
    if (d_ang>d0) then
       cos_cylindre=0.d0
    else
       cos_cylindre=0.5d0*(1.d0+cos(PI*d_ang/d0))
    endif
!----------------------------------------------------------------------
  end function cos_cylindre
!----------------------------------------------------------------------

!-----------------------------------------------------------------
  end function moho_filtre
!-----------------------------------------------------------------

!-----------------------------------------------------------------
  subroutine get_indexloc(phi,theta,i,j,dr,LARG,ii,jj,phir,thetar)
!-----------------------------------------------------------------
    
    use model_crust_berkeley_par

    implicit none
    real, intent(in) :: theta,phi,dr
    integer, intent(in) :: i,j,LARG
    doubleprecision :: t,p,eps
    real, intent(out) :: thetar,phir
    integer, intent(out) :: ii,jj
    eps=1.d-8
!       
    p  =phi+(i-1-LARG/2)*dr/real(LARG/2)
    t  =theta+(j-1-LARG/2)*dr/real(LARG/2)
    if (p<0.d0-eps) p=p+360.d0
    if (p>=360.d0-eps) p=p-360.d0
    if (t>180.d0-eps) then
       t=t-180.d0
       p=360.d0-p
    else if (t< 0.d0-eps) then
       t=180.d0+t
       p=360.d0-p
    endif
    if (p<0.d0-eps) p=p+360.d0
    if (p>=360.d0-eps) p=p-360.d0
!
    ii=nint(p/drin)+1
    if (ii>NBP) ii=NBP
    jj=nint(t/drin)+1
    if (jj>NBT) jj=NBT
    thetar=t
    phir  =p
!-----------------------------------------------------------------
  end subroutine get_indexloc
!-----------------------------------------------------------------

!-----------------------------------------------------------------
  subroutine crust_bilinear_variable(theta,phi,gll_ind,rho_cr,vp_cr,vsv_cr,vsh_cr)
!theta phi en radians
!dr en degre
!reponse en metre
!-----------------------------------------------------------------
    use model_crust_berkeley_par
    implicit none
!    real, dimension(2) :: crust_array
    integer :: gll_ind
    doubleprecision :: theta,phi,dx,dy,t,p,D_DEG
    doubleprecision, intent(out) :: rho_cr,vp_cr,vsv_cr,vsh_cr
    integer :: broj_lokacija,j,m
    real, parameter :: pi=3.141592653589793,deg2rad=pi/180.,rad2deg=180./pi

    ! sfrench 20110103 model grid spacing (degrees)
    D_DEG = 2.d0

    t = 90 - theta/deg2rad

    ! sfrench 20110103 new bounds on latitude in accord with new model grid spacing
    !if(t>89.d0)  t = 89.d0
    !if(t<-89.d0) t = -89.d0
    if(t>90.d0 - D_DEG) t = 90.d0 - D_DEG
    if(t<D_DEG - 90.d0) t = D_DEG - 90.d0

    p = phi/deg2rad
    if(p>180.d0) p = p - 360.d0
    broj_lokacija = size(crust_array, DIM = 1)

    m = 0

    rho_cr = 0.d0 
    vp_cr = 0.d0 
    vsv_cr = 0.d0
    vsh_cr = 0.d0

    do j = 1,broj_lokacija

       if(m>=4) then
          exit
       else

          if( gll_ind == nint(crust_array(j,3) + 1) ) then

             dy = dabs(crust_array(j,1)-t)        

             dy = dy / D_DEG ! sfrench 20110103 : normalized to model grid spacing

             if(dy<1.d0) then
             
                dx = dabs(crust_array(j,2)-p)

                if (dx>180.0) dx = 360.0-dx

                dx = dx / D_DEG ! sfrench 20110103 : normalized to model grid spacing

                if(dabs(dx)<1.d0) then
                   ! Increment number of found locations. Must be <= 4 
                   m = m + 1
                   rho_cr = rho_cr + (1.d0-dx)*(1.d0-dy)*crust_array(j,4) 
                   vp_cr  = vp_cr  + (1.d0-dx)*(1.d0-dy)*crust_array(j,5) 
                   vsv_cr = vsv_cr + (1.d0-dx)*(1.d0-dy)*crust_array(j,6) 
                   vsh_cr = vsh_cr + (1.d0-dx)*(1.d0-dy)*crust_array(j,7)
                  !  print*,'m= ',m,'dx= ',dx,' dy= ',dy, 'Vsv= ',vsv_cr


                endif
  
             endif
             
          endif
       
       endif

    enddo
    
!-----------------------------------------------------------------
  end subroutine crust_bilinear_variable
!-----------------------------------------------------------------

!-----------------------------------------------------------------
  subroutine read_crust_smooth_variable(unit)
!-----------------------------------------------------------------
    use model_crust_berkeley_par
    implicit none
    integer, intent(in) :: unit
    integer :: j,l,broj_lokacija
!
    
    read(unit,*) broj_lokacija

    ! sfrench 20110103 note: crust_array is <t> <p> <gll_ind> <rho> <vp> <vsv> <vsh>
    allocate(crust_array(broj_lokacija,7))
    crust_array(:,:)=-1000.
    do j=1,broj_lokacija
       read(unit,*) (crust_array(j,l),l=1,7)
    enddo
    print*,'read_crust_smooth_variable: I have read ',broj_lokacija,' crustal inputs!'
    
!-----------------------------------------------------------------
  end subroutine read_crust_smooth_variable
!-----------------------------------------------------------------

!-----------------------------------------------------------------
  subroutine read_crustmoho_filtre(unit)
!-----------------------------------------------------------------
    use model_crust_berkeley_par
    implicit none
    integer, intent(in) :: unit
    integer :: i,j
!
    
    read(unit,*) NBP,NBT,drin
    if (drin/=2.) STOP 'read_crust_filtre: dr muste be ==2'    
    allocate(moho_start(NBP,NBT))
   moho_start(:,:)=-1000.
    do j=1,NBP
       do i=1,NBT
          read(unit,*) moho_start(j,i)
       enddo
    enddo
    moho_start(:,:)=moho_start(:,:)/1000.
    NBT=NBT-1 
!-----------------------------------------------------------------
  end subroutine read_crustmoho_filtre
!-----------------------------------------------------------------

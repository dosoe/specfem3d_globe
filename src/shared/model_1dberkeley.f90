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

!-------------------------------
!
! 1D Berkeley model 
!
! Add infos... 
!
!-------------------------------

module model_1dberkeley_par
  ! Added by <FM> Feb. 2022
  use constants, only: A3d_folder

  ! number of layers in model1Dberkeley.dat
  integer :: NR_REF_BERKELEY
  integer :: NR_inner_core_berk
  integer :: NR_outer_core_berk
  integer :: NR_water_berk
  integer :: ifanis_berk
  integer :: tref_berk
  integer :: ifdeck_berk

  ! model_1dberkeley_variables
  double precision, dimension(:), allocatable :: &
    Mref_V_radius_berkeley,                      &
    Mref_V_density_berkeley,                     &
    Mref_V_vpv_berkeley,                         &
    Mref_V_vph_berkeley,                         &
    Mref_V_vsv_berkeley,                         &
    Mref_V_vsh_berkeley,                         &
    Mref_V_eta_berkeley,                         &
    Mref_V_Qkappa_berkeley,                      &
    Mref_V_Qmu_berkeley

  ! Utpal Kumar, Feb, 2022
  ! define the berkeley 1D model 
  character (len=100) :: berk_model1D = trim(A3d_folder)//'model1D.dat'
  integer :: modemohoberk = -1

end module model_1dberkeley_par


  !--------------------------------------------
  subroutine model_1dberkeley_broadcast()
  !--------------------------------------------
  !
  ! reads and broadcasts berkeley 1D model
  !
  use constants
  use model_1dberkeley_par
  !
  implicit none
  !
  integer :: i, icode !, ier
  integer, parameter :: lunit = 54
  character (len=100) :: filename, title
  !
  ! define the berkeley 1D model 
  ! Utpal Kumar, Feb, 2022
  filename = berk_model1D
  !
  ! root nodes read header
  !
  if(myrank==0)then
    !
    open(lunit,file=trim(filename),status='old')
    !
    read(lunit,100,iostat=icode) title
    read(lunit,*  ,iostat=icode) ifanis_berk,&
                                 tref_berk,  &
                                 ifdeck_berk
    !                             
    read(lunit,*  ,iostat=icode) NR_REF_BERKELEY,    &
                                 NR_inner_core_berk, &
                                 NR_outer_core_berk, &
                                 NR_water_berk
    !
  endif
  !
  ! broadcast header values
  !
  !call MPI_BCAST(ifanis_berk       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(tref_berk         ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(ifdeck_berk       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(NR_REF_BERKELEY   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(NR_inner_core_berk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(NR_outer_core_berk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !call MPI_BCAST(NR_water_berk     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !
  call BCAST_ALL_SINGLEI(ifanis_berk       )
  call BCAST_ALL_SINGLEI(tref_berk         )
  call BCAST_ALL_SINGLEI(ifdeck_berk       )
  call BCAST_ALL_SINGLEI(NR_REF_BERKELEY   )
  call BCAST_ALL_SINGLEI(NR_inner_core_berk)
  call BCAST_ALL_SINGLEI(NR_outer_core_berk)
  call BCAST_ALL_SINGLEI(NR_water_berk     )
  !
  ! allocate arrays
  !
  allocate(Mref_V_radius_berkeley(NR_REF_BERKELEY),  &
           Mref_V_density_berkeley(NR_REF_BERKELEY), &
           Mref_V_vpv_berkeley(NR_REF_BERKELEY),     &
           Mref_V_vsv_berkeley(NR_REF_BERKELEY),     &
           Mref_V_Qkappa_berkeley(NR_REF_BERKELEY),  &
           Mref_V_Qmu_berkeley(NR_REF_BERKELEY),     &
           Mref_V_vph_berkeley(NR_REF_BERKELEY),     &
           Mref_V_vsh_berkeley(NR_REF_BERKELEY),     &
           Mref_V_eta_berkeley(NR_REF_BERKELEY)      )
  !
  ! root proc reads data
  !
  if(myrank==0)then
  !
  do i = 1,NR_REF_BERKELEY
    read(lunit,*)Mref_V_radius_berkeley(i),  &
                 Mref_V_density_berkeley(i), &
                 Mref_V_vpv_berkeley(i),     &
                 Mref_V_vsv_berkeley(i),     &
                 Mref_V_Qkappa_berkeley(i),  &
                 Mref_V_Qmu_berkeley(i),     &
                 Mref_V_vph_berkeley(i),     &
                 Mref_V_vsh_berkeley(i),     &
                 Mref_V_eta_berkeley(i)     
  enddo
  !
  close(lunit)
  !
  endif
  !
  ! broadcast data
  !
  call BCAST_ALL_DP(Mref_V_radius_berkeley ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_density_berkeley,NR_REF_BERKELEY) 
  call BCAST_ALL_DP(Mref_V_vpv_berkeley    ,NR_REF_BERKELEY)     
  call BCAST_ALL_DP(Mref_V_vph_berkeley    ,NR_REF_BERKELEY)     
  call BCAST_ALL_DP(Mref_V_vsv_berkeley    ,NR_REF_BERKELEY)     
  call BCAST_ALL_DP(Mref_V_vsh_berkeley    ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_eta_berkeley    ,NR_REF_BERKELEY)  
  call BCAST_ALL_DP(Mref_V_Qkappa_berkeley ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_Qmu_berkeley    ,NR_REF_BERKELEY)
  !
  ! reading formats
  !
100 format(a80)
105 format(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)   
  !
  !----------------------------------------
  end subroutine model_1dberkeley_broadcast
  !----------------------------------------
  



!
!----------------------------------------------
!

  subroutine model_1dberkeley(x,rho,vpv,vph,vsv,vsh,eta,Qkappa,Qmu,iregion_code,CRUSTAL)

  use constants
  use model_1dberkeley_par

  implicit none

! model_1dref_variables

! input:
! dimensionless radius x

! output: non-dimensionalized
!
! mass density             : rho
! compressional wave speed : vpv
! compressional wave speed : vph
! shear wave speed         : vsv
! shear wave speed         : vsh
! dimensionless parameter  : eta
! shear quality factor     : Qmu
! bulk quality factor      : Qkappa

  double precision :: x,rho,vpv,vph,vsv,vsh,eta,Qmu,Qkappa
  integer :: iregion_code
  logical :: CRUSTAL

  ! local parameters
  double precision :: r,frac,scaleval
  integer :: i, mohonodeval
  logical, parameter :: mimic_native_specfem = .true. 

  ! compute real physical radius in meters
  r = x * EARTH_R

  i = 1
  do while(r >= Mref_V_radius_berkeley(i) .and. i /= NR_REF_BERKELEY)
    i = i + 1
  enddo

! make sure we stay in the right region
  if(mimic_native_specfem .and. iregion_code == IREGION_INNER_CORE .and. i > NR_inner_core_berk) i = NR_inner_core_berk

  if(mimic_native_specfem .and. iregion_code == IREGION_OUTER_CORE .and. i < NR_inner_core_berk+2) i = NR_inner_core_berk+2
  if(mimic_native_specfem .and. iregion_code == IREGION_OUTER_CORE .and. i > NR_outer_core_berk) i = NR_outer_core_berk

  if(mimic_native_specfem .and. iregion_code == IREGION_CRUST_MANTLE .and. i < NR_outer_core_berk+2) i = NR_outer_core_berk+2

  ! if crustal model is used, mantle gets expanded up to surface
  ! for any depth less than 24.4 km, values from mantle below moho are taken
  !! Utpal Kumar, Feb, 2022
  if (modemohoberk < 0) then
    call est_moho_node(mohonodeval)
    modemohoberk = mohonodeval
  end if

  if(mimic_native_specfem .and. CRUSTAL .and. i > modemohoberk) then
      i = modemohoberk ! Warining : may need to be changed if file is modified !
  end if
  !
  if(i == 1) then
    ! first layer in inner core
    rho    = Mref_V_density_berkeley(i)
    vpv    = Mref_V_vpv_berkeley(i)
    vph    = Mref_V_vph_berkeley(i)
    vsv    = Mref_V_vsv_berkeley(i)
    vsh    = Mref_V_vsh_berkeley(i)
    eta    = Mref_V_eta_berkeley(i)
    Qkappa = Mref_V_Qkappa_berkeley(i)
    Qmu    = Mref_V_Qmu_berkeley(i)
  else
    ! interpolates between one layer below to actual radius layer,
    ! that is from radius_ref(i-1) to r using the values at i-1 and i
    frac = (r-Mref_V_radius_berkeley(i-1))/(Mref_V_radius_berkeley(i)-Mref_V_radius_berkeley(i-1))
    ! interpolated model parameters
    rho = Mref_V_density_berkeley(i-1)  + frac * (Mref_V_density_berkeley(i)- Mref_V_density_berkeley(i-1))
    vpv = Mref_V_vpv_berkeley(i-1)      + frac * (Mref_V_vpv_berkeley(i)    - Mref_V_vpv_berkeley(i-1)    )
    vph = Mref_V_vph_berkeley(i-1)      + frac * (Mref_V_vph_berkeley(i)    - Mref_V_vph_berkeley(i-1)    )
    vsv = Mref_V_vsv_berkeley(i-1)      + frac * (Mref_V_vsv_berkeley(i)    - Mref_V_vsv_berkeley(i-1)    )
    vsh = Mref_V_vsh_berkeley(i-1)      + frac * (Mref_V_vsh_berkeley(i)    - Mref_V_vsh_berkeley(i-1)    )
    eta = Mref_V_eta_berkeley(i-1)      + frac * (Mref_V_eta_berkeley(i)    - Mref_V_eta_berkeley(i-1)    )
    Qkappa = Mref_V_Qkappa_berkeley(i-1)+ frac * (Mref_V_Qkappa_berkeley(i) - Mref_V_Qkappa_berkeley(i-1) )
    Qmu = Mref_V_Qmu_berkeley(i-1)      + frac * (Mref_V_Qmu_berkeley(i)    - Mref_V_Qmu_berkeley(i-1)    )
  endif

  ! make sure Vs is zero in the outer core even if roundoff errors on depth
  ! also set fictitious attenuation to a very high value (attenuation is not used in the fluid)
  if(mimic_native_specfem .and. iregion_code == IREGION_OUTER_CORE) then
    vsv = 0.d0
    vsh = 0.d0
    Qkappa = 3000.d0
    Qmu = 3000.d0
  endif

  ! non-dimensionalize
  ! time scaling (s^{-1}) is done with scaleval
  scaleval=dsqrt(PI*GRAV*EARTH_RHOAV)
  rho=rho/EARTH_RHOAV
  vpv=vpv/(EARTH_R*scaleval)
  vph=vph/(EARTH_R*scaleval)
  vsv=vsv/(EARTH_R*scaleval)
  vsh=vsh/(EARTH_R*scaleval)

  end subroutine model_1dberkeley

  !! Utpal Kumar, Feb, 2022
  !! Subroutine to decide whether the moho node has already been computed
  subroutine est_moho_node(estmohonode)
    use model_1dberkeley_par
    implicit none
    integer :: estmohonode
    
    if (modemohoberk < 0) then
      call det_moho_node(berk_model1D, estmohonode)
      ! print*,"Determining Moho node ",estmohonode
    else
      estmohonode = modemohoberk
    end if
    
    end subroutine est_moho_node
  

!! Utpal Kumar, Feb 2022
!! subroutine to get the moho node number
subroutine det_moho_node(model_file, mohonodebrk)
  implicit none


  !! Declare vars
  integer :: FID = 10
  character (len=100), intent(in) :: model_file
  character (len=100) :: CTMP
  real (kind=8), allocatable :: radius(:),density(:),vpv(:),vsv(:),qkappa(:),qmu(:),vph(:),vsh(:),eta(:)
  real (kind=8), allocatable :: derivdensity(:), derivVp(:), derivVs(:)
  real (kind=8) :: tol = 2.0d0**(-5), maxderivdensity = 0.
  integer :: i = 0, IERR = 0
  integer :: totalheadersval = 3
  integer, intent(inout) :: mohonodebrk
  integer :: num_lines, totalDiscont
  integer :: j

  ! model_file = berk_model1D
  mohonodebrk = 1 !initialize moho node
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

  j = 1
  do i = 1, num_lines-1
      if (abs(radius(i+1)-radius(i)) < tol) then
          if ((abs(vsv(i)) > tol) .and. (abs(vsv(i+1)) > tol) .and. (abs(derivVs(i)) > tol) &
          .and. (abs(radius(i)-6371000.) < 90000.)) then
              
              if (abs(derivdensity(i)) > maxderivdensity) then
                  maxderivdensity = abs(derivdensity(i))
                  mohonodebrk = i
              end if
          end if
          j = j + 1
      end if
  end do
  
  close(FID)
  deallocate(radius,density,vpv,vsv,qkappa,qmu,vph,vsh,eta)
  deallocate(derivdensity, derivVp, derivVs)


end subroutine det_moho_node

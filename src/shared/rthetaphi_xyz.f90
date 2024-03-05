!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

  subroutine xyz_2_rthetaphi(x,y,z,r,theta,phi)

! convert x y z to r theta phi, single precision call

  use constants, only: CUSTOM_REAL,SMALL_VAL_ANGLE,ZERO

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: x,y,z
  real(kind=CUSTOM_REAL), intent(out) :: r,theta,phi

  double precision xmesh,ymesh,zmesh

  xmesh = dble(x)
  ymesh = dble(y)
  zmesh = dble(z)

  if (zmesh > -SMALL_VAL_ANGLE .and. zmesh <= ZERO) zmesh = -SMALL_VAL_ANGLE
  if (zmesh < SMALL_VAL_ANGLE .and. zmesh >= ZERO) zmesh = SMALL_VAL_ANGLE

  theta = real(datan2(dsqrt(xmesh*xmesh+ymesh*ymesh),zmesh), kind=CUSTOM_REAL)

  if (xmesh > -SMALL_VAL_ANGLE .and. xmesh <= ZERO) xmesh = -SMALL_VAL_ANGLE
  if (xmesh < SMALL_VAL_ANGLE .and. xmesh >= ZERO) xmesh = SMALL_VAL_ANGLE

  phi = real(datan2(ymesh,xmesh), kind=CUSTOM_REAL)

  r = real(dsqrt(xmesh*xmesh + ymesh*ymesh + zmesh*zmesh), kind=CUSTOM_REAL)

  end subroutine xyz_2_rthetaphi

!
!------------------------------------------------------------------------------
!

  subroutine xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

! convert x y z to r theta phi, double precision call

  use constants, only: SMALL_VAL_ANGLE,ZERO

  implicit none

  double precision, intent(in) :: x,y,z
  double precision, intent(out) :: r,theta,phi

  double precision :: xmesh,ymesh,zmesh

  xmesh = x
  ymesh = y
  zmesh = z

  if (zmesh > -SMALL_VAL_ANGLE .and. zmesh <= ZERO) zmesh = -SMALL_VAL_ANGLE
  if (zmesh < SMALL_VAL_ANGLE .and. zmesh >= ZERO) zmesh = SMALL_VAL_ANGLE

  theta = datan2(dsqrt(xmesh*xmesh+ymesh*ymesh),zmesh)

  if (xmesh > -SMALL_VAL_ANGLE .and. xmesh <= ZERO) xmesh = -SMALL_VAL_ANGLE
  if (xmesh < SMALL_VAL_ANGLE .and. xmesh >= ZERO) xmesh = SMALL_VAL_ANGLE

  phi = datan2(ymesh,xmesh)

  r = dsqrt(xmesh*xmesh + ymesh*ymesh + zmesh*zmesh)

  end subroutine xyz_2_rthetaphi_dble

!
!------------------------------------------------------------------------------
!

  subroutine rthetaphi_2_xyz(x,y,z,r,theta,phi)

! convert r theta phi to x y z

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: r,theta,phi
  real(kind=CUSTOM_REAL),intent(out) :: x,y,z

  x = r * sin(theta) * cos(phi)
  y = r * sin(theta) * sin(phi)
  z = r * cos(theta)

  end subroutine rthetaphi_2_xyz


!
!------------------------------------------------------------------------------
!

  subroutine rthetaphi_2_xyz_dble(x,y,z,r,theta,phi)

! convert r theta phi to x y z

  implicit none

  double precision,intent(in) :: r,theta,phi
  double precision,intent(out) :: x,y,z

  x = r * sin(theta) * cos(phi)
  y = r * sin(theta) * sin(phi)
  z = r * cos(theta)

  end subroutine rthetaphi_2_xyz_dble

!
!------------------------------------------------------------------------------
!

  subroutine latlon_2_thetaphi_dble(lat,lon,theta,phi)

! converts lat/lon (in degree) to theta/phi (in rad), double precision call

  use constants, only: DEGREES_TO_RADIANS

  implicit none

  double precision, intent(in) :: lat,lon
  double precision, intent(out) :: theta,phi

  ! converts to theta/phi
  theta = (90.d0 - lat) * DEGREES_TO_RADIANS
  phi = lon * DEGREES_TO_RADIANS

  ! reduces range for colatitude to 0 and PI, for longitude to 0 and 2*PI
  call reduce(theta,phi)

  end subroutine latlon_2_thetaphi_dble

!
!------------------------------------------------------------------------------
!

  subroutine latlon_2_geocentric_thetaphi_dble(lat,lon,theta,phi,is_mesh_elliptical)

! returns latitude/longitude in radians for geocentric location

  use constants, only: DEGREES_TO_RADIANS

  implicit none

  double precision,intent(in) :: lat,lon
  double precision,intent(out) :: theta,phi

  logical,intent(in) :: is_mesh_elliptical

  ! converts geographic (lat) to geocentric latitude and converts to co-latitude (theta)
  call lat_2_geocentric_colat_dble(lat,theta,is_mesh_elliptical)

  ! longitude (in rad)
  phi = lon * DEGREES_TO_RADIANS

  ! reduces range for colatitude to 0 and PI, for longitude to 0 and 2*PI
  call reduce(theta,phi)

  end subroutine latlon_2_geocentric_thetaphi_dble

!
!------------------------------------------------------------------------------
!

  subroutine lat_2_geocentric_colat_dble(lat_prime,theta,is_mesh_elliptical)

! converts geographic latitude (lat_prime) (in degrees) to geocentric colatitude (theta) (in radians)

  use constants, only: PI_OVER_TWO,DEGREES_TO_RADIANS,ASSUME_PERFECT_SPHERE
  use shared_parameters, only: ONE_MINUS_F_SQUARED

  implicit none

  ! latitude (in degrees)
  double precision,intent(in) :: lat_prime
  ! co-latitude (in radians)
  double precision,intent(inout) :: theta
  ! mesh ellipticity
  logical, intent(in) :: is_mesh_elliptical

  ! local parameter
  double precision :: val

  if (.not. ASSUME_PERFECT_SPHERE) then
    ! converts geographic (lat_prime) to geocentric latitude and converts to co-latitude (theta)
    if (is_mesh_elliptical) then
      ! geocentric colatitude, accounting for ellipticity factor
      val = PI_OVER_TWO - atan( ONE_MINUS_F_SQUARED*dtan(lat_prime * DEGREES_TO_RADIANS) )
    else
      ! for perfect sphere, geocentric and geographic latitudes are the same
      ! converts geographic latitude (in degrees) to geocentric co-latitude (in radians)
      val = PI_OVER_TWO - lat_prime * DEGREES_TO_RADIANS
    endif
  else
    ! for perfect sphere, geocentric and geographic latitudes are the same
    ! converts latitude (in degrees) to co-latitude (in radians)
    val = PI_OVER_TWO - lat_prime * DEGREES_TO_RADIANS
  endif

  ! return result
  theta = val

  end subroutine lat_2_geocentric_colat_dble

!
!------------------------------------------------------------------------------
!

  subroutine geographic_2_geocentric_lat_dble(lat_prime,lat,is_mesh_elliptical)

! converts geographic latitude (lat_prime) (in degrees) to geocentric colatitude (theta) (in radians)

  use constants, only: PI_OVER_TWO,DEGREES_TO_RADIANS,RADIANS_TO_DEGREES,ASSUME_PERFECT_SPHERE
  use shared_parameters, only: ONE_MINUS_F_SQUARED

  implicit none

  ! latitude (in degrees) (geographic/geodetic)
  double precision,intent(in) :: lat_prime
  ! latitude (in radians) (geocentric)
  double precision,intent(inout) :: lat
  ! mesh ellipticity
  logical, intent(in) :: is_mesh_elliptical

  ! local parameter
  double precision :: val

  if (.not. ASSUME_PERFECT_SPHERE) then
    ! converts geographic (lat_prime) to geocentric latitude (lat)
    if (is_mesh_elliptical) then
      ! geocentric latitude, accounting for ellipticity factor
      val = datan( ONE_MINUS_F_SQUARED * dtan(lat_prime * DEGREES_TO_RADIANS) ) * RADIANS_TO_DEGREES
    else
      ! for perfect sphere, geocentric and geographic latitudes are the same
      val = lat_prime
    endif
  else
    ! for perfect sphere, geocentric and geographic latitudes are the same
    val = lat_prime
  endif

  ! return result
  lat = val

  end subroutine geographic_2_geocentric_lat_dble

!
!------------------------------------------------------------------------------
!

! not used yet...
!
!  subroutine geographic_2_geocentric_colat_dble(theta_prime,theta)
!
!! converts geographic colatitude (theta_prime) to geocentric colatitude (theta) (in radians)
!
!  use constants, only: PI_OVER_TWO
!  use shared_parameters, only: ONE_MINUS_F_SQUARED
!
!  implicit none
!
!  double precision,intent(in) :: theta_prime
!  double precision,intent(inout) :: theta
!
!  ! temporary variable for theta
!  ! (in case the input and output variable is the same, the compiler might do some weird things otherwise)
!  double precision :: val,lat_prime
!
!  if (.not. ASSUME_PERFECT_SPHERE) then
!    ! mesh is elliptical
!    ! geocentric colatitude, accounting for ellipticity factor
!    ! latitude in rad
!    lat_prime = PI_OVER_TWO - theta_prime
!    ! range from atan is [-PI/2,PI/2], thus theta should always be within [0,PI]
!    val = PI_OVER_TWO - datan( ONE_MINUS_F_SQUARED * dtan(lat_prime) )
!  else
!    ! mesh is spherical, thus geocentric and geographic colatitudes are identical
!    val = theta_prime
!  endif
!
!  ! return result
!  theta = val
!
!  end subroutine geographic_2_geocentric_colat_dble
!
!!
!------------------------------------------------------------------------------
!

  subroutine thetaphi_2_geographic_latlon_dble(theta_in,phi_in,lat,lon,is_mesh_elliptical)

! convert geocentric theta/phi (in rad) to geographic lat/lon (in degrees)

  use constants, only: DEGREES_TO_RADIANS,PI_OVER_TWO,RADIANS_TO_DEGREES

  implicit none

  double precision,intent(in) :: theta_in,phi_in
  double precision,intent(out) :: lat,lon

  logical,intent(in) :: is_mesh_elliptical

  ! local parameters
  double precision :: theta,phi,theta_prime

  ! initializes
  theta = theta_in
  phi = phi_in

  ! reduces range for colatitude to 0 and PI, for longitude to 0 and 2*PI
  call reduce(theta,phi)

  ! converts geocentric to geographic colatitude
  if (is_mesh_elliptical) then
    ! note: for example, the source/station information is given in geographic latitude/longitude
    !       (lat/lon given with respect to a reference ellipsoid).
    !       we need to convert the geocentric mesh positions (theta,phi) to geographic ones (lat/lon),
    !       correcting geocentric latitude for ellipticity if the mesh is elliptical.
    !
    !       for spherical meshes, the conversion from x/y/z to r/theta/phi to get colatitude/longitude is enough
    !       and won't need the correction factor of ellipticity.
    call geocentric_2_geographic_colat_dble(theta,theta_prime)
    ! sets colatitude to geographic colatitude
    theta = theta_prime
  endif

  ! gets geographic latitude and longitude in degrees
  lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
  lon = phi * RADIANS_TO_DEGREES

  end subroutine thetaphi_2_geographic_latlon_dble

!
!------------------------------------------------------------------------------
!

  subroutine xyz_2_rlatlon_dble(x,y,z,r,lat,lon,is_mesh_elliptical)

! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)

  use constants, only: RADIANS_TO_DEGREES,PI_OVER_TWO

  implicit none

  double precision,intent(in) :: x,y,z
  double precision,intent(out) :: r,lat,lon
  logical, intent(in) :: is_mesh_elliptical

  ! local parameters
  double precision :: theta,phi,theta_prime

  ! converts location to (geocentric) position radius/colatitude/longitude
  call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

  ! reduces range for colatitude to 0 and PI, for longitude to 0 and 2*PI
  call reduce(theta,phi)

  ! converts geocentric to geographic colatitude
  if (is_mesh_elliptical) then
    ! note: for example, the source/station information is given in geographic latitude/longitude
    !       (lat/lon given with respect to a reference ellipsoid).
    !       we need to convert the geocentric mesh positions (theta,phi) to geographic ones (lat/lon),
    !       correcting geocentric latitude for ellipticity if the mesh is elliptical.
    !
    !       for spherical meshes, the conversion from x/y/z to r/theta/phi to get colatitude/longitude is enough
    !       and won't need the correction factor of ellipticity.
    call geocentric_2_geographic_colat_dble(theta,theta_prime)
    ! sets colatitude to geographic colatitude
    theta = theta_prime
  endif

  ! gets geographic latitude and longitude in degrees
  lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
  lon = phi * RADIANS_TO_DEGREES

  end subroutine xyz_2_rlatlon_dble

!
!------------------------------------------------------------------------------
!

  subroutine geocentric_2_geographic_colat_dble(theta,theta_prime)

! converts geocentric colatitude (theta) to geographic colatitude (theta_prime) (in radians)

! see the book of Dahlen and Tromp, 1998: chapter 14.1.5 Geographic versus geocentric colatitude

! (see also in file geographic_versus_geocentric_coordinates_from_Dahlen_Tromp_1998.pdf in the "doc" directory of the code)

!
! The location of a seismic station situated upon the Earth's surface is conventionally specified
! by giving its elevation (e) above the geoid, its longitude (phi) with respect to Greenwich,
! and its geographical colatitude (theta_prime), which is the angle between the normal n to the reference ellipsoid
! and the normal z to the equatorial plane. (...)
! It is necessary to convert the geographic colatitude (theta_prime) to the corresponding geocentric colatitude (theta),
! which is the angle between the radius vector r and z.
!
!
! @jeroen:
! "
! The conversion is given in the book of Dahlen and Tromp, 1998 (equation 14.32),
! and using the observed ellipticity given by (14.24) this gives 1.00670466 (not 1.006760466).
! If you use the hydrostatic value given by (14.23),
! you would get 1.0066711.
! One could argue for either one, but I prefer to use the best fitting observed value.
!
! With regards to the conversion from geographic to geocentric coordinates:
! all station and event locations are reported by bulletins in map coordinates,
! that is, geographic coordinates. Thus the assumption should be that
! locations given in the STATIONS and CMTSOLUTION files are geographic
! coordinates. When you are reading a topographic map, again the
! assumption is that we are dealing with geographic coordinates,
! and those coordinates need to be converted to geocentric coordinates.
! "
!
! for WGS84 ellipsoid, eccentricity squared e**2 = 0.00669437999014
!
! used here:           e**2 = 0.0066600069180170474
!
! the geocentric latitude (lat) and geographic latitude (lat_prime) derives as
!   tan(lat) = (1 - e**2) * tan(lat_prime)           (1)
!
! from (1) it is clear that the geographic latitude can be written as
!   lat_prime = atan( 1/(1 - e**2) * tan(lat) )
!
! since
!   latitude lat = PI/2 - theta
! with theta being colatitude and
!   tan( PI/2 - theta ) = cot( theta ) = 1/tan(theta) = cos(theta)/sin(theta)
!
! we find that the geographic colatitude (theta_prime) and geocentric colatitude (theta) relate as
!   tan(theta_prime) = ( 1 - e**2 ) tan(theta)
!
! or using the latitude expressions from above, we have the geographic colatitude (theta_prime) as
!   theta_prime = PI/2 - atan( 1/(1 - e**2) * 1/tan(theta) )
!               = PI/2 - atan( 1/(1 - e**2) * cos(theta)/sin(theta) )
!               = PI/2 - atan( 1.00670466  * cos(theta)/sin(theta) )

! -----------------------------------------------------------------------------------

! Email from Daniel Peter, May 2014:

! First, these are some excerpts from: Dahlen & Tromp, 1998: chapter
! 14.1.5 Geographic versus geocentric colatitude

! The location of a seismic station situated upon the Earth's surface is conventionally specified
! by giving its elevation (e) above the geoid, its longitude (phi) with respect to Greenwich,
! and its geographical colatitude (theta_prime), which is the angle
! between the normal n to the reference ellipsoid
! and the normal z to the equatorial plane. (...)
! it is necessary to convert the geographic colatitude (theta_prime) to
! the corresponding geocentric colatitude (theta),
! which is the angle between the radius vector r and z.

! Second, a while ago, in our email exchange with Jeroen about this, he states:

! "
! The conversion is given by D & T (14.32), and using the observed ellipticity given by (14.24)
! this gives 1.00670466 (not 1.006760466).
! If you use the hydrostatic value given by (14.23),
! you would get 1.0066711.
! One could argue for either one, but I prefer to use the best fitting observed value.
!
! With regards to the conversion from geographic to geocentric coordinates:
! all station and event locations are reported by bulletins in map coordinates,
! that is, geographic coordinates. Thus the assumption should be that
! locations given in the STATIONS and CMTSOLUTION files are geographic
! coordinates. When you are reading a topographic map, again the
! assumption is that we are dealing with geographic coordinates,
! and those coordinates need to be converted to geocentric coordinates.
! "

! It is a tricky thing to have a correct way of doing this, but I think we
! are coherent when dealing with it, also in case ellipticity is turned
! on. here some points to consider:

! mesher:

! - topography:
! topography is always added on top of a spherical mesh
! -> the geocentric mesh positions are converted to geographic ones to
! find the corresponding elevation

! - ellipticity:
! is added at the very end of the compute_element_** routines, i.e. after
! assigning mantle properties and topography.
! it stretches out all mesh nodes in radial direction.
! -> the ellipticity factor is based on the geocentric colatitude of the
! initial mesh node location.

! I think in the mesher the conversions are all fine and we don't need to
! check for ELLIPTICITY in the conversion routines. topography data is in
! general given in geographic positions, so we need to convert from
! geocentric to geographic locations, no matter what ellipticity flag we
! have set. in that sense, the flag ASSUME_PERFECT_SPHERE is overriding
! this assumption, i.e. then we assume that topography data would be given
! at geocentric positions.

! The tricky part comes in the solver:

! solver:

! - moment tensor:
! the moment tensor of the CMTSOLUTION is rotated from a spherical to a
! Cartesian representation, using the geocentric colatitude
! - vertical component:
! the vertical component at source/receiver locations is calculated with
! respect to geocentric colatitude
! - target position in the mesh:
! to find the best source/receiver position in the mesh, the target
! location is using the geocentric colatitude and corrects the radius for
! topography and ellipticity when calculating the x/y/z Cartesian
! position. We then locate this x/y/z location in the mesh.

! Note that we always convert the geographic latitude to geocentric
! colatitude for source/station positions, unless ASSUME_PERFECT_SPHERE is
! set. if that flag is set, we assume that source/station locations are
! given in geocentric positions already.

! I think the confusing part is for the elliptical mesh. in case of a
! spherical mesh, we would locate the geocentric position in the spherical
! mesh, which should be fine. however, for the elliptical case where
! ELLIPTICITY is turn on, the mesh would be stretched out and x/y/z
! positions of the mesh nodes are then what? geocentric or geographic now?
! and what is the Cartesian position of a location given in geographic
! coordinates? since we map topography onto geocentric colatitude and
! stretch out the ellipticity in radial direction, i think we always have
! to convert first the source/receiver positions to geocentric ones to
! then locate the best Cartesian x/y/z location in the mesh. so, i think,
! we are doing it the right way.

! Not sure if i made this point clear, and i would welcome some more
! thoughts about the rotations. the rotations for components and moment
! tensors are using a vertical based on the geocentric colatitude? i'm not
! sure if this is correct. so, my question would rather be: what is
! vertical for a seismometer?

! I would think that when you install one, you would use the plumb line to
! orient it. this would mean that you are orienting it along the local
! vertical, which is always perpendicular to the geoid. i guess you know
! that the geoid looks like a potato and not like a perfect sphere, so
! vertical on a seismometer seems not to be the same like the radial
! direction in a geocentric framework...

! Anyway, I'm sorry that I'm not giving a clear answer, but hopefully this
! gives some thoughts.

  use constants, only: PI_OVER_TWO,TINYVAL,ASSUME_PERFECT_SPHERE,USE_OLD_VERSION_5_1_5_FORMAT
  use shared_parameters, only: ONE_MINUS_F_SQUARED

  implicit none

  double precision,intent(in) :: theta
  double precision,intent(inout) :: theta_prime

  ! note: september, 2014
  ! factor: 1/(1 - e^2) = 1/(1 - (1 - (1-f)^2)) = 1/( (1-f)^2 )
  !         with eccentricity e^2 = 1 - (1-f)^2
  ! see about Earth flattening in constants.h: flattening factor changed to 1/299.8
  !                                            f = 1/299.8 -> 1/( (1-f)^2 ) = 1.0067046409645724
  double precision :: FACTOR_TAN

  ! temporary variable for theta_prime
  ! (in case the input and output variable is the same, the compiler might do some weird things otherwise)
  double precision :: val

  ! factor
  FACTOR_TAN = 1.d0 / ONE_MINUS_F_SQUARED

  ! note: instead of 1/tan(theta) we take cos(theta)/sin(theta) and avoid division by zero

  if (.not. ASSUME_PERFECT_SPHERE) then
    ! mesh is elliptical
    if (USE_OLD_VERSION_5_1_5_FORMAT) then
      val = PI_OVER_TWO - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
    else
      ! converts geocentric colatitude theta to geographic colatitude theta_prime
      val = PI_OVER_TWO - datan(FACTOR_TAN*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
    endif
  else
    ! mesh is spherical, thus geocentric and geographic colatitudes are identical
    val = theta
  endif

  ! return result
  ! range from atan is [-PI/2,PI/2], thus theta_prime should always be within [0,PI]
  theta_prime = val

  end subroutine geocentric_2_geographic_colat_dble

!
!------------------------------------------------------------------------------
!

  subroutine xyz_2_rlatlon_cr(x_in,y_in,z_in,r_out,lat_out,lon_out,is_mesh_elliptical)

! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)

  use constants, only: CUSTOM_REAL,RADIANS_TO_DEGREES,PI_OVER_TWO

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: x_in,y_in,z_in
  real(kind=CUSTOM_REAL),intent(out) :: r_out,lat_out,lon_out
  logical, intent(in) :: is_mesh_elliptical

  ! local parameters
  double precision :: r,theta,phi,theta_prime,lat,lon
  double precision :: x,y,z

  ! converts to double
  x = x_in
  y = y_in
  z = z_in

  ! converts location to (geocentric) position radius/colatitude/longitude
  call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

  ! reduces range for colatitude to 0 and PI, for longitude to 0 and 2*PI
  call reduce(theta,phi)

  ! converts geocentric to geographic colatitude
  if (is_mesh_elliptical) then
    ! note: for example, the source/station information is given in geographic latitude/longitude
    !       (lat/lon given with respect to a reference ellipsoid).
    !       we need to convert the geocentric mesh positions (theta,phi) to geographic ones (lat/lon),
    !       correcting geocentric latitude for ellipticity if the mesh is elliptical.
    !
    !       for spherical meshes, the conversion from x/y/z to r/theta/phi to get colatitude/longitude is enough
    !       and won't need the correction factor of ellipticity.
    call geocentric_2_geographic_colat_dble(theta,theta_prime)
    ! sets colatitude to geographic colatitude
    theta = theta_prime
  endif

  ! gets geographic latitude and longitude in degrees
  lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
  lon = phi * RADIANS_TO_DEGREES

  ! converts to custom real
  r_out = real(r,kind=CUSTOM_REAL)
  lat_out = real(lat,kind=CUSTOM_REAL)
  lon_out = real(lon,kind=CUSTOM_REAL)

  end subroutine xyz_2_rlatlon_cr

!
!------------------------------------------------------------------------------
!

! not used yet...
!
!  subroutine geocentric_2_geographic_colat_cr(theta,theta_prime)
!
!! converts geocentric colatitude (theta) to geographic colatitude (theta_prime) (in radians)
!
!  use constants, only: CUSTOM_REAL
!
!  implicit none
!
!  real(kind=CUSTOM_REAL),intent(in) :: theta
!  real(kind=CUSTOM_REAL),intent(inout) :: theta_prime
!
!  ! local parameters
!  double precision :: dtheta,dtheta_prime
!
!  ! gets double precision value
!  dtheta = dble(theta)
!
!  ! converts geocentric to geographic colatitude
!  call geocentric_2_geographic_dble(dtheta,dtheta_prime)
!
!  ! gets custom-real value
!  theta_prime = real(dtheta_prime, kind=CUSTOM_REAL)
!
!  end subroutine geocentric_2_geographic_colat_cr
!
!
!------------------------------------------------------------------------------
!

  subroutine xyz_2_latlon_minmax(nspec,nglob,ibool,xstore,ystore,zstore, &
                                 lat_min,lat_max,lon_min,lon_max,is_mesh_elliptical)

! returns minimum and maximum values of latitude/longitude of given mesh points;
! latitude in degree between [-90,90], longitude in degree between [0,360]

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,HUGEVAL

  implicit none

  integer,intent(in) :: nspec,nglob

  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nspec)

  ! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore,ystore,zstore

  double precision,intent(out) :: lat_min,lat_max,lon_min,lon_max
  logical, intent(in) :: is_mesh_elliptical

  ! local parameters
  double precision :: x,y,z
  double precision :: r,lat,lon
  !double precision :: r_min,r_max
  integer :: ispec,i,j,k,iglob

  ! initializes
  lat_min = HUGEVAL
  lat_max = -HUGEVAL
  lon_min = HUGEVAL
  lon_max = -HUGEVAL
  !r_min = HUGEVAL
  !r_max = -HUGEVAL

  ! loops over all elements
  do ispec = 1,nspec

    ! loops only over corners
    do k = 1,NGLLZ,NGLLZ-1
      do j = 1,NGLLY,NGLLY-1
        do i = 1,NGLLX,NGLLX-1

          ! gets x/y/z coordinates
          iglob = ibool(i,j,k,ispec)
          x = xstore(iglob)
          y = ystore(iglob)
          z = zstore(iglob)

          ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
          call xyz_2_rlatlon_dble(x,y,z,r,lat,lon,is_mesh_elliptical)

          ! stores min/max
          if (lat < lat_min ) lat_min = lat
          if (lat > lat_max ) lat_max = lat

          if (lon < lon_min ) lon_min = lon
          if (lon > lon_max ) lon_max = lon

          !if (r < r_min ) r_min = r
          !if (r > r_max ) r_max = r
        enddo
      enddo
    enddo
  enddo

  ! limits latitude to [-90.0,90.0]
  if (lat_min < -90.d0 ) lat_min = -90.d0
  if (lat_max > 90.d0 ) lat_max = 90.d0

  ! limits longitude to [0.0,360.0]
  if (lon_min < 0.d0 ) lon_min = lon_min + 360.d0
  if (lon_min > 360.d0 ) lon_min = lon_min - 360.d0
  if (lon_max < 0.d0 ) lon_max = lon_max + 360.d0
  if (lon_max > 360.d0 ) lon_max = lon_max - 360.d0

  end subroutine xyz_2_latlon_minmax

!
!-------------------------------------------------------------------------------------------
!

  subroutine get_greatcircle_distance(theta1,phi1,theta2,phi2,dist)

! returns great-circle distance (by haversine formula) in radians

  implicit none

  double precision, intent(in) :: theta1,phi1,theta2,phi2    ! positions colat/phi in rad
  double precision, intent(out) :: dist
  ! local parameters
  double precision :: mid_theta,mid_phi,sintheta,sinphi
  double precision :: a

! note: the haversine formula calculates the great-circle distance for positions on a sphere.
!       It is slightly more accurate than the law of cosines formula:
!          d = acos( sin(theta1) * sin(theta2) + cos(theta1) * cos(theta2) * cos(phi1 - phi2) )
!       Both are valid for distances on a sphere.
!
!       We could correct the geodetic lat/lon positions to get geocentric lat/lon and then calculate the distance.
!       This would be more accurate, but adds additional computational costs.

  ! using latitudes:
  !  lat1 = PI_OVER_TWO - theta1    ! in rad
  !  lat2 = PI_OVER_TWO - theta2
  !  ! haversine formula
  !  mid_lat = 0.5d0 * (lat1 - lat2)
  !  mid_lon = 0.5d0 * (phi1 - phi2)
  !  sinlat = dsin(mid_lat)
  !  sinlon = dsin(mid_lon)
  !  a = dsqrt( sinlat*sinlat + dcos(lat1) * dcos(lat2) * sinlon*sinlon )

  ! here, we're using colatitudes:
  !  -> lat1 - lat2 = pi/2 - theta1 - (pi/2 - theta2)
  !                 = theta2 - theta1
  !     it doesn't matter though if we use (theta2 - theta1) or (theta1 - theta2) below as it will be used
  !     together with sin^2.
  !
  !  -> cos(lat1)   = cos(pi/2 - theta1)
  !                 = sin(theta1)
  !     this is the main difference wrt/ to the formula for latitudes

  ! haversine formula
  mid_theta = 0.5d0 * (theta1 - theta2)
  mid_phi = 0.5d0 * (phi1 - phi2)

  sintheta = dsin(mid_theta)
  sinphi = dsin(mid_phi)

  a = dsqrt( sintheta*sintheta + dsin(theta1) * dsin(theta2) * sinphi*sinphi )

  if (abs(a) > 1.d0) a = sign(1.d0,a) * 1.d0 ! limit to +/- 1 numerical accuracy for asin(..)

  ! great-circle (epicentral) distance in rad
  dist = 2.d0 * dasin(a)

  end subroutine get_greatcircle_distance

!
!-------------------------------------------------------------------------------------------
!

! not used yet...
!
!  subroutine get_greatcircle_distance_cosines(theta1,phi1,theta2,phi2,dist)
!
!! returns great-circle distance (by law of cosines formula) in radians
!
!  use constants, only: RADIANS_TO_DEGREES
!
!  implicit none
!
!  double precision, intent(in) :: theta1,phi1,theta2,phi2    ! positions colat/phi in rad
!  double precision, intent(out) :: dist
!  ! local parameters
!  double precision :: a
!
!  ! (formula by law of cosines)
!  a = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2)
!
!  if (abs(a) > 1.d0) a = sign(1.d0,a) * 1.d0 ! limit to +/- 1 numerical accuracy for acos(..)
!
!  ! great-circle (epicentral) distance in rad
!  dist = acos(a)
!
!  end subroutine get_greatcircle_distance_cosines
!
!
!
!-------------------------------------------------------------------------------------------
!

! not used yet...
!
!  subroutine get_greatcircle_distance_latlon_haversine(lat1,lon1,lat2,lon2,dist)
!
!! returns great-circle distance (by haversine formula) in degrees
!
!  use constants, only: DEGREES_TO_RADIANS, RADIANS_TO_DEGREES, TWO_PI
!
!  implicit none
!
!  double precision, intent(in) :: lat1,lon1,lat2,lon2
!  double precision, intent(out) :: dist
!
!  ! local parameters
!  double precision :: lat1_rad,lon1_rad,lat2_rad,lon2_rad
!  double precision :: mid_lat,mid_lon,sinlat,sinlon
!  double precision :: a,c
!
!! note: the haversine formula calculates the great-circle distance for lat/lon positions on a sphere.
!!       It is slightly more accurate than the law of cosines formula:
!!          d = acos( sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2) )
!!       Both are valid for distances on a sphere.
!!
!!       We could correct the geodetic lat/lon positions to get geocentric lat/lon and then calculate the distance.
!!       This would be more accurate, but adds additional computational costs.
!!
!!       Given we only use epicentral distances for output infos, we reduce here the costs and
!!       use the simple haversine formula.
!
!  ! converts to radians
!  lat1_rad = lat1 * DEGREES_TO_RADIANS
!  lon1_rad = lon1 * DEGREES_TO_RADIANS
!
!  lat2_rad = lat2 * DEGREES_TO_RADIANS
!  lon2_rad = lon2 * DEGREES_TO_RADIANS
!
!  mid_lat = 0.5d0 * (lat1_rad - lat2_rad)
!  mid_lon = 0.5d0 * (lon1_rad - lon2_rad)
!
!  sinlat = dsin(mid_lat)
!  sinlon = dsin(mid_lon)
!
!  a = dsqrt( sinlat*sinlat + dcos(lat1_rad) * dcos(lat2_rad) * sinlon*sinlon)
!
!  ! great-circle (epicentral) distance in rad
!  c = 2.d0 * asin(a)
!
!  ! distance in degrees
!  dist = c * RADIANS_TO_DEGREES
!
!  end subroutine get_greatcircle_distance_latlon_haversine



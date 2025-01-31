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

! compute the Euler angles and the associated rotation matrix

  subroutine euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

  use constants, only: DEGREES_TO_RADIANS

  implicit none

  double precision,intent(out) :: rotation_matrix(3,3)
  double precision,intent(in) :: CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH

  ! local parameters
  double precision :: alpha,beta,gamma
  double precision :: sina,cosa,sinb,cosb,sing,cosg

! compute colatitude and longitude and convert to radians
  alpha = CENTER_LONGITUDE_IN_DEGREES * DEGREES_TO_RADIANS
  beta = (90.0d0 - CENTER_LATITUDE_IN_DEGREES) * DEGREES_TO_RADIANS
  gamma = GAMMA_ROTATION_AZIMUTH * DEGREES_TO_RADIANS

  sina = dsin(alpha)
  cosa = dcos(alpha)
  sinb = dsin(beta)
  cosb = dcos(beta)
  sing = dsin(gamma)
  cosg = dcos(gamma)

! define rotation matrix
  rotation_matrix(1,1) = cosg*cosb*cosa-sing*sina
  rotation_matrix(1,2) = -sing*cosb*cosa-cosg*sina
  rotation_matrix(1,3) = sinb*cosa
  rotation_matrix(2,1) = cosg*cosb*sina+sing*cosa
  rotation_matrix(2,2) = -sing*cosb*sina+cosg*cosa
  rotation_matrix(2,3) = sinb*sina
  rotation_matrix(3,1) = -cosg*sinb
  rotation_matrix(3,2) = sing*sinb
  rotation_matrix(3,3) = cosb

  end subroutine euler_angles

!
!-------------------------------------------------------------------------------------------
!

  subroutine determine_chunk_corners_latlon(CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH, &
                                            ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
                                            corners_lat,corners_lon)

  use constants, only: DEGREES_TO_RADIANS,RADIANS_TO_DEGREES,ONE,PI,TWO_PI,PI_OVER_TWO,R_UNIT_SPHERE

  implicit none

  double precision,intent(in) :: CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH
  double precision,intent(in) :: ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES
  double precision,intent(out) :: corners_lat(4),corners_lon(4)

  ! local parameters
  ! rotation matrix from Euler angles
  integer :: i,j,ix,iy,icorner
  double precision :: rotation_matrix(3,3)
  double precision :: vector_ori(3),vector_rotated(3)
  double precision :: r_corner,theta_corner,phi_corner,lat,long,colat_corner
  double precision :: x,y,gamma,rgt,xi,eta
  double precision :: x_top,y_top,z_top
  double precision :: ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

  ! initializes
  corners_lat(:) = 0.d0
  corners_lon(:) = 0.d0

  ! compute rotation matrix from Euler angles
  call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

  ! convert width to radians
  ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * DEGREES_TO_RADIANS
  ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES * DEGREES_TO_RADIANS

  ! loop on the four corners of the chunk to display their coordinates
  icorner = 0
  do iy = 0,1
    do ix = 0,1

    icorner = icorner + 1

    xi  = - ANGULAR_WIDTH_XI_RAD/2.  + dble(ix)*ANGULAR_WIDTH_XI_RAD
    eta = - ANGULAR_WIDTH_ETA_RAD/2. + dble(iy)*ANGULAR_WIDTH_ETA_RAD

    x = dtan(xi)
    y = dtan(eta)

    gamma = ONE/dsqrt(ONE+x*x+y*y)
    rgt = R_UNIT_SPHERE*gamma

    ! define the mesh points at the top surface
    x_top = -y*rgt
    y_top = x*rgt
    z_top = rgt

    ! rotate top
    vector_ori(1) = x_top
    vector_ori(2) = y_top
    vector_ori(3) = z_top
    do i=1,3
      vector_rotated(i) = 0.0d0
      do j=1,3
        vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
      enddo
    enddo
    x_top = vector_rotated(1)
    y_top = vector_rotated(2)
    z_top = vector_rotated(3)

    ! convert to latitude and longitude
    call xyz_2_rthetaphi_dble(x_top,y_top,z_top,r_corner,theta_corner,phi_corner)
    call reduce(theta_corner,phi_corner)

    ! convert geocentric to geographic colatitude
    call geocentric_2_geographic_dble(theta_corner,colat_corner)

    if (phi_corner > PI) phi_corner = phi_corner - TWO_PI

    ! compute real position of the source
    lat = (PI_OVER_TWO-colat_corner) * RADIANS_TO_DEGREES
    long = phi_corner * RADIANS_TO_DEGREES

    corners_lat(icorner) = lat
    corners_lon(icorner) = long

    enddo
  enddo

  end subroutine determine_chunk_corners_latlon

!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, October 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine get_jacobian_boundaries(myrank,iboun,nspec,xstore,ystore,zstore, &
    dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
              jacobian2D_xmin,jacobian2D_xmax, &
              jacobian2D_ymin,jacobian2D_ymax, &
              jacobian2D_bottom,jacobian2D_top, &
              normal_xmin,normal_xmax, &
              normal_ymin,normal_ymax, &
              normal_bottom,normal_top, &
              NSPEC2D_BOTTOM,NSPEC2D_TOP, &
              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

  implicit none

  include "constants.h"

  integer nspec,myrank,dummy_var,ispec_tmp
  integer NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX

  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer ibelm_bottom(NSPEC2D_BOTTOM),ibelm_top(NSPEC2D_TOP)

  logical iboun(6,nspec)

  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  real(kind=CUSTOM_REAL) jacobian2D_xmin(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) jacobian2D_xmax(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) jacobian2D_ymin(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) jacobian2D_ymax(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM)
  real(kind=CUSTOM_REAL) jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP)

  real(kind=CUSTOM_REAL) normal_xmin(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) normal_xmax(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) normal_ymin(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) normal_ymax(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM)
  real(kind=CUSTOM_REAL) normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP)

  double precision dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ)
  double precision dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ)
  double precision dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY)
  double precision dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY)

! global element numbering
  integer ispec

! counters to keep track of number of elements on each of the boundaries
  integer ispecb1,ispecb2,ispecb3,ispecb4,ispecb5,ispecb6

  double precision xelm(NGNOD2D),yelm(NGNOD2D),zelm(NGNOD2D)

! arrays for sorting routine
  integer, dimension(:), allocatable :: ind,ninseg,iglob,locval,iwork
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: work
  double precision, dimension(:), allocatable :: xstore_selected,ystore_selected,zstore_selected

! check that the parameter file is correct
  if(NGNOD /= 27) call exit_MPI(myrank,'elements should have 27 control nodes')
  if(NGNOD2D /= 9) call exit_MPI(myrank,'surface elements should have 9 control nodes')

  ispecb1 = 0
  ispecb2 = 0
  ispecb3 = 0
  ispecb4 = 0
  ispecb5 = 0
  ispecb6 = 0

  do ispec=1,nspec

! determine if the element falls on a boundary

! on boundary: xmin

  if(iboun(1,ispec)) then

    ispecb1=ispecb1+1
    ibelm_xmin(ispecb1)=ispec

!   specify the 9 nodes for the 2-D boundary element
    xelm(1)=xstore(1,1,1,ispec)
    yelm(1)=ystore(1,1,1,ispec)
    zelm(1)=zstore(1,1,1,ispec)
    xelm(2)=xstore(1,NGLLY,1,ispec)
    yelm(2)=ystore(1,NGLLY,1,ispec)
    zelm(2)=zstore(1,NGLLY,1,ispec)
    xelm(3)=xstore(1,NGLLY,NGLLZ,ispec)
    yelm(3)=ystore(1,NGLLY,NGLLZ,ispec)
    zelm(3)=zstore(1,NGLLY,NGLLZ,ispec)
    xelm(4)=xstore(1,1,NGLLZ,ispec)
    yelm(4)=ystore(1,1,NGLLZ,ispec)
    zelm(4)=zstore(1,1,NGLLZ,ispec)
    xelm(5)=xstore(1,(NGLLY+1)/2,1,ispec)
    yelm(5)=ystore(1,(NGLLY+1)/2,1,ispec)
    zelm(5)=zstore(1,(NGLLY+1)/2,1,ispec)
    xelm(6)=xstore(1,NGLLY,(NGLLZ+1)/2,ispec)
    yelm(6)=ystore(1,NGLLY,(NGLLZ+1)/2,ispec)
    zelm(6)=zstore(1,NGLLY,(NGLLZ+1)/2,ispec)
    xelm(7)=xstore(1,(NGLLY+1)/2,NGLLZ,ispec)
    yelm(7)=ystore(1,(NGLLY+1)/2,NGLLZ,ispec)
    zelm(7)=zstore(1,(NGLLY+1)/2,NGLLZ,ispec)
    xelm(8)=xstore(1,1,(NGLLZ+1)/2,ispec)
    yelm(8)=ystore(1,1,(NGLLZ+1)/2,ispec)
    zelm(8)=zstore(1,1,(NGLLZ+1)/2,ispec)
    xelm(9)=xstore(1,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
    yelm(9)=ystore(1,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
    zelm(9)=zstore(1,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)

    call compute_jacobian_2D(myrank,ispecb1,xelm,yelm,zelm,dershape2D_x, &
                  jacobian2D_xmin,normal_xmin,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
    !! test reord
    allocate (xstore_selected(ispecb1))
    allocate (ystore_selected(ispecb1))
    allocate (zstore_selected(ispecb1))
    allocate(ind(ispecb1))
    allocate(ninseg(ispecb1))
    allocate(iglob(ispecb1))
    allocate(locval(ispecb1))
    allocate(ifseg(ispecb1))
    allocate(iwork(ispecb1))
    allocate(work(ispecb1))

    do ispec_tmp=1,ispecb1
      xstore_selected(ispec_tmp) = xstore(1,1,1,ibelm_xmin(ispec_tmp))
      ystore_selected(ispec_tmp) = ystore(1,1,1,ibelm_xmin(ispec_tmp))
      zstore_selected(ispec_tmp) = zstore(1,1,1,ibelm_xmin(ispec_tmp))
    enddo

    call sort_array_coordinates(ispecb1,xstore_selected,ystore_selected,zstore_selected, &
            ibelm_xmin,iglob,locval,ifseg,dummy_var,ind,ninseg,iwork,work)

    deallocate (xstore_selected)
    deallocate (ystore_selected)
    deallocate (zstore_selected)
    deallocate(ind)
    deallocate(ninseg)
    deallocate(iglob)
    deallocate(locval)
    deallocate(ifseg)
    deallocate(iwork)
    deallocate(work)
    !! test reord
  endif

! on boundary: xmax

  if(iboun(2,ispec)) then

    ispecb2=ispecb2+1
    ibelm_xmax(ispecb2)=ispec

!   specify the 9 nodes for the 2-D boundary element
    xelm(1)=xstore(NGLLX,1,1,ispec)
    yelm(1)=ystore(NGLLX,1,1,ispec)
    zelm(1)=zstore(NGLLX,1,1,ispec)
    xelm(2)=xstore(NGLLX,NGLLY,1,ispec)
    yelm(2)=ystore(NGLLX,NGLLY,1,ispec)
    zelm(2)=zstore(NGLLX,NGLLY,1,ispec)
    xelm(3)=xstore(NGLLX,NGLLY,NGLLZ,ispec)
    yelm(3)=ystore(NGLLX,NGLLY,NGLLZ,ispec)
    zelm(3)=zstore(NGLLX,NGLLY,NGLLZ,ispec)
    xelm(4)=xstore(NGLLX,1,NGLLZ,ispec)
    yelm(4)=ystore(NGLLX,1,NGLLZ,ispec)
    zelm(4)=zstore(NGLLX,1,NGLLZ,ispec)
    xelm(5)=xstore(NGLLX,(NGLLY+1)/2,1,ispec)
    yelm(5)=ystore(NGLLX,(NGLLY+1)/2,1,ispec)
    zelm(5)=zstore(NGLLX,(NGLLY+1)/2,1,ispec)
    xelm(6)=xstore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
    yelm(6)=ystore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
    zelm(6)=zstore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
    xelm(7)=xstore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
    yelm(7)=ystore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
    zelm(7)=zstore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
    xelm(8)=xstore(NGLLX,1,(NGLLZ+1)/2,ispec)
    yelm(8)=ystore(NGLLX,1,(NGLLZ+1)/2,ispec)
    zelm(8)=zstore(NGLLX,1,(NGLLZ+1)/2,ispec)
    xelm(9)=xstore(NGLLX,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
    yelm(9)=ystore(NGLLX,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
    zelm(9)=zstore(NGLLX,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)

    call compute_jacobian_2D(myrank,ispecb2,xelm,yelm,zelm,dershape2D_x, &
                  jacobian2D_xmax,normal_xmax,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
    !! test reord
    allocate (xstore_selected(ispecb2))
    allocate (ystore_selected(ispecb2))
    allocate (zstore_selected(ispecb2))
    allocate(ind(ispecb2))
    allocate(ninseg(ispecb2))
    allocate(iglob(ispecb2))
    allocate(locval(ispecb2))
    allocate(ifseg(ispecb2))
    allocate(iwork(ispecb2))
    allocate(work(ispecb2))

    do ispec_tmp=1,ispecb2
      xstore_selected(ispec_tmp) = xstore(1,1,1,ibelm_xmax(ispec_tmp))
      ystore_selected(ispec_tmp) = ystore(1,1,1,ibelm_xmax(ispec_tmp))
      zstore_selected(ispec_tmp) = zstore(1,1,1,ibelm_xmax(ispec_tmp))
    enddo

    call sort_array_coordinates(ispecb2,xstore_selected,ystore_selected,zstore_selected, &
            ibelm_xmax,iglob,locval,ifseg,dummy_var,ind,ninseg,iwork,work)

    deallocate (xstore_selected)
    deallocate (ystore_selected)
    deallocate (zstore_selected)
    deallocate(ind)
    deallocate(ninseg)
    deallocate(iglob)
    deallocate(locval)
    deallocate(ifseg)
    deallocate(iwork)
    deallocate(work)
    !! test reord
  endif

! on boundary: ymin

  if(iboun(3,ispec)) then

    ispecb3=ispecb3+1
    ibelm_ymin(ispecb3)=ispec

!   specify the 9 nodes for the 2-D boundary element
    xelm(1)=xstore(1,1,1,ispec)
    yelm(1)=ystore(1,1,1,ispec)
    zelm(1)=zstore(1,1,1,ispec)
    xelm(2)=xstore(NGLLX,1,1,ispec)
    yelm(2)=ystore(NGLLX,1,1,ispec)
    zelm(2)=zstore(NGLLX,1,1,ispec)
    xelm(3)=xstore(NGLLX,1,NGLLZ,ispec)
    yelm(3)=ystore(NGLLX,1,NGLLZ,ispec)
    zelm(3)=zstore(NGLLX,1,NGLLZ,ispec)
    xelm(4)=xstore(1,1,NGLLZ,ispec)
    yelm(4)=ystore(1,1,NGLLZ,ispec)
    zelm(4)=zstore(1,1,NGLLZ,ispec)
    xelm(5)=xstore((NGLLX+1)/2,1,1,ispec)
    yelm(5)=ystore((NGLLX+1)/2,1,1,ispec)
    zelm(5)=zstore((NGLLX+1)/2,1,1,ispec)
    xelm(6)=xstore(NGLLX,1,(NGLLZ+1)/2,ispec)
    yelm(6)=ystore(NGLLX,1,(NGLLZ+1)/2,ispec)
    zelm(6)=zstore(NGLLX,1,(NGLLZ+1)/2,ispec)
    xelm(7)=xstore((NGLLX+1)/2,1,NGLLZ,ispec)
    yelm(7)=ystore((NGLLX+1)/2,1,NGLLZ,ispec)
    zelm(7)=zstore((NGLLX+1)/2,1,NGLLZ,ispec)
    xelm(8)=xstore(1,1,(NGLLZ+1)/2,ispec)
    yelm(8)=ystore(1,1,(NGLLZ+1)/2,ispec)
    zelm(8)=zstore(1,1,(NGLLZ+1)/2,ispec)
    xelm(9)=xstore((NGLLX+1)/2,1,(NGLLZ+1)/2,ispec)
    yelm(9)=ystore((NGLLX+1)/2,1,(NGLLZ+1)/2,ispec)
    zelm(9)=zstore((NGLLX+1)/2,1,(NGLLZ+1)/2,ispec)

    call compute_jacobian_2D(myrank,ispecb3,xelm,yelm,zelm,dershape2D_y, &
                  jacobian2D_ymin,normal_ymin,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
    !! test reord
    allocate (xstore_selected(ispecb3))
    allocate (ystore_selected(ispecb3))
    allocate (zstore_selected(ispecb3))
    allocate(ind(ispecb3))
    allocate(ninseg(ispecb3))
    allocate(iglob(ispecb3))
    allocate(locval(ispecb3))
    allocate(ifseg(ispecb3))
    allocate(iwork(ispecb3))
    allocate(work(ispecb3))

    do ispec_tmp=1,ispecb3
      xstore_selected(ispec_tmp) = xstore(1,1,1,ibelm_ymin(ispec_tmp))
      ystore_selected(ispec_tmp) = ystore(1,1,1,ibelm_ymin(ispec_tmp))
      zstore_selected(ispec_tmp) = zstore(1,1,1,ibelm_ymin(ispec_tmp))
    enddo

    call sort_array_coordinates(ispecb3,xstore_selected,ystore_selected,zstore_selected, &
            ibelm_ymin,iglob,locval,ifseg,dummy_var,ind,ninseg,iwork,work)

    deallocate (xstore_selected)
    deallocate (ystore_selected)
    deallocate (zstore_selected)
    deallocate(ind)
    deallocate(ninseg)
    deallocate(iglob)
    deallocate(locval)
    deallocate(ifseg)
    deallocate(iwork)
    deallocate(work)
    !! test reord
  endif

! on boundary: ymax

  if(iboun(4,ispec)) then

    ispecb4=ispecb4+1
    ibelm_ymax(ispecb4)=ispec

!   specify the 9 nodes for the 2-D boundary element
    xelm(1)=xstore(1,NGLLY,1,ispec)
    yelm(1)=ystore(1,NGLLY,1,ispec)
    zelm(1)=zstore(1,NGLLY,1,ispec)
    xelm(2)=xstore(NGLLX,NGLLY,1,ispec)
    yelm(2)=ystore(NGLLX,NGLLY,1,ispec)
    zelm(2)=zstore(NGLLX,NGLLY,1,ispec)
    xelm(3)=xstore(NGLLX,NGLLY,NGLLZ,ispec)
    yelm(3)=ystore(NGLLX,NGLLY,NGLLZ,ispec)
    zelm(3)=zstore(NGLLX,NGLLY,NGLLZ,ispec)
    xelm(4)=xstore(1,NGLLY,NGLLZ,ispec)
    yelm(4)=ystore(1,NGLLY,NGLLZ,ispec)
    zelm(4)=zstore(1,NGLLY,NGLLZ,ispec)
    xelm(5)=xstore((NGLLX+1)/2,NGLLY,1,ispec)
    yelm(5)=ystore((NGLLX+1)/2,NGLLY,1,ispec)
    zelm(5)=zstore((NGLLX+1)/2,NGLLY,1,ispec)
    xelm(6)=xstore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
    yelm(6)=ystore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
    zelm(6)=zstore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
    xelm(7)=xstore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
    yelm(7)=ystore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
    zelm(7)=zstore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
    xelm(8)=xstore(1,NGLLY,(NGLLZ+1)/2,ispec)
    yelm(8)=ystore(1,NGLLY,(NGLLZ+1)/2,ispec)
    zelm(8)=zstore(1,NGLLY,(NGLLZ+1)/2,ispec)
    xelm(9)=xstore((NGLLX+1)/2,NGLLY,(NGLLZ+1)/2,ispec)
    yelm(9)=ystore((NGLLX+1)/2,NGLLY,(NGLLZ+1)/2,ispec)
    zelm(9)=zstore((NGLLX+1)/2,NGLLY,(NGLLZ+1)/2,ispec)

    call compute_jacobian_2D(myrank,ispecb4,xelm,yelm,zelm,dershape2D_y, &
                  jacobian2D_ymax,normal_ymax,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)

    !! test reord
    allocate (xstore_selected(ispecb4))
    allocate (ystore_selected(ispecb4))
    allocate (zstore_selected(ispecb4))
    allocate(ind(ispecb4))
    allocate(ninseg(ispecb4))
    allocate(iglob(ispecb4))
    allocate(locval(ispecb4))
    allocate(ifseg(ispecb4))
    allocate(iwork(ispecb4))
    allocate(work(ispecb4))

    do ispec_tmp=1,ispecb4
      xstore_selected(ispec_tmp) = xstore(1,1,1,ibelm_ymax(ispec_tmp))
      ystore_selected(ispec_tmp) = ystore(1,1,1,ibelm_ymax(ispec_tmp))
      zstore_selected(ispec_tmp) = zstore(1,1,1,ibelm_ymax(ispec_tmp))
    enddo

    call sort_array_coordinates(ispecb4,xstore_selected,ystore_selected,zstore_selected, &
            ibelm_ymax,iglob,locval,ifseg,dummy_var,ind,ninseg,iwork,work)

    deallocate (xstore_selected)
    deallocate (ystore_selected)
    deallocate (zstore_selected)
    deallocate(ind)
    deallocate(ninseg)
    deallocate(iglob)
    deallocate(locval)
    deallocate(ifseg)
    deallocate(iwork)
    deallocate(work)
    !! test reord

  endif

! on boundary: bottom

  if(iboun(5,ispec)) then

    ispecb5=ispecb5+1
    ibelm_bottom(ispecb5)=ispec

    xelm(1)=xstore(1,1,1,ispec)
    yelm(1)=ystore(1,1,1,ispec)
    zelm(1)=zstore(1,1,1,ispec)
    xelm(2)=xstore(NGLLX,1,1,ispec)
    yelm(2)=ystore(NGLLX,1,1,ispec)
    zelm(2)=zstore(NGLLX,1,1,ispec)
    xelm(3)=xstore(NGLLX,NGLLY,1,ispec)
    yelm(3)=ystore(NGLLX,NGLLY,1,ispec)
    zelm(3)=zstore(NGLLX,NGLLY,1,ispec)
    xelm(4)=xstore(1,NGLLY,1,ispec)
    yelm(4)=ystore(1,NGLLY,1,ispec)
    zelm(4)=zstore(1,NGLLY,1,ispec)
    xelm(5)=xstore((NGLLX+1)/2,1,1,ispec)
    yelm(5)=ystore((NGLLX+1)/2,1,1,ispec)
    zelm(5)=zstore((NGLLX+1)/2,1,1,ispec)
    xelm(6)=xstore(NGLLX,(NGLLY+1)/2,1,ispec)
    yelm(6)=ystore(NGLLX,(NGLLY+1)/2,1,ispec)
    zelm(6)=zstore(NGLLX,(NGLLY+1)/2,1,ispec)
    xelm(7)=xstore((NGLLX+1)/2,NGLLY,1,ispec)
    yelm(7)=ystore((NGLLX+1)/2,NGLLY,1,ispec)
    zelm(7)=zstore((NGLLX+1)/2,NGLLY,1,ispec)
    xelm(8)=xstore(1,(NGLLY+1)/2,1,ispec)
    yelm(8)=ystore(1,(NGLLY+1)/2,1,ispec)
    zelm(8)=zstore(1,(NGLLY+1)/2,1,ispec)
    xelm(9)=xstore((NGLLX+1)/2,(NGLLY+1)/2,1,ispec)
    yelm(9)=ystore((NGLLX+1)/2,(NGLLY+1)/2,1,ispec)
    zelm(9)=zstore((NGLLX+1)/2,(NGLLY+1)/2,1,ispec)

    call compute_jacobian_2D(myrank,ispecb5,xelm,yelm,zelm,dershape2D_bottom, &
                  jacobian2D_bottom,normal_bottom,NGLLX,NGLLY,NSPEC2D_BOTTOM)

    !! test reord
    allocate (xstore_selected(ispecb5))
    allocate (ystore_selected(ispecb5))
    allocate (zstore_selected(ispecb5))
    allocate(ind(ispecb5))
    allocate(ninseg(ispecb5))
    allocate(iglob(ispecb5))
    allocate(locval(ispecb5))
    allocate(ifseg(ispecb5))
    allocate(iwork(ispecb5))
    allocate(work(ispecb5))

    do ispec_tmp=1,ispecb5
      xstore_selected(ispec_tmp) = xstore(1,1,1,ibelm_bottom(ispec_tmp))
      ystore_selected(ispec_tmp) = ystore(1,1,1,ibelm_bottom(ispec_tmp))
      zstore_selected(ispec_tmp) = zstore(1,1,1,ibelm_bottom(ispec_tmp))
    enddo

    call sort_array_coordinates(ispecb5,xstore_selected,ystore_selected,zstore_selected, &
            ibelm_bottom,iglob,locval,ifseg,dummy_var,ind,ninseg,iwork,work)

    deallocate (xstore_selected)
    deallocate (ystore_selected)
    deallocate (zstore_selected)
    deallocate(ind)
    deallocate(ninseg)
    deallocate(iglob)
    deallocate(locval)
    deallocate(ifseg)
    deallocate(iwork)
    deallocate(work)
    !! test reord

  endif

! on boundary: top

  if(iboun(6,ispec)) then

    ispecb6=ispecb6+1
    ibelm_top(ispecb6)=ispec

    xelm(1)=xstore(1,1,NGLLZ,ispec)
    yelm(1)=ystore(1,1,NGLLZ,ispec)
    zelm(1)=zstore(1,1,NGLLZ,ispec)
    xelm(2)=xstore(NGLLX,1,NGLLZ,ispec)
    yelm(2)=ystore(NGLLX,1,NGLLZ,ispec)
    zelm(2)=zstore(NGLLX,1,NGLLZ,ispec)
    xelm(3)=xstore(NGLLX,NGLLY,NGLLZ,ispec)
    yelm(3)=ystore(NGLLX,NGLLY,NGLLZ,ispec)
    zelm(3)=zstore(NGLLX,NGLLY,NGLLZ,ispec)
    xelm(4)=xstore(1,NGLLY,NGLLZ,ispec)
    yelm(4)=ystore(1,NGLLY,NGLLZ,ispec)
    zelm(4)=zstore(1,NGLLY,NGLLZ,ispec)
    xelm(5)=xstore((NGLLX+1)/2,1,NGLLZ,ispec)
    yelm(5)=ystore((NGLLX+1)/2,1,NGLLZ,ispec)
    zelm(5)=zstore((NGLLX+1)/2,1,NGLLZ,ispec)
    xelm(6)=xstore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
    yelm(6)=ystore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
    zelm(6)=zstore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
    xelm(7)=xstore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
    yelm(7)=ystore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
    zelm(7)=zstore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
    xelm(8)=xstore(1,(NGLLY+1)/2,NGLLZ,ispec)
    yelm(8)=ystore(1,(NGLLY+1)/2,NGLLZ,ispec)
    zelm(8)=zstore(1,(NGLLY+1)/2,NGLLZ,ispec)
    xelm(9)=xstore((NGLLX+1)/2,(NGLLY+1)/2,NGLLZ,ispec)
    yelm(9)=ystore((NGLLX+1)/2,(NGLLY+1)/2,NGLLZ,ispec)
    zelm(9)=zstore((NGLLX+1)/2,(NGLLY+1)/2,NGLLZ,ispec)

    call compute_jacobian_2D(myrank,ispecb6,xelm,yelm,zelm,dershape2D_top, &
                  jacobian2D_top,normal_top,NGLLX,NGLLY,NSPEC2D_TOP)

    !! test reord
    allocate (xstore_selected(ispecb6))
    allocate (ystore_selected(ispecb6))
    allocate (zstore_selected(ispecb6))
    allocate(ind(ispecb6))
    allocate(ninseg(ispecb6))
    allocate(iglob(ispecb6))
    allocate(locval(ispecb6))
    allocate(ifseg(ispecb6))
    allocate(iwork(ispecb6))
    allocate(work(ispecb6))

    do ispec_tmp=1,ispecb6
      xstore_selected(ispec_tmp) = xstore(1,1,1,ibelm_top(ispec_tmp))
      ystore_selected(ispec_tmp) = ystore(1,1,1,ibelm_top(ispec_tmp))
      zstore_selected(ispec_tmp) = zstore(1,1,1,ibelm_top(ispec_tmp))
    enddo

    call sort_array_coordinates(ispecb6,xstore_selected,ystore_selected,zstore_selected, &
            ibelm_top,iglob,locval,ifseg,dummy_var,ind,ninseg,iwork,work)

    deallocate (xstore_selected)
    deallocate (ystore_selected)
    deallocate (zstore_selected)
    deallocate(ind)
    deallocate(ninseg)
    deallocate(iglob)
    deallocate(locval)
    deallocate(ifseg)
    deallocate(iwork)
    deallocate(work)
    !! test reord

  endif

  enddo

! check theoretical value of elements at the bottom
  if(ispecb5 /= NSPEC2D_BOTTOM) then
    call exit_MPI(myrank,'ispecb5 should equal NSPEC2D_BOTTOM')
  endif

! check theoretical value of elements at the top
  if(ispecb6 /= NSPEC2D_TOP) call exit_MPI(myrank,'ispecb6 should equal NSPEC2D_TOP')

  nspec2D_xmin = ispecb1
  nspec2D_xmax = ispecb2
  nspec2D_ymin = ispecb3
  nspec2D_ymax = ispecb4

  end subroutine get_jacobian_boundaries

! -------------------------------------------------------

  subroutine compute_jacobian_2D(myrank,ispecb,xelm,yelm,zelm,dershape2D,jacobian2D,normal,NGLLA,NGLLB,NSPEC2DMAX_AB)

  implicit none

  include "constants.h"

! generic routine that accepts any polynomial degree in each direction

  integer ispecb,NGLLA,NGLLB,NSPEC2DMAX_AB,myrank

  double precision xelm(NGNOD2D),yelm(NGNOD2D),zelm(NGNOD2D)
  double precision dershape2D(NDIM2D,NGNOD2D,NGLLA,NGLLB)

  real(kind=CUSTOM_REAL) jacobian2D(NGLLA,NGLLB,NSPEC2DMAX_AB)
  real(kind=CUSTOM_REAL) normal(3,NGLLA,NGLLB,NSPEC2DMAX_AB)

  integer i,j,ia
  double precision xxi,xeta,yxi,yeta,zxi,zeta
  double precision unx,uny,unz,jacobian

  do j=1,NGLLB
    do i=1,NGLLA

    xxi=ZERO
    xeta=ZERO
    yxi=ZERO
    yeta=ZERO
    zxi=ZERO
    zeta=ZERO
    do ia=1,NGNOD2D
      xxi=xxi+dershape2D(1,ia,i,j)*xelm(ia)
      xeta=xeta+dershape2D(2,ia,i,j)*xelm(ia)
      yxi=yxi+dershape2D(1,ia,i,j)*yelm(ia)
      yeta=yeta+dershape2D(2,ia,i,j)*yelm(ia)
      zxi=zxi+dershape2D(1,ia,i,j)*zelm(ia)
      zeta=zeta+dershape2D(2,ia,i,j)*zelm(ia)
    enddo

!   calculate the unnormalized normal to the boundary
    unx=yxi*zeta-yeta*zxi
    uny=zxi*xeta-zeta*xxi
    unz=xxi*yeta-xeta*yxi
    jacobian=dsqrt(unx**2+uny**2+unz**2)
    if(jacobian == ZERO) call exit_MPI(myrank,'2D Jacobian undefined')

!   normalize normal vector and store surface jacobian

! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      jacobian2D(i,j,ispecb)=sngl(jacobian)
      normal(1,i,j,ispecb)=sngl(unx/jacobian)
      normal(2,i,j,ispecb)=sngl(uny/jacobian)
      normal(3,i,j,ispecb)=sngl(unz/jacobian)
    else
      jacobian2D(i,j,ispecb)=jacobian
      normal(1,i,j,ispecb)=unx/jacobian
      normal(2,i,j,ispecb)=uny/jacobian
      normal(3,i,j,ispecb)=unz/jacobian
    endif

    enddo
  enddo

  end subroutine compute_jacobian_2D


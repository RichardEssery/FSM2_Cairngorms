!-----------------------------------------------------------------------
! Read a (lon,lat,time) netCDF data file
!-----------------------------------------------------------------------
subroutine READNC(Ncols,Nrows,Ntime,varname,var)

use netcdf

implicit none

character(len=*), intent(in) :: &
  varname             ! Variable name

integer, intent(in) :: &
  Ncols,             &! Number of columns in grid
  Nrows,             &! Number of rows in grid 
  Ntime               ! Number of timesteps

real, intent(out) :: &
  var(Ncols,Nrows,Ntime) ! Variable on grid

integer :: &
  ncid,              &! NetCDF dataset ID
  varid,             &! NetCDF variable ID
  status              ! NetCDF error status

status = nf90_open(varname//'.nc', nf90_nowrite, ncid)
status = nf90_inq_varid(ncid, varname, varid)
status = nf90_get_var(ncid, varid, var)
status = nf90_close(ncid)

end subroutine READNC

!-----------------------------------------------------------------------
! Read and flatten a (lon,lat,time) netCDF data file
!-----------------------------------------------------------------------
subroutine READNC_LAND(lat0,lon0,Nlat,Nlon,Nland,Ntime,ilat,ilon,      &
                       varname,var)

use netcdf

implicit none

character(len=*), intent(in) :: &
  varname             ! Variable name

integer, intent(in) :: &
  lat0,              &! First latitude point
  lon0,              &! First longitude point
  Nlat,              &! Number of latitude points
  Nlon,              &! Number of longitude points
  Nland,             &! Number of land points
  Ntime,             &! Number of timesteps
  ilat(Nland),       &! Latitude indices of land points
  ilon(Nland)         ! Longitude indices of land points

real, intent(out) :: &
  var(Nland,Ntime)    ! Variable on land points

integer :: &
  l,                 &! Land point counter
  ncid,              &! NetCDF dataset ID
  varid,             &! NetCDF variable ID
  status              ! NetCDF error status

real :: &
  var2d(Nlon,Nlat,Ntime) ! Variable on grid

print*,'Reading ',varname
status = nf90_open(varname//'.nc',nf90_nowrite,ncid)
status = nf90_inq_varid(ncid,varname,varid)
status = nf90_get_var(ncid,varid,var2d,start=(/lon0,lat0,1/),  &
                                 count=(/Nlon,Nlat,Ntime/))
do l = 1, Nland
  var(l,:) = var2d(ilon(l),ilat(l),:)
end do
status = nf90_close(ncid)

end subroutine READNC_LAND


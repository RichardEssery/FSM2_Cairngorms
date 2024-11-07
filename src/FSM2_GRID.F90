!----------------------------------------------------------------------!
! Flexible Snow Model (FSM version 2.1.0) with netCDF grid driving     !
!                                                                      !
! Richard Essery                                                       !
! School of GeoSciences                                                !
! University of Edinburgh                                              !
!----------------------------------------------------------------------!
program FSM2_GRID

#include "OPTS.h"

use netcdf

use CONSTANTS, only: &
  hcon_air,          &! Thermal conductivity of air (W/m/K)
  hcon_clay,         &! Thermal conductivity of clay (W/m/K)
  hcon_sand           ! Thermal conductivity of sand (W/m/K)

use LAYERS, only: &
  Dzsnow,            &! Minimum snow layer thicknesses (m)
  Dzsoil,            &! Soil layer thicknesses (m)
  fvg1,              &! Fraction of vegetation in upper canopy layer
  Ncnpy,             &! Number of canopy layers
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  zsub                ! Subcanopy wind speed diagnostic height (m)

use PARAMETERS, only: &
  fcly,              &! Soil clay fraction
  fsnd,              &! Soil sand fraction
  rgr0                ! Fresh snow grain radius (m)

use SOILPROPS, only: &
  b,                 &! Clapp-Hornberger exponent
  hcap_soil,         &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil,         &! Thermal conductivity of dry soil (W/m/K)
  sathh,             &! Saturated soil water pressure (m)
  Vcrit,             &! Volumetric soil moisture at critical point
  Vsat                ! Volumetric soil moisture at saturation

implicit none

! Grid dimensions
integer :: &
  Ncols,             &! Number of columns in grid
  Nrows               ! Number of rows in grid

! Site characteristics
real :: &
  lat,               &! Latitude (radians)
  noon                ! Time of solar noon (hour)

! Meteorological driving data
integer :: &
  year,              &! Year
  month,             &! Month of year
  day                 ! Day of month
real :: &
  dt,                &! Timestep (s)
  azim,              &! Solar azimuth (radians)
  elev,              &! Solar elevation (radians)
  hour,              &! Hour of day
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind speed measurement height (m)
real, allocatable :: &
  LW(:,:,:),         &! Incoming longwave radiation (W/m^2)
  Ps(:,:,:),         &! Surface pressure (Pa)
  Qa(:,:,:),         &! Specific humidity (kg/kg)
  Rf(:,:,:),         &! Rainfall rate (kg/m^2/s)
  Sdif(:,:,:),       &! Diffuse shortwave radiation (W/m^2)
  Sdir(:,:,:),       &! Direct-beam shortwave radiation (W/m^2)
  Sf(:,:,:),         &! Snowfall rate (kg/m^2/s)
  Shad(:,:,:),       &! Shadow fraction
  Ta(:,:,:),         &! Air temperature (K)
  trans(:,:),        &! Wind-blown snow transport rate (kg/m^2/s)
  Ua(:,:,:)           ! Wind speed (m/s)

! Model state variables
integer, allocatable :: &
  Nsnow(:,:)          ! Number of snow layers
real, allocatable :: &
  albs(:,:),         &! Snow albedo
  Tsrf(:,:),         &! Snow/ground surface temperature (K)
  Dsnw(:,:,:),       &! Snow layer thicknesses (m)
  Qcan(:,:,:),       &! Canopy air space humidities
  Rgrn(:,:,:),       &! Snow layer grain radii (m)
  Sice(:,:,:),       &! Ice content of snow layers (kg/m^2)
  Sliq(:,:,:),       &! Liquid content of snow layers (kg/m^2)
  Sveg(:,:,:),       &! Snow mass on vegetation layers (kg/m^2)
  Tcan(:,:,:),       &! Canopy air space temperatures (K)
  Tsnow(:,:,:),      &! Snow layer temperatures (K)
  Tsoil(:,:,:),      &! Soil layer temperatures (K)
  Tveg(:,:,:),       &! Vegetation layer temperatures (K)
  Vsmc(:,:,:)         ! Volumetric moisture content of soil layers
logical :: &
  start_file          ! True if start file exists

! Diagnostics
real, allocatable :: &
  fsnow(:,:),        &! Ground snowcover fraction
  H(:,:),            &! Sensible heat flux to the atmosphere (W/m^2)
  LE(:,:),           &! Latent heat flux to the atmosphere (W/m^2)
  LWout(:,:),        &! Outgoing LW radiation (W/m^2)
  LWsub(:,:),        &! Subcanopy downward LW radiation (W/m^2)
  Melt(:,:),         &! Surface melt rate (kg/m^2/s)
  Roff(:,:),         &! Runoff from snow (kg/m^2/s)
  snd(:,:),          &! Snow depth (m)
  snw(:,:),          &! Total snow mass on ground (kg/m^2) 
  subl(:,:),         &! Sublimation rate (kg/m^2/s)
  svg(:,:),          &! Total snow mass on vegetation (kg/m^2)
  SWincum(:,:),      &! Cumulated incoming SW radiation (W/m^2)  
  SWout(:,:),        &! Outgoing SW radiation (W/m^2)
  SWsub(:,:),        &! Subcanopy downward SW radiation (W/m^2)
  SWoutcum(:,:),     &! Cumulated outgoing SW radiation (W/m^2)  
  Tsub(:,:),         &! Subcanopy air temperature (K)
  Usub(:,:),         &! Subcanopy wind speed (m/s)
  Wflx(:,:,:)         ! Water flux into snow layer (kg/m^2/s)

! Topography
real, allocatable :: &
  aspct(:,:),        &! Slope aspect (radians)
  fsky(:,:),         &! Sky view fraction
  slope(:,:)          ! Slope angle (radians)

! Vegetation characteristics
real, allocatable :: &
  alb0(:,:),         &! Snow-free ground albedo
  vegh(:,:),         &! Canopy height (m)
  VAI(:,:)            ! Vegetation area index

! NetCDF variables
character(len=33) :: &
  units               ! Time units
integer :: &
  dimids(3),         &! x, y, t dimension IDs
  grid,              &! Latlon (0) or UTM (1) grid
  Ntime,             &! Number of timesteps
  ncid,              &! Dataset ID
  rec,               &! Record number
  status,            &! Error status
  time,time0,        &! Time (hours)
  varid,             &! Input variable ID
  varids(8)           ! Output variable IDs
real, allocatable :: &
  x(:),              &! x coordinates
  y(:)                ! y coordinates

! Counters
integer :: &
  i,                 &! Grid row counter
  j,                 &! Grid column counter
  n                   ! Timestep counter   

namelist /drive/ lat,noon

call FSM2_PARAMS

! Canopy, snow and soil layers
fvg1 = 0.5
zsub = 1.5
#if CANMOD == 1
Ncnpy = 1
#endif
#if CANMOD == 2
Ncnpy = 2
#endif
Nsmax = 3
Nsoil = 4
allocate(Dzsnow(Nsmax))
allocate(Dzsoil(Nsoil))
Dzsnow = (/0.1, 0.2, 0.4/)
Dzsoil = (/0.1, 0.2, 0.4, 0.8/)

! Grid dimensions and vegetation characteristics
status = nf90_open('ancil.nc', nf90_nowrite, ncid)
grid = 0
status = nf90_inq_dimid(ncid,'lon',dimids(1))
status = nf90_inq_dimid(ncid,'lat',dimids(2))
if (status /= nf90_noerr) then
  grid = 1
  status = nf90_inq_dimid(ncid,'x',dimids(1))
  status = nf90_inq_dimid(ncid,'y',dimids(2))
end if
status = nf90_inquire_dimension(ncid,dimids(1),len=Ncols)
status = nf90_inquire_dimension(ncid,dimids(2),len=Nrows)
allocate(x(Ncols))
allocate(y(Nrows))
if (grid == 0) then
  status = nf90_inq_varid(ncid,'lon',varid)
  status = nf90_get_var(ncid,varid,x)
  status = nf90_inq_varid(ncid,'lat',varid)
  status = nf90_get_var(ncid,varid,y)
else
  status = nf90_inq_varid(ncid,'x',varid)
  status = nf90_get_var(ncid,varid,x)
  status = nf90_inq_varid(ncid,'y',varid)
  status = nf90_get_var(ncid,varid,y)
end if
allocate(alb0(Ncols,Nrows))
allocate(VAI(Ncols,Nrows))
allocate(vegh(Ncols,Nrows))
status = nf90_inq_varid(ncid,'alb0',varid)
status = nf90_get_var(ncid, varid,alb0)
status = nf90_inq_varid(ncid,'VAI',varid)
status = nf90_get_var(ncid,varid,VAI)
status = nf90_inq_varid(ncid,'vegh',varid)
status = nf90_get_var(ncid,varid,vegh)
status = nf90_close(ncid)

! Soil properties
b = 3.1 + 15.7*fcly - 0.3*fsnd
hcap_soil = (2.128*fcly + 2.385*fsnd)*1e6 / (fcly + fsnd)
sathh = 10**(0.17 - 0.63*fcly - 1.58*fsnd)
Vsat = 0.505 - 0.037*fcly - 0.142*fsnd
Vcrit = Vsat*(sathh/3.364)**(1/b)
hcon_soil = (hcon_air**Vsat) * ((hcon_clay**fcly)*(hcon_sand**(1 - fcly))**(1 - Vsat))

! Allocate and initialize state variables
allocate(albs(Ncols,Nrows))
allocate(Nsnow(Ncols,Nrows))
allocate(Tsrf(Ncols,Nrows))
allocate(Dsnw(Nsmax,Ncols,Nrows))
allocate(Qcan(Ncnpy,Ncols,Nrows))
allocate(Rgrn(Nsmax,Ncols,Nrows))
allocate(Sice(Nsmax,Ncols,Nrows))
allocate(Sliq(Nsmax,Ncols,Nrows))
allocate(Sveg(Ncnpy,Ncols,Nrows))
allocate(Tcan(Ncnpy,Ncols,Nrows))
allocate(Tsnow(Nsmax,Ncols,Nrows))
allocate(Tsoil(Nsoil,Ncols,Nrows))
allocate(Tveg(Ncnpy,Ncols,Nrows))
allocate(Vsmc(Nsoil,Ncols,Nrows))
inquire(file='start.txt',exist=start_file)
if (start_file) then  ! Initialize from file start if it exists
  open(8,file='start.txt')
  read(8,*) albs
  read(8,*) Dsnw
  read(8,*) Nsnow
  read(8,*) Qcan
  read(8,*) Rgrn
  read(8,*) Sice
  read(8,*) Sliq
  read(8,*) Sveg
  read(8,*) Tcan
  read(8,*) Tsnow
  read(8,*) Tsoil
  read(8,*) Tsrf
  read(8,*) Tveg
  read(8,*) Vsmc
  close(8)
else                  ! Cold start
  albs = 0.8
  Dsnw = 0
  Nsnow = 0
  Qcan = 0
  Rgrn = rgr0
  Sice = 0
  Sliq = 0
  Sveg = 0
  Tcan = 285
  Tsnow = 273
  Tsoil = 285
  Tsrf = 285
  Tveg = 285
  Vsmc = 0.5*Vsat
end if

! Allocate diagnostic output arrays
allocate(fsnow(Ncols,Nrows))
allocate(H(Ncols,Nrows))
allocate(LE(Ncols,Nrows))
allocate(LWout(Ncols,Nrows))
allocate(LWsub(Ncols,Nrows))
allocate(Melt(Ncols,Nrows))
allocate(Roff(Ncols,Nrows))
allocate(snd(Ncols,Nrows))
allocate(snw(Ncols,Nrows))
allocate(subl(Ncols,Nrows))
allocate(svg(Ncols,Nrows))
allocate(SWout(Ncols,Nrows))
allocate(SWsub(Ncols,Nrows))
allocate(Tsub(Ncols,Nrows))
allocate(Usub(Ncols,Nrows))
allocate(Wflx(Nsmax,Ncols,Nrows))

! Driving data characteristics
dt = 3600
zT = 2
zU = 10
status = nf90_open('LWdown.nc', nf90_nowrite, ncid)
status = nf90_inq_dimid(ncid,'time',dimids(3))
status = nf90_inquire_dimension(ncid,dimids(3),len=Ntime)
status = nf90_inq_varid(ncid,'time',varids(3))
status = nf90_get_var(ncid,varids(3),time0,start=(/1/))
status = nf90_get_att(ncid,varids(3),'units',units)
read(units(13:16),*) year
read(units(18:19),*) month
read(units(21:22),*) day

! Allocate and read driving data
allocate(LW(Ncols,Nrows,Ntime))
allocate(Ps(Ncols,Nrows,Ntime))
allocate(Qa(Ncols,Nrows,Ntime))
allocate(Rf(Ncols,Nrows,Ntime))
allocate(Sf(Ncols,Nrows,Ntime))
allocate(Sdif(Ncols,Nrows,Ntime))
allocate(Sdir(Ncols,Nrows,Ntime))
allocate(Shad(Ncols,Nrows,Ntime))
allocate(Ta(Ncols,Nrows,Ntime))
allocate(Ua(Ncols,Nrows,Ntime))
allocate(trans(Ncols,Nrows))
call READNC(Ncols,Nrows,Ntime,'LWdown',LW)
call READNC(Ncols,Nrows,Ntime,'Psurf',Ps)
call READNC(Ncols,Nrows,Ntime,'Qair',Qa)
call READNC(Ncols,Nrows,Ntime,'Rainf',Rf)
call READNC(Ncols,Nrows,Ntime,'Snowf',Sf)
#if SLOPES == 1
call READNC(Ncols,Nrows,Ntime,'SWdif',Sdif)
call READNC(Ncols,Nrows,Ntime,'SWdir',Sdir)
#else
call READNC(Ncols,Nrows,Ntime,'SWdown',Sdif)
Sdir(:,:,:) = 0
#endif
call READNC(Ncols,Nrows,Ntime,'Tair',Ta)
call READNC(Ncols,Nrows,Ntime,'Wind',Ua)
trans(:,:) = 0

#if SLOPES == 1
! Topography
status = nf90_open('topog.nc', nf90_nowrite, ncid)
allocate(aspct(Ncols,Nrows))
allocate(fsky(Ncols,Nrows))
allocate(slope(Ncols,Nrows))
status = nf90_inq_varid(ncid,'aspect',varid)
status = nf90_get_var(ncid,varid,aspct)
status = nf90_inq_varid(ncid,'sky',varid)
status = nf90_get_var(ncid,varid,fsky)
status = nf90_inq_varid(ncid,'slope',varid)
status = nf90_get_var(ncid,varid,slope)
status = nf90_close(ncid)
call READNC(Ncols,Nrows,Ntime,'shadowfrac',Shad)
read(5,drive)
lat = (3.14159/180)*lat
#endif

! Prepare NetCDF output file
status = nf90_create('FSM2out.nc',nf90_clobber,ncid)
if (grid == 0) then
  status = nf90_def_dim(ncid,'lon',Ncols,dimids(1))
  status = nf90_def_var(ncid,'lon',NF90_REAL,dimids(1),varids(1))
  status = nf90_def_dim(ncid,'lat',Nrows,dimids(2))
  status = nf90_def_var(ncid,'lat',NF90_REAL,dimids(2),varids(2))
else
  status = nf90_def_dim(ncid,'x',Ncols,dimids(1))
  status = nf90_def_var(ncid,'x',NF90_REAL,dimids(1),varids(1))
  status = nf90_def_dim(ncid,'y',Nrows,dimids(2))
  status = nf90_def_var(ncid,'y',NF90_REAL,dimids(2),varids(2))
end if
status = nf90_def_dim(ncid,'time',NF90_UNLIMITED,dimids(3))
status = nf90_def_var(ncid,'time',NF90_REAL,dimids(3),varids(3))
status = nf90_def_var(ncid,'snc',NF90_REAL,dimids,varids(4))
status = nf90_def_var(ncid,'snw',NF90_REAL,dimids,varids(5))
status = nf90_def_var(ncid,'alb',NF90_REAL,dimids,varids(6))
status = nf90_def_var(ncid,'rgrn',NF90_REAL,dimids,varids(7))
status = nf90_def_var(ncid,'lwc',NF90_REAL,dimids,varids(8))
status = nf90_put_att(ncid,varids(3),'units',units)
status = nf90_put_att(ncid,varids(4),'long_name','snowcover fraction')
status = nf90_put_att(ncid,varids(4),'units','-')
status = nf90_put_att(ncid,varids(5),'long_name','mass of snowpack')
status = nf90_put_att(ncid,varids(5),'units','kg m-2')
status = nf90_put_att(ncid,varids(6),'long_name','albedo')
status = nf90_put_att(ncid,varids(6),'units','-')
status = nf90_put_att(ncid,varids(7),'long_name','snow grain size')
status = nf90_put_att(ncid,varids(7),'units','m')
status = nf90_put_att(ncid,varids(8),'long_name','liquid water content of snow')
status = nf90_put_att(ncid,varids(8),'units','kg m-2')
status = nf90_enddef(ncid)
status = nf90_put_var(ncid,varids(1),x)
status = nf90_put_var(ncid,varids(2),y) 

! Run the model and write daily output
rec = 1
allocate(SWincum(Ncols,Nrows))
allocate(SWoutcum(Ncols,Nrows))
SWincum(:,:) = 0
SWoutcum(:,:) = 0
do n = 1, Ntime
  write(6,fmt='(a,i3,a)',advance='no') achar(13),int(100*n/float(Ntime)),'% complete'
  hour = mod(n,24) 
  time = time0 + n 
#if SLOPES == 1
  call SOLARPOS(year,month,day,hour,lat,noon,azim,elev)
  if (elev <= 0) then
    Sdif(:,:,n) = Sdif(:,:,n) + Sdir(:,:,n) 
    Sdir(:,:,n) = 0
  end if
#endif
  do i = 1, Nrows
  do j = 1, Ncols
#if SLOPES == 1
    Sdif(j,i,n) = fsky(j,i)*Sdif(j,i,n)
    Sdir(j,i,n) = shad(j,i,n)*Sdir(j,i,n)*(cos(slope(j,i)) +  & 
                       sin(slope(j,i))*cos(azim - aspct(j,i))/tan(elev))
#endif
    call FSM2_TIMESTEP(                                                &
                       ! Driving variables                             &
                       dt,elev,zT,zU,LW(j,i,n),Ps(j,i,n),Qa(j,i,n),    &
                       Rf(j,i,n),Sdif(j,i,n),Sdir(j,i,n),Sf(j,i,n),    &
                       Ta(j,i,n),trans(j,i),Ua(j,i,n),                 &
                       ! Vegetation characteristics                    &
                       alb0(j,i),vegh(j,i),VAI(j,i),                   &
                       ! State variables                               &
                       albs(j,i),Tsrf(j,i),Dsnw(:,j,i),Nsnow(j,i),     &
                       Qcan(:,j,i),Rgrn(:,j,i),Sice(:,j,i),            &
                       Sliq(:,j,i),Sveg(:,j,i),Tcan(:,j,i),            &
                       Tsnow(:,j,i),Tsoil(:,j,i),Tveg(:,j,i),          &
                       Vsmc(:,j,i),                                    &
                       ! Diagnostics                                   &
                       fsnow(j,i),H(j,i),LE(j,i),LWout(j,i),           &
                       LWsub(j,i),Melt(j,i),Roff(j,i),snd(j,i),        &
                       snw(j,i),subl(j,i),svg(j,i),SWout(j,i),         &
                       SWsub(j,i),Tsub(j,i),Usub(j,i),Wflx(:,j,i)      )  
  end do
  end do
  SWincum(:,:) = SWincum(:,:) + Sdif(:,:,n) + Sdir(:,:,n)
  SWoutcum(:,:) = SWoutcum(:,:) + SWout(:,:)
  if (hour==12) then
    status = nf90_put_var(ncid,varids(3),time,start=(/rec/))
    status = nf90_put_var(ncid,varids(4),fsnow,start=(/1,1,rec/),  &
                          count=(/Ncols,Nrows,1/))
    status = nf90_put_var(ncid,varids(5),snw,start=(/1,1,rec/),    &
                          count=(/Ncols,Nrows,1/))
    status = nf90_put_var(ncid,varids(6),SWoutcum/SWincum,start=(/1,1,rec/),    &
                          count=(/Ncols,Nrows,1/))
    status = nf90_put_var(ncid,varids(7),Rgrn(1,:,:),start=(/1,1,rec/),    &
                          count=(/Ncols,Nrows,1/))
    status = nf90_put_var(ncid,varids(8),Sliq(1,:,:),start=(/1,1,rec/),    &
                          count=(/Ncols,Nrows,1/))
    SWincum(:,:) = 0
    SWoutcum(:,:) = 0
    rec = rec + 1
  end if
end do
status=nf90_close(ncid) 
print*,minval(snw),maxval(snw)
write(6,*)

! Write out state variables at end of run
open(8,file='dump.txt')
write(8,*) albs
write(8,*) Dsnw
write(8,*) Nsnow
write(8,*) Qcan
write(8,*) Rgrn
write(8,*) Sice
write(8,*) Sliq
write(8,*) Sveg
write(8,*) Tcan
write(8,*) Tsnow
write(8,*) Tsoil
write(8,*) Tsrf
write(8,*) Tveg
write(8,*) Vsmc
close(8)

end program FSM2_GRID

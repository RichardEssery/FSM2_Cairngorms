########################################################################
# Compile and run FSM2 with UKV or WFDE5 driving
#
# Richard Essery
# School of GeoSciences
# University of Edinburgh
########################################################################
WDIR=`pwd`
NDIR=/usr/include/
SDIR=src
METM=WFDE5
MDIR=$METM'_Cairngorms/'
ODIR=$METM'_Output/'

mkdir -p $ODIR
cd $SDIR

cat > OPTS.h << EOF
/* Process options                                  : Possible values */
#define ALBEDO 2   /* snow albedo                   : 1, 2            */
#define CANINT 1   /* canopy interception of snow   : 1, 2            */
#define CANMOD 1   /* forest canopy layers          : 1, 2            */
#define CANRAD 1   /* canopy radiative properties   : 1, 2            */
#define CANUNL 1   /* unloading of canopy           : 1, 2            */
#define CONDCT 1   /* snow thermal conductivity     : 0, 1            */
#define DENSTY 1   /* snow density                  : 0, 1, 2         */
#define EXCHNG 1   /* turbulent exchange            : 0, 1            */
#define HYDROL 1   /* snow hydraulics               : 0, 1, 2         */
#define SGRAIN 1   /* snow grain growth             : 1, 2            */
#define SNFRAC 1   /* snow cover fraction           : 1, 2            */
/* Driving data options                             : Possible values */
#define SLOPES 0   /* slope shading                 : 0, 1            */
#define ZOFFST 1   /* measurement height offset     : 0, 1            */
EOF
cp OPTS.h $WDIR/$ODIR'/OPTS.h'

echo 'Compiling FSM2'
gfortran -cpp -O3 -o $WDIR/FSM2 -I$NDIR FSM2_MODULES.F90               \
FSM2_PARAMS.F90 FSM2_GRID.F90 FSM2_TIMESTEP.F90 READNC.F90             \
CANOPY.F90 INTERCEPT.F90 LUDCMP.F90 PSIMH.F90 QSAT.F90 SNOW.F90        \
SOIL.F90 SOLARPOS.F90 SRFEBAL.F90 SWRAD.F90 THERMAL.F90 TRIDIAG.F90    \
TWOSTREAM.F90 -lnetcdff -lnetcdf 
rm *.mod
cd $WDIR

rm -f start.txt
cp $MDIR'ancil.nc' ancil.nc
years=(2016 2017 2018 2019 2020 2021)
months=(01 02 03 04 05 06 07 08 09 10 11 12)
vars=(LWdown Psurf Qair Rainf Snowf SWdown Tair Wind)
for y in ${years[*]}
do
  for m in ${months[*]}
  do
    if test -f $MDIR'LWdown'$y$m'.nc'; then
      echo 'Running FSM2 open for' $y $m
      for v in ${vars[*]}
      do
        ln -s -f $MDIR$v$y$m'.nc' $v.nc
      done
      time ./FSM2 < nlst
      mv dump.txt start.txt
      mv FSM2out.nc $ODIR'FSM2out'$y$m'.nc'
    fi
  done
done
rm -f start.txt
rm *.nc


#!/bin/sh

homedir=/home/mbueti/wind_stress

name=$1
startdate=$2
stormid=$3

yyyy=`echo $startdate | cut -b 1-4`

echo "Running WindOnly for $name starting $startdate"
cd $homedir/output

grep ' '$stormid $homedir/input/vitals/syndat_tcvitals.$yyyy | $homedir/scripts/gfdl_pre_sortvit.sh > track
echo "000  000 000000    000000 0000 000  000 000  00  000 0000  000 00  00 0000 0000 0000 0000 0000 0000 0000 0000" >> track

echo "Track PreSorting finished. Running WindOnly"

cp $homedir/input/parameters.${name}.${startdate} parameters.inp

ln -s -f  parameters.inp       fort.10
ln -s -f  track                fort.15

$homedir/exe/windonly > windonly.out
rm fort.*

echo "WindOnly is finished. Packing output."

tar -czf $name.$startdate.tar.gz TXY.* WSURF.* windonly.out parameters.inp track
rm -rf TXY.* WSURF.* windonly.out parameters.inp track

echo "Output packed and available in $name.$startdate.tar.gz"

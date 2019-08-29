#!/bin/sh
# the script below cleans up log.lammps to get thermo out
file=outfiles/log.lammps
if test -s $file
then
awk '
/Step/,/Loop/ { print }
' $file > x.dat
awk '
!($0 ~ /Step/) { print }
' x.dat > y.dat
awk '
!($0 ~ /Loop/) { print }
' y.dat > outfiles/thermo.dat
rm x.dat
rm y.dat
echo "Generated thermo.dat"
echo "Done!"
else
echo "where is log.lammps?"
fi

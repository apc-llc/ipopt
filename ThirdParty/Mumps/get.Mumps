#!/bin/sh

# Set the following to the latest MUMPS version.
#  THERE MUST BE NO SPACE BEFORE AND AFTER THE EQUAL (=) OPERATOR.
mumps_ver=4.10.0

set -e

wgetcmd=wget
wgetcount=`which wget 2>/dev/null | wc -l`
if test ! $wgetcount = 1; then
  echo "Utility wget not found in your PATH."
  if test `uname` = Darwin; then
    wgetcmd=ftp
    echo "Using ftp command instead."
  elif test `uname` = FreeBSD; then
    wgetcmd=fetch
    echo "Using fetch command instead."
  else
    exit -1
  fi
fi

echo " "
echo "Running script for preparing the source code for MUMPS"
echo " "

if [ ! -e mumps.fetch ]; then
  echo "Downloading the source code from ..."
  if ! ( $wgetcmd http://mumps.enseeiht.fr/MUMPS_${mumps_ver}.tar.gz && touch mumps.fetch ) ;
  then
    exit 1;
  fi
fi

echo "Unpacking the source code..."
tar -xf MUMPS_${mumps_ver}.tar.gz

rm -rf MUMPS
mv MUMPS_${mumps_ver} MUMPS

echo "Apply a patch to work around a bug in MUMPS."
echo " "

patch -p0 < mumps.patch

echo " "
echo "Verify that there are no error message in the output above."
echo " "
echo "Done preparing the source code for MUMPS."
echo " "


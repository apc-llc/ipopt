#!/bin/sh

set -e

wgetcmd="wget"
wgetcount=`which wget 2>/dev/null | wc -l`
if test ! $wgetcount = 1; then
  echo "Utility wget not found in your PATH."
  if test `uname` = Darwin; then
    wgetcmd="curl -L -k -O"
    echo "Using ftp command instead."
  elif test `uname` = FreeBSD; then
    wgetcmd=fetch
    echo "Using fetch command instead."
  else
    exit -1
  fi
fi

coinasl=solvers-20130815.tgz

echo " "
echo "Running script for preparing the source code for the ASL"
echo " "

if [ ! -e asl.fetch ]; then
  echo "Downloading the source code from projects.coin-or.org..."
  if $wgetcmd http://www.coin-or.org/BuildTools/ASL/$coinasl && touch asl.fetch ;
  then
    mv $coinasl solvers.tgz
  else
    echo
    echo "Failed, try downloading the source code from www.ampl.com..."
    if ! ( $wgetcmd http://www.ampl.com/netlib/ampl/solvers.tgz && touch asl.fetch ) ;
    then
      exit 1;
    fi
  fi
fi

rm -rf solvers

echo "Unpacking the source code..."
tar -xf solvers.tgz

echo "Adding No_dtoa to CFLAGS..."
cd solvers
mv makefile.u makefile.u.orig
sed -e 's/CFLAGS = /CFLAGS = -DNo_dtoa /g' makefile.u.orig > makefile.u 
cd ..

echo " "
echo "Done preparing the source code for ASL."
echo " "


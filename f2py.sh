#!/bin/bash

FILES='read_ramses_fortran.f90'
F2PY=f2py
FORT=gfortran
BASEDIR=$(dirname "$0")

for f in $FILES
do
    bn=$(basename "$f" .f90)
    if [[ $FORT == "gfortran" ]]; then
        # OMP on
        $FORT -x f95-cpp-input -c -fopenmp $f
        $F2PY -lgomp --f90exec=$FORT --f77exec=$FORT --f90flags='-fopenmp -O3 -x f95-cpp-input'  -c $f -m $bn
        # OMP off
        # $FORT -x f95-cpp-input -c $f
        # $F2PY -c -lgomp --f90exec=$FORT --f77exec=$FORT $f -m $bn --opt='-O3 -x f95-cpp-input'
    fi

    if [[ $FORT == "ifort" ]]; then
        # OMP on
        $FORT -O3 -foptimize-sibling-calls -c $f
        $F2PY -c $f -m $bn --fcompiler=intelem --opt='-O3 -heap-arrays 500000000 -foptimize-sibling-calls -fpp -m64 -free -fopenmp' -liomp5
        # OMP off
        # $FORT -O3 -foptimize-sibling-calls -c $f
        # $F2PY -c $f -m $bn --fcompiler=intelem --opt='-O3 -heap-arrays -foptimize-sibling-calls -fpp -m64 -free' -liomp5
    fi
done
rm *.o *.mod

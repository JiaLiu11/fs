#! /usr/bin/env bash

# C compiler
CC=`which g++`; CFLAGS=" -O3 -Wall -lgsl -lgslcblas"


#$CC Freestreaming.cpp testfs.cpp arsenal.cpp mistools.cpp gauss_quadrature.cpp -o fs.e  $CFLAGS
$CC LdMatching.cpp Freestreaming.cpp CellData.cpp lmtest.cpp EOS.cpp arsenal.cpp mistools.cpp gauss_quadrature.cpp -o lm.e  $CFLAGS

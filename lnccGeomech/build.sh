#!/bin/bash

obj_folder="/home/rsilva/projects/lnccProgram/bin"
filename="lnccGeomech"
fullpath=$obj_folder/$filename

#  set "flags_DEBUG=/c /warn:all /nologo /Z7 /Od /traceback /gen-interfaces /warn:interfaces /check /Qtrapuv /fpe:0"
flags_DEBUG="-c -Wall -Wconversion -g -fbacktrace -fbounds-check -ffpe-trap=zero,overflow -ffree-line-length-512 -Og"

if  not exist $obj_folder mkdir -p $obj_folder
pushd $obj_folder
gfortran $flags_DEBUG ..\fonte\variaveisGlobais.F90
gfortran $flags_DEBUG ..\fonte\malha.F90
gfortran $flags_DEBUG ..\fonte\leituraEscrita.F90
gfortran $flags_DEBUG ..\fonte\funcoesDeForma.F90
gfortran $flags_DEBUG ..\fonte\utilitarios.F90
gfortran $flags_DEBUG ..\fonte\mInputReader.F90
gfortran $flags_DEBUG ..\fonte\propGeoFisica.F90
gfortran $flags_DEBUG ..\fonte\solverGaussSkyline.F90
gfortran $flags_DEBUG ..\fonte\leituraEscritaSimHidroGeoMec.F90
gfortran $flags_DEBUG ..\fonte\geomecanica.F90
gfortran $flags_DEBUG ..\fonte\hidrodinamica.F90
gfortran $flags_DEBUG ..\fonte\driverGeocreep.F90

gfortran -g variaveisGlobais.o malha.o leituraEscrita.o funcoesDeForma.o utilitarios.o mInputReader.o propGeoFisica.o solverGaussSkyline.o leituraEscritaSimHidroGeoMec.o geomecanica.o hidrodinamica.o driverGeocreep.o -o $fullpath
popd

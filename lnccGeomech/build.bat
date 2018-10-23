@echo off

if not defined vsOn (
    call "C:\Program Files (x86)\Intel\Composer XE 2015\bin\compilervars.bat" intel64 vs2010
   set vsOn=true
)

setlocal
SET me=%~n0

set obj_folder="lnccGeomech\bin"

set "flags_all=/diag-disable"
set "flags=/c /warn:all /nologo"
set "flags_DEBUG=/Z7 /Od /traceback /warn:interfaces /check /fpe:0"

if not exist %obj_folder% mkdir %obj_folder%
pushd %obj_folder%
ifort %flags% %flags_DEBUG% ..\fonte\variaveisGlobais.F90
ifort %flags% %flags_DEBUG% ..\fonte\malha.F90
ifort %flags% %flags_DEBUG% ..\fonte\leituraEscrita.F90
ifort %flags% %flags_DEBUG% ..\fonte\funcoesDeForma.F90
ifort %flags% %flags_DEBUG% ..\fonte\utilitarios.F90
ifort %flags% %flags_DEBUG% ..\fonte\mInputReader.F90
ifort %flags% %flags_DEBUG% ..\fonte\propGeoFisica.F90
ifort %flags% %flags_DEBUG% ..\fonte\solverGaussSkyline.F90
ifort %flags% %flags_DEBUG% ..\fonte\leituraEscritaSimHidroGeoMec.F90
ifort %flags% %flags_DEBUG% ..\fonte\geomecanica.F90
ifort %flags% %flags_DEBUG% ..\fonte\hidrodinamica.F90
ifort %flags% %flags_DEBUG% ..\fonte\driverGeocreep.F90

ifort /nologo /debug:all /exe:lnccGeomech.exe variaveisGlobais.obj malha.obj leituraEscrita.obj funcoesDeForma.obj utilitarios.obj mInputReader.obj propGeoFisica.obj solverGaussSkyline.obj leituraEscritaSimHidroGeoMec.obj geomecanica.obj hidrodinamica.obj driverGeocreep.obj
popd

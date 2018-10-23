@echo off

if not defined vsOn (
    call "C:\Program Files (x86)\Intel\Composer XE 2015\bin\compilervars.bat" intel64 vs2010
   set vsOn=true
)

cd %cwd%

devenv /debugexe "..\..\bin\lnccGeomech.exe"

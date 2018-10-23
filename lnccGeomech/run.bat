@echo off

set "dirExp=data\%1"

pushd %dirExp%
start "" "..\..\bin\lnccGeomech.exe"
popd

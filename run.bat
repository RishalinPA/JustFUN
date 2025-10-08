@echo off
echo  Running ...
gfortran -O5 -fopenmp -o kolT kolT.f95  && .\kolT.exe
echo DONE

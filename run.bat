@echo off
echo  Running ...
gfortran -O5 -fopenmp -o kol kol.f95  && .\kol.exe
echo DONE
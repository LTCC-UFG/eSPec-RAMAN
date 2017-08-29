# Raman COMPILATION MAKEFILE
#
# Vinicius Vaz da Cruz
#

PROG=raman
CC=icc
FC=ifort -nofor_main
#FC=gfortran
#ifort flags
#FFLAGS=-traceback -fast
#gfortran flags
#FFLAGS=-O3
#CFLAGS=-O3 
#triolith flags
CFLAGS=-O0 -I/software/apps/fftw/3.3.2/i1214/include -L/software/apps/fftw/3.3.2/i1214/lib -lfftw3
FFLAGS=-traceback -O0 -mcmodel large -shared-intel -L/software/apps/fftw/3.3.2/i1214/lib -lfftw3 
#-check all -check noarg_temp_created -g

all: main

main: src/main.o src/rdinput.o src/readallspl.o src/splinesurf.o src/fourier.o
	$(FC) $(FFLAGS) -o $(PROG) src/main.o src/rdinput.o src/readallspl.o src/splinesurf.o src/fourier.o -lfftw3

main.o: src/main.c
	$(CC) $(CFLAGS) -c src/main.c

rdinput.o: src/rdinput.c
	$(CC) $(CFLAGS) -c src/rdinput.c

splinesurf.o: src/splinesurf.f src/splinesurf.h
	$(FC) $(FFLAGS) -c src/splinesurf.f

fourier.o: src/fourier.c src/fourier.h
	$(CC) $(CFLAGS) -c src/fourier.c 

elk:
	cp paral_espec_raman.sh paral_run_espec_raman_elk.sh
	sed -i 's/#%header2%//g' paral_run_espec_raman_elk.sh
	sed -i 's/#%header1%/#Elk environment/g' paral_run_espec_raman_elk.sh
	sed -i 's/#%header3%//g' paral_run_espec_raman_elk.sh
	sed -i 's/%especpath%/\/home\/vinicius\/programming\/eSPec\/espec_v07.x/g' paral_run_espec_raman_elk.sh
	sed -i 's/%ramanpath%/\/home\/vinicius\/programming\/eSPec-RAMAN\/raman/g' paral_run_espec_raman_elk.sh
	sed -i 's/%fcorrelpath%/\/home\/vinicius\/programming\/eSPec-RAMAN\/fcorrel\/correl/g' paral_run_espec_raman_elk.sh
	sed -i 's/%pythonpath%/\/home\/vinicius\/programming\/eSPec-RAMAN\/functions.py/g' paral_run_espec_raman_elk.sh	
	sed -i 's/%pyversion%/python/g' paral_run_espec_raman_elk.sh	
	chmod +x paral_run_espec_raman_elk.sh
triolith:
	cp paral_espec_raman.sh paral_run_espec_raman_triolith.sh
	sed -i 's/#%header2%/module load buildenv-intel\/2015-1/g' paral_run_espec_raman_triolith.sh
	sed -i 's/#%header1%/#Triolith environment:/g' paral_run_espec_raman_triolith.sh
	sed -i 's/#%header3%/export LD_LIBRARY_PATH=\/software\/apps\/intel\/composer_xe_2015.1.133\/compiler\/lib\/intel64/g' paral_run_espec_raman_triolith.sh
	sed -i 's/%especpath%/\/proj\/xramp2015\/progs\/eSPec-latest\/eSPec_v0.7\/espec_v07.x/g' paral_run_espec_raman_triolith.sh
	sed -i 's/%ramanpath%/\/proj\/xramp2015\/progs\/eSPec-latest\/eSPec-RAMAN\/raman/g' paral_run_espec_raman_triolith.sh
	sed -i 's/%fcorrelpath%/\/proj\/xramp2015\/progs\/eSPec-latest\/eSPec-RAMAN\/fcorrel\/correl/g' paral_run_espec_raman_triolith.sh
	sed -i 's/%pythonpath%/\/proj\/xramp2015\/progs\/eSPec-latest\/eSPec-RAMAN\/functions.py/g' paral_run_espec_raman_triolith.sh
	sed -i 's/#%triolithpython%/module add python\/3.5.2/g' paral_run_espec_raman_triolith.sh
	sed -i 's/%pyversion%/python3/g' paral_run_espec_raman_triolith.sh 
	chmod +x paral_run_espec_raman_triolith.sh

clean:
	rm *~ src/*~ src/*.o

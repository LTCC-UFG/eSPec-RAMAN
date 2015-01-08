# general makefile. MAKE SURE that you use GNU make.
# pawsa@theochem.kth.se
# NOTE: refer to gnumake.dvi in difficulties.
# distribution contains only F90 files, soft F file links are created 
# on demand.
system=IFC
MAIN=ram
sfx=${OSTYPE}-${MACHTYPE}
FC=f77
FFLAGS=-O

# IBM AIX SETTINGS *************************************************
ifeq ($(system),AIX)
F90=xlf90
FOPTS= -O3 -Q
# -WF,-DUSE_ESSL
#-O3 -Q -WF,-DUSE_BLAS
LOPTS=-L/r1/pawsa/lib -llapack -lspline -lblas -lfft
EXT=F
endif

# GNU FORTRAN SETTINGS *************************************************
ifeq ($(system),GNU)
F90=gfortran
FC=gfortran
FOPTS= -O3 
# -WF,-DUSE_ESSL
#-O3 -Q -WF,-DUSE_BLAS
LOPTS= -L. -llapack -lspline -lblas -lfft -static
EXT=F90
endif

# PDC Strindberg SETTINGS *******************************************
ifeq ($(system),PDC)
F90=mpxlf90
FOPTS= -O3 -Q -WF,-DUSE_BLAS
# -WF,-DUSE_ESSL
#-O3 -Q -WF,-DUSE_BLAS
LOPTS=-L/usr/local/numlib/lib -L$(HOME)/lib -llapackcom -lspline -lblas -lfft
EXT=F
MAIN=ramp
endif

# LINUX PG settings ************************************************
ifeq ($(system),PGI)
F90=pgf90
FOPTS= -fast -Mfree -Mdalign -Mnoframe
LOPTS=-L/home/pawsa/lib liblapack.a -lspline -lblas -lfft -lpgftnrtl -lg2c
EXT=F
endif

# LINUX PG/64 settings *********************************************
ifeq ($(system),PGI64)
F90=pgf90
FOPTS= -fastsse -Mvect=idiom -Mfree
LOPTS=-Lspline -lspline -lacml
EXT=F90
endif

# LINUX IFC settings ***********************************************
# Inter fortran compiler
ifeq ($(system),IFC)
FC =ifort
F90=ifort
FOPTS= -i_dynamic -O -pad -xW 
FFLAGS = $(FOPTS)
LOPTS=-xW -Lfft -Lspline -lspline -lfft -llapack -lblas
#,--rpath=/opt/intel_fc_80/lib
EXT=F90
endif

# LINUX PG settings ************************************************
ifeq ($(system),ABSOFT)
F90=f90
FOPTS= -O
LOPTS=-L/home/pawsa/lib -llapack -lspline -lblas -lfft
EXT=F90
endif

# LINUX PG settings ************************************************
ifeq ($(system),BEOWOLF)
F90=pgf90
FOPTS= -fast -Mfree -Mdalign -Mnoframe -I/usr/local/mpi/include
LOPTS=-L$(HOME)/lib -llapack -lspline -lblas -lfft -lpgftnrtl -lmpich \
	-L/usr/local/vol/mpich-1.1.2/build/LINUX/ch_p4/lib 
EXT=F
MAIN=ramp
endif

# CRAY T3E SETTINGS ************************************************
ifeq ($(system),CRAY)
F90=f90
FOPTS=-DSYS_CRAY -DUSE_BLAS -O -apad
#FOPTS= -O3 -apad -DSYS_CRAY -DUSE_BLAS
LOPTS=-L$(HOME)/lib -lspline
EXT=F90
MAIN=ramp
endif

# SUN SETTINGS *****************************************************
ifeq ($(system),NAG)
F90=f95
FOPTS=-O4 -target=native -Oblock=32 -ieee=nonstd -DUSE_BLAS
#-O4 -DUSE_BLAS
LOPTS=-L$(HOME)/lib -lspline -xlic_lib=sunperf -L/pkg/math/lib -lslatec
EXT=F90
endif

# SUN WORKSHOP 5 SETTINGS ******************************************
ifeq ($(system),SUN)
F90=f90
FOPTS=-g -DUSE_BLAS
#-fast -native -fns=yes 
# up to -O4
LOPTS=-L$(HOME)/lib -lspline -xlic_lib=sunperf  -L/pkg/math/lib  -lslatec
EXT=F90
endif

# ALPHA SETTINGS *************************************************
ifeq ($(system),ALPHA)
F90=f90
EXT=F90
FOPTS= -fast -free -DUSE_BLAS -pipeline
LOPTS=-L. -ldxml -lfft -lspline
endif

# ANGUS (Linux/Alpha) SETTINGS ***********************************
ifeq ($(system),ANGUS)
F90=fort
EXT=F90
FOPTS=  -DUSE_BLAS
LOPTS=-L$(HOME)/lib -lcxml -lcpml -lfft -lspline
endif

# GENERAL SETTINGS *************************************************
OBJSRAM= lapint.o lanczos.o ioproc.o pfsolv.o 

ram: $(OBJSRAM) $(MAIN).o libfft.a libspline.a
	$(F90) $^ -o $@ $(LOPTS)

dist: Ram.tgz
	@echo "Distribution placed in Ram.tgz"

Ram.tgz: $(OBJSRAM:.o=.F90) ram.F90 ramp.F90 Makefile
	@echo "Packing " $^
	rm -f $@
	tar cvf - $^ | gzip> $@

libfft.a:
	$(MAKE) -C fft all FC=$(FC) FFLAGS="$(FFLAGS)"

libspline.a:
	$(MAKE) -C spline all FC=$(FC) FFLAGS="$(FFLAGS)"

# pattern implicit rules
%.o: %.$(EXT)
	$(F90) -c $(FOPTS) $<

# below, the dependices are described. The dependices should refer rather to
# module files but the module files have different extensions or don't
# exist at all (as on Cray).

ram.o: ram.$(EXT) lapint.o ioproc.o pfsolv.o 

lapint.o: lapint.$(EXT)

lanczos.o: lanczos.$(EXT) lapint.o

ioproc.o: ioproc.$(EXT) lapint.o lanczos.o

pfsolv.o: pfsolv.$(EXT) lanczos.o ioproc.o lapint.o

ramp.o: ramp.$(EXT)

clean: 
	-rm -f *.o ram *.lst *.mod *.M work.pc* *.d *.a
	$(MAKE) -C fft    clean
	$(MAKE) -C spline clean

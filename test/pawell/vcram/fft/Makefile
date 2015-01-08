# the makefile for the scientific libriaries
# Pawel Salek, pawsa@theochem.kth.se

sfx=${OSTYPE}-${MACHTYPE}

AR= ar
FC    = ifort
FFLAGS= -O3 -xW
FFTOBJS= zfftb.o zfftf.o zffti.o zheev.o
OBJS=$(FFTOBJS) 

#.f.o:
#	$(F77) -c $(FFLAGS) $<

all: libfft.a 
	cp -p libfft.a ..

libfft.a : $(FFTOBJS)
	-$(RM) $@
	$(AR) cvq $@ $(FFTOBJS)

clean : 
	-$(RM) -f *.o *~ *.d work.pc* libfft.a

dist:
	tar czf Fft.tgz $(FFTOBJS:o=f) Makefile

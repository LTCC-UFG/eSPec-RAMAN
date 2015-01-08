OBJS = spline.o seval.o ttidbs.o
FC     = ifort
FFLAGS = -O3 -xW
sfx=${OSTYPE}-${MACHTYPE}


all: libspline.a
	cp -p libspline.a ..

clean:
	$(RM) *.o *~ *.d work.pc* libspline.a

libspline.a: $(OBJS)
	-$(RM) $@
	$(AR) cvq $@ $(OBJS)


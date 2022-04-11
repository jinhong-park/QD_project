include ../makefile_settings

sources=arpack_zn.h arpack_zn.cpp arpack_debug.h arpack_debug.cpp

objs=arpack_debug.o arpack_zn.o

%.o : %.cpp %.h
	$(CPPCOMPILER) $(CPPCOMPILEFLAGS) -I$(VLASTNE) $< -o $@


all : arpack.a


arpack.a: $(objs)
	$(AR) $(ARF) arpack.a $(objs)

clean:
	rm -f $(objs) arpack.a

backup:
	tar -cvzf arpack.tar.gz $(sources) makefile
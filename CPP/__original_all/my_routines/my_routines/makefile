include ../makefile_settings

sources=auxiliaries.cpp auxiliaries.h

objs=auxiliaries.o

%.o : %.cpp $(headers)
	$(CPPCOMPILER) $(CPPCOMPILEFLAGS)  $< -o $@

auxiliaries.a: $(objs)
	$(AR) $(ARF) auxiliaries.a $(objs)

all: auxiliaries.a

backup: 
	tar -cvzf my_routines.tar.gz $(sources) makefile

clean:
	rm -f auxiliaries.a $(objs)


include ../makefile_settings

region2d_version=region2d
#region2d_version=region2d_new


#unused main files: main_two_electronMartin.cpp 
sources2=main_graphene.cpp main_QPC.cpp main_impurities.cpp main_laser.cpp main_other.cpp main_majorana.cpp main_transport.cpp main_1e.cpp main_1e-Leon.cpp

#unused sources: 
sources=main.h energy.cpp test.cpp dens_mat.cpp fit.cpp $(region2d_version).cpp classes.cpp hamilt.cpp $(region2d_version).h classes.h ph_funcs.cpp two_electron.cpp two_electron2.cpp two_electron2.h laser.cpp laser.h transport.cpp transport.h

#for 1d skyrmion: simplified bundle
sources_skr=

%.o : %.cpp main.h classes.h
	$(CPPCOMPILER) $(CPPCOMPILEFLAGS) -I$(VLASTNE) -I$(ARPACK) $< -o $@

#unused targets: two_electron
targets=auxs.a graphene QPC two_electron impurities laser one_electron majorana transport skyrmion skyrmionP skyrmionC

#unused objs: 
objs=energy.o test.o hamilt.o fit.o $(region2d_version).o ph_funcs.o dens_mat.o classes.o two_electron2.o two_electron.o laser.o transport.o

objs_skr=classes.o

all: $(targets)

auxs.a: $(objs)
	ar r auxs.a $(objs)

two_electron: $(sources) auxs.a main_two_electronMartin.o main_two_electronMartin.cpp
	$(CPPCOMPILER)  main_two_electronMartin.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o two_electron 

one_electron-other: $(sources) auxs.a main_other.o main_other.cpp
	$(CPPCOMPILER)  main_other.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o one_electron-other

one_electron: $(sources) auxs.a main_1e-Leon.o main_1e-Leon.cpp
	$(CPPCOMPILER)  main_1e-Leon.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o one_electron

QPC: $(sources) auxs.a main_QPC.o main_QPC.cpp
	$(CPPCOMPILER)  main_QPC.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o QPC 

impurities: $(sources) auxs.a main_impurities.o main_impurities.cpp
	$(CPPCOMPILER)  main_impurities.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o impurities

graphene: $(sources) auxs.a main_graphene.o main_graphene.cpp
	$(CPPCOMPILER)  main_graphene.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o graphene

laser: $(sources) auxs.a main_laser.o main_laser.cpp
	$(CPPCOMPILER)  main_laser.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o laser

majorana: $(sources) auxs.a main_majorana.o main_majorana.cpp
	$(CPPCOMPILER)  main_majorana.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o majorana 

transport: $(sources) auxs.a main_transport.o main_transport.cpp
	$(CPPCOMPILER)  main_transport.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o transport

skyrmion: classes.h classes.cpp region1d.cpp region1d.o
	$(CPPCOMPILER)  region1d.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o skyrmion

skyrmionP: classes.h classes.cpp region1d-profile.cpp region1d-profile.o
	$(CPPCOMPILER)  region1d-profile.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o skyrmionP

skyrmionC: classes.h classes.cpp region1d-conv.cpp region1d-conv.o
	$(CPPCOMPILER)  region1d-conv.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o skyrmionC


clean:
	rm -f $(targets) $(objs) main_two_electron.o main_laser.o main_graphene.o main_impurities.o main_QPC.o main_other.o main_transport.o region1d.o

backup:
	tar -cvf working.tar.gz --gzip --ignore-failed-read $(sources) $(sources2) makefile data/gnuplot*wf


include ../makefile_settings

region2d_version=region2d
#region2d_version=region2d_new


#unused main files: main_two_electronMartin.cpp 
sources2=main_graphene.cpp main_QPC.cpp main_impurities.cpp main_laser.cpp main_other.cpp main_majorana.cpp main_transport.cpp main_1e.cpp main_1e-Leon.cpp

#unused sources: 
sources=main.h energy.cpp test.cpp dens_mat.cpp fit.cpp $(region2d_version).cpp classes.cpp hamilt.cpp $(region2d_version).h classes.h ph_funcs.cpp two_electron.cpp two_electron2.cpp two_electron2.h laser.cpp laser.h transport.cpp transport.h

#for Jin: simplified bundle
sources-Jin=main.h energy.cpp test.cpp dens_mat.cpp fit.cpp region2d.cpp classes.cpp hamilt.cpp region2d.h classes.h ph_funcs.cpp two_electron2.h laser.h transport.h two_electron.cpp two_electron2.cpp

%.o : %.cpp main.h classes.h
	$(CPPCOMPILER) $(CPPCOMPILEFLAGS) -I$(VLASTNE) -I$(ARPACK) $< -o $@

#unused targets: two_electron
targets=auxs.a graphene two_electron impurities laser one_electron majorana transport skyrmion skyrmion-check one_electron-Jin

#unused objs: 
objs=energy.o test.o hamilt.o fit.o $(region2d_version).o ph_funcs.o dens_mat.o classes.o two_electron2.o two_electron.o laser.o transport.o

objs-Jin=energy.o test.o hamilt.o fit.o region2d.o ph_funcs.o dens_mat.o classes.o two_electron2.o two_electron.o

all: one_electron-Jin

auxs.a: $(objs)
	ar r auxs.a $(objs)

auxs-Jin.a: $(objs-Jin)
	ar r auxs-Jin.a $(objs-Jin)

region1d.o : region1d.cpp region1d.h main.h classes.h
	$(CPPCOMPILER) $(CPPCOMPILEFLAGS) -I$(VLASTNE) -I$(ARPACK) region1d.cpp -o region1d.o

other: $(sources-Jin) auxs-Jin.a main_other.o main_other.cpp
	$(CPPCOMPILER)  main_other.o auxs-Jin.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o other

two_electron: $(sources) auxs.a main_two_electron.o main_two_electron.cpp
	$(CPPCOMPILER)  main_two_electron.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o two_electron 

one_electron-other: $(sources) auxs.a main_other.o main_other.cpp
	$(CPPCOMPILER)  main_other.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o one_electron-other

one_electron: $(sources) auxs.a main_1e-Leon.o main_1e-Leon.cpp
	$(CPPCOMPILER)  main_1e-Leon.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o one_electron

one_electron-Jin: $(sources-Jin) auxs-Jin.a main_1e-Jin.o main_1e-Jin.cpp
	$(CPPCOMPILER)  main_1e-Jin.o auxs-Jin.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o one_electron

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

skyrmion: classes.h classes.cpp region1d.cpp region1d.h region1d.o main_skyrmion.cpp main_skyrmion.o
	$(CPPCOMPILER)  main_skyrmion.o region1d.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o skyrmion

skyrmion-check: classes.h classes.cpp region1d.cpp region1d.h main_skyrmion-check.cpp main_skyrmion-check.o
	$(CPPCOMPILER)  main_skyrmion-check.o region1d.o auxs.a $(ARPACKL) $(VLASTNEL) $(CPPLINKFLAGS) -o skyrmionC



clean:
	rm -f $(targets) $(objs) main_two_electron.o main_laser.o main_graphene.o main_impurities.o main_QPC.o main_other.o main_transport.o region1d.o

backup:
	tar -cvf working.tar.gz --gzip --ignore-failed-read $(sources) $(sources2) makefile data/gnuplot*wf

backup-Jin:
	tar -cvf working-Jin.tar.gz --gzip --ignore-failed-read $(sources-Jin) makefile






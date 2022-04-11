#include "main.h"

//!parameters - all energies in meV
double T=0.1;                                                        //temperature in Kelvins
double Bstar=8, Estar;                                               //energy resolution in Tesla
double muF=8;                                                        //Fermi energy
double muS=muF/1e3;                                                  //0.1% spin accumulation
double alpha=0.05;                                                   //lever arm
double Vg=0;                                                         //gate voltage in milli volts
double B=0, Ez;                                                      //external field and Zeeman energy
double Gextra=0.0;                                                   //additional conductance

int Fder,which;
double cf;

void reset(double Bstar_=0, double B_=0) {
  T=0.1;
  Bstar=Bstar_;
  if (Bstar==0) Bstar=8; 
  Estar=Bstar*2*(0.39/2)*units.hbar*units.e/(2*units.m_e)/units.meV;
  B=B_;
  Ez=B*(0.39/2)*units.hbar*units.e/(2*units.m_e)/units.meV;
  Vg=0;
}

double F(double E)
{
  double ratio=E/cf;
  if (ratio<-100) ratio=-100;
  if (ratio>100) ratio=100;
  double x=exp(ratio);
  switch (Fder) {
    case 0 : return(1/(x+1));
    case 1 : return(-(x/cf)/((x+1)*(x+1)));
    case 2 : return(-1/(cf*cf)*x/(1+x)/(1+x)*(1-2*x/(1+x)));
  }
  return(0);
}

complex<double> underintegral(double E)
{
  double res;
  cf=T*units.kB/units.meV;
  Fder=1;
  res=-F(E);

  cf=-Estar/(2*M_PI);
  E-=alpha*Vg;
  switch (which) {
    case (1) : {
      Fder=0;
      return(res*(F(E-Ez)+F(E+Ez))/2);
    }
    case (2) : {
      Fder=0;
      return(res*(F(E+Ez)-F(E-Ez))/2);
    }
    case (3) : {
      Fder=1;
      return(res*(-1)*(F(E+Ez)+F(E-Ez))/2);
    }
    case (4) : {
      Fder=1;
      return(res*(-1)*(F(E+Ez)-F(E-Ez))/2);
    }
    case (5) : {
      Fder=2;
      return(res*(F(E+Ez)+F(E-Ez))/2);
    }
  }
  return(0);
}

void rangefix(double& E0)
{
  double g0=abs(underintegral(0));
  double E=E0;
  int step=0;
  do {
    double g=max(abs(underintegral(E)), abs(underintegral(-E)));
    if (g/g0 < 1e-10) E/=1.1;
    else if (g/g0 > 1e-8) E*=1.1;
    else break;
    step++;
  } while (step<1000);
  if (step<1000) E0=E;
}

double G(int which_)
{
  which=which_;
  double range=1*units.kB/units.meV*20;
  //rangefix(range);
  
  double err=1e-8;
  int maxsteps=14;
  complex<double> res=qromb(underintegral, -range, range, maxsteps, 4, err);
  return(res.real());
}

void calculate(double &Isig, double &Vsig, double& G1, double &G2, double& G3, double &G4)
{
  fprintf(logg, "Vg=%e [mV], T=%e [K], Estar=%e [meV] (Bstar=%e), B=%e [T], Ez=%e [meV], muS=%e [meV]:\n", Vg, T, Estar, Bstar, B, Ez, muS);
  G1=G(1);G2=G(2);G3=G(3);G4=G(4);
  
  //current for zero voltage
  Isig=-2*G2*muS - G3*muS*muS;
  //to units
  Isig*=units.e/(2*M_PI*units.hbar) * units.meV;
  //to picoamps
  Isig*=1e12;
  
  //voltage for zero current
  Vsig=(2*G2*muS+G3*muS*muS)/(2*(G1+Gextra)+muS*G4);
  //to units
  Vsig*=units.meV/units.e;
  //to nanoovolts
  Vsig*=1e9;
  
  fprintf(logg, "\tG1=%e [2e/h], G2=%e [2e/h], G3=%e [e/h meV], G4=%e [e/h meV], Isig=%e [pA], Vsig=%e[nV]\n", G1, G2, G3, G4, Isig, Vsig);
}

int main()
{

  const char* name[]={"QPC-gate", "QPC-temp", "QPC-resol", "QPC-B"};
  output_class output(4,"data/",name,"QPC-log");
  
  //output.set_appendix("");
  output.clear_files();                     //clear the output files
  output.open_files();                     //clear the output files

  //-------------------------------------------------
  //!I-sig and V-sig as a function of the gate
  for (prec p1=-10;p1<10.1;p1+=0.1) {

    reset();
    Vg=p1;

    prec Isig, Vsig, G1, G2, G3, G4;
    calculate(Isig, Vsig, G1, G2, G3, G4);
    prec G5=G(5)*(0.39/2)*units.hbar*units.e/(2*units.m_e)/units.meV*units.e/(2*M_PI*units.hbar);
    G5*=units.e*1e12/1e6;//pA / \mu V B
    fprintf(output.file(0), "%e %e %e %e %e %e\n", -Vg, G1, G3, Isig, Vsig, G5);
  }
  
  //-------------------------------------------------
  //!I-sig and V-sig as a function of the temperature
  for (prec p1=0.01;p1<11;p1*=1.1) {

    reset();
    T=p1;

    prec Isig, Vsig, G1, G2, G3, G4;
    calculate(Isig, Vsig, G1, G2, G3, G4);
    fprintf(output.file(1), "%e %e %e %e %e\n", T, G1, G3, Isig, Vsig);
  }

  //-------------------------------------------------
  //!I-sig and V-sig as a function of the Bstar
  for (prec p1=0.1;p1<110.01;p1*=1.1) {
    
    reset(p1,0);

    prec Isig, Vsig, G1, G2, G3, G4;
    calculate(Isig, Vsig, G1, G2, G3, G4);
    fprintf(output.file(2), "%e %e %e %e %e\n", Bstar, G1, G3, Isig, Vsig);
  }
  
  //-------------------------------------------------
  //!I-sig and V-sig as a function of the B
  for (prec p1=-0.5;p1<0.51;p1+=0.01) {
    
    reset(8,p1);
    
    prec Isig, Vsig, G1, G2, G3, G4;
    calculate(Isig, Vsig, G1, G2, G3, G4);
    fprintf(output.file(3), "%e %e %e %e %e %e %e\n", B, G1, G3, Isig, Vsig, G2, G4);
  }
            
  return(0);
}


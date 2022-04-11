#include <iostream>
#include <fstream>

#include "main.h"


void check_ipr()
{
  FILE* file_in=fopen("data/laser3/lasing_modes-int-conv","r");
  FILE* file_out=fopen("data/laser3/lasing_modes-int-conv-ipr","w");
  
  double ipr[20], D, psi2sum[20], psi4sum[20], x, dx=0;

  do {
    for (int state=0;state<20;state++) psi2sum[state]=psi4sum[state]=0;

    for (int line=0;line<999;line++) {
      if (fscanf(file_in,"%le %le", &D, &x)==EOF) {fclose(file_out);exit(0);}
      if (line==0) dx=x;
      if (line==1) dx=x-dx;
      for (int state=0;state<20;state++) {
	double re, im, psi2;
	int res=fscanf(file_in,"%le %le %le", &re, &im, &psi2);
	psi2sum[state]+=psi2;
	psi4sum[state]+=psi2*psi2;
      }
    }
    fprintf(file_out,"%e ",D);
    for (int state=0;state<20;state++) if (psi2sum[state]!=0) fprintf(file_out,"%e ",psi2sum[state]*psi2sum[state]/psi4sum[state]*dx); else fprintf(file_out,"0 ");
    fprintf(file_out,"\n");
    printf("ipr got at pump %e\n", D);
  } while(true);
}


void average(int length)
{
  FILE* file_in=fopen("data/laser3/lasing_modes-int-conv","r");
  FILE* file_out=fopen("data/laser3/lasing_modes-int-conv-ave2","w");
  
  double psi[1000*20], x[1000], D;

  do {
    for (int line=0;line<999;line++) {
      if (fscanf(file_in,"%le %le", &D, x+line)==EOF) {fclose(file_out);exit(0);}
      if (line>length) fprintf(file_out,"%e %e ", D, x[line-length/2]);

      for (int state=0;state<20;state++) {
	int pntr=state*1000+line;
	double re, im;
	int res=fscanf(file_in,"%le %le %le", &re, &im, psi+pntr);
	if (line>length) {
	  double sum=0;
	  for (int i=0;i<length;i++) sum+=psi[pntr-i];
	  fprintf(file_out,"%e ", sum/length);
	}
      }
      if (line>length) fprintf(file_out,"\n");
    }
    fprintf(file_out,"\n");
    printf("averaged at pump %e from length %i\n", D, length);
  } while(true);
}

void normalize()
{
  FILE* file_in=fopen("data/laser3/lasing_modes-int","r");
  FILE* file_out=fopen("data/laser3/lasing_modes-int-conv","w");
  
  double im[1000*20], re[1000*20], psi[1000*20], x[1000], D, sum[20];

  do {
    for (int state=0;state<20;state++) sum[state]=0;

    for (int line=0;line<999;line++) {
      if (fscanf(file_in,"%le %le", &D, x+line)==EOF) {fclose(file_out);exit(0);}
      for (int state=0;state<20;state++) {
	int pntr=line*20+state;
	int res=fscanf(file_in,"%le %le %le", re+pntr, im+pntr, psi+pntr);
	sum[state]+=psi[pntr];
      }
    }

    for (int line=0;line<999;line++) {
      fprintf(file_out,"%e %e ", D, x[line]);
      for (int state=0;state<20;state++) {
	int pntr=line*20+state;
	if (sum[state]>0) fprintf(file_out,"%e %e %e ", re[pntr]/sum[state], im[pntr]/sum[state], psi[pntr]/sum[state]);
	else fprintf(file_out,"%e %e %e ", 0.0, 0.0, 0.0);
      }
      fprintf(file_out,"\n");
    }
    fprintf(file_out,"\n");
    printf("converted at pump %e: sum of the lowest: %e\n",D,sum[0]);
  } while(true);
}

void convert_for_xmgrace(int which)
{
  FILE* file_in=fopen("data/laser3/lasing_modes-non-conv","r");
  double psi[1000*1000], x[1000], D[1000];
  int pump=0;
  bool further=true;

  do {
    for (int line=0;line<999;line++) {
      if (fscanf(file_in,"%le %le", D+pump, x+line)==EOF) {further=false;break;}
      for (int state=0;state<20;state++) {
	double aux1, aux2, aux3;
	int res=fscanf(file_in,"%le %le %le", &aux1, &aux2, &aux3);
	if (state==which) psi[pump*1000+line]=aux3;
      }
    }
    pump++;
  } while (further);
  fclose(file_in);
  
  char name[1000], app[10]; 
  itoa(which,app);
  strcpy(name, "data/laser3/lasing_modes-non-conv-xm");
  strcat(name, app);
  FILE* file_out=fopen(name,"w");

  for (int line=0;line<999;line++) {
    fprintf(file_out,"%e ", x[line]);
    for (int i=0;i<pump;i++) fprintf(file_out,"%e ", psi[i*1000+line]);
    fprintf(file_out,"\n");
  }
  fclose(file_out);
  file_out=fopen("data/laser3/lasing_modes-non-conv-xmD","w");
  for (int i=0;i<pump;i++) fprintf(file_out,"%i %e\n", i, D[i]);
  fprintf(file_out,"\n");
  fclose(file_out);
}



void give_fermi_contact_coupling()
{
  double p[]={0.3, 0.2, 0.5,1}, g[]={1.35, 1.70, 0.96,0.5553*2}, eta[]={2700,2700,4500,186}, ave=0;
  for (int i=0;i<3;i++) ave+=p[i]*g[i]*eta[i];
  double betaG=2*units.mi/(3.0)*units.muB*units.muN*ave;
  double AG=betaG/(pow(units.aGaAs,3)/8);
  double BG=AG*(3/2.0)/(0.44/2)/units.muB;
  double betaS=2*units.mi/(3.0)*units.muB*units.muN*p[3]*g[3]*eta[3];
  double AS=betaS/(pow(units.aSi,3)/8);
  double BS=AS*(1/2.0)/(2.0/2)/units.muB;
  AG/=(1e-3*units.meV);
  AS/=(1e-3*units.meV);
  printf("A(GaAs) = %e [mueV], A(Si)= %e [mueV] times 29Si content, naturaly 4,7perc gives %e [mueV]\n",AG,AS,AS*0.047 );
  betaG/=(1e-3*units.meV)*pow(units.nm,3);
  betaS/=(1e-3*units.meV)*pow(units.nm,3);
  printf("beta(GaAs) = %e [mueV nm^3], beta(Si)= %e [mueV nm^3] times 29Si content, naturaly 4,7perc gives %e [mueV nm^3]\n",betaG,betaS,betaS*0.047 );
  printf("fully polarized nuclei correspond to field B(GaAs)=%e and B(Si natural)=%e\n",BG,BS*0.047);
}

void give_res_gap_renorm()
{
  double lso=1e-6/1;
  double m=units.m_e*(0.067+0.024*0);
  double vf=units.hbar/(2*m*lso);
  double Ef=m*vf*vf/2/units.meV;
  double vs=vf, vc=vf, Kc=0.8, Ks=1;
  double Delta=1;
  double K2=(vs/Ks+vc*Kc)/(vc/Kc+vs*Ks);
  double K=sqrt(K2);
  double Ed=M_PI*units.hbar*vf/units.aGaAs/units.meV;
  double ratio=Ed/Delta;
  double expon=(1-K)/(2-K);
  double res=pow(ratio, expon);
  printf("factor K=%e, exponent=%e, E-num=%e, ratio=%e, the enhancement factor is %e, Ef=%e\n",sqrt(K2), expon, Ed, ratio, res, Ef);
}

void give_res_min_acc()
{
  double T=1.0;
  double mus=sqrt( 2*max(180*units.meV/1000, 2*M_PI*units.kB*T)*sqrt(units.hbar*1e+4*units.kB*T));
  printf("detection limit on spin accumulation is %e [meV]\n",mus/units.meV);
}
  
void give_res_el_nuc_spins()
{
  //electron
  double B=1, l=30e-9, ct=2480, Xi=1.4e+9*units.e, sigma=10*units.e, m=0.067, beta=4e-6*units.e*1e-27, rho=5300, T=1, w=8e-9, lso=1e-6, I=1.5, J=0.5;
  double EzN=1.2*units.muN*B;
  double Ql=EzN/units.hbar/ct*l;
  double hw=units.hbar*units.hbar/(2*units.m_e*m*l*l);
  double alphapz=Ql*Xi*l*2/hw;
  double alphadf=sigma*Ql*Ql/hw;
  printf("electron parameters: Ql=%e, hw=%e [meV], alphapz=%e, alphadf=%e\n", Ql, hw/units.meV, alphapz, alphadf);
  
  double Epz=sqrt( pow(units.hbar,7)*pow(ct,5)*rho )/ (Xi*beta*units.m_e*m);
  double gamma1=units.kB*T*EzN*EzN/(Epz*Epz)*(l*l*l*l)/(w*w*lso*lso)*9*I*I/(2*pow(M_PI,3))/units.hbar;
  printf("Epz=%e [meV], Gamma(1)=%e [1/s]\n", Epz/units.meV, gamma1);

  
  l=4.19e-9, ct=2358, Xi=3.4e+8*units.e, sigma=5*units.e, m=1.0/4.5, w=3.24e-9, I=5.0/2, J=3.0/2, rho=5650, beta=1.0/3*units.e*pow(0.61e-9,3)/4, T=10;
  double bV=beta/(4*M_PI*l*l*w/3), lambda=0.05, lambdap=0.15, Delta=100*units.meV;
  Ql=bV*J/units.hbar/ct*l;
  hw=units.hbar*units.hbar/(2*units.m_e*m*l*l);
  alphapz=Ql*Xi*l*2/hw;
  alphadf=sigma*Ql*Ql/hw;
  printf("hole parameters: beta/V=%e [meV], Ql=%e, hw=%e [meV], alphapz=%e, alphadf=%e\n", bV/units.meV, Ql, hw/units.meV, alphapz, alphadf);
  
  Epz=sqrt( pow(units.hbar,7)*pow(ct,5)*rho )/ (Xi*beta*units.m_e*m);
  gamma1=pow(3.0,5)/pow(2.0,5)/pow(M_PI,4)*I*I*J*J*units.kB*T/(Epz*Epz)*(lambda*lambda*beta*beta)/(l*l*w*w*w*w)/units.hbar;
  printf("Epz=%e [meV], Gamma(1)=%e [1/s]\n", Epz/units.meV, gamma1);
  
  double gammadp=pow(J,2)*pow(I,4)*pow(bV,6)*pow(sigma,2)/( 2*M_PI*pow(units.hbar,4)*pow(ct,5)*rho*pow(Delta,4))*pow(lambdap,2)*units.kB*T;
  printf("Gamma(DP)=%e [1/s]\n", gammadp);
}

void give_res_helical1()
{
  double ct=2480, rho=5300, QDR=2e-29, piezo=1.4e9, I=1.5, J=4e-6*units.e*pow(units.nm,3);

  //basic transition energy characteristics
  double energy=30e-3*units.kB/units.e;
  double omega=energy*units.e/units.hbar;
  double Q=omega/ct;
  printf("transition energy (corresponding to 30 mK) %e [eV]\ncorresponding angular frequency %e [1/s]\ntransversal phonon wavevector %e [1/nm] and wavelength %e [nm]\nequivalent magnetic field %e [T]\n", energy, omega, Q*units.nm, 2*M_PI/Q/units.nm, energy*units.e/units.muN);
  
  //dipole-dipole interaction strength for nearest neighbors in GaAs (a/sqrt(2) apart)
  double res, aux1, aux2, aux3;
  res=units.mi/(4*M_PI)*pow(units.muN,2)/pow(units.aGaAs/sqrt(2),3)*I*I/units.e;
  aux1=res*units.e/units.hbar;
  aux2=res*units.e/(units.muN*I);
  printf("dipole-dipole interaction strength is %e eV corresponsing to frequency %e and magnetic field %e\n",res,aux1, aux2);

  res*=units.e/units.hbar/M_PI/pow(50.0,2);
  printf("the resulting energy diffusion time over the wire radius distance is %e\n",1/res);
  
  //dipole-dipole phonon emission by recoil
  res=1/(2*M_PI*units.hbar*pow(ct,3)*rho)*pow(I,4);
  aux1=3/pow(units.aGaAs/sqrt(2),4)*units.mi/(4*M_PI)*pow(units.muN,2);
  res*=omega*aux1*aux1;
  printf("relaxation rate for dipole dipole and phonon emission is %e (aux=%e)\n",res,aux1);
  
  //quadrupole relaxation
  res=pow(units.e*QDR,2)*pow(omega,5)/(2*M_PI*units.hbar*pow(ct,7)*rho)*pow(piezo,2)*pow(I/4/(2*I-1),2);
  printf("relaxation rate for quadrupole and piezoelectric field is %e\n",res);
  
  res=units.m_e*0.067/(M_PI*pow(units.hbar,3))*pow(J/pow(50*units.nm,2),2)*energy*units.e/(8*units.meV);
  printf("relaxation rate electron flip-flop on shell excitation is %e\n",res);
  
}

double C_func(double g)
{
  return(sin(M_PI*g)/2*tgamma(1-g)*tgamma(1-g)*pow(2*M_PI,2*g-4));
}
  

void give_res_helical2()
{
  
  double Kc=0.5, g=0.75, gp=0.67, gpz=0.33;
  double vf=2e+5, A0=90e-3*units.meV, m=units.m_e*0.067;
  double deltaa=units.hbar*vf/units.aGaAs;
  double kf=m*vf/units.hbar, Ef=m*vf*vf/2;
  
  double T=30e-3, L=10e-6, I=3.0/2;
  double Nperp=2500, Npar=L/units.aGaAs;
  
  double J2kf=A0*A0/(deltaa)*C_func(g)*pow(deltaa/(units.kB*T),2-2*g)*pow(tgamma(g/2)/tgamma(1-g/2),2);
  //double J2kf_un=A0*A0/(deltaa)*C_func(1)*pow(tgamma(1.0/2)/tgamma(0/2),2);
  //double cost=J2kf/Nperp*I;
  //printf("J2kf=%e (without renormalization %e) [meV].\n", J2kf/units.meV, J2kf_un/units.meV);
  double lambdaT=units.hbar*vf/(units.kB*T);
  printf("lambdaT= %e mu m, temperature = %e [meV]\n",lambdaT*1e6, units.kB*T/units.meV);
  double E_RKKY=I*I*J2kf;

  printf("RKKY energy per particle: %e [meV]. Constituents: constant C=%e, J(2kf)=%e meV\n", E_RKKY/units.meV, C_func(g), J2kf/units.meV);

  double rho=8/pow(units.aGaAs,3);
  double Jdd=units.mi*units.muN*units.muN*rho;
  double E_dd=Jdd*I*I*0.25*Nperp;
  printf("dip-dip energy per particle: %e [meV]. Constituents: constant Jdd=%e meV\n", E_dd/units.meV,  Jdd/units.meV);

  double B=0.001;
  double EZ=pow(units.muN*B,2)/(2*J2kf/Nperp/Nperp);
  double Bcrit=2*I*J2kf/Nperp/units.muN;
  printf("external field energy per particle: %e [meV], critical field %e. Constituents: Zeeman energy %e meV supressed by Zeeman/J2kf ratio %e\n", EZ/units.meV, Bcrit, units.muN*B/units.meV, units.muN*B/J2kf*Nperp*Nperp);
  
  //printf("constant C=%e, J(2kf)=%e meV, energy cost of nuclear flip %e [meV] (temp %e, field %e)\n", C, J2kf/units.meV, cost/units.meV, cost/units.kB, cost/(1.3*units.muN));
  
  //double Jint=2*Npar/Nperp*I*A0*A0/2*pow(kf*units.aGaAs,2*g-1)/(4*M_PI*2*Ef)*tgamma(g-0.5)/(tgamma(g)*sqrt(M_PI))*tgamma(2*g)*sin(M_PI/2*(2*g-1));
  //printf("energy cost of nuclear flip from summing J_ij %e [meV] (temp %e, field %e)\n", Jint/units.meV, Jint/units.kB, Jint/(1.3*units.muN));
 
}

#define BIRTH 1

class birthdays_list
{
  int N;			//number of people
  int D;			//number of days
  int * birthday;		//list of datums (length N)
  FILE* logg;
  int entries;			//number of candidates
  char **name;			//list of names
  float *est;			//list of estimates
    
  public:
  birthdays_list(int D, int N);
  ~birthdays_list();
  void generate();		//generate new set of N integeres from 1 to D and sort ascendingly
  void shift();			//shift all numbers (modulo D) if first and the last are occupied 
  bool evaluate(int length);	//true if at least length are on consequtive days
  void logg_list();		//list data into the logg
  void flush_stat(int hits, int pntr);//output the table differences to estimates
};

birthdays_list::birthdays_list(int D_, int N_)
{
  logg=fopen("data/log","w+");
  if (logg==0) {printf("opening the log file failed\n"); exit(1);}
  D=D_;
  N=N_;
  birthday=new int[N];
  FILE* file=fopen("data/table.txt","r");
  int res=fscanf(file,"%d\n",&entries);
  name=new char*[entries];
  est=new float[entries];
  for (int i=0;i<entries;i++) {
    name[i]=new char[100];
    int res=fscanf(file,"%s\t%f\n",name[i],est+i);
  }
  #if BIRTH > 0 
  fprintf(logg,"class birthday_list created with days %i, people %i and %i estimates\n",D,N,entries);
  for (int i=0;i<entries;i++) {
    fprintf(logg,"%s\t%f\n",name[i],est[i]);
  }
#endif
}

birthdays_list::~birthdays_list() 
{
  #if BIRTH > 1 
  fprintf(logg,"class birthday_list deleted\n");
#endif

  delete birthday; 
  fclose(logg);
};

//generate new set of N integeres from 1 to D
void birthdays_list::generate()
{
  for (int i=0;i<N;i++) {
    int b=(int) ceil( D*((double) rand())/RAND_MAX );
    if (b==0) {printf("rnd failed, b=%i\n",b);b++;}
    if (b==D+1) {printf("rnd failed, b=%i\n",b);b--;}
    birthday[i]=b;
  }
#if BIRTH > 1 
  fprintf(logg,"birthdays generated randomnly with output:\n");
  logg_list();
#endif
  picsrt_inplace(N,birthday,'a');
#if BIRTH > 1 
  fprintf(logg,"after sorting this reads:\n");
  logg_list();
#endif
}

//shift all numbers (modulo D) if first and the last are occupied
void birthdays_list::shift()
{
  int shifts=0;
  do {
    if (birthday[0]==1 && birthday[N-1]==D) {//shift
#if BIRTH > 1 
  fprintf(logg,"shift necessary, since b[0]=%i and b[%i]=%i\n",birthday[0],N-1,birthday[N-1]);
#endif
      int hmn=0;
      for (int i=0;i<N;i++) if (++birthday[i]==D+1) {birthday[i]=1;hmn++;}
      for (int i=N-1;i>hmn;i--) birthday[i]=birthday[i-hmn];
      for (int i=0;i<hmn;i++) birthday[i]=1;
      shifts++;
    }
    else break;
  } while (true); 
#if BIRTH > 1
  fprintf(logg,"after %i shiftings the birthdays are:\n",shifts);
  logg_list();
#endif
}

//true if at least length are on consequtive days
bool birthdays_list::evaluate(int length)
{
  int act_length=1;
  for (int i=0;i<N-1;i++) {
    int dif=birthday[i+1]-birthday[i];
#if BIRTH > 1 
  fprintf(logg,"evaluation at i=%i with dif=%i and act_length=%i and length=%i\n",i,dif,act_length,length);
#endif
    //if (dif==0) continue;
    //else if (dif==1) act_length++;
    if (dif==0) act_length++;
    else {
      act_length=1;
#if BIRTH > 1 
   fprintf(logg,"act_length reset to 1\n");
#endif
    }
    if (act_length==length) return(true);
  }
  return(false);
}

void birthdays_list::logg_list()
{
  for (int i=0;i<N;i++) fprintf(logg,"%i ", birthday[i]);
  fprintf(logg,"\n");
}

void birthdays_list::flush_stat(int hits, int pntr)
{

  int res=system("clear");
  double act_est=((double) hits)/ ((double) pntr);
  printf("\nname\t\tresult\t\t\tsuccess probability\n");
  for (int i=0;i<entries;i++) {
    double prob=1-erf(sqrt(2*pntr)*abs(act_est-est[i]));
    printf("%s\t\t%f%%\t\t%f%%\n", name[i], est[i]*100, prob*100);
  }
  double prob=((double) hits)/ ((double) pntr);
  double rel=1/sqrt(pntr);
  int digits=floor(-log10(rel));
  printf("%i out of %i: actual estimate %f%%, %i digits precision", hits, pntr, (float) (prob*100), digits); 
   
  #if BIRTH > 0
  //fprintf(logg,"from hits=%i and tries=%i the probability is %e and rel. precision %e (that is %i valid digits)\n",hits, pntr, prob, rel, digits);
#endif

}

void dices()
{
  
  timeb cas;
  birthdays_list list(365,28);
  int pntr=0;
  int hits=0;
  ftime(&cas);
  long int last=cas.time*1000+cas.millitm;
  do {
    list.generate();
    list.shift();
    hits+=list.evaluate(3);
    pntr++; 
    if (pntr<0) {pntr--;break;}
    ftime(&cas);
    long int now=cas.time*1000+cas.millitm;
    if (now-last>200) {list.flush_stat(hits,pntr);last=now;}
    if (pntr<10) usleep(1000*1000);
    else if (pntr<10) usleep(100*1000);
    else if (pntr<500) usleep(10*1000);
  }
  while (hits<1000000000);
}


#include <gsl/gsl_multimin.h>

double fitting_error(const gsl_vector *v, void *params);

int main_for_fitting_error()
{
  randomize(99);
  randomize(999);
  double sets_vals[]={-1.343, 1.512, 1.705, 0.6526, 0.8942, 0.3333, 0.5485, 0.8166, 0.3449,0.3412,0.6740,1.127,0.4774, -0.4514, 1.634, 0.8357, 0.4736, 0.25, 162.8, 162.1, 0.25, 0.8992, 0.5371, 0.25, 163.6, 162.8, 0.25};
  size_t sets=9, params=3;
  FILE* iter_file=fopen("data/minim_status","a+");
  
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  
  int status;
  double size;
  size_t iter=0, pntr=0;
  
  minex_func.n = params;
  minex_func.f = fitting_error;
  //minex_func.params = &iter;
  minex_func.params = &pntr;
  
  s = gsl_multimin_fminimizer_alloc (T, params);
  
  do {
        
    iter = 0;
    
    fprintf (iter_file,"\nrestarting at pntr=%i\n\n", (int) pntr);

    /* Starting point */
    x = gsl_vector_alloc (params);
    ss = gsl_vector_alloc (params);

    if (sets>pntr) {//predefined values as the starting point
      gsl_vector_set (x, 0, sets_vals[pntr*3+0]);
      gsl_vector_set (x, 1, sets_vals[pntr*3+1]);
      gsl_vector_set (x, 2, sets_vals[pntr*3+2]);
      
      /* Set initial step sizes */
      gsl_vector_set (ss, 0, 0.05);
      gsl_vector_set (ss, 1, 0.05);
      gsl_vector_set (ss, 2, 0.05);
    }
    else {
      gsl_vector_set (x, 0, generuj());
      gsl_vector_set (x, 1, generuj());
      gsl_vector_set (x, 2, generuj());
      
      /* Set initial step sizes */
      gsl_vector_set (ss, 0, (0.5-generuj())/3);
      gsl_vector_set (ss, 1, (0.5-generuj())/3);
      gsl_vector_set (ss, 2, (0.5-generuj())/3);
    }
    
    /* Initialize method */
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    
    do /* iterate */
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status) break;
      
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-3);
      
      if (status == GSL_SUCCESS) fprintf (iter_file,"converged to minimum at\n");
      
      fprintf (iter_file,"%5d %10.3e %10.3e %10.3e f() = %7.3f size = %.3f\n",(int) iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), gsl_vector_get (s->x, 2), s->fval, size);
      
    } while (status == GSL_CONTINUE && iter < 10 && false);
    
  } while (pntr++<10);
    
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  fclose(iter_file);
  return(status);
}



//int main()
double fitting_error(const gsl_vector *v, void *params)
{
  
  //!GaAs parameters
  parameters.set("eps0",12.9,"-");          //static dielectric constant
  parameters.set("density",5300,"kg/m3");   //density
  parameters.set("cl",5290,"m/s");          //longitudinal velocity
  parameters.set("ct",2480,"m/s");          //transversal velocity
  parameters.set("def",7,"eV");             //deformation potential
  parameters.set("piezo",1.4E9,"eV/m");     //piezoelectric constant
  parameters.set("m_eff",0.067,"m_e");      //effective mass (bulk: 0.067).
    
  parameters.set("gx",-0.364,"-");           //g factors
  parameters.set("gy",-0.364,"-");
  parameters.set("gz",-0.364,"-");
  
  parameters.set("br",1.24/2*3.79,"meVA");       //bychkov-rashba
  parameters.set("dress",1.26*7.164,"meVA");    //dresselhaus
  parameters.set("dress3",27.5,"eVA^3");  //dresselhaus cubic
  parameters.set("T",0.25,"K");        //temperature
  
  parameters.set("exch",0.0,"meV");
  parameters.set("xMn",0,"-");

  //!quantum dot setup
  parameters.set("lv",34,"nm");             //confinement length l_0
  
  parameters.set("B",6.5,"Tesla");         //magnetic field (earth's magn field ~5E-5T)
  parameters.set("theta_B",0.5,"pi");        //0->perpendicular field; 0.5->in-plane field
  parameters.set("phi_B",0.25,"pi");          //in-plane orientation with respect to crystallographic x axis

  parameters.set("ALHO",0,"-");
  parameters.set("d",76.5,"nm");  		        //half of the distance between quantum dot minima
  parameters.set("phi_d",0.25,"pi");         //in-plane orientation with respect to crystallographic x axis

  //!computational parameters
  int J = 4;                               //[./.]
  parameters.set("eig",J,"-");              //[J] number of single electron states (including spin) to obtain in diagonalization. Must be large enough because we will kick states when sorting.
  parameters.set("eigmax",2*J,"-");         //[4*eig] max. number of states in sorting for the previous - must be large enough to accomodate all crossings. Must be at least 3*eig.

  //!grids
  parameters.set("extent",5,"-");		        //[5] Don't change!
  parameters.set("grid_step",0.2,"-"); 	    //[0.2/0.15] determines grid dimension and single electron functions precision. (For comparison: dim = 2*extent/grid_step (approx) => dim = 50 @ grid_step=0.2).
  parameters.set("CE_exp",5,"-");          //[10/20] Sets precision of the relaxation integrals.
  parameters.set("sets",2,"-");             //=2S+1, where S is the spin of the single electron case (sets=2 => electrons with spin 1/2, which is needed for single-e relax)).
                          
  //!sorting
  parameters.set("threshold",8E-1,"-");		  //[1E-5] tolerated unprecision in sorting accoring to the symmetries.
  parameters.set("sorting_method",-4,"-"); 	//0->no sorting; 1->I; 2->Ix,Iy; 3->lz; 4->distance; 5->FD-states; *-1->adds spin. NOTE: 1e and b states MUST BE fine! f states need only to be fine up to the states needed.
                                            //example: use -1 for double dot (vs magnetic field) spectrum.


  const char* name[]={"measured-relax.dat","measured-1e.dat","measured-fitted.dat","measured-fitting-results.dat"};
  output_class output(4,"data/",name,"log");

  //output.clear_files();                     //clear the output files
  
  prec error_sum;
  
  
  double sets_vals[]={\
  -1.343, 1.512, 1.705,\
  0.6526, 0.8942, 0.3333,\
  0.5485, 0.8166, 0.3449,\
  0.3412,0.6740,1.127,\
  0.4774,-0.4514, 1.634, 
  0.8357,0.4736,0.25,\
  162.8, 162.1, 0.25,\
  0.8992, 0.5371, 0.25,\
  163.6, 162.8, 0.25};
  
  //-------------------------------------------------
  //!the outer loop
  for (prec p1=0;p1<10;p1++) 
  //for (prec p1=0;p1<10.01;p1=max(0.1,1.5*p1)) 
  {
    //parameters.set("phi_B",p1,"pi");
    prec alpha=sqrt(5.0/4/(1+p1*p1));
    prec beta=p1*alpha;
    //parameters.set("br",1.24*beta,"meVA");       //bychkov-rashba
    //parameters.set("dress",1.26*alpha,"meVA");    //dresselhaus

    parameters.set("br",gsl_vector_get(v,0),"meVA");       //bychkov-rashba
    parameters.set("dress",gsl_vector_get(v,1),"meVA");   //dresselhaus    
    //parameters.set("phi_B",gsl_vector_get(v,2),"pi");      //orientation
    //parameters.set("phi_d",gsl_vector_get(v,2),"pi");

    parameters.set("br",sets_vals[(int) p1*3+0],"meVA");       //bychkov-rashba
    parameters.set("dress",sets_vals[(int) p1*3+1],"meVA");   //dresselhaus    
    parameters.set("phi_B",sets_vals[(int) p1*3+2],"pi");      //orientation
    parameters.set("phi_d",sets_vals[(int) p1*3+2],"pi");

    
    FILE* file_exp=fopen("data/experiment.png.2.dat","r");
    float next_e,next_G;
    do {
      fscanf(file_exp,"%f %f\n",&next_e,&next_G);
      if (next_e>0) break;
    } while (true);
    bool stop_fitting=false;
    error_sum=0;
    int pntr=0;

    int app = ((int *) params)[0]; 
    
    //output.set_appendix("th",(int) round(p1*100));      //appendix to the output files
    //output.set_appendix("iter",app);      //appendix to the output files
    output.set_appendix("set",(int) p1);      //appendix to the output files
    output.clear_files();                          //clear the output files
    output.open_files();    

    state_info *states=0, *states_previous=0;

    //!initialize the sorting classes for single and double electron states
    sorting_class sorting1e;
    sorting1e.init((int) parameters.read("eigmax","-"), (int) parameters.read("eig","-"),3);

    //!the inner loop
    prec step_max=100, step=0;
    for (prec p2=0.0;p2<350 && !stop_fitting;p2+=step) {
      //parameters.set("d",p2,"nm");
      //parameters.set("B",p2,"Tesla");
      //if (p2<134.01 || p2>140) step=2; else step=0.2;
      
      step=min(step_max, next_e-p2+1e-5);
      
      prec Ebias=p2*1e-3*units.meV/(2*units.e*parameters.read("d","nm")*units.nm);
      fprintf(logg,"bias energy of %f mueV converted to field %.2e V/m\n",p2,Ebias);
      parameters.set("Ebias",Ebias,"V/m");

      sprintf(buffer,"parameters set: p1=%f, p2=%f\n",p1,p2);
      splachni(logg,buffer,1+4);

      //!the single electron diagonalization
      parameters.recompute();
      states=diagonalize();

      message(logg,"sorting...",1);
      sorting1e.inspect_new_1e_states(states);
      sorting1e.sort(p2, states->sorted2un, states->unsorted2s);
      message(logg,"...done\n",1);

      //some info about the single electron states into log and an output file
      fprintf(output.file(1),"%.14e ",p2);
      states->show_set("clear");
      states->show_set("energy");
      statistics_multi(*states, 0, 2+8, output.file(1),2);
  
      states->show_set("Ix");states->show_set("Iy");
      states->show_set("I");
      states->show_set("s");states->show_set("lz");
      //states->show_set("sx");states->show_set("sy");
      //states->show_set("sz");
      statistics_multi(*states, 0, 1+2+4+8+16, logg,4);
  
      //!!! relax rates
      //the looping parameter into output file
      sprintf(buffer,"%.14e ",p2);
      splachni(output.file(0),buffer,2);
      
      //for the lowest spin down state, first excited spin up to the ground state - phonon/interaction type
      content* temp=new content[states->r1->Nw*states->r1->sets];
      prec fitted_result=0;
      for (int from=0; from<4; from++ ) {//initial state
        for (int to=0; to<4; to++) {    //final state
          prec rate_sum=0;
          for (int which=0; which<3; which++) {//channel
            //#0: Def-LA (def=Xi_d,def2=Xi_u); #1,#2: Piezo(LA,TA); #3: Def-TA (def2=Xi_u)
            prec rate=0;
            if (from!=to) rate=transition_rate(*states,states->sorted2un[from],states->sorted2un[to],which);
            fprintf(output.file(0),"%e ",rate);
            rate_sum+=rate;
          }
          if (from==2 && to!=3) fitted_result+=rate_sum;
          //dipole matrix elements
          prec phi_d=parameters.read("phi_d","pi")*M_PI;
          content* ket=states->vysl->eigenvecs+states->r1->Nw*states->r1->sets*states->sorted2un[to];
          content* bra=states->vysl->eigenvecs+states->r1->Nw*states->r1->sets*states->sorted2un[from];
          states->r1->op_ket(ket, temp, states->r1->r_x, states->r1->set);
          content x=0;
          for (int i=0;i<states->r1->sets;i++) x+=states->r1->braket(bra,i,temp,i); 
          states->r1->op_ket(ket, temp, states->r1->r_y, states->r1->set);
          content y=0;
          for (int i=0;i<states->r1->sets;i++) y+=states->r1->braket(bra,i,temp,i);
          content rd=x*cos(phi_d)+y*sin(phi_d), rcd=x*sin(phi_d)-y*cos(phi_d);
          prec c=units.length/units.nm;
          fprintf(output.file(0),"%e %e %e ",rate_sum, abs(rd)*c, abs(rcd)*c);
        }
      }
      delete temp;
      fprintf(output.file(0),"%e \n", fitted_result);
      
      if (p2>next_e && !stop_fitting) {//if we fit the experimental data
        //calculate the error
        prec error=pow(log10(fitted_result)-log10(next_G),2);
        error_sum+=error;
        fprintf(output.file(2),"%i %e %e %e %e %e %e\n", pntr, p2, next_e, fitted_result, (double)  next_G, error, error_sum);
        pntr++;
        //read the next entry from the data table        
        if (fscanf(file_exp,"%f %f\n",&next_e,&next_G)!=2) {
          fclose(file_exp);
          fprintf(output.file(3),"%e %i %e \n", p1, pntr, error_sum);
          stop_fitting=true; 
        }
        step=min(step_max, next_e-p2+1e-5);
      }

      //step=sorting1e.adjust_step();
      if (states_previous!=0) clean_up(states_previous);
      states_previous=states;
      fflush(0);
    } //!second parameter (p2) loop ends here
    clean_up(states);
    sorting1e.deallocate();
  } //!first parameter (p1) loop ends here 

  return(error_sum);
}

#include <string>

char coding_aux (char in, int shift)
{
  const char* keyboard[3]={"qwertyuiop","asdfghjkl","zxcvbnm"};
  int length[]={10,9,7};
  bool found=false;
  int ii,jj;
  for (int i=0;i<3;i++) {
    for (int j=0;j<length[i];j++) {
      if (in==keyboard[i][j]) {found=true;ii=i;jj=j;break;}
    }
  }
  //printf("%c found:%i at ii=%i, jj=%i\n", in, found, ii, jj);
  if (!found) return(in);
  jj+=shift;
  while (jj<0) {jj+=length[ii];}
  while (jj>=length[ii]) {jj-=length[ii];}
  //printf("returning %c\n", keyboard[ii][jj]);
  return(keyboard[ii][jj]);
}

void coding (char* message, int shift)
{
  for (int i=0;i<strlen(message);i++) {
    printf("%c",coding_aux(message[i],shift));
    shift*=-1;
  }
  printf("\n");
}


int main()
{
  /*coding("ponorit sa, alebo radsej plavat?",1);
  coding("plavat ako lod, ktora sa klze po tvojej vonavej a vlhkej lesklej hladine.",2);
  coding("ale musime pritom uzavriet elektricky okruh, medzi styrmi rukami a styrmi ocami.",3);
  coding("len vtedy sa mozeme do seba spravne zacvaknut v podpalubi.",1);
  coding("len vtedy sa lod spravne knise a vzdychaju jej pritom vsetky plachty.",2);
  coding("napinaju a uvolnuju sa lana a vydavaju pritom zvuky, ktore sa podobaju na starodavnu piesen o modrych ociach.",3);
  coding("pribeh o tom ako sa lod vynara a zase ponara, stale tam a spat, von a dnu, do tvojho zvlneneho oceanu.",1);*/
  
  const char* name[]={"noname"};
  output_class output(0,"data/",name,"log");
  output.clear_files();
  output.open_files();
  
  int N=101;
  double k=N/2.0, theta=0, delta=1, t=1;
  
  lapack_dstev_prototype lapack(N);
  double *d=new double[N];
  double *e=new double[N];
  double *vecs=new double[N*N];
  double *vals=new double[N];
  
  progress_bar_prototype progress_bar("lapack dstev");
  progress_bar.start(true);
    for (int i=0;i<N;i++) {d[i]=delta*cos(M_PI/N*k*i+theta); e[i]=t;}
    lapack.diagonalize(d, e, vals, vecs, true);
  for (int i=0;i<N;i++) printf("eigval nr %i is %e \n",i,vals[i]);
  progress_bar.finished(true);
  
  delete d;
  delete e;
  delete vecs;
  delete vals;
  exit(0);
  
  //check_de_noise();
  //exit(0);
  //printf("hello\n");
  //dices();
  //give_fermi_contact_coupling();
  //give_res_gap_renorm();
  //exit(0);
  //check_ipr();
  //average(10); exit(0);
  //for (int i=0;i<20;i++) convert_for_xmgrace(i); exit(0);
  //normalize();
  
  //printf("units.length=%e\n",units.length);
  //!aux: units calculation
  //give_res_helical2();
  //exit(0);
  return(0);
}












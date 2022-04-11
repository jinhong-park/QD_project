#include "main.h"
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>

#define CASE 1
#define MAJO 3
#define MAJOT 3
//#define MAJOSCAN 2

double transmutation_gap(double& zlow, double& zhigh);

void constructCNT(content* C)
{
  content ii=content(0,1);
#if CASE == 1
  content aux[]={
    1.0,ii,1.0,ii,1.0,-ii,1.0,-ii,\
    ii,1.0,ii,1.0,-ii,1.0,-ii,1.0,\
    1.0,-ii,1.0,-ii,1.0,ii,1.0,ii,\
    -ii,1.0,-ii,1.0,ii,1.0,ii,1.0};
#endif
#if CASE == 2
  content aux[]={
    1.0,ii,1.0,ii,1.0,-ii,1.0,-ii,\
    ii,1.0,ii,1.0,-ii,1.0,-ii,1.0,\
    1.0,-ii,1.0,-ii,1.0,ii,1.0,ii,\
    -ii,1.0,-ii,1.0,ii,1.0,ii,1.0};
#endif
    for (int i=0;i<8*4;i++) C[i]=aux[i];
}

class Hamiltonian {

public:
  content *M;
  void constructH0();
  void addHk(content k);
  void mE();
  
  Hamiltonian() {M=new content[8*8];};
  ~Hamiltonian() {delete M;};
};

struct pars {

  double s, s2;
  double z, z2;
  double n, n2;
  double theta;
  double mu, mu2;
  double E, E2;
  
  class Hamiltonian H;
  
  content *u, *vt, *vecs8, *vecs4, *C, *k;
  double *sv;
  lapack_zgesvd_prototype* svd8, *svd4;
  
  void init();
  inline void refresh();
  pars() {};
  ~pars();
} p;

void pars::init() {

  svd8=new lapack_zgesvd_prototype(8,8);
  svd4=new lapack_zgesvd_prototype(4,4);

  k=new content[4];
  vt=new content[8*8];
  u=new content[8*8];
  C=new content[8*4];
  vecs8=new content[8*8];
  vecs4=new content[4*4];
  sv=new double[8];
  
  constructCNT(C); 
}

pars::~pars() {
   delete svd8;
   delete svd4;
   delete k;
   delete vt;
   delete u;
   delete C;
   delete vecs8;
   delete vecs4;
   delete sv;
}

inline void pars::refresh()
{
  p.s2=p.s*p.s;
  p.n2=p.n*p.n;
  p.mu2=p.mu*p.mu;
  p.z2=p.z*p.z;
  p.E2=p.E*p.E;
}

//return different solutions for wavevectors with non-negative imaginary part and their number
int veck(content* k)
{
  p.refresh();
  double s1=p.E2 - p.n2 - p.s2 + p.mu2;
  double d1=p.n2*p.s2+(p.E2-p.s2)*p.mu2;
  if (abs(d1)<1e-12) d1=0;
  k[0]=sqrt(s1-2.0*sqrt((content) d1));
  k[1]=sqrt(s1+2.0*sqrt((content) d1));
  if (abs(s1)<1e-12) k[1]=-k[0];    //tweak by hand
  if (abs(d1)<1e-12) k[1]=k[0];    //tweak by hand
  double s2=p.E2 - p.z2 - p.s2 + p.mu2;
  double d2=p.z2*p.s2+(p.E2-p.s2)*p.mu2;
  if (abs(d2)<1e-12) d2=0;
  k[2]=sqrt(s2-2.0*sqrt((content) d2));
  k[3]=sqrt(s2+2.0*sqrt((content) d2));
  if (abs(s2)<1e-12) k[3]=-k[2];    //tweak by hand
  if (abs(d2)<1e-12) k[3]=k[2];    //tweak by hand
#if MAJO>1
  fprintf(logg,"\t\tveck gives the following 4 k vectors (out of s1=%e, d1=%e, s2=%e, d2=%e)\n", s1,d1,s2,d2);
  for (int i=0;i<4;i++) fprintf(logg,"\t\t\tk[%i]=%e%+ei\n",i,k[i].real(),k[i].imag());
#endif
  
  for (int i=0;i<4;i++) if (k[i].imag()<0) k[i]*=-1;
  bool redundant[4]={0,0,0,0};
  for (int i=1;i<4;i++) {
    for (int j=0;j<i;j++) {
#if MAJO>3
  fprintf(logg,"\t\tredundancy check:abs(k[%i]-k[%i])=%e\n", i,j,abs(k[i]-k[j]));
#endif
      if (redundant[j]) continue;
      if (abs(k[i]-k[j])<1e-8 || abs(k[i].imag())<1e-12) redundant[i]=true;
    }
  }
  int l=0;
  for (int i=0;i<4;i++) if (!redundant[i]) k[l++]=k[i];
  
#if MAJO>1
  fprintf(logg,"\t\tafter redundeness check the following %i stayed\n", l);
  for (int i=0;i<l;i++) fprintf(logg,"\t\t\tk[%i]=%e%+ei\n",i,k[i].real(),k[i].imag());
#endif
  return(l);
}

//give the gap as the minimum of the E^2 as a function of the K^2
double gap()
{
  p.refresh();
  double E2min, Emin, Kmin;
  if (p.mu==0) Emin=abs(abs(p.n)-abs(p.s)); 
  else {
    double Kminnumer=p.mu2*p.mu2-(p.mu2+p.s2)*p.n2;
    if (Kminnumer>0) E2min=p.s2*(1.0-p.n2/p.mu2);
    else E2min=p.n2+p.s2+p.mu2-2*sqrt(p.n2*(p.s2+p.mu2));
    if (E2min<0) Emin=0; else Emin=sqrt(E2min);
  }
  double branch1=Emin;
  if (p.mu==0) Emin=abs(abs(p.z)-abs(p.s)); 
  else {
    double Kminnumer=p.mu2*p.mu2-(p.mu2+p.s2)*p.z2;
    if (Kminnumer>0) E2min=p.s2*(1.0-p.z2/p.mu2);
    else E2min=p.z2+p.s2+p.mu2-2*sqrt(p.z2*(p.s2+p.mu2));
    if (E2min<0) Emin=0; else Emin=sqrt(E2min);
  }
  double branch2=Emin;
#if MAJO>1
  fprintf(logg,"gap in branch 1 is %e, in branch 2 is %e\n",branch1, branch2); 
#endif
  return(min(branch1,branch2));
}

//construct the Hamiltonian at k=0
void Hamiltonian::constructH0()
{
  double nr=p.n*cos(p.theta);
  double ni=p.n*sin(p.theta);
  
  for (int i=0;i<8*8;i++) M[i]=0;
  
  //block B3
  M[0*8+4]=content(0,-p.s-p.z);
  M[0*8+5]=content(0,-p.mu);
  M[1*8+4]=content(0,-p.mu);
  M[1*8+5]=content(0,p.s-p.z);
  
  //block B4
  M[2*8+6]=content(0,-nr-p.s);
  M[2*8+7]=content(0,ni-p.mu);
  M[3*8+6]=content(0,-ni-p.mu);
  M[3*8+7]=content(0,p.s-nr);
  
  //block B1
  M[4*8+0]=content(0,p.s+p.z);
  M[4*8+1]=content(0,p.mu);
  M[5*8+0]=content(0,p.mu);
  M[5*8+1]=content(0,-p.s+p.z);
  
  //block B2
  M[6*8+2]=content(0,nr+p.s);
  M[6*8+3]=content(0,p.mu+ni);
  M[7*8+2]=content(0,p.mu-ni);
  M[7*8+3]=content(0,nr-p.s);
}

//add the k-dependent part into the Hamiltonian
void Hamiltonian::addHk(content k)
{
  M[0*8+4]+=k;
  M[1*8+5]-=k;
  M[2*8+6]-=k;
  M[3*8+7]+=k;
  M[4*8+0]+=k;
  M[5*8+1]-=k;
  M[6*8+2]-=k;
  M[7*8+3]+=k;
}

//subtract the energy times identity
void Hamiltonian::mE() {for (int i=0;i<8;i++) M[i*8+i]-=p.E;}


//phase at given parameters. Output: nullspace dimension + singular values
int phase(double* sv, int debug=0, FILE* file=0)
{
  if (file==0) file=logg;

#if MAJO > 2
  if (debug==0) debug=2;
#endif

  int l=veck(p.k), pntr=0;
  
  for (int i=0;i<l;i++) {
    p.H.constructH0();
    p.H.addHk(p.k[i]);
    p.H.mE();

  if (debug>3) {
    fprintf(file, "for k#%i the H-E was constructed as:\n",i);
    for (int i=0;i<8;i++) {
      for (int j=0;j<8;j++) fprintf(file,"%.4e%+.4e ",p.H.M[i*8+j].real(),p.H.M[i*8+j].imag());
      fprintf(file,"\n");
    }
  }
    
    //nullspace of the hamiltonian for the wavevector solution
    p.svd8->decompose(p.H.M, p.sv, p.vt, p.u);
    int hmn;
    for (hmn=0;hmn<8;hmn++) if (p.sv[7-hmn]>1e-12*p.sv[0]) break;//find number of zero singular values and copy nullvectors into vecs
  if (debug>1) {
    fprintf(file,"\t\tSVD decomposition of H-E gave %i zero values in the following set\n\t\t", hmn);
    for (int i=0;i<8;i++) fprintf(file,"%e, ", p.sv[i]);
    fprintf(file,"\n");
  }
    for (int j=0;j<hmn;j++) {
      for (int k=0;k<8;k++) p.vecs8[pntr*8+k]=conj(p.vt[(7-j)*8+k]);//a nullvector=(row of VT)^*
      pntr++;							   
    }
  }
#if MAJO >0
  if (pntr<4) fprintf(file,"too few (%i) nullvectors of HmE found at E=%e\n", pntr, p.E);
#endif
  if (pntr>4) {
    fprintf(file,"WARNING: too many (%i) nullvectors of HmE found at E=%e... ignoring them\n", pntr, p.E);
    FILE* file=fopen("data/majorana-too-many","a+");
    fprintf(file,"too many nullvectors encountered at s=%e, n=%e, z=%e, mu=%e, theta=%e, E=%e\n", p.s, p.n, p.z, p.mu, p.theta, p.E);
    fprintf(file,"\tfrom %i independent k-vectors:\n", l);
    for (int i=0;i<l;i++) {
      fprintf(file,"\t\t\tk[%i]=%e%+ei gave the following singular values:\n",i,p.k[i].real(),p.k[i].imag());
      p.H.constructH0();
      p.H.addHk(p.k[i]);
      p.H.mE();
      p.svd8->decompose(p.H.M, p.sv, p.vt, p.u);
      for (int i=0;i<8;i++) fprintf(file,"%e, ", p.sv[i]);
      fprintf(file,"\n");
    }
    fclose(file);
    pntr=4;
  }
 
  //contract the nullvectors into a 4x4 matrix
  for (int i=0;i<4;i++) {
    for (int j=0;j<pntr;j++) {
      p.vecs4[i*4+j]=0;
      for (int k=0;k<8;k++) p.vecs4[i*4+j]+=p.C[i*8+k]*p.vecs8[j*8+k];
    }
  }

  //svd of the contracted matrix
  if (pntr>0) p.svd4->decompose(p.vecs4,p.sv,0,0,false,4,pntr);///!!!4 and pntr were swapped???

  //swap the svd values
  for (int i=0;i<pntr;i++) sv[i]=p.sv[pntr-i-1];
  
  //identify the nullspace rank
  int hmn;
  for (hmn=0;hmn<pntr;hmn++) if (sv[hmn]>1e-6*sv[pntr-1]) break;
  if (debug>1) {
    fprintf(file,"\t\tSVD decomposition of nullspace (matrix dim %i) at E=%e gave %i zero values from the following set\n\t\t",pntr, p.E, hmn);
    for (int i=0;i<pntr;i++) fprintf(file,"%e, ", sv[i]);
    fprintf(file,"\n");
  }
  return(hmn);
}

//returns number of local minima on the sv(E) curve away from E=0
int prescan(double *mins, double& step, int what=0)
{
  double Egap=gap(), limit;
  if (what==0) limit=Egap; else limit=2*M_PI;
  step=(limit*(1-1e-6))/100;

#if MAJO>0
  fprintf(logg,"\t\tprescan (what=%i) started at parameters s=%.3e, z=%.3e, n=%.3e, mu=%.3e, theta=%.3e with gap %e, upper limit %e, step %e\n", what, p.s, p.z, p.n, p.mu, p.theta, Egap, limit, step);
#endif
  //!!!change this lower limit on the gap
  if (Egap<1e-6) return(0);
  //if (step<1e-4) step=1e-4;

  int pntr=0;
  double * sv=new double[4], dprev=0, dnow, now, prev=0, param;
  do {
    if (what==0) {
      param=p.E+=step;
      phase(sv);
      now=sv[0];
    } else {
      double zlow=0, zhigh=2;
      now=transmutation_gap(zlow, zhigh);
      param=p.theta+=step;
    }
    dnow=now-prev;
    if (dprev<-step/1e+8 && dnow>step/1e+8 && (what==0 || now<100)) mins[pntr++]=param-step; //local minimum found
    if (pntr==10) {
      fprintf(logg,"WARNING: too many minima, most probably some corruption, everything will be trashed\n");
      pntr=0;
      break;
    }
#if MAJO>1
    fprintf(logg,"\t\tprescan step: at param=%e, watched value=%e, previous=%e, dnow=%e, dprev=%e. Number of minima: %i\n", param, now, prev, dnow, dprev, pntr);
#endif
    prev=now;
    dprev=dnow;
  } while((abs(dnow)>1e-12 || what!=0) && param<limit-step);
  delete sv;
  return(pntr);
}

double sv_aux[4], coor[10];
complex<double> val[10];
double minimize_func (double x, void * params)
{
  int what=*(int *) params;
  if (what==0) {
    p.E=x;
    phase(sv_aux);
    return(sv_aux[0]);
  } 
  if (what==1) { 
    p.theta=x;
    double zlow=0, zhigh=2;
    return(transmutation_gap(zlow, zhigh));
  }
  if (what==2) {
    double err;
    complex<double> res=polint(coor, val,  3, x, err);
    return(res.real());
  }
  return(0);
}

double minimize(double a, double b, double &m, int error, int what=0)
{
  gsl_set_error_handler_off();
#if MAJO>0
  fprintf(logg,"\t\tgsl-minimize started at parameters s=%.3e, z=%.3e, n=%.3e, mu=%.3e, theta=%.3e with interval a=%e, b=%e, minimum estimated at E=%e\n", p.s, p.z, p.n, p.mu, p.theta, a, b, m);
#endif
  int status, iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  gsl_function F;
  
  F.function = &minimize_func;
  F.params = (void *) &what;

  T = gsl_min_fminimizer_goldensection;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, a, b);

#if MAJO>1
  fprintf (logg,"using %s method\n", gsl_min_fminimizer_name (s));
  fprintf (logg,"%5s [%9s, %9s] %9s %10s\n","iter", "lower", "upper", "min","err");
  fprintf (logg,"%5d [%.4e, %.4e] %.4e %+.4e\n", iter, a, b, m, b - a);
#endif


  do {
    iter++;
    status = gsl_min_fminimizer_iterate (s);
    m = gsl_min_fminimizer_x_minimum (s);
    a = gsl_min_fminimizer_x_lower (s);
    b = gsl_min_fminimizer_x_upper (s);
    status = gsl_min_test_interval (a, b, 1e-9, 0.0);
#if MAJO>1
    if (status == GSL_SUCCESS) fprintf (logg,"Converged:\n");
    fprintf (logg,"%5d [%.4e, %.4e] %.4e %+.4e\n", iter, a, b, m, b - a);
#endif
  } while (status == GSL_CONTINUE && iter < max_iter);
  gsl_min_fminimizer_free (s);
  gsl_set_error_handler (NULL);
  error=status;
  if (error==0) return(minimize_func (m, F.params)); else return(0);
}

//returns the rank of the localized solution, the corresponding energy and the tag
int give_rank(double *sv, double& E, char& tag)
{
  int ranks[4]={0,0,0,0};
  tag='N';
  //first check the energy zero for possible majoranas
  p.E=E=0;
  double* sv0=new double[4];
  int rank0=phase(sv0);
  ranks[rank0]++;
  if (rank0==1) tag='M';
  if (rank0==2) tag='F';
  
  //check possible E nonzero minima
  //first: find the minima within the gap
  double* mins=new double[10], step;
  int hmn=prescan(mins, step), rank;
#if MAJO > 0
  fprintf(logg,"at zero energy the rank is %i, from sv[0]=%e, sv[1]=%e. In addition, prescan gave %i minima\n", rank0, sv0[0], sv0[1], hmn);
#endif

  //check each minimum
  for (int i=0;i<hmn;i++) {
    double m=mins[i];
    int error=0;
    minimize(mins[i]-step,mins[i]+step,mins[i],error,0);
#if MAJO > 0
    if (error==0) fprintf(logg,"%i-th local minimum, estimated at %e, nailed down: yes, at %e\n",i,m,mins[i]);
    else fprintf(logg,"%i-th local minimum, estimated at %e, nailed down:no, WARNING: gsl_error code=%i\n",i,m,error);
#endif
    if (error==0) {
      rank=phase(sv);
      ranks[rank]++;
      if (rank>0) {
	E=p.E;
	tag='f';  
	if (rank>1) {
	  fprintf(logg,"WARNING: new phase found, doubt it!\n");
          FILE* file=fopen("data/majorana-new-phase","a+");
          fprintf(file,"more than two fermions inside the gap above E=0 suspected at s=%.14e, n=%.14e, z=%.14e, mu=%.14e, theta=%.14e\n", p.s, p.n, p.z, p.mu, p.theta);
          fprintf(file,"\tfrom ranks=(%i,%i,%i), the %i minima are at energies: ", ranks[0], ranks[1], ranks[2], hmn);
          for (int i=0;i<hmn;i++) {
	    fprintf(file,"%e, ",mins[i]);
            p.E=mins[i]=0;
            phase(sv,3,file);
	  }
          fprintf(file,"zero energy: rank=%i, sv[0]=%e, sv[1]=%e\n",rank0, sv0[0],sv0[1]);
          fclose(file);
          //exit(1);
	}
      }
#if MAJO > 0
    fprintf(logg,"\tthe two smallest sv are %e and %e and the rank is %i\n", sv[0], sv[1], rank);
#endif
    }
  }
#if MAJO > 0
  fprintf(logg,"zero and non-zero minima gave ranks: 0x%i times, 1x%i times, 2x%i times\n", ranks[0], ranks[1], ranks[2]);
#endif
  if (ranks[1]*ranks[2]>0) {
    fprintf(logg,"WARNING: multiple solutions found!\n");
    FILE* file=fopen("data/majorana-multiple","a+");
    fprintf(file,"multiple solutions suspected at s=%.14e, n=%.14e, z=%.14e, mu=%.14e, theta=%.14e\n", p.s, p.n, p.z, p.mu, p.theta);
    fprintf(file,"\tfrom ranks=(%i,%i,%i), the %i minima are at energies: ", ranks[0], ranks[1], ranks[2], hmn);
    for (int i=0;i<hmn;i++) fprintf(file,"%e, ",mins[i]);
    fprintf(file,"zero energy: rank=%i, sv[0]=%e, sv[1]=%e\n",rank0, sv0[0],sv0[1]);
    p.E=0;
    phase(sv,3,file);
    fclose(file);
    //exit(1);
  }
  for (int r=2;r>-1;r--) if (ranks[r]>0) return(r);
  delete sv0;
  delete mins;
  return(0);
}

double transmutation_gap(double& zlow, double& zhigh)
{
  char tag;
  //first the range, remmember highest occurence of majorana, lowest of frac. charge
  double F=-100, M=-100, E, * sv=new double[4];
  p.z=zlow;
  do {
    int rank=give_rank(sv, E, tag);
#if MAJOT > 1
  fprintf(logg,"trans step 1: at z=%e the returned tag was %c: (actual values of M=%e and F=%e will be updated accrdingly)\n", p.z, tag, M, F);
#endif
    if (tag=='M') M=p.z;
    if (tag=='F' || tag=='f') {F=p.z;break;}
    p.z+=(zhigh-zlow)/100;
  } while (p.z<zhigh);

#if MAJOT > 0
  fprintf(logg,"trans. gap initial scan gave M=%e and F=%e\n", M, F);
#endif

  if (F==-100 || M==-100) return(100); //gap not defined as one of the phases was not seen
  
  //identify end of M
  zlow=M;zhigh=F;
  do {
    p.z=(zlow+zhigh)/2;
    int rank=give_rank(sv, E, tag);
    if (tag=='M') zlow=p.z; else zhigh=p.z;
#if MAJOT > 1
  fprintf(logg,"trans step 2: at z=%e the returned tag was %c: (actual values of zlow=%e and zhigh=%e)\n", p.z, tag, zlow, zhigh);
#endif
  } while (zhigh-zlow>1e-8);
  M=p.z;
#if MAJOT > 0
  fprintf(logg,"trans. gap nailed down end of M region at M=%e\n", M);
#endif


  //identify beginning of F
  zlow=M;zhigh=F;
  do {
    p.z=(zlow+zhigh)/2;
    int rank=give_rank(sv, E, tag);
    if (tag=='F' || tag=='f') zhigh=p.z; else zlow=p.z;
#if MAJOT > 1
  fprintf(logg,"trans step 3: at z=%e the returned tag was %c: (actual values of zlow=%e and zhigh=%e)\n", p.z, tag, zlow, zhigh);
#endif
  } while (zhigh-zlow>1e-8);
  F=p.z;
#if MAJOT > 0
  fprintf(logg,"trans. gap nailed down start of F region at F=%e\ntrans gap ended with M=%e, F=%e and the gap=%e\n", F, M, F, F-M);
#endif

  zlow=M;
  zhigh=F;
  delete sv;
  return(F-M);
}

//minimize the M to F gap wrt the angle theta, input - theta estimate, output theta min, and the corresponding gap
double transmute(double & theta)
{
  int step=0;
  double theta_step=M_PI/100;
  p.theta=theta-theta_step*5;
  double gpprev[100], thetas[100];
  int inc=0, dec=0;
  bool long_dec=false, long_inc=false, success=false;
  do {//find the local minimum
      double zlow=0,zhigh=2;
      double gp=transmutation_gap(zlow, zhigh);
      if (step>0) {
	if (gp<gpprev[step-1]) {inc=0; dec++;} else {inc++; dec=0;}
      }
      if (dec>3) long_dec=true;
      if (inc>3) long_inc=true;
#if MAJOT>0
      fprintf(logg,"transmute step %i, theta=%e, gap=%e (zlow=%e, zhigh=%e), resulting in inc=%i, dec=%i\n", step, p.theta, gp, zlow, zhigh, inc, dec);
#endif
      if (long_dec && long_inc) {success=true;break;}
      gpprev[step]=gp;
      thetas[step]=p.theta;
      p.theta+=theta_step;
  } while (step++<50);

  //minimize the polynomial interpolation - the minimum appeared inc spots ago
  if (success) {
#if MAJOT > 0
      fprintf(logg,"local trans minimum will be searched for from the following data set (step=%i, inc=%i)\n",  step, inc);
#endif
    for (int i=0;i<3;i++) {
      coor[i]=thetas[step-inc+i-1];
      val[i]=gpprev[step-inc+i-1];
#if MAJOT > 0
      fprintf(logg,"\ti=%i, coor=%e, val=%e\n", i, coor[i], val[i].real());
#endif
    }
    theta=thetas[step-inc];
    int error=0;
    double gp_min=minimize(thetas[step-inc-1],thetas[step-inc+1],theta,error,2);
    if (error==0) {
#if MAJOT > 0
      fprintf(logg,"local trans minimum, estimated at %e, nailed down: yes, at %e\n", thetas[step-inc],theta);
#endif
      return(gp_min);
    }
    else {
#if MAJOT > 0
      fprintf(logg,"local trans minimum, estimated at %e, nailed down:no, WARNING: gsl_error code=%i\n", theta, error);
#endif
      return(-1);
    }
  }
#if MAJOT > 0
  fprintf(logg,"no local trans minimum found\n");
#endif
  return(-1);
}

int main()
{
  const char* app="_final_ntheta_chosen";
  const char* name[]={"majorana4"};
  output_class output(1,"data/",name,"log");
  //output.set_appendix(app);

  output.clear_files();                    //clear the output files
  output.open_files();                     //clear the output files
    
  p.init(); 
  
  { 
    p.s=0.5;p.n=1.5;p.z=-2;p.mu=1.75;p.theta=0;p.E=0;
    double* sv=new double[4], E;
    char tag;
    phase(sv,3);
    //int rank=give_rank(sv, E, tag);
    delete sv;
    exit(1);
  }

  progress_bar_prototype progress_bar(app);
  progress_bar.start(true);
  //progress_bar.add((int) (1.5*2*50*(4/0.05+1)),0);
  progress_bar.add((int) (401*2001),0);

  
  //p.s=1;p.n=0.0;p.z=1.0;p.mu=0.0;
  //double theta=M_PI;
  //double gap_min=transmute(theta);
  //exit(1);

  /*double theta_in=M_PI, theta_out=0;
  for (p.mu=0;p.mu<1.5;p.mu+=0.01) {
    theta_out=theta_in;
    double gap_min=transmute(theta_out);
    fprintf(output.file(0),"%e %e %e\n", p.mu, theta_out/M_PI, gap_min);
    if (gap_min!=-1) theta_in=theta_out;
    progress_bar.add(0,1);
    fflush(0);
  }
  exit(1);
*/
  
  p.s=0.5+1e-4;p.z=1.0+1e-5;p.mu=0.2+1e-5;

  int tagN;
  char tag;
  double * sv=new double[4], E;
  //for (p.z=-1;p.z<1.001;p.z+=0.2) {
    //for (p.mu=-1;p.mu<1.001;p.mu+=0.2) {
      for (p.n=-2;p.n<2.001;p.n+=0.01) {
        for (p.theta=0;p.theta<2.001*M_PI;p.theta+=M_PI*0.001) {
	  int rank=give_rank(sv, E, tag);
	  if (tag=='M') tagN=2; else if (tag=='f') tagN=3; else if (tag=='F') tagN=1; else tagN=0;
	  fprintf(output.file(0),"%e %e %e %e %e %i %e %e %e %c %i\n", p.s, p.z, p.n, p.mu, p.theta, rank, E, gap(), sv[0], tag, tagN);
	  phase(sv);
	  progress_bar.add(0,1);
	}
        fprintf(output.file(0),"\n");
      }
      fprintf(output.file(0),"\n\n");
    //}
  //}
  delete sv;
  progress_bar.finished();
  return(0);
}

  /*
  p.E=0;
  p.mu=1.51;
  p.s=0.5;
  p.z=0.0;
  p.n=-1.5;
  p.theta=0;
   
  double* sv=new double[4];
  double* mins=new double[10], step, res;
  for (int i=0;i<1000;i++) {
    p.z=(-1+generuj()*2);
    p.n=(-1+generuj()*2);
    p.theta=generuj()*2*M_PI;
    int hmn=prescan(mins, step);
    int success=minimize(mins[0]-2*step,mins[0],mins[0]-step);
    fprintf(logg,"minimization ended with E=%e\n",p.E);
    int rank=phase(sv);
    double Eformula=abs(p.n*p.z*sin(p.theta))/sqrt(p.z*p.z+p.n*p.n-2*p.z*p.n*cos(p.theta));
    fprintf(output.file(0),"%e %e %e %e %e %e\n", p.z, p.n, p.theta, p.E, sv[3], abs(p.E-Eformula));
    printf("step %i\n",i);
  }*/


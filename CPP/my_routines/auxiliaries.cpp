#include "auxiliaries.h"

#define AUXIL 0

char buffer[1000];
//FILE* logg=fopen("data/log","a+");
FILE* logg=0;


void splachni(FILE* f,const char* co,int kam)
{
   bool ciel[3];
   int bit_mask=2*2*2/2;    /*najvyssi bit*/
   
   for (int i=2;i>=0;i--)
   {
      if (kam/bit_mask) ciel[i]=true;   /*najvyssi bit je plny*/
         else ciel[i]=false;
      kam=kam-bit_mask*ciel[i];
      bit_mask=bit_mask/2;             /*dalsi nizsi bit*/
   }

   if (ciel[0]) printf("%s",co);          /*na obrazovku*/
   if (ciel[1]) {
     if (f==0) {
       printf("file is not opened!\n");
       exit(1);
     }
     fprintf(f,"%s",co);       /*do s�boru predaneho argumentom*/
   }
   
   if (ciel[2]) {
     if (logg==0) {
       printf("log file is not opened!\n");
       exit(1);
     }
     fprintf(logg,"%s",co);    /*do suboru log*/
   }
}


void message(const char* message)
{
  sprintf(buffer,"%s",message);
  splachni(logg,buffer,where_to_print_critical);
}

void message(const char* message, int where_to)
{
  sprintf(buffer,"%s",message);
  splachni(logg,buffer,where_to);
}

void message(FILE* file,const char* message, int where_to)
{
  sprintf(buffer,"%s",message);
  splachni(file,buffer,where_to);
}

void other_name(const char * name) {
  if (logg!=0) fclose(logg);
  logg=fopen(name,"a+");
  if (logg==0) {
    printf("log file not opened successfully!\n");
    exit(1);
  }
}

void clear_file(const char * name)
{
  FILE* file=fopen(name, "w+");
  if (file!=0) fclose(file);
}

output_class::output_class(int NF_, const char* path, const char** name, const char* logname)
{
  NF=NF_;
  filenamebase=new char*[NF];
  filename=new char*[NF];
  files=new FILE*[NF];
  for (int i=0;i<NF;i++) {
    files[i]=0;
    filenamebase[i]=new char[strlen(name[i])+strlen(path)+1];
    filename[i]=new char[strlen(name[i])+strlen(path)+100];
    strcpy(filenamebase[i],path);
    strcat(filenamebase[i],name[i]);
    strcpy(filename[i],filenamebase[i]);
  }
  logfilenamebase=new char[strlen(path)+strlen(logname)+1];
  logfilename=new char[strlen(path)+strlen(logname)+100];
  strcpy(logfilenamebase, path);
  strcat(logfilenamebase, logname);
  strcpy(logfilename, logfilenamebase);
}

output_class::~output_class()
{
  close_files();
  for (int i=0;i<NF;i++) {
    delete filenamebase[i];
    delete filename[i];
  }
  delete files;
  delete filenamebase;
  delete filename;
  delete logfilenamebase;
  delete logfilename;
}

void output_class::clear_files(int which)
{
  if (which!=-1) clear_file(filename[which]); else {
    for (int i=0;i<NF;i++) clear_file(filename[i]);
    clear_file(logfilename);
  }
}

void output_class::close_files()
{
  for (int i=0;i<NF;i++) {
    if (files[i]!=0) {
      //printf("file %i pointer %p will be closed\n",i,files[i]);
      fclose(files[i]);
      files[i]=0;
    }
    //else printf("file %i pointer %p is already closed\n",i, files[i]);
  }
}

void output_class::open_files()
{
  close_files();
  for (int i=0;i<NF;i++) {
    char string[1000];
  loop:
    files[i]=fopen(filename[i],"a");
    if (files[i]==0) {
      printf("file %s not opened successfully for append\n",filename[i]);
      exit(1);
      //system("pwd");
      //sprintf(string,"ls %s\n",filename[i]); system(string);
      //files[i]=fopen(filename[i],"r");
      //if (files[i]==0) printf("file %s not opened successfully for reading\n",filename[i]);
      //else fclose(files[i]);
      //system("sleep 1");
      //goto loop;
    }
    //else if (logg!=0) fprintf(logg,"file %s opened successfully pointing to %p\n",filename[i], files[i]);
  }
  other_name(logfilename);
  for (int i=0;i<NF;i++) fprintf(logg,"file(%i): %s was opened\n",i,filename[i]);
}

void output_class::set_appendix(const char* tag, int p)
{
  char aux[100];
  if (p<0) aux[0]=0; else itoa(p,aux);
  for (int i=0;i<NF;i++) {
    strcpy(filename[i],filenamebase[i]);
    strcat(filename[i],tag);
    strcat(filename[i],aux);
  }
  strcpy(logfilename,logfilenamebase);
  strcat(logfilename,tag);
  strcat(logfilename,aux);
}

FILE* output_class::file(int i)
{
  if (i<0 || i>NF-1) {printf("asking for a non-existing file %i (out of %i)\n",i, NF-1);exit(1);}
  if (files[i]==0) printf("warning: file %s not open\n",filename[i]);
  return(files[i]);
}


//!************************************************

/*gama funkcia s (polo)ciselnym argumentom
gama(i,j)=gama funkcia argumentu (i/j); j={1,2}*/
double gama(int i, int j)
{
   double c=1;
   while (i>j) {
      i-=j;
      c*=(i+0.0)/j;
   }
   if (i!=j) c=c*sqrt(M_PI); /*ostala gama (1/2)*/
   return(c);
}

double erfi(double x)
{
  double term=2/sqrt(M_PI)*x;
  double sum=term;
  int n=0;
  do {
    n++;
    term*=x*x/n;
    sum+=term/(2*n+1);
    if (term/sum<1e-6 || n>100) break;
  } while (true);
  //fprintf(logg,"check of erfi: erfi(%e)=%e\n",x,sum);
  return(sum);
}


double fac(int n) {return(gama(n+1,1));}

double pown(double x, int n) 
{
  double res=1;
  for (int i=0;i<n;i++) res*=x;
  return(res); 
}

//converts integer into ascii
void itoa(int i,char *num2)
{
  bool sign=false;
  int count=0,j;
  char num[16];

  if (i<0) {
    i=-i;
    sign=true;
    num2[0]='-';
  }

  do {
    j=i%10;
    i/=10;
    num[count++]='0'+j;
  } while(i>0);
  num[count]=0;
  j=strlen(num);
  for (int k=0;k<j;k++) num2[k+sign]=num[j-k-1];
  num2[j+sign]=0;
}

int quadratic(double p, double q, double& x1, double& x2)
{
  double D=p*p-4*q;
  if (D<0) return(-1);
  D=sqrt(D);
  x1=(-p+D)/2;
  x2=(-p-D)/2;
  #if AUXIL > 1
  double res1=x1*x1+p*x1+q, res2=x2*x2+p*x2+q;
  fprintf(logg,"real roots of the quadratic equation %.5e ,%.5ei) give residua %e %e\n",x1,x2,res1,res2);
  #endif
  if (abs(x1)>abs(x2)) {D=x1;x1=x2;x2=D;}
  return(0);
}

int complex_bilinear(complex<double>* coef, complex<double>& x)
{
  //scale coefs
  int which=0;
  for (int i=1;i<3;i++) if (abs(coef[i])>abs(coef[which])) which=i;
  for (int i=0;i<3;i++) coef[i]/=abs(coef[which]);
  coef[1]=conj(coef[1]);
  
  double ap=(coef[0]+coef[1]).real();
  double bp=(coef[0]+coef[1]).imag();
  double am=(coef[0]-coef[1]).real();
  double bm=(coef[0]-coef[1]).imag();
  double a3=coef[2].real();
  double b3=coef[2].imag();
  
  double D1=-bm*bp-am*ap;
  double D2=bm*a3-ap*b3;
  
  double eps=1e-14;
  int info=0;
  double i=0,r=0;
  if (abs(D1)<eps) {
    if (abs(D2)<eps) {//infinity solutions
      info=1;
      if (abs(bm)>0) {i=0;r=-b3/bm;}
      else if (abs(am)>0) {r=0;i=-b3/am;}
      else {i=r=0;}
    }
    else {//no solution
      i=r=0;
      info=-1;
    }
  } else {
    r=(bp*b3+am*a3)/D1;
    i=(ap*b3-bm*a3)/D1;
    info=0;
  }
  x=complex<double> (r,i);
  #if AUXIL > 1
  double res=abs(coef[0]*x+conj(coef[1])*conj(x)+coef[2]);
  fprintf(logg,"complex root of the complex bilinear equation (%.5e%+.5ei) gives residuum %e (info=%i)\n",r,i,res,info);
  #endif
  return(info);
}

//finds roots of quadratic equation with complex coefficients c0 |x|^2 + c1 x + c2 x^* + c3 =0
//info = -1 if no solution, info = 0 if ok, info = 1 if infinity solutions
int complex_quadratic(complex<double>* coef, complex<double>& x1, complex<double>& x2)
{
  //fprintf(logg,"complex quadratic: abs(c[0])=%.5e, abs(c[1])=%.5e, abs(c[2])=%.5e, abs(c[3])=%.5e\n", abs(coef[0]), abs(coef[1]), abs(coef[2]), abs(coef[3]));

  double eps=1e-14;
  int info;
  //check if really quadratic equation, rescale coefficients
  int which=0;
  for (int i=1;i<4;i++) if (abs(coef[i])>abs(coef[which])) which=i;
  for (int i=0;i<4;i++) coef[i]/=abs(coef[which]);
  if (abs(coef[0])<eps) {info=complex_bilinear(coef+1,x1); x2=x1; return(info);}
  for (int i=1;i<4;i++) coef[i]/=coef[0];
  coef[0]=1.0;
  
  //fprintf(logg,"complex quadratic: abs(c[0])=%.5e, abs(c[1])=%.5e, abs(c[2])=%.5e, abs(c[3])=%.5e\n", abs(coef[0]), abs(coef[1]), abs(coef[2]), abs(coef[3]));

  //new variables 
  double ap=(coef[1]+conj(coef[2])).real();
  double bp=(coef[1]+conj(coef[2])).imag();
  double am=(coef[1]-conj(coef[2])).real();
  double bm=(coef[1]-conj(coef[2])).imag();
  double a3=coef[3].real();
  double b3=coef[3].imag();
  double p, q, aux, r1=0, r2=0, i1=0, i2=0;
   
  //fprintf(logg,"complex quadratic: ap=%.5e, am=%.5e, bp=%.5e, bm=%.5e, a3=%.5e, b3=%.5e\n",ap, am, bp, bm, a3, b3);

  //decide upon the imaginary part of the equation
  if (abs(am)<eps && abs(bm)<eps) {
    if (abs(b3)<eps) {//no condition from imaginary part
      aux=ap*ap/4+bp*bp/4-a3;
      if (aux<0) info=-1; //no solution
      else {
        if (abs(aux)<eps) info=0; else info=1;
        i1=i2=bp/2;
        r1=-ap/2+sqrt(aux);
        r2=-ap/2-sqrt(aux);
      }
      fprintf(logg,"c.q. branched to single equation : aux=%e, r1=%e, r2=%e, i1=%e, i2=%e, info=%i\n",aux, r1, r2, i1, i2, info);
    }
    else {
      info=-1;//no solution
      fprintf(logg,"c.q. branched to no solution : b3=%e, info=%i\n",b3, info);
    }
  }
  else {
    if (abs(am)>abs(bm)) {//I=-(b3+bm R)/am
      p=2*b3*bm/(am*am)+ap+bp*bm/am;
      q=b3*b3/(am*am)+bp*b3/am+a3;
      aux=1+bm*bm/(am*am);
      info=quadratic(p/aux, q/aux, r1, r2);
      i1=-(b3+bm*r1)/am;
      i2=-(b3+bm*r2)/am;
      fprintf(logg,"c.q. branched to larger am : p=%e, q=%e, aux=%e, r1=%e, r2=%e, i1=%e, i2=%e, info=%i\n",p, q, aux, r1, r2, i1, i2, info);
    }
    else {//R=-(b3+am I)/bm
      p=2*am*b3/(bm*bm)-ap*am/bm-bp;
      q=b3*b3/(bm*bm)-ap*b3/bm+a3;
      aux=1+am*am/(bm*bm);
      info=quadratic(p/aux, q/aux, i1, i2);
      r1=-(b3+am*i1)/bm;
      r2=-(b3+am*i2)/bm;
      fprintf(logg,"c.q. branched to larger bm : p=%e, q=%e, aux=%e, r1=%e, r2=%e, i1=%e, i2=%e, info=%i\n",p, q, aux, r1, r2, i1, i2, info);
    }
  }

  x1=complex<double>(r1,i1);
  x2=complex<double>(r2,i2);
  if (r1*r1+i1*i1>r2*r2+i2*i2) {complex<double> aux2=x1;x1=x2;x2=aux2;}
  if (info==-1) {x1=x2=r1=r2=i1=i2=0;}
  #if AUXIL > 1
  double res1=abs(coef[0]*x1*conj(x1)+coef[1]*x1+coef[2]*conj(x1)+coef[3]);
  double res2=abs(coef[0]*x2*conj(x2)+coef[1]*x2+coef[2]*conj(x2)+coef[3]);
  fprintf(logg," complex roots of the complex quadratic equation (%.5e%+.5ei), (%.5e%+.5ei) give residua %e %e (info=%i)\n",r1,i1,r2,i2,res1,res2,info);
  #endif
  return(info);
}

extern "C" void F77NAME(zgeev) (char* JOBVL, char* JOBVR, int* N, complex<double>* A, int* LDA, complex<double>* W, complex<double>* VL,  int* LDVL, complex<double>* VR, int* LDVR, complex<double>* WORK, int* LWORK, double* RWORK, int* INFO );

int cubic(double* coef, complex<double>* root)
{
  //rescale the coefficients
  if (abs(coef[0])==0) {printf("zero coefficient in cubic equation\n");return(0);}
  double c0=coef[3]/coef[0];
  double c1=coef[2]/coef[0];
  double c2=coef[1]/coef[0];
  
  complex<double> CM[9];//construct the companion matrix
  for (int i=0;i<9;i++) CM[i]=0;
  CM[0+2]=-c0;
  CM[3+0]=1.0; CM[3+2]=-c1;
  CM[6+1]=1.0; CM[6+2]=-c2;
  
  char job='N';
  int N=3, lwork=18, info;
  complex<double> vl[9], vr[9], work[18];
  double rwork[6];
  F77NAME(zgeev)(&job, &job,  &N, CM, &N, root, vl, &N, vr, &N, work, &lwork, rwork, &info);
  if (info!=0) {printf("cubic equation was not solved, info=%i\n", info);exit(1);}
  complex<double> aux;//sort the roots
  if (abs(root[0].imag())>abs(root[1].imag())) {aux=root[1]; root[1]=root[0];root[0]=aux;}
  if (abs(root[1].imag())>abs(root[2].imag())) {aux=root[2]; root[2]=root[1];root[1]=aux;}
  if (abs(root[0].imag())>abs(root[1].imag())) {aux=root[1]; root[1]=root[0];root[0]=aux;}
  N=0;//number of real solutions
  for (int i=0;i<3;i++) if (root[i].imag()/root[i].real()<1e-8) N++;

  //check
#if AUXIL > 1
  complex<double> y[4]={0}, tot=0;
  for (int i=0;i<3;i++) {
    for (int j=0;j<4;j++) y[i]+=coef[j]*pow(root[i],3-j);
    tot+=y[i];
  }
  if (abs(tot)>1e-14) {
    fprintf(logg,"cubic equation roots total residuum %.2e (residua %.2e %.2e %.2e)\n",abs(tot), abs(y[0]), abs(y[1]), abs(y[2]) );
  }
#endif
  return(N);
}


int cubic2(double* coef, complex<double>* root)
{
  double const eps=1e-15;
  double ext=abs(coef[0]);
  for (int i=1;i<4;i++) if (abs(coef[i])>ext) ext=abs(coef[i]);
  double a=coef[0]/ext;
  double b=coef[1]/ext;
  double c=coef[2]/ext;
  double d=coef[3]/ext;

  if (abs(a)<eps) {fprintf(logg,"zero coefficient a in cubic equation\n");exit(1);}

  double aux1=2*b*b*b-9*a*b*c+27*a*a*d;
  double aux2=b*b-3*a*c;

  double Q;
  if (abs(aux2)<eps) Q=aux1; else {
    Q=aux1*aux1-4*aux2*aux2*aux2;
    if (Q>0) Q=sqrt(Q); else Q=0;
  }
  double C=pow( (Q+aux1)/2 ,1.0/3 );

  int N=0;
  if (abs(Q)<eps && abs(aux2)<eps) {
    root[0]=root[1]=root[2]=-b/(3*a);
    N=3;
  } else if (abs(Q)<eps) {
    root[0]=-(9*a*a*d-4*a*b*c+b*b*b)/(a*aux2);
    root[1]=root[2]=-(b*c-9*a*d)/(2*aux2);
    N=3;
  } else {
    complex<double> aux3=complex<double>(1,sqrt(3));
    root[0]=-b/(3*a)-C/(3*a)-aux2/(C*3*a);
    root[1]=-b/(3*a)+C*aux3/(6*a)+conj(aux3)*aux2/(6*a*C);
    root[2]=conj(root[1]);
    N=1;
  }
  #if AUXIL > 1
  complex<double> y[4]={0}, tot=0;
  for (int i=0;i<3;i++) {
    for (int j=0;j<4;j++) y[i]+=coef[j]*pow(root[i],3-j);
    tot+=y[i];
  }
  if (abs(tot)>1e-14) {
    fprintf(logg,"cubic equation roots total residuum %.2e (residua %.2e %.2e %.2e)\n",abs(tot), abs(y[0]), abs(y[1]), abs(y[2]) );
  }
  #endif
  return(N);
}

//returns x coth(x) as (1+u)/(1+d)
void xcothx(double x, double &u, double &d)
{
  if (abs(x)>0.1) {d=0;u=x/tanh(x)-1;return;}
  double x2=x*x;
  double up=u=x2/2;
  double dp=d=x2/6;
  int fac=3;  
  do {
    up*=x2/(fac*(fac+1));
    u+=up;
    dp*=x2/((fac+1)*(fac+2));
    d+=dp;
    fac+=2;
  } while (up>1e-16);
}
    
//Brillouin function
double brillouin(int J, double x)
{
  if (x==0) return(0);
  double u1,d1,u2,d2;
  xcothx((2*J+1)/(2*J)*x, u1, d1);
  xcothx(1/(2*J)*x, u2, d2);
  double res=u1+d2+u1*d2-u2-d1-u2*d1;
  res/=x*(1+d1)*(1+d2);
  return(res);
}

/*
double norm(const complex<double>& x)
{
  double x1=x.real();
  double x2=x.imag();
  return(x1*x1+x2*x2);
}

double abs(const complex<double>& x)
{
  return(sqrt(norm(x)));
}
*/


//gives i mod n such that the result is from region
//<-N/2,N/2> or <-N/2,N/2+1> for N even and N+1 odd
int mod_pm(int i,int n)
{
  int j=pos_mod(i,n);
  if (j>(n+1)/2) j-=n;
  return(j);
}

//the nearest integer number
int my_round(double x) {
  double up,down;
  up=ceil(x);
  down=floor(x);
  if (up-x<x-down) return((int) up); else return((int) down);
}


#define  A1  (-3.969683028665376e+01)
#define  A2   2.209460984245205e+02
#define  A3  (-2.759285104469687e+02)
#define  A4   1.383577518672690e+02
#define  A5  (-3.066479806614716e+01)
#define  A6   2.506628277459239e+00

#define  B1  (-5.447609879822406e+01)
#define  B2   1.615858368580409e+02
#define  B3  (-1.556989798598866e+02)
#define  B4   6.680131188771972e+01
#define  B5  (-1.328068155288572e+01)

#define  C1  (-7.784894002430293e-03)
#define  C2  (-3.223964580411365e-01)
#define  C3  (-2.400758277161838e+00)
#define  C4  (-2.549732539343734e+00)
#define  C5   4.374664141464968e+00
#define  C6   2.938163982698783e+00

#define  D1   7.784695709041462e-03
#define  D2   3.224671290700398e-01
#define  D3   2.445134137142996e+00
#define  D4   3.754408661907416e+00

#define P_LOW   0.02425
/* P_high = 1 - p_low*/
#define P_HIGH  0.97575

double normsinv(double p)
{
  double x=0, q, r, u, e;
  if ((0 < p )  && (p < P_LOW)){
     q = sqrt(-2*log(p));
     x = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
  }
  else{
        if ((P_LOW <= p) && (p <= P_HIGH)){
           q = p - 0.5;
           r = q*q;
           x = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
        }
        else{
                if ((P_HIGH < p)&&(p < 1)){
                   q = sqrt(-2*log(1-p));
                   x = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
                }
        }
  }

// If you are compiling this under UNIX OR LINUX, you may uncomment this block for better accuracy.
  if(( 0 < p)&&(p < 1)){
     e = 0.5 * erfc(-x/sqrt(2)) - p;
     u = e * sqrt(2*M_PI) * exp(x*x/2);
     x = x - u/(1 + x*u/2);
  }
  return(x);
}


double erfinv(double p)
{
  const double a=8.0*(M_PI-3.0)/(3.0*M_PI*(4.0-M_PI));
  if (p==1) return(1e300);
  if (p==-1) return(-1e300);
  double a1=log(1-p*p);
  double a2=(2/M_PI/a)+a1/2;
  if (p>0) return(sqrt(sqrt(a2*a2-a1/a)-a2)); else return(-(sqrt(a2*a2-a1/a)-a2));
}

//find first n solutions of equation J_l(x)=J_{l+1}(x)
void bessel_solver(int l, int n, double* res)
{
  if (l<0) l=-l;
  if (n<1) n=1; else if (n>100) n=100;
  double x=0, x_prev, step=0.01, dif=0, dif_prev;
  int pntr=0;
  
  do {
    dif_prev=dif;
    x_prev=x;
    x+=step;
    dif=gsl_sf_bessel_Jn(l,x)-gsl_sf_bessel_Jn(l+1,x);
    if (dif*dif_prev<0) {//solution identified, start the interval halving
      double xl=x_prev, xr=x, xm;
      do {
        xm=(xr+xl)/2;
        double dif_loc=gsl_sf_bessel_Jn(l,x)-gsl_sf_bessel_Jn(l+1,x);
        if (dif_loc*dif_prev<0) xr=xm; else xl=xm;
        if (xr-xl<1e-12) break;
      } while(true);
      res[pntr++]=x=xm;
      dif=0;
      if (pntr>2) step=(res[pntr-1]-res[pntr-2])/4; else step=res[0]/4;
    }
    if (pntr==n) break;
  } while(true);
}
          
//sort according to x in 'a'/'d'e scending order
void picsrt_inplace(int N, int* x, char order)
{
  int i;
  for (int j=1;j<N;j++) {
    int xaux=x[j];
    for (i=j-1;i>-1;i--) {
      if (order=='a' && x[i]<xaux) break;
      if (order=='d' && x[i]>xaux) break;
      x[i+1]=x[i];
    }
    x[i+1]=xaux;
  }
}

//sort according to x in 'a'/'d'e scending order
void picsrt_inplace(int N, double* x, char order)
{
  int i;
  for (int j=1;j<N;j++) {
    double xaux=x[j];
    for (i=j-1;i>-1;i--) {
      if (order=='a' && x[i]<xaux) break;
      if (order=='d' && x[i]>xaux) break;
      x[i+1]=x[i];
    }
    x[i+1]=xaux;
  }
}

//sort according to x in 'a'/'d'e scending order
void picsrt_inplace(int N, complex<double>* x, char method, char order)
{
  int i;
  for (int j=1;j<N;j++) {
    complex<double> xaux=x[j];
    for (i=j-1;i>-1;i--) {
      bool test=false;
      switch (method) {
        case 'R' : {test=(abs(x[i].real())<abs(xaux.real())); break;}
        case 'r' : {test=(x[i].real()<xaux.real()); break;}
        case 'I' : {test=(abs(x[i].imag())<abs(xaux.imag())); break;}
        case 'i' : {test=(x[i].imag()<xaux.imag()); break;}
        case 'a' : {test=(abs(x[i])<abs(xaux)); break;}        
      }
      if (order=='a' && test) break;
      if (order=='d' && !test) break;
      x[i+1]=x[i];
    }
    x[i+1]=xaux;
  }
}

//sort according to x in 'a'/'d'e scending order
//create the index and rank arrays, the input array is not touched
void picsrt_index(int N, double* x, int * index, int*rank, char order)
{
  for (int i=0;i<N;i++) index[i]=i;
  int i, iaux;
  for (int j=1;j<N;j++) {
    iaux=index[j];
    double xaux=x[iaux];
    for (i=j-1;i>-1;i--) {
      if (order=='a' && x[index[i]]<xaux) break;
      if (order=='d' && x[index[i]]>xaux) break;
      index[i+1]=index[i];
    }
    index[i+1]=iaux;
  }
  if (rank!=0) for (int j=0;j<N;j++) rank[index[j]]=j;
}

void picsrt_index(int N, complex<double>* x, int * index, int* rank, char method, char order)
{
  for (int i=0;i<N;i++) index[i]=i;
  int i, iaux;
  for (int j=1;j<N;j++) {
    iaux=index[j];
    complex<double> xaux=x[iaux];
    for (i=j-1;i>-1;i--) {
      bool test=false;
      switch (method) {
        case 'R' : {test=(abs(x[index[i]].real())<abs(xaux.real())); break;}
        case 'r' : {test=(x[index[i]].real()<xaux.real()); break;}
        case 'I' : {test=(abs(x[index[i]].imag())<abs(xaux.imag())); break;}
        case 'i' : {test=(x[index[i]].imag()<xaux.imag()); break;}
        case 'a' : {test=(abs(x[index[i]])<abs(xaux)); break;}        
      }
      if (order=='a' && test) break;
      if (order=='d' && !test) break;
      index[i+1]=index[i];
    }
    index[i+1]=iaux;
  }
  if (rank!=0) for (int j=0;j<N;j++) rank[index[j]]=j;
}

void sort_from_index(int N, double * x, int * index)
{
  double* aux=new double[N];
  for (int i=0;i<N;i++) aux[i]=x[i];
  for (int i=0;i<N;i++) x[i]=aux[index[i]];
  delete[] aux;
}

//returns a matrix (vector in case of single derivatives when asked for by "vec"=true)
//representing a two dimensional differential operator
//dx^i dy^j using 2l+1 by 2l+1 points (including the central one)
//apart from the coefficient 1/ (hx^i hy^j)
//coefficients are stored in the matrix such that the desired operator is proportional to
//DF=sum_{i,j=-l,l} c[i*(2*l+1)+j] * f(i,j)
void differential_operator(int i, int j, int l, double coef, double* out) 
{
#define MAX_L 4

  double result[(2*MAX_L+1)*(2*MAX_L+1)];  
  for (int k=0;k<(2*MAX_L+1)*(2*MAX_L+1);k++) result[k]=0;

  if (l<1 || i<0 || j<0) {printf("invalid choice in differential operator: %i, %i, %i\n",i,j,l);exit(1);}
  if (l>4) {printf("precision %i not implemented in differential operator\n",l);exit(1);}
  if (i+j>3) {printf("derivative %i, %i not implemented in differential operator\n",i,j);exit(1);}

  bool swap=false;
  if (j>i) {int k; k=i; i=j; j=k; swap=true;}
  
  if (i==0 && j==0) result[l*(2*l+1)+l]=1;

  if (i==1 && j==0) {
     double vektor[MAX_L][2*MAX_L+1]={\
  			{-1,0,1},\
			{-1,8,0,-8,1},\
			{-1,9,-45,0,45,-9,1},\
			{-1,32.0/3,-56,224,0,-224,56,-32.0/3,1}};
     double koefs[MAX_L]={1.0,-1.0/12,1.0/60,-1.0/280};
     for (int x=0;x<2*l+1;x++) {
       result[x*(2*l+1)+l]=vektor[l-1][x]*koefs[l-1];
     }
  }

  if (i==2 && j==0) {
     double vektor[MAX_L][2*MAX_L+1]={\
			{1,-2,1},\
			{1,-16,30,-16,1},\
			{2,-27,270,-245*2,270,-27,2},\
			{-9,128,-1008,8064,-7175*2,8064,-1008,128,-9}};
     double koefs[MAX_L]={1.0,-1.0/12,1.0/180,1.0/(9*560)};
     for (int x=0;x<2*l+1;x++) {
       result[x*(2*l+1)+l]=vektor[l-1][x]*koefs[l-1];
     }
  }

  if (i==3 && j==0) {
     double vektor[MAX_L][2*MAX_L+1]={\
			{0,0,0},\
			{-0.5,1,0,-1,0.5},\
			{1,-8,13,0,-13,8,-1},\
			{-7,72,-338,488,0,-488,338,-72,7}};

     double koefs[MAX_L]={1.0,1.0,1.0/8,1.0/240};
     for (int x=0;x<2*l+1;x++) {
       if (swap) result[l*(2*l+1)+x]=vektor[l-1][x]*koefs[l-1];
       else result[x*(2*l+1)+l]=vektor[l-1][x]*koefs[l-1];
     }
  }
  
  //in the next matrix is written such that index y goes backwards
  //what is accounted for as the result is written
  if (i==1 && j==1) {
    double matrix[MAX_L][(2*MAX_L+1)*(2*MAX_L+1)]={\
                          {-1,   0,   1,\
	                    0,   0,   0,\
	                    1,   0,  -1},\
	\
	             {-1,   0,   0,   0,   1,\
	               0,  16,   0, -16,   0,\
		       0,   0,   0,   0,   0,\
		       0, -16,   0,  16,   0,\
		       1,   0,   0,   0,  -1},\
	 \
	        {-2,   0,   0,   0,   0,   0,   2,\
	          0,  27,   0,   0,   0, -27,   0,\
		  0,   0,-270,   0, 270,   0,   0,\
		  0,   0,   0,   0,   0,   0,   0,\
		  0,   0, 270,   0,-270,   0,   0,\
		  0, -27,   0,   0,   0,  27,   0,\
		  2,   0,   0,   0,   0,   0, -2},\
		  \
     {-1/112.0,   0,   0,   0,   0,   0,   0,   0,1/112.0,\
	     0,8.0/63, 0,   0,   0,   0,   0,-8.0/63,0,\
	     0,   0,  -1,   0,   0,   0,   1,   0,   0,\
	     0,   0,   0,   8,   0,  -8,   0,   0,   0,\
	     0,   0,   0,   0,   0,   0,   0,   0,   0,\
	     0,   0,   0,  -8,   0,   8,   0,   0,   0,\
	     0,   0,   1,   0,   0,   0,  -1,   0,   0,\
	     0,-8.0/63,0,   0,   0,   0,   0,8.0/63, 0,\
       1/112.0,   0,   0,   0,   0,   0,   0,   0,-1/112.0}};

    double koefs[MAX_L]={1.0,-1.0/48,1.0/720,-1.0/20};

     for (int x=0;x<2*l+1;x++) {
       for (int y=0;y<2*l+1;y++) {
         result[x*(2*l+1)+l]=matrix[l-1][y*(2*l+1)+(2*l-x)]*koefs[l-1];
       }
     }
  }
  if (i==2 && j==1) {
    double matrix[MAX_L][(2*MAX_L+1)*(2*MAX_L+1)]={\
                          { 1,  -2,   1,\
	                    0,   0,   0,\
	                   -1,   2,  -1},\
	\
	              {1,   0,  -2,   0,   1,\
	               0, -32,  64, -32,   0,\
		       0,   0,   0,   0,   0,\
		       0,  32, -64,  32,   0,\
		      -1,   0,   2,   0,  -1},\
	 \
	         {4,   0,   0,  -8,   0,   0,   4,\
	          0, -81,   0, 162,   0, -81,  0,\
		  0,   0,1620,-3240,1620, 0,   0,\
		  0,   0,   0,   0,   0,   0,   0,\
		  0,   0,-1620,3240,-1620, 0,   0,\
		  0,  81,   0,-162,   0,  81,   0,\
		 -4,   0,   0,   8,   0,   0,  -4},\
		  \
     {1.0/224,   0,   0,   0,-2.0/224,0,   0,   0,1.0/224,\
	     0,-16.0/189,0, 0,32.0/189,0, 0,-16.0/189,0,\
	     0,   0,   1,   0,  -2,   0,   1,   0,   0,\
	     0,   0,   0, -16,  32, -16,   0,   0,   0,\
	     0,   0,   0,   0,   0,   0,   0,   0,   0,\
	     0,   0,   0,  16, -32,  16,   0,   0,   0,\
	     0,   0,  -1,   0,   2,   0,  -1,   0,   0,\
	     0,16.0/189,0,  0,-32.0/189,0,  0,16.0/189,0,\
      -1.0/224,   0,   0,   0,2.0/224,0,   0,   0,-1.0/224}};

    double koefs[MAX_L]={1.0,-1.0/48,1.0/2160,-1.0/20};

     for (int x=0;x<2*l+1;x++) {
       for (int y=0;y<2*l+1;y++) {
         result[x*(2*l+1)+y]=matrix[l-1][y*(2*l+1)+(2*l-x)]*koefs[l-1];
       }
     }
  }
  
  //if swap - swap along +x +y diagonal
     for (int x=0;x<2*l+1;x++) {
       for (int y=0;y<2*l+1;y++) {
         if (swap) out[x*(2*l+1)+y]=result[y*(2*l+1)+x]*coef;
         else out[x*(2*l+1)+y]=result[x*(2*l+1)+y]*coef;
       }
     }
 /*    if (swap) printf("check: matrix %i, %i is\n",j,i);
     else printf("check: matrix %i, %i is\n",i,j);

     for (int j=2*l;j>-1;j--) {
       for (int i=0;i<2*l+1;i++) {
         printf("%f\t",out[i*(2*l+1)+j]);
       }
       printf("\n");
     }
     printf("\n\n"); */
}

//returns coefficients defining the differential operator d^d/dx^d at the zero coordinate
//input: d (previous - the derivative), n = number of available points (=2b+1 previously)
//vector of coordinates of the available points (in internal length units)
void construct_differential_operator(int d, int n, double* coors, double* coefs)
{
  //construct the matrix
  gsl_matrix *m = gsl_matrix_alloc(n,n);
  for (int i=0;i<n;i++) {//Taylor expansion terms (derivatives)
    for (int j=0;j<n;j++) {//point = column
      if (i==0) gsl_matrix_set(m, i, j, 1);
      else {
        double aux=coors[j];
        for (int k=1;k<i;k++) aux*=coors[j]/(k+1);
        //fprintf(logg,"at i=%i, j=%i, setting %e\n",i,j,aux);
        gsl_matrix_set(m, i, j, aux);
      }
    }
  }
  double* m_aux=new double[n*n];
  for (int i=0;i<n;i++) for (int j=0;j<n;j++) m_aux[i*n+j]=gsl_matrix_get(m,i,j);

  //construct the RHS
  gsl_vector *b = gsl_vector_calloc (n);
  gsl_vector_set(b,d,1);
  //solve the set of equations
  gsl_permutation * p = gsl_permutation_alloc (n);
  gsl_vector *x=gsl_vector_alloc(n);
  int s;
  int error=gsl_linalg_LU_decomp (m, p, &s);
  if (error) {printf("gsl_linalg_LU_decomp returned with error=%i in construct_differential_operator\n",error);exit(1);}
  error=gsl_linalg_LU_solve (m, p, b, x);
  if (error) {printf("gsl_linalg_LU_solve returned with error=%i in construct_differential_operator\n",error);exit(1);}
  //copy results
  for (int i=0;i<n;i++) coefs[i]=gsl_vector_get(x,i);
#if AUXIL > 1
  fprintf(logg,"differential operator d^%i/dx^%i in %i-point scheme:\ncoordinates:\t",d,d,n);
  for (int i=0;i<n;i++) fprintf(logg,"%i %+.2e ",i, coors[i]);
  fprintf(logg,"\ncoefficients:\t");
  for (int i=0;i<n;i++) fprintf(logg,"%i %+.2e ",i, coefs[i]);
  fprintf(logg,"\n");
  double co=abs(coors[1]-coors[0]); 
  fprintf(logg,"coordinates/coefficients rescaled by %+.3e/%+.3e:\t",co,1.0/co);
  for (int i=0;i<n;i++) fprintf(logg,"%i %+.2e %+.2e ",i, coors[i]/co, coefs[i]*co);
  fprintf(logg,"\n");
  double sum=0;
  for (int i=0;i<n;i++) {
    double res=0;
    for (int j=0;j<n;j++) {
      res+=m_aux[i*n+j]*coefs[j];
      //fprintf(logg,"m[%i, %i]=%+.3e x[%i] =%+.3e ",i,j,m_aux[i*n+j],j,coefs[j]);
    }
    //fprintf(logg,"\nline %i: M x =%e vs vector b=%e\n",i,res,gsl_vector_get(b,i));
    sum+=abs(res-gsl_vector_get(b,i));
  }
  if (sum>1e-12) fprintf(logg,"residual error of the linear set of equations is %e in construct_differential_operator for d=%i, n=%i\n", sum, d, n);
#endif
  //fprintf(logg,"construct differential operator (for d=%i, n=%i) returns:\n",d,n);
  //for (int i=0;i<n;i++) fprintf(logg,"coors[%i]=%+e->coefs[%i]=%e\n",i,coors[i], i,coefs[i]);
  
  //release memory
  gsl_permutation_free (p);
  gsl_vector_free (x);
  gsl_vector_free(b);
  gsl_matrix_free(m);
  delete[] m_aux;
}


//!**************************************************

 void stopky::pusti()
{
  prvok *novy;
  novy=(prvok*) malloc(sizeof(prvok));
  if (bezi) posledny->dalsi=novy;
  else prvy=novy;
  novy->dalsi=0;
  posledny=novy;
  ftime(&(posledny->cas));
  bezi++;
}

long int stopky::doba(int ktora)
{
  if ((ktora<0) || (ktora>=bezi)) return(0);
  prvok ktory;
  if (ktora==0) ktory.dalsi=posledny;
  else {
    ktory.dalsi=prvy;
    for (int i=1;i<ktora;i++) ktory.dalsi=ktory.dalsi->dalsi;
  }
  ftime(&(ktory.cas));
  long int x;
  x=(ktory.cas.time-ktory.dalsi->cas.time)*1000+ktory.cas.millitm-ktory.dalsi->cas.millitm;
  return(x);
}

long int stopky::zastav(int ktora)
{
  if ((ktora<0) || (ktora>bezi) || (bezi==0)) return(0);
  prvok predch, nasl, tento; //vsetky tri budu obsahovat v ukazovateli dalsi adresu prvku svojho mena, nie dalsieho po nom!!!
  if (ktora==0) ktora=bezi;
  if (ktora==1) {
    predch.dalsi=0;
    tento.dalsi=prvy;
  }
  else {
    predch.dalsi=prvy;
    for (int i=1;i<ktora-1;i++) predch.dalsi=predch.dalsi->dalsi;
    tento.dalsi=predch.dalsi->dalsi;
  }
  if (ktora==bezi) nasl.dalsi=0;
  else nasl.dalsi=tento.dalsi->dalsi;	//oba pripady by mali dat to iste
  ftime(&(tento.cas));
  long int x;
  x=(tento.cas.time-tento.dalsi->cas.time)*1000+tento.cas.millitm-tento.dalsi->cas.millitm;
  if (ktora==1) prvy=nasl.dalsi;
  else predch.dalsi->dalsi=nasl.dalsi;
  if (ktora==bezi) posledny=predch.dalsi;
  free(tento.dalsi);
  bezi--;
  return(x);
}


//!***************************************

 int u_getch(void)
    {
	static signed int returned = (-1), fd;
	static struct termios nueva, antes;
	fd = fileno(stdin);
	tcgetattr(fd, &antes);
	nueva = antes;
	nueva.c_lflag &= ~(ICANON|ECHO);
	tcsetattr(fd, TCSANOW, &nueva);
	returned = getchar();
	tcsetattr(fd, TCSANOW, &antes);
	return returned;
    }

int stop()
{
  int znak=u_getch();
  if (znak==27) exit(1);
  return(znak);
}

 double generuj()  /*nahodne cislo z intervalu <0,1> */
{
  double x=0,y;
  for (int i=0;i<3;i++)
  {
    y=rand();
    x=x/(RAND_MAX+0.0)+y;
  }

  return(x/(RAND_MAX+0.0));
}

//naplni generator aktualnou milisekundou a vrati ju
int randomize(int seed)
{
  timeb cas;
  ftime(&cas);
  if (seed==-1) seed=cas.millitm;
  srand(seed);
  return(seed);
}


//linear kongruential without shuffling
double ran1(int i)
{
int a=16807;
int m=2147483647;		//2^31-1
int p=2836;
int q=127773;
double rm=(1.0/m);
double EPS=1.2e-7;
double RNMX=(1-EPS);

  static int in,out;
  double output;

  if (i!=0) in=i;		//seed
  else if (in==0) in=12345;

  int k=in/q;

  out=a*(in-k*q)-p*k;	//out=in*a mod m without overflows
  if (out<0) out+=m;

  in=out;
  if ((output=rm*out)>RNMX) return (RNMX);
  else return(output);		//return value from <0 to 1)
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!                         progress bar
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

progress_bar_prototype::progress_bar_prototype(const char* title_)
{
  //file_desc_stdin = fileno(stdin);
  //tcgetattr(file_desc_stdin, &preserve);
  //temp=preserve;
  //temp.c_lflag &= ~(ICANON | ECHO | ECHOE | ECHOK | ECHONL | ECHOPRT | ECHOKE | ICRNL | ISIG);
  //temp.c_cc[VMIN]=1;
  //temp.c_cc[VTIME] = 0;
  //tcsetattr(file_desc_stdin, TCSANOW, &temp);

  length=0;
  todo=0;
  done=0;
  last_update=-21;
  title=new char[100];
  strcpy(title, title_);
  timer=new stopky();
  aux=new char[100];
  aux2=new char[100];
  mess=new char[100];
}

progress_bar_prototype::~progress_bar_prototype()
{
  
  delete title;
  delete aux;
  delete aux2;
  delete mess;
  
  delete timer;
  //tcsetattr(file_desc_stdin, TCSANOW, &preserve);

};

void progress_bar_prototype::reset(const char* title_)
{
  print_message(false);
  strcpy(title,title_);
  timer->zastav();
  last_update=-21;
  todo=0;
  done=0;
}

void progress_bar_prototype::start(bool on)
{ 
  timer->pusti();
  print_message(on);
}

void progress_bar_prototype::print_message(bool on)
{
  long int now=timer->doba();
  //printf("pg called: now=%li, last_update=%li asked to print: %i\n",now,last_update,on);
  if ((now-last_update)<100 && on) return;

  //update eta
  if (done==0 || todo<done) eta=0;
  else eta=now*(todo-done)/done;

  //delete actual message
  //tcsetattr(file_desc_stdin, TCSANOW, &temp);
  //for (int i=0;i<length;i++) putchar('\b');//backspace
  for (int i=0;i<length;i++) fputc('\b',stderr);//backspace
  //tcsetattr(file_desc_stdin, TCSANOW, &preserve);
  length=0;
  
  //construct and output message
  if (!on) return;
  strcpy(mess, title);
  sprintf(aux," %i/%i ", (int) my_round(done), (int) my_round(max(todo,done)));
  strcat(mess,aux);
  int l=strlen(mess);
  for (int i=0;i<30-l;i++) strcat(mess," ");
  strcat(mess, "  ETA ");
  time2string((long int) round(eta), aux);
  strcat(mess,aux);
  length=strlen(mess);
  //tcsetattr(file_desc_stdin, TCSANOW, &temp);
  //for (int i=0;i<length;i++) putchar(mess[i]);
  fprintf(stderr,"%s",mess);
  //tcsetattr(file_desc_stdin, TCSANOW, &preserve);
  last_update=now;
}	

void progress_bar_prototype::add(double todo_, double done_, bool print)
{
  todo+=todo_;
  done+=done_;
  if (print) print_message(print);
}

long int progress_bar_prototype::finished(bool on)
{

  //delete actual message
  //tcsetattr(file_desc_stdin, TCSANOW, &temp);
  //for (int i=0;i<length;i++) putchar('\b');//backspace
  for (int i=0;i<length;i++) fputc('\b',stderr);//backspace
  //tcsetattr(file_desc_stdin, TCSANOW, &preserve);
  length=0;
  
  if (!on) return(timer->zastav());

  //construct and output message
  strcpy(mess, title);
  sprintf(aux," %i/%i ", (int) my_round(done), (int) my_round(max(todo,done)));
  strcat(mess,aux);
  int l=strlen(mess);
  for (int i=0;i<30-l;i++) strcat(mess," ");
  strcat(mess, " took ");
  time2string(timer->doba(), aux);
  strcat(mess,aux);
  length=strlen(mess);
  //for (int i=0;i<length;i++) putchar(mess[i]);
  //putchar('\n');
  fprintf(stderr,"%s\n",mess);
  length=0;
  return(timer->zastav());
}

void progress_bar_prototype::time2string(long int time, char* res)
{ 
  int hrs=time / (60*60*1000);
  time-=hrs*60*60*1000;
  int min=time / (60*1000);
  time-=min*60*1000;
  int sec=time/1000;
  time-=sec*1000;

  res[0]='\0';
  if (hrs>0) {
    sprintf(aux2, "%ih:", hrs);
    strcat(res, aux2);
  }
  sprintf(aux2, "%02i:%02i.%03li s               ", min,sec,time);
  strcat(res, aux2);

}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!                         Numerical Recipes Routines
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//#define NR1
//#define NR2

//#define NRCHECK

complex<double> polint(double* x, complex<double>* y,  int N, double x0, double& err)
{
  static complex<double>* c=0;
  static complex<double>* d=0;
  static int N_old=0;
  if (N<1) {printf("wrong order (%i) in polint\n",N);exit(1);}

  if (N!=N_old) {
    if (c!=0) delete c;
    if (d!=0) delete d;
    c=new complex<double>[N];
    d=new complex<double>[N];
    N_old=N;
  }
  
  int ns=0;
  double difmin=abs(x0-x[0]);

  //initialize and find the closest point
  for (int i=0;i<N;i++) {
    double dif=abs(x0-x[i]);
    if (dif<difmin) {
      difmin=dif;
      ns=i;
    }
    c[i]=d[i]=y[i];
#ifdef NR2
    fprintf(logg,"polint: source data (%i) y[%e]=%e%+ei\n",i, x[i], y[i].real(), y[i].imag());
#endif
  }
#ifdef NR2
  fprintf(logg,"polint: interpolating at x0=%e\n",x0);
#endif

  complex<double> y0=y[ns];
  complex<double> dy;

  //iterate recursively
  for (int m=1;m<N;m++) {//m is the actual order (table column)
    for (int i=0;i<N-m;i++) {//i is the table row
      double ho=x[i]-x0;
      double hp=x[i+m]-x0;
      complex<double> w=c[i+1]-d[i];
      double den=ho-hp;
      if (den==0) {printf("x values (%i and %i) the same in polint\n",i,i+m);exit(1);}
      w/=den;
      d[i]=hp*w;
      c[i]=ho*w;
    }
    if (2*ns<N-m) dy=c[ns];
    else {dy=d[ns-1];ns--;}
    y0+=dy;
  }
  err=abs(dy);
#ifdef NR2
  fprintf(logg,"resulting in: %e%+ei with error estimate %e\n",y0.real(), y0.imag(),err);
#endif
  return(y0);
}

void polint_check(int N)
{
  complex<double>* cn=new complex<double>[N];
  complex<double>* y=new complex<double>[N+1];
  double* x=new double[N+1];

  for (int i=0;i<N;i++) cn[i]=complex<double> (generuj(), generuj());
  for (int i=0;i<N+1;i++) {
    if (i>0) x[i]=x[i-1]-generuj()/10; else x[i]=generuj()/10;
    y[i]=0;
    for (int j=0;j<N;j++) y[i]+=cn[j]*pow(x[i],j);
  }

  double x0=(x[N+1]-x[0])*generuj()*0+x[0]+generuj();
  complex<double> res=0, respol;
  double err;
  for (int j=0;j<N;j++) res+=cn[j]*pow(x0,j);
  respol=polint(x,y,N-1,x0,err);
  printf("%i order polint (of polynome of %i order): error estimate %e vs exact %e\n",N-1,N,err,abs(respol-res));
  respol=polint(x,y,N,x0,err);
  printf("%i order polint (of polynome of %i order): error estimate %e vs exact %e\n",N,N,err,abs(respol-res));
  respol=polint(x,y,N+1,x0,err);
  printf("%i order polint (of polynome of %i order): error estimate %e vs exact %e\n",N+1,N,err,abs(respol-res));

  delete[] cn;
  delete[] y;
  delete[] x;

}


double trapzd( complex<double> f(double) , double a, double b, complex<double>& s, int n)
{
  if (n==1) {
    s=(f(a)+f(b))*(b-a)/2.0;
    return(0);
  }
  else {
    int it=1;
    for (int i=1;i<n-1;i++) it*=2;
    double del=(b-a)/it;
    double x=a+del/2;
    complex<double> sum=0;
    for (int j=0;j<it;j++) {
      sum+=f(x);
      x+=del;
    }
    sum*=del;
    s=(s+sum)/2.0;
    return(abs(s-sum));
  }
}

double midpnt( complex<double> f(double) , double a, double b, complex<double>& s, int n)
{
  if (n==1) {
    s=(b-a)*f((a+b)/2);
    return(0);
  }
  else {
    int it=1;
    for (int i=1;i<n-1;i++) it*=3;
    double del=(b-a)/(3*it);
    double ddel=2*del;
    double x=a+del/2;
    complex<double> sum=0;
    for (int j=0;j<it;j++) {
      sum+=f(x);
      x+=ddel;
      sum+=f(x);
      x+=del;
    }
    sum*=del;
    s=s/3.0+sum;
    return(abs(2.0*s-3.0*sum));
  }
}


void trapzd_check(complex<double> f(double), complex<double> F(double))
{
  double a=0, b=0.25;
  complex<double> s1,s2;
  complex<double> exact=F(b)-F(a);
  for (int n=1;n<11;n++) {
    double err=trapzd(f,a,b,s1,n);
    printf("trapzd check: step %i: error estimate %e vs exact %e (int:%e)\n", n, err, abs(s1-exact), abs(exact));
    err=midpnt(f,a,b,s2,n);
    printf("midpnt check: step %i: error estimate %e vs exact %e (int:%e)\n", n, err, abs(s2-exact), abs(exact));
  }
}

complex<double> qromb(complex<double> f(double), double a, double b, int& maxsteps, int order, double& err)
{
  if (maxsteps<order) {printf("not enough steps allowed (%i) for polint of order %i in qromb\n",maxsteps,order);exit(1);}
  complex<double>* s=new complex<double>[maxsteps+1];
  double* h=new double[maxsteps+1];
#ifdef RELAXCHECK
  if (h==0 || s==0) {printf("memory alloc failed in qromb (maybe maxsteps too big? (%i)\n",maxsteps);exit(1);}
#endif
  h[0]=1;
  complex<double> est;
  double error;
  for (int j=0;j<maxsteps;j++) {
    trapzd(f, a, b, s[j],j+1);
    if (j+1>=order) {
      est=polint(h+j+1-order, s+j+1-order, order, 0.0, error);
#ifdef NR1
      fprintf(logg,"qromb: interpolation at step %i: error %e vs estimate %e and err %e\n",j+1,error, abs(est), err);
#endif
      if (error<=err*abs(est)) {
        err=error;
        maxsteps=j+1;
        delete[] s;
        delete[] h;
        return(est);
      }
    }
    s[j+1]=s[j];
    h[j+1]=h[j]/4;
  }
  err=error;
  delete[] s;
  delete[] h;
  return(est);
}

complex<double> qromo(complex<double> f(double), double a, double b, int& maxsteps, int order, double& err)
{
  if (maxsteps<order) {printf("not enough steps allowed (%i) for polint of order %i in qromo\n",maxsteps,order);exit(1);}
  complex<double>* s=new complex<double>[maxsteps+1];
  double* h=new double[maxsteps+1];
  h[0]=1;
  complex<double> est;
  double error;
  for (int j=0;j<maxsteps;j++) {
    midpnt(f, a, b, s[j],j+1);
    if (j+1>=order) {
      est=polint(h+j+1-order, s+j+1-order, order, 0.0, error);
#ifdef NR1
      fprintf(logg,"qromo: interpolation at step %i: error %e vs estimate %e and err %e\n",j+1,error, abs(est), err);
#endif
      if (error<=err*abs(est)) {
        err=error;
        maxsteps=j+1;
        delete[] s;
        delete[] h;
        return(est);
      }
    }
    s[j+1]=s[j];
    h[j+1]=h[j]/9;
  }
  err=error;
  delete[] s;
  delete[] h;
  return(est);
}

void qromb_check(complex<double> f(double), complex<double> F(double))
{
  double a=0, b=0.25;
  complex<double> exact=F(b)-F(a);
  int order=4;
  for (int i=order;i<12;i++) {
    int maxsteps=i;
    double error=0;
    complex<double> est=qromb(f,a,b,maxsteps,order,error);
    printf("qromb check: after %i steps estimated error %e vs exact error %e\n", maxsteps, error,abs(est-exact));
    est=qromo(f,a,b,maxsteps,order,error);
    printf("qromo check: after %i steps estimated error %e vs exact error %e\n", maxsteps, error,abs(est-exact));
  }
}

complex<double> auxf(double x) {return(complex<double> (sin(2*M_PI*x),cos(2*M_PI*x)));}
complex<double> auxF(double x) {return(complex<double> (-cos(2*M_PI*x)/(2*M_PI),sin(2*M_PI*x)/(2*M_PI)));}
complex<double> auxf2(double x) {return(complex<double> (x*exp(-x*x),-x*exp(-x*x)));}
complex<double> auxF2(double x) {return(complex<double> (-exp(-x*x)/2,exp(-x*x)/2));}

void nr_check1() 
{
  for (int i=2;i<10;i++) polint_check(i);
  trapzd_check(auxf, auxF);
  qromb_check(auxf, auxF);
  trapzd_check(auxf2, auxF2);
  qromb_check(auxf2, auxF2);
}


/*
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!                         banded matrix prototype
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void banded_matrix_prototype::cart2diag(int i, int j, int &d, int &seq)
{
  if (i<0 || i>=N || j<0 || j>=N) 
    {printf("banded_matrix_prototype::cart2diag wrong input: i=%i, j=%i)\n",i,j); exit(1);} 
  d=i-j;
  if (d<0) d=-1-2*d; else d*=2;
  seq=j;
  if (d<0 || d>=2*N-1 || seq<0 || seq>=N) 
    {printf("banded_matrix_prototype::cart2diag wrong output: d=%i, seq=%i (from i=%i, j=%i))\n",d,seq,i,j); exit(1);} 
}

void banded_matrix_prototype::diag2cart(int d, int seq, int &i, int &j)
{
  if (d<0 || d>=2*N-1 || seq<0 || seq>=N) 
    {printf("banded_matrix_prototype::diag2cart wrong input: d=%i, seq=%i)\n",d,seq); exit(1);} 
  j=seq;
  if (d%2==1) d=(-1-d)/2; else d/=2;
  i=d+j;
  if (i<0 || i>=N || j<0 || j>=N) 
    {printf("banded_matrix_prototype::diag2cart wrong output: i=%i, j=%i from d=%i, seq=%i))\n",i,j,d, seq); exit(1);} 
}

void banded_matrix_prototype::add(int i, int j, complex<double> x)
{
  if (x.real()==0 && x.imag()==0) return;
  int d, seq;
  cart2diag(i,j,d,seq);
  complex<double>* & p=diagonal[d];
  if (p==0) {
    p=new complex<double> [N];
    for (int k=0;k<N;k++) p[k]=0;
  }
  if (p[seq].real()==0 && p[seq].imag()==0) pntr[d]++;
  p[seq]+=x;
}

complex<double> banded_matrix_prototype::read(int i, int j)
{
  int d, seq;
  cart2diag(i,j,d,seq);
  complex<double>* & p=diagonal[d];
  if (p==0) return(0);
  return(p[seq]);
}

  
banded_matrix_prototype::banded_matrix_prototype(int N)
{
  this->N=N;
  if (N<0) {printf("banded_matrix_prototype constructor wrong input: N=%i)\n",N); exit(1);}
  pntr=new int[2*N];
  diagonal=new complex<double>* [2*N];
  for (int i=0;i<2*N-1;i++) {
    pntr[i]=0;
    diagonal[i]=0;
  }
  AB=new complex <double>[0];
  IPIV=new int[N];
}

banded_matrix_prototype::~banded_matrix_prototype()
{
  delete AB;
  delete IPIV;
  delete pntr;
  for (int i=0;i<2*N-1;i++) if (diagonal[i]!=0) delete diagonal[i];
  delete diagonal; 
}


void banded_matrix_prototype::matrix_info(bool out)
{
  Nl=Nu=KL=KU=0;
  if (out) fprintf(logg,"nonzero lower diagonals: (#diagonal offset x number of nonzero entries)\n");
  for (int k=1;k<N;k++) if (diagonal[2*k]!=0) {
    Nl++;
    KL=k;
    if (out) fprintf(logg,"#%3ix%3i ",k, pntr[2*k]);
  }
  if (out) fprintf(logg,"\nnonzero upper diagonals: (#diagonal offset x number of nonzero entries)\n");

  for (int k=1;k<N;k++) if (diagonal[2*k-1]!=0) {
    Nu++;
    KU=k;
    if (out) fprintf(logg,"#%ix%i ",k, pntr[2*k-1]);
  }
  if (out) fprintf(logg,"\nwith N=%i, number of upper diagonals=%i (last one is %i), number of lower diagonals=%i (last one is %i)\n", N, Nu, KU, Nl, KL);
  if (KU>=N || KL>=N) {printf("wrong number of diagonals in matrix_info: N=%i, KL=%i, KU=%i\n",N,KL,KU); exit(1);}
}

void banded_matrix_prototype::matrix2lapack()
{
  matrix_info(false);
  fprintf(logg,"converting banded matrix into lapack form: N=%i, KU=%i, KL=%i\n",N, KU, KL);
  LDAB=2*KL+KU+1;
  delete AB;
  AB=new complex <double>[N*LDAB];
  if (AB==0) {printf("failed to allocate %e MB of memory in matrix2lapack\n", N*LDAB*sizeof(complex<double>)/1e+6 );exit(1);}
  for (int i=0;i<N*LDAB;i++) AB[i]=0;
  for (int k=1;k<=KU;k++) {
    int i=KL+KU-k;
    complex<double>* & p=diagonal[2*k-1];
    if (p==0) continue;
    for (int j=0;j<N;j++) AB[j*LDAB+i]=p[j];//in AB the transpose is for fortran
  }
  for (int k=0;k<=KL;k++) {//includes the main diagonal being k=0
    int i=KL+KU+k;
    complex<double>* & p=diagonal[2*k];
    if (p==0) continue;
    for (int j=0;j<N;j++) AB[j*LDAB+i]=p[j];//in AB the transpose is for fortran
  }
  //alternative way:
  //for (int i=0;i<N*LDAB;i++) AB[i]=0;
  //for (int i=0; i<N;i++) for (int j=0;j<N;j++) {
  //  if (i-j>=-KU && i-j<=KL) AB[KL+KU+i-j+LDAB*j]=read(i,j);
  //} 
}


void banded_matrix_prototype::LU_factorize()
{
  int INFO;
  F77NAME(zgbtrf)(&N,&N,&KL,&KU,AB,&LDAB,IPIV,&INFO);
  if (INFO!=0) {printf("banded matrix LU factorization failed: zgbtrf INFO=%i\n",INFO);exit(1);}  
}

//solve M x = y, do not overwrite y
void banded_matrix_prototype::LU_solve(complex<double>* x, complex<double> *y)
{
  for (int i=0;i<N;i++) x[i]=y[i];
  char TRANS='N';
  int NRHS=1, INFO=0;
  F77NAME(zgbtrs)(&TRANS, &N, &KL,&KU, &NRHS, AB, &LDAB, IPIV, x, &N, &INFO);
  if (INFO!=0) {printf("banded LU-factorized operator could not be inverted: zgbtrs INFO=%i\n",INFO);exit(1);}
}*/

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!                         sparse matrix prototype
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

extern "C"
{
  void F77NAME(zgbtrf)( int* , int* , int* , int* , complex<double>* , int* , int* , int*);
  void F77NAME(zgbtrs)( char*, int* , int* , int*, int*, complex<double>* , int* , int* , complex<double>*, int*, int*);
}

sparse_matrix_prototype::sparse_matrix_prototype(int N, int max_length)
{
  this->max_length=max_length;
  this->N=N;
  if (N<0 || max_length<0 || max_length>N) {printf("sparse_matrix_prototype constructor wrong input: N=%i, max_length=%i)\n",N,max_length); exit(1);}
  length=new int[N];
  label=new int*[N];
  value=new complex<double>*[N];
  diagonal=new int[2*N];
  for (int i=0;i<N;i++) {
    diagonal[i]=diagonal[N+i]=length[i]=0;
    label[i]=new int[max_length];
    value[i]=new complex<double>[max_length];
  }
  AB=new complex <double>[0];
  IPIV=new int[N];  
}

sparse_matrix_prototype::~sparse_matrix_prototype()
{
  for (int i=0;i<N;i++) {
    length[i]=0;
    delete label[i];
    delete value[i];
  }
  delete value;
  delete label;
  delete length;
  delete AB;
  delete IPIV;
  delete diagonal;
}

void sparse_matrix_prototype::add(int i, int j, complex<double> x)
{
  //first check whether this entry is already non zero
  for (int k=0;k<length[i];k++) if (label[i][k]==j) {value[i][k]+=x; return;}
  //increase number of entries, 
  if (++length[i]>=max_length) {//if we reached limit, enlarge the array dimension
    complex<double>* value_old=value[i];
    int* label_old=label[i];
    value[i]=new complex<double>[length[i]];
    label[i]=new int[length[i]];
    for (int k=0;k<length[i]-1;k++) {
      value[i][k]=value_old[k];
      label[i][k]=label_old[k];
    }
    delete value_old;
    delete label_old;
  }
  //put the new entry
  value[i][length[i]-1]=x;
  label[i][length[i]-1]=j;
}

complex<double> sparse_matrix_prototype::read(int i, int j)
{
  for (int k=0;k<length[i];k++) if (label[i][k]==j) return(value[i][k]);
  return(0);
}

void sparse_matrix_prototype::operate(complex<double>*x, complex<double>*y) //y = M x
{
  for (int i=0;i<N;i++) {
    y[i]=0;
    for (int k=0;k<length[i];k++) y[i]+=x[label[i][k]]*value[i][k];
  }
}

void sparse_matrix_prototype::matrix_structure(bool out)
{
  //first find out the farthest occupied diagonals
  int ave=0, nonzero=0, maximal=0, minimal=N;
  for (int i=0;i<N;i++) {
    int l=length[i];
    if (minimal>l) minimal=l;
    if (l==0) continue;
    ave+=l;
    if (l>maximal) maximal=l;
    for (int j=0;j<l;j++) {
      int d, seq;
      cart2diag(i,label[i][j],d,seq);
      diagonal[d]++;
    }
  }
  if (out) fprintf(logg,"number of line entries is from %i to %i with average %e\n", minimal, maximal, ((double) ave)/N);
  Nl=Nu=KL=KU=0;
  if (out) fprintf(logg,"nonzero lower diagonals: (#diagonal offset x number of nonzero entries)\n");
  for (int k=1;k<N;k++) if (diagonal[2*k]>0) {
    Nl++;
    KL=k;
    if (out) fprintf(logg,"#%ix%i ", k, diagonal[2*k]);
  }
  if (out) fprintf(logg,"\nnonzero upper diagonals: (#diagonal offset x number of nonzero entries)\n");

  for (int k=1;k<N;k++) if (diagonal[2*k-1]>0) {
    Nu++;
    KU=k;
    if (out) fprintf(logg,"#%ix%i ",k, diagonal[2*k-1]);
  }
  if (out) fprintf(logg,"\nwith N=%i, number of upper diagonals=%i (last one is %i), number of lower diagonals=%i (last one is %i)\n", N, Nu, KU, Nl, KL);
  if (KU>=N || KL>=N) {printf("wrong number of diagonals in matrix_info: N=%i, KL=%i, KU=%i\n",N,KL,KU); exit(1);}
}
  
void sparse_matrix_prototype::matrix2banded()
{
  fprintf(logg,"converting sparse matrix into banded lapack form: N=%i, KU=%i, KL=%i\n",N, KU, KL);
  LDAB=2*KL+KU+1;
  delete AB;
  AB=new complex <double>[N*LDAB];
  if (AB==0) {printf("failed to allocate %e MB of memory in matrix2lapack\n", N*LDAB*sizeof(complex<double>)/1e+6 );exit(1);}
  for (int i=0;i<N*LDAB;i++) AB[i]=0;

  complex<double>* AB_offset=AB+KL+KU;
  int LDAB1=LDAB-1;
  for (int i=0;i<N;i++) {
    int l=length[i];
    for (int k=0;k<l;k++) {
      int j=label[i][k];
      AB_offset[i+LDAB1*j]=value[i][k];
    }
  }
  //the formula is:
  /*
  for (int i=0; i<N;i++) for (int j=0;j<N;j++) {
    if (i-j>=-KU && i-j<=KL) AB[KL+KU+i-j+LDAB*j]=read(i,j);
  } */ 
}
void sparse_matrix_prototype::LU_banded_factorize()
{
  int INFO;
  F77NAME(zgbtrf)(&N,&N,&KL,&KU,AB,&LDAB,IPIV,&INFO);
  if (INFO!=0) {printf("banded matrix LU factorization failed: zgbtrf INFO=%i\n",INFO);exit(1);}  
}

//solve M x = y, do not overwrite y
void sparse_matrix_prototype::LU_banded_solve(complex<double>* x, complex<double> *y)
{
  for (int i=0;i<N;i++) x[i]=y[i];
  char TRANS='N';
  int NRHS=1, INFO=0;
  F77NAME(zgbtrs)(&TRANS, &N, &KL,&KU, &NRHS, AB, &LDAB, IPIV, x, &N, &INFO);
  if (INFO!=0) {printf("banded LU-factorized operator could not be inverted: zgbtrs INFO=%i\n",INFO);exit(1);}
}




void FourierTransform(double *data, int nn, int isign)
{
  
  {  //check input
    int n=nn;
    while (n>1) if (n%2==1) {printf("length (nn=%i) not a power of 2 in FourierTransform",nn);exit(1);} else n/=2;
    if (isign!=1 && isign!=-1) {printf("sign (isign=%i) not +-1 in FourierTransform",isign); exit(1);}
    if (data==0) {printf("array pointer zero in FourierTransform"); exit(1);}
  }
  
  //first the bit reversal of the data
  int n=2*nn, j=1;
  for (int i=1;i<=n;i+=2) {
    if (j>i) {//swap the two
      double temp=data[j-1];
      data[j-1]=data[i-1];data[i-1]=temp;
      temp=data[j];
      data[j]=data[i];data[i]=temp;
    }
    int m=nn;
    while (m>=2 && j>m) {j-=m;m/=2;}
    j+=m;
  }
  
  //now the Danielson-Lanczos part
  int mmax=2;     //current length of the transform
  while (n>mmax) {//outer loop
    int istep=2*mmax;
    double theta=2*M_PI/(isign*mmax);
    double wpr=-2*sin(theta/2)*sin(theta/2);
    double wpi=sin(theta);
    double wr=1;
    double wi=0;
    for (int m=1;m<=mmax;m+=2) {
      for (int i=m;i<=n;i+=istep) {
        int j=i+mmax;
        double tempr=wr*data[j-1]-wi*data[j];
        double tempi=wr*data[j]+wi*data[j-1];
        data[j-1]=data[i-1]-tempr;
        data[j]=data[i]-tempi;
        data[i-1]+=tempr;
        data[i]+=tempi;
      }
      double wtemp=wr;
      wr+=wr*wpr-wi*wpi;
      wi+=wi*wpr+wtemp*wpi;
    }
    mmax=istep;
  }
}
  
void RealFourierTransform(double *data, int n, int isign)
{
  
  {//check input
    int nn=n;
    while (nn>1) if (nn%2==1) {printf("length (n=%i) not a power of 2 in RealFourierTransform",n);exit(1);} else nn/=2;
    if (isign!=1 && isign!=-1) {printf("sign (isign=%i) not +-1 in RealFourierTransform",isign); exit(1);}
    if (data==0) {printf("array pointer zero in RealFourierTransform"); exit(1);}
  }
  
  //reshuffle
  double theta=M_PI/(n/2.0);
  double c1=0.5, c2;
  if (isign==1) {c2=-0.5;FourierTransform(data,n/2,1);}
  else {c2=0.5;theta=-theta;}
  double wpr=-2*sin(theta/2)*sin(theta/2);
  double wpi=sin(theta);
  double wr=1+wpr;
  double wi=wpi;
  int n2p3=n+3;
  for (int i=2;i<=n/4;i++) {
    int i1=2*i-1;
    int i2=i1+1;
    int i3=n2p3-i2;
    int i4=i3+1;
    double h1r=c1*(data[i1-1]+data[i3-1]);
    double h1i=c1*(data[i2-1]-data[i4-1]);
    double h2r=-c2*(data[i2-1]+data[i4-1]);
    double h2i=c2*(data[i1-1]-data[i3-1]);
    data[i1-1]=h1r+wr*h2r-wi*h2i;
    data[i2-1]=h1i+wr*h2i+wi*h2r;
    data[i3-1]=h1r-wr*h2r+wi*h2i;
    data[i4-1]=-h1i+wr*h2i+wi*h2r;
    double wtemp=wr;
    wr+=wr*wpr-wi*wpi;
    wi+=wi*wpr+wtemp*wpi;
  }
  if (isign==1) {
    double h1r=data[0];
    data[0]=h1r+data[1];
    data[1]=h1r-data[1];
  } else {
    double h1r=data[0];
    data[0]=c1*(h1r+data[1]);
    data[1]=c1*(h1r-data[1]);
    FourierTransform(data,n/2,-1);
  }
}

/*void DeNoise(double delta, double v, double l, double t, int length, double* data, double &average, double& min, double& max, double &rms, double & maxd, FILE* logfile)
{
  
  //log input
  if (logfile!=0) fprintf(logfile, "de-noise called with params:\ndelta=%e, v=%e, l=%e, t=%e, length=%i\n",delta, v, l, t, length);
  
  //check input
  if (delta<0 || v<0 || l<0 || length<2 || data==0 || length>1000000) {
    if (logfile!=0) fprintf(logfile,"wrong input, exiting\n");
    exit(1);
  }

  //find the power of 2 and pad by zeros
  int N=2;
  while (N<length) N*=2;
  double* work=new double[N];
  if (work==0) {
    if (logfile!=0) fprintf(logfile,"allocation of work failed, exiting\n");
    exit(1);
  }
  for (int i=length;i<N;i++) work[i]=0;

  //copy to work removing the linear trend
  double a=(data[length-1]-data[0])/(length-1);
  double b=data[0];
  for (int i=0;i<length;i++) work[i]=data[i]-(a*i+b);

  //forward FFT
  RealFourierTransform(work, N, 1);
  
  //apply filter [smoothing e.g. 0->f(i-kc) with f(0)=1, f(x>>1)=0]
  int kc=(int) floor(N*delta*v/l);
  if (logfile!=0) fprintf(logfile,"based on characteristic dimensions the border index is kc=%i vs length=%i\n",kc,N);
  for (int i=2*kc+1;i<N;i++) work[i]=0;
  if (kc<N/2) work[1]=0;//reshuffling format of Real Data FFT
  
  //backward FFT
  RealFourierTransform(work, N, -1);
  
  //statistics
  double sumd;
  average=sumd=rms=0;
  for (int i=0;i<length;i++) {
    double we=2*work[i]/N+(a*i+b); //finish the backward FFT scale, add back the linear trend

    average+=we;  //statistics of estimate
    if (i==0 || we<min) min=we;
    if (i==0 || we>max) max=we;

    double dif=we-data[i]; //statistics of differences
    sumd+=dif;
    rms+=dif*dif;
    if (i==0 || abs(dif)>maxd) maxd=abs(dif);
    data[i]=we;
  }
  delete work;
  
  //finish the statistics
  average/=length;
  sumd/=length;
  rms=sqrt(rms/length-sumd*sumd);
}*/
  
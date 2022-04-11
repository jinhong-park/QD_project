#include "main.h"
#include "gsl/gsl_linalg.h"

struct rates {
  int n;           //dimension of the density matrix
  prec *Gamma;     //transition rates
  content *Omega;  //oscillating field matrix elements (tot)
  content *OmegaB;  //oscillating field matrix elements (mag)
  content *OmegaE;  //oscillating field matrix elements (el)
  content **Omega_detailed;
                   //field matrix elements: a field of 5 values for each couple: Bx, By, Bz, Ex, Ey (x axis defined by d_0_angle, y by d_0_angle+M_PI/2 )
  prec *e;         //energies of the states
  prec *rho;       //diagonal density matrix elements
  int a,b;         //the indexes of the resonant states
  int pinned;      //if one of the res. states is choosen in advance
  prec gamma;      //decoherence rate of the resonant states
  prec Omega2;         //square of the disturbing field matrix element
                   //between the res. states
  prec omega;      //fequency of the field
  prec domega;     //difference of the field and  
                   //omega of the resonant states
  prec I;         //field-induced rate
};


void list_m(prec *input, int n)
{
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      sprintf(buffer,"%+.1e ",input[i*n+j]);
      splachni(logg,buffer,2);
    }
    message(logg,"\n",2);
  }
  message(logg,"\n\n",2);
}

void list_m(content *input, int n)
{
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      sprintf(buffer,"%+.1e%+.1ei ",input[i*n+j].real(),input[i*n+j].imag());
      splachni(logg,buffer,2);
    }
    message(logg,"\n",2);
  }
  message(logg,"\n\n",2);
}

void list_v(prec *input, int n)
{
  for (int i=0;i<n;i++) {
    sprintf(buffer,"%+.5e ",input[i]);
    splachni(logg,buffer,2);
  }
  message(logg,"\n\n",2);
}

void list_rates(rates *data)
{
  message(logg,"transition rates:\n",2);
  list_m(data->Gamma,data->n);
  message(logg,"field matrix elements of Bx:\n",2);
  list_m(data->Omega_detailed[0],data->n);
  message(logg,"field matrix elements of By:\n",2);
  list_m(data->Omega_detailed[1],data->n);
  message(logg,"field matrix elements of Bz:\n",2);
  list_m(data->Omega_detailed[2],data->n);
  message(logg,"for the set direction (B_osc):\n",2);
  list_m(data->OmegaB,data->n);
  
  message(logg,"field matrix elements of Ex:\n",2);
  list_m(data->Omega_detailed[3],data->n);
  message(logg,"field matrix elements of Ey:\n",2);
  list_m(data->Omega_detailed[4],data->n);
  message(logg,"for the set direction (E_osc):\n",2);
  list_m(data->OmegaE,data->n);

  message(logg,"total field matrix elements:\n",2);
  list_m(data->Omega,data->n);
 
  message(logg,"energies:\n\n",2);
  list_v(data->e,data->n);
}

void list_abs(rates *data)
{
  sprintf(buffer,"resonant states: %i %i\n",data->a,data->b);
  splachni(logg, buffer, 2);
  sprintf(buffer,"domega=%.4e\t(field freq.=%.4e)\n",data->domega,data->omega);
  splachni(logg, buffer, 2);
  sprintf(buffer,"field induced rate=%.4e\n(decoh.=%.4e, field m.el. =%.4e)\n\n",data->I,data->gamma,sqrt(data->Omega2));
  splachni(logg, buffer, 2);
  sprintf(buffer,"state a=%i,\tenergy=%e, rho(a)=%e\n",data->a,data->e[data->a],data->rho[data->a]);
  splachni(logg,buffer,2);
  sprintf(buffer,"state b=%i,\tenergy=%e, rho(b)=%e\n",data->b,data->e[data->b],data->rho[data->b]);
  splachni(logg,buffer,2);
  sprintf(buffer,"difference:\tenergy=%e, rho=%e\n\n",data->e[data->b]-data->e[data->a],data->rho[data->a]-data->rho[data->b]);
  splachni(logg,buffer,2);
}

void free_rates(rates* r)
{
  delete r->Gamma;
  delete r->Omega;
  delete r->OmegaB;
  delete r->OmegaE;
  for (int i=0;i<5;i++) delete r->Omega_detailed[i];
  delete r->Omega_detailed;
  delete r->e;
  delete r->rho;
  delete r;
}

void overlap(ret_from_arpack_zn& vysl, reg1& r1, int i, int j, content* oF) {
  
  if (r1.sets!=2) {printf("max sets %i not implemented in overlap\n",r1.sets);exit(1);}
  int dim=r1.Nw*r1.sets;
  content *tmp=new content[dim];

  //magnetic field
  r1.op_ket(vysl.eigenvecs+dim*j, tmp, reg1::sigmax,reg1::set);
  oF[0]=r1.braket(vysl.eigenvecs+dim*i,0,tmp,0);
  oF[0]+=r1.braket(vysl.eigenvecs+dim*i,1,tmp,1);
  r1.op_ket(vysl.eigenvecs+dim*j, tmp, reg1::sigmay,reg1::set);
  oF[1]=r1.braket(vysl.eigenvecs+dim*i,0,tmp,0);
  oF[1]+=r1.braket(vysl.eigenvecs+dim*i,1,tmp,1);
  r1.op_ket(vysl.eigenvecs+dim*j, tmp, reg1::sigmaz,reg1::set);
  oF[2]=r1.braket(vysl.eigenvecs+dim*i,0,tmp,0);
  oF[2]+=r1.braket(vysl.eigenvecs+dim*i,1,tmp,1);
  
  //electric field
  r1.operate(vysl.eigenvecs+dim*j, tmp,0,0, reg1::lin_com, content(1,0), xc ,reg1::set);
  r1.operate(vysl.eigenvecs+dim*j, tmp,1,1, reg1::lin_com, content(1,0), xc ,reg1::set);
  oF[3]=r1.braket(vysl.eigenvecs+dim*i,0,tmp,0);
  oF[3]+=r1.braket(vysl.eigenvecs+dim*i,1,tmp,1);
  
  r1.operate(vysl.eigenvecs+dim*j, tmp,0,0, reg1::lin_com, content(1,0), yc ,reg1::set);
  r1.operate(vysl.eigenvecs+dim*j, tmp,1,1, reg1::lin_com, content(1,0), yc ,reg1::set);
  oF[4]=r1.braket(vysl.eigenvecs+dim*i,0,tmp,0);
  oF[4]+=r1.braket(vysl.eigenvecs+dim*i,1,tmp,1);
  
  delete[] tmp;
}

struct rates *fill_rates(state_info& states, int n)
#define LOG_FILL1
{
  
  prec B=parameters.values[parameters_class::B];
  prec theta_B=parameters.values[parameters_class::theta_B];
  prec phi_B=parameters.values[parameters_class::phi_B];
  prec gx=parameters.values[parameters_class::gx];
  prec gy=parameters.values[parameters_class::gy];
  prec gz=parameters.values[parameters_class::gz];
  prec d=parameters.values[parameters_class::d];
  prec phi_d=parameters.values[parameters_class::phi_d];
  prec Eosc=parameters.values[parameters_class::Eosc];
  prec phi_Eosc=parameters.values[parameters_class::phi_Eosc];
  prec Bosc=parameters.values[parameters_class::Bosc];
  prec phi_Bosc=parameters.values[parameters_class::phi_Bosc];
  prec theta_Bosc=parameters.values[parameters_class::theta_Bosc];
  
  rates* r=new rates[1];
  r->Gamma=new prec[n*n];
  r->Omega=new content[n*n];
  r->OmegaB=new content[n*n];
  r->OmegaE=new content[n*n];
  r->Omega_detailed=new content*[5];
  for (int i=0;i<5;i++) r->Omega_detailed[i]=new content[n*n];
  r->e=new prec[n];
  r->rho=new prec[n];
  r->n=n;
  r->pinned=-1;

  for (int s1=0;s1<n;s1++) {
    int u1=states.sorted2un[s1];
//    sprintf(buffer,"sorted #%i is unsorted #%i\n",is,i);
//    splachni(logg,buffer,4);
    for (int s2=0;s2<n;s2++) {
      int u2=states.sorted2un[s2];
//      sprintf(buffer,"sorted #%i is unsorted #%i\n",js,j);
//      splachni(logg,buffer,4);
        
      int index=s1*n+s2;
      if (u1 == u2) r->Gamma[s1*n+s2]=0;
      else {
        r->Gamma[index]=transition_rate(states,u1,u2,0);
        r->Gamma[index]+=transition_rate(states,u1,u2,1);
        r->Gamma[index]+=transition_rate(states,u1,u2,2);
      }
      
      content oF[5];
      overlap(*states.vysl,*states.r1, u1, u2, oF);
      
      r->OmegaB[index]=(cos(theta_Bosc)*oF[2]*gz+sin(theta_Bosc)*(oF[0]*cos(phi_Bosc)*gx+oF[1]*sin(phi_Bosc)*gy))*units.omega*Bosc;
      r->OmegaE[index]=(oF[3]*cos(phi_Eosc)+oF[4]*sin(phi_Eosc))*units.omega*Eosc;

      r->Omega[index]=r->OmegaB[index]+r->OmegaE[index];

      r->Omega_detailed[0][index]=(oF[0]*gx*cos(phi_d)+oF[1]*gy*sin(phi_d))*units.omega*Bosc;
      r->Omega_detailed[1][index]=(-oF[0]*gx*sin(phi_d)+oF[1]*gy*cos(phi_d))*units.omega*Bosc;
      r->Omega_detailed[2][index]=oF[2]*gz*units.omega*Bosc;
      r->Omega_detailed[3][index]=(oF[3]*cos(phi_d)+oF[4]*sin(phi_d))*units.omega*Eosc;
      r->Omega_detailed[4][index]=(oF[4]*cos(phi_d)-oF[3]*sin(phi_d))*units.omega*Eosc;
    }
    r->e[s1]=states.vysl->eigenvals[u1].real()*units.omega;
  }
#ifdef LOG_FILL1
  message(logg,"FILL_RATES:\nrates filled with result:\n",2);
  list_rates(r);
#endif
  return(r);
}

void find_resonant_states(struct rates *r)
{
//#define LOG_FILL2
  prec omega=r->omega;
  r->domega=1e+100;
  int imin=0,jmin=0;
  
  for (int i=0;i<r->n;i++) {
    for (int j=0;j<r->n;j++) {
      if (i==j) continue;
      if (r->pinned>-1) if (r->pinned!=i && r->pinned!=j) continue;
      prec odif=omega-fabs(r->e[i]-r->e[j]);
      if (fabs(odif)<fabs(r->domega)) {
        r->domega=odif;
        imin=i;
        jmin=j;
      }
    }
  }
  r->a=imin;
  r->b=jmin;
  
  if (r->e[imin]>r->e[jmin]) {
    r->a=jmin;
    r->b=imin;
  }
  r->domega=r->omega-(r->e[jmin]-r->e[imin]);
  
#ifdef LOG_FILL2
  message(logg,"FIND_RESONANT_STATES:states identified with result:\n",2);
  list_rates(r);
#endif
}

void matrix_eq(rates *r, prec *Ares, prec *bres)
{
//#define LOG_LINEQ

  int n=r->n;
  prec *A=new prec[n*n];
  
  for (int i=0;i<n*n;i++) A[i]=0;

  for (int i=0;i<n;i++) {
    for (int k=0;k<n;k++) {
      A[i*n+k]+=r->Gamma[k*n+i];
      A[i*n+i]-=r->Gamma[i*n+k];
    }
    if (i==r->a) {
      A[i*n+r->a]-=r->I;
      A[i*n+r->b]+=r->I;
    }
    if (i==r->b) {
      A[i*n+r->b]-=r->I;
      A[i*n+r->a]+=r->I;
    }
#ifdef LOG_LINEQ
    sprintf(buffer,"equation #%i\nline of matrix A[n x n]:\n",i);
    splachni(logg,buffer,2);
    list_v(A+i*n,n);
#endif
    for (int k=0;k<n-1;k++) A[i*n+k]-=A[i*n+n-1];
#ifdef LOG_LINEQ
    message(logg,"line of matrix Ares[n-1 x n-1] and RHS (-bres[n-1])\n",2);
    list_v(A+i*n,n);
#endif
  }

  
  for (int i=0;i<n-1;i++) {
    for (int j=0;j<n-1;j++) Ares[i*(n-1)+j]=A[i*n+j];
    bres[i]=-A[i*n+n-1];
  }
#ifdef LOG_LINEQ
  message(logg,"matrix equation bulid:\nLHS:\n",2);
  list_m(Ares,r->n-1);
  message(logg,"RHS:\n",2);
  list_v(bres,r->n-1);
#endif
}

void solve(int n, prec *A, prec *B, prec *res)
{ 
//#define LOG_LINEQ_R
  gsl_matrix_view m=gsl_matrix_view_array (A, n, n);
  gsl_vector_view b=gsl_vector_view_array (B, n);
  gsl_vector *x = gsl_vector_alloc (n);
  //int s;
  gsl_permutation * p = gsl_permutation_alloc (n);
  //gsl_linalg_LU_decomp (&m.matrix, p, &s);
  //gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
  
  gsl_vector *tau=gsl_vector_alloc(n);
  gsl_vector *norm=gsl_vector_alloc(n);
  int signum;
  gsl_linalg_QRPT_decomp (&m.matrix, tau, p, &signum, norm);
  gsl_linalg_QRPT_solve (&m.matrix, tau, p, &b.vector, x);
  gsl_vector_free(tau);
  gsl_vector_free(norm);      
      
  prec sum=0;
  for (int i=0;i<n;i++) {
    res[i]=(prec) gsl_vector_get(x,i);
    sum+=res[i];
  }
  res[n]=1-sum;
  gsl_permutation_free(p);
  gsl_vector_free(x);
#ifdef LOG_LINEQ_R
  message(logg,"linear equations solved:\n",2);
  list_v(res,n+1);
#endif
}

void compute_I(rates *r) {
  r->gamma=0;  
  for (int i=0;i<r->n;i++) r->gamma+=r->Gamma[r->a*r->n+i]+r->Gamma[r->b*r->n+i];
  r->Omega2=(r->Omega[r->a*r->n+r->b]*conj(r->Omega[r->a*r->n+r->b])).real();
  r->domega=r->omega-(r->e[r->b]-r->e[r->a]);

  if (r->gamma>0) r->I=r->gamma/(r->domega*r->domega+r->gamma*r->gamma)*r->Omega2/4;
  else r->I=0;
}

prec absorption(struct rates *r, prec omega, int mode)
{
  
//#define LOG_ABS
  prec *A=new prec[(r->n-1)*(r->n-1)];
  prec *b=new prec[r->n-1];
  
  r->omega=omega;
  if (mode!=0) find_resonant_states(r);
  compute_I(r);
  matrix_eq(r,A,b);
  solve(r->n-1,A,b,r->rho);
  prec absor=r->I*2*(r->e[r->b]-r->e[r->a])*(r->rho[r->a]-r->rho[r->b]);
  absor*=1.054e-34/1.602e-19;
#ifdef LOG_ABS
  sprintf(buffer,"ABSORPTION:\nof %e with the following parameters being used:\n",absor);
  splachni(logg,buffer,2);
  list_abs(r);
#endif
  delete[] A;
  delete[] b;
  return(absor);
}

prec find_value(struct rates *r, prec val, prec step, prec eps, int what)
{
//#define LOG_FIND1
//#define LOG_FIND2

  const int max=100;
  bool success=true;
  prec omega0=r->e[r->b]-r->e[r->a];
  prec omega;
#ifdef LOG_FIND1
  sprintf(buffer,"FIND_VALUE:\nlooking for value %.4e with resonance at %.4e with step %.4e\n",val,omega0,step);
  splachni(logg,buffer,2);
#endif
  prec val_act;
  int count=0;
  do {//first find a smaller value
    count++;
    omega=omega0+step;
    step*=2;
    val_act=absorption(r,omega,0);
    if (what==2) val_act=r->rho[r->b];
#ifdef LOG_FIND2
    fprintf(logg,"find 1 step:%i, omega=%e, actual value=%e\n",count,omega,val_act);
#endif
    if (count>max) {
      message(logg,"maximal steps reached!!!",2);
      success=false;
    } 
  } while(count<max && val_act>val); 

#ifdef LOG_FIND2    
  sprintf(buffer,"after %i steps found smaller value (%.4e)\n",count,val_act);
  splachni(logg,buffer,2);
#endif

  prec left,right; 
  if (omega0>omega) {
    left=omega;
    right=omega0;
  } else {
    left=omega0;
    right=omega;
  }
  prec m,del;
  count=0;
    
  do {//by half-interval method find the desired value
    count++;
    m=(right+left)/2;
    del=absorption(r,m,0)-val;
    if (what==2) del=r->rho[r->b]-val;
    if (del*step>0) left=m; else right=m;
#ifdef LOG_FIND2  
    fprintf(logg,"find 2 step:%i, del=%e, val=%e, l=%e, ri-l=%e \n",count, del,val,left,right-left); 
#endif
    if (count>max) {
      message(logg,"maximal steps reached!!!",2);
      success=false;
    } 
  } while (fabs(del)>eps*val && fabs(right-left)>eps*fabs(left) && success) ;
#ifdef LOG_FIND1 
  sprintf(buffer,"by half-interval method found value %.4e at freq. %.4e (eps=%.4e) after %i steps\nFIND_VALUE END\n\n",del+val,m,eps,count);
  splachni(logg,buffer,2);
#endif
  if (success) return(m); else return(-1);
}
    

void scan(struct rates *r, int a, int b, prec* scanrates, FILE* file, int where, int mode, prec eps)
{
  prec temperature=parameters.values[parameters_class::T];
  
  const int max=100;
#define LOG_SCAN
  if (r->e[a]>r->e[b]) {int aux=a;a=b;b=aux;}
  prec omega0=r->e[b]-r->e[a];
  r->pinned=a;
  r->a=a; r->b=b;
  prec a0=absorption(r,omega0,0);
  prec Ir=r->I;
  prec Ir0;
  if (r->Gamma[b*r->n+a]==0) Ir0=0; else Ir0=Ir/r->Gamma[b*r->n+a];
  prec rhobb0=r->rho[r->b];
  prec rhoaa0=r->rho[r->a];
  prec tau;
  if (temperature>0) tau=exp(-(r->e[r->b]-r->e[r->a])/units.omega/temperature);
  else
  tau=0;
  prec left=0, right=0;
  
  if (a0>0) {
    prec step=2*r->gamma*sqrt((1+tau+2*Ir0)/(1+tau));
    right=find_value(r,a0/2,step,eps,1);
    left=find_value(r,a0/2,-step,eps,1);
  }
  prec a_width=right-left;
  //if (a_width<0) a_width=0;
  
    
#ifdef LOG_SCAN
  sprintf(buffer, "SCAN: steps computed from\ntau=%.4e, Ir0=%.4e, r->gamma=%.4e\n",tau, Ir0, r->gamma);
  splachni(logg,buffer,2);
  
  sprintf(buffer, "resonant freq. of states %i and %i scanned\nhalf values (absorption) found at %.4e and %.4e (width=%.4e)\n",a,b,left,right,right-left);
  splachni(logg,buffer,2);
#endif

  left=right=0;
  if (rhobb0>0 && rhobb0<1 && Ir0>0) {
    prec step=4*r->gamma*sqrt(2*Ir0*(0.5+Ir0)+tau*(1+3*tau+3*Ir0));
    step/=sqrt(Ir0-tau*(1+tau+3*Ir0));
    right=find_value(r,rhobb0/2,step,eps,2);
    left=find_value(r,rhobb0/2,-step,eps,2);
  }
  prec rhobb_width=right-left;
  //if (rhobb_width<0) rhobb_width=0;

  #ifdef LOG_SCAN
  sprintf(buffer, "resonant freq. of states %i and %i scanned\nhalf values (rhobb) found at %.4e and %.4e (width=%.4e)\n\n",a,b,left,right,right-left);
  splachni(logg,buffer,2);
#endif
  
  scanrates[0]=r->Gamma[r->a*r->n+r->b];
  scanrates[1]=r->Gamma[r->b*r->n+r->a];
  scanrates[2]=r->gamma;
  for (int i=0;i<3;i++) scanrates[3+i]=abs(r->Omega_detailed[i][r->a*r->n+r->b]);
  scanrates[6]=abs(r->OmegaB[r->a*r->n+r->b]);
  for (int i=0;i<2;i++) scanrates[7+i]=abs(r->Omega_detailed[i+3][r->a*r->n+r->b]);
  scanrates[9]=abs(r->OmegaE[r->a*r->n+r->b]);
  scanrates[10]=sqrt(r->Omega2);
  scanrates[11]=omega0;
  scanrates[12]=a0;
  scanrates[13]=a_width;
  scanrates[14]=Ir;
  scanrates[15]=Ir0;
  scanrates[16]=rhoaa0;
  scanrates[17]=rhobb0;
  scanrates[18]=rhobb_width;
  
  if (mode==0) return;
    
  if (a0==0) return;
  prec step0=(right-left)/mode*2;
  int count=0;
  for (int sign=-1;sign<2;sign+=2) {
    prec omega=omega0,absor;
    do {
      count++;
      absor=absorption(r,omega,0);
      sprintf(buffer,"%e\t%e\t%e\n",omega-omega0,omega,absor);
      splachni(file,buffer,2);
      omega+=sign*step0;
    } while (absor/a0>eps && count<max); 
  }
#ifdef LOG_SCAN
  sprintf(buffer,"resonant freq. scanned to both sides up to relative value (eps=)%.4e\ntotally %i points written into the file\nSCAN END\n\n",eps,count);
  splachni(logg,buffer,2);
#endif
}
/*
class rho_w2special 
{
  int i,j;        // the two special states where the electron is be injected
  int n;          //number of considered states
  content* A;     //matrix Lambda
  int N;          //dimension of matrix Lambda = n*n
  vysl_from_arpack_zn vysl;   //eigensystem of Lambda
  content* S;     //matrix elements of the sigma matrixes
  state_info* states;
  prec theta, phi;//direction along which the spin will be injected into the two states
  content* rho0; //initial density matrix
  content* rho;  //density matrix at time t

  int two2one(int i, int j);
  void one2two(int n, int &i, int &j);
  content matrix_element(int i, int j, reg1::mean_value_of_operator op);
  void fill_matrix_elements();
  
  void fill_lambda();
  void diagonalize_lambda();
  void compute_rho(prec t);
  prec give_spin(int which);

public:
  rho_w2special(int i_in, int j_in, int n_in, state_info* states_in);
  ~rho_w2special();
  
  void inject_spin(prec theta, prec phi);
  content give_rho(int i, int j);
  void spin_at_t(prec &sx, prec &sy, prec &sz, prec t);
};*/

#define RHO_W2SPECIAL
//#define RHO_W2SPECIAL2
/*  
void multiply_by_lambda(int n, void* M, content* x, content *y) 
{
  content* A=(content*) M;
  for (int i=0;i<n;i++) {
    y[i]=content(0,0);
    for (int j=0;j<n;j++) {
      y[i]+=A[i*n+j]*x[j];
    }
  }
}*/
rho_w2special::rho_w2special(int i_in, int j_in, int n_in, state_info* states_in) 
{
  states=states_in;
  i=i_in;
  j=j_in;
  n=n_in;
  N=n*n;
  A=new content[N*N];
  U=new content[N*N];
  Um1=new content[N*N];
  d=new content[N];
  S=new content[3*n*n];
  rho0=new content[n*n];
  rho=new content[n*n];
  //vysl.eigenvals=new content[N];
  //vysl.eigenvecs=new content[N*N];
  
#ifdef RHO_W2SPECIAL
  fprintf(logg,"\nobject of class rho_w2special created:\nbasis comprises %i states, the two special states are %i and %i\n",n,i,j);
#endif
}

rho_w2special::~rho_w2special() 
{
  delete A;
  delete U;
  delete Um1;
  delete S;
  delete d;
  delete rho0;
  delete rho;
  //delete vysl.eigenvals;
  //delete vysl.eigenvecs;
}

int rho_w2special::two2one(int i, int j) {return(i*n+j);}

void rho_w2special::one2two(int a, int &i, int& j) {i=a/n;j=a-i*n;}

void rho_w2special::fill_lambda()
{
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      int a=two2one(i,j);
      for (int k=0;k<n;k++) {
        for (int l=0;l<n;l++) {
          int b=two2one(k,l);
          content& y=A[a*N+b];
          y=content(0,0);
          if (i==k && j==l) y=content(0,-states->omega(i,j))*1e-3*units.e/units.hbar;
          for (int plus=-1;plus<2;plus+=2) {
            for (int lambda=0;lambda<3;lambda++) {
              for (int m=0;m<n;m++) {
                if (j==l) y-=M_PI*spectral_density(*states,i,m,m,k,plus,lambda,-plus*states->omega(m,k));
                if (k==i) y-=M_PI*spectral_density(*states,l,m,m,j,plus,lambda,+plus*states->omega(l,m));
              }
              y+=M_PI*spectral_density(*states,i,k,l,j,plus,lambda,+plus*states->omega(l,j));
              y+=M_PI*spectral_density(*states,i,k,l,j,plus,lambda,-plus*states->omega(i,k));
            }
          }
        }
      }
    }
  }
}

/*
void rho_w2special::diagonalize_lambda()
{
  fill_lambda();
  arpack_zn arp(N,n*n,-1e-12,100000,true,(char*)"SM");
  arp.go_for_it(N,multiply_by_lambda,(void*) A,&vysl,false);
  arp.show_stats(&vysl,logg,4,true,true);
  arp.check_eigenset(n*n,multiply_by_lambda,N, (void*) A, vysl.eigenvals, vysl.eigenvecs,4);
}*/
  
prec rho_w2special::compute_rho(prec t)
{
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      int a=two2one(i,j);
      rho[a]=0;
      for (int b=0;b<N;b++) {
        for (int c=0;c<N;c++) {
          rho[a]+=U[b*N+a]*Um1[c*N+b]*rho0[c]*exp(d[b]*t);//!Fortran to C = transpose
        }
      }
    }
  }

  content normH=0, trace=0;
  for (int i=0;i<n;i++) {
    for (int j=i;j<n;j++) {
      if (i==j) trace+=rho[i*n+j];
      normH+=rho[i*n+j]-conj(rho[j*n+i]);
    }
  } 
  trace-=content(1,0);
#ifdef RHO_W2SPECIAL2
  sprintf(buffer,"residual norms: hermitivity:%e, trace:%e\n",abs(normH),abs(trace));
  splachni(logg,buffer,4);
#endif
  return(max(abs(normH),abs(trace)));
}

content rho_w2special::matrix_element(int i, int j, reg1::mean_value_of_operator op)
{
  content* temp=new content[states->r1->Nw*states->r1->sets];
  content* pntri=states->vysl->eigenvecs+states->r1->Nw*states->r1->sets*i;
  content*  pntrj=states->vysl->eigenvecs+states->r1->Nw*states->r1->sets*j;
  states->r1->op_ket(pntri, temp,op,reg1::set);
  content res=0;
  for (int i=0;i<states->r1->sets;i++) res+=states->r1->braket(pntrj,i,temp,i);
  delete[] temp;
  return(res);
}

void rho_w2special::fill_matrix_elements()
{
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      S[i*n+j]=matrix_element(i,j, reg1::sigmax);
      S[n*n+i*n+j]=matrix_element(i,j, reg1::sigmay);
      S[2*n*n+i*n+j]=matrix_element(i,j, reg1::sigmaz);
    }
  }
}

void rho_w2special::inject_spin(prec theta, prec phi)
{
  prec sii, sjj, sij;
  content sijc;
  
  sii=(sin(theta)*cos(phi)*S[i*n+i]+sin(theta)*sin(phi)*S[n*n+i*n+i]+cos(theta)*S[2*n*n+i*n+i]).real();
  sjj=(sin(theta)*cos(phi)*S[j*n+j]+sin(theta)*sin(phi)*S[n*n+j*n+j]+cos(theta)*S[2*n*n+j*n+j]).real();
  sijc=sin(theta)*cos(phi)*S[i*n+j]+sin(theta)*sin(phi)*S[n*n+i*n+j]+cos(theta)*S[2*n*n+i*n+j];
  sij=abs(sijc);
  
  prec aux=(sii-sjj)/sqrt(4*sij*sij+(sii-sjj)*(sii-sjj));
  prec alpha2=(1+aux)/2;
  content alpha=sqrt(alpha2);
  content beta=sqrt(1-alpha2)*exp(content(0,-arg(sijc)));
  
#ifdef RHO_W2SPECIAL
  fprintf(logg,"from matrix elements sii=%.4e, sjj=%.4e, sij=%.4e%+.4ei,\nthe coeficients are alpha=%.4e%+.4ei, beta=%.4e%+.4ei:\n",sii,sjj, sijc.real(), sijc.imag(), alpha.real(), alpha.imag(), beta.real(),beta.imag());
#endif

  for (int k=0;k<n*n;k++) rho0[k]=content(0,0);

  rho0[i*n+i]=alpha*conj(alpha);
  rho0[i*n+j]=alpha*conj(beta);
  rho0[j*n+i]=conj(alpha)*beta;
  rho0[j*n+j]=beta*conj(beta);

#ifdef RHO_W2SPECIAL
  fprintf(logg,"spin injected along (%.3e,%.3e)PI, the density matrix is:\n",theta/M_PI, phi/M_PI);
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) fprintf(logg,"%+e%+ei ",rho0[i*n+j].real(),rho0[i*n+j].imag());
      fprintf(logg,"\n");
  }
#endif
}

content rho_w2special::give_rho(int i, int j) {return(rho[i*n+j]);}

prec rho_w2special::give_spin(int which)
{ 
  if (which<0 || which>2) {printf("only three spin directions are possible in rho_w2special::spin_at_t\n"); exit(1);}
  content x=0;
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      x+=S[which*n*n+i*n+j]*rho[j*n+i];
    }
  }
  return(x.real());
}

void rho_w2special::spin_at_t(prec &sx, prec &sy, prec &sz, prec t)
{
  prec norm=compute_rho(t);
  if (norm>1e-5) printf("density matrix not Hermitian or of trace 1 !!! (residual norm:%e)\n",norm);
#ifdef RHO_W2SPECIAL2
  fprintf(logg,"rho at time %.4e:\n",t);
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) fprintf(logg,"%+.4e%+.4ei ",rho[i*n+j].real(),rho[i*n+j].imag());
    fprintf(logg,"\n");
  }
#endif
  sx=give_spin(0);
  sy=give_spin(1);
  sz=give_spin(2);
}
  

#define F77NAME(x) x ## _

extern "C" 
{
  void F77NAME(zgeev) (char* JOBVL,char*  JOBVR,int* N, content*  A, int* LDA, content*  W, content * VL, int* LDVL, content* VR, int*  LDVR, content*  WORK, int* LWORK, prec* RWORK, int* INFO );

}
  
void rho_w2special::diagonalize_lambda()
{
  char jobvl='N', jobvr='V';
  complex<double>* a=new complex<double>[N*N];
  for (int b=0;b<N;b++) for (int c=0;c<N;c++) a[b*N+c]=A[c*N+b];//!Fortran to C = transpose
  int lda=N;
  complex<double>* vl=new complex<double>[N*N];
  int ldvl=N;
  int ldvr=N;
  complex<double> work_aux;
  int lwork=-1;
  double* rwork=new double[2*N];
  int info;
  
  F77NAME(zgeev) (&jobvl, &jobvr, &N, a, &lda, d, vl, &ldvl, U, &ldvr, &work_aux, &lwork, rwork, &info);
  
  if (info<0) {printf("illegal parameter #%i in zgeev called from rho_w2special::diagonalize_lambda\n",info);exit(1);}
  if (info>0) {printf("only %i eigenvectors computed successfully in zgeev called from rho_w2special::diagonalize_lambda\n",info);exit(1);}
  
  lwork=(int) round(work_aux.real());
#ifdef RHO_W2SPECIAL
  fprintf(logg,"zgeev asks for work of dim %i in rho_w2special::diagonalize_lambda\n",lwork);
#endif
  if (lwork<1) exit(1);
  complex<double>* work=new content[lwork];
  
  F77NAME(zgeev) (&jobvl, &jobvr, &N, a, &lda, d, vl, &ldvl, U, &ldvr, work, &lwork, rwork, &info);


  if (info<0) {printf("illegal parameter #%i in zgeev called from rho_w2special::diagonalize_lambda\n",info);exit(1);}
  if (info>0) {printf("only %i eigenvectors computed successfully in zgeev called from rho_w2special::diagonalize_lambda\n",info);exit(1);}
  
  delete[] a;
  delete[] vl;
  delete[] work;
  delete[] rwork;
}

prec rho_w2special::check_eigensystem(int where)
{
  prec max=0;
  message(logg,"checking eigenvectors\nvec\teigenvalue\t\t\tnorm\t\tresidual norm\n",where);
  for (int b=0;b<N;b++) {
    prec norm=0, normdif=0;
    for (int a=0;a<N;a++) {
      content aux=0;
      for (int c=0;c<N;c++) aux+=A[a*N+c]*U[b*N+c];//!Fortran to C = transpose 
      aux-=d[b]*U[b*N+a];//!Fortran to C = transpose
      normdif+=(aux*conj(aux)).real();
      norm+=(U[b*N+a]*conj(U[b*N+a])).real();//!Fortran to C = transpose
    }
    sprintf(buffer,"%2i\t%+e%+ei\t%.4e\t%.4e\n",b,d[b].real(), d[b].imag(),norm,normdif);
    splachni(logg,buffer,where);
    if (max<normdif) max=normdif;
  }
  return(max);
} 

prec rho_w2special::check_inversion(int where)
{
  prec max=0;
  message(logg,"checking inversion\nvec\tnorm\t\tresidual norm\n",where);
  for (int b=0;b<N;b++) {
    prec norm=0, normdif=0;
    for (int a=0;a<N;a++) {
      content aux=0;
      for (int c=0;c<N;c++) aux+=Um1[c*N+a]*U[b*N+c];//!Fortran to C = transpose
      if (a==b) aux-=content(1,0);
      normdif+=(aux*conj(aux)).real();
      norm+=(U[b*N+a]*conj(U[b*N+a])).real();//!Fortran to C = transpose
    }
    sprintf(buffer,"%2i\t%.4e\t%.4e\n",b,norm,normdif);
    splachni(logg,buffer,where);
    if (max<normdif) max=normdif;
  }
  return(max);
} 

void rho_w2special::solve()
{
  int where=0;
#ifdef RHO_W2SPECIAL
  where=4;
#endif
  
  message(logg,"filling matrix elements...",1+4);
  fill_matrix_elements();
  message(logg,"...done\n",1+4);
#ifdef RHO_W2SPECIAL
  for (int which=0;which<3;which++) {
    fprintf(logg,"matrix elements of ");
    if (which==0) fprintf(logg,"sx:\n");
    if (which==1) fprintf(logg,"sy:\n");
    if (which==2) fprintf(logg,"sz:\n");
    for (int i=0;i<n;i++) {
      for (int j=0;j<n;j++) fprintf(logg,"%+e%+ei ",S[which*n*n+i*n+j].real(),S[which*n*n+i*n+j].imag());
      fprintf(logg,"\n");
    }
  }
#endif
  
  message(logg,"filling matrix Lambda...",1+4);
  fill_lambda(); 
  message(logg,"...done\n",1+4);
#ifdef RHO_W2SPECIAL
  fprintf(logg,"matrix Lambda:\n");
  for (int a=0;a<N;a++) {
    for (int b=0;b<N;b++) fprintf(logg,"%+.2e%+.2ei ",A[a*N+b].real(),A[a*N+b].imag());
    fprintf(logg,"\n");
  }
#endif
  
  message(logg,"diagonalizing Lambda...",1+4);
  diagonalize_lambda();
  prec max=check_eigensystem(where);
  if (max>1e-5) {
    message(logg,"eigenvectors wrong? !!!:\n",1+4);
    check_eigensystem(4);
  } 
  else fprintf(logg,"eigenvectors ok, maximal residual norm %.4e",max);
  message(logg,"...done\n",1+4);
#ifdef RHO_W2SPECIAL
  fprintf(logg,"matrix U:\n");
  for (int a=0;a<N;a++) {
    for (int b=0;b<N;b++) fprintf(logg,"%+.2e%+.2ei ",U[b*N+a].real(),U[b*N+a].imag());//!Fortran to C = transpose
    fprintf(logg,"\n");
  }
#endif
  
  message(logg,"inverting matrix U...",1+4);
  invert(Um1,U,N);
  max=check_inversion(where);
  if (max>1e-5) {
    message(logg,"inversion wrong? !!!:\n",1+4);
    check_inversion(4);
  }
  else fprintf(logg,"inverted matrix ok, maximal residual norm %.4e",max);
  message(logg,"...done\n",1+4);
#ifdef RHO_W2SPECIAL
  fprintf(logg,"matrix Um1:\n");
  for (int a=0;a<N;a++) {
    for (int b=0;b<N;b++) fprintf(logg,"%+.2e%+.2ei ",Um1[b*N+a].real(),Um1[b*N+a].imag());//!Fortran to C = transpose
    fprintf(logg,"\n");
  }
#endif
}


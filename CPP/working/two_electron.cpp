#include "main.h"
#include "fftw3.h"

//degree of polynomial approximation (up to three)
#define LEVEL 3   
//points to the left (and right) for the derivatives
#define MAX_L 4

void cons2four(int i, int& n, int N)
{
  int max=N/2+N%2;
  int min=-N/2;
  n=i;
  if (i>max) n-=N;
  if (n<min || n>max) {
    n=0;printf("coordinate off in cons2four\n");exit(1);
  }
}

bool four2cons(int n, int& i, int N)
{
  int max=N/2+N%2;
  int min=-N/2; 
  bool off=false;
  if (n<min) {n+=N;off=true;}
  if (n>max) {n-=N;off=true;}
  i=n;
  if (n<0) i+=N;
  if (i<0 || i>=N) {
    i=0;printf("coordinate off in four2cons\n");exit(1);
  }
  return(off);
}
  

void update(content* coefs, prec x, int px, int py, int n, int m, int Nx, int Ny)
{
  
  prec matrix[(2*MAX_L+1)*(2*MAX_L+1)];
  differential_operator(px,py,MAX_L,1,matrix);
  //load_matrix(matrix,px,py);
  
  for (int i=0;i<2*MAX_L+1;i++) {
    for (int j=0;j<2*MAX_L+1;j++) {
      prec c=matrix[i*(2*MAX_L+1)+j];
      if (c==0) continue;
      int ni,mj;
      if (four2cons(n+i-MAX_L,ni,Nx)) continue;
      if (four2cons(m+j-MAX_L,mj,Ny)) continue;
      
      coefs[ni*Ny+mj]+=x*c;
    }
  }
}

prec integral_aux(int i, int j, prec x, prec y) 
{
  prec result=0;
  prec r=sqrt(x*x+y*y);
  
  if (i==0 && j==0) result=y*log(x+r)+x*log(y+r);
  
  if (i==1 && j==0) result=y/2*r+x*x/2*log(y+r);
  if (i==0 && j==1) result=x/2*r+y*y/2*log(x+r);
  
  if (i==2 && j==0) result=x*y/6*r+x*x*x/3*log(y+r)-y*y*y/6*log(x+r);
  if (i==1 && j==1) result=r*r*r/3;
  if (i==0 && j==2) result=x*y/6*r+y*y*y/3*log(x+r)-x*x*x/6*log(y+r);
  
  if (i==3 && j==0) result=x*x*y*r/12-y*y*y*r/6+x*x*x*x*log(y+r)/4;
  if (i==2 && j==1) result=x*(2*x*x+y*y)*r/8-y*y*y*y*log(x+r)/8;
  if (i==1 && j==2) result=x*x*y*r/8+y*y*y*r/4-x*x*x*x*log(y+r)/8;
  if (i==0 && j==3) result=x*(y*y-2*x*x)*r/12+y*y*y*y*log(x+r)/4;
  
  if (i==4 && j==0) result=(2*x*x*x*y-3*x*y*y*y)*r/40+3*y*y*y*y*y*log(x+r)/40+x*x*x*x*x*log(y+r)/5;
  if (i==3 && j==1) result=(3*x*x-2*y*y)*r*r*r/15;
  if (i==2 && j==2) result=x*y*r*r*r/10-y*y*y*y*y*log(x+r)/10-x*x*x*x*x*log(y+r)/10;
  if (i==1 && j==3) result=(3*y*y-2*x*x)*r*r*r/15;
  if (i==0 && j==4) result=(2*x*y*y*y-3*y*x*x*x)*r/40+3*x*x*x*x*x*log(y+r)/40+y*y*y*y*y*log(x+r)/5;

  if (i+j>4 || i<0 || j<0) {
    printf("wrong choice in integral_aux\n");
    exit(1);
  }
  return(result);
}


prec integral(int i, int j, int n, int m,prec dx, prec dy)
{
  prec result=integral_aux(i,j,dx*n+dx/2,dy*m+dy/2);
  result+=integral_aux(i,j,dx*n-dx/2,dy*m-dy/2);
  result-=integral_aux(i,j,dx*n+dx/2,dy*m-dy/2);
  result-=integral_aux(i,j,dx*n-dx/2,dy*m+dy/2);
  return(result);
}

prec attenuation(prec theta) 
{
  prec res=1;
  if (theta!=0) res=(3.0-4*cos(theta)+cos(2*theta))*(6+theta*theta)/(3*theta*theta*theta*theta);
  return(res);
}

content* ComputeCoefficients(int Nx, int Ny, prec dx, prec dy, int level)
{
  content* coefs=new content[Nx*Ny*2];
  
  for (int i=0;i<Nx*Ny;i++) coefs[i]=0;
  
  for (int n=0;n<Nx;n++) {
    int nf;
    cons2four(n,nf,Nx);
    for (int m=0;m<Ny;m++) {
      int mf;
      cons2four(m,mf,Ny);
      
      prec qn=nf*dx;
      prec qm=mf*dy;
/*      
      update(coefs,integral(0,0,nf,mf,dx,dy),0,0,nf,mf,Nx,Ny);
      
      update(coefs,integral(1,0,nf,mf,dx,dy)/dx,1,0,nf,mf,Nx,Ny);
      update(coefs,integral(0,0,nf,mf,dx,dy)*(-qn)/dx,1,0,nf,mf,Nx,Ny);
      update(coefs,integral(0,1,nf,mf,dx,dy)/dy,0,1,nf,mf,Nx,Ny);
      update(coefs,integral(0,0,nf,mf,dx,dy)*(-qm)/dy,0,1,nf,mf,Nx,Ny);
*/
      
      for (int k=0;k<=level;k++) {
        for (int l=0;l<=level-k;l++) {
          for (int i=0;i<=k;i++) {
            for (int j=0;j<=l;j++) {
              prec f=pown(-qn,k-i)*pown(-qm,l-j)/(fac(i)*fac(k-i)*fac(j)*fac(l-j));
              f/=pown(dx,k)*pown(dy,l);
              update(coefs,f*integral(i,j,nf,mf,dx,dy),k,l,nf,mf,Nx,Ny);
            }
          }
        }
      }
    }
  }

  for (int i=0;i<Nx;i++) {
    int n;
    cons2four(i,n,Nx);
    prec ax=attenuation(((prec) n)/( (prec) Nx)*2*M_PI);
    for (int j=0;j<Ny;j++) {
      int m;
      cons2four(j,m,Ny);
      prec ay=attenuation(((prec) m)/( (prec) Ny)*2*M_PI);
      coefs[i*Ny+j]*=ax*ay*ax*ay;
    }
  }
  
  return(coefs);
}

void expandtable(content* from,int Nx, int Ny,content* to ,int Nxa,int Nya)
{
  if (Nxa<Nx || Nya<Ny) {
    printf("wrong dimensions in expandtable!!!(Nx:%i, Nxa:%i, Ny:%i, Nya:%i)\n", Nx, Nxa, Ny, Nya);
    exit(1);
  }
  for (int i=0;i<Nxa;i++) {
    for (int j=0;j<Nya;j++) {
      if (i<Nx && j<Ny) to[i*Nya+j]=from[i*Ny+j]; else to[i*Nya+j]=0;
    }
  }
} 

/*
class CE_prototype {
  state_info* I;                           //the corresponding states and region structures
  int N, Ne;                               //number of states/entries [the latter is N*(N+1)/2]
  int Nx, Ny, Nxa, Nya, blowup;            //dimensions of the non-expanded and expanded tables
  int length_t;                            //length of the expanded table
  prec hx, hy;                             //grid steps
  prec eps0;

  complex<double>** FT;                     //Fourier Transforms of state pairs - array of pointers
  bool *FT_computed;                       //whether the specific pair is computed

  complex<double>* coefs;                  //coefficients for the FT pair integration
  complex<double>* coefs2;                  //the same, retyped
  bool CoefsInitialized;                   //whether they are calculated

  complex<double>* aux_region;             //auxiliary space - the region
  complex<double>* aux_table;              //auxiliary space - the non-expanded table
  complex<double>* aux_exp_table;          //auxiliary space - the expanded table (zeros added to the previous one)

  int ij2seq(int i, int j, bool swap);     //converts i, j into sequential label, indicates whether they have been swapped
  void create_table(int i, int j);         //create the non-expanded and expanded table for state pair {i, j}
  complex<double>* compute_FT(int i, int j, bool& cc);//returns pointer to the pair FT,computes if not yet, indicates whether complex conjugation (and inversion of the k vector) is needed
public:
  CE_prototype(state_info *I, int N, int blowup);//alocate memory except for the FT pairs
  ~CE_prototype();                         //deallocate all memory
  void deallocate();                       //deallocate memory for the FT pairs (check if allocated)
  void reset(state_info* I);               //delete all "computed" flags, renew the state info (do not deallocate FT memory)
  complex<double> CoulombElement(int ui, int uj, int uk, int ul);//Coulomb element of (unsorted) states C_ijkl in meV
};*/

CE_prototype::CE_prototype(state_info*I, int N, int blowup, double MB)
{
  if (N<1 || blowup<0) {printf("wrong parameters in CE_prototype constructor N=%i, blowup=%i\n", N, blowup);exit(1);}

  CE_prototype::I=I;
  CE_prototype::N=N;
  CE_prototype::blowup=blowup;

  I->r1->give_par_fft_2(Nx,Ny);
  Nxa=Nx*(1+blowup);Nya=(1+blowup)*Ny;
  length_t=Nxa*Nya;
  I->r1->give_par(hx,hy);
  aux_region=new complex<double>[I->r1->Nw];
  aux_table=new complex<double>[Nx*Ny];
  kinv=new int[length_t];
  CoefsInitialized=false;

  Ne=N*(N+1)/2;
  FT=new complex<double>*[Ne];
  FT_computed=new bool[Ne];
  FT_to_compute=new bool[Ne];
  for (int i=0;i<Ne;i++) {FT[i]=0; FT_computed[i]=FT_to_compute[i]=0;}
  eps0=parameters.read("eps0","-");

  double mem=sizeof(complex<double>)*length_t/1e+6;
  fprintf(logg, "Initializing CE_prototype. %i (selected) states gives %i pairs. Grid %i by %i will be expanded %i times. A single FT pair consumes %f MB of memory, with %i pairs gives %.1f MB requirement. ", N, Ne, Nx, Ny, blowup, mem , Ne, mem*Ne);
  if (mem*Ne>MB) remmember=false; else remmember=true;
  const char * msg []={"processor", "memory"};
  fprintf(logg, "Running in %s extensive regime\n", msg[remmember]);
  reset(I);
  if (remmember) {
    FT[0]=new content[Ne*length_t];
    if (FT[0]==0) {printf("memory allocation failed in CE_prototype contructor\n"); exit(0);}
    for (int i=0;i<Ne;i++) FT[i]=FT[0]+i*length_t;
  } else FT[0]=new content[2*length_t];
}

CE_prototype::~CE_prototype()
{
  delete FT[0];
  delete FT;
  delete FT_computed;
  delete FT_to_compute;
  if (CoefsInitialized) delete coefs;
  delete aux_region;
  delete aux_table;
  delete kinv;
}


int CE_prototype::ij2seq(int ui, int uj, bool& swap)
{
  int s=I->unsorted2s[ui];
  int l=I->unsorted2s[uj];
  if (s<0 || l<0 || s>=N || l>=N) {printf("wrong sorted labels in ij2seq: ui=%i, s=%i, uj=%i, l=%i\n",ui,s,uj,l); exit(0);}
  swap=false;
  if (s>l) {int aux=l; l=s; s=aux; swap=true;}
  int sl=(s*N-(s*(s-1))/2+l-s);
  if (sl<0 || sl>=Ne) {printf("wrong index in ij2seq: s=%i, l=%i, sl=%i, Ne=%i, ui=%i, uj=%i, N=%i\n", s, l, sl, Ne, ui, uj, N); exit(0);}
  return(s*N-(s*(s-1))/2+l-s);
}

void CE_prototype::deallocate()
{
  int hmn=0;
  for (int i=0;i<Ne;i++) if (FT[i]!=0) {hmn++; delete FT[i]; FT[i]=0; FT_computed[i]=false;}
  fprintf(logg, "Deallocated %i pairs\n", hmn);
}

void CE_prototype::reset(state_info* I)
{
  CE_prototype::I=I;
  for (int i=0;i<Ne;i++) FT_computed[i]=false;
  if (CoefsInitialized) delete coefs;
  coefs=ComputeCoefficients(Nxa,Nya,2*M_PI/(Nxa*hx),2*M_PI/(Nya*hy),LEVEL);
  CoefsInitialized=true;
  aux_exp_table=I->Fourier.address_in(Nxa,Nya);
  for (int i=0;i<length_t;i++) aux_exp_table[i]=0;
  
  int im, jm;
  for (int i=0;i<Nxa;i++) {
    if (i==0) im=0; else im=Nxa-i;
    for (int j=0;j<Nya;j++) {
      if (j==0) jm=0; else jm=Nya-j;
      int k=i*Nya+j;
      int ki=im*Nya+jm;
      kinv[k]=ki;
      coefs[length_t+k]=coefs[ki];
    }
  }
}

void CE_prototype::create_table(int ui, int uj)
{
  complex<double>* vecs=I->vysl->eigenvecs;
  int length=I->r1->Nw*I->r1->sets;

  I->r1->braket_vec(vecs+ui*length,vecs+uj*length,aux_region,reg1::set);
  I->r1->region2table(aux_region,0,aux_table,Nx,Ny);

  //aux_exp_table=I->Fourier.address_in(Nxa,Nya);
  //for (int i=0;i<length_t;i++) aux_exp_table[i]=0;
  for (int i=0;i<Nx;i++) for (int j=0;j<Ny;j++) aux_exp_table[i*Nya+j]=aux_table[i*Ny+j];
  //expandtable(aux_table,Nx,Ny,aux_exp_table,Nxa,Nya);

  //if (ui==1 && uj==16) for (int i=0; i<Nxa*Nya;i++) fprintf(logg,"class: %i %e%+ei\n", i, aux_exp_table[i].real(), aux_exp_table[i].imag());
}

complex<double>* CE_prototype::compute_FT(int ui, int uj, bool& cc)
{
  int sl=ij2seq(ui,uj,cc);
  //fprintf(logg,"class-compute FT called with ui=%i uj=%i, resulting in sl=%i and cc=%i, computed=%i\n",ui, uj, sl, cc, FT_computed[sl])
  if (remmember) {
    if (! FT_computed[sl]) {
      if (cc) create_table(uj,ui); else create_table(ui,uj);
      //if (ui==1 && uj==16) for (int i=0; i<Nxa*Nya;i++) fprintf(logg,"class-fin: %i %e%+ei\n", i, aux_exp_table[i].real(), aux_exp_table[i].imag());
      complex<double>* aux1=FT[sl];
      complex<double>* aux2=I->Fourier.execute();
      for (int i=0;i<length_t;i++) aux1[i]=aux2[i];
    //if (ui==1 && uj==16) for (int i=0; i<Nxa*Nya;i++) fprintf(logg,"class-fout: %i %e%+ei\n", i, aux1[i].real(), aux1[i].imag());
      FT_computed[sl]=true; 
    }
    return(FT[sl]);
  }
  if (cc) create_table(uj,ui); else create_table(ui,uj);
  return(I->Fourier.execute());
}

complex<double> CE_prototype::CoulombElement(int ui, int uj, int uk, int ul)
{
  bool cc1, cc2, kinv1, kinv2;
  complex<double>* input1=compute_FT(ui,uk,cc1);
  if (!remmember) {
    complex<double>* aux=FT[0]+length_t;
    for (int i=0;i<length_t;i++) aux[i]=input1[i];
    input1=aux;
  }
  complex<double>* input2=compute_FT(uj,ul,cc2);
  if (cc1) kinv1=false; else kinv1=true;
  if (cc2) kinv2=true; else kinv2=false;

  int im, jm;
  complex<double> x1, x2, result=0;
  for (int i=0;i<Nxa;i++) {
    if (i>0) im=Nxa-i; else im=0;
    for (int j=0;j<Nya;j++) {
      if (j>0) jm=Nya-j; else jm=0;
      int k=i*Nya+j;
      int ki=im*Nya+jm;
      if (kinv[k]!=ki) {printf("kinv[k=%i]=%i, kinv=%i\n",k,kinv[k],ki);exit(1);}
      if (kinv1) x1=input1[ki]; else x1=input1[k];
      if (kinv2) x2=input2[ki]; else x2=input2[k];
      if (cc1) x1=conj(x1);
      if (cc2) x2=conj(x2);
      result+=x1*x2*coefs[k];
    }
  }
  
  if (!cc1 && cc2) {
    complex<double> result2=0;
    for (int k=0;k<length_t;k++) result2+=input1[kinv[k]]*conj(input2[kinv[k]])*coefs[k];
    if (abs(result-result2)>1e-14) {printf("the two differ by %e\n",abs(result-result2));exit(1);}
  }
          //norm     //to Joule                                    //to meV
  result*=1/(2*M_PI)*units.wavevector*units.e2prime/eps0/(1e-3*units.e);
  return(result);
}

complex<double> CE_prototype::CoulombElementTweak(int ui, int uj, int uk, int ul)
{
  if (eps0==0) return(0);
  bool cc1, cc2;
  complex<double>* input1=compute_FT(ui,uk,cc1);
  if (!remmember) {
    complex<double>* aux=FT[0]+length_t;
    for (int i=0;i<length_t;i++) aux[i]=input1[i];
    input1=aux;
  }
  complex<double>* input2=compute_FT(uj,ul,cc2);
  complex<double>* coefsinv=coefs+length_t;
  complex<double> result=0;
  if (cc1 && !cc2) for (int k=0;k<length_t;k++) result+=conj(input1[k])*input2[k]*coefs[k];
  else if (!cc1 && cc2) for (int k=0;k<length_t;k++) result+=input1[k]*conj(input2[k])*coefsinv[k];
  else if (cc1 && cc2) for (int k=0;k<length_t;k++) result+=conj(input1[k])*conj(input2[kinv[k]])*coefs[k];
  else if (!cc1 && !cc2) for (int k=0;k<length_t;k++) result+=input1[kinv[k]]*input2[k]*coefs[k];

          //norm     //to Joule                                    //to meV
  result*=1/(2*M_PI)*units.wavevector*units.e2prime/eps0/(1e-3*units.e);
  return(result);
}


//computes a Coulomb element <ij| C | kl> using Fourier transform
//init = 1 initialized coefficient
//       0 used previously initialized coefficients
//      -1 deallocate
content CoulombElement(state_info &I, int ui, int uj, int uk, int ul, int blowup, int initialize)
{
#define CE0
  static bool CoefsInitialized=false;
  static content* coefs;
  static int Nx_old=0;
  static int Ny_old=0;

  int Nx,Ny,Nxa,Nya;   
  I.r1->give_par_fft_2(Nx,Ny);
  Nxa=Nx*(1+blowup);Nya=(1+blowup)*Ny;
  //Nxa=Nx+blowup;Nya=Ny+blowup;
  int length_r=I.r1->sets*I.r1->Nw;
  int length_t=Nxa*Nya;
  content* vecs=I.vysl->eigenvecs;
  prec hx,hy;
  I.r1->give_par(hx,hy);
  
  if ((Nxa!=Nx_old || Nya!=Ny_old || !CoefsInitialized ) && initialize==0) {
    message(logg,"dimensions changed or not initialized - initialization forced in Coulomb element\n",4);
    initialize=1;
    printf("CE dims: Nxa=%i, Nya=%i\n",Nxa, Nya);
  }
  
  //destruct everything
  if (initialize!=0 && CoefsInitialized) {
     delete coefs;
     CoefsInitialized=false;
  }
  
  //construct everything
  if (initialize==1) {
    coefs=ComputeCoefficients(Nxa,Nya,2*M_PI/(Nxa*hx),2*M_PI/(Nya*hy),LEVEL);
    CoefsInitialized=true;
    Nx_old=Nxa; Ny_old=Nya;
      
#ifdef CE0
    sprintf(buffer,"Coulomb element initialized with geometry parameters:\nNx=%i\tNy=%i,\tNxa=%i\tNya=%i\nhx=%f\thy=%f\n\n",Nx,Ny,Nxa,Nya,hx,hy);
    splachni(logg,buffer,4);
#endif
  }
  
  if (initialize==-1) return(0);
  
  //check
  if (! CoefsInitialized) {
    printf("Coulomb element : not initialized!!!\n");
    exit(1);
  }
  if (Nxa!=Nx_old || Nya!=Ny_old) {
    printf("dimensions changed in Coulomb element without initialization!!!\n");
    exit(1);
  }
  
  
  //ok
  content* aux_region=new content[length_r];
  content* aux_table=new content[Nx*Ny];
  content* fourier_out, *fourier_in;
  content* fourier=new content[length_t];
  
  fourier_in=I.Fourier.address_in(Nxa,Nya);
  
  I.r1->braket_vec(vecs+ui*length_r,vecs+uk*length_r,aux_region,reg1::set);
  I.r1->region2table(aux_region,0,aux_table,Nx,Ny);
  
  /*
  prec aux1=0;
  for (int i=0;i<Nx*Ny;i++) aux1+=abs(aux_table[i])*abs(aux_table[i]);
  sprintf(buffer,"sum of squares for the first function is %e\n",aux1);
  splachni(logg,buffer,4);
  */
  
  expandtable(aux_table,Nx,Ny,fourier_in,Nxa,Nya);
  //if (ui==1 && uj==4 && uk==16 && ul==1) for (int i=0; i<Nxa*Nya;i++) fprintf(logg,"CE   : %i %e%+ei\n", i, fourier_in[i].real(), fourier_in[i].imag());
  /*prec aux2=0;
  for (int i=0;i<length_t;i++) aux2+=abs(fftw_in[i])*abs(fftw_in[i]);
  sprintf(buffer,"sum of squares for the first function after expanding is %e\n",aux2);
  splachni(logg,buffer,4);*/

  //if (ui==1 && uj==4 && uk==16 && ul==1) for (int i=0; i<Nxa*Nya;i++) fprintf(logg,"CE   -fin: %i %e%+ei\n", i, fourier_in[i].real(), fourier_in[i].imag());
  fourier_out=I.Fourier.execute();
  //if (ui==1 && uj==4 && uk==16 && ul==1) for (int i=0; i<Nxa*Nya;i++) fprintf(logg,"CE   -fout: %i %e%+ei\n", i, fourier_out[i].real(), fourier_out[i].imag());

  for (int i=0;i<Nxa;i++) {
    int im=0;
    if (i>0) im=Nxa-i;
    for (int j=0;j<Nya;j++) {
      int jm=0;
      if (j>0) jm=Nya-j;
      fourier[i*Nya+j]=fourier_out[im*Nya+jm];
      //if (ui==1 && uj==4 && uk==16 && ul==1) fprintf(logg,"CE   -inv: %i %i %i %i %e%+ei %e%+ei\n", i, j, i*Nya+j, im*Nya+jm, fourier[i*Nya+j].real(), fourier[i*Nya+j].imag(), fourier_out[i*Nya+j].real(), fourier_out[i*Nya+j].imag());
    }
  }
  
  /*
  prec aux3=0;
  for (int i=0;i<length_t;i++) aux3+=abs(fourier[i])*abs(fourier[i]);
  sprintf(buffer,"sum of squares for the first function after fourier is %e, value at (0,0) is %e , (0,1) is %e\n",aux3, abs(fourier[0]),abs(fourier[1]));
  splachni(logg,buffer,4);
  */

  
  I.r1->braket_vec(vecs+uj*length_r,vecs+ul*length_r,aux_region,reg1::set);
  I.r1->region2table(aux_region,0,aux_table,Nx,Ny);
  /*prec aux4=0;
  for (int i=0;i<Nx*Ny;i++) aux4+=abs(aux_table[i])*abs(aux_table[i]);
  sprintf(buffer,"sum of squares for the second function is %e\n",aux4);
  splachni(logg,buffer,4);*/
  
  expandtable(aux_table,Nx,Ny,fourier_in,Nxa,Nya);
  /*prec aux5=0;
  for (int i=0;i<length_t;i++) aux5+=abs(fftw_in[i])*abs(fftw_in[i]);
  sprintf(buffer,"sum of squares for the second function after expanding is %e\n",aux5);
  splachni(logg,buffer,4);*/
  
  I.Fourier.execute();
  /*prec aux6=0;
  for (int i=0;i<length_t;i++) aux6+=abs(fftw_out[i])*abs(fftw_out[i]);
  sprintf(buffer,"sum of squares for the second function after fourier is %e, value at (0,0) is %e, (0,1) is %e\n",aux6, abs(fftw_out[0]), abs(fftw_out[1]));
  splachni(logg,buffer,4);*/
    
  content result=0;
  for (int i=0;i<length_t;i++) {
    //content aux=fftw_out[i]*fourier[i]*coefs[i];
    result+=fourier_out[i]*fourier[i]*coefs[i];
    //result+=coefs[i];
    //sprintf(buffer,"%i %e %e %e %e\n",i,aux.real(), aux.imag(), fftw_in[i].real(), fftw_in[i].imag());
    //splachni(logg,buffer,2);
  }
  //if (ui==1 && uj==4 && uk==16 && ul==1) fprintf(logg,"CE   -res: %e%+ei\n", result.real(), result.imag());
  
  /*prec aux7=0;
  for (int i=0;i<length_t;i++) aux7+=abs(coefs[i])*abs(coefs[i]);
  sprintf(buffer,"sum of squares of the coefficients is %e, value at (0,0) is %e, (0,1) is %e\n",aux7, abs(coefs[0]), abs(coefs[1]));
  splachni(logg,buffer,4);*/

  
  delete[] aux_region;
  delete[] aux_table;
  delete[] fourier;
  
  prec eps0=parameters.read("eps0","-");
  
          //norm     //to Joule                                    //to meV
  result*=1/(2*M_PI)*units.wavevector*units.e2prime/eps0/(1e-3*units.e);
  
  return(result);

}

void CoulombElementCalibrate(state_info &I, int i, int j, int k, int l, prec precision, int& blowup, int mx)
{
#define COULOMBELEMENT1
  int Nx, Ny;
  I.r1->give_par_fft_2(Nx,Ny);
  prec mistake0=0;
  blowup=0;
  content result0=CoulombElement(I, i, j, k, l, blowup,1);
  sprintf(buffer,"%i %e\n",0,abs(result0));
  splachni(logg,buffer,4);

  do {
    blowup++;
    content result1=CoulombElement(I, i, j, k, l, blowup,1);
    prec mistake1=abs(result1-result0)/(abs(result1)+abs(result0));
#ifdef COULOMBELEMENT1
    sprintf(buffer, "with grid enlargement factor %i the (%i,%i,%i,%i) C. el. gave %e(abs)\nresulting in the actual error estimate of %e (previous was %e)\n",blowup, i,j,k,l, abs(result1), mistake1, mistake0);
    //sprintf(buffer,"%i %e\n",plus,abs(result1));
    splachni(logg,buffer,4);
#endif
    if (mistake1<min(mistake0,precision)) break;
    mistake0=mistake1;
    result0=result1;
    if (blowup>mx) {
      sprintf(buffer,"no decaying error under precision %e reached after %i steps; using the maximal step\n",precision, blowup); 
      splachni(logg,buffer,1+2);
      break;
    }
  } while (true);
}


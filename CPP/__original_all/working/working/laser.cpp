//optionally: remove "inverse" from braket calls and the "r->grid->volumes" factor multiplication

//insert proper sorting of the CF bases inside the CF pool class (will be necessary beyond FFA)

//collective rather than sequential state updates?

//iterate each state until convergence then move to another?

//update total profile before frequency reset?


#include "laser.h"

//relative ratio of the difference between the lasing mode and TM eigenfrequency to see whether they correspond
#define C1 1e-2 
//relative modal integral of a lasing mode to be descarded
#define C2 1e-7
//how much above should the pump be enlarged after convergence
#define C3 1e-2
//not used
#define C4 0e-0
//threshold for the relative change of the norm and phase in the lasing mode amplitude
#define C5 1e-5
//threshold for the relative precision of the zero imaginary part of the frequency in identify root
#define C6 1e-6
//error goal in scan_TM
#define C7a 1e-3
//max error tolerance in scan_TM
#define C7b 1e-2
//allowed extent of k_step in scan_TM
#define C8 10
//minimal value of LP.Nbasis
#define C9 25




laser_parameter_prototype LP;

void CF_operate(int N, void *A, content *x, content *y)
{
  region* pregion=(region*) A;
  prec hx, h2n2, shift;
  content bc;
  pregion->give_par(hx,hx);
  pregion->give_par3(bc,shift);
  h2n2=-hx*hx*LP.n2;
  //bc=content(0,1)*0.0;
  content phase=2.0-exp(bc*hx)-(1.0+bc*hx/2.0)/(1.0-bc*hx/2.0)*0.0;
#ifdef OPEN_LEFT
  y[0]=(x[1]-x[0]*phase)/(h2n2*LP.n2x[0])+x[0]*shift;//outward flow
#else
  //y[0]=(x[1]-x[0])/(h2n2*LP.n2x[0])+x[0]*shift;  //zero derivative
  y[0]=(x[1]-2.0*x[0])/(h2n2*LP.n2x[0])+x[0]*shift;  //zero wavefunction
#endif
  for (int i=1;i<N-1;i++) {
    y[i]=(x[i-1]-2.0*x[i]+x[i+1])/(h2n2*LP.n2x[i])+x[i]*shift;
  }
#ifdef OPEN_RIGHT
  y[N-1]=(x[N-2]-x[N-1]*phase)/(h2n2*LP.n2x[N-1])+x[N-1]*shift;//outward flow
#else
  y[N-1]=(x[N-2]-2.0*x[N-1])/(h2n2*LP.n2x[N-1])+x[N-1]*shift;//zero wavefunction
#endif
}
  
void CF_operate2(int N, void *A, content *x, content *y)
{
  region* pregion=(region*) A;
  pregion->operate_tot(x,y);
  //content* aux=new content[N];
  //CF_operate2(N, A, x, aux);
  //for (int i=0;i<N;i++) fprintf(logg, "%i : |y-y_alt| = %e\n",i,abs(y[i]-aux[i]));
  //delete aux;
  //exit(0);
}

void set_matrix(int N, content* H)
{
  content *x=new content[N];
  for (int i=0;i<N;i++) x[i]=0;
  for (int i=0;i<N;i++) {    
    x[i]=1.0;
    if (i>0) x[i-1]=0;
    CF_operate(N, LP.r, x, H+i*N);
  }
}

//auxiliary for the "Hamiltonian" build
content refractive_indexm2(prec x, prec y)
{
  int k;
  LP.r->coordinate2s(x,y,k);
  return(1.0/LP.n2x[k]);
}

#define DBG(N,fmt, ...) ( if (LASER>N) {for (int i=0;i<N;i++) fprintf(logg,"\t"); fprintf(logg, fmt, ##__VA_ARGS__);} )


void compute_overlaps(int N, content* v1, content * v2, content *O)
{
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      O[i*N+j]=0;
      for (int k=0;k<N;k++) O[i*N+j]+=conj(v1[i*N+k])*v2[j*N+k];
      //if (abs(O[i*N+j])>0.5) fprintf(logg,"overlap match (%e) : i=%i, j=%i\n",abs(O[i*N+j]),i,j);
    }
  }
}
    


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SIMPLE SORTING


simple_sorting_prototype::simple_sorting_prototype(int N, int entries) 
{
  if (N<0 || entries<0 || N>1000 || entries>5) {printf("input N=%i, entries=%i invalid in simple_sorting constructor\n",N, entries);exit(1);}
  simple_sorting_prototype::N=N;
  simple_sorting_prototype::entries=entries;
  value=new content[entries*N];
  coor=new double[entries];
  coor_at_zero=new double[N];
  value_at_zero=new double[N];
  valid_entries=new int[N];
  for (int i=0;i<N;i++) valid_entries[i]=0;
  estimate=new content[N];
  position=new int[N];
  error1=new double[3*N];
  error2=error1+N;  
  error3=error2+N;
  crosses_zero=new bool[N];
  index=new int[N];
  vecR=new content[N*N];
  O=new content[N*N];
}


simple_sorting_prototype::~simple_sorting_prototype() {
  delete value;
  delete coor;
  delete valid_entries;
  delete position;
  delete estimate;
  delete error1;
  delete coor_at_zero;
  delete value_at_zero;
  delete crosses_zero;
  delete index;
  delete vecR;
  delete O;
}

void simple_sorting_prototype::debug_info()
{
  fprintf(logg, "something wrong in sorting, the actual state is:\n");
  fprintf(logg, "N=%i, entries=%i, zeros=%i", N, entries, zeros);
  for (int j=0;j<entries;j++) fprintf(logg, "c[%i]=%e ",j,coor[j]);
  fprintf(logg,"\n");
  for (int i=0;i<N;i++) {
    fprintf(logg, "state %i: valid entries=%i, error1=%e, error2=%e, error3=%e, coor_at_zero=%e, value_at_zero=%e, crosses_zero=%i, position=%i, index=%i\n", i, valid_entries[i], error1[i], error2[i], error3[i], coor_at_zero[i], value_at_zero[i], crosses_zero[i], position[i], index[i]);
    for (int j=0;j<entries;j++) fprintf(logg, "v[%i]=%e%+ei ",j,value[i*entries+j].real(), value[i*entries+j].imag());
    fprintf(logg,"estimate=%e%+ei\n", estimate[i].real(), estimate[i].imag());
  }
}

content simple_sorting_prototype::give_value(int i, int offset, double* c) 
{
#if SORTING > 0 
  fprintf(logg,"\t\tsorting class: asking for value of sorted state %i (valid_entries %i) with offset %i and coordinate %p\n", i, valid_entries[i], offset, c);
#endif
  if (valid_entries[i]>offset) {
    if (c!=0) *c=coor[offset];
    return(value[i*entries+offset]); 
  }
  if (c!=0) *c=0;
  return(0);
}

void simple_sorting_prototype::copy(simple_sorting_prototype* orig)
{
#if SORTING > 0
  fprintf(logg,"\tsorting class: copying an instance from %p\n", orig);
#endif  
  for (int i=0;i<entries;i++) coor[i]=orig->coor[i];
  for (int i=0;i<N;i++) {
    position[i]=orig->position[i];
    valid_entries[i]=orig->valid_entries[i];
    for (int j=0;j<entries;j++) value[i*entries+j]=orig->value[i*entries+j];
#if SORTING > 0
    fprintf(logg,"\t\tsorting class: copy: pos[i]=%i, coor[i mod entries]=%e, ve[i]=%i val[i,0].real()=%e\n",position[i], coor[i%entries], valid_entries[i], value[i*entries+0].real());
#endif
    for (int j=0;j<N;j++) vecR[i*N+j]=orig->vecR[i*N+j];
  }
  zeros=-1;
}

int simple_sorting_prototype::zeros_crossed(double cl, double cr) 
{
#if SORTING > 0
  fprintf(logg,"\t\tsorting class: zeros crossed called with cl=%e, cr=%e\n",cl,cr);
#endif
  if (cl==cr) {cl=coor[1]; cr=coor[0];}
  zeros=0;
  int which=0;
  for (int i=0;i<N;i++) {
    coor_at_zero[i]=value_at_zero[i]=crosses_zero[i]=false;
    if (valid_entries[i]<2) continue;
    if (valid_entries[i]>entries) {printf("can not happen in zeros_crossed:i=%i, ve=%i en=%i\n",i,valid_entries[i], entries); exit(0);}
    if (coor[0]==coor[1] || (valid_entries[i]>2 && (coor[0]==coor[2] || coor[1]==coor[2]))) fprintf(logg, "polint will crash in zeros crossed: coor={%e, %e, %e}, while trying state %i for zero crossing\n",coor[0], coor[1], coor[2],i);
    content vm;
    double vl, vr, err, cm;
    if (cl==coor[1]) vl=value[i*entries+1].imag(); else vl=polint(coor, value+i*entries, valid_entries[i], cl, err).imag();
    if (cr==coor[0]) vr=value[i*entries+0].imag(); else vr=polint(coor, value+i*entries, valid_entries[i], cr, err).imag();
#if SORTING > 0
    fprintf(logg,"\t\tsorting class: trying state %i for zero crossing with vl=%e vr=%e\n",i,vl,vr);
#endif    
    if (vl*vr>0) continue;
    double cl_aux=cl, cr_aux=cr;
    do {
      cm=(cl_aux+cr_aux)/2;
      vm=polint(coor, value+i*entries, valid_entries[i], cm, err);
      if (vm.imag()*vl>0) {vl=vm.imag(); cl_aux=cm;} else {vr=vm.imag(); cr_aux=cm;}
    } while (abs((cl_aux-cr_aux)/(cl_aux+cr_aux))>1e-15 && abs(vm.imag()/vm.real())>1e-12);
    coor_at_zero[i]=cm;
    value_at_zero[i]=vm.real();
    crosses_zero[i]=true;
    zeros++;
    if (abs(value_at_zero[i])>abs(value_at_zero[which])) which=i;
  }
  #if SORTING > 0
  for (int i=0;i<N;i++) if (crosses_zero[i]) fprintf(logg,"\t\tsorting class: sorted eigenvalue %i valid entries %i crosses zero %i at coordinate %e with real value %e\n", i, valid_entries[i], crosses_zero[i], coor_at_zero[i], value_at_zero[i]);
  #endif
  #if LASER > 3
  if (zeros>0) fprintf(logg,"\t\t\t\t%i sorted values crossing(s) spotted: largest is %e at coordinate %e\n", zeros, value_at_zero[which], coor_at_zero[which]);
  #endif
  return(zeros);
}

bool simple_sorting_prototype::zero_estimate(int i, double &coor, double& value)
{
#if SORTING > 0
  fprintf(logg,"\t\tsorting class: zero estimate of state %i called for (number of zeros within the window: %i)\n",i, zeros);
#endif      
  if (zeros==-1) zeros_crossed();
  value=value_at_zero[i]; 
  coor=coor_at_zero[i];
  return(crosses_zero[i]); 
}

bool simple_sorting_prototype::extrapolate(int i, double &coordinate, double target, int from_last)
{
#if SORTING > 0
  fprintf(logg,"\t\tsorting class: extrapolate to for state %i called for target value %e\n",i, target);
#endif
  if (valid_entries[i]<2) return(false);
  double* vals=new double[entries];
  for (int j=0;j<valid_entries[i];j++) vals[j]=value[i*entries+j].imag();
  if (from_last==-1) coordinate=extrapolate_to_value(target, coor, vals, min(valid_entries[i]-1,2));
  else coordinate=extrapolate_to_value(target, coor, vals, 1);
  delete vals;
  if (coordinate==0) return(false);
  return(true); 
}

//eigvals=true according to eigenvalues, otherwise according to eigenvectors
void simple_sorting_prototype::update(content* new_value, double new_coor, content* der, int special, bool eigvals)
{
  //which set to discard: most past value, if not special decision is made
  int which=-1;
  if (special!=-1) {
    //discard the one which has the largest imaginary part of the same sign as the new value
    double v=new_value[position[special]].imag(), ext=0;
    for (int k=0;k<valid_entries[special];k++) {
      double act=value[special*entries+k].imag();
      if (act*v>0) if (which==-1 || abs(act)>abs(ext)) {which=k;ext=act;}
    }
  }
  //otherwise discard the most past value
  if (which==-1) which=entries-1;
  
  //update the groups
  int invalid=0;
  double error_tot=0;
  for (int j=which;j>0;j--) coor[j]=coor[j-1];
  coor[0]=new_coor;
  for (int i=0;i<N;i++) {
    //update the data sets
    for (int j=which;j>0;j--) value[i*entries+j]=value[i*entries+j-1];
    value[i*entries+0]=new_value[position[i]];
    if (valid_entries[i]-1<which) valid_entries[i]++;    
    //discard the history of unreliably predicted data
    if (eigvals && error2[i]>10) {valid_entries[i]=1;error2[i]=0;invalid++;}
    if (!eigvals && error3[i]>0.5) {valid_entries[i]=1;error3[i]=0;invalid++;}
    //calculate the total error
    if (eigvals) error_tot+=error2[i]; else error_tot+=error3[i];
    //fill in the derivatives
    if (der!=0) {if (valid_entries[i]>1 && (error2[i]<1 || !eigvals)) der[i]=(value[i*entries+0]-value[i*entries+1])/(coor[0]-coor[1]); else der[i]=0;}
  }
  
  //check whether the coordinates are not the same
  for (int i=0;i<N;i++) {
    int j;
    for (j=1;j<valid_entries[i];j++) if (coor[0]==coor[j]) break;
    if (j==valid_entries[i] || valid_entries[i]==0) continue;
    fprintf(logg,"WARNING: new coordinate (%e) of sorted state %i is the same as %i-th one in valid history. Decreasing the valid history.", coor[0], i, j);
    valid_entries[i]=j;
  }

#if SORTING > 0
  fprintf(logg,"\tsorting class: after update the data are (special point %i resulted in replacement of set %i):\n", special, which);
  for (int i=0;i<N;i++) {
    fprintf(logg,"\t\tgroup %i, estimate, eigenval and overlap(method=%i) errors %e, %e and %e, valid entries reset to %i\nentries: ", i, eigvals, error1[i], error2[i], error3[i], valid_entries[i]);
    for (int j=0;j<valid_entries[i];j++) fprintf(logg, "#%i=%+.5e%+.5ei ", j, value[i*entries+j].real(), value[i*entries+j].imag());
    fprintf(logg,"\n");
    if (der!=0) {
      content der_prev=0;
      if (valid_entries[i]>2) der_prev=(value[i*entries+1]-value[i*entries+2])/(coor[1]-coor[2]);
      double der_ch=abs(der[i])+abs(der_prev);
      if (der_ch>0) der_ch=abs(der[i]-der_prev)/der_ch;
      fprintf(logg,"\tderivative is %.5e%+.5e vs previous value %.5e%+.5e (rel. change %e)\n", der[i].real(), der[i].imag(), der_prev.real(), der_prev.imag(), der_ch);
    }
  }
#endif 
  
#if LASER > 3
  fprintf(logg,"\t\t\tsorting finished with total error %e (%i values failed to be predicted)\n", error_tot, invalid);
#endif
}

void simple_sorting_prototype::sort_mother(int methods, content *new_vecL,content *new_vecR, content* new_value, double new_coor, double* errors, content* der, int special)
{
  bool eigvals=false, step1=true, step2=true;
  if (methods%2) eigvals=true;methods/=2;
  if (methods%2) step1=false; methods/=2;
  if (methods%2) step2=false; methods/=2;
  
  if (step1) {
    int errorsN=0;
    if (errors!=0) errorsN=(int) (errors[0]);
    zeros=-1;
    if (valid_entries[0]==0) {
      for (int i=0;i<N;i++) {position[i]=i;error3[i]=error2[i]=error1[i]=0;}
      if (errors!=0) errors[0]=errors[1]=errors[2]=0;
    }
    else {
      for (int i=0;i<N;i++) {
        position[i]=-1;
        index[i]=i;
        estimate[i]=polint(coor, value+i*entries, valid_entries[i], new_coor, error1[i]);
      }
      //pre-sort new values according to the real part magnitude
      if (eigvals || errors!=0) picsrt_index(N, new_value, index, 0, 'R', 'd');
      if (!eigvals) compute_overlaps(N, new_vecL, vecR, O);
      
      //sort one by one according to a chosen method
      for (int i=0;i<N;i++) {
        int which=-1;
        double diff_min=0, diff=0;
        for (int j=0;j<N;j++) {
          if (position[j]!=-1) continue;
          if (eigvals) diff=abs(estimate[j]-new_value[index[i]]);
          else diff=abs(1-(conj(O[index[i]*N+j])*O[index[i]*N+j]).real());
          if (which==-1 || diff<diff_min) {diff_min=diff;which=j;}
        }
        position[which]=index[i];
        if (!eigvals) error3[which]=diff_min; else error3[which]=0;
        error2[which]=abs(estimate[which]-new_value[index[i]])/abs(estimate[which]+new_value[index[i]]);
        
  #if SORTING > 0
        fprintf(logg,"\t\tsorting class: pre-sorted new state #%i (index %i) with value %e%+ei is recognized as sorted group #%i with %i valid entries. Resulting in estimation errors: error1=%e, error2=%e, and error3=%e\n", i, index[i], new_value[index[i]].real(), new_value[index[i]].imag(), which, valid_entries[which], error1[which], error2[which], error3[which]);
  #endif
        if (i<errorsN) {//take into account only errors from few states of highest real part of the eigenvalue
          if (i==0 || errors[2]<error1[which]) errors[2]=error1[which];
          if (i==0 || errors[1]<error2[which]) errors[1]=error2[which];
          if (i==0 || errors[0]<error3[which]) errors[0]=error3[which];
        }
      }
    }
  }
  if (step2) {
    //update history
    update(new_value, new_coor, der, special, eigvals); 
    
    //store the matrix of eigenvectors
    if (!eigvals) for (int i=0;i<N;i++) for (int j=0;j<N;j++) vecR[i*N+j]=new_vecR[position[i]*N+j];
  }
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LASER PARAMETER


void laser_parameter_prototype::init() 
{
  //gamma_perp=1e+10;
  //gamma_par=1e+8;
  //g=25*units.Debye;
  //n2=1.5;
  //na=0.0025*units.Na;
  //D0=na;
  //omegaa=2*M_PI*units.c/(640*units.nm);
  //gamma_perp=omegaa/50;
  //LP.FFA=true;

  #if LASER > 0
  fprintf(logg,"laser parameters initialized: electric field units: %e [m^2/V^2]\npopulation inversion units: %e [1]\nnatural frequency %e vs wavevector*speed of light %e \n", A(), B(), omegaa, units.wavevector*units.c);
  #endif
  
  inputs=new content*[10];
  inverse=new bool[10];
  conjugate=new bool[10];
  
  Dinc=1;
}


void laser_parameter_prototype::deallocate() 
{
  delete inputs;
  delete inverse;
  delete conjugate;
}



//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CF BASIS


CF_basis_prototype::CF_basis_prototype(int dim, int N, bool allocate)
{
  if (dim<0 || N<0 || N>dim) {printf("dim=%i, N=%i is an invalid option in CF_basis_prototype constructor\n",dim,N);exit(1);}
  CF_basis_prototype::N=N;
  CF_basis_prototype::dim=dim;
  psi=new content*[N];
  psi[0]=0;
  km2=0;
  pos=new int[N];
  if (allocate) {
    km2=new content[N+1];
    psi[0]=new content[N*dim];
  }
}


CF_basis_prototype::~CF_basis_prototype()
{
  if (psi[0]!=0) delete psi[0];
  delete psi;
  if (km2!=0) delete km2;
  delete pos;
}


int CF_basis_prototype::load(ret_from_arpack_zn& vysl, double k, double k1, double k2)
{
  if (k1==k2) k1=k2=k;
  CF_basis_prototype::k=k;
  if (km2!=0 && km2!=vysl.eigenvals) fprintf(logg,"km differs in CF_basis.load at k=%e, pointer lost???\n", LP.x(k));
  km2=vysl.eigenvals;
  double* dist=new double[N];
  int hmn=0;
  for (int i=0;i<N;i++) {
    psi[i]=vysl.eigenvecs+dim*i;
    if (LP.with_shifts) km2[i]-=LP.shift_H;
    double distr=max( max(0, -km2[i].real()+k1*k1), max(0, -k2*k2+km2[i].real()) );
    if (distr==0) hmn++;
    dist[i]=distr*distr+km2[i].imag()*km2[i].imag();
    fprintf(logg,"eigenvector %i has eigenvalue km2 (%.3e%+.3ei), with k^2=%e added. From k1-k2=(%e-%e) the distance on real axis %e and total %e. States within=%i\n", i, km2[i].real(), km2[i].imag(), LP.shift_H, k1, k2, distr, dist[i], hmn);
  }
  picsrt_index(N, dist, pos, 0, 'a');
  if (k1==k2) hmn=N;
  if (hmn<C9) hmn=C9;
  if (hmn>N) hmn=N;
  for (int i=0;i<N;i++) {
    fprintf(logg,"after sorting: position %i eigenvector %i has eigenvalue km2 (%.3e%+.3ei) with distance %e (included %i)\n", i, pos[i], km2[pos[i]].real(), km2[pos[i]].imag(), dist[pos[i]], (i<hmn));
  }
  delete dist;
  return(hmn);
}


void CF_basis_prototype::normalize()
{
  prec hx,hy;
  LP.r->give_par(hx,hy);
  content* aux=LP.aux;
  for (int i=0;i<N;i++) {
    for (int k=0;k<dim;k++) aux[k]=psi[i][k]*LP.n2x[k]*LP.n2;
    content eta=1.0/sqrt(LP.r->braket(psi[i],0,aux,0,1));
    for (int k=0;k<dim;k++) psi[i][k]*=eta;
  }
  
  #ifdef LASER_CHECK
  double diff=0;        //orthonormality
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) {
      LP.inputs[0]=psi[i];
      LP.inputs[1]=psi[j];
      LP.inputs[2]=LP.n2x;
      content me=LP.n2*LP.r->braketN(3,LP.inputs, LP.inverse, LP.conjugate);
      if (i==j) me-=1.0;
      if (abs(me)>N*1e-11) fprintf(logg,"CF basis <%i %i> scalar product %e%+e difference from %i\n", i, j, me.real(), me.imag(), (i==j));
      diff+=abs(me);
    }
  }
  if (diff>N*N*1e-11) fprintf(logg,"\t\tWARNING: CF basis at frequency %e is ONO within %e\n", LP.x(k) ,diff);
  
  content* delta=new content[2*dim*dim]; //completeness
  content* deltap=delta+dim*dim;
  for (int i=0;i<dim;i++) {
    for (int j=0;j<dim;j++) {
      delta[i*dim+j]=0;
      for (int k=0;k<N;k++) delta[i*dim+j]+=psi[k][i]*psi[k][j]*LP.n2x[i]*LP.n2x[j]*LP.n2*LP.n2;
      deltap[j*dim+i]=delta[i*dim+j];
    }
  }
  content tot=0, tot2=0;
  for (int i=0;i<N;i++) {
    for (int j=0;j<dim;j++) {
      content res=LP.r->braket(delta+j*dim,0,psi[i],0,1)-LP.n2x[j]*LP.n2*psi[i][j];
      tot+=conj(res)*res;
      content res2=LP.r->braket(deltap+j*dim,0,psi[i],0,1)-LP.n2x[j]*LP.n2*psi[i][j];
      tot2+=conj(res2)*res2;
    }
  }
  double res=sqrt(tot.real())/(dim*N);
  double res2=sqrt(tot2.real())/(dim*N);
  if (res+res2>1e-14) fprintf(logg,"WARNING: delta functions integrals average residuum over first index %e over second index %e\n", res, res2);
  delete delta;
  #endif
}

void CF_basis_prototype::plot(FILE * file)
{
  for (int k=0;k<dim;k++) {
    prec x,y;
    LP.r->s2coordinate(k,x,y);
    fprintf(file, "%e %e ",x*units.length/units.nm, (LP.n2x[k]*LP.n2).real());
    for (int i=0;i<N;i++)  fprintf(file, "%e %e %e ", psi[pos[i]][k].real(), psi[pos[i]][k].imag(), pow(abs(psi[pos[i]][k]),2));
    fprintf(file, "\n");
  }
  //double k=this->k;
  //for (int i=0;i<N;i++) fprintf(logg,"CF state %i (unsorted %i) has eigenvalue2 (%.3e%+.3ei) vs k2=%.3e (diff %.3e)\n", i, pos[i], km2[pos[i]].real(), km2[pos[i]].imag(), k*k, abs(km2[pos[i]]-k*k)/(k*k));
}

//returns the correlation of basis functions i and j shifted by dx grid points
content CF_basis_prototype::correlation1(int i, int j, int dx)
{
  //shift with periodic boundary conditions
  for (int k=0;k<dim;k++) LP.aux[k]=psi[pos[j]][(k+dx) % dim];
  //scalar product
  LP.clear_braket_flags();
  LP.inputs[0]=LP.aux;
  LP.inputs[1]=psi[pos[i]];
  LP.inputs[2]=LP.n2x;
  return(LP.r->braketN(3,LP.inputs,LP.inverse, LP.conjugate)*LP.n2);  
}

void CF_basis_prototype::filter(content* f, content *fcf)
{
  LP.clear_braket_flags();
  LP.inputs[0]=LP.n2x;
  LP.inputs[1]=f;
  for (int i=0;i<LP.dim;i++) fcf[i]=0;
  for (int i=0;i<N;i++) {
    LP.inputs[2]=psi[pos[i]];
    content a=LP.r->braketN(3,LP.inputs, LP.inverse, LP.conjugate)*LP.n2;
    for (int j=0;j<LP.dim;j++) fcf[j]+=a*psi[pos[i]][j];
  }
#ifdef LASER_CHECK_EPS
  double res=0;
  for (int i=0;i<LP.dim;i++) res+=abs(fcf[i]-f[i]);
  if (res>1e-0*LP.dim) fprintf(logg,"WARNING: CF filter gave residuum %e\n",res/LP.dim);
#endif
}

content CF_basis_prototype::correlation2(int i, int j, CF_basis_prototype* CF2)
{
  LP.clear_braket_flags();
  LP.inputs[0]=psi[pos[i]];
  LP.inputs[1]=CF2->psi[CF2->pos[j]];
  LP.inputs[2]=LP.n2x;
  return(LP.r->braketN(3,LP.inputs,LP.inverse, LP.conjugate)*LP.n2);
}

content CF_basis_prototype::correlation3(int i, int j, CF_basis_prototype* CF2)
{
  double ni=LP.r->braket(psi[pos[i]],0,psi[pos[i]],0).real();
  double nj=LP.r->braket(CF2->psi[CF2->pos[j]],0,CF2->psi[CF2->pos[j]],0).real();
  content oij=LP.r->braket(psi[pos[i]],0,CF2->psi[CF2->pos[j]],0);
  return(sqrt((oij*conj(oij)).real()/(ni*nj)));
}


void CF_basis_prototype::compare_with(CF_basis_prototype* CF2)
{
  //compare eigenvalues
  for (int i=0;i<N;i++) {
    content e1=km2[i];
    content e2=CF2->km2[CF2->pos[i]];
    fprintf(logg,"index %i first basis eigenvalue %.5ei%+.5ei second basis eigenvalue %.5ei%+.5ei (abs(diff)=%.5e)\n",i,e1.real(), e1.imag(), e2.real(), e2.imag(), abs(e1-e2));
  }
  
  //calculate the overlap matrix
  double* O=new double[N*N];
  content *aux=LP.aux;
  for (int i=0;i<N;i++) {
    for (int k=0;k<dim;k++) aux[k]=psi[i][k]*LP.n2x[k]*LP.n2;
    for (int j=0;j<N;j++) {
      //!!EXPERIMENTAL
      //content braket=LP.r->braket(aux,0,CF2->psi[CF2->pos[j]], 0 ,1);
      content braket=LP.r->braket(psi[i],0,CF2->psi[CF2->pos[j]],0)/sqrt(LP.r->braket(psi[i],0,psi[i],0)*LP.r->braket(CF2->psi[CF2->pos[j]],0,CF2->psi[CF2->pos[j]],0));
      if (abs(braket)>0.5) fprintf(logg,"overlap: <1.basis %i (pos %i) | 2.basis %i (pos %i) = %.5ei%+.5ei\n", i, pos[i], j, CF2->pos[j], braket.real(), braket.imag());
      O[i*N+j]=(braket*conj(braket)).real();
    }
  }
  
  //sort 
  int method=2;				//according to the overlaps
  for (int i=0;i<N;i++) pos[i]=-1;//initially none is taken
  for (int round=1;round<3;round++) { //in the first round assign only very clear pairs
    for (int i=0;i<N;i++) {		//for each "new" state
      bool already=false;
      for (int j=0;j<N;j++) if (pos[j]==i) already=true; //if this was already assigne in the first round, skip it
      if (already) continue;
      int which=-1;
      double diff_min=0, diff=0;
      for (int j=0;j<N;j++) {		//look at all "old" states
	if (pos[j]!=-1) continue;	//which were not yet taken
	if (method==1) diff=abs(CF2->km2[CF2->pos[j]]-km2[i]);
	else diff=1-O[i*N+j];
	if (which==-1 || diff<diff_min) {diff_min=diff;which=j;}
      }
      if (diff_min<1e-2 || round==2) pos[which]=i;

      fprintf(logg,"\t\tcompare with sorting round %i. New state #%i with value %e%+ei is recognized as reference state #%i based on diff=%e (|dE|=%e) pos[%i]=%i\n",round, i, km2[i].real(), km2[i].imag(), which, diff_min, abs(km2[i]-CF2->km2[which]),which,pos[which]);
    }
  }
  delete O;
}

void CF_basis_prototype::spit_km2(FILE * file, double l2)
{
  int *index=new int[N];
  picsrt_index(N, km2, index, 0, 'r', 'a');
  for (int i=0;i<N;i++) fprintf(file, "%i %e %e\n",i,l2*km2[index[i]].real(), l2*km2[index[i]].imag());
  printf("average Re[km2] distance is %e\n", l2*(km2[N-1]-km2[0]).real()/N);
  delete index;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CF_pool_prototype
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CF_pool_prototype::CF_pool_prototype(double grain, double maxmem) 
{
  this->grain=grain; 
  maxitems=(int) floor(maxmem/(sizeof(content)*(LP.Ndiag+1)*LP.dim));
  items=0;
  k=new double[maxitems];
  k[0]=0;
  undel=new int[maxitems];
  CF_basis=new CF_basis_prototype*[maxitems];
#if LASER > 0
  fprintf(logg,"CF basis pool created with grain %e and maximum %i items filling %e Bytes of memory\n",grain, maxitems, maxmem);
#endif  
}

CF_pool_prototype::~CF_pool_prototype()
{
  delete k; 
  for (int i=0;i<items;i++) delete CF_basis[i];
  delete CF_basis;
  delete undel;
}

int CF_pool_prototype::closest(double k0)
{
  if (items<2) return(0);
  if (k0<k[0]) return(0);
  if (k0>k[items-1]) return(items-1);
  int l=0, r=items-1, m;
  do {
    m=(r+l)/2;
    if (k0<k[m]) r=m; else l=m;
  } while (r-l>1);
  if (r!=l+1) {printf("fail in CF_pool_prototype::closest\n");exit(1);}
  if (k0-k[l]<k[r]-k0) return(l); else return(r);
}



CF_basis_prototype* CF_pool_prototype::give_CF_basis(double k0, int undeletable)
{
  CF_basis_prototype* res;

  int which=closest(k0);
  double dk=(k0-k[which]), scale=units.wavevector*units.c/LP.gamma_perp;
  bool suitable=((grain>=0 && grain*abs(dk)*scale<1) || abs(dk)<=1e-12*abs(k0));
  if (items==0) suitable=false;
#if LASER > 3
  fprintf(logg,"\t\t\tCF basis pool asked for frequency %e. The closest basis #%i with frequency %e is within tolerance (dk %e vs 1/grain %e):%i\n", LP.x(k0), which, LP.x(k[which]), dk*scale, 1/grain, suitable );
#endif  
  //if there is a close enough basis, return that, with a frequency reset to k0
  if (suitable) {
    res=CF_basis[which];
    res->k=k0;
    if (undel[which]==0) undel[which]=undeletable;
    return(res);
  }
  //if got here, it means there is no basis available and a new one will be created
  res=new CF_basis_prototype(LP.dim, LP.Ndiag);
  ret_from_arpack_zn vysl;
  vysl.eigenvals=new content[LP.Ndiag+1];
  vysl.eigenvecs=new content[LP.Ndiag*LP.dim];
  if (LP.with_shifts) {
    LP.r->operate_tot_init_shift(-k0*k0);
    LP.shift_H=-k0*k0;
  }
  LP.r->operate_tot_init_boundary(k0);
  if (!LP.arpack->go_for_it(LP.dim, CF_operate,(void*) LP.r, &vysl, false)) {//diag fail
    delete vysl.eigenvals;
    delete vysl.eigenvecs;
    return(0);
  }
  int N_basis_aux=res->load(vysl, k0, k1, k2);
  res->normalize();
  //if this is not the first basis, order basis states according to an existing one
  if (items!=0) res->compare_with(CF_basis[which]); 
  //if this is the first basis, set the basis dimension according to "load"
  else {LP.Nbasis=N_basis_aux; if (LP.Nbasis<C9) LP.Nbasis=C9;}
  //insert the basis in the proper place
  if (dk>0 && items!=0) which++;
  items++;
  for (int i=items;i>which;i--) {
    CF_basis[i]=CF_basis[i-1];
    k[i]=k[i-1];
    undel[i]=undel[i-1];
  }
  CF_basis[which]=res;
  k[which]=k0;
  undel[which]=undeletable;
#if LASER > 3
  fprintf(logg,"\t\t\tnew CF_basis wit frequency %e added into the pool on position %i undeletable:%i\n", LP.x(k[which]), which, undel[which]);
#endif
  
  //if there are too many entries, delete one
  int which2=0;
  if (items==maxitems) {
    for (int i=0;i<items-1;i++) 
      if (k[i+1]-k[i]<k[which2+1]-k[which2] && which2!=which && (undel[i]==0)) {which2=i;}
    if (undel[which2]) {printf("have not found a deletable CF_basis in the pool (which2=%i)\n", which2);exit(1);}
#if LASER > 3
    fprintf(logg,"\t\t\tCF_basis %i with frequency %e will be deleted from the pool\n", which2, LP.x(k[which2]));
#endif    
    delete CF_basis[which2];
    for (int i=which2;i<items;i++) {
      CF_basis[i]=CF_basis[i+1];
      k[i]=k[i+1];
      undel[i]=undel[i+1];
    }
    items--;
  }
  return(res);
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!THRESHOLD MATRIX


TM_prototype::TM_prototype(int N) {
  TM_prototype::N=N;
  M=new content[N*N];
  Tau=new content[N*N];
  eig=new content[2*N];
  eig_der=new content[N];
  vecL=new content[2*N*N];
  vecR=new content[2*N*N];
  pos=new int[N];
  lapack=new lapack_zgee_prototype(N);
}


TM_prototype::~TM_prototype() {
  delete M;
  delete vecL;
  delete vecR;
  delete eig;
  delete eig_der;
  delete pos;
  delete lapack;
  delete Tau;
}

void TM_prototype::copy(TM_prototype& TM_source, int* index)
{
  for (int i=0;i<N;i++) {
    int j=i;
    if (index!=0) j=index[i];
    eig[i]=TM_source.eig[j]; eig[i+N]=TM_source.eig[j+N];
    eig_der[i]=TM_source.eig_der[j];
    pos[i]=TM_source.pos[j];
    for (int k=0; k<N; k++) {
      M[i*N+k]=TM_source.M[j*N+k];
      Tau[i*N+k]=TM_source.Tau[j*N+k];
      vecL[i*N+k]=TM_source.vecL[j*N+k]; vecL[i*N+k+N*N]=TM_source.vecL[j*N+k+N*N];
      vecR[i*N+k]=TM_source.vecR[j*N+k]; vecR[i*N+k+N*N]=TM_source.vecR[j*N+k+N*N];
    }
  }
}


void TM_prototype::diagonalize(bool skipvecs) {
  if (skipvecs) lapack->diagonalize(M, eig+N); else lapack->diagonalize(M,eig+N,vecL+N*N,vecR+N*N,true);
  for (int i=0;i<N;i++) pos[i]=i;
  for (int i=0;i<N;i++) for (int j=0;j<N-1;j++) {
    bool swap=false;
    if (LP.D!=0) {if (abs(eig[pos[j]+N]-1/LP.D)>abs(eig[pos[j+1]+N]-1/LP.D)) swap=true;}
    else if (eig[pos[j]+N].real()<eig[pos[j+1]+N].real()) swap=true;
    if (swap) {int aux1=pos[j];pos[j]=pos[j+1]; pos[j+1]=aux1;}
  }
  for (int i=0;i<N;i++) {
    eig[i]=eig[pos[i]+N];
    if (!skipvecs) {
      content norm=0;
      for (int j=0;j<N;j++) norm+=conj(vecL[pos[i]*N+j+N*N])*vecR[pos[i]*N+j+N*N];
      for (int j=0;j<N;j++) {
        vecL[i*N+j]=vecL[pos[i]*N+j+N*N];
        vecR[i*N+j]=vecR[pos[i]*N+j+N*N]/norm;
      } 
    }
  }
  #if LASER > 3
  fprintf(logg,"\t\t\t\tTM matrix was diagonalized with the following eigenvalues:\n");
  for (int i=0;i<N;i++) fprintf(logg,"\t\t\t#%2i. eigenvalue: %+.5e%+.5ei (abs=%e)\n",i, eig[i].real(), eig[i].imag(), abs(eig[i]));
  #endif
  #ifdef LASER_CHECK
  if (!skipvecs) {
    //double res=is_eigensys(N,N,M, eig, vecL,vecR,false);
    //if (abs(res)>1e-14) fprintf(logg,"\t\t\tWARNING: TM eigensystem precision %+.5e\n", res);
    //res=is_orthonormal(N,N,vecL,vecR,false);
    //if (abs(res)>1e-14) fprintf(logg,"\t\t\tWARNING: TM eigenvectors ONO precision %+.5e\n", res);
  }
  #endif
}

void TM_prototype::update_from_Tau(double k)
{
  for (int i=0;i<N;i++) {
    content L=LP.Lambda(k, CF_basis->km2[CF_basis->pos[i]]);
    for (int j=0;j<N;j++) M[i*N+j]=L*Tau[i*N+j];
  }
}

void TM_prototype::construct(CF_basis_prototype* CF_basis, bool Tau_only) 
{
  TM_prototype::CF_basis=CF_basis;
  #if LASER > 3
  if (Tau_only) fprintf(logg,"\t\t\tconstructing the Tau matrix from CF_basis %p\n",CF_basis);
  else fprintf(logg,"\t\t\tconstructing the TM matrix from CF_basis at frequency %e\n",LP.x(CF_basis->k));
  #endif
  content* aux=LP.aux;
  int dim=LP.dim;
  for (int i=0;i<N;i++) {
    for (int k=0;k<dim;k++) aux[k]=LP.d[k]/LP.total_profile[k]*CF_basis->psi[CF_basis->pos[i]][k];
    for (int j=0;j<=i;j++) { 
      Tau[i*N+j]=LP.r->braket(aux,0,CF_basis->psi[CF_basis->pos[j]],0,1);
      Tau[j*N+i]=Tau[i*N+j];
    }
  }
  #ifdef LASER_CHECK_EPS
  content *Tau_inv=new content[N*N];
  for (int i=0;i<N;i++) {
    for (int k=0;k<dim;k++) aux[k]=LP.total_profile[k]/LP.d[k]*LP.n2*LP.n2*LP.n2x[k]*LP.n2x[k]*CF_basis->psi[CF_basis->pos[i]][k];
    for (int j=0;j<=i;j++) { 
      Tau_inv[i*N+j]=LP.r->braket(aux,0,CF_basis->psi[CF_basis->pos[j]],0,1);
      Tau_inv[j*N+i]=Tau_inv[i*N+j];
    }
  }
  content res=0.0;
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) {
      if (i==j) res+=1.0;
      for (int k=0;k<N;k++) res-=Tau[i*N+k]*Tau_inv[k*N+j];
    }
  }
  if (abs(res)>1e-1*N*N) fprintf(logg, "WARNING: Tau * Tau_inv is one with residuum %e \n", abs(res));
  #endif
  for (int i=0;i<N;i++) {
    content L=LP.Lambda(CF_basis->k, CF_basis->km2[CF_basis->pos[i]]);
    for (int j=0;j<N;j++) {
      M[i*N+j]=L*Tau[i*N+j];
#ifdef LASER_CHECK_EPS
      Tau_inv[j*N+i]/=L;
#endif
#if LASER > 4
      fprintf(logg,"\t\t\t\tTM elements M[%i, %i] = (%+.5e%+.5ei)\n\t\t\t\t(Lambda=(%+.5e%+.5ei), k=%e, km2[%i]=(%+.5e%+.5ei)\n", i, j, Tau[i*N+j].real(), Tau[i*N+j].imag(), L.real(), L.imag(), CF_basis->k, i, CF_basis->km2[CF_basis->pos[i]].real(),   CF_basis->km2[CF_basis->pos[i]].imag());
#endif
    }
  }
#ifdef LASER_CHECK_EPS
  res=0.0;
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) {
      if (i==j) res+=1.0;
      for (int k=0;k<N;k++) res-=M[i*N+k]*Tau_inv[k*N+j];
    }
  }
  if (abs(res)>1e-1*N*N) fprintf(logg, "WARNING: M * M_inv is one with residuum %e \n", abs(res));
  delete Tau_inv;
#endif
  

  
#ifdef LASER_CHECK_EPS
  content * M_CF=new content[N*N];
  content * M_invCF=new content[N*N];
  content * aux2=new content[dim];
  content * aux_cf=new content[dim];
  for (int k=0;k<dim;k++) aux2[k]=LP.n2x[k]*LP.total_profile[k]/LP.d[k];
  for (int i=0;i<N;i++) {//inverse TM run through CF filter
    for (int k=0;k<dim;k++) aux[k]=CF_basis->psi[CF_basis->pos[i]][k]*aux2[k];
    CF_basis->filter(aux, aux_cf);
    LP.inputs[0]=aux_cf;
    LP.inputs[0]=aux;
    LP.inputs[1]=LP.n2x;
    for (int j=0;j<N;j++) {
      LP.inputs[2]=CF_basis->psi[CF_basis->pos[j]];
      M_invCF[i*N+j]=LP.r->braketN(3,LP.inputs,LP.inverse, LP.conjugate)*LP.n2*LP.n2;
      //fprintf(logg,"M_invCF[%i,%i]=%e%+ei vs M_inv[%i,%i]=%e%+ei\n", i, j, M_invCF[i*N+j].real(), M_inv[i*N+j].imag(), i, j,  M_inv[i*N+j].real(), M_inv[i*N+j].imag());
    }
  }
  for (int i=0;i<N;i++) {//TM run through CF filter
    for (int k=0;k<dim;k++) aux[k]=CF_basis->psi[CF_basis->pos[i]][k]/aux2[k];
    CF_basis->filter(aux, aux_cf);
    LP.inputs[0]=aux_cf;
    LP.inputs[0]=aux;
    LP.inputs[1]=LP.n2x;
    for (int j=0;j<N;j++) {
      LP.inputs[2]=CF_basis->psi[CF_basis->pos[j]];
      M_CF[j*N+i]=LP.r->braketN(3,LP.inputs,LP.inverse, LP.conjugate);
      //fprintf(logg,"M_CF[%i,%i]=%e%+ei vs M[%i,%i]=%e%+ei\n", i, j, M_CF[i*N+j].real(), M_CF[i*N+j].imag(), i, j,  M[i*N+j].real(), M[i*N+j].imag());
    }
  } //check they give unity
  res=0.0;
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) {
      if (i==j) res+=1.0;
      for (int k=0;k<N;k++) res-=M_invCF[i*N+k]*M_CF[k*N+j];
    }
  }
  if (abs(res)>1e-1*N*N) fprintf(logg, "WARNING: M_invCF * M_CF is one with residuum %e \n", abs(res));
  delete aux2;
  delete aux_cf;
  delete M_CF;
  delete M_invCF;
#endif
  
}



//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LASING MODE


lasing_mode_prototype::lasing_mode_prototype(int dim, int N, double k) {
  lasing_mode_prototype::k=k;
  k_last_update=0;
  A=new content[N];
  A_save=new content[N];
  profile=new double[dim];
  profile_save=new double[dim];
  amplitude=new content[dim]; 
  c=new content[N*N];
  d=new content[N*N];
  e=new content[20];
  TM=new TM_prototype(N);
  lasing_mode_prototype::CF_basis=LP.CF_pool->give_CF_basis(k);
}


lasing_mode_prototype::~lasing_mode_prototype() {
  delete A;
  delete profile;
  delete amplitude;
  delete c;
  delete d;
  delete e;
  delete TM;
}


void lasing_mode_prototype::init_guess(double lambda)
{
  //TM->construct(lasing_mode_prototype::CF_basis);
  TM->construct(LP.CF_pool->give_CF_basis(lasing_mode_prototype::k));
  TM->diagonalize(false);
  int hmn=4, N=TM->N;
  evaluate_e(hmn);
  //double lambda=TM->eig[0].real();
  double gamma=LP.Gamma(k);
  double omega=sqrt((lambda-1/LP.D)/(gamma*e[1].real()));
  #if LASER > 1
  fprintf(logg,"\tnew lasing state estimate for amplitude is %.5e, from Re(lambda)=%.5e, 1/D=%.5e, gamma=%.5e, e=%.5e%+.5ei (e[0] equals lambda with residuum %.5e)\n", omega, lambda, 1/LP.D, gamma, e[1].real(), e[1].imag(), abs(e[0]-lambda));
  #endif
  if (lambda<1/LP.D) {printf("inconsistency: lambda(%.5e) - 1/D(%.5e) < 0 (%.5e)\n", lambda, 1/LP.D, lambda-1/LP.D); exit(1);}
  #if LASER > 4
  for (int j=0;j<hmn;j++) fprintf(logg,"\t\t\t\tcoefficient e[%i]=(%.5e%+.5ei)\n", j, e[j].real(), e[j].imag());
  #endif
  for (int i=0;i<N;i++) A[i]=omega*TM->vecR[0+i];
  update_profile();
}


void lasing_mode_prototype::update_profile() {
  double gamma=LP.Gamma(k);
  double sum_old=0, sum_new=0;
  norm=0;
  for (int k=0;k<LP.dim;k++) {
    content aux=0;
    for (int j=0;j<LP.Nbasis;j++) aux+=A[j]*CF_basis->psi[CF_basis->pos[j]][k];
    amplitude[k]=aux;
    double aux2=(aux*conj(aux)).real();
    norm+=aux2;
    sum_old+=profile[k];
    profile[k]=gamma*aux2;
    sum_new+=profile[k];
  }
  norm=sqrt(norm);
  #if LASER > 3
  fprintf(logg,"\t\t\tupdating profile of lasing mode at frequency %e with gamma %e to %e from %e (rel. change by %.5e)\n", LP.x(this->k), gamma, sum_new, sum_old, 1-sum_old/sum_new);
  #endif
}

void lasing_mode_prototype::update_modal_power()
{
  double old=modal_power, modal_integral2=0;
  modal_integral=0;
  for (int k=0;k<LP.dim;k++) {
    modal_integral+=profile[k]*LP.d[k].real()/LP.total_profile[k];
    modal_integral2+=profile[k]*LP.n2x[k].imag();
  }
  modal_integral*=LP.D*LP.Gamma(lasing_mode_prototype::k);
  modal_integral2*=LP.n2;
  modal_power=modal_integral-modal_integral2;
  
  //alternative expression in 1D : directly the outward flux
  //at the rhs
  prec outfluxr=-((amplitude[LP.dim-1]+amplitude[LP.dim-1])*conj(amplitude[LP.dim-1]-amplitude[LP.dim-2])).imag()/2.0;
  //at the lhs
  prec outfluxl=-((amplitude[0]+amplitude[0])*conj(amplitude[0]-amplitude[1])).imag()/2.0;
  prec hx,hy;
  LP.r->give_par(hx,hy);
  prec factor=lasing_mode_prototype::k*hx*units.wavevector*units.length;
  prec modal_power_alt=(outfluxl*0+outfluxr)/(factor*factor);
  //alternative: using the boundary condition, the outward derivative is equal to i k Psi, so that the modal power is proportional to |psi|^2
  prec modal_power_alt2=(amplitude[LP.dim-1]*conj(amplitude[LP.dim-1])).real()/factor;

#if LASER > 3
  fprintf(logg,"\t\t\tupdating modal power of lasing mode at frequency %e to %e from the old %e (rel. change by %.5e) (mode integral %e model_integral2=%e, power from integral %e)\n", LP.x(lasing_mode_prototype::k), modal_power_alt2, old, 1-old/modal_power, modal_integral, modal_integral2, modal_power);
  
  fprintf(logg,"\t\t\talternative modal power: from lhs %e (ignored) and rhs %e divided by the square of the frequency factor %e it follows as %e; yet alternative %e\n", outfluxl, outfluxr, factor, modal_power_alt, modal_power_alt2);
#endif
  //if the integral gives nonsense, accept the value based on derivative
  if (modal_power<0 && modal_power_alt>0) modal_power=modal_power_alt; 
  modal_power=modal_power_alt2;
}


double lasing_mode_prototype::update_frequency()
{
  int N=LP.Nbasis;
  int which=0;//find largest amplitude component
  for (int i=1;i<N;i++) if (abs(A[i])>abs(A[which])) which=i;

  content x=0;//calculate the product Tau * A
  for (int i=0;i<N;i++) x+=TM->Tau[which*N+i]*A[i];
  x/=A[which];
    
  //coefficients for the cubic equation
  content k2=lasing_mode_prototype::CF_basis->km2[lasing_mode_prototype::CF_basis->pos[which]];
  double R=k2.real(), I=k2.imag(), a=x.real(), b=x.imag();
  double omega=LP.omegaa/units.c/units.wavevector;
  double gamma=LP.gamma_perp/units.c/units.wavevector;
  double coef[4];
  coef[0]=-b;
  coef[1]=b*omega+a*gamma;
  coef[2]=b*R-a*I;
  coef[3]=-b*R*omega-b*I*gamma+a*I*omega-R*a*gamma;
  content root[3];
  int n=cubic(coef, root);
  
  which=0;//choose the closest root & check it
  double k=CF_basis->k;
  for (int i=1;i<n;i++) if (abs(root[i].real()-k)<abs(root[which].real()-k)) which=i;
#if LASER > 3
  fprintf(logg, "\t\t\tlasing mode frequency updated from %.5e to %.5e (accepting root %i)\n", LP.x(k), LP.x(root[which].real()), which);
  fprintf(logg, "\t\t\tcubic equation %+.2e x^3 %+.2e x^2 %+.2e x %+.2e gave %i real roots (all roots: (%+.2e%+.2ei,%+.2e%+.2ei,%+.2e%+.2ei)\n",coef[0], coef[1], coef[2], coef[3], n, root[0].real(), root[0].imag(), root[1].real(), root[1].imag(), root[2].real(), root[2].imag() );
  k=root[which].real();
  content y=coef[0]*k*k*k+coef[1]*k*k+coef[2]*k+coef[3];
  fprintf(logg, "\t\t\tauxilliary coefs: r=%.5e, i=%.5e, a=%.5e, b=%.5e, omega=%.5e, gamma=%.5e, chosen root residuum %e\n", R, I, a, b, omega, gamma, abs(y));
#endif
  return(root[which].real());
}

void lasing_mode_prototype::evaluate_e(int hmn)
{
  int dim=LP.dim;
  content *amp=new content[4*dim];
  content *ampi=amp+dim;
  content *ampL=amp+2*dim;
  content *hole=amp+3*dim;
  for (int i=0;i<dim;i++) amp[i]=ampL[i]=0;
  for (int j=0;j<TM->N;j++) {//components of an eigenvector
    for (int k=0;k<dim;k++) {//over grid
      amp[k]+=TM->vecR[j]*CF_basis->psi[CF_basis->pos[j]][k];
      ampL[k]+=conj(TM->vecL[j])*CF_basis->psi[CF_basis->pos[j]][k]*LP.Lambda(this->k,CF_basis->km2[CF_basis->pos[j]]);
    }
  }
  for (int i=0;i<dim;i++) {
    ampL[i]*=LP.d[i]*amp[i];
    ampi[i]=1.0;
    hole[i]=1.0/LP.total_profile[i];
  }

  LP.clear_braket_flags();
  LP.inputs[0]=ampL;
  LP.inputs[1]=hole;
  LP.inputs[2]=ampi;
  for (int i=0;i<hmn;i++) {
    e[i]=LP.r->braketN(3,LP.inputs, LP.inverse, LP.conjugate);
    for (int j=0;j<dim;j++) {
      ampi[j]*=amp[j]*conj(amp[j]);
      hole[j]/=LP.total_profile[j];
    }
  }
  delete amp;
  
  #if LASER > 4
  for (int i=0;i<hmn;i++) fprintf(logg,"\t\t\t\te coefficients: e[%i]=%.5e%+.5ei\n", i, e[i].real(), e[i].imag());
  #endif
}


void lasing_mode_prototype::evaluate_cd(int hmn)
{
  int N=TM->N;
  if (hmn==0) hmn=N;
  //find the "zeroth" mode
  int which=0;
  for (int i=1;i<N;i++) if (abs(TM->eig[i]-1/LP.D)<abs(TM->eig[which]-1/LP.D)) which=i;
  if (which!=0) {printf("closest to 1/D is not number zero (but %i) in evaluate_cd\n",which);exit(0);}

  //compute the auxiliary amplitudes
  int dim=LP.dim;
  content *amp=new content[(2*hmn+1)*dim];
  content *ampL=amp+hmn*dim;
  content *da0p=amp+2*hmn*dim;
  for (int i=0;i<hmn*dim;i++) amp[i]=ampL[i]=0;
  for (int i=0;i<hmn;i++) {//eigenvectors of threshold matrix
    for (int j=0;j<N;j++) {//components of an eigenvector
      for (int k=0;k<dim;k++) {//over grid
        amp[i*dim+k]+=TM->vecR[i*N+j]*CF_basis->psi[CF_basis->pos[j]][k];
        ampL[i*dim+k]+=conj(TM->vecL[i*N+j])*CF_basis->psi[CF_basis->pos[j]][k]*LP.Lambda(this->k,CF_basis->km2[CF_basis->pos[j]]);
      }
    }
  }
  for (int i=0;i<dim;i++) da0p[i]=LP.d[i]*amp[which*dim+i]/(LP.total_profile[i]*LP.total_profile[i]);
  
  LP.clear_braket_flags();
  
  for (int j=0;j<hmn;j++) {
    LP.inputs[0]=da0p;
    LP.inputs[1]=ampL+j*dim;
    for (int k=0;k<hmn;k++) {
      LP.inputs[2]=amplitude;
      LP.inputs[3]=amp+k*dim;
      LP.conjugate[3]=true;
      c[j*N+k]=LP.r->braketN(4,LP.inputs, LP.inverse, LP.conjugate);
      LP.conjugate[2]=true;
      LP.conjugate[3]=false;  
      d[j*N+k]=LP.r->braketN(4,LP.inputs, LP.inverse, LP.conjugate);
      
      #if LASER > 4
      fprintf(logg,"\t\t\t\tcd coefficients: c[%i,%i]=%.5e%+.5ei\td[%i,%i]=%.5e%+.5ei\n", j, k, c[j*N+k].real(), c[j*N+k].imag(), j, k, d[j*N+k].real(), d[j*N+k].imag());
      #endif
    }
  }
#ifdef LASER_CHECK
  content *aux=LP.aux;
  for (int i=0;i<dim;i++) aux[i]=LP.d[i]/LP.total_profile[i];
  double tot=0;
  for (int i=0;i<hmn;i++) {
    for (int j=0;j<hmn;j++) {
      LP.inputs[0]=aux;
      LP.inputs[1]=ampL+i*dim;
      LP.inputs[2]=amp+j*dim;
      LP.conjugate[2]=LP.conjugate[3]=false;
      content check1=LP.r->braketN(3,LP.inputs, LP.inverse, LP.conjugate);
      content check2=0;
      for (int k=0;k<N;k++) {
        content aux2=0;
        for (int l=0;l<N;l++) aux2+=TM->M[k*N+l]*TM->vecR[j*N+l];
        check2+=conj(TM->vecL[i*N+k])*aux2;
      }
      #if LASER > 4
      fprintf(logg,"\t\t\t\tcd coefficients check: <a[%i] | TM | a[%i]>: 1.way %e 2.way %e (rel. prec. %e)\n", i, j, abs(check1), abs(check2), abs(check1-check2));
      #endif
      tot+=abs(check1-check2);
    }
  }
  if (tot>hmn*hmn*1e-14) fprintf(logg,"\t\tWARNING: cd coefficients check: 1.way vs 2.way are the same with precision %e\n", tot);
  #endif
  delete amp;
}


void lasing_mode_prototype::perturbative_amplitude_update(content *A_new, int hmn)
{
  int N=TM->N;
  double gamma=LP.Gamma(k);
  double D=LP.D;
  int which=0;
  for (int i=1;i<N;i++) if (abs(TM->eig[i]-1/D)<abs(TM->eig[which]-1/D)) which=i;
  if (which!=0) {printf("closest eigenvalue to 1/D is not zero (but %i) in update_amplitude\n",which);exit(0);}
  if (hmn==0) hmn=N;

  //express old amplitudes in TM ONO basis
  content *omega=new content[N*2];
  content *dm=omega+N;
  for (int i=0;i<N;i++) {
    omega[i]=dm[i]=0;
    for (int j=0;j<N;j++) omega[i]+=conj(TM->vecL[i*N+j])*A[j];
  }

  #if LASER > 2
  fprintf(logg,"\t\t\tlasing mode amplitude perturbative update: at Gamma=%e, inverse pump 1/D=%e, the closest TM eigenvalue (%i) is lambda-1/D=%.5e%+.5ei\n", gamma, 1/D, which, TM->eig[which].real()-1/D, TM->eig[which].imag());
  #endif

  //estimate the frequency at which the eigenvalue is real
  content lambda=TM->eig[which], dlambda=TM->eig_der[which];
  content* eig_ren=new content[N];
  double dk=-lambda.imag()/dlambda.imag();
  for (int i=0;i<N;i++) eig_ren[i]=TM->eig[i]+dk*TM->eig_der[i];
  bool estimate=true;//whether to estimate frequency to eigenvalue real
  if (dlambda.imag()==0 || abs(dk/k)>1e-1) estimate=false;
  #if LASER > 2
  fprintf(logg,"\t\t\t\tfrom zeroth eigenvalue %.5e%+.5ei and its derivative %.5e%+.5ei, the frequency adjustement is %.3e (relative %.3e estimates will be used: %i)\n", lambda.real(), lambda.imag(), dlambda.real(), dlambda.imag(), LP.x(dk), abs(dk/k), estimate );
  if (estimate) {
    fprintf(logg, "\t\t\t\testimated eigenvalues: \n\t\t\t\tindex\teigenvalue\t\t\tincrement\t\t\trelative\n");
    for (int i=0;i<TM->N*0+2;i++) fprintf(logg, "\t\t\t\t%i\t%.5e%+.5ei\t%.5e%+.5ei \t%.5e\n", i, TM->eig[i].real(), TM->eig[i].imag(), eig_ren[i].real(), eig_ren[i].imag(), abs(eig_ren[i]/TM->eig[i])-1 );
  }
  #endif

  int hmn2=2*hmn;

  double* M=new double[2*hmn2*hmn2];
  double* Mc=M+hmn2*hmn2;
  double* V=new double[2*hmn2];
  double* Vc=V+hmn2;
  for (int j=0;j<hmn;j++) {//equation index
    for (int k=0;k<hmn;k++) {//unknown index
      content aux=d[j*N+k];
      if (j==k && j!=which) aux+=(1/LP.D-TM->eig[j])/(gamma*omega[which]);
      content p=aux+conj(c[j*N+k]);
      content m=aux-conj(c[j*N+k]);
      M[(0+j)*hmn2+(0+k)]=p.real();//real part of equation, multiplicator at real part of k-th unknown
      M[(0+j)*hmn2+(hmn+k)]=-p.imag();//real part of equation, multiplicator at imag part of k-th unknown
      M[(hmn+j)*hmn2+(0+k)]=m.real();//imag part of equation, multiplicator at real part of k-th unknown
      M[(hmn+j)*hmn2+(hmn+k)]=m.imag();//imag part of equation, multiplicator at imag part of k-th unknown
    }
    content rhs;
    if (estimate) rhs=omega[j]*(eig_ren[j]-1/LP.D)/(omega[which]*gamma);
    else { 
      if (j!=which) rhs=omega[j]*(TM->eig[j]-1/LP.D)/(omega[which]*gamma);
      else rhs=omega[j]*(TM->eig[j].real()-1/LP.D)/(omega[which]*gamma);
    }
    V[0+j]=rhs.real();
    V[hmn+j]=rhs.imag();
  }
  //copy the matrix and the rhs vector into auxilliary arrays Mc (transpose for fortran!!), Vc
  for (int i=0;i<hmn2;i++) {
    //fprintf(logg,"equation #%i: lhs: ",i);
    for (int j=0;j<hmn2;j++) {
      Mc[i*hmn2+j]=M[j*hmn2+i];//transpose for fortran
      //fprintf(logg," %+.3e x[%i] ",Mc[i*hmn2+j],j);
    }
    Vc[i]=V[i];
    //fprintf(logg,"rhs=%.3e\n",Vc[i]);
  }
  int *IPIV=new int[hmn2];
  int INFO,NRHS=1;
  if (hmn>0) dgesv_(&hmn2, &NRHS, Mc, &hmn2, IPIV, Vc, &hmn2, &INFO);
  double tot=0, tot2=0;
  for (int i=0;i<hmn2;i++) {//check residuum
    double rhs=V[i];
    tot2+=abs(rhs);
    for (int j=0;j<hmn2;j++) rhs-=M[i*hmn2+j]*Vc[j];
    tot+=abs(rhs);
  }
  if (tot/tot2>1e-14*N || INFO!=0) fprintf(logg,"WARNING: linear equations not solved? (INFO=%i, residuum=%e, tot=%e, tot2=%e)\n", INFO, tot/tot2, tot, tot2);
  #if LASER > 4
  fprintf(logg,"\t\t\t\tlinear equations solved with flag INFO=%i, fulfilled with residuum %e\n",INFO,tot);
  #endif  
  for (int i=0;i<hmn;i++) dm[i]=content(Vc[i], Vc[i+hmn]);

  #if LASER > 3
  fprintf(logg,"\t\t\tperturbative update of lasing mode amplitude gave the following changes in TM ONO basis:\n\t\t\tindex\told value\t\t\tincrement\t\t\trelative change\n");
  for (int i=0;i<hmn;i++) {
    double diff=abs(dm[i]/omega[i]);
    fprintf(logg,"\t\t\t%i\t%.5e%+.5ei\t%.5e%+.5ei\t%.5e\n", i, omega[i].real(), omega[i].imag(), dm[i].real(), dm[i].imag(), diff);
  }
  #endif
  delete IPIV;
  delete M;
  delete V;
  //write down the resulting amplitudes
  for (int i=0;i<N;i++) {
    A_new[i]=0;
    for (int j=0;j<N;j++) A_new[i]+=(omega[j]+dm[j])*TM->vecR[j*N+i];
  }
  delete omega;
  delete eig_ren;
}


content lasing_mode_prototype::iterative_step(content *A_new)
{
  int N=TM->N;
  //check if the closest to real eigenvalue is also the largest one
  bool ok=true;
  for (int i=0;i<N;i++) {
    A_new[i]=0;
    if (abs(TM->eig[i])>abs(TM->eig[0])) ok=false;
  }
  
  if (ok) {//if yes, multiply the vector by the matrix
    for (int i=0;i<N;i++) {
      content* M=TM->M+i*N, aux=0;
      for (int j=0;j<N;j++) aux+=M[j]*A[j];
      A_new[i]=aux*LP.D;
    }
  } else {//if not, exclude the potentially problematic TM eigenvectors
    TM->diagonalize(false);
    for (int i=0;i<N;i++) {//calculate projections of vector A with TM eigenvectors
      content omega=0;
      if (abs(TM->eig[i])>abs(TM->eig[0])) continue;
      for (int j=0;j<N;j++) omega+=conj(TM->vecL[i*N+j])*A[j];
      omega*=TM->eig[i]*LP.D;
      content * M=TM->vecR+i*N;
      for (int j=0;j<N;j++) A_new[j]+=omega*M[j];
    }
  }
  
  content braket=0, norm=0;
  for (int i=0;i<N;i++) {//overlap of the original and updated vector A
    braket+=conj(A[i])*A_new[i];
    norm+=conj(A[i])*A[i];
  }
  
  #if LASER > 3
  fprintf(logg,"\t\t\titerative update of lasing mode amplitude gave the following changes in TM ONO basis:\n\t\t\tindex\told value\t\t\tincrement\t\t\trelative change\n");
  if (ok) TM->diagonalize(false);
  for (int i=0;i<N;i++) {
    content omega=0, dm=0;
    for (int j=0;j<N;j++) {
      omega+=conj(TM->vecL[i*N+j])*A[j];
      dm+=conj(TM->vecL[i*N+j])*(A_new[j]-A[j]);
    }
    double diff=abs(dm/omega);
    fprintf(logg,"\t\t\t%i\t%.5e%+.5ei\t%.5e%+.5ei\t%.5e\n", i, omega.real(), omega.imag(), dm.real(), dm.imag(), diff);
  }
  #endif

  return(braket/norm);
}


content lasing_mode_prototype::perturbative_step(content *A_new, int hmn)
{
  TM->diagonalize(false);
  evaluate_cd(hmn);
  perturbative_amplitude_update(A_new, hmn);
  content braket=0, norm=0;
  for (int j=0;j<TM->N;j++) {
    braket+=conj(A[j])*A_new[j];
    norm+=conj(A[j])*A[j];
  }
  return(braket/norm);
}


void lasing_mode_prototype::amplitude2A()
{
  int N=TM->N;
  content *A_new=A;
  content *aux=LP.aux;
  #if LASER > 3
  A_new=new content[N];
  #endif
  for (int i=0;i<N;i++) {
    for (int k=0;k<LP.dim;k++) aux[k]=CF_basis->psi[CF_basis->pos[i]][k]*LP.n2x[k]*LP.n2;
    A_new[i]=LP.r->braket(aux,0,amplitude,0,1);
  }
  #if LASER > 3
    for (int i=0;i<N;i++) {
    fprintf(logg,"\t\t\t\tindex %i CF_basis position %i old amplitude %.5ei%+.5ei new amplitude %.5ei%+.5ei\n", i, CF_basis->pos[i], A[i].real(), A[i].imag(), A_new[i].real(), A_new[i].imag());
    A[i]=A_new[i];
  }
  delete A_new;
  #endif
  
  #ifdef LASER_CHECK
  content *amp_old=amplitude;
  amplitude=new content[LP.dim];
  update_profile();

  content diff=0, norm=0;
  for (int i=0;i<LP.dim;i++) {
    diff+=amplitude[i]-amp_old[i];
    norm+=amp_old[i]*conj(amp_old[i]);
  }
  if (abs(diff/norm)>1e-8) fprintf(logg,"\t\tWARNING: CF basis redefinition resulted in (unwanted) amplitude change of %e\n", abs(diff/norm));
  delete amp_old;
  #endif
}

void lasing_mode_prototype::omega_check(int i, const char * where)
{
  content omega0;
  double sum=0;
  for (int a=0;a<LP.Nbasis;a++) {
    content omegaa=0;
    for (int b=0;b<LP.Nbasis;b++) omegaa+=conj(TM->vecL[a*LP.Nbasis+b])*A[b];
    if (a==0) omega0=omegaa; else sum+=abs(omegaa);
  }
  if (sum>abs(omega0)*1e-6) {
    fprintf(logg,"WARNING: A into omega check: %s",where);
    fprintf(logg,"lasing mode %i: vector A: omega0=%e%+ei, sum of absolute values of the rest is %e (relative)\n", i, omega0.real(), omega0.imag(), sum/abs(omega0));
  }
}

void lasing_mode_prototype::step2(flag_prototype& flag, int step)
{
  content* A_new=new content[TM->N];
  TM->construct(LP.CF_pool->give_CF_basis(lasing_mode_prototype::k));
  
  flag.iterative=false;
  if (flag.bans>10) flag.iterative=true;   //too many perturbative fails for this state
  if (flag.freq_reset) flag.iterative=true;//frequency was reset, we need to check the state first
  flag.freq_reset=false;		//clear the flag
  if (flag.ban_step!=-1 && step<flag.ban_step+10) flag.iterative=true; //a recent fail
  if (modal_integral/LP.total_modal_integral<C2*(1-Dth/LP.D)) flag.iterative=true;//this mode is zero anyway

  content braket;
  if (flag.iterative) { //simple iterative update;
      TM->diagonalize(true); //diagonalize the threshold matrix without eigenvectors
      braket=iterative_step(A_new);
  } 
  else { //try an perturbative update and check its stability
    for (int i=0;i<LP.dim;i++) profile_save[i]=profile[i];
    for (int i=0;i<TM->N;i++) A_save[i]=A[i];
    TM->diagonalize(false);
    content diff_before=LP.D*TM->eig[0]-1.0;
    braket=perturbative_step(A_new, 1);
    
    //construct the resulting threshold matrix
    for (int i=0;i<TM->N;i++) A[i]=A_new[i];
    update_profile();
    for (int i=0;i<LP.dim;i++) LP.total_profile[i]+=profile[i]-profile_save[i];
    TM_prototype TM_trial(TM->N);
    TM_trial.construct(LP.CF_pool->give_CF_basis(lasing_mode_prototype::k));
    TM_trial.diagonalize(true);
    content diff_after=LP.D*TM_trial.eig[0]-1.0;
    bool passed=true;
    if (abs(diff_after)>abs(diff_before)) passed=false;
    if (flag.updates++>200) passed=false;//!EXPERIMENTAL - too many updates, take iterative
    #if LASER > 2
    fprintf(logg,"\t\tperturbative update tried: diff before %+.3e/%+.3e vs after %+.3e/%+.3e, braket-1 after %+.3e/%+.3e meaning passed %i (updates so far: %i)\n", diff_before.real(), diff_before.imag(), diff_after.real(), diff_after.imag(), braket.real()-1.0, braket.imag(), passed, flag.updates);
    #endif
    
    //decide whether the perturbative step is in a good direction
    if (!passed) {//nope - fold back
      flag.ban_step=step;
      flag.bans++;
      for (int i=0;i<TM->N;i++) A[i]=A_save[i];
      //LP.total_profile[i]-=profile[i]-profile_save[i];
    }
    else TM->copy(TM_trial);//accept the new TM matrix
    //in any case, use the diagonalized TM matrix for iterative update
    braket=iterative_step(A_new);
  }

  //accept the last values
  for (int i=0;i<TM->N;i++) A[i]=A_new[i];
  delete A_new;
  update_profile();
 
  flag.braket=braket;
  flag.diff=LP.D*TM->eig[0]-1.0;
  flag.last_step=step;
  
  #if LASER > 2
  const char* msg[] = {"perturbative", "iterative"};
  fprintf(logg,"\t\tupdating lasing mode (%s: bans=%i, ban_step=%i) at freq %.2e: resulting bracket of real-1/imag %+.3e/%+.3e at diff=%+.3e/%+.3e, rel. modal integral %e\n", msg[flag.iterative], flag.bans, flag.ban_step, LP.x(lasing_mode_prototype::k), braket.real()-1.0, braket.imag(),  flag.diff.real(), flag.diff.imag(),modal_integral/LP.total_modal_integral);
  #endif
}
    


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FLAG PROTOTYPE

void flag_prototype::reset()
{
  bans=updates=0;
  last_step=ban_step=-1;
  diff=braket=0;
  converged=fixed=false;
  iterative=freq_reset=true;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SNAPSHOT PROTOTYPE

snapshot_prototype::snapshot_prototype() 
{
  D=-1;
  Nlm=0;
  ks=Dths=0;
  amplitudes=0;
}

snapshot_prototype::~snapshot_prototype() 
{
  if (ks!=0) {
    delete ks;
    delete Dths;
    delete amplitudes;
  }
}

void snapshot_prototype::save(laser_prototype& laser)
{
  D=LP.D;
  Nlm=laser.Nlm;
  if (ks!=0) {
    delete ks;
    delete Dths;
    delete amplitudes;
  }
  ks=new double[Nlm];
  Dths=new double[Nlm];
  amplitudes=new content[Nlm*LP.dim];
  
  for (int i=0;i<Nlm;i++) {
    ks[i]=laser.lasing_mode[i]->k;
    Dths[i]=laser.lasing_mode[i]->Dth;
    for (int k=0;k<LP.dim;k++) amplitudes[i*LP.dim+k]=laser.lasing_mode[i]->amplitude[k];
  }
}

double snapshot_prototype::load(laser_prototype& laser)
{
      
  fprintf(logg, "loading from snapshot: ");
  if (laser.Nlm>Nlm) {
    fprintf(logg, "more actual (%i) than saved lasing modes (%i): deleting\n", Nlm, laser.Nlm);
    for (int i=Nlm;i<laser.Nlm;i++) delete laser.lasing_mode[i];
    laser.Nlm=Nlm;
  }
  else if (laser.Nlm<Nlm) {
    fprintf(logg, "less actual (%i) than saved lasing modes (%i): creating\n", Nlm, laser.Nlm);
    for (int i=laser.Nlm;i<Nlm;i++) {
      laser.lasing_mode[i]=new lasing_mode_prototype(LP.dim, LP.Nbasis, ks[i]);
      laser.lasing_mode[i]->Dth=Dths[i];
    }
    laser.Nlm=Nlm;
  } else fprintf(logg, "the same number of modes (%i)\n", Nlm);

  for (int i=0;i<Nlm;i++) {
    laser.lasing_mode[i]->k=ks[i];
    laser.lasing_mode[i]->CF_basis=LP.CF_pool->give_CF_basis(ks[i]);
    for (int k=0;k<LP.dim;k++) laser.lasing_mode[i]->amplitude[k]=amplitudes[i*LP.dim+k];
    laser.lasing_mode[i]->amplitude2A();
    laser.lasing_mode[i]->update_profile();
  }
  laser.update_total_profile();
  return(D);
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LASER PROTOTYPE


laser_prototype::laser_prototype(int N)
{
  //LP.Nbasis=LP.Ndiag=laser_prototype::N=N;
  LP.Nbasis=LP.Ndiag=N;
  Nlm=0;
  N_real=0;
  TM_real=new TM_prototype*[N_real];
  real2lm=new int[N_real];
  TM_real_which=new int[N_real];
  TM_eigreal=new content[N];
  TM_k=new double[N];
  lasing_mode=new lasing_mode_prototype*[100];
  snapshot=new snapshot_prototype;
}


laser_prototype::~laser_prototype()
{
  for (int i=0;i<N_real;i++) delete TM_real[i];
  delete TM_real;
  delete TM_k;
  delete TM_eigreal;
  for (int i=0;i<Nlm;i++) delete lasing_mode[i];
  delete lasing_mode;
  destruct_grid();
  delete CF_pool;
  delete snapshot;
}

void laser_prototype::destruct_grid()
{
  delete LP.d;
  delete LP.nx;  
  delete LP.n2x;
  delete LP.total_profile;
  delete LP.r;
  delete LP.arpack;
  delete LP.aux;
}

void laser_prototype::construct_grid()
{
  #if LASER>0
  fprintf(logg,"constructing the grid:\nreading parameters...");
  #endif
  parameters.recompute();
  parameters.list(logg,1+4);

  int dim_x=(int) parameters.read("dim_x","-");
  int dim_y=(int) parameters.read("dim_y","-");
  prec length_x=parameters.read("length_x","nm")*units.nm/units.length;
  prec length_y=parameters.read("length_y","nm")*units.nm/units.length;
  if (dim_y!=2) {printf("only 1D laser implemented but dim_y=%i. Decrease parameter flat_y\n",dim_y); exit(1);}
  
  #if LASER>0
  fprintf(logg,"creating region...");
  #endif
  LP.r=new region(dim_x,-length_x,2*length_x,'n',dim_y,-length_y,2*length_y,'n',inside,1);
  region* r=LP.r;
  r->act_prec_set((region::operators_precision) int(parameters.read("precision","-")));
  r->symmetrize_ops(false, false);
  r->set_inv_angle(0);
  //r->ant(potential, units.length/units.nm, units.energy/units.meV);
  prec hx,hy;
  LP.r->give_par(hx,hy);
  prec l=parameters.read("length_x","nm")-hx*units.length/units.nm;
  printf("length(r)=%e (from points:%e)\n",2*l,(LP.r->Nw-1)*hx*units.length/units.nm);
  printf("length(g)=%e (from points:%e)\n",2*l+hx*units.length/units.nm,(LP.r->Nw)*hx*units.length/units.nm);
  printf("length returned by region:%e\n",LP.r->true_length()*units.length/units.nm);

  #if LASER>0
  fprintf(logg,"creating arpack...");
  #endif
  LP.arpack=new arpack_zn(r->Nw, LP.Ndiag, -1e-12, 1000, true, (char*) "SM");
  
  LP.total_profile=new double[r->Nw];
  LP.d=new content[r->Nw];
  LP.aux=new content[r->Nw];
  LP.dim=r->Nw;
  LP.nx=new content[r->Nw];
  LP.n2x=new content[r->Nw];
}

void laser_prototype::construct_Hamiltonian(const char* filename, bool save)
{
  LP.r->operate(0,LP.d,0,0,region::lin_com,content(1,0), pump_envelope, region::set);

  if (filename==0 || save) {//read data from function
    LP.r->operate(0,LP.nx,0,0,region::lin_com,content(1,0),refractive_index_envelope,region::set);
    for (int i=0;i<LP.r->Nw;i++) LP.n2x[i]=LP.nx[i]*LP.nx[i];
  }
  
  if (filename!=0) {//load or save data
    FILE* file=0;
    if (save) file=fopen(filename,"w");
    if (!save) file=fopen(filename,"r");
    if (file==0) {printf("opening file %s failed\n",filename); exit(1);}
    double x,y, nr, ni;
    content nx;
    for (int i=0;i<LP.r->Nw;i++) {
      LP.r->s2coordinate(i,x,y);
      if (save) {
        nx=LP.nx[i]*sqrt(LP.n2);
        fprintf(file,"%e %e %e\n", x*units.length/units.nm, nx.real(), nx.imag());
      } else {
        fscanf(file,"%le %le %le\n", &x, &nr, &ni);
        nx=content(nr,ni)/sqrt(LP.n2);
        LP.nx[i]=nx;
        LP.n2x[i]=nx*nx;
      }
    }
    fclose(file);
  }

  //check
  FILE* file=fopen("data/laser/refractive_index_check.txt","w");
  if (file==0) {printf("opening file rf_i_check failed\n"); exit(1);}
  double x, y, nr, ni;
  for (int i=0;i<LP.r->Nw;i++) {
    LP.r->s2coordinate(i,x,y);
    content n=LP.nx[i]*sqrt(LP.n2);
    fprintf(file,"%.8e %.8e %.8e\n", x*units.length/units.nm, n.real(), n.imag());
  }
  fclose(file);

  #if LASER>0
  fprintf(logg,"building the hamiltonian\n");
  #endif
  LP.r->operate_tot_init();
  //r->ant(potential, units.length/units.nm, units.energy/units.meV);
  
  content u=-1.0/(LP.n2); //units are: wavevector ^ 2
  LP.r->operate_tot_init(0, 0, reg1::dxdx, u, refractive_indexm2, exp_integral);
  LP.r->operate_tot_init_shift(0);
  LP.shift_H=0;
  //r->ant(potential, units.length/units.nm, units.energy/units.meV);
  LP.r->operate_tot_init_boundary();
  LP.r->operate_tot_init_boundary(0);
  //r->ant(potential, units.length/units.nm, units.energy/units.meV);
}



bool laser_prototype::rescale_CF_grid(double scan_extent_new)
{
  for (int i=0;i<CF_pool->items;i++) CF_pool->undel[i]=0;
#if LASER > 0
  fprintf(logg,"expanding/shrinking the CF_basis grid with scan_step %e: frequency interval scan extent %e with %i terms will be change into ", scan_step, scan_extent, terms);
#endif
  scan_extent=round(2*scan_extent_new/scan_step)*scan_step;
  CF_pool->k1=k1=(LP.omegaa-LP.gamma_perp*scan_extent)/units.c/units.wavevector;
  CF_pool->k2=k2=(LP.omegaa+LP.gamma_perp*scan_extent)/units.c/units.wavevector;
  terms=(int) (2*scan_extent/scan_step)+1;

#if LASER > 0
  fprintf(logg,"extent %e with %i terms\n", scan_extent, terms);
#endif
  
  progress_bar_prototype progress_bar("expanding CF grid");
  progress_bar.start();
  progress_bar.add(terms, 0);
  for (int i=0;i<terms;i++) {
    double k=k1+i*k_step;
#if LASER > 1
    fprintf(logg,"\tCF_basis #%i at frequency %.5e:\n",i, LP.x(k));
#endif
    if (CF_pool->give_CF_basis(k, 1)==0) {progress_bar.finished();return(false);}
    progress_bar.add(0,1);
  }
  progress_bar.finished();
  return(true);
}

//basis_resolution is the limit distance in frequency [in units of gamma_perp] below which a basis is interpolated rather than created (-1 means no interpolation, 0 means FFA)
//scan_extent_init is the initial extent for the frequency range for the TM scans [in units of gamma_perp]
//scan_step the step in frequency [in units of gamma_perp] for TM scans
bool laser_prototype::construct_CF_grid(double scan_extent_init, double scan_step, double basis_resolution, FILE* file)
{
  //if (1/scan_step>1000) {printf("too large/small step (%e) in construct_CF_grid\n",scan_step);exit(1);}
  scan_extent=terms=0;
  
  laser_prototype::scan_step=scan_step;
  k_step=scan_step*LP.gamma_perp/units.c/units.wavevector;
  
  double grain=1/basis_resolution;
  if (basis_resolution==0) grain=0;
  LP.CF_pool=CF_pool=new CF_pool_prototype(grain,100e+6);
  bool success=rescale_CF_grid(scan_extent_init);
  
  //check whether I cover the scan interval
  for (int i=0;i<LP.CF_pool->items;i++) {
    content* km2=LP.CF_pool->CF_basis[i]->km2;
    double km2min=km2[0].real(), km2max=km2[0].real();
    for (int j=1;j<LP.Ndiag;j++) {
      if (km2[j].real()<km2min) km2min=km2[j].real();
      if (km2[j].real()>km2max) km2max=km2[j].real();
    }
    bool ok=false;
    if (km2min<k1*k1 && km2max>k2*k2) ok=true;
    fprintf(logg, "CF basis #%i span check: k2min~%e vs k1^2=%e and k2max~%e vs k2^2=%e : ok=%i\n", i, km2min, k1*k1, km2max, k2*k2, ok);
    if (!ok) printf("WARNING: CF basis %i does not span the scan interval: ratio of limits (should be above 1): lower %e  upper %e: increase number of basis states!\n", i, (k1*k1)/km2min, km2max/(k2*k2) );
  }
 
  if (file==0 || !success) return(success);

  for (int j=terms/2;j<terms/2+1;j++) {
    double k=k1+j*k_step;
    CF_basis_prototype* CF_basis=CF_pool->give_CF_basis(k);
    for (int i=0;i<LP.Ndiag;i++) {
      int p=CF_basis->pos[i];
      double ipr, x0;
      LP.r->ipr_and_com(CF_basis->psi[p], ipr, x0);
      fprintf(file, "%i %i %e %e %e %e %e\n", i, p, CF_basis->km2[p].real(), CF_basis->km2[p].imag(), abs(CF_basis->km2[p]), ipr*units.length/units.nm, x0*units.length/units.nm);
    }
  }
    fflush(0);
  //output the imaginary part of the CF frequencies as a function of the real part for the CF basis in the middle
  /*fprintf(file, "\n\n");
  CF_basis_prototype* CF_basis=CF_pool->give_CF_basis((k1+k2)/2);
  prec length=LP.r->true_length(), k=CF_basis->k, n=sqrt(LP.n2);
  prec scale=units.wavevector*length*units.length;
  fprintf(file, "boundary condition (over i) k=%e, length=%e (scale=%e), product=%e\n", k, length, scale, k*length);
  for (int j=0;j<N;j++) {
    content km=sqrt(CF_basis->km2[CF_basis->pos[j]]);
    if (km.real()<0) km=-km;
    content lhs=tan(n*km*length), rhs=-content(0,1)*n*km/k;
    content res=(lhs-rhs)/(lhs+rhs);
    fprintf(file, "%e %e %e\n", km.real()*scale, km.imag()*scale, abs(res));
    //content a=atan(rhs)/(sqrt(LP.n2)*km);
    //fprintf(file, "eigenvalue %i: implicit equation resuduum %e%+ei (virtual length %e%+ei)\n",res.real(), res.imag(), a.real(), a.imag());    
    //fprintf(file, "eigenvalue %i: implicit equation lhs=%e%+ei, rhs=%e%+ei, relative difference %e%+ei\n",j,lhs.real(), lhs.imag(), rhs.real(), rhs.imag(), res.real(), res.imag());
    content psiN=CF_basis[i]->psi[CF_basis[i]->pos[j]][dim-1], psiN1=CF_basis[i]->psi[CF_basis[i]->pos[j]][dim-2], psiN2=CF_basis[i]->psi[CF_basis[i]->pos[j]][dim-3] ;
    double hx,hy;
    LP.r->give_par(hx,hy);
    content logder=(psiN-psiN1)/(psiN+psiN1)*2.0/hx/content(0,1);
    fprintf(file, "eigenfuncton %i: log der at the right end (over i) %e%+ei\n",j,logder.real(), logder.imag());
    content logder2=(psiN-2.0*psiN1+psiN2)/(hx*hx)/(-psiN1*km*km*n*n);
    fprintf(file, "eigenfuncton %i: log der2 at the right end (over -km2 n2) %e%+ei\n",j,logder2.real(), logder2.imag());
    if (j==0) {
      content* psi=CF_basis[i]->psi[CF_basis[i]->pos[j]];
      for (int k=0;k<dim;k++) {
        content theory=sin((k+1.0)*hx*km*n);
        content ratio=theory/psi[k];
        fprintf(file, "%e %e%+ei %e%+ei %e %e\n", (k+1.0)/dim, theory.real(), theory.imag(), psi[k].real(), psi[k].imag(), ratio.real(), ratio.imag());
      }
    }
  }
*/
  return(success);
}


void laser_prototype::scan_CF(double R1, double R2, int terms, FILE* file)
{
  if (terms<2 || terms>1000) {printf("too little/many terms (%i) in construct_CF_grid\n",terms);exit(1);}
  double R_step=(R2-R1)/(terms-1);

#if LASER > 0
  fprintf(logg,"scanning CF bases (%i terms in the length interval <%e,%e> resulting in step %e)\n",terms, R1, R2, R_step);
#endif

  progress_bar_prototype progress_bar("CF scan");
  progress_bar.start();
  progress_bar.add(terms, 0);
  for (int i=0;i<terms;i++) {
    double R=R1+i*R_step;
#if LASER > 1
    fprintf(logg,"\tCF_basis #%i at length %.5e:\n",i, R);
#endif
    
    parameters.set("flat_x", R, "nm");
    parameters.set("lv", R/50, "nm");    
    parameters.recompute();
    construct_grid();    
    CF_basis_prototype* CF_basis=LP.CF_pool->give_CF_basis(LP.omegaa/units.c/units.wavevector);
    
    if (file!=0) {
      //prec scale=units.wavevector*units.nm;
      prec scale=units.wavevector/(2*M_PI/(R*units.nm));
      fprintf(file, "%e ",R);
      for (int j=0;j<LP.Ndiag;j++) fprintf(file, "%e %e ", (CF_basis->km2[CF_basis->pos[j]]).real()*scale*scale, (CF_basis->km2[CF_basis->pos[j]]).imag()*scale*scale);
      fprintf(file, "\n");
    }
    progress_bar.add(0,1);
  }
  progress_bar.finished();
}


void laser_prototype::update_total_profile() 
{
  double tot=0;
  for (int i=0;i<LP.dim;i++) {
    LP.total_profile[i]=1;
    for (int j=0;j<Nlm;j++) LP.total_profile[i]+=lasing_mode[j]->profile[i];
    tot+=LP.total_profile[i]-1;
  }
  LP.total_modal_integral=0;
  for (int i=0;i<Nlm;i++) {
    lasing_mode[i]->update_modal_power();
    LP.total_modal_integral+=lasing_mode[i]->modal_integral;
  }
  #if LASER > 2
  fprintf(logg,"\t\ttotal profile updated, total modal integral=%e\n\t\trelative modal integrals:", LP.total_modal_integral);
  for (int i=0;i<Nlm;i++) fprintf(logg,"#%i=%.2e ", i, lasing_mode[i]->modal_integral/LP.total_modal_integral);
  fprintf(logg,"\n");
  #endif
}


void laser_prototype::scan_TM(int hmn, FILE* file, FILE* file2) 
{
  #if LASER > 0
  fprintf(logg,"starting TM scan at pumping %e with %i lasing modes\n", LP.D, Nlm);
  #endif

  int N_sorting=0, N_zeros=0;
  simple_sorting_prototype sorting_aux(LP.Nbasis, 3);
  simple_sorting_prototype* sorting[1000];
  double errors[6];
  TM_prototype TM(LP.Nbasis);
  if (LP.FFA) TM.construct(LP.CF_pool->give_CF_basis(k1), true);

  progress_bar_prototype progress_bar("scan TM");
  progress_bar.start();
  progress_bar.add(terms, 0);
  //STEP 1: identify all zero crossings within the grid
  double k=k1, k_step_dyn=k_step;
  int i=0;
  do {
    if (LP.FFA) TM.update_from_Tau(k); else TM.construct(CF_pool->give_CF_basis(k));
#if SORTING > 0
      fprintf(logg,"\t\tTM scan sort: starting at freq %e with k1=%e, k2=%e, k_step=%e\n",LP.x(k), k1,k2,k_step);
#endif
    TM.diagonalize(LP.TMsorting);
    errors[0]=min(LP.Nbasis,max(10,2*Nlm));
    sorting_aux.sort_mother(LP.TMsorting+4,TM.vecL,TM.vecR,TM.eig, k, errors, TM.eig_der);
#if LASER>2
    fprintf(logg,"\t\tTM scan sort: at frequency %e the sort_mother step 1 returned error %e with goal %e and maximum tolerable %e, at freq. step to minimal step ratio %e\n", k, errors[LP.TMsorting], C7a, C7b, k_step_dyn/k_step);
#endif
    if (errors[LP.TMsorting]!=0) {//if the error is unknown, continue with the previous step
      //if the error is not acceptable, restart
      if (errors[LP.TMsorting]>C7b && k_step_dyn>k_step) {
        k-=k_step_dyn-k_step;//go almost back
        k_step_dyn=k_step;
        continue;
      } 
      //otherwise accept and recalculate the step in frequency to achieve the error
      k_step_dyn/=errors[LP.TMsorting]/C7a;
      if (k_step_dyn<k_step) k_step_dyn=k_step;
      if (k_step_dyn>k_step*C8) k_step_dyn=k_step*C8;
    }
    sorting_aux.sort_mother(LP.TMsorting+2,TM.vecL,TM.vecR,TM.eig, k, errors, TM.eig_der);

    if (file!=0) {
      fprintf(file, "%e ",LP.x(k));
      for (int j=0;j<LP.Nbasis;j++) {
        content val=sorting_aux.give_value(j);
        fprintf(file, "%e %e ",val.real(), val.imag());
      }
      fprintf(file, "%e %e %e\n", errors[0], errors[1], errors[2]);
    }
    int inc=sorting_aux.zeros_crossed();
    N_zeros+=inc;  
    #if LASER > 2
    if (inc>0) fprintf(logg,"\t\t%i new TM real eigenvalues at frequency %e (total %i now)\n", inc, LP.x(k), N_zeros);
    #endif
    //if zero(s) were crossed, make a copy of the sorting class - for each zero one copy
    while (inc-->0) {
      sorting[N_sorting]=new simple_sorting_prototype(LP.Nbasis, 3);
      sorting[N_sorting]->copy(&sorting_aux);
      if (N_sorting++==1000) {printf("too many crossings identified in scan_TM, probably wrong completely\n"); exit(0);}
    }
    progress_bar.add(0, 1);
    k+=k_step_dyn;
  } while (k<k2);
  progress_bar.finished();
  if (file!=0) fprintf(file,"\n");

  //STEP 2: identify state which crosses zero and its eigenvalue
  if (N_zeros==0) {printf("no crossings of zero found in scan_TM\n"); exit(1);}
  int pntr=0;
  int* is=new int[N_zeros*3];
  int* js=is+N_zeros, *pos=is+N_zeros*2;
  double* k_est=new double[N_zeros*2];
  double* v_est=k_est+N_zeros, v;
  for (int i=0;i<N_sorting;i++) {
    int pntr_sorting=0;
    for (int j=0;j<LP.Nbasis;j++) {
      if (!sorting[i]->zero_estimate(j, k, v )) continue;
      pos[pntr]=pntr;
      js[pntr]=j;
      is[pntr]=i;
      k_est[pntr]=k;
      v_est[pntr]=v;
      pntr++;
      if (pntr_sorting++>0) i++;//if this sorting has more than one zero crossings, move to the next copy
      #if LASER > 2
      fprintf(logg,"\t\t%i-th TM real eigenvalue:  %i. in sorting %i at frequecy is %e with value estimate %e\n", pntr-1, js[pntr-1], is[pntr-1], LP.x(k_est[pntr-1]), v_est[pntr-1]);
      #endif
    }
  }
  if (pntr!=N_zeros) {printf("pntr (%i) is not equal N_zeros (%i) in scan_TM (N_sorting=%i)\n",pntr, N_zeros, N_sorting);exit(1);}
  //sort crossed zeros according to estimated real part of the eigenvalue
  for (int i=0;i<N_zeros;i++) 
    for (int j=0;j<N_zeros-i-1;j++) if (v_est[pos[j]]<v_est[pos[j+1]]) { int pos_aux=pos[j]; pos[j]=pos[j+1]; pos[j+1]=pos_aux; }
  //de+allocate space
  for (int i=0;i<N_real;i++) delete TM_real[i];
  delete TM_real;
  delete TM_real_which;
  TM_real=new TM_prototype*[N_zeros];
  TM_real_which=new int[N_zeros];

  //STEP 3: run the bracketing until hmn largest values found
  progress_bar.reset("scan TM: bracketing");
  progress_bar.add(N_zeros,0);
  progress_bar.start();
  N_real=0;
  for (int i=0;i<N_zeros;i++) {

    double k;
    content v;
    int label=pos[i];  
    bool success=identify_root2(js[label], sorting[is[label]], k, v, k_step/2);
    
    #if LASER > 2
    fprintf(logg,"\t\t%i-th TM real eigenvalue searched for(success:%i): difference to estimated frequency %e and value %.5e is %e and %.5e%+.5ei\n", N_real, success, LP.x(k_est[label]), v_est[label], 1-k/k_est[label], (v-v_est[label]).real(), (v-v_est[label]).imag());
    #endif
    if (!success) continue;
    TM_real[N_real]=new TM_prototype(LP.Nbasis);
    TM_real[N_real]->construct(CF_pool->give_CF_basis(k));
    TM_real[N_real]->diagonalize(false);
    TM_real_which[N_real]=sorting[is[label]]->position[js[label]];
    #if LASER > 2
    fprintf(logg,"\t\t%i-th TM real matrix: CF frequency %e and %i-th eigenvalue %e vs k=%e v=%e\n", N_real, k, TM_real_which[N_real],  TM_real[N_real]->eig[TM_real_which[N_real]].real(), k, v.real());
    #endif
    TM_eigreal[N_real]=v;
    TM_k[N_real]=k;
    progress_bar.add(0,1);
    if (++N_real==hmn) break;
  }
  if (N_real<hmn) fprintf(logg,"WARNING: not enough real eigenvalues found (%i vs %i asked for)\n",N_real,hmn);
  progress_bar.finished();

  delete is;
  delete k_est;
  for (int i=0;i<N_sorting;i++) delete sorting[i];
  //for (int i=0;i<Nlm;i++) lasing_mode[i]->omega_check(i,"end of scan_TM\n");
  
  //output of normalized amplitudes of real frequency eigenmodes
  if (file2!=0) {
    double* norm=new double[hmn];
    for (int i=0;i<hmn;i++) {
      norm[i]=0;
      if (i<Nlm) {norm[i]=LP.D;continue;}
      if (i<N_real) {
        CF_basis_prototype* CF_basis=TM_real[i]->CF_basis;
        for (int k=0;k<LP.dim;k++) {
          content val=0;
          for (int j=0;j<LP.Nbasis;j++) val+=CF_basis->psi[CF_basis->pos[j]][k]*TM_real[i]->vecR[TM_real_which[i]*LP.Nbasis+j];
          norm[i]+=(val*conj(val)).real();
        }
        continue;
      }
      norm[i]=1;
    }
    
    for (int k=0;k<LP.dim;k++) {
      prec x,y;
      LP.r->s2coordinate(k,x,y);
      fprintf(file2, "%e %e ",x*units.length/units.nm, LP.D*LP.d[k].real()/LP.total_profile[k]);
      for (int i=0;i<hmn;i++)  {
        content val=0;
        if (i<Nlm) val=lasing_mode[i]->amplitude[k];
        else if (i<N_real) {
          CF_basis_prototype* CF_basis=TM_real[i]->CF_basis;
          for (int j=0;j<LP.Nbasis;j++) val+=CF_basis->psi[CF_basis->pos[j]][k]*TM_real[i]->vecR[TM_real_which[i]*LP.Nbasis+j];
        }
        else val=0;
        fprintf(file2, "%e %e %e ", val.real(), val.imag(), pow(abs(val),2)/norm[i]);
      }
      fprintf(file2, "\n");
    }
    delete norm;
  }
}

bool laser_prototype::identify_root2(int i, simple_sorting_prototype * sorting, double &k, content& v, double k_last_update, int history)
{
  if (k_last_update<0) k_last_update*=-1;
  bool success=false; 
  content vleft=0, vright=0, vsave=sorting->give_value(i, 0, &k);
  double kleft, kright, k_est, ksave=k;
  
  if (k_last_update<1e-6) k_last_update=1e-6;
  //if k_last_update estimate is provided, use it, if not, use grid points
  if (k_last_update!=0) {
    kleft=k-abs(k_last_update);
    kright=k+abs(k_last_update);
  } else {
    kleft=k1+k_step*floor((k-k1)/k_step);
    kright=kleft+k_step;
    k_last_update=k_step;
  }
  
  //calculate eigenvalues in the neighborhood
  TM_prototype TM(LP.Nbasis);
  if (LP.FFA) TM.construct(LP.CF_pool->give_CF_basis(k), true);
  for (int j=0;j<2;j++) {
    if (j==0) k=kleft; else k=kright;  
    if (LP.FFA) TM.update_from_Tau(k); else TM.construct(CF_pool->give_CF_basis(k));
    TM.diagonalize(LP.TMsorting);
    sorting->sort_mother(LP.TMsorting,TM.vecL,TM.vecR, TM.eig, k, 0, TM.eig_der);
    v=sorting->give_value(i, 0, &k);
    if (j==0) vleft=v; else vright=v;
  }
  #if LASER > 2
  fprintf(logg,"IL3 started with: i=%i, k=%e, v=%.5e%+.5ei, k_last_update=%e, history=%i\nresulting in kleft=%e (vleft=%.5e%+.5ei) and kright=%e (vright=%.5e%+.5ei)\n", i, LP.x(ksave), vsave.real(), vsave.imag(), k_last_update, history, LP.x(kleft), vleft.real(), vleft.imag(), LP.x(kright), vright.real(), vright.imag() );
  #endif
  
  int step=0;
  do {  
    if (!sorting->extrapolate(i,k_est,0, history)) {
        fprintf(logg,"WARNING: do not see any zero crossing in identifying limits3 (history=%i)\n", history);
        break;
    }
    
    if ((k_est-k)>k_last_update*3) k+=(k_last_update*=3);
    else if ((k_est-k)<-k_last_update*3) k-=(k_last_update*=3);
    else k=k_est;
    if (k<k1 || k>k2) {fprintf(logg,"WARNING: gone out of frequency window <%e,%e> in identifying limits3 k=%e\n", LP.x(k1), LP.x(k2), LP.x(k));break;}
 
    if (LP.FFA) TM.update_from_Tau(k); else TM.construct(CF_pool->give_CF_basis(k));
    TM.diagonalize(LP.TMsorting);
    sorting->sort_mother(LP.TMsorting,TM.vecL,TM.vecR,TM.eig, k, 0, TM.eig_der, i*(history-2) + (history-3));
    v=sorting->give_value(i);
    if (abs(v.imag()/v.real())<C6) {success=true;break;}
    
    #if LASER > 2
    fprintf(logg,"\t\t\tidentifying limits3 step %i: at frequency k=%e (est. from %i with last update=%e), closest eigenvalue to (%.5e,%+.5e) is (%.5e,%+.5e)\n", step, LP.x(k), history, k_last_update, vsave.real(), vsave.imag(), v.real(), v.imag());
    #endif
  } while (step++<100);
  #if LASER > 2
  fprintf(logg,"\t\tidentifying limits3 converged:%i after %i steps with the estimated zero at %e with eigenvalue (%.5e%+.5e)\n", success, step, LP.x(k), v.real(), v.imag());
  #endif
  if (!success) {k=ksave;v=vsave;}
  return(success);
}

bool laser_prototype::update_pump(double &D, double Dnext, double Dnextnext)
{
  D=LP.D;
  if (Dnext<0) {fprintf(logg,"negative pump value asked for (%e) in update_pump\n", Dnext);exit(1);}
  #if LASER > 1
  fprintf(logg,"\tsetting a new value for the pump: actual=%.5e, next mode at %.5e next next mode at %.5e\n", D, Dnext, Dnextnext);
  #endif

  if (Nlm==0) {
    D=min(Dnext*(1+C3*LP.Dinc), Dnextnext);
    #if LASER > 0
    fprintf(logg,"first lasing mode: setting the pump at %.3e = (1 + %.2e) / largest eigenvalue of TM\n", D, C3*LP.Dinc);
    #endif
    return(true);
  }
  
  if (Dnext<D) {
    #if LASER > 0
    fprintf(logg,"the pump converged: the pump exceeded the goal (%.5e vs %.5e), actual value will be kept\n", D, Dnext);
    #endif
    return(true);
  }
  D=min(D*(1+C3*LP.Dinc/sqrt(Nlm)),Dnextnext);
  #if LASER > 1
  fprintf(logg,"\tthe pump has not converged, will be enlarged to %.5e\n", D);
  #endif
  return(false);
}


int laser_prototype::next_lasing_mode_pump(double& knext, double& Dnext, double& Dnextnext, FILE* file)
{
  #if LASER > 1
  fprintf(logg,"\tchecking the real TM eigenvalues (%i found at %i lasing modes)\n", N_real, Nlm);
  #endif
  if (N_real<=Nlm) {fprintf(logg, "number of real TM eigenvalues (%i) not more than lasing modes (%i)\n", N_real,Nlm); exit(1);}
  //do all lasing modes correspond to a real eigenvalue of 1/D?
  delete real2lm;
  real2lm=new int[N_real];
  for (int i=0;i<N_real;i++) real2lm[i]=-1;
  bool* identified=new bool[N_real];
  for (int i=0;i<N_real;i++) identified[i]=false;
  for (int i=0;i<Nlm;i++) {
    double k=lasing_mode[i]->k;
    double ext=0;
    int which=-1;
    for (int j=0;j<N_real;j++) {//!EXPERIMENTAL - TRYING TO FIGHT REPLICAS
      //if (identified[j]) continue;
      double aux=abs(k-TM_k[j]);
      if (which==-1 || aux<ext) {ext=aux;which=j;}
    }
    bool found=(ext<C1*abs(k));
    #if LASER > 2
    fprintf(logg,"\t\t\t%i lasing mode at frequency %.5e: real mode %i (with D*lambda-1=%.3e%+.3ei) is closest in frequency (being %.5e) at rel. freq. distance %.5e: meaning was found=%i\n", i, LP.x(k), which, LP.D*TM_eigreal[which].real()-1.0, TM_eigreal[which].imag(),  LP.x(TM_k[which]), ext, found);
    #endif
    if (found) {real2lm[which]=i; identified[which]=true;}
    else if (lasing_mode[i]->modal_integral>C2*LP.total_modal_integral*(1-lasing_mode[i]->Dth/LP.D)) {
      fprintf(logg,"%i lasing mode does not correspond to a TM eigenvalue [ext=%e vs C1*abs(k)=%e]\n", i, ext, C1*abs(k));
      exit(1);
    }
    //check whether a lasing mode can be discarded
    if (lasing_mode[i]->modal_integral<C2*LP.total_modal_integral*(1-lasing_mode[i]->Dth/LP.D)) {
      #if LASER > 1
      fprintf(logg,"\t%i lasing mode died off, will be descarded (Dth/D=%e, relative_modal_integral=%e)\n", i, lasing_mode[i]->Dth/LP.D, lasing_mode[i]->modal_integral/LP.total_modal_integral);
      #endif
      delete lasing_mode[i];
      for (int j=i;j<Nlm-1;j++) lasing_mode[j]=lasing_mode[j+1];
      Nlm--;
      for (int j=0;j<N_real;j++) if (real2lm[j]==i) real2lm[j]=-1;
      i--;
      //printf("deleting lasing mode: stopping\n");
      //exit(0);
    } 
  }
  double* v=new double[2*N_real];
  double* k=v+N_real;
  int* new2real=new int[N_real];
  int pntr=0;
  for (int i=0;i<N_real;i++) {
    if (identified[i]) continue;
    v[pntr]=TM_eigreal[i].real();
    new2real[pntr]=i;
    k[pntr++]=TM_k[i];
  }
  int* index=new int[pntr];
  picsrt_index(pntr, v, index);
  sort_from_index(pntr, v, index);
  sort_from_index(pntr, k, index);
  int i2return=new2real[index[0]];
  delete index;
  //picsrt_inplace(pntr, v, k);//!!!!! EXPERIMENTAL
  if (file!=0) {
    fprintf(file, "%e ", LP.D);
    for (int i=0;i<Nlm;i++) fprintf(file, "%e %e ", LP.x(lasing_mode[i]->k), 1/abs(lasing_mode[i]->TM->eig[0]));
    for (int i=0;i<pntr;i++) fprintf(file, "%e %e ", LP.x(k[i]), 1/v[i]);
    for (int i=pntr;i<8-Nlm;i++) fprintf(file, "0 0 ");
    fprintf(file, "\n");
  }
  if (pntr<2) {fprintf(logg,"no two unassigned real values of TM matrix were found: pntr=%i\n", pntr); exit(1);}
  #if LASER > 1
  fprintf(logg,"\tthe sorted unassigned real TM eigenvalues:\n");
  for (int i=0;i<pntr; i++) fprintf(logg,"\t%i-th value %.5e frequency %e (grid spot %i) corresponding to pump %.5e\n", i, v[i], LP.x(k[i]), (int) round((k[i]-k1)/k_step), 1/v[i]);
  #endif
  Dnext=1/v[0];
  Dnextnext=1/v[1];
  knext=k[0];
  delete identified;
  delete v;
  //for (int i=0;i<Nlm;i++) lasing_mode[i]->omega_check(i, "end of next lasing mode pump\n");
  delete new2real;
  return(i2return);
}

content laser_prototype::reset_mode_frequency(int mode)
{
  simple_sorting_prototype sorting(LP.Nbasis,3);
  TM_prototype& TM=*lasing_mode[mode]->TM;
  TM.diagonalize(LP.TMsorting);
  sorting.sort_mother(LP.TMsorting,TM.vecL,TM.vecR,TM.eig, lasing_mode[mode]->k);
  
  //keep the old values in case the sorting fails
  double knew=lasing_mode[mode]->k;
  content v=LP.D*lasing_mode[mode]->TM->eig[0]-1.0;
  //if (identify_limits(0,&sorting)) identify_root(&sorting,0,knew,v);
  //if (identify_limits2(0,&sorting,lasing_mode[mode]->k_last_update)) identify_root(&sorting,0,knew,v);//!EXPERIMENTAL
  if (!identify_root2(0,&sorting, knew, v, lasing_mode[mode]->k_last_update,3)) identify_root2(0,&sorting, knew, v, lasing_mode[mode]->k_last_update,2);//!EXPERIMENTAL
  
  lasing_mode[mode]->CF_basis=CF_pool->give_CF_basis(knew, 1);
  double kold=lasing_mode[mode]->k;
  lasing_mode[mode]->k=knew;
  lasing_mode[mode]->amplitude2A();
  double kalt=lasing_mode[mode]->update_frequency();
#if LASER > 2
  fprintf(logg,"\t\tfrequency update: from %e by reset %e vs update %e\n", kold, knew-kold, kalt-kold);
#endif
  lasing_mode[mode]->k_last_update=knew-kold;
  //for (int i=0;i<Nlm;i++) lasing_mode[i]->omega_check(i, "end of reset mode frequency\n");
  return(LP.D*v-1.0);
}

int laser_prototype::iterate_at_pump(flag_prototype* flag, int step)
{
  int which=-1;

  //run all in the beginning to estimate diff, or if some is fixed, iterate only it
  for (int i=0;i<Nlm;i++) if (flag[i].last_step==-1 || flag[i].fixed) {which=i;break;}
  
  if (which==-1) {
    for (int i=0;i<Nlm;i++) {
      if (flag[i].converged) continue;
      if (which==-1) which=i;
      #if LASER > 2
      fprintf(logg,"\t\titeration at pump: considering mode %i (freq. %e) diff=%.3e%+.3ei last_step=%i\n", i, LP.x(lasing_mode[i]->k), flag[i].diff.real(), flag[i].diff.imag(), flag[i].last_step);
      #endif
      switch (LP.strategy) {
	case 0 : {//take the one with largest error
	  if (abs(flag[i].diff)>abs(flag[which].diff)) which=i;
	  break;
	}
	case 2 : {//take the one with largest error times the state weight
	  if (abs(flag[i].diff)*lasing_mode[i]->modal_integral>abs(flag[which].diff)*lasing_mode[which]->modal_integral) which=i;
	  break;
	}
	case 1 : {//fix the one with largest norm
	  if (lasing_mode[i]->modal_integral>lasing_mode[which]->modal_integral) which=i;
	  break;
	}
      }
      //from time to time iterate each, or after freqency was reset
      if (LP.strategy!=1) if (step-flag[i].last_step>Nlm*2 || flag[i].freq_reset) {which=i;break;}
    }
    if (LP.strategy==1) flag[which].fixed=true;
  }
  
  #if LASER > 2
  fprintf(logg,"\t\titeration at pump: chosen mode %i (with abs(diff)=%e)\n", which,abs(flag[which].diff));
  #endif

  //take the step
  lasing_mode[which]->step2(flag[which], step);
  
  //decide if frequency is needed to reset
  if (abs(flag[which].diff.imag())>abs(flag[which].diff.real())) {
    update_total_profile();//!EXPERIMENTAL
    flag[which].diff=reset_mode_frequency(which);
    flag[which].freq_reset=true;
  }
  
  bool converged=false;
  //if the mode is irrelevant (too small and not new), declare convergence
  if (lasing_mode[which]->modal_integral<C2*LP.total_modal_integral*(1-lasing_mode[which]->Dth/LP.D)) converged=true;

  //if the error is below threshold, declare convergence
  if (abs(flag[which].diff)<C5) converged=true;

  //if we converged, release the flags
  if (converged) {
    flag[which].converged=true;
    flag[which].fixed=false;
    flag[which].freq_reset=false;
  } //if not, reset all modes as not converged
  else for (int i=0;i<Nlm;i++) flag[i].converged=false;

  #if LASER > 2
  fprintf(logg,"\t\tresulting diff (real/imag)=%+.3e/%+.3e frequency reset: %i gives converged=%i\n", flag[which].diff.real(), flag[which].diff.imag(), flag[which].freq_reset,  flag[which].converged);
  #endif
  return(which);
}

//calculates the unnormalized overlap between two modes, two types
double overlap_aux(int dim, content* a, content *b, int type)
{
  content sum=0;
  if (type==0) for (int k=0;k<dim;k++) sum+=abs(a[k])*abs(b[k]);
  else if (type==1) for (int k=0;k<dim;k++) sum+=a[k]*b[k]*LP.n2x[k];       
  return((sum*conj(sum)).real());
  return(0);
}

int laser_prototype::select_one_mode()
{
  if (Nlm==0) return(-1);
  int mode=0;
  double wmin=1;
  for (int i=0;i<Nlm;i++) {
    double weight=lasing_mode[mode]->modal_integral/(LP.total_modal_integral*(1-lasing_mode[mode]->Dth/LP.D));
    if (weight<wmin || i==0) {wmin=weight; mode=i;}
  }
  //if (wmin>1e-3) mode=-1;
  fprintf(logg,"a single mode selected for convergence: %i\n",mode);
  return(mode);
}

void laser_prototype::converge_one_mode(int mode)
{
  if (mode==-1) return;
  fprintf(logg,"a single mode convergence started\n");
  int step=0;
  flag_prototype flag;
  flag.reset();
  do {
    if (step++<20) flag.freq_reset=true;
    lasing_mode[mode]->step2(flag, step);
    lasing_mode[mode]->update_profile();
    update_total_profile();
    if (abs(flag.diff.imag())>abs(flag.diff.real())) {
      flag.diff=reset_mode_frequency(mode);
      flag.freq_reset=true;
    }
    if (abs(flag.diff)<C5 || lasing_mode[mode]->modal_integral<C2*LP.total_modal_integral*(1-lasing_mode[mode]->Dth/LP.D)) break;
  } while (step<1000);
  fprintf(logg,"a single mode convergence finished after %i steps with diff=%.3e%+.3e\n", step, flag.diff.real(), flag.diff.imag());
}

void laser_prototype::plot_lasing_modes(FILE * file)
{
  for (int k=0;k<LP.dim;k++) {
    prec x,y;
    LP.r->s2coordinate(k,x,y);
    fprintf(file, "%e %e ",x*units.length/units.nm, (LP.n2x[k]*LP.n2).real());
    for (int m=0;m<Nlm;m++) {
      content aux=lasing_mode[m]->amplitude[k]/lasing_mode[m]->norm;
      fprintf(file, "%e %e %e ", aux.real(), aux.imag(), (aux*conj(aux)).real());
    }
    fprintf(file, "\n");
  }
}

bool laser_prototype::check_replicas(bool del)
{
  bool replicas=false;
  double* norm=new double[Nlm]; 
  for (int m=0;m<Nlm;m++) norm[m]=sqrt(overlap_aux(LP.dim,lasing_mode[m]->amplitude,lasing_mode[m]->amplitude, 0));
  for (int i=0;i<Nlm;i++) {
    for (int j=i+1;j<Nlm;j++) {
      prec dk=abs(lasing_mode[i]->k-lasing_mode[j]->k)/(lasing_mode[i]->k+lasing_mode[j]->k);
      prec O=overlap_aux(LP.dim,lasing_mode[i]->amplitude,lasing_mode[j]->amplitude, 0)/(norm[i]*norm[j]);
      fprintf(logg,"checking for replicas: lasing mode pair %i-%i has frequency distance %e and overlap %e\n", i, j, dk, O);
      if ((dk<1e-4 && abs(O-1)<1e-3) || dk<1e-8) {//I have found a replica
        replicas=true;
        fprintf(logg,"WARNING: I found a replica\n");
	if (!del) continue;
        fprintf(logg,"state %i will be added to state %i and deleted after that (norm[%i]=%e, norm[%i]=%e\n", j, i, j, norm[j], i, norm[i]);
	if (norm[j]<norm[i]) {
	  double fac=sqrt(1+norm[j]*norm[j]/(norm[i]*norm[i]));
	  for (int k=0;k<LP.dim;k++) lasing_mode[i]->amplitude[k]=lasing_mode[i]->amplitude[k]*fac;
	} else {
	  double fac=sqrt(1+norm[i]*norm[i]/(norm[j]*norm[j]));
	  for (int k=0;k<LP.dim;k++) lasing_mode[i]->amplitude[k]=lasing_mode[j]->amplitude[k]*fac;
	}
	//if (norm[i]<norm[j]) for (int k=0;k<LP.dim;k++) lasing_mode[i]->amplitude[k]=lasing_mode[j]->amplitude[k];
	delete lasing_mode[j];
	for (int k=j;k<Nlm-1;k++) lasing_mode[k]=lasing_mode[k+1];
	Nlm--;
	j--;
        lasing_mode[i]->amplitude2A();
        update_total_profile();

	//converge the merged mode first by itself
	converge_one_mode(i);
      }
    }
  }
  delete norm;
  if (replicas && del) update_total_profile();
  return(replicas);
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TLM_APPROX

TLM_approx_prototype::TLM_approx_prototype(laser_prototype& laser)
{
  TLM_approx_prototype::laser=&laser;
  Nl=laser.Nlm;
  Nr=laser.N_real;
  if (Nl<0 || Nr<0) {printf("wrong parameters in TLM_approx constructor Nl=%i, Nr=%i\n",Nl,Nr); exit(0);} 
  I=new double[Nr];
  dI=new content[Nr];
  lambda=new content[Nr];
  xi=new content[Nr*Nr];
  xi_inv=new content[Nr*Nr];
  Dth=new double[Nr];
  ipr=new double[Nr];
  x0=new double[Nr];
  real2lm=new int[Nr];
  lm2real=new int[Nr];
  o[0]=new double[Nr*2];
  o[1]=o[0]+Nr;
  oi[0]=new double[Nr*2];
  oi[1]=oi[0]+Nr;
  #if LASER > 0
  fprintf(logg,"TLM_approx starts with %i lasing modes and %i real TM eigenvalues at pump %e\n", Nl, Nr, LP.D);
  #endif
}


TLM_approx_prototype::~TLM_approx_prototype()
{
  delete lambda;
  delete xi;
  delete xi_inv;
  delete I;
  delete dI;
  delete Dth;
  delete ipr;
  delete x0;
  delete real2lm;
  delete lm2real;
  delete o[0];
  delete oi[0];
}

void TLM_approx_prototype::initialize()
{
  //fill in the input values - calculate the auxiliary amplitudes psi.aL, psi.aR, psi.Lambda^-1*aR for both lasing and non-lasing states
  int dim=LP.dim;
  content* PsiaL=new content[dim*Nr];
  content* PsiaR=new content[dim*Nr];
  content* PsiAaR=new content[dim*Nr];
  for (int j=0;j<dim*Nr;j++) PsiaL[j]=PsiaR[j]=PsiAaR[j]=0;
  double* k=new double[Nr];
  for (int i=0;i<Nl;i++) lm2real[i]=-1;
  
  int pntr=0;
  for (int m=0;m<Nr;m++) {
    bool lasing=false;
    if (laser->real2lm[m]!=-1) lasing=true;
    int which;
    TM_prototype* TM;
    content* A, omega;
    
    if (lasing) {
      TM=laser->lasing_mode[laser->real2lm[m]]->TM;
      TM->diagonalize(false);
      which=0;
      laser->lasing_mode[laser->real2lm[m]]->amplitude2A();//!trying to get ipr/x0 the same in TLM and direct call in main_laser
      A=laser->lasing_mode[laser->real2lm[m]]->A;
      k[m]=laser->lasing_mode[laser->real2lm[m]]->k;
    }
    else {
      TM=laser->TM_real[m];
      which=laser->TM_real_which[m];
      A=TM->vecR+which*LP.Nbasis;
      k[m]=laser->TM_k[m];
    }
    lambda[m]=TM->eig[which];
    content* aL=TM->vecL+which*LP.Nbasis;
    content* aR=TM->vecR+which*LP.Nbasis;
    
    for (int i=0;i<LP.Nbasis;i++) {
      content* psi=TM->CF_basis->psi[TM->CF_basis->pos[i]];
      content km2=TM->CF_basis->km2[TM->CF_basis->pos[i]];
      for (int j=0;j<dim;j++) {
        PsiaL[m*dim+j]+=conj(aL[i])*psi[j];
        PsiaR[m*dim+j]+=aR[i]*psi[j];
        PsiAaR[m*dim+j]+=aR[i]*psi[j]/LP.Lambda(k[m],km2);
      }
    }
#ifdef LASER_CHECK
    content res=0, res2=0;//check that A is the TM eigenvector with eigenvalue lambda
    for (int a=0;a<LP.Nbasis;a++) {
      res+=lambda[m]*A[a];
      res2+=TM->eig[which]*TM->vecR[which*LP.Nbasis+a];
      for (int b=0;b<LP.Nbasis;b++) res-=TM->M[a*LP.Nbasis+b]*A[b];
      for (int b=0;b<LP.Nbasis;b++) res2-=TM->M[a*LP.Nbasis+b]*TM->vecR[which*LP.Nbasis+b];
    }
    if (abs(res)>1e-6*LP.Nbasis*LP.Nbasis) fprintf(logg,"WARNING: |lambda A - TM*A|=%e (alt=%e) for state %i (real index)\n", abs(res), abs(res2), m);
#endif
    //calculate ipr, com
    LP.r->ipr_and_com(PsiaR+m*dim, ipr[m], x0[m]);
    ipr[m]*=units.length/units.nm;
    x0[m]*=units.length/units.nm;

    dI[m]=I[m]=0;
    if (lasing) {//this mode is already lasing at the beginning
      real2lm[m]=pntr;
      Dth[m]=laser->lasing_mode[laser->real2lm[m]]->Dth;
      lm2real[pntr]=m;
      omega=0;
      for (int b=0;b<LP.Nbasis;b++) omega+=conj(aL[b])*A[b];
      I[m]=(omega*conj(omega)).real();
      if (++pntr>=Nr) {printf("pntr %i reached Nr %i in TLM_approx initialize\n", pntr, Nr); exit(1);}
    } else {                    //this mode is not yet lasing
      Dth[m]=-1;
      real2lm[m]=-1;
    }
#if LASER > 1
    fprintf(logg, "\tTLM_approx %i-th real mode eigenvalue %e%+ei (vs 1/pump=%e) %i-th in the corresponding TM_real\n", m, lambda[m].real(), lambda[m].imag(), 1/LP.D, which);
    fprintf(logg, "\tlasing:%i (laser->real2lm[%i]=%i) real2lm[%i]=%i, lm2real[%i]=%i, Dth[%i]=%e, I[%i]=%e\n", lasing, m, laser->real2lm[m], m, real2lm[m], pntr-lasing, lm2real[pntr-lasing], m, Dth[m], m, I[m]);
#endif
#ifdef LASER_CHECK
    if (lasing) {
      double res=0;
      for (int j=0;j<dim;j++) res+=abs(omega*PsiaR[m*dim+j]-laser->lasing_mode[laser->real2lm[m]]->amplitude[j]);
      fprintf(logg, "\tTLM_approx lasing mode %i alternative amplitude residuum %e\n", m, res);
    }
#endif
  }
  if (pntr<Nl) {
    fprintf(logg, "less lasing modes (pntr=%i vs Nlm=%i) identified from real modes themself in TLM_approx initialize\n", pntr, Nl);
    Nl=pntr;
  }
    
  
  //claculate the maximal overlaps of two types. Consider only lasing modes as candidates for "closest" mode
  double *norm=new double[Nr*2];
  double *oM=new double[Nr*Nr*2];
  for (int type=0;type<2;type++) {
    //calculate norms
    for (int m=0;m<Nr;m++) norm[type*Nr+m]=sqrt(overlap_aux(dim,PsiaR+m*dim,PsiaR+m*dim, type));
    
    //calculate overlap matrix for <any | lasing>
    for (int m=0;m<Nr;m++) {
      for (int l=0;l<Nr;l++) {
	double aux=0;
	if (laser->real2lm[l]!=-1 && l!=m) aux=overlap_aux(dim,PsiaR+m*dim, PsiaR+l*dim, type)/(norm[type*Nr+m]*norm[type*Nr+l]);
        oM[type*Nr*Nr+m*Nr+l]=aux;
	//fprintf(logg,"calculating overlap (type %i): <real=%i|real=%i> (the latter lasing index %i) is %e with norms %e and %e\n", type, m, l, laser->real2lm[l], aux, norm[type*Nr+m], norm[type*Nr+l]);
      }
    }
    //find the maximum of each line
    for (int m=0;m<Nr;m++) {
      int lmax=m;
      for (int l=0;l<Nr;l++) if (oM[type*Nr*Nr+m*Nr+lmax]<oM[type*Nr*Nr+m*Nr+l]) lmax=l;
      oi[type][m]=lmax;
      o[type][m]=oM[type*Nr*Nr+m*Nr+lmax];
      	//fprintf(logg,"maximal overlap of type %i for real state %i is with real state %i lasing index %i being %e\n", type, m, lmax, laser->real2lm[lmax], o[type][m]);
    }
  }
  delete norm;
  delete oM;
  
  //calculate the xi matrix elements
  content* n4d=new content[dim];
  content* tp=new content[dim];
  for (int i=0;i<dim;i++) {
    n4d[i]=LP.n2*LP.n2*LP.n2x[i]*LP.n2x[i]/LP.d[i];  
    tp[i]=LP.total_profile[i];
  }

  for (int a=0;a<Nr;a++) {
    for (int b=0;b<Nr;b++) {
      LP.clear_braket_flags();
      LP.inputs[0]=PsiaL+a*dim;
      LP.inputs[1]=n4d;
      LP.inputs[2]=PsiaR+b*dim;
      LP.inputs[3]=PsiaR+b*dim; LP.conjugate[3]=true;
      LP.inputs[4]=PsiAaR+a*dim;
      xi[a*Nr+b]=LP.r->braketN(5, LP.inputs, LP.inverse, LP.conjugate)*LP.Gamma(k[b]);
#if LASER > 3
      fprintf(logg,"\t\t\t\txi matrix element: (a=%i, b=%i) is %.5e%+.5e\n",  a, b, xi[a*Nr+b].real(), xi[a*Nr+b].imag());
#endif
    }
//#ifdef LASER_CHECK
    //comparing/replacing the lambda by the TM of "the kernel inverse" "eigenvalue"
    LP.clear_braket_flags();
    LP.inputs[0]=PsiaL+a*dim;
    LP.inputs[1]=n4d;
    LP.inputs[2]=PsiAaR+a*dim;
    LP.inputs[3]=tp;
    content aux=LP.r->braketN(4,LP.inputs, LP.inverse, LP.conjugate);
#if LASER > 2
    fprintf(logg, "\t\t\t\tthe inverse of TM %.5e%+.5ei will replace the inverse of eigenvalue %.5e%+.5ei (residuum %.5e%+.5ei)\n", aux.real(), aux.imag(), (1.0/lambda[a]).real(), (1.0/lambda[a]).imag(), aux.real()-(1.0/lambda[a]).real(), aux.imag()-(1.0/lambda[a]).imag());
    lambda[a]=1.0/aux;
#endif
  }

  delete n4d;
  delete tp;
  delete PsiaL;
  delete PsiaR;
  delete PsiAaR;
  delete k;
  
  //check labels if some modes were not assigned, remove them from considerations
  for (int a=0;a<Nl;a++) if (lm2real[a]==-1) {
    for (int b=a;b<Nl-1;b++) lm2real[b]=lm2real[b+1];
    Nl--;
  }
  
}

void TLM_approx_prototype::inverse_xi(int i)
{
  if (Nl==0 && i==-1) return;
  int Nt=Nl;
  if (i!=-1) Nt++;
  for (int a=0;a<Nl;a++) if (lm2real[a]<0 || lm2real[a]>=Nr) {
    for (int b=0;b<Nl;b++) fprintf(logg,"ERROR: lm2real[%i]=%i\n",b,lm2real[b]);
    for (int b=0;b<Nr;b++) fprintf(logg,"ERROR: real2lm[%i]=%i\n",b,real2lm[b]);
    printf("corrupt labels in inverse_xi, exiting...\n");
    exit(1);
  }
  content* xi_aux=new content[Nt*Nt];
  for (int a=0;a<Nl;a++) for (int b=0;b<Nl;b++) xi_aux[a*Nt+b]=xi[lm2real[a]*Nr+lm2real[b]];  
  if (i!=-1) {
    for (int a=0;a<Nl;a++) {
      xi_aux[a*Nt+Nl]=xi[lm2real[a]*Nr+i]; 
      xi_aux[Nl*Nt+a]=xi[i*Nr+lm2real[a]];
    }
    xi_aux[Nl*Nt+Nl]=xi[i*Nr+i];
  }
  //fprintf(logg,"TLM will call inverse with Nt=%i\n",Nt);
  //for (int i=0;i<Nt;i++) fprintf(logg,"xi_aux[%i, %i]=%e
  //fflush(0);
  invert(xi_inv, xi_aux, Nt);
  #ifdef LASER_CHECK
  double res=0;
  for (int a=0;a<Nt;a++) {
    for (int b=0;b<Nt;b++) {
      content aux=0;
      if (a==b) aux=-1.0;
      for (int k=0;k<Nt;k++) aux+=xi_inv[a*Nt+k]*xi_aux[k*Nt+b];
      res+=abs(aux);
    }
  }
  if (res>1e-12*Nt*Nt) fprintf(logg,"WARNING: inverse check in inverse_xi gave residuum %e\n", res);
  #endif
  delete xi_aux;
}

void TLM_approx_prototype::intensities(double D)
{
  inverse_xi(-1);
  for (int a=0;a<Nl;a++) {
    dI[a]=0;
    for (int b=0;b<Nl;b++) {
      dI[a]+=xi_inv[a*Nl+b]*(D-1.0/lambda[lm2real[b]]);
#if LASER > 3
      fprintf(logg,"\t\t\t\txi_inv matrix element: (a=%i, b=%i) is %.5e%+.5e\n",  a, b, xi_inv[a*Nl+b].real(), xi_inv[a*Nl+b].imag());
#endif
    }
    #if LASER > 2
    fprintf(logg, "\t\t\tat pump %e, the intensity increment for lasing state a=%i (%i real) is %e (neglected imaginary part %e)\n", D, a, lm2real[a], dI[a].real(), dI[a].imag());
    for (int b=0;b<Nl;b++)  fprintf(logg, "\t\t\tthe contribution from lasing state %i (real %i) which is %e above the pump is through xi^-1 matrix element %e%+ei\n", b, lm2real[b], D-1.0/lambda[lm2real[b]].real(), xi_inv[a*Nl+b].real(), xi_inv[a*Nl+b].imag() );
    #endif
  }
#ifdef LASER_CHECK
  for (int a=0;a<Nl;a++) {
    content aux=D-1.0/lambda[lm2real[a]];
    for (int b=0;b<Nl;b++) aux-=xi[lm2real[a]*Nr+lm2real[b]]*dI[b];
    if (abs(aux)>1e-12) fprintf(logg,"WARNING: intensity residuum for lasing state %i (real %i) is %e\n",a,lm2real[a],abs(aux));
  }
#endif   
}

int TLM_approx_prototype::next_pump_threshold(double& D)
{
  if (Nl==Nr) {printf("no more non-lasing states in next_pump_threshold\n"); exit(1);}
  double ext=0;
  int which=-1;
  for (int i=0;i<Nr;i++) {
    if (real2lm[i]!=-1) continue;//this one is already lasing
    
    inverse_xi(-1);
    content lhs=1.0;
    for (int a=0;a<Nl;a++) for (int b=0;b<Nl;b++) lhs-=xi[i*Nr+lm2real[a]]*xi_inv[a*Nl+b];
    
    content rhs=1.0/lambda[i];
    for (int a=0;a<Nl;a++) for (int b=0;b<Nl;b++) rhs-=xi[i*Nr+lm2real[a]]*xi_inv[a*Nl+b]/lambda[lm2real[b]];
    
    double act=(rhs/lhs).real();
#if LASER > 2
    fprintf(logg, "\t\t\ttrying %i-th non-lasing mode (Dth=%e): lhs=%.5e%+.5ei and rhs=%.5e%+.5ei gives D=%e\n", i, Dth[i], lhs.real(), lhs.imag(), rhs.real(), rhs.imag(), act);
#endif
    if (which==-1 || (act<ext && act>0)) {which=i; ext=act;}
  }
#if LASER > 2
  fprintf(logg, "\t\tat pump %e the smallest positive pump threshold is %e for the next mode %i\n", D, ext, which);
#endif
  D=ext;
  return(which);
}

int TLM_approx_prototype::predict(lasing_state_char_prototype& lsc, FILE* file)
{
  double D=LP.D;//initial pump
  initialize();//fill in parameters
  int pntr=0;
  do {
    pntr++;
    intensities(D);
    double Dnext=D;
    int mode=-1;
    if (Nl<Nr) mode=next_pump_threshold(Dnext);
    if (Dnext==0 && D==0) {printf("something wrong, D=Dnext=0 in TLM_approx_predict\n"); exit(1);}
    
    if (file!=0) {
      fprintf(file, "%e %i %i %e ", D, Nl, mode, Dnext); 
      for (int i=0;i<Nl;i++) fprintf(file, "%e ", I[lm2real[i]]+dI[lm2real[i]].real());
      for (int i=0;i<Nr-Nl;i++) fprintf(file, "%e ", 0.0);
      fprintf(file, "\n");
    }
    
    if (lsc.Dth[0]==0) {
      if (mode==-1 || Dnext<D) break;
      D=Dnext;
    } else {
      if (D>lsc.Dth[0]) break;
      else D=min(D+lsc.Dth[1],Dnext);
    }
    if (D==Dnext) {
      Dth[mode]=Dnext;
      real2lm[mode]=Nl;
      lm2real[Nl]=mode;
      Nl++;
    }    
  } while (true);
  int* index=new int[Nl];
  for (int i=0;i<Nl;i++) {
    lsc.Dth[i]=Dth[lm2real[i]];
    lsc.k[i]=laser->TM_k[lm2real[i]];
    lsc.ipr[i]=ipr[lm2real[i]];
    lsc.x0[i]=x0[lm2real[i]];
    lsc.o[0][i]=o[0][lm2real[i]];
    lsc.o[1][i]=o[1][lm2real[i]];
    lsc.oi[0][i]=oi[0][lm2real[i]];
    lsc.oi[1][i]=oi[1][lm2real[i]];    
  }
  picsrt_index(Nl,lsc.Dth,index,0,'a');
  sort_from_index(Nl,lsc.Dth,index);
  sort_from_index(Nl,lsc.k,index);
  sort_from_index(Nl,lsc.ipr,index);
  sort_from_index(Nl,lsc.x0,index);
  sort_from_index(Nl,lsc.o[0],index);
  sort_from_index(Nl,lsc.o[1],index);
  sort_from_index(Nl,lsc.oi[0],index);
  sort_from_index(Nl,lsc.oi[1],index);
  //picsrt_inplace(Nl,resD,'a');
  delete index;
  return(Nl);
}




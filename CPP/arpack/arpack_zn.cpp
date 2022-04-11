// replace all "z"  by "c" ("z" for double) and save under new name!!!
#include "arpack_zn.h"
#define ARPACK 1
#define DIAG 2

#define EPS 1e-8

/*what is not variable for setting:
bmat - generalized or standard
iparam[7-1] - mode 1: A x = a x
HowMny - how many of eigenvec to return
*/
arpack_zn::arpack_zn(int dim, int howmany_eval,arpack_zn::precision tolerance,int max_iter_,bool eigen_vec, const char* which_eig,int explic_ncv)
{
#if ARPACK > 0
  fprintf(logg,"starting arpack with input: dim=%i, eig=%i, tol=%e, max_iter=%i, eigenvecs=%i, which=%s ncv=%i\n",dim, howmany_eval, tolerance, max_iter_, eigen_vec, which_eig, explic_ncv);
#endif

  nev=howmany_eval;
  ncv=max(nev*2+10,20);		//# of Lanczos vectors used
  ncv=min(ncv,dim);
  if (explic_ncv) ncv=explic_ncv;
  bmat='I';			//standard('I') or generalized ('G')
  which=new char[3];
  strcpy(which,which_eig);	//what kind of eigenvalues to arpack_zn::compute (f.e. "SM" = smallest)
  tol=tolerance;		//relative accuracy

//parameters for ..aupd

  max_iter=max_iter_;

  iparam=new arpack_zn::integer[11];

  n=dim;			//dimension of the problem = N
  resid=new arpack_zn::content[dim];		//vector with lenght N
  V=new arpack_zn::content[dim*ncv];		// array N times ncv
  ldv=dim;			//leading dimension of V
  ipntr=new arpack_zn::integer[14];	//lotsa OUTPUT

  workd=new arpack_zn::content[3*dim];	//array 3*N
  lworkl=3*ncv*ncv+5*ncv;	//=ncv*(ncv+8)
  workl=new arpack_zn::content[lworkl];	//array of dim lworkl
  rwork=new arpack_zn::precision[ncv];	//REAL array of dim ncv
  info=0;			//communication flag


//parameters for ..eupd

  rvec=eigen_vec;		//true = eigenvector also
  HowMny='A';			//'A' =all (=nev)     'S' -some, specified by select

  //!!!!!!!!!!!CHANGE: z and d will be now allocated outside!!!!!!!
  //d=new arpack_zn::content[nev+1];		//eigenvalues
  //z=new arpack_zn::content[dim*nev];		//array of dim N times nev (can be set to V)
  
  ldz=dim;			//The leading dimension of the array Z. In any case,  LDZ .ge. 1.
  sigma=0;			//shift (not significant for mode 1)

  selekt=new arpack_zn::logical[ncv];	// array of dim ncv : select[i] = true - eigenvector i wanted
  workev=new arpack_zn::content[2*ncv];	//workspace

  with_progress_bar=false;	//default - no progress bar
  
  #if ARPACK > 0
  fprintf(logg,"arpack arrays allocated with: nev=%i, ncv=%i, ldv=%i, lworkl=%i, ldz=%i, HowMny=%c\n",nev, ncv, ldv, lworkl, ldz, HowMny);
  #endif
  
};

arpack_zn::~arpack_zn()
{
  delete  resid;
  delete  V;
  
  //!!!CHANGE z and d will be now deleted outside
  //delete  d;
  //delete  z;
  
  delete  workd;
  delete  workl;
  delete  ipntr;
  delete  selekt;
  delete  iparam;
  delete  which;
  delete  workev;
  delete  rwork;
}

bool arpack_zn::make_sure(int i, int iter)
{
  if (i!=0) {
    char msg[100];
    if (iter<0) sprintf(msg,"zneupd"); else sprintf(msg,"znaupd after %i iterations",iter);
    if (i>0 || (i==-14 && iter<0)) {
      sprintf(buffer,"Failed to converge in arpack routine %s: info = %i \n\n returning...\n",msg,i);
      splachni(logg,buffer,1+4);
    }
    else {
      sprintf(buffer,"Input error occured in arpack routine %s: info = %i \n\n stopping...\n",msg,i);
      splachni(logg,buffer,1+4);
      exit(1);      
    }
  }
  return(i==0);
}


int arpack_zn::check_eigenset(int for_hmn,arpack_zn::pfAxy pf,int n, void *A, arpack_zn::content* d, arpack_zn::content *z, int where)
{
   arpack_zn::content *vec=new arpack_zn::content[n];

   //first check orthonormality
   int hmn_o=0;
   for (int i=0;i<for_hmn;i++) {
     for (int j=i;j<for_hmn;j++) {
       content braket=0;
       for (int k=0;k<n;k++) braket+=conj(z[i*n+k])*z[j*n+k];
       if (i==j) braket-=1.0;
       double res=abs(braket);
       if (res>EPS) {
         fprintf(logg,"vectors non-orthonormality above tolerance(%e):|<#%i|#%i>-%i|=%e\n",EPS,i,j,(i==j),res);
         hmn_o++;
       }
     }
   }
   //then eigen-ality
   int hmn_e=0;
   for (int i=0;i<for_hmn;i++) {
     (*pf)(n,A,&z[n*i],vec);       
     content braHmEket=0;
     for (int k=0;k<n;k++) braHmEket+=vec[k]-d[i]*z[i*n+k];
     double res=abs(braHmEket);
     if (res>EPS) {
       fprintf(logg,"vectors non-eigen-ality above tolerance: abs(H-E|%i>)=%e (corresponding eigenval: %e%+ei)\n",i,res, d[i].real(), d[i].imag());
       hmn_e++;
     }
   }
   delete[] vec;
   return(hmn_e+hmn_o);
 }

bool arpack_zn::go_for_it(int n, arpack_zn::pfAxy pf, void *A,ret_from_arpack_zn *output,bool check)
{
  bool failed=false;
  ido=info=0;
  
  iparam[1-1]=1;    //exact shifts
  iparam[3-1]=max_iter;   //max # of iterations
  iparam[4-1]=1;    //must be
  iparam[7-1]=1;    //mode 1 - A x = a x

  stopky s;

    /*info=1;
    srand(1);
    randomize();
    for (int i=0;i<n;i++) resid[i]=(rand()+0.0)/RAND_MAX;*/
  //sprintf(buffer, "calling ARPACK with parameters:\nn=%i, nev=%i, ncv=%i, ldv=%i\n",n,nev,ncv,ldv);
  //splachni(logg,buffer,1+4);

  output->times.main=output->times.after=0;
  output->times.user=output->times.checking=0;
  int step=0;
  do {
    step++;
    s.pusti();
    /*c  call ..aupd
     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
       IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO )*/
    F77NAME(znaupd)	(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, V, \
    			&ldv,iparam, ipntr, workd, workl, &lworkl,rwork, &info);
    output->times.main+=s.zastav();
    if (with_progress_bar) {
      progress_bar->add(0,1);
      progress_bar->print_message(true);
    }
    if (ido==1 || ido==-1) {
      s.pusti();
      (*pf)(n,A,&workd[ipntr[0]-1],&workd[ipntr[1]-1]);
      output->times.user+=s.zastav();
      // y= A x
      // x = workd[ipntr[1]], y = workd[ipntr[2]]
      // what means *pf is called with (&A,&x,&y)
    }
    else break;
  } while (true);

  if (!make_sure(info, step)) failed=true;

  z=output->eigenvecs;
  d=output->eigenvals;

  s.pusti();
  /*call cneupd
     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, WORKEV, BMAT,
       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD,
       WORKL, LWORKL, RWORK, INFO )*/
    if (!failed) {
        F77NAME(zneupd)	(&rvec, &HowMny, selekt, d, z, &ldz, &sigma, \
  		workev, &bmat, &n, which, &nev, &tol, resid, &ncv, V, &ldv, \
		iparam, ipntr, workd, workl, &lworkl,rwork, &info);
        //printf("ncv=%i, nev=%i\n",ncv, nev);
        if (!make_sure(info, -1)) failed=true;
    }
  output->times.after=s.zastav();

  //make sure, that all eigenvecs called for are correct
  s.pusti();
  int hmn_bad=0;
  if (rvec && check) hmn_bad=check_eigenset(nev,pf,n,A,d,z,0);
  if (hmn_bad>0) fprintf(logg,"number of failed checks on eigenvectors at EPS=%e: %i\n",EPS,hmn_bad);
  output->times.checking=s.zastav();

  output->perform.iter=iparam[3-1];
  output->perform.converged=nev-hmn_bad;
  output->perform.opx=iparam[9-1];
  output->perform.reort=iparam[11-1];

  return(!failed);
};

void arpack_zn::show_stats(ret_from_arpack_zn *vysl,FILE* file, int where_to, bool newline, bool header)
{
  if (header) {
    message(file,"timings & performance in arpack_zn routine:\n",where_to);
    message(file,"(1) dimension of the vector,\ntimes: (2)main routine ..aupd, (3) postprocessing routine ..eupd\n     (4) user provided matrix multiplication, (5) checking the validity\n",where_to);
    message(file,"perfomance: (6) # of Arnoldi iterations, (7) # of matrix OP*x multiplications, (8) # of reorthogonalizations (9) # of converged values\n",where_to);
    message(file,"    dim    main   after  user   check   iter   opx    reort converged\n",where_to);
  }
  #define mac(x) {\
  		if (x>9999999) sprintf(buffer,"%12i ",x);\
		else if (x>99999) sprintf(buffer,"%8i ",x);\
		else sprintf(buffer,"%6i ",x);\
		splachni(file,buffer,where_to);}
  mac(n)
  mac(vysl->times.main)
  mac(vysl->times.after)
  mac(vysl->times.user)
  mac(vysl->times.checking)
  mac(vysl->perform.iter)
  mac(vysl->perform.opx)
  mac(vysl->perform.reort)
  mac(vysl->perform.converged)

  #undef mac
  if (newline) message(file,"\n",where_to);
}

void arpack_zn::set_progress_bar(progress_bar_prototype* pb)
{
  progress_bar=pb;
  if (pb!=0) with_progress_bar=true; else with_progress_bar=false;
}


void arpack_zn::set_parameter(int ktory, void* vysledok)
{
  switch (ktory) {
    case 1: {//ido - arpack_zn::integer
      int *medzip=(int*) vysledok;
      ido=*medzip;
      break;
    }
    case 2: {//ncv - arpack_zn::integer
      int *medzip=(int*) vysledok;
      ncv=*medzip;
      break;
    }
    case 3: {//tolerance - arpack_zn::double
      double *medzip=(double*) vysledok;
      tol=*medzip;
      break;
    }
    default: {
      sprintf(buffer,"neplatny argument v arpack_zn.set_parameter:%i",ktory);
      splachni(logg,buffer,where_to_print_critical);
    }
  }
}
void arpack_zn::read_parameter(int ktory, void* vysledok)
{
  switch (ktory) {
    case 1: {//ido - arpack_zn::integer
      int *medzip=(int*) vysledok;
      *medzip=ido;
      break;
    }
    case 2: {//ncv - arpack_zn::integer
      int *medzip=(int*) vysledok;
      *medzip=ncv;
      break;
    }
    default: {
      sprintf(buffer,"neplatny argument v arpack_zn.read_parameter:%i",ktory);
      splachni(logg,buffer,where_to_print_critical);
    }
  }
}

int* (select_aux)(complex<double> *x) {int* i=new int; *i=true; printf("SELECT called???\n"); exit(0); return(i);}

lapack_zgee_prototype::lapack_zgee_prototype(int N)
{
  if (N<1 || N>10000) {printf("lapack_zgee_prototype input wrong N=%i\n",N); exit(1);}
  lapack_zgee_prototype::N=N;
  jobvl=jobvr=jobvs='V';
  sort='N';
  ldvr=ldvl=ldvs=lda=N;
  vl=vs=new complex<double>[ldvs*N];
  vr=new complex<double>[ldvr*N];
  rwork=new double[2*N];
  w=new complex<double>[N];
  A=new complex<double>[lda*N];
  work=new complex<double>[2*N];
  bwork=new int[N];
  
  //calculate the needed workspace
  lwork=-1;
  F77NAME(zgees)(& jobvs, &sort, select_aux, &N, A, &lda, &sdim, w, vs, &ldvs, work, &lwork, rwork, bwork, &info);
  int lwork1=(int) round(work[0].real());
  lwork=-1;
  F77NAME(zgeev)(& jobvl, &jobvr, &N, A, &lda, w, vl, &ldvl, vr, &ldvl, work, &lwork, rwork, &info);
  int lwork2=(int) round(work[0].real());
  lwork=max(lwork1,lwork2);
  #if DIAG > 0
  fprintf(logg,"lapack zgee class allocated with N=%i, lwork=%i (max of zgees=%i zgeev=%i)\n",N,lwork,lwork1,lwork2);
  #endif
  if (lwork<1 || info!=0) {printf("lapack_zgee_prototype error upon lwork calculation call: info=%i, lwork=%i\n",info, lwork); exit(1);}
  delete w;
  delete work;
  work=new complex<double>[lwork];
  jobvs='N';
}

lapack_zgee_prototype::~lapack_zgee_prototype()
{
  delete A;
  delete vs;
  delete vr;
  delete work;
  delete rwork;
  delete bwork;
}
int lapack_zgee_prototype::diagonalize(complex<double> *H, complex<double> *val, complex<double> *vecL, complex<double> *vecR, bool check)
{
  //FILE* file=fopen("zgeex_aux","w");
  //fprintf(file,"zgeex diag called: H=%p val=%p vecL=%p vecR=%p check=%i\nlda=%i N=%i\n", H, val, vecL, vecR, check, lda, N);
  //for (int i=0;i<lda;i++) {
    //for (int j=0;j<N;j++) fprintf(file,"%.14e%+.14ei ", H[i*N+j].real(), H[i*N+j].imag());
    //fprintf(file,"\n");
  //}
  //fclose(file);
  if (vecL==0 && vecR==0) check=false;
  for (int i=0;i<lda;i++) for (int j=0;j<N;j++) A[j*lda+i]=H[i*N+j];//!transpose for fortran
  
  if (vecL==0 && vecR==0) F77NAME(zgees)(& jobvs, &sort, select_aux, &N, A, &lda, &sdim, val, vs, &ldvs, work, &lwork, rwork, bwork, &info); 
  else F77NAME(zgeev)(&jobvl, &jobvr, &N, A, &lda, val, vecL, &ldvl, vecR, &ldvr, work, &lwork, rwork, &info);
  
  if (info<0) {printf("lapack_zgee_prototype error upon diagonalization call: info=%i\n",info); exit(1);}
  if (info!=0) {
    fprintf(logg,"lapack_zgee_prototype was unable to calculate all eigenvalues: info=%i\n",info);
    if (info>=N) info=0;
  }
  
  if (check) {
    double resL=0, resR=0;
    complex<double> auxL, auxR;
    for (int i=0;i<N;i++) {
      rwork[i]=rwork[i+N]=0;
      for (int j=0;j<N;j++) {
        if (vecL!=0) {
          auxL=conj(vecL[i*N+j])*val[i];
          for (int k=0;k<N;k++) auxL-=conj(vecL[i*N+k])*H[k*N+j];
          rwork[i]+=abs(auxL);
        }
        if (vecR!=0) {
          auxR=vecR[i*N+j]*val[i];
          for (int k=0;k<N;k++) auxR-=H[j*N+k]*vecR[i*N+k];
          rwork[i+N]+=abs(auxR);
        }
      }
      resL+=rwork[i];
      resR+=rwork[i+N];
    }
    if (max(resL,resR)>N*1e-10) fprintf(logg, "WARNING: lapack zgee l/r-eigensystems exact with precisions %e/%e\n",resL, resR);
    #if DIAG > 0
    fprintf(logg, "lapack zgee l/r-eigensystems exact with precisions %e/%e\n",resL, resR);
    if (resL>N*1e-10) for (int i=0;i<N;i++) fprintf(logg, "l-eigenvector #%i residuum %e\n", i, rwork[i]);
    if (resR>N*1e-10) for (int i=0;i<N;i++) fprintf(logg, "r-eigenvector #%i residuum %e\n", i, rwork[i+N]);
    #endif
  } 
  return(info);
}

lapack_zhee_prototype::lapack_zhee_prototype(int N)
{
  if (N<1 || N>10000) {printf("lapack_zhee_prototype input wrong N=%i\n",N); exit(1);}
  lapack_zhee_prototype::N=lda=N;
  jobz='V';
  uplo='U';
  rwork=new double[3*N-2];
  w=new double[N];
  A=new complex<double>[lda*N];
  work=new complex<double>[2*N-1];
  
  //calculate the needed workspace
  lwork=-1;
  F77NAME(zheev)(& jobz, &uplo, &N, A, &lda, w, work, &lwork, rwork, &info);
  int lwork=(int) round(work[0].real());
  #if DIAG > 0
  fprintf(logg,"lapack zhee class allocated with N=%i, lwork=%i\n",N,lwork);
  #endif
  if (lwork<1 || info!=0) {printf("lapack_zhee_prototype error upon lwork calculation call: info=%i, lwork=%i\n",info, lwork); exit(1);}
  delete w;
  delete work;
  work=new complex<double>[lwork];
}

lapack_zhee_prototype::~lapack_zhee_prototype()
{
  delete A;
  delete work;
  delete rwork;
}
int lapack_zhee_prototype::diagonalize(complex<double> *H, double *val, complex<double> *vec, bool check)
{
  if (vec==0) {jobz='N';check=false;} else jobz='V';
  for (int i=0;i<lda;i++) for (int j=0;j<N;j++) A[j*lda+i]=H[i*N+j];//!transpose for fortran
  F77NAME(zheev)(& jobz, &uplo, &N, A, &lda, val, work, &lwork, rwork, &info); 

  if (info<0) {printf("lapack_zhee_prototype error upon diagonalization call: info=%i\n",info); exit(1);}
  if (info!=0) {
    fprintf(logg,"lapack_zhee_prototype was unable to calculate all eigenvalues: info=%i\n",info);
    if (info>=N) info=0;
  }

  if (check) {
    double res=0;
    for (int i=0;i<N;i++) {
      rwork[i]=0;
      for (int j=0;j<N;j++) {
        complex<double> aux=conj(vec[i*N+j])*val[i];
        for (int k=0;k<N;k++) aux-=conj(vec[i*N+k])*H[k*N+j];
        rwork[i]+=abs(aux);
      }
      res+=rwork[i];
    }
    if (res>N*1e-14) fprintf(logg, "WARNING: lapack zhee eigensystem exact with precisions %e\n",res);
    #if DIAG > 0
    fprintf(logg, "lapack zhee eigensystem exact with precisions %e\n",res);
    if (res>N*1e-14) for (int i=0;i<N;i++) fprintf(logg, "eigenvector #%i residuum %e\n", i, rwork[i]);
    #endif
  } 
  return(info);
}


lapack_dstev_prototype::lapack_dstev_prototype(int N)
{
  if (N<1 || N>10000) {printf("lapack_dstev_prototype input wrong N=%i\n",N); exit(1);}
  lapack_dstev_prototype::N=ldz=N;
  work=new double[max(2*N-2,1)];
  E=new double[N];
  D=new double[N];
#if DIAG > 0
  fprintf(logg,"lapack dstev class allocated with N=%i\n",N);
#endif
}

lapack_dstev_prototype::~lapack_dstev_prototype()
{
  delete work;
  delete E;
  delete D;
}

int lapack_dstev_prototype::diagonalize(double *d, double *e, double *vals, double *vecs, bool check)
{
  if (vecs==0) {jobz='N';check=false;} else jobz='V';
  Z=vecs;
  for (int i=0;i<ldz;i++) {D[i]=d[i];E[i]=e[i];}//!copy the input matrix
#if DIAG > 1
  fprintf(logg,"dstev called with input:\n");
  for (int i=0;i<N;i++) fprintf(logg,"d[%i]=%e, e[%i]=%e\n",i,D[i],i,E[i]);
#endif
  F77NAME(dstev)(& jobz, &N, D, E, Z, &ldz, work, &info);
  for (int i=0;i<ldz;i++) vals[i]=D[i];
#if DIAG > 1
  fprintf(logg,"dstev returned output:\n");
  for (int i=0;i<N;i++) {
    if (jobz=='V') fprintf(logg,"vecs[%i,0]=%e ",i,vecs[i*N+0]);
    fprintf(logg,"vals[%i]=%e\n",i,vals[i]);
  }
#endif
  
  if (info<0) {printf("lapack_dstev_prototype error upon diagonalization call: info=%i\n",info); exit(1);}
  if (info!=0) {
    fprintf(logg,"lapack_dstev_prototype was unable to calculate all eigenvalues: info=%i\n",info);
  }
  
  if (check) {
    double res=0;
    for (int i=0;i<N;i++) {
      work[i]=0;
      for (int j=0;j<N;j++) {
        double aux=vecs[i*N+j]*(vals[i]-d[j]);
        if (j>0) aux-=vecs[i*N+j-1]*e[j-1];
        if (j<N-1) aux-=vecs[i*N+j+1]*e[j];
        work[i]+=abs(aux);
      }
      res+=work[i];
    }
    if (res>N*1e-14) fprintf(logg, "WARNING: lapack dstev eigensystem exact with precisions %e\n",res);
#if DIAG > 0
    fprintf(logg, "lapack dstev eigensystem exact with precisions %e\n",res);
    if (res>N*1e-14) for (int i=0;i<N;i++) fprintf(logg, "eigenvector #%i residuum %e\n", i, work[i]);
#endif
  }
  return(info);
}



lapack_zgesvd_prototype::lapack_zgesvd_prototype(int M, int N)
{
  if (N*M<1 || N>10000 || M>10000) {printf("lapack_zgesvd_prototype input wrong N=%i, M=%i\n",N,M); exit(1);}
  lapack_zgesvd_prototype::N=N;
  lapack_zgesvd_prototype::M=M;
  
  jobu=jobvt='A';
  lda=M;
  ldu=M;
  ldvt=N;
  rwork=new double[5*min(N,M)];
  S=new double[min(N,M)];
  A=new complex<double>[M*N];
  U=new complex<double>[M*M];
  VT=new complex<double>[N*N];
  work=new complex<double>[3*max(N,M)];
    
  //calculate the needed workspace
  lwork=-1;
  F77NAME(zgesvd)(&jobu, &jobvt, &M, &N, A, &lda, S, U, &ldu, VT, &ldvt, work, &lwork, rwork, &info);
  lwork=(int) round(work[0].real());
  #if DIAG > 0
  fprintf(logg,"lapack zgesvd class allocated with M=%i, N=%i, lwork=%i\n", M, N, lwork);
  #endif
  if (lwork<1 || info!=0) {printf("lapack_zgesvd_prototype error upon lwork calculation call: info=%i, lwork=%i\n",info, lwork); exit(1);}
  delete S;
  delete work;
  work=new complex<double>[lwork];
}

lapack_zgesvd_prototype::~lapack_zgesvd_prototype()
{
  delete A;
  delete U;
  delete VT;
  delete work;
  delete rwork;
}

int lapack_zgesvd_prototype::decompose(complex<double> *a, double *s, complex<double> *vt, complex<double>* u, bool check, int M, int N)
{
  if (N==-1) N=this->N;
  if (M==-1) M=this->M;
  if (N<1 || N>this->N || M<1 || M>this->M) {printf("lapack_zgesvd_prototype wrong input: N=%i (this->N=%i), M=%i (this->M=%i)\n",N,this->N,M,this->M); exit(1);}
    
  if (u==0) {jobu='N';check=false;} else jobu='A';  
  if (vt==0) {jobvt='N';check=false;} else jobvt='A';
  
  for (int i=0;i<M;i++) for (int j=0;j<N;j++) A[i+j*M]=a[i*N+j];//!transpose for fortran

  F77NAME(zgesvd)(&jobu, &jobvt, &M, &N, A, &M, s, U, &M, VT, &N, work, &lwork, rwork, &info);

  if (info<0) {printf("lapack_zgesvd_prototype error upon decomposition call: info=%i\n",info); exit(1);}
  if (info!=0) fprintf(logg,"lapack_zgesvd_prototype was unable to calculate all singular values: info=%i\n",info);
  
#if DIAG > 1
  fprintf(logg, "lapack zgesvd returned the following:\n");
  fprintf(logg, "U matrix:\n");
  for (int i=0;i<M;i++) {
    for (int j=0;j<M;j++) fprintf(logg,"%.4e%+.4e ",U[i*M+j].real(),U[i*M+j].imag());
    fprintf(logg,"\n");
  }
  fprintf(logg, "s values:\n");
  for (int j=0;j<min(M,N);j++) fprintf(logg,"%.4e ",s[j]);
  fprintf(logg,"\n");

  fprintf(logg, "VT matrix:\n");
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) fprintf(logg,"%.4e%+.4e ",VT[i*N+j].real(),VT[i*N+j].imag());
    fprintf(logg,"\n");
  }
  fprintf(logg, "for the original a matrix:\n");
  for (int i=0;i<M;i++) {
    for (int j=0;j<N;j++) fprintf(logg,"%.4e%+.4e ",a[i*N+j].real(),a[i*N+j].imag());
    fprintf(logg,"\n");
  }
#endif
  
  if (jobu!='N') for (int i=0;i<M;i++) for (int j=0;j<M;j++) u[i*M+j]=U[j*M+i];   //transpose back
  if (jobvt!='N') for (int i=0;i<N;i++) for (int j=0;j<N;j++) vt[i*N+j]=VT[j*N+i];//transpose back
  
  if (check) {
    double res=0;
    for (int i=0;i<M;i++) {
      for (int j=0;j<N;j++) {
	complex<double> aux=a[i*N+j];
        for (int k=0;k<min(N,M);k++) aux-=u[i*M+k]*s[k]*vt[k*N+j];
        res+=abs(aux);
#if DIAG > 1
    fprintf(logg, "lapack zgesvd eigensystem check element %i-%i: residuum=%e\n",i,j,abs(aux));
#endif

      }
    }
    if (res>max(N,M)*1e-12) fprintf(logg, "WARNING: lapack zgesvd eigensystem exact with precision %e\n",res);
    #if DIAG > 0
    fprintf(logg, "lapack zgesvd eigensystem exact with precision %e\n",res);
    #endif
  } 
  return(info);
}

/*int lapack_zgesvd_prototype::decompose2(complex<double> *a, double *s, complex<double> *vt, complex<double>* u, bool check)
{
  if (u==0) {jobu='N';u=U;} else jobu='A';  
  if (vt==0) {jobvt='N';vt=VT;} else jobvt='A';

  if (check) for (int i=0;i<M;i++) for (int j=0;j<N;j++) A[i*N+j]=a[i*N+j];//input copy

#if DIAG < 2
  if (!check) for (int i=0;i<M;i++) for (int j=0;j<N;j++) A[i*N+j]=a[i*N+j];//input copy
#endif

  F77NAME(zgesvd)(&jobu, &jobvt, &M, &N, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, &info);

  if (info<0) {printf("lapack_zgesvd_prototype error upon decomposition call: info=%i\n",info); exit(1);}
  if (info!=0) {
    fprintf(logg,"lapack_zgesvd_prototype was unable to calculate all singular values: info=%i\n",info);
    if (info>=N) info=0;
  }
  
#if DIAG > 1
  fprintf(logg, "lapack zgesvd returned the following:\n");
  fprintf(logg, "u matrix:\n");
  for (int i=0;i<M;i++) {
    for (int j=0;j<M;j++) fprintf(logg,"%.4e%+.4e ",u[i*M+j].real(),u[i*M+j].imag());
    fprintf(logg,"\n");
  }
  fprintf(logg, "s values:\n");
  for (int j=0;j<max(M,N);j++) fprintf(logg,"%.4e ",s[j]);
  fprintf(logg,"\n");

  fprintf(logg, "vt matrix:\n");
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) fprintf(logg,"%.4e%+.4e ",vt[i*N+j].real(),vt[i*N+j].imag());
    fprintf(logg,"\n");
  }
  fprintf(logg, "for the original a matrix:\n");
  for (int i=0;i<M;i++) {
    for (int j=0;j<N;j++) fprintf(logg,"%.4e%+.4e ",A[i*N+j].real(),A[i*N+j].imag());
    fprintf(logg,"\n");
  }
#endif
  
  if (check) {
    double res=0;
    for (int i=0;i<M;i++) {
      for (int j=0;j<N;j++) {
	complex<double> aux=A[i*N+j];
        for (int k=0;k<min(N,M);k++) aux-=u[i*M+k]*s[k]*vt[k*N+j];
        res+=abs(aux);
#if DIAG > 1
    fprintf(logg, "lapack zgesvd eigensystem check element %i-%i: residuum=%e\n",i,j,abs(aux));
#endif

      }
    }
    if (res>max(N,M)*1e-12) fprintf(logg, "WARNING: lapack zgesvd eigensystem exact with precision %e\n",res);
    #if DIAG > 0
    fprintf(logg, "lapack zgesvd eigensystem exact with precision %e\n",res);
    #endif
  } 
  return(info);
}*/


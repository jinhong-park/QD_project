#include "main.h"

//extern void hamiltonian1(int, void *, content* ,content *);

//richardson interpolation
//N - number of measured points, E - measured values, h - values of the parameter
//outcome[4] - interp from 3 point, from 2 point (best), from 2 point (mean), single value
//!!! allocated outside!!!
//return value - the lowest existing interpol in the previuos row
//at lowest parameter, if 2 point interp exists otherwise 0
//returning value the best interp. extracted
//tol - the intrap. value will be exepted, if its distance from the best single value is less than tol*(max single valu- min)
//e - exponent of the first error term, s - step in exponents in error terms
prec richardson(int N,int e, int s, prec *E, prec *h, prec tol, prec* outcome)
{
//#define FIT1
  prec *r2=0;if (N>1) r2=new prec[N-1];	//fit z dvojic
  prec *r3=0;if (N>2) r3=new prec[N-2];	//fit z trojic
  bool *r2acc=0; if (N>1) r2acc=new bool[N-1];	//accepted for final outcome?
  bool *r3acc=0; if (N>2) r3acc=new bool[N-2];	//accepted for final outcome?

  bool este=true, first=true;
  int vyh=0;
  prec val,best=0;

  //zisti rozptyl
  prec min=E[0],max=E[0];
  for (int i=1;i<N;i++) {
    if (E[i]<min) min=E[i];
    if (E[i]>max) max=E[i];
  }
  prec diff=max-min;

  while (este) {
  //usporiadaj ich vzostupne podla parametra
  for (int i=0;i<N-1;i++) {
    for (int j=0;j<N-1-i;j++) {
      if (h[j]>h[j+1]) {	//flip it
        val=h[j]; h[j]=h[j+1]; h[j+1]=val;
	val=E[j]; E[j]=E[j+1]; E[j+1]=val;
      }
    }
  }

  if (first) best=E[0];		//best guess - least parameter value
  first=false;

  //kontrola - vyhod vsetky nemonotonne
  prec znam=E[0]-E[1];
  este=false;
  for (int i=1;i<N-1;i++) {
    if ((E[i]-E[i+1])*znam<0) {
      E[i+1]=E[N-1];
      h[i+1]=h[N-1];
      N--;
      vyh++;
      este=true;
      break;
    }
  }
  }

  //mozes hadat - najprv po dvojiciach, podla h^4
  prec q,qe,r2m=0,r2m2=0,r2mc=0,r2m2c=0,weight2=0,w,we;
  int hmn2=0;
  for (int i=0;i<N-1;i++) {
    q=h[i]/h[i+1];
    qe=1; for (int j=0;j<e;j++) qe*=q;
    r2[i]=(E[i]-qe*E[i+1])/(1-qe);
    if (fabs((r2[i]-best))<tol*diff) {
      hmn2++;
      w=2/(h[i]+h[i+1]);
      we=1; for (int j=0;j<e;j++) we*=w;
      r2mc+=r2[i]*we;
      r2m2c+=r2[i]*r2[i]*we;
      r2acc[i]=true;
      weight2+=we;
    }
    else r2acc[i]=false;
    r2m+=r2[i];
    r2m2+=r2[i]*r2[i];
  }
  if (N>1) {r2m/=N-1;r2m2/=N-1;}
  prec err2=0;
  if (hmn2>0) {
    r2mc/=weight2;r2m2c/=weight2;
    if (hmn2>1) err2=sqrt((r2m2c-r2mc*r2mc)/weight2);
  }


  //teraz po trojiciach
  prec p,qs,ps,pe,r3m=0,r3m2=0,r3mc=0,r3m2c=0,weight3=0,we_s;
  int hmn3=0;
  for (int i=0;i<N-2;i++) {
    q=h[i]/h[i+1];
    qs=1;for (int j=0;j<s;j++,qs*=q);
    qe=1;for (int j=0;j<e;j++,qe*=q);
    p=h[i+1]/h[i+2];
    ps=1;for (int j=0;j<s;j++,ps*=p);
    pe=1;for (int j=0;j<e;j++,pe*=p);
    r3[i]=(E[i]-qe*E[i+1])*(1-1/ps)-qe*(E[i+1]-pe*E[i+2])*(qs-1);
    r3[i]/=(1-qe)*(1-1/ps)+(1-pe)*(1-qs)*qe;
    if (fabs((r3[i]-best))<tol*diff) {
      w=3/(h[i]+h[i+1]+h[i+2]);
      we_s=1; for(int j=0;j<e+s;j++,we_s*=w);
      weight3+=we_s;
      hmn3++;
      r3mc+=r3[i]*we_s;
      r3m2c+=r3[i]*r3[i]*we_s;
      r3acc[i]=true;
    }
    else r3acc[i]=false;
    r3m+=r3[i];
    r3m2+=r3[i]*r3[i];
  }
  if (N>2) { r3m/=N-2;r3m2/=N-2;}
  prec err3=0;
  if (hmn3>0) {
    r3mc/=weight3;r3m2c/=weight3;
    if (hmn3>1) err3=sqrt((r3m2c-r3mc*r3mc)/weight3);
  }

#ifdef FIT1
  //tlac vysledkov
  sprintf(buffer,"Given %i points, used %i points, with results\n",N+vyh, N);
  splachni(logg,buffer,4);

  sprintf(buffer,"With tolerance %f, best single value %e and variance %e\naccepted %i couple guesses and %i tripple guesses\n",tol,best,diff,hmn2,hmn3);
  splachni(logg,buffer,4);

  if (hmn3>0) {
    sprintf(buffer,"thus getting mean accepted tripple guess at E=%.8e +- %e\n",r3mc,err3);
    splachni(logg,buffer,4);
  }
  else message(file,"no accepted tripple guesses :(\n",4);

  if (hmn2>0) {
    sprintf(buffer,"and getting mean accepted couple guess at E=%.8e +- %e\n",r2mc,err2);
    splachni(logg,buffer,4);
  }
  else message(file,"no accepted couple guesses :(\n",4);


  sprintf(buffer,"without worrying about tolerance:\n tripple:%.8e +- %e\n couple:%.8e +- %e\n",r3m,sqrt((r3m2-r3m*r3m)/(N-2)) ,r2m,sqrt((r2m2-r2m*r2m)/(N-1)));
  splachni(logg,buffer,4);

  message(file,"parameter, value, guess from couple, accepted?, from tripple, accepted?\n",4);
  for (int i=0;i<N+vyh;i++) {
    sprintf(buffer,"%e\t%e\t",h[i],E[i]);
    splachni(logg, buffer,4);
    if (i<N-1) {
      sprintf(buffer,"%e ",r2[i]);
      splachni(logg, buffer,4);
      if (r2acc[i]) message(file," accepted ",4);
      else message(file," not acc. ",4);
    }
   if (i<N-2) {
      sprintf(buffer,"%e",r3[i]);
      splachni(logg, buffer,4);
      if (r3acc[i]) message(file," accepted\n",4);
      else message(file," not acc.\n",4);
    }
    else message(file,"\n",4);
  }
#endif

  //koniec - give the best guess

  for (int i=0;i<4;i++) outcome[i]=0;
  if (hmn3>0) {
    val=r3[0];
    outcome[0]=r3[0];
  }
  else if (hmn2>0) val=r2[0];
      else val=best;
  if (hmn2>0) {
    outcome[1]=r2[0];
    outcome[2]=r2mc;
  }
  outcome[3]=best;

  if (N>1) delete r2;
  if (N>2) delete r3;
  if (N>1) delete r2acc;
  if (N>2) delete r3acc;
#ifdef FIT1
  sprintf(buffer,"Finally returning: highes achieved interp., 3 point, 2 point, 2 point mean, single\n%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",val,outcome[0],outcome[1],outcome[2],outcome[3]);
  splachni(logg,buffer,4);
#endif
  return(val);
}

//N number of measurements, eig # of eigenvecs, e,s exponents of error for richardson
//dim, stepdim- dim for measurements, vysl[5] - from outside with some space allocated
//here (0,..,4) are for the: best achieved interp, 3 point int, 2 point, 2 point mean, single
//returns total time of computation [ms]
int fitvalues(int eig, int e, int s, ret_from_arpack_zn* vysl, int N, int stepdim)
{

  if (N<0) {message("wrong N in fitvalues\n");exit(1);}
  
#define FIT0
  prec dim_x=parameters.values[parameters_class::dim_x];
  prec dim_y=parameters.values[parameters_class::dim_y];
  prec length_x=parameters.values[parameters_class::length_x];
  prec length_y=parameters.values[parameters_class::length_y];
  prec dLx=parameters.values[parameters_class::dLx];
  prec dLy=parameters.values[parameters_class::dLy];
  parameters.list(logg,4); 

  int dim_x0=(int) dim_x, dim_y0=(int) dim_y;

  prec *E=new prec[N*eig];
  prec *h=new prec[N*eig];
  int size,time=0;

  for (int i=0;i<N;i++) {
#ifdef FIT0
    sprintf(buffer,"(%i out of %i):extracting at dim=%i,%i...",i+1,N,dim_x0,dim_y0);
    splachni(logg,buffer,4);
#endif
    reg1 r1(dim_x0,dLx-length_x,2*length_x,'d',dim_y0,dLy-length_y,2*length_y,'d',inside,2);

    size=r1.sets*r1.Nw;
    arpack_zn arp1(size,eig,-1e-12,100000,true,(char*) "SM");
    hamiltonian_init(r1);
    arp1.go_for_it(size,hamiltonian,(void*) &r1,vysl,true);
    time+=vysl->times.main+vysl->times.user;
    arp1.show_stats(vysl,logg,4,true,true);
    sort_0(size,eig,vysl,'a','a');
#ifdef FIT0
    message(logg,"done\n",4);
#endif

    for (int j=0;j<eig;j++) {
      E[j*N+i]=vysl->eigenvals[j].real();
      h[j*N+i]=1.0/dim_x0;
    }
    dim_x0+=stepdim;
    dim_y0+=stepdim;
  }
  prec vals[4];
  for (int i=0;i<eig;i++) {
    vysl->eigenvals[i]=richardson(N,e,s,E+i*N,h+i*N,100.0,vals);
  }
  delete E;
  delete h;

  return(time);
}
void clean_up(state_info* states) 
{
  CoulombElement(*states,0,0,0,0,0,-1);
  delete states->r1;
  delete states->vysl->eigenvals;
  delete states->vysl->eigenvecs;
  if (states->vysl->eigenvecs_plaq) delete states->vysl->eigenvecs_plaq;
  if (states->vysl->eigenvecs_orig) delete states->vysl->eigenvecs_orig;
  delete states->vysl;
  delete states;
}

extern "C"
{
  void F77NAME(zgetrf)( int* , int* , content* , int* , int* , int*);
  void F77NAME(zgbtrf)( int* , int* , int* , int* , content* , int* , int* , int*);
}

state_info* diagonalize()
{
  message(logg,"reading parameters...",1+4);
  parameters.recompute();
  parameters.list(logg,1+4);
  int dim_x=(int) parameters.values[parameters_class::dim_x];
  int dim_y=(int) parameters.values[parameters_class::dim_y];
  prec length_x=parameters.values[parameters_class::length_x];
  prec length_y=parameters.values[parameters_class::length_y];
  prec dLx=parameters.values[parameters_class::dLx];
  prec dLy=parameters.values[parameters_class::dLy];
  int eig=(int) parameters.values[parameters_class::eig];
  int eigmax=(int) parameters.values[parameters_class::eigmax];
  int geometry=(int) parameters.values[parameters_class::geometry];
  int sets=(int) parameters.values[parameters_class::sets];
  
  if (sets!=dof.length) {printf("inconsistency in parameters: sets=%i, while dof.length=%i; check dof class in classes.h\n",sets, dof.length);exit(1);}
  
  message(logg,"creating region...",1+4);
  reg1* r1;
  
  if (geometry!=1) {
    //2D potential
   r1=new reg1(dim_x,dLx-length_x,2*length_x,'d',dim_y,dLy-length_y,2*length_y,'d',inside,sets);
  }
  else {    //1D ring - periodic boundary conditions in x
    r1=new reg1(dim_x,dLx-length_x,2*length_x,'p',dim_y,dLy-length_y,2*length_y,'d',inside,sets);
  }

  r1->act_prec_set((region::operators_precision) (int) parameters.read("precision","-"));
  r1->symmetrize_ops((bool) parameters.read("peierls","-"), true); 
  r1->set_inv_angle(parameters.values[parameters_class::phi_d]);
  //test1(16,16,logg,1+4);
  
  message(logg,"allocating space...",1+4);
  ret_from_arpack_zn* vysl=new ret_from_arpack_zn;
  vysl->eigenvecs=new content[eig*r1->sets*r1->Nw];
  vysl->eigenvals=new content[eig+1];
  vysl->eigenvecs_plaq=vysl->eigenvecs_orig=0;
  
  //!!!nuclear spins
  if (parameters.read("exch","meV")*parameters.read("xMn","-")*parameters.read("nuclear_state","-")!=0) spin_impurities.reinitialize(r1,(int) parameters.read("nuclear_state","-") );

  if ((int) parameters.read("system","-") == 'g') {
    #define BRUTE_FORCExx
    //fprintf(logg,"arpack_dim=%i from sets=%i and Nw_g=%i\n",r1->g_length,r1->sets,r1->Nw_g);
    message(logg,"creating arpack...",1+4);
    arpack_zn arp(r1->g_length,eig,-1e-12,100000,true,(char*) "LM");
    message(logg,"building hamiltonian...\n",1+4);
    progress_bar_prototype progress_bar("");
    sprintf(buffer,"linear dimension:%i, ",r1->g_length);
    splachni(logg,buffer,1+4);

    #ifdef BRUTE_FORCE
    if (r1->HG!=0) delete r1->HG; 
    r1->HG=new content[r1->g_length*r1->g_length];
    if (r1->M!=0) delete r1->M; 
    r1->M=new content[r1->g_length*r1->g_length];
    for (int i=0;i<r1->g_length*r1->g_length;i++) r1->HG[i]=r1->M[i]=0;
    sprintf(buffer,"hamiltonian:%.0f MB, eigenvectors:%.0f MB\n",r1->g_length*r1->g_length*3*sizeof(content)/1e6, eig*r1->g_length*sizeof(content)/1e6);
    splachni(logg,buffer,1+4);
    #endif
    delete r1->HG_sparse; r1->HG_sparse=new sparse_matrix_prototype(r1->g_length,20);
    delete r1->M_sparse; r1->M_sparse=new sparse_matrix_prototype(r1->g_length,4);
    hamiltonian_init_graphene(*r1);
    //check_for_hermitivity(*r1,logg,1+4);

    //r1->ant(mass_gr, units.length/units.nm, units.energy/units.meV); 

    #ifdef BRUTE_FORCE
    r1->convertM2locvec();
    if (r1->LU!=0) delete r1->LU; 
    r1->LU=new content[r1->g_length*r1->g_length];
    if (r1->IPIV!=0) delete r1->IPIV; 
    r1->IPIV=new int[r1->g_length];
    for (int n=0;n<r1->g_length*r1->g_length;n++) r1->LU[n]=r1->HG[n];
    int INFO;
    progress_bar.reset("LU-factorization");
    progress_bar.start();
    F77NAME(zgetrf)(&r1->g_length,&r1->g_length,r1->LU,&r1->g_length,r1->IPIV,&INFO);
    progress_bar.finished(true);
    if (INFO!=0) {printf("graphene hamiltonian could not be LU-factorized: zgetrf INFO=%i\n",INFO);exit(1);}
    #endif
    
    r1->HG_sparse->matrix_structure(false);
    r1->HG_sparse->matrix2banded();
    sprintf(buffer,"hamiltonian-b:%.0f MB, eigenvectors:%.0f MB\n",r1->g_length*r1->HG_sparse->LDAB*sizeof(content)/1e6, 3*eig*r1->g_length*sizeof(content)/1e6);
    splachni(logg,buffer,1+4);
    //r1->ant(mass_gr, units.length/units.nm, units.energy/units.meV); 
    
    progress_bar.reset("LU-b/factorization");
    progress_bar.start();
    r1->HG_sparse->LU_banded_factorize();
    progress_bar.finished(true);

    progress_bar.reset("arpack iterations");
    progress_bar.add(0,0);
    progress_bar.start();
    arp.set_progress_bar(&progress_bar);
    #ifdef BRUTE_FORCE
    arp.go_for_it(r1->g_length, hamiltonian_graphene_with_checking,(void*) r1,vysl,true);
    #else
    arp.go_for_it(r1->g_length, hamiltonian_graphene,(void*) r1,vysl,true);
    #endif
    progress_bar.finished(true);
    arp.show_stats(vysl,logg,4,true,true);
    #ifdef BRUTE_FORCE
    r1->graphene_direct_check(10, vysl->eigenvecs, vysl->eigenvals);
    #endif
    
    //the eigenvalues
    for (int i=0;i<eig;i++) {
      vysl->eigenvals[i]=1.0/vysl->eigenvals[i];
      if (vysl->eigenvals[i].imag()>1e-8)
      fprintf(logg,"warning: non-zero imaginary part of energy: state %i E=%e%+ei\n",i,vysl->eigenvals[i].real(), vysl->eigenvals[i].imag());
    }
    
    message(logg,"presorting...",1+4);
    sort_0(r1->g_length,eig,vysl,'r','a');
    
    message(logg,"converting eigensolutions...",1+4);
    //convert and normalize
    //original to plaquette
    vysl->eigenvecs_plaq=new content[eig*r1->g_length];
    for (int i=0;i<eig;i++) r1->M_sparse->operate(vysl->eigenvecs+i*r1->g_length, vysl->eigenvecs_plaq+i*r1->g_length);
    normalize(r1->g_length, eig, vysl->eigenvecs_plaq);
    
    sprintf(buffer,"plaq-eigenvector basis ONO with residuum %e\n",is_orthonormal(r1->g_length,eig,vysl->eigenvecs_plaq,0,false));
    splachni(logg,buffer,1+4);

    //original to standard
    vysl->eigenvecs_orig=vysl->eigenvecs;
    vysl->eigenvecs=new content[eig*r1->sets*r1->Nw];
    for (int i=0;i<eig;i++) {//eigenvector
      for (int j=0;j<r1->sets/2;j++) {// (grid+sigma) components
        r1->psip2psig(vysl->eigenvecs_orig+i*r1->g_length, j, vysl->eigenvecs+(2*i*r1->sets/2+2*j)*r1->Nw, vysl->eigenvecs+(2*i*r1->sets/2+2*j+1)*r1->Nw);
      }
    }
    normalize(r1->sets*r1->Nw, eig, vysl->eigenvecs);
  }
  
  else  {
    message(logg,"creating arpack...",1+4);
    arpack_zn arp(r1->sets*r1->Nw,eig,-1e-12,100000,true,(char*)"SM");
    //r1->ant(potential, units.length/units.nm, units.energy/units.meV);
    message(logg,"building hamiltonian...\n",1+4);
    hamiltonian_init(*r1);
    //check_for_hermitivity(*r1,logg,1+4);
    
    progress_bar_prototype progress_bar("diagonalizing");
    progress_bar.add((int) round(0.4*pow(r1->sets*r1->Nw,1.8/2)),0);
    progress_bar.start();
    arp.set_progress_bar(&progress_bar);
    arp.go_for_it(r1->sets*r1->Nw,hamiltonian,(void*) r1,vysl,true);
    progress_bar.finished(true);
    arp.show_stats(vysl,logg,4,true,true);
    
    message(logg,"presorting...",1+4);
    sort_0(r1->sets*r1->Nw,eig,vysl,'a','a');
  }

  //message(logg,"fitting phase...",1+4);
  //fit_the_phase(vysl_aux,*r1,eig,3,1e-3);

  message(logg,"creating state_info...done\n",1+4);
  state_info* states=new state_info(r1,vysl,eig,eigmax);
  
  return(states);
}

void spot_crossing(state_info& I, prec *param_to_change, FILE* file, int where_to)
{
  static prec last_energies[2][1000];
  static prec last_param=0;
  static int filled[1000]={0};
  
  prec& param=*param_to_change;
  prec param_old=param;
  
  for (int u1=0;u1<I.unsorted;u1++) {
    int s1=I.unsorted2s[u1];
    for (int u2=u1+1;u2<min(I.unsorted,u1+4);u2++) {
      int s2=I.unsorted2s[u2];
      if (filled[s1]<2 || filled[s2]<2) continue;
      prec val1=last_energies[1][s1]-last_energies[1][s2];
      prec val2=last_energies[0][s1]-last_energies[0][s2];
      prec val3=(I.vysl->eigenvals[u1]-I.vysl->eigenvals[u2]).real();
      if (val2*val3<0 || (fabs(val1)>fabs(val2) && fabs(val3)>fabs(val2))) {
        //crossing found
#ifdef MINI0        
        sprintf(buffer,"(anti)crossing of states %i-%i of width %e probably occurs at param=%e to %e",u1,u2,last_param,param);
        splachni(file[1],buffer,4);
#endif

        prec step=param-last_param;
        param=last_param;
        prec min=minimize(u1,u2,I.unsorted,&param,val2,step);
        param=param_old;
      }
    }
  }
  for (int s=0;s<I.sorted;s++) {
    int u=I.sorted2un[s];
    last_energies[1][s]=last_energies[0][s];
    if (u!=-1) {
      last_energies[0][s]=I.vysl->eigenvals[u].real();
      filled[s]+=1;
    }
    else filled[s]=0;
  }
}

prec minimize(int a, int b, int eig, prec* x, prec start_val, prec step, int dimstep, int N) 
{
  prec x_val=*x,val;
  ret_from_arpack_zn vysl;
  step/=2;
  int i;
  
  for (i=1;i>0;i++) {
    
    *x+=step;
    state_info* states=diagonalize();
    val=(states->vysl->eigenvals[a]-states->vysl->eigenvals[b]).real();
    clean_up(states);

#ifdef MINI1    
    sprintf(buffer,"at param %e, step %e value %e found\n(states: %i, %i)\n",*x,step,val,a,b);
    splachni(logg,buffer,4);
#endif

    if (fabs(val)>start_val) {
      *x-=step;
      step/=2;
#ifdef MINI1    
    message(logg,"stepping back\n",4);
#endif

    }
    else {
      start_val=val;
#ifdef MINI1    
      message(logg,"stepping ahead\n",4);
#endif
    }
    if (step<1e-10) break;
  }
#ifdef MINI1    
    sprintf(buffer,"after %i steps, restoring parameter at %e\nreturning value of %e found at param %e\n",i,x_val,val,*x);
    splachni(logg,buffer,4);
#endif

  *x=x_val;
  return(val);
}


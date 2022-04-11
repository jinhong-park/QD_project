#include "laser.h"

//D(x) / D_0
content pump_envelope(prec x, prec y) 
{
#ifdef HALF_PUMPED
  double d=0.01;
  if (x>0) d=1;
#else
  double l0=100*units.nm/units.length;
  double x0=-500*units.nm/units.length;
  double d=1/(1+((x-x0)/l0)*((x-x0)/l0));
  if (x<-800*units.nm/units.length || x>-200*units.nm/units.length) d=0.01; else d=1;
#endif
  return(content(d,0));
}


 
//n(x) (only the envelope, the overall scale squared is in LP.n2)
content refractive_index_envelope(prec x, prec y)
{
  //uniform
  //return(content(1,0));
  
  //smooth
  //double l0=-500*units.nm/units.length;
  //n=1.0/(1+(x/l0)*(x/l0));
  
  //two regions
  //content n=1.0;
  //if (x<0) n=1.0;
  
  static double xi=8*units.nm/units.length;//correlation length of the disorder
  static double nMin=2, nMax=10;//the refractive index is randomly distributed from nMin to nMax
  static bool init=false;//is it already initialized
  static double* values=0;//randomly generated values for refractive index
  static double x0=0;  //the coordinate of the left end of the cavity
  static int length=0;//number of entries in the array
  if (!init) {
    if (values!=0) delete values;
    length=(int) (parameters.read("length_x","nm")*units.nm/units.length*2/xi)+1;
    values=new double[length];
    for (int i=0;i<length;i++) values[i]=((nMax-nMin)*generuj()+nMin)/LP.n2;
#ifdef SYMMETRIC
    int half=(length-1)/2;
    for (int i=0;i<=half;i++) values[length-1-i]=values[i];
#endif
    init=true;
    x0=x;
  }
#ifdef SYMMETRIC
  if (x<0) x=-x;
#endif
  int index=(int) floor((x-x0)/xi);
  double res=(x-x0)/xi-index;
  //res=0;
  if (index>=length-1) index=length-2;
  if (index<0) index=0;
  return(content(values[index]*(1-res)+values[index+1]*res) + 0.001*content(0,1));
  //interacting: xi=8, Min=2, Max=10, im=0.4*0.25*i=0.1i, eig=80, lv=2, flat_x=1000, LP.n2=1.5*1.5, omegaa=60
}

//n(x) (only the envelope, the overall scale squared is in LP.n2)
content refractive_index_envelope2(prec x, prec y)
{
  static double xi=20*units.nm/units.length;//correlation length of the disorder
  static double n2Min=2, n2Max=10;//the refractive index is randomly distributed from nMin to nMax
  static bool init=false;//is it already initialized
  static double* values=0;//randomly generated values for refractive index
  static double x0=0;  //the coordinate of the left end of the cavity
  static int length=0;//number of entries in the array
  if (!init) {
    if (values!=0) delete values;
    length=(int) (parameters.read("length_x","nm")*units.nm/units.length*2/xi)+1;
    values=new double[length];
    for (int i=0;i<length;i++) values[i]=((n2Max-n2Min)*generuj()+n2Min);
    init=true;
    x0=x;
  }
  int index=(int) floor((x-x0)/xi);
  double res=(x-x0)/xi-index;
  //res=0;
  if (index>=length-1) index=length-2;
  if (index<0) index=0;
  return((content(values[index]*(1-res)+values[index+1]*res) + 0.001*content(0,1))/sqrt(LP.n2)/sqrt(LP.n2));
}


int main()
{
  parameters.set("precision",0,"-"); 
  parameters.set("grid_step",1,"-");
  parameters.set("lv",2.0,"nm");
  parameters.set("phi_d",0.0,"pi");
  parameters.set("d",0.0,"nm");
  
  parameters.set("eig",60,"-");   //number of states in CF bases for diagonalization

  parameters.set("geometry",2.0,"-");
  parameters.set("flat_x",1000,"nm");
  parameters.set("flat_y",parameters.read("lv","nm"),"nm");
  parameters.set("sets",1,"-");
  
  LP.gamma_perp=1e+14;
  LP.gamma_par=1e+8;
  LP.g=25*units.Debye;
  LP.n2=1.5*1.5;
  LP.na=0.0025*units.Na;
  LP.D0=LP.na;
  LP.omegaa=2*M_PI*units.c/(640*units.nm);
  //LP.gamma_perp=LP.omegaa;
  LP.FFA=false;
  LP.with_shifts=true;
  LP.TMsorting=0;
  
  const char* name[]={"lasing_map", "TM_scan", "TM_vals", "lasing_mode_iter", "CF_kms", "CFs", "TLMs", "Dhx", "lasing_states"};
  output_class output(9,"data/laser/",name,"log-laser");
  FILE* stat=fopen("data/laser/lasing_map_all","w");
  //FILE* Dhx=fopen("data/laser/effective_pump_all","w");


  for (int i=10;i>=0;i--) {
    output.set_appendix("Nlm",i);
    output.clear_files();                            //clear the output files
  }
  output.open_files();

//for (prec p1=50;p1<50.1;p1+=1) {
  //parameters.set("eig",p1,"-");
  //parameters.set("flat_y",parameters.read("lv","nm"),"nm");
  LP.init();
  laser_prototype laser((int) parameters.read("eig","-"));
    
  //laser.scan_CF(100, 5000, 200, output.file(4));
  //exit(0);
  
  //initialize
  #if LASER > 0 
  fprintf(logg,"constructing the spatial grid and CF basis set\n");
  #endif
  laser.construct_grid();
  double length=LP.r->true_length();
  //LP.omegaa=50/length/units.length*units.c;
  fprintf(logg, "semiclassics regime: k R (2*flat_x) = %e\n", LP.omegaa/units.c*length*units.length);
  printf("omegaa=%e, gamma_perp=%e\n",LP.omegaa, LP.gamma_perp);
  //original values
  LP.gamma_perp=4.0*units.c/(units.length*length);  
  LP.gamma_par=0.001*units.c/(units.length*length);
  LP.omegaa=60*units.c/(units.length*length);
  //tuned values
  LP.gamma_perp=0.4*units.c/(units.length*length);
  LP.gamma_par=0.001*units.c/(units.length*length);
  LP.omegaa=(60-4*0.2)*units.c/(units.length*length);
  printf("omegaa=%e, gamma_perp=%e\n",LP.omegaa, LP.gamma_perp);
  //exit(0);
  
  //xi=p1*units.nm/units.length;
  //init=0;
  
  /*
  //comparison of the arpack and lapack diagonalization speed
  content *vec=new content[LP.dim*LP.dim*2];
  content *val=new content[LP.dim+1];
  ret_from_arpack_zn vysl;
  vysl.eigenvals=val;
  vysl.eigenvecs=vec;
  int eig=parameters.read("eig","-");
  for (int i=2*eig;i<2*eig+1;i+=2) {
    arpack_zn arpack(LP.dim, eig , -1e-12, 100000, true, (char*) "SM",i);
    char *msg=new char[100];
    sprintf(msg,"arpack nev=%i, ncv=%i",eig,i);
    progress_bar_prototype progress_bar(msg);
    progress_bar.start(true);
    arpack.go_for_it(LP.dim, CF_operate, (void*) LP.r, &vysl, false);
    progress_bar.finished(true);
  }
  
  lapack_zgee_prototype lapack(LP.dim);
  content *H=new content[LP.dim*LP.dim];
  set_matrix(LP.dim, H);
  progress_bar_prototype progress_bar("lapack zheev");
  progress_bar.start(true);
  lapack.diagonalize(H, val, vec, vec+LP.dim*LP.dim);
  progress_bar.finished(true);
  exit(0);*/
  //randomize(40);
for (int run=0;run<1;run++) {  
  double scan_extent=15.0, scan_step=1.0/200, basis_resolution=0.0/1000;
  LP.FFA=true;
  LP.TMsorting=0;
  do {//if diagonalization fails, take another disorder configuration
    laser.construct_Hamiltonian("data/laser/refractive_index_saved.txt");
    if (laser.construct_CF_grid(scan_extent, scan_step, basis_resolution, output.file(4))) break;
    delete LP.CF_pool;
    fprintf(logg,"diagonalization failed, taking another disorder configuration\n");
  } while (true);
  laser.Nlm=0;
  laser.update_total_profile();
  for (int i=0;i<LP.dim;i++) {
    double x,y;
    LP.r->s2coordinate(i,x,y);
    content n=sqrt(LP.n2) * LP.nx[i];
    //fprintf(output.file(7),"x=%+.8e [nm], n(x)=%+.8e%+.8ei\n", x*units.length/units.nm, n.real(), n.imag() );
    fprintf(output.file(7),"%+.8e %+.8e %+.8e\n", x*units.length/units.nm, n.real(), n.imag() );
  }
  CF_basis_prototype* CF_basis=laser.CF_pool->give_CF_basis(LP.omegaa/units.c/units.wavevector);
  CF_basis->plot(output.file(5));
  double ipr, x0, iprsum=0;
  for (int i=0;i<LP.Nbasis;i++) {
    LP.r->ipr_and_com(CF_basis->psi[i], ipr, x0);
    iprsum+=ipr;
  }
  printf("mean CF ipr:%e\n", iprsum/LP.Nbasis*units.length/units.nm); 
  //fprintf(output.file(0),"%e %e\n", p1, iprsum/LP.N*units.length/units.nm); 
  //laser.scan_TM(8, output.file(1), output.file(6));
  
  //fflush(0);

  //laser.destruct_grid();
  //continue;
  //exit(0);
//}
  /*for (int b=0;b<100;b++) {
    content self=0, off=0;
    prec k=laser.CF_basis[b]->k;
    fprintf(output.file(5),"%e ",k*units.wavevector*parameters.read("length_x","nm")*units.nm);
    for (int i=0;i<LP.N;i++) {
      for (int j=0;j<LP.N;j++) {
        content c=abs(laser.CF_basis[49]->correlation2(i,j,laser.CF_basis[b]));
        if (i==j) self+=c; else off+=c;
        //if (i==0) fprintf(output.file(5),"%e ",abs(c));
      }
    }
    //fprintf(output.file(5),"\n");
    self/=LP.N;
    off/=LP.N*(LP.N-1);
    fprintf(output.file(5),"%e %e\n", abs(self), abs(off));
  }
  exit(0);
  
  prec hx;
  LP.r->give_par(hx,hx);
  prec* maxc=new prec[LP.N*LP.N];
  for (int dx=0;dx<LP.dim/2;dx++) {
    content self=0, off=0;
    fprintf(output.file(5),"%e ",dx*hx*units.length/units.nm);
    for (int i=0;i<LP.N;i++) {
      for (int j=0;j<LP.N;j++) {
        content c=laser.CF_basis[0]->correlation1(i,j,dx);
        if (i==j) self+=c; else off+=c;
        if (abs(c)>maxc[i*LP.N+j]) maxc[i*LP.N+j]=abs(c);
        if (i==0) fprintf(output.file(5),"%e ",abs(c));
      }
    }
    fprintf(output.file(5),"\n");
    self/=LP.N;
    off/=LP.N*(LP.N-1);
    //fprintf(output.file(5),"%e %e\n", abs(self), abs(off));
  }
  //for (int i=0;i<LP.N;i++) fprintf(output.file(5),"%i %e\n", i, maxc[i*LP.N+0]);  
  delete maxc;
  exit(0);
  */
  
  bool plot_at_pump_flag[5]={1,1,1,1,1};
  double plot_at_pump[5]={0.01, 0.1, 1, 10, 100};
  
  //main iteration cycle
  int step1=0, step2;
  flag_prototype* flag=0;
  bool new_state=true;
  do {
    printf("main loop: at pump %e with %i lasing modes\n", LP.D, laser.Nlm);
    step1++;
    #if LASER > 0 
    fprintf(logg,"(re)starting the main cycle with %i lasing modes: step %i\n", laser.Nlm, step1);
    #endif
    
    for (int i=0;i<5;i++) {
      if (LP.D>plot_at_pump[i] && plot_at_pump_flag[i]) {    
        fprintf(output.file(8),"\n\n");
        laser.plot_lasing_modes(output.file(8)); 
        plot_at_pump_flag[i]=false;
      }
    }
    
    if (flag!=0) delete flag;
    flag=new flag_prototype[laser.Nlm];
    for (int i=0;i<laser.Nlm;i++) flag[i].reset();
    
    step2=0;
    int hmn;
    double k0[20]={0}, p0[20]={0};
    progress_bar_prototype progress_bar("converging at pump");
    progress_bar.start();
    LP.perturb_fail=0;
    do {//iterate at new pump value
      if (laser.Nlm==0) break;
      int which=laser.iterate_at_pump(flag, step2);
      hmn=0;
      for (int i=0;i<laser.Nlm;i++)  if (flag[i].converged) hmn++;

      #if LASER > 1
      fprintf(logg,"\titeration at a given pump value taken: step %i, %i/%i converged.\nMode number-frequency-modal power-norm: ", step2, hmn, laser.Nlm);
      for (int i=0;i<laser.Nlm;i++) fprintf(logg,"(%i,%+.2f,%.2e,%.2e) ", i, LP.x(laser.lasing_mode[i]->k), laser.lasing_mode[i]->modal_power, laser.lasing_mode[i]->norm);
      fprintf(logg,"\n");
      #endif
      laser.update_total_profile();

      //aux output
      bool headers=true;
      const char* string; 
      char c='p';
      if (flag[which].iterative) c='i';
      if (headers) string="step=%3i conv=%i state=%i method=%c diff=%e "; else string="%i %i %i %c %e ";
      if (step2!=0) fprintf(output.file(3),string, step2, hmn, which, c, abs(flag[which].diff));
      for (int i=0;i<laser.Nlm;i++) {
        double k=laser.lasing_mode[i]->k;
        double p=laser.lasing_mode[i]->modal_power;
        if (step2!=0) {
          if (headers) string="dk=%+.2e k=%+.2e dp=%+.2e p=%+.2e diff.R=%+.2e diff.I=%+.2e "; else string="%e %e %e %e %e %e ";
          if (i==which) fprintf(output.file(3),string, (k-k0[i])/k0[i], LP.x(k), (p-p0[i])/p0[i], p, flag[i].diff.real(), flag[i].diff.imag());
          else fprintf(output.file(3),string, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        k0[i]=k;p0[i]=p;
      }
      if (step2!=0) {fprintf(output.file(3),"\n"); fflush(0);}
      //end of aux output
      progress_bar.add(1,1);
      //if (step2 % 1000 ==0) {laser.scan_TM(8, output.file(1), output.file(6));}

    } while (hmn!=laser.Nlm && step2++<2000);
    int time=progress_bar.finished();
    if (step2 % 300==0 && step2>0) laser.converge_one_mode(laser.select_one_mode());
    #if LASER > 0 
    double tot=0;
    for (int i=0;i<laser.Nlm;i++) tot+=abs(flag[i].diff);
    fprintf(logg,"main cycle: converged at pump after %i steps with total residuum %.5e\n", step2, tot);
    #endif
    if ((laser.check_replicas() || step2>1999)) {
      if (LP.Dinc>0.01) {
	if (++LP.strategy==3) LP.strategy=0;
	fprintf(logg,"replica(s) found  (or not converged) at pump %e, changing strategy and restarting from snapshot\n", LP.D);
	LP.D=(LP.D+laser.snapshot->load(laser))/2;
	LP.Dinc/=2;
	continue;
      }
      else if (step2<=1999) laser.check_replicas(true); //continue anyway, there is no other option left
    }

    if (new_state) laser.scan_TM(max(8,laser.Nlm+8), output.file(1), output.file(6));
    else laser.scan_TM((laser.Nlm+8));
    double D, knew, Dnext, Dnextnext;
    int next_mode=laser.next_lasing_mode_pump(knew, Dnext, Dnextnext, output.file(2));
    if (Dnextnext<LP.D || (Dnext<laser.snapshot->D && LP.Dinc>0.01)) { //we failed, restart with smaller pump increase
      fprintf(logg,"modes shaked at pump %e (Dnextnext=%e, Dnext=%e, snapshot.D=%e, restarting from snapshot (turned off at the moment))\n", LP.D, Dnextnext, Dnext, laser.snapshot->D);
      //LP.D=(LP.D+laser.snapshot->load(laser))/2;
      LP.Dinc/=2;
      //continue;
    }
    
    new_state=laser.update_pump(D, Dnext, Dnextnext);
    if (new_state) {//new lasing state found
      LP.D=D;
      #if LASER > 0 
      fprintf(logg,"new lasing state (%i in TM_real) found at frequency %e and pump %e\n",next_mode, LP.x(knew), Dnext);
      #endif
      laser.lasing_mode[laser.Nlm]=new lasing_mode_prototype(LP.dim, LP.Nbasis, knew);
      laser.lasing_mode[laser.Nlm]->Dth=Dnext;
      laser.lasing_mode[laser.Nlm]->init_guess(laser.TM_eigreal[next_mode].real());
      laser.Nlm++;
      laser.update_total_profile();
      
       #if LASER > 1
      laser.lasing_mode[laser.Nlm-1]->TM->construct(laser.lasing_mode[laser.Nlm-1]->CF_basis);
      int N=LP.Nbasis;
      content *A=laser.lasing_mode[laser.Nlm-1]->A;
      content res=0, norm=0;
      for (int i=0;i<N;i++) {
        content rhs=A[i]/LP.D;
        norm+=A[i]*conj(A[i]);
        for (int j=0;j<N;j++) rhs-=laser.lasing_mode[laser.Nlm-1]->TM->M[i*N+j]*A[j];
        res+=rhs*A[i];
      }
      fprintf(logg,"\t\t\tcheck of iterative estimate: rel. residuum is %.3e\n", abs(res/norm));
      #endif
      output.set_appendix("Nlm",laser.Nlm);
      output.open_files();
    } else {
      laser.snapshot->save(laser);
      if (step2<500) LP.Dinc*=1.1; else LP.Dinc/=1.1;
      if (LP.Dinc>10) LP.Dinc=10;
      
      //closeness to the CFs
      for (int i=0;i<laser.Nlm;i++) {
	int jmax=0;
	prec asum=0, amax=0;
	for (int j=0;j<LP.Nbasis;j++) {
	  double a=((laser.lasing_mode[i]->A[j])*conj(laser.lasing_mode[i]->A[j])).real();
	  asum+=a;
	  if (a>amax) {amax=a;jmax=j;}
	}
	fprintf(logg,"lasing state %i: |a_max|^2 / sum |a|^2 = %e (largest CF #%i)\n",i,amax/asum, jmax);
      }
          
      //smoothening out of the refractive index?
      double n=0, n2=0, h=0, h2=0, nh=0;
      for (int i=0;i<LP.dim;i++) {
	double naux=LP.n2x[i].real();
	double haux=LP.D*LP.d[i].real()/LP.total_profile[i];
	n+=naux; n2+=naux*naux;
	h+=haux; h2+=haux*haux;
	nh+=naux*haux;
      fprintf(output.file(7),"%.4e %4i %e\n", LP.D, i, haux); 
      }
      double nvar =LP.dim*n2/(n*n)-1;
      double hvar =LP.dim*h2/(h*h)-1;
      double covar=LP.dim*nh/(n*h)-1;
	
      //estimate lasing modes by TLM approx
      lasing_state_char_prototype lsc(20);
      TLM_approx_prototype TLM_approx(laser);
      lsc.Dth[0]=0;
      int Nlm_estimate=TLM_approx.predict(lsc);
      fprintf(logg,"the TLM approx predicts %i lasing states with thresholds:",Nlm_estimate);
      //sprintf(buffer,"%e %i %i %i %i %e %e ", LP.D, step2, time, laser.Nlm, Nlm_estimate, hvar, covar);
      sprintf(buffer,"%e %e %e %i %i %e %e ", LP.D, h/LP.dim, nvar, laser.Nlm, Nlm_estimate, hvar, covar);
      fprintf(output.file(0),"%s",buffer);
      fprintf(stat,"%s",buffer);
      
      for (int i=0;i<16;i++) {
	if (i<Nlm_estimate) {
          fprintf(logg, "%e ",lsc.Dth[i]);
          fprintf(output.file(0), "%e ", lsc.Dth[i]);
          fprintf(stat, "%e ", lsc.Dth[i]);
          if (i<laser.Nlm) {
	    //the mode is lasing - one more free spot at the end
	    sprintf(buffer,"%e %e %e %e %e %i %e %i %e ", LP.x(lsc.k[i]), lsc.ipr[i], lsc.x0[i], laser.lasing_mode[i]->modal_power, lsc.o[0][i], (int) lsc.oi[0][i], lsc.o[1][i], (int) lsc.oi[1][i], 0.0);
	    fprintf(output.file(0), "%s", buffer);
	    fprintf(stat, "%s", buffer);
	    
	    double ipr, x0;
	    LP.r->ipr_and_com(laser.lasing_mode[i]->amplitude, ipr, x0);
	    ipr*=units.length/units.nm;
	    x0*=units.length/units.nm;
	    if (abs(ipr-lsc.ipr[i])+abs(x0-lsc.x0[i]) > 1e-8) fprintf(logg,"WARNING: lm %i: TLM (ipr=%e, x0=%e) vs direct (ipr=%e, x0=%e) -> diff=(%e, %e)\n", i, lsc.ipr[i], lsc.x0[i], ipr, x0, lsc.ipr[i]-ipr, lsc.x0[i]-x0);
	  }
          else {
	    sprintf(buffer,"%e %e %e %e %e %i %e %i %e ",0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0.0);
	    fprintf(output.file(0), "%s", buffer);
	    fprintf(stat, "%s", buffer);
	  }
	}
	else {
	  sprintf(buffer,"%e %e %e %e %e %e %i %e %i %e ",0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0.0);
	  fprintf(output.file(0), "%s", buffer);
	  fprintf(stat, "%s", buffer);
	}
      }
      fprintf(logg, "\n");
      fprintf(output.file(0), "\n");
      fprintf(stat, "\n");
      LP.D=D;
      //exit(0);
    }
    laser.CF_pool->refresh_undeletable();
    
#if LASER > 0 
    fprintf(logg,"main cycle: pump set to %e\n",LP.D);
#endif

  } while (LP.D<110);
  
}
}


//}


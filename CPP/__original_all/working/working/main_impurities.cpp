#include "main.h"

//temperature tweak: under in classes.cpp
//multisequence iterations min ->1 instead of 3

//FINAL VERSION: 
//change the g-factor to almost zero, check its sign for ZnTe


prec V;   //volume of the grid cell
prec v0;  //density of cations
prec Eexch;//the exchnage energy with which the discretized exchange hamiltonian can be written as 
//H_kl = \delta_kl exch \sigma_z * I_z / I
//exch = (5/2) (3/2) (beta/3) (1/V) 
prec B2E[4]; //external magnetic field coupled to spin operator of impurity in energy units (g_Mn muB B): B2E[0-2] is the vector, B2E[3] is the vector magnitude
prec kBT; //temperature in energy units
prec e2T; //energy of a fully aligned impurity spin divided by temperature (5/2 3/2 beta/3 1/V) / kBT
two_electron_state_prototype* state2e=0;
region* r1;
int* s2un2e, * un2s2e;
output_class* output;
double ES, ET0, ETmaxshift; //the unperturbed energy, unperturbed singlet, triplet, maximal T shift
double non_mon_thresh=1e-3; //threshold to declare false if error increases
int multisequence_iterations_min=3, multisequence_iterations_max=20;//limits for number of sequences
int iterationsmax=100;	    //limit on number of iterations within a single sequence
int sector, state_in_sector=0;//identifications of the state we are interested in (S=0/0, T+-=+-1/0, T0=0/1)
prec R, width; //doping radius: Mn spin is modulated by a factor 1/(1+exp((r-R)/width))
               //width>0 amounts to overdoped, width<0 amounts to underdoped

bool left_only=false, homog=true, bistable=false, B_in_impurities=false;
int patterns_output=0;
int init_pattern=1, init_pattern_min=1, init_pattern_max=1;//how to initialize the patterns

//returns the lowest spin polarized / unpolarized f-state depending on subsector =  0, +1, -1
int subsector_ground_state(int subsector, int which)
{
  int Nf=parameters.read("eig2e","-");
  int* pos=new int[Nf];
  for (int i=0;i<Nf;i++) pos[i]=i;
  for (int i=0;i<Nf;i++) for (int j=0;j<Nf-1;j++) {
    if (state2e->read('f',pos[j],"energy")>state2e->read('f',pos[j+1],"energy")) {
      int aux=pos[j]; pos[j]=pos[j+1]; pos[j+1]=aux;
    }
  }  
  
  prec S2, SB;
  int state=0, pntr=-1;
  do {
    S2=state2e->read('f',pos[state], "S2");
    SB=state2e->read('f',pos[state], "SB");
    if (abs(SB-subsector)<0.1) pntr++;
    fprintf(logg, "f-state %i with S2=%e SB=%e results in subsector count %i (pntr)\n", state, S2, SB, pntr);
    if (++state>=(int) parameters.read("eig2e","-")) {printf("2e state in subsector %i not found\n", subsector); exit(1);}
  } while (pntr<which);
  fprintf(logg, "state %i in subsector %i is the f-state %i\n", pntr, subsector, pos[state-1]);
  return(pos[state-1]);
}


double E2l0(double E)
{
  double x = units.hbar*units.hbar/(E*units.meV*parameters.read("m_eff","m_e")*units.m_e);
  double l0=sqrt(x)/units.nm;
  //printf("conversion from energy %e meV to confinement length %e nm\n",E, l0);
  return(l0);
}

struct disp_stat {
  prec *values;
  prec sum, sum2, ave, disp, error;
  int N, length;
  disp_stat(int length_) {
    length=length_;
    values=new prec[length];
    sum=sum2=N=0;
  }
  ~disp_stat() {delete values;}
  void new_value(prec x) {
    for (int i=length-1;i>0;i--) values[i]=values[i-1];
    values[0]=x;
    if (N<length) N++;
    sum=sum2=0;
    for (int i=0;i<N;i++) {sum+=values[i]; sum2+=values[i]*values[i];}
    ave=sum/N;
    disp=max(0,sum2/N-ave*ave);
    disp=sqrt(disp);
    error=abs(disp/ave);
    if (ave==0) error=0;
  }
};

extern "C" void F77NAME(dgels) (char*, int*, int*, int *, double*, int*, double*, int*, double*, int*,int*);

struct least_square_prototype {
  prec T0, estimate, error;
  prec * x, *y, *A, *b, *WORK;
  int N, occurences;
  
  least_square_prototype(prec T0_) {
    T0=T0_;
    estimate=error=N=occurences=0;
    x=new prec[10*100];
    y=x+100;
    b=y+100;
    A=b+100;
    WORK=A+300;
  }
  ~least_square_prototype() {delete x;}
  
  void new_value(prec x_new, prec y_new, int occurences_now) {
    occurences+=occurences_now;
    x[N]=x_new;
    y[N++]=y_new;
    for (int i=0;i<N;i++) for (int k=0;k<3;k++) A[i*3+k]=pow(x[i],k);
    for (int i=0;i<N;i++) b[i]=y[i];
    char TRANS='T';
    int M=3, NRHS=1, LWORK=300, INFO, LDB=max(M,N);
    F77NAME(dgels)(&TRANS, &M, &N, &NRHS, A, &M, b, &LDB, WORK, &LWORK, &INFO);
    if (INFO!=0) {printf("LLS error: routine dgels returned %i\n",INFO); exit(1);}
    prec new_estimate=0;
    fprintf(logg,"LLS predicts koefs: b[0]=%e, b[1]=%e, b[2]=%e\n",b[0],b[1],b[2]);
    for (int k=0;k<3;k++) new_estimate+=b[k]*pow(T0,k);
    if (new_estimate==0) error=1; else error=abs(1-estimate/new_estimate);
    if (N==1) error=0;
    estimate=new_estimate;
  }
};
  

void print_pattern(FILE* file)
{
  for (int n=0;n<r1->Nw;n++) {
    prec x, y, s[4];
    spin_impurities.position(n,x,y);
    state2e->spinatn(s, n);
    fprintf(file,"%.4e %.4e %.4e %.4e %.4e %.4e\n", x*units.length/units.nm, y*units.length/units.nm, spin_impurities.I[2][n], s[3]*e2T, s[3]*(1.0/2)/V, s[0]/V);
  }
  fprintf(file,"\n\n");
}

void update_constants()
{
  prec hx,hy,width_z=parameters.read("width_z","nm");
  r1->give_par(hx,hy);
  hx*=units.length/units.nm;
  hy*=units.length/units.nm;
  V=hx*hy*width_z;
  v0=pow(units.aGaAs/units.nm,3)/4;
  //fprintf(logg,"V=%e v0=%e V/v0=%e\n",V, v0, V/v0);
  
  Eexch=-(5.0/2)*(3.0/2)/V*(parameters.read("exch","meV")*units.meV*v0/3)/units.energy;
  B2E[3]=2*units.muB*parameters.read("B","Tesla")/units.energy;//g(Mn)=2 (minus sign absorbed into muB), no 5/2!!!
  prec phi_B=parameters.read("phi_B", "pi")*M_PI;
  prec theta_B=parameters.read("theta_B", "pi")*M_PI;
  B2E[0]=B2E[3]*cos(phi_B)*sin(theta_B);
  B2E[1]=B2E[3]*sin(phi_B)*sin(theta_B);
  B2E[2]=B2E[3]*cos(theta_B);
  kBT=parameters.read("T","K")*units.kB/units.energy;
  e2T=abs(Eexch/kBT);
  if (kBT==0) e2T=0;
  ETmaxshift=parameters.read("exch","meV")*parameters.read("xMn", "-")*(5.0/2)*(1.0/2)*2;
}

void electron_1e_diagonalization(state_info* &states)
{
  for (int i=0;i<2;i++) {
    if (i==0) spin_impurities.reinitialize(0);
    if (i==1) spin_impurities.reinitialize(states->r1,6);
    parameters.recompute();
    if (states!=0) clean_up(states);
    states=diagonalize();

    //some info about the single electron states into log and an output file
    states->show_set("clear");
    states->show_set("energy");
    statistics_multi(*states, 0, 2+8, output->file(1),2);

    states->show_set("Ix");states->show_set("Iy");
    states->show_set("I");
    states->show_set("s");states->show_set("lz");
    states->show_set("sx");states->show_set("sy");
    states->show_set("sz");
    statistics_multi(*states, 0, 1+2+4+8+16, logg,4);
  }
}

void electron_diagonalization(state_info* &states, sorting_class& sorting2e)
{
  //!the single electron diagonalization
  prec B=parameters.read("B","Tesla");
  if (parameters.read("d","nm")==0 && B==0) parameters.set("B",0.01, "Tesla");
  randomize();
  spin_impurities.reinitialize(0);
  parameters.recompute();
  if (states!=0) clean_up(states);
  states=diagonalize();
  parameters.set("B",B, "Tesla");
  r1=states->r1;
  
  //message(logg,"sorting...",1);
  //sorting1e.inspect_new_1e_states(states);
  //sorting1e.sort(0, states->sorted2un, states->unsorted2s);
  //message(logg,"...done\n",1);
  
  //some info about the single electron states into log and an output file
  states->show_set("clear");
  states->show_set("energy");
  statistics_multi(*states, 0, 2+8, output->file(1),2);
  
  states->show_set("Ix");states->show_set("Iy");
  states->show_set("I");
  states->show_set("s");states->show_set("lz");
  //states->show_set("sx");states->show_set("sy");
  //states->show_set("sz");
  statistics_multi(*states, 0, 1+2+4+8+16, logg,4);
  
  if (B==0) parameters.set("B",0.01, "Tesla");
  if (state2e!=0) delete state2e;
  state2e=new two_electron_state_prototype(states, &sorting2e);
  state2e->diagonalize(0, s2un2e, un2s2e);
  parameters.set("B",B, "Tesla");
  state2e->statistics("b",1+2+4,logg,4,20);
  state2e->statistics("f",1+2+4,logg,4);
  state2e->statistics("f",1+2+4,output->file(2),2);
  ES =state2e->read('f',subsector_ground_state(0,0),"energy");
  ET0=state2e->read('f',subsector_ground_state(0,1),"energy");
  prec ETp=state2e->read('f',subsector_ground_state(1,0),"energy");
  prec ETm=state2e->read('f',subsector_ground_state(-1,0),"energy");
  fprintf(logg, "lowest state energies, singlet: %e triplets: '+'=%e,'0'=%e,'-'=%e giving exchange %e\n", ES, ETp, ET0, ETm, ET0-ES);
  state2e->reducedDM(0,0);
}

//update the impurity spins according to the hole spin density
void align_impurities(state_info* state=0)
{
  prec norm2e=0;
  for (int n=0;n<r1->Nw;n++) {
    prec s[4];
    prec x,y;
    r1->s2coordinate(n,x,y);
    prec r=sqrt(x*x+y*y)*units.length/units.nm;
    if (state==0) state2e->spinatn(s,n); else {
      spinor psi=spinor(state->vysl->eigenvecs[n],state->vysl->eigenvecs[n+r1->Nw]);
      for (int i=0;i<4;i++) s[i]=(psi.braket(psi.SigmaTimesMe(pauli+i*4))).real();
    }
    s[1]=s[2]=0;
    prec modul=1;
    if (!homog) {
      if (R!=0) {
	if (width!=0) modul/=1+exp(min((r-R)/width,100));
	else {if (r>R) modul=0; else modul=1;}
      }
      if (left_only && x>0*units.nm/units.length) modul=0;
    }
    for (int i=0;i<3;i++) s[i+1]=Eexch*s[i+1]+B2E[i]*B_in_impurities;
    spin_impurities.set_impurity(n,s+1,kBT,modul);         
    norm2e+=s[0];
  }
  if (abs(norm2e-2)>1e-8) fprintf(logg,"2particle state norm=%e\n", norm2e);
}

prec iteration_step()
{
  //diagonalize the two hole system
  state2e->diagonalize(0, s2un2e, un2s2e,6);
  state2e->statistics("b",1+2+4,logg,4,20);
  state2e->statistics("f",1+2+4,logg,4);        
  int state=subsector_ground_state(sector, state_in_sector);
  state2e->reducedDM(state,state);
  state2e->distribution_characteristics(5);
  
  //update mpurities, calculate and output ensemble characteristics
  align_impurities();
  spin_impurities.distribution_characteristics(state2e);
  if (patterns_output>1) print_pattern(output->file(5));
  
  //return the energy shift
  double Eunperturbed;
  if (sector==0) Eunperturbed=ES; else Eunperturbed=ET0;
  return(state2e->read('f',state,"energy")-Eunperturbed);
}

bool iteration_sequence(int& stepmax, prec &E)
{
  bool success=false;
  if (++init_pattern>init_pattern_max) init_pattern=init_pattern_min;
  if (sector==0) spin_impurities.reinitialize(r1,init_pattern,R);
  if (sector==1) spin_impurities.reinitialize(r1,8,R);
  if (sector==-1) spin_impurities.reinitialize(r1,9,R);
  update_constants();
  
  //output statistics from distributions at the iteration beginning
  spin_impurities.distribution_characteristics(state2e);
  int state=subsector_ground_state(sector, state_in_sector);
  state2e->reducedDM(state,state);
  state2e->distribution_characteristics(5);
  fprintf(logg,"iter_seq: Rmax=%e\n", state2e->Rmax);
  fprintf(output->file(4),"%2i %+e %+e %+e %+e %+e %+e %2i %+e %+e %+e %+e %2i %2i\n", 0, spin_impurities.ave[2], spin_impurities.d[2], spin_impurities.dphi[2], state2e->Rmax, abs(state2e->M[2]/state2e->M[0]), state2e->Qsph, spin_impurities.borders, 1.0, 0.0, 0.0, 0.0, 0, 0);
  fflush(output->file(4));
  if (patterns_output>1) print_pattern(output->file(5));
  
  if (parameters.read("xMn","-")==0) {
    E=0;
    stepmax=1;
    return(true);
  }
  
  int step=1, decreases_error=0, decreases_update=0, steady_pattern=0;
  disp_stat disp(5);
  prec prev_error=0, valuem1=0, valuem2=0; 
  do {//the loop
  
    //make one iteration
    E=iteration_step();
    disp.new_value(E);
    if (disp.error<prev_error) decreases_error++; else decreases_error=0;
    prev_error=disp.error;
    double prev_update=valuem1-valuem2, update=E-valuem1;
    valuem2=valuem1;valuem1=E;
    if (abs(prev_update)>abs(update) || prev_update*update>0) decreases_update++; else decreases_update=0;
    
    //output the result of the iteration step
    fprintf(output->file(4),"%2i %+e %+e %+e %+e %+e %+e %2i %+e %+e %+e %+e %2i %2i\n", step, spin_impurities.ave[2], spin_impurities.d[2], spin_impurities.dphi[2], state2e->Rmax, abs(state2e->M[2]/state2e->M[0]), state2e->Qsph, spin_impurities.borders, 1-spin_impurities.overlap, E, disp.disp, disp.error, decreases_error, decreases_update);
    fflush(output->file(4));

    //the inner loop end conditions
    
    //if we become non-monotonic after several steps and there is no underrelaxation, declare fail
    if ((decreases_error==0 || decreases_update==0) && (step>8) && spin_impurities.under==0)
	if (disp.error>non_mon_thresh && (decreases_error+decreases_update)<4) break;
    
    //if the spin impurities patter did not changed for the last 1 step(s)
    if (abs(spin_impurities.overlap-1)<1e-20) steady_pattern++; else steady_pattern=0;
    if (steady_pattern==1) {success=true;break;}
    
    //if the dispersion error is below a small threshold and there were at least few steps
    if (disp.error<1e-5 && step>5) {success=true;break;}
    
    //if the dispersion error is below a larger threshold but there is a monotonous trend in updates
    if (disp.error<1e-4 && decreases_update>10) {success=true;break;}
    
    //if the value dropped very small, accept zero
    if (abs(E)<1e-6 && decreases_update>10) {E=0; success=true;break;}
      
    //if you reached the end or there is a very long monotonous trend
    if (decreases_update>90 || step==stepmax) {
      //if the energy is much larger than the pattern change, accept the energy
      if (abs(spin_impurities.overlap-1)/abs(E)<1) success=true;
      //if the energy is much smaller than he pattern change, accept zero
      if (abs(spin_impurities.overlap-1)/abs(E)>1) {E=0;success=true;}
    }
  } while (step++<stepmax);
  stepmax=step;
  fprintf(output->file(4),"\n");
  return(success);
}

//do several sequences, if all of the are convergent, declare "convergent"
//returns the number of taken sequences and the value statistics
bool multisequence (disp_stat** disp)
{
  int sequence=1;
  bool success=true;
  
  do {
    prec E;
    int iterations=iterationsmax;
    if (!iteration_sequence(iterations, E)) success=false;
    
    //statistics of the energy shift and pattern at the sequences end 
    disp[spin_impurities.pattern]->new_value(E);
    fprintf(logg,"multiseq: Rmax=%e\n", state2e->Rmax);
    fprintf(output->file(3),"%2i/%2i %+e %+e %+e %+e %+e %+e %2i %+e %+e %+e %2i\n", sequence, iterations, spin_impurities.ave[2], spin_impurities.d[2], spin_impurities.dphi[2], state2e->Rmax, abs(state2e->M[2]/state2e->M[0]), state2e->Qsph, spin_impurities.borders, 1-spin_impurities.overlap, E, disp[spin_impurities.pattern]->error, spin_impurities.pattern);
    fflush(0);
    //fprintf(output->file(4),"multisequence conditions: N=%i, min=%i, error=%e, sequence=%i, multi_max=%i, success=%i\n",disp[spin_impurities.pattern]->N, multisequence_iterations_min, disp[spin_impurities.pattern]->error, sequence, multisequence_iterations_max, success );
    if ((disp[spin_impurities.pattern]->N>=multisequence_iterations_min && disp[spin_impurities.pattern]->error<1e-3) || sequence==multisequence_iterations_max || success==false) break;
  } while (sequence++>-1);
  return(success);
}

void LLS_extrapolation(least_square_prototype** least_square)
{
  //increased temperature extrapolation
  prec T0=parameters.read("T","K");
  prec T=T0, step1=0.01, step2=0.001, prev=0;
  int valid=0, invalid=0, lowest=0, most_frequent=0, step=1;
  bool convergent=false;
  do {
    disp_stat* disp[4];
    for (int i=0;i<4;i++) disp[i]=new disp_stat(multisequence_iterations_max);
    bool convergent_now=multisequence(disp);
    for (int i=0;i<4;i++) {
      if (disp[i]->N>disp[most_frequent]->N) most_frequent=i;
      if (disp[i]->ave<disp[lowest]->ave) lowest=i;
    }
    if (convergent_now) {
      if (!convergent && patterns_output>0) print_pattern(output->file(6));//first convergent pattern
      convergent=true;
      if (step!=1 && abs(disp[most_frequent]->ave)>abs(prev)) fprintf(output->file(7),"refusing convergent value due to wrong T-behavior: prev=%e actual=%e\n", prev, disp[most_frequent]->ave);
      else {
        valid++;
        for (int i=0;i<4;i++) if (disp[i]->N>0) least_square[i]->new_value(T,disp[i]->ave, disp[i]->N);
        fprintf(output->file(7),"patterns spotted following number of times: no %ix, uniform %ix, dipole %ix, core-halo %ix\n", disp[0]->N, disp[1]->N, disp[2]->N, disp[3]->N);
        fprintf(output->file(7),"new value(s) into LLS (%e +- %e for most frequent one, being %i): actual estimate is %e (relative change from previous %e)\n", disp[most_frequent]->ave, disp[most_frequent]->disp, most_frequent, least_square[most_frequent]->estimate, least_square[most_frequent]->error);
      }
      prev=disp[most_frequent]->ave;
      //step==1 means it is convergent at initial temperature
      if (step==1 || (valid>5 && least_square[most_frequent]->error<1e-4) || valid>30) {
        if (step==1) {
          least_square[most_frequent]->error=disp[most_frequent]->disp/disp[most_frequent]->ave;
          if (disp[most_frequent]->ave==0) least_square[most_frequent]->error=0;
        }
        for (int i=0;i<4;i++) delete disp[i];
        break;
      }
    }
    for (int i=0;i<4;i++) delete disp[i];
    if (convergent && !convergent_now) invalid++; 
    if (convergent) T+=step2; else T+=step1;
    parameters.set("T",T,"K");
    for (int i=1;i<7;i++) fprintf(output->file(i),"temp rescale T=%f\n\n", T);
    fprintf(output->file(7),"convergence loop: step %i at T=%e convergent %i.  From so far valid=%i invalid=%i LLS predicts value %e (reliability %e)\n", step, T, convergent, valid, invalid, least_square[most_frequent]->estimate, least_square[most_frequent]->error);
  } while (step++<100);
  if (step==1) fprintf(output->file(7), "convergent at T0, no temperature fitting, error taken as the multisequence dispersion divided by the average\n"); 
  fprintf(output->file(7),"after %i steps and %i convergent values LLS predicts value %e (reliability %e)\n", step, valid, least_square[most_frequent]->estimate, least_square[most_frequent]->error);
  parameters.set("T",T0,"K");

}

void single_hole_part(prec* res, int to)
{
  if (to==0) return;
  parameters.set("sets",2, "-");
  dof.length=2;
  
  prec B=parameters.read("B","Tesla"), xMn=parameters.read("xMn","-");
  if (B==0) parameters.set("B",0.001, "Tesla");
  parameters.recompute();
  state_info* states=0;

  //sorting_class sorting1e;
  //sorting1e.init((int) parameters.read("eigmax","-"), (int) parameters.read("eig","-"),3);
 
  //!the single electron diagonalization
  for (int round=0;round<to;round++) {
    
    switch (round) {
      case 0 : {//retain the spin pattern
	states=diagonalize();
        res[0]=states->read(0,"energy");
        res[1]=states->read(1,"energy");
	res[2]=states->read(2,"energy");
        break;
      } 
      case 1 : {
        spin_impurities.reinitialize(r1,8-3,0);
	prec e=0, diff;
	int step=0;
	do { 
	  states=diagonalize();
	  align_impurities(states);
	  prec e_new=states->read(0,"energy");
	  diff=e-e_new;
          fprintf(logg,"single e diag loop step %i : e=%e, e_new=%e, diff=%e\n", step, e, e_new, diff);
	  e=e_new;
	  if (abs(diff)<1e-12 || step++==20) break;
	  clean_up(states);
	} while (true);
        fprintf(logg,"single e diag loop end: steps %i : e=%e, diff=%e\n", step, e, diff);
        res[3]=states->read(0,"energy");
	break;
      }
      case 2 : {
	parameters.set("xMn",0,"-");
        states=diagonalize();
        res[4]=states->read(0,"energy");
	break;
      }
    }
 
    //message(logg,"sorting...",1);
    //sorting1e.inspect_new_1e_states(states);
    //sorting1e.sort(0, states->sorted2un, states->unsorted2s);
    //message(logg,"...done\n",1);

    //some info about the single electron states into log and an output file
    states->show_set("clear");
    states->show_set("energy");
    //statistics_multi(*states1e, 0, 2+8, output.file(1),2);

    states->show_set("Ix");states->show_set("Iy");
    states->show_set("I");
    states->show_set("s");states->show_set("lz");
    statistics_multi(*states, 0, 1+2+4+8+16, logg,4);

    //the looping parameter into output file
    //sprintf(buffer,"%.14e ",p2);
    //splachni(output.file(0),buffer,2);
    if (states!=0) clean_up(states);
  }
      
  //step=sorting1e.adjust_step();
  parameters.set("sets",1, "-");
  dof.length=1;
  parameters.set("B",B, "Tesla");
  parameters.set("xMn",xMn,"-");

}


int main()
{
  //!GaAs parameters
  parameters.set("m_eff",0.2*0+0.106*0+0.21,"m_e");             //effective mass - used to be 1.0/(3.8+0.7)
  parameters.set("eps0",9.4+10.6*0,"-");               //static dielectric constant - used to be 12.9

  
  parameters.set("gx",-0.44/1e3,"-");               //g factors - used to be -0.44
  parameters.set("gy",-0.44/1e3,"-");
  parameters.set("gz",-0.44/1e3,"-");
  
  
  parameters.set("exch",-1050-1321.7*0,"meV");		//-1 eV (beta * N_0)
  parameters.set("xMn",0.0+0.00544*0,"-");       // x_Mn
  parameters.set("T", 0.0+17*0, "K");         //temperature
  parameters.set("width_z",1.83+1.0*0,"nm");     //width of the 2DEG	- used to be 3 nm
  
  //!quantum dot setup
  parameters.set("lv",E2l0(25*0+25.6711*0+37.7936*0+61.5921*0+30),"nm");			//confinement length - used to be 20 meV
  parameters.set("ALHO",0,"-");
  
  spin_impurities.under=1e-10;

  parameters.set("B",0.0,"Tesla");   		   //magnetic field
  parameters.set("phi_B",0.0,"pi");
  parameters.set("theta_B",0.0,"pi");

  parameters.set("d",0,"nm");   		//interdot distance
  parameters.set("phi_d",0.0,"pi");  
  parameters.set("ring",0,"nm"); 
  
  //!computational parameters
  parameters.set("eig",60/2,"-");    		//number of single electron states (including spin) to obtain in diagonalization
  parameters.set("eigmax",150,"-");      	//max. number of states in sorting for the previous - must be large enough to accomodate all crossings
  parameters.set("selected2e",21,"-");		//number of single electron states (no spin) to built the two electron basis; due to the spin it must be at least eig >= 2*selected, but if the states are sorted, it should be larger to accomodate all crossings
  parameters.set("eig2e",100-80*0,"-");    		//goal - number of the two electron states to get

  //!grids
  parameters.set("extent",5.8,"-");		//leave as is
  parameters.set("grid_step",0.1,"-"); 	//determines grid dimension and single electron functions precision
  parameters.set("CE_exp",4,"-");		//find out the rule, then I optimize
  parameters.set("sets",1,"-");     //2S+1
  parameters.set("precision",3,"-");
  parameters.set("system",'A',"-"); //GaAs

  dof.length=1;
  
  //!sorting
  parameters.set("threshold",0.2,"-");		//tolerated unprecision in sorting accoring the symmetries, value is ad-hoc, we will see
  parameters.set("sorting_method",2,"-");  	//(0 no sorting; 1 I; 2 Ix,Iy; 3 lz; 4 distance;  5 FD states; *-1 spin)

  const char* name[]={"measured-imp0","measured-1e","measured-f", "measured-imp1", "measured-imp2", "plot-imp", "plot-imp-mins","log-LLS"};

  int NF=8;
  output=new output_class(NF,"data/final/addons/",name,"log");     //SPOT
  //output->set_appendix("P");       //appendix to the output files
  
  output->clear_files();                            //clear the output files
  //output->open_files();
  //for (int i=0;i<7;i++) fprintf(output->file(i), "pars NEW: d set to 0.3 lv ");
  //output->close_files();
  
  output->open_files();    
  state_info *states=0;
  units.list(logg,4);
  
  //!initialize the sorting classes for single and double electron states
  sorting_class sorting2e;
  int Ns=(int) parameters.read("selected2e","-");
  s2un2e=new int [Ns*Ns];
  un2s2e=new int [Ns*Ns];
  for (int a=0;a<Ns*Ns;a++) un2s2e[a]=s2un2e[a]=a;
  sorting2e.init(Ns*Ns,Ns*Ns,1);
  bool diag=true;
  //TraceOn('z', 3, 2, 1, 1);
  
  //parameters.set("sets",2,"-");
  //electron_1e_diagonalization(states);
  //exit(0);

//SPOT-start   
  //!the outer loop
  prec step=5e-2;
  for (prec p1=1e-6;p1<1.1e-2;p1+=step) {
    //if (p1>10) {
      //parameters.set("lv",E2l0(61.5921),"nm");
      //parameters.set("eps0",9.4,"-");
    //}
    //if (p1>1.49 && p1<1.6) step=0.01; else step=0.1;
      //p1=0.0;
      //parameters.set("d",parameters.read("d","nm")+0.1*parameters.read("lv","nm"),"nm");diag=true;
    //}//step=10.1;

    //if (p1>10) step=10;
    //if (p1<3.5 || p1>4.5) iterationsmax=50; else iterationsmax=100;
    
    //parameters.set("lv",E2l0(p1),"nm");
    parameters.set("B",p1,"T");
    //parameters.set("d",(p1+0.0)*parameters.read("lv","nm"),"nm");
    diag=true;
    //parameters.set("T",p1*0+2,"K");
    //parameters.set("eps0",12.9/p1,"-");               //static dielectric constant
    //if (p1==0) parameters.set("eps0",0,"-");               //static dielectric constant
    //R=p1;//*parameters.read("lv","nm");
    //width=0.0;//*parameters.read("lv","nm");//minus = underdoped
    homog=true;
    left_only=false;
    patterns_output=1;
    bistable=true;
    B_in_impurities=false;
    init_pattern_min=5;
    init_pattern=init_pattern_max=7;
    
    fprintf(output->file(0),"\n");

    //!the inner loop
    for (prec p2=0.0;p2<0.011;p2+=10.002) {
    //for (prec p2=0.02;p2<0.021;p2+=0.001) {
      //if (p1==3.1 && p2<0.0165) continue;
 //SPOT-end
      
      parameters.set("xMn",p2,"-");      
      //parameters.set("T",p2,"K");      
      sprintf(buffer,"parameters set: p1=%f, p2=%f\n",p1,p2);
      splachni(logg,buffer,1+4);
      
      if (diag) {electron_diagonalization(states, sorting2e);diag=false;}
      update_constants();

      //the looping parameter into all output files
      //for (int i=1;i<NF;i++) fprintf(output->file(i),"\npars R=%6f (width=%3f)[nm], d=0.0 lv, x=%f T=%6f ",p1, width, p2, parameters.read("T","K"));
      for (int i=1;i<NF;i++) fprintf(output->file(i),"\npars B=%6f x=%f T=%6f ", p1, p2, parameters.read("T","K"));                            //SPOT
      //fprintf(output->file(0),"%e %e %e ",p1, p2, parameters.read("T","K"));
      //for (int i=1;i<NF;i++) fprintf(output->file(i),"\npars d=%6f lv, x=%f T=%6f ",p1, 0.01, parameters.read("T","K"));                            //SPOT
      fprintf(output->file(0),"%e %e %e ", p1, p2, parameters.read("T","K"));
      
      prec T0=parameters.read("T","K");
      least_square_prototype* least_square[8];
      for (int i=0;i<8;i++) least_square[i]=new least_square_prototype(T0);
      sector=1;
      //multisequence_iterations_min=6;
      multisequence_iterations_min=1;
      if ((!homog || T0!=0) || false) {             //SPOT
        for (int i=1;i<NF;i++) fprintf(output->file(i), "\n\n  TRIPLET\n");
        LLS_extrapolation(least_square);
      } else {
        least_square[1]->new_value(T0,ETmaxshift,1);
        int state=subsector_ground_state(sector, state_in_sector);
        state2e->reducedDM(state,state);
        state2e->distribution_characteristics(5);
      }
      prec resS[5], resT[5];
      for (int i=0;i<5;i++) resS[i]=resT[i]=0;
      single_hole_part(resT,0);

      int state=subsector_ground_state(sector, state_in_sector);
      content auxRR=state2e->rdotr(state,state);
      prec RR1T=auxRR.real(), RR2T=auxRR.imag()*pow(units.length/units.nm,2);
      prec RmaxT=state2e->Rmax;
      prec RSmaxT=state2e->RSmax;
      sector=0;
      for (int i=1;i<NF;i++) fprintf(output->file(i), "\n  SINGLET\n");
      if (bistable) multisequence_iterations_min=6;
      if (parameters.read("xMn","-")!=0) LLS_extrapolation(least_square+4);   //SPOT
      else {
        least_square[4+0]->new_value(T0,0,1);
        int state=subsector_ground_state(sector, state_in_sector);
        state2e->reducedDM(state,state);
        state2e->distribution_characteristics(5);
      }
      single_hole_part(resS,0);

      state=subsector_ground_state(sector, state_in_sector);
      auxRR=state2e->rdotr(state,state);
      prec RR1S=auxRR.real(), RR2S=auxRR.imag()*pow(units.length/units.nm,2);
      prec RmaxS=state2e->Rmax;
      prec RSmaxS=state2e->RSmax;

      //output results
      double Exch=ET0-ES, ExchMn=0;
      fprintf(output->file(0),"%e %e %e ", ES, ET0, ET0-ES);
      int most_frequent=0, lowest=0;
      
      for (int state=0;state<2;state++) {//0-T, 1-S
        for (int pattern=0;pattern<4;pattern++) {
          if (least_square[state*4+pattern]->occurences>least_square[state*4+most_frequent]->occurences) most_frequent=pattern;
          if (least_square[state*4+pattern]->estimate<=least_square[state*4+lowest]->estimate && least_square[state*4+pattern]->occurences>0) lowest=pattern;
        }
        int N1=least_square[state*4+most_frequent]->occurences, N2=least_square[state*4+lowest]->occurences;
        if (N1<1 || N2<1) {printf("no entry anywhere? (N=%i in most frequent %i and N=%i in lowest %i)\n", N1, most_frequent, N2, lowest); exit(1);}
        double T=least_square[state*4+lowest]->x[least_square[state*4+lowest]->N-1];
        double shift=least_square[state*4+lowest]->estimate;
        double error=least_square[state*4+lowest]->error;
        fprintf(output->file(0),"%e %2i %i %e %e ", T , N2, lowest, shift, error);
        if (state==0) ExchMn+=shift; else ExchMn-=shift;//shift(T)-shift(S)
      }
      int ground;
      if (Exch+ExchMn>0) ground=1; else ground=0;
      fprintf(output->file(0),"%e %e %i ", ETmaxshift, ExchMn, ground);
      if (bistable) for (int i=0;i<4;i++) {
        //if (i==lowest) continue;
        if (least_square[4+i]->N==0) fprintf(output->file(0),"%i/%i %i ",i,0,0); else
          fprintf(output->file(0),"%i/%i %e ",i,least_square[4+i]->occurences,least_square[4+i]->estimate);
      }
      fprintf(output->file(0),"%e %e %e %e %e %e %e %e ", RmaxT, RmaxS, RR1T, RR2T, RR1S, RR2S, RSmaxT, RSmaxS);
      fprintf(output->file(0),"%e %e %e %e %e %e %e\n", resT[0], resT[1], resT[2], resS[0], resS[1], resS[2], resS[3]);
      
      fflush(0);
      for (int i=0;i<4;i++) delete least_square[i];

    } //!second parameter (p2) loop ends here

  } //!first parameter (p1) loop ends here 
  
  delete state2e;
  delete s2un2e;
  delete un2s2e; 
  clean_up(states);
  sorting2e.deallocate();
  
  delete output;
  return(0);
}


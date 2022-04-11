#include "main.h"
#include "fftw3.h"

#define RELAXCHECK
#define RELAX1
//#define RELAX2

struct exptable_prototype {//<i exp(ikr) j> <k exp(-ikr) l> = <i exp(ikr) j> <l exp(ikr) k>^*
  bool sq;       //i==l and j==k
  bool sq_rev;   //i==k and j==l
  content *d1,*d2;//the two Fourier transforms (for i-j, and k-l)
  content *dataxy;//their product
  int nx,ny;      //dimension of the expanded table
  prec kxo,kyo;   //steps of the expanded table (in internal units of wavevector)

  //gives product of two fourier transforms, attenuated and expanded on a bigger grid
  //i-l are unsorted labels, expand is how many times to expand the grid (1=unchaged)
  void prepare(state_info& states, int i, int j, int k, int l, int expand);
};

struct inputxyz_prototype {
  //phonon wavevector input
  prec kx,ky,kz;	//components in cartesian coors
  prec phi,k,K;		//inplane radial coordinates and total magnitude
  //input for the relaxation:
  prec energy;		//total energy in internal units (equal to hbar K velocity)
  prec c1,c2;	    //auxiliary prefactors for |M(k)| for the deformation (1, S_u/S_d) or (S_d/S_u,1) or (0,0) depending which out of Sigma_u, Sigma_d is nonzero
  int Ix;		//what to integrate (previous sum.ix): I0 - Fourier over whole space, I1 - Fourier, I2 - previous divided by kz, I3 - previous times mf, I4 - previous times overz
  //auxiliary for the integration:
  int interaction;	//which interaction [0-deformation 1-piezo long 2-piezo trans 3-deformation transversal (Si)]
  prec k_max;	//upper limit for integral over k
  prec width_z_pi;      //prefactor of the argument of z-overlap, tu be multiplied by kz(in internal units)
  //parts of the underintegral function:
  content overxy;	//Fourier
  prec mf,overz;	//the direction dependent couplings, overlap along z axis
  //numerics parameters
  int polint_order;	//degree of the polynomial interpolation for the Fourier grid
  int qromx_order; 	//which order to use in qromb/qromo integrations over phi and k
  int qromx_maxsteps;	//how many evaluations in qromb (#eval = 2^maxsteps)
  prec qromx_err; //aimed at relative error

  prec factor_M();//geometrical couplings |M(K)|^2 of the electron-phonon interaction
  void prepare(exptable_prototype& t);//

};

prec overlap_z(prec p);//z-part of the overlap
prec occupation(prec energy, int plus);//thermal factor [1-e(-E/kB T)]^-1+1
prec prefactor(prec energy, int interaction);//the electron-phonon coupling strength - constant part
content value_polint(prec kx, prec ky, exptable_prototype& table, int order);//returns the value of Fourier at (kx, ky) making a polynomial interpolation of "order" order

inputxyz_prototype* inputxyz_global;
exptable_prototype* table_global;

//returns the integrated function evaluated at phi, k is already set in inputxyz
content at_phi(double phi)
{
  inputxyz_global->kx=inputxyz_global->k*cos(phi);
  inputxyz_global->ky=inputxyz_global->k*sin(phi);
  //inputxyz_global->phi=phi;
  inputxyz_global->overxy=value_polint(inputxyz_global->kx, inputxyz_global->ky, *table_global, inputxyz_global->polint_order);
  if (inputxyz_global->Ix<3) return(inputxyz_global->overxy); 
  inputxyz_global->mf=inputxyz_global->factor_M();
  return(inputxyz_global->overxy*inputxyz_global->mf);
}

//the intergal over phi is done
//if overkz is included, the input is kz, otherwise it is k
content at_k(double x)
{
#ifdef RELAXCHECK
  if (inputxyz_global->K<x && inputxyz_global->Ix!=0) {printf("argument of at_k (x=%e) is larger than the maximal vector K=%e\n",x, inputxyz_global->K);exit(1);}
#endif

  if (inputxyz_global->Ix<2) {//x is k 
    inputxyz_global->k=x;
    inputxyz_global->kz=sqrt(inputxyz_global->K*inputxyz_global->K-x*x);
  } else {//x is kz
    inputxyz_global->kz=x;
    inputxyz_global->k=sqrt(inputxyz_global->K*inputxyz_global->K-x*x);
  }

  int steps=inputxyz_global->qromx_maxsteps;
  prec err=inputxyz_global->qromx_err;

  content result=qromb(at_phi, 0, 2*M_PI, steps, inputxyz_global->qromx_order, err);

#ifdef RELAX2
  fprintf(logg,"at k=%e (x=%e) [1/nm] the integral over phi converged in %i steps (qromb, order=%i, allowed steps=%i) with value=(%.5e%+.5e), |value|=%e and error estimate=%e (vs goal %e [relative])\n",inputxyz_global->k*units.wavevector*units.nm, x*units.wavevector*units.nm, steps, inputxyz_global->qromx_order, inputxyz_global->qromx_maxsteps, result.real(), result.imag(), abs(result), err, inputxyz_global->qromx_err);
#endif

  if (inputxyz_global->Ix<2) result*=inputxyz_global->k;
  if (inputxyz_global->Ix>3) {
    inputxyz_global->overz=overlap_z(inputxyz_global->kz*inputxyz_global->width_z_pi);
    result*=inputxyz_global->overz;
  }
  return(result);
}

//integrate over the radial coordinate
content I_over_k()
{
  int steps=inputxyz_global->qromx_maxsteps;
  prec err=inputxyz_global->qromx_err;
  prec xmin=0, xmax=inputxyz_global->k_max;
  if (inputxyz_global->Ix>=2) {
    xmin=sqrt(pow(inputxyz_global->K,2)-pow(inputxyz_global->k_max,2));
    xmax=inputxyz_global->K;
  }
  content result=qromb(at_k, xmin, xmax, steps, inputxyz_global->qromx_order, err);
#ifdef RELAX2
  fprintf(logg,"the integral over k (%e to %e [1/nm]) converged in %i steps (qromb, order=%i, allowed steps=%i) with value=(%.5e%+.5e), |value|=%e and error estimate=%e (vs goal %e)\n", xmin*units.wavevector*units.nm, xmax*units.wavevector*units.nm, steps, inputxyz_global->qromx_order, inputxyz_global->qromx_maxsteps, result.real(), result.imag(), abs(result), err, inputxyz_global->qromx_err);
#endif
  return(result);
}

//the electron-phonon coupling strength - constant part
prec prefactor(prec energy, int interaction)
{
  prec pref;
  //into SI units
  energy*=units.energy;
  prec def=parameters.read("def","eV")*units.e;
  prec def2=parameters.read("def2","eV")*units.e;
  prec piezo=parameters.read("piezo","eV/m")*units.e;
  prec rho=parameters.read("density","kg/m3");
  prec cl=parameters.read("cl","m/s");
  prec ct=parameters.read("ct","m/s");

  switch (interaction) {
    case (0) : {//deformation longitudinal
      //it was scaled by the larger one in factor_M()
      pref=pow(energy,2)*pow(max(abs(def),abs(def2)),2)/(4*M_PI*M_PI*rho*pow(cl,4)*pow(units.hbar,3));
      break;
    }
    case (1) : {//piezoelectric longitudinal
      pref=pow(piezo,2)/(4*M_PI*M_PI*rho*cl*cl*units.hbar);
      break;
    }
    case (2) : {//piezoelectric transversal
      pref=pow(piezo,2)/(4*M_PI*M_PI*rho*ct*ct*units.hbar);
      break;
    }
    case (3) : {//deformation shear transversal
      pref=pow(energy,2)*pow(def2,2)/(4*M_PI*M_PI*rho*pow(ct,4)*pow(units.hbar,3));
      break;
    }
    default: {
      printf("wrong interaction (%i) in prefactor\n",interaction);
      pref=0;
      exit(1);
    }
  }
  return(pref);
}

//returns electron-phonon coupling (squared) depending on the phonon wavevector direction ("Mp"-factor)
//input: phonon k-vector (internal units) and the interaction (0-4)
prec inputxyz_prototype::factor_M()
{
  if (interaction==0) {//logitudinal deformation - only here the interaction strength is inserted, in others only in "prefactor"
    prec val=c1;
    if (K!=0) val+=c2*kz*kz/(K*K);
    return(val*val);
  }
#ifdef RELAXCHECK
  if (fabs(K*K-kz*kz-k*k)>1e-12*K*K || fabs(k*k-kx*kx-ky*ky)>1e-12*k*k) {
    printf ("not correct k/K in factor_M: difference, K:  %e %e, k: %e %e\n",K*K-kz*kz-ky*ky-kx*kx,K*K,k*k-kx*kx-ky*ky,k*k);
    exit(1);
  }
#endif

  if (k==0) return(0);//all futher are zero if the phonon is out of plane
  switch (interaction) {
    case (1): {//piezoelectric longitudinal
      prec val=6*kz*kx*ky/(K*K*K);
      return(val*val);
    }
    case (2): {//piezoelectric transversal
      prec val1=2*kz*(kx*kx-ky*ky)/(k*K*K);
      prec val2=2*kx*ky*(2*kz*kz-k*k)/(K*K*K*k);
      return(val1*val1+val2*val2);
    }
    case (3): {//transversal deformation (Si - shear)
      prec val=kz*k/(K*K);
      return(val*val);
    }
    default : {
      printf("wrong interaction tag (%i) in factor_M\n",interaction);
      exit(1);
    }
  }
}

//gives z-part of the overlap (dimensionless) for a infinite square well potential
//input - p = k_z * l_z / pi (dimensionless)
prec overlap_z(prec p)
{
  if (p<=0) return(1);
  prec p2=p*p;
  prec val=32*(1-cos(p*M_PI))/(M_PI*M_PI*p2*(p2-4)*(p2-4));
  return(val);
}

//thermal occupation factor; the input energy is absolute, the sign is in plus
//if plus is +1, add 1 to the exponential, if -1, nothing
prec occupation(prec energy, int plus)
{
  if (energy<0) {printf ("only non-negative energies allowed in occupation (energy:%e)\n",energy);exit(1);}
  prec temperature=parameters.read("T","K");
  prec n=0,value=0;
  if (temperature!=0) value=n=1/(exp(energy*units.energy/units.kB/temperature)-1);
  if (plus==1) value+=1;
  else if (plus!=-1) {printf("not allowed value of plus in occupation\n");exit(1);}
  fprintf(logg,"phonon occupation factor (at temperature %e Kelvin): n=%e, energy difference=%e meV, thermal factor=%e\n", temperature, n, energy*plus*units.energy/units.meV, value);
  return(value);
} 


//computes the fourier transform (NOT ATTENUATED -- THIS MUST BE DONE AFTER) 
//into the table with dimensions nx, ny corresponding to two states i, j (unsorted index)
//The spin part is contracted.
content* plane_wave_overlap(state_info& states, int nx, int ny, int i, int j)
{
  int Nx,Ny;
  states.r1->give_par_fft_2(Nx,Ny);
  content* fourier_in=states.Fourier.address_in(nx,ny);
  content* vec=new content[states.r1->Nw];
  content* table=new content[Nx*Ny];

  int dim=states.r1->Nw*states.r1->sets;
  states.r1->braket_vec(states.vysl->eigenvecs+dim*i,states.vysl->eigenvecs+dim*j,vec,reg1::set);
  states.r1->region2table(vec,0,table,Nx,Ny);
  delete[] vec;
  expandtable(table,Nx,Ny,fourier_in,nx,ny);
  delete[] table;
#ifdef RELAX2
  fprintf(logg,"computing Fourier transform of eigenstate product Psi_%i^+ Psi_%i (on a grid of %i by %i points expanded to %i by %i)\n", i, j, Nx, Ny,nx,ny);
#endif
  return(states.Fourier.execute());
}

void exptable_prototype::prepare(state_info& states, int i, int j, int k, int l, int expand)
{
  if (expand<1) {printf("too small expand factor (%i) in exptable_prototype::prepare\n",expand);exit(1);}
  states.r1->give_par_fft_2(nx,ny);
  nx*=expand;
  ny*=expand;

  states.r1->give_par_fft(kxo,kyo);
  kxo/=expand;
  kyo/=expand;
  dataxy=new content[nx*ny];
  sq=sq_rev=false;
  if (i==k && j==l) sq_rev=true;
  if (i==l && j==k) sq=true;

#ifdef RELAX1
  fprintf(logg,"filling the table: wavefunctions (%i, %i, %i, %i) (giving flags sq=%i, sq_rev=%i). Grid will be expanded to a %i by %i corresponding to Fourier wavevector grid steps %e, %e [1/nm]\n", i, j, k, l, sq, sq_rev, nx, ny, kxo*units.wavevector*units.nm, kyo*units.wavevector*units.nm);
#endif

  d1=d2=plane_wave_overlap(states, nx, ny, i,j);
  if (!sq && !sq_rev) {
    d1=new content[nx*ny];
    for (int i=0;i<nx*ny;i++) d1[i]=d2[i];
    d2=plane_wave_overlap(states,nx,ny,l,k);
  }

  int pntr_rev=nx*ny-1;
  int pntr=0;
  for (int i=0;i<nx;i++) {
    prec ax=states.Fourier.attenuation_factors[i];
    for (int j=0;j<ny;j++) {
      prec ay=states.Fourier.attenuation_factors[nx+j];
      if (!sq_rev) dataxy[pntr]=d1[pntr]*conj(d2[pntr])*ax*ax*ay*ay;
      else {
        dataxy[pntr]=d1[pntr]*d2[pntr_rev]*ax*ax*ay*ay;
        pntr_rev--;
      }
      pntr++;
    }
  }
  if (!sq && !sq_rev) delete d1;
}

//sets the limit for the integration (maximal k-vector)
void inputxyz_prototype::prepare(exptable_prototype& t)
{
  //control flags defaults
  qromx_order=4;
  qromx_maxsteps=14;
  qromx_err=1e-5;
  polint_order=5;
  Ix=(int) parameters.read("Ix","-");
  //a control flag overwritten
  /*switch ((int) parameters.read("unused3","-")) {
    case (1) : {qromx_order=(int) parameters.read("unused4","-"); break;}
    case (2) : {qromx_maxsteps=(int) parameters.read("unused4","-"); break;}
    case (3) : {qromx_err=parameters.read("unused4","-"); break;}
    case (4) : {polint_order=(int) parameters.read("unused4","-"); break;}
  }*/

  //interaction strengths
  c1=parameters.read("def","eV")*units.e;
  c2=parameters.read("def2","eV")*units.e;
  prec c3=max(abs(c1),abs(c2));
  if (c3!=0) {c1/=c3; c2/=c3;}

  //phonon wavevector
  K=fabs(energy)*units.energy/units.hbar;
  if (interaction<2) K/=parameters.read("cl","m/s"); else K/=parameters.read("ct","m/s");
  K/=units.wavevector;

  //maximal k in the integration
  k_max=min(t.nx*t.kxo/2, t.ny*t.kyo/2);
  if (Ix!=0) k_max=min(K,k_max);

#ifdef RELAX2
  fprintf(logg,"setting integration limits: the energy of %e meV corresponds for %i interaction to phonon wavevector %e [1/nm] versus the Fourier table steps of %e, %e and maxima %e, %e [1/nm]\n", energy*units.energy/units.meV, interaction, K*units.wavevector*units.nm, t.kxo*units.wavevector*units.nm, t.kyo*units.wavevector*units.nm, t.nx*t.kxo/2*units.wavevector*units.nm, t.ny*t.kyo/2*units.wavevector*units.nm);
#endif

}

//given real momentum returns value from table by 2D polynomial interpolation
//order means number of points to be used (at least 2)
content value_polint(prec kx, prec ky, exptable_prototype& table, int order)
{
  if (order<2) {printf("too small order (%i) in value_polint\n",order);exit(1);}
  static int previous_order=0;
  static content* datax=0;
  static content* datay=0;
  static prec* coor=0;
  static int* coori=0;

  if (previous_order!=order) {
    if (datax!=0) delete datax;
    if (datay!=0) delete datay;
    if (coor!=0) delete coor;
    if (coori!=0) delete coori;
    datax=new content[order];
    coor=new prec[order];
    coori=new int[order];
    datay=new content[order];
#ifdef RELAX2
    fprintf(logg, "value_polint: order changed from %i to %i, initializing\n",previous_order, order);
#endif
    for (int i=0;i<order;i++) {
      if ((i+1) % 2) coori[i]=-(i+1)/2; else coori[i]=(i+1)/2;
      coor[i]=coori[i];
#ifdef RELAX2
    fprintf(logg, "coordinate shift %i at position %i\n",coori[i],i);
#endif
    }
    previous_order=order;
  }

  int cx,cy;     //coordinates of the closest grid point
  prec dx,dy;     //normalized length from the closest grid point
  
  cx=(int) floor(kx/table.kxo);
  dx=kx/table.kxo-cx;
  cx=pos_mod(cx,table.nx);

  cy=(int) floor(ky/table.kyo);
  dy=ky/table.kyo-cy;
  cy=pos_mod(cy,table.ny);

#ifdef RELAX3
  fprintf(logg, "x-axis: for the input momentum %e and the table step %e the closest grid point to the left is at %e with coordinate %i and distance %e\n", kx, table.kxo, cx*table.kxo, cx, dx);
  fprintf(logg, "y-axis: for the input momentum %e and the table step %e the closest grid point to the left is at %e with coordinate %i and distance %e\n", ky, table.kyo, cy*table.kyo, cy, dy);
#endif

  prec err;
  for (int j=0;j<order;j++) {
    int lcy=pos_mod(cy+coori[j],table.ny);
    for (int i=0;i<order;i++) {
      int lcx=pos_mod(cx+coori[i],table.nx);
      datax[i]=table.dataxy[lcx*table.ny+lcy];
    }
    datay[j]=polint(coor,datax,order,dx,err);
  }
  content result=polint(coor,datay,order,dy,err);
#ifdef RELAX3
  fprintf(logg, "returning value (%e%+ei), as refinement of the closest point to the left with value (%e%+ei)\n", result.real(), result.imag(), table.dataxy[cx*table.ny+cy].real(),table.dataxy[cx*table.ny+cy].imag());
#endif
  return(result);
}


content spectral_density(state_info& states, int i, int j, int k, int l, int plus, int interaction, prec energy)
{
  energy*=units.meV/units.energy;

  if (energy<0) {
#ifdef RELAX1
    fprintf(logg,"spectral density for states (%i,%i,%i,%i), plus=%i, interaction=%i, omega=%e: 0 (energy<0)\n",i, j, k, l, plus, interaction, energy);
#endif
    return(content(0,0));
  }
  if (occupation(energy,plus)==0) {
#ifdef RELAX1
    fprintf(logg,"spectral density for states (%i,%i,%i,%i), plus=%i, interaction=%i, energy=%e: 0 (zero T)\n",i, j, k, l, plus, interaction, energy);
#endif
    return(content(0,0));
  }
  
  //compute the Fourier transform on an expanded grid
  exptable_prototype t; 
  table_global=&t;
  t.prepare(states, i, j, k, l, (int) parameters.read("CE_exp","-"));
 
  //input structure
  inputxyz_prototype in;
  inputxyz_global=&in;
  in.energy=energy;
  in.interaction=interaction;
  in.width_z_pi=parameters.read("width_z","nm")*units.nm/units.length/M_PI;
  in.prepare(t);
 
  //do the integral and normalize
  content res=I_over_k();
  res*=occupation(energy,plus);//thermal factor
  res*=prefactor(energy,interaction);//units
  //normalization
  res*=pow(units.wavevector,2);//for I0-I1 (k dk)
  if (in.Ix>1) res/=units.wavevector;//for 1/kz 
  res/=M_PI;//to spectral density

#ifdef RELAX1
  fprintf(logg,"spectral density for states (%i,%i,%i,%i), plus=%i, interaction=%i, underintegral=%i, energy=%e [meV]: %.6e%+.6ei\n",i, j, k, l, plus, interaction, (int) parameters.read("Ix","-"), energy*units.energy/units.meV, res.real(), res.imag());
#endif
  
  delete t.dataxy;
  return(res);
}

prec transition_rate(state_info& states,int from, int to, int interaction)
{
  prec energy=(states.vysl->eigenvals[from]-states.vysl->eigenvals[to]).real()*units.energy/units.meV;
  //prec energy=states.read(from,"energy")-states.read(to,"energy");
#ifdef RELAX1
  fprintf(logg,"transition_rate: energies: %e(%i) -> %e(%i)\n",states.vysl->eigenvals[from].real()*units.energy/units.meV, from, states.vysl->eigenvals[to].real()*units.energy/units.meV, to);
#endif
  int plus=1;
  if (energy<=0) {plus=-1;energy*=-1;}
  content res=spectral_density(states, from, to , to, from, plus, interaction, energy);
  return(M_PI*res.real());
}

//"unsorted" states must be already sorted according to the energy!!!!!!!!!!!!!!!!1
void transition_rates(state_info& states) 
{
  for (int from=states.unsorted-1;from>-1;from--) {
    for (int to=0;to<states.unsorted;to++) {
      prec rate[4];
      for (int i=0;i<3;i++) {
        prec energy=(states.vysl->eigenvals[from]-states.vysl->eigenvals[to]).real();
        if (from>to) {
          if (energy<0) {
            message("unsorted states are not sorted according to energy!\nin transition_rates\n");
            exit(1);
          }
          rate[i]=transition_rate(states,from,to,i);
        }
        else rate[i]=states.rates[4*states.unsorted*to+4*from+i]*occupation(-energy,-1)/(1+occupation(-energy,-1));
      }
      rate[3]=rate[0]+rate[1]+rate[2];
      for (int i=0;i<4;i++) states.rates[from*4*states.unsorted+to*4+i]=rate[i];
    }
  }
    
#ifdef RELAX1
    sprintf(buffer,"\n\ntable of transition rates:\n");
    splachni(logg,buffer,4);
    
    message(logg,"deformation potential\n",4);
    for (int from=0;from<states.unsorted;from++) {
      for (int to=0;to<states.unsorted;to++) {
        sprintf(buffer,"%e ",states.rates[from*4*states.unsorted+to*4+0]);
        splachni(logg,buffer,4);
      }
      message(logg,"\n",4);
    }
    
    message(logg,"\n\npiezoelectric (total)\n",4);
    for (int from=0;from<states.unsorted;from++) {
      for (int to=0;to<states.unsorted;to++) {
        sprintf(buffer,"%e ",states.rates[from*4*states.unsorted+to*4+1]+states.rates[from*4*states.unsorted+to*4+2]);
        splachni(logg,buffer,4);
      }
      message(logg,"\n",4);
    }
#endif
}

/*
void two_el_rate(state_info& states, int length, int* indexes, content* coefs, prec energy, prec* result)
//#define TWO_EL_RATE0
//#define TWO_EL_RATE1
{
  //check input parameters
  int eig=(int) parameters.read("eig","-");
  //fprintf(logg,"I AM HERE\n");
#ifdef TWO_EL_RATE0
  sprintf(buffer,"two_el_rate:\n\nlength of the index array: %i\n",length);
  splachni(logg,buffer,4);
#endif
  if (length<1) {fprintf(logg,"no input data in two_el_rate\n");exit(1);}
#ifdef TWO_EL_RATE1
  sprintf(buffer,"indexes:\n");splachni(logg,buffer,4);
#endif
  for (int i=0;i<length;i++) { 
    //int aux=states.sorted2un[indexes[i]];
    int aux=indexes[i];
#ifdef TWO_EL_RATE1
    sprintf(buffer,"#%2i:%2i\n",i,aux);splachni(logg,buffer,4);
#endif
    if (aux<1 || aux>eig) {fprintf(logg,"wrong index in two_el_rate:%i\n",aux);exit(1);}
  }
#ifdef TWO_EL_RATE0
  sprintf(buffer,"indexes look ok in two_el_rate\n");splachni(logg,buffer,4);
#endif
#ifdef TWO_EL_RATE1
  sprintf(buffer,"coefs:\n");splachni(logg,buffer,4);
#endif
  content trace=0;
  for (int i=0;i<length;i++) {
    trace+=coefs[i*length+i];
    for (int j=0;j<length;j++) {
      content aux=coefs[i*length+j];
#ifdef TWO_EL_RATE1
      sprintf(buffer,"%e%+e ",aux.real(),aux.imag());splachni(logg,buffer,4);
#endif
    }
  }
#ifdef TWO_EL_RATE0
  sprintf(buffer,"coefs look ok in two_el_rate\nthe matrix trace is (%e%+ei)\n",trace.real(), trace.imag());splachni(logg,buffer,4);
#endif  
#ifdef TWO_EL_RATE1
  sprintf(buffer,"energy:%e (being %e meV)\n",energy, energy*units.energy/units.meV);splachni(logg,buffer,4);
#endif

  //norm_of_difference(length, indexes, coefs);

  //construct the "wavefunctions"
  int wfl=states.r1->Nw*states.r1->sets;
  content cst=content(1,0);//content(1.0/sqrt(wfl),0);
  content* wf=new content[wfl*2];
  for (int n=0;n<wfl;n++) {
    wf[n]=content(0,0);
    wf[wfl+n]=cst;
    for (int i=0;i<length;i++) {
      int offseti=wfl*(indexes[i]-1);
      for (int j=0;j<length;j++) {
        int offsetj=wfl*(indexes[j]-1);
        wf[n]+=coefs[i*length+j]*conj(states.vysl->eigenvecs[offseti+n])*states.vysl->eigenvecs[offsetj+n];
      }
    }
  }

  //prepare the space for auxiliary "wavefunctions"
  //two wavefunctions will be needed
  content en0=states.vysl->eigenvals[0];
  content en1=states.vysl->eigenvals[1];
  content* wf_old=states.vysl->eigenvecs;
  
  states.vysl->eigenvecs=wf;
  states.vysl->eigenvals[0]=content(0,0);
  states.vysl->eigenvals[1]=content(energy,0);
  
  for (int i=0;i<3;i++) {
#ifdef TWO_EL_RATE0
  sprintf(buffer,"computing rate #%i\t",i);splachni(logg,buffer,4);
#endif  
    result[i]=transition_rate(states,1,0,i);
#ifdef TWO_EL_RATE0
  sprintf(buffer,"with result:%e\n",result[i]);splachni(logg,buffer,4);
#endif  
  }
 
  //delete the auxiliary wavefunctions
  delete wf;
  states.vysl->eigenvecs=wf_old;
  states.vysl->eigenvals[0]=en0;
  states.vysl->eigenvals[1]=en1;
  return;
}*/

/*struct sum_prototype { //to delete
  content i0,i1,i2,i3,i4;
  prec norm, pref;
};*/


/*struct stat_prototype {//to delete
  struct inputxyz_prototype max[2];
  struct inputxyz_prototype min[2];
  struct sum_prototype sum;
  int pntr;

  void init(); 
  void put(inputxyz_prototype& in);
  void normalize(inputxyz_prototype& in, int interaction, int plus);
  content puke(exptable_prototype& t, int log);
};*/


/*
void stat_prototype::init() {
  pntr=0;
  for (int i=0;i<2;i++) {
    max[i].kz=max[i].overz=max[i].mf=0;
    min[i].kz=min[i].overz=min[i].mf=1e3;
  }
  sum.i1=sum.i2=sum.i3=sum.i4=0;
}


void stat_prototype::normalize(inputxyz_prototype& in,int interaction, int plus)
{
  //printf("normalization2: %e\n",M_PI*kx_max*ky_max/pntr*t.norm));
  sum.norm=in.stepkx*in.stepky*pow(units.wavevector,2);//integral over kx and ky
  sum.norm/=units.wavevector;//for 1/kz
  sum.norm/=2*M_PI;//originally relaxation rate to spectral density now

  sum.pref=occupation(in.energy,plus);
  sum.pref*=prefactor(in.energy,in.interaction);
}

void stat_prototype::put(inputxyz_prototype& value){

  pntr++;
  if (min[0].kz>value.kz) min[0].kz=value.kz;
  if (max[0].kz<value.kz) max[0].kz=value.kz;
  if (min[0].mf>value.mf) min[0].mf=value.mf;
  if (max[0].mf<value.mf) max[0].mf=value.mf;
  if (min[0].overz>value.overz) min[0].overz=value.overz;
  if (max[0].overz<value.overz) max[0].overz=value.overz;
  
  if (abs(value.overxy)>max(abs(sum.i0)/(STEPS*STEPS)*1e-3,1e-12)) {
    if (min[1].kz>value.kz) min[1].kz=value.kz;
    if (max[1].kz<value.kz) max[1].kz=value.kz;
    if (min[1].mf>value.mf) min[1].mf=value.mf;
    if (max[1].mf<value.mf) max[1].mf=value.mf;
    if (min[1].overz>value.overz) min[1].overz=value.overz;
    if (max[1].overz<value.overz) max[1].overz=value.overz;
  }
  content aux=value.overxy;
  sum.i1+=aux;
  aux/=value.kz;
  sum.i2+=aux;
  aux*=value.mf;
  sum.i3+=aux;
  aux*=value.overz;
  sum.i4+=aux;

  //if (abs(aux)<0 || abs(aux)>1e5) {
    //sprintf(buffer,"too large value in stat_prototype(transition rate): %e\n",abs(aux));
    //splachni(logg,buffer,1+4);
  //}
#ifdef RELAXCHECK
  if (abs(value.overxy)<0) message(logg,"overxy is negative\n",1+4);
  if (value.kz<0) message(logg,"kz is negative\n",1+4);
  if (value.mf<0) message(logg,"mf is negative\n",1+4);
  if (value.overz<0) message(logg,"overz is negative\n",1+4);
#endif
}

content stat_prototype::puke(exptable_prototype& t, int where){

  sum.i0=0;
  for (int i=0;i<t.nx*t.ny;i++) sum.i0+=abs(t.dataxy[i]);
  
  sum.i0*=t.kxo*t.kyo*pow(units.wavevector,2)/(2.0*M_PI);
  sum.i1*=sum.norm;
  sum.i2*=sum.norm;
  sum.i3*=sum.norm;
  sum.i4*=sum.norm;

//write summary:
#ifdef RELAX2
  sprintf(buffer,"\nintegrals (with normalization %e):\ni0(overlap xy over all points):(%e%+ei) \ni1(overlap xy over the actual area):(%e%+ei) \ni2(previous with 1/kz):(%e%+ei) \ni3(previous with m_factor):(%e%+ei)\ni4(previous with overlap z):(%e%+ei)\n (i1-i4 computed in %i points)\n",sum.norm,sum.i0.real(),sum.i0.imag(),sum.i1.real(),sum.i1.imag(),sum.i2.real(),sum.i2.imag(),sum.i3.real(),sum.i3.imag(),sum.i4.real(),sum.i4.imag(),pntr);
  splachni(logg,buffer,where);
  
  sprintf(buffer,"\nspan of parameters:\n1/k_z: absolute (%.4e,%.4e), relevant (%.4e,%.4e)\n",1/max[0].kz, 1/min[0].kz,1/max[1].kz,1/min[1].kz);
  splachni(logg,buffer,where);
  
  sprintf(buffer,"m_factor: absolute (%.4e,%.4e), relevant (%.4e,%.4e)\n",min[0].mf,max[0].mf,min[1].mf,max[1].mf);
  splachni(logg,buffer,where);
  
  sprintf(buffer,"overlap z: absolute (%.4e,%.4e), relevant (%.4e,%.4e)\n",min[0].overz,max[0].overz,min[1].overz,max[1].overz);
  splachni(logg,buffer,where);
#endif
  //return(sum.i0); 
  //return(sum.i1*units.wavevector); 
  //return(sum.i2); 
  return(sum.i4*sum.pref); 
}

void sum_through_the_grid(inputxyz_prototype& in, stat_prototype& stat, exptable_prototype& t)
{
  prec width_z_pi=parameters.read("width_z","nm")*units.nm/units.length/M_PI;
  in.kx=-in.kx_max;
  while (in.kx<in.kx_max) {
    in.ky=-in.ky_max;
    while (in.ky<in.ky_max) {
      in.kz=delta(in);
      if (in.kz!=0) {
  content aux=in.overxy=four_spline(in.kx,in.ky,t);
  in.overxy=value_polint(in.kx, in.ky,t,4);
  //if (abs(in.overxy-aux)>1e-12) { printf("wrong polint of order 2: spline (%e%+ei) vs polint (%e%+ei)\n",in.overxy.real(), in.overxy.imag(), aux.real(), aux.imag());exit(1);}
  in.overz=overlap_z(in.kz*width_z_pi);
  in.mf=factor_M(in);
  stat.put(in);
      }
      in.ky+=in.stepky;
    }
    in.kx+=in.stepkx;
  }
}

//given real momentum returns value from table by bilinear spline
content four_spline(prec kx, prec ky, exptable_prototype& t)
{
  int lcx,lcy,rcx,rcy;    //coordinates of the left, right grid point
  prec dx,dy;     //normalized length from left grid point
  
  lcx=(int) floor(kx/t.kxo);
  dx=kx/t.kxo-lcx;
  rcx=pos_mod(lcx+1,t.nx);
  lcx=pos_mod(lcx,t.nx);
  
  lcy=(int) floor(ky/t.kyo);
  dy=ky/t.kyo-lcy;
  rcy=pos_mod(lcy+1,t.ny);
  lcy=pos_mod(lcy,t.ny);
  
  content f00=t.dataxy[lcx*t.ny+lcy];
  content f10=t.dataxy[rcx*t.ny+lcy];
  content f01=t.dataxy[lcx*t.ny+rcy];
  content f11=t.dataxy[rcx*t.ny+rcy];
  
  content result=(1-dx)*(1-dy)*f00+dx*(1-dy)*f10+dy*(1-dx)*f01+dx*dy*f11;
  //if (abs(result)<0) message(logg, "four spline: result negative\n",1+4);
  return(result);
}

*/


/*

//#include "gsl/gsl_errno.h"
//#include "gsl/gsl_spline.h"

//void fourier(ret_from_arpack_zn& vysl, reg1& r1,int from, int to, int mx, int my, exptable_prototype& output, int command);
void finer_grid(prec* t1, int nx, int ny, prec* t2,int mx, int my);
void finer_grid(content* t1, int nx, int ny, content* t2, int mx, int my);
prec norm(prec* t,int n);
prec norm2(content* t,int n);
void draw_table(content* t, int nx, int ny, FILE* file, int what);
void shift_table(content* tab,int tab_nx,int tab_ny,int shift_x, int shift_y);


//sum of absolute values
prec norm(prec* t,int n)
{
  prec val=0;
  for (int i=0;i<n;i++)
    val+=fabs(t[i]);
  return(val);
  
}

//sum of absolute values squared
prec norm2(content* t,int n)
{
  prec val=0;
  for (int i=0;i<n;i++)
    val+=(t[i]*conj(t[i])).real();
  return(val);
  
}

//write a complex table (nx times ny) into file
//what - 1 real part, 2 - imag part, 4 - absolute value
void draw_table(content* t, int nx, int ny, FILE* file, int what)
{
  prec val=0;
  for (int k=0;k<3;k++) {
    if (what%2==0) {
      what/=2;
      continue; 
    }  
    what/=2;
    
    for (int i=0;i<nx;i++) {
      for (int j=0;j<ny;j++) {
        switch(k) {
          case (0): {val=t[i*ny+j].real();break;}
          case (1): {val=t[i*ny+j].imag();break;}
          case (2): {val=sqrt((t[i*ny+j]*conj(t[i*ny+j])).real());break;}
        }
        sprintf(buffer,"%e ",val);
        splachni(file,buffer,2);
      }
      fprintf(file,"\n");
    }
  }
}

//interpolate a complex function (t1) defined on a grid points (nx times ny), result in t2, where between two old points inserts mx or my new points on the corresponding axes
void finer_grid(content* t1, int nx, int ny, content* t2, int mx, int my)
{
  //dimensions of the new table
  int tab_nx=nx+mx*(nx-1);
  int tab_ny=ny+my*(ny-1);

  //auxiliary for input
  prec *aux1_real=new prec[nx*ny];
  prec *aux1_imag=new prec[nx*ny];
   
  //auxiliary for output
  prec *aux2_real=new prec[tab_nx*tab_ny];
  prec *aux2_imag=new prec[tab_nx*tab_ny];
  
  for (int i=0;i<nx*ny;i++) {
    aux1_real[i]=t1[i].real();
    aux1_imag[i]=t1[i].imag();
  }
  
  finer_grid(aux1_real,nx,ny,aux2_real,mx,my);
  finer_grid(aux1_imag,nx,ny,aux2_imag,mx,my);
  

  //from auxiliary output into output
  for (int i=0;i<tab_nx*tab_ny;i++) {
    t2[i]=content(aux2_real[i],aux2_imag[i]);
  }

  delete aux1_real;
  delete aux2_real;
  delete aux1_imag;
  delete aux2_imag;

}

//interpolate the function on a grid into a finer grid
//in between each couple of points in the old grid inserts mx (my) new points
void finer_grid(prec* t1, int nx, int ny, prec* t2,int mx, int my)
{
  //type of interpolation (cspline, linear,...)
  const gsl_interp_type* type=gsl_interp_cspline;
  //dimensions of the new table
  int tab_nx=nx+mx*(nx-1);
  int tab_ny=ny+my*(ny-1);
  
  //first make interpolation on each x-line in table t1
  prec* auxx=new prec[ny];		//values of argument
  for (int k=0;k<ny;k++) auxx[k]=k;
  gsl_spline *spline=gsl_spline_alloc(type, ny);
  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  
  for (int i=0;i<nx;i++) {
    
     //compute spline
    gsl_spline_init(spline, auxx , t1+i*ny, ny);
      
     //save interpolated data in the new table
    int index_x=i*(mx+1);	//x coordinate in the new table
    for (int j=0;j<tab_ny;j++) {
      prec length=(double) j/ (double) (my+1);
      t2[index_x*tab_ny+j]=gsl_spline_eval (spline, length , acc);
    }
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  delete auxx;
   
   
   //now along each y-line in the table t2 
  spline=gsl_spline_alloc(type, nx);
  acc=gsl_interp_accel_alloc();
  auxx=new prec[nx];
  prec *auxy=new prec[nx];
  for (int k=0;k<nx;k++) auxx[k]=k;
  
  for (int j=0;j<tab_ny;j++) {
    
     //first create a vector to feed into spline
    for (int k=0;k<nx;k++) auxy[k]=t2[k*(mx+1)*tab_ny+j];
     
     //compute spline
    gsl_spline_init(spline, auxx , auxy, nx);
      
     //save interpolated data in the new table
    for (int i=0;i<tab_nx;i++) {
      prec length=(double) i/ (double) (mx+1);
      t2[i*tab_ny+j]=gsl_spline_eval (spline, length , acc);
    } 
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  delete auxx;
  delete auxy;
}  

//move the values in the table by a constant "distance" respecting periodicity
void shift_table(content* tab,int tab_nx,int tab_ny,int shift_x, int shift_y)
{

  content* taux=new content[tab_nx*tab_ny];
  for (int i=0;i<tab_nx*tab_ny;i++) taux[i]=tab[i];
  
  for (int i=0;i<tab_nx;i++) {
    int iloc=pos_mod(i+shift_x,tab_nx);
    for (int j=0;j<tab_ny;j++) {
      int jloc=pos_mod(j+shift_y,tab_ny);
      tab[iloc*tab_ny+jloc]=taux[i*tab_ny+j];
    }
  }
  delete taux;
}

void spin_and_charge_rate(state_info& states, int eig, int interaction, int log, prec* spin, prec* charge)
{
  ret_from_arpack_zn vysl(*(states.vysl));
  reg1 r1(*(states.r1));
//#define LOG_SCH
  //najprv najdi najnizsiu spin down hladinu
  int state1,state2=-1;
  prec spin11=0,spin22=0;
  for (state1=0;state1<eig;state1++) {
    spin11=(r1.mean_value(vysl.eigenvecs+r1.Nw*r1.sets*state1,reg1::sigmaz)).real();
    if (spin11<0.1) break;
    if (state1==eig-1) {
      sprintf(buffer,"no spin down state in spin_and_charge_rate!!!\nreturning zeros\n\n");
      splachni(logg,buffer,1+2);
      spin[0]=charge[0]=spin[1]=charge[1]=0;
      return;
    }
  }
  //a k nej najblizsiu spin up okrem ground state
  prec vzd=100*(vysl.eigenvals[0].real());
  prec vzd_loc;
  for (int j=1;j<eig;j++) {
    prec s=r1.mean_value(vysl.eigenvecs+r1.Nw*r1.sets*j,reg1::sigmaz).real();
    vzd_loc=fabs((vysl.eigenvals[j]-vysl.eigenvals[state1]).real());
    if (s>-0.1 && vzd_loc<vzd) {
      state2=j;
      spin22=s;
      vzd=vzd_loc;
    }
    if (j==eig-1 && state2==-1) {
      sprintf(buffer,"no spin up state in spin_and_charge_rate\nreturning zeros\n\n");
      splachni(logg,buffer,1+2);
      spin[0]=charge[0]=spin[1]=charge[1]=0;
      return;
    }
  }
  
  //zdiagonalizuj sigmaz
  content spin12=r1.braket(vysl.eigenvecs+r1.Nw*r1.sets*state1, 0, vysl.eigenvecs+r1.Nw*r1.sets*state2,0)-r1.braket(vysl.eigenvecs+r1.Nw*r1.sets*state1, 1, vysl.eigenvecs+r1.Nw*r1.sets*state2,1);
  prec aux1=sqrt((spin11-spin22)*(spin11-spin22)+4*(spin12*conj(spin12)).real());
  content aux2=conj(spin12)*content(2,0);
  content v1=(spin11-spin22-aux1)/aux2;   //spin down - v1*vec1+vec2
  content v2=(spin11-spin22+aux1)/aux2;   //spin up        - || -
  prec norm1=1/sqrt((v1*conj(v1)).real()+1);
  prec norm2=1/sqrt((v2*conj(v2)).real()+1);
  
    //create new eigenfunctions space
  content* aux=new content[r1.Nw*r1.sets*3];

  //ground state
  r1.operate(vysl.eigenvecs+r1.Nw*r1.sets*0,aux+r1.Nw*r1.sets*0,0,0, reg1::lin_com, content(1,0), 0 ,reg1::set);
  r1.operate(vysl.eigenvecs+r1.Nw*r1.sets*0,aux+r1.Nw*r1.sets*0,1,1, reg1::lin_com, content(1,0), 0 ,reg1::set);

  //'spin' (spin down)
  r1.operate(vysl.eigenvecs+r1.Nw*r1.sets*state2,aux+r1.Nw*r1.sets*1,0,0, reg1::lin_com, content(norm1,0), 0 ,reg1::set);
  r1.operate(vysl.eigenvecs+r1.Nw*r1.sets*state1,aux+r1.Nw*r1.sets*1,0,0, reg1::lin_com, v1*norm1, 0 ,reg1::add);
  
  r1.operate(vysl.eigenvecs+r1.Nw*r1.sets*state2,aux+r1.Nw*r1.sets*1,1,1, reg1::lin_com, content(norm1,0), 0 ,reg1::set);
  r1.operate(vysl.eigenvecs+r1.Nw*r1.sets*state1,aux+r1.Nw*r1.sets*1,1,1, reg1::lin_com, v1*norm1, 0 ,reg1::add);
  
  //'charge' (spin up)
  r1.operate(vysl.eigenvecs+r1.Nw*r1.sets*state2,aux+r1.Nw*r1.sets*2,0,0, reg1::lin_com, content(norm2,0), 0 ,reg1::set);
  r1.operate(vysl.eigenvecs+r1.Nw*r1.sets*state1,aux+r1.Nw*r1.sets*2,0,0, reg1::lin_com, v2*norm2, 0 ,reg1::add);
  
  r1.operate(vysl.eigenvecs+r1.Nw*r1.sets*state2,aux+r1.Nw*r1.sets*2,1,1, reg1::lin_com, content(norm2,0), 0 ,reg1::set);
  r1.operate(vysl.eigenvecs+r1.Nw*r1.sets*state1,aux+r1.Nw*r1.sets*2,1,1, reg1::lin_com, v2*norm2, 0 ,reg1::add);
  
  //check
  prec s1=r1.mean_value(aux+r1.Nw*r1.sets*1,reg1::sigmaz).real();
  prec s2=r1.mean_value(aux+r1.Nw*r1.sets*2,reg1::sigmaz).real();
#ifdef LOG_SCH
  sprintf(buffer,"spin down and up level: %i, %i\n",state1,state2); 
  splachni(logg,buffer,log);
  sprintf(buffer,"spins of these states: %e, %e, 'spin overlap':%e\n",spin11,spin22,norm(spin12));
  splachni(logg,buffer,log);
  sprintf(buffer,"spin up and down after diagonalization: %e, %e\n",s1,s2);
  splachni(logg,buffer,log);
#endif

  //prepare structure vysl
  content aux0[3];
  content* vecs_old=vysl.eigenvecs;
  content* vals_old=vysl.eigenvals;
  
  aux0[0]=vysl.eigenvals[0];
  aux0[1]=vysl.eigenvals[state1];
  aux0[2]=vysl.eigenvals[state2];

  vysl.eigenvecs=aux;
  vysl.eigenvals=aux0;
  
  //everything ready, compute transition rate
#ifdef LOG_SCH
  message(logg,"\ncoherent spin rate:---------------------\n",log);
#endif
  spin[0]=transition_rate(states, 1, 0, interaction);
#ifdef LOG_SCH
  message(logg,"coherent charge rate:-------------------\n",log);  
#endif
  charge[0]=transition_rate(states, 2, 0,interaction);
  
  vysl.eigenvecs=vecs_old;
  vysl.eigenvals=vals_old;
  delete aux;
  
  //non-coherent transition rate
  prec t1=transition_rate(states,state1,0,interaction);
  prec t2=transition_rate(states,state2,0,interaction);
  spin[1]=norm1*norm1*(t2+(v1*conj(v1)).real()*t1);
  charge[1]=norm2*norm2*(t2+(v2*conj(v2)).real()*t1);
  
#ifdef LOG_SCH
  message(logg,"\nnon coherent spin rate:-------------------\n",log);
  sprintf(buffer,"coefficients(squared): %e (state %i), %e (state %i)\nrate:%e\n",(v1*conj(v1)).real()*norm1*norm1,state1,norm1*norm1,state2,spin[1]); 
  splachni(logg,buffer,log);
  message(logg,"non coherent charge rate:-------------------\n",log);  
  sprintf(buffer,"coefficients(squared): %e (state %i), %e (state %i)\nrate:%e\n",(v2*conj(v2)).real()*norm2*norm2,state1,norm2*norm2,state2,charge[1]); 
  splachni(logg,buffer,log);
#endif

  return;
}

prec norm_of_difference(int length, int* indexes, content* coefs)
{
  prec sum=0, tot=0;
  bool compare=true;
  int length_old; 
  int* indexes_old=0;
  content* coefs_old=0;
  FILE* file=fopen("data/aux_file_C","r");
  if (file==0) {
    fprintf(logg,"no data saved\n");
    compare=false;
  }
  else {
    //read the data
    fscanf(file,"%d",&length_old);
    if (length_old<1 || length_old>100) {
      fprintf(logg,"saved length wrong: %i\n",length_old);
      return(0);
    }
    indexes_old=new int[length_old];
    coefs_old=new content[length_old*length_old];
    for (int j=0;j<length_old;j++) fscanf(file,"%d",indexes_old+j);
    for (int j=0;j<length_old;j++) for (int k=0;k<length_old;k++) {
      float im, re;
      fscanf(file,"%e %e", &re, &im);
      coefs_old[j*length_old+k]=content(re, im);
    }
    if (length!=length_old) {
      fprintf(logg,"number of states changed\n");
      compare=false;
    }
    fclose(file);
  }
  
  //save data
  file=fopen("data/aux_file_C","w");    
  fprintf(file,"%i\n",length);
  for (int j=0;j<length;j++) fprintf(file,"%i ",indexes[j]);
  fprintf(file,"\n\n");
  for (int j=0;j<length;j++) for (int k=0;k<length;k++) fprintf(file,"%e %e\n", coefs[j*length+k].real(), coefs[j*length+k].imag());
  fclose(file);

  
  if (compare) { 
    //compute the two matrices 
    content* c=new content[100*100];
    content* c_old=new content[100*100];
    for (int i=0;i<10000;i++) c[i]=c_old[i]=0;
    for (int j=0;j<length_old;j++) for (int k=0;k<length_old;k++) c_old[indexes_old[j]*length_old+indexes_old[k]]+=coefs_old[j*length_old+k];
    for (int j=0;j<length;j++) for (int k=0;k<length;k++) c[indexes[j]*length+indexes[k]]+=coefs[j*length+k];
    //and the difference
    for (int i=0;i<10000;i++) {sum+=abs(c[i]-c_old[i]);tot+=abs(c[i]);}
    delete c;
    delete c_old;
    fprintf(logg, "at length %i, the matrix change from the previous run is: %e (relative:%e)\n",length,sum, sum/tot);
  }
  if (indexes_old!=0) delete indexes_old;
  if (coefs_old!=0) delete coefs_old;
  
  return(sum);
}*/


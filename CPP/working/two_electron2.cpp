#include "main.h"

//!change group eigenvalues into real instead of complex

#define TWOE 2
#define TWOE_CHECK

prec is_hermitian(int N, content* H, bool print)
{
  prec res=0, norm=0;
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) {
      if (print) fprintf(logg,"%.5e%+.5e ",H[i*N+j].real(),H[i*N+j].imag());
      if (i<=j) res+=abs(H[i*N+j]-conj(H[j*N+i]));
      norm+=abs(H[i*N+j]);
    }
    if (print) fprintf(logg,"\n");
  }
  if (norm>0) return(res/norm); else return(0);
}

prec is_orthonormal(int N, int Neig, content* vecL, content* vecR, bool print)
{
  if (vecR==0) vecR=vecL;
  prec res=0;
  content norm=0;
  for (int i=0;i<Neig;i++) {
    for (int j=0;j<Neig;j++) {
      content braket=0;
      for (int k=0;k<N;k++) braket+=conj(vecL[i*N+k])*vecR[j*N+k];
      if (i==j) norm+=braket;
      if (print) fprintf(logg,"scalar product of vectors |< #%i | #%i >| = %.5e\n",i,j,abs(braket));
      if (i==j) braket-=1.0;
      res+=abs(braket);
    }
  }
  if (abs(norm)>0) return(abs(res/norm)); else return(0);
}


prec is_eigensys(int N, int Neig, content* H, content* val, content* vecL, content* vecR, bool print)
{
  if (vecR==0) vecR=vecL;
  prec res=0;
  content norm=0;
  for (int i=0;i<Neig;i++) {
    for (int j=0;j<Neig;j++) {
      content braHmEket=0;
      for (int k=0;k<N;k++) {
        if (i==j) norm+=conj(vecL[i*N+k])*vecR[j*N+k];
        content aux=0;
        for (int l=0;l<N;l++) aux+=H[k*N+l]*vecR[j*N+l];
        braHmEket+=conj(vecL[i*N+k])*(aux-val[j]*vecR[j*N+k]);
      }
      if (print) fprintf(logg,"scalar product |< #%i | H - E_%i  #%i >| = %.5e\n",i,j,j,abs(braHmEket));
      res+=abs(braHmEket);
    }
  }
  if (abs(norm)>0) return(res/abs(norm)); else return(0);
}


/*prec is_eigensys2(int N, int Neig, content* H, content* val, content* vecL, content* vecR, bool print) 
{
  //compute VL^dagger . (H-eigval 1) . VR
  content* HmE=new content[2*N*Neig];
  content* res=HmE+N*Neig;
  for (int i=0;i<Neig;i++) {
    for (int k=0;k<N;k++) HmE[i*N+k]=H[i];
    HmE[i*N+i]-=val[i];
  }
  //first HmE times VR^T both 
  char TRANSA='T', TRANSB='N';
  ZGEMM(&TRANSA, &TRANSB, &Neig, &N, &Neig, &alpha, H,  
  */

bool symmetry_prototype::symmetry1e(int u, int& id)
{
  int Ii, Ixi, Iyi;
  prec I=states->read(u,"I");
  if (abs(1.0-abs(I))<tolerance) Ii=(int) my_round(I); else Ii=0;
  prec Ix=states->read(u,"Ix");
  if (abs(1.0-abs(Ix))<tolerance) Ixi=(int) my_round(Ix); else Ixi=0;
  prec Iy=states->read(u,"Iy");
  if (abs(1.0-abs(Iy))<tolerance) Iyi=(int) my_round(Iy); else Iyi=0;
  prec lz=states->read(u,"lz");

  switch (type) {
 
    case 0 : {id=1;return(true);}
    case 1 : {
      if (Ii==0) return(false);
      if (Ii==1) id=1; else id=2;
      return(true);
    }
    case 2 : {
      if (Ixi*Iyi==0) return(false);
      IxIy2id(Ixi, Iyi, id);
      return(true);
    }
    case 3 : {
      if (abs(my_round(lz)-lz)<tolerance) {
      id=(int) my_round(lz);
        //fprintf(logg,"from lz=%e, symmetry %i identified\n",lz,id);
	return(true);
      }
      else return(false);
    }
  }
  return(false);
}

void symmetry_prototype::id2IxIy(int id, int& Ix, int& Iy)
{
  switch (id) {
    case 1 : {Ix=1; Iy=1; break;}
    case 2 : {Ix=-1; Iy=1; break;}
    case 3 : {Ix=-1; Iy=-1; break;}
    case 4 : {Ix=1; Iy=-1; break;}
    default :  {
      sprintf(buffer,"wrong symmetry identificator in id2IxIy: %i\n",id);
      splachni(logg,buffer,1+4);
      exit(1);
    }
  }
}

void symmetry_prototype::IxIy2id(int Ix, int Iy, int& id)
{
  if (Ix==1 && Iy==1) {id=1;return;}
  if (Ix==-1 && Iy==1) {id=2;return;}
  if (Ix==-1 && Iy==-1) {id=3;return;}
  if (Ix==1 && Iy==-1) {id=4;return;}
  sprintf(buffer,"wrong symmetries IxIy2id: Ix=%i, Iy=%i\n",Ix, Iy);
  splachni(logg,buffer,1+4);
  exit(1);
}

int symmetry_prototype::symmetry2e(int u1, int u2, int p)
{
  int id1, id2;
  bool sharp1=symmetry1e(u1, id1);
  bool sharp2=symmetry1e(u2, id2);
  if (!sharp1 || !sharp2) return(p);
  switch (type) {
    case 0 : return(p);
    case 1 : {
      if (id1==id2) return(p);
      return(2*p);
    }
    case 2 : {
      int Ix1, Ix2, Iy1, Iy2, id;
      id2IxIy(id1,Ix1,Iy1);
      id2IxIy(id2,Ix2,Iy2);
      IxIy2id(Ix1*Ix2,Iy1*Iy2,id);
      return(id*p);
    }
    case 3 : return(2*(id1+id2)+(p==1));
  }
  return(0);
}

int symmetry_prototype::symmetry2einv(int id)
{

  //printf("symm inv called with id=%i\n",id);
  switch (type) {
    case 0 : ;
    case 1 : ;
    case 2 : return((int) max(id,-id));
    case 3 : return(id/2);
  }
  return(0);
}


void symmetry_prototype::initialize(state_info* states_, int type_)
{
  states=states_;
  type=type_;
  if (type<0 || type>3) {
    sprintf(buffer,"wrong symmetry type in symmetry initialization: %i\n",type_);
    splachni(logg,buffer,1+4);
    exit(1);
  }
  tolerance=parameters.read("threshold","-");
  if (tolerance<0){
    sprintf(buffer,"warning: tolerance negative in symmetry initialization: %e\n",tolerance);
    splachni(logg, buffer,4);
  }
}


void group_prototype::initialize(int N)
{
  if (N<1) {
    sprintf(buffer, "the length of arrays not positive (N=%i) in group_prototype::initialize\n",N);
    splachni(logg,buffer,1+4);
    exit(1);
  }
  this->N=N;
  k=new int[N];
  H=new content[N*N];
  eigenvecs=new content[N*N];
  eigenvals=new content[N];
}

void group_prototype::release()
{
  if (k!=0) delete k; else {sprintf(buffer, "k unallocated in group.release\n");splachni(logg,buffer,1+4);exit(1);}
  if (H!=0) delete H; else {sprintf(buffer, "H unallocated in group.release\n");splachni(logg,buffer,1+4);exit(1);}
  if (eigenvecs!=0) delete eigenvecs; else {sprintf(buffer, "eigenvecs unallocated in group.release\n");splachni(logg,buffer,1+4);exit(1);}
  if (eigenvals!=0) delete eigenvals; else {sprintf(buffer, "eigenvals unallocated in group.release\n");splachni(logg,buffer,1+4);exit(1);}
}


two_electron_state_prototype::two_electron_state_prototype(state_info* states_, sorting_class * sorting_)
{
  selected=0;
  statek=0;
  group=0;
  CE=0;
  statea=0;
  stateb=0;
  Hb=0;
  statef=0;
  statef_vecs_begin=0;
  rho=0;
  mat_ele.derivs=mat_ele.orbital=mat_ele.orbital2e=mat_ele.spinor=mat_ele.spin_imp=mat_ele.spin_imp2e=0;
  
  sorting=sorting_;
  states=states_;
  Ns=(int) parameters.read("selected2e","-");
  if (Ns<1) {
    sprintf(buffer,"number of selected 1e states:%i\n",Ns);
    splachni(logg,buffer,1+4);
    exit(1);
  }
  Nf=(int) parameters.read("eig2e","-");
  if (Nf<1) {
    fprintf(logg, "number of 2e eigenstates to get too small(%i) -> reset to one\n",Nf);
    Nf=1;
  }
  progress_bar = new progress_bar_prototype("generic");
  timing.c_elems=0;
  timing.final_diag=0;
}

two_electron_state_prototype::~two_electron_state_prototype() 
{
  if (selected!=0) delete selected;
  if (statek!=0) delete statek;
  if (group!=0) {
    for (int g=0;g<Ng;g++) group[g].release();
    delete group;
  }
  if (statea!=0) {
    if (statea[0].aij!=0) delete statea[0].aij;
    delete statea;
  }
  
  if (stateb!=0) delete stateb;
  if (Hb!=0) delete Hb;
  if (statef!=0) {
    if (statef[0].cfkS!=0) delete statef[0].cfkS;
    if (statef[0].prob_g!=0) delete statef[0].prob_g;
    delete statef;
  }
  if (statef_vecs_begin!=0) delete statef_vecs_begin;
  delete progress_bar;
  if (mat_ele.derivs!=0) delete mat_ele.derivs;
  if (mat_ele.orbital!=0) delete mat_ele.orbital;  
  if (mat_ele.orbital2e!=0) delete mat_ele.orbital2e;
  if (mat_ele.spinor!=0) delete mat_ele.spinor;
  if (mat_ele.spin_imp!=0) delete mat_ele.spin_imp;  
  if (mat_ele.spin_imp2e!=0) delete mat_ele.spin_imp2e;
  
  if (rho!=0) delete rho;
}

void two_electron_state_prototype::diagonalize(prec par, int* sorted2un, int *unsorted2s, int step)
{
  if (step<0 || step>7) {printf("wrong step in state2e.diagonalize: %i\n",step); exit(1);}
  switch (step) {
    case 0 : 
      //printf("selected\n");
      construct_selected();
      fflush(logg);
    case 1 : 
      if (selected==0) {printf("selected states do not exist in state2e.diagonalize step-1\n"); exit(1);}
      //printf("statek\n");
      construct_statek();
      fflush(logg);
    case 2 : 
      if (statek==0) {printf("k-states do not exist in state2e.diagonalize step-2\n"); exit(1);}
      //printf("groups const\n");
      construct_groups();
      fflush(logg);
      diagonalize_groups();
      fflush(logg);
      //printf("groups diag\n");
    case 3 : 
      if (group==0) {printf("groups do not exist in state2e.diagonalize step-3\n"); exit(1);}
      //printf("statea\n");
      construct_statea(par, sorted2un, unsorted2s);
      fflush(logg);
    case 4 : 
      if (statea==0) {printf("a-states do not exist in state2e.diagonalize step-4\n"); exit(1);}
      //printf("stateb\n");
      construct_stateb();
      fflush(logg);
    case 5 : 
      if (selected==0) {printf("selected states do not exist in state2e.diagonalize step-5\n"); exit(1);}
      //printf("single-e matrix elements-1\n");
      mat_ele_fill_derivs();	//derivatives
      fflush(logg);
      mat_ele_fill_orbital();	//orbital part
      mat_ele_fill_spinor();	//spinor part
      fflush(logg);
      if (mat_ele.orbital2e!=0) delete mat_ele.orbital2e;
      mat_ele.orbital2e=new content[Na*Na*3];
      mat_ele_convert(mat_ele.orbital, mat_ele.orbital2e, "SOC m.e. 2e");
      fflush(logg);
    case 6 :  {
      if (selected==0) {printf("selected states do not exist in state2e.diagonalize step-6\n"); exit(1);}
      //printf("single-e matrix elements-2\n");
      prec Arho=-parameters.read("exch","meV")*units.meV/units.energy*parameters.read("xMn", "-")*(5.0/2)*(1.0/2);
      if (spin_impurities.I[0]==0 || spin_impurities.N==0) Arho=0;
      if (Arho!=0) {
	mat_ele_fill_spin_imp();//spin impurity part
        if (mat_ele.spin_imp2e!=0) delete mat_ele.spin_imp2e;
        mat_ele.spin_imp2e=new content[Na*Na*3];
        mat_ele_convert(mat_ele.spin_imp, mat_ele.spin_imp2e, "nuclear m.e. 2e");
        fflush(logg);
      }
    }
    case 7 :
      if (stateb==0) {printf("b-states do not exist in state2e.diagonalize step-5\n"); exit(1);}
      //printf("statef\n");
      construct_statef();
      fflush(logg);
      //printf("finished\n");
  }
}


void two_electron_state_prototype::construct_selected()
{
  int type=3;//try all symmetries staring with the most special
  prec tolerance=parameters.read("threshold","-");
  //construct selected states
  if (selected!=0) delete selected;
  selected=new selected_prototype[Ns];
  int assigned;
restart:
  int id2=0, idsum=0;
  symmetry.initialize(states,type);
  assigned=0;
  for (int s=0;s<states->sorted;s++) {
    int u=states->sorted2un[s];
    //the state must be visible
    if (u==-1) {
      sprintf(buffer,"warning: sorted state #%i no more visible\n",s);
      splachni(logg, buffer,4);
      continue;
    }
    //the state must be spin up
    if (states->read(u,"s")<1-tolerance) continue;

    int id;
    if (!symmetry.symmetry1e(u,id)) {
      //if the state does not have a sharp symmetry, reset the symmetry type and restart
      sprintf(buffer,"WARNING: state (sorted:%i, unsorted:%i) without sharp symmetry of type %i, restarting with lower orbital symmetry...\n",s,u,type);
      splachni(logg, buffer,4);
      if (type==0) {
        sprintf(buffer,"symmetry cannot be lowered anymore in construct_selected\n");
        splachni(logg,buffer,1+4);
        exit(1);
      }
      type--;
      goto restart;
    }
    id2+=id*id;
    idsum+=id;

    selected[assigned].u=u;
    selected[assigned].valley=1;
    assigned++;
    if (assigned==Ns) break;
  }
  if (assigned<Ns) {
    sprintf(buffer, "not enough states to fill all %i selected (only %i suitable found\n",Ns,assigned);
    splachni(logg,buffer,1+4);
    exit(1);
  }
  if (id2==idsum*idsum && type==3) {//all states the same symmetry with lz
    type--;
    goto restart;
  }

#if TWOE > 0
  fprintf(logg,"at symmetry type %i, selected were the following %i states:\n",type, Ns);
  for (int i=0;i<Ns;i++) {
    int id;
    symmetry.symmetry1e(selected[i].u,id);
    fprintf(logg, " selected #%2i: unsorted=%2i, valley=%i, id=%+2i\n", i, selected[i].u, selected[i].valley, id);
  }
  fprintf(logg,"\n");
#endif
}

void two_electron_state_prototype::construct_statek()
{
  prec muB=0;
  int sets=(int) parameters.read("sets","-");
  if (sets==1) sprintf(buffer,"no");
  else if (sets==2) {muB=zeeman();sprintf(buffer,"%e [meV] of",muB*units.energy/units.meV);}
  else {printf("sets=%i not implemented in two_electron_state_prototype\n",sets);exit(1);}
#if TWOE > 0
  fprintf(logg,"%s Zeeman energy will be subtracted from 1e spectrum (muB=%e vs Zeeman=%e)\n\n",buffer,muB,zeeman());
#endif
  Nk=Ns*Ns;
  if (statek!=0) delete statek;
  statek=new statek_prototype[Nk];
  int k=0;
  for (int i=0;i<Ns;i++) {
    for (int j=i;j<Ns;j++) {
      for (int p=1;p>-2;p-=2) {
        if (i==j && p==-1) continue;
        statek[k].i=i;
        statek[k].j=j;
        statek[k].p=p;
        statek[k].id=symmetry.symmetry2e(selected[i].u,selected[j].u,p);
        if (i==j) statek[k].norm=1.0/2; else statek[k].norm=1.0/sqrt(2);
        statek[k].energy=(states->read(selected[i].u,"energy")+states->read(selected[j].u,"energy"))*units.meV/units.energy-muB;
#if TWOE > 1
        fprintf(logg," k-state #%2i added: i=%2i, j=%2i, p=%+i, id=%+i, norm=%f energy=%f [meV]\n",k,statek[k].i, statek[k].j, statek[k].p, statek[k].id, statek[k].norm, statek[k].energy*units.energy/units.meV);
#endif
        k++;
      }
    }
  }
#if TWOE > 0
  fprintf(logg,"k-state array filled by %i states\n\n",k);
#endif
  if (k!=Ns*Ns) exit(1); 
}

int two_electron_state_prototype::less_than_E_in_group(prec energy, int g, int where)
{
  int sum=0;
  switch (where) {
    case 0 : { 
      for (int i=0;i<group[g].N;i++) 
      if (statek[group[g].k[i]].energy<=energy) sum++;
      break;
    }
    case 1 : {
      for (int i=0;i<group[g].Neig;i++) 
      if (abs(group[g].eigenvals[i])<=energy) sum++;
      break;
    }
    default : {
      sprintf(buffer, "wrong choice in less_than_E_in_group: where=%i\n",where);
      splachni(logg,buffer,1+4);
      exit(1);
    }
  }
  return(sum);
}

void two_electron_state_prototype::construct_groups(int Neigtot)
{

  //find minimal and maximal symmetry identificators and energies for k-states
  prec energy_min, energy_max;
  energy_min=energy_max=statek[0].energy;
  symmetry.min=symmetry.max=statek[0].id;
  for (int k=1;k<Nk;k++) {
    if (statek[k].id>symmetry.max) symmetry.max=statek[k].id;
    if (statek[k].id<symmetry.min) symmetry.min=statek[k].id;
    if (statek[k].energy<energy_min) energy_min=statek[k].energy;
    if (statek[k].energy>energy_max) energy_max=statek[k].energy;
  }
#if TWOE > 0
  fprintf(logg,"k-state symmetry (2e) identificators range: (%i,%i) and energy range: (%.5e,%.5e)\n",symmetry.min,symmetry.max, energy_min, energy_max);
#endif


  //allocate groups - more then needed, but it is ok
  if (group!=0) {
    for (int g=0;g<Ng;g++) group[g].release();
    delete group;
  }
  group=new group_prototype[symmetry.max-symmetry.min+1];
  Ng=0;

  //find out number of members in each symmetry group - add the group to the pool only if it has nonzero number of memebrs
  int id=symmetry.min;
  do {
    int sum=0;
    for (int k=0;k<Nk;k++) if (statek[k].id==id) sum++;
    if (sum>0) {
      group[Ng].initialize(sum);
      group[Ng].id=id;
#if TWOE > 0
      fprintf(logg,"group #%i added: symmetry identificator:%+i, number of items:%2i\n",Ng,id,sum);
#endif
      int kg=0;
      for (int k=0;k<Nk;k++) if (statek[k].id==id) {group[Ng].k[kg++]=k;statek[k].g=Ng;}
#if TWOE > 1
      fprintf(logg," the k-states in the group:\n");
      for (int kg=0;kg<group[Ng].N;kg++) {
        fprintf(logg,"  state#%i=kg: k-label=%2i, id=%+i\n",kg,group[Ng].k[kg], statek[group[Ng].k[kg]].id);
      }
      fprintf(logg,"\n");
#endif
      Ng++;
    }
    id++;
  } while (id<=symmetry.max);

#if TWOE > 0
    {int sum=0;
    for (int g=0;g<Ng;g++) sum+=group[g].N;
    fprintf(logg,"total number of states in all groups is %i\n\n", sum);
    if (sum!=Nk) exit(1);}
#endif

  //set the energy cutoff to ask for Neigtot states
  //sanity check: Neigtot should be from 1 to Nk
  if (Neigtot<1 || Neigtot>Nk) {fprintf(logg,"Neigtot was reset from %i to Nk in construct_groups\n",Neigtot);Neigtot=Nk;}

  int step=0;
  prec energy;
  do {
    int sum=0;
    energy=(energy_min+energy_max)/2;
    for (int g=0;g<Ng;g++) sum+=less_than_E_in_group(energy,g,0);
    if (sum<Neigtot) energy_min=energy;
    if (sum>Neigtot) energy_max=energy;
    if (sum==Neigtot) break;
    if (step++==20) {
      energy=energy_max;
      break;
    }
  } while (true);
  int sum=0;
  for (int g=0;g<Ng;g++) {
    sum+=group[g].Neig=less_than_E_in_group(energy,g,0);
#if TWOE > 0
    fprintf(logg," in group #%i asking for %i eigenvalues\n",g, group[g].Neig);
#endif
  }
#if TWOE > 0
  fprintf(logg,"totally %i eigenstates are being asked for, vs %i planned\n\n",sum, Neigtot);
#endif

}

int diagonalize_arp_expand_by;
void diagonalize_arp_aux(int N, void *A, content* x, content* y)
{
  content* H= (content *) A;
  N-=diagonalize_arp_expand_by;
  for (int i=0;i<N;i++) {
    y[i]=0;
    for (int j=0;j<N;j++) y[i]+=H[i*N+j]*x[j];
  }
  for (int i=0;i<diagonalize_arp_expand_by;i++) y[N+i]=x[N+i]*(1.0+i)*1e+6;
} 


void diagonalize_arp(int N, int Neig, content *H, content *vecs, content * vals, const char* arpack_type, progress_bar_prototype* progress_bar)
{
#if TWOE > 1
  fprintf(logg,"diagonalization called, matrix dimension:%i, eigenvectors asked for:%i\n", N, Neig);
#endif
  diagonalize_arp_expand_by=max(Neig-(N-2),0); 
  int N_arp=N+diagonalize_arp_expand_by;
#if TWOE > 1
  fprintf(logg,"expanding by:%i to dimension:%i\n", diagonalize_arp_expand_by, N_arp);
  fprintf(logg,"the matrix is hermitian with relative precision %.5e\n",is_hermitian(N,H));
#endif

  arpack_zn arp(N_arp,Neig,-1e-12,1000,true, arpack_type);
  arp.set_progress_bar(progress_bar);
  ret_from_arpack_zn vysl;
  vysl.eigenvals=new content[Neig+1];
  vysl.eigenvecs=new content[N_arp*Neig];
  prec thr=-1e-12;
  do {
    arp.set_parameter(3,&thr);
    if (arp.go_for_it(N_arp,diagonalize_arp_aux,(void*) H,&vysl,true)) break;
    if (thr<0) thr=-thr; else thr*=1000;
    printf("arpack failed to converge, resetting precision to %e\n",thr);
  } while (true);
#if TWOE > 1
  arp.show_stats(&vysl,logg,4,true,true);
#endif
 
  //sort in ascending order
  int* pos=new int[Neig];
  for (int i=0;i<Neig;i++) pos[i]=i;

  char what=arpack_type[1];
  for (int i=0;i<Neig;i++) {
    for (int j=0;j<Neig-i-1;j++) {
      bool swap=false;
      switch (what) {
        case 'M' : { if (abs(vysl.eigenvals[j])>abs(vysl.eigenvals[j+1])) swap=true;break;}
        case 'R' : { if (abs(vysl.eigenvals[j].real())>abs(vysl.eigenvals[j+1].real())) swap=true;break;}
        case 'I' : { if (abs(vysl.eigenvals[j].imag())>abs(vysl.eigenvals[j+1].imag())) swap=true;break;}
        default: {fprintf(logg,"wrong sorting type (%c) in diagonalize_arp\n",what);exit(0);}
      }
      if (swap) {
        content aux=vysl.eigenvals[j]; vysl.eigenvals[j]=vysl.eigenvals[j+1]; vysl.eigenvals[j+1]=aux;
        int aux2=pos[j]; pos[j]=pos[j+1]; pos[j+1]=aux2;
      }
    }
  }

  for (int i=0;i<Neig;i++) {
    if (vecs!=0) for (int j=0;j<N;j++) vecs[i*N+j]=vysl.eigenvecs[pos[i]*N_arp+j];
    vals[i]=vysl.eigenvals[i];
  }

  delete[] pos;
  delete vysl.eigenvecs;
  delete vysl.eigenvals;
}

void two_electron_state_prototype::diagonalize_groups()
{
  //compute the Coulomb elements
  
  //find the total number of them
  int Ngmax=0, CE_tot=0;
  for (int g=0;g<Ng;g++) {
    CE_tot+=(group[g].N*(group[g].N+1))/2;
    if (group[g].N>Ngmax) Ngmax=group[g].N;
  }

  //prepare the coefficients for the integral
  int CE_exp=(int) parameters.read("CE_exp","-");
  progress_bar->reset("Coulomb m.e. - grid");
  progress_bar->start();
  //if (CE!=0) delete CE;
  CE_prototype* CE=new CE_prototype(states, Ns, CE_exp, 8000);//CE class
  //CoulombElement(*states, 0, 0, 0, 0, CE_exp, 1);
  progress_bar->finished(true);

  progress_bar->reset("Coulomb m.e.");
  progress_bar->add(CE_tot,0);
  progress_bar->start();

  //compute the Coulomb matrix elements
  for (int g=0;g<Ng;g++) {//for each group
    for (int kg1=0;kg1<group[g].N;kg1++) {//can really start now
      for (int kg2=kg1;kg2<group[g].N;kg2++) {
        statek_prototype* k1=statek+group[g].k[kg1];
        statek_prototype* k2=statek+group[g].k[kg2];
        int ui=selected[k1->i].u;
        int uj=selected[k1->j].u;
        int uk=selected[k2->i].u;
        int ul=selected[k2->j].u;
        int p1=k1->p;
        int p2=k2->p;
        if (p1!=p2) {sprintf(buffer, "exchange symmetry of the two k states differ, but Coulomb element asked for\n");splachni(logg,buffer,1+4);exit(1);}
        
        //!start-new method of CE calculation
        //content calt=CE->CoulombElement(ui,uj,uk,ul);
        //if (ui==uj || uk==ul) calt*=2; else calt+=((double) p1)*CE->CoulombElement(ui,uj,ul,uk);
	content c=CE->CoulombElementTweak(ui,uj,uk,ul);
	if (ui==uj || uk==ul) c*=2; else c+=((double) p1)*CE->CoulombElementTweak(ui,uj,ul,uk);
	//if (abs(c-calt)>1e-15) fprintf(logg,"c and calt differ by %e for ui=%i, uj=%i, uk=%i, ul=%i\n",abs(c-calt), ui, uj, ul, uk);
	//else fprintf(logg,"c and calt ok for ui=%i, uj=%i, uk=%i, ul=%i\n", ui, uj, ul, uk);
        //!end-new method
	
        //content c=CoulombElement(*states,ui,uj,uk,ul, CE_exp, 0);
        //if (ui==uj || uk==ul) c*=2; else c+=((double) p1)*CoulombElement(*states,ui,uj,ul,uk, CE_exp,0);

#if IMP_CHECK
        /*if (generuj()*1000>990) {
          content calt1=CoulombElement(*states, ui, uj, uk, ul, CE_exp, 0);
          content calt2=CoulombElement(*states, ui, uj, ul, uk, CE_exp, 0);
          content c1=CE.CoulombElement(ui,uj,uk,ul);
          content c2=CE.CoulombElement(ui,uj,ul,uk);
          if (abs(c1-calt1)/abs(c1+calt1)>1e-10) fprintf(logg, "Coulomb elements differ at 1: ui=%i, uj=%i, uk=%i, ul=%i, c=%.3e%+.3ei, calt=%.3e%+.3ei (ratio-1: %.3e%+.3ei)\n", ui, uj, uk, ul, c1.real(), c1.imag(), calt1.real(), calt1.imag(), (c1/calt1-1.0).real(), (c1/calt1-1.0).imag());
          if (abs(c1-calt1)/abs(c1+calt1)>1e-10) fprintf(logg, "Coulomb elements differ at 2: ui=%i, uj=%i, ul=%i, uk=%i, c=%.3e%+.3ei, calt=%.3e%+.3ei (ratio-1: %.3e%+.3ei)\n", ui, uj, ul, uk, c2.real(), c2.imag(), calt2.real(), calt2.imag(), (c2/calt2-1.0).real(), (c2/calt2-1.0).imag() );
        }*/
#endif
        

        progress_bar->add(0,1);
        progress_bar->print_message(true);
        
        c*=2*(k1->norm)*(k2->norm)*units.meV/units.energy;
        group[g].H[kg1*group[g].N+kg2]=c;
        if (kg1!=kg2) group[g].H[kg2*group[g].N+kg1]=conj(c);
        else group[g].H[kg1*group[g].N+kg2]+=k1->energy;
#if TWOE > 2
        fprintf(logg," Coulomb element: group=%i, kg1=%i->k=%i->(ui=%i, uj=%i, p=%i) kg2=%i->k=%i->(ui=%i, uj=%i, p=%i) : result=(%.5e%+.5ei)\n", g, kg1, group[g].k[kg1], ui, uj, p1, kg2, group[g].k[kg2], uk, ul, p2, c.real(), c.imag());
#endif
      }
    }
#ifdef TWOE_CHECK
  fprintf(logg,"group %i matrix is hermitian with relative precision %.5e\n\n",g, is_hermitian(group[g].N,group[g].H));
#endif
  }
  delete CE;

  timing.c_elems=progress_bar->finished(true);
  progress_bar->reset("Groups diag");
  //for (int g=0;g<Ng;g++) progress_bar->add((int) round(0.4*pow(group[g].N,1.8/2)),0);
  progress_bar->add(Ng,0);
  progress_bar->start();

  //diagonalize all the groups, check if TWOE1
  content* auxval=new content[Ngmax];
  content* auxvec=new content[Ngmax*Ngmax];
  for (int g=0;g<Ng;g++) {
    progress_bar->print_message(true);
#ifdef TWOE_CHECK
    fprintf(logg,"calling the diagonalization for group %i with dimension=%i, eig=%i, matrix hermitian up to %.5e:\n",g, group[g].N, group[g].Neig,is_hermitian(group[g].N,group[g].H));
#endif
    //!new method for diagonalization
    lapack_zgee_prototype lapack(group[g].N);
    lapack.diagonalize(group[g].H, auxval, auxvec, group[g].eigenvecs, false);
    for (int i=0;i<group[g].N;i++) group[g].eigenvals[i]=auxval[i];
    progress_bar->add(0,1);

    //diagonalize_arp(group[g].N, group[g].Neig, group[g].H, group[g].eigenvecs, group[g].eigenvals, (char *) "SM", progress_bar);
    //!end new method
#ifdef TWOE_CHECK
    fprintf(logg,"the obtained eigenvectors are orthonormal with precision:%.5e:\n",is_orthonormal(group[g].N,group[g].N,group[g].eigenvecs));
    fprintf(logg,"the obtained system is the Hamiltonian eigensystem with precision:%.5e:\n",is_eigensys(group[g].N,group[g].N,group[g].H, group[g].eigenvals, group[g].eigenvecs));
#endif
  }
  delete[] auxval;
  delete[] auxvec;
  progress_bar->finished(true);
}

void two_electron_state_prototype::construct_statea(prec par, int *sorted2un, int *unsorted2s)
{
  Na=0;
  for (int g=0;g<Ng;g++) Na+=group[g].Neig;
#if TWOE > 0
  fprintf(logg,"constructing the a-states: %i states will be filled\n",Na);
#endif
  if (statea!=0) {
    if (statea[0].aij!=0) delete statea[0].aij;
    delete statea;
  }
  statea=new statea_prototype[Na];

  int a=0;
  for (int g=0;g<Ng;g++) {
    for (int ag=0;ag<group[g].Neig;ag++) {
      statea[a].ag=ag;
      statea[a].g=g;
#if TWOE > 1
  prec w_max=0;
  int k_max=0;
  for (int kg=0;kg<group[g].N;kg++) {
    prec w=abs(group[g].eigenvecs[ag*group[g].N+kg]);
    if (w>w_max) {w_max=w;k_max=group[g].k[kg];}
  }
  fprintf(logg," a-state #%2i added: group %i, g-eigenvector#%2i=ag (energy=%f [meV]) k-state %i with weight %e\n",a,g,ag,(group[g].eigenvals[ag]).real()*units.energy/units.meV, k_max, w_max);
#endif
      a++;

    }
  }
#if TWOE > 1
  fprintf(logg,"\n");
#endif

  if (a!=Na) {
    sprintf(buffer,"wrong number of a-states in construction: %i instead of %i\n",a,Na);
    splachni(logg,buffer,1+4);
    exit(1);
  }
  
  //identify the lowest singlet and triplet state
  prec minS=0, minT=0;
  int aS=0, aT=0;
  bool firstS=true, firstT=true;
  for (int a=0;a<Na;a++) {
    int g=statea[a].g;
    int ag=statea[a].ag;
    switch (statek[group[g].k[0]].p) {
      case (1) : {
        if (firstS) {
          firstS=false;
          minS=(group[g].eigenvals[ag]).real();
          aS=a;
        }
        if ((group[g].eigenvals[ag]).real()<minS) {
          minS=(group[g].eigenvals[ag]).real();
          aS=a;
        }
        break;
      }
      case (-1) : {
        if (firstT) {
          firstT=false;
          minT=(group[g].eigenvals[ag]).real();
          aT=a;
        }
        if ((group[g].eigenvals[ag]).real()<minT) {
          minT=(group[g].eigenvals[ag]).real();
          aT=a;
        }
        break;
      }
      default: {printf("cannot get here: lowest singlet/triplet\n");exit(1);}
    }
  }

  if (firstS || firstT) {
    fprintf(logg, "lowest singlet and/or triplet not found; unperturbed exchange unquantified\n");
    exchange_un=0;
  }
  else exchange_un=minT-minS;

#if TWOE > 0
  fprintf(logg,"the lowest singlet (a=%i)/triplet (a=%i) found at %e / %e [meV] resulting in unperturbed exchange of %e [meV]\n",aS,aT,minS*units.energy/units.meV, minT*units.energy/units.meV, exchange_un*units.energy/units.meV);
#endif

  //sort the a-states such that lowest S/T is #0/1
  //first swap the found lowest pair
  statea_prototype aux;
  aux=statea[aS];  statea[aS]=statea[0];  statea[0]=aux;
  if (aT==0) aT=aS;
  aux=statea[aT];  statea[aT]=statea[1];  statea[1]=aux;
  
  //sort the rest according to the energy
  for (int a1=2;a1<Na;a1++) {
    for (int a2=2;a2<Na-(a1-2)-1;a2++) {
      prec e1=(group[statea[a2].g].eigenvals[statea[a2].ag]).real();
      prec e2=(group[statea[a2+1].g].eigenvals[statea[a2+1].ag]).real();
      if (e1>e2) {//swap
        aux=statea[a2];  statea[a2]=statea[a2+1];  statea[a2+1]=aux;
      }
    }
  }

#if TWOE > 1
  fprintf(logg," a-state were sorted: Singlet, Triplet, then according to the energy\n");
  for (int a=0;a<Na;a++) {
    int g=statea[a].g;
    int ag=statea[a].ag;
    prec w_max=0;
    int k_max=0;
    for (int kg=0;kg<group[g].N;kg++) {
      prec w=abs(group[g].eigenvecs[ag*group[g].N+kg]);
      if (w>w_max) {w_max=w;k_max=group[g].k[kg];}
    }
    fprintf(logg," a-state #%2i: group %i, g-eigenvector#%2i=ag (energy=%f [meV]) k-state %i with weight %e\n",a,g,ag,(group[g].eigenvals[ag]).real()*units.energy/units.meV, k_max, w_max);
  }
  fprintf(logg,"\n");
#endif

  //sort according to parity
  int type_aux=symmetry.type;
  int type=(int) abs(parameters.read("sorting_method","-"));
  if (type==5) type=3;
  if (type==4) type=0;
  symmetry.type=type;
  level_u* levels=new level_u[Na];
  for (int a=0;a<Na;a++) {
    int g=statea[a].g;
    int ag=statea[a].ag;
    int k=group[g].k[0];
    int i=statek[k].i;
    int j=statek[k].j;
    int p=statek[k].p;
    int u1=selected[i].u;
    int u2=selected[j].u;
    levels[a].parity=symmetry.symmetry2e(u1,u2,p);
    levels[a].value=(group[g].eigenvals[ag]).real();
#if TWOE > 1
    fprintf(logg," a-state %i: (symmetry %i (of type %i), energy %f [meV]) from: group %i, g-eigenvector %2i, k-state %i, (i=%i->ui=%i, j=%i->uj=%i, p=%i)\n", a,levels[a].parity, symmetry.type, levels[a].value*units.energy/units.meV, g, ag, k, i, u1, j, u2, p);
#endif
  }
  symmetry.type=type_aux;
  sorting->inspect_new_2e_states(levels);
  delete[] levels;
  sorting->sort(par, sorted2un, unsorted2s);
#if TWOE > 1
  fprintf(logg," unsorted2s/sorted2un arrays:\n");
  for (int a=0;a<Na;a++) fprintf(logg,"sorted2un[%i]=%i, unsorted2s[%i]=%i\n",a, sorted2un[a], a, unsorted2s[a]);
#endif
  statea_prototype* statea_aux=new statea_prototype[Na];
  for (int a=0;a<Na;a++) statea_aux[a]=statea[a];
  for (int a=0;a<Na;a++) statea[unsorted2s[a]]=statea_aux[a];
  delete[] statea_aux;

#if TWOE > 1
  fprintf(logg," a-state were sorted according to parity:\n");
  for (int a=0;a<Na;a++) {
    int g=statea[a].g;
    int ag=statea[a].ag;
    prec w_max=0;
    int k_max=0;
    for (int kg=0;kg<group[g].N;kg++) {
      prec w=abs(group[g].eigenvecs[ag*group[g].N+kg]);
      if (w>w_max) {w_max=w;k_max=group[g].k[kg];}
    }
    fprintf(logg," a-state #%2i: group %i, g-eigenvector#%2i=ag (energy=%f [meV]) k-state %i with weight %e\n",a,g,ag,(group[g].eigenvals[ag]).real()*units.energy/units.meV, k_max, w_max);
  }
  fprintf(logg,"\n");
#endif

  statea[0].aij=new content[Na*Ns*Ns];
  for (int a=1;a<Na;a++) {
    statea[a].aij=statea[a-1].aij+Ns*Ns;
    for (int i=0;i<Ns*Ns;i++) statea[a].aij[i]=0;
  }

  for (int a=0;a<Na;a++) {
    int g=statea[a].g;
    int ag=statea[a].ag;
    int N=group[g].N;
    for (int kg=0;kg<N;kg++) {
      statek_prototype& sk=statek[group[g].k[kg]];
      content aux=group[g].eigenvecs[ag*N+kg]*sk.norm;
      statea[a].aij[sk.i*Ns+sk.j]+=aux;
      statea[a].aij[sk.j*Ns+sk.i]+=aux*((double) sk.p);
    }
  }
#ifdef TWOE_CHECK
  fprintf(logg,"a-states were expressed in |i x j> basis. They are orthonormal with precision:%.5e:\n",is_orthonormal(Na,Ns*Ns,statea[0].aij));
#endif
  
  double* auxij=new double[Ns*Ns];
  int* index=new int[Ns*Ns];
  for (int i=0;i<Ns;i++) for (int j=0;j<Ns;j++) auxij[i*Ns+j]=(statea[0].aij[i*Ns+j]*conj(statea[0].aij[i*Ns+j])).real();
  picsrt_index(Ns*Ns,auxij,index);
  fprintf(logg,"singlet components sorted in descending order:\n|state > x |state>\t\tweigth (|f|^2)\n");
  for (int k=0;k<20;k++) fprintf(logg,"|%i> x |%i> \t\t\t %e\n", index[k]/Ns, index[k] % Ns, auxij[index[k]]);
	  
}

void two_electron_state_prototype::construct_stateb()
{
  prec muB=zeeman();
  if (parameters.read("gz","-")<0) muB*=-1;
  //find out how many there will be
  Nb=0;
  for (int a=0;a<Na;a++) if (statek[group[statea[a].g].k[0]].p==+1) Nb+=1; else Nb+=3;
#if TWOE > 0
  fprintf(logg,"constructing the b-states: %i states will be added\n",Nb);
#endif
  if (stateb!=0) delete stateb;
  stateb=new stateb_prototype[Nb];
  int b=0;
  for (int a=0;a<Na;a++) {
    if (statek[group[statea[a].g].k[0]].p==+1) {
      stateb[b].energy=(group[statea[a].g].eigenvals[statea[a].ag]).real();
      stateb[b].a=a;
      stateb[b++].Sigma=0;
#if TWOE > 1
      fprintf(logg," b-state #%i added: a-state %i, spinor=%i, energy=%f [meV]\n",b-1,a,stateb[b-1].Sigma, stateb[b-1].energy*units.energy/units.meV);
#endif
    }
    else {
      for (int Sigma=1;Sigma<4;Sigma++) {
        stateb[b].energy=(group[statea[a].g].eigenvals[statea[a].ag]).real();
        stateb[b].a=a;
        stateb[b].Sigma=Sigma;
        stateb[b++].energy+=-muB*2*(Sigma-2);
#if TWOE > 1
        fprintf(logg," b-state #%i added: a-state %i, spinor=%i, energy=%f [meV]\n",b-1,a,stateb[b-1].Sigma, stateb[b-1].energy*units.energy/units.meV);
#endif
      }
    }
  }
#if TWOE > 1
    fprintf(logg,"\n");
#endif
  if (b!=Nb) {
    sprintf(buffer,"wrong number of b-states in construction: %i instead of %i\n",b,Nb);
    splachni(logg,buffer,1+4);
    exit(1);
  }
}

content r2c(prec x, prec y) { return(content(x*x+y*y,0));}

void two_electron_state_prototype::mat_ele_fill_derivs()
{
  progress_bar->reset("SOC m.e.");
  progress_bar->start();
  int N=mat_ele.derivs_length=8;
  progress_bar->add(N*Ns*Ns,0);

  int dim=states->r1->Nw*states->r1->sets;
  content* aux=new content[dim];
  if (mat_ele.derivs!=0) delete mat_ele.derivs;
  mat_ele.derivs=new content[N*Ns*Ns];

  for (int i=0;i<N;i++) {
    states->r1->operate_tot_init();
    for (int j=0;j<states->r1->sets;j++) {
      content im=content(0,1);
      switch (i) {
        case 0 : {states->r1->operate_tot_init(j,j,region::dx,-im,0,exp_integral);break;}
        case 1 : {states->r1->operate_tot_init(j,j,region::dy,-im,0,exp_integral);break;}
        case 2 : {states->r1->operate_tot_init(j,j,region::dxdxdy,im,0,exp_integral);break;}
        case 3 : {states->r1->operate_tot_init(j,j,region::dxdydy,im,0,exp_integral);break;}
        case 4 : {states->r1->operate_tot_init(j,j,region::lin_com,1.0,xc,exp_integral);break;}
        case 5 : {states->r1->operate_tot_init(j,j,region::lin_com,1.0,yc,exp_integral);break;}
        case 6 : {//-i(x dy - y dx)
          states->r1->operate_tot_init(j,j,reg1::dy,-im,xc,exp_integral);
          states->r1->operate_tot_init(j,j,reg1::dx,im,yc,exp_integral);
          break;
        }
        case 7 : {states->r1->operate_tot_init(j,j,region::lin_com,1.0,r2c,exp_integral);break;}
      }
    }
    for (int s1=0;s1<Ns;s1++) {
      for (int s2=0;s2<Ns;s2++) {
        progress_bar->add(0,1);
        progress_bar->print_message(true);

        content* left=states->vysl->eigenvecs+dim*selected[s1].u;
        content* right=states->vysl->eigenvecs+dim*selected[s2].u;
        states->r1->operate_tot(right,aux);
        mat_ele.derivs[(s1*Ns+s2)*N+i]=0;
        for (int j=0;j<states->r1->sets;j++) mat_ele.derivs[(s1*Ns+s2)*N+i]+=states->r1->braket(left, j, aux,j);
      }
    }
  }
  delete[] aux;
  progress_bar->finished(true);

#if TWOE > 1
  prec u=parameters.read("lv","nm")*units.nm/units.length;
  for (int s1=0;s1<Ns;s1++) {
    for (int s2=0;s2<Ns;s2++) {
        fprintf(logg," matrix element of derivatives: <ui=%2i | op |uj=%2i>:\ndx=%5e, dy=%5e, dxdxdy=%5e, dxdydy=%5e, x=%5e, y=%5e, lz=%5e [1/lv, 1/lv^3, lv, 1], r2=%5e [lv]\n",selected[s1].u, selected[s2].u, abs(mat_ele.derivs[(s1*Ns+s2)*N+0])*u, abs(mat_ele.derivs[(s1*Ns+s2)*N+1])*u, abs(mat_ele.derivs[(s1*Ns+s2)*N+2])*u*u*u, abs(mat_ele.derivs[(s1*Ns+s2)*N+3])*u*u*u, abs(mat_ele.derivs[(s1*Ns+s2)*N+4])/u, abs(mat_ele.derivs[(s1*Ns+s2)*N+5])/u, abs(mat_ele.derivs[(s1*Ns+s2)*N+6]), abs(mat_ele.derivs[(s1*Ns+s2)*N+7])/(u*u));
    }
  }
#endif

#ifdef TWOE_CHECK
  aux=new content[Ns*Ns];
  for (int i=0;i<N;i++) {
    for (int s1=0;s1<Ns;s1++) {
      for (int s2=0;s2<Ns;s2++) aux[s1*Ns+s2]=mat_ele.derivs[(s1*Ns+s2)*N+i];
    }
    fprintf(logg,"orbital matrix elements: for operator %i the matrix is hermitian with relative precision %.5e\n",i, is_hermitian(Ns,aux));
  }
  delete[] aux;
#endif

}

void two_electron_state_prototype::mat_ele_fill_orbital()
{
  if (mat_ele.orbital!=0) delete mat_ele.orbital;
  mat_ele.orbital=new content[Ns*Ns*3];

  prec br=parameters.read("br","meVA")*units.meV*1e-10*units.wavevector/units.energy;
  prec dress=parameters.read("dress","meVA")*units.meV*1e-10*units.wavevector/units.energy;
  prec dress3=parameters.read("dress3","eVA^3")*units.e*1e-30*pow(units.wavevector,3)/units.energy;

  for (int s1=0;s1<Ns;s1++) {
    for (int s2=0;s2<Ns;s2++) {
      int label=s1*Ns+s2;

      content kx=mat_ele.derivs[label*mat_ele.derivs_length+0];
      content ky=mat_ele.derivs[label*mat_ele.derivs_length+1];
      content kxkxky=mat_ele.derivs[label*mat_ele.derivs_length+2];
      content kxkyky=mat_ele.derivs[label*mat_ele.derivs_length+3];

      //bychkov-rashba term
      mat_ele.orbital[label*3+0]=ky*br;
      mat_ele.orbital[label*3+1]=-kx*br;
      mat_ele.orbital[label*3+2]=0;

      //if (s1==0 && s2==2) fprintf(logg," br at sx: <ui=%i| w |uj=%i>= %5e\n",selected[0].u,selected[2].u,abs(dy)*br*units.energy/units.meV); 

      //dresselhaus term
      mat_ele.orbital[label*3+0]+=-kx*dress;
      mat_ele.orbital[label*3+1]+=ky*dress;

      //if (s1==0 && s2==2) fprintf(logg," d at sx: <ui=%i| w |uj=%i>= %5e\n",selected[0].u,selected[2].u,abs(dx)*dress*units.energy/units.meV); 

      //cubic dresselhaus term
      mat_ele.orbital[label*3+0]+=kxkyky*dress3;
      mat_ele.orbital[label*3+1]+=-kxkxky*dress3;

      //if (s1==0 && s2==2) fprintf(logg," d3 at sx: <ui=%i| w |uj=%i>= %5e\n",selected[0].u,selected[2].u,abs(dxdydy)*dress3*units.energy/units.meV); 

#if TWOE > 1
    fprintf(logg," matrix element - orbital vector: <ui=%2i | w |uj=%2i>=(%5e, %5e, %5e) [meV]\n",selected[s1].u, selected[s2].u, abs(mat_ele.orbital[label*3+0])*units.energy/units.meV, abs(mat_ele.orbital[label*3+1])*units.energy/units.meV, abs(mat_ele.orbital[label*3+2])*units.energy/units.meV);
#endif
    }
  }

#ifdef TWOE_CHECK
  content* aux=new content[Ns*Ns];
  for (int i=0;i<3;i++) {
    for (int s1=0;s1<Ns;s1++) {
      for (int s2=0;s2<Ns;s2++) {
        aux[s1*Ns+s2]=mat_ele.orbital[(s1*Ns+s2)*3+i];
      }
    }
    fprintf(logg,"orbital matrix elements: for coordinate %i the matrix is hermitian with relative precision %.5e\n",i, is_hermitian(Ns,aux));
  }
  delete[] aux;
#endif

}

//matix element of <Sigma1 | sigma1 \otimes sigma2 | Sigma2>
//where Sigma = 0, 1, 2, 3 ... S, T+, T0, T-
//and pauli   = 0, 1, 2, 3 ... 1, sx, sy, sz
//quantization direction of constituent spins 1/2 is along B
content spinme(int Sigma1, int Sigma2, int pauli1, int pauli2)
{
  prec phi_B=parameters.read("phi_B","pi")*M_PI;
  prec theta_B=parameters.read("theta_B","pi")*M_PI;
  prec c=cos(theta_B/2);
  prec s=sin(theta_B/2);
  content e=exp(content(0,phi_B));

  spinor spin[2];
  spin[0].up=c;
  spin[0].down=s*e;

  spin[1].up=-s;
  spin[1].down=c*e;

  //spinors in basis S=sum_ab C[a*2+b] |a b>, a,b=0/1 
  //singlet, triplets -1, 0, 1
  prec Sigma[4*2*2]={0, 1/sqrt(2), -1/sqrt(2), 0,\
		     1, 0, 0, 0,\
                     0, 1/sqrt(2), 1/sqrt(2), 0,\
                     0, 0, 0, 1};
  //Pauli matrices 1, x, y, z
		     //!EXPERIMENTAL
  /*content pauli[4*4]={1.0, 0, 0, 1.0,\
		      0, 1.0, 1.0, 0,\
		      0, -content(0,1), content(0,1), 0,\
                      1.0, 0, 0, -1.0};*/

  content res=0;
  for (int s1=0;s1<2;s1++) {
    for (int s2=0;s2<2;s2++) {
      prec v1=Sigma[Sigma1*4+s1*2+s2];
      if (v1==0) continue;
      for (int sp1=0;sp1<2;sp1++) {
        for (int sp2=0;sp2<2;sp2++) {
          prec v2=Sigma[Sigma2*4+sp1*2+sp2];
          if (v2==0) continue;
          content aux1=spin[s1].braket(spin[sp1].SigmaTimesMe(pauli+pauli1*4));
          content aux2=spin[s2].braket(spin[sp2].SigmaTimesMe(pauli+pauli2*4));
	  res+=aux1*aux2*v1*v2;
          //fprintf(logg,"adding term for <Sigma1=%+ix%+i | sp=%i | Sigma2=%+ix%+i> with v1=%.3e and v2=%.3e resulting in (%.3e%+.3e) for label=%i\n",i,j,p,k,l,v1,v2,aux.real(), aux.imag(),label);
        }
      }
    }
  }
  return(res);
}

void two_electron_state_prototype::mat_ele_fill_spinor()
{
  if (mat_ele.spinor!=0) delete mat_ele.spinor;
  mat_ele.spinor=new content[4*4*3*2];

  for (int Sigma1=0;Sigma1<4;Sigma1++) {
    for (int Sigma2=0;Sigma2<4;Sigma2++) {
      for (int pauli=0;pauli<3;pauli++) {
	content* me=mat_ele.spinor+(Sigma1*4+Sigma2)*3*2;
        me[pauli]=spinme(Sigma1,Sigma2,pauli+1,0);
        me[3+pauli]=spinme(Sigma1,Sigma2,0,pauli+1);
#if TWOE > 1
        fprintf(logg," matrix element - Pauli matrix vector: <Sigma=%i | sigma x 1|Sigma=%i>\n\tsx=(%+5e%+5ei)\n\tsy=(%+5e%+5ei)\n\tsz=(%+5e%+5ei)\n",Sigma1, Sigma2, me[0].real(), me[0].imag(), me[1].real(), me[1].imag(), me[2].real(), me[2].imag());
        fprintf(logg," matrix element - Pauli matrix vector: <Sigma=%i | 1 x sigma |Sigma=%i>\n\tsx=(%+5e%+5ei)\n\tsy=(%+5e%+5ei)\n\tsz=(%+5e%+5ei)\n",Sigma1, Sigma2, me[0+3].real(), me[0+3].imag(), me[1+3].real(), me[1+3].imag(), me[2+3].real(), me[2+3].imag());
#endif
      }
    }
  }
#ifdef TWOE_CHECK
  content aux[16*2];
  for (int pauli=0;pauli<3;pauli++) {
    for (int Sigma1=0;Sigma1<4;Sigma1++) {
      for (int Sigma2=0;Sigma2<4;Sigma2++) {
        aux[Sigma1*4+Sigma2]=mat_ele.spinor[(Sigma1*4+Sigma2)*3*2+pauli];
        aux[Sigma1*4+Sigma2+16]=mat_ele.spinor[(Sigma1*4+Sigma2)*3*2+pauli+3];
      }
    }
    fprintf(logg,"spinor matrix elements: pauli matrix x 1 %i is hermitian with relative precision %.5e\n",pauli, is_hermitian(4,aux));
    fprintf(logg,"spinor matrix elements: 1 x pauli matrix %i is hermitian with relative precision %.5e\n",pauli, is_hermitian(4,aux+16));
  }
#endif
}


void two_electron_state_prototype::mat_ele_fill_spin_imp()
{
  progress_bar->reset("nuclear m.e.");
  progress_bar->start();
  progress_bar->add(Ns*Ns,0);
  
  if (mat_ele.spin_imp!=0) delete mat_ele.spin_imp;
  mat_ele.spin_imp=new content[Ns*Ns*3];
  int N=states->r1->Nw;
  int length=states->r1->Nw*states->r1->sets;
  for (int i=0;i<Ns;i++) {
    content* psi1=states->vysl->eigenvecs+selected[i].u*length;
    if (selected[i].u>states->unsorted) 
      {printf("wrong coversion: s=%i to u=%i vs eig=%i\n",i,selected[i].u,states->unsorted);exit(1);}
    for (int j=0;j<Ns;j++) {
      content* psi2=states->vysl->eigenvecs+selected[j].u*length;
      if (selected[j].u>states->unsorted) 
      {printf("wrong coversion: s=%i to u=%i vs eig=%i\n",j,selected[j].u,states->unsorted);exit(1);}
      content overlap[3]={0};
      for (int n=0; n<N; n++) {
        content aux=conj(psi1[n])*psi2[n];
        for (int x=0;x<3;x++) overlap[x]+=aux*spin_impurities.I[x][n];
      }
      for (int x=0;x<3;x++) {
        mat_ele.spin_imp[(i*Ns+j)*3+x]=overlap[x];
#if TWOE > 2
        fprintf(logg,"selected states matrix element due to impurities: <Psi(sel=%i, u=%i) | I[%i] | Psi(sel=%i, u=%i) > = %.3e%+.3e\n",i, selected[i].u, x, j,selected[j].u, overlap[x].real(), overlap[x].imag());
#endif        
      }
      progress_bar->add(0,1);
    }
  }
  progress_bar->finished(true);
}

void two_electron_state_prototype::mat_ele_convert(content *me1e, content *me2e, const char* mess)
{
  
  progress_bar->reset(mess);
  progress_bar->start();
  progress_bar->add(Na*Na,0);
  
  //int pntr=0;
  
  for (int a1=0;a1<Na;a1++) {
    for (int a2=0;a2<Na;a2++) {
      content* res=me2e+(a1*Na+a2)*3;
      res[0]=res[1]=res[2]=0;
      for (int i=0;i<Ns;i++) {
	for (int j=0;j<Ns;j++) {
	  content c1=conj(statea[a1].aij[i*Ns+j]);
	  if (c1.real()==0 && c1.imag()==0) continue;//{pntr++;continue;}
	  for (int ip=0;ip<Ns;ip++) {
	    content c1c2=c1*statea[a2].aij[ip*Ns+j];
	    content* me=me1e+(i*Ns+ip)*3;	      
	    res[0]+=c1c2*me[0];
	    res[1]+=c1c2*me[1];
	    res[2]+=c1c2*me[2];
	  }
	}
      }
      progress_bar->add(0,1);
    }
  }
  progress_bar->finished(true);
  //printf("perc. of zeros %e\n",(100.0*pntr)/(Na*Na*Ns*Ns));
}

content two_electron_state_prototype::mat_ele_from2e(int b1, int b2, content* me2e)
{
  content *orb=me2e+(stateb[b1].a*Na+stateb[b2].a)*3;
  content *spin=mat_ele.spinor+(stateb[b1].Sigma*4+stateb[b2].Sigma)*3*2;
  return(2.0*(orb[0]*spin[0]+orb[1]*spin[1]+orb[2]*spin[2]));
}

content two_electron_state_prototype::mat_ele_from1e(int b1, int b2, content* me1e)
{
  int a1=stateb[b1].a;
  int a2=stateb[b2].a;
  int ag1=statea[a1].ag;
  int ag2=statea[a2].ag;
  int g1=statea[a1].g;
  int g2=statea[a2].g;
  int N1=group[g1].N;
  int N2=group[g2].N;
  
  content res[3]={0,0,0};
  for (int kg1=0;kg1<N1;kg1++) {
    int k1=group[g1].k[kg1];
    content ck1p=conj(group[g1].eigenvecs[ag1*N1+kg1])*statek[k1].norm;
    int i1=statek[k1].i, j1=statek[k1].j;
    prec p1=statek[k1].p;
    for (int kg2=0;kg2<N2;kg2++) {
      int k2=group[g2].k[kg2];
      content ck2=group[g2].eigenvecs[ag2*N2+kg2]*statek[k2].norm;
      int i2=statek[k2].i, j2=statek[k2].j;
      prec p2=statek[k2].p;
      if (i1!=i2 && i1!=j2 && j1!=i2 && j1!=j2) continue;
      for (int x=0;x<3;x++) {
        content aux=0;
        if (j1==j2) aux+=me1e[(i1*Ns+i2)*3+x];
        if (i1==j2) aux+=me1e[(j1*Ns+i2)*3+x]*p1;
        if (j1==i2) aux+=me1e[(i1*Ns+j2)*3+x]*p2;
        if (i1==i2) aux+=me1e[(j1*Ns+j2)*3+x]*p1*p2;
        res[x]+=aux*ck1p*ck2;
      }
    }
  }
  content sum=0;
  for (int x=0;x<3;x++) sum+=res[x]*2.0*mat_ele.spinor[(stateb[b1].Sigma*4+stateb[b2].Sigma)*3*2+x];
  return(sum);
}

/*
content two_electron_state_prototype::mat_ele_spin_imp(int b1, int b2)
{
  int a1=stateb[b1].a;
  int a2=stateb[b2].a;
  int ag1=statea[a1].ag;
  int ag2=statea[a2].ag;
  int g1=statea[a1].g;
  int g2=statea[a2].g;
  int N1=group[g1].N;
  int N2=group[g2].N;
  
  content res_alt[3]={0,0,0};
  for (int kg1=0;kg1<N1;kg1++) {
    int k1=group[g1].k[kg1];
    content ck1p=conj(group[g1].eigenvecs[ag1*N1+kg1])*statek[k1].norm;
    int i1=statek[k1].i, j1=statek[k1].j;
    prec p1=statek[k1].p;
    for (int kg2=0;kg2<N2;kg2++) {
      int k2=group[g2].k[kg2];
      content ck2=group[g2].eigenvecs[ag2*N2+kg2]*statek[k2].norm;
      int i2=statek[k2].i, j2=statek[k2].j;
      prec p2=statek[k2].p;
      if (i1!=i2 && i1!=j2 && j1!=i2 && j1!=j2) continue;
      for (int x=0;x<3;x++) {
        if (j1==j2) res_alt[x]+=mat_ele.spin_imp[(i1*Ns+i2)*3+x]*ck1p*ck2;
        if (i1==j2) res_alt[x]+=mat_ele.spin_imp[(j1*Ns+i2)*3+x]*ck1p*ck2*p1;
        if (j1==i2) res_alt[x]+=mat_ele.spin_imp[(i1*Ns+j2)*3+x]*ck1p*ck2*p2;
        if (i1==i2) res_alt[x]+=mat_ele.spin_imp[(j1*Ns+j2)*3+x]*ck1p*ck2*p1*p2;
      }
    }
  }
  content sum=0;
  for (int x=0;x<3;x++) {
   //sum+=res_alt[x]*2.0*spinme(stateb[b1].Sigma, stateb[b2].Sigma, x+1, 0);
   sum+=res_alt[x]*2.0*mat_ele.spinor[(stateb[b1].Sigma*4+stateb[b2].Sigma)*3*2+x];
  }
#if TWOE > 2
  if (b1==b2 || abs(sum)>1e-6) fprintf(logg,"b-state matrix element due to impurities: <a=%i, S=%i | I | a=%i, S=%i > = %.3e%+.3e\n",a1, stateb[b1].Sigma, a2, stateb[b2].Sigma, sum.real(), sum.imag());
#endif
  return(sum);
}*/


void two_electron_state_prototype::construct_statef()
{
  //prec muB=zeeman();

  //allocate space for the final Hamiltonian
  if (Hb!=0) delete Hb;
  Hb=new content[Nb*Nb];

  progress_bar->reset("final matrix");
  progress_bar->start();
  progress_bar->add(Nb*(Nb+1)/2,0);
  //5/2 for the impurity spin, 1/2 for the hole spin
  prec Arho=-parameters.read("exch","meV")*units.meV/units.energy*parameters.read("xMn", "-")*(5.0/2)*(1.0/2);
  if (spin_impurities.I[0]==0 || spin_impurities.N==0) {
    fprintf(logg,"no spin impurities present\n");
    Arho=0;
  } else fprintf(logg,"spin impurities present\n");

  //build the Hamiltonian: 
  for (int b1=0;b1<Nb;b1++) {
    for (int b2=0;b2<Nb;b2++) {
      int label=b1*Nb+b2;
      if (b2<b1) {//Hermitian
        Hb[label]=conj(Hb[b2*Nb+b1]);
        continue;
      }
      else {
        Hb[label]=mat_ele_from2e(b1,b2,mat_ele.orbital2e); //spin-orbit matrix element
	//content aux=mat_ele_soc(b1,b2); //spin-orbit matrix element
	//if (abs(aux-Hb[label])>1e-15) fprintf(logg," |aux-Hb[label]| = %e for b1=%i and b2=%i\n", abs(aux-Hb[label]), b1, b2);

        if (Arho!=0) {
	  //content m1=mat_ele_spin_imp(b1,b2)*Arho;//spin-impurity matrix element
	  //content m2=mat_ele_from2e(b1,b2,mat_ele.spin_imp2e)*Arho;
	  //if (abs(m1-m2)>1e-20) fprintf(logg," |m1-m2| = %e for b1=%i and b2=%i\n", abs(m1-m2), b1, b2);
	  Hb[label]+=mat_ele_from2e(b1,b2,mat_ele.spin_imp2e)*Arho;
	}
        if (b1==b2) {//the diagonal: a-state energy + Zeeman energy
          Hb[label]+=stateb[b1].energy;
        }
        progress_bar->add(0,1);
        progress_bar->print_message(true);
      }
#if TWOE > 2
  fprintf(logg,"final matrix element for b states : | <b1=%i | H |b2=%i> | =%e [meV]\n", b1, b2, abs(Hb[label])*units.energy/units.meV);
#endif
    }
  }
#ifdef TWOE_CHECK
  fprintf(logg,"the final matrix is hermitian with relative precision %.5e\n", is_hermitian(Nb,Hb));
#endif
  timing.final_matrix=progress_bar->finished(true);

  //sanity check: Nf should be from 1 to Nb
  if (Nf<1)  {fprintf(logg,"Nf was reset from %i to 1 in construct_statef\n",Nf);Nf=1;}
  if (Nf>Nb) {fprintf(logg,"Nf was reset from %i to %i (number of b-states) in construct_statef\n",Nf,Nb);Nf=Nb;}
  //degeneracy check - in that case enlarge Nf
  int Nf_reset=Nf;
  prec diag=0, offdiag=0;
  do {
    if (Nf_reset==Nb || Nf_reset-Nf>20) break;
    offdiag=abs(Hb[(Nf_reset-1)*Nb+Nf_reset]);
    diag=abs(Hb[(Nf_reset-1)*Nb+Nf_reset-1]-Hb[(Nf_reset)*Nb+Nf_reset]);
    if (diag>1e-4 || offdiag>1) break;
    Nf_reset++;
#if TWOE > 0
    fprintf(logg,"\nexpanding Nf by one to %i due to possible degeneracy issue: (distance to the next diagonal %.3e, offdiagonal element %.3e)\n",Nf_reset,diag, offdiag);
#endif  
  } while (true);
#if TWOE > 0
  if (Nf_reset!=Nf) fprintf(logg,"the number of Nf states will be enlarged by %i to %i due to possible degeneracy issue: (distance to the next diagonal %.3e, offdiagonal element %.3e)\n",Nf_reset-Nf,Nf_reset,diag, offdiag);
  else fprintf(logg,"the number of Nf states (%i) looks safe wrt possible degeneracy: distance to the next diagonal %.3e\n",Nf,diag);
#endif  
  Nf=Nf_reset;

  //allocate space for eigenvectors
  content* vecs=new content[Nf*Nb];
  content* vals=new content[Nf];

  //diagonalize Hb by arpack
  progress_bar->reset("final diag");
  progress_bar->add((int) round(0.4*pow(Nb,1.8/2)),0);
  progress_bar->start();

  //TraceOn('z',1,1,1,1);
  diagonalize_arp(Nb, Nf, Hb, vecs, vals, (char* ) "SM",progress_bar);
  timing.final_diag=progress_bar->finished(true);
  //TraceOff();


#ifdef TWOE_CHECK
  fprintf(logg,"The final matrix was diagonalized. The obtained eigenvectors are orthonormal with precision:%.5e:\n",is_orthonormal(Nb,Nf,vecs));
  fprintf(logg,"the obtained system is the Hamiltonian eigensystem with precision:%.5e:\n",is_eigensys(Nb,Nf,Hb,vals,vecs));
#endif

  //progress_bar->reset("computing statistics");
  //progress_bar->start();
  if (statef!=0) {
    delete statef[0].cfkS;
    delete statef[0].prob_g;
    delete statef;
  }
  statef=new statef_prototype[Nf];

  //sort f-states according to b-states of maximal weight
  //find a maximal weight for each f-state
  for (int f=0;f<Nf;f++) {
    prec wmax=0;
    for (int b=0;b<Nb;b++) {
      prec w=(vecs[f*Nb+b]*conj(vecs[f*Nb+b])).real();
      if (w>wmax) {
        bool already=false;
        for (int f2=0;f2<f;f2++) if (statef[f2].b==b) already=true;
        if (!already) {wmax=w;statef[f].b=b;}
      }
    }
#if TWOE > 1
    fprintf(logg,"f-state %i is closest to b-state %i with weight %e\n",f,statef[f].b,wmax);
#endif
  }
  //sort acsending according to the max weight labels
  int* pos=new int[Nf];
  for (int f=0;f<Nf;f++) pos[f]=f;

  for (int f1=0;f1<Nf;f1++) {
    for (int f2=0;f2<Nf-f1-1;f2++) {
      if (statef[f2].b>statef[f2+1].b) {//swap
        int aux=pos[f2]; pos[f2]=pos[f2+1]; pos[f2+1]=aux;
        aux=statef[f2].b; statef[f2].b=statef[f2+1].b; statef[f2+1].b=aux;
      }
    }
  }
#if TWOE > 1
    fprintf(logg,"f-states will be placed according to i->pos[i]\n");
    for (int f=0;f<Nf;f++) fprintf(logg,"pos[%i]=%i\n",f,pos[f]);
    fprintf(logg,"\n");
#endif

  //fill in f-states
  //statef[0].cfb=new content[Nf*Nb];
  for (int f=0;f<Nf;f++) {
    statef[f].cfb=vecs+pos[f]*Nb;
    statef[f].energy=vals[pos[f]].real();
  }
  if (statef_vecs_begin!=0) delete statef_vecs_begin;
  statef_vecs_begin=vecs;
  delete[] pos;
  delete[] vals;

  //express f-states in the k-state x Sigma basis
  statef[0].cfkS=new content[Nf*Nk*4];
  for (int i=0;i<Nf*Nk*4;i++) statef[0].cfkS[i]=0;

  for (int f=0;f<Nf;f++) {
    statef[f].cfkS=statef[0].cfkS+Nk*4*f;
    for (int b=0;b<Nb;b++) {
      content cfb=statef[f].cfb[b];
      int a=stateb[b].a;
      int Sigma=stateb[b].Sigma;
      int g=statea[a].g;
      int ag=statea[a].ag;
      int N=group[g].N;
      for (int kg=0;kg<N;kg++) {
        int k=group[g].k[kg];
        content cak=group[g].eigenvecs[ag*N+kg];
        statef[f].cfkS[k*4+Sigma]+=cak*cfb;
        if ((statek[k].p==-1 && Sigma==0) || (statek[k].p==1 && Sigma>0)) {
          sprintf(buffer, "wrong exchange symmetry combination in statef: Sigma=%i, p=%i (k-state k=%i, group g=%i, eigenstate kg=%i)\n",Sigma, statek[k].p, k, g, kg);
          splachni(logg,buffer,1+4);
          exit(1);
        }
      }
    }
  }

  //compute statistics - mean values of spin and spectral decomposition of symmetry operators
  statef[0].prob_g=new prec[Ng*Nf];
  for (int f=0;f<Nf;f++) {
    statef[f].prob_g=statef[0].prob_g+f*Ng;
    for (int g=0;g<Ng;g++) statef[f].prob_g[g]=0;
    statef[f].S2=0;
    statef[f].SB=0;

    for (int k=0;k<Nk;k++) {
      int g=statek[k].g;
      for (int Sigma=0;Sigma<4;Sigma++) {
        prec prob=(statef[f].cfkS[k*4+Sigma]*conj(statef[f].cfkS[k*4+Sigma])).real();
        statef[f].prob_g[g]+=prob;
        if (statek[k].p==-1) statef[f].S2+=prob;
        if (Sigma>0) statef[f].SB+=-prob*(Sigma-2);
      }
    }
    statef[f].p_max=statef[f].prob_g[0];
    statef[f].g_max=0;
    for (int g=1;g<Ng;g++) 
      if (statef[f].prob_g[g]>statef[f].p_max) {
        statef[f].p_max=statef[f].prob_g[g];
        statef[f].g_max=g;
      }
  }

#ifdef TWOE_CHECK
  //check the orthonormality of f-states in both bases
  fprintf(logg,"The f-states are orthonormal with precision %e in b-state basis and %e in k-state basis\n", is_orthonormal(Nb,Nf,statef_vecs_begin), is_orthonormal(Nk*4,Nf,statef[0].cfkS));
#endif
  //progress_bar->finished(true);
}

//print the energies and symmetry info for f-states:
//file, where: 1 - stdout, 2 - file, 4 - log
//what - 1 energies - 2 characteristics 4-heading (means also a new line, otherwise not)
void two_electron_state_prototype::statistics(const char* st, int what, FILE* file, int where, int Nmax)
{
  bool energies=(what % 2);
  what/=2;
  bool chars=(what % 2);
  what/=2;
  bool heading= (what % 2);

  int choice;
  if (strcmp(st,"b")==0) choice=0;
  else if (strcmp(st,"f")==0) choice=1;
  else if (strcmp(st,"f-b")==0) choice=2;
  else {
    sprintf(buffer, "wrong choice in two electron statistics: %s\n",st);
    splachni(logg,buffer,1+4);
    exit(1);
  }
  
  const char* Snames[]={"S ", "T+", "T0", "T-"};
  switch (choice) {
    case 0 : {
      if (heading) {
	sprintf(buffer, "b-state\t");
	splachni(file, buffer, where);
	if (energies) {
	  sprintf(buffer, "energy[meV]\t\t");
	  splachni(file, buffer, where);
	}
	if (chars) {
   	  sprintf(buffer, "Sigma\torbital symmetry\n");
	  splachni(file, buffer, where);
        }
      }
      int N=Nb;
      if (Nmax!=-1) N=min(Nb,Nmax);
	
      for (int b=0; b<N; b++) {
	if (heading) {
	  sprintf(buffer, "%2i\t", b);
	  splachni(file, buffer, where);
	}
	if (energies) {
	  sprintf(buffer, "%.14e\t", stateb[b].energy*units.energy/units.meV);
	  splachni(file, buffer, where);
	}
	if (chars) {
          int a=stateb[b].a;
          int g=statea[a].g;
	  sprintf(buffer, "%s\t%i(group %i)\t", Snames[stateb[b].Sigma],  symmetry.symmetry2einv(group[g].id),g);
	  splachni(file, buffer, where);
	}
	if (heading) {
	  sprintf(buffer, "\n");
	  splachni(file, buffer, where);
	}
      }
      break;
    }
    case 2 : ;
    case 1 : {
      if (heading) {
	sprintf(buffer, "f-state\t");
	splachni(file, buffer, where);
	if (energies) {
	  sprintf(buffer, "energy[meV]\t\t");
	  splachni(file, buffer, where);
	}
	if (chars) {
   	  sprintf(buffer,"Sigma^2[hbar^2]\t\tSigma.B[hbar]\t\torb. symm.\t its precision\t\t~b-state\tsoc shift\n");
	  splachni(file, buffer, where);
        }
      }

      int N=Nf;
      if (Nmax!=-1) N=min(Nf,Nmax);

      for (int f=0; f<N; f++) {
	if (heading) {
	  sprintf(buffer, "%2i\t", f);
	  splachni(file, buffer, where);
	}
	if (energies) {
          prec energy=statef[f].energy;
          if (choice==2) energy-=stateb[statef[f].b].energy;
	  sprintf(buffer, "%.14e\t", energy*units.energy/units.meV);
	  splachni(file, buffer, where);
	}
	if (chars) {
	  sprintf(buffer, "%.14e\t%.14e\t%i(group %i)\t%.14e\t%s (%i)\t%.14e", statef[f].S2, statef[f].SB, symmetry.symmetry2einv(group[statef[f].g_max].id),statef[f].g_max, statef[f].p_max, Snames[stateb[statef[f].b].Sigma], statef[f].b,  (statef[f].energy-stateb[statef[f].b].energy)*units.energy/units.meV);
	  splachni(file, buffer, where);
	}
	if (heading) {
	  sprintf(buffer, "\n");
	  splachni(file, buffer, where);
	}
      }
      break;
    }
  }
}

prec two_electron_state_prototype::read(char states, int s, const char* tag)
{
  int quantity;
  if (strcmp(tag,"energy")==0) quantity=0;
  else if (strcmp(tag,"S2")==0) quantity=1;
  else if (strcmp(tag,"SB")==0) quantity=2;
  else {
    printf("wrong quatity asked for in two electron read: %s\n",tag);
    exit(1);
  }
  switch (states) {
    case ('b') : {
      if (s<0 || s>=Nb) {printf("no b-state #%i in two electron read\n",s);exit(1);}
      if (quantity==0) return(stateb[s].energy*units.energy/units.meV);
      if (quantity==1) return(1-stateb[s].Sigma==0);
      if (quantity==2) {
        if (stateb[s].Sigma==0) return(0);
        return(-(stateb[s].Sigma-2));
      }
      break;
    }
    case ('f') : {
      if (s<0 || s>=Nf) {printf("no f-state #%i in two electron read\n",s);exit(1);}
      if (quantity==0) return(statef[s].energy*units.energy/units.meV);
      if (quantity==1) return(statef[s].S2);
      if (quantity==2) return(statef[s].SB);      
      break;
    }
    default : {printf("no %c states in two electron read\n",states);exit(1);}
  }
  return(0);
}

content two_electron_state_prototype::read(char states, int s1, int s2, int particle, const char* tag)
{
  int quantity;
  if (strcmp(tag,"energy")==0) return(read(states, s1, tag));
  if (strcmp(tag,"dipolex")==0) quantity=4;
  else if (strcmp(tag,"dipoley")==0) quantity=5;
  else if (strcmp(tag,"Lz")==0) quantity=6;
  else {
    printf("wrong quatity asked for in two electron read: %s\n",tag);
    exit(1);
  }
  switch (states) {
    case ('b') : {
      if (min(s1,s2)<0 || max(s1,s2)>=Nb) {printf("no b-states #%i and #%i in two electron read\n",s1,s2);exit(1);}
      //if (stateb[s1].Sigma!=stateb[s2].Sigma) return(0);
      content sum=0;
      int g1=statea[stateb[s1].a].g, g2=statea[stateb[s2].a].g;
      int ag1=statea[stateb[s1].a].ag, ag2=statea[stateb[s2].a].ag;
      content* c1=group[g1].eigenvecs+ag1*group[g1].N;
      for (int kg1=0;kg1<group[g1].N;kg1++) {
        int k1=group[g1].k[kg1];
        int i1=statek[k1].i, j1=statek[k1].j, p1=statek[k1].p;
        content* c2=group[g2].eigenvecs+ag2*group[g2].N;
        for (int kg2=0;kg2<group[g2].N;kg2++) {
          content aux=0;
          int k2=group[g2].k[kg2];
          int i2=statek[k2].i, j2=statek[k2].j, p2=statek[k2].p;
          if (j1==j2) aux+=mat_ele.derivs[(i1*Ns+i2)*mat_ele.derivs_length+quantity];
          if (i1==j2) aux+=mat_ele.derivs[(j1*Ns+i2)*mat_ele.derivs_length+quantity]*double(p1);
          if (j1==i2) aux+=mat_ele.derivs[(i1*Ns+j2)*mat_ele.derivs_length+quantity]*double(p2);
          if (i1==i2) aux+=mat_ele.derivs[(j1*Ns+j2)*mat_ele.derivs_length+quantity]*double(p1*p2);
          sum+=aux*statek[k1].norm*statek[k2].norm*conj(c1[kg1])*c2[kg2];
        }
      }
      if (particle==2) sum*=double(statek[group[g1].k[0]].p*statek[group[g2].k[0]].p);
      if (quantity==4 || quantity==5) sum*=units.length/units.nm;
      return(sum);
      break;
    }
    case ('f') : {
      if (min(s1,s2)<0 || max(s1,s2)>=Nf) {printf("no f-states #%i and #%i in two electron read\n",s1,s2);exit(1);}
      content sum=0;
      int g1=statea[stateb[s1].a].g, g2=statea[stateb[s2].a].g;
      int ag1=statea[stateb[s1].a].ag, ag2=statea[stateb[s2].a].ag;
      for (int k1=0;k1<Nk;k1++) {
        int i1=statek[k1].i, j1=statek[k1].j, p1=statek[k1].p;
        for (int Sigma1=0;Sigma1<4;Sigma1++) {
          content c1=statef[s1].cfkS[k1*4+Sigma1];
          for (int k2=0;k2<Nk;k2++) {
            int i2=statek[k2].i, j2=statek[k2].j, p2=statek[k2].p;
            int Sigma2=Sigma1;
            content c2=statef[s2].cfkS[k2*4+Sigma2];
            content aux=0;
            if (j1==j2) aux+=mat_ele.derivs[(i1*Ns+i2)*mat_ele.derivs_length+quantity];
            if (i1==j2) aux+=mat_ele.derivs[(j1*Ns+i2)*mat_ele.derivs_length+quantity]*double(p1);
            if (j1==i2) aux+=mat_ele.derivs[(i1*Ns+j2)*mat_ele.derivs_length+quantity]*double(p2);
            if (i1==i2) aux+=mat_ele.derivs[(j1*Ns+j2)*mat_ele.derivs_length+quantity]*double(p1*p2);
            sum+=aux*statek[k1].norm*statek[k2].norm*conj(c1)*c2;
          }
        }
      }
      if (quantity==4 || quantity==5) sum*=units.length/units.nm;
      return(sum);
      break;
    }
    default : {printf("no %c states in two electron read\n",states);exit(1);}
  }
  return(0);
}


void two_electron_state_prototype::spinatn(prec *s, int n)
{
  if (n<0 || n>=states->r1->Nw) {printf("wrong seq label (%i) in spinatn (Nw=%i)\n",n,states->r1->Nw); exit(1);}
  /*content SS[16*4];
  for (int Sigma1=0;Sigma1<4;Sigma1++) 
    for (int Sigma2=0;Sigma2<4;Sigma2++) {
      for (int i=0;i<4;i++) SS[(Sigma1*4+Sigma2)*4+i]=spinme(Sigma1, Sigma2, 0, i);//times 1/2 for hole spin for components 1-3 (1 for component 0)
      double delta;
      if (Sigma1==Sigma2) delta=1; else delta=0;
      if (abs(SS[(Sigma1*4+Sigma2)*4+0]-delta)>1e-14) {printf("difference 1 %e\n", delta);exit(1);}
      for (int i=1;i<4;i++) {
        if (abs(SS[(Sigma1*4+Sigma2)*4+i]-mat_ele.spinor[(Sigma1*4+Sigma2)*6+3+i-1])>1e-14) {printf("difference 2 i=%i\n",i);exit(1);}
      }
    }*/

  content sum[4]={0,0,0,0}, sum_alt[4]={0,0,0,0};
  int length=states->r1->Nw*states->r1->sets;
  for (int i1=0;i1<Ns;i1++) {
    content psi1=states->vysl->eigenvecs[selected[i1].u*length+n];
    for (int i2=0;i2<Ns;i2++) {
      content psi2=states->vysl->eigenvecs[selected[i2].u*length+n];
      for (int Sigma1=0;Sigma1<4;Sigma1++) {
        for (int Sigma2=0;Sigma2<4;Sigma2++) {
          content rhoij=rho[(i1*Ns+i2)*16+Sigma1*4+Sigma2];
          //for (int i=0;i<4;i++) {
            //sum[i]+=2.0*rhoij*psi1*conj(psi2)*SS[(Sigma2*4+Sigma1)*4+i];
          //}
          content x=rhoij*psi1*conj(psi2);
          if (Sigma1==Sigma2) sum_alt[0]+=x;
          content* aux=mat_ele.spinor+(Sigma2*4+Sigma1)*6+3-1;
          sum_alt[1]+=x*aux[1];
          sum_alt[2]+=x*aux[2];
          sum_alt[3]+=x*aux[3];
        }
      }
    }
  }
  for (int i=0;i<4;i++) {
    s[i]=sum_alt[i].real()*2;
    //if (abs(s[i]-sum[i])>1e-14) {printf("difference 3 %i\n", i);exit(1);}
  }
  #ifdef TWOE_CHECK
  for (int i=0;i<4;i++) if (abs(sum[i].imag())>1e-8*abs(s[i])) fprintf(logg,"warning: imaginary part of local spin %e too large wrt real part %e in spinatn()\n",sum[i].imag(), s[i]);
  #endif
}

//rho is the density, nmax is the initial guess, steps is the number of tries
double two_electron_state_prototype::MC_max_search(content* rho, int guess, int steps)
{  
  //Monte carlo search
  //find the interpolated maximum
  content ** derivs=new content*[6];
  states->r1->value_at_init(derivs,rho);
  //the crude position of the maximum
  prec xmax0, ymax0, hx, hy;
  states->r1->s2coordinate(guess, xmax0, ymax0);
  states->r1->give_par(hx,hy);
  //do Monte Carlo search in the vicinity of crude maximum
  prec xmax=xmax0;
  prec ymax=ymax0;
  prec fmax=states->r1->value_at(xmax0,ymax0,derivs).real();
  for (int i=0;i<steps;i++) {
    prec x=xmax+hx*(2*generuj()-1)/5;
    prec y=ymax+hy*(2*generuj()-1)/5;
    prec f=states->r1->value_at(x,y,derivs).real();
    //fprintf(logg,"searching: x=%e, y=%e, f=%e, fmax=%e, cmax=%e, xmax=%e, ymax=%e\n",x,y,f,fmax,cmax, xmax, ymax);
    if (f>fmax) {
      fmax=f;
      xmax=x;
      ymax=y;
    }
  }
  prec Rmax=sqrt(xmax*xmax+ymax*ymax)*units.length/units.nm;
  fprintf(logg,"searching result: xmax=%e, ymax=%e, Rmax=%e, fmax=%e, guess=%i\n",xmax,ymax,Rmax,fmax,guess);
  states->r1->value_at_release(derivs);
  delete[] derivs;
  return(Rmax);
}

void two_electron_state_prototype::distribution_characteristics(int I)
{
  prec s[4], x, y, cmax=0, Qxx=0, Qyy=0, Qxy=0;
  prec dx=0, dy=0;
  int nmax=0, nSmax=0;
  for (int i=0;i<I;i++) M[i]=0;
  content * rho=new content[states->r1->Nw*2];
  content * rhoS=rho+states->r1->Nw;
  for (int n=0;n<states->r1->Nw;n++) {
    spinatn(s,n);
    states->r1->s2coordinate(n,x,y);
    prec r=sqrt(x*x+y*y);
    prec phi=atan(y/x);
    if (x==0 && y==0) phi=0;
    if (x<0) phi+=M_PI;
    content phase=exp(content(0,phi));
    content factor=1.0;
    prec c=s[0];
    rho[n]=c;
    if (c>cmax) {cmax=c;nmax=n;Rmax=r*units.length/units.nm;}
    rhoS[n]=s[3];
    if (rhoS[n].real()>rhoS[nSmax].real()) nSmax=n;
    //fprintf(logg,"x=%e, y=%e, n=%i, c=%e, cmax=%e, nmax=%i\n",x,y,n,c,cmax,nmax);  
    for (int i=0;i<I;i++) {
      M[i]+=c*factor;
      factor*=phase;
    }
    dx+=c*x;
    dy+=c*y;
    //if (abs(x-r*phase.real())+abs(y-r*phase.imag())>1e-12) fprintf(logg,"x=%e vs Re(rip)=%e, y=%e vs Im(rip)=%e\n",x,r*phase.real(), y, r*phase.imag());
    Qxx+=x*x*c;
    Qxy+=x*y*c;
    Qyy+=y*y*c;
  }
  //for (int i=1;i<I;i++) M[i]/=M[0];
  //fprintf(logg,"abs(M[1])=%e vs d=%e\n",abs(M[1]), sqrt(dx*dx+dy*dy)/abs(M[0]));
  Qsph=sqrt((Qxx-Qyy)*(Qxx-Qyy)+4*Qxy*Qxy) / (Qxx+Qyy);
  
  Rmax=MC_max_search(rho, nmax, 100000);
  RSmax=MC_max_search(rhoS, nSmax, 100000);
  delete[] rho;
  
  /*if (true) {//Monte carlo search
    //find the interpolated maximum
    content ** derivs=new content*[6];
    states->r1->value_at_init(derivs,rho);
    //the crude position of the maximum
    prec xmax0, ymax0, hx, hy;
    states->r1->s2coordinate(nmax, xmax0, ymax0);
    states->r1->give_par(hx,hy);
    //do Monte Carlo search in the vicinity of crude maximum
    prec xmax=xmax0;
    prec ymax=ymax0;
    prec fmax=states->r1->value_at(xmax0,ymax0,derivs).real();
    for (int i=0;i<100000;i++) {
      prec x=xmax+hx*(2*generuj()-1)/5;
      prec y=ymax+hy*(2*generuj()-1)/5;
      prec f=states->r1->value_at(x,y,derivs).real();
      //fprintf(logg,"searching: x=%e, y=%e, f=%e, fmax=%e, cmax=%e, xmax=%e, ymax=%e\n",x,y,f,fmax,cmax, xmax, ymax);
      if (f>fmax) {
        fmax=f;
        xmax=x;
        ymax=y;
      }
    }
    Rmax=sqrt(xmax*xmax+ymax*ymax)*units.length/units.nm;
    fprintf(logg,"searching result: xmax=%e, ymax=%e, Rmax=%e, fmax=%e, cmax=%e\n",xmax,ymax,Rmax,fmax,cmax);
    states->r1->value_at_release(derivs);
    delete derivs;
    delete rho;
    //exit(0);
  }*/
}

void two_electron_state_prototype::reducedDM(int f1, int f2)
{
  if (f1>=Nf || f2>=Nf || f1<0 || f2<0) {
    printf("reduced density matrix asked for non-existing f-state (f1=%i, f2=%i, Nf=%i)\n", f1, f2, Nf);
    exit(1);
  }
  #if TWOE > 0
  fprintf(logg,"reduced density matrix asked for states <f#%i | f#%i>\n", f1, f2);
  #endif
  
  if (rho==0) rho=new content[Ns*Ns*16];
  for (int i=0;i<Ns*Ns*16;i++) rho[i]=0;
#ifdef TWOE_CHECK
  content* rho2=new content[Ns*Ns*16];
  for (int i=0;i<Ns*Ns*16;i++) rho2[i]=0;
#endif

  for (int k1=0;k1<Nk;k1++) {
    int i1=statek[k1].i, j1=statek[k1].j;
    prec p1=statek[k1].p;
    for (int k2=0;k2<Nk;k2++) {
      int i2=statek[k2].i, j2=statek[k2].j;
      prec p2=statek[k2].p;
      content aux=0;
      for (int Sigma1=0;Sigma1<4;Sigma1++) {
        for (int Sigma2=0;Sigma2<4;Sigma2++) { 
          content aux=statef[f1].cfkS[k1*4+Sigma1]*conj(statef[f2].cfkS[k2*4+Sigma2]);
          aux*=statek[k1].norm*statek[k2].norm;

          if (i1==i2) rho[(j1*Ns+j2)*16+Sigma2*4+Sigma1]+=aux;
          if (j1==i2) rho[(i1*Ns+j2)*16+Sigma2*4+Sigma1]+=aux*p1;
          if (i1==j2) rho[(j1*Ns+i2)*16+Sigma2*4+Sigma1]+=aux*p2;
          if (j1==j2) rho[(i1*Ns+i2)*16+Sigma2*4+Sigma1]+=aux*p1*p2;
          #ifdef TWOE_CHECK
          if (j1==j2) rho2[(i1*Ns+i2)*16+Sigma2*4+Sigma1]+=aux;
          if (i1==j2) rho2[(j1*Ns+i2)*16+Sigma2*4+Sigma1]+=aux*p1;
          if (j1==i2) rho2[(i1*Ns+j2)*16+Sigma2*4+Sigma1]+=aux*p2;
          if (i1==i2) rho2[(j1*Ns+j2)*16+Sigma2*4+Sigma1]+=aux*p1*p2;
          #endif
        }
      }
    }
  }

#ifdef TWOE_CHECK
  //check the trace and equivalence
  content tr=0, diff=0, norm=0;
  for (int s1=0;s1<Ns;s1++) {
    for (int s2=0;s2<Ns;s2++) {
      for (int Sigma1=0;Sigma1<4;Sigma1++) {
        for (int Sigma2=0;Sigma2<4;Sigma2++) {
          content aux1=rho[(s1*Ns+s2)*16+Sigma1*4+Sigma2];
          content aux2=rho2[(s1*Ns+s2)*16+Sigma1*4+Sigma2];
          if (Sigma1==0) aux2*=-1.0;
          if (Sigma2==0) aux2*=-1.0;
          diff+=aux1-aux2;
          norm+=aux1+aux2;
          if (Sigma1==Sigma2 && s1==s2) tr+=aux1;
        }
      }
    }
  }
  if (f1==f2) tr-=1.0;
  fprintf(logg,"reduced DM: trace residuum %e, alternative residuum %e\n", abs(tr), abs(diff/norm));
  delete[] rho2;
#endif
}

content two_electron_state_prototype::rdotr(int f1, int f2)
{
  if (f1>=Nf || f2>=Nf || f1<0 || f2<0) {
    printf("r.r asked for non-existing f-state (f1=%i, f2=%i, Nf=%i)\n", f1, f2, Nf);
    exit(1);
  }
  int N=mat_ele.derivs_length;
  
  content rr=0, r2=0, unity=0;
  for (int k1=0;k1<Nk;k1++) {
    int i1=statek[k1].i, j1=statek[k1].j;
    prec p1=statek[k1].p;
    for (int k2=0;k2<Nk;k2++) {
      int i2=statek[k2].i, j2=statek[k2].j;
      prec p2=statek[k2].p;

      content aux_rr=0, aux_r2=0;
      for (int xy=4;xy<6;xy++) {
	aux_rr+=mat_ele.derivs[(i1*Ns+i2)*N+xy]*mat_ele.derivs[(j1*Ns+j2)*N+xy];
	aux_rr+=p1*mat_ele.derivs[(j1*Ns+i2)*N+xy]*mat_ele.derivs[(i1*Ns+j2)*N+xy];
	aux_rr+=p2*mat_ele.derivs[(i1*Ns+j2)*N+xy]*mat_ele.derivs[(j1*Ns+i2)*N+xy];
	aux_rr+=p1*p2*mat_ele.derivs[(j1*Ns+j2)*N+xy]*mat_ele.derivs[(i1*Ns+i2)*N+xy];
      }

      if (j1==j2) aux_r2+=mat_ele.derivs[(i1*Ns+i2)*N+7];
      if (i1==j2) aux_r2+=p1*mat_ele.derivs[(j1*Ns+i2)*N+7];
      if (j1==i2) aux_r2+=p2*mat_ele.derivs[(i1*Ns+j2)*N+7];
      if (i1==i2) aux_r2+=p1*p2*mat_ele.derivs[(j1*Ns+j2)*N+7];

      //if (abs(aux_rr)>1e-7) fprintf(logg,"contribution to rr: k1=%i, k2=%i, aux_rr=%.3e%+.3ei\n", k1, k2, aux_rr.real(), aux_rr.imag());
      //if (abs(aux_r2)>1e-7) fprintf(logg,"contribution to r2: k1=%i, k2=%i, aux_r2=%.3e%+.3ei\n", k1, k2, aux_r2.real(), aux_r2.imag());

      for (int Sigma=0;Sigma<4;Sigma++) {
        content amp=conj(statef[f1].cfkS[k1*4+Sigma])*statef[f2].cfkS[k2*4+Sigma];
	if (k1==k2) unity+=amp;
        amp*=statek[k1].norm*statek[k2].norm;
      	rr+=amp*aux_rr;
	r2+=amp*aux_r2;
      }
    }
  }
  #if TWOE > 1
  fprintf(logg,"r.r asked for states <f#%i | f#%i>, getting rr=%.3e%+.3ei, r2=%.3e%+.3ei, rr/r2=%.3e%+.3ei, unity=%e (normalization)\n", f1, f2, rr.real(), rr.imag(), r2.real(), r2.imag(), (rr/r2).real(), (rr/r2).imag(), unity.real());
  #endif
  
  return(content( (rr/r2).real(), r2.real() ) );
}


//phonon induced relaxation f1 -> f2
void two_electron_state_prototype::relaxation(int f1, int f2, prec* result)
{
  if (f1>=Nf || f2>=Nf || f1<0 || f2<0) {
    sprintf(buffer, "relaxation rate asked for non-existing f-state (f1=%i, f2=%i, Nf=%i)\n", f1, f2, Nf);
    splachni(logg, buffer, 1+4);
    return;
  }
  
#if TWOE > 0
  sprintf(buffer,"relaxation rates (f#%i->f#%i) at energy difference %e [meV]\n", f1, f2, (statef[f1].energy-statef[f2].energy)*units.energy/units.meV);
  splachni(logg,buffer,4);
#endif  

  //first build matrices g and h
  content* g12=new content[Ns*Ns];
  for (int i=0;i<Ns*Ns;i++) g12[i]=0;
#ifdef TWOE_CHECK
  content* h12=new content[Ns*Ns];
  for (int i=0;i<Ns*Ns;i++) h12[i]=0;
#endif

  for (int k1=0;k1<Nk;k1++) {
    int i1=statek[k1].i, j1=statek[k1].j;
    prec p1=statek[k1].p;
    for (int k2=0;k2<Nk;k2++) {
      int i2=statek[k2].i, j2=statek[k2].j;
      prec p2=statek[k2].p;
      content aux=0;
      for (int Sigma=0;Sigma<4;Sigma++) aux+=conj(statef[f1].cfkS[k1*4+Sigma])*statef[f2].cfkS[k2*4+Sigma];
      aux*=statek[k1].norm*statek[k2].norm;

      if (j1==j2) g12[i1*Ns+i2]+=aux;
      if (i1==j2) g12[j1*Ns+i2]+=aux*p1;
      if (j1==i2) g12[i1*Ns+j2]+=aux*p2;
      if (i1==i2) g12[j1*Ns+j2]+=aux*p1*p2;
#ifdef TWOE_CHECK
      if (j1==j2) h12[i1*Ns+i2]+=aux*p1*p2;
      if (i1==j2) h12[j1*Ns+i2]+=aux*p2;
      if (j1==i2) h12[i1*Ns+j2]+=aux*p1;
      if (i1==i2) h12[j1*Ns+j2]+=aux;
#endif
    }
  }

#ifdef TWOE_CHECK
  //check the trace and equivalence
  content trg=0, trh=0, diff=0, norm=0;
  for (int s=0;s<Ns;s++) {
    trg+=g12[s*Ns+s];
    trh+=h12[s*Ns+s];
    diff+=g12[s*Ns+s]-h12[s*Ns+s];
    norm+=abs(g12[s*Ns+s])+abs(h12[s*Ns+s]);
  }
  if (f1==f2) {trg-=1.0;trh-=1.0;}
  fprintf(logg,"trace of g=%e, trace of h=%e, they are equal with precison %e\n", abs(trg), abs(trh), abs(diff/norm));
  delete[] h12;
#endif

  //construct the "wavefunctions"
  int wfl=states->r1->Nw*states->r1->sets;
  content cst=content(1,0);//content(1.0/sqrt(wfl),0);
  content* wf=new content[wfl*2];
  for (int n=0;n<wfl;n++) {
    wf[n]=0;
    wf[wfl+n]=cst;
    for (int s1=0;s1<Ns;s1++) {
      int offset1=wfl*selected[s1].u;
      for (int s2=0;s2<Ns;s2++) {
        int offset2=wfl*selected[s2].u;
        wf[n]+=(2.0*g12[s1*Ns+s2])*conj(states->vysl->eigenvecs[offset1+n])*states->vysl->eigenvecs[offset2+n];
      }
    }
  }
  delete[] g12;

  //prepare the space for auxiliary "wavefunctions"
  //two wavefunctions will be needed
  content en0=states->vysl->eigenvals[0];
  content en1=states->vysl->eigenvals[1];
  content* wf_old=states->vysl->eigenvecs;
  
  states->vysl->eigenvecs=wf;
  states->vysl->eigenvals[0]=content(statef[f2].energy,0);
  states->vysl->eigenvals[1]=content(statef[f1].energy,0);
  
  for (int i=0;i<4;i++) {
#if TWOE > 0
  sprintf(buffer,"computing rate #%i\t",i);splachni(logg,buffer,4);
#endif  
    result[i]=transition_rate(*states,1,0,i);
#if TWOE > 0
  sprintf(buffer,"with result:%e\n",result[i]);splachni(logg,buffer,4);
#endif  
  }
 
  //delete the auxiliary wavefunctions
  delete[] wf;
  states->vysl->eigenvecs=wf_old;
  states->vysl->eigenvals[0]=en0;
  states->vysl->eigenvals[1]=en1;
  return;
}

















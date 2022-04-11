#include "main.h"

inline void nl2npnm(int n, int l, int & np, int &nm)
{
  nm=np=n;
  if (l>0) np+=l; else nm-=l;
}

inline void npnm2nl(int &n, int &l, int  np, int nm)
{
  l=np-nm;
  n=(np+nm-abs(l))/2;
}


prec energy(int n, int l, int sigma)
{
  prec B=parameters.read("B","Tesla");
  prec theta_B=parameters.read("theta_B","pi")*M_PI;
  prec lB=parameters.read("lB","nm")*units.nm;
  prec gz=parameters.read("gz","-");
  prec m_eff=parameters.read("m_eff","m_e");
  
  prec aux1=pow(units.hbar,2)/(units.m_e*m_eff*lB*lB);
  prec aux2=units.hbar*units.e/(2*units.m_e);
  //in SI units
  prec res=aux1*(2*n+abs(l)+1)+B*cos(theta_B)*(l/m_eff+gz*sigma/2)*aux2;
  //in dimensionless units
  res/=units.energy;
  //sprintf(buffer, "energy: quantum numbers n=%+i, l=%+i, s=%+i, aux1=%e, aux2=%e,  energy[meV]=%+e\n",n,l,sigma,aux1*units.energy/units.meV,aux2*units.energy/units.meV,res*units.energy/units.meV);splachni(logg,buffer,4);
  return(res);
}

prec energy(prec n, prec l, int sigma)
{
  prec B=parameters.read("B","Tesla");
  prec theta_B=parameters.read("theta_B","pi")*M_PI;
  prec lB=parameters.read("lB","nm")*units.nm;
  prec gz=parameters.read("gz","-");
  prec m_eff=parameters.read("m_eff","m_e");
  
  prec aux1=pow(units.hbar,2)/(units.m_e*m_eff*lB*lB);
  prec aux2=units.hbar*units.e/(2*units.m_e);
  //in SI units
  prec res=aux1*(2*n+abs(l)+1)+B*cos(theta_B)*(l/m_eff+gz*sigma/2)*aux2;
  //in dimensionless units
  res/=units.energy;
  //sprintf(buffer, "energy: quantum numbers n=%+.2f, l=%+.2f, s=%+i, energy=%+e\n",n,l,sigma,res*units.energy/units.meV);splachni(logg,buffer,4);
  return(res);
}

prec zeeman()
{
  prec B=parameters.read("B","Tesla");
  prec theta_B=parameters.read("theta_B","pi")*M_PI;
  prec phi_B=parameters.read("phi_B","pi")*M_PI;
  prec Bx=sin(theta_B)*cos(phi_B)*B;
  prec By=sin(theta_B)*sin(phi_B)*B;
  prec Bz=cos(theta_B)*B;
  prec gx=parameters.read("gx","-");
  prec gy=parameters.read("gy","-");
  prec gz=parameters.read("gz","-");
  prec res=sqrt(gx*Bx*gx*Bx+gy*By*gy*By+gz*Bz*gz*Bz)*units.e*units.hbar/units.m_e/4;
  return(res/units.energy);
}


//correction to the state |n,l,s> according to corr
//1 : BR-BR
//2 : D - D
//4 : BR-D3
//8 : D -D3
//16: D3-D3
//31: all of them
prec correction(int corr,prec n, prec l, int s)
{
  prec B=parameters.read("B","Tesla");
  prec theta_B=parameters.read("theta_B","pi")*M_PI;
  prec lB=parameters.read("lB","nm")*units.nm;
  prec lv=parameters.read("lv","nm")*units.nm;
  prec gz=parameters.read("gz","-");
  prec m_eff=parameters.read("m_eff","m_e");
  prec br=parameters.read("br","meVA")*units.meV*1e-10;
  prec dress=parameters.read("dress","meVA")*units.meV*1e-10;
  prec dress3=parameters.read("dress3","eVA^3")*units.e*1e-30;
  
  prec ad=dress/pow(units.hbar,2)*2*units.m_e*m_eff*lv;
  prec abr=br/pow(units.hbar,2)*2*units.m_e*m_eff*lv;
  prec ad3=dress3*units.m_e*m_eff/pow(units.hbar,2)/lv;
  prec E0=pow(units.hbar,2)/(2*units.m_e*m_eff*lv*lv);

  prec np,nm;
  prec theta,tp1,tm1,tp2,tm2,az2t,tot=0,err;

  nm=np=n;
  if (l>0) np+=l; else nm-=l;

  theta=B*cos(theta_B)*units.e*lB*lB/(2*units.hbar);//theta=B*cos(theta_B)/(2*lB*lB);
  tp1=(theta+1);
  tm1=(theta-1);
  tp2=tp1*tp1;
  tm2=tm1*tm1;
  az2t=2*theta*gz*m_eff/2;//az2t=2*theta*gz;


  if (corr%2==1) {//BR-BR
    if (s==1) {
      err=(1+np)*tp2/(-tp1+az2t);
      err+=nm*tm2/(-tm1+az2t);
    }
    else {
      err=-(1+nm)*tm2/(-tm1+az2t);
      err+=-np*tp2/(-tp1+az2t);
    }
    err/=2;
    tot+=abr*abr*err;
  }
  corr/=2;
  if (corr%2==1) {//D-D
    if (s==1) {
      err=(1+nm)*tm2/(tm1+az2t);
      err+=np*tp2/(tp1+az2t);
    }
    else {
      err=-(1+np)*tp2/(tp1+az2t);
      err+=-nm*tm2/(tm1+az2t);
    }
    err/=2;
    tot+=ad*ad*err;
  }

  corr/=2;
  if (corr%2==1) {//BR-D3
    err=0;
    tot+=abr*ad3*err;
  }
  corr/=2;
  if (corr%2==1) {//D-D3
    if (s==1) {
      err=-(1+nm)*tm1*(nm*tm2+2*(1+np*tp2+theta*theta))/((tm1+az2t)*tp1);
      err-=np*tp1*(np*tp2+(1+2*nm)*tm2)/((tp1+az2t)*tm1);
    }
    else {
       err=(1+np)*tp1*(2*nm*tm2+np*tp2+2*(1+theta*theta))/((tp1+az2t)*tm1);
      err+=nm*tm1*(nm*tm2+(1+2*np)*tp2)/((tm1+az2t)*tp1);
    }
    err*=(-1+theta*theta)/(2*pow(lB/lv,2));
    tot+=ad3*ad*err;
  }
  corr/=2;
  if (corr%2==1) {//D3-D3
    if (s==1) {

      err=nm*(np-1)*np*tm2*tp2*tp2/(3+theta+az2t);
      err+=np*tp2*(np*tp2+tm2*(1+2*nm))*(np*tp2+tm2*(1+2*nm))/(1+theta+az2t);
      err+=nm*(nm-1)*(nm-2)*tm2*tm2*tm2/(3-3*theta+az2t);
       err+=(1+nm)*tm2*(nm*tm2+2*(1+np*tp2+theta*theta))*(nm*tm2+2*(1+np*tp2+theta*theta))/(-1+theta+az2t);
      err+=9*(nm-1)*nm*(np+1)*tm2*tm2*tp2/(1-3*theta+az2t);
      err+=(nm+1)*(nm+2)*(np+1)*tm2*tm2*tp2/(-3+theta+az2t);
      err+=9*nm*(np+1)*(np+2)*tm2*tp2*tp2/(-1-3*theta+az2t);
      err+=(np+1)*(np+2)*(np+3)*tp2*tp2*tp2/(-3-3*theta+az2t);

      err/=8*pow(lB,4);
    }
    else {

      err=-np*(np-1)*(np-2)*tp2*tp2*tp2/(-3-3*theta+az2t);
      err+=-9*np*(np-1)*(nm+1)*tm2*tp2*tp2/(-1-3*theta+az2t);
      err+=-nm*np*(nm-1)*tm2*tm2*tp2/(-3+theta+az2t);
      err+=-9*(nm+1)*(nm+2)*np*tm2*tm2*tp2/(1-3*theta+az2t);
      err+=-nm*tm2*(nm*tm2+(1+2*np)*tp2)*(nm*tm2+(1+2*np)*tp2)/(-1+theta+az2t);
      err+=-(nm+1)*(nm+2)*(nm+3)*tm2*tm2*tm2/(3-3*theta+az2t);
      err+=-(1+np)*tp2*(np*tp2+2*(nm*tm2+1+theta*theta))*(np*tp2+2*(nm*tm2+1+theta*theta))/(1+theta+az2t);
      err+=-(nm+1)*(np+2)*(np+1)*tm2*tp2*tp2/(3+theta+az2t);
      err/=8*pow(lB,4);
    }
    tot+=ad3*ad3*err;
  }
  //sprintf(buffer, "correction returns %e\n",tot);splachni(logg,buffer,4);
  return(tot*E0);
}

prec level(int which, int corr, int& n, int &l, int &s, bool wc)
{
  #define _how_many 10000
  static prec mins[_how_many];
  static int arr_n[_how_many];
  static int arr_l[_how_many];
  static int arr_s[_how_many];

  if (which<0) { //compute all values
    //in n,l, are maximal values up which to go
    int count=0;
    for (int i=0;i<=n;i++) {
      for (int j=-l;j<=l;j++) {
        for (int k=-1;k<2;k+=2) {
	  arr_n[count]=i;
	  arr_l[count]=j;
	  arr_s[count]=k;
	  mins[count]=energy(i,j,k);
          if (wc) mins[count]+=correction(corr,i,j,k);
	  count++;
	}
      }
    }

    prec tempp;
    int tempi;
    for (int i=1;i<count;i++) {
      for (int j=0;j<count-i;j++) {
        if (mins[j]<=mins[j+1]) continue;
	tempp=mins[j];mins[j]=mins[j+1];mins[j+1]=tempp;
        tempi=arr_n[j];arr_n[j]=arr_n[j+1];arr_n[j+1]=tempi;
        tempi=arr_l[j];arr_l[j]=arr_l[j+1];arr_l[j+1]=tempi;
        tempi=arr_s[j];arr_s[j]=arr_s[j+1];arr_s[j+1]=tempi;
      }
    }
    if (! wc) {
      for (int i=0;i<count;i++) mins[i]+=correction(corr,arr_n[i],arr_l[i],arr_s[i]);
    }
    //sprintf(buffer, "level returns zero\n");splachni(logg,buffer,4);
    return(0);
  }
  else {//now read the values previously computed values
    n=arr_n[which-1];
    l=arr_l[which-1];
    s=arr_s[which-1];
    //sprintf(buffer, "level returns %e\n",mins[which-1]);splachni(logg,buffer,4);
    return(mins[which-1]);
  }
}


void levels(FILE* file, int where_to,int corr, int hmn,int n, int l, int numbers, bool wc)
{
  prec val;
  int np,nm;

  bool n_out=false, l_out=false,np_out=false, nm_out=false, s_out=false;
  if (numbers%2==1) n_out=true;
  numbers/=2;
  if (numbers%2==1) l_out=true;
  numbers/=2;
  if (numbers%2==1) np_out=true;
  numbers/=2;
  if (numbers%2==1) nm_out=true;
  numbers/=2;
  if (numbers%2==1) s_out=true;
  numbers/=2;


  int s;
  level(-1,corr,n,l,s,wc);

    for (int i=0;i<hmn;i++) {
      val=level(i+1,corr,n,l,s,wc);
      nl2npnm(n,l,np,nm);
      sprintf(buffer,"%.10e ",val);
      splachni(file, buffer,where_to);
      if (n_out) {
        sprintf(buffer,"%i ",n);
        splachni(file, buffer,where_to);
      }
      if (l_out) {
        sprintf(buffer,"%i ",l);
        splachni(file, buffer,where_to);
      }
      if (np_out) {
        sprintf(buffer,"%i ",np);
        splachni(file, buffer,where_to);
      }
      if (nm_out) {
        sprintf(buffer,"%i ",nm);
        splachni(file, buffer,where_to);
      }
      if (s_out) {
        sprintf(buffer,"%i ",s);
        splachni(file, buffer,where_to);
      }
    }
    message(file,"\n",where_to);
}

//computes correction to energy in the case of a ring
//which = 0 magnetic field independent part
//which = 1 field dependent part
prec ring_correction(prec lz, int which)
{
  prec B=parameters.values[parameters_class::B];
  prec theta_B=parameters.values[parameters_class::theta_B];
  prec lv=parameters.values[parameters_class::lv];
  prec ring=parameters.values[parameters_class::ring];
  prec m_eff=parameters.values[parameters_class::m_eff];
  
  prec res;
  switch (which) {
    case (0) : {
      //transversal LHO ground state
      res=1/pow(lv,2);
      //field independent part of the kinetic energy (first order pert.)
      res+=lz*lz*(0.5*pow(lv/ring,2)+0.75*pow(lv/ring,4))/pow(ring,2);
      break;
    }
    case (1) : {
      //field dependent part of the kinetic energy (first order pert.)
      prec flux=M_PI*ring*ring*cos(theta_B)*B/(2*M_PI);
      res=pow(flux/ring,2)*1.5*pow(lv/ring,2);
      break;
    }
    default : {printf("wrong choice 'which'=%i in ring_correction\n",which);exit(1);}
  }
  //to proper units
  prec mult=units.hbar*units.hbar/(2*units.m_e*m_eff*units.length*units.length)/units.energy;
  return(res*mult);
}

int glob_n,glob_l,glob_al;
prec nrm,rho20m;
prec x0;
content phase;

content anal_aux(prec x, prec y)
{
  prec B=parameters.values[parameters_class::B];
  prec theta_B=parameters.values[parameters_class::theta_B];

  x-=x0;
  int& n=glob_n;
  int& l=glob_l;
  int& a=glob_al;

  prec res;
  content res2,res3=content(1,0);
  prec rho2=(x*x+y*y)*rho20m;

  if (n==0) res=1;
  else if (n==1) res=1+a-rho2;
  else if (n==2) res=1.0/2*(2+3*a+a*a-4*rho2-2*a*rho2+rho2*rho2);
  else {
    message("too high n in laguerre in anal_aux1()!!!");
    exit(1);
  }
  res3=content(res*exp(-rho2/2)*nrm,0);
  if (l>0) res2=content(x,y); else res2=content(x,-y);
  for (int i=0;i<a;i++) res3*=res2;
  if (B*cos(theta_B)!=0) res3*=exp(y*phase);
  return(res3);
}

void norm_function(reg1 &r1,content *output)
{
  prec nrm=0;
  for (int i=0;i<r1.Nw*r1.sets;i++) nrm+=norm(output[i]);
  nrm=1/sqrt(nrm);
  for (int i=0;i<r1.Nw*r1.sets;i++) output[i]*=nrm;
}


void anal_function(reg1& r1,content *output,int n, int l, int s, prec d_0,prec sign, reg1::setting_or_adding soa)
{
  prec B=parameters.values[parameters_class::B];
  prec theta_B=parameters.values[parameters_class::theta_B];
  prec lB=parameters.values[parameters_class::lB];

  x0=d_0;
  glob_n=n;
  glob_l=l;
  glob_al=abs(l);
  if (s==-1) s=1; else s=0;

  prec hx,hy;
  r1.give_par(hx,hy);

  rho20m=1/(lB*lB);
  nrm=M_PI/rho20m;
  for (int i=n+1;i<n+glob_al;i++) nrm*=i;
  nrm=sign/sqrt(nrm*hx*hy);

  phase=content(0,B*cos(theta_B)/2/rho20m*d_0);

  r1.operate(0,output,0,s,reg1::lin_com,content(1,0),anal_aux,soa);
  r1.operate(0,output,0,1-s,reg1::lin_com,content(0,0),0,soa);
}


//zapise do vysl(eigenfuctions aj eigenvals) eig najnizsich rieseni single dot problemu
void anal_spectrumSD(ret_from_arpack_zn & vysl,reg1& r1,int eig, prec d_0)
{
  level(-1,0,eig,eig,eig,false);
  int n,l,s;

  for (int i=0;i<eig;i++) {
    vysl.eigenvals[i]=level(i+1,0,n,l,s,false);
    anal_function(r1,vysl.eigenvecs+i*r1.Nw*r1.sets,n,l,s,d_0,1,reg1::set);
    norm_function(r1,vysl.eigenvecs+i*r1.Nw*r1.sets);
  }
}

bool g_function(reg1& r1, content *output,int g,int n,int l,int spin,prec d_0)
{
  int s[]={1,1,1,1,1,-1,1,-1,1,1,-1,-1,1,-1,-1,1};
  if ((l==0) && (g==2 || g==3)) return(false);

  int pnt=(g-1)*4;
  if (abs(l)%2==1) pnt=(4-g)*4;

  anal_function(r1,output,n,l,spin,d_0,s[pnt],reg1::set);
  anal_function(r1,output,n,l,spin,-d_0,s[1+pnt],reg1::add);
  anal_function(r1,output,n,-l,spin,-d_0,s[2+pnt],reg1::add);
  anal_function(r1,output,n,-l,spin,d_0,s[3+pnt],reg1::add);

  norm_function(r1,output);

  //for (int i=0;i<4;i++) printf("sign %i",s[pnt+i]);
  //printf("\nand normed - group %i\n",g);

return(true);
}

//zapise do vysl(eigenfuctions) eig rieseni double dot problemu usporiadanych podla poradia v d->\inf
void anal_spectrumDD(ret_from_arpack_zn & vysl,reg1& r1,int eig, prec d_0)
{
  level(-1,0,eig,eig,eig,false);
  int n,l,s;
  int pnt=0,g=4,hmn=0;
  content *pntr=vysl.eigenvecs;

  do {
    if (g==4) do level((++pnt),0,n,l,s,false); while(l>0);
    //printf("from level at pnt=%i extracted n=%i,l=%i,s=%i\n",pnt,n,l,s);
    do {g++;if (g==5) g=1;} while (! g_function(r1,pntr,g,n,l,s,d_0));
    //printf("function written at g=%i\n",g);
    pntr+=r1.Nw*r1.sets;
    hmn++;
  } while (hmn<eig);
}

int give_parity(state_info* states, int u, prec eps, int method)
{
  char system=parameters.values[parameters_class::system];
  
  //content *pntr=vysl.eigenvecs+r1.sets*r1.Nw*eig;
  int par=0;
  bool spin;
  if (method>=1) spin=false; else {spin=true;method*=-1;}
  
  prec p1=0,p2=0,p3=0;
  if (method==1 || method==2) {
    p1=states->read(u,"Ix");
    p2=states->read(u,"Iy");
    p3=states->read(u,"I");
    if (system=='g') p3=states->read(u,"Isigmaz");
  }
  
  switch (method) {
    case (0): {//nothing
      par=0;
      break;  
    }
      
    case (2): {//all inversions
      if (fabs(1-p1)<eps && fabs(1-p2)<eps && fabs(1-p3)<eps) par=1;
      else if (fabs(-1-p1)<eps && fabs(1-p2)<eps && fabs(-1-p3)<eps) par=2;
      else if (fabs(-1-p1)<eps && fabs(-1-p2)<eps && fabs(1-p3)<eps) par=3;
      else if (fabs(1-p1)<eps && fabs(-1-p2)<eps && (-1-p3)<eps) par=4;
      break;
    }
    case (1): {//only inv
      if (fabs(1-p3)<eps) par=1;
      else if (fabs(-1-p3)<eps) par=2;
      break;
    }
    case (5): {//np, nm
      prec np=states->read(u,"np");
      prec nm=states->read(u,"nm");
      
      if (fabs(my_round(np)-np)<eps || fabs(my_round(nm)-nm)<eps) 
        par=(int) (my_round(np)+10)*100+(my_round(nm)+10);
      else par=0;
      break;
    }
    case (3): ;//lz
    case (6): {////graphene total jz
      prec L;
      if (method==3) L=states->read(u,"lz"); else L=states->read(u,"g-jz");
      par=(int) my_round(L);
      if (fabs(L-par)>eps) par=0;//no definite parity
      else {//definite parity found
        if (par>=0) par=par*2+1; else par*=-2;
      }
      break;
    }
    case (4): {//nothing (maybe only spin) 
      par=1;
      break;
    }
      
    default: {printf("wrong method in parity choosen!!!\n");}
  }
  if (spin) {
    prec s=states->read(u,"s");    
    if (fabs(1-s)<eps) par*=1;
    else if (fabs(-1-s)<eps) par*=-1;
    else par=0;
  }
  //for graphene only - discriminate positive and negative energy states
  if (system=='g' && par!=0) if(states->read(u,"energy")<0) par+=1000;

  return(par);
}

//write info about a state (label sorted - s)
//format: +1 - header +2-maximal precision (9) +4 new line at the end
//what: defined in state info, added by tags by show_set(tag)
//show_set(clear) = clear all flags
void statistics(state_info& I, int s, int format,  FILE* file, int where_to)
{
  bool header=false;
  bool maxprec=false;
  bool newline=false;
  char bufe[10]="%+.14e ";

  format/=1;if (format%2) header=true;
  format/=2;if (format%2) maxprec=true;
  format/=2;if (format%2) newline=true;
  if (!maxprec) bufe[3]='0';
  
  const char* tag;
  const char* name;
  const char* unit;
  char str[100];
  char buffer[10000];
  
  int u=I.sorted2un[s];
  if (header) sprintf(buffer,"%2i.sorted (%+2i. unsorted) ",s,u);
  else buffer[0]=0;

  for (int t=0;t<I.length;t++) {
    if (! I.show[t]) continue;
    I.tag_info(t,tag,name,unit);
    if (header) {
      strcat(buffer,tag);strcat(buffer,"[");strcat(buffer,unit);strcat(buffer,"]=");
    }
    prec val=0;
    if (u!=-1) val=I.read(u,t);
    sprintf(str,bufe,val);
    strcat(buffer,str);
  }
  if (newline) strcat(buffer,"\n");
  splachni(file,buffer,where_to);
}

//writes for all states info defined by tags in state_info 
//format: +1 - header +2-maximal precision +4 new line after every eigenvecs +8 newline at the end +16 only visible states
//what_global (encoded binary): 
//1 zeeman energy 2 tunneling energy 4 geff
//8 times.main 16 times.user 32 perform.iter 64 perform.opx
void statistics_multi(state_info& I, int what_global, int format,  FILE* file, int where_to)
{  
  prec B=parameters.read("B","Tesla");
  prec m_eff=parameters.read("m_eff","m_e");
  
  char bufe[10]="+.14e ";
  bool only_visible=false, header=false;
  
  if (format%2) header=true;
  if ((format/2)%2) bufe[3]='0';
  if ((format/16)%2) only_visible=true;
  
  prec globals[7];
  globals[0]=I.read(I.sorted2un[1],"energy")-I.read(I.sorted2un[0],"energy");
  globals[1]=I.read(I.sorted2un[2],"energy")-I.read(I.sorted2un[0],"energy");
  globals[2]=globals[0]*units.meV/(B*units.hbar*units.e/(2*units.m_e));
  globals[3]=I.vysl->times.main;
  globals[4]=I.vysl->times.user;
  globals[5]=I.vysl->perform.iter;
  globals[6]=I.vysl->perform.opx;
  
  char buf[100];
  if ( header) {
    const char* tag;
    const char* name;
    const char* unit;
    
    message(file,"##. ",where_to);
    for (int t=0;t<I.length;t++) {
      if (!I.show[t]) continue;
      I.tag_info(t,tag,name,unit);
      strcpy(buf,tag);strcat(buf,"[");strcat(buf,unit);strcat(buf,"] ");
      sprintf(buffer,"%11s ",buf);
      splachni(file, buffer,where_to);
    }
    message(file,"\n",where_to);
  }

  for (int s=0;s<I.sorted;s++) {
    if (I.sorted2un[s]==-1 && only_visible) continue; 
    statistics(I,s,format, file, where_to);
  }
  buf[0]=buffer[0]=0;
  
  for (int i=0;i<7;i++) {
    if (what_global%2) {
      sprintf(buf,bufe,globals[i]);
      strcat(buffer,buf);
      what_global/=2;
    }
  }
  if ((format/2)%2) strcat(buffer,"\n");
  splachni(file, buffer,where_to);
}

//#define EXTRAPOLATE
void extrapolate_coefs(prec *coors, prec *vals, int method, prec *coefs)
{
  //y=ax^2+bx+c;
  prec a, b, c;
  switch (method) {
    case 1: {
      b=(vals[0]-vals[1])/(coors[0]-coors[1]);
      c=vals[0]-b*coors[0];
      if (coors[0]==coors[1]) {
        message(logg,"\nzero distance in coordinates in extrapolate_coefs(1)\n",1+4);
        exit(1);
      }
      coefs[0]=c;coefs[1]=b;
      break;
    }
    case 2: {
      prec b1=(vals[0]-vals[1])/(coors[0]-coors[1]);
      prec b2=(vals[1]-vals[2])/(coors[1]-coors[2]);
      a=(b1-b2)/(coors[0]-coors[2]);
      b=b1-a*(coors[0]+coors[1]);
      c=vals[0]-a*coors[0]*coors[0]-b*coors[0];
      if (coors[0]==coors[1] || coors[0]==coors[2] ||coors[1]==coors[2]) {
        sprintf(buffer,"\nzero distance in coordinates in extrapolate_coefs(2):\n%e %e %e\n",coors[0],coors[1],coors[2]);
        splachni(logg,buffer,1+4);
        exit(1);
      }
      coefs[0]=c;coefs[1]=b;coefs[2]=a;
      break;
    }
    default: {
      message(logg,"wrong method in extrapolate_coefs\n",1+4);
      exit(1);
    }
  }
  #ifdef EXTRAPOLATE
  sprintf(buffer,"coeficients extrapolated by method %i\n",method);
  splachni(logg,buffer,4);
  for (int i=0;i<method+1;i++) {
    sprintf(buffer,"par[%i]=%e\tvalue[%i]=%e\tcoef[%i]=%e\n",i,coors[i],i,vals[i],i,coefs[i]);
    splachni(logg,buffer,4);
  }
  #endif
}


prec extrapolate_to_value(prec value, prec *coors, prec *vals, int method)
{

  prec coefs[3];
  extrapolate_coefs(coors,vals,method,coefs);
  prec a=coefs[2],b=coefs[1],c=coefs[0]-value;
  
  prec res=0;
  switch (method) {
    case 1 : {//linear
      if (b==0) res=0; else res=-c/b;
      break;
    }
    case 2 : {//quadratic
      if (a!=0) {
        prec d=b*b-4*a*c;
        if (d>0) {
          prec res1=(-b-sqrt(d))/(2*a);
          prec res2=(-b+sqrt(d))/(2*a);
          if (fabs(res1-coors[2])<fabs(res2-coors[2])) res=res1;
          else res=res2;
        }
        else res=0;
      }
      else {
        if (b==0) res=0; else res=-c/b;
      }
      break;
    }
    default: {
      message("wrong method in extrapolate_to_0!!!\n");
      exit(1);
    }
  }

#ifdef EXTRAPOLATE
  sprintf(buffer,"EXTRAP: by %i method a zero at %.4e computed\nfrom coefficients a=%.4e, b=%.4e, c=%.4e\nthose got from input data:\ncoordinates:",method, res,a,b,c);
  splachni(logg,buffer,4);
  for (int i=0;i<method+1;i++) {
    sprintf(buffer,"%.4e\t",coors[i]);
    splachni(logg,buffer,4);
  }
  message(logg,"\nvalues:\t",4);
  for (int i=0;i<method+1;i++) {
    sprintf(buffer,"%.4e\t",vals[i]);
    splachni(logg,buffer,4);
  }
  message(logg,"\n\n",4);
#endif
  return(res);
}
  
prec extrapolate_to_acr(prec *coors, prec *vals)
{
  prec d20=vals[0]-vals[2];
  prec p20=coors[0]-coors[2];
  prec d10=vals[1]-vals[2];
  prec p10=coors[1]-coors[2];
  
  if (p10*p20==0) {printf("zero distance in coordinates in extrapolate_to_acr\n");exit(1);}
  prec num=(d20*p10*p10-d10*p20*p20);
  prec den=(d20*p10-d10*p20);
  if (den==0) {printf("some shit in extrapolate_to_acr\nreturning zero\n"); return(0);}
  prec res=coors[2]+num/(den*2);

#ifdef EXTRAPOLATE
  sprintf(buffer,"EXTRAP: zero at %.4e (in %.4e) computed\nfrom input data:\ncoordinates:",res,num/(2*den));
  splachni(logg,buffer,4);
  for (int i=0;i<3;i++) {
    sprintf(buffer,"%.4e\t",coors[i]);
    splachni(logg,buffer,4);
  }
  message(logg,"\nvalues:\t",4);
  for (int i=0;i<3;i++) {
    sprintf(buffer,"%.4e\t",vals[i]);
    splachni(logg,buffer,4);
  }
  message(logg,"\n\n",4);
#endif
  return(res);
}

  
void predict(ret_from_arpack_zn& vysl, int eig, prec param,prec at_param, int method,prec* result, bool reset)
{
  static int pntr=0;
  static prec values[10][2];
  static prec params[2];

  if (reset) pntr=0;

  if (pntr<=method) {
    params[pntr]=param;
    for (int i=0;i<eig;i++) {
      values[i][pntr]=vysl.eigenvals[i].real();
      result[i]=0;
    }
    pntr++;
    return;
  }
  pntr--;
  for (int i=0;i<pntr;i++) params[i]=params[i+1];
  params[pntr]=param;
  for (int i=0;i<eig;i++) {
    for (int j=0;j<pntr;j++) {
      values[i][j]=values[i][j+1];
    }
    values[i][pntr]=vysl.eigenvals[i].real();
  }

  for (int i=0;i<eig;i++) {
    
    prec coefs[3];
    coefs[2]=0;
    extrapolate_coefs(values[i],params,method,coefs);
    prec a=coefs[2],b=coefs[1],c=coefs[0];
    result[i]=a*at_param*at_param+b*at_param+c;
  }
  
  pntr+=2;
  if (pntr>method+1) pntr--;
}


/*
#ifdef ENERGY_PERT
    level(-1,0,eig,eig,eig,0);
    for (int i=0;i<eig;i++) {
      int n,l,s;
      level(i+1,0,n,l,s,0);
      prec en=energy(n,l,s);
      if (ENERGY_PERT==1) {
        en+=correction(1,n,l,s);
        en+=correction(2,n,l,s);
      }
      vysl.eigenvals[i]=content(en,0);
    }
#endif //ENERGY_PERT
  
#ifdef FUNCTION_PERT
    level(-1,0,eig,eig,eig,0);
    content* aux=new content[r1.Nw*r1.sets];
    for (int i=0;i<eig;i++) {
      content *pntr=vysl.eigenvecs+r1.Nw*r1.sets*i;
      int n,l,s,spin=0;
      if (s==-1) spin=1;
      level(i+1,0,n,l,s,0);
      anal_function(r1,aux,n,l,s,d_0, 1, reg1::set);
      r1.operate(aux,pntr,0,0,reg1::lin_com,content(1,0),0,reg1::set);
      r1.operate(aux,pntr,1,1,reg1::lin_com,content(1,0),0,reg1::set);
      if (FUNCTION_PERT==1) {
        r1.operate(aux,pntr,1,0,reg1::lin_com,content(ab/2,ad/2),xc,reg1::add);
        r1.operate(aux,pntr,0,1,reg1::lin_com,content(-ab/2,ad/2),xc,reg1::add);
        r1.operate(aux,pntr,1,0,reg1::lin_com,content(-ad/2,-ab/2),yc,reg1::add);
        r1.operate(aux,pntr,0,1,reg1::lin_com,content(ad/2,-ab/2),yc,reg1::add);
      }
      norm_function(r1,pntr);
    }
#endif //FUNCTION_PERT

*/
void vec_dot_sigma(content* vec, content* Sigma)
{
  Sigma[0*2+0]=vec[2]+vec[3];
  Sigma[1*2+1]=-vec[2]+vec[3];  
  Sigma[0*2+1]=vec[0]-content(0,1)*vec[1];
  Sigma[1*2+0]=vec[0]+content(0,1)*vec[1];
} 

void vector2angles(double* in, double* out)
{
  for (int i=0;i<3;i++) out[i]=0;
  for (int i=0;i<3;i++) out[0]+=in[i]*in[i];
  out[0]=sqrt(out[0]);
  if (out[0]!=0) {
    out[1]=acos(in[2]/out[0]);
    //if (out[1]<0) out[1]+=2*M_PI;
    //if (out[1]>M_PI) out[1]-=M_PI;
    if (sin(out[1])!=0) {
      out[2]=acos(in[0]/(out[0]*sin(out[1])));
      //if (out[2]<0) out[2]*=-1;
      if (in[1]*sin(out[2])<0) out[2]=2*M_PI-out[2];
    }
  }
  prec res=in[2]-out[0]*cos(out[1]);
  res+=in[1]-out[0]*sin(out[1])*sin(out[2]);
  res+=in[0]-out[0]*sin(out[1])*cos(out[2]);
  if (abs(res)>1e-8) {
    printf("wrong vector2angle conversion:\n");
    printf("input: (%.6e, %.6e, %.6e) resulted in angles (%.6e, %.6e, %.6e) with residuum %e\n",in[0],in[1],in[2],out[0],out[1],out[2],res);
    exit(1);
  }
}

void spinor2angles(content* xi, double& theta, double& phi)
{
  if (xi[1]==content(0,0)) {theta=0;phi=0;return;}
  content ratio=xi[0]/xi[1];
  double r=abs(ratio);
  phi=-arg(ratio);
  if (phi<0) phi+=2*M_PI;
  theta=2*asin(sqrt(1/(1+r*r)));
  double res=abs(ratio-cos(theta/2)/sin(theta/2)*exp(content(0,-phi)));
  if (res>1e-8) {printf("wrong conversion in spinor2angles\n");exit(1);}
  //fprintf(logg,"input spinor: (%.2e%+.2ei,%.2e%+.2ei) results in angles (%.3e,%.3e) [pi]\n",xi[0].real(),xi[0].imag(),xi[1].real(),xi[1].imag(),theta,phi);
  return;
}

void rotate_inplane(content* kin, double angle, content* kout)
{
  static double angle_old=0;
  static double c=1,s=0;
  if (angle!=angle_old) {
    c=cos(angle);
    s=sin(angle);
    angle_old=angle;
  }
  kout[0]=c*kin[0]-s*kin[1];
  kout[1]=c*kin[1]+s*kin[0];
}

void rotate_inplane(double x_in, double y_in, double angle, double& x_out, double& y_out)
{
  static double angle_old=0;
  static double c=1,s=0;
  if (angle!=angle_old) {
    c=cos(angle);
    s=sin(angle);
    angle_old=angle;
  }
  x_out=c*x_in-s*y_in;
  y_out=c*y_in+s*x_in;
}


void waveguide2dot(prec xw, prec yw, prec& x, prec& y)
{
  prec wg_x=parameters.values[parameters_class::wg_x];
  prec wg_y=parameters.values[parameters_class::wg_y];
  prec wg_width=parameters.values[parameters_class::wg_width];
  prec wg_phi=parameters.values[parameters_class::wg_phi];
  
  xw+=wg_x;
  yw+=wg_y-wg_width/2;
  rotate_inplane(xw,yw,wg_phi,x,y);

  //x=xw*cos(wg_phi)-yw*sin(wg_phi);
  //y=yw*cos(wg_phi)+xw*sin(wg_phi);
}

bool dot2waveguide(prec x, prec y, prec& xw, prec& yw)
{
  prec wg_x=parameters.values[parameters_class::wg_x];
  prec wg_y=parameters.values[parameters_class::wg_y];
  prec wg_width=parameters.values[parameters_class::wg_width];
  prec wg_phi=parameters.values[parameters_class::wg_phi];

  rotate_inplane(x,y,-wg_phi,xw,yw);
  //xw=x*cos(-wg_phi)-y*sin(-wg_phi);
  //yw=y*cos(-wg_phi)+x*sin(-wg_phi);
  
  xw-=wg_x;
  yw-=wg_y-wg_width/2;
  if (xw>-wg_x && yw>0 && yw<wg_width) return(true);//inside waveguide means from middle of the dot towards the waveguide
  return(false);
}

void UnitaryTransformation(prec x, prec y, spinor& xi)
{
  prec ld=parameters.read("ldress","nm")*units.nm/units.length;
  prec lbr=parameters.read("lbr","nm")*units.nm/units.length;

  prec ldm1=0,lbrm1=0;
  if (ld!=0) ldm1=1/ld;
  if (lbr!=0) lbrm1=1/lbr;

  //spin-orbit vector
  prec n[3]={x*ldm1-y*lbrm1,x*lbrm1-y*ldm1,0};

  //unitary transformation
  prec n2=0;
  for (int i=0;i<3;i++) n2+=n[i]*n[i];
  prec vn2=abs(sqrt(n2));
  content aux[4]={0};
  aux[3]=cos(vn2/2);
  for (int i=0;i<3;i++) if (vn2>0) aux[i]=sin(vn2/2)*content(0,1)*n[i]/vn2; else aux[i]=content(0);
  content sigma[4];
  vec_dot_sigma(aux,sigma);
  xi=xi.SigmaTimesMe(sigma);
}


content psi_ground(prec x, prec y) {
  prec lB=parameters.read("lB","nm")*units.nm/units.length;
  prec res=exp(-(x*x+y*y)/(2*lB*lB));
  res/=sqrt(M_PI)*lB;
  return(content(res,0));
};


int psi_ground_spin,psi_ground_component,psi_ground_power;
content psi_ground_app(prec x, prec y) { 
  prec theta_B=parameters.read("theta_B","pi")*M_PI;
  prec phi_B=parameters.read("phi_B","pi")*M_PI;
  //prec ld=parameters.read("ldress","nm")*units.nm/units.length;
  //prec lbr=parameters.read("lbr","nm")*units.nm/units.length;

  //prec ldm1=0,lbrm1=0;
  //if (ld!=0) ldm1=1/ld;
  //if (lbr!=0) lbrm1=1/lbr;

  //magnetic field vector
  //prec b[3]={sin(theta_B)*cos(phi_B),sin(theta_B)*sin(phi_B),cos(theta_B)};

  //spin-orbit vector
  //prec n[3]={x*ldm1-y*lbrm1,x*lbrm1-y*ldm1,0};

  //spin-orbit magnetic field
  //prec b_so[3]={n[1]*b[2]-n[2]*b[1],n[2]*b[0]-n[0]*b[2],n[0]*b[1]-n[1]*b[0]};

  //total field
  //prec b_tot[3];
  //for (int i=0;i<3;i++) b_tot[i]=b[i]+b_so[i]*0;

  //locally aligned spinor
  prec b_tot_angles[3]={1,theta_B,phi_B};
  //vector2angles(b_tot,b_tot_angles);
  spinor xi;
  if (psi_ground_spin==1) {
    xi.up=cos(b_tot_angles[1]/2);
    xi.down=sin(b_tot_angles[1]/2)*exp(content(0,b_tot_angles[2]));
  }
  else {
    xi.up=-sin(b_tot_angles[1]/2);
    xi.down=cos(b_tot_angles[1]/2)*exp(content(0,b_tot_angles[2]));
  }

  //unitary transformation
  UnitaryTransformation(x,y,xi);
  
  //final
  content res=psi_ground(x,y);
  if (psi_ground_component==0) res*=xi.up;
  if (psi_ground_component==1) res*=xi.down;
  return(pow(res,psi_ground_power));
};



#include "main.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_eigen.h"

#define TRANS1
//#define TRANS2
//#define TRANS3
//#define TRANSCHECK

//to do once - the geometry dependent (energy independent) tasks
void transport_class::prepare_geometry(bool full)
{
  method_full=full;

  //regions();
  //!!!before prepare_regions() was called here, now I have split it to have spin impurities possible
  //region->ant(potential,units.length/units.nm, units.energy/units.meV);
  
  int i,j;
  region->sec2gridlimits(0,i,j);
  xlines=j-i+1;
  if (xlines<1) {printf("prepare geometry: no points inside the conductor (imin=%i, imax=%i)\n",i,j);exit(1);}
#ifdef TRANS2
  fprintf(logg,"conductor within grid range: (%i-%i): total %i xlines\n",i,j,xlines);
#endif

  dim=new int[xlines];
  for (int si=0;si<xlines;si++) {
    dim[si]=region->sec2dim(0,si);
#ifdef TRANS2
    fprintf(logg,"conductor (sector 0) line %i contains %i points\n",si,dim[si]);
#endif
  }

  message(logg,"allocating space...",1+4);
  allocate_g_function();
  
  self_array=new content*[xlines*xlines];
  self_filled=new bool[xlines*xlines];
  for (int i=0;i<xlines*xlines;i++) self_filled[i]=false;

  for (int l=1;l<Nl;l++) {
    int lmin,lmax,imin,imax;
    region->sec2gridlimits(l,lmin,lmax);
    region->sec2gridlimits(0,imin,imax);
    int ifrom=max(lmin-b,imin)-imin;
    int ito=min(lmax+b,imax)-imin;
    for (int i=ifrom;i<=ito;i++) {
      for (int iy=0;iy<dim[i];iy++) {
        if (!region->secnearbysec(0,i,iy,l)) continue;
        for (int j=ifrom;j<=ito;j++) {
          for (int jy=0;jy<dim[j];jy++) {
            if (!region->secnearbysec(0,j,jy,l)) continue;
            self_filled[i*xlines+j]=true;
            break;
          }
        }
      }
    }
  }
  int tot=0;
  for (int i=0;i<xlines;i++) for (int j=0;j<xlines;j++) if (self_filled[i*xlines+j]) tot+=dim[i]*dim[j]*4;
  self_array_base=new content[tot];
#ifdef TRANS1
    fprintf(logg,"the following self energy submatrices are allocated:");
#endif
  int pntr=0;
  for (int i=0;i<xlines;i++) {
    for (int j=0;j<xlines;j++) {
      if (self_filled[i*xlines+j]) {
        self_filled[i*xlines+j]=false;
        self_array[labelM(i,j)]=self_array_base+pntr;
        pntr+=dim[i]*dim[j]*4;
#ifdef TRANS1
        fprintf(logg,"%i-%i ",i, j);
#endif
      }
    }
  }

  int mx=0;
  for (int i=0;i<xlines;i++) if (mx<dim[i]) mx=dim[i];
  auxM=new content[mx*mx*4];

  s_matrix.s=new content*[Nl*Nl];
  s_matrix.M=new int[Nl];
  
  message(logg,"...done\n",1+4);

  message(logg,"building the transversal-only-dependent Hamiltonian in the leads and computing lead transverse modes...",1+4);
  leads_geometry_init();
  message(logg,"...done\n",1+4);
  
  message(logg,"building the Hamiltonian including leads...",1+4);
  hamiltonian_init(*region);

  //no spin-orbit terms, no orbital effects of the magnetic field
  prec br=parameters.read("br","meVA");
  prec dress=parameters.read("dress","meVA");
  prec dress3=parameters.read("dress3","eVA^3");
  bool peierls=(bool) parameters.read("peierls","-");
  parameters.set("br",0,"meVA");
  parameters.set("dress",0,"meVA");
  parameters.set("dress3",0,"eVA^3");
  parameters.set("peierls",0,"-");
  hamiltonian_init(*region2);
  parameters.set("br",br,"meVA");
  parameters.set("dress",dress,"meVA");
  parameters.set("dress3",dress3,"eVA^3");
  parameters.set("peierls",peierls,"-");

  message(logg,"...done\n",1+4);
}

//geometry dependent tasks for leads
void transport_class::leads_geometry_init()
{
  
  for (int l=1;l<Nl;l++) {
    hamiltonian_for_lead_init(*region, lead[l].heading);
    int i0=0,iy0=0,iadd=0,iyadd=0,imin=0,imax=0,N=0;
    region->sec2gridlimits(l,imin,imax);
    switch (lead[l].heading) {
      case 1 : {i0=imax-imin;iy0=0;iadd=0;iyadd=1;N=region->sec2dim(l,i0);break;}
      case 2 : {i0=0;iy0=0;iadd=0;iyadd=1;N=region->sec2dim(l,i0);break;}
      case 3 : {i0=0;iy0=region->sec2dim(l,i0)-1;iadd=1;iyadd=0;N=imax-imin+1;break;}
      case 4 : {i0=0;iy0=0;iadd=1;iyadd=0;N=imax-imin+1;break;}
    }
    //fprintf(logg,"building H for lead %i: i0=%i, iy0=%i, iadd=%i, iyadd=%i, N=%i\n",l,i0,iy0,iadd,iyadd,N);
    lead[l].points=N;
    content* ham=new content [N*2*N*2];
    for (int p1=0;p1<N;p1++) {
      for (int s1=0;s1<2;s1++) {
        for (int p2=0;p2<N;p2++) {
          for (int s2=0;s2<2;s2++) {
            ham[(p1+s1*N)*N*2+(p2+s2*N)]=region->Hsecsec(l,i0+p1*iadd,iy0+p1*iyadd,s1,l,i0+p2*iadd,iy0+p2*iyadd,s2);
            //fprintf(logg,"H entry: p1=%i, s1=%i, p2=%i, s2=%i: H=(%e,%e)\n",p1,s1,p2,s2,ham[(p1+s1*N)*N*2+(p2+s2*N)].real(),ham[(p1+s1*N)*N*2+(p2+s2*N)].imag());
          }
        }
      }
    }

    lead[l].mode=new mode_prototype[N*2];
    content* res=new content[N*2*N*2];
    prec* en=new prec[N*2];
#ifdef TRANS2 
    fprintf(logg,"transversal hamiltonian [nat. units * hy^2]:\n");
    for (int i=0;i<N*2;i++) {
      for (int j=0;j<N*2;j++) {
        fprintf(logg,"%.2e%+.2ei ",ham[i*2*N+j].real()*hy*hy,ham[i*2*N+j].imag()*hy*hy);
      }
      fprintf(logg,"\n");
    }
#endif
    diag(ham,res,en,N*2);
#ifdef TRANS2
    content aux=0;
    for (int i=0;i<2*N;i++) {
      for (int j=0;j<2*N;j++) {
        for (int k=0;k<2*N;k++) aux+=ham[i*2*N+k]*res[k+2*N*j];
        aux-=en[i]*res[i*2*N+j];
      }
    }
    fprintf(logg,"transversal wavefunctions are eigensystem with precision %e\n", abs(aux));
#endif

    for (int m=0;m<N*2;m++) {
      lead[l].mode[m].wavefunction=res+N*2*m;
      lead[l].mode[m].energy=en[m]-en[0]*0;
      //identify the spinor of this mode
      prec mx=0;
      spinor largest, act;
      for (int i=0;i<N;i++) {
        act.up=lead[l].mode[m].wavefunction[i];
        act.down=lead[l].mode[m].wavefunction[N+i];
        prec ac=(act.braket(act)).real();
        if (ac>mx) {mx=ac; largest=act;}
      }
      largest.up/=sqrt(mx);
      largest.down/=sqrt(mx);
      lead[l].mode[m].normalized_spinor=largest;

#ifdef TRANS2
      content aux=0;
      for (int y=0;y<2*N;y++) {
        aux+=en[m]*res[m*2*N+y];
        for (int yp=0;yp<2*N;yp++) aux-=ham[y*2*N+yp]*res[m*2*N+yp];//lead[l].mode[m].wavefunction[yp];
      }
      fprintf(logg,"lead %i mode %i wavefunction (energy:%e [meV]) is eigenstate with precision %e; values:\n\t",l,m,lead[l].mode[m].energy*units.energy/units.meV, abs(aux));
      for (int y=0;y<N*2;y++) fprintf(logg,"%.2e%+.2ei ",lead[l].mode[m].wavefunction[y].real(), lead[l].mode[m].wavefunction[y].imag());
      fprintf(logg,"\n");
      prec theta, phi;
      lead[l].mode[m].normalized_spinor.angles(theta, phi);
      fprintf(logg,"the spinor part of the eigenmode is: (%e%+ei, %e%+ei) [quantization axis along angles (%e, %e)]\n",lead[l].mode[m].normalized_spinor.up.real(), lead[l].mode[m].normalized_spinor.up.imag(), lead[l].mode[m].normalized_spinor.down.real(), lead[l].mode[m].normalized_spinor.down.imag(),theta, phi);
#endif
    }
    delete en;
    delete ham;
  }
}

//construct the region including leads
void transport_class::prepare_regions()
{
  message(logg,"parameters recomputed and listed in transport_class::prepare_regions()\n",1+4);
  parameters.recompute();
  parameters.list(logg,1+4);
  message(logg,"creating regions with leads...",1+4);

  b=(int) parameters.read("precision","-")+1;
  Nl=(int) parameters.read("leads","-")+1;
  int dim_x=(int) parameters.read("dim_x","-");
  int dim_y=(int) parameters.read("dim_y","-");
  prec length_x=parameters.read("length_x","nm")*units.nm/units.length;
  prec length_y=parameters.read("length_y","nm")*units.nm/units.length;
  prec dLx=parameters.read("dLx","nm")*units.nm/units.length;
  prec dLy=parameters.read("dLy","nm")*units.nm/units.length;
  int geometry=(int) parameters.read("geometry","-");
  int leads=(int) parameters.read("leads","-");
  int sets=(int) parameters.read("sets","-");
  if (leads<1) {printf("no active leads (#%i) in transport\n",leads);exit(1);}
  if (geometry!=1) {
    //2D potential
    region=new reg1(dim_x,dLx-length_x,2*length_x,'d',dim_y,dLy-length_y,2*length_y,'d',inside,sets);
    region2=new reg1(dim_x,dLx-length_x,2*length_x,'d',dim_y,dLy-length_y,2*length_y,'d',inside,sets);
  }
    else {
    //1D 
    region=new reg1(dim_x,dLx-length_x,2*length_x,'d',dim_y,dLy-length_y,2*length_y,'d',inside,sets);
    region2=new reg1(dim_x,dLx-length_x,2*length_x,'d',dim_y,dLy-length_y,2*length_y,'d',inside,sets);
  }
  region->give_par(hx,hy);
  //region->ant(potential,units.length/units.nm, units.energy/units.meV);
  region->sectorize(inside_sector,Nl);
  //region->ant(potential,units.length/units.nm, units.energy/units.meV);
  region2->sectorize(inside_sector,Nl);
  //region->ant(potential,units.length/units.nm, units.energy/units.meV);

  //!!!!
  int i,j;
  region->sec2gridlimits(0,i,j);
  xlines=j-i+1;
  //for each lead: 
  //1. check the lead is rectangular
  //2. check it has at least b points in the heading direction
  for (int l=1;l<Nl;l++) {
    int imin,imax;
    region->sec2gridlimits(l,imin,imax);
    int idim0=region->sec2dim(l,0);
    for (int i=imin+1;i<=imax;i++) {
      if (idim0!=region->sec2dim(l,i-imin)) {
        printf("lead %i is not rectangular\n",l);
        printf("\t line 0: %i points, line %i: %i points\n",idim0,i,region->sec2dim(l,i-imin));
        exit(1);
      }
    }
    int length=0;
    switch (lead[l].heading) {
      case 1 : ;
      case 2 : {length=imax-imin+1;break;}
      case 3 : ;
      case 4 : {length=idim0;break;}
    }
    if (length<b) {
      printf("not enough points (%i vs b=%i) in lead %i\n", length, b, l);
      exit(1);
    }
#ifdef TRANS1
    fprintf(logg, "lead %i heading:%i grid xlines:%i-%i, number of points on each %i\n",l, lead[l].heading, imin, imax, idim0);
#endif
  }
#ifdef TRANS1
  fprintf(logg, "grid dimensions: hx=%.4e hy=%.4e [nm]; %i xlines inside conductor\n",hx*units.length/units.nm, hy*units.length/units.nm, xlines);
#endif
}

void transport_class::allocate_g_function()
{
  if (method_full) {
    g_array=new content*[2*xlines*xlines];
    int tot=0;
    for (int i=0;i<xlines;i++) for (int j=0;j<xlines;j++) tot+=dim[i]*dim[j]*4;
    g_array_base=new content[2*tot];  

    int pntr=0;
    for (int i=0;i<xlines;i++) {
      for (int j=0;j<xlines;j++) {
        g_array[labelg(0,i,j)]=g_array_base+pntr;
        g_array[labelg(1,i,j)]=g_array_base+pntr+tot;
        pntr+=dim[i]*dim[j]*4;
      }
    }
  }
  else {
    g_array=new content*[(b+1)*(b+1)*xlines];
    int tot=0;
    for (int up=0;up<xlines;up++) {
      for (int i=0;i<=up;i++) {
        for (int j=0;j<=up;j++) {
          tot+=dim[i]*dim[j]*4;
          cycle(up,j);
        }
        cycle(up,i);
      }
    }
    g_array_base=new content[tot];  

    int pntr=0;
    for (int up=0;up<xlines;up++) {
      for (int i=0;i<=up;i++) {
        for (int j=0;j<=up;j++) {
#ifdef TRANS3
          fprintf(logg,"g_array allocated: up=%i, i=%i, j=%i\n",up,i,j);
#endif
          g_array[labelg(up,i,j)]=g_array_base+pntr;
          pntr+=dim[i]*dim[j]*4;
          cycle(up,j);
        }
        cycle(up,i);
      }
    }
  }
}
  

//deallocate space from geometry dependent tasks
void transport_class::clean_geometry()
 {
  for (int l=1;l<Nl;l++) {
    delete lead[l].mode[0].wavefunction;
    delete lead[l].mode;
  }
  delete dim;

  delete self_array_base;
  delete self_array;
  delete g_array_base;
  delete g_array;

  delete self_filled;

  delete region;
  delete region2;

  delete auxM;

  delete s_matrix.s;
  delete s_matrix.M;
}

//the energy dependent tasks
void transport_class::prepare_with_E()
{
  for (int i=0;i<xlines;i++) for (int j=0;j<xlines;j++) self_filled[labelM(i,j)]=false;

  message(logg,"computing lead green function...\n",1+4);
  leads_with_E_init();
  message(logg,"...done\n",1+4);

  message(logg,"computing self energy...",1+4);
  compute_self();
  message(logg,"...done\n",1+4);
}

//energy dependent tasks for leads
void transport_class::leads_with_E_init()
{
  prec Ef=parameters.read("fermi_e","meV")*units.meV/units.energy;
  prec m_eff=parameters.read("m_eff","m_e");
  for (int l=1;l<Nl;l++) {
    lead[l].propagating_modes=lead[l].points*2;
    lead[l].orbital_modes=0;
    for (int m=0;m<lead[l].points*2;m++) {
      prec e=Ef-lead[l].potential-lead[l].mode[m].energy;
#ifdef TRANS2
      fprintf(logg,"lead %i mode %i, with energy %e (energy E0 %e, potential %e, fermi Ef %e) [meV]\n", l, m, e*units.energy/units.meV, lead[l].mode[m].energy*units.energy/units.meV, lead[l].potential*units.energy/units.meV, Ef*units.energy/units.meV);
#endif
      if (e<0) {
        lead[l].propagating_modes=m;
        break;
      }
      lead[l].mode[m].orb=-1;
      lead[l].mode[m].spin=-1;
      lead[l].mode[m].wavevector=sqrt(2*units.m_e*m_eff*e*units.energy)/units.hbar/units.wavevector;
      lead[l].mode[m].velocity=region->wave_velocity(lead[l].mode[m].wavevector)*units.hbar/units.length/(units.m_e*m_eff)/units.velocity;
      lead[l].mode[m].longitudinal=new content[b*b];
      for (int j=0;j<b;j++) {
        for (int k=0;k<b;k++) {
          lead[l].mode[m].longitudinal[j*b+k]=leadg_long(j+1,k+1,lead[l].mode[m].wavevector,lead[l].mode[m].velocity);
        }
      }
#ifdef TRANS2
      fprintf(logg,"\twave vector %e [1/nm], velocity %e [m/s]\n",lead[l].mode[m].wavevector*units.wavevector*units.nm, lead[l].mode[m].velocity*units.velocity);
      fprintf(logg,"longitudinal greens function:\n\t\t");
      for (int j=0;j<b;j++) {
        for (int k=0;k<b;k++) {
          fprintf(logg,"(%i-%i): %.2e%+.2ei ",j,k, lead[l].mode[m].longitudinal[j*b+k].real(), lead[l].mode[m].longitudinal[j*b+k].imag());
        }
      }
      fprintf(logg,"\n");
#endif
    }
#ifdef TRANS1
    fprintf(logg,"total %i propagating modes in lead %i with energies:\n",lead[l].propagating_modes,l);
    for (int m=0;m<lead[l].propagating_modes;m++) fprintf(logg, "%e ",lead[l].mode[m].energy*units.energy/units.meV);
    fprintf(logg,"\n");
    if (lead[l].propagating_modes<lead[l].points*2) {
      fprintf(logg,"lowest pinched-off mode energy:%e vs Ef=%e minus lead offset=%e \n",lead[l].mode[lead[l].propagating_modes].energy*units.energy/units.meV, 
Ef*units.energy/units.meV, lead[l].potential*units.energy/units.meV);
    }
#endif
  }

  pair_modes();

  for (int i=1;i<Nl;i++) {
    for (int j=1;j<Nl;j++) {
      s_matrix.s[i*Nl+j]=new content[lead[i].propagating_modes*lead[j].propagating_modes];
    }
    s_matrix.M[i]=lead[i].propagating_modes;
  }
}

//computes self energy
void transport_class::compute_self()
{
  //Sigma_{I J} = sum_{leads} sum_{IP} H_{I IP} sum_{JP} G(lead)_{IP JP} H_{JP J}
  //where multiindexes are: I=(i yi si), J=(j yj sj), Ip=(ip yip sip), Jp=(jp yjp sjp)
  //and xlines (i,j,ip,jp) values refer to the no-lead region regiono 
  //I and J are inside the region, IP and JP are inside the lead 
  //all labels correspond to region with leads labels!!!

  for (int l=1;l<Nl;l++) {

    int lmin,lmax,imin,imax;
    region->sec2gridlimits(l,lmin,lmax);
    region->sec2gridlimits(0,imin,imax);
    int ifrom=max(lmin-b,imin)-imin;
    int ito=min(lmax+b,imax)-imin;
    int modes=lead[l].propagating_modes;
#ifdef TRANS2
    fprintf(logg,"self energy for lead %i, heading %i with %i propagating modes, sector 0 xlines:%i-%i\n",l,lead[l].heading,modes,ifrom,ito);
#endif

    //sum through I and J
    for (int i=ifrom;i<=ito;i++) {
      for (int iy=0;iy<dim[i];iy++) {
        if (!region->secnearbysec(0,i,iy,l)) continue;
        for (int j=ifrom;j<=ito;j++) {
          for (int jy=0;jy<dim[j];jy++) {
            if (!region->secnearbysec(0,j,jy,l)) continue;

//fprintf(logg,"point inside i=%i,j=%i, yi=%i, yj=%i\n",i,j,yi,yj);
            //sum through IP and the H_{I,IP} value 
            int lfromip=max(i+imin-b,lmin)-lmin;
            int ltoip=min(i+imin+b,lmax)-lmin;
            int lfromjp=max(j+imin-b,lmin)-lmin;
            int ltojp=min(j+imin+b,lmax)-lmin;
            for (int is=0;is<2;is++) {
              for (int js=0;js<2;js++) {

                content SIJ=0;
                for (int ip=lfromip;ip<=ltoip;ip++) {
                  int ipdim=region->grid2secdim(l,ip+lmin);
                  for (int ipy=0;ipy<ipdim;ipy++) {
    //fprintf(logg,"point outside ip=%i, yip=%i\n",ip,yip);
                    //printf("calling iszeroH from self with: i=%i, ip=%i, yi=%i, yip=%i\n",i,ip,yi,yip);
                    if (iszeroH(0,i,iy,l,ip,ipy)) continue;
    //fprintf(logg,"coupled point outside ip=%i, yip=%i\n",ip,yip);
                    for (int ips=0;ips<2;ips++) {
                      content HIIP=region2->Hsecsec(0,i,iy,is,l,ip,ipy,ips);
  
                      //sum through JP and the H_{JP,J} value
                      content GHIPJ=0;
                      for (int jp=lfromjp;jp<=ltojp;jp++) {
                        int jpdim=region->grid2secdim(l,jp+lmin);
                        for (int jpy=0;jpy<jpdim;jpy++) {
                          //fprintf(logg,"checking JP jp=%i, yjp=%i\n",jp,yjp);
                          //fprintf(logg,"is in\n");
                          if (iszeroH(l,jp,jpy,0,j,jy)) continue;
                          //fprintf(logg,"is coupled\n");
                          for (int jps=0;jps<2;jps++) {
                              content HJPJ=region2->Hsecsec(l,jp,jpy,jps,0,j,jy,js);

                            //the lead greens function - sum through the modes
                            //li/lj: the perpendicular distance from the boudary minus one
                            //wi/wj: the transversal wavefunction label (longitudinal distance)
                            int li=0,lj=0,wi=0,wj=0;
                            switch (lead[l].heading) {
                              case 1 : {li=lmax-(ip+lmin);lj=lmax-(jp+lmin);wi=ipy;wj=jpy;break;}
                              case 2 : {li=ip;lj=jp;wi=ipy;wj=jpy;break;}
                              case 3 : {li=ipdim-ipy-1; lj=jpdim-jpy-1; wi=ip;wj=jp;break;}
                              case 4 : {li=ipy; lj=jpy; wi=ip;wj=jp;break;}
                            }
                            //fprintf(logg,"through %i modes\n",modes);

                            content GIPJP=0;
                            for (int m=0;m<modes;m++) {
                              content aux=lead[l].mode[m].longitudinal[li*b+lj]*lead[l].mode[m].wavefunction[wi+ips*lead[l].points]*conj(lead[l].mode[m].wavefunction[wj+jps*lead[l].points]);
                              GIPJP+=aux;
#ifdef TRANS3
    if (abs(HIIP*aux*HJPJ)>1e-12) {
      fprintf(logg,"self energy entry:\n");
      fprintf(logg,"H(i=%i, iy=%i, is=%i, ip=%i, ipy=%i, ips=%i)=%e%+ei\n",i, iy, is, ip, ipy, ips, HIIP.real(), HIIP.imag());
      fprintf(logg,"G(IP, JP) from lead %i mode %i li=%i, lj=%i, wi=%i, wj=%i, long=%.2e%+.2ei, |chi(IP)|=%.2e |chi(JP)|=%e\n",l, m, li, lj, wi, wj, lead[l].mode[m].longitudinal[li*b+lj].real(), lead[l].mode[m].longitudinal[li*b+lj].imag(), abs(lead[l].mode[m].wavefunction[wi+ips*lead[l].points]), abs(lead[l].mode[m].wavefunction[wj+jps*lead[l].points]));
      fprintf(logg,"H(jp=%i, jpy=%i, jps=%i, j=%i, jy=%i, js=%i)=%e%+ei\n",jp, jpy, jps, j, jy, js, HJPJ.real(), HJPJ.imag());
    }
#endif
                            }
                            GHIPJ+=GIPJP*HJPJ;
                          }
                        }
                      }
                      SIJ+=HIIP*GHIPJ;
                    }
                  }
                }
                if (!self_filled[labelM(i,j)]) self_filled[labelM(i,j)]=true;
                self_array[labelM(i,j)][label(i,j,iy,jy,is,js)]=SIJ;//*hx*hy;
                //hy is for the transversal wavefunctions - twice square root of hy is compensated
                //by the normalization of the Greens function on 1 instead of delta(y-y') similarly as in glonditudinal
                //the dimension of G is 1/energy
              }
            }
          }
        }
      }
    }
  }
}

//deallocate space from energy dependent tasks
void transport_class::clean_with_E()
{

  for (int l=1;l<Nl;l++) {
    for (int m=0;m<lead[l].propagating_modes;m++) {
      delete lead[l].mode[m].longitudinal;
    }
  }
  for (int i=1;i<Nl;i++) for (int j=1;j<Nl;j++) delete s_matrix.s[i*Nl+j];
}

//compute the green function between the first and the last xline using recursive method
//full -- compute the full inverse matrix (true) or just the four corners (false)
void transport_class::recursive_green(bool full,  progress_bar_prototype* progress_bar)
{
  if (full!=method_full) {//the method has changed - change the allocation space for the g-function
    method_full=full;
    delete g_array;
    delete g_array_base;
    allocate_g_function();
  }

  bool with_progress_bar=true;
  if (progress_bar==0) with_progress_bar=false;
  sprintf(buffer,"starting the recursive method(full:%i)\nstep 0 ",method_full);
  splachni(logg,buffer,1+4);
  prec Ef=parameters.read("fermi_e","meV")*units.meV/units.energy;
  

  //content* aux=new content[n0*n0];
  
  fill_EmHmS(auxM,0,0,Ef);
  //sprintf(buffer,"\tmatrix EmHmS(0,0) loaded (norm:%e)\n",normM(aux,0,0));splachni(logg,buffer,4);  
  invert(g(0,0,0),auxM,dim[0]*2);
  //g_filled[labelM(0,0)]=true;
  //check_invert(g(0,0),aux,0);
  //delete aux;
  //sprintf(buffer,"\tmatrix g(0,0) computed (norm:%e). Starting the loop\n",normM(g(0,0),0,0));splachni(logg,buffer,4);
  //check_it(0);
  //print_matrixes();
  if (with_progress_bar) progress_bar->add(xlines,1);

  for (int i=1;i<xlines;i++) {
    //if (i%(max(xlines/10,1))==0) {sprintf(buffer,"%i ",i);splachni(logg,buffer,1+4);}
    //content* aux=new content[ni*ni];
    //printf("round: i=%i\n",i);
    //last down corner: g(i,i)
    fill_EmHmS(auxM,i,i,Ef);
    VGV(auxM,i,i,true,-1);
    invert(g(i,i,i),auxM,dim[i]*2);
    //g_filled[labelM(i,i)]=true;
    //check_invert(g(i,i),aux,i);
    //delete aux;
    //sprintf(buffer,"\tg(%i,%i) computed (norm:%e)\n",i,i,normM(g(i,i),i,i));splachni(logg,buffer,4);

    //last row and last column: g(i,j) and g(j,i)
    for (int j=0;j<i;j++) {
      //int nj=dim[j]*2;
      //content* aux=new content[ni*nj];
      VG(auxM,i,j,false,1);
      //sprintf(buffer,"\t\t(norm:%e) auxilary: VG (%i,%i) computed\n",normM(aux,i,j),i,j);splachni(logg,buffer,1+4);
      MM(g(i,i,j),g(i,i,i),auxM,i,i,i,j);
      //g_filled[labelM(i,j)]=true;
      //sprintf(buffer,"\tg(%i,%i) computed (norm:%e)\n",i,j,normM(g(i,j),i,j));splachni(logg,buffer,4);
      GV(auxM,j,i,false,1);
      //sprintf(buffer,"\t\t(norm:%e) auxilary: GV (%i,%i) computed\n",normM(aux,j,i),j,i);splachni(logg,buffer,1+4);
      MM(g(i,j,i),auxM,g(i,i,i),j,i,i,i);
      //g_filled[labelM(j,i)]=true;
      //delete aux;
      //sprintf(buffer,"\tg(%i,%i) computed (norm:%e)\n",j,i,normM(g(j,i),j,i));splachni(logg,buffer,4);
      if (!method_full) cycle(i,j);
    }

    //inside: g(j,k)
    for (int j=0;j<i;j++) {
      //int nj=dim[j]*2;
      //content* aux=new content[nj*ni];
      GV(auxM,j,i,false,1);
      for (int k=0;k<i;k++) {
        //int nk=dim[k]*2;
        MM(g(i,j,k),auxM,g(i,i,k),j,i,i,k);
        add(g(i,j,k),g(i-1,j,k),j,k,1);
        //g_aux_filled[labelM(j,k)]=true;
        //sprintf(buffer,"\t(norm:%e) g_aux(%i,%i) computed\n",normM(g_aux(j,k),j,k),j,k);splachni(logg,buffer,1+4);
        if (!method_full) cycle(i,k);
      }
      //delete aux;
      if (!method_full) cycle(i,j);
    }
  
    //swap g_aux and g
    /*for (int j=0;j<i;j++) {
      for (int k=0;k<i;k++) {
        if (!g_filled[labelM(j,k)]) {printf("g(j,k) empty???\n");exit(1);}
        content* auxp=g_array[labelM(j,k)];
        g_array[labelM(j,k)]=g_aux(j,k);
        g_aux_array[labelM(j,k)]=auxp;
        g_aux_filled[labelM(j,k)]=false;
        //sprintf(buffer,"\tg_aux(%i,%i) copied to g(%i,%i) (norm:%e) \n",j,k,j,k,normM(g(j,k),j,k));splachni(logg,buffer,4);
        if (!method_full) cycle(i,k);
      }
      if (!method_full) cycle(i,j);
    }*/
#ifdef TRANSCHECK
    if ((i==xlines-1) && method_full) {fprintf(logg,"checking at step i=%i\n",i);check_it(i);}
#endif
    if (with_progress_bar) progress_bar->add(0,1);

  }
  message(logg,"\n",1+4);
  //print_matrixes();
}

//returns the longitudinal part of the green function of a semiinfinite lead
content transport_class::leadg_long(int i, int j, prec k, prec v)
{
  prec ka=k*hx;
  content res=exp(content(0,ka*abs(i-j)))-exp(content(0,ka*abs(i+j)));
  res*=content(0,-1)/(v*units.velocity*units.hbar); //amplitude of the 1d greens function
  //printf("lead_log called with i=%i, j=%i, ka=%f, result=%f%+fi\n",i,j,ka,res.real(),res.imag());
  res*=units.length*hx;
                            //in the algorithm, glead is an inverse of a lead hamiltonian matrix, H(i,j).G(j,k)=1_ik
                            //whereas the upper is a discretized greens function, defined H(x')g(x,x')=delta(x-x')
  res*=units.energy;        //the result is in units of the inverse energy
  return(res);
}

//initialize the lead
void transport_class::set_lead(int l, int heading, double x, double y, double width, double V)
{
  int leads=(int) parameters.read("leads","-");
  if (l<1 || l>min(leads,4)) {
    printf("initializing a lead (#%i) out of range (0-%i)\n",l,min(leads,4));
    exit(1);
  }
  if (heading<1 || heading>4) {
    printf("initializing a lead (#%i): heading (%i) out of range (0-3)\n",l,heading);
    exit(1);
  }
  lead[l].heading=heading;
  lead[l].x=x*units.nm/units.length;
  lead[l].y=y*units.nm/units.length;
  lead[l].width=width*units.nm/units.length;
  lead[l].potential=V*units.meV/units.energy;
}

void transport_class::cycle(int i, int &j) 
{
  //if (i==xlines-1) j=i-1;
  if (j==0) j=max(0,i-b);
}


void transport_class::amplitude_aux1(int l, int& ifrom, int &ito) 
{
  int lmin,lmax,imin,imax;
  region->sec2gridlimits(l,lmin,lmax);
  region->sec2gridlimits(0,imin,imax);
  switch (lead[l].heading) {
    case 1 : {ifrom=ito=lmax+1;break;}
    case 2 : {ifrom=ito=lmin-1;break;}
    case 3 : ;
    case 4 : {ifrom=lmin;ito=lmax;break;}
  }
  ifrom-=imin;
  ito-=imin;
}

void transport_class::amplitude_aux2(int l, int i, int& yifrom, int &yito) 
{
  switch (lead[l].heading) {
    case 1 : ;
    case 2 : {yifrom=0;yito=dim[i]-1;break;}
    case 3 : {yifrom=yito=0;break;}
    case 4 : {yifrom=yito=dim[i]-1;break;}
  }
}

//returns transversal label of a point in the lead neighbouring the input grid point
//if no, returns -1 - assumes no holes in the lead
int transport_class::jcor(int l, int si, int siy)
{
  int gi,gj;
  region->sec2grid(0,si,siy,gi,gj);

  //printf("neighbouring point in lead:\ninput: lead=%i, xline i=%i, yi=%i (regiono)\n",l,i,yi);
  switch (lead[l].heading) {
    case 1 : {gi--;break;}
    case 2 : {gi++;break;}
    case 3 : {gj--;break;}
    case 4 : {gj++;break;}
  }
  int i,j,sec;
  region->grid2sec(gi,gj, sec, i,j);
  if (sec!=l) return(-1);
  switch (lead[l].heading) {
    case 1 : ;
    case 2 : return(j);
    case 3 : ;
    case 4 : return(i);
  }
  return(j);
}


//tranmission amplitude from lead, mode into a second combination of these two
content transport_class::amplitude(int lead_from, int mode_from, int lead_to, int mode_to)
{

#ifdef TRANS3
  fprintf(logg,"computing amplitude: from lead=%i mode=%i to lead=%i mode=%i\n", lead_from, mode_from, lead_to, mode_to);
#endif
  
  prec Ef=parameters.read("fermi_e","meV")*units.meV/units.energy;

  if (lead_from<1 || lead_from>=Nl || lead_to<1 || lead_to>=Nl ) {
    printf("invalid lead (%i/%i) in amplitude!\n",lead_from,lead_to); exit(1);}

  if (lead[lead_from].propagating_modes<=mode_from) {
    //fprintf(logg,"incoming lead mode not propagating, returning transmission zero in amplitude\n");
    return(0);
  }
  if (lead[lead_to].propagating_modes<=mode_to) {
    //fprintf(logg,"outgoing lead mode not propagating, returning transmission zero in amplitude\n");
    return(0);
  }

  //sprintf(buffer,"transmission:\nusing green function (%i,%i) (norm:%e) \n",i,j,normM(g(i,j),i,j));splachni(logg,buffer,1+4);
  
  content res=0;
#ifdef TRANS3
  int aux1, aux2;
  amplitude_aux1(lead_from, aux1, aux2);
  fprintf(logg,"amplitude loop i: bounds (%i, %i)\n", aux1, aux2);
  amplitude_aux2(lead_from, aux1, aux1, aux2);
  fprintf(logg,"\tthe first xline bounds (%i, %i)\n", aux1, aux2);
  amplitude_aux1(lead_to, aux1, aux2);
  fprintf(logg,"amplitude loop j: bounds (%i, %i)\n", aux1, aux2);
  amplitude_aux2(lead_to, aux1, aux1, aux2);
  fprintf(logg,"\tthe first xline bounds (%i, %i)\n", aux1, aux2);
#endif

  int ifrom=0, ito=0;
  amplitude_aux1(lead_from, ifrom, ito);
  for (int i=ifrom;i<=ito;i++) {
    int iyfrom=0, iyto=0;
    amplitude_aux2(lead_from, i, iyfrom, iyto);
    for (int iy=iyfrom;iy<=iyto;iy++) {
      int wi=jcor(lead_from,i,iy);
      if (wi==-1) continue;
      int jfrom=0, jto=0;
      amplitude_aux1(lead_to, jfrom, jto);
      for (int j=jfrom;j<=jto;j++) {
        int jyfrom=0, jyto=0;
        amplitude_aux2(lead_to, j, jyfrom, jyto);
        for (int jy=jyfrom;jy<=jyto;jy++) {
          int wj=jcor(lead_to,j,jy);
          if (wj==-1) continue;
          for (int is=0;is<2;is++) {
            for (int js=0;js<2;js++) {
              content x=g(xlines-1,j,i,jy,iy,js,is);
              x*=lead[lead_from].mode[mode_from].wavefunction[wi+is*lead[lead_from].points];
              x*=conj(lead[lead_to].mode[mode_to].wavefunction[wj+js*lead[lead_to].points]);
              res+=x;
#ifdef TRANS3
    if (abs(x)>1e-12) {
      fprintf(logg,"amplitude entry:\n");
      fprintf(logg,"G(i=%i, yi=%i, si=%i, j=%i, yj=%i, sj=%i)=%e%+ei\n",i, iy, is, j, jy, js, g(xlines-1,j,i,jy,iy,js,is).real(), g(xlines-1,j,i,jy,iy,js,is).imag());
      fprintf(logg,"transversal wavefunction of from-mode %i at wi=%i : %.2e%+.2ei\n", mode_from, wi, lead[lead_from].mode[mode_from].wavefunction[wi+is*lead[lead_from].points].real(),  lead[lead_from].mode[mode_from].wavefunction[wi+is*lead[lead_from].points].imag());
      fprintf(logg,"transversal wavefunction of to-mode %i at wj=%i : %.2e%+.2ei\n", mode_to, wj, lead[lead_to].mode[mode_to].wavefunction[wj+js*lead[lead_to].points].real(), lead[lead_to].mode[mode_to].wavefunction[wj+js*lead[lead_to].points].imag());
    }
#endif
            }
          }
        }
      }
    }
  }
  
  prec v1=lead[lead_from].mode[mode_from].velocity;
  prec v2=lead[lead_to].mode[mode_to].velocity;
  //res*=content(0,1)*2.0*(hy)*sqrt(v1*v2)/(hx*hy);
  res*=content(0,1)*units.hbar;
  res*=sqrt(v1*v2)*units.velocity;
  res*=1/(units.energy*units.length*hx*units.length*hy);    //for the greens function - inverse of the Hamiltonian normalized to a delta function along both coordinates
  res*=1/(units.length*hy);                //for the two transversal wavefunctions
  res*=pow(units.length*hy,2);//for the two integrations along transversal direction

  if (lead_from==lead_to && mode_from==mode_to) res-=1.0;
  return(res);
}

//fills in the scattering matrix
content** transport_class::scattering_matrix(int* M)
{
#ifdef TRANS1
  fprintf(logg,"transmission probability:\n\tfrom lead\tfrom mode\tto lead\tto mode\n");
#endif
  for (int lf=1;lf<Nl;lf++) {
    M[lf]=lead[lf].propagating_modes;
    for (int lt=1;lt<Nl;lt++) {
      for (int mf=0;mf<lead[lf].propagating_modes;mf++) {
        for (int mt=0;mt<lead[lt].propagating_modes;mt++) {
          content x=0;
          s_matrix.s[lt*Nl+lf][mt*lead[lf].propagating_modes+mf]=x=amplitude(lf,mf,lt,mt);
          prec T=(x*conj(x)).real();
#ifdef TRANS1
          fprintf(logg,"%e\t%i\t\t%i\t\t%i\t\t%i\n",T,lf,mf,lt,mt);
#endif
        }
      }
    }
  }
  //check_unitarity();
  return(s_matrix.s);
}

//multiplies 2by2 matrices res = a x b
void M2M_aux(content* a, content* b, content* res )
{
  content aux[4];
  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++) {
      aux[i*2+j]=0;
      for (int k=0;k<2;k++) {
        aux[i*2+j]+=a[i*2+k]*b[k*2+j];
      }
    }
  }
  for (int i=0;i<4;i++) res[i]=aux[i];
}

//sums all mode to mode transmission according to the trace formula
//T_mn= Tr[t_mn^dagger sigma_alpha t_mn sigma_beta]
//indexed by lead (from, to) and sigma = 0,1,2,3 = id, x,y,z are Pauli matrices
prec transport_class::lead2lead_spinresolved(int lt, int alpha, int lf, int beta)
{
  content pauli[4*4]={1.0, 0, 0, 1.0,\
                      0, 1.0, 1.0, 0,\
                      0, -content(0,1), content(0,1), 0,\
                      1.0, 0, 0, -1.0};
  if (lf<1 || lf>=Nl || lt<1 || lt>=Nl ) {
    printf("invalid lead (from=%i/to=%i) in lead2lead_spinresolved!\n",lf,lt); exit(1);}
  if (alpha<0 || alpha>3 || beta<0 || beta>3 ) {
    printf("invalid sigma matrix (alpha=%i/beta=%i) in lead2lead_spinresolved!\n",alpha,beta); exit(1);}

  content S[4], Sdag[4], aux[4];
  prec res=0;
  //sum through all orbital modes
  for (int mof=0;mof<lead[lf].orbital_modes;mof++) {
    for (int mot=0;mot<lead[lt].orbital_modes;mot++) {
      scattering_spin_submatrix(lt, mot, lf, mof,S);//2x2 matrix
      for (int i=0;i<2;i++) for (int j=0;j<2;j++) Sdag[i*2+j]=conj(S[j*2+i]);//matrix conjugate
      M2M_aux(Sdag,pauli+4*alpha,aux);
      M2M_aux(aux,S,aux);
      M2M_aux(aux,pauli+4*beta,aux);
      res+=aux[0*2+0].real()+aux[1*2+1].real();//trace
    }
  }
  return(res);
}

//fills 2x2 scattering submatrix indexed by lead and orbital(!) mode
//if the time reversal pair is not propagating, the corresponding transmission amplitudes are zero
//S is allocated outside
content transport_class::scattering_spin_submatrix(int lt, int mot, int lf, int mof, content *S)
{
  if (lf<1 || lf>=Nl || lt<1 || lt>=Nl ) {
  printf("invalid lead (%i/%i) in scattering_spin_submatrix!\n",lf,lt); exit(1);}
  for (int i=0;i<4;i++) S[i]=0;

  //find propagating modes labels
  int mf[2]={-1,-1}, mt[2]={-1,-1};
  for (int m=0;m<lead[lf].propagating_modes;m++) if (lead[lf].mode[m].orb==mof) {mf[0]=m;break;}
  for (int m=0;m<lead[lt].propagating_modes;m++) if (lead[lt].mode[m].orb==mot) {mt[0]=m;break;}
  if (mf[0]==-1 || mt[0]==-1) return(0);//at least one is not propagating

  //the time reversal pair modes
  mf[1]=lead[lf].mode[mf[0]].trs;
  mt[1]=lead[lt].mode[mt[0]].trs;

  //S[0*2+0]=s_matrix.s[lt*Nl+lf][mt[0]*s_matrix.M[lf]+mf[0]];
  //if (mf[1]!=-1) S[0*2+1]=s_matrix.s[lt*Nl+lf][mt[0]*s_matrix.M[lf]+mf[1]];
  //if (mt[1]!=-1) S[1*2+0]=s_matrix.s[lt*Nl+lf][mt[1]*s_matrix.M[lf]+mf[0]];
  //if (mf[1]!=-1 && mt[1]!=-1) S[1*2+1]=s_matrix.s[lt*Nl+lf][mt[1]*s_matrix.M[lf]+mf[1]];
  for (int i=0;i<2;i++) {
    if (mf[i]==-1) continue;
    for (int j=0;j<2;j++) {
      if (mt[j]==-1) continue;
      content x=s_matrix.s[lt*Nl+lf][mt[j]*s_matrix.M[lf]+mf[i]];
      S[0*2+0]+=x*conj(lead[lf].mode[mf[i]].normalized_spinor.up)*lead[lt].mode[mt[j]].normalized_spinor.up;
      S[0*2+1]+=x*conj(lead[lf].mode[mf[i]].normalized_spinor.down)*lead[lt].mode[mt[j]].normalized_spinor.up;
      S[1*2+0]+=x*conj(lead[lf].mode[mf[i]].normalized_spinor.up)*lead[lt].mode[mt[j]].normalized_spinor.down;
      S[1*2+1]+=x*conj(lead[lf].mode[mf[i]].normalized_spinor.down)*lead[lt].mode[mt[j]].normalized_spinor.down;
    }
  }
  return(lead[lt].mode[mt[0]].trs_phase/lead[lf].mode[mf[0]].trs_phase);
}

//matrix element of the Hamitlonian plus the self energy
inline content transport_class::HpS(int i, int j, int iy, int jy, int is, int js)
{
#ifdef TRANSCHECK
  if (iy>=dim[i] || jy>=dim[j]) printf("out of range in HpS\n");
#endif
  content res;
  int ai,aj,aiy,ajy;
  region->sec2active(0,i,iy,ai,aiy);
  region->sec2active(0,j,jy,aj,ajy);
  res=region->Hsecsec(0,i,iy,is,0,j,jy,js);
  if (self_filled[labelM(i,j)]) res+=self(i,j,iy,jy,is,js);
  return(res);
}

//fills the matrix representing the E-H-Sigma of the region without leads part (i,j)
inline void transport_class::fill_EmHmS(content* res, int i, int j, prec Ef)
{
  for (int iy=0;iy<dim[i];iy++) {
    for (int jy=0;jy<dim[j];jy++) {
      int lab=label(i,j,iy,jy,0,0);
      res[lab]=-HpS(i,j,iy,jy,0,0);
      if (i==j && iy==jy) res[lab]+=Ef;
      lab=label(i,j,iy,jy,0,1);
      res[lab]=-HpS(i,j,iy,jy,0,1);
      lab=label(i,j,iy,jy,1,0);
      res[lab]=-HpS(i,j,iy,jy,1,0);
      lab=label(i,j,iy,jy,1,1);
      res[lab]=-HpS(i,j,iy,jy,1,1);
      if (i==j && iy==jy) res[lab]+=Ef;
    }
  }
}

//multiply off diagonal part of the Hamiltonian conecting the last added line i with the green function: 
//result=matrix (i,j)=sum_ip V_iip G_ipj, ip<i
//g[up=i-1] is used
void transport_class::VG(content* res,int i,int j,bool update,content coef)
{
//  if (i==16) fprintf(logg,"VG called: i=%i, j=%i\n",i,j);
  for (int iy=0;iy<dim[i];iy++) {
    for (int is=0;is<2;is++) {
      for (int jy=0;jy<dim[j];jy++) {
        for (int js=0;js<2;js++) {
          int lab=label(i,j,iy,jy,is,js);
          if (!update) res[lab]=0;
          for (int ip=0;ip<i;ip++) {
            if (iszeroh(i,ip)) continue;
            for (int ipy=0;ipy<dim[ip];ipy++) {
              if (iszerov(i,ip,iy,ipy)) continue;
  //if (i==16) fprintf(logg,"VG entry: i=%i, yi=%i, j=%i, yj=%i, ip=%i, yip=%i\n",i,yi,j,yj,ip,yip);
              //for (int ips=0;ips<2;ips++) {
                //printf("from VG\n");
                //res[lab]+=HpS(i,ip,iy,ipy,is,ips)*g(ip,j,ipy,jy,ips,js)*coef;
                res[lab]+=HpS(i,ip,iy,ipy,is,0)*g(i-1,ip,j,ipy,jy,0,js)*coef;//sum through ips
                res[lab]+=HpS(i,ip,iy,ipy,is,1)*g(i-1,ip,j,ipy,jy,1,js)*coef;
              //}
            }
          }
        }
      }
    }
  }
}

//multiply off diagonal part of the Hamiltonian conecting the last added line j with the green function: 
//result=matrix (i,j)=sum_jp G_ijp V_jpj, jp<j
void transport_class::GV(content* res,int i,int j,bool update,content coef)
{
  for (int iy=0;iy<dim[i];iy++) {
    for (int jy=0;jy<dim[j];jy++) {
      //is=0, js=0
      int lab=label(i,j,iy,jy,0,0);
      if (!update) res[lab]=0;
      for (int jp=0;jp<j;jp++) { if (iszeroh(jp,j)) continue;
        for (int jpy=0;jpy<dim[jp];jpy++) { if (iszerov(jp,j,jpy,jy)) continue;
            res[lab]+=g(j-1,i,jp,iy,jpy,0,0)*HpS(jp,j,jpy,jy,0,0)*coef;//sum over jps
            res[lab]+=g(j-1,i,jp,iy,jpy,0,1)*HpS(jp,j,jpy,jy,1,0)*coef;
        }
      }
      //is=0, js=1
      lab=label(i,j,iy,jy,0,1);
      if (!update) res[lab]=0;
      for (int jp=0;jp<j;jp++) { if (iszeroh(jp,j)) continue;
        for (int jpy=0;jpy<dim[jp];jpy++) { if (iszerov(jp,j,jpy,jy)) continue;
            res[lab]+=g(j-1,i,jp,iy,jpy,0,0)*HpS(jp,j,jpy,jy,0,1)*coef;//sum over jps
            res[lab]+=g(j-1,i,jp,iy,jpy,0,1)*HpS(jp,j,jpy,jy,1,1)*coef;
        }
      }
      //is=1, js=0
      lab=label(i,j,iy,jy,1,0);
      if (!update) res[lab]=0;
      for (int jp=0;jp<j;jp++) { if (iszeroh(jp,j)) continue;
        for (int jpy=0;jpy<dim[jp];jpy++) { if (iszerov(jp,j,jpy,jy)) continue;
            res[lab]+=g(j-1,i,jp,iy,jpy,1,0)*HpS(jp,j,jpy,jy,0,0)*coef;//sum over jps
            res[lab]+=g(j-1,i,jp,iy,jpy,1,1)*HpS(jp,j,jpy,jy,1,0)*coef;
        }
      }
      //is=1, js=1
      lab=label(i,j,iy,jy,1,1);
      if (!update) res[lab]=0;
      for (int jp=0;jp<j;jp++) { if (iszeroh(jp,j)) continue;
        for (int jpy=0;jpy<dim[jp];jpy++) { if (iszerov(jp,j,jpy,jy)) continue;
            res[lab]+=g(j-1,i,jp,iy,jpy,1,0)*HpS(jp,j,jpy,jy,0,1)*coef;//sum over jps
            res[lab]+=g(j-1,i,jp,iy,jpy,1,1)*HpS(jp,j,jpy,jy,1,1)*coef;
        }
      }
    }
  }
}
//multiply off diagonal part of the Hamiltonian conecting the last added line i with the green function and the hamiltonian again: 
//result=matrix (i,j)=sum_ipjp V_iip G_ipjp V_jpj
//note it must be i=j
void transport_class::VGV(content* res,int i,int j,bool update,content coef)
{
  //printf("VGV called with i=%i, j=%i\n",i,j);
  for (int iy=0;iy<dim[i];iy++) {
    for (int is=0;is<2;is++) {
      for (int jy=0;jy<dim[j];jy++) {
        for (int js=0;js<2;js++) {
          int lab=label(i,j,iy,jy,is,js);
          if (!update) res[lab]=0;
          for (int ip=0;ip<i;ip++) { 
            if (iszeroh(i,ip)) continue;
            for (int ipy=0;ipy<dim[ip];ipy++) {
              if (iszerov(i,ip,iy,ipy)) continue;
              for (int ips=0;ips<2;ips++) {
                content GVIPJ=0;
                for (int jp=0;jp<j;jp++) {
                  if (iszeroh(jp,j)) continue;
                  for (int jpy=0;jpy<dim[jp];jpy++) {
                    if (iszerov(jp,j,jpy,jy)) continue;
                    //for (int jps=0;jps<2;jps++) {
                      //printf("from VGV\n");
                      //GVIPJ+=g(ip,jp,ipy,jpy,ips,jps)*HpS(jp,j,jpy,jy,jps,js);
                      GVIPJ+=g(i-1,ip,jp,ipy,jpy,ips,0)*HpS(jp,j,jpy,jy,0,js);//sum through jps
                      GVIPJ+=g(i-1,ip,jp,ipy,jpy,ips,1)*HpS(jp,j,jpy,jy,1,js);
                    //}
                  }
                }
                res[lab]+=HpS(i,ip,iy,ipy,is,ips)*GVIPJ*coef;
              }
            }
          }
        }
      }
    }
  }
}




prec transport_class::lead2lead(int from, int to)
{
  prec x=0;
  for (int mf=0;mf<lead[from].propagating_modes;mf++) {
    for (int mt=0;mt<lead[to].propagating_modes;mt++) {
      content aux=amplitude(from,mf,to,mt);
      x+=(aux*conj(aux)).real();
    }
  }
  return(x);
}

prec transport_class::trace(int from, int to)
{
  content sum=0;
  for (int i=0;i<xlines;i++) {
    for (int yi=0;yi<dim[i];yi++) {
      for (int si=0;si<2;si++) {
        for (int j=0;j<xlines;j++) {
          for (int yj=0;yj<dim[j];yj++) {
            for (int sj=0;sj<2;sj++) {
              for (int k=0;k<xlines;k++) {
                for (int yk=0;yk<dim[k];yk++) {
                  for (int sk=0;sk<2;sk++) {
                    for (int l=0;l<xlines;l++) {
                      for (int yl=0;yl<dim[l];yl++) {
                        for (int sl=0;sl<2;sl++) {
                          if (to==0) {if (i>b-1 || j>b-1) continue;}
                          else if (i<xlines-b || j<xlines-b) continue;
                          if (from==0) {if (k>b-1 || l>b-1) continue;} 
                          else if (k<xlines-b || l<xlines-b) continue;
                          content st=self(i,j,yi,yj,si,sj);
                          content stt=self(j,i,yj,yi,sj,si);
                          content gr=g(xlines-1,j,k,yj,yk,sj,sk);
                          content sf=self(k,l,yk,yl,sk,sl);
                          content sft=self(l,k,yl,yk,sl,sk);
                          content ga=g(xlines-1,i,l,yi,yl,si,sl);
                          content res=(st-conj(stt))*gr*(sf-conj(sft))*conj(ga)*(-1.0);
                          sum+=res;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (abs(sum.imag())>1e-5) printf("trace result not real: %f%+fi\n",sum.real(),sum.imag());
  return(sum.real());
}

extern content potentialc(prec, prec);

void transport_class::draw_potential(FILE* file, int where)
{
  region->draw_function(potentialc, -1,1, file, where);
}

//transition and reflection amplitudes for a low(w)-high(wp)-low(l) region with wavevectors k/kp for low/high region
//the incoming wave has amplitude 1 at the beginning of the first low region
void one_barrier(content &t, content &r, content k, content kp, prec w, prec wp, prec l) 
{
  content i=content(0,1);
  content den=(k-kp)*(k-kp)*exp(2.0*i*kp*wp)-(k+kp)*(k+kp);
  r=exp(i*2.0*k*w)*(k*k-kp*kp)*(exp(2.0*i*kp*wp)-1.0)/den;
  t=-4.0*exp(-i*(k-kp)*wp)*k*kp/den;
  t*=exp(i*k*(w+wp+l));
}

//input is a transfer matrix M giving (most far) right amplitudes in terms of (most far) left amplitudes
//call updates this matrix adding an additional barrier step of the one_barrier type
void update(content* M,content k, content kp, prec w, prec wp, prec l)
{
  content t,r;
  one_barrier(t,r,k,kp,w,wp,l);
  content Mp[4];
  content T[4];
  T[0*2+0]=-r; T[0*2+1]=1.0;
  T[1*2+0]=-r*r+t*t; T[1*2+1]=r;
  
  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++) {
      Mp[i*2+j]=0.0;
      for (int k=0;k<2;k++) Mp[i*2+j]+=T[i*2+k]*M[k*2+j]/t;
    }
  }
  for (int i=0;i<4;i++) M[i]=Mp[i];
}

//reflection and transmission amplitude of a region of a total length x, composed of step-like barriers: low(w)-high(wp)
//defined by energies Ef=particle energy in low region, Eb=barrier height
//everything in internal units (energy, length)
void rectangular_barriers(content &t, content &r, prec Ef, prec Eb, prec w, prec wp, prec x)
{  
  prec conv=sqrt(2*units.energy*units.m_e*parameters.read("m_eff","m_e")/pow(units.hbar*units.wavevector,2));
  content k=conv*sqrt(content(Ef,0));
  content kp=conv*sqrt(content(Ef-Eb,0));
  
  sprintf(buffer,"rectangular barrier\n\tw %f, wp %f, x %f [nm], Ef %f, Eb %f [meV], k %f%+fi, kp %f%+fi [1/nm]\n",w*units.length/units.nm, wp*units.length/units.nm, x*units.length/units.nm, Ef*units.energy/units.meV, Eb*units.energy/units.meV, k.real()*units.wavevector*units.nm, k.imag()*units.wavevector*units.nm, kp.real()*units.wavevector*units.nm, kp.imag()*units.wavevector*units.nm);
  splachni(logg,buffer,4);
  if (Ef<0) {t=0;r=1.0;return;}
  
  if (x<0) {printf("x can not be negative in rectangular barriers\n");exit(1);}
  
  int n=(int) floor(x/(w+wp));
  x-=n*(w+wp);
  sprintf(buffer, "number of w+wp regions (scatterers) is %i, residual length is %f [nm]\n",n,x*units.length/units.nm);splachni(logg,buffer,1+4);
  content M[4];
  M[0*2+0]=M[1*2+1]=content(1,0);
  M[0*2+1]=M[1*2+0]=content(0,0);
  for (int i=0;i<n;i++) {
    update(M,k,kp,w/2,wp,w/2);
    content aux;
    aux=M[0*2+0];M[0*2+0]=M[1*2+0];M[1*2+0]=aux;
    aux=M[0*2+1];M[0*2+1]=M[1*2+1];M[1*2+1]=aux;
  }

  if (x<-1e-10 or x>w+wp+1e-10) {printf("x can not be %e here in rectangular barriers\n",x);exit(1);}
  
  if (x>w/2) {
    x-=w/2;
    if (x>wp) {
      x-=wp;
      if (x>w/2) {printf("can not be!\n");exit(1);}
      update(M,k,kp,w/2,wp,x);
    }
    else update(M,k,kp,w/2,x,0);
  }
  else update(M,k,kp,x,0,0);
  t=-(M[0*2+0]*M[1*2+1]-M[0*2+1]*M[1*2+0])/M[0*2+1];
  r=-M[0*2+0]/M[0*2+1];

  sprintf(buffer,"rectangular barrier result: |t|^2=%e, |r|^2=%e\n",abs(t)*abs(t), abs(r)*abs(r));
  splachni(logg,buffer,4);

  return;
}

//step barriers
void transport_class::check_GF(FILE* file, bool full)
{
  extern prec GF_check_Eb, GF_check_w, GF_check_wp;
  extern bool GF_check_flag;
  bool temp=GF_check_flag;
  GF_check_flag=true;
  prec Ef=parameters.read("fermi_e","meV");
  prec ll=lead[1].x*units.length/units.nm;
  prec rl=lead[2].x*units.length/units.nm;
  GF_check_Eb=1;
  GF_check_w=30;
  GF_check_wp=60-30;
  prec x=(rl-ll);
  printf("series of step barriers (low %f - high %f of total distance %f [nm]) of height %f for Ef=%f meV\n", GF_check_w, GF_check_wp, x, GF_check_Eb, Ef);

  Ef*=units.meV/units.energy;
  GF_check_Eb*=units.meV/units.energy;
  GF_check_w*=units.nm/units.length;
  GF_check_wp*=units.nm/units.length;
  x*=units.nm/units.length;

  prepare_geometry(full);
  prepare_with_E();
  recursive_green(full);

  printf("1. rectangular barrier\n");
  prec muB=units.e*units.hbar/units.m_e/4.0*parameters.read("B","Tesla")*parameters.read("gz","-")/units.energy;
  content ta1, ra1, ta2, ra2;
  rectangular_barriers(ta1,ra1,Ef-lead[1].mode[0].energy,GF_check_Eb,GF_check_w,GF_check_wp,x);
  rectangular_barriers(ta2,ra2,Ef+2*muB-lead[1].mode[0].energy,GF_check_Eb,GF_check_w,GF_check_wp,x);
  
  //numerics
  printf("2. numerics using Fisher-Lee\n");
  content tn1=amplitude(1,0,2,0);
  content rn1=amplitude(1,0,1,0);
  content tn2=amplitude(1,1,2,1);
  content rn2=amplitude(1,1,1,1);

  //1D formula
  printf("3. 1D formula - meaninful only if the space is 1 point along the y axis\n");
  content t1D,r1D;
  oneDamplitude(t1D, r1D);
  
  //numerics
  printf("4. numerics using trace - takes too long, uncomment if... maybe works\n");
  prec Tn=0;//trace(0,1);
  prec Rn=0;//trace(0,0);
  //sprintf(buffer, "(0->0) %f, (0->1) %f, (1->0) %f, (1->1) %f\n",Rn,Tn,trace(1,0),trace(1,1));
  //splachni(logg,buffer,1+4); 
  
  sprintf(buffer, "%e mode 0: anal:%e vs num:%e mode 1: anal:%e vs num:%e\n", Ef/units.meV*units.energy, abs(ta1*conj(ta1)), abs(tn1*conj(tn1)), abs(ta2*conj(ta2)), abs(tn2*conj(tn2)));
  splachni(file,buffer,2);
  
  printf("\n\n\t\ttransmission\t\t\t\treflection\n");
  printf("1.(up)\t\t%f%+fi (sq:%.3e)\t%f%+fi (sq:%.3e)\n",ta1.real(),ta1.imag(),abs(ta1)*abs(ta1),ra1.real(),ra1.imag(),abs(ra1)*abs(ra1));
  printf("2.(up)\t\t%f%+fi (sq:%.3e)\t%f%+fi (sq:%.3e)\n",tn1.real(),tn1.imag(),abs(tn1)*abs(tn1),rn1.real(),rn1.imag(),abs(rn1)*abs(rn1));
  printf("3.\t\t%f%+fi (sq:%.3e)\t%f%+fi (sq:%.3e)\n",t1D.real(),t1D.imag(),abs(t1D)*abs(t1D),r1D.real(),r1D.imag(),abs(r1D)*abs(r1D));
  printf("1.(down)\t%f%+fi (sq:%.3e)\t%f%+fi (sq:%.3e)\n",ta2.real(),ta2.imag(),abs(ta2)*abs(ta2),ra2.real(),ra2.imag(),abs(ra2)*abs(ra2));
  printf("2.(down)\t%f%+fi (sq:%.3e)\t%f%+fi (sq:%.3e)\n",tn2.real(),tn2.imag(),abs(tn2)*abs(tn2),rn2.real(),rn2.imag(),abs(rn2)*abs(rn2));

  int M[2];
  scattering_matrix(M);
  check_trs();

  //printf("4.\t\t\t\t (sq:%.3e)\t\t\t\t (sq:%.3e)\n",Tn, Rn);
  clean_with_E();
  clean_geometry();
  GF_check_flag=temp;
}

//extract transmission and reflection amplitudes for a one dimensional case (all xlines have length (dim) 1)
//(note: valid only for spin along z)
void transport_class::oneDamplitude(content &t, content &r)
{
  content gn1=g(xlines-1,xlines-1,0,0,0,0,0);
  content g11=g(xlines-1,0,0,0,0,0,0);
  prec v1=lead[1].mode[0].velocity;
  prec v2=lead[2].mode[0].velocity;
  prec k=lead[1].mode[0].wavevector;
  prec conv=units.hbar*units.velocity/(units.energy*units.length)/hx;
  t=conv*content(0,1)*sqrt(v1*v2)*gn1;//*exp(content(0,-k*hx*(xlines-2)));
  r=conv*content(0,1)*sqrt(v1*v2)*g11-1.0;
  printf("\t1D amplitudes:\n\t\tgreen function values\n\tGN1=%f%+fi, G11=%f%+fi\n",gn1.real(),gn1.imag(),g11.real(),g11.imag());
  printf("\t\twave velocities %f, %f, wave vector %f, grid step %f\n",v1,v2,k,hx);
}

void transport_class::oneDspindensity(FILE* file)
{
  //prec vup=lead[1].mode[0].velocity;
  //prec vdown=lead[1].mode[1].velocity;
  content g11up=g(xlines-1,0,0,0,0,0,0)*0.0+1.0;
  content g11down=g(xlines-1,0,0,0,0,1,1)*0.0+1.0;
  for (int i=0;i<xlines;i++) {
    content gi1up=g(xlines-1,i,0,0,0,0,0);
    content gi1down=g(xlines-1,i,0,0,0,1,1);
    double sz=(gi1up/g11up*conj(gi1up/g11up)-gi1down/g11down*conj(gi1down/g11down)).real();
    fprintf(file, "%i %e\n",i,sz);
  }
}


//makes a time reversal of a mode of lead and returns the mode to which it is transformed
//returns the overlap <mode_out | T mode_in>
//if TRS mode is not in the lead (not propagating), returns -1 in mode out and zero for the coefficient
content transport_class::mode_trs(int l, int mode_in, int& mode_out)
{
  if (l<1 || l>=Nl || mode_in<0 ||mode_in>=lead[l].propagating_modes) {printf("mode_trs: wrong lead(%i) or its mode (%i)\n",l,mode_in);exit(0);}
  int N=lead[l].points;
  //make time reversal of a mode
  content* trs_wf=new content[N*2];
  for (int i=0;i<N;i++) {
    for (int is=0;is<2;is++) {
      trs_wf[i+N*is]=conj(lead[l].mode[mode_in].wavefunction[i+N*(1-is)])*(1.0*(is==0)-1.0*(is==1));
    }
  }
#ifdef TRANS3
  fprintf(logg,"time reversal of the lead %i mode %i:\n",l,mode_in);
  for (int is=0;is<2;is++) {
    if (is==0) fprintf(logg,"orig-up:\t"); else fprintf(logg,"orig-down:\t");
    for (int i=0;i<N;i++) fprintf(logg,"%.3e%+.3e ",lead[l].mode[mode_in].wavefunction[i+N*is].real(), lead[l].mode[mode_in].wavefunction[i+N*is].imag());
    fprintf(logg,"\n");
  }
  for (int is=0;is<2;is++) {
    if (is==0) fprintf(logg,"trs-up:\t"); else fprintf(logg,"trs-down:\t");
    for (int i=0;i<N;i++) fprintf(logg,"%.3e%+.3e ",trs_wf[i+N*is].real(),trs_wf[i+N*is].imag());
    fprintf(logg,"\n");
  }
#endif

  //find to which it is closest
  content omax=0, olap;
  for (int m=0;m<lead[l].propagating_modes;m++) {
    olap=0;
    for (int i=0;i<2*N;i++) olap+=trs_wf[i]*conj(lead[l].mode[m].wavefunction[i]);
    if (abs(olap)>abs(omax)) {omax=olap;mode_out=m;}
#ifdef TRANS2
    fprintf(logg,"time reversed mode overlap (abs) with mode %i: %e:\n",m,abs(olap));
#endif
  }
  //check that the overlap is large, otherwise finish with result "not found"
  if (abs(omax)<0.5) {mode_out=-1;olap=0;}
#ifdef TRANS2
    fprintf(logg,"time reversed mode %i is closest to mode %i with fidelity %e:\n",mode_in,mode_out,abs(omax));
#endif
  delete trs_wf;
  return(omax);
}

//assign modes according to spin pair - using the trs symmetry
void transport_class::pair_modes()
{
  for (int l=1;l<Nl;l++) {
    int orb=0;
    for (int m=0;m<lead[l].propagating_modes;m++) {
      if (lead[l].mode[m].orb!=-1) continue; //already assigned - skip
      int mt,mo;
      content phase=mode_trs(l,m,mt);
      lead[l].mode[m].orb=orb;
      lead[l].mode[m].spin=0;
      lead[l].mode[m].trs=-1;
      lead[l].mode[m].trs_phase=phase;

      if (mt!=-1) {//a spin flipped comrade is found
        lead[l].mode[m].trs=mt;
        lead[l].mode[mt].trs=m;
        lead[l].mode[mt].orb=orb;
        lead[l].mode[mt].spin=1;
        //check: T(T(m)) should be the original mode m
        content phase2=mode_trs(l,mt,mo);
        if (mo!=m) {printf("pair_modes: lead %i mode %i: trs(%i)=%i, but trs(%i)=%i!!!\n",l,m,m,mt,mt,mo);exit(1);}
        lead[l].mode[mt].trs_phase=phase2;
      }
      orb++;
    }
    lead[l].orbital_modes=orb;
  }
#ifdef TRANS2
  fprintf(logg,"time reversed pairs labels of modes:\n");
  for (int l=1;l<Nl;l++) {
    fprintf(logg,"\tlead %i with %i orbital (%i total) propagating modes\n",l,lead[l].orbital_modes, lead[l].propagating_modes);
    for (int m=0;m<lead[l].propagating_modes;m++) {
      fprintf(logg,"\t\tpropagating mode %i: orb=%i, spin=%i, trs pair=%i, trs phase=(%.3e%+.3ei)\n", m, lead[l].mode[m].orb, lead[l].mode[m].spin, lead[l].mode[m].trs, lead[l].mode[m].trs_phase.real(), lead[l].mode[m].trs_phase.imag());
    }
  }
#endif
}


//check for self-duality (TRS symmetry) of the scattering matrix (logs residuals)
prec transport_class::check_trs()
{
  //construct the scattering matrix from its self-dual property
  content** dual_s=new content*[Nl*Nl];
  for (int l1=1;l1<Nl;l1++) {
    int dim1=lead[l1].propagating_modes;
    for (int l2=1;l2<Nl;l2++) {
      int dim2=lead[l2].propagating_modes;
      int lab=l1*Nl+l2;
      int labt=l2*Nl+l1;
      dual_s[lab]=new content[max(dim1*dim2,1)];
      for (int m1=0;m1<dim1;m1++) {
        int mt1=lead[l1].mode[m1].trs;
        for (int m2=0;m2<dim2;m2++) {
          int mt2=lead[l2].mode[m2].trs;
          if (mt1==-1 || mt2==-1) dual_s[lab][m1*dim2+m2]=0;
          else dual_s[lab][m1*dim2+m2]=s_matrix.s[labt][mt2*dim1+mt1]* 
                 lead[l2].mode[mt2].trs_phase/lead[l1].mode[mt1].trs_phase;
        }
      }
    }
  }
  prec res=0;
  //compute norm of the difference
  for (int l1=1;l1<Nl;l1++) {
    int dim1=lead[l1].propagating_modes;
    for (int l2=1;l2<Nl;l2++) {
      int dim2=lead[l2].propagating_modes;
      int lab=l1*Nl+l2;
      for (int m1=0;m1<dim1;m1++) {
        for (int m2=0;m2<dim2;m2++) {
          res+=abs(dual_s[lab][m1*dim2+m2]-s_matrix.s[lab][m1*dim2+m2]);
#ifdef TRANS3
          fprintf(logg,"scattering matrix: element (l1=%i, m1=%i, l2=%i, m2=%i): self dual (%.3e%+.3ei) vs original (%.3e%+.3ei) abs(difference)=%e\n",l1,m1,l2,m2, dual_s[lab][m1*dim2+m2].real(), dual_s[lab][m1*dim2+m2].imag(), s_matrix.s[lab][m1*dim2+m2].real(), s_matrix.s[lab][m1*dim2+m2].imag() ,abs(dual_s[lab][m1*dim2+m2]-s_matrix.s[lab][m1*dim2+m2]));
#endif
        }
      }
    }
  }
  sprintf(buffer,"the scattering matrix is self-dual (TRS) with tolerance %e\n",res);
  splachni(logg,buffer,1+4);
  for (int l1=1;l1<Nl;l1++) for (int l2=1;l2<Nl;l2++) delete dual_s[l1*Nl+l2];
  delete dual_s;

  //check the 2x2 spin submatrices
  content S[4], St[4], aux[4];
  content sy[4]={0, -content(0,1), content(0,1), 0};

  prec tot=0;
  for (int l1=1;l1<Nl;l1++) {
    for (int l2=1;l2<Nl;l2++) {
      for (int mo1=0;mo1<lead[l1].orbital_modes;mo1++) {
        for (int mo2=0;mo2<lead[l2].orbital_modes;mo2++) {
          scattering_spin_submatrix(l2, mo2, l1, mo1, S);
          content phase=1.0;
          scattering_spin_submatrix(l1, mo1, l2, mo2, St);
          content x=St[0*2+1];
          St[0*2+1]=St[1*2+0];
          St[1*2+0]=x;
          M2M_aux(sy,St,aux);
          M2M_aux(aux,sy,aux);
          prec res=0;
#ifdef TRANS3
          fprintf(logg,"\tthe scattering spin submatrix (lead %i mode %i -> lead %i mode %i):\n", l1, mo1, l2, mo2);
          for (int i=0;i<4;i++) fprintf(logg,"(%+.3e%+.3ei) ",S[i].real(),S[i].imag());
          fprintf(logg,"\n");
          fprintf(logg,"\tits dual brother (lead %i mode %i -> lead %i mode %i) times the phase (%+.3e%+.3ei):\n", l2, mo2, l1, mo1, phase.real(),phase.imag());
          for (int i=0;i<4;i++) fprintf(logg,"(%+.3e%+.3ei) ",(phase*aux[i]).real(),(phase*aux[i]).imag());
          fprintf(logg,"\n");
#endif
          for (int i=0;i<4;i++) res+=abs(S[i]-aux[i]*phase);
#ifdef TRANS2
          fprintf(logg,"\tthe scattering spin submatrix (lead %i mode %i -> lead %i mode %i) is self-dual (TRS) with tolerance %e\n", l1, mo1, l2, mo2, res);
#endif
          tot+=res;
        }
      }
    }
  }
  sprintf(buffer,"the scattering spin submatrices are self-dual (TRS) with tolerance sum %e\n",tot);

  return(res);
}

//check for unitarity of the scattering matrix
prec transport_class::check_unitarity()
{
  //checking res rules (unitarity of S)
  prec res=0;
  for (int lf=1;lf<Nl;lf++) {
    for (int lt=1;lt<Nl;lt++) {
      for (int mf=0;mf<lead[lf].propagating_modes;mf++) {
        for (int mt=0;mt<lead[lt].propagating_modes;mt++) {
          content aux=0;
          if (mt==mf && lf==lt) aux=-1.0;
          for (int ls=1;ls<Nl;ls++) {
            for (int ms=0;ms<lead[ls].propagating_modes;ms++) {
              aux+=conj(s_matrix.s[ls*Nl+lt][ms*lead[lt].propagating_modes+mt])*s_matrix.s[ls*Nl+lf][ms*lead[lf].propagating_modes+mf];
            }
          }
#ifdef TRANS3
          fprintf(logg,"unitarity: entry (S^+S)_(%i-%i,%i-%i)=%e\n",lf,mf,lt,mt,abs(aux));
#endif
          res+=abs(aux);
        }
      }
    }
  }
#ifdef TRANS1
  sprintf(buffer,"the scattering matrix is unitary with unprecision %e\n",res);
  splachni(logg,buffer,1+2);
#endif
  return(res);
}

//converts two multiindexes into a label into a matrix stored as a row of numbers (vector)
inline int transport_class::label(int i, int j, int iy, int jy, int is, int js) 
{
  int ni=dim[i], nj=dim[j];
  //printf("label: yi=%i, si=%i, yj=%i, sj=%i, ni=%i, nj=%i\n",yi,si,yj,sj,ni,nj);
#ifdef TRANSCHECK
  if (iy>=ni || is>1 || jy>=nj || js>1 || is<0 || js<0 ||iy<0 || jy<0)
    {printf("label out of range:yi=%i, si=%i, yj=%i, sj=%i, ni=%i, nj=%i\n",iy,is,jy,js,ni,nj);exit(1);} 
#endif
  return((is*ni+iy)*nj*2+js*nj+jy);
}

//returns label in the matrixes space
inline int transport_class::labelM(int i, int j) 
{
  return(i*xlines+j);
}

inline int transport_class::labelg(int up, int i, int j) 
{
  if (method_full) return((up%2)*xlines*xlines+labelM(i,j));
  if (i!=0) i-=up-b;
  if (j!=0) j-=up-b;
  return(up*(b+1)*(b+1)+i*(b+1)+j);
}

inline content* transport_class::g(int up, int i, int j) 
{
  return(g_array[labelg(up,i,j)]);
}

inline content transport_class::g(int up, int i, int j, int iy, int jy, int is, int js) 
{
#ifdef TRANSCHECK
  if (iy>=dim[i] || jy>=dim[j]) {printf("out of range in g: i=%i, yi=%i, j=%i, yj=%i\n",i,iy,j,jy);exit(1);}
#endif
  return(g(up,i,j)[label(i,j,iy,jy,is,js)]);
}

inline content* transport_class::self(int i, int j) 
{
  return(self_array[labelM(i,j)]);
}
inline content transport_class::self(int i, int j, int iy, int jy, int is, int js) 
{
#ifdef TRANSCHECK
  if (iy>=dim[i] || jy>=dim[j]) {printf("out of range in self: i=%i, yi=%i, j=%i, yj=%i\n",i,iy,j,jy);exit(1);}
#endif
  return(self(i,j)[label(i,j,iy,jy,is,js)]);
}

//if the two points are further apart then b in any direction returns true
//input are sector coordinates
inline bool transport_class::iszeroH(int sec1, int si1, int sj1, int sec2, int si2, int sj2)
{
  int i1,j1,i2,j2;
  region->sec2grid(sec1,si1,sj1,i1,j1);
  region->sec2grid(sec2,si2,sj2,i2,j2);
  if (abs(i1-i2)<=b && abs(j1-j2)<=b) return(false);
  return(true);
}

//true if both the hamiltonian and self energy are zero for these two sector 0 coordinates
inline bool transport_class::iszeroh(int i, int j)
{
  if (abs(i-j)<=b) return(false);    //hamiltonian is nonzero
  if (self_filled[labelM(i,j)]) return(false);
  return(true);
}

//true if both the hamiltonian and self energy are zero for these two sector 0 coordinates
inline bool transport_class::iszerov(int i, int j, int iy, int jy)
{
#ifdef TRANSCHECK
  if (iy>=dim[i] || jy>=dim[j]) {printf("out of range in iszerov\n");exit(1);}
#endif
  if (!iszeroH(0,i,iy,0,j,jy)) return(false);//hamiltonian is nonzero
  for (int l=0;l<Nl;l++) {
    if (region->secnearbysec(0,i,iy,l) && region->secnearbysec(0,j,jy,l)) return(false);
  }
  return(true);
}

//multiplication of two matrixes with given multiindexes
void transport_class::MM(content* res,content* a, content* b,int i1,int j1, int i2, int j2)
{
#ifdef TRANSCHECK
  if (j1!=i2) {printf("dimensions of matrixes for multiplication do not correspond in MM(...)\n");exit(1);}
#endif
  int ni1=dim[i1];
  int nj1=dim[j1];
  int nj2=dim[j2];
  for (int i1y=0;i1y<ni1;i1y++) {
    for (int j2y=0;j2y<nj2;j2y++) {
      int lab=label(i1,j2,i1y,j2y,0,0);
      res[lab]=0;
      for (int j1y=0;j1y<nj1;j1y++) {
        res[lab]+=a[label(i1,j1,i1y,j1y,0,0)]*b[label(j1,j2,j1y,j2y,0,0)];
        res[lab]+=a[label(i1,j1,i1y,j1y,0,1)]*b[label(j1,j2,j1y,j2y,1,0)];
      }
      lab=label(i1,j2,i1y,j2y,0,1);
      res[lab]=0;
      for (int j1y=0;j1y<nj1;j1y++) {
        res[lab]+=a[label(i1,j1,i1y,j1y,0,0)]*b[label(j1,j2,j1y,j2y,0,1)];
        res[lab]+=a[label(i1,j1,i1y,j1y,0,1)]*b[label(j1,j2,j1y,j2y,1,1)];
      }
      lab=label(i1,j2,i1y,j2y,1,0);
      res[lab]=0;
      for (int j1y=0;j1y<nj1;j1y++) {
        res[lab]+=a[label(i1,j1,i1y,j1y,1,0)]*b[label(j1,j2,j1y,j2y,0,0)];
        res[lab]+=a[label(i1,j1,i1y,j1y,1,1)]*b[label(j1,j2,j1y,j2y,1,0)];
      }
      lab=label(i1,j2,i1y,j2y,1,1);
      res[lab]=0;
      for (int j1y=0;j1y<nj1;j1y++) {
        res[lab]+=a[label(i1,j1,i1y,j1y,1,0)]*b[label(j1,j2,j1y,j2y,0,1)];
        res[lab]+=a[label(i1,j1,i1y,j1y,1,1)]*b[label(j1,j2,j1y,j2y,1,1)];
      }
    }
  }
}


//add two matrixes res=res+in with multiindexes i,j
void transport_class::add(content* res, content* in, int i, int j, prec coef)
{
  int ni=dim[i];
  int nj=dim[j];
  for (int iy=0;iy<ni;iy++) {
    for (int is=0;is<2;is++) { 
      for (int jy=0;jy<nj;jy++) {
        for (int js=0;js<2;js++) { 
          int lab=label(i,j,iy,jy,is,js);
          res[lab]+=in[lab]*coef;
        }
      }
    }
  }
}

//add a constant to a diagonal matrix with multiindexes i,j res=res+x*id 
void transport_class::addC(content* res, content x, int i, int j)
{
#ifdef TRANSCHECK
  if (i!=j) {printf("not a diagonal matrix in addC\n");exit(1);}
#endif
  int ni=dim[i];
  for (int iy=0;iy<ni;iy++) {
    for (int is=0;is<2;is++) { 
      res[label(i,i,iy,iy,is,is)]+=x;
    }
  }
}


//compute norm of a matrix in
prec transport_class::normM(content* in, int i, int j)
{
  int ni=dim[i];
  int nj=dim[j];
  prec rez=0;
  for (int iy=0;iy<ni;iy++) {
    for (int is=0;is<2;is++) { 
      for (int jy=0;jy<nj;jy++) {
        for (int js=0;js<2;js++) { 
          content x=in[label(i,j,iy,jy,is,js)];
          rez+=(x*conj(x)).real();
        }
      }
    }
  }
  return(sqrt(rez));
}

//checks that the greens function is the inverse of E-H-Sigma cut after line i. Computes the norm of
//[G (E-H-S)]_{jk} - delta(j,l) for each j,k from 0 to i (including)
//result is the list of residual norms into the log
void transport_class::check_it(int i) 
{
  prec Ef=parameters.read("fermi_e","meV")*units.meV/units.energy;
  sprintf(buffer,"checking up to xline %i :\n",i);
  splachni(logg,buffer,4);
  prec* norms=new prec [(i+1)*(i+1)];
  prec rmax=0;
  for (int j=0;j<i+1;j++) {
    int nj=dim[j]*2;
    for (int k=0;k<i+1;k++) {
      int nk=dim[k]*2;
      content* res=new content[nj*nk];
      content* aux=new content[nj*nk];
      for (int a=0;a<nj*nk;a++) res[a]=0;
      
      for (int jp=0;jp<i+1;jp++) {
        int njp=dim[jp]*2;
        content* aux0=new content[njp*nk];
        fill_EmHmS(aux0,jp,k,Ef);
        MM(aux,g(xlines-1,j,jp),aux0,j,jp,jp,k);
        add(res,aux,j,k,1);
        delete aux0;
      }
      if (j==k) addC(res,-1,j,k);
      rmax=max(norms[j*(i+1)+k]=normM(res,j,k),rmax);
      delete res;
      delete aux;
    }
  } 
  sprintf(buffer,"maximal norm of a residual: %e\n",rmax);splachni(logg,buffer,4);
  if (rmax>1e-12) {
    for (int j=0;j<i+1;j++) {
      for (int k=0;k<i+1;k++) {
        sprintf(buffer,"\t element (%i,%i) has residual norm %e\n",j,k,norms[j*(i+1)+k]);splachni(logg,buffer,4);
      }
    }
  }
  delete norms;
}

//print the upper left entry of the two matrixes
void transport_class::print_matrixes()
{
  prec Ef=parameters.read("fermi_e","meV")*units.meV/units.energy;
  for (int i=0;i<xlines;i++) {
    fprintf(logg,"E-H-S%i ",i);
    for (int j=0;j<xlines;j++) {
      content h=-HpS(i,j,0,0,0,0);
      if (i==j) h+=Ef;
      fprintf(logg,"%e%+ei ",h.real(),h.imag());
    }
    fprintf(logg,"\n");
  } 
  fprintf(logg,"\n\n");  
  for (int i=0;i<xlines;i++) {
    fprintf(logg,"S%i ",i);
    for (int j=0;j<xlines;j++) {
      content h=0;
      if (self_filled[labelM(i,j)]) {
        h=self(i,j,0,0,0,0);
        fprintf(logg,"%e%+ei ",h.real(),h.imag());
      }
      else fprintf(logg,"N.A. ");
    }
    fprintf(logg,"\n");
  } 
}



void transport_class::check_invert(content* a, content *inv, int i)
{
  int ni=dim[i]*2;
  content *aux=new content[ni*ni];
  MM(aux,a,inv,i,i,i,i);
  addC(aux,-1.0,i,i);
  sprintf(buffer,"inverted matrix checked: residual norm:%e \n",normM(aux,i,i));
  splachni(logg,buffer,4);
  delete aux;
}

//invert matrix of dimension n by n
void invert(content* res,content* in,int n)
{
  gsl_matrix_complex * aux=gsl_matrix_complex_alloc(n, n);
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      gsl_complex x;
      GSL_SET_COMPLEX(&x,in[i*n+j].real(),in[i*n+j].imag());
       gsl_matrix_complex_set(aux, i, j, x);
    }
  }
  gsl_permutation* p=gsl_permutation_alloc(n);
  int signum;
  int ans=gsl_linalg_complex_LU_decomp (aux, p, &signum);
  if (ans!=0) {sprintf(buffer,"return value from LU_decomp:%i!!!\n",ans);splachni(logg,buffer,1+4);}
  gsl_matrix_complex * inverse=gsl_matrix_complex_alloc(n, n);
  ans=gsl_linalg_complex_LU_invert (aux, p, inverse);
  if (ans!=0) {sprintf(buffer,"return value from LU_invert:%i!!!\n",ans);splachni(logg,buffer,1+4);}
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      gsl_complex x=gsl_matrix_complex_get(inverse, i, j);
      res[i*n+j]=content(GSL_REAL(x),GSL_IMAG(x));
    }
  }
  gsl_permutation_free(p);
  gsl_matrix_complex_free(aux);
  gsl_matrix_complex_free(inverse);
}

//diagonalize hermitian matrix in of dimension n by n, return vectors in w, energies in e
void diag(content *in, content *w, prec* e, int n)
{
  gsl_eigen_hermv_workspace* sp=gsl_eigen_hermv_alloc(n);
  gsl_matrix_complex* A=gsl_matrix_complex_alloc(n,n);
  gsl_matrix_complex *evec=gsl_matrix_complex_alloc(n, n);
  gsl_vector* eval=gsl_vector_alloc(n);
  
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      gsl_complex x;
      GSL_SET_COMPLEX(&x,in[i*n+j].real(),in[i*n+j].imag());
       gsl_matrix_complex_set(A, i, j, x);
    }
  }
  int ans=gsl_eigen_hermv(A, eval,evec,sp);
  if (ans!=0) printf("return value from eigen_hermv:%i!\n",ans);
  for (int i=0;i<n;i++) e[i]=gsl_vector_get(eval,i);
  
  int* pos=new int[n];
  for (int i=0;i<n;i++) pos[i]=i;

  for (int i=0;i<n;i++) {
    for (int j=0;j<n-i-1;j++) {
      if (e[j]<=e[j+1]) continue;
      prec aux;
      aux=e[j];e[j]=e[j+1];e[j+1]=aux;
      int auxi;
      auxi=pos[j];pos[j]=pos[j+1];pos[j+1]=auxi;
    }
  }

  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      gsl_complex x=gsl_matrix_complex_get(evec, i, pos[j]);
      w[i+j*n]=content(GSL_REAL(x),GSL_IMAG(x));
    }
  }
  delete pos;

  gsl_vector_free(eval);
  gsl_matrix_complex_free(A);
  gsl_matrix_complex_free(evec);
  gsl_eigen_hermv_free(sp);
}

#define TRANSW1
#define TRANSW2
#define TRANSWCHECK


//constructor
transport_window_prototype::transport_window_prototype(prec E1_, prec E2_, int NP_)
{
  E1=E1_;
  E2=E2_;
  NP=NP_;
  Nl=(int) parameters.read("leads","-");
  step=(E2-E1)/(NP-1);
  
  t=new prec*[Nl*Nl*16];
  t[0]=new prec[Nl*Nl*16*NP];
  for (int i=1;i<Nl*Nl*16;i++) t[i]=t[0]+i*NP;
  N=new int*[Nl*2];
  N[0]=new int[Nl*2*NP];
  for (int i=1;i<Nl*2;i++) N[i]=N[0]+i*NP;
  mu=new prec[Nl*2];
#ifdef TRANSWCHECK
  if (E2<E1) {printf("wrong energy interval (%e,%e) in transport_window_prototype constructor\n",E1, E2);exit(1);}
  if (NP<2) {printf("too few points for interpolation (%i) in transport_window_prototype constructor\n",NP);exit(1);}
  if (Nl<1) {printf("no leads (%i) in transport_window_prototype constructor\n",Nl);exit(1);}
#endif
#ifdef TRANSW1
  fprintf(logg,"transport window class constructed:\n\tenergy window (%e,%e) [meV] split into %i points\n",E1, E2, NP);
#endif
}

//destructor
transport_window_prototype::~transport_window_prototype()
{
  delete t[0];
  delete t;
  delete N[0];
  delete N;
  delete mu;
}

void transport_window_prototype::set_spin_axis(int a)
{
  spin_axis=a;
#ifdef TRANSWCHECK
  if (a<1 || a>3) {printf("wrong spin quantization axis (%i) in transport_window_prototype::set_spin_axis\n",a);exit(1);}
#endif
}

//transform the ordered set of data into equally spaced set
void transport_window_prototype::input(int NP_, prec* e_, int** N_, prec** t_)
{
#ifdef TRANSWCHECK
  if (NP<1) {printf("too few points (%i) in transport_window_prototype::input\n",NP_);exit(1);}
#endif

  int pntr=0;
  if (e_[0]>E1 || e_[NP_-1]<E2) fprintf(logg, "warning: data range (%e,%e) in transport_window_prototype::input do not cover the whole interval (%e,%e) - putting zeros everywhere outside\n",e_[0],e_[NP_-1],E1,E2);

  for (int ie=0;ie<NP;ie++) {
    prec act=E1+step*ie;
    if (act<e_[0] || act>e_[NP-1]) {//outside the range - set all transmissions to zero
      for (int i=0;i<Nl*Nl*16;i++) t[i][ie]=0;
      for (int i=0;i<Nl;i++) {N[i*2][ie]=N[i*2+1][ie]=0;}
      continue;
    }
    if (act==e_[pntr]) {
      for (int i=0;i<Nl;i++) for (int j=0;j<Nl;j++) for (int a=0;a<4;a++) for (int b=0;b<4;b++) 
        t[(i*Nl+j)*16+a*4+b][ie]=t_[(i*Nl+j)*16+a*4+b][pntr];
      for (int i=0;i<Nl;i++) {N[i*2][ie]=N_[i*2][pntr];N[i*2+1][ie]=0;}
      continue;
    }
    while (! (e_[pntr]<act && e_[pntr+1]>act) ) {
      pntr++;
      if (pntr>NP_-2) {printf("the loop can not got here: pntr=%i, NP_=%i\n",pntr, NP_);exit(1);}
    }
    if (e_[pntr+1]<=e_[pntr]) {printf("this can not happen: e_[pntr+1]=%e, e[pntr]=%e, act=%e\n",e_[pntr+1], e_[pntr], act);exit(1);}

    bool change=false;
    for (int i=0;i<Nl;i++) {
      N[i*2][ie]=N_[i*2][pntr];
      N[i*2+1][ie]=0;
      if (N_[i*2][pntr]!=N_[i*2][pntr+1]) change=true;
    }

    for (int i=0;i<Nl;i++) for (int j=0;j<Nl;j++) for (int a=0;a<4;a++) for (int b=0;b<4;b++) {
      int label=(i*Nl+j)*16+a*4+b;
      prec x=(act-e_[pntr])/(e_[pntr+1]-e_[pntr]);
      if (change) x=0;
      t[label][ie]=t_[label][pntr]+(t_[label][pntr+1]-t_[label][pntr])*x;
    } 
  }
}

void transport_window_prototype::set_mu(int l, int z, prec mul)
{
#ifdef TRANSWCHECK
  if (z<0 || z>1 || l<1 || l>Nl) {printf("wrong input (l=%i, z=%i) in transport_window_prototype::set_mu\n",l,z);exit(1);}
#endif
  mu[l*2+z]=mul;
}

void transport_class::derivative(prec *V, prec Bstep, prec* results)
{
  prec current_unit=units.e*units.e/(2*M_PI*units.hbar)*1e-3*1e+9;
  //fprintf(logg,"current unit is %e nA /mV\n",current_unit);
  prec T[4*4*16];
  prec I[6*4];
  prec Vqpc[5];
  prec tQPCder=0;
  prec Bold=parameters.read("B","Tesla");
  prec phiold=parameters.read("phi_B","pi");
  prec aux=0;

  //prec nso_phi=M_PI/2;
  //if (parameters.read("dress","meVA")!=0) nso_phi=atan(parameters.read("br","meVA")/parameters.read("dress","meVA"));
  prec nso_phi=0;
  //if (parameters.read("br","meVA")!=0) nso_phi=atan(parameters.read("dress","meVA")/parameters.read("br","meVA"));
  //nso_phi/=M_PI;

  for (int pntr=0;pntr<3;pntr++) {
    prec B=0, phi=0;
    switch (pntr) {
      case 0 : {B=0;phi=nso_phi;break;}
      case 1 : {B=Bstep;phi=nso_phi;break;}
      case 2 : {B=Bstep;phi=nso_phi+0.5;break;}
      case 3 : {B=-Bstep;phi=nso_phi;break;}
      case 4 : {B=-Bstep;phi=nso_phi+0.5;break;}
    }
     
    parameters.set("B",B,"Tesla");
    parameters.set("phi_B",phi,"pi");

    int M[4];
    prepare_geometry(false);
    if (Nl!=4) {printf("wrong number of leads (%i) in transport_class::derivative\n",Nl);exit(1);}
    prepare_with_E();
    recursive_green(false);
    scattering_matrix(M);
    //check_unitarity();
    //check_trs();
    for (int i=1;i<Nl;i++) for (int j=1;j<Nl;j++) for (int a=0;a<4;a++) for (int b=0;b<4;b++) {
      T[(i*Nl+j)*16+a*4+b]=lead2lead_spinresolved(i,a,j,b);
      //fprintf(logg,"T_{%i%i}^{%i%i}=%e\n",i,j,a,b,T[(i*Nl+j)*16+a*4+b]);
    }

    Vqpc[pntr]=(T[(3*Nl+1)*16+0*4+0]*V[1]+T[(3*Nl+2)*16+0*4+0]*V[2])/(M[3]-T[(3*Nl+3)*16+0*4+0]);
    fprintf(logg, "QPC potential is %e [mV] (and V[1]=%e, V[2]=%e, T_31^00=%e, T_32^00=%e, M[3]=%i)\n",Vqpc[pntr], V[1], V[2], T[(3*Nl+1)*16+0*4+0], T[(3*Nl+2)*16+0*4+0],M[3]);

    //set QPC current to zero
    if (pntr==0) {V[3]=Vqpc[pntr];aux=M[3]-T[(3*Nl+3)*16+0*4+0];}

    //compute the spin and charge current
    for (int a=0;a<4;a++) {
      I[pntr*4+a]=0;
      if (a==0) I[pntr*4+a]=M[3]*V[3];
      for (int i=1;i<4;i++) I[pntr*4+a]-=T[(3*Nl+i)*16+a*4+0]*V[i];
    }
    for (int a=0;a<4;a++) I[pntr*4+a]*=current_unit;
    fprintf(logg, "currents through QPC: I[0]=%e, I[1]=%e, I[2]=%e, I[3]=%e\n", I[pntr*4+0], I[pntr*4+1], I[pntr*4+2], I[pntr*4+3]);

    clean_with_E();
    clean_geometry();
  }

  results[0]=(I[1*4+0]-I[0*4+0])/(1*Bstep);
  results[1]=(I[2*4+0]-I[0*4+0])/(1*Bstep);
  results[2]=cos(nso_phi*M_PI)*I[0*4+1]+sin(nso_phi*M_PI)*I[0*4+2];
  results[3]=-sin(nso_phi*M_PI)*I[0*4+1]+cos(nso_phi*M_PI)*I[0*4+2];
  results[4]=I[0*4+3];
  results[5]=Vqpc[0];
  results[6]=(Vqpc[1]-Vqpc[0])/(1*Bstep);
  results[7]=(Vqpc[2]-Vqpc[0])/(1*Bstep);

  parameters.set("B",Bold,"Tesla");
  parameters.set("phi_B",phiold,"pi");
  fprintf(logg, "QPC voltage derivatives are (1:%e,2:%e), with the part due to current derivative only (1:%e, 2:%e)\n", results[6], results[7], results[0]/current_unit/aux, results[1]/current_unit/aux);

}

prec transport_class::QPC_derivative(prec Estep, prec& T)
{
  extern prec GF_check_Eb;
  prec fermi_e=parameters.read("fermi_e","meV");
  int M[4];
  prepare_geometry(false);
  prepare_with_E();
  recursive_green(false);
  scattering_matrix(M);
  T=M[3]-lead2lead_spinresolved(3,0,3,0);
  clean_with_E();

  parameters.set("fermi_e",fermi_e+Estep,"meV");
  prepare_with_E();
  recursive_green(false);
  scattering_matrix(M);
  prec der=(M[3]-lead2lead_spinresolved(3,0,3,0)-T)/Estep;
  clean_with_E();
  clean_geometry();
  parameters.set("fermi_e",fermi_e,"meV");
  sprintf(buffer, "computing energy derivative at GF=%e and fermi_e=%e giving der=%e at T=%e\n",GF_check_Eb,fermi_e, der, T);
  splachni(logg,buffer,1+4);
  return(der);
}


prec transport_class::set_QPC_sensitive(prec Estep, prec Vstep1, prec threshold1, int& maxstep1, prec Vstep2, prec threshold2, int& maxstep2, prec estimate, prec& T_res)
{

  extern prec GF_check_Eb;
  prec der, prev, T;
  int pntr1=0, pntr2=0;
  int M[4];

  //first - make sure the QPC is closed
  do {
    der=QPC_derivative(Estep,T);
    pntr1++;
    sprintf(buffer, "QPC sensitive stage 1: step %i, der=%e at T=%e at QPC offset %e\n", pntr1, der, T, GF_check_Eb);
    splachni(logg, buffer, 1+4);
    GF_check_Eb+=Vstep1;
  } while ((T>threshold1 || der<0) && pntr1<maxstep1);
  GF_check_Eb-=Vstep1;

  //compute the QPC transmission derivative wrt Fermi energy
  der=0;
  do {
    pntr2++;
    prev=der;
    T_res=T;
    der=QPC_derivative(Estep,T);
    sprintf(buffer, "stage 2: step %i, QPC energy derivative at QPC offset %e is %e (previous %e, transmission %e)\n", pntr2, GF_check_Eb, der, prev, T);
    splachni(logg,buffer,1+4);
    GF_check_Eb-=Vstep2;
  } while ((der>prev || T<threshold2) && pntr2<maxstep2);
  GF_check_Eb+=2*Vstep2;

  if (pntr2==maxstep2) {//failed
    GF_check_Eb=estimate;
    prev=QPC_derivative(Estep,T_res);
  }
  maxstep1=pntr1;
  maxstep2=pntr2;
  sprintf(buffer, "stage 3: finished with GF=%e at T=%e and the derivative der=%e after %i+%i steps\n", GF_check_Eb, T_res, prev, maxstep1, maxstep2);
  splachni(logg,buffer,1+4);
  return(prev);
}

#include "main.h"

const bool vulcano_potential_ring = false; //this option is important only if ring!=0

struct lead_prototype lead[4+1];

/*
content aux_x(prec x, prec y) {return(content(x,0));}
content aux_y(prec x, prec y) {return(content(y,0));}

prec vektorx(prec x, prec y) {
  prec B_, theta_B_;
  parameters.read("B_",B_);
  parameters.read("theta_B_",theta_B_);
  return(-B_*cos(theta_B_)/2*y); 
}
 
prec vektory(prec x, prec y) {
  prec B_, theta_B_;
  parameters.read("B_",B_);
  parameters.read("theta_B_",theta_B_);
  return(parameters.B_*cos(parameters.theta_B_)/2*x); 
}
content vektorcx(prec x, prec y) { return(content(vektorx(x,y),0));}
content vektorcy(prec x, prec y) { return(content(vektory(x,y),0));}

prec vektor2(prec x,prec y) {
  prec B_, theta_B_;
  parameters.read("B_",B_);
  parameters.read("theta_B_",theta_B_);
  return(B_*B_*cos(theta_B_)*cos(theta_B_)*(x*x+y*y)/4);
}

content vektor2c(prec x,prec y) {
  prec B_, theta_B_;
  parameters.read("B_",B_);
  parameters.read("theta_B_",theta_B_);
  return(content(B_*B_*cos(theta_B_)*cos(theta_B_)*(x*x+y*y)/4,0));
}

content r2c(prec x, prec y) {return(content(x*x+y*y,0));}
*/
    
content exp_integral(prec xi, prec yi, prec xf, prec yf)
 {
   prec B=parameters.values[parameters_class::B];
   prec theta_B=parameters.values[parameters_class::theta_B];
   
   prec ring=parameters.values[parameters_class::ring];
   prec geometry=(int) parameters.values[parameters_class::geometry];
   
   prec value=0;
   
   if (geometry==1) value=B*cos(theta_B)/2*ring*(xf-xi);
   else value=B*cos(theta_B)/2*(xi*(yf-yi)-yi*(xf-xi));
   if (value==0) return(content(1,0));
   return(exp(content(0,value)));
 }

 content potentialc(prec x, prec y) { return(content(potential(x,y),0));}
 content potential_grc(prec x, prec y) { return(potential_gr(x,y));}
 content mass_grc(prec x, prec y) { return(mass_gr(x,y));}

 content spin_potential(prec x, prec y){  return(spin_impurities.potential(x,y)); }
 content xc(prec x, prec y) {return(content(x,0));}
 content yc(prec x, prec y) {return(content(y,0));}

 int ring_eiphi_exp;
 content ring_eiphi(prec x, prec y)
 {
   prec ring=parameters.values[parameters_class::ring];
   return(exp(content(0,x/ring*ring_eiphi_exp)));
 }

int soc_modulation_aux;//0-modulation function, 1-partial x derivative, 2-partial y derivative
content soc_modulation(prec x, prec y)
{
  int leads=(int) parameters.values[parameters_class::leads];
  prec pinned=parameters.values[parameters_class::flat_x]*0;
  pinned+=100000*units.nm/units.length;
  prec ramp=0*units.nm/units.length;
  int region=0;
  if (x<pinned-ramp) region=-1;
  else if (x>pinned) region=+1;
  switch (soc_modulation_aux) {
    case 0 : {
      if (region==-1) return(1.0);
      if (region==+1) return(0.0);
      return((pinned-x)/ramp);
      break;
    }
    case 1 : {
      if (region==0) return(-1/ramp);
      return(0.0);
      break;
    }
    case 2 : {
      return(0.0);
      break;
    }
    default : {
      printf("wrong control flag (%i) in soc_modulation\n",soc_modulation_aux);
      exit(1);
    }
  }
}
  
 /*

content h_00(prec x, prec y)
{ return(content(vektor2(x,y)+potential(x,y)+parameters.azz*parameters.B_*cos(parameters.theta_B_),0));}
content h_01(prec x, prec y)
{ return(content(parameters.ab*vektory(x,y)-parameters.ad*vektorx(x,y), parameters.ab*vektorx(x,y)-parameters.ad*vektory(x,y))); }
content h_10(prec x, prec y)
{ return(content(parameters.ab*vektory(x,y)-parameters.ad*vektorx(x,y), -parameters.ab*vektorx(x,y)+parameters.ad*vektory(x,y))); }
content h_11(prec x, prec y)
{ return(content(vektor2(x,y)+potential(x,y)-parameters.azz*parameters.B_*cos(parameters.theta_B_),0)); }

content hd3_01_1(prec x, prec y) {
return(content(-4*vektory(x,y),-4*vektorx(x,y)));
}

content hd3_01_2(prec x, prec y) {
return(content(4*vektory(x,y)*vektorx(x,y),-2*vektory(x,y)*vektory(x,y)));
}

content hd3_01_3(prec x, prec y) {
return(content(2*vektorx(x,y)*vektorx(x,y),-4*vektorx(x,y)*vektory(x,y)));
}

content hd3_01_4(prec x, prec y) {
return(content(2*vektorx(x,y)*vektory(x,y)*vektory(x,y),2*vektorx(x,y)*vektorx(x,y)*vektory(x,y)));
}

content hd3_10_1(prec x, prec y) {
return(content(-4*vektory(x,y),4*vektorx(x,y)));
}

content hd3_10_2(prec x, prec y) {
return(content(-4*vektory(x,y)*vektorx(x,y),-2*vektory(x,y)*vektory(x,y)));
}

content hd3_10_3(prec x, prec y) {
return(content(-2*vektorx(x,y)*vektorx(x,y),-4*vektorx(x,y)*vektory(x,y)));
}

content hd3_10_4(prec x, prec y) {
return(content(2*vektorx(x,y)*vektory(x,y)*vektory(x,y),-2*vektorx(x,y)*vektorx(x,y)*vektory(x,y)));
}
*/

void hamiltonian_for_lead_init(reg1& region, int heading)
{
   int precision=(int) parameters.read("precision","-");

   prec m_eff=parameters.read("m_eff","m_e");
   prec B=parameters.read("B","Tesla");
   prec theta_B=parameters.read("theta_B","pi")*M_PI;
   prec phi_B=parameters.read("phi_B","pi")*M_PI;
   prec gx=parameters.read("gx","-");
   prec gy=parameters.read("gy","-");
   prec gz=parameters.read("gz","-");
   /*prec br=parameters.read("br","meVA")*units.meV*1e-10;
   prec dress=parameters.read("dress","meVA")*units.meV*1e-10;
   prec dress3=parameters.read("dress3","eVA^3")*units.e*1e-30;
   int geometry=(int) parameters.read("geometry","-");
   int growth=(int) parameters.read("growth","-");
   if (growth!=1 && growth!=110) {
     sprintf(buffer,"wrong growth (%i) in hamiltonian_init\nchoosing growth 1\n",growth);
     splachni(logg,buffer,1+4);
     growth=1;
   }*/
  region.act_prec_set((reg1::operators_precision) precision);
  region.symmetrize_ops((bool) parameters.read("peierls","-"), true);
  
  bool op_rmmbr[reg1::all_op]={false};
  for (int i=0;i<reg1::all_op;i++) region.remember_set((reg1::operators) i,op_rmmbr[i]);

  //only kinetic energy along y direction -- no magnetic field!!!!
  content u=pow(units.hbar,2)/(2*units.m_e*m_eff)/pow(units.length,2)/units.energy; //units
  //parameters.set("B",0,"Tesla");
  region.operate_tot_init();
  if (heading==1 || heading==2) {
    region.operate_tot_init(0,0,reg1::dydy,-u,0,exp_integral);
    region.operate_tot_init(1,1,reg1::dydy,-u,0,exp_integral);
  }
  if (heading==3 || heading==4) {
    region.operate_tot_init(0,0,reg1::dxdx,-u,0,exp_integral);
    region.operate_tot_init(1,1,reg1::dxdx,-u,0,exp_integral);
  }
  //parameters.set("B",B,"Tesla");

  
  u=units.hbar*units.e/(4*units.m_e)/units.energy;
  //Zeeman "z"
  prec Bz=B*cos(theta_B);
  region.operate_tot_init(0,0,reg1::lin_com,Bz*gz*u,0,exp_integral);//here was aux_Z
  region.operate_tot_init(1,1,reg1::lin_com,-Bz*gz*u,0,exp_integral);//here was aux_Z
  
  //Zeeman "x, y"
  prec Bx=B*sin(theta_B)*cos(phi_B);
  prec By=B*sin(theta_B)*sin(phi_B);
  region.operate_tot_init(1,0,reg1::lin_com,content(gx*Bx,-gy*By)*u,0,exp_integral);
  region.operate_tot_init(0,1,reg1::lin_com,content(gx*Bx,gy*By)*u,0,exp_integral);

}

void hamiltonian_init_kinetic_potential(reg1& region, int sets)
{
  int geometry=(int) parameters.read("geometry","-");
  //!potential
  for (int i=0;i<sets;i++) region.operate_tot_init(i,i,reg1::lin_com, content(1,0),potentialc,exp_integral);
  
  //!kinetic energy (if 1D problem, only in x variable)
  prec m_eff=parameters.read("m_eff","m_e");
  double u=pow(units.hbar,2)/(2*units.m_e*m_eff)/pow(units.length,2)/units.energy; //units
  for (int i=0;i<sets;i++) {
    region.operate_tot_init(i,i,reg1::dxdx,-u,0,exp_integral);
    if (geometry!=1) region.operate_tot_init(i,i,reg1::dydy,-u,0,exp_integral);
  }
}

void hamiltonian_init_spin_impurities(reg1& region, int sets)
{
  //!spin impurities
  prec Arho=-parameters.read("exch","meV")*parameters.read("xMn","-")*(5.0/2)*(1.0/2);
  if (Arho!=0 && spin_impurities.N!=0) {
    if (dof_prototype::spin<dof.Ndof) {
      printf("spin impurities are present, but the particle has no spin (%i vs Ndof %i)\n",dof_prototype::spin, dof.Ndof);
      exit(1);
    }
    dof.reset("indexes");
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
      spin_impurities.set_indexes(i,j);
      dof.indexC[dof_prototype::spin]=1-2*i;
      dof.convert("C2I");
      int I=dof.index;
      dof.indexC[dof_prototype::spin]=1-2*j;
      dof.convert("C2I");
      int J=dof.index;
      region.operate_tot_init(I,J,reg1::lin_com, content(1,0), spin_potential, exp_integral);
    }
  }
}

void hamiltonian_init_Zeeman(reg1& region)
{
  prec B=parameters.read("B","Tesla");
  prec theta_B=parameters.read("theta_B","pi")*M_PI;
  prec phi_B=parameters.read("phi_B","pi")*M_PI;
  prec gx=parameters.read("gx","-");
  prec gy=parameters.read("gy","-");
  prec gz=parameters.read("gz","-");
  content u=units.hbar*units.e/(4*units.m_e)/units.energy;
  
  dof.reset("all");
  dof.operatorC[dof_prototype::spin]=dof_prototype::paulix;
  dof.convert("OP2M",units.muB*parameters.read("gx","-")/2*B*sin(theta_B)*cos(phi_B)/units.energy);
  dof.operatorC[dof_prototype::spin]=dof_prototype::pauliy;
  dof.convert("OP2M",units.muB*parameters.read("gy","-")/2*B*sin(theta_B)*sin(phi_B)/units.energy);
  dof.operatorC[dof_prototype::spin]=dof_prototype::pauliz;
  dof.convert("OP2M",units.muB*parameters.read("gz","-")/2*B*cos(theta_B)/units.energy);
  for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++)
    if (abs(dof.M_sum[i*dof.length+j])!=0) region.operate_tot_init(j,i,reg1::lin_com,dof.M_sum[i*dof.length+j],0,exp_integral);
}

void hamiltonian_init_SOI_2D_100_lin(reg1& region, double br, double dress, double u)
{
  
  dof.reset("all");
  dof.operatorC[dof_prototype::spin]=dof_prototype::paulix;
  dof.convert("OP2M",-dress*content(0,-1)*u);
  dof.operatorC[dof_prototype::spin]=dof_prototype::pauliy;
  dof.convert("OP2M",-br*content(0,-1)*u);
  for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++)
    if (abs(dof.M_sum[i*dof.length+j])!=0) region.operate_tot_init(j,i,reg1::dx,dof.M_sum[i*dof.length+j],0,exp_integral);
  
  dof.reset("all");
  dof.operatorC[dof_prototype::spin]=dof_prototype::paulix;
  dof.convert("OP2M",br*content(0,-1)*u);
  dof.operatorC[dof_prototype::spin]=dof_prototype::pauliy;
  dof.convert("OP2M",dress*content(0,-1)*u);
  for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++)
    if (abs(dof.M_sum[i*dof.length+j])!=0) region.operate_tot_init(j,i,reg1::dy,dof.M_sum[i*dof.length+j],0,exp_integral);
  
  
  //!spatially dependent spin-orbit interaction strengths
  soc_modulation_aux=0;
  /*region.operate_tot_init(1,0,reg1::dy, content(-dress,-br)*u,soc_modulation,exp_integral);
   region.operate_tot_init(0,1,reg1::dy, content(dress,-br)*u,soc_modulation,exp_integral);
   region.operate_tot_init(1,0,reg1::dx, content(br,dress)*u,soc_modulation,exp_integral);
   region.operate_tot_init(0,1,reg1::dx, content(-br,dress)*u,soc_modulation,exp_integral);*/
  /*soc_modulation_aux=2;
   region.operate_tot_init(1,0,reg1::lin_com, content(-dress,-br)*u/2.0,soc_modulation,exp_integral);
   region.operate_tot_init(0,1,reg1::lin_com, content(dress,-br)*u/2.0,soc_modulation,exp_integral);
   soc_modulation_aux=1;
   region.operate_tot_init(1,0,reg1::lin_com, content(br,dress)*u/2.0,soc_modulation,exp_integral);
   region.operate_tot_init(0,1,reg1::lin_com, content(-br,dress)*u/2.0,soc_modulation,exp_integral);*/
}

void hamiltonian_init_SOI_2D_110_lin(reg1& region, double br, double dress, double u)
{
  region.operate_tot_init(1,0,reg1::dy, content(0,-br)*u,0,exp_integral);
  region.operate_tot_init(1,0,reg1::dx, content(br,0)*u,0,exp_integral);
  region.operate_tot_init(0,1,reg1::dy, content(0,-br)*u,0,exp_integral);
  region.operate_tot_init(0,1,reg1::dx, content(-br,0)*u,0,exp_integral);
  region.operate_tot_init(0,0,reg1::dx, content(0,-dress/2)*u,0,exp_integral);
  region.operate_tot_init(1,1,reg1::dx, content(0,dress/2)*u,0,exp_integral);
}

void hamiltonian_init_SOI_2D_100_cubic(reg1& region, double dress3, double u)
{
  dof.reset("all");
  dof.operatorC[dof_prototype::spin]=dof_prototype::paulix;
  dof.convert("OP2M",dress3*content(0,1)*u);
  for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++)
    if (abs(dof.M_sum[i*dof.length+j])!=0) region.operate_tot_init(j,i,reg1::dxdydy,dof.M_sum[i*dof.length+j],0,exp_integral);
  
  dof.reset("all");
  dof.operatorC[dof_prototype::spin]=dof_prototype::pauliy;
  dof.convert("OP2M",-dress3*content(0,1)*u);
  for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++)
    if (abs(dof.M_sum[i*dof.length+j])!=0) region.operate_tot_init(j,i,reg1::dxdxdy,dof.M_sum[i*dof.length+j],0,exp_integral);
}

void hamiltonian_init_SOI_2D_110_cubic(reg1& region, double dress3, double u)
{
  region.operate_tot_init(0,0,reg1::dxdxdx, content(dress3/2,0)*u,0,exp_integral);
  region.operate_tot_init(1,1,reg1::dxdxdx, content(-dress3/2,0)*u,0,exp_integral);
  
  region.operate_tot_init(0,0,reg1::dxdydy, content(-dress3,0)*u,0,exp_integral);
  region.operate_tot_init(1,1,reg1::dxdydy, content(dress3,0)*u,0,exp_integral);
}

void hamiltonian_init_SOI_1D_100_lin(reg1& region, double br, double dress, double ring)
{
  //the coupling strength is (all tilda) hbar^2/(4 m l_d R)
  //while the dimensionless constant dress corresponds to (all tilda) hbar^2/2 m l_d l_0
  //therefore to convert to ring means (dress/2R) with R dimensionless (Rtilda/l_0)
  //the differential operator \partial_\phi = (L/2\pi) \partial_x = R \partial_x
  //sigma_\pm = (sigma_x \pm i\sigma_y)/2 => sigma_+ = (0,1), sigma_-=(1,0)
  double ring_eiphi_exp=-1;
  content u=dress*content(0,1)/(2.0*ring)/units.energy;
  region.operate_tot_init(1,0,reg1::dx,content(0,-ring)/units.length*u,ring_eiphi,exp_integral);
  region.operate_tot_init(1,0,reg1::lin_com,u/(-2.0),ring_eiphi,exp_integral);
  ring_eiphi_exp=1;
  region.operate_tot_init(0,1,reg1::dx,content(0,ring)/units.length*u,ring_eiphi,exp_integral);
  region.operate_tot_init(0,1,reg1::lin_com,u/(-2.0),ring_eiphi,exp_integral);
  ring_eiphi_exp=-1;
  u=br*content(1,0)/(2.0*ring)/units.energy;
  region.operate_tot_init(0,1,reg1::dx,content(0,-ring)/units.length*u,ring_eiphi,exp_integral);
  region.operate_tot_init(0,1,reg1::lin_com,u/(-2.0),ring_eiphi,exp_integral);
  ring_eiphi_exp=1;
  region.operate_tot_init(1,0,reg1::dx,content(0,-ring)/units.length*u,ring_eiphi,exp_integral);
  region.operate_tot_init(1,0,reg1::lin_com,u/(2.0),ring_eiphi,exp_integral);
}

void hamiltonian_init_SOI_1D_100_cubic(reg1& region, double dress3, double ring, double lv)
{
  //the coupling strength is (all tilda) (gamma_c/2) i/16w^2R
  //while the dimensionless constant dress3 corresponds to (all tilda) (gamma_c/2)/l_0^3
  //therefore to convert to ring means (i dress3/16Rw^2) with R,w dimensionless
  //parameter lv is used as the ring transversal width
  content u=dress3*content(0,1)/(16.0*ring*lv*lv);
  ring_eiphi_exp=1;
  region.operate_tot_init(0,1,reg1::dx,content(0,-ring)/units.length*u,ring_eiphi,exp_integral);
  region.operate_tot_init(0,1,reg1::lin_com,u/2.0,ring_eiphi,exp_integral);
  ring_eiphi_exp=-1;
  region.operate_tot_init(1,0,reg1::dx,content(0,ring)/units.length*u,ring_eiphi,exp_integral);
  region.operate_tot_init(1,0,reg1::lin_com,u/2.0,ring_eiphi,exp_integral);
  ring_eiphi_exp=-3;
  region.operate_tot_init(0,1,reg1::dx,content(0,-3*ring)/units.length*u,ring_eiphi,exp_integral);
  region.operate_tot_init(0,1,reg1::lin_com,-9.0/2.0*u,ring_eiphi,exp_integral);
  ring_eiphi_exp=3;
  region.operate_tot_init(1,0,reg1::dx,content(0,3*ring)/units.length*u,ring_eiphi,exp_integral);
  region.operate_tot_init(1,0,reg1::lin_com,-9.0/2.0*u,ring_eiphi,exp_integral);
}

void  hamiltonian_init_SOI(reg1& region)
{
  int geometry=(int) parameters.read("geometry","-");
  int sets=(int) parameters.read("sets","-");
  int growth=(int) parameters.read("growth","-");
  if (growth!=1 && growth!=110) {
    sprintf(buffer,"wrong growth (%i) in hamiltonian_init_SOI\nchoosing growth 1\n",growth);
    splachni(logg,buffer,1+4);
    growth=1;
  }
  
  prec br=parameters.read("br","meVA")*units.meV*1e-10;
  prec dress=parameters.read("dress","meVA")*units.meV*1e-10;
  prec dress3=parameters.read("dress3","eVA^3")*units.e*1e-30;
  prec u=1/units.length/units.energy;
  if (geometry!=1 && dof_prototype::spin<dof.Ndof) {//2D problem
    if (br!=0 || dress!=0) {
      if (growth==1) hamiltonian_init_SOI_2D_100_lin(region, br, dress, u);
      if (growth==110) {
        if (sets!=2) {
          fprintf(logg,"growth direction 110 not implemented for dof class and linear SOC\n");
          exit(1);
        }
        hamiltonian_init_SOI_2D_110_lin(region, br, dress, u);
      }
    }
    if (dress3!=0) {
      u=1/pow(units.length,3)/units.energy;
      if (growth==1) hamiltonian_init_SOI_2D_100_cubic(region, dress3, u);
      if (growth==110) {
        if (sets!=2) {
          fprintf(logg,"growth direction 110 not implemented for dof class and cubic SOC\n");
          exit(1);
        }
        hamiltonian_init_SOI_2D_110_cubic(region, dress3, u);
      }
    }
  }
  if (geometry==1 && dof_prototype::spin<dof.Ndof) { //1D ring
    prec ring=parameters.read("ring","nm")*units.nm;
    prec lv=parameters.read("lv","nm")*units.nm;
    
    if (br!=0 || dress!=0) {
      if (growth==1) {
        if (sets!=2) {
          fprintf(logg,"linear SOC not implemented for 1D ring and dof class\n");
          exit(1);
        }
        hamiltonian_init_SOI_1D_100_lin(region, br, dress, ring);
      }
      if (growth==110) {
        //region.operate_tot_init(1,0,reg1::dy, content(0,-br),0,exp_integral);
        printf("spin-orbit D&BR along 110 not implemented in 1d ring\n");
        exit(1);
      }
    }
    if (dress3!=0) { //Dresselhaus cubic
      if (growth==1) {
        if (sets!=2) {
          fprintf(logg,"linear SOC not implemented for 1D ring and dof class\n");
          exit(1);
        }
        hamiltonian_init_SOI_1D_100_cubic(region, dress3, ring, lv);
      }
      if (growth==110) {
        //region.operate_tot_init(0,0,reg1::dxdxdx, content(2*dress3/2,0),0,exp_integral);
        printf("spin-orbit D3 along 110 not implemented in 1d ring\n");exit(1);
      }
    }
  }
}

//z-energies
double Ezi(int i) {
    return(i+0.5);
}

//matrix elements of the z operator
complex<double> zij(int i, int j) {
    if (i==j-1) return(sqrt(j/2.0));
    if (i==j+1) return(sqrt(i/2.0));
    return(0);
}

void  hamiltonian_init_finite_width(reg1& region)
{
  if (dof.length<2) {printf("dof.length=%i too small in hamiltonian_init_finite_width\n",dof.length); exit(1);}
  
  int length=dof.length/2;
  complex<double>* M_in=new complex<double>[length*length];
  complex<double>* M_out=new complex<double>[dof.length*dof.length];
  
  //the z-coordinate dependent part
  double m=parameters.read("m_eff","m_e")*units.m_e;
  double lz=parameters.read("width_z","nm")*units.nm;
  double B=parameters.read("B","Tesla");
  double theta_B=parameters.read("theta_B","pi")*M_PI;
  double phi_B=parameters.read("phi_B","pi")*M_PI;
  double Bin=B*sin(theta_B), Bz=B*cos(theta_B);
  double Bx=Bin*cos(phi_B), By=Bin*sin(phi_B);
  double lzB=pow( pow(lz,-4) + pow(units.e*Bin/units.hbar,2), -0.25);

  //the diagonal energies
  double E=units.hbar*units.hbar/(m*lzB*lzB)/units.energy;
  for (int i=0;i<dof.length;i++) region.operate_tot_init(i,i,reg1::lin_com,E*Ezi(i/2),0,0);
  
  //the matrix of z-coordinate elements
  for (int i=0;i<length;i++) for (int j=0;j<length;j++) {
    M_in[i*length+j]=zij(i,j);
    fprintf(logg,"z-coordinate mat-element [%i,%i]=%e%+ei\n",i,j,M_in[i*length+j].real(),M_in[i*length+j].imag());
  }
  //the spin-independent terms
  
  //expand by the spin dof
  dof.add_one_tensor(dof_prototype::unity, M_in, dof.Ndof-1, M_out);
  for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++) {
    fprintf(logg,"expanded z-coordinate mat-element [%i,%i]=%e%+ei\n",i,j,M_out[i*dof.length+j].real(), M_out[i*dof.length+j].imag());
  }
  
  //momentum terms
  complex<double> ii=complex<double>(0,1.0);
  E=(lzB*units.e*units.hbar/m/units.length)/units.energy;
  fprintf(logg,"scale of momentum dependent z-terms: %e\n",E*units.energy/units.meV);
  for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++)
    if (abs(M_out[i*dof.length+j])!=0) region.operate_tot_init(j,i,reg1::dx,-ii*By*E*M_out[i*dof.length+j],0,0);
  for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++)
    if (abs(M_out[i*dof.length+j])!=0) region.operate_tot_init(j,i,reg1::dy,ii*Bx*E*M_out[i*dof.length+j],0,0);
  //confinement terms
  E=(lzB*units.e*units.e/(4*m)*Bz*units.length)/units.energy;
  fprintf(logg,"scale of position dependent z-terms: %e\n",E*units.energy/units.meV);
  for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++)
    if (abs(M_out[i*dof.length+j])!=0) region.operate_tot_init(j,i,reg1::lin_com,-Bx*E*M_out[i*dof.length+j],xc,0);
  for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++)
    if (abs(M_out[i*dof.length+j])!=0) region.operate_tot_init(j,i,reg1::lin_com,-By*E*M_out[i*dof.length+j],yc,0);

  delete[] M_in;
  delete[] M_out;
}


void hamiltonian_init(reg1& region)
{
   int sets=(int) parameters.read("sets","-");
   if (sets!=dof.length) fprintf(logg,"sets=%i and dof.length=%i do not match (may be intentional for 2e diagonalization)\n",sets, dof.length);

  region.act_prec_set((reg1::operators_precision) parameters.read("precision","-"));
  region.symmetrize_ops((bool) parameters.read("peierls","-"), true);
  
  bool op_rmmbr[reg1::all_op]={false};
  for (int i=0;i<reg1::all_op;i++) region.remember_set((reg1::operators) i,op_rmmbr[i]);
  region.operate_tot_init();

  //kinetic+potential energies
  hamiltonian_init_kinetic_potential(region, sets);
  //nuclear spins
  if (dof_prototype::spin<dof.Ndof) hamiltonian_init_spin_impurities(region, sets);
  //Zeeman
  if (dof_prototype::spin<dof.Ndof) hamiltonian_init_Zeeman(region);
  //SOI
  if (dof_prototype::spin<dof.Ndof) hamiltonian_init_SOI(region);
  //z dependent part
  region.symmetrize_ops(false, true);
  if (dof.Ndof-dof_prototype::spin>1) hamiltonian_init_finite_width(region);
}

void sort_0(int N,int eig, ret_from_arpack_zn* vysl, char method, char order)
{
  int* index=new int[eig];
  picsrt_index(eig, vysl->eigenvals, index,0, method, order);
  
  complex<double>* vecs_new=new complex<double>[N*eig];
  complex<double>* vals_new=new complex<double>[eig+1];
  for (int i=0;i<eig;i++) {
    vals_new[i]=vysl->eigenvals[index[i]];
    for (int j=0;j<N;j++) vecs_new[i*N+j]=vysl->eigenvecs[index[i]*N+j];
  }
  delete[] index;
  delete vysl->eigenvecs;
  delete vysl->eigenvals;
  vysl->eigenvecs=vecs_new;
  vysl->eigenvals=vals_new;
}

void normalize(int N, int eig, complex<double>* vecs)
{
  for (int i=0;i<eig;i++) {
    double norm=0;
    complex<double>* offset=vecs+i*N;
    for (int k=0;k<N;k++) norm+=(offset[k]*conj(offset[k])).real();
    norm=1/sqrt(norm);
    for (int k=0;k<N;k++) offset[k]*=norm;
  }
}



/*
content find_phase(reg1& r1, ret_from_arpack_zn& vysl,int eig, int wch, prec lev)
{
  content *input=vysl.eigenvecs+eig*r1.sets*r1.Nw;
  content trial=0,s=0,s2=0;
  int i,j,n=0,nsucc=0,wch2,k;
  prec x,y;
  bool ins;
  int where=0;
  prec d_0=parameters.values[parameters_class::d];
  prec length_x=parameters.values[parameters_class::length_x];
  

  for (k=1;k<5;k++) {
   x=-d_0+length_x/5*k;
   r1.nearest_LD_point(x,0,i,j,ins);
   int label=r1.ac2sl(i,j);
   if (! ins) continue;
   for (wch2=0;wch2<2;wch2++) {
     if (wch2!=wch && wch!=3) continue;
     trial=input[label+r1.Nw*wch2];
     if (abs(trial)<lev) continue;
     trial/=abs(trial);
     k=100;
     break;
   }
   if (k==100) break;
  }
  if (k<100) trial=content(1,0);
#ifdef PHASE
  sprintf(buffer,"k,x,phase: %i,%e,%e%+ei\n",k,x,trial.real(),trial.imag());
  splachni(logg,buffer,4);
#endif
  return(trial);
}

void fit_the_phase(ret_from_arpack_zn& vysl,reg1 &r1, int eig, int wch, prec lev)
{
  for (int i=0;i<eig;i++) {
    content phase=find_phase(r1, vysl ,i,wch,lev);
    content *pntr=vysl.eigenvecs+i*r1.Nw*r1.sets;
    for (int j=0;j<r1.Nw*r1.sets;j++) pntr[j]*=conj(phase);
  }
}
*/
prec kappa3,field2,grid;

void hmt(int dim,void *A, content* x,content *y){
  prec& a=grid;
  prec pot;
  for (int i=0;i<dim;i++) {
    if (a*i<0.1) pot=1000; else pot=kappa3*i*a+field2*i*i*a*a;
/*    y[i]=(pot+15/(6*a*a))*x[i];
    if (i>0) y[i]-=4.0*x[i-1]/(3*a*a);
    if (i<dim-1) y[i]-=4.0*x[i+1]/(3*a*a);
    if (i>1) y[i]+=x[i-2]/(12*a*a);
    if (i<dim-2) y[i]+=x[i+2]/(12*a*a);*/
    y[i]=(pot+2/(a*a))*x[i];
    if (i>0) y[i]-=x[i-1]/(a*a);
    if (i<dim-1) y[i]-=x[i+1]/(a*a);
  }
}

void dresselhaus(prec kappa)
{
  prec length=5;
  int dim=200;
  int eig=2;
  
  //define constants
  prec B=parameters.values[parameters_class::B];
  prec theta_B=parameters.values[parameters_class::theta_B];
  grid=length/(dim-1);
  kappa3=kappa*kappa*kappa;
  field2=B*B*sin(theta_B)*sin(theta_B);
  
  //diagonalize
  sprintf(buffer,"diagonalizing for dresselhaus:\n");
  splachni(logg,buffer,1+2);
  ret_from_arpack_zn vysl;
  vysl.eigenvals=new content[eig+1];
  vysl.eigenvecs=new content[eig*dim];
  arpack_zn arp1(dim,eig,1e-10,10000,true,(char*)"SM");
  arp1.go_for_it(dim,hmt,(void*) 0,&vysl,true);
  sort_0(dim,eig,&vysl,true,1);
  
  //FILE* file=fopen("data/moment","a");
  /*for (int i=0;i<dim;i++) {
    //length and potential
    fprintf(file,"%e\t%e\t",i*grid,kappa3*i*grid+field2*i*i*grid*grid);
    for (int j=0;j<eig;j++) {
      //eigenvecs and eigenvals
      fprintf(file,"%e\t%e\t",vysl.eigenvecs[i+j*dim].real(),vysl.eigenvals[j].real());
    }
    fprintf(file,"\n");
  }
  fprintf(file,"\n");*/
  
    
  //find mean momentum squared and length
  content mom2=0,l=0;
  for (int i=0;i<dim;i++) {
    vysl.eigenvecs[i+dim]=2.0*vysl.eigenvecs[i];
    if (i>0) vysl.eigenvecs[i+dim]-=vysl.eigenvecs[i-1];
    if (i<dim-1) vysl.eigenvecs[i+dim]-=vysl.eigenvecs[i+1];
    mom2+=conj(vysl.eigenvecs[i])*vysl.eigenvecs[i+dim]/(grid*grid);
    l+=conj(vysl.eigenvecs[i])*vysl.eigenvecs[i]*(grid*i);
  }
  
  prec dress=(3.1e-3*mom2).real();
  prec width_z=l.real();
  sprintf(buffer,"with result:%e%+ei at mean length %e\n",mom2.real()*3.1e-3,mom2.imag()*3.1e-3,l.real());
  splachni(logg,buffer,1+2);
  
  //fprintf(file,"%e %e %e %e %e %e\n",kappa,B_inplane,mom2.real(),l.real(),ad,vysl.eigenvals[0].real());
  //fclose(file);
  parameters.values[parameters_class::width_z]=width_z;
  parameters.values[parameters_class::dress]=dress;
  delete vysl.eigenvecs;
  delete vysl.eigenvals;

}


/*
bool inside(prec x, prec y)
{
  if (((x-d_0x)*(x-d_0x)/(l_0B*l_0B*length_0_x*length_0_x)+ (y-d_0y)*(y-d_0y)/(l_0B*l_0B*length_0_y*length_0_y))<1) return(true);
  if (((x+d_0x)*(x+d_0x)/(l_0B*l_0B*length_0_x*length_0_x)+ (y+d_0y)*(y+d_0y)/(l_0B*l_0B*length_0_y*length_0_y))<1) return(true);
  //if ((x>-d_0) && (x<d_0) && (y>-length_0_y*l_0B) && (y<length_0_y*l_0B)) return(true);
  return(false);
}*/
/*
//RING!!!
bool inside(prec x, prec y)
{
  prec length_x=parameters.values[parameters_class::length_x];
  prec length_y=parameters.values[parameters_class::length_y];
  if ((x*x/(length_x*length_x)+ y*y/(length_y*length_y))<1) return(true);
  return(false);
}
 */
 //************quartic potential*****************
 /*
 prec potential(prec x, prec y) {
   prec xx=0;
   if (x>0 && y>0) xx=x*y*0;
   x=x*x-d_0*d_0;
   return(av*(y*y+x*x/(4*d_0*d_0)+xx));
 }
*/


//**********exponential potential***************
/*
 prec potential(prec x, prec y) {
   prec vo=8*av;
   prec al=av;
   if (x<0 || Asym_) x+=d_0; else x-=d_0;
   return(vo*(1-exp(-al*(y*y+x*x))));
 }
*/

prec GF_check_Eb, GF_check_w, GF_check_wp;
bool GF_check_flag=false;
//potential for a series of step barriers -- uncomment line in potential(x,y) in case geometry 2
content potential_flat_GF_check(prec x, prec y) 
{ x+=0.25*units.nm/units.length;
  prec ll=lead[1].x;
  x-=ll;
  if (x<GF_check_w/2) return(content(0,0));
  x-=GF_check_w/2;
  int box=(int) floor(x/(GF_check_w+GF_check_wp));
  x-=box*(GF_check_w+GF_check_wp);
  if (x<GF_check_wp) return(content(GF_check_Eb,0));
  else return(content(0,0));
}

prec potential_flat(prec x, prec y) 
{ 
  prec m_eff=parameters.values[parameters_class::m_eff];
  prec d_x=parameters.values[parameters_class::d_x];
  prec d_y=parameters.values[parameters_class::d_y];
  prec d=parameters.values[parameters_class::d];

  prec prefac=pow(units.length,2)*pow(units.hbar,2)/(2*units.m_e*m_eff*pow(units.length,4));
  prefac/=units.energy;

  prec U=0;

  //if (d>0) U=100*units.meV/units.energy*exp(-(x*x+y*y)/(d*d));

  //QPC neck
  prec Etop=GF_check_Eb*units.meV/units.energy;
  prec cap=GF_check_Eb*units.meV/units.energy;
  //prec cap=7.0*units.meV/units.energy;

  //!quadratic along x and quadratic along y, limited from below: ver.2
  U=(y*y/pow(d_y,4)-(x*x)/pow(d_x,4))*prefac;
  if (U<-cap) U=-cap;
  U+=(GF_check_Eb-31.6)*units.meV/units.energy;

  //!y confinement length dependent on x
  //d_y*=(1+pow(abs(x)/d_x,2));
  //U=Etop*0+y*y/pow(d_y,4)*prefac;
  //U*=sqrt(1-x*x/lp/lp);

  int QPC=3;
  if ((abs(y-lead[QPC].y)<lead[QPC].width/2) && (x>parameters.values[parameters_class::flat_x])) {
    //inside the QPC:
    cap=10*units.meV/units.energy;
    prec l=(lead[QPC].x-parameters.values[parameters_class::flat_x])/2;
    x-=parameters.values[parameters_class::flat_x]+l;
    y-=lead[QPC].y;

    //!quadratic along x and quadratic along y, limited from below
    //U=(y*y/pow(d_y,4)-(x*x)/pow(d_x,4))*prefac;
    //if (U<-cap) U=-cap;
    //if (d_x==0) U=0;
    //U+=(GF_check_Eb*0-31.6)*units.meV/units.energy;
  }
  else {
    //lead 1,2 to the lead 3 shielding
    if (parameters.values[parameters_class::leads]==3) {
      prec r1=pow(x-(lead[1].x+min(lead[3].x,parameters.values[parameters_class::flat_x]))/2,2)+pow(y-(lead[1].y+lead[3].y)/2,2);
      prec r2=pow(x-(lead[2].x+min(lead[3].x,parameters.values[parameters_class::flat_x]))/2,2)+pow(y-(lead[2].y+lead[3].y)/2,2);
      prec r3=x*x+pow(y-lead[3].y*0,2);
      if (d>0) U=100*units.meV/units.energy*(exp(-r1*r1/(d*d))*0+exp(-r2*r2/(d*d))*0+exp(-r3*r3/(d*d)));
      //U+=31.6*units.meV/units.energy;
      //prec lx=(parameters.values[parameters_class::flat_x]-x)*units.length/units.nm;
      //prec ly=abs(y-lead[3].y)*units.length/units.nm;
      //prec fac=(1.0-min(lx,50.0)/50.0)*(1.0-min(ly,20.0)/20.0);
      //U+=-31.6*units.meV/units.energy*min(fac,1);
    }
  }

  return(U);
}



//*************standard ring/single dot/double dot****************//

//inside any of the leads - returns its number, otherwise zero (does not mean it is inside the conductor)
int inside_lead(prec x, prec y) 
{
  int leads=(int) parameters.values[parameters_class::leads];
  for (int l=1;l<=abs(leads);l++) {
    switch (lead[l].heading) {
      case 1 : {if (x<lead[l].x && abs(y-lead[l].y)<lead[l].width/2) return(l);break;}
      case 2 : {if (x>lead[l].x && abs(y-lead[l].y)<lead[l].width/2) return(l);break;}
      case 3 : {if (y<lead[l].y && abs(x-lead[l].x)<lead[l].width/2) return(l);break;}
      case 4 : {if (y>lead[l].y && abs(x-lead[l].x)<lead[l].width/2) return(l);break;}
    }
  }
  return(0);
}

//on the axis of one of the leads
bool on_a_lead_axis(prec x, prec y)
{
  int leads=(int) parameters.values[parameters_class::leads];
  for (int l=1;l<=leads;l++) {//make sure the lead corridor is inside
    switch (lead[l].heading) {
      case 1 : ;
      case 2 : {if ((abs(y-lead[l].y)<lead[l].width/2) && (x*lead[l].x>0)) return(true);break;}
      case 3 : ;
      case 4 : {if ((abs(x-lead[l].x)<lead[l].width/2) && (y*lead[l].y>0)) return(true);break;}
    }
  }
  return(false);
}

//beyond the edge (longitudinaly) of one of the leads
bool beyond_leads(prec x, prec y) 
{
  int leads=(int) parameters.values[parameters_class::leads];
  bool ans=true;
  for (int l=1;l<=abs(leads);l++) {
    switch (lead[l].heading) {
      case 1 : {if (x>lead[l].x) ans=false;break;}
      case 2 : {if (x<lead[l].x) ans=false;break;}
      case 3 : {if (y>lead[l].y) ans=false;break;}
      case 4 : {if (y<lead[l].y) ans=false;break;}
    }
  }
  return(ans);
}

//the sector this point belongs to (0:conductor, 1,...:leads, -1:outside
int inside_sector(prec x, prec y) 
{
  //fprintf(logg,"inside_sector: (x=%.3e, y=%.3e)\n",x,y);
  int l=inside_lead(x,y);
  //fprintf(logg,"inside_lead=%i, on_a_lead_axis=%i, beyond_leads=%i, inside=%i\n", inside_lead(x,y), on_a_lead_axis(x,y),  beyond_leads(x,y), inside(x,y));
  if (l>0) return(l);
  if (inside(x,y)) return(0);
  return(-1);
}

//graphene define by Fermi-Dirac like mass profile
//#define FDlike

//all active points including the leads
bool inside(prec x, prec y)
{
  bool ans=false;
  //ans=true;

  prec d_x=parameters.values[parameters_class::d_x];
  prec d_y=parameters.values[parameters_class::d_y];  
  prec dLx=parameters.values[parameters_class::dLx];
  prec dLy=parameters.values[parameters_class::dLy];
  prec lB=parameters.values[parameters_class::lB];
  prec extent=parameters.values[parameters_class::extent];
  prec ring=parameters.values[parameters_class::ring];
  int geometry=(int) parameters.values[parameters_class::geometry];
  int leads=(int) parameters.values[parameters_class::leads];
  bool ALHO=(bool) parameters.values[parameters_class::ALHO];
  
  prec m_eff=parameters.values[parameters_class::m_eff];
  prec lv=parameters.values[parameters_class::lv];
  prec prefac=pow(units.length,2)*pow(units.hbar,2)/(2*units.m_e*m_eff*pow(lv*units.length,4));
  prefac/=units.energy;
  
  int system=(int) parameters.values[parameters_class::system];
  if (system=='g') {
#ifdef FDlike
    if (sqrt(x*x+y*y)<parameters.values[parameters_class::gr_w]) ans=true;
#else
    if (mass_gr(x,y)/parameters.values[parameters_class::gr_V0]<exp(2*extent)) ans=true;
#endif
  }
  else {
    switch (geometry) {
      case (3) :
      case (0) : {  //single dot/double dot/ring
	x-=dLx;
	y-=dLy;
	prec r2right=(x-d_x)*(x-d_x)+(y-d_y)*(y-d_y);
	prec r2left=(x+d_x)*(x+d_x)+(y+d_y)*(y+d_y);
	prec dist=lB*extent;
	if (fabs(sqrt(r2right)-ring)<dist || fabs(sqrt(r2left)-ring)<dist) ans=true;
	if (vulcano_potential_ring && ring>0) { 
	  if (prefac*(pow(ring,4)/min(r2right,r2left) + min(r2right,r2left) - 2*pow(ring,2))/4.>5E3) ans=false; //cutting the grid outside the ring using extent, inside using a potential threshold
	  }
	if (ALHO) {//anisotropic LHO lx*=(1+d/lv), ly/=(1+d/lv)
	  prec d=parameters.values[parameters_class::d];
	  prec phi_d=parameters.values[parameters_class::phi_d];
	  prec d_x=cos(phi_d), d_y=sin(phi_d);
	  prec rd2=(x*d_x+y*d_y)*(x*d_x+y*d_y);
	  prec rcd2=(x*d_y-y*d_x)*(x*d_y-y*d_x);
	  ans=false;
	  prec fac=1+d/lB;
	  //prec dist2=(lB+d)*extent;
	  if (rcd2/(dist/fac*dist/fac) + rd2/(dist*fac*dist*fac)<1) ans=true;  
	}
	if (geometry==3) {
	  prec xw, yw;
	  if (dot2waveguide(x,y,xw,yw)) ans=true;
	}
	break;
      }
      case (1) : {return(true);break;}
      case (2) : {
	prec flat_x=parameters.values[parameters_class::flat_x];
	prec flat_y=parameters.values[parameters_class::flat_y];
	if (abs(x)<flat_x && abs(y)<flat_y) ans=true;
	if (on_a_lead_axis(x,y)) ans=true;
	//if (x*x+y*y<=pow(12.1*units.nm/units.length,2)) ans=false;
	break;
      }
    }
  }
  if (leads>0) {
    if (beyond_leads(x,y)) ans=false;
    if (inside_lead(x,y)) ans=true;
  }
  return(ans);
}


prec potential(prec x, prec y) {

  prec m_eff=parameters.values[parameters_class::m_eff];
  prec d_x=parameters.values[parameters_class::d_x];
  prec d_y=parameters.values[parameters_class::d_y];
  prec phi_d=parameters.values[parameters_class::phi_d];
  prec d=parameters.values[parameters_class::d];
  prec lv=parameters.values[parameters_class::lv];
  prec Ebias=parameters.values[parameters_class::Ebias];
  prec ring=parameters.values[parameters_class::ring];
  int geometry=(int) parameters.values[parameters_class::geometry];
  int leads=(int) parameters.values[parameters_class::leads];
  bool ALHO=(bool) parameters.values[parameters_class::ALHO];

  prec r2right=(x-d_x)*(x-d_x)+(y-d_y)*(y-d_y);
  prec r2left=(x+d_x)*(x+d_x)+(y+d_y)*(y+d_y);
  //prec r2=x*x+y*y;
  prec res=0;

  prec prefac=pow(units.length,2)*pow(units.hbar,2)/(2*units.m_e*m_eff*pow(lv*units.length,4));
  prefac/=units.energy;
  
  //!anisotropic LHO!!! d is the enlargement of the confinement lengths, added to lv
  if (ALHO) {
    prec d_x=cos(phi_d), d_y=sin(phi_d);
    prec rd2=(x*d_x+y*d_y)*(x*d_x+y*d_y);
    prec rcd2=(x*d_y-y*d_x)*(x*d_y-y*d_x);
    r2left=r2right=rd2/pow(1+d/lv,4)+rcd2*pow(1+d/lv,4);
  }
  
  //!exponential potential
  //p1=p2=1-exp(-pow(abs(ydot)/4/lv,1)-pow(abs(xdot)/4/(lv+d),1));
  //prefac=10*units.meV/units.energy;
  
  //printf("prefac in potential: %e\n",prefac);
  //exit(1);
  
  switch (geometry) {
    case (3) :
    case (0) : {
      if (ring==0) res=min(r2left,r2right)*prefac;                 //bi-quadratic
      //if (ring==0) res=((rd2-d2)*(rd2-d2)/(4*d2)+rcd2)*prefac;   //quartic
      else {
        prec vleft=r2left-2*sqrt(r2left)*ring+ring*ring;
        prec vright=r2right-2*sqrt(r2right)*ring+ring*ring;
        res=min(vleft,vright)*prefac;
        if(vulcano_potential_ring && ring>0) { //vulcano 
          if(min(r2left,r2right)==0) res=prefac*1E10;//removes the divergence at the center of the ring 
          else { 
            res=prefac*min((pow(ring,4)/r2left + r2left - 2*pow(ring,2))/4 , (pow(ring,4)/r2right + r2right - 2*pow(ring,2))/4); 
          }
        } 
      }
      if (Ebias!=0) res+=-units.e*Ebias*(x*cos(phi_d)+y*sin(phi_d))*units.length*units.el_intensity/units.energy;
      if (geometry==3) {
        prec xw, yw;
        if (dot2waveguide(x,y,xw,yw)) {
          double barrier=parameters.values[parameters_class::wg_V];
          if (res>barrier) res=barrier;
        }
      }
      break;
    }
    case (1) : {
      prec phi=x/ring;
      break;
    }
    case (2) : {
      if (GF_check_flag) res=(potential_flat_GF_check(x,y)).real();
      else res=potential_flat(x,y);
      break;
    }
  }
  
  if (leads>0) switch (leads) {
    for (int l=1;l<=abs(leads);l++) {
      switch (lead[l].heading) {
        case 1 : {if (x<lead[l].x && abs(y-lead[l].y)<lead[l].width/2) res=lead[l].potential;break;}
        case 2 : {if (x>lead[l].x && abs(y-lead[l].y)<lead[l].width/2) res=lead[l].potential;break;}
        case 3 : {if (y<lead[l].y && abs(x-lead[l].x)<lead[l].width/2) res=lead[l].potential;break;}
        case 4 : {if (y>lead[l].y && abs(x-lead[l].x)<lead[l].width/2) res=lead[l].potential;break;}
      }
    }
  }
  res+=impurities.potential(x,y);
  return(res);
} 



prec potential_waveguide(prec x)
{
  double length=parameters.values[parameters_class::wg_length];
  double V=parameters.values[parameters_class::wg_V]+parameters.values[parameters_class::dot_offset];
  if (x<=0) return(V);
  if (x>length) return(0); 
  return(V*(1-x/length));
}

prec tunnel_barrier(prec x)
{
  double length=parameters.values[parameters_class::wg_length];
  if (x<=0 || x>=length) return(0);
  double l1=30*units.nm/units.length, l2=8*units.nm/units.length;
  //double ext1=l1*5, ext2=l2*5;
  double V=parameters.values[parameters_class::wg_V];
  //prec x1=(x-ext1);
  //prec x2=(x-(length-ext2));  
  //prec res=V*(exp(-x1*x1/(l1*l1))+exp(-x2*x2/(l2*l2)));
  //if (x>l1 && x<length-l2) res=max(res,7*units.meV/units.energy);
  if (x<l1 || x>length-l2) return(V); else return(7*units.meV/units.energy);
}


//*****************flat region************************
/*bool inside(prec x, prec y) 
{
  //prec lw=parameters.values[parameters_class::L_lead_width];
  //prec rw=parameters.values[parameters_class::R_lead_width];
  //prec lp=parameters.values[parameters_class::L_lead_pos];
  //prec rp=parameters.values[parameters_class::R_lead_pos];
  
  if (abs(y)<max(lw,rw)/2) return(true);
  return(false); 
}


prec potential(prec x, prec y) {
  //prec a=(aux_Z(x,y)).real();
  //return(a*0.001*E0*1e-3*units.e/units.energy);
  //prec E=parameters.values[parameters_class::def];
  //if (abs(x)<25/units.length*1e-9) return(bar);
  return(0);
} 
*/

content u_const;
void hamiltonian(int N, void *A, content *x, content *y)
{
  reg1* pregion=(reg1*) A;
  pregion->operate_tot(x,y);
}

void hamiltonian_graphene(int N, void *A, content *x, content *y)
{
  region* pregion=(region*) A;
  pregion->inverseLUtimesM(x,y);
}

void hamiltonian_graphene_with_checking(int N, void *A, content *x, content *y)
{
  region* pregion=(region*) A;
  pregion->inverseLUtimesM_with_checking(x,y);
}


//dof_prototype dof;

//quadratic 'potential' which reaches 1 meV at 20 nm from the center
/*prec potentialg(prec x, prec y) 
{
  prec r=sqrt(x*x+y*y)*units.length/units.nm;
  prec res=(r/20)*(r/20)+pow(r/40,4);
  res*=1*units.meV/units.energy;
  return(res);
} */

//tunable double dot 'mass'
prec mass_gr(prec x, prec y) 
{
  prec dd=parameters.values[parameters_class::gr_dd];
  prec w=parameters.values[parameters_class::gr_w];
  prec r0=parameters.values[parameters_class::gr_r0];
  
  prec phi_d=parameters.values[parameters_class::phi_d];
  prec d=parameters.values[parameters_class::d];

#ifdef FDlike
  prec r=sqrt(x*x+y*y);
  prec res=1.0/(1.0+exp(-(r-r0)/dd));
#else
  prec xd=x*cos(phi_d)+y*sin(phi_d);
  prec yd=y*cos(phi_d)-x*sin(phi_d);
  prec d1=sqrt((xd+d)*(xd+d)+yd*yd);  
  prec d2=sqrt((xd-d)*(xd-d)+yd*yd);
  
  prec v1=exp((d1-r0)/dd);
  prec v2=exp((d2-r0)/dd);
  prec vbr=exp((abs(xd)-d)/dd)+exp((abs(yd)-w)/dd);
  
  prec res=min(vbr,min(v1,v2));
#endif
  
  res*=parameters.values[parameters_class::gr_V0];
  return(res);
} 

//tunable double dot 'potential'
prec potential_gr(prec x, prec y) 
{
  prec r0=parameters.values[parameters_class::gr_r0];
  prec phi_d=parameters.values[parameters_class::phi_d];
  prec d=parameters.values[parameters_class::d];

  prec xd=x*cos(phi_d)+y*sin(phi_d);
  prec yd=y*cos(phi_d)-x*sin(phi_d);

  prec res=0;
  if (xd>-d+r0 && xd<d-r0) res=1;

  res*=parameters.values[parameters_class::gr_U];
  return(res);
} 

void hamiltonian_init_graphene(reg1& region)
{
  int N=dof.length;	 
  if (N!=(int) parameters.read("sets","-")) {
    printf("length(%i) vs sets(%i) in hamiltonian_init_graphene\n", N, (int) parameters.read("sets","-"));
  }
  //number of (sigma + orbital) sets is N/2
  //sigma is supposed to be the first dof, then:
  //sigma = (I mod 2) while (I/2) labels the "subblocks"

  //kinetic energy
  content u=parameters.read("vf","m/s")*units.hbar/units.length/units.energy; //units
  dof.reset("all");
  dof.operatorC[dof_prototype::sigma]=dof_prototype::splus;
  dof.convert("OP2M",u);
  for (int i=0;i<N;i++) for (int j=0;j<N;j++) {
    //fprintf(logg,"M(%i,%i)=%.3e%+.3ei\n",i,j,dof.M_sum[i*N+j].real(),dof.M_sum[i*N+j].imag());
    if (abs(dof.M_sum[i*N+j])!=0) {
      region.update_g_matrix(i%2,j%2,"k-",i/2,j/2,dof.M_sum[i*N+j],0,exp_integral, region.HG, region.HG_sparse);
    }
  }
  dof.reset("all");
  dof.operatorC[dof_prototype::sigma]=dof_prototype::sminus;
  dof.convert("OP2M",u);
  for (int i=0;i<N;i++) for (int j=0;j<N;j++) {
    //fprintf(logg,"M(%i,%i)=%.3e%+.3ei\n",i,j,dof.M_sum[i*N+j].real(),dof.M_sum[i*N+j].imag());
    if (abs(dof.M_sum[i*N+j])!=0) {
      region.update_g_matrix(i%2,j%2,"k+",i/2,j/2,dof.M_sum[i*N+j],0,exp_integral,region.HG, region.HG_sparse);
    }
  }  
  dof.reset("all");
  dof.operatorC[dof_prototype::sigma]=dof_prototype::unity;
  dof.convert("OP2M",1.0);
  for (int i=0;i<N;i++) for (int j=0;j<N;j++) {
    //fprintf(logg,"M(%i,%i)=%.3e%+.3ei\n",i,j,dof.M_sum[i*N+j].real(),dof.M_sum[i*N+j].imag());
    if (abs(dof.M_sum[i*N+j])!=0) {
      region.update_g_matrix(i%2,j%2,"unity",i/2,j/2,dof.M_sum[i*N+j],0,exp_integral, region.M, region.M_sparse);
    }
  }
  //fprintf(logg,"\nend of kinetic energy (g_length=%i)\n",region.g_length);
  
  //"extrinsic" gap (=mass) to confine it
  dof.reset("all");
  dof.operatorC[dof_prototype::sigma]=dof_prototype::pauliz;
  dof.convert("OP2M",1.0);
  for (int i=0;i<N;i++) for (int j=0;j<N;j++) {
    //fprintf(logg,"M(%i,%i)=%.3e%+.3ei\n",i,j,dof.M_sum[i*N+j].real(),dof.M_sum[i*N+j].imag());
    if (abs(dof.M_sum[i*N+j])!=0) {
      region.update_g_matrix(i%2,j%2,"unity",i/2,j/2,dof.M_sum[i*N+j],mass_grc,exp_integral, region.HG, region.HG_sparse);
    }
  }
  
  //electrostatic potential
  dof.reset("all");
  dof.operatorC[dof_prototype::sigma]=dof_prototype::unity;
  dof.convert("OP2M",1.0);
  for (int i=0;i<N;i++) for (int j=0;j<N;j++) {
    //fprintf(logg,"M(%i,%i)=%.3e%+.3ei\n",i,j,dof.M_sum[i*N+j].real(),dof.M_sum[i*N+j].imag());
    if (abs(dof.M_sum[i*N+j])!=0) {
      region.update_g_matrix(i%2,j%2,"unity",i/2,j/2,dof.M_sum[i*N+j],potential_grc,exp_integral, region.HG, region.HG_sparse);
    }
  }
  
  //the rest - no orbital operators
  dof.reset("all");
  int tau=1;//only "up" of the kappa subspace
  
  //gap
  dof.reset("operators");
  dof.operatorC[dof_prototype::sigma]=dof_prototype::pauliz;
  dof.convert("OP2M",tau*parameters.read("Dg","meV")*units.meV/units.energy);
  
  //the rest only if spin dof is considered
  if (dof_prototype::spin<dof.Ndof) {

    //Zeeman
    prec B=parameters.read("B","Tesla");
    prec theta_B=parameters.read("theta_B","pi")*M_PI;
    prec phi_B=parameters.read("phi_B","pi")*M_PI;
    dof.reset("operators");
    dof.operatorC[dof_prototype::spin]=dof_prototype::paulix;
    dof.convert("OP2M",units.muB*parameters.read("gx","-")/2*B*sin(theta_B)*cos(phi_B)/units.energy);
    dof.operatorC[dof_prototype::spin]=dof_prototype::pauliy;
    dof.convert("OP2M",units.muB*parameters.read("gy","-")/2*B*sin(theta_B)*sin(phi_B)/units.energy);
    dof.operatorC[dof_prototype::spin]=dof_prototype::pauliz;
    dof.convert("OP2M",units.muB*parameters.read("gz","-")/2*B*cos(theta_B)/units.energy);


    //intrinsic spin-orbit
    dof.reset("operators");
    dof.operatorC[dof_prototype::sigma]=dof_prototype::pauliz;
    dof.operatorC[dof_prototype::spin]=dof_prototype::pauliz;
    dof.convert("OP2M",tau*parameters.read("Di","meV")*units.meV/units.energy);

    //Rashba spin-orbit
    dof.reset("operators");
    dof.operatorC[dof_prototype::sigma]=dof_prototype::paulix;
    dof.operatorC[dof_prototype::spin]=dof_prototype::pauliy;
    dof.convert("OP2M",tau*parameters.read("Dr","meV")*units.meV/units.energy);
    dof.operatorC[dof_prototype::sigma]=dof_prototype::pauliy;
    dof.operatorC[dof_prototype::spin]=dof_prototype::paulix;
    dof.convert("OP2M",-parameters.read("Dr","meV")*units.meV/units.energy);

    for (int i=0;i<N;i++) for (int j=0;j<N;j++) {
      //fprintf(logg,"M(%i,%i)=%.3e%+.3ei\n",i,j,dof.M_sum[i*N+j].real(),dof.M_sum[i*N+j].imag());
      if (abs(dof.M_sum[i*N+j])!=0) region.update_g_matrix(i%2,j%2,"unity",i/2,j/2,dof.M_sum[i*N+j],0,exp_integral, region.HG, region.HG_sparse);
    }
  }
}


/*void hamiltonian_init_graphene2(reg1& region)
{
  region.act_prec_set((reg1::operators_precision) (int) parameters.read("precision","-"));
  region.symmetrize_ops((bool) parameters.read("peierls","-"), true);
  bool op_rmmbr[reg1::all_op]={false};
  for (int i=0;i<reg1::all_op;i++) region.remember_set((reg1::operators) i,op_rmmbr[i]);
  region.operate_tot_init();

  int N=dof.length;	 //number of sets (length of the consecutive index)
  if (N!=(int) parameters.read("sets","-")) 
    {printf("length(%i) vs sets(%i) in hamiltonian_init_graphene\n", N, (int) parameters.read("sets","-"));}
  //content u=content(0,-1)*parameters.read("vf","m/s")*units.hbar/units.length/units.energy; //units
  //region.operate_tot_init(0,1,reg1::dx,u,0,exp_integral);
  //region.operate_tot_init(1,0,reg1::dx,u,0,exp_integral);
  //region.operate_tot_init(0,1,reg1::dy,content(0,1)*u,0,exp_integral);
  //region.operate_tot_init(1,0,reg1::dy,content(0,-1)*u,0,exp_integral);
  //return;

  //no potential - confinement is fixed by the boundary?
  //for (int i=0;i<N;i++) region.operate_tot_init(i,i,reg1::lin_com, content(1,0),potentialgc,exp_integral);
  dof.reset("all");
  dof.operatorC[dof_prototype::sigma]=dof_prototype::pauliz;
  //dof.convert("OP2M",100.0*units.meV/units.energy);   
  dof.convert("OP2M",1.0);  
  for (int i=0;i<N;i++) for (int j=0;j<N;j++) {
    fprintf(logg,"M(%i,%i)=%.3e%+.3ei\n",i,j,dof.M_sum[i*N+j].real(),dof.M_sum[i*N+j].imag());
    if (abs(dof.M_sum[i*N+j])!=0) region.operate_tot_init(j,i,reg1::lin_com,dof.M_sum[i*N+j],potentialc,exp_integral);
  }
  prec hx, hy;
  region.give_par(hx,hy);
  prec coors[2*1+1], coefs[2*1+1];
  for (int i=-1;i<=+1;i++) coors[i+1]=(i-0.5)*hx;
  construct_differential_operator(1,2*1+1,coors,coefs);

  //kinetic energy is linear in k
  content u=content(0,-1)*parameters.read("vf","m/s")*units.hbar/units.length/units.energy; //units
  u_const=u*coefs[1]*coefs[2];
  //printf("correction factor is %e%+ei\n",u_const.real(), u_const.imag());
  dof.reset("all");
  dof.operatorC[dof_prototype::sigma]=dof_prototype::paulix;
  dof.convert("OP2M",u);
  for (int i=0;i<N;i++) for (int j=0;j<N;j++) {
    fprintf(logg,"M(%i,%i)=%.3e%+.3ei\n",i,j,dof.M_sum[i*N+j].real(),dof.M_sum[i*N+j].imag());
    if (abs(dof.M_sum[i*N+j])!=0) {
      if (i==j) region.operate_tot_init(j,i,reg1::dx,dof.M_sum[i*N+j],0,exp_integral);
      if (i>j) region.operate_tot_init(j,i,reg1::dxp,dof.M_sum[i*N+j],0,exp_integral);
      if (i<j) region.operate_tot_init(j,i,reg1::dxp,dof.M_sum[i*N+j],0,exp_integral);
    }
  }
  dof.reset("all");
  dof.operatorC[dof_prototype::sigma]=dof_prototype::pauliy;
  dof.convert("OP2M",u);
  for (int i=0;i<N;i++) for (int j=0;j<N;j++) {
    fprintf(logg,"M(%i,%i)=%.3e%+.3ei\n",i,j,dof.M_sum[i*N+j].real(),dof.M_sum[i*N+j].imag());
    if (abs(dof.M_sum[i*N+j])!=0) {
      if (i==j) region.operate_tot_init(j,i,reg1::dy,dof.M_sum[i*N+j],0,exp_integral);
      if (i>j) region.operate_tot_init(j,i,reg1::dyp,dof.M_sum[i*N+j],0,exp_integral);
      if (i<j) region.operate_tot_init(j,i,reg1::dyp,dof.M_sum[i*N+j],0,exp_integral);
    }
  }
  fprintf(logg,"\nend of orbital energy\n");

  //the rest - no orbital operators
  dof.reset("all");
  int tau=1;//only "up" of the kappa subspace
  
  //Zeeman (if spin dof is considered)
  if (dof_prototype::spin<dof.Ndof) {
    prec B=parameters.read("B","Tesla");
    prec theta_B=parameters.read("theta_B","pi")*M_PI;
    prec phi_B=parameters.read("phi_B","pi")*M_PI;
    dof.reset("operators");
    dof.operatorC[dof_prototype::spin]=dof_prototype::paulix;
    dof.convert("OP2M",units.muB*parameters.read("gx","-")/2*B*sin(theta_B)*cos(phi_B)/units.energy);
    dof.operatorC[dof_prototype::spin]=dof_prototype::pauliy;
    dof.convert("OP2M",units.muB*parameters.read("gy","-")/2*B*sin(theta_B)*sin(phi_B)/units.energy);
    dof.operatorC[dof_prototype::spin]=dof_prototype::pauliz;
    dof.convert("OP2M",units.muB*parameters.read("gz","-")/2*B*cos(theta_B)/units.energy);
  }
  
  //gap
  dof.reset("operators");
  dof.operatorC[dof_prototype::sigma]=dof_prototype::pauliz;
  dof.convert("OP2M",tau*parameters.read("Dg","meV")*units.meV/units.energy);
  
  //intrinsic spin-orbit
  if (parameters.read("Di","meV")!=0) {
    dof.reset("operators");
    dof.operatorC[dof_prototype::sigma]=dof_prototype::pauliz;
    dof.operatorC[dof_prototype::spin]=dof_prototype::pauliz;
    dof.convert("OP2M",tau*parameters.read("Di","meV")*units.meV/units.energy);
  }

  //Rashba spin-orbit
  if (parameters.read("Dr","meV")!=0) {
    dof.reset("operators");
    dof.operatorC[dof_prototype::sigma]=dof_prototype::paulix;
    dof.operatorC[dof_prototype::spin]=dof_prototype::pauliy;
    dof.convert("OP2M",tau*parameters.read("Dr","meV")*units.meV/units.energy);
    dof.operatorC[dof_prototype::sigma]=dof_prototype::pauliy;
    dof.operatorC[dof_prototype::spin]=dof_prototype::paulix;
    dof.convert("OP2M",-parameters.read("Dr","meV")*units.meV/units.energy);
  }
  
  for (int i=0;i<N;i++) for (int j=0;j<N;j++) {
    fprintf(logg,"M(%i,%i)=%.3e%+.3ei\n",i,j,dof.M_sum[i*N+j].real(),dof.M_sum[i*N+j].imag());
    if (abs(dof.M_sum[i*N+j])!=0) region.operate_tot_init(j,i,reg1::lin_com,dof.M_sum[i*N+j],0,exp_integral);
  }
  fflush(logg);
}*/

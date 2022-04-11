#include "main.h" 

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!two component spinor
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //Pauli matrices 1, x, y, z
  content pauli[]=   {1.0, 0, 0, 1.0,\
		      0, 1.0, 1.0, 0,\
		      0, -content(0,1), content(0,1), 0,\
                      1.0, 0, 0, -1.0};

void spinor::construct(double theta, double phi)
{
  up=cos(theta/2);
  down=sin(theta/2)*exp(content(0,phi));
}

void spinor::angles(double& theta, double& phi)
{
  if (down==content(0,0)) {theta=0;phi=0;return;}
  content ratio=up/down;
  double r=abs(ratio);
  phi=-arg(ratio);
  if (phi<0) phi+=2*M_PI;
  theta=2*asin(sqrt(1/(1+r*r)));
  double res=abs(ratio-cos(theta/2)/sin(theta/2)*exp(content(0,-phi)));
  if (res>1e-8) {printf("wrong conversion in spinor::angles\n");exit(1);}
  //fprintf(logg,"input spinor: (%.2e%+.2ei,%.2e%+.2ei) results in angles (%.3e,%.3e) [pi]\n",xi[0].real(),xi[0].imag(),xi[1].real(),xi[1].imag(),theta,phi);
  return;
}

spinor spinor::SigmaTimesMe(content* Sigma)
{
  spinor res;
  res.up=Sigma[0*2+0]*up+Sigma[0*2+1]*down;
  res.down=Sigma[1*2+0]*up+Sigma[1*2+1]*down;
  return(res);
}

spinor spinor::add(double coef, spinor s)
{
  spinor res;
  res.up=up+coef*s.up;
  res.down=down+coef*s.down;
  return(res);
}

spinor spinor::multiply(double coef)
{
  spinor res;
  res.up=coef*up;
  res.down=coef*down;
  return(res);
}

spinor spinor::conjugate()
{
  spinor res;
  res.up=conj(up);
  res.down=conj(down);
  return(res);
}

content spinor::braket(spinor ket)
{
  content res=conj(up)*ket.up+conj(down)*ket.down;
  return(res);
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!***********************************acr_stack
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//used for adjusting the step
//calling another_in with consecutive number creates (if not existing) stacks keeping track of the values at previous steps
//length: 2 - linear extrapolation, two values enough
//        3 - quadratic extrapolation
//kind:   0 - normal data
//        1 - need derivative before exptrapolating - length needs to be the method + 1
//adjust step goes through all stacks, those filled give step
//the smallest in the range (min, max) is taken
//a step is set shorter_step times shortest step to get to the predicted zero

#define ACR_STACK0
//#define ACR_STACK1

acr_stack::acr_stack() {hmn_exist=0;}

acr_stack::~acr_stack()
{
  for (int i=0;i<hmn_exist;i++) {
    stack* second=first->next;
    delete first;
    first=second;
  }
}

stack* acr_stack::give_stack(int which)
{
  stack *wished=first;
  for (int i=0;i<which;i++) wished=wished->next;
  return(wished);
}


void acr_stack::another_in(int which, prec coor, prec value,int length, int kind,char* tag_in)
{

  if (which>hmn_exist+1) {
    message("ACR_STACK: stacks can be created only consecutively\n",1+4);
    exit(1);
  }
  if (which==hmn_exist) {//create a new stack
#ifdef ACR_STACK1
    sprintf(buffer,"ACR_STACK: stack #%i with tag '%s' of length %i and kind %i will be created\n", which, tag_in, length, kind);
    splachni(logg,buffer,4);
#endif
    stack* neu=new stack(length, kind, tag_in);
    if (hmn_exist==0) {//the first one

    first=last=neu;
    }
    //next one
    last->next=neu;
    last=neu;
    hmn_exist++;
  }

  give_stack(which)->another_in(coor, value, tag_in);
}

prec acr_stack::adjust_step(prec position, prec min, prec max, prec shorter_step)
{
  int pntr=0;
  prec step=max;
  int w=-1;
  for (int i=0;i<hmn_exist;i++) {
    stack* x=give_stack(i);

#ifdef ACR_STACK1
  sprintf(buffer,"ACR_STACK: for adjustind step stack#%i '%s' is tried: filled=%i\n",i,x->tag,x->filled);
  splachni(logg,buffer,4);
#endif
    
    if (x->filled<x->length) continue;
    pntr++;
    prec zero_at=give_stack(i)->coordinate(0);
    prec distance=zero_at-position;
    prec step_i;
    if (distance>0) step_i=distance/shorter_step;
    else step_i=-distance/shorter_step*2;
    if (step_i<step) {
      step=step_i;
      w=i;
    }
#ifdef ACR_STACK1
    sprintf(buffer,"ACR_STACK: zero prediction lead to distance %.4e giving step %.4e\n\n",distance, step_i);
    splachni(logg,buffer,4);
#endif
  }
  if (step<min || pntr==0) step=min;
#ifdef ACR_STACK0
  if (w==-1) {
  sprintf(buffer,"ACR_STACK: no stack predicted lower step as maximum\n");
  splachni(logg,buffer,4);
  }
  else {
    stack* lowest=give_stack(w);
    sprintf(buffer,"ACR_STACK: stack#%i '%s' predicted the lowest step from:\ncoordinates:\t",w,give_stack(w)->tag);
    splachni(logg,buffer,4);
    for (int i=0;i<lowest->length;i++) {
      sprintf(buffer,"%.4e\t",lowest->coors[i]);
      splachni(logg,buffer,4);
    }
    message(logg,"\nvalues:\t",4);
    for (int i=0;i<lowest->length;i++) {
      sprintf(buffer,"%.4e\t",lowest->vals[i]);
      splachni(logg,buffer,4);
    }
    message(logg,"\n\n",4);
  }
  sprintf(buffer,"ACR_STACK: step of %.4e was set\n\n",step);
  splachni(logg,buffer,4);
#endif
  return(step);
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!***********************************stack
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//used to extrapolate
//length of values with coordinates are stored
//and used to extrapolate a value at given coordinate
//or approximate a coordinate for a given value

#define STACK0
#define STACK1

stack::stack(int l, int k, char* tag_in)
{
  filled=0;
  vals=new prec[l];
  coors=new prec[l];
  length=l;
  kind=k;
  pointer_to_empty=0;
  strcpy(tag, tag_in);
  
#ifdef STACK0
  sprintf(buffer,"stack with tag '%s' of kind %i and length %i created\n\n",tag,k,l);
  splachni(logg,buffer,4);
#endif
}

stack::~stack()
{
  delete vals;
  delete coors;
#ifdef STACK0
  sprintf(buffer,"stack (with tag '%s') deleted\n",tag);
  splachni(logg,buffer,4);
#endif
}

void stack::another_in(prec coor, prec value, char* tag_in)
{
  for (int i=length-1;i>0;i--) {
    vals[i]=vals[i-1];
    coors[i]=coors[i-1];
  }
  vals[0]=value;
  coors[0]=coor;
  if (pointer_to_empty<length) pointer_to_empty++;
  filled=pointer_to_empty;

  strcpy(tag,tag_in);
  
#ifdef STACK1
  sprintf(buffer,"a couple (%.2e,%.2e) added to stack with message '%s'\nresulting status: filled=%i,\ncoordinates:\t",coor,value,tag,filled);
  splachni(logg,buffer,4);
  for (int i=0;i<length;i++) {
    sprintf(buffer,"%.4e\t",coors[i]);
    splachni(logg,buffer,4);
  }
  message(logg,"\nvalues:\t",4);
  for (int i=0;i<length;i++) {
    sprintf(buffer,"%.4e\t",vals[i]);
    splachni(logg,buffer,4);
  }
  message(logg,"\n\n",4);

#endif

}

void stack::derivative(prec* coors_out, prec* vals_out)
{
  for (int i=0;i<length-1;i++) {
    coors_out[i]=(coors[i+1]+coors[i])/2;
    vals_out[i]=(vals[i+1]-vals[i])/(coors[i+1]-coors[i]);
  }
}

prec stack::coordinate(prec value)
{
  prec coordinate;
  
  switch (kind) {
    case 1 : {//normal
      coordinate=extrapolate_to_value(value,coors,vals,length-1);
      break;
    }
    case 2 : {//do derivative
      prec *coors_d=new prec[length-1];
      prec *vals_d=new prec[length-1];
      derivative(coors_d, vals_d);
      coordinate=extrapolate_to_value(value,coors_d,vals_d,length-2);
      delete coors_d;
      delete vals_d;
      break;
    }
    case 3: {//fit to anticrossing
      //only if value is zero
      if (value!=0) { message("anticrossing-like fitting only to zero energy distance in stack::coordinate\n",1+4);exit(1);}
      if (length<3) { message("anticorssing stack must be at least of length 3\n",1+4); exit(1);}
      coordinate=extrapolate_to_acr(coors, vals);
      break;
    }
    default: {
      message("wrong kind of stack in coordinate\n",1+4);
      exit(1);
    }
  }
#ifdef STACK1
  sprintf(buffer,"stack '%s'\n predicted value %e at %e\nfrom:\ncoordinates:\t",tag,value, coordinate);
  splachni(logg,buffer,4);
  for (int i=0;i<length;i++) {
    sprintf(buffer,"%.4e\t",coors[i]);
    splachni(logg,buffer,4);
  }
  message(logg,"\nvalues:\t",4);
  for (int i=0;i<length;i++) {
    sprintf(buffer,"%.4e\t",vals[i]);
    splachni(logg,buffer,4);
  }
  message(logg,"\n\n",4);
#endif

  return(coordinate);
}
prec stack::value(prec coordinate)
{
  
  if (kind==2) {
    message("prediction not implemented for anticrossing stack in stack::value\n",1+4); exit(1);}
  
  prec coefs[3];
  extrapolate_coefs(coors,vals,length-1,coefs);
  //prec a=coefs[2],b=coefs[1],c=coefs[0];
  
  prec res=0;
  prec cn=1;
  for (int i=0;i<length;i++) {
    res+=coefs[i]*cn;
    cn*=coordinate;
  }
#ifdef STACK1
  sprintf(buffer,"stack '%s'\n predicted at coordinate %e a value %e\nfrom:\ncoordinates:\t",tag,coordinate, res);
  splachni(logg,buffer,4);
  for (int i=0;i<length;i++) {
    sprintf(buffer,"%.4e\t",coors[i]);
    splachni(logg,buffer,4);
  }
  message(logg,"\nvalues:\t",4);
  for (int i=0;i<length;i++) {
    sprintf(buffer,"%.4e\t",vals[i]);
    splachni(logg,buffer,4);
  }
  message(logg,"\n\n",4);
#endif

  return(res);
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!***********************************units
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
class units_class units;

units_class::units_class()
{
  m_e=9.10938188e-31;
  m_p=1.672621777e-27;
  hbar=1.05457148e-34;
  e=1.60217646e-19;
  eps=8.854187817e-12;
  mi=4*M_PI*1e-7;
  e2prime=e*e/(4*M_PI*eps);
  kB=1.3806503e-23;
  meV=e*1e-3;
  nm=1e-9;
  Na=6.022141e+23;
  c=299792458.0;
  Debye=1/c*1e-21;
  aGaAs=5.65e-10;
  muB=e*hbar/m_e/2;
  muN=e*hbar/m_p/2;
  aSi=5.43e-10;

  double scale=1000;

  mass=m_e;
  length=scale*hbar*hbar/(mass*e2prime);
  energy=hbar*hbar/(mass*length*length);
  omega=energy/hbar;
  B=hbar/(e*length*length);
  A=length*B;
  lC=hbar*hbar/(e2prime*mass);

  wavevector=1/length;
  momentum=hbar/length;
  el_intensity=energy/(e*length);
  magneton=energy/B;
  
  velocity=momentum/mass;
  density=mass/(length*length*length);
  temperature=energy/kB;
  
  so_strength=energy*length;
  so3_strength=energy*length*length*length;
  g_factor=(4*m_e)/mass;
}

void units_class::list(FILE* file, int where_to)
{
  sprintf(buffer,"\n\nnatural constants values used (SI units):\n\n");
  sprintf(buffer,"hbar=%e\nelectron mass=%e\nelectron charge=%e\n",hbar,m_e,e);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"vacuum permitivitty=%e\nvacuum permeabilitty=%e\nBoltzmann constant=%e\n",eps,mi,kB);
  splachni(logg,buffer,where_to);

  sprintf(buffer,"\nnatural dimensionfull constants: given by 1 parameter - length set by hand\nother are electron mass, el. charge, hbar\n\n");
  splachni(logg,buffer,where_to);
  sprintf(buffer,"mass=%e [kg]\n",mass);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"length=%e [m]\n",length);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"energy=%e [J], %e [meV]\n",energy,energy/meV);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"omega=%e [1/s]\n",omega);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"B=%e [T]\n",B);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"A=%e [T m]\n",A);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"wavevector=%e [1/m]\n",wavevector);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"velocity=%e [m/s]\n",velocity);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"density=%e [kg/m^3]\n",density);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"temperature=%e [K]\n",temperature);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"momentum=%e [kg m/s]\n",momentum);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"el_intensity=%e [V/m]\n",el_intensity);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"magneton=%e [J/T], %e [meV/T]\n",magneton,magneton/meV);
  splachni(logg,buffer,where_to);
  sprintf(buffer,"so_strength=%e [J m], %e [meV A]\n",so_strength,so_strength/(meV*1e-10));
  splachni(logg,buffer,where_to);
  sprintf(buffer,"g-factor=%e [1]\n",g_factor);
  splachni(logg,buffer,where_to);
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!***********************************parameters
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
class parameters_class parameters;


//physical parameters: names
#define N_p 62
#define N_c 28

const char* N_physical[N_p]={\
  "B","phi_B","theta_B","gx","gy",\
  "gz","lv","lB","growth","br",\
  "lbr","dress","ldress","dress3","d",\
  "d_x","d_y","phi_d","width_z","T",\
  "Bosc","phi_Bosc","theta_Bosc","Eosc","phi_Eosc",\
  "ring","eps0","density","cl","ct",\
  "def","piezo","flat_x","flat_y","exch",\
  "xMn","ALHO","nuclear_state","fermi_e","Ebias",\
  "m_eff","imp_ext","imp_str","imp_num", "wg_x",\
  "wg_y","wg_width","wg_phi","wg_V","wg_length",\
  "dot_offset","def2","vf","Dg","Di",\
  "Dr", "gr_dd", "gr_w", "gr_r0", "gr_V0",\
  "gr_U","pairing"\
};

//physical parameters: description
const char* D_physical[N_p]={\
  "strenght of magnetic field","direction of magnetic field: angle phi",\
  "direction of magnetic field: angle theta","g-factor along x",\
  "g-factor along y","g-factor along z",\
  "confining length","effective length (derived)",\
  "growth direction","strength of Bychkov-Rashba",\
  "spin orbit length (BR) (derived)","strength of linear Dresselhaus",\
  "spin orbit length (D) (derived)","strength of cubic Dresselhaus",\
  "interdot distance", "interdot distance along x (derived)",\
  "interdot distance along y (derived)", "interdot orientation: angle phi",\
  "width of the well (relaxation)", "temperature (relaxation) do not use - is always zero",\
  "strength of oscillating magnetic field","orientation of oscillating magnetic field: angle phi",\
  "orientation of oscillating magnetic field: angle theta","strength of oscillating electic field",\
  "orientation of oscillating electric field: angle phi","ring radius",\
  "static dielectric constant","material density",\
  "longitudinal sound velocity","transversal sound velocity",\
  "deformation potential", "piezoelectric potential",\
  "region extent along x (from -value to value) if geometry==2", "region extent along y (from -value to value) if geometry==2",\
  "coupling to impurity spins = beta * N_0","impurity spins concentration",\
  "anisotropic LHO potential? (y/n=1/0)","type of nuclear ensemble for 1e diag",\
  "fermi energy", "additional potential bias electric field",\
  "effective mass", "impurities: extent",\
  "impurities: strength", "impurities: number",\
  "waveguide shift x", "waveguide shift y",\
  "waveguide width", "waveguide azimuth (heading out)",\
  "barrier height", "barrier length",\
  "potential offset of the dot wrt the lead","deformation shear transversal phonon potential",\
  "fermi velocity in graphene", "graphene gap",\
  "graphene intrinsic soc", "graphene rashba soc",\
  "graphene mass dd", "graphene mass w",\
  "graphene mass r0", "graphene mass V0",\
  "graphene potential U", "pairing gap Delta0"\
};

//physical parameters: unit names
const char* U_physical[N_p]={\
  "Tesla", "pi", "pi", "-", "-", "-", "nm", "nm", "-", "meVA",\
  "nm", "meVA", "nm", "eVA^3", "nm", "nm", "nm", "pi", "nm", "K",\
  "Tesla", "pi", "pi", "V/m", "pi", "nm", "-", "kg/m3", "m/s", "m/s",\
  "eV", "eV/m", "nm", "nm", "meV", "-", "-", "-", "meV", "V/m",\
  "m_e", "nm", "meV", "-", "nm", "nm", "nm", "pi", "meV", "nm",\
  "meV","eV","m/s","meV","meV","meV","nm","nm","nm","meV",\
  "meV", "meV"
  };

//computational parameters: names
const char* N_computational[N_c]={\
  "extent","grid_step","length_x","length_y","dLx0",\
  "dLy0","dLx","dLy","dim_x","dim_y",\
  "dim","N","dimstep","precision","peierls",\
  "eig","eigmax","geometry","leads","threshold",\
  "sorting_method","eig2e","selected2e","wg_extent","CE_exp",\
  "Ix","sets","system",\
};

//computational parameters: descriptions
const char* D_computational[N_c]={\
  "expand region [lB]","grid step [lB]",\
  "range for x axis [nm] (derived)","range for y axes [nm] (derived)",\
  "x shift of lower left corner (asked for)[nm]","y shift of lower left corner (asked for) [nm]",\
  "x shift of lower left corner (final) [nm] (derived)","y shift of lower left corner (final) [nm] (derived)",\
  "number of grid points on x axis plus 1 (derived)","number of grid points on y axis plus 1 (derived)",\
  "resolution (derived)","number of steps in Richardson extrapolation",\
  "step in Richardson","precision of differential operators",\
  "are operators symmetrized by Peierls phase?","number of eigenstates to get",\
  "maximal number of states in sorting routines", "(0/1/2/3 = 2D/1D ring/flat/waveguide)",\
  "how to treat leads include(# of)/exclude/ignore (>0,0,-1)","threshold in assigning parity (allowed unprecision)",\
  "sorting method for 1e states (0 no sorting; 1 I; 2 Ix,Iy; 3 lz; 4 distance; 5 FD states; *-1 spin)", "number of two electron eigenvectors",\
  "number of single electron states to build 2e basis", "extent of the waveguide",\
  "expansion in the Fourier Transforms (zero padding)", "what to integrate in spectral density (0-3 for debug, 4-default)",\
  "number of copies (=2S+1) where S is spin","system: g/A/S graphene/GaAs/Silicon",\
};

const char* U_computational[N_c]={\
  "-", "-", "nm", "nm", "nm", "nm", "nm", "nm", "-", "-",\
  "-", "-", "-", "-", "-", "-", "-", "-", "-", "-",\
  "-", "-","-","nm","-","-","-","-",\
  };

const char* noname="-";

parameters_class::parameters_class() {//constructor

  for (int i=0;i<N_p;i++) {
    names[i]=N_physical[i];
    descriptions[i]=D_physical[i];
    unit_names[i]=U_physical[i];
    unit_values[i]=1;
    values[i]=0;
  }
  for (int i=0;i<N_c;i++) {
    names[i+N_p]=N_computational[i];
    descriptions[i+N_p]=D_computational[i];
    unit_names[i+N_p]=U_computational[i];
    unit_values[i+N_p]=1;
    values[i+N_p]=0;
  }

  unit_values[B]=units.B;
  unit_values[theta_B]=unit_values[phi_B]=1/M_PI;
  unit_values[gx]=unit_values[gy]=unit_values[gz]=units.g_factor;
  unit_values[6]=unit_values[7]=units.length/1e-9;
  unit_values[8]=1;
  unit_values[9]=unit_values[11]=units.so_strength/(1e-3*units.e*1e-10);
  unit_values[10]=unit_values[12]=units.length/1e-9;
  unit_values[13]=units.so3_strength/(units.e*1e-30);
  unit_values[14]=unit_values[15]=unit_values[16]=units.length/1e-9;
  unit_values[17]=1/M_PI;
  unit_values[18]=units.length/1e-9;
  unit_values[20]=units.B;
  unit_values[21]=unit_values[22]=1/M_PI;
  unit_values[Eosc]=units.el_intensity;
  unit_values[24]=1/M_PI;
  unit_values[25]=units.length/1e-9;
  unit_values[27]=units.mass/(units.length*units.length*units.length);
  unit_values[28]=unit_values[29]=units.velocity;
  unit_values[30]=unit_values[31]=units.energy/units.e;
  unit_values[32]=unit_values[33]=units.length/1e-9;
  unit_values[exch]=units.energy/units.meV;
  unit_values[xMn]=1;
  unit_values[ALHO]=1;
  unit_values[nuclear_state]=1;
  unit_values[37]=unit_values[38]=units.energy/(units.e*1e-3);
  unit_values[Ebias]=units.el_intensity;
  unit_values[40]=1;
  unit_values[imp_ext]=units.length/units.nm;
  unit_values[imp_str]=units.energy/units.meV;
  unit_values[43]=1;
  unit_values[wg_x]=units.length/units.nm;
  unit_values[wg_y]=units.length/units.nm;
  unit_values[wg_width]=units.length/units.nm;
  unit_values[wg_phi]=1/M_PI;
  unit_values[wg_V]=units.energy/units.meV;
  unit_values[wg_length]=units.length/units.nm;
  unit_values[dot_offset]=units.energy/units.meV;  
  unit_values[def2]=units.energy/units.e;
  unit_values[vf]=1;
  unit_values[Dg]=unit_values[Di]=unit_values[Dr]=units.energy/units.meV;
  unit_values[gr_dd]=unit_values[gr_w]=unit_values[gr_r0]=units.length/units.nm;
  unit_values[gr_V0]=unit_values[gr_U]=unit_values[pairing]=units.energy/units.meV;

  unit_values[extent]=1;
  unit_values[grid_step]=1;  
  unit_values[length_x]=unit_values[length_y]=units.length/1e-9;
  unit_values[dLx0]=unit_values[dLy0]=unit_values[dLx]=unit_values[dLy]=units.length/1e-9;
  unit_values[wg_extent]=units.length/1e-9;
  unit_values[CE_exp]=1;
  unit_values[Ix]=1;
  unit_values[sets]=1;
  unit_values[system]=1;
  
  values[wg_extent]=0;
  values[growth]=1;
  values[precision]=reg1::eight;
  values[peierls]=true;

  //default values
  set("width_z", 11, "nm");
  set("N", 1, "-");
  set("dimstep", 1, "-");
  set("geometry", 0, "-");  
  set("ALHO", 0, "-");

  set("eps0",12.9,"-");               //static dielectric constant
  set("density",5300,"kg/m3");        //density
  set("cl",5290,"m/s");               //longitudinal velocity
  set("ct",2480,"m/s");               //transversal velocity
  set("def",7,"eV");                  //deformation potential
  set("def2",0,"eV");                 //deformation shear transversal phonon potential
  set("piezo",1.4e9,"eV/m");          //piezoelectric constant
  set("m_eff",0.067,"m_e");           //effective mass
  
  set("gx",-0.44,"-");               //g factors
  set("gy",-0.44,"-");
  set("gz",-0.44,"-");
  
  set("br",3.3*0,"meVA");		//bychkov-rashba
  set("dress",4.5*0,"meVA");		//dresselhaus
  set("dress3",27.5*0,"eVA^3");	//dresselhaus cubic 
  
  set("Bosc",0,"Tesla");	        //oscillating B
  set("phi_Bosc",0,"pi");
  set("theta_Bosc",0,"pi");
  
  set("Eosc",0,"V/m");		//oscillating E
  set("phi_Eosc",0,"pi");
  
  set("leads",0,"-"); //no leads
  
  set("dim_x",0.0,"-");
  set("dim_y",0.0,"-");
  set("geometry",0.0,"-");
  set("ring",0.0,"nm");
  set("Ebias",0.0,"V/m");

  set("phi_B",0,"pi");
  set("theta_B",0,"pi");

  set("phi_d",0.0,"pi");  
  
  set("grid_step",0.25,"-");
  set("extent",5,"-");	
  set("dot_offset",0,"meV");
  
  set("nuclear_state",0,"-");

  set("CE_exp",1,"-");  
  set("Ix",4,"-");
  set("sets",2,"-");
  set("system",'A',"-");

}

double parameters_class::read(const char* name, const char* unit)
{
  int index=name2index(name);
  double value=values[index]*unit_values[index];
  if (strcmp(unit_names[index],unit)!=0) printf("wrong unit (%s vs correct: %s) in parameter %s\n",unit,unit_names[index],name);
  return(value);
}

void parameters_class::set(const char* name, const double value,const char* unit)
{
  int index=name2index(name);
  values[index]=value/unit_values[index];
  if (unit==0) return;
  if (strcmp(unit_names[index],unit)!=0) printf("wrong unit (%s vs correct: %s) in parameter %s\n",unit,unit_names[index], name);
  return;
}

int parameters_class::name2index(const char* name) 
{
  int i;
  for (i=0;i<N_p+N_c;i++) {
    if (strcmp(name, names[i])==0) break;
  }
  if (i==N_p+N_c) {
    printf("unknown parameter %s used in parameters_class::name2index\n",name);
    exit(1);
  }
  return(i);
}

void parameters_class::recompute()
{
  prec B=parameters.read("B","Tesla");
  prec theta_B=parameters.read("theta_B","pi")*M_PI;
  prec lv=parameters.read("lv","nm")*units.nm;

  //effective length - in nanometers
  prec lB=1/sqrt(sqrt(pow(B*cos(theta_B)*units.e/2/units.hbar,2)+1/pow(lv,4)))/units.nm;
  parameters.set("lB",lB,"nm");

  //spin-orbit lengths
  prec m_eff=parameters.read("m_eff","m_e");
  prec br=parameters.read("br","meVA")*(units.meV*1e-10);
  prec dress=parameters.read("dress","meVA")*(units.meV*1e-10);
  prec lbr=pow(units.hbar,2)/(2*units.m_e*m_eff*br);
  if (br==0) lbr=0;
  prec ldress=pow(units.hbar,2)/(2*units.m_e*m_eff*dress);
  if (dress==0) ldress=0;
  set("lbr",lbr/units.nm,"nm");
  set("ldress",ldress/units.nm,"nm");

  //eigestates/geometry/system parameters sanity check
  int geometry=(int) parameters.read("geometry","-");
  int leads=(int) parameters.read("leads","-");
  if (geometry!=2) {//if not transport -- eigenstates will be computed
    int eig=(int) parameters.read("eig","-");
    int eigmax=(int) parameters.read("eigmax","-");
    if (eigmax<eig) {
      message(logg,"value of eigmax must be at least eig\neigmax reset to eig\n",4);
      parameters.set("eigmax",eig,"-");
    }
  }
  char system=(char) parameters.read("system","-");
  fprintf(logg,"apart from spatial %i degree(s) of freedom considered: ",dof.Ndof);
  if (dof_prototype::spin<dof.Ndof) fprintf(logg,"spin ");
  if (dof_prototype::sigma<dof.Ndof) fprintf(logg,"sigma ");
  if (dof_prototype::tau<dof.Ndof) fprintf(logg,"tau ");
  if (dof_prototype::eta<dof.Ndof) fprintf(logg,"eta ");
  if (dof_prototype::kappa<dof.Ndof) fprintf(logg,"kappa ");
  fprintf(logg,". The system tag is '%c' and number of components is %i.\n",system, (int) parameters.read("sets","-"));
  //if (system=='A') if (dof_prototype::sigma<=dof.Ndof) {printf("sigma (%i vs end=%i) not defined in GaAs [in recompute()]: edit classes.h dof_prototype\n", dof_prototype::sigma, dof.Ndof);exit(1);}
  
  //set up the geometry
  if (leads>0) lB=lv/units.nm;//!!!NO grid adjustements for the magnetic field length if there are leads present
  prec grid_step=parameters.read("grid_step","-")*lB;
  prec extent=parameters.read("extent","-")*lB;
  
  if (system=='g') {
    extent=(2*parameters.read("extent","-")+3)*parameters.read("gr_dd","nm");
    grid_step=parameters.read("grid_step","-")*parameters.read("lv","nm");
  }
  
  prec length_x, length_y;
  int dim_x, dim_y;
 
  switch (geometry) {//specify grid step and number of points along x and y axis
    case (3) : //the waveguide is present in a 2d potential
    case (0) : {//2d double dot or ring
      bool ALHO=(bool) parameters.read("ALHO","-");
      prec d_x=0,d_y=0;
      prec d=fabs(parameters.read("d","nm"));
      prec phi_d=parameters.read("phi_d","pi")*M_PI;
      prec ring=fabs(parameters.read("ring","nm"));
      if (ALHO && ring!=0) {printf("ring with ALHO not implemented [in recompute()]\n"); exit(1);}
      if (ALHO) {
        d_x=d_y=0;
        prec epar=extent*(1+d/lv*units.nm);
        prec eper=extent/(1+d/lv*units.nm);
        length_x=max(fabs(epar*cos(phi_d)), fabs(eper*sin(phi_d)));
        length_y=max(fabs(epar*sin(phi_d)), fabs(eper*cos(phi_d)));
        //fprintf(logg,"ring pars: d=%e, epar=%e, eper=%e, l_x=%e, l_y=%e\n",d, epar,eper,length_x, length_y);
      } else {
        d_x=cos(phi_d)*d;
        d_y=sin(phi_d)*d;
	if (system=='g') ring+=fabs(parameters.read("gr_r0","nm"));
        length_x=extent+fabs(d_x)+fabs(ring);
        length_y=extent+fabs(d_y)+fabs(ring);
      }
      parameters.set("d_x",d_x,"nm");
      parameters.set("d_y",d_y,"nm");


      if (geometry==3) {//expand the space to have the waveguide inside
        prec wg_phi=parameters.read("wg_phi","pi")*M_PI;
        prec wg_extent=parameters.read("wg_extent","nm");
        prec wg_width=parameters.read("wg_width","nm");
        prec x1,y1,x2,y2;
        waveguide2dot(wg_extent*units.nm/units.length,0,x1,y1);//one farthers corner
        waveguide2dot(wg_extent*units.nm/units.length,wg_width*units.nm/units.length,x2,y2);//second farthers corner
        prec xm=max(abs(x1),abs(x2))/units.nm*units.length;
        prec ym=max(abs(y1),abs(y2))/units.nm*units.length;
        length_x=max(length_x,xm);
        length_y=max(length_y,ym);
      }

      //adjust to integer number of points
      if (grid_step<=0) {message("grid_step zero / negative in parameters.recompute()\n",1+4);exit(1);}
      dim_x=(int) my_round(2*length_x/grid_step);
      dim_y=(int) my_round(2*length_y/grid_step);

      if (d>0 || true) {//always hit the middle potential peak in a double dot
        if (dim_y%2!=0) dim_y++;
        if (dim_x%2!=0) dim_x++;
      }
      length_x=grid_step*dim_x/2.0;
      length_y=grid_step*dim_y/2.0;
      break;
    }
    case (2) : {//flat region with leads
      //dim_x=(int) my_round(2*extent/grid_step);
      //dim_y=(int) my_round(2*extent/grid_step);
      dim_x=max((int) my_round(2*parameters.read("flat_x","nm")/grid_step),1);
      dim_y=max((int) my_round(2*parameters.read("flat_y","nm")/grid_step),1);
      length_x=grid_step*dim_x/2.0;
      length_y=grid_step*dim_y/2.0;

      prec leads=abs(parameters.read("leads","-"));
      prec xexp=0, yexp=0;
      for (int l=1;l<=leads;l++) {
        switch (lead[l].heading) {
          case 1 : {xexp=max(xexp,grid_step*4-lead[l].x*units.length/units.nm-length_x);break;}
          case 2 : {xexp=max(xexp,grid_step*4+lead[l].x*units.length/units.nm-length_x);break;}
          case 3 : {yexp=max(yexp,grid_step*4-lead[l].y*units.length/units.nm-length_y);break;}
          case 4 : {yexp=max(yexp,grid_step*4+lead[l].y*units.length/units.nm-length_y);break;}
        }
      }
      printf("adding %e and %e [nm] to x and y\n",xexp,yexp);
      dim_x+=(int) my_round(2*xexp/grid_step);
      dim_y+=(int) my_round(2*yexp/grid_step);

      if (dim_y%2!=0) dim_y++;
      if (dim_x%2!=0) dim_x++;
      length_x=grid_step*dim_x/2.0;
      length_y=grid_step*dim_y/2.0;
      break;
    }
    case (1) : {//1d ring (circle with radius ring represented by periodic boundary conditions)
      length_x=M_PI*ring;
      dim_x=(int) my_round(ring/grid_step);
      grid_step=2*length_x/dim_x;
      dim_y=2;
      length_y=dim_y*grid_step/2;
      break;
    }
    default : {printf("wrong geometry(=%i) in parameters.recompute()\n",geometry);exit(1);}
  }

  parameters.set("dim_x",dim_x,"-");
  parameters.set("dim_y",dim_y,"-");
  parameters.set("length_x",length_x,"nm");
  parameters.set("length_y",length_y,"nm");
  parameters.set("grid_step",2*length_x/dim_x/lB,"-");
}

void parameters_class::list(FILE* file, int where_to)
{
  message(file,"\n\nchosen parameters:\n",where_to);

  message(file,"physical parameters:\n",where_to);
  message(file,"tag\t\tvalue\t\tdimensionless\tunit\tdescription\n",where_to);
  for (int i=0;i<N_p;i++) {
    if (strlen(names[i])>7) 
      sprintf(buffer,"%s\t%.4e\t%.4e\t%s\t%s\n",names[i],values[i]*unit_values[i],values[i],unit_names[i],descriptions[i]);
    else sprintf(buffer,"%s\t\t%.4e\t%.4e\t%s\t%s\n",names[i],values[i]*unit_values[i],values[i],unit_names[i],descriptions[i]);
    splachni(file,buffer,where_to);
  }

  message(file,"\ncomputational parameters:\n",where_to);
  message(file,"tag\t\tvalue\t\tdimensionless\tdescription\n",where_to);
  for (int i=N_p;i<N_p+N_c;i++) {
    if (strlen(names[i])>7) 
      sprintf(buffer,"%s\t%.4e\t%.4e\t%s\n",names[i],values[i]*unit_values[i],values[i],descriptions[i]);
    else 
      sprintf(buffer,"%s\t\t%.4e\t%.4e\t%s\n",names[i],values[i]*unit_values[i],values[i],descriptions[i]);
    splachni(file,buffer,where_to);
  }
  message(file,"\n\n",where_to);
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!***********************************state_info
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//collected info about the states
//after diagonalizing in region r1 
//and writing eigenvalues and vectors into vysl
//this class computes, stores and gives various parameters of a state

const char* CH_units[]={\
  "hbar/2","hbar/2","hbar/2","hbar","nm^2",\
  "nm^-2","-","-","-","hbar/2",\
  "-","-","hbar","meV/T","-",\
  "-","meV","meV","meV","meV",\
  "hbar/2","meV","meV","meV","meV",\
  "meV","hbar/2","hbar/2","-","nm",\
  "nm"\
};
  

const char* CH_tags[]={\
  "sx","sy","sz","lz","r2","dd","Ix","Iy","I","s",\
  "np","nm","Lz","mi","l","n","E0","E","corr0","corr",\
  "S","energy","energyI","corr_r0","corr_r1","aux","g-jz","sigmaz","Isigmaz","r_x",\
  "r_y"\
};

const char* CH_descriptions[]={\
  "spin along x","spin along y","spin along z",\
  "orbital momentum along z","r squared","laplace",\
  "inversion along x","inversion along y","complete inversion",\
  "spin along B","Fock label(+)","Fock label(-)",\
  "orbital momentum along z (kinematic)","magnetic moment along z",\
  "Fock label(l)","Fock label(n)","Fock energy(0)",\
  "Fock energy(1)","so correction(0)","so correction(1)",\
  "spin(+-1)","energy","energy-imaginary part","ring correction(0)",\
  "ring correction(1)","auxiliary","total \"J\" for graphene","sigmaz for graphene","sigmaz tensor times I", "dipole moment along x",\
  "dipole moment along y"\
};


state_info::state_info(reg1* r1_in, ret_from_arpack_zn* vysl_in, int eig, int N)
{
  r1=r1_in;
  vysl=vysl_in;
  unsorted=eig;
  sorted=N;
  length=31;
  
  x=new prec*[unsorted];
  actual=new bool*[unsorted];
  
  unsorted2s=new int[unsorted];
  sorted2un=new int[sorted];
    
  rates=new prec[4*eig*eig];
  
  for (int s=0;s<sorted;s++) sorted2un[s]=-1;
  
  for (int u=0;u<unsorted;u++) {
    unsorted2s[u]=sorted2un[u]=u;
    x[u]=new prec[length];
    actual[u]=new bool[length];
    for (int t=0;t<length;t++) {
      actual[u][t]=false;
    }
  }
  
  tags=new const char*[length];
  names=new const char*[length];
  unit_names=new const char*[length];
  show=new bool[length];
  unit_values=new prec[length];


  for (int t=0;t<length;t++) {
    tags[t]=CH_tags[t];
    names[t]=CH_descriptions[t];
    unit_names[t]=CH_units[t];
    show[t]=false;
    unit_values[t]=1;
  }
  unit_values[0]=unit_values[1]=unit_values[2]=2/units.hbar;
  unit_values[3]=1/units.hbar;
  unit_values[4]=1/pow(units.nm,2);
  unit_values[5]=pow(units.nm,2);
  unit_values[9]=2/units.hbar;
  unit_values[12]=1/units.hbar;
  unit_values[13]=1/units.meV; //meV / Tesla
  unit_values[16]=unit_values[17]=unit_values[18]=unit_values[19]=1/units.meV;
  unit_values[20]=2/units.hbar;
  unit_values[21]=unit_values[22]=unit_values[23]=unit_values[24]=unit_values[25]=1/units.meV;
  unit_values[26]=2/units.hbar;
  unit_values[27]=2/units.hbar;
  unit_values[29]=unit_values[30]=1/units.nm;
}

state_info::~state_info()
{
  for (int u=0;u<unsorted;u++) {
    delete x[u];
    delete actual[u];
}
  
  delete x;
  delete actual;
  delete unsorted2s;
  delete sorted2un;
  delete rates;
  
  delete tags;
  delete unit_names;
  delete names;
  delete show;
  delete unit_values;
}

void state_info::show_set(const char* what)
{
  if (strcmp(what,"clear")==0) {
    for (int t=0;t<length;t++) show[t]=false;
    return;
}
  show[find_tag(what)]=true;
}
  

prec state_info::read(int u, const char *tag) {
  int t=find_tag(tag);
  return(read(u,t));
}

prec state_info::read(int u, int t) {
  if (u<0) return(-1);
  compute(u,t);
  return(x[u][t]*unit_values[t]);
}

void state_info::tag_info(int t, const char* &tag, const char* & name, const char* &unit)
{
  tag=tags[t];
  name=names[t];
  unit=unit_names[t];
}

void state_info::list_tags(FILE* file, int where_to)
{
  for (int t=0;t<length;t++) {
    sprintf(buffer,"%2i %s\t%s\t%s\n",t,tags[t],unit_names[t],names[t]);
    splachni(file,buffer, where_to);
  }
}

int state_info::find_tag(const char* tag)
{
  int t;
  for (t=0;t<length;t++) {
    if (strcmp(tag,tags[t])==0) break;
    if (t==length-1) {
      sprintf(buffer,"tag %s not found in state_info::find_tag\n",tag);
      splachni(logg,buffer,1+4);
      exit(1);
}
}
  return(t);
}
//return energy difference in units of frequency, that is (E_i-E_j)/hbar
prec state_info::omega(int i, int j)
{
  if (i==j) return(0);
  int si=sorted2un[i];
  int sj=sorted2un[j];
  if (si==-1 || sj==-1) return(0);
  int t=find_tag("energy");
  prec ei=read(si,t);
  prec ej=read(sj,t);
  prec omega=(ei-ej)/unit_values[t];
  return(omega);
}

void state_info::new_states() 
{
  for (int u=0;u<unsorted;u++) 
    for (int t=0;t<length;t++) actual[u][t]=false;
}
    
void state_info::compute(int u, int t)
{
  if (u<0) return;
  if (actual[u][t]) return;
  
  actual[u][t]=true;
  
  content* pntr=vysl->eigenvecs+u*r1->Nw*r1->sets;
  
  prec B=parameters.read("B","Tesla");
  prec theta_B=parameters.read("theta_B","pi")*M_PI;
  prec phi_B=parameters.read("phi_B","pi")*M_PI;
  prec lB=parameters.read("lB","nm")*units.nm;
  prec gz=parameters.read("gz","-");  
  prec m_eff=parameters.read("m_eff","m_e");
  int geometry=(int) parameters.read("geometry","-");
  
  //first compute dependent parameters, if any, then the value (if not already)
  switch (t) {
    case 0: {//sx
      if (dof_prototype::spin>=dof.Ndof) x[u][0]=cos(phi_B)*sin(theta_B);
      else {
	//PREVIOUS
	//prec prev=r1->mean_value(pntr,reg1::sigmax).real();
	dof.reset("operators");
	dof.operatorC[dof_prototype::spin]=dof_prototype::paulix;
	dof.convert("OP2M");
	content aux=0;
	for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++)
          if (abs(dof.M[i*dof.length+j])!=0) aux+=r1->braket(pntr,i,pntr,j)*dof.M[i*dof.length+j];
        x[u][0]=aux.real();
	//if (abs(x[u][0]-prev)>-1e-8) fprintf(logg,"difference of old (%e) and dof(%e) %e in sigmax exp. value\n",prev, x[u][0], abs(x[u][0]-prev));
      }
      x[u][0]*=units.hbar/2;
      break;
    }
    case 1: {//sy
      if (dof_prototype::spin>=dof.Ndof) x[u][1]=sin(phi_B)*sin(theta_B);
      else {
	//PREVIOUS
	//prec prev=r1->mean_value(pntr,reg1::sigmay).real();
	dof.reset("operators");
	dof.operatorC[dof_prototype::spin]=dof_prototype::pauliy;
	dof.convert("OP2M");
	content aux=0;
	for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++) 
	  if (abs(dof.M[i*dof.length+j])!=0) aux+=r1->braket(pntr,i,pntr,j)*dof.M[i*dof.length+j];
	x[u][1]=aux.real();
	//if (abs(x[u][1]-prev)>-1e-8) fprintf(logg,"difference of old (%e) and dof(%e) %e in sigmay exp. value\n",prev, x[u][1],abs(x[u][1]-prev));
      }
      x[u][1]*=units.hbar/2;
      break;
    }
    case 2: {//sz
      if (dof_prototype::spin>=dof.Ndof) x[u][2]=cos(theta_B);
      else {
	//PREVIOUS
	//prec prev=r1->mean_value(pntr,reg1::sigmaz).real();
        dof.reset("operators");
	dof.operatorC[dof_prototype::spin]=dof_prototype::pauliz;
	dof.convert("OP2M");
	content aux=0;
	for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++) 
          if (abs(dof.M[i*dof.length+j])!=0) aux+=r1->braket(pntr,i,pntr,j)*dof.M[i*dof.length+j];
        x[u][2]=aux.real();
	//if (abs(x[u][2]-prev)>1e-8) fprintf(logg,"difference of old (%e) and dof(%e) %e in sigmaz exp. value\n",prev, x[u][2],abs(x[u][2]-prev));
      }
      x[u][2]*=units.hbar/2;
      break;
    }
    case 3: {//lz
      if (geometry==1) x[u][3]=r1->mean_value(pntr,reg1::lz1Dx).real();
      else x[u][3]=r1->mean_value(pntr,reg1::lz).real();
      x[u][3]*=units.hbar;
      break;}
    case 4: {//r2
      x[u][4]=r1->mean_value(pntr,reg1::r2).real()*pow(units.length,2);break;}
    case 5: {//dd
      x[u][5]=r1->mean_value(pntr,reg1::dd_mv).real()/pow(units.length,2);break;}
    case 6: {//Ix
      x[u][6]=r1->mean_value(pntr,reg1::invx_mv).real();break;}
    case 7: {//Iy
      x[u][7]=r1->mean_value(pntr,reg1::invy_mv).real();break;}
    case 8: {//I
      x[u][8]=r1->mean_value(pntr,reg1::inv_mv).real();break;}
    case 9: {//s
      compute(u,0);compute(u,1);compute(u,2);
      x[u][9]=cos(theta_B)*x[u][2]+sin(theta_B)*(x[u][0]*cos(phi_B)+x[u][1]*sin(phi_B));
      break;
    }
    case 10: {//np
      compute(u,3);compute(u,4);compute(u,5);
      x[u][10]=(x[u][4]/(lB*lB)-lB*lB*x[u][5]-2.0+2.0*x[u][3]/units.hbar)/4.0;
      break;
    }
    case 11: {//nm
      compute(u,3);compute(u,4);compute(u,5);
      x[u][11]=(x[u][4]/(lB*lB)-lB*lB*x[u][5]-2.0-2.0*x[u][3]/units.hbar)/4.0;
      break;
    }
    case 12: {//Lz
      compute(u,3);compute(u,4);
      x[u][12]=x[u][3]+B*cos(theta_B)/2*x[u][4]*units.e;
      break;
    } 
    case 13: {//mi
      compute(u,3);compute(u,2);
      x[u][13]=(-x[u][3]/m_eff-gz*x[u][2])*units.e/(2*units.m_e);
      break;
    }
    case 14: {//l
      compute(u,10);compute(u,11);
      x[u][14]=x[u][10]-x[u][11];
      break;
    }
    case 15: {//n
      compute(u,10);compute(u,11);
      x[u][15]=min(x[u][10],x[u][11]);
      break;
    }
    case 16: {//E0
      compute(u,14);compute(u,15);compute(u,20);
      x[u][16]=energy(x[u][15],x[u][14],(int) (2*x[u][20]/units.hbar))*units.energy;
      break;
    }
    case 17: {//E
      compute(u,14);compute(u,15);compute(u,20);
      x[u][17]=energy((int) rint(x[u][15]),(int) rint(x[u][14]),(int) (2*x[u][20]/units.hbar))*units.energy;
      break;
    }
    case 18: {//corr0
      compute(u,14);compute(u,15);compute(u,20);
      x[u][18]=correction(31,rint(x[u][15]),rint(x[u][14]),(int) (x[u][20]));
      break;
    }
    case 19: {//corr
      compute(u,14);compute(u,15);compute(u,20);
      x[u][19]=correction(31,x[u][15],x[u][14],(int) (x[u][20]));
      break;
    }
    case 20: {//S
      compute(u,9);
      if (x[u][9]>0) x[u][20]=units.hbar/2; else x[u][20]=-units.hbar/2;
      break;
    }
    case 21: {//energy
      x[u][21]=vysl->eigenvals[u].real()*units.energy;
      break;
    }
    case 22: {//energy-imaginary part
      x[u][22]=vysl->eigenvals[u].imag()*units.energy;
      break;
    }
    case 23: {//ring correction(0)
      compute(u,3);
      x[u][23]=ring_correction(x[u][3],0)*units.energy;
      break;
    }
    case 24: {//ring correction(1)
      x[u][24]=ring_correction(0,1)*units.energy;
      break;
    }
    case 25: {//aux: E+ring corr(0)+ring corr(1)
      compute(u,21);compute(u,23);compute(u,24);
      x[u][25]=(x[u][21]+x[u][23]+x[u][24])*units.energy;
      break;
    }   
    case 26: {//total jz - graphene: lz+sigmaz +sz (the first in hbar, the latter in hbar/2)
      compute(u,3);compute(u,2);compute(u,27);
      x[u][26]=x[u][3]+x[u][2]+x[u][27];
      break;
    }
    case 27: {//sigmaz - graphene: sigmaz
      if (dof_prototype::sigma>=dof.Ndof) {
	fprintf(logg,"sigma(%i) dof is not among considered ones(%i)\n",dof_prototype::sigma,dof.Ndof); exit(1);
      }
      dof.reset("operators");
      dof.operatorC[dof_prototype::sigma]=dof_prototype::pauliz;
      dof.convert("OP2M");
      content aux=0;
      for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++) 
	if (abs(dof.M[i*dof.length+j])!=0) aux+=r1->braket(pntr,i,pntr,j)*dof.M[i*dof.length+j];
      x[u][27]=aux.real()*units.hbar/2;
      break;
    }
    case 28: {//sigmaz tensor product total Inversion - graphene: 
      if (dof_prototype::sigma>=dof.Ndof) {
	fprintf(logg,"sigma(%i) dof is not among considered ones(%i)\n",dof_prototype::sigma,dof.Ndof); exit(1);
      }
      dof.reset("operators");
      dof.operatorC[dof_prototype::sigma]=dof_prototype::pauliz;
      dof.convert("OP2M");
      content aux=0;
      content* temp_mv=new content[r1->Nw*r1->sets];
      r1->op_ket(pntr,temp_mv,reg1::inv_mv,reg1::set);
      for (int i=0;i<dof.length;i++) for (int j=0;j<dof.length;j++) 
	if (abs(dof.M[i*dof.length+j])!=0) aux+=r1->braket(pntr,i,temp_mv,j)*dof.M[i*dof.length+j];
      x[u][28]=aux.real();
      delete temp_mv;
      break;
    }
    case 29: {//r_x
      x[u][29]=r1->mean_value(pntr,reg1::r_x).real()*units.length;break;}
    case 30: {//r_y
      x[u][30]=r1->mean_value(pntr,reg1::r_y).real()*units.length;break;}
    default : {
      message("wrong tag in state_info::compute\n",1+4);
      exit(1);
    }
  }
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!Impurities
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
class impurities_class impurities;

impurities_class::impurities_class() 
{
  N=0;
}

impurities_class::~impurities_class() 
{
  free_arrays();
}
void impurities_class::free_arrays()
{
  if (N<1) return;
  delete strength;
  delete extent;
  delete loc_x;
  delete loc_y;
}

void impurities_class::generate()
{
  int N_new=(int) parameters.read("imp_num","-");
  if (N!=N_new) {
    free_arrays();
    N=N_new;
    strength=new prec[N];
    extent=new prec[N];
    loc_x=new prec[N];
    loc_y=new prec[N];
  }
  int geometry=(int) parameters.read("geometry","-");
  prec x_low, x_len, y_low, y_len;

  if (geometry==2) {
    x_low=-parameters.read("flat_x","nm");
    x_len=-2*x_low;
    y_low=-parameters.read("flat_y","nm");
    y_len=-2*y_low;
  }
  else {
    x_low=parameters.read("dLx","nm")-parameters.read("length_x","nm");
    x_len=2*parameters.read("length_x","nm");
    y_low=parameters.read("dLy","nm")-parameters.read("length_y","nm");
    y_len=2*parameters.read("length_y","nm");
  }

  prec E=parameters.read("imp_str","meV")*units.meV/units.energy;
  prec l=parameters.read("imp_ext","nm")*units.nm/units.length;
  average=0;

  for (int i=0;i<N;i++) {
    strength[i]=(generuj()*0+1-0.5)*E;
    extent[i]=(generuj()*0+1-0.5)*l;
    average+=2*M_PI*extent[i]*extent[i]*strength[i];
    loc_x[i]=(x_low+generuj()*x_len)*units.nm/units.length;
    loc_y[i]=(y_low+generuj()*y_len)*units.nm/units.length;
  }
  //loc_x[0]=loc_y[1]=40*units.nm/units.length;
  //loc_x[1]=loc_y[0]=0;
}
  prec impurities_class::potential(prec x, prec y)
{
  prec pot=0;
  //printf("total %i impurities:\n",N);
  for (int i=0;i<N;i++) {
    prec dist=pow(loc_x[i]-x,2)+pow(loc_y[i]-y,2);
    prec act=strength[i]*exp(-pow(dist/extent[i],2)/2);
    pot+=act;
    //printf("impurity %i: loc(%e, %e), strength %e, extent %e\n\at distance to actual point %e the strength is %e\n", i, loc_x[i], loc_y[i], strength[i], extent[i], dist, act);
  }
  //exit(1);
  return(pot-average);
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!spin impurities
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
class spin_impurities_class spin_impurities;
#define IMP 2

spin_impurities_class::spin_impurities_class()
{
  N=0;
  for (int i=0;i<3;i++) {
    I[i]=new prec[N];
    I_old[i]=new prec[N];
  }
}

void spin_impurities_class::allocate()
{
  for (int i=0;i<3;i++) {
    delete I_old[i];
    I_old[i]=new prec[N];
    delete I[i];
    I[i]=new prec[N];
    for (int n=0;n<N;n++) I[i][n]=0;
  }
}


spin_impurities_class::~spin_impurities_class()
{
  for (int i=0;i<3;i++) {
    delete I_old[i];
    delete I[i];
  }
}

void spin_impurities_class::reinitialize(class region* r1, int type, prec aux)
{
  double phi_B=parameters.read("phi_B","pi")*M_PI;
  double theta_B=parameters.read("theta_B","pi")*M_PI;
  double B[3]={sin(theta_B)*cos(phi_B),sin(theta_B)*sin(phi_B),cos(theta_B)};
  
  Arho=-parameters.read("exch","meV")*units.meV/units.energy*parameters.read("xMn", "-")*(5.0/2)*(1.0/2);
  spin_impurities_class::r1=r1;
  if (r1==0) N=0; else N=r1->Nw;
  allocate();

  for (int n=0;n<N;n++) {
    //in case you would need the coordinates
    prec x,y;
    r1->s2coordinate(n,x,y);
    prec norm=0;
    for (int i=0;i<3;i++) {
      switch (type) {
        case 0 : {I[i][n]=0; break;}						//put to zero
        case 1 : ;
        case 5 : {if (i==2) I[i][n]=2*(generuj()>0.5)-1; else I[i][n]=0; break;}//a random unit vector along z
        case 2 : {if (i==2) I[i][n]=1; else I[i][n]=0; break;}			//a unit vector along z
        case 3 : {I[i][n]=generuj()*2-1; break;}				// a random unit vector
        case 4 : {								//gaussian distribution with unit dispersion
          I[i][n]=sqrt(2)*erfinv(generuj())/sqrt(3);
          if (generuj()<0.5) I[i][n]*=-1;
          break;
        }        
        case 6 : {if (i==2) {//dipole along x axis
	  if (abs(x)>1e-8) {if (x>0) I[i][n]=+1; else I[i][n]=-1;} else I[i][n]=0; break;}
	}
        case 7 : {//core-halo with boundary at aux
          if (i==2) {
            double r=sqrt(x*x+y*y)*units.length/units.nm;
            if (r>aux) I[i][n]=+1; else I[i][n]=-1;
          } 
          break;
        }
        case 8 : {I[i][n]=-B[i];break;}//uniform antialigned with the field
        case 9 : {I[i][n]=+B[i];break;}//uniform aligned with the field
        case 10 : {//smooth dipole aligned with the field
	  double x0=20*units.nm/units.length; 
	  double scale=100*units.nm/units.length;
	  double profile=1/(1+(x-x0)*(x-x0)/(scale*scale)) -1/(1+(x+x0)*(x+x0)/(scale*scale));
	  if (-x>0) profile=1; else profile=-1;
	  if (x<x0 && x>-x0) profile=1; else profile=0;
	  I[i][n]=+B[i]*profile;
	  break;
	}
        default :{fprintf(logg,"wrong type (%i) in spin_impurities_reinitialize\n",type);exit(1);}
      }
      norm+=I[i][n]*I[i][n];
    }
    if (type==3) for (int i=0;i<3;i++) I[i][n]/=sqrt(norm);
  } 
  double sum2=0, sum[3]={0,0,0};
  for (int n=0;n<N;n++) {
    for (int i=0;i<3;i++) sum[i]+=I[i][n];
    sum2+=I[0][n]*I[0][n]+I[1][n]*I[1][n]+I[2][n]*I[2][n];
  }
  double ave=0;
  for (int i=0;i<3;i++) {
    sum[i]/=max(N,1);
    ave+=sum[i]*sum[i];
  }
  sum2/=max(N,1);
  #if IMP>0
  const char* mess[]={"zero", "unit random along z", "unit vector along z", "unit random", "gaussian", "unit random along z", "perfect dipole along x","core halo","uniform antialigned", "uniform aligned","smooth dipole along the field"};
  fprintf(logg,"%i impurity spins were initialized to %s (average [%e,%e,%e], dispersion:%e)\n",N, mess[type], sum[0], sum[1], sum[2], sum2-ave*ave);
  #endif
}

void spin_impurities_class::position(int n, prec& x, prec& y)
{
  r1->s2coordinate(n,x,y);
}

void spin_impurities_class::update_overlap()
{
  prec norm_new=0, norm_old=0;
  overlap=0;
  for (int i=0;i<3;i++) for (int n=0;n<N;n++) {
    norm_old+=I_old[i][n]*I_old[i][n];
    norm_new+=I[i][n]*I[i][n];
    overlap+=I[i][n]*I_old[i][n];
  }
  if (norm_new*norm_old!=0) overlap/=sqrt(norm_new*norm_old); else overlap=1.0;
}


void spin_impurities_class::set_impurity(int n, prec * v, prec kBT, prec x)
{
  prec sum=0;
  for (int i=0;i<3;i++) sum+=v[i]*v[i];
  if (sum>0) for (int i=0;i<3;i++) v[i]/=sqrt(sum); else x=0;
  if (kBT!=0 && sum!=0) x*=brillouin(5.0/2, sqrt(sum)/kBT);

  prec overlap=0;
  for (int i=0;i<3;i++) overlap+=-v[i]*I[i][n]*x;
  #if IMP>2
  if (overlap<0) fprintf(logg,"impurity %i reset antialigned to (%.2e, %.2e, %.2e) with brillouin factor %e from (%.2e, %.2e, %.2e) resulting in <old|new>=%e\n", n, v[0], v[1], v[2], x, I[0][n], I[1][n], I[2][n], overlap);
  #endif
  for (int i=0;i<3;i++) {
    I_old[i][n]=I[i][n];
    I[i][n]=-v[i]*x;
    prec ratio=I[i][n]/I_old[i][n];
    if (I_old[i][n]==0) ratio=0;
    if (ratio>0.9 && ratio<1) I[i][n]*=pow(ratio,under);
    if (I[i][n]<-1) I[i][n]=-1;
    if (I[i][n]>1) I[i][n]=1;    
    //I[i][n]=-v[i]*x*(1-under)+ under*I_old[i][n];
  }
}

void spin_impurities_class::set_indexes(int si, int sj)
{
  spin_impurities_class::si=si;
  spin_impurities_class::sj=sj;
}

content spin_impurities_class::potential(prec x, prec y)
{
  if (N==0) return(content(0,0));
  int n;
  bool in=r1->coordinate2s(x,y,n);
  #if IMP>2
  fprintf(logg,"coors (%e, %e) give seq label %i, corresponding impurity spin is (%+.2e, %+.2e, %+.2e) (sigma indexes si=%i sj=%i and Arho=%e)\n", x, y, n, I[0][n], I[1][n], I[2][n], si, sj, Arho);
  #endif
  if (!in) {
    printf("point (%e, %e) claimed not inside by region->coordinate2s - can not happen here in spin impurity potential\n",x,y); 
    exit(1);
    return(0);
  }
  if (n>=N) {
    printf("seq label %i larger than or equal to N=%i in spin impurities potential\n", n, N); 
    exit(1);
    return(0);
  }
  if (si==0 && sj==0) return(I[2][n]*Arho);
  if (si==1 && sj==1) return(-I[2][n]*Arho);
  if (si==0 && sj==1) return(I[0][n]*Arho+content(0,1)*I[1][n]*Arho);
  if (si==1 && sj==0) return(I[0][n]*Arho-content(0,1)*I[1][n]*Arho);
  printf("can not get here in spin impurity potential\n"); 
  exit(1);
}

int spin_impurities_class::circle(prec r, prec h, int i)
{
  prec phi=0;
  int sign=0, candidate=0;
  //int pntr=0;
  do {
    int n;
    if (r1->coordinate2s(r*cos(phi), r*sin(phi), n)) {
      if (I[i][n]>0) sign=+1; else sign=-1;
      if (candidate==0) candidate=sign;
      if (sign!=candidate) break;
    }
    phi+=h/r;
    //fprintf(logg,"circle step %i: phi=%e, sign=%i, candidate=%i\n",pntr++, phi, sign, candidate);
  }  while (phi>0 && phi<2*M_PI);
  if (sign==candidate) return(sign); else return(0);
}

void spin_impurities_class::distribution_characteristics(two_electron_state_prototype* state2e)
{
  prec dx[3], dy[3], Qxx[3], Qyy[3], Qxy[3], tot0=0, d0=0, Q0s=0, Q0a=0, Qa0[3], avenon=0;
  
  for (int i=0;i<3;i++) dx[i]=dy[i]=Qxx[i]=Qxy[i]=Qyy[i]=ave[i]=0;
  double* w=new double[N*4];
  int nmax=0;
  for (int n=0;n<N;n++)  {
    if (state2e!=0) state2e->spinatn(w+4*n, n); else w[4*n+0]=1;
    if (w[4*n]>w[4*nmax]) nmax=n;
  }
  prec x,y;
  r1->s2coordinate(nmax,x,y);
  Rmax=sqrt(x*x+y*y)*units.length/units.nm;

  for (int n=0;n<N;n++) {
    prec x,y;
    r1->s2coordinate(n,x,y);
    for (int i=0;i<3;i++) {
      prec s=I[i][n]*w[4*n];
      ave[i]+=s;
      dx[i]+=x*s;
      dy[i]+=y*s;
      Qxx[i]+=x*x*s;
      Qxy[i]+=x*y*s;
      Qyy[i]+=y*y*s;
      if (i==2) {
        avenon+=I[i][n];
        tot0+=abs(s);
        Q0a+=abs(s*x*y);
        Q0s+=abs(s)*(x*x+y*y)/2;
      }
    }
  }
  for (int i=0;i<3;i++) {
    ave[i]/=tot0;
    d[i]=sqrt(dx[i]*dx[i]+dy[i]*dy[i]);
    dphi[i]=atan(dx[i]/dy[i]);
    prec Qsaux=(Qxx[i]+Qyy[i])/2;
    prec Qaaux=(Qxx[i]-Qyy[i])/2;
    Qs[i]=Qsaux/Q0s;
    Qa[i]=sqrt(Qaaux*Qaaux+Qxy[i]*Qxy[i])/Q0a;
    Qphi[i]=atan(Qxy[i]/Qaaux);
    prec Q1=Qsaux+Qa[i]*Q0a;
    prec Q2=Qsaux-Qa[i]*Q0a;
    Qsph[i]=abs(Qa[i]*Q0a/Qsaux);
    if (Qsaux==0) Qsph[i]=0;
#if IMP > 1
    if (i==2) {
      fprintf(logg,"distributon characteristics:\naverage: ave=%e (non-weigthed=%e), N=%i\ndipole: dx=%e, dy=%e, d=%e dphi=%e\n",ave[i], avenon/N,N,dx[i],dy[i],d[i],dphi[i]);
      fprintf(logg,"Quadrupole: Qxx=%e, Qxy=%e, Qyy=%e, (Qxx+Qyy)/2=%e vs max=%e, Qs=%e, (Qxx-Qyy)/2=%e vs max=%e, Qa=%e\n", Qxx[i], Qxy[i], Qyy[i], Qsaux, Q0s, Qs[i], Qaaux, Q0a, Qa[i]);
      fprintf(logg,"Quadrupole - eigenvalues: Q1=%e, Q2=%e (sphericity: %e)\n", Q1, Q2, Qsph[i]);
    }
#endif
  }
  for (int n=0;n<N;n++) {
    prec x,y;
    r1->s2coordinate(n,x,y);
    d0+=abs(I[2][n]*w[4*n])*abs(x*dx[2]+y*dy[2])/d[2];
  }
#if IMP > 1
  fprintf(logg,"d0=%e\n",d0);
#endif
  d[2]/=d0;
  delete w;
  
  //rotationally symmetric borders -- only for z-component
  prec hx, hy, r=0;
  prec length_x=parameters.read("length_x","nm")*units.nm/units.length;
  r1->give_par(hx, hy);
  borders=0;
  int s_prev=0;
  do {
    r+=hx;
    if (!inside(r,0) || r>length_x) break;
    int s=circle(r,hx,2);
    if (s==0) continue;
    if (s_prev==0) {s_prev=s;continue;}
    if (s==s_prev) continue;
    borders++;
#if IMP > 1
    fprintf(logg, "borders increased to %i at r=%e (s=%i, s_prev=%i)\n", borders, r, s, s_prev);
#endif
    s_prev=s;
  } while (true);
  
  //decide about the pattern
  pattern=0;
  if (abs(abs(ave[2])-1)<2e-1) pattern=1;
  else if (abs(d[2]-1)<2e-1) pattern=2;
  if (borders==2 && abs(d[2])<2e-1) pattern=3;
  
#if IMP > 1
  fprintf(logg,"resulting pattern: %i\n", pattern);  
#endif
  
  update_overlap();
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Fourier in 2D
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


//constructor
Fourier2D::Fourier2D()
{
Initialized=false;
}

//destructor
Fourier2D::~Fourier2D()
{
  if (Initialized) {
    fftw_destroy_plan(plan);
    delete fftw_in;
    delete fftw_out;
    delete attenuation_factors;
    Initialized=false;
  }
}

content* Fourier2D::address_in(int Nx_now, int Ny_now)
{
  
  if ((Nx!=Nx_now || Ny!=Ny_now) && Initialized) {
    //dimensions changed - reinitialize
    fftw_free(fftw_in);
    fftw_free(fftw_out);
    delete attenuation_factors;
    fftw_destroy_plan(plan);
    Initialized=false;
  }
  
  if (!Initialized) {
    Nx=Nx_now;Ny=Ny_now;
    fftw_in=(content*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny);
    fftw_out=(content*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny);
    plan=fftw_plan_dft_2d(Nx,Ny,(fftw_complex*) fftw_in,(fftw_complex*) fftw_out, FFTW_FORWARD, FFTW_ESTIMATE);
    attenuation_factors=new prec[Nx+Ny];
    fill_attenuators();
    Initialized=true;
  }

  return(fftw_in);
}

content* Fourier2D::execute()
{
  fftw_execute(plan);
  return(fftw_out);
}

//fourier to consecutive coordinate (along which=1,2 axis)
int Fourier2D::four2cons(int n, int which)
{
  int N;
  switch (which) {
    case 1: N=Nx;break;
    case 2: N=Ny;break;
    default: {message("wrong choice in Fourier2D::four2cons\n",1+4);exit(1);}
  }
  //int max=N/2+N%2;
  //int min=-N/2; 
  //if (n<min) {n+=N;}
  //if (n>max) {n-=N;}
  if (n<0) n+=N;
  if (n<0 || n>=N) {
    n=0;message("coordinate off in Fourier2D::four2cons\n",1+4);exit(1);
  }
  return(n);
}

//consecutive to fourier coordinate (along which=1,2 axis)
int Fourier2D::cons2four(int i, int which)
{
  int N;
  switch (which) {
    case 1: N=Nx;break;
    case 2: N=Ny;break;
    default: {message("wrong choice in Fourier2D::cons2four\n",1+4);exit(1);}
  }
  int max=N/2+N%2;
  int min=-N/2; 

  if (i>max) i-=N;
  if (i<min || i>max) {
    i=0;message("coordinate off in Fourier2D::cons2four\n",1+4);exit(1);
  }
  return(i);
}

void Fourier2D::fill_attenuators()
{
  for (int i=0;i<Nx;i++) {
    int n=cons2four(i,1);
    attenuation_factors[i]=attenuation(((prec) n)/( (prec) Nx)*2*M_PI);
  }
  for (int j=0;j<Ny;j++) {
    int m=cons2four(j,2);
    attenuation_factors[j+Nx]=attenuation(((prec) m)/( (prec) Ny)*2*M_PI);
  }
}

prec Fourier2D::attenuation(prec theta) 
{
  prec res=1;
  if (theta!=0) res=(3.0-4*cos(theta)+cos(2*theta))*(6+theta*theta)/(3*theta*theta*theta*theta);
  return(res);
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!                         parity
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define PARITY0
//#define PARITY1

//sorting_class sorting;

void sorting_class::init(int s, int u, int mve)
{
  if (mve<1 || mve>3) {printf("invalid parameter: max_valid_entries (%i) in sorting_class::init\n",mve);mve=1;}
  max_valid_entries=mve;
  deallocate();
  
  sorted=s;
  unsorted=u;
  
#ifdef PARITY0
  sprintf(buffer,"allocating and initializing sorting class:\n%i unsorted states will be sorted to %i positions\n",unsorted, sorted);
  message(logg,buffer,4);
#endif

  
  level=new level_s[sorted];
  new_state=new level_u[unsorted];
  
  gap_value=new prec[sorted*sorted];
  gap_param=new prec[sorted*sorted];
  
  initialized=false;
  too_many=false;
  for (int s=0;s<sorted;s++) {
    level[s].valid_entries=0;
    level[s].visible=false;
    level[s].s2u=-1;
    level[s].parity=0;
    for (int i=0;i<3;i++) level[s].vals[i]=0;
    for (int s2=0;s2<sorted;s2++) {
      int label=s*sorted+s2;
      gap_value[label]=gap_param[label]=-1;
    }
  }
  dog_tag=new int[sorted*sorted*2];
  allocated=true;
}

void sorting_class::save(FILE* file)
{
  if (!allocated) {
    message("sorting_class asked to save, but is not allocated yet\n",1+4);
    exit(1);
  }
  
  fprintf(file,"seen,initialized,too_many,sorted,unsorted\n");
  fprintf(file,"%i\t%i\t%i\t%i\t%i\n",seen,initialized,too_many,sorted,unsorted);
  
  fprintf(file,"pars:\n");
  fprintf(file,"%e\t%e\t%e\n",pars[0],pars[1],pars[2]);
  
  for (int s=0;s<sorted;s++) {
    fprintf(file,"sorted_level_%i\n",s);
    fprintf(file,"parity,valid_entries,visible,s2u\n");
    fprintf(file,"%i\t%i\t%i\t%i\n",level[s].parity,level[s].valid_entries,level[s].visible,level[s].s2u);
    fprintf(file,"energy_values:\n");
    for (int i=0;i<level[s].valid_entries;i++)  fprintf(file,"%e\t",level[s].vals[i]);
    if (level[s].valid_entries>0) fprintf(file,"\n");
  }
}

void sorting_class::load(FILE* file)
{
  if (!allocated) {
    message("sorting_class asked to load, but is not allocated yet\n",1+4);
    exit(1);
  }
  
  char* aux_string=new char[1000];
  
  fscanf(file,"%s",aux_string);
  #ifdef PARITY0
  sprintf(buffer,"reading global info: %s\n",aux_string);splachni(logg,buffer,4);
  #endif  
  int seenr,initializedr,too_manyr;
  int sortedr,unsortedr,methodr;
  fscanf(file,"%d\t%d\t%d\t%d\t%d\t%d\n",&seenr,&initializedr,&too_manyr,&sortedr,&unsortedr,&methodr);
  if (sorted!=sortedr) {printf("saved and actual parameter differ: sorted (%i,%i)\n",sortedr,sorted);exit(1);} 
  if (unsorted!=unsortedr) {printf("saved and actual parameter differ: unsorted (%i,%i)\n",unsortedr,unsorted);exit(1);} 
  
  seen=seenr;
  initialized=initializedr;
  too_many=too_manyr;
  #ifdef PARITY0
  sprintf(buffer,"seen set to %i\n",seen);splachni(logg,buffer,4);
  sprintf(buffer,"initialized set to %i\n",initialized);splachni(logg,buffer,4);
  sprintf(buffer,"too_many set to %i\n",too_many);splachni(logg,buffer,4);
  #endif
  float parsr[3];
  fscanf(file,"%s",aux_string);
  #ifdef PARITY0
  sprintf(buffer,"reading parameters: %s\n",aux_string);splachni(logg,buffer,4);
  #endif
  fscanf(file,"%e\t%e\t%e\n",parsr,parsr+1,parsr+2);
  for (int i=0;i<3;i++) {
    pars[i]=parsr[i];
    #ifdef PARITY0
    sprintf(buffer,"par[%i]=%e\t",i,pars[i]);splachni(logg,buffer,4);
    #endif
  }
  #ifdef PARITY0
  message(logg,"\n",4);
  #endif
  
  
  int parityr,valid_entriesr,visibler,s2ur;
  for (int s=0;s<sorted;s++) {
    fscanf(file,"%s",aux_string);
    #ifdef PARITY1
    sprintf(buffer,"reading state info: %s\n",aux_string);splachni(logg,buffer,4);
    #endif
    fscanf(file,"%s",aux_string);
    #ifdef PARITY1
    sprintf(buffer,"reading: %s\n",aux_string);splachni(logg,buffer,4);
    #endif
    fscanf(file,"%d\t%d\t%d\t%d\n",&parityr,&valid_entriesr,&visibler,&s2ur);
    
    #ifdef PARITY1
    sprintf(buffer,"state set: parity:%i,valid_entries:%i,visible:%i,s2u:%i\n",parityr,valid_entriesr,visibler,s2ur); splachni(logg,buffer,4);
    #endif
    level[s].parity=parityr;
    level[s].valid_entries=valid_entriesr;
    level[s].visible=visibler;
    level[s].s2u=s2ur;
    
    fscanf(file,"%s",aux_string);
    #ifdef PARITY1
    sprintf(buffer,"reading: %s\n",aux_string);splachni(logg,buffer,4);
    #endif    
    float valsr;
    for (int i=0;i<level[s].valid_entries;i++) {
      fscanf(file,"%e\t",&valsr);
      level[s].vals[i]=valsr;
      #ifdef PARITY1
      sprintf(buffer,"state set: vals[%i]:%e\n",i,level[s].vals[i]); splachni(logg,buffer,4);
      #endif
    }
  }
}

void sorting_class::deallocate()
{

#ifdef PARITY0
  sprintf(buffer,"deallocating sorting class (allocated:%i)\n",allocated);
  message(logg,buffer,4);
#endif

 if (allocated) {
   delete new_state;
   delete level;
   dog.reset();
   delete dog_tag;
   allocated=false;
   delete gap_value;
   delete gap_param;
 }
} 
  
void sorting_class::update_history()
{
  for (int s=0;s<sorted;s++) {
    if (level[s].visible) {
      for (int i=2;i>0;i--) level[s].vals[i]=level[s].vals[i-1];
      level[s].vals[0]=new_state[level[s].s2u].value;
      level[s].valid_entries=min(max_valid_entries,level[s].valid_entries+1);
    }
  }
  pars[2]=pars[1];pars[1]=pars[0];pars[0]=act_parameter;
}

void sorting_class::forecast()
{

  for (int s=0;s<sorted;s++) {
    if (level[s].visible==false) continue;
    prec extrap;
    prec coefs[3];
    #ifdef PARITY1
    sprintf(buffer,"forecasting position of state %i...",s);
    message(logg,buffer,4);
    #endif
    if (level[s].valid_entries>1) {//if at least two values
      extrapolate_coefs(pars,level[s].vals,level[s].valid_entries-1,coefs);
      extrap=coefs[1]*act_parameter+coefs[0];
      if (level[s].valid_entries>2) extrap+=coefs[2]*act_parameter*act_parameter;
    }
    else {//only one value - that is used as approximation
      extrap=level[s].vals[0];
    }
    level[s].forecast=extrap;
    #ifdef PARITY1
    sprintf(buffer,"at parameter %e\n",extrap);
    message(logg,buffer,4);
    #endif
  }
}


void sorting_class::inspect_new_1e_states(state_info* states)
{
   prec threshold=parameters.read("threshold","-");
   int method=(int) parameters.read("sorting_method","-");
   all_definite=true;                 //definite parity for each state is assigned
   //states=s;
   for (int u=0;u<unsorted;u++) {
     new_state[u].value=states->read(u,"energy");
     new_state[u].parity=give_parity(states,u,threshold,method);
     if (new_state[u].parity==0) all_definite=false;
     new_state[u].assigned=false;
   }
   if (method==0) is_method_zero=true; else is_method_zero=false; 
}

void sorting_class::inspect_new_2e_states(level_u* levels)
{
   int method=(int) parameters.read("sorting_method","-");
   all_definite=true;                 //definite parity for each state is assigned
   for (int u=0;u<unsorted;u++) {
     new_state[u].value=levels[u].value;
     new_state[u].parity=levels[u].parity;
     //if (new_state[u].parity==0) all_definite=false;
     new_state[u].assigned=false;
   }
   if (method==0) is_method_zero=true; else is_method_zero=false; 
}


bool sorting_class::do_I_sort()
{
  if (too_many) {
    message(logg,"number of levels exceeded maximum - keeping previous order!",1+2);
    return(false);
  }
  if (is_method_zero) {return(false);}

  //we initialize the levels if not yet inicialized and not ended due to small eigmax
  if (! initialized && ! too_many) {
    initialized=true;
    seen=unsorted-1;
    
    for (int s=0;s<unsorted;s++) {
      level[s].visible=true;
      level[s].s2u=s;
      new_state[s].u2s=s;
      level[s].valid_entries=0;
      level[s].parity=new_state[s].parity;
    }
    if (! all_definite) { 
      message(logg,"no definite parity -> not inicialized!\n",1+2);
      initialized=false;
    }
    update_history();
    return(false);
  }
  else if (! all_definite) message("no definite parity -> sorting according to the distance in this round!\n");
  for (int s=0;s<sorted;s++) level[s].assigned=false;
  return(true);
}

int sorting_class::closest_sorted2unsorted(int u, bool ignore_parity) 
{
  prec min=-1;
  int min_s=-1;
  bool found=false;
  for (int s=0;s<sorted;s++) {
    if (! level[s].visible || level[s].assigned) continue;
    if (new_state[u].parity!=0 && !ignore_parity) 
      if (level[s].parity!=new_state[u].parity) continue;
    prec dist=abs(new_state[u].value-level[s].forecast);
    #ifdef PARITY1
    sprintf(buffer,"%i. unsorted a free position (%i. sorted) with the same parity (%i) found ... trying (act. distance and minimum:%.4e %.4e)\n", u,s,level[s].parity,dist,min);
    splachni(logg,buffer,4);
    #endif
    if (dist<min || !found) {
      min=dist;
      min_s=s;
      found=true;
    }
  }
  #ifdef PARITY1
  if (found) {
    sprintf(buffer,"for unsorted %i (value=%.4e),\nclosest value %.4e (distance=%.4e) found\n", u,new_state[u].value,level[min_s].forecast,min);
    splachni(logg,buffer,4);
  }
  #endif
  return(min_s);
}

void sorting_class::state_found(int u, int s) 
{
  if (s!=-1) {
    #ifdef PARITY1
    sprintf(buffer,"%i. unsorted assigned %i. sorted position\n",u,s);
    splachni(logg,buffer,4);
    #endif
    level[s].assigned=true;
    level[s].s2u=u;
  }
  else {
    #ifdef PARITY1
    message(logg,"a new state found\n",4);
    #endif
    seen++; 
    s=seen;
    if (s>=sorted) {
      message("parity.seen exceeds maximum in parity!!! updating positions stopped\n");
      too_many=true;
    }
    else {
      level[s].parity=new_state[u].parity;
      level[s].s2u=u;
      level[s].assigned=true;
      level[s].visible=true;
    }
  }
  new_state[u].u2s=s;
  new_state[u].assigned=true;
}

void sorting_class::copy_results(int * s2u, int* u2s)
{
  for (int s=0;s<sorted;s++) {
    if (!level[s].assigned) {
      level[s].visible=false;
      level[s].valid_entries=0;
      s2u[s]=-1;
    }
    else s2u[s]=level[s].s2u;
  }
   if (! too_many) for (int u=0;u<unsorted;u++) u2s[u]=new_state[u].u2s;
  
  #ifdef PARITY0
  char num[100000],vis[100000],pos[100000],parch[100000],aux[10];

  strcpy(num,"sorted  :");
  strcpy(vis,"visible :");
  strcpy(pos,"unsorted:");
  strcpy(parch,"parity  :");


  for (int s=0;s<min(seen+1,sorted);s++) {
    sprintf(aux,"%3i ",s);
    strcat(num,aux);
    sprintf(aux,"%3i ",level[s].visible);
    strcat(vis,aux);
    sprintf(aux,"%3i ",level[s].s2u);
    strcat(pos,aux);
    sprintf(aux,"%3i ",level[s].parity);
    strcat(parch,aux);
  }
  splachni(logg,num,4);splachni(logg,"\n",4);
  splachni(logg,vis,4);splachni(logg,"\n",4);
  splachni(logg,pos,4);splachni(logg,"\n",4);
  splachni(logg,parch,4);splachni(logg,"\n",4);
  #endif //PARITY0
}

void sorting_class::sort(prec p, int* s2u, int* u2s)
{
   act_parameter=p;

  #ifdef PARITY0
    sprintf(buffer,"all have definite parity:%i\n",all_definite);
    message(logg,buffer,4);
  #endif
  
  #ifdef PARITY0
    sprintf(buffer,"deciding whether to sort...");
    message(logg,buffer,4);
  #endif
  bool sorting=do_I_sort();
  #ifdef PARITY0
    sprintf(buffer,"%i\n",sorting);
    message(logg,buffer,4);
  #endif

  if (sorting) {
  
    #ifdef PARITY0
      sprintf(buffer,"forecasting states positions\n");
      message(logg,buffer,4);
    #endif
    forecast();
    for (int u=0;u<unsorted;u++) {
      #ifdef PARITY1
      sprintf(buffer,"unsorted state %i assigning at...",u);
      message(logg,buffer,4);
      #endif
      int s=closest_sorted2unsorted(u,!all_definite);
      #ifdef PARITY1
      sprintf(buffer,"position %i\n",s);
      message(logg,buffer,4);
      #endif
      state_found(u,s);
    }
    #ifdef PARITY0
    sprintf(buffer,"copying results to states\n");
    message(logg,buffer,4);
    #endif
    copy_results(s2u, u2s);
    #ifdef PARITY0
    sprintf(buffer,"updating history\n");
    message(logg,buffer,4);
    #endif
    update_history();
    update_gaps(p);
  }
}

prec sorting_class::dog_init(prec s_max, prec s_min, prec f)
{
  step_max=s_max; step_min=s_min; expand=f;
  dog_N=0;
  for (int j=0;j<2*sorted*sorted;j++) dog_tag[j]=-1;
  dog.reset();
  return(step_min);
}

void sorting_class::dogs_fill(prec p)
{
  //go through all pairs of visible states and update stack for each
  for (int s1=0;s1<sorted;s1++) {
    int u1=level[s1].s2u;
    for (int s2=s1+1;s2<sorted;s2++) {
      int u2=level[s2].s2u;
      int label=2*(s1*sorted+s2);

      //if any of the two is not visible, discard the stack
      if (u1==-1 || u2==-1) {
        if (dog_tag[label]>-1) dog.reset(dog_tag[label]);
        dog_tag[label]=-1;
        if (dog_tag[label+1]>-1) dog.reset(dog_tag[label+1]);
        dog_tag[label+1]=-1;
        continue;
      }

      //both are visible -- if a stack does not exist, create one
      prec val=new_state[u1].value-new_state[u2].value;
      if (dog_tag[label]==-1) {
        dog_tag[label]=dog_N;
        dog_N++;
      }
      if (dog_tag[label+1]==-1) {
        dog_tag[label+1]=dog_N;
        dog_N++;
      }

      //update the stack
      char* tag=new char[200];
      sprintf(tag,"states %2i and %2i (sorted)",s1,s2);
      dog.another_in(dog_tag[label],p,val,3,1,strcat(tag,"  - closest (quadratic approx)"));
      dog.another_in(dog_tag[label+1],p,val,2,1,strcat(tag,"  - closest (linear)"));
      delete tag;
    }
  }
}

void sorting_class::update_gaps(prec p)
{
  //go through all pairs of visible states and update gap for each
  for (int s1=0;s1<sorted;s1++) {
    int u1=level[s1].s2u;
    if (u1==-1) continue;
    for (int s2=0;s2<sorted;s2++) {
      int u2=level[s2].s2u;
      if (u2==-1) continue;
      
      //both are visible -- if a stack does not exist, create one
      int label=s1*sorted+s2;
      prec val=fabs(new_state[u1].value-new_state[u2].value);
      
      if (val<gap_value[label] || gap_value[label]<0) {
        gap_value[label]=val;
        gap_param[label]=p;
      }
    }
  }
}

void sorting_class::gap(prec& val, prec& p, int s1, int s2) 
{
  int label=s1*sorted+s2;
  val=gap_value[label];
  p=gap_param[label];
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!                         dof prototype
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

dof_prototype dof;

dof_prototype::dof_prototype() 			//constructor
{
  Ndof=Ndof_max;
  length=1;
  for (int i=0;i<Ndof;i++) length*=2;
  indexC=new int[Ndof+1];
  operatorC=new int[Ndof+1];
  M=new complex<double>[length*length];
  M_sum=new complex<double>[length*length];
  braC=new int[Ndof+1];
  ketC=new int[Ndof+1];
  pntr=new int[Ndof+1];
}

dof_prototype::~dof_prototype() 			//destructor
{
  delete indexC;
  delete operatorC;
  delete M;
  delete M_sum;
  delete braC;
  delete ketC;
  delete pntr;
}

void dof_prototype::reset_Ndof(int Ndof_new)
{
    if (logg!=0) fprintf(logg,"changing the value of components in dof, from %i to %i (max %i)",Ndof, Ndof_new, Ndof_max);
    if (Ndof_new<Ndof_max) Ndof=Ndof_new;
    else {printf("failed in dof{prototype::reset_Ndof\n");exit(1);}
    length=1;
    for (int i=0;i<Ndof;i++) length*=2;
}

void dof_prototype::reset(const char* what)
{
  const char* choices[]={"operators", "matrix", "indexes", "all"};
  int choice;
  for (choice=0;choice<4;choice++) if (strcmp(what,choices[choice])==0) break;
  switch (choice) {
    case 0 : {for (int i=0;i<Ndof+1;i++) operatorC[i]=0;break;}
    case 1 : {for (int i=0;i<length*length;i++) M_sum[i]=0;break;}
    case 2 : {for (int i=0;i<Ndof+1;i++) indexC[i]=1;break;}
    case 3 : {reset("operators");reset("matrix");reset("indexes");break;}
    default : {printf("what %s not found in dof_prototype::reset\n",what); exit(1);}
  }
}

//sigma matrix 1,x,y,z,plus,minus can be written as  c1 |ket1> <bra1| + c2|ket2><bra2|
//where ci are complex coefficients, and ket/bra=1/-1 for up/down spinor
complex<double> dof_prototype::sigma2ketbra(int sigma, int term, int& ket, int& bra)
{
  if ((term<0 || term>1) || (sigma<0 || sigma>5)) {
    printf("sigma=%i/term=%i invalid in sigma2ketbra\n",sigma, term); exit(1);
  }
  complex<double> ii=complex<double>(0,1);
  switch (sigma) {
    case 0 : {if (term==0) {ket=1;bra=1;return(1.0);} else {ket=-1;bra=-1;return(1.0);} break;}
    case 1 : {if (term==0) {ket=1;bra=-1;return(1.0);} else {ket=-1;bra=1;return(1.0);} break;}
    case 2 : {if (term==0) {ket=1;bra=-1;return(-ii);} else {ket=-1;bra=1;return(ii);} break;}
    case 3 : {if (term==0) {ket=1;bra=1;return(1.0);} else {ket=-1;bra=-1;return(-1.0);} break;}
    case 4 : {if (term==1) {ket=1;bra=-1;return(1.0);} else {ket=1;bra=-1;return(0.0);} break;}
    case 5 : {if (term==1) {ket=-1;bra=1;return(1.0);} else {ket=-1;bra=1;return(0.0);} break;}
    default : return(0);
  }
}
void dof_prototype::convert_log(const char* mess, const char* what, complex<double> c)
{
  fprintf(logg,"convert stat (%s) with what=%s, coef=%e%+ei\n",mess, what, c.real(), c.imag());
  for (int i=0; i<Ndof+1; i++) fprintf(logg,"indexC[%i]=%i ",i,indexC[i]);
  fprintf(logg,"index=%i\n",index);
  for (int i=0; i<Ndof+1; i++) fprintf(logg,"operatorC[%i]=%i ",i,operatorC[i]);
  fprintf(logg,"\n");
  for (int i=0; i<length; i++) for (int j=0; j<length; j++) fprintf(logg,"element(%i,%i): M=%+e%+ei M_sum=%+e%+ei\n", i, j, M[i*length+j].real(), M[i*length+j].imag(), M_sum[i*length+j].real(),M_sum[i*length+j].imag());
}

void dof_prototype::convert(const char* what, complex<double> coef)
{
  //convert_log("start",what,coef);
  const char* choices[]={"C2I", "I2C", "OP2M"};
  int choice;
  for (choice=0;choice<3;choice++) if (strcmp(what,choices[choice])==0) break;
  switch (choice) {
    case 0 : {
      index=0;
      for (int i=0, p=1;i<Ndof+1;i++, p*=2) if (indexC[i]==-1) index+=p; else if (indexC[i]!=1) {printf("component index(%i)=%i out of scope in dof_prototype::convert\n",i,indexC[i]); exit(1);}
      break;
    }
    case 1 : {
      if (index<0 || index>=length) {printf("consequtive index %i out of scope(0,%i) in dof_vector_prototype::setI\n",index,length); exit(1);}
      for (int i=0;i<Ndof+1;i++,index/=2) if (index % 2) indexC[i]=-1; else indexC[i]=1;
      break;
    }
    case 2 : {
      if (Ndof<0) {printf("Ndof=%i in dof_prototype::convert\n",Ndof);exit(1);}
      for (int i=0;i<length*length;i++) M[i]=0;
      for (int i=0;i<Ndof+1;i++) pntr[i]=0;
      do {
        complex<double> c=1;
        for (int i=0;i<Ndof+1;i++) c*=sigma2ketbra(operatorC[i],pntr[i],ketC[i],braC[i]);
	      int ket=0, bra=0;
        for (int i=0, p=1;i<Ndof+1;i++, p*=2) {
          if (ketC[i]==-1) ket+=p;
          if (braC[i]==-1) bra+=p;
        }
        M[ket*length+bra]=c;
        pntr[0]++;
        for (int i=0;i<Ndof;i++) if (pntr[i]==2) {pntr[i]=0;pntr[i+1]++;}
      } while (pntr[Ndof]<2);
      if (abs(coef)!=0) for (int i=0;i<length*length;i++) M_sum[i]+=coef*M[i];
      break;
    }
    default : {printf("what %s not found in dof_prototype::convert\n",what); exit(1);}
  }
  //convert_log("end",what,coef);
}



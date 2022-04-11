#include "region2d.h"
//#define DEBUG1
//#define DEBUG2
//#define DEBUG3
#define REGION_CHECK
#define GRPH_CENTRALx

typedef double prec;
typedef complex<prec> content;


region::region(int Nx, prec Lx0, prec Lx, char bc_x, int Ny, prec Ly0, prec Ly, char bc_y, function_inside f, int sets)
{
  fprintf(logg,"starting with region...\n");
  region::sets=sets;

  func=f;
  act_prec=two;
  operators_sym1=false;
  operators_sym2=true;

  b=all_prec;
  if (bc_x=='p') pbc_x=true; else pbc_x=false;
  if (bc_y=='p') pbc_y=true; else pbc_y=false;
  if (pbc_x || pbc_y) {printf("periodic bc not implemented in region\n");exit(1);}
  
  neumann=dirichlet=periodic=false;
  if (bc_x=='d' && bc_y=='d') dirichlet=true;
  else if (bc_x=='n' && bc_y=='n') neumann=true;
  else {printf("boundary conditions x=%c y=%c not implemented in region\n",bc_x, bc_y);exit(1);}
  
  hx=Lx/Nx;	
  hy=Ly/Ny;
  
  if (abs(hx-hy)>1e-12*abs(hx+hy)) {printf("hx=%e and hy=%e different in region2d. This is probably not desired\n",hx, hy); exit(1);}
  
  kx=2*M_PI/((Nx+1)*hx);ky=2*M_PI/((Ny+1)*hy);
  tab_nx=Nx+1;tab_ny=Ny+1;
  
  Nx+=2*b-2;Ny+=2*b-2;
  if (pbc_x) Nx++;
  if (pbc_y) Ny++;
  
  region::Nx=Nx;
  Nx1=Nx+1; Nxb=Nx-b; Nxb1=Nxb+1;

  region::Ny=Ny;
  Ny1=Ny+1; Nyb=Ny-b;Nyb1=Nyb+1;

  Np=Nx1*Ny1;
  
  Nr_max=100;
  Nr=N_btp=0;
  relative_point_method='g';
  
  Nl=1;
  sectorized=false;

  region::Lx0=Lx0;
  region::Ly0=Ly0;
  region::Lx=Lx;
  region::Ly=Ly;

  inicialize_diferential_operators();
  //inicialize_diferential_operators_test();
  //exit(0);

  inicialize();
  fill_arrays();

  for (int i=0;i<all_op;i++) {
    remember[i]=false;
    for (int j=0;j<sets;j++) dif_op_filled[i][j]=false;
  }

  fprintf(logg,"prepared for work...\n");
  if (Nw<1) message("no points inside!\n",1+4);
}
//sets the remember flag for this operator
void region::remember_set(operators op,bool val)
{
  remember[op]=val;
}

//sets the symetrization on /works only with operate_tot_init/
void region::symmetrize_ops(bool val1, bool val2)
{
  operators_sym1=val1;
  operators_sym2=val2;
}

//sets the actual precision to be used
void region::act_prec_set(operators_precision val)
{
  act_prec=val;
}

region::~region()
{
  clean();
  fprintf(logg,"done with region...\n\n");
}

//allocate memory for all arrays
void region::inicialize1()
{
  xline=new xline_prototype[Nx1];
  yline=new yline_prototype[Ny1];

  for (int j=0;j<Ny1;j++) yline[j].fb=hy*(j-b+1-0.5*(pbc_y==true))+Ly0;
  for (int i=0;i<Nx1;i++) {
    xline[i].fb=hx*(i-b+1-0.5*(pbc_x==true))+Lx0;
    xline[i].np=0;
  }

  int hmn=0;			//points total
  for (int i=b;i<Nxb1;i++) {
    int hmn_line=0;		//points on this xline
    for (int j=b;j<Nyb1;j++) {
      //!!!remove points which have 3 outside neighbours - necessary for graphene
      //if(func(xd(i),yd(j))) {
      if (inside_with_neighbors(i,j)) {
        hmn++;hmn_line++;
      }
    }
    xline[i].np=hmn_line;
    if (hmn_line>0) {
      xline[i].active_line_point=new active_line_point_prototype[hmn_line];
    }
    else {
      xline[i].active_line_point=0;
    }
  }
  for (int i=b;i<Nxb1;i++) if (xline[i].np>0) {first_active_xline=i;break;}
  for (int i=Nxb1;i>=b;i--) if (xline[i].np>0) {last_active_xline=i;break;}
  lines=last_active_xline-first_active_xline+1;

  Nw=hmn;
}  
void region::inicialize2()
{
  for (int i=0;i<all_op;i++) {
    dif_op_vals[i]=new content*[sets];
    dif_op_filled[i]=new bool[sets];
  }
}

void region::inicialize3()
{
  sector=new sector_prototype[1];
  sector[0].imin=first_active_xline;
  sector[0].imax=last_active_xline;
  sector[0].dim=new int[lines];
  sector[0].offset=new int[lines];
  for (int i=sector[0].imin;i<=sector[0].imax;i++) {
    sector[0].dim[i-sector[0].imin]=xline[i].np;
    sector[0].offset[i-sector[0].imin]=0;
  }
}

void region::inicialize4()
{
  grid_point=new grid_point_prototype[Np];
  active_point=new active_point_prototype[Nw];
  sl_btp=new int[Nw];
  others=new content[Nw];
  temp_mv=new content[sets*Nw];
}
void region::inicialize5()
{
  multiples=new multiple**[sets*sets];
  for (int m=0;m<sets*sets;m++) {
    multiples[m]=new multiple*[Nw];
    for (int i=0;i<Nw;i++) multiples[m][i]=new multiple[8*b+1];
  }
}
void region::inicialize6()
{
  int aux=8*b+2;
  local_vector_prototype::max_length=aux;
  local_vector_labels=new int[sets*sets*Nw*aux*3];
  local_vector_coefs=new content[sets*sets*Nw*aux*3];
  local_vector=new local_vector_prototype*[sets*sets];
  for (int m=0;m<sets*sets;m++) {
    local_vector[m]=new local_vector_prototype[Nw];
    for (int i=0;i<Nw;i++) {
      local_vector[m][i].sl=local_vector_labels+((m*Nw+i)*3+0)*local_vector_prototype::max_length;
      local_vector[m][i].sl_r=local_vector_labels+((m*Nw+i)*3+1)*local_vector_prototype::max_length;
      local_vector[m][i].sl_s=local_vector_labels+((m*Nw+i)*3+2)*local_vector_prototype::max_length;
      local_vector[m][i].coef=local_vector_coefs+((m*Nw+i)*3+0)*local_vector_prototype::max_length;
      local_vector[m][i].coef_r=local_vector_coefs+((m*Nw+i)*3+1)*local_vector_prototype::max_length;
      local_vector[m][i].coef_s=local_vector_coefs+((m*Nw+i)*3+2)*local_vector_prototype::max_length;
      local_vector[m][i].length=0;
      local_vector[m][i].length_r=0;
      local_vector[m][i].length_s=0;
    }
  }
}
void region::inicialize7()
{  
  int aux=8*b+2;
  relative_point=new relative_point_prototype[Nr_max];
  int lc=aux*2+2*sets, ll=aux*2;
  relative_point_coefs=new content[Nr_max*lc];
  relative_point_labels=new int[Nr_max*ll];
  for (int i=0;i<Nr_max;i++) {
    relative_point[i].coef=relative_point_coefs+i*lc;
    relative_point[i].sl=relative_point_labels+i*ll;
    relative_point[i].coef2=relative_point_coefs+i*lc+aux;
    relative_point[i].sl2=relative_point_labels+i*ll+aux;
    relative_point[i].fe=relative_point_coefs+i*lc+2*aux;
    relative_point[i].fi=relative_point_coefs+i*lc+2*aux+sets;
  }
}
void region::inicialize8()
{
  M_sparse=new sparse_matrix_prototype(0,0);
  HG_sparse=new sparse_matrix_prototype(0,0);
  
  HG=M=LU=0;
  IPIV=0;
  M_length=0;
  M_locvecval=0;
  M_locvecind=0;
  M_length_max=18;
}

void region::inicialize()
{
  inicialize1();
  inicialize2();
  inicialize3();
  inicialize4();
  inicialize5();
  inicialize6();
  inicialize7();
  inicialize8();
}

void region::clean_multiples()
{
  if (multiples[0][0]==0) message(logg,"looks as multiples are not existing!!\n",where_to_print_critical);
  else
    for (int m=0;m<sets*sets;m++) {
      for (int i=0;i<Nw;i++) {
        for(int k=0;k<8*b+1;k++) {
          multiples[m][i][k].upsl.sl=-1;
          multiples[m][i][k].val=0;
         }
      }
    }
  for (int m=0;m<sets*sets;m++) for (int i=0;i<Nw;i++) {
    local_vector[m][i].length=local_vector[m][i].length_r=local_vector[m][i].length_s=0;
  }
}


//deallocate memory
void region::clean()
{
  for (int i=0;i<all_prec;i++) {
    for (int j=0;j<all_op;j++) {
      delete p_op_mat[j][i];
    }
  }
  if (sectorized) for (int l=0;l<Nl;l++) {
    delete sector[l].dim; 
    delete sector[l].offset;
  }
  else {
    delete sector[0].dim; 
    delete sector[0].offset;
  }
  delete sector;

  for (int i=0;i<Nx1;i++) {
    if (xline[i].np>0) {
      if (sectorized) {
        for (int aj=0;aj<xline[i].np;aj++) delete xline[i].active_line_point[aj].nearby_sector;
      }
      delete xline[i].active_line_point;
    }
  }

  delete grid_point;
  delete active_point;
  delete relative_point;
  delete relative_point_coefs;
  delete relative_point_labels;
  delete sl_btp;
  delete xline;
  delete yline;
  delete others;
  delete temp_mv;
  
  new_values();
  for (int i=0;i<all_op;i++) {
    delete dif_op_vals[i];
    delete dif_op_filled[i];
  }

  for (int m=0;m<sets*sets;m++) {
    for (int i=0;i<Nw;i++) delete multiples[m][i];
    delete multiples[m];
  }
  delete multiples;
  
  delete local_vector_coefs;
  delete local_vector_labels;
  for (int m=0;m<sets*sets;m++) delete local_vector[m];
  delete local_vector;
  
  //graphene
  deallocate_graphene();
}

bool region::inside_with_neighbors(int i, int j)
{
  bool in=func(xd(i), yd(j));
  if (!in) return(false);
  //!!!the following line brings it to the previous version - before the graphene problem (three neighbors out) appeared
  //if (in) return(true);
  int number_outside_neighbors=0;
  if (!func(xd(i-1), yd(j))) number_outside_neighbors++;
  if (!func(xd(i+1), yd(j))) number_outside_neighbors++;
  if (!func(xd(i), yd(j-1))) number_outside_neighbors++;
  if (!func(xd(i), yd(j+1))) number_outside_neighbors++;
  if (number_outside_neighbors>2) return(false); else return(true);
}


//label points with sequential labels, account for PBC by defining equivalent couples such, that b && Nb lines are IDENTICAL!!!
// i\in<0,b-1> : [i, Nxb-b+i], [Nxb+i,b+i] for x direction
void region::fill_arrays()
{
  int label=0;
  bool in;

  //grid points
  for (int i=0;i<Nx1;i++) {
    int sl_line=0;
    for (int j=0;j<Ny1;j++) {
      grid_point_prototype &gp=grid_point[i*Ny1+j];
      gp.sl=gp.sl_g[0]=gp.sl_g[1]=-1;
      gp.active_coordinate.i=-1;
      gp.active_coordinate.j=-1;
      gp.sector=-1;
      gp.type=0;
      gp.type_graphene=0;
      gp.plaq=-1;
      gp.b_sum_i=gp.b_sum_j=0;
      //!!!remove points which have 3 outside neighbours - necessary for graphene
      //in=func(xd(i),yd(j));
      in=inside_with_neighbors(i,j);
      if (((i<b) || (i>Nxb) || (j<b) || (j>Nyb))) in=false;
      //fprintf(logg,"fill_arrays: point grid (i=%i, j=%i): in:%i, inside(x=%.3e,y=%.3e)=%i\n",i,j,in,xd(i),yd(j),func(xd(i),yd(j)));
      if (in) {
        gp.sl=label++;
        gp.active_coordinate.i=i-first_active_xline;
        gp.active_coordinate.j=sl_line++;
        gp.sector=0;
        gp.type=1;
      }
    }
  }
  //identify and label active components for graphene wavefunction within the pseudospin subspace
  Nw_g=0;
  #ifdef GRPH_CENTRAL
  int central_sl;
  if (!coordinate2s(0,0,central_sl)) fprintf(logg,"central point out of the grid: graphene needs modification\n");
  #endif
  for (int i=0;i<Nx1;i++) {
    for (int j=0;j<Ny1;j++) {
      grid_point_prototype &gp=grid_point[i*Ny1+j];
      if (gp.type!=1) continue;
      #ifdef GRPH_CENTRAL
      if (gp.sl==central_sl) {gp.type_graphene=3;continue;}
      #else
      if (gp.sl==0 || gp.sl==Nw-1) {gp.type_graphene=4;continue;}
      #endif
      bool border=false;
      for (int di=-1;di<2;di++) for (int dj=-1;dj<2;dj++) if (grid_point[(i+di)*Ny1+j+dj].type!=1) {
        border=true;
        gp.b_sum_i+=di;
        gp.b_sum_j+=dj;        
      }
      if (border) gp.type_graphene=2; else gp.type_graphene=1;
      gp.sl_g[0]=Nw_g++;
      if (!border) gp.sl_g[1]=Nw_g++; else gp.sl_g[1]=gp.sl_g[0];
    }
  }

  //identify and label plaquettes for graphene hamiltonian for the pseudospin subspace
  Nw_plaq=0;
  for (int i=0;i<Nx1-1;i++) {
    for (int j=0;j<Ny1-1;j++) {
      grid_point_prototype &gp=grid_point[i*Ny1+j];
      if (gp.type!=1) continue;
      bool plaq=true;
      for (int di=0;di<2;di++) for (int dj=0;dj<2;dj++) if (grid_point[(i+di)*Ny1+j+dj].type!=1) plaq=false;
      if (plaq) gp.plaq=Nw_plaq++;
    }
  }
  if (Nw_g!=2*Nw_plaq) fprintf(logg,"number of plaquettes (%i) not equal to twice the number of wavefunction components (%i) for graphene\n",Nw_plaq,Nw_g);
  g_length=Nw_g*sets/2;
  fprintf(logg,"graphene: Nw_plaq=%i, Nw_g=%i, sets/2=%i, g_length=%i\n",Nw_plaq, Nw_g, sets/2, g_length);


/*  //periodic boundary conditions
  #define aux(a,b,c,d) \
    {grid_point[(a)*Ny1+(b)].sl=-grid_point[(c)*Ny1+(d)].sl;}
  if (pbc_y) {
    for (int i=b;i<Nxb1;i++) {
      for (int j=0;j<b;j++) {
        aux(i,j,i,Nyb-b+j+1)
	aux(i,Nyb1+j,i,b+j)
      }
    }
  }
  if (pbc_x) {
    for (int j=b;j<Nyb1;j++) {
      for( int i=0;i<b;i++) {
        aux(i,j,Nxb-b+i+1,j)
	aux(Nxb1+i,j,b+i,j)
      }
    }
  }
  if (pbc_x && pbc_y) {
    for (int i=0;i<b;i++) {
      for (int j=0;j<b;j++) {
        aux(i,j,Nxb-b+i,Nyb-b+j)
	aux(i,Ny-j,Nxb-b+i,b+b-j)
	aux(Nx-i,j,b+b-i,Nyb-b+j)
	aux(Nx-i,Ny-j,b+b-i,b+b-j)
      }
    }
  }
#undef aux
*/

  //xlines
  for (int i=0;i<Nx1;i++) {
    int count=0;
    xline[i].fpc=0;
    xline[i].fpsl=-1;
    bool first=true;
    for (int j=b;j<Nyb1;j++) {
      if (grid_point[i*Ny1+j].type!=1) continue;
      int label=c2sl(i,j);
      active_point[label].coordinate.i=i;
      active_point[label].coordinate.j=j;
      active_point[label].volume=1;
      if (first) {
        first=false;
        xline[i].fpc=j;
        xline[i].fpsl=label;
      }
      xline[i].active_line_point[count].j=j;
      count++;
    }
    if (xline[i].np!=count) {printf("number of points on xline %i different (now:%i previously:%i)!\n",i,count,xline[i].np); exit(1);}
  }
  if (neumann && relative_point_method=='r') active_point[0].volume=active_point[Nw-1].volume=0.5;//!!!correct only for 1D and NN
}


void region::sectorize(int (inside_sector)(prec, prec), int Nl_)
{
  sectorized=true;
  Nl=Nl_;

  //sector limits allocate
  delete sector[0].dim;
  delete sector[0].offset;
  delete sector;
  sector=new sector_prototype[Nl];
  for (int l=0;l<Nl;l++) {
    sector[l].imin=Nx1;
    sector[l].imax=-1;
  }

  //first tag all grid points by sector labels
  for (int i=0;i<Nx1;i++) {
    //identify sector label for each active point
    for (int j=b;j<Nyb1;j++) {
      int sec=inside_sector(xd(i),yd(j));
      //fprintf(logg,"sectorize: sec %i (previously %i) for grid point (i=%i, j=%i)->(x=%.3e, y=%.3e)\n",sec,grid_point[i*Ny1+j].sector,i,j,xd(i),yd(j));
      if (sec==-1 || grid_point[i*Ny1+j].type==0) {//the point is outside active region
        if (grid_point[i*Ny1+j].sector!=-1) {printf("discrepancy of previous sector %i (seq.label:%i) and actual sector %i in sectorize; grid point (i=%i, j=%i)->(x=%.3e, y=%.3e)\n", grid_point[i*Ny1+j].sector, grid_point[i*Ny1+j].type, sec, i, j, xd(i), yd(j));exit(1);}
        continue;
      }
      if (i<sector[sec].imin) sector[sec].imin=i;
      if (i>sector[sec].imax) sector[sec].imax=i;
      grid_point[i*Ny1+j].sector=sec;
      //fprintf(logg,"sec %i point added grid(i=%i, j=%i)\n",sec,i,j);
    }
  }
  for (int l=0;l<Nl;l++) {
    if (sector[l].imin==Nx1 || sector[l].imax==-1) {printf("region::sectorize: no points for sector %i found\n",l);exit(1);}
  }

  for (int l=0;l<Nl;l++) {
    int length=sector[l].imax-sector[l].imin+1;
    sector[l].dim=new int[length];
    sector[l].offset=new int[length];
    for (int si=0;si<length;si++) sector[l].dim[si]=sector[l].offset[si]=0;
  }

  //xlines offsets and active points neighbours
  bool* seen=new bool[Nl];
  for (int i=first_active_xline;i<=last_active_xline;i++) {
    for (int l=0;l<Nl;l++) seen[l]=false;
    int aj=0;
    for (int j=b;j<Nyb1;j++) {
      if (grid_point[i*Ny1+j].type!=1) continue;//if not, it is an active point
      if (aj>=xline[i].np) {printf("beyond active points: aj=%i, xline[i].np=%i; grid_point[%i,%i].type=%i\n",aj, xline[i].np,i,j,grid_point[i*Ny1+j].type);exit(1);}

      //check for sector offsets
      int sec=grid_point[i*Ny1+j].sector;
      sector[sec].dim[i-sector[sec].imin]++;
//printf("adding point of sector %i on sec. line %i (grid line %i)\n",sec,i-sector[sec].imin,i);
      if (!seen[sec]) {//first point iniside the sector l
        seen[sec]=true;
        sector[sec].offset[i-sector[sec].imin]=aj;
      }
      //check the neighbourhood for other sectors
      xline[i].active_line_point[aj].nearby_sector=new bool[Nl];
      xline[i].active_line_point[aj].nearby_any_sector=false;
      for (int l=0;l<Nl;l++) xline[i].active_line_point[aj].nearby_sector[l]=false;
      int length=(int) act_prec+1;
      for (int ri=-length;ri<length+1;ri++) {
        for (int rj=-length;rj<length+1;rj++) {
          int nearby_sector=grid_point[(i+ri)*Ny1+(j+rj)].sector;
          if (nearby_sector!=sec && nearby_sector!=-1) xline[i].active_line_point[aj].nearby_sector[nearby_sector]=true;
        }
      }
      for (int l=0;l<Nl;l++) if (xline[i].active_line_point[aj].nearby_sector[l]) xline[i].active_line_point[aj].nearby_any_sector=true;
      
      //if (xline[i].active_line_point[aj].nearby_any_sector) {
        //printf("point nearby found: i=%i, aj=%i; nearby_sector:",i,aj);
        //for (int l=0;l<Nl;l++) printf("%i ",xline[i].active_line_point[aj].nearby_sector[l]);
      //}
      aj++;
    }
  }
  delete seen;

}

int region::local_vector_prototype::max_length;
void region::local_vector_prototype::update(int sl_new, content coef_new, int type) 
{
  int *l;       //pointer to the length
  int *label;   //array of corresponding labels
  content* c;   //array of coefficients
  switch (type) {
    case 1 : {  //inside points
      l=&length;
      label=sl;
      c=coef;
      break;
    }
    case 2 : {  //relative points
      l=&length_r;
      label=sl_r;
      c=coef_r;
      break;
    }
    case 3 : {  //special inside points
      l=&length_s;
      label=sl_s;
      c=coef_s;
      break;
    }
    default : {
      printf("wrong type %i in local_vector_prototype::update\n",type);
      exit(1);
    }
  }
    
  int i;
  for (i=0;i<*l;i++) if (label[i]==sl_new) break;
  if (i<*l) c[i]+=coef_new;
  else {
    (*l)++;
    if ((*l)==max_length) {
      printf("local_vector has grown beyond maximal length (%i)\n",max_length);
      exit(1);
    }
    c[i]=coef_new;
    label[i]=sl_new;
  }
}

//input, output are indexed by label in <0,Nw-1>
void region::dn_point_operation(content * input, content *output, int nl, int nh,int ml, int mh, prec* koeficients)
{
  int xcoor[(2*all_prec+1)*(2*all_prec+1)];	//space for relative coordinates
  int ycoor[(2*all_prec+1)*(2*all_prec+1)];
  prec koef[(2*all_prec+1)*(2*all_prec+1)];	//and koeficients

  int label,next,index,count;
  int i,j,k,l;
  prec x;

  count=0;
  for (k=nl;k<=nh;k++) {
    for (l=ml;l<=mh;l++) {
      index=(k-nl)+(l-ml)*(nh-nl+1);
      x=koeficients[index];
      if (x==0) continue;
      xcoor[count]=k;
      ycoor[count]=l;
      koef[count]=x;
      count++;
    }
  }

  for (i=b;i<Nxb1;i++) {
    for (j=b;j<Nyb1;j++) {
      if (grid_point[i*Ny1+j].type!=1) continue;
      label=c2sl(i,j);
      output[label]=0;
      for(k=0;k<count;k++) {
        if (grid_point[(i+xcoor[k])*Ny1+(j+ycoor[k])].type!=1) continue;
        next=ac2sl(i+xcoor[k],j+ycoor[k]);
        output[label]+=input[next]*koef[k];
#ifdef DEBUG3
        fprintf(logg,"dn_point_operation: point with label %3i\n  neighbour label %3i  value %.18e+i%.18e\n  multiplied by koef %+.18e\n actualized output %.18e+i%.18e\n ",label, next, input[next].real(), input[next].imag(),koef[k],output[label].real(),output[label].imag());
#endif
      }
    }
  }
}

void region::dn_point_operation_tot(int w_f_x, int w_f_y, int nl, int nh, int ml, int mh, prec* koef, content f(prec,prec), content kf, content exp_integral(prec, prec, prec, prec))
{
  int index_aux[(2*all_prec+1)*(2*all_prec+1)];
  int index2_aux[(2*all_prec+1)*(2*all_prec+1)];
  int label_aux=0;
  for (int k=-b;k<=b;k++) {
    for (int l=-b;l<=b;l++) {
      index_aux[label_aux]=relative2multiple(k,l);
      index2_aux[label_aux]=relative2matrix(nl,nh,ml,mh,k,l);
      label_aux++;
    }
  }

  int label,neighbour,index,index2;
  content value;

  //multiples and local vector are doing the same, with the latter compressing the matrix into a set of points used
  //multiples are kept only because of the routine Hsecsec for the recursive GF, which is faster through the multiples than through the local vector
  for (int i=b;i<Nxb1;i++) {
    for (int j=b;j<Nyb1;j++) {
      if (grid_point[i*Ny1+j].type!=1) continue;
      label=c2sl(i,j);
      multiple* &pmult=multiples[w_f_x*sets+w_f_y][label];
      label_aux=-1;
      for (int k=-b;k<=b;k++) {
        for (int l=-b;l<=b;l++) {
          label_aux++;
          if ((index=index_aux[label_aux])==-1) continue;
          if ((index2=index2_aux[label_aux])==-1) continue;
          multiple &mult=pmult[index];
          if (grid_point[(i+k)*Ny1+j+l].type!=1) continue;
          neighbour=ac2sl(i+k,j+l);
          value=koef[index2]*kf;
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SYMMETRIZED HAMILTONIAN???
          //if (f!=0) value*=(f(xd(i),yd(j))+f(xd(i+k),yd(j+l)))/2.0;
          if (f!=0) value*=f(xd(i),yd(j));
          if (operators_sym1) value*=exp_integral(xd(i),yd(j),xd(i+k),yd(j+l));
          if ((value.real()==0 && value.imag()==0)) continue;
          if (mult.upsl.sl>0 && mult.upsl.sl!=neighbour) {
          sprintf(buffer, "WRONG: point at (%3i, %3i) with label %3i considering relative at (%3i, %3i)\nmultiple index %3i, matrix index %3i, old neighbour label %3i new label %3i\n",i,j,label,k,l,index, index2,mult.upsl.sl,neighbour);
          splachni(logg, buffer,where_to_print_critical);}
          mult.upsl.sl=neighbour;
          mult.val+=value;
        }
      }
    }
  }
  
  for (int a=0;a<Nw;a++) {//for all active points
    int i=active_point[a].coordinate.i;
    int j=active_point[a].coordinate.j;
    int label_aux=-1, matrix_label;
    for (int k=-b;k<=b;k++) {//the local matrix relative coordinates
      for (int l=-b;l<=b;l++) {
        label_aux++;
        if ((matrix_label=index2_aux[label_aux])==-1) continue;
        int type=grid_point[(i+k)*Ny1+(j+l)].type;
        if (type==0 && dirichlet) continue;
        content value=koef[matrix_label]*kf;
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SYMMETRIZED HAMILTONIAN???
        //if (f!=0) value*=(f(xd(i),yd(j))+f(xd(i+k),yd(j+l)))/2.0;
        if (f!=0) value*=f(xd(i),yd(j));
        if (operators_sym1) value*=exp_integral(xd(i),yd(j),xd(i+k),yd(j+l));
        if ((value.real()==0 && value.imag()==0)) continue;
        //if (a==0 && neighbour==1) fprintf(logg,"loc. mat before update(wfx=%i wfy=%i):length=%i value[0]=%e\n",w_f_x,w_f_y,local_vector[w_f_x*sets+w_f_y][a].length, local_vector[w_f_x*sets+w_f_y][a].coef[0].real()); 
        if (type==1) local_vector[w_f_x*sets+w_f_y][a].update(grid_point[(i+k)*Ny1+(j+l)].sl,value,1);
        else if (neumann) {
          if (type==0) {//new active relative point
            relative_point[Nr].sl0=a;//temporary - will be changed later to "nearest" inside point
            relative_point[Nr].coordinate.i=i+k;
            relative_point[Nr].coordinate.j=j+l;
            relative_point[Nr].solved=false;
            relative_point[Nr].length=0;
            grid_point[(i+k)*Ny1+(j+l)].sl=Nr;
            grid_point[(i+k)*Ny1+(j+l)].type=2;
            Nr++;
            if (Nr_max-1==Nr) {printf("maximal active relative points reached\n");exit(1);}
            #ifdef DEBUG2
            fprintf(logg,"new relative point - number of relative points updated to %i\n",Nr);
            #endif
          }
          local_vector[w_f_x*sets+w_f_y][a].update(grid_point[(i+k)*Ny1+(j+l)].sl,value,2);
          #ifdef DEBUG2
          fprintf(logg,"relative point: (at coors %i-%i) reached from active point %i (at coors %i-%i)\n", i+k, j+l, a, i, j);
          #endif
        }
        if (periodic) {printf("periodic bc not implemented in dn_point_operation_tot\n"); exit(1);}
        //if (a==0 && neighbour==1) fprintf(logg,"loc. mat after update(wfx=%i wfy=%i):length=%i value[0]=%e\n",w_f_x,w_f_y,local_vector[w_f_x*sets+w_f_y][a].length, local_vector[w_f_x*sets+w_f_y][a].coef[0].real()); 
      }
    }
  }

  
#ifdef REGION_CHECK
  for (int i=b;i<Nxb1;i++) {
    for (int j=b;j<Nyb1;j++) {
    if (grid_point[i*Ny1+j].type!=1) continue;
    label=c2sl(i,j);
    //fprintf(logg,"multiples: wfx=%i wfy=%i nl=%i nh=%i ml=%i mh=%i for point with label %3i at (%3i, %3i):\n",w_f_x, w_f_y,nl, nh, ml, mh,label,i,j);
    local_vector_prototype &lv=local_vector[w_f_x*sets+w_f_y][label];
    int neighbours=0,l;
    for (int k=0;k<8*b+1;k++) {
      multiple &mult=multiples[w_f_x*sets+w_f_y][label][k];
      neighbour=mult.upsl.sl;
      value=mult.val;
      //fprintf(logg,"neighbour %3i value %.6e+i%.6e:\n",neighbour,value.real(),value.imag());
      if (abs(value)==0) continue;
      neighbours++;
      for (l=0;l<lv.length;l++) if (neighbour==lv.sl[l]) break;
      if (l==lv.length) fprintf(logg, "active point %i : neighbour %i not found in local_vector\n",label, neighbour);
      else if (abs(value-lv.coef[l])>1e-15) fprintf(logg, "active point %i : neighbour %i multiple value=%.6e+i%.6e versus local_vector value %.6e+i%.6e (diff=%e)\n",label, neighbour, value.real(), value.imag(), lv.coef[l].real(), lv.coef[l].imag(),abs(value-lv.coef[l]));
    }
    //for (int k=0;k<lm.length;k++) {
      //sprintf(buffer,"loc. matrix: neighbour %3i value %.6e+i%.6e:\n",lm.sl[k],lm.coef[k].real(),lm.coef[k].imag());
      //splachni(logg,buffer,2);
    //}
    if (neighbours!=lv.length) fprintf(logg,"active point %i : neighbours in multiple %i versus local_vector length %i!!!\n",label,neighbours,lv.length);
    }
  }
#endif
}

//updates the boundary condition using k [internal units - wavevector]
void region::operate_tot_init_boundary(content k)
{
  bc_coefficient=content(0,1)*k;
}

//updates the overall shift [internal units - energy]
void region::operate_tot_init_shift(double shift)
{
  overall_shift=shift;
}


//convert the local vector into boundary point into local vector into inside points plus special point
//update the special points array
void region::operate_tot_init_boundary()
{
  if (dirichlet) return;
  if (periodic) {printf("periodic bc not implemented in operate_tot_init_boundary\n"); exit(1);}
  for (int from=0;from<sets;from++) {
    for (int to=0;to<sets;to++) {
      for (int a=0;a<Nw;a++) {
        local_vector_prototype &lv=local_vector[from*sets+to][a];
        if (lv.length_r==0) continue;
        #ifdef DEBUG2
        fprintf(logg,"local vector (from=%i, to=%i, a=%i) with %i relative points encountered\n",from, to, a,lv.length_r);
        #endif
        //the point reaches out to the boundary
        for (int i=0;i<lv.length_r;i++) {
          if (relative_point_method=='r') {
            content c=lv.coef_r[i];
            relative_point_prototype &rp=relative_point[lv.sl_r[i]];
            #ifdef DEBUG2
            fprintf(logg,"relative point %i (solved:%i) is cared for (relates to %i active points and special point #%i)\n",lv.sl_r[i],rp.solved,rp.length,rp.sl0);
            #endif
            if (!rp.solved) express_relative_point(lv.sl_r[i]);
            for (int j=0;j<rp.length;j++) lv.update(rp.sl[j],c*rp.coef[j],1);
            lv.update(rp.sl0, rp.c1*c,3);
          }
          else interpolate_relative_point(lv.sl_r[i]);
        }

        //insert the point inside the list of boundary touching points
        int i;
        for (i=0;i<N_btp;i++) if (a==sl_btp[i]) break;
        if (i==N_btp) {
          N_btp++;
          sl_btp[i]=a;
          #ifdef DEBUG2
          fprintf(logg,"the active point %i not found among btp's - added - actual number of btp's %i\n",a, N_btp);
          #endif
        }
      }
    }
  }
}

int region::relative_point_prototype::max_length;
void region::relative_point_prototype::update(int sl_new, content coef_new) 
{
  int i;
  for (i=0;i<length;i++) if (sl[i]==sl_new) break;
  if (i<length) coef[i]+=coef_new;
  else {
    length++;
    if (length==max_length) {
      printf("local_vector has grown beyond maximal length (%i)\n",max_length);
      exit(1);
    }
    coef[i]=coef_new;
    sl[i]=sl_new;
  }
}

//fills in the data, calls itself recursively
void region::express_relative_point(int sl)
{
  relative_point_prototype &rp=relative_point[sl];//myself
#ifdef DEBUG2
  fprintf(logg,"expressing relative point %i through inside points (is already done=%i)\n",sl,rp.solved);
  fprintf(logg,"relative point: length %i, actual special point %i (coef %.3e%+.3e), list of active points:\n",rp.length,rp.sl0,rp.c1.real(), rp.c1.imag());
  for (int i=0;i<rp.length;i++) fprintf(logg,"\tactive point %i (coef %.3e%+.3e)\n",rp.sl[i],rp.coef[i].real(), rp.coef[i].imag());
#endif
  if (rp.solved) return;
  
  coordinate_prototype pos;  //position yourself on the actual special point and look towards me
  pos.i=active_point[rp.sl0].coordinate.i;
  pos.j=active_point[rp.sl0].coordinate.j;
  int di=rp.coordinate.i-pos.i;
  int dj=rp.coordinate.j-pos.j;
  int si=0, sj=0;
  if (di<0) si=-1; else if (di>0) si=1;
  if (dj<0) sj=-1; else if (dj>0) sj=1;

  //walk towards me
  do {
    if (grid_point[pos.i*Ny1+pos.j].type==1) //if inside, update my special point
    rp.sl0=grid_point[pos.i*Ny1+pos.j].sl;
    else if (grid_point[pos.i*Ny1+pos.j].type==2) //if outside, make sure it is solved for
    express_relative_point(grid_point[pos.i*Ny1+pos.j].sl);
    else printf("point of wrong type %i in express_relative_point called for sl=%i\n",grid_point[pos.i*Ny1+pos.j].type, sl);
    pos.i+=si;
    pos.j+=sj;
  } while (pos.i!=rp.coordinate.i || pos.j!=rp.coordinate.j);
  #ifdef DEBUG2
  fprintf(logg,"all points towards the grid should be ok; the special point was updated to %i\n",rp.sl0);
  #endif
  //all points are ok
  int order=act_prec+1;
  prec *coors=new prec[2*order+1];
  prec *d=new prec[2*order+1];
  prec h;    //grid step
  int k;
  if (si==0 && sj*sj==1) {h=hy;k=rp.coordinate.j-active_point[rp.sl0].coordinate.j;}
  else if (sj==0 && si*si==1) {h=hx;k=rp.coordinate.i-active_point[rp.sl0].coordinate.i;}
  else {printf("other than dx/dy not implemented in express_relative_point\n");exit(1);}
  if (k<0) k*=-1;//the outward oriented derivative
  for (int i=0;i<2*order+1;i++) coors[i]=(k-i)*h;
  construct_differential_operator(1,2*order+1,coors,d);
  //walk away from me, updating the coefficients
  int i=0;//pointer into d's
  rp.c1=1/d[0];//the special point coefficient to start with (if no boundary points except me, that's it)
  //fprintf(logg,"updating special point to %e\n",1/d[0]);
  do {
    pos.i-=si;
    pos.j-=sj;
    i++;
    grid_point_prototype &gp=grid_point[pos.i*Ny1+pos.j];
    if (gp.type==1) rp.update(gp.sl,-d[i]/d[0]);//the point is inside
    else if (gp.type==2) {//the point is outside
      relative_point_prototype &rp_aux=relative_point[gp.sl];
      if (!rp_aux.solved) {printf("relative point (%i) should be solved for now, but is not in express_relative_point\n",gp.sl);exit(1);}//make sure the point is solved for
      if (rp_aux.sl0!=rp.sl0) {printf("the special points of relative points %i and %i are not the same (%i vs %i) in express_relative_point\n", sl, gp.sl, rp_aux.sl0, rp.sl0);exit(1);}//make sure special points are the same 
      //all inside points that the relative point is expressed through 
      for (int j=0;j<rp_aux.length;j++) rp.update(rp_aux.sl[j],-d[i]/d[0]*rp_aux.coef[j]);
      //finally update the special point coefficient
      rp.c1-=d[i]*rp_aux.c1/d[0];
      //fprintf(logg,"updating special point by %e * %e to %e\n",rp_aux.c1.real(), -d[i]/d[0], rp.c1.real()); 
    }
  } while (i<2*order);
  rp.solved=true;
  #ifdef DEBUG2
  fprintf(logg,"after the update: length %i, actual special point %i (coef %.3e%+.3e), list of active points:\n",rp.length,rp.sl0,rp.c1.real(), rp.c1.imag());
  for (int i=0;i<rp.length;i++) fprintf(logg,"\tactive point %i (coef %.3e%+.3e)\n",rp.sl[i],rp.coef[i].real(), rp.coef[i].imag());
  #endif
  delete coors;
  delete d;
}

//fill in length and arrays sl2, coef2 in relative point sl
void region::interpolate_relative_point(int sl)
{
  relative_point_prototype &rp=relative_point[sl];//myself
  rp.re=sqrt(xd(rp.coordinate.i)*xd(rp.coordinate.i)+yd(rp.coordinate.j)*yd(rp.coordinate.j));
  coordinate_prototype pos=rp.coordinate, pos_min=rp.coordinate;
  //find closest point inside
  do {
    double r_min=0;
    for (int di=-1;di<2;di++) {
      int i=pos.i+di;
      if (i<0 || i>Nx) continue;
      for (int dj=-1;dj<2;dj++) {
        int j=pos.j+dj;
        if (j<0 || j>Ny) continue;
        double r=sqrt(xd(i)*xd(i)+yd(j)*yd(j));
        if (r_min==0 || r_min>r) {pos_min.i=i;pos_min.j=j;r_min=r;}
      }
    }
    pos=pos_min;
    if (grid_point[pos.i*Ny1+pos.j].type==1) break;        
  } while(true);
  rp.ri=sqrt(xd(pos.i)*xd(pos.i)+yd(pos.j)*yd(pos.j));

  //1D only: interpolation is simply a copy of the nearest inside point
  rp.length2=1;
  rp.sl2[0]=grid_point[pos.i*Ny1+pos.j].sl;
  rp.coef2[0]=1;
  rp.solved=true;
#ifdef DEBUG2
  fprintf(logg,"relative point %i interpolated:\nclosest active point at (%i,%i) of type %i with seq label %i\n",sl,pos.i,pos.j,grid_point[pos.i*Ny1+pos.j].type,rp.sl2[0]);
#endif
  if (rp.sl2[0]<0 || rp.sl2[0]>=Nw) exit(0);
}

void region::extrapolate_relative_points(content *input)
{
  //general:
  /*for (int r=0;r<Nr;r++) {
    relative_point_prototype& rp=relative_point[r];
    int length2=rp.length2;
    for (int set=0;set<max_sets;set++) {
      rp.fi[set]=0;
      for (int i=0;i<length2;i++) rp.fi[set]+=input[set*Nw+rp.sl2[i]]*rp.coef[i];
      rp.fe[set]=rp.fi[set]*exp(content(0,k)*(rp.re-rp.ri));
    }
  }*/
  //simplified:1D and 1 set only
  for (int r=0;r<Nr;r++) {
    relative_point_prototype& rp=relative_point[r];
    if (act_prec==0) {
      prec x=xd(rp.coordinate.i);
      rp.fi[0]=input[rp.sl2[0]]*rp.coef2[0];
      //right end is open
      if (x>0) {
        rp.fe[0]=rp.fi[0]*exp(bc_coefficient*(rp.re-rp.ri)); //finite derivative
        //rp.fe[0]=rp.fi[0];                                 //zero derivative
        //rp.fe[0]=0;                                          //zero wavefunction
      }
      //at the left end the cavity is closed or open
      if (x<0) {
        rp.fe[0]=rp.fi[0]*exp(bc_coefficient*(rp.re-rp.ri)); 
        //rp.fe[0]=rp.fi[0];
        //rp.fe[0]=0;
      }
      continue;
    }
    if (act_prec!=1) {fprintf(logg,"precision %i not implemented in extrapolate\n",act_prec); exit(1);}
    rp.fe[0]=input[rp.sl2[0]]*exp(bc_coefficient*(rp.re-rp.ri));
    //int dist=rp.coordinate.i-active_point[rp.sl2[0]].coordinate.i;
    //if (abs(dist)==1) rp.fe[0]=input[rp.sl2[0]]*exp(bc_coefficient*3.0*(rp.re-rp.ri));
    //else {
      //rp.fe[0]=27.0*(exp(bc_coefficient*(rp.re-rp.ri))-exp(-bc_coefficient*(rp.re-rp.ri)))*input[rp.sl2[0]];
      //rp.fe[0]+=exp(1.5*bc_coefficient*(rp.re-rp.ri))*input[rp.sl2[0]-dist/2];
    //}
    //fprintf(logg, "extrapolating relative point %i: input point %i (%.2e%+.2e) gives fi=(%.2e%+.2e) and fe=(%.2e%+.2e)\n", r, rp.sl2[0], input[rp.sl2[0]].real(), input[rp.sl2[0]].imag(), rp.fi[0].real(), rp.fi[0].imag(), rp.fe[0].real(), rp.fe[0].imag());
  }  
}


void rot2D(prec& x, prec& y, prec c, prec s)
{
  prec aux_x=c*x-s*y;
  prec aux_y=c*y+s*x;
  x=aux_x;
  y=aux_y;
}

void region::operate(content *input, content* output,int w_f_x, int w_f_y, operators op, content koef, content f(prec,prec), setting_or_adding soa)
{
  content * akt_in;
  bool only_koef=false;
  if (f==0) only_koef=true;
  output+=w_f_y*Nw;		//to be indexed by label in <1,Nw>
  input+=w_f_x*Nw;

  if (op<all_op) {//is it diferential operator?
    int x=op_mat_dim[op][act_prec].x;
    int y=op_mat_dim[op][act_prec].y;

    if (remember[op]) {//is this one under remember verwendung?
      akt_in=dif_op_vals[op][w_f_x];
      if (! dif_op_filled[op][w_f_x]) {//is it already computed?
        //it is ! - compute it!
        dif_op_vals[op][w_f_x]=akt_in=new content[Nw];
        dif_op_filled[op][w_f_x]=true;
        dn_point_operation(input,akt_in,-x,x,-y,y,p_op_mat[op][act_prec]);
      }
    }
    else { //no remember bothers
      akt_in=others;
      dn_point_operation(input,akt_in,-x,x,-y,y,p_op_mat[op][act_prec]);
    }
    //everything is prepared for final linear combination
  }

  else {//it is ! a diferential operator
    if (op==lin_com) {
        if (input-w_f_x*Nw==0) {
          input=others;
	        for (int i=0;i<Nw;i++) others[i]=1;
        }
        akt_in=input;
    }
    else if (op>=invx && op<=inv) {
      int xi=1; if (op==invx || op==inv) xi=-1;
      int yi=1; if (op==invy || op==inv) yi=-1;
      prec c=cos(inv_angle);
      prec s=sin(inv_angle);
      content* f[6];
      value_at_init(f, input);
      for (int i=0;i<Nw;i++) {
        int j,k;
        sl2c(i,j,k);
        prec x=xd(j), y=yd(k);
        rot2D(x,y,c,-s);
        x*=xi;y*=yi;
        rot2D(x,y,c,s);
        others[i]=value_at(x, y, f);
      }
      value_at_release(f);
      akt_in=others;
    }
    else {
      printf("non - valid value in operate\n");
      akt_in=0; exit(1);
    }
  }


  if (only_koef) {
    if (soa==add) for (int i=0;i<Nw;i++) {
      output[i]+=koef*akt_in[i];
#ifdef DEBUG3
      sprintf(buffer,"operate/add: point with (label) %3i\n input value[label] %.18e+i%.18e\n multiplied by(koef) %+.18e+i%+.18e\noutput %.18e+i%.18e\n ",i+1 ,akt_in[i].real(),akt_in[i].imag(),koef.real(),koef.imag(),output[i].real(),output[i].imag());
      splachni(logg,buffer,2);
#endif
    }
    else
    for (int i=0;i<Nw;i++) {
      output[i]=koef*akt_in[i];
#ifdef DEBUG3
      sprintf(buffer,"operate/set: point with label %3i\n input value[i] %.18e+i%.18e\n multiplied by(koef) %+.18e+i%+.18e\noutput %.18e+i%.18e\n ",i+1 ,akt_in[i].real(),akt_in[i].imag(),koef.real(),koef.imag(),output[i].real(),output[i].imag());
      splachni(logg,buffer,2);
#endif
    }
  }
  else{
    int label,i,j;
    if (soa==add) {
      for (i=b;i<Nxb1;i++) {
        for (j=b;j<Nyb1;j++) {
          if (grid_point[i*Ny1+j].type!=1) continue;
          label=c2sl(i,j);
          output[label]+=koef*f(xd(i),yd(j))*akt_in[label];
#ifdef DEBUG3
          sprintf(buffer,"operate/add: point with label %3i\n input value[label] %.18e+i%.18e\n multiplied by(koef*function) %+.18e+i%+.18e\noutput %.18e+i%.18e\n ",label ,akt_in[label].real(),akt_in[label].imag(),(koef*f(xd(i),yd(j))).real(),(koef*f(xd(i),yd(j))).imag(),output[label].real(),output[label].imag());
          splachni(logg,buffer,2);
#endif
        }
      }
    }
    else {
      for (i=b;i<Nxb1;i++) {
        for (j=b;j<Nyb1;j++) {
          if (grid_point[i*Ny1+j].type!=1) continue;
          label=c2sl(i,j);
          output[label]=koef*f(xd(i),yd(j))*akt_in[label];
#ifdef DEBUG3
          sprintf(buffer,"operate/set: point with label %3i\n input value[label] %.18e+i%.18e\n multiplied by(koef*function) %+.18e+i%+.18e\noutput %.18e+i%.18e\n ",label ,akt_in[label].real(),akt_in[label].imag(),(koef*f(xd(i),yd(j))).real(),(koef*f(xd(i),yd(j))).imag(),output[label].real(),output[label].imag());
          splachni(logg,buffer,2);
#endif
        }
      }
    }
  }
}

double graphene_center_aux[5]={-1.0/12, 1.0/3, 0, 1.0/3, -1.0/12};

content region::graphene_boundary_aux(int i, int j, int sig)
{
  if (sig==0) return(content(1.0));
  double x=xd(i), y=yd(j);
  //grid_point_prototype& gp=grid_point[i*Ny1+j];
  //x=gp.b_sum_i;
  //y=gp.b_sum_j;
  bool inv=false;
  if (x<0) {inv=true; x=-x; y=-y;}    
  double phi;
  if (y==0) phi=0; else phi=atan(x/y);
  if (inv) phi+=M_PI;
  return(content(0,1)*exp(content(0,1)*phi));
}

//update graphene Hamiltonian by applying |sigi><sigj| op (unity/k+/k-) to i-j-th subsector of the Hamiltonian, which  is a matrix with linear dimension of g_length = N * Nw_g and N=sets/2 (N will be 2 if spin and sigma dof are considered)
//with possible functional, constant and peierls spatial dependence
void region::update_g_matrix(int sigi, int sigj, const char* op, int ioffset, int joffset, content kf, content f(prec,prec), content integral(prec, prec, prec, prec), content *H, sparse_matrix_prototype* H_sparse)
{
  //sanity and input
  if (sigi<0 || sigj<0 || sigi >1 || sigj >1) {printf("sigi=%i, sigj=%i wrong in region::update_g_matrix\n",sigi,sigj);exit(1);}
  if (ioffset<0 || joffset<0 || ioffset>=sets/2 || joffset>=sets/2) {printf("ioffset=%i, joffset=%i wrong in region::update_g_matrix\n",ioffset,joffset);exit(1);}
  const char* choices[]={"unity","k+","k-"};
  int choice;
  for (choice=0;choice<3;choice++) if (strcmp(choices[choice],op)==0) break;
  if (choice==3) {printf("operator=%s wrong in region::update_g_matrix\n",op);exit(1);}
  bool peierls_symmetrization=operators_sym1 && (choice!=0);

  //printf("peierls is on: %i with phases: x-(%e%+ei) y-(%e%+ei)\n",peierls_symmetrization, phasex.real(), phasex.imag(), phasey.real(), phasey.imag());
  
  //input/output offset
  //content* Hoffset=H+g_length*ioffset*Nw_g+joffset*Nw_g;
  
  //coefficients of the plaquette operator
  content c[2][2], ii=content(0,1);
  
  //one equation for every plaquette
  for (int i=0;i<Nx1;i++) {
    for (int j=0;j<Ny1;j++) {
      grid_point_prototype& gp_output=grid_point[i*Ny1+j];
      if (gp_output.plaq==-1) continue;
      //int output_index=gp_output.plaq*2+sigi;
      int output_index=(gp_output.plaq*2+sigi)*sets/2+ioffset;
      
      //if (peierls_symmetrization) vf_final*=integral(xd(i)+hx/2,yd(j)+hy/2,xd(i+di),yd(j+dj));            
      content phasex0=1.0, phasey0=1.0, phasex1=1.0, phasey1=1.0;

      if (peierls_symmetrization) {
        phasex0=integral(xd(i), yd(j), xd(i+1), yd(j));
        phasex1=integral(xd(i), yd(j+1), xd(i+1), yd(j+1));
        phasey0=integral(xd(i), yd(j), xd(i), yd(j+1));
        phasey1=integral(xd(i+1), yd(j), xd(i+1), yd(j+1));
      }
      switch (choice) {
        case 0 : {//unity
          c[0][0]=c[1][0]=c[0][1]=c[1][1]=0.25;
          break;
        }
        case 1 : {//k+ = (-i) dx + i(-i)dy = -i dx + dy
          c[0][0]= (-ii*(-1.0/hx)         - 1.0/hy)/2.0;
          c[1][0]= (-ii*(+1.0/hx)*phasex0 - 1.0/hy)/2.0;
          c[0][1]= (-ii*(-1.0/hx)         + 1.0/hy*phasey0)/2.0;
          c[1][1]= (-ii*(+1.0/hx)*phasex1 + 1.0/hy*phasey1)/2.0;
          break;
        }
        case 2 : {//k- = (-i) dx - i(-i)dy = -i dx - dy
          c[0][0]= (-ii*(-1.0/hx)         + 1.0/hy)/2.0;
          c[1][0]= (-ii*(+1.0/hx)*phasex0 + 1.0/hy)/2.0;
          c[0][1]= (-ii*(-1.0/hx)         - 1.0/hy*phasey0)/2.0;
          c[1][1]= (-ii*(+1.0/hx)*phasex1 - 1.0/hy*phasey1)/2.0;
          break;
        }
        default : c[0][0]=c[1][0]=c[0][1]=c[1][1]=0;
      }
      
      
      for (int di=0;di<2;di++) {
        for (int dj=0;dj<2;dj++) {
          grid_point_prototype& gp_input=grid_point[(i+di)*Ny1+(j+dj)];
          int input_index=gp_input.sl_g[sigj]*sets/2+joffset;
          content vf_final=kf;
          if (f!=0) vf_final*=f(xd(i)+hx/2,yd(j)+hy/2);
          switch (gp_input.type_graphene) {
            case 1 : {
              content aux=c[di][dj]*vf_final;
              if (H!=0) H[output_index*g_length+input_index]+=aux;
              if (H_sparse!=0) H_sparse->add(output_index,input_index, aux);
              break;
            }
            case 2 : {//border point - apply boundary condition
              vf_final*=graphene_boundary_aux(i+di,j+dj,sigj);
              content aux=c[di][dj]*vf_final;
              if (H!=0) H[output_index*g_length+input_index]+=aux;
              if (H_sparse!=0) H_sparse->add(output_index,input_index, aux);
              //PREVIOUS: Hoffset[output_index*g_length+input_index]+=c[di][dj]*vf_final;
              break;
            }
            case 3 : {//central point - is interpolated by its neighbors - must not be also boundary!!
              for (int p=-2;p<3;p++) {
                if (p==0) continue;
                input_index=grid_point[(i+di+p)*Ny1+(j+dj)].sl_g[sigj]*sets/2+joffset;
                //PREVIOUS Hoffset[output_index*g_length+input_index]+=c[di][dj]*graphene_center_aux[p+2]*vf_final;
                content aux=c[di][dj]*graphene_center_aux[p+2]*vf_final;
                if (H!=0) H[output_index*g_length+input_index]+=aux;
                if (H_sparse!=0) H_sparse->add(output_index,input_index, aux);
                input_index=grid_point[(i+di)*Ny1+(j+dj+p)].sl_g[sigj]*sets/2+joffset;
                 //PREVIOUS Hoffset[output_index*g_length+input_index]+=c[di][dj]*graphene_center_aux[p+2]*vf_final;
                if (H!=0) H[output_index*g_length+input_index]+=aux;
                if (H_sparse!=0) H_sparse->add(output_index,input_index, aux);
              }
              break;
            }
	    case 4 : {//first / last point. It is assumed zero.
	      break;
	    }
            default : {
             printf("wrong gp.type_g=%i at (%i,%i) in region::update_g_matrix\n", gp_input.type_graphene, i+di,j+dj);exit(1);
            }
          }
        }
      }
    }
  }
}

void region::convertM2locvec()
{
  if (M_length!=0) delete M_length;
  M_length=new int[g_length];     //number of non-zero entries on a given line of the M matrix
  if (M_locvecval!=0) delete M_locvecval;
  M_locvecval=new content[g_length*M_length_max];   //values of non-zero entries
  if (M_locvecind!=0) delete M_locvecind;
  M_locvecind=new int[g_length*M_length_max];   //indexes of non-zero entries
  for (int i=0;i<g_length;i++) {
    M_length[i]=0;
    content* M_i=M+i*g_length;
    content* M_locvecval_i=M_locvecval+i*M_length_max;
    int* M_locvecind_i=M_locvecind+i*M_length_max;
    for (int j=0;j<g_length;j++) {
      if (abs(M_i[j])!=0) {
        M_locvecval_i[M_length[i]]=M_i[j];
        M_locvecind_i[M_length[i]++]=j;
        if (M_length[i]>M_length_max) {
          printf("loc vec in convertM2locvec grown beyond limits (%i)\n",M_length_max); exit(1);}
      }
    }
  }
}



void region::inverseLUtimesM(content* x, content* y)
{
  M_sparse->operate(x,y);
  HG_sparse->LU_banded_solve(y,y);
}

extern "C"
{
  void F77NAME(zgetrs)( char*, int* , int* , content* , int* , int* , content*, int*, int*);
}

void region::inverseLUtimesM_with_checking(content* x, content* y)
{
  //fprintf(logg,"entering inverseLUtimesM (with g_length %i)\n", g_length);
  //first z = M x
  for (int i=0; i<g_length; i++) {
    y[i]=0;
    content* M_locvecval_i=M_locvecval+i*M_length_max;
    int* M_locvecind_i=M_locvecind+i*M_length_max;
    int length=M_length[i];
    for (int j=0;j<length;j++) y[i]+=M_locvecval_i[j]*x[M_locvecind_i[j]];
    //fprintf(logg,"line %i: input (%.3e%+.3ei) through locvec of length %i results in (%.3e%+.3ei)\n",i,x[i].real(), x[i].imag(), length, y[i].real(),y[i].imag());
  }
  complex<double>* y_alt=new complex<double>[g_length];
  M_sparse->operate(x,y_alt);
  double res=0;
  for (int i=0; i<g_length; i++) res+=abs(y[i]-y_alt[i]);
  if (res>-1e-15) fprintf(logg,"LU inverse check: y and y_alt (M result) difference %e\n",res); 
  
  //then y=A^-1 z 
  content *z=new content[g_length*2];
  content *z2=z+g_length;
  for (int i=0; i<g_length; i++) z[i]=y[i];
  //for (int i=0; i<g_length; i++) fprintf(logg,"y[%i] before zgetrs call = (%e%+ei)\n", i, y[i].real(), y[i].imag());
  //for (int i=0; i<g_length; i++) for (int j=0; j<g_length; j++) fprintf(logg,"HG[%i,%i]=(%e%+ei) and LU[%i,%i]= (%e%+ei)\n", i, j, HG[i*g_length+j].real(), HG[i*g_length+j].imag(), i, j, LU[i*g_length+j].real(), LU[i*g_length+j].imag());
  char TRANS='T';
  int NRHS=1, INFO;
  F77NAME(zgetrs)(&TRANS, &g_length, &NRHS, LU, &g_length, IPIV, y, &g_length, &INFO);
  if (INFO!=0) {printf("LU-factorized operator could not be inverted: zgetrs INFO=%i\n",INFO);exit(1);}
  res=0;
  for (int i=0;i<g_length;i++) {
    z2[i]=0;
    for (int j=0;j<g_length;j++) z2[i]+=HG[i*g_length+j]*y[j];
    res+=abs(z2[i]-z[i]);
  }
  if (res>1e-12) fprintf(logg,"check on the LU inverse: residuum of the difference of the HG * HG^-1 * M x = %e\n",res);
  
  HG_sparse->LU_banded_solve(y_alt,y_alt);
  res=0;
  for (int i=0; i<g_length; i++) res+=abs(y[i]-y_alt[i]);
  if (res>-1e-15) fprintf(logg,"LU inverse check: y and y_alt difference (HG^-1 *M result) %e\n",res); 
  delete y_alt;
  delete z;
}

//check that for each eigenvector 
//1.    A x lambda = M x
//2.    A x = lambda^-1 M x
void region::graphene_direct_check(int eig, content* vecs, content *vals)
{
  content* z=new content[g_length*3];
  content* Mx=z+g_length;
  content* Mxp=z+2*g_length;

  //eigenvectors
  for (int n=0;n<eig;n++) {
    content* vec=vecs+n*g_length;
    content val=vals[n];

    //0.   lambda x = OP x
    inverseLUtimesM(vec,z);
    double res=0;
    for (int i=0;i<g_length;i++) res+=abs(val*vec[i]-z[i]);
    fprintf(logg,"check #0 for eigenvector %i: lambda x - OP x = %e\n",n,res);
    
    //1.   lambda x = LU^-1 M x
    res=0;
    for (int i=0;i<g_length;i++) {
      Mx[i]=0;
      Mxp[i]=0;
      for (int j=0;j<g_length;j++) {
        Mx[i]+=M[i*g_length+j]*vec[j];
        //if (abs(M[i*g_length+j])>0) fprintf(logg,"nonzero element M(%i,%i)=(%e%+ei)\n", i, j, M[i*g_length+j].real(), M[i*g_length+j].imag());
      }
      for (int j=0;j<M_length[i];j++) {
        Mxp[i]+=M_locvecval[i*M_length_max+j]*vec[M_locvecind[i*M_length_max+j]];      
        //fprintf(logg,"nonzero loc vec element locvec(%i,%i)->%i = (%e%+ei)\n", i, j, M_locvecind[i*M_length_max+j], M_locvecval[i*M_length_max+j].real(), M_locvecval[i*M_length_max+j].imag());
      }
      z[i]=Mx[i];
      res+=abs(Mx[i]-Mxp[i]);
    }
    fprintf(logg,"Mx and Mxp are equal up to %e\n",res);
    char TRANS='T';
    int NRHS=1, INFO;
    F77NAME(zgetrs)(&TRANS, &g_length, &NRHS, LU, &g_length, IPIV, z, &g_length, &INFO);
    res=0;
    for (int i=0;i<g_length;i++) res+=abs(val*vec[i]-z[i]);
    fprintf(logg,"check #1 for eigenvector %i: lambda x - LU^-1 M x = %e\n",n,res);

    //2.    HG x = lambda^-1 M x    
    for (int i=0;i<g_length;i++) {
      z[i]=0;
      for (int j=0;j<g_length;j++) z[i]+=HG[i*g_length+j]*vec[j];
    }
    res=0;
    for (int i=0;i<g_length;i++) res+=abs(z[i]-Mx[i]/val);
    fprintf(logg,"check #2 for eigenvector %i: HG x - lambda^-1 M x = %e\n",n,res);
    
    //scalar product
    for (int m=0;m<eig;m++) {
      content* My=Mxp;
      content* y=vecs+m*g_length, *x=vec;
      for (int i=0;i<g_length;i++) {
        My[i]=0;
        for (int j=0;j<g_length;j++) My[i]+=M[i*g_length+j]*y[j];
      }
      //x, y, Mx, My
      content r;
      //r=0; for (int i=0;i<g_length;i++) r+=conj(x[i])*y[i];
      //fprintf(logg,"< v[%i] | v[%i] > = %e%+ei\n",n,m,r.real(), r.imag());
      //r=0; for (int i=0;i<g_length;i++) r+=conj(x[i])*My[i];
      //fprintf(logg,"< v[%i] | M v[%i] > = %e%+ei\n",n,m,r.real(), r.imag());
      //r=0; for (int i=0;i<g_length;i++) r+=conj(Mx[i])*y[i];
      //fprintf(logg,"< Mv[%i] | v[%i] > = %e%+ei\n",n,m,r.real(), r.imag());
      r=0; for (int i=0;i<g_length;i++) r+=conj(Mx[i])*My[i];
      fprintf(logg,"< Mv[%i] | Mv[%i] > = %e%+ei\n",n,m,r.real(), r.imag());
    }
  }
  delete z;
  
  //matrixes
  //HG and HG_sparse
  double res=0, tot=0, res_alt=0;
  for (int i=0;i<g_length;i++) {
    for (int j=0;j<g_length;j++) {
      content h1=HG[i*g_length+j], h2=HG_sparse->read(i,j), h3=0;
      if (i-j>=-HG_sparse->KU && i-j<=HG_sparse->KL) 
	h3=HG_sparse->AB[HG_sparse->KL+HG_sparse->KU+i-j+HG_sparse->LDAB*j]; 
      if (abs(h1-h2)>1e-15) fprintf(logg,"HG (%e%+ei) and HG banded (%e%+ei) differ at (i=%i, j=%i) by %e\n",h1.real(), h1.imag(),h2.real(), h2.imag(), i,j,abs(h1-h2));
      //if (abs(h1-h3)>1e-15) fprintf(logg,"HG (%e%+ei) and AB (%e%+ei) differ at (i=%i, j=%i) by %e\n",h1.real(), h1.imag(),h3.real(), h3.imag(), i,j,abs(h1-h3));
      res+=abs(h1-h2);
      res_alt+=abs(h1-h3);
      tot+=abs(h1+h2);
    }
  }
  fprintf(logg,"check: HG differs from HG sparse by %e out of (%e)\n", res, tot);
  //fprintf(logg,"check: HG differs from AB banded by %e out of (%e)\n", res_alt, tot);
  
  //M and M_sparse
  res=tot=0;
  for (int i=0;i<g_length;i++) {
    for (int j=0;j<g_length;j++) {
      content h1=M[i*g_length+j], h2=M_sparse->read(i,j);
      if (abs(h1-h2)>1e-15) fprintf(logg,"M (%e%+ei) and M sparse (%e%+ei) differ at (i=%i, j=%i) by %e\n",h1.real(), h1.imag(),h2.real(), h2.imag(), i,j,abs(h1-h2));
      res+=abs(h1-h2);
      tot+=abs(h1+h2);
    }
  }
  fprintf(logg,"check: M differs from M sparse by %e out of (%e)\n", res, tot);  
}
    

//converts the sigma subspace wavefunction into grid wavefunction inserting values using boundary conditions; the first argument has length Nw_g*sets/2, the other two have lengths Nw
//offset is the dof internal index/2 (cutting away the sigma dof bit)
void region::psip2psig(content *psip, int offset, content* psig_up, content* psig_down)
{
  psip+=offset;
  for (int i=0;i<Nx1;i++) {
    for (int j=0;j<Ny1;j++) {
      grid_point_prototype& gp=grid_point[i*Ny1+j];
      switch (gp.type_graphene) {
        case 0 : continue;
        case 1 : {
          psig_up[gp.sl]=psip[gp.sl_g[0]*sets/2];
          psig_down[gp.sl]=psip[gp.sl_g[1]*sets/2];
          break;
        }
        case 2 : {
          psig_up[gp.sl]=psip[gp.sl_g[0]*sets/2]*graphene_boundary_aux(i,j,0);
          psig_down[gp.sl]=psip[gp.sl_g[1]*sets/2]*graphene_boundary_aux(i,j,1);
          break;
        }
        case 3 : {
          psig_up[gp.sl]=psig_down[gp.sl]=0;
          for (int p=-2;p<3;p++) {
            if (p==0) continue;
            int* sl_g=grid_point[(i+p)*Ny1+(j)].sl_g;
            psig_up[gp.sl]+=psip[sl_g[0]*sets/2]*graphene_center_aux[p+2];
            psig_down[gp.sl]+=psip[sl_g[1]*sets/2]*graphene_center_aux[p+2];
            sl_g=grid_point[(i)*Ny1+(j+p)].sl_g;
            psig_up[gp.sl]+=psip[sl_g[0]*sets/2]*graphene_center_aux[p+2];
            psig_down[gp.sl]+=psip[sl_g[1]*sets/2]*graphene_center_aux[p+2];
          }
          break;
        }        
	case 4 : {
	  psig_up[gp.sl]=psig_down[gp.sl]=0;
	  break;
	}
      }
    }
  }
}



void region::deallocate_graphene()
{
  delete HG_sparse;
  delete M_sparse;
  if (M!=0) delete M;
  if (HG!=0) delete HG;
  if (LU!=0) delete LU;
  if (IPIV!=0) delete IPIV;
  if (M_length!=0) delete M_length;
  if (M_locvecval!=0) delete M_locvecval;
  if (M_locvecind!=0) delete M_locvecind;
}


void region::operate_tot_init()
{
  clean_multiples();
  Nr=N_btp=0;
  overall_shift=0;
}


void region::operate_tot_init(int w_f_x, int w_f_y, operators op, content koef, content f(prec,prec),content exp_integral(prec, prec, prec, prec))
{
  if (op<all_op) {//is it diference operator?
    int x=op_mat_dim[op][act_prec].x;
    int y=op_mat_dim[op][act_prec].y;

    dn_point_operation_tot(w_f_x,w_f_y, -x,x,-y,y,p_op_mat[op][act_prec],f,koef,exp_integral);
    return;
  }

  else {//it is ! a diference operator
    switch (op) {
      case(lin_com): {
        prec aux=1;
        dn_point_operation_tot(w_f_x,w_f_y,0,0,0,0,&aux,f,koef,exp_integral);
        break;
      }
      default: {fprintf(logg,"non - valid value in operate\n");exit(1);}
    }
  }
}


void region::operate_tot(content *input, content* output)
{
  if (relative_point_method=='g') extrapolate_relative_points(input);
  for (int a=0;a<Nw*sets;a++) output[a]=input[a]*overall_shift;
  content *input_aux,*output_aux;
  for (int from=0;from<sets;from++) {
    input_aux=input+Nw*from;
    for (int to=0;to<sets;to++) {
      output_aux=output+Nw*to;
      //local_vector_prototype *plm=local_vector[from*sets+to];
      for (int a=0;a<Nw;a++) {//through all active points
        //local_vector_prototype &lm=plm[a];
        local_vector_prototype &lv=local_vector[from*sets+to][a];
        //int length=lm.length;
        for (int l=0;l<lv.length;l++) *output_aux+=lv.coef[l]*input_aux[lv.sl[l]];
        output_aux++;
      }
      //through all boundary touching points
      for (int i=0;i<N_btp;i++) {
        local_vector_prototype &lv=local_vector[from*sets+to][sl_btp[i]];
        output_aux=output+Nw*to+sl_btp[i];

        if (relative_point_method=='r') {
          #ifdef DEBUG3
          fprintf(logg,"boundary condition point#%i (from % i to %i) output point %i: number of related special points %i, their list:\n",i, from, to, sl_btp[i], lv.length_s);
          for (int j=0;j<lv.length_s;j++) fprintf(logg,"\tspecial active point %i (coef %.3e%+.3e, bc_coef %.3e%+.3e )\n",lv.sl_s[j],lv.coef_s[j].real(), lv.coef_s[j].imag(), bc_coefficient.real(), bc_coefficient.imag());
          #endif
          int gi=active_point[sl_btp[i]].coordinate.i;
          prec x=xd(gi);
          //if (x>0) 
          for (int j=0;j<lv.length_s;j++) *output_aux+=lv.coef_s[j]*input_aux[lv.sl_s[j]]*bc_coefficient;
        }
        else {
#ifdef DEBUG3
            fprintf(logg,"boundary touching point#%i (from % i to %i) output point %i: number of relative points %i, their list:\n",i, from, to, sl_btp[i], lv.length_r);
            for (int r=0;r<lv.length_r;r++) fprintf(logg,"\trelative point %i (coef %.3e%+.3e, fe %.3e%+.3e )\n",lv.sl_r[r],lv.coef_r[r].real(), lv.coef_r[r].imag(), relative_point[lv.sl_r[r]].fe[from].real(), relative_point[lv.sl_r[r]].fe[from].imag());
#endif       
          for (int r=0;r<lv.length_r;r++) *output_aux+=lv.coef_r[r]*relative_point[lv.sl_r[r]].fe[from];//it must be to==from
        }
      }
    }
  }
}

prec region::wave_velocity(prec k)
{
  bool ok;
  for (int label=0;label<Nw;label++) {
    int i,j;
    sl2c(label,i,j);
    if (i<b || i>Nxb1 || j<b || j>Nyb1) {printf("wave_velocity: i,j out of range\n");exit(1);}
    ok=true;
    for (int l=-act_prec-1;l<=act_prec+1;l++) if (grid_point[(i+l)*Ny1+j].type!=1) ok=false;
    if (ok) {
      content *aux_in=new content[Nw];
      content *aux_out=new content[Nw];
      for (int l=-act_prec-1;l<=act_prec+1;l++) {
        if (grid_point[(i+l)*Ny1+j].type!=1) {printf("wave_velocity: l1 out of range\n");exit(1);}
        int l1=c2sl(i+l,j);
        aux_in[l1]=exp(content(0,k*hx*l));
      }
      operate(aux_in,aux_out,0,0,dx,content(1,0),0,set);
      prec res=(conj(aux_in[label])*aux_out[label]).imag();
      //printf("velocity:\n\tderivative at label: %f%+fi, wavefunction: %f%+fi\n",aux_out[label].real(),aux_out[label].imag(),aux_in[label].real(),aux_in[label].imag());
      delete aux_in;
      delete aux_out;
      return(res);
    }
  }
  printf("velocity could not be computed - not enough points inside\nusing crude estimate (k itself)\n");
  return(k);
}



content region::mv_lz_1(prec x, prec y)
{
  return(content(x,0));
}


content region::mv_lz_2(prec x, prec y)
{
  return(content(y,0));
}



content region::mv_r2(prec x, prec y)
{
  return(content(x*x+y*y,0));
}

content region::mv_r4(prec x, prec y)
{ prec w=x*x+y*y;
  return(content(w*w,0));
}

content region::mv_x(prec x, prec y)
{
  return(content(x,0));
}

content region::mv_y(prec x, prec y)
{
  return(content(y,0));
}



void region::op_ket(content *input, content *output, mean_value_of_operator op, setting_or_adding soa)
{

  if ((op==sigmax || op==sigmay || op==sigmaz) && sets!=2) {
    printf("operator %i not defined for sets=%i!=2\n",op, sets);
    exit(1);
  }

  switch (op) {
    case (unit): {
      for (int set=0;set<sets;set++) operate(input,output,set,set,lin_com,content(1,0),0,soa);
      break;
    }
    case (sigmax): {
      operate(input,output,1,0,lin_com,content(1,0),0,soa);
      operate(input,output,0,1,lin_com,content(1,0),0,soa);
      break;
    }
    case (sigmay): {
      operate(input,output,1,0,lin_com,content(0,-1),0,soa);
      operate(input,output,0,1,lin_com,content(0,1),0,soa);
      break;
    }
    case (sigmaz): {
      operate(input,output,0,0,lin_com,content(1,0),0,soa);
      operate(input,output,1,1,lin_com,content(-1,0),0,soa);
      break;
    }
    case (lz): {
      for (int set=0;set<sets;set++) operate(input,output,set,set,dy,content(0,-1),mv_lz_1,soa);
      for (int set=0;set<sets;set++) operate(input,output,set,set,dx,content(0,1),mv_lz_2,add);
      break;
    }
    case (lz1Dx): {
      prec points=Nx-2*b+1;
      for (int set=0;set<sets;set++) operate(input,output,set,set,dx,content(0,-1)*points/(2*M_PI)*hx,0,soa);
      break;
    }
    case (r2): {
      for (int set=0;set<sets;set++) operate(input,output,set,set,lin_com,content(1,0),mv_r2,soa);
      break;
    }
    case (r4): {
      for (int set=0;set<sets;set++) operate(input,output,set,set,lin_com,content(1,0),mv_r4,soa);
      break;
    }
    case (r_x): {
      for (int set=0;set<sets;set++) operate(input,output,set,set,lin_com,content(1,0),mv_x,soa);
      break;
    }
    case (r_y): {
      for (int set=0;set<sets;set++) operate(input,output,set,set,lin_com,content(1,0),mv_y,soa);
      break;
    }
    case (dd_mv): {
      for (int set=0;set<sets;set++) operate(input,output,set,set,dydy,content(1,0),0,soa);
      for (int set=0;set<sets;set++) operate(input,output,set,set,dxdx,content(1,0),0,add);
      break;
    }
    case (invx_mv): {
      for (int set=0;set<sets;set++) operate(input,output,set,set,invx,content(1,0),0,soa);
      break;
    }
    case (invy_mv): {
      for (int set=0;set<sets;set++) operate(input,output,set,set,invy,content(1,0),0,soa);
      break;
    }
    case (inv_mv): {
      for (int set=0;set<sets;set++) operate(input,output,set,set,inv,content(1,0),0,soa);
      break;
    }
    default: printf("non-valid operator in mean_value\n");
  }
}


content region::mean_value(content *input, mean_value_of_operator op, int type)
{
    op_ket(input,temp_mv,op,set);
    content result=0;
    for (int set=0;set<sets;set++) result+=braket(input,set,temp_mv,set,type);
    return(result);
}


content region::braket(content *input1, int w_f_x, content *input2, int w_f_y, int type)
{
  input1+=w_f_x*Nw;
  input2+=w_f_y*Nw;
  content result=0;
  switch (type) {
    case 0 : {for (int i=0;i<Nw;i++) result+=conj(input1[i])*input2[i]*active_point[i].volume;break;}
    case 1 : {for (int i=0;i<Nw;i++) result+=input1[i]*input2[i]*active_point[i].volume;break;}
    default : {printf("wrong type %i in region::braket\n",type);exit(1);}
  }
  return(result);
}

content region::braketN(int N, content **input, bool* inverse, bool* conjugate)
{
  content result=0;
  content aux[10];
  if (N<1 || N>10) {printf("wrong number of vectors %i in region::braket\n",N);exit(1);}
  for (int i=0;i<Nw;i++) {
    for (int j=0;j<N;j++) {
      aux[j]=input[j][i];
      if (inverse[j]) aux[j]=1.0/aux[j];
      if (conjugate[j]) aux[j]=conj(aux[j]);
    }
    for (int j=1;j<N;j++) aux[0]*=aux[j];
    result+=aux[0]*active_point[i].volume;
  }
  return(result);
}

void region::ipr_and_com(content *input, double& ipr, double& x0)
{ 
    double sum2=0, sum4=0, sum2x=0;
    for (int k=0;k<Nw;k++) {
      double aux=(input[k]*conj(input[k])).real();
      double x,y;
      s2coordinate(k,x,y);
      sum2+=aux;
      sum2x+=aux*x;
      sum4+=aux*aux;
    }
    ipr=sum2*sum2/sum4*hx;	//inverse participation ratio
    x0=sum2x/sum2;		//center of mass
}

void region::renormalize_to_volume(content *input, int w_f_x)
{
  input+=w_f_x*Nw;
  int imax=0;
  prec max=0, sum=0;
  for (int i=0;i<Nw;i++) {
    prec val=(conj(input[i])*input[i]).real()*active_point[i].volume;
    sum+=val;
    if (val>max) {val=max;imax=i;}
  }
  content phase=input[imax]/abs(input[imax])*sqrt(sum);
  phase=1.0/phase;
  if ((input[Nw-1]*phase).real()<0) phase*=-1;
  for (int i=0;i<Nw;i++) input[i]*=phase;
}

void region::braket_vec(content *input1, content *input2,content *output, setting_or_adding soa)
{
  content value;
  for (int w=0;w<sets;w++) {
    input1+=w*Nw;
    input2+=w*Nw;
    for (int i=0;i<Nw;i++) {
      value=conj(input1[i])*input2[i]*active_point[i].volume;
      if (soa==add || w!=0) output[i]+=value; else output[i]=value;
    }
  }
}

void region::braket_vec(double *input1, content *input2,content *output, setting_or_adding soa)
{
  content value;
  for (int w=0;w<sets;w++) {
    input1+=w*Nw;
    input2+=w*Nw;
    for (int i=0;i<Nw;i++) {
      value=input1[i]*input2[i]*active_point[i].volume;
      if (soa==add || w!=0) output[i]+=value; else output[i]=value;
    }
  }
}


void region::new_values()
{

  for (int i=0;i<all_op;i++) {
    for (int j=0;j<sets;j++) {
      if (dif_op_filled[i][j]) {
        delete dif_op_vals[i][j];
        dif_op_filled[i][j]=false;
      }
    }
  }
}


void region::constants(prec l2nm)
{
  printf("Lx0:%e\tLx:%e\thx:%e\tkx:%e\n",Lx0*l2nm,Lx*l2nm,hx*l2nm,kx/l2nm);
  printf("b:%i\tNxb:%i\tNxb1:%i\n",b,Nxb,Nxb1);
  printf("Nx:%i\tNx1:%i\tpbc_x:%i\ttab_nx:%i\n",Nx,Nx1,pbc_x,tab_nx);

  printf("Ly0:%e\tLy:%e\thy:%e\tky:%e\n",Ly0*l2nm,Ly*l2nm,hy*l2nm,ky/l2nm);
  printf("b:%i\tNyb:%i\tNyb1:%i\n",b,Nyb,Nyb1);
  printf("Ny:%i\tNy1:%i\tpbc_y:%i\ttab_ny:%i\n",Ny,Ny1,pbc_y,tab_ny);

  printf("Np:%i\tNw:%i\n",Np,Nw);

}


//walk through the created net
extern bool dot2waveguide(prec x, prec y, prec& xw, prec& yw);
//extern content aux_Z(prec x, prec y);
extern prec potential_waveguide(prec x);
void region::ant(function_potential f, double l2nm, double e2meV)
{
  bool inversions=false;
  bool conversions=true;
  bool waveguide=false;
  bool local_vectors=false;
  bool graphene=false;
  int sets=0;
  int sigma=0;

  coordinate_prototype pos={1,1};
  char buf[2]=" ";
  char znak=1;
  int label,lbl=0;
  prec x=0,y=0;
  coordinate_prototype posIx, posIy, posI;
  int lblIx=0,lblIy=0,lblI=0;

  //bool active_inside;

  do {
    printf("Up/Down/Left/Right / Grid info/Inv/Conv/Waveguide/locVec / graphenE / Quit\n");
    switch (znak) {
      case('u'): {if (pos.j<Ny) pos.j++;break;}
      case('d'): {if (pos.j>0) pos.j--;break;}
      case('l'): {if (pos.i>0) pos.i--;break;}
      case('r'): {if (pos.i<Nx) pos.i++;break;}
      case('g'): {constants(l2nm);break;}
      case('i'): {inversions=1-inversions;break;}
      case('e'): {if (!graphene) graphene=true; else
        {sigma++;if (sigma==2) {sigma=0;graphene=false;}}
        break;}
      case('c'): {conversions=1-conversions;break;}
      case('w'): {waveguide=1-waveguide;break;}
      case('v'): {if (!local_vectors) local_vectors=true; 
      else {sets++; if (sets>=sets*sets) {sets=0;local_vectors=false;}}
        break;}
      case('q'): {printf("\n");return;}
      default: continue;
    };
    invert(pos.i,pos.j,posIx.i, posIx.j,invx);
    invert(pos.i,pos.j,posIy.i, posIy.j,invy);
    invert(pos.i,pos.j,posI.i, posI.j,inv);
    
    coordinate_prototype active=grid_point[pos.i*Ny1+pos.j].active_coordinate;
    
    printf("\n");
    
    for(int j=Ny;j>=0;j--) {
      for(int i=0;i<Nx1;i++) {
        label=c2sl(i,j);
        int type=grid_point[i*Ny1+j].type;
        int sector=grid_point[i*Ny1+j].sector;
        buf[0]='\0';
        if (type==0) buf[0]='*';
        if (type==2) buf[0]='-';
        if (type==1) buf[0]='0'+sector;
        if (i==posIx.i && j==posIx.j) {buf[0]='X';lblIx=label;} //Ix
        if (i==posIy.i && j==posIy.j) {buf[0]='Y';lblIy=label;} //Iy
        if (i==posI.i && j==posI.j) {buf[0]='I';lblI=label;} //I
        if (i==pos.i && j==pos.j) {buf[0]='M';lbl=label;} //me
        message(buf,1);
      }
      printf("\n");
    }
    grid_point_prototype& me=grid_point[pos.i*Ny1+pos.j];
    //point info
    printf("point: sequential label=%i ",lbl);
    prec pot=f(xd(pos.i),yd(pos.j));
    printf("grid coord:(%i,%i)=(%.2e nm,%.2e nm) (potential:%.4e [meV])\n",pos.i,pos.j,xd(pos.i)*l2nm,yd(pos.j)*l2nm,pot*e2meV);
    
    printf("xline[%i] (active xlines:%i-%i, tot:%i): fb=%.4e, fpc=%i, fpsl=%i, np=%i\n",pos.i, first_active_xline, last_active_xline, lines, xline[pos.i].fb, xline[pos.i].fpc,xline[pos.i].fpsl,xline[pos.i].np);
    
    if (me.type==1) printf("point type %i (volume %.2e), sector %i:",me.type,active_point[me.sl].volume, me.sector);
    else printf("point type %i, sector %i:",me.type, me.sector);
    int sec=me.sector;
    if (sec>-1) printf("grid limits:(%i-%i), sector on grid line[%i]:offset=%i,dim=%i ", sector[sec].imin, sector[sec].imax, pos.i, sector[sec].offset[pos.i-sector[sec].imin], sector[sec].dim[pos.i-sector[sec].imin]);
    else printf("not inside ");
    printf("nearby sectors: ");
    if (sectorized) {
      if (me.type==1) {//an active point
        for (int l=0;l<Nl;l++) {if (xline[pos.i].active_line_point[active.j].nearby_sector[l]) printf("%i ",l);}
        printf("nearby any:%i\n",xline[pos.i].active_line_point[active.j].nearby_any_sector);
      }
      else printf ("not an active point\n");
    }
    else printf("not initialized\n");
    
    if (conversions) {
      int ai,aj,si,sj,auxs,auxi,auxj;
      grid2active(pos.i, pos.j, ai, aj);
      grid2sec(pos.i,pos.j,auxs,si,sj);
      printf("grid/sec2grid/active2grid: (%i,%i)", pos.i, pos.j);
      if (sec==-1) printf("/not an active point\n");
      else {
        sec2grid(sec,si,sj,auxi,auxj);printf("/(%i,%i)",auxi,auxj);
        active2grid(ai,aj,auxi,auxj);printf("/(%i,%i)\n",auxi,auxj);
      }
      printf("sec/grid2sec/active2sec: (%i,%i:%i)",si, sj, sec);
      grid2sec(pos.i,pos.j,auxs,auxi,auxj);printf("/(%i,%i:%i)",auxi,auxj,auxs);
      if (sec==-1) printf("/not an active point\n");
      else {
        active2sec(ai,aj,auxs,auxi,auxj);printf("/(%i,%i:%i)\n",auxi,auxj,auxs);
      }
      printf("active/grid2active/sec2active: (%i,%i)",active.i, active.j);
      grid2active(pos.i,pos.j,auxi,auxj);printf("/(%i,%i)",auxi,auxj);
      if (sec==-1) printf("/not an active point\n");
      else {
        sec2active(sec,si,sj,auxi,auxj);printf("/(%i,%i)\n",auxi,auxj);
      }
    }
    if (inversions) {
      printf("1 =[(%2i,%2i)(%+.2e,%+.2e)]",pos.i,pos.j,xd(pos.i)*l2nm,yd(pos.j)*l2nm);
      if (lblIx==0) printf("(bb)");
      printf(" ");
      printf("I =[(%2i,%2i)(%+.2e,%+.2e)]",posI.i,posI.j,xd(posI.i)*l2nm,yd(posI.j)*l2nm);
      if (lblI==0) printf("(bb)");
      printf("\n");
      printf("Ix=[(%2i,%2i)(%+.2e,%+.2e)]",posIx.i,posIx.j,xd(posIx.i)*l2nm,yd(posIx.j)*l2nm);
      if (lblIx==0) printf("(bb)");
      printf(" ");
      printf("Iy=[(%2i,%2i)(%+.2e,%+.2e)]",posIy.i,posIy.j,xd(posIy.i)*l2nm,yd(posIy.j)*l2nm);
      if (lblIy==0) printf("(bb)");
      printf("\n");
    }
    if (waveguide) {
      prec xw,yw;
      bool wg_ins=dot2waveguide(xd(pos.i),yd(pos.j),xw,yw);
      printf("waveguide coord: (%.2e,%.2e) (inside:%i, wg_potential %e)\n", xw*l2nm, yw*l2nm, wg_ins, potential_waveguide(xw)*e2meV);
    }
    if (local_vectors) {
      printf("local vector:");
      if (me.type==0) printf("point out\n");
      if (me.type==1) {
        bool btp=false;
        for (int i=0;i<N_btp;i++) if (sl_btp[i]==me.sl) btp=true;
        printf("active point [set->set=%i] (among btp's:%i):\n",sets, btp);
        local_vector_prototype lv=local_vector[sets][me.sl];
        printf("%i related active points:",lv.length);
        for (int i=0;i<lv.length;i++) printf("%i (c=%.2e) ",lv.sl[i],lv.coef[i].real()); printf("\n");
        printf("%i related relative points:",lv.length_r);
        for (int i=0;i<lv.length_r;i++) printf("%i (c=%.2e) ",lv.sl_r[i],lv.coef_r[i].real()); printf("\n");
        printf("%i related special points:",lv.length_s);
        for (int i=0;i<lv.length_s;i++) printf("%i (c=%.2e) ",lv.sl_s[i],lv.coef_s[i].real()); printf("\n");
      }
      if (me.type==2) {
        relative_point_prototype rp=relative_point[me.sl];
        printf("relative point %i (method %c):\ncoors (%i-%i),", me.sl, relative_point_method, rp.coordinate.i, rp.coordinate.j);
        if (relative_point_method=='r') {
          printf("solved for=%i, special point %i (c=%.2e), %i related active points\n", rp.solved, rp.sl0, rp.c1.real(),rp.length);
          for (int i=0;i<rp.length;i++) printf("%i (c=%.2e) ",rp.sl[i],rp.coef[i].real()); printf("\n");
        }
        if (relative_point_method=='g') {
          printf("re=%.2e ri=%.2e [nm], interpolated=%i by %i active points:\n", rp.re*l2nm, rp.ri*l2nm, rp.solved, rp.length2);
          for (int i=0;i<rp.length2;i++) printf("%i (c=%.2e) ",rp.sl2[i],rp.coef2[i].real()); printf("\n");
        }      
      }
    }
    if (graphene) {
      const char* mess[]={"out","inside","border","central","first/last"};
      if (me.type_graphene!=2) printf("graphene tags: type=%i(%s), sl_up=%i, sl_down=%i, plaq=%i (Nw_g=%i, Nw_plaq=%i)\n", me.type_graphene, mess[me.type_graphene], me.sl_g[0], me.sl_g[1], me.plaq, Nw_g, Nw_plaq);
      else printf("graphene tags: type=%i(%s:boundary vector (%i,%i)), sl_up=%i, sl_down=%i, plaq=%i (Nw_g=%i, Nw_plaq=%i)\n", me.type_graphene, mess[me.type_graphene], me.b_sum_i, me.b_sum_j, me.sl_g[0], me.sl_g[1], me.plaq, Nw_g, Nw_plaq);
      if (me.plaq!=-1) {
	printf("sig-out=%i subspace:\n",sigma);
	int output=me.plaq*2+sigma;
	if (M_length==0) printf("matrix M undefined");
	else {
	  printf("matrix M entries:\n ");
	  for (int i=0;i<M_length[output];i++) printf("(%i)-(%.2e%+.2ei) ",M_locvecind[output*M_length_max+i], M_locvecval[output*M_length_max+i].real(), M_locvecval[output*M_length_max+i].imag());
	}
	if (HG==0) printf("\nmatrix HG undefined");
	else {
	  printf("\nmatrix HG entries:\n: ");
	  for (int i=0;i<g_length;i++) {
	    content HGij=HG[g_length*output+i];
	    if (abs(HGij)>0) printf("(%i)-(%.2e%+.2ei) ",i,HGij.real(), HGij.imag());
	  }
	}
	if (M_sparse->N==0) printf("\nmatrix M_sparse undefined");
	else {
	  printf("\nmatrix M_sparse (this line length %i) entries:\n ",M_sparse->length[output]);
	  for (int i=0;i<M_sparse->length[output];i++) printf("(%i)-(%.2e%+.2ei) ",M_sparse->label[output][i], M_sparse->value[output][i].real(), M_sparse->value[output][i].imag());
	}
	if (HG_sparse->N==0) printf("\nmatrix HG_sparse undefined");
	else {
	  printf("\nmatrix HG_sparse (this line length %i) entries:\n ",HG_sparse->length[output]);
	  for (int i=0;i<HG_sparse->length[output];i++) printf("(%i)-(%.2e%+.2ei) ",HG_sparse->label[output][i], HG_sparse->value[output][i].real(), HG_sparse->value[output][i].imag());
	}
	printf("\n");
      }
    }
  }
  while((znak=u_getch())!='q');
}

//test diferential operator
//output[15]: norm+max of function, norm+max of derivation, norm+max of measured der.
//norm+max of (der-mder), norm+max of (der-mder)/der
void region::test(operators op, content func(prec,prec), content der(prec, prec), prec *output)
{
  content* buf1=new content[4*Nw];
  content* buf2=buf1+Nw;
  content* buf3=buf1+2*Nw;
  content* buf4=buf1+3*Nw;

  operate(0,buf1,0,0,lin_com,content(1,0),func,set);
  operate(0,buf2,0,0,lin_com,content(1,0),der,set);

  if (operators_sym1) {
    operators_sym1=false;
    operate_tot_init();
    operate_tot_init(0,0,op,content(1,0),0,0);
    operate_tot(buf1,buf3);
    operate_tot(buf1,buf4);
    operate(buf2,buf4,0,0,lin_com,content(-1,0),0,add);
    operators_sym1=true;
  }
  else {
    operate(buf1,buf3,0,0,op,content(1,0),0,set);
    operate(buf1,buf4,0,0,op,content(1,0),0,set);
    operate(buf2,buf4,0,0,lin_com,content(-1,0),0,add);
  }

  int count=0,countz=0;
  content val[5];
  content mx[5]={0};
  content norm[5]={0};

  //fprintf(logg,"test: range(b=%i, Nxb=%i, Nyb=%i\n",b,Nxb,Nyb);

  for (int i=b+b;i<Nxb-b;i++) {
    for (int j=b+b;j<Nyb-b;j++) {
      if (grid_point[i*Ny1+j].type!=1) continue;
      int index=c2sl(i,j);
      //fprintf(logg,"test: (i=%i/j=%i)->index=%i\n",i,j,index);
      val[0]=buf1[index];		//function
      val[1]=buf2[index];		//der
      val[2]=buf3[index];	//mder
      val[3]=buf4[index];	//der-mder
      if (fabs(val[1].real())>1e-10 && fabs(val[1].imag())>1e-10) val[4]=content(val[3].real()/val[1].real(),val[3].imag()/val[1].imag());
      else { val[4]=0; countz++;}

      for (int k=0;k<5;k++) {
        norm[k]+=val[k]*conj(val[k]);
        if (fabs(val[k].real())>fabs(mx[k].real())) mx[k]=content(val[k].real(),mx[k].imag());
        if (fabs(val[k].imag())>fabs(mx[k].imag())) mx[k]=content(mx[k].real(),val[k].imag());
      }
      count++;
    }
  }
  for (int k=0;k<5;k++) {
    output[k*3+0]=sqrt(norm[k].real())/(count-countz*(k==4));
    output[k*3+1]=mx[k].real();
    output[k*3+2]=mx[k].imag();
  }

  delete buf1;
}

extern void spinor2angles(content* xi, double* angles);

void region::draw(complex<double> *input,int w_f, int what, content zero, int only_inside, FILE* file, int where_to)
{ 
  int how=1;
  if (what<0) {what*=-1;how=-1;}
  input+=w_f*Nw;
  int label;
  complex<double> value;
  double out;
  what*=2;

  if (how==1) {//as a table
    for (int k=0;k<3;k++) {
      what/=2;
      if (what%2==0) continue;
      for (int j=0;j<Ny1;j++) {
        for (int i=0;i<Nx1;i++) {
          label=ac2sl(i,j);
          if (grid_point[i*Ny1+j].type==0) {
            if (only_inside) continue;
            else value=zero;
          }
          else value=input[label];
          if (k==0) out=value.real();
          else if (k==1) out=value.imag();
          else out=sqrt(value.real()*value.real()+value.imag()*value.imag());
          sprintf(buffer,"%e ",out);
          splachni(file,buffer,where_to);
        }
        message(file,"\n",where_to);
      }
    }
  }
  else {//as a column with coordinates
  //message(file,"x\t\ty\t\tRe(Psi)\t\tIm(Psi)\n\n\n",where_to);
    for (int j=b;j<Nyb1;j++) {
      for (int i=b;i<Nxb1;i++) {
        int label=ac2sl(i,j);
        if (grid_point[i*Ny1+j].type==0) {
          if (only_inside) continue;
          else value=zero;
        }
        else value=input[label];
//        if (abs(value)<1e-4) continue;
        //double theta,phi;
        //content xi[2];
        //if (label==0) xi[0]=xi[1]=content(0,0);
        //else if (w_f==0) {xi[0]=value;xi[1]=input[label+Nw];}
        //else {xi[1]=value;xi[0]=input[label-Nw];}
        //spinor2angles(xi,theta,phi);
        sprintf(buffer,"%+.6e\t%+.6e\t%+.6e\t%+.6e\t%+.6e\n",xd(i),yd(j),value.real(),value.imag(),abs(value));
        //sprintf(buffer,"%+.6e\t%+.6e\t%+.6e\t%+.6e\n",xd(i),yd(j),theta/M_PI,phi/M_PI);
        splachni(file, buffer, where_to);
      }
      message(file,"\n",where_to);
    }
  }
}

void region::draw_function(content (f)(prec,prec),int what, int only_inside, FILE* file, int where)
{
  prec max=0;
  for (int n=0;n<Nw;n++) {
    int i,j;
    sl2c(n,i,j);
    prec loc=(f(xd(i),yd(j))).real();
    if (loc>max) max=loc;
  }
  
  content* aux=new content[Nw];
  operate(0,aux,0,0,lin_com,content(1,0),f,set);
  draw(aux,0,what,max,only_inside,file,where);
  delete aux;
}



//converts a linear vector(input, w_f) of points inside region into a 2d (nx times ny = tab_nx times tab_ny) table (output)
//the table is indexed following: point from region with coordinates ~(i*hx,j*hy) is stored in the table at point (i,j) -> i*ny+j  
//on the boundary the one row of zeros is stored
void region::region2table(content *input, int w_f, content *output, int nx, int ny)
{
  input+=Nw*w_f;
  if (nx*ny!=tab_ny*tab_nx) {
    printf("Wrong dimensions in region2table!!!\n");
    exit(1);
  }
  for (int i=0;i<nx;i++) {
    for (int j=0;j<ny;j++) {
      output[i*ny+j]=0;
    }
  }
  int offset=b-1;
  for (int i=offset;i<Nxb+2;i++) {
    for (int j=offset;j<Nyb+2;j++) {
      if (grid_point[i*Ny1+j].type==0) continue;
      int label=ac2sl(i,j);
      output[(i-offset)*ny+(j-offset)]=input[label];
    }
  }
}

//converts 2d (nx times ny = tab_nx times tab_ny) table (input) into a linear vector representing the points inside the region (output)
//the boundary in the table is ignored
void region::table2region(content *input, content *output, int w_f,int nx, int ny)
{
  output+=Nw*w_f;
  if (nx*ny!=tab_ny*tab_nx) {
    printf("Wrong dimensions in table2region!!!\n");
    exit(1);
  }

  int i,j;
  int offset=b-1;
  for (int label=0;label<Nw;label++) {
    sl2c(label, i,j);
    output[label]=input[(i-offset)*ny+j-offset];
  }
}

//for coordinates i=0,...,tab_nx-1, j=0,...,tab_ny-1 gives proper
//vector (omega_x,omega_y) which is from the region (-some value,+some value)
void region::give_omega(int i, int j, prec& omega_x, prec& omega_y)
{
  if (i<0 || i>=tab_nx || j<0 || j>=tab_ny) 
     printf("strange coordinates in give_omega:i=%i, j=%i, tab_nx=%i, tab_ny=%i\n",i,j,tab_nx,tab_ny);
  i=mod_pm(i,tab_nx);
  j=mod_pm(j,tab_ny);
  omega_x=i*kx;
  omega_y=j*ky;
}

void region::value_at_init(content** f, content* input)
{
  for (int i=1;i<6;i++) {
    f[i]=new content[Nw];
    operators op=dx;
    switch (i) {
        case 1 : {op=dx;break;}
        case 2 : {op=dy;break;}
        case 3 : {op=dxdy;break;}
        case 4 : {op=dxdx;break;}
        case 5 : {op=dydy;break;}
    }
    dn_point_operation(input,f[i],-op_mat_dim[op][act_prec].x,op_mat_dim[op][act_prec].x,-op_mat_dim[op][act_prec].y,op_mat_dim[op][act_prec].y,p_op_mat[op][act_prec]);
  }
  f[0]=input;
}

void region::value_at_release(content **f)
{
  for (int i=1;i<6;i++) delete f[i];
}

content region::value_at(prec x, prec y, content** f)
{
  int i,j;
  bool inside=nearest_grid_point(x,y,i,j);
  if (!inside) return(content(0));
  int xi,yi;
  
  content res=0;
  prec dx=x-xd(i), dy=y-yd(j);
  if (grid_point[i*Ny1+j].type==0) return(0);
  int label=ac2sl(i,j);
  
  res=f[0][label]+dx*f[1][label]+dy*f[2][label];
  res+=f[3][label]*dx*dy+f[4][label]*dx*dx/2.0+f[5][label]*dy*dy/2.0;

  return(res);
}

//returns grid coordinates of the nearest point from point defined by coordinates x, y
//if the grid point is outside the region, returns inside=false, sl=0
bool region::nearest_grid_point(prec x, prec y, int &i, int &j)
{
  if (x<xd(0) || x>xd(Nx) ||y<yd(0) || y>yd(Ny)) {//point is out of region
    i=j=0;
    return(false);
  }
  
  i=(int) floor((x-xd(0))/hx);
  j=(int) floor((y-yd(0))/hy);
  if (x-xd(i)>hx/2) i++;
  if (y-yd(j)>hy/2) j++;
  
  prec dx=x-xd(i), dy=y-yd(j);
  if ((abs(dx)-hx/2)>1e-12*hx || (abs(dy)-hy/2)>1e-12*hy) {
    printf("wrong nearest_grid_point:\nx=%f xd(i)=%f i=%i dx=%f hx=%f, abs(dx)-hx/2=%e\n",x,xd(i),i,dx,hx,abs(dx)-hx/2);
    printf("y=%f yd(j)=%f j=%i dy=%f hy=%f, abs(dy)-hy/2=%e\n",y,yd(j),j,dy,hy,abs(dy)-hy/2);
    exit(1);
  }
  return(true);
  //printf("nearest_grid_point:x=%f xd(%i)=%f dx=%f hx=%f, abs(dx)-hx/2=%e\n",x,i,xd(i),dx,hx,abs(dx)-hx/2);
  //printf("y=%f yd(%i)=%f dy=%f hy=%f, abs(dy)-hy/2=%e\n\n",y,j,yd(j),dy,hy,abs(dy)-hy/2);
}

int region::invertconseq(int n, operators op)
{
  if (n<0 || n>Nw) return(-1);
  int i,j, ii, ij;
  sl2c(n,i,j);
  invert(i,j,ii,ij,op);
  return(c2sl(ii,ij));
}

//returns coodinates of closest point corresponding to input point with inversion op 
void region::invert(int i, int j, int& ii, int &ij, operators op)
{
  int xi=1,yi=1;
  switch (op) {
    case (invx) : {xi=-1;break;}
    case (invy) : {yi=-1;break;}
    case (inv) : {xi=yi=-1;break;}
    default : {message("wrong operator in invert\n");}
  }
  bool inside=nearest_grid_point(xd(i)*xi,yd(j)*yi,ii,ij);
  if (!inside) {message("inverted point out of range!\n");ii=ij=0;return;}
  if (xd(i)*xi-xd(ii)>hx/2) ii++;
  if (yd(j)*yi-yd(ij)>hy/2) ij++;
  
  if (fabs(xline[i].fb-xline[ii].fb*xi)>1e-10 || fabs(yline[j].fb-yline[ij].fb*yi)>1e-10) {
    message("coordinates not opposite in inversion!!!\n");
    sprintf(buffer,"coords (%f,%f)-->(%f,%f) with inversion %i\n",xline[i].fb,yline[j].fb,xline[ii].fb,yline[ij].fb,op);
    splachni(logg,buffer,1+4);
    sprintf(buffer,"index coords (%i,%i)-->(%i,%i) \n",i,j,ii,ij);
    splachni(logg,buffer,1+4);
  }
  return;
}




//fill matrix koeficients for differential operators
void region::inicialize_diferential_operators()
{
#define aux(i,x,from,to)\
  { to=new prec[i]; for (int j=0;j<i;j++) to[j]=x*from[j];}

  //dx, dy
  {
    prec vektor[all_prec][2*all_prec+1]={\
        {-1,0,1},\
      {-1,8,0,-8,1},\
      {-1,9,-45,0,45,-9,1},\
      {-1,32.0/3,-56,224,0,-224,56,-32.0/3,1}};
    prec koefs[all_prec]={1.0,-1.0/6,1.0/30,-1.0/140};

    for (int l=0;l<all_prec;l++) {
      aux((2*l+3),1/(2*hx)*koefs[l],vektor[l],p_op_mat[dx][l])
      op_mat_dim[dx][l].x=l+1;        //2 is from definition of F
      op_mat_dim[dx][l].y=0;        //in terms of f (1/4)

      aux(2*l+3,1/(2*hy)*koefs[l],vektor[l],p_op_mat[dy][l])
      op_mat_dim[dy][l].x=0;
      op_mat_dim[dy][l].y=l+1; 
    }
  }
  
  //dxp dxm dyp dym
  {
    prec coors[2*all_prec+1], coefs[2*all_prec+1];
    for (int l=0;l<all_prec;l++) {
      for (int i=-l-1;i<=l+1;i++) coors[i+l+1]=(i-0.5)*hx;
      p_op_mat[dxp][l]=new prec[2*l+3];
      construct_differential_operator(1,2*l+3,coors,p_op_mat[dxp][l]);
      op_mat_dim[dxp][l].x=l+1;
      op_mat_dim[dxp][l].y=0;
      for (int i=-l-1;i<=l+1;i++) coors[i+l+1]=(i-0.5)*hy;
      p_op_mat[dyp][l]=new prec[2*l+3];
      construct_differential_operator(1,2*l+3,coors,p_op_mat[dyp][l]);
      op_mat_dim[dyp][l].x=0;
      op_mat_dim[dyp][l].y=l+1;
      for (int i=-l-1;i<=l+1;i++) coors[i+l+1]=(i+0.5)*hx;
      p_op_mat[dxm][l]=new prec[2*l+3];
      construct_differential_operator(1,2*l+3,coors,p_op_mat[dxm][l]);
      op_mat_dim[dxm][l].x=l+1;
      op_mat_dim[dxm][l].y=0;
      for (int i=-l-1;i<=l+1;i++) coors[i+l+1]=(i+0.5)*hy;
      p_op_mat[dym][l]=new prec[2*l+3];
      construct_differential_operator(1,2*l+3,coors,p_op_mat[dym][l]);
      op_mat_dim[dym][l].x=0;
      op_mat_dim[dym][l].y=l+1;
    }
  }

  //dxdx, dydy
  {
    prec vektor[all_prec][2*all_prec+1]={\
          {1,-2,1},\
      {1,-16,30,-16,1},\
      {2,-27,270,-245*2,270,-27,2},\
      {-9,128,-1008,8064,-7175*2,8064,-1008,128,-9}};
    prec koefs[all_prec]={1.0,-1.0/12,1.0/180,1.0/(9*560)};

    for (int l=0;l<all_prec;l++) {
      aux(2*l+3,2/(2*hx*hx)*koefs[l],vektor[l],p_op_mat[dxdx][l])
      op_mat_dim[dxdx][l].x=l+1;
      op_mat_dim[dxdx][l].y=0;

      aux(2*l+3,2/(2*hy*hy)*koefs[l],vektor[l],p_op_mat[dydy][l])
      op_mat_dim[dydy][l].x=0;
      op_mat_dim[dydy][l].y=l+1;
   }
  }
  
  //dxdxdx, dydydy
  {
    prec vektor[all_prec][2*all_prec+1]={\
      {0,0,0},\
      {-0.5,1,0,-1,0.5},\
      {1,-8,13,0,-13,8,-1},\
      {-7,72,-338,488,0,-488,338,-72,7}};
      prec koefs[all_prec]={1.0,1.0,1.0/8,1.0/240};

      for (int l=0;l<all_prec;l++) {
        aux(2*l+3,2/(2*hx*hx*hx)*koefs[l],vektor[l],p_op_mat[dxdxdx][l])
        op_mat_dim[dxdxdx][l].x=l+1;
        op_mat_dim[dxdxdx][l].y=0;

        aux(2*l+3,2/(2*hy*hy*hy)*koefs[l],vektor[l],p_op_mat[dydydy][l])
        op_mat_dim[dydydy][l].x=0;
        op_mat_dim[dydydy][l].y=l+1;
      }
  }

  //laplace
  for (int l=0;l<all_prec;l++){
    prec matrix[(2*all_prec+1)*(2*all_prec+1)]={0};
    for (int i=0;i<2*l+3;i++) {
      for (int j=0;j<2*l+3;j++) {
        int index=i+j*(2*l+3);
  if (i==l+1) matrix[index]+=p_op_mat[dydy][l][j];
  if (j==l+1) matrix[index]+=p_op_mat[dxdx][l][i];
      }
    }
    aux((2*l+3)*(2*l+3),1,matrix,p_op_mat[dd][l])
    op_mat_dim[dd][l].x=l+1;
    op_mat_dim[dd][l].y=l+1;
  }

  { //dxdy
    //beware!!! matrix elements are being saved in the opposite direction
    //in y as are written here, && in the same direction in the x sence
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //this can be account by (-1)^(#of dy) && it is so done here
    prec matrix[all_prec][(2*all_prec+1)*(2*all_prec+1)]={\
                          {-1,   0,   1,\
                      0,   0,   0,\
                      1,   0,  -1},\
  \
               {-1,   0,   0,   0,   1,\
                 0,  16,   0, -16,   0,\
           0,   0,   0,   0,   0,\
           0, -16,   0,  16,   0,\
           1,   0,   0,   0,  -1},\
   \
          {-2,   0,   0,   0,   0,   0,   2,\
            0,  27,   0,   0,   0, -27,   0,\
      0,   0,-270,   0, 270,   0,   0,\
      0,   0,   0,   0,   0,   0,   0,\
      0,   0, 270,   0,-270,   0,   0,\
      0, -27,   0,   0,   0,  27,   0,\
      2,   0,   0,   0,   0,   0, -2},\
      \
     {-1/112.0,   0,   0,   0,   0,   0,   0,   0,1/112.0,\
       0,8.0/63, 0,   0,   0,   0,   0,-8.0/63,0,\
       0,   0,  -1,   0,   0,   0,   1,   0,   0,\
       0,   0,   0,   8,   0,  -8,   0,   0,   0,\
       0,   0,   0,   0,   0,   0,   0,   0,   0,\
       0,   0,   0,  -8,   0,   8,   0,   0,   0,\
       0,   0,   1,   0,   0,   0,  -1,   0,   0,\
       0,-8.0/63,0,   0,   0,   0,   0,8.0/63, 0,\
       1/112.0,   0,   0,   0,   0,   0,   0,   0,-1/112.0}};

    prec koefs[all_prec]={1.0,-1.0/12,1.0/180,-1.0/5};

    for (int l=0;l<all_prec;l++) {
      aux((2*l+3)*(2*l+3),-1/(4*hx*hy)*koefs[l],matrix[l],p_op_mat[dxdy][l])
      op_mat_dim[dxdy][l].x=l+1;
      op_mat_dim[dxdy][l].y=l+1;
     }
   }
  { //dxdxdy, dxdydy
    prec matrix[all_prec][(2*all_prec+1)*(2*all_prec+1)]={\
                          { 1,  -2,   1,\
                      0,   0,   0,\
                     -1,   2,  -1},\
  \
                {1,   0,  -2,   0,   1,\
                 0, -32,  64, -32,   0,\
           0,   0,   0,   0,   0,\
           0,  32, -64,  32,   0,\
          -1,   0,   2,   0,  -1},\
   \
           {4,   0,   0,  -8,   0,   0,   4,\
            0, -81,   0, 162,   0, -81,  0,\
      0,   0,1620,-3240,1620, 0,   0,\
      0,   0,   0,   0,   0,   0,   0,\
      0,   0,-1620,3240,-1620, 0,   0,\
      0,  81,   0,-162,   0,  81,   0,\
     -4,   0,   0,   8,   0,   0,  -4},\
      \
     {1.0/224,   0,   0,   0,-2.0/224,0,   0,   0,1.0/224,\
       0,-16.0/189,0, 0,32.0/189,0, 0,-16.0/189,0,\
       0,   0,   1,   0,  -2,   0,   1,   0,   0,\
       0,   0,   0, -16,  32, -16,   0,   0,   0,\
       0,   0,   0,   0,   0,   0,   0,   0,   0,\
       0,   0,   0,  16, -32,  16,   0,   0,   0,\
       0,   0,  -1,   0,   2,   0,  -1,   0,   0,\
       0,16.0/189,0,  0,-32.0/189,0,  0,16.0/189,0,\
      -1.0/224,   0,   0,   0,2.0/224,0,   0,   0,-1.0/224}};

    prec koefs[all_prec]={1.0,-1.0/24,1.0/1080,-1.0/10};

    for (int l=0;l<all_prec;l++) {
      aux((2*l+3)*(2*l+3),-2/(4*hx*hx*hy)*koefs[l],matrix[l],p_op_mat[dxdxdy][l])
      op_mat_dim[dxdxdy][l].x=l+1;
      op_mat_dim[dxdxdy][l].y=l+1;
      for (int i=0;i<2*l+3-1;i++) { //transpose matrix ,but along (+x,+y) quadrant axis!
        for (int j=0;j<(2*l+3-i-1);j++) {
        prec val=matrix[l][i+(2*l+3)*j];
        matrix[l][i+(2*l+3)*j]=matrix[l][(2*l+2-j)+(2*l+2-i)*(2*l+3)];
        matrix[l][(2*l+2-j)+(2*l+2-i)*(2*l+3)]=val;
        }
      }
      aux((2*l+3)*(2*l+3),2/(4*hx*hy*hy)*koefs[l],matrix[l],p_op_mat[dxdydy][l])
      op_mat_dim[dxdydy][l].x=l+1;
      op_mat_dim[dxdydy][l].y=l+1;
    }
  }
}
#undef aux

void region::listop(operators op, operators_precision l) 
{
  for (int i=0;i<2*op_mat_dim[op][l].x+1;i++) {
    for (int j=0;j<2*op_mat_dim[op][l].y+1;j++) {
      printf("%+f\t",p_op_mat[op][l][(2*op_mat_dim[op][l].y+1)*i+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void region::inicialize_diferential_operators_test()
{
  for (int l=0;l<all_prec;l++) {
    prec* res=new prec[(2*l+3)*(2*l+3)];
    printf("\noperator dx, precision %i\ncandidate:\n",l);
    
    differential_operator(1,0,l+1,1/hx,res);
    for (int i=0;i<(2*l+3);i++) { 
      for (int j=0;j<(2*l+3);j++) { 
        printf("%+f\t",res[i*(2*l+3)+j]);
      }
      printf("\n");
    }
    printf("\n\n correct:\n");
    listop(dx,(operators_precision) l);
    printf("\n");
    
    printf("\noperator dy, precision %i\ncandidate:\n",l);
    differential_operator(0,1,l+1,1/hy,res);
    for (int i=0;i<(2*l+3);i++) { 
      for (int j=0;j<(2*l+3);j++) { 
        printf("%+f\t",res[i*(2*l+3)+j]);
      }
      printf("\n");
    }
    printf("\n\n correct:\n");
    listop(dy,(operators_precision) l);
    printf("\n");
  }
}
    
/*
      
    
    differential_operator(0,1,l+1,1/hy,p_op_mat[dy][l]);
    differential_operator(2,0,l+1,1/(hx*hx),p_op_mat[dxdx][l]);
    differential_operator(0,2,l+1,1/(hy*hy),p_op_mat[dydy][l]);
    differential_operator(3,0,l+1,1/(hx*hx*hx),p_op_mat[dxdxdx][l]);
    differential_operator(0,3,l+1,1/(hy*hy*hy),p_op_mat[dydydy][l]);
    differential_operator(1,1,l+1,1/(hx*hy),p_op_mat[dxdy][l]);
    differential_operator(2,1,l+1,1/(hx*hx*hy),p_op_mat[dxdxdy][l]);
    differential_operator(1,2,l+1,1/(hx*hy*hy),p_op_mat[dxdydy][l]);
    
    for (int i=0;i<(2*l+3)*(2*l+3);i++) p_op_mat[dd][l][i]=p_op_mat[dydy][l][i]+p_op_mat[dxdx][l][i];
  }
}*/


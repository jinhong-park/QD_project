#ifndef R2d
#define R2d

//using namespace std;

#include "auxiliaries.h"

class region {

typedef double prec;
typedef complex<prec> content;

typedef bool (function_inside)(prec,prec);
typedef prec (function_potential)(prec, prec);
typedef content (function_potentialc)(prec, prec);


public:
enum direction{ left, right, down, up, stay};
enum operators { dx, dy, dxp, dyp, dxm, dym, dd, dxdy, dxdx, dydy, dxdxdy, dxdydy, dxdxdx, dydydy, all_op, lin_com, invx, invy, inv};
enum operators_precision{two, four, six, eight, all_prec};
enum setting_or_adding {set, add};
enum del_cr_app {del, create, append};
enum mean_value_of_operator {unit,sigmax,sigmay,sigmaz,lz,lz1Dx,r2,r4,r_x,r_y,dd_mv,invx_mv, invy_mv,inv_mv,all_mv};

prec inv_angle;                 //angle phi_d along which inversions Ix, Iy are evaluated
void set_inv_angle(prec angle) {inv_angle=angle;}

private:
struct coordinate_prototype {		//(xline, yline) coordinates of point
  int i;			//xline coordinate
  int j;			//yline coordinate
};

struct active_line_point_prototype{//info for each active point on xline
  bool *nearby_sector;//true if nearby sector #i (false for the sector it is in itself)
  bool nearby_any_sector;//other than it is in
  int j;      //y grid coordinate (init:fill_arrays())
};

struct xline_prototype{
  int np;     //# of active points on this line (init:initialize())
  int fpc;      //y grid coordinate of the first active point (init:fill_arrays())
  int fpsl;     //sequential label of the first active point (init:fill_arrays())
  prec fb;      //real length (init:initialize())
  active_line_point_prototype *active_line_point;//active points indexed sequentially (partial init:fill_arrays())
};
xline_prototype *xline;   //indexed by [0..Nx1-1]

struct yline_prototype{
  prec fb;      //real length (init:initialize())
};
yline_prototype *yline;   //ylines, indexed by [0..Ny1-1]

struct sector_prototype {//xline sector limits and yline offsets and dimensions
  int imin, imax; //grid labels
  int* dim;     //number of active points of sector on xline (sector coordinate)
  int* offset;  //sequential label of the first sector point on xline (sector coordinate)
};
sector_prototype *sector;//indexed [0..0] (if not sectorized), or [0..Nl-1] (if sectorized)


struct grid_point_prototype{//info for each grid point
  int type;                      //0-out, 1-active, 2-connected to active through the boundary condition 
  int sl;                        //type 1: sequential label of point (init:fill_arrays()) [0..Nw-1]
                                 //type 2: [0..Nr-1] (filled while building the Hamiltonian) - each time a new one is encountered
  int type_graphene;		 //0-out, 1-inside, 2-boundary, 3-central
  int sl_g[2];	 		//sequential labels (into Nw_plag) of the up and down sigma components (-1 if type_g=0/3), the same if type_g=2
  int plaq;			 //sequential label for the corresponding plaquette (-1 if not completely inside)
  int b_sum_i, b_sum_j;          //sum of the vectors towards all neareast points outside the grid
  
  coordinate_prototype active_coordinate; //active coordinates (sequential) (init:fill_arrays())
  int sector;                    //which sector it belongs to: -1:not active, 0:conductor, 1,...,Nl:lead, 
                                 //(preliminary init[0/-1]:fill_arrays() final init:sectorize())
  
};
grid_point_prototype *grid_point;   //all the cros overs of the net lines indexed by [0..Nx*(Ny+1)+0..Ny]

struct active_point_prototype {//info on sequential point
  coordinate_prototype coordinate; //(init:fill_arrays())
  double volume;               //volume corresponding to the cell (1 for dirichlet, <0,1> for neumann - implemented only for 1D!!!
};
active_point_prototype *active_point;//indexed by [0..Nw-1]

int *sl_btp, N_btp;             //labels of active points which reach inside through boundary conditions, their number
content bc_coefficient;         //boundary condition coefficients (works only for 1D), otherwise each active point has its own
  
struct relative_point_prototype {//is expressed through c0 + c1*sl0 + sum_i=1..length coef[i]*sl[i]
  static int max_length;        //maximal length of the matrix
  coordinate_prototype coordinate;  //grid coordinate of this point
  //first alternative - ('recursive solution')
  bool solved;                 //is already expressed through active points only
  int sl0;                     //special active point label
  content c0, c1;
  int length;
  int *sl;
  content *coef;
                                //inserts a new point - if not there yet, increase the length; 
  void update(int sl_new, content coef_new);

  //alternative extrapolation/interpolation ('global parameters'):
  prec ri, re;                  //radii of the interpolated and extrapolated point
  content* fe, *fi;             //fe = exp[ik(re-ri)] fi (for each set)
  int length2;                  //fi = sum_{i=0..length2-1} f[sl2[i]]*coef2[i]
  int *sl2;
  content *coef2;
};
relative_point_prototype *relative_point; //[0..Nr-1] //Nr is number of active relative points
int Nr, Nr_max;                           //number / maximal number of active relative points (100)
void express_relative_point(int sl);//expresses this relative point through active points and the special point, calls itself recursively; input = sequential label into relative points array
content *relative_point_coefs;//array to store coefs
int* relative_point_labels;
void interpolate_relative_point(int sl);//fill in length and arrays sl2, coef2 in relative point sl
void extrapolate_relative_points(content *input); //refresh fe's and fi's using the input and boundary condition
char relative_point_method;//'g'/'r' using 'global parameters' or 'recursive solving' for relative points?
  
struct local_vector_prototype { //operator acting locally = linear combination of nearby points
  static int max_length;        //maximal length of the matrix
  int length, length_r, length_s;//number of nearby active     inside / relative  / special points 
  int *sl, *sl_r, *sl_s;         //sequential labels of nearby inside / relative  / special points
  content *coef, *coef_r, *coef_s;//coefficients of the linear combination
                                //inserts a new point - if not there yet, increase the length; 
                                //type: 1-inside, 2-relative, 3-special
  void update(int sl_new, content coef_new, int type);
};
local_vector_prototype **local_vector;//indexed [set_from*sets+set_to][0..Nw-1(active point label)]
int *local_vector_labels;       //array to store labels
content *local_vector_coefs;    //array to store coefs

//ARRAYS
content* others;		//buffer for other operators
				//with size # of points inside
				// *(NOT times sets!!!)
content* temp_mv;		//buffer for mean value operations
				// sets*(size # of points inside)

//CONSTANTS

prec Lx0,Ly0;		//left down corner (point 0,0)
prec Lx,Ly;		//dimension of the region
int Nx,Nx1,Ny,Ny1;	//# of xlines, +1, y lines, +1
int Nxb,Nxb1,Nyb,Nyb1;  //border of the inside of the region
			//Nxb=Nx-b(last line of the active region) , Nxb1=Nx-(b-1) (first line of the upper border) 
prec hx,hy;		//grid distance in x,y direction hx=Lx/(Nx-1)

			//FFT info:
int tab_nx,tab_ny;	//number of xlines, ylines in the table for FFT
prec kx,ky;		//grid distance in reciprocal space

function_inside* func;//'inside' function - definition of the region
bool inside_with_neighbors(int i, int j);//calls the previous - inside only if also the neighbors are inside (smoothens the boundary, no single standing points

operators_precision act_prec; //actual precision used
bool pbc_x;   //periodic boundary conditions?
bool pbc_y;
bool neumann, dirichlet, periodic;//boundary conditions

			//diferential operators:
prec* p_op_mat[all_op][all_prec];
			//difference operators defined by matrices


struct {int x,y;} op_mat_dim[all_op][all_prec];
			//dimension of the corresponding matrix
			//(-x,x)x(-y,y)

content **dif_op_vals[all_op];
			//value of difference operator stored
bool *dif_op_filled[all_op];
			//exist it = memory allocated and filled by values?
bool remember[all_op];
			//should I remember this operator when computing?
bool operators_sym1, operators_sym2;
			//1: will be operators symetrized? - means will function 'integral' be used?
      //2: will the functional coefficients symmetrized? - means will the operators forced to be hermitian?

struct multiple {
        struct grid_point_prototype upsl;	//label of the neighbour
        content val;}; 		//value for the multiplication

multiple*** multiples;
			//for each inside point a table with neighbour labels
			//and values for multiplication - this all for 4 part of the
			//hamiltonian

int b;		        //boundary - # of lines on boundary

int Np;			//# of points in the net: Np=(Nx+1) * (Ny+1)
int Nl;     //number of sectors (number of leads plus one for the conductor)
bool sectorized;//sector labels initialized?

public:
int Nw;			//# of points inside the region
int sets; //number of copies of the grid points in the wavefunction (2S+1, where S is the spin)
int lines;              //# of active xlines = last-first+1
int first_active_xline, last_active_xline;  //first and last nonempty xline

//FUNCTION PROTOTYPES

//CONSTRUCTOR
//Nx - # of xlines, Lx0 left side(lower bound) x-coordinate, Lx-x-length of the region, 
//pbc_x- periodic boundary conition in along the x - coordinate
//boundary - # of boundary lines (>= maximal stretch of matrix operators)
//f - function defining the region
region(int Nx, prec Lx0, prec Lx, char bc_x, int Ny, prec Ly0, prec Ly, char bc_y, function_inside f, int sets);
//DESTRUCTOR
~region();

//sets the remember flag for this operator
void remember_set(operators op,bool val);

//sets the actual precision to be used
void act_prec_set(operators_precision val);

//sets the symetrization on /works only with operate_tot_init/
void symmetrize_ops(bool s1, bool s2=true);

//fill in arrays according to sectors (conductor->0, leads->l+1, outside->-1)
//Nl = number of leads
void sectorize(int (inside)(prec x, prec y), int Nl);

//computed operator over the field is set(added) onto field of output with coeff. koef (and possibly with functional koeficient)
//if input==0 input is reset to point on a vector of 1
void operate(content *input, content* output,int which_field_x, int which_field_y, operators op, content koeficient, content func_koef(prec,prec) ,setting_or_adding);

//clears the table multiples
void operate_tot_init();
  
//sets the overall shift [internal units - energy]
void operate_tot_init_shift(double shift);

//convert the local vector into boundary point into local vector into inside points plus special point
//update the special points array
void operate_tot_init_boundary();

//update the boundary condition - works only for 1D as of now - the outgoing (real for CF states) wavevector in internal units
void operate_tot_init_boundary(content k);

//adds values into the table multiples according to operator
void operate_tot_init(int w_f_x, int w_f_y, operators op, content koef, content f(prec,prec),content integral(prec, prec, prec, prec));

private:
double overall_shift;           //a constant shift used in operate_tot: output + = shift * input
public:
    
//provides the multiplication output=local_vectors \otimes input
void operate_tot(content *input, content* output);

//computes temp_mv=(op)|\Psi> 
void op_ket(content *input,content *output, mean_value_of_operator op,setting_or_adding soa);

//computes mean value of an operator = sum <\Psi|op\Psi>
content mean_value(content *input, mean_value_of_operator op, int type=0);

//computes scalar product 
//type = 0 : <input1|input2>
//type = 1 : <input1^*|input2>
content braket(content *input1, int which_field_x, content *input2, int which_field_y, int type=0);

//multiplies N vectors defined by input, possibly inverse and / or conjugate each of them when multiplying
content braketN(int N, content **input, bool* inverse, bool* conjugate);


//normalizes the vectors to integral with cell weighted by their relative volume (if neumann, boundary cells are smaller - so far works only for 1D); fits also the phase
void renormalize_to_volume(content *input, int w_f_x);

//compute the inverse participation ratio and the center of mass (both in internal length units)
void ipr_and_com(content *input, double& ipr, double& x0);

/*computes scalar product: sets/adds(acc. to soa) output(+)= input1^* * input 2 */
void braket_vec(content *input1, content *input2,content *output, setting_or_adding soa);

/*computes scalar product: sets/adds(acc. to soa) output(+)= input1^* * input 2 */
void braket_vec(double *input1, content *input2,content *output, setting_or_adding soa);

//returns || f - R_which(f)|| / 2 ||f|| for
//R_{1,2,3}= {inv x, inv y, inv}
//prec parity(content *input, int which);

//clear all (flags of) computed field - new data arrived
//and deallocate memory
void new_values();

//output values into a table with coordinates (i,-j) (beginning in the left upper corner)
//what : 1 real part 2 imag part 4 abs
//if what is negative, a column with coordinates is printed:  (x,y,Re Psi, Im Psi, abs(Psi) \n)
void draw(content* input, int w_f, int what, content zero, int only_inside, FILE* file, int where_to);

//draws function f using the above
void draw_function(content (f)(prec,prec), int what, int only_inside, FILE* file, int where);

//converts a linear vector(input, w_f) of points inside region into a 2d (nx times ny) table (output)
//the table should be Ny+2 times Nx+2 large to save also the boundary zeros
void region2table(content *input, int w_f, content *output, int nx, int ny);
//converts 2d (nx times ny) table (input) into a linear vector representing the points inside the region (output)
//the table should be at least Ny+2 times Nx+2 
void table2region(content *input, content *output, int w_f,int nx, int ny);

//for coordinates i=0,...,tab_nx-1, j=0,...,tab_ny-1 gives proper
//vector (omega_x,omega_y) which is from the region (-some value,+some value)
void give_omega(int i, int j, prec& omega_x, prec& omega_y);


//walk through the created net
void ant(function_potential, double l2nm, double e2meV);

//test diferential operator
//output[10]: norm+max of function, norm+max of derivation, norm+max of measured der.
//	  norm+max of (der-mder), norm+max of (der-mder)/der
void test(operators op, content func(prec,prec), content der(prec, prec), prec *output);

int invertconseq(int n, operators op);
void invert(int i, int j, int& ii, int &ij, operators op);

//interpolate function using second order taylor at real space point
//the function and derivatives are in f: f[0], f[1],..., f[5] = f, fx, fy, fxy, fxx, fyy
content value_at(prec x, prec y, content** f);
//compute the derivatives for the above routine, allocate memory
void value_at_init(content** f, content* input);
//release memory
void value_at_release(content **f);


void give_par(prec &h_x, prec &h_y) {h_x=hx;h_y=hy;}
void give_par_2(prec &L_x0, prec &L_y0) {L_x0=Lx0;L_y0=Ly0;}
void give_par_fft(prec &k_x, prec& k_y) {k_x=kx;k_y=ky;}
void give_par_fft_2(int &n_x, int& n_y) {n_x=tab_nx;n_y=tab_ny;}
double true_length() {if (relative_point_method=='r') return((Nw-1)*hx); else return((Nw+0.5)*hx);}
void give_par3(content &bc, prec& shift) {bc=bc_coefficient; shift=overall_shift;}

//converts consecutive label into real space coordinates (in dimless units)
//void sl2xy(int label, double& x, double &y);
//returns real space coordinates of the maximum of xi_bra*xi_ket (two component spinors!)
void maxofbracket(content* bra, content* ket, double& x, double &y);

private:

//converts coordinate to sequential label - can account for periodic boundary conditions
//inline int c2sl(coordinate);
inline int c2sl(int i, int j);
inline int ac2sl(int i, int j);
inline void sl2c(int , int&, int &);


//auxiliary functions for computing lz operator
static content mv_lz_1(prec x, prec y);
static content mv_lz_2(prec x, prec y);
static content mv_r2(prec x, prec y);
static content mv_r4(prec x, prec y);
static content mv_x(prec x, prec y);
static content mv_y(prec x, prec y);

//converts relative coordinate into index in multiple table
inline int relative2multiple(int, int);

//converts relative coordinate into the index of a matrix <nl,nh>x<ml,mh>
inline int relative2matrix(int, int,int,int,int,int);


//does matrix multiplication (linear operator in a point defined by matrix koeficients
//matrix is defined as (nl,..,nh)x(ml,..,mh) with indexes being relative coordinates to the actual point
void dn_point_operation(content * input, content *output, int nl, int nh,int ml, int mh, prec* koeficients);

//similarily as previous, but: it does not write into output, but actualises the table multiples and if operators_sym1 is on, it adds *exp(integral) to each value in multiple
void dn_point_operation_tot(int w_f_x, int w_f_y, int nl, int nh,int ml, int mh, prec* koeficients, content f(prec,prec), content kf, content integral(prec, prec, prec, prec));


//returns distance in real (natural) units of length
inline prec xd(int i);
inline prec yd(int j);

//fill matrix koeficients for differential operators
void inicialize_diferential_operators();

//auxillary
void listop(operators op, operators_precision l); 
void inicialize_diferential_operators_test();

//allocate memory for all arrays
void inicialize();
  void inicialize1();
  void inicialize2();
  void inicialize3();
  void inicialize4();
  void inicialize5();
  void inicialize6();
  void inicialize7();
  void inicialize8();
  
//deallocate memory
void clean();

//set all values in multiples to 0
void clean_multiples();

//label points with sequential labels, account for PBC by defining equivalent couples such, that b and Nb lines are IDENTICAL!!!
// i\in<0,b-1> : [i, Nxb-b+i], [Nxb+i,b+i] for x direction
void fill_arrays();

public:

//!graphene
sparse_matrix_prototype * M_sparse;//energy stored as a sparse matrix
sparse_matrix_prototype * HG_sparse;//Hamiltonian stored as a sparse matrix

content*  M, *HG, *LU;  //Dirac equation matrixes: energy(RHS), Hamiltonian(LHS) and LU decomposition
  //matrix M converted to a local vector by lines: each line is compacted listing
int* M_length;          //number of nonzero components (usually four, 8 for points around center)
int M_length_max;       //maximal length - 9 times other dof than sigma
content* M_locvecval;   //their values
int* M_locvecind;       //their indexes (from 0 to g_length)
int Nw_g;		    //number of complex components of the wavefunction, somewhat smaller than 2 * Nw
int Nw_plaq;		//number of plaquettes on which Dirac equation is considered, must equal to Nw_g/2
int g_length;   //length of the input/output vector and the dimension of the matrices M, LU, HG
                //g_length = Nw_g * (sets/2)
int* IPIV;      //auxiliary for LU decomposition
  
//apply the operator for the ARPACK: [ H^-1 * M ] 
void inverseLUtimesM(content *x, content* yd); //y = LU^-1 ( M ( x ) ) 
void inverseLUtimesM_with_checking(content *x, content* yd);

			//converts the sigma subspace wavefunction into grid wavefunction inserting values using boundary conditions; the first argument has length g_length, the other two have lengths Nw; offset selects the component other that sigma
void psip2psig(content *psip, int offset, content* psig_up, content* psig_down);
			//update graphene Hamiltonian by applying |sigi><sigj| op (unity/k+/k-) to i-j-th subsector of the Hamiltonian, which  is a matrix with linear dimension of g_length = N * Nw_g and N=sets/2 (N will be 2 if spin and sigma dof are considered)
      //with possible constant, and functional, and peierls spatially dependent multiplication
void update_g_matrix(int sigi, int sigj, const char* op, int ioffset, int joffset, content kf, content f(prec,prec), content integral(prec, prec, prec, prec), content *H, sparse_matrix_prototype * H_sparse);
void convertM2locvec();//convert the explixitly constructed M matrix into a compressed form
void graphene_direct_check(int eig, content* vecs, content *vals);
void deallocate_graphene();
  //auxiliary routine defining boundary outward normal derivative azimuthal angle phi
  //returns ii * exp(ii phi)
content graphene_boundary_aux(int i, int j, int sig);


//print out constants
void constants(prec l2nm);

//converts sequential point label [0,...,Nw-1] to length coordinates
void s2coordinate(int sl, prec& x, prec &y);
//converts length coordinates to nearest sequential point
bool coordinate2s(prec x, prec y, int& sl);


//active2grid
void active2grid(int ai, int aj, int &i, int &j);
//grid to active (all sectors) coordinates; returns -1 if not inside
bool grid2active(int i, int j, int &ai, int &aj);
//converts the sector into grid coordinates
void sec2grid(int sec, int si, int sj, int &i, int &j);
//active2sector
void active2sec(int ai, int aj, int &sec, int &si, int &sj);
//converts the grid into sector coordinates
bool grid2sec(int i, int j, int& sec, int& si, int& sj);
//active2sector
void sec2active(int sec, int si, int sj, int& ai, int& aj);
//number of points in sector on the grid xline
int grid2secdim(int sec, int i);
//xlines range for this sector (grid coordinates)
void sec2gridlimits(int sec, int &i, int &j);
//is the point nearby sector sec2? input:sec coordinates
bool secnearbysec(int sec, int si, int sj, int sec2);
//dimension of the sector xline
int sec2dim(int sec, int i);
//returns grid coodinates of the nearest point
//if the grid point is outside the region, returns false, others zeros
bool nearest_grid_point(prec x, prec y, int &i, int &j);

//hamiltonian matrix element, input are two points with sector coordinates 
content Hsecsec(int sec1, int i, int iy, int is, int sec2, int j, int jy, int js);

//velocity of a wave e^(i k x) along x direction
prec wave_velocity(prec k);

};


inline double region::xd(int i)
{
  return(xline[i].fb);
}


inline double region::yd(int j)
{
  return(yline[j].fb);
}

inline int region::ac2sl(int i, int j)
{
  int value=grid_point[i*Ny1+j].sl;
  if (value<0) return(-value);
  else return(value);
}

inline int region::c2sl(int i, int j)
{
  if (i<0 && j<0 && i>Nx && j>Ny) {printf("c2sl:out of range\n");exit(1);}
  return(grid_point[i*Ny1+j].sl);
}

inline void region::sl2c(int label, int& i, int &j)
{
  if (label<0 || label>=Nw) {printf("sl2c:label %i outside range (%i-%i)\n",label,0,Nw-1);exit(1);}
  for (i=label/(Ny+1-2*b)+b;i<Nxb1;i++)       //initial guess
    if (xline[i].np>0 && xline[i].fpsl+xline[i].np>label) break; //I have the right line
    j=label-xline[i].fpsl+xline[i].fpc;//initial guess
    while (c2sl(i,j)!=label) {
      j++;  //right point
      if (j>Ny) {printf("searching for label %i failed in sl2c\n", label);exit(1);}
    }
}

inline int region::relative2multiple(int i, int j)
{
  if (i<-b || i>b || j<-b || j>b) return(-1);
  if (i==0 && j==0) return(0);
  if (i==0) return(j+b+1*(j<0));
  if (j==0) return(i+3*b+1*(i<0));
  if (j==i) return(j+5*b+1*(j<0));
  if (j==-i) return(j+7*b+1*(j<0));
  return(-1);
}


inline int region::relative2matrix(int nl, int nh,int ml,int mh, int k,int l)
{
  if ((nl>k) || (nh<k) || (ml>l) || (mh<l)) return(-1);
  return((k-nl)+(l-ml)*(nh-nl+1));
}



//is the point nearby sector sec2? input:sec coordinates
inline bool region::secnearbysec(int sec, int si, int sj, int sec2)
{
  int ai,aj,i,j;
  sec2active(sec,si,sj,ai,aj);
  sec2grid(sec,si,sj,i,j);
  return(xline[i].active_line_point[aj].nearby_sector[sec2]);
}

inline void region::s2coordinate(int sl, prec& x, prec &y)
{
  x=xd(active_point[sl].coordinate.i);
  y=yd(active_point[sl].coordinate.j);
}

inline bool region::coordinate2s(prec x, prec y, int& sl)
{
  sl=-1;
  int i,j;
  bool in=nearest_grid_point(x,y,i,j);
  if (!in) return(false);
  if (grid_point[i*Ny1+j].type!=1) return(false);
  sl=grid_point[i*Ny1+j].sl;
  return(true);
}

//active2grid
inline void region::active2grid(int ai, int aj, int &i, int &j)
{
  i=ai+first_active_xline;
  j=xline[i].active_line_point[aj].j;
}

//grid to active (all sectors) coordinates; returns -1 if not inside
inline bool region::grid2active(int i, int j, int &ai, int &aj)
{
  ai=grid_point[i*Ny1+j].active_coordinate.i;
  aj=grid_point[i*Ny1+j].active_coordinate.j;
  if (ai==-1 || aj==-1) return(false);
  return(true);
}

//converts the sector into grid coordinates
inline void region::sec2grid(int sec, int si, int sj, int &i, int &j)
{
#ifdef DEBUG
  if (sec<0 || sec>Nl) {printf("sec2grid:wrong sector id=%i, with Nl=%i\n",sec,Nl);exit(1);}
#endif
  i=si+sector[sec].imin;
  j=xline[i].active_line_point[sj+sector[sec].offset[si]].j;
}

//active2sector
inline void region::active2sec(int ai, int aj, int &sec, int &si, int &sj)
{
  int i,j;
  active2grid(ai,aj,i,j);
  sec=grid_point[i*Ny1+j].sector;
#ifdef DEBUG
  if (sec==-1) {printf("active2sec:wrong sector id=%i, of the point grid:(%i,%i), active:(%i,%i)\n",sec,i,j,ai,aj);exit(1);}
#endif
  si=i-sector[sec].imin;
  sj=aj-sector[sec].offset[si];
}

//converts the sector into grid coordinates
inline bool region::grid2sec(int i, int j, int& sec, int& si, int& sj)
{
  int ai, aj;
  if (!grid2active(i,j,ai,aj)) {sec=si=sj=-1;return(false);}
  active2sec(ai,aj,sec,si,sj);
  return(true);
}

//sector2active 
inline void region::sec2active(int sec, int si, int sj, int& ai, int& aj)
{
  int i,j;
  sec2grid(sec,si,sj,i,j);
  grid2active(i,j,ai,aj);
}

//grid to active (all sectors) coordinates; returns -1 if not inside
inline int region::grid2secdim(int sec, int i)
{
#ifdef DEBUG
  if (i<0 || i>=Nx1) {printf ("grid2dim: xline %i out of range (%i-%i)\n",i,0,Nx1);exit(1);}
#endif
  return(sector[sec].dim[i-sector[sec].imin]);
}

//dimension of the sector xline
inline int region::sec2dim(int sec, int i)
{
  return(sector[sec].dim[i]);
}

//xlines range for this sector (grid coordinates)
inline void region::sec2gridlimits(int sec, int &i, int &j)
{
  i=sector[sec].imin;
  j=sector[sec].imax;
}

//!THIS NEEDS TO BE REWRITTEN TO BE FASTER - CONSTRUCT MATRIX COEFFICIENTS (multiples) WITHOUT ANY COMPRESSION (that is relative2multiple function)
//!multiples are stored using w_f_x (input) w_f_y (output) for the spinor subspace, which is opposite to the convention
//!in the (output | H | input) 
inline complex<double> region::Hsecsec(int sec1, int i, int iy, int is, int sec2, int j, int jy, int js) 
//inline complex<double> region::Hsecsec(int sec2, int j, int jy, int js, int sec1, int i, int iy, int is) 
{
  int gi, gj, giy, gjy;
  sec2grid(sec1,i,iy,gi,giy);
  sec2grid(sec2,j,jy,gj,gjy);
  int ir=gj-gi;
  int jr=gjy-giy;
  if (abs(ir)>b || abs(jr)>b) return(0);
  int label_multiple=relative2multiple(ir,jr);
  if (label_multiple==-1) return(0);
  int label=c2sl(gi,giy);
  //int label2=c2sl(gj,gjy);
  //complex<double> val=0;
  //local_vector_prototype & lm=local_vector[js*sets+is][label-1];
  //int l,l2;
  //for (l=0;l<lm.length;l++) if (lm.sl[l]==label2) {val=lm.coef[l];break;}//!here is <-> js are swapped
  complex<double> val=multiples[js*sets+is][label][label_multiple].val;//!here is <-> js are swapped
#ifdef DEBUG
  int label2=c2sl(gj,gjy);
  int label_multiple2=relative2multiple(-ir,-jr);
  if (label_multiple2*label_multiple<0) {
    printf("multiple labels fckdp(\n");
    printf("ir:%i, jr:%i, ml:%i, ml2:%i\n",ir,jr,label_multiple,label_multiple2);
    exit(1);
  }
  if (label_multiple2==-1) return(0);
  //complex<double> val2=0;
  //local_vector_prototype & lm2=local_vector[is*sets+js][label2-1];
  //for (l2=0;l2<lm2.length;l2++) if (lm2.sl[l2]==label1) {val2=lm2.coef[l2];break;}
  complex<double> val2=multiples[is*sets+js][label2][label_multiple2].val;//!here also
  if (abs(val-conj(val2))>=1e-10) {
    printf("they are not conjugated: %e%+ei vs %e%+ei\n",val.real(),val.imag(),val2.real(),val2.imag());
    printf("\tgi=%i, giy=%i, is=%i is (label %i) / gj=%i, gjy=%i, js=%i (label %i)\n",gi,giy,is,label,gj,gjy,js,label2);
    printf("\tlabel_multiple,label_multiple2 %i,%i\n",label_multiple,label_multiple2);
    printf("\tfirst (label=%i) points to neighbour label %i, second (label=%i) points to neighbour label %i\n",label,multiples[js*sets+is][label][label_multiple].upsl.sl,label2,multiples[js*sets+is][label2][label_multiple2].upsl.sl);
  }
#endif
  return(val);
}






#endif

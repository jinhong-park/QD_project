#include "main.h"

//#define HALF_PUMPED
//#define SYMMETRIC
//#define OPEN_LEFT
#define OPEN_RIGHT

#define LASER 2
#define SORTING 0
#define LASER_CHECK
//#define LASER_CHECK_EPS

class simple_sorting_prototype {
  int N;                              //length of the array
  int entries;                        //max number of entries (past) = order of polynomial approximation
  int * valid_entries;
  content * value;                    //value [ (i=0,..N-1) * entries + past] past = 0, ... entries-1
  content * estimate;                 //estimates at the new coordinate
  double * coor;                      //coordinates
  int *index;                         //auxiliary for sorting
  content *vecR;		      //N by N matrix of right eigenvectors (indexed by unsorted labels) so that new arrived eigenvectorL is compared to vecR[sorted2unsorted[i]]
  content *O;			      //N by N matrix of overlaps
public:
  int *position;                      //position [old] = new (sorted2unsorted)
private:
  double *error1, *error2, *error3;    //relative error of the estimate, relative error of the estimate vs assigned value, overlap difference from 1 
  double *coor_at_zero;               //estimated coordinates where the imaginary part of the value reaches zero
  double *value_at_zero;              //values at zero
  int zeros;                          //number of zero crossings (-1 means uninitialized)
  bool *crosses_zero;                 //whether this one crosses zero
  
public:
  simple_sorting_prototype(int N, int entries); //constructor
  void copy(simple_sorting_prototype* orig);    //a copy of data of an existing object
  ~simple_sorting_prototype();        //destructor
  bool zero_was_crossed(content& left, content &right); //did some level(s) cross zero? if yes, l/r equal Re(largest) at boundaries
  
  //updates the history. If special!=-1 the index special is used to discard one set from the history
  void update(content* new_value, double new_coor, content* der, int special, bool eigvals);
  
  void sort_mother(int methods, content *vecL, content *vecR, content* new_value, double new_coor, double* errors=0, content* der=0, int special=-1); //sort the non-hermitian matrix eigenvectors according to their eigenvalues (!methods bit 0) or overlaps (methods bit 0) possibly with a given special index. Evaluate error[0-2] as 0:error of the estimate, 1:the difference of the estimate to the identified new value, 2:overlap difference from 1. Evaluate eigenvalues derivatives if der!=0. methods bit 1 = skip step 1 of the algorithm (try), methods bit 2 = skip step 2 (accept)
  content give_value(int i, int offset=0, double* coor=0); //return the sorted value #i from the actual set (backwards by offset, also the corresponding coordinate
  int zeros_crossed(double cl=0, double cr=0);//returns how many of sorted states cross zero between cl (coor[0] by default) and cr (coor[1])
  bool zero_estimate(int i, double& coor, double& value);//returns position and real part of the i-th state (returns whether it crosses zero)
  bool extrapolate(int i, double &coordinate, double value, int from_last=-1);//extrapolates valid values to value and returns the corresponding coordinate; returns success
  void debug_info();//flush out the state
};

struct laser_parameter_prototype {

  //CF basis are psi(x) : = psi^Num (x) / units.length^(dim/2)
  //electric field Psi (x) = c_m psi_m^Num(x) => c_m are in Volt / meter
  double gamma_perp;              //polarization decay (decoherence)        [frequency] 10^12-13 1/s Optics Express 2008
  double gamma_par;               //population inversion decay (relaxation) [frequency] 10^8-9 1/s Optics Express 2008
  double g;                       //charge dipole element                   [charge length] 25 Debye Spectrochimica Acta 22, 2121 (1966)
  double n2;                      //dielectric constant                     [1]
  double D0,D,Dinc;               //population inversion (~pump strength)   [1/volume], [1], increase (relative,~ C3)
  double omegaa;                  //active atom frequency                   [frequency] 30 gamma_perp Science 2008
                                                                            //630 rhodamin640 dye in methanol
  double na;                      //active atom density                     [1/volume] 0.0025 M Nature 368, 436
  
  //global parameters
  region* r;                      //the grid
  content *aux;                   //workspace of dimension r->Nw
  int dim;                        //r->Nw
  int Nbasis,Ndiag;               //number of CF states to use as basis, to ask for in diagonalization
  content* d;                     //pump profile D(x)/D
  content* nx;                    //refractive index envelope (n(x) = sqrt(LP.n2) * LP.nx[x])
  content* n2x;                   //refractive index envelope (n^2(x) = LP.n2 * LP.n2x[x])
  double *total_profile;          //total profile = 1 + sum of weighted profiles = 1+sum_mu Gamma(mu) |Psi_mu(x)|^2
  double total_modal_integral;    //sum of modal integrals of lasing modes
  double length;                  //length of the active region - depends on method g/r
  arpack_zn* arpack;              //arpack
  double shift_H;                 //shift in the Hamiltonian (- k*k )
  bool FFA;                       //fixed frequency approximation - CF basis are freq. independent
  bool TMsorting;		  //1 - according eigenvalues, 0 - according eigenvectors
  bool with_shifts;               //if shift of the spectrum at -k*k should be applied
  int perturb_fail;		  //if perturbative starts to fail, do not do it for some number of steps
  //class CF_basis_prototype* CF_basis_FFA;//the only CF basis is stored here
  class CF_pool_prototype* CF_pool;//the CF basis pool
  
  content ** inputs;               //auxilliary for region2d -> braket(N,...)
  bool * inverse, * conjugate;
  void clear_braket_flags() {for (int i=0;i<10;i++) conjugate[i]=inverse[i]=false;}//set all inverse and conjugate to false

  
  inline double x(double k) {return(k*units.wavevector*units.c-omegaa)/gamma_perp;}  ////conversion of frequency/wavevector to dimless units
  inline content x(content k) {return(k*units.wavevector*units.c-omegaa)/gamma_perp;}  ////conversion of frequency/wavevector to dimless units
  double Gamma(double k) {return(1/(1+x(k)*x(k)));}                         //lorenzian frequency/wavevector filter   [1]
  double A() {return(4*g*g/(units.hbar*units.hbar*gamma_par*gamma_perp));}  //electric field^2 "units"   [meter^2/Volt^2]
  double B() {return(D0*g*g/(units.eps*units.hbar*gamma_perp));}            //population inversion "units"    [1]
  content Lambda(double k, content kp2) {return( -k*k/(k*k-kp2) / (x(k)+content(0,1)) );} //threshold matrix diagonal [1]
  void init();                    //constructor; initialize parameters
  void deallocate();              //destructor
  
  int strategy;			  //update strategy
};

class CF_basis_prototype {
public:
  int N;          //number of terms
  int dim;        //dimension of the grid (single term length)
  double k;       //boundary condition
  content *km2;    //complex eigenvalues (km^2)
  content **psi;  //complex eigenfunctions
  int *pos;       //positions of sorted eigenstates
  
  CF_basis_prototype(int dim, int N, bool allocate=false); //constructor - eigenvecs/eigenvals allocated *optionally*
  ~CF_basis_prototype();              //destructor - deallocate *all* space
  int load(ret_from_arpack_zn& vysl, double k, double k1=0, double k2=0);//copy the states and possibly sort ascending according to "distance" to the line (k1, k2)
  void normalize();//normalize to "volume"
  void filter(content* f, content *fcf);//projects function f on the span of the CF basis
  void plot(FILE * file);//plot CFs into file
  content correlation1(int i, int j, int dx);//returns the correlation of basis functions i and j shifted by dx grid points  
  content correlation2(int i, int j, CF_basis_prototype* CF2);//returns the correlation of basis functions i from this and j from a different CF basis  
  content correlation3(int i, int j, CF_basis_prototype* CF2);//returns the normalized scalar product i from this and j from a different CF basis  
  void compare_with(CF_basis_prototype* CF2);//compare how close are two CF bases and sort "this" according to the overlaps with CF2 in the same order
  void spit_km2(FILE* file, double L2);//list eigenvalues [ sqrt(km2) ] sorted according to the real part
};

class CF_pool_prototype {
public:
  int items, maxitems;          //number of stored CF_bases and the limit
  CF_basis_prototype** CF_basis;//stored CF_bases
  double *k, k1, k2;            //corresponding frequencies of the CF_bases and covered interval (copy of laser_prototype k1, k2)
  int *undel;                   //whether the table can be deleted if limit is reached: -1 never, 0 yes,1 after the pump is enlarged
  
  double grain;                 //if an existing k is within (gamma_perp/grain) from the asked for, the basis is accepted as the valid approximation. Otherwise the basis is created. Grain=0 is FFA, grain<0 means always a new basis is created
  
  int closest(double k0);       //return the pointer to the existing basis with closest frequency to k0
  CF_basis_prototype* give_CF_basis(double k0, int undeletable=0);    //returns the closest if within tolerance, otherwise construct one, including the allocation of the memory and inserting the new basis into the ordered stock. If limit is reached, one basis is deleted
  CF_pool_prototype(double grain, double maxmem);
  ~CF_pool_prototype();
  void refresh_undeletable() {for (int i=0;i<items;i++) if (undel[i]==1) undel[i]=0;};
};


struct TM_prototype {
  int N;                             //linear dimension of the matrix
  content *M, *Tau;                  //matrix elements
  content *eig, *eig_der;            //eigenvalues and their derivative wrt k
  content *vecL, *vecR;              //l- and r-eigenvectors (optional)
  int* pos;                          //positions during sorting

  lapack_zgee_prototype* lapack;     //class for lapack diagonalization
  CF_basis_prototype* CF_basis;      //corresponding CF_basis
  
  TM_prototype(int N);                //allocate space
  ~TM_prototype();                    //deallocate space
  void copy(TM_prototype& TM_source, int* index=0);  //copy the matrix, with index table
  void diagonalize(bool skipvecs=true);  //diagonalize, optionally with eigenvectors; sorting is ascending according to distance to 1/D; if D==0 it is descending according to the real part
  void construct(CF_basis_prototype* CF_basis, bool Tau_only=false);  //construct TM (Tau_only is not used currently, the whole matrix is always constructed)
  void update_from_Tau(double k);     //construct TM as L(k) times Tau, where Tau is already computed
};

struct flag_prototype {
  content braket, diff;			//braket and diff of the previous(input)/actual(output) step
  int last_step;			//step of the last update (begins with -1)
  int ban_step;				//step when perturbative update was banned (no ban=-1)
  int bans, updates;			//how many bans (begins with zero) and updates
  bool converged, iterative, freq_reset;//if converged, if iterative is to be / was taken, if frequency is to be reset
  bool fixed;				//single mode fixed until convergence
  void reset();				//set all flags/counters at the iteration start
};

class lasing_mode_prototype {
public:
  CF_basis_prototype* CF_basis;     //corresponding CF basis
  double k, k_last_update;          //lasing frequency and the last update in freqency reset
  TM_prototype* TM;                 //threshold matrix for the actual iteration
  content* A, *A_save;              //amplitudes of basis states
  content *amplitude;               //Psi(x) = A_m psi_m(x)
  double* profile, *profile_save;   //mode profile Gamma(k) |Psi(x)|^2 = Gamma(k) |sum_m A_m psi_m(x)|^2 - spatial function 
  double modal_integral, modal_power;//modal power (in units of 2 epsilon_0 E_0^2 Omega_\mu)
  double norm;                       //sqrt(sum of |amplitude|^2)
  content *c, *d, *e;               //parameters for the perturbative step
  double Dth;                       //lasing threshold pump

  lasing_mode_prototype(int dim, int N, double k); //constructor - allocate space
  ~lasing_mode_prototype();         //destructor - deallocate space
  void init_guess(double lambda);   //set the initial amplitude computing e coefficients and TM matrix
  void update_profile();            //recompute profile, including the constant Gamma and amplitude Psi
  void update_modal_power();        //recompute modal power (total profile should be calculated first)
  double update_frequency();        //reset frequency such that the largest amplitude does not change phase upon iteration - does not work
                                    //solve for the amplitude (defined by amplitudes of CF basis states) perturbatively, considering first hmn CF_states (if zero, all of them)
  void perturbative_amplitude_update(content* A_new, int hmn=0);
  
  //evaluates coefficients c d e for the perturbative solution for amplitude, assumes TM was already diagonalized 
  void evaluate_cd(int hmn=0);
  void evaluate_e(int hmn=2);
  
  content iterative_step(content *A_new);//amplitude update using the direct map A_new = D TM . A
  content perturbative_step(content *A_new, int hmn=0);//amplitude update using the perturbative solution of direct map
  bool step(content& braket, content& diff, bool iterative_force=false, int Nlm=0); //take step - returns whether iterative (or perturbative) was taken, together with resulting change of amplitude vector and difference of TM eigenvalue from 1/D
  void step2(flag_prototype& flag, int step);
  void amplitude2A();              //recompute the amplitude expansion coefficients in (most probably new) CF_basis 
  void omega_check(int i, const char * where);//computes omega[i] = <aL[i] | A> and prints omega0 and omega0/sum_i omegai -1 into log with additional message "where"
};

class laser_prototype {
public:
  int Nlm;                            //number of lasing modes
  //int N;                              //how many states to consider in constant flux bases
  lasing_mode_prototype** lasing_mode; //lasing modes

  int N_real;                         //number of real eigenvalues
  TM_prototype** TM_real;             //threshold matrices for all real eigenvalues
  int* real2lm;                       //corresponding index into lasing modes (if -1, it is not lasing)
  int* TM_real_which;                 //index of the TM eigenstate crossing the zero axis

  double* TM_k;                       //real frequencies at which TM has real eigenvalues
  content* TM_eigreal;                //corresponding eigenvalues
  
  //CF_basis_prototype** CF_basis;      //CF_basis grid
  CF_pool_prototype* CF_pool;//CF basis pool
  double k1,k2,k_step;                //the covered frequency interval (in internal frequency units)
  double scan_step, scan_extent;      //the covered frequency interval (in units of gamma_perp)
  int terms;                          //number of "grid" points - this may be changed dynamically by redefining scan_extent;
  
  laser_prototype(int N);             //constructor 
  ~laser_prototype();                 //destructor
  void construct_grid();              //construct the grid and arpack, and allocate all space
  void destruct_grid();               //deallocate the previous
  void construct_Hamiltonian(const char* filename=0, bool save=false);//construct the Hamiltonian (load from/save to file?)
        //compute the CF in equally spaced intervals to cover the frequency window
  bool construct_CF_grid(double scan_extent_init, double scan_step, double basis_resolution, FILE* file=0);
  bool rescale_CF_grid(double scan_extent_new);//change the interval for TM scans (input in units of gamma_perp)
  void scan_CF(double R1, double R2, int terms, FILE* file);//auxiliary

  void update_total_profile();        //update the total profile
        //locate frequencies within the interval at which eigenvalues are real, possibly output to files
        //consider only the largest hmn values
  void scan_TM(int hmn, FILE* file1=0, FILE* file2=0);

  //from sorting with one set of values, identify limits within the CF grid where the i-th eigenvalue gets real
       //find the frequency where the i-th sorted value reaches zero by search using polynomial fitting
       //returns frequency and complex eigenvalue
  bool identify_root2(int i, simple_sorting_prototype * sorting, double &k, content& v, double k_last_update=0, int history=3);

         //correct for the mode frequency if it became too imaginary, returns residuum of largest eigenvalue
  content reset_mode_frequency(int mode);
        //returns the frequency such that RHS in iteration is real (does not work curently)
  double update_frequency(content* a_old, content* a_new, double k);

        //change the pump intensity within limits, returns true if the next lasing mode is below actual pump
        //Dnext and Dnextnext are estimates for the next and second next lasing mode thresholds
  bool update_pump(double &D, double Dnext, double Dnextnext);
        //Dnext is the pump at which the next lasing mode appears (and corresponding frequency)
        //Dnextnext is the second next threshold estimate. Returns TM_real index of the new state 
  int next_lasing_mode_pump(double &knext, double& Dnext, double& Dnextnext, FILE* file=0);

  int iterate_at_pump(flag_prototype* flag, int step);//make one iteration for update of active lasing modes
  bool check_replicas(bool del=false);//check for replicas, optionally delete them
  int select_one_mode();//find a very small mode
  void converge_one_mode(int mode); //converge just one mode (probably decaying)
  void plot_lasing_modes(FILE * file);//plot lasing modes

  void main_cycle();                  //main iteration cycle 
  
  class snapshot_prototype* snapshot;
};

class snapshot_prototype {
  int Nlm;
  double* ks, *Dths;
  content *amplitudes;
public:
  double D;
  snapshot_prototype();//constructor
  ~snapshot_prototype();//destructor
  void save(laser_prototype& laser);//save the snapshot at the actual pump
  prec load(laser_prototype& laser);//load snapshot, returns the pump
};


struct lasing_state_char_prototype {
  double* k;			//lasing frequency
  double* Dth;			//threshold
  double* ipr;			//inverse participation ratio (in nm)
  double* x0;			//center of mass (in nm)
  double* o[2];			//maximal overlap with some other state
  double* oi[2];		//index of state giving maximal overlap
  
  lasing_state_char_prototype(int N) {
    k=new double[8*N];
    Dth=k+N;
    ipr=Dth+N;
    x0=ipr+N;
    o[0]=x0+N;
    o[1]=o[0]+N;
    oi[0]=o[1]+N;
    oi[1]=oi[0]+N;
  }
  ~lasing_state_char_prototype() {delete k;}
};


class TLM_approx_prototype {

  laser_prototype *laser;
  int Nl, Nr;                    //number of lasing modes and real TM eigenvalues
  double *I;                     //initial values of intensities of lasing modes
  content *dI;                   //increments of intensities of lasing modes
  content* lambda;               //eigenvalues of TM
  content* xi;                   //xi coefficients
  
  double *Dth;                   //running pump, pump threshold for each state (-1 means the mode is not yet lasing)
  double *ipr, *x0, *o[2], *oi[2];
  int* real2lm, *lm2real;        //conversion between real (up to Nr) and lasing mode (up to Nl) indexes
  
  content* xi_inv;               //the inverse of xi - dimension Nl x Nl
  void inverse_xi(int i);        //calculates the inverse of chi for the set of lasing modes, possibly expanded by state i
  void intensities(double D);    //calculates intensities of the lasing modes by inverting the xi matrix
 
public:
  TLM_approx_prototype(laser_prototype& laser); //constructor
  ~TLM_approx_prototype();       //destructor
  void initialize();             //fills in parameters from the laser, constructs matrix xi
          //finds the next threshold for pump as minimum from all not-yet lasing modes;D is the initial/threshold I/O pump
          //returns the index of the corresponding state
  int next_pump_threshold(double& D); 
  int predict(lasing_state_char_prototype& lasing_state_char, FILE* file=0);//returns max number of lasing modes and their thresholds, possibly construct the intensity vs pump diagram into file
};


//!following are defined in laser.cpp
extern laser_parameter_prototype LP;                              //parameters
void CF_operate(int N, void *A, content *x, content *y);          //"Hamiltonian" for diagonalization by arpack 
content pump_envelope(prec x, prec y);                            //D(x) / D_0
content refractive_index_n(prec x, prec y);                       //n(x)^power


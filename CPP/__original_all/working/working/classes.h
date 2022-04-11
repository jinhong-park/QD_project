#include <fftw3.h>
#include "arpack_zn.h"
#include <stdio.h>
#include "auxiliaries.h"

typedef double prec;
typedef complex<double> content;
class region;
typedef region reg1;

//universal constants and "units" derived out of them -
// that is universal constants combined to form the appropriate dimensions
class units_class {
  
  public:
    units_class();
    ~units_class(){};
  
    double m_e, hbar, e, eps, mi;
    double e2prime,kB,lC,meV,nm;
    double Debye, Na, c, aGaAs, muB;
    double m_p, muN, aSi;

    double mass,length,energy,omega,B;
    double A,wavevector,momentum,el_intensity,magneton;
    double velocity,density,temperature;
    double so_strength, so3_strength, g_factor;

    void list(FILE* ,int);

};
 
//physical and computational parameters
//say a parameter (such as magnetic field) has in SI units value x (such as 1 [Tesla])
//the unit corresponding to this parameter is "B" in units_class (such as B=2.64 [Tesla])
//stored in unit_values=2.64 
//parameter is written into this class through set("B",1,"T") whereby a value
//1/2.64 (=0.38) is written into the array values
//by reading through read("B",B) the value 0.38*2.64 (=1) is written into B
//if a parameter is not in SI units (such as energy in meV), then the 
//unit_values should contain the appropriate unit divided by the conversion factor
//for example, for the energy unit_values=units.energy/(1e-3*units.e)
//therefore, if such non SI parameter is used in formulas, it must be MULTIPLIED by the conversion factor
class parameters_class {
  const char* names[100];
  const char* descriptions[100];
  const char* unit_names[100];
  double unit_values[100];

  int name2index(const char* name);

  public:
    double values[100];

    parameters_class();
    ~parameters_class(){};

    enum {
      B,phi_B,theta_B,gx,gy,\
      gz,lv,lB,growth,br,\
      lbr,dress,ldress,dress3,d,\
      d_x,d_y,phi_d,width_z,T,\
      Bosc,phi_Bosc,theta_Bosc,Eosc,phi_Eosc,\
      ring,eps0,density,cl,ct,\
      def,piezo,flat_x,flat_y,exch,\
      xMn,ALHO,nuclear_state,fermi_e,Ebias,\
      m_eff,imp_ext,imp_str,imp_num,wg_x,\
      wg_y,wg_width,wg_phi,wg_V,wg_length,\
      dot_offset,def2,vf,Dg,Di,\
      Dr,gr_dd,gr_w,gr_r0,gr_V0,\
      gr_U,pairing,\
      extent,grid_step,length_x,length_y,dLx0,\
      dLy0,dLx,dLy,dim_x,dim_y,\
      dim,N,dimstep,precision,peierls,\
      eig,eigmax,geometry,leads,threshold,\
      sorting_method,eig2e,selected2e,wg_extent,CE_exp,\
      Ix,sets,system,\
    };
    
    double read(const char* name, const char* unit);
    void set(const char* name, const double value,const char* unit);
    void recompute();
    void list(FILE* file, int where_to);

};




//struct lead

//class leads_class {

//impurities with gaussian potential of random strength and extension, both drawn from uniform distributions
//offset is a total integral of the impurities potential, and should be thus subtracted if zero average is required
//generate() generate a new set of impurities
//potential returns the impurities potential in the plane 
class impurities_class {
  int N;
  prec* strength;
  prec* extent;
  prec* loc_x;
  prec* loc_y;
  prec average;

  void free_arrays();

public:
  impurities_class();
  ~impurities_class();
  void generate();
  prec potential(prec x, prec y);
};

class spin_impurities_class {
  prec *I_old[3];    //expectation values of the spin components (previous)
  prec Arho;         //coupling A [volume x energy] divided by the volume corresponding to a single spin
  class region* r1;//pointer to the region
  int si, sj;           //sigma matrix indexes for the potential
  int circle(prec r, prec h, int i);//returns sign of component i of spin on circle of radius r. If not all the same, returns 0. h is the grid step;
  void update_overlap();    //calculates the overlap <I_old | I> over the whole distribution
  
public:
  prec under;         //underrelaxation coefficient I_new = I_new*(1-under) + I_old*under
  prec *I[3];         //expectation values of the spin components
  prec ave[3], d[3], dphi[3], Qs[3], Qa[3], Qphi[3], Qsph[3], Rmax, overlap; //average, dipole moment and orientation, Quadrupole and orientation, sphericity, radius of maximal density, and overlap with the previous set
  int borders, pattern;        //number of rotationally symmetric borders and resulting pattern class (0 - nothing, 1 - monopole, 2 - dipole, 3 - core and halo)
  spin_impurities_class();
  ~spin_impurities_class();
  void allocate();        //allocate and possible deallocate the I, I old arrays
  void reinitialize(class region* r1, int type=0, prec aux=0);   //set region and initialize spins to : zero (0), random unit vector (1:z component only, 2:all components), 3: gaussian vector with unit dispersion <v \cdot v>=1, 5:same as 1, 6: dipole along x, 7: core+halo with boundary at aux
  int N;         //number of impurity spins equal to the number of active points
  void position(int n, prec& x, prec& y);       //give position of impurity #n in length units
    //set {\bf I}_n along {bf v} (in energy units) using the Brillouin function B(5/2, |v|/kBT) as the multiplicator; mudolation factor is for inhomogeneous profile x_Mn
  void set_impurity(int n, prec * v, prec kBT, prec modul=1.0);
  void set_indexes(int si, int sj);//set the indexes for the potential
  content potential(prec x, prec y);//potential at point (x,y) with sigma matrix indexes si, sj
  void distribution_characteristics(class two_electron_state_prototype* state2e=0);//calculates the spin distribution characteristics possibly using the hole density as a weight
  
};

//used to extrapolate
//length of values with coordinates are stored
//and used to extrapolate a value at given parameter
//or approximate a parameter for a given value
class stack {
  int pointer_to_empty;             //pointer to the first empty storage place
  public:
    stack(int l, int k, char* tag); //constructor
    ~stack();                       //desctructor
    stack* next;                    //pointer to next 
    prec *vals;                     //array of values
    prec *coors;                    //array of corresponding coordinates
    char tag[100];                  //message for debug
    int length;                     //number of storage bins
    int kind;                       //0 - normal 1 - derivative
    int filled;                     //number of valid entries
    prec coordinate(prec value);    //extrapolate coordinate where the value will be reached
    prec value(prec coordinate);    //predict a value at a given coordinate
    void another_in(prec coor, prec value, char* tag);
                                    //fill in another pair od data
    void derivative(prec *coors, prec *vals);    //derivate the data row
    void empty() {filled=0;};
};

//used for adjusting the step
//calling another_in with consecutive number creates (if not existing) stacks keeping track of the values at previous steps
//length: 2 - linear extrapolation, two values enough
//        3 - quadratic extrapolation
//kind:   0 - normal data
//        1 - need derivative before exptrapolating - length needs to be the method + 1
//adjust step goes through all stacks, those filled give step
//the smallest in the range (min, max) is taken
//a step is set shorter_step times shorter to get to the predicted zero
class acr_stack {
  class stack* first;
  class stack* last;
  int hmn_exist;
  stack* give_stack(int which);
  public:
    acr_stack();
    ~acr_stack();
    void another_in(int which, prec coor, prec value,int length, int kind,char* name);
    prec adjust_step(prec position, prec min, prec max, prec shorter_step);
    void reset() {for (int i=0;i<hmn_exist;i++) give_stack(i)->empty();};
    void reset(int i) {if (i<hmn_exist) give_stack(i)->empty(); else {("corrupt in acr_stack\n");exit(1);}};
};

//class for Fourier transforms -- once created, can be called by 
//address_in with the dimensions of the table which returns the address where the input is to be written
//if the dimensions change or the class is not initialized, it initializes
//if the dimensions do not change, it just does the transform without new memory allocations
class Fourier2D {
  private:
    content* fftw_in;
    content* fftw_out;
    int Nx;
    int Ny;
    fftw_plan plan;
    bool Initialized;
    void fill_attenuators();
    prec attenuation(prec theta);
  public:
    Fourier2D();
    ~Fourier2D();
    content* execute();
    content* address_in(int Nx, int Ny);
    prec* attenuation_factors;
    int cons2four(int i, int which);
    int four2cons(int n, int which);
};

struct level_s {
   int parity;
   bool visible;
   prec vals[3];
//   prec pars[3];
   int valid_entries;
   prec forecast;
   int s2u;
   bool assigned;
};

struct level_u {
  prec value;
  int parity;
  int u2s;
  bool assigned;
};

class  sorting_class {
      //concerning parity and sorting

    level_s *level;     //structure level_s for each of the sorted levels
    int seen;           //highest seen level
    bool initialized;    //are the parities of the states already fixed?
    bool allocated;     //is memory allocated?
    bool too_many;      //has seen grown bigger than eigmax?
    int sorted;         //number of sorted states (eigmax)
    int unsorted;       //number of unsorted states (eig)
    
    bool all_definite;  //do all states have definite parity?
    //prec threshold;   //threshold to assign the parity
    bool is_method_zero;//is yes, no sorting what so ever

    int max_valid_entries;//how many values used to extrapolate: 1 only the last, 2/3:linear/quadratic extrapolation 
    
    prec act_parameter;
    //struct state_info* states;
    
    level_u* new_state; //new states being considered;
    
    
    prec* gap_value;         //minimal gap for all the couples
    prec* gap_param;         //minimal gap occured at param
    
    prec pars[3];
    
    acr_stack dog;
    int dog_N;
    int* dog_tag;
    prec step_max,step_min,expand;
    
  public:
    sorting_class() {initialized=allocated=false;}; //constructor
    ~sorting_class(){deallocate();} //destructor
    void init(int sorted, int unsorted, int max_valid_entries_=3);
    void deallocate();
    void save(FILE* file);
    void load(FILE* file);
    
  private:
    void update_history();          //add entries into the arrays of values for visible levels
    void forecast();                //extrapolate values for each level to given parameter
    bool do_I_sort(); 		    //do I sort (true) or skip (false) ?
    int closest_sorted2unsorted(int u, bool ignore_parity);//find closest state according to parity (if not zero for the unsorted) and distance from forecast 
    void copy_results(int* s2u, int* u2s);
    void state_found(int u, int s);
    void update_gaps(prec p);
    
  public:
    // fill in needed values from 1e states
    void inspect_new_1e_states(class state_info* states);
    void inspect_new_2e_states(level_u* levels);
    void sort(prec p, int* s2u, int* u2s);
    
    prec dog_init(prec step_max, prec step_min, prec expand);
    void dogs_fill(prec p);
    prec adjust_step(){return(dog.adjust_step(act_parameter,step_min,step_max,expand));};
    void gap(prec& value, prec& param, int s1, int s2);
};

//collected info about the states
//after diagonalizing in region r1 
//and writing eigenvalues and vectors into vysl
//this class computes, stores and gives various parameters of a state
class state_info {

  public:  
    reg1* r1;
    ret_from_arpack_zn* vysl;
    prec *rates;       //4* unsorted x unsorted table of transition rates:
                       //rate which from state "from" to state "to" is in 
                       //[from*unsorted*4+to*4+which]
                       //where which: {0,1,2,3}={df, pz1, pz2, sum of the three}]
                       //temperature is taken into account (if energy difference is negative, rate is zero
    int unsorted;      //number of unsorted states (eig - smaller)
    int sorted;        //number of sorted states   (eigmax - larger)
    int *unsorted2s;      //unsorted (smaller set) to sorted (larger)
    int *sorted2un;        //sorted (larger set) to unsorted (smaller)
                         //(-1 means the state is not visible any more)
    int length;        //number of characteristics
    bool* show;        //what will be shown by statistics
    class Fourier2D Fourier; //"static" class for repeated Fourier transforms
  
  private:
    
    //concerning characteristics
    prec** x;            //characteristics of the state (such as spin)
  //prec* rates;         //transition rates

    const char **tags;         //tags (such as sz)
    const char **names;        //names (such as spin along z)
    const char **unit_names;        //unit labels
    prec *unit_values;   //units to convert from internal units into SI
    bool **actual;       //is the value actual?
  
    void compute(int i, int t); //compute a value with a given tag label
    int find_tag(const char *tag); //find tag identifier
    
  public:

    prec omega(int i, int j);//(E_i-E_j)/hbar in intertnal units
  
    void show_set(const char* what);
    void list_tags(FILE*, int);//list for each tag: number, tag, unit, value
  
    void new_states(); //new states -> set all flags on tags on all states to 'not computed yet'
  
    prec read(int u, int t);    //read a value with a given tag label
    prec read(int u, const char *tag);    //read a value with a given tag
  
    void tag_info(int t, const char* &tag, const char* & name, const char* &unit);
                          //find info about a tag
    state_info(reg1*, ret_from_arpack_zn*, int eig, int N);   
                          //constructor
    ~state_info();          //destructor
    
};

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!density matrix for coherently injected spin
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

class rho_w2special 
{
  int i,j;        // the two special states where the electron is be injected
  int n;          //number of considered states
  content* A;     //matrix Lambda
  content* U;     //right eigenvectors of Lambda
  content* Um1;   //the inverse of U
  content* d;     //eigenvalues of Lambda
  int N;          //dimension of matrix Lambda = n*n
  //ret_from_arpack_zn vysl;   //eigensystem of Lambda
  content* S;     //matrix elements of the sigma matrixes
  state_info* states;
  content* rho0; //initial density matrix
  content* rho;  //density matrix at time t

  int two2one(int i, int j);
  void one2two(int n, int &i, int &j);
  content matrix_element(int i, int j, reg1::mean_value_of_operator op);
  void fill_matrix_elements();
  
  void fill_lambda();
  void diagonalize_lambda();    
  prec check_eigensystem(int where);  
  prec check_inversion(int where);
  prec compute_rho(prec t);
  prec give_spin(int which);

  public:
    rho_w2special(int i_in, int j_in, int n_in, state_info* states_in);
    ~rho_w2special();
  
    void inject_spin(prec theta, prec phi);
    void solve();
    content give_rho(int i, int j);
    void spin_at_t(prec &sx, prec &sy, prec &sz, prec t);
};

/*
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h" 
#include "gsl/gsl_complex.h" 
#include "gsl/gsl_complex_math.h" 
*/

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!extended and evanescent waves in 2D with spin orbit interaction
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

extern content pauli[];  //Pauli matrices 1, x, y, z - classes.cpp


struct spinor {
  spinor(content up=0, content down=0) {this->up=up;this->down=down;};
  content up, down;
  void angles(double& theta, double& phi);
  spinor SigmaTimesMe(content* Sigma);
  spinor add(double coef, spinor s);
  spinor multiply(double coef);
  spinor conjugate();
  content braket(spinor ket);
  void construct(double theta, double phi);
};

class so_wave_class
{
  friend class so_waveguide_class;
   //all "human" inputs/outputs are in SI units (meV, nm, ...). units.length, units. energy are used to convert to dimensionless units
   //everything is stored in internal units
                             //k vectors (2D, along the waveguide axis): inicidet wave -real, reflected +real, transmitted -imag
   content k0[2],k[2],kpert[2]; //initial guess, final result and perturbative result for a complex vector of a wave
   double e;                 //E-V
   double spin;              //spin aligned/antialigned to the "field" (+1, -1)
   spinor xi;                //spinor components
   double theta, phi;        //spinor angles
   double b[3];              //magnetic field
   double alpha, beta, gamma;//so strengths
   double csg;		     //for 1D integration
   content pi_field[3];      //r, theta, phi of the complex spin-orbit field [in energy, 1, 1]
   content espin_ref;        //spin complex energy set for reference at the beginning of iteration
   content espin;            //spin complex energy - iterated
                             //which object to use to call Hamiltonian for the 1D integration
   static so_wave_class* derivs_object;
   double energy1D;         //energy for the 1D integration
   double (*V)(double x);    //potential for the integration
   spinor* value;            //values of the numerically integrated wavefunction
   double* coor;             //coordinates of the numerically integrated wavefunction
   int buffer_length, filled;//number of allocated/filled entries in the above arrays (the last entry in the array has index buffer_length/filled-1)


public:
   //constructor
   so_wave_class();
   //destructor
   ~so_wave_class();
                             //initialize variables, convert to dimensionless units
   void initialize(double e_, double spin_, content* k0);
   void iterate_for_k(content esigma_ref_opposite);     //find the exact wavevector
                             //read out the result by "name" and "unit", convert to dimensionfull units as you go
   content result(const char* name, const char* unit);

private:
                             //input did not pass through the check - exit
   void wrong_input(const char* name, double value);
   void construct_spinor();  //compute the two spinor components and angles of the spin up direction
   void extract_vectors();   //compute angles/moduli of spin-orbit complex "magnetic field" vectors
                             //spin "energy"
   content esigma(content* k); 
                             //spin-orbit induced complex magnetic field
   void pi(content* k,content* result);
                             //velocity - absolute value, depends on the spin, defined through v=i/hbar [H,x] = \partial_p H
   prec velocity();
                             //Hamiltonian for 1D integration: give ddpsidxdx based on psi and dpsidx at position x [internal]
   void Hamiltonian(spinor& dddpsidxdxdx, spinor& ddpsidxdx, spinor& dpsidx, spinor& psi, double x);
                             //integrate with initial condition, all inputs in internal units
   void integrate1D(double energy, double (V)(double x), struct spinor& ddpsidxdx, struct spinor& dpsidx, struct spinor& psi);
                             //interface for the gsl routine
   static int derivs(double x, const double* y, double* f, void *params);
   spinor value_at(double x); //value of the numerical solution in the intermediate region x in internal units
};

class so_waveguide_class
{
                             //total energy, transversal mode energy [internal]
  double energy_tot, energy_trans;
  int channel;               //transversal quantum number [wavefunction_y = sin(n pi / w)]

public:
                             //incoming, reflected and transmitted waves [spin up/down of the particular wave]
  so_wave_class wave_i[2], wave_t[2], wave_r[2];
                             //transmission and reflection amplitudes (transmitted/reflected to incident)
                             //indexed: [incident][transmitted/reflected]
                             //flux amplitude(s1 -> s2) = amplitude(s1 -> s2) * sqrt(velocity2/velocity1)
  content r[2][2], t[2][2], r_flux[2][2], t_flux[2][2];
                             //transmission and reflection amplitudes: (incident/reflected to transmitted)
                             //indexed: [transmitted][incident/reflected]
  content c_i[2][2], c_r[2][2];

private:
  double (*V)(double x);    //potential for the integration, V(0)=barrier+dot offset and V(length)=0 input/output in internal units
  bool timesdV_sc;          //multiply the scattering wave by the potential difference?
  int spin_sc, component_sc;//spin and the component of the scattering wave
  static so_waveguide_class* object_sc;
  content Sigma_sc[4];      //operator for the mean_spin (equal to unit matrix, or \sigma \cdot {\bf n} )

public:
                            //initialize waveguide mode - chanel from 1,..., returns true if the lower in energy spin is propagating
                            //energy is in meV, V(x) is in internal both x and returned value
  bool set_channel(double energy_, int channel_, double (V)(double x));
                            //integrate for both spins, compute transmision amplitudes
                            //returns the total (flux) transmission probability
  void propagate_through_barrier();
                            //construct scattering wave for: incident spin, spin component of the wave (0, 1) and multiply by the potential difference *(V(0)-V(x)) if true
  void wavefunction_set(int spin, int component, bool timesdV);
                            //for the outside region_2d routines, all in internal
			    //x y are dot coordinates
  static content wavefunction(double x, double y);
                            //returns overlap of the state wavefunction #unsorted with the scattered wave of "spin" times the potential difference
                            //because of the normalization of the scattered wave, the dimension of the result is inverse length times energy squared -- it should be multiplied by 1/units.length * units.energy^2
  content wg_overlap(state_info* states, int unsorted, int spin);
			    //tunneling rate in 1/s for unsorted state and incident spin
  double tunneling_rate(state_info* states, int unsorted, int spin);
                            //spin quanitzation direction for (the unnormalized spinor) the scattering wave of incident spin s at position x, wavefunction (<psi(x)|psi(x)>) probability is returned
  double spin_expectation(double x, int s, double& theta, double& phi);
                            //constructs the scattering submatrix from the transmission coefficients
                            //according the formula t_ab= sum_xy <a|transmitted> <incident|b> t_{trans inc}
                            //and a,b = up/down along the axis (spinor coefficients) defined by (theta, phi)
  void scattering_submatrix(content* T, prec theta, prec phi);
                            //construct scattering states for integration of the spin (op=0, theta, phi define axis) or probability (op=1) mean value
  void mean_spin_set(bool op, prec theta, prec phi);
                            //returns the mean value of the previously set quantity
  static content mean_spin(prec x);
  prec mean_spin_in_channel(char channel, int spin);

                            //returns the spin polarization and magnitude of the outgoing particle flux (for unit incoming flux) and transmission matrix elements
  //these are defined by Gs = <scattering state| {v_[at out going lead],sigma}/2 | scattering state> /C
  // G = <scattering state| v_[at out going lead]| scattering state> / C
  // where C=<scattering state| v_[at incoming lead] | scattering state>
  void current_polarization(prec &G, prec &Gs, prec* t2);

private:
                            //extract transmision amplitudes from the value of the spinor at x=wg_length
                            //this final value originated from transmitted wave of "spin"
  void extract_transmissions(spinor psi, spinor dpsidx, int spin);
                            //convert transmision coefs from all/transmitted to all/incident (c_i, c_r -> r, t)
  void convert_transmissions();
                            //construct spinor corresponding to incident spin wave scattered, at position x, y, possibly times the potential difference, x, y are waveguide coordinates
  spinor scattered_wave(double x, double y, int spin, bool timesdV);
  spinor scattered_wave_x(double x, int spin);
  double scattered_wave_y(double y);

  //computes current at the interface
  void current();
                            //check solution of linear set of equations A[NxN].root[N]=b[N]
  void check_solution(int N, int INFO, content* A, content* root, content* b);

};

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!transport: recursive green's function
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

struct s_matrix_prototype {
  int* M;               //number of channels in each lead
  content** s;          //s[i,j]=submatrix for the transmission from the lead j to the lead i
                        //s_matrix.s[lt*Nl+lf][mt*M[lf]+mf]=amplitude(lf,mf,lt,mt);
};


struct mode_prototype {
  //dimension equals to the number of points of the corresponding lead times 2 for the spin
  content * wavefunction;//wavefunction for each of the mode 
  prec energy;          //energy of the mode
  prec velocity;        //longitudinal velocity of the mode
  prec wavevector;      //longitudinal wavevector
  content *longitudinal;//the longitudinal part of the lead.mode green function, the dimension is b*b
  int trs, orb, spin;   //labels distributing the modes into spin opposite couples, according to trs symmetry
  content trs_phase;    //overlap of this mode time reversed and the trs pair: <m_trs|T m_this> 
  spinor normalized_spinor;//spinor part of the mode
};

struct lead_prototype {
  int heading;          //0/1/2/3 towards -x/x/-y/y axis
  prec x,y;             //coordinates where the lead axis touches the inside fo the region(in internal length units)
  prec width;           //width (from lower towards the positive axis values; in internal length units)
  prec potential;       //potential of the lead
  int points;           //number of modes in the lead (equal to number of points in perpendicular direction)
                        //number of propagating modes and propagating orbital modes (not counting the spin)
  int propagating_modes, orbital_modes;
  mode_prototype* mode; //dimension equals modes, but only the propagating ones are initialized
};

struct transport_class {
  int b;                //# of nearest neighbours
  int Nl;               //number of leads
  prec hx,hy;           //grid step
  reg1* region;         //region with the leads included
  reg1* region2;        //region with the Hamiltonian without spin-orbit (for the self energy)

  int* dim;             //number of active points on xlines
  int xlines;           //number of xlines with active points

  bool method_full;     //whether the copmlete greens function is computed or just the first-to-last point
                      
  //wavefunction is indexed by (xline, sequential label, spin) 
  //operators (Hamiltonian,...) are then indexed by two such (multi)indexes
  // and stored as a matrix of matrices:

  //xxx(i,j)= matrix for points on xline(i) \times xline(j), that is
  //[i*xlines+j] is pointer to matrix[dim(i)*2*dim(j)*2]
  content **self_array;       //self energy
  content **g_array,**g_aux_array;  //green functions

  content *self_array_base, *g_array_base, *g_aux_array_base;
  bool *self_filled, *g_filled, *g_aux_filled;

  content* auxM;//auxiliary matrix
  
  //struct lead_prototype lead[4];         //leads
  struct s_matrix_prototype s_matrix;    //the s-matrix

//!internal routines: 
                     //!matrix manipulations

  //converts two multiindexes into a label into a matrix stored as a row of numbers (vector)
  int label(int i, int j, int yi, int yj, int si, int sj);
  //returns label in the matrixes space
  int labelM(int i, int j);  
  int labelg(int up, int i, int j);
  //multiplication of two matrixes with given multiindexes
  void MM(content* res,content* a, content* b,int i1,int j1, int i2, int j2);
  //add two matrixes res=res+in with multiindexes i,j
  void add(content* res, content* in, int i, int j, prec coef);
  //add a constant to a diagonal matrix with multiindexes i,j res=res+x*id 
  void addC(content* res, content x, int i, int j);
  //compute norm of a matrix in
  prec normM(content* in, int i, int j);

                     //!greens functions, hamiltonian, self-energy
  //functions which return the pointer to the greens function / GF aux / self energy matrix corresponding to xlines i,j
  //and the values of this matrix specified by the y and spin coordinates
  //content* g(int i, int j);
  content* g(int up, int i, int j);
  //content g(int i, int j, int yi, int yj, int si, int sj);
  content g(int up, int i, int j, int yi, int yj, int si, int sj);
  //content* g_aux(int i, int j);
  content* self(int i, int j);
  content self(int i, int j, int yi, int yj, int si, int sj); 
  //true if both the hamiltonian and self energy are zero for these two xline indexes
  bool iszeroh(int i, int j);
  //if the two points are further apart then b in any direction returns true
  //input are sector coordinates
  bool iszeroH(int s1, int si1, int sj1, int s2, int si2, int sj2);
  //true if both the hamiltonian and self energy are zero for these two indexes
  bool iszerov(int i, int j, int yi, int yj);
  //matrix element of the Hamitlonian plus the self energy
  content HpS(int i, int j, int yi, int yj, int si, int sj);
  //fills the matrix representing the E-H-Sigma of the region without leads part (i,j)
  void fill_EmHmS(content* res, int i, int j, prec Ef);
  //returns the longitudinal part of the green function of a semiinfinite lead
  content leadg_long(int i, int j, prec k, prec v);
  //computes self energy
  void compute_self();
  //multiply off diagonal part of the Hamiltonian conecting the last added line i with the green function: 
  //result=matrix (i,j)=sum_ip V_iip G_ipj, ip<i
  void VG(content* res,int i,int j,bool update,content coef);
  //multiply off diagonal part of the Hamiltonian conecting the last added line j with the green function: 
  //result=matrix (i,j)=sum_jp G_ijp V_jpj, jp<j
  void GV(content* res,int i,int j,bool update,content coef);
  //multiply off diagonal part of the Hamiltonian conecting the last added line i with the green function and the hamiltonian again: 
  //result=matrix (i,j)=sum_ipjp V_iip G_ipjp V_jpj
  //note it must be i=j
  void VGV(content* res,int i,int j,bool update,content coef);


                    //!allocation, initialization, clean-up
  //construct the region including leads
  void prepare_regions();
  //allocate space for the g-function, according to the method
  void allocate_g_function();
  //to do once - the geometry dependent (energy independent) tasks - initialize the method (full/linear)
  void prepare_geometry(bool full);
  //the energy dependent tasks
  void prepare_with_E();
  //deallocate space from geometry dependent tasks
  void clean_geometry();
  //deallocate space from energy dependent tasks
  void clean_with_E();
  //geometry dependent tasks for leads
  void leads_geometry_init();
  //energy dependent tasks for leads
  void leads_with_E_init();
  //initialize the lead
  void set_lead(int l, int heading, double x, double y, double width, double V);

                  //!auxiliaries
  void check_invert(content* a, content *inv, int i);
  //checks that the greens function is the inverse of E-H-Sigma cut after line i. Computes the norm of
  //[G (E-H-S)]_{jk} - delta(j,l) for each j,k from 0 to i (including)
  //result is the list of residual norms into the log
  void check_it(int i) ;
  //print the upper left entry of the three matrixes
  void print_matrixes();
  void cycle(int i, int &j); 
  //returns y label of a point in the lead neighbouring with point y on xline inside
  //if no, returns -1 - assumes no holes in the lead
  int jcor(int l, int xline, int y);
  //true if the point is nearby (within b steps from) the lead l; input are r_wo indexes
  bool point_nearby_lead(int l, int ai, int aj);
  void draw_potential(FILE* file, int where);
  //set limits for region without leads i-labels neighbouring a lead heading to heading
  void amplitude_aux1(int heading, int& ifrom, int &ito); 
  //set limits for region without leads yi-labels neighbouring a lead heading to heading on line i
  void amplitude_aux2(int heading, int i, int& yifrom, int &yito); 

  //makes a time reversal of a mode of lead and returns the mode to whihc it is transformed
  //returns the overlap <mode_out | T mode_in>
  //if TRS mode is not in the lead (not porpagating), returns -1 in mode out and zero for the coefficient
  content mode_trs(int lead, int mode_in, int& mode_out);
  //check for self-duality (TRS symmetry) of the scattering matrix (logs residuals)
  prec check_trs();
  //assign modes according to spin pair - using the trs symmetry
  void pair_modes();


  //checks for the unitarity of the scattering matrix, returns the residual
  prec check_unitarity();

                    //!high level routines
  //compute the green function between the first and the last xline using recursive method
  //full -- compute the full inverse matrix (true) or just the four corners (false)
  void recursive_green(bool full, progress_bar_prototype* progress_bar=0);
  //tranmission amplitude from lead, mode, spin into a second combination of these three
  content amplitude(int lead_from, int mode_from, int lead_to, int mode_to);
  //fills 2x2 scattering submatrix (S allocated outside) indexed by lead and orbital(!) mode
  //if the time reversal pair is not propagating, the corresponding transmission amplitudes are zero
  //returns the phase = phase_to/phase_from
  content scattering_spin_submatrix(int lf, int mof, int lt, int mot, content *S);
  //sums all mode to mode transmission according to the trace formula
  //T_mn= Tr[t_mn^dagger sigma_alpha t_mn sigma_beta]
  //indexed by lead (from, to) and sigma = 0,1,2,3 = id, x,y,z are Pauli matrices
  prec lead2lead_spinresolved(int lt, int alpha, int lf, int beta);
  //check the algorithm on a series of step barriers - requires change of the potential in hamilt.cpp
  void check_GF(FILE* file, bool full);
  //extract transmission and reflection amplitudes for a one dimensional case (all xlines have length (dim) 1)
  //(note: valid only for spin along z)
  void oneDamplitude(content &t, content &r);
  //outputs spin density of the scattering state as a function of position
  void oneDspindensity(FILE* file);
  //fills in the scattering matrix
  content** scattering_matrix(int* M);
  prec lead2lead(int from, int to);
  prec trace(int from, int to);

  //sets charge current through lead 3 to zero, computes the spin current through lead 3, 
  //in addition, computes the derivative of charge current with respect to Bx and By (numerically applying field Bstep)
  //input: voltages V[1], V[2]
  //output (all concern lead 3 at B=0): results[0-4]={dI/dBx, dI/dBy, I_spinx, I_spiny, I_spinz}
  void derivative(prec *V, prec Bstep, prec* results);

  //tunes GF_check_Eb such that the total transmission through the QPC is maximal
  //in two stages:
  //1. steps GF up by Vstep1 until the transmission is below threshold1 (in maximal maxsteps1)
  //2. steps GF down by Vstep2 until the transmission derivative starts to drop and the transmission itself is at least threshold 2 (in maximal maxsteps2)
  //returns the transmission derivative and the actual number of steps
  //estimate is a value to set for the GF if the search fails
  //T is transmission at the maximum derivative
  prec set_QPC_sensitive(prec Estep, prec Vstep1, prec threshold1, int& maxstep1, prec Vstep2, prec threshold2, int& maxstep2, prec estimate, prec& T);

  //computes the derivative wrt energy of the total transmission through the lead 3 using difference with the step Estep
  //T is the transmission
  //result is in 1/meV units
  prec QPC_derivative(prec Estep, prec& T);

};

class transport_window_prototype {
  prec E1, E2;                    //interval for which the transmissions are given (in meV)
  int NP;                         //number of point for the interpolation
  int Nl;                         //number of leads
  prec step;                      //step in energy (E2-E1)/NP-1 => E(i) = E1 + step*i; i=0,...,NP-1 (in meV)
  prec **t;                       //transmissions (interpolated): T_{ij}^{ab}(e) -> t[(i*Nl+j)*16+a*4+b][i] (dimless)
  int **N;                        //number of open channels for each lead and spin: N_i^a(e) -> N[i*4+a][i]
  prec *mu;                       //spin resolved chemical potentials of the leads mu_i^z -> mu[i*2+z] (in meV)
  int spin_axis;                  //spin quantization axis for the spin current formulas (meaning of "z", 1/2/3)

public:
                          //constructor/destructor
  transport_window_prototype(prec E1_, prec E2_, int NP_);
  ~transport_window_prototype();
  void set_spin_axis(int a);                         //a=1,2,3, defining the direction meant by the "z"=1
  void input(int np, prec* e, int** N, prec** t);      //transform the ordered set of data into equally spaced set
  void set_mu(int l, int z, prec mul);               //set the chemical potential of the lead
  prec zero_current(int l, int z);                   //finds the (spin) chemical potential such that the (spin) current is zero
  prec current(int l, int z);                        //(spin) current through lead l in units of e^2/h
};

//!experimental
class CE_prototype {
  state_info* I;                           //the corresponding states and region structures
  int N, Ne;                               //number of states/entries [the latter is N*(N+1)/2]
  int Nx, Ny, Nxa, Nya, blowup;            //dimensions of the non-expanded and expanded tables
  int length_t;                            //length of the expanded table
  prec hx, hy;                             //grid steps
  prec eps0;

  complex<double>** FT;                    //Fourier Transforms of state pairs - array of pointers
  bool *FT_computed, *FT_to_compute;       //whether the specific pair is / will be computed
  bool remmember;                          //whether accumulate FT pairs in memory

  complex<double>* coefs;                  //coefficients for the FT pair integration
  bool CoefsInitialized;                   //whether they are calculated

  complex<double>* aux_region;             //auxiliary space - the region
  complex<double>* aux_table;              //auxiliary space - the non-expanded table
  complex<double>* aux_exp_table;          //auxiliary space - the expanded table (zeros added to the previous one)
  int *kinv;				   //table of inverse momenta kinv[k]=-k

  int ij2seq(int i, int j, bool& swap);     //converts i, j into sequential label, indicates whether they have been swapped
  void create_table(int i, int j);         //create the expanded table for state pair i, j
  complex<double>* compute_FT(int i, int j, bool& cc);//returns pointer to the pair FT,computes if not yet, indicates whether complex conjugation (and inversion of the k vector) is needed
public:
  CE_prototype(state_info *I, int N, int blowup, double MB);//alocate memory except for the FT pairs
  ~CE_prototype();                         //deallocate all memory
  void deallocate();                       //deallocate memory for the FT pairs (check if allocated)
  void reset(state_info* I);               //delete all "computed" flags, renew the state info (do not deallocate FT memory)
  int dry_run(bool reset, int ui, int uj, int uk, int ul);//set flags for those that will be required to be computed
  bool memory(double MB);                  //decide whether accumulate in memory or not, MB is the limit
  complex<double> CoulombElement(int ui, int uj, int uk, int ul);//Coulomb element of (unsorted) states C_ijkl in meV  
  complex<double> CoulombElementTweak(int ui, int uj, int uk, int ul);//Coulomb element of (unsorted) states C_ijkl in meV
};

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!graphene: generalization of wavefunction degrees of freedom
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

struct dof_prototype {
  int Nc, length;						//number of components / maximal value of the consecutive index
  //component names; end sets the boundary
  //!for graphene sigma dof is special and has to be the first!!!
  //enum name {spin,end,sigma};		//GaAs
  //enum name {end,spin,sigma};			//laser and two electron case
  //enum name {sigma,end,spin,endx,xspin};			//graphene
    //enum name {spin,tau,end,sigma};			//skyrmion
  enum name {end,spin,tau,sigma};			//skyrmion check1
  enum op {unity, paulix, pauliy, pauliz, splus, sminus};			//operator basis within a single dof
  int *indexC;							//values of the component vector (+1/-1 for each of Nc)
  int index;							  //consecutive index
  void convert(const char* what, complex<double> coef=0);	//conversion of: "C2I/I2C" : indexC into index and back, "OP2M" : operatorC into matrix M - here if c!=0, update M_sum, too
  int* operatorC;						//tensor product of an operator (one of four) for each dof
  complex<double>* M, *M_sum;					//matrix corresponding to the operator term, sum of matrixes
					//the matrix M is defined as M = M_IJ |I><J|; M_IJ=M[i*length+J]
  //void addM(complex<double> coef)				//add the actual matrix to the sum
    //{for (int i=0;i<length*length;i++) M_sum[i]+=c*M[i];};
  void reset(const char* what);					//zero the matrix sum
  dof_prototype();						//constructor
  ~dof_prototype();						//destructor

  int* braC, *ketC, *pntr;					//auxilliary
  complex<double> sigma2ketbra(int sigma, int term, int& ket, int& bra);
  void convert_log(const char* mess, const char* what, complex<double> coef=0);
};


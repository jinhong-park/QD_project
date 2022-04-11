//#include <stdio.h>
//using namespace std;


#include "region2d.h"
#include "arpack_debug.h"
#include "arpack_zn.h"
//#include "auxiliaries.h"

typedef region reg1;
typedef double prec;
typedef complex<prec> content;

#include "classes.h"
#include "two_electron2.h"

extern prec E0;

bool (inside)(prec, prec);
int (inside_sector)(prec, prec);
prec (potential)(prec, prec);
prec (potential_gr)(prec, prec);
prec (mass_gr)(prec, prec);
content (potential_grc)(prec ,prec );
content (mass_grc)(prec ,prec );

extern class parameters_class parameters;
extern class units_class units;
//extern class sorting_class sorting;
extern class impurities_class impurities;
extern class spin_impurities_class spin_impurities;
extern struct lead_prototype lead[5];
extern struct laser_parameter_prototype LP;
extern struct dof_prototype dof;

content xc(prec x, prec y); 
content yc(prec x, prec y);
 
//!************** hamilt.cpp **************!//

content exp_integral(prec, prec, prec, prec);

//x,y su pointery na vektory o dlzke 2*N=(N*spin up,N*spin down)
//A je ukazovatel na reg1 - region, ktory je volany na vykonavanie jednotlivych algebr. operacii na vektore x
//vysledkom je prepisany vektor y=H x
void hamiltonian1(int N, void *A, content *x, content *y);
void hamiltonian1_init(reg1&);

void hamiltonian(int N, void *A, content *x, content *y);
void hamiltonian_graphene(int N, void *A, content *x, content *y);
void hamiltonian_graphene_with_checking(int N, void *A, content *x, content *y);
void hamiltonian_init(reg1&);
void hamiltonian_init_graphene(reg1& region);
void hamiltonian_for_lead_init(reg1&, int heading);

//podla clanku SS vyrata konstanty:
//ab, ad, ak E<>0; av, omega_u ak lv<>0; length vzdy
//print==true - aj ich vsetky vypise
void recompute_constants();


//prehodi obsah N prvkovych vektorov x,y
void flip_0(int N, content*x, content *y);

//usporiada vlastne cisla a vektory(o dlzke N, pocet eig) adresovane vo vysl podla (method='a/r/i/R/I' : real, imag, abs value) (v)zostupne (order='a/d')
void sort_0(int N,int eig, ret_from_arpack_zn* vysl, char method, char order);

//normalize a set of eig complex vectors of lentgh N such that sum_k |vec[i*N+k]|^2 = 1 for each i
void normalize(int N, int eig, complex<double>* vecs);

//for 'eig' eigenvectors uses function find_phase to find a common phase
//then divides by this phase
//wch, lev - arguments for phase
void fit_the_phase(ret_from_arpack_zn& vysl,reg1 &r1, int eig, int wch, prec lev);

//picks a point on the x line at some distance from -d_0
//find a value of function, if bigger than lev, this is the phase
//if not moves certain distance to the right
//wch - (0,1) - check only spin (1,-1);  (3) - check both
content find_phase(reg1& r1, ret_from_arpack_zn& vysl,int eig, int wch, prec lev);

//update linear dresselhaus constant
void dresselhaus(prec kappa);

//potential of the waveguide barrier - input/output in internal units
double potential_waveguide(double x);
double tunnel_barrier(double x);

//!************ energy.cpp ***********!//

//energia analytickeho vl. stavu s kvant. cislami n,l,s
prec energy(int n, int l, int s);

//energia analytickeho vl. stavu s kvant. cislami n,l,s
prec energy(prec n, prec l, int s);

//Zeeman energy in interal units
prec zeeman();

//correction to the state |n,l,s> according to corr
//1 : BR-BR
//2 : D - D
//4 : BR-D3
//8 : D -D3
//16: D3-D3
prec correction(int corr, prec n, prec l, int s);

//if which <0 : from ranges n <0,n>, l <-l,l>, s <-1,1>
//sort values of energy+correction according to energy+(wc==true)*correction
//if which >0 return previously computed values and quantum numbers in n,l,s
//energiu pritom rata s korekciami podla corr
prec level(int which,int corr, int& n, int &l, int &s, bool wc);

//graf najnizsich (hmn) anal. hladin s vystupom do file/where_to
//n,l -maximum of n,l up to which to go
//numbers 1+2+4: print out also the quantum numbers /1-n, 2-l,4-np,8-nm,16-s/
//wc,corr - ako v predchadzajucom
void levels(FILE* file, int where_to,int corr, int hmn, int n, int l,int numbers, bool wc);

//computes correction to energy in the case of a ring
//which = 0 magnetic field independent part
//which = 1 field dependent part
prec ring_correction(prec lz, int which);

//zapise do output Fock-Darwin funkciu
void anal_function(reg1& r1,content *output,int n, int l, int s, prec d_0,prec sign, reg1::setting_or_adding soa);
    
//prenasobi funkciu aby norma bola 1
void norm_function(reg1 &r1,content *output);
    
//zapise do vysl(eigenfuctions aj eigenvals) eig najnizsich rieseni single dot  problemu
void anal_spectrumSD(ret_from_arpack_zn & vysl,reg1& r1,int eig, prec d_0);

//zapise do vysl(eigenfuctions) eig rieseni double dot problemu usporiadanych podla poradia v single dot
void anal_spectrumDD(ret_from_arpack_zn & vysl,reg1& r1,int eig, prec d_0);


//returns 'parity' according to eigenvalues of operators inv_x, inv_y, spin
// method==1: +1 +1 ->1, +1 -1 ->2, -1 -1->3, -1 +1 ->4
// method==2: only inv:   +1 -> 1, -1->2
// method==5: np, nm
// method==3: lz
// method==4: all are one (spin can be considered)
//if method<0: consider also spin: method*=spin
// where eigenvalue is decided to be definite, if differs from whole number by less then eps
//if not definite parity=0
int give_parity(state_info* states,int state, prec eps, int method);

//find a position of a given state
//invesrion to position from sort_after_parity
//int give_pos(int j,int N,int* position);

//fit values according to method (1-linear, 2- kvadratic) and predict value at parameter at_param
//reset is reset
//OUTPUT - result[eig] - first method times is set to 0
//INPUT - values in vysl at parameter param
//number of states - eig
//void predict(ret_from_arpack_zn& vysl, int eig, prec param,prec at_param, int method ,prec* result, bool reset=false);

//extrapolates a row of values with coordinates by linear (1)
//or quadratic(2) method and gives two (c,b) or three (c,b,a)
//coefficients of the quadratic fucntion ax^2+bx+c
void extrapolate_coefs(prec *coors, prec *vals, int method, prec *coefs);

//extrapolates values with coefficients to reach zero by
//(1) linear or (2) quadratic method (gives the root closer to last coordinate or zero if there are no roots)
prec extrapolate_to_value(prec value, prec *coors, prec *vals, int method);

//extrapolates the value of the parameter at which there is a complete anticrossing (that is unperturbed energies are exactly equal)
//the stack must be of a minimal length 3
//there should be one half of differencies in the energy insereted into vals at the corresponding value of the parameter stored in coors
prec extrapolate_to_acr(prec *coors, prec *vals);

//write info about a state
//format: +1 - header +2-maximal precision (9) +4 new line at the end
//what: defined in state info, added by tags by show_set(tag)
//show_set(clear) = clear all flags
void statistics(state_info& I, int state, int format,  FILE* file, int where_to);

//writes for all states info defined by tags in state_info 
//format: +1 - header +2-maximal precision +4 new line after every eigenvecs +8 newline at the end
//+16 only visible states
//what_global (encoded binary): 
//1 zeeman energy 2 tunneling energy 4 geff
//8 times.main 16 times.user 32 perform.iter 64 perform.opx
void statistics_multi(state_info& I, int what_global, int format,  FILE* file, int where_to);


//fills two by two complex matrix Sigma using the four vector vec according to
//Sigma= sum_i vec_i . Pauli_i, where Pauli_3 = identity
void vec_dot_sigma(content* vec, content* Sigma);

//converts (x,y,z) to (r,theta,phi)
void vector2angles(double* in, double* out);

//convert a two component spinor (spin up, spion down) into (theta,phi) along the spinor quantization axis
void spinor2angles(content* xi, double& theta, double& phi);

//rotate a two component complex vector by "angle"
void rotate_inplane(content* kin, double angle, content* kout);

//conversion between the dot and waveguide coordinates
void waveguide2dot(prec xw, prec yw, prec& x, prec& y);
bool dot2waveguide(prec x, prec y, prec& xw, prec& yw); //returns true if the point is inside the waveguide

//rotates the spinor according to transformation removing the linear spin-orbit terms
void UnitaryTransformation(prec x, prec y, spinor& xi);

//the orbital ground state in harmonic potential - the (real) Gaussian
content psi_ground(prec x, prec y);

//unitary transformation applied on the ground state gaussian of spin "psi_ground_spin"
//the "psi_ground_component" component of the spinor is returned
//the final complex number is raised to power "psi_ground_power"
extern int psi_ground_spin, psi_ground_component,psi_ground_power;
content psi_ground_app(prec x, prec y);

//!************* fit.cpp ***************//

//richardson interpolation
//N - number of measured points, E - measured values, h - values of the parameter
//outcome[4] - interp from 3 point, from 2 point (best), from 2 point (mean), single value
// allocated outside!!!
//return value - the lowest existing interpol in the previuos row
//at lowest parameter, if 2 point interp exists otherwise 0
//returning value the best interp. extracted
//tol - the intrap. value will be exepted, if its distance from the best single value is less than tol*(max single valu- min)
//e - exponent of the first error term, s - step in exponents in error terms
prec richardson(int N,int e, int s, prec *E, prec *h, prec tol, prec* outcome);


//N-pocet nameranych udajov, E-N hodnot energie(budu sa fitovat), h - hodnoty parametra, tol - ak odhad < tol* rozpatie hodnot E, odhad sa prijme
//ak bol prijaty odhad z trojic, vrati ho a v guess2 je odhad z dvojic, ak nebol prijaty z 3,len z 2, vrati z 2 a guess2=0, inak vrati hodnotu nameranu pri najmensom parametri
prec fitit(int N, prec *E, prec *h, prec tol, prec& guess2,FILE* file, int where_to);

//eig kolko eigenvalues, N kolko merani, e exponent prveho clena chyby, s krok v exponentoch pre dalsich clenov, dim pri akej dimenzii zacat, dimstep  aky je krok, vysl kam zapisat vysledky (musi byt vonku alokovane miesto! - ako pri poslednej (najvyssej) dimenzii=dim+(N-1)*dimstep -> region -> size=2*r.Nw ->
//vysl.eigenvals[eig], vysl.eigenvecs[eig*size]
int fitvalues(int eig, int e, int s, ret_from_arpack_zn *vysl, int N=1, int stepdim=1);

//create adequate region
//and allocate space into vysl
//recompute_constants should be called just before

//deletes eigenvalues and eigenvectors, vysl, region
void clean_up(state_info* states);

//creates region, vysl, and state_info
//diagonalizes hamiltonian, results are in vysl
state_info* diagonalize();

//find a minimal value for a difference between energies of states a, b
//parameter is given by an address - has to be of type prec
//start, step - for paramter
//start_value - difference found at start
//routine supposes monotonous behaviour
//dimstep, N - params for richardson
prec minimize(int a, int b, int eig, prec* x, prec start_val, prec step, int N=1, int dimstep=1);

//if energy differences between states indicate a possible crossing,
//measure energy difference at the crossing
void spot_crossing(state_info& I, prec *param_to_change, FILE* file, int where_to);


//!*********** test.cpp **************//

//run 10 checks for 2 random vectors <v |H| w> ?= conj(<w |H| v>)
void check_for_hermitivity(reg1& r1, FILE*, int);

//test individual operator
void test2(int dimx, int dimy, reg1::operators op,bool symetrized, prec *outcome, FILE* file, int where_to);

//run test for each differential operator for each precision
void test1(int dimx, int dimy, FILE* file, int where_to);

//writes all important physical/computational parameters int file
//void info(FILE* file, int where_to, int eig, int N, int dimstep);

//prints scalar products of eigenvectors
//size = r1.Nw and not r1.sets*r1.Nw
void ONO_table(ret_from_arpack_zn* vysl,int size, int eig, FILE* file, int where_to);

//compares the output and timings from APRACK and LAPACK for random matrices diag
void check_diag(FILE* file);            

//!*********** ph_func.cpp **************//

//gives transition rate of found spectrum between two levels - from(initial) and to(final)
//interaction: 0 - deformation [taking into account both longitudinal (def), and shear (def2) phonons]
//interaction: 1 - piezoelectric longitudinal
//interaction: 2 - piezoelectric transversal
//interaction: 3 - deformation transversal [shear (def2)]
prec transition_rate(state_info& states,int from, int to, int interaction);

//result =  8 transition rates for each state up to eigmax
//4 - 3 types of rates to the ground state + their sum
//4 - sum of 3 types to all other states (including ground) + their sum
void transition_rates(state_info& states);

//'spin' and 'charge' rate for the lowests excited state 
void spin_and_charge_rate(ret_from_arpack_zn& vysl, reg1& r1, int eig, int which, prec* spin, prec* charge);

//returns the spectral density, that is a sum over phonons (a, lambda) of the product of 
    //overlap xy:  <i|exp(\ii q_a.r)|j> <k|exp(-\ii q_a.r)|l>(=<l|exp(\ii q_a.r)|k>^*)
    //overlap z 
    //occupation: n^plus(omega_\alpha)
    //phonon couplings |M_a|^2
    //where a=(interaction, q), interaction=0,1,2,3 for def, piez long, piez trans, def trans (Si)
//energy is in meV
content spectral_density(state_info& states, int i, int j, int k, int l, int plus, int interaction, prec energy);

//two electron relaxation rate
//input: length of the index array
//array of indexes (refering to unsorted states)
//array of complex coefficients (of length at least length*length)
//energy difference of the states (in internal energy units)
//result - allocated outside, array of three real numbers (double)
//void two_el_rate(state_info& states, int length, int* indexes, content* coefs, prec energy, prec* result);


//!********************************dens_mat.cpp

struct rates; //info needed to find absorption

//find needed info 
struct rates *fill_rates(state_info& states, int n);

//find needed info dependent on omega
void find_resonant_states(struct rates *r, prec omega, int pinned);

//build a matrix (LHS) and vector (RHS of a set of (n-1) linear eq.)
void matrix_eq(rates *r, prec *Ares, prec *bres);

//solve these equations
void solve(int n, prec *A, prec *B, prec *res);

//give absorption
prec absorption(struct rates *r, prec omega);

//by half-interval method finds a frequancy with a given absorption (val)
//it is supposed that it is nearby resonance of states r->a and r->b and it must be lower as the resonant frequancy of these states
//the relative precision to use is eps
prec find_value(struct rates *r, prec val, prec step, prec eps);

//find larger and lower freq. around resonance of states a,b, where absorp. is one half of the resonant absorption
//if (mode>0) scan around the resonant frequency (with step such that there are mode points in the boundaries given by the previous) until the absorption drops under given relative value
//left - relative precision in giving the half-value freq.
//right - the epsilon when to stop if mode>0
//on return these two contain the left and right half-vlaue frequencies 
void scan(struct rates *r, int a, int b, prec* scanrates, FILE* file, int where, int mode, prec eps);


//!*********** dens_mat2.cpp **************//

//data needed to compute the derivative
struct derivs;

//computes the vector of derivatives for the components of the density matrix according to data from derivs, time, omega
// input and output have dimension n^2
void derivative(derivs *data, prec t, prec omega, content *input, content *output, int wherelog);

//computes time evolution of the density matrix components up to time tmax, or if the change of the diagonals in one step is lower than 1e-10
//all elements are written into file at each time step
void find_equilibrium(derivs *data, prec omega, content *result, prec tmax, FILE * file, int wherelog);

//one step ahead according to runge kuta method
void rk_step(derivs *data,content *values, prec omega, prec time, prec step, int wherelog);

//gives absorption of energy from the perturbing field
prec absorption(derivs *data, content *rho);

//fill information into derivs
derivs *fill_derivs(ret_from_arpack_zn& vysl, reg1& r1, int n);

//fill frequencies in the rotating wave approx.
//method==0 => considering all freq.
//method==1 => considering only the lowest freq.
void fill_derivs2(struct derivs *data, prec omega, int method, int where);

//!******************two electron************************//

//copy first table to the second. Starting in the lower left conrer, padding with zeros if the second table is larger 
void expandtable(content* from,int Nx, int Ny,content* to ,int Nxa,int Nya);

//computes a Coulomb element <ij| C | kl> using Fourier transform
//blowup - how many times to enlarge the space grid before Fourier is executed
//init = 1 initialize everything
//       0 use previously initialized FFTW and coefficients
//      -1 mop&lobby
content CoulombElement(state_info &I, int i, int j, int k, int l, int blowup, int initialize);

//calibrates C.e. by enlarging consecutively the grid and stopping after mistake(now) (=differences of two consecutive results; relative) is smaller than mistake(-1) and smaller than the precision 
//max - max number of steps
//result is in blowup
void CoulombElementCalibrate(state_info &I, int i, int j, int k, int l, prec precision, int& blowup, int max);

//diagonalize dense complex (in general non-Hermitian) matrix using arpack
//N - dimension of the matrix
//Neig - number of eigenvalues to return
//vecs, vals - output, allocated outside, if vecs=0, eigenvecsots are not returned
//arpack_type = 'S/L'+'M/R/I' smallest/largest+magnitude/real/imaginary 
void diagonalize_arp(int N, int Neig, content *H, content *vecs, content * vals, const char* arpack_type, progress_bar_prototype* progress_bar=0);

//check hermitivity/orthonormality/eigenality of matrix(es). type == 0 means bra corresponds to conjugation, type !=0 bra corresponds to transposition
prec is_hermitian(int N, content* H, bool print=false);
prec is_orthonormal(int N, int Neig, content* vecL, content* vecR=0, bool print=false);
prec is_eigensys(int N, int Neig, content* H, content* val, content* vecL, content* vecR=0, bool print=false);


//!******************transport************************//

//invert matrix of dimension n by n
void invert(content* res,content* in,int n);

void clean_geometry();
void clean_with_E();
void prepare_with_E();
void prepare_geometry();
void leads_geometry_init();
void leads_with_E_init();
void compute_self();

//compute the green function between the first and the last xline using recursive method
void recursive_green(bool full);

//tranmission amplitude from lead, mode, spin into a second combination of these three
content amplitude(int lead_from, int mode_from, int lead_to, int mode_to);

//total tranmission probability from lead to lead
prec lead2lead(int from, int to);

//transmission amplitudes through a series of rectangular barriers
void rectangular_barriers(content &t, content &r, prec Ef, prec Eb, prec w, prec wp, prec l);

//transmission/reflection amplitudes for one dimensional case
void oneDamplitude(content &t, content &r);

prec trace(int from, int to);

void draw_potential(FILE* file, int where);

//invert general complex matrix n by n
void invert(content* res,content* in,int n);

//diagonalize hermitian matrix in of dimension n by n, return vectors in w, energies in e
void diag(content *in, content *w, prec* e, int n);

//checks the algorithm by comparing analytical and numerical tramissions through
//a series of step functions -- needs modification of potential in hamilt.cpp
void check_GF(FILE* file, bool full);

//fills in the scattering matrix, returns pointer to arrays (matrices) s_ij for ij=leads
//and fills array with number of channels
content** scattering_matrix(int* M);

//!******************tunneling************************//

//compute the tunnel rates for the aligned_state aligned to the fermi energy (dot is raised by an offset) in all propagating channels. The rates for channels are returned in the array allocated outside
//the sum of all channel rates is returned
double tunneling_rate(state_info* states, int aligned_state, int unsorted, int& channel, double* channel_rate, double fit2value);


//!******************laser************************//
//D(x) / D_0; input in internal lengths
content pump_envelope(prec x, prec y);
//n(x) / n_0; input in internal lengths - this function is called only once
content refractive_index_envelope(prec x, prec y);

void set_matrix(int N, content *H);

void CF_operate(int N, void *A, content *x, content *y);


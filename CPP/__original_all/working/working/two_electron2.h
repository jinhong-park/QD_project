//#include "main.h"

struct symmetry_prototype {
  class state_info* states;
  int type;			//0/1/2/3 - No/I/Ix+Iy/Lz
  int min, max;			//interval for 2e symmetry ids
  prec tolerance;		//tolerance for symmetry sorting
				//returns tue if symmetry is below tolerance, 1e symmetry id for unsorted state
  bool symmetry1e(int u, int& id);
				//returns 2e symmetry id for k state
  int symmetry2e(int u1, int u2, int p);
				//converts symmetry id into the orbital quantum number (trashes p)
  int symmetry2einv(int id);

  void initialize(state_info* states_, int type_);

				//conversions between id and {Ix, Iy} for type 2
  void id2IxIy(int id, int& Ix, int& Iy);
  void IxIy2id(int Ix, int Iy, int& id);
};


struct selected_prototype {
  int u;			//unsorted label
  int valley;			//valley
};

struct statek_prototype {
  int i,j;			//selected labels
  int p;			//exchange symmetry
  int id;			//symmetry identificator
  int g;			//to which group it belongs
  prec norm;			//normalization: |k> = norm * ( |i>|j>+p|j>|i> ) 
  prec energy;			//energy (in internal units)
};

struct group_prototype {
  int id;			//2e symmetry identificator
  int N;			//number of k-states in this group
  int Neig;			//number of eigenstates asked for
  int* k;			//labels of k-states in this group
  content* H;			//Hamiltonian [N*N]
  content* eigenvecs;		//eignevectors [N*N]
  content* eigenvals;		//eigenvalues [N]

  void initialize(int N);	//allocate space
  void release();		//deallocate space
};

struct statea_prototype {
  int g;			//which group it is from
  int ag;			//which eigenstate from the group
  content* aij;			//coefficient in the |i x j> basis (selected labels)
};

struct stateb_prototype {
  int a;			//a-state
  int Sigma;			//Sigma = 0, 1, 2, 3 ... S, T+, T0, T-
  prec energy;			//orbital plus Zeeman (2 mu B S_z)
};

struct statef_prototype {
  content* cfb;			//coefficients in the b-states basis
  prec energy;			//energy
  content* cfkS;		//coefficients in the k-Sigma-state basis 
  prec S2;			//mean value of S2 operator
  prec SB;			//mean value for S \times B operator
  int g_max;			//group with maximal probability
  prec p_max;			//the weight of the maximal group (max of prob_g)
  prec *prob_g;			//probabilities for all groups (symmetries)
  int b;			//b-state to which this f-state is most similar
};

struct mat_ele_prototype {	//matrices for all pairs of selected states:
  int derivs_length;		//length of the "derivs" (times Ns*Ns entries) 
  content* derivs;		//derivatives: kx, ky, kxkxky, kxkyky, x, y, x ky - y kx-> 7*Ns*Ns
  content* orbital;		//vector to multiply Pauli matrix vector -> 3*Ns*Ns (selected states)
  content* orbital2e;		//vector to multiply Pauli matrix vector -> 3*Na*Na (a-states)
  content* spinor;		//Pauli matrix vector -> 2*3*4*4
  content* spin_imp;		//spin impurities (selected states)
  content* spin_imp2e;		//spin impurities (a-states)
};

struct timing_prototype {	//times the diagonalization took [ms]
  prec c_elems;			//coulomb elements
  prec final_diag;		//final diagonalization
  prec final_matrix;		//final matrix construction
};

class two_electron_state_prototype {

  state_info * states;		//single electron states

  int Ns;			//number of selected-states (all spin up)

  selected_prototype* selected;	//selected states

  symmetry_prototype symmetry;	//information about the symmetry

  int Nk;			//number of k-states
			
  statek_prototype* statek;	//k-states

  int Ng;			//number of symmetry groups

  group_prototype* group;	//symmetry groups
  CE_prototype* CE;       //Coulomb element class

  int Na;			//number of a-states = sum_g group[g].Neig

  statea_prototype* statea;	//a-states

  sorting_class* sorting;	//class which sorts the a-states

  int Nb;			//number of b-states

  stateb_prototype* stateb;	//b-states

  content* Hb;			//the final Hamiltonian to diagonalize

  int Nf;			//number of f-states

  statef_prototype* statef;	//f-states
  content* statef_vecs_begin;	//beginning of the eigenvecstor array

  mat_ele_prototype mat_ele;	//single particle matrix elements

  prec exchange_un;		//unperturbed exchange

				//progress bar
  progress_bar_prototype* progress_bar;
  
  content* rho;                 //reduced (single electron) density matrix

  //auxilliary for distribution characteristics - find the point of maximum of rho by interpolation
  //rho is the density, guess is the initial guess (sequential label), steps is the number of tries
  double MC_max_search(content* rho, int guess, int steps);
  
public:
  
  content M[10];                 //charge distributon moments
  double Rmax, Qsph, RSmax;      //position of the maximal density, distribution sphericity (0=perfect circle, infinity=line), RSmax position of the minimal spin density
  void distribution_characteristics(int I);//calculate the above from rho
  
				//constructor
  two_electron_state_prototype(class state_info* states_, sorting_class* sorting_);
				//destructor
  ~two_electron_state_prototype();


  //run the diagonalization (skipping first "step" steps)
  void diagonalize(prec par, int *sorted2un, int *unsorted2s, int step=0);

  //phonon induced relaxation f1 -> f2
  //result is an outside allocated vector of length 3 (def, pz l, pz t)
  void relaxation(int f1, int f2, prec* result);

  //print the energies and symmetry info for states:
  //st="b" / "f" / "f-b"
  //file, where: 1 - stdout, 2 - file, 4 - log
  //what - 1 energies - 2 symmetries - 4 heading (mean also a new line, otherwise not)
  //maximal number of states to show
  void statistics(const char* st, int what, FILE* file, int where, int Nmax=-1);

  //expectation value of "tag": of a state number s from group st
  //so far: st=b/f/f-b, tag=energy [meV]
  prec read(char st, int s, const char* tag);
  
  //matrix element of "tag" (single electron operator with index "particle"): 
  //between states number s1, s2 from group st
  //so far: st=b (spinor part is ignored), tag=dipolex, dipoley [nm] 
  content read(char st, int s1, int s2, int particle, const char *tag);

  struct timing_prototype timing;
  
  //computes the reduced density matrix rho=<f1 | f2>_1_orb where the overlap is taken only over the first particle orbital dof
  //the result is rho = rho_ij_Sigma1Sigma2 |i><j| |Sigma1><Sigma2|
  //where i,j = 0,..., Ns-1 and Sigma = 0,1,2,3
  //idexing is rho_ij_Sigma1Sigma2 = rho[(i*Ns+j)*16+Sigma1*4+Sigma2]
  //rho allocated inside
  void reducedDM(int f1, int f2);

  //computes the matrix element of <f1 | xx+yy | f2>
  content rdotr(int f1, int f2);
  
  //return the expectation value of sx,sy,sz at grid point with seq label n
  //that is 2 Tr[rho sigma delta(r-Rn) ]
  void spinatn(prec *s, int n);

private:
				//fill the selected-state array (type of symmetry)
  void construct_selected();


  void construct_statek();	//fill the k-state array

				//number of states with energy less than "energy" in a group
				//where = 0/1 - before/after the diagonalization
  int less_than_E_in_group(prec energy, int g, int where);

				//fill the groups, asking for Neigtot states from the diagonalization (if Neigtot=0, use all k-states)
  void construct_groups(int Neigtot=0);

  void diagonalize_groups();	//build hamiltonians, and diagonalize them, give time estimates and progress

  void construct_statea(prec par, int* sorted2un, int *unsorted2s);	//fill in the a-states, par for sorting, sorted2un and unsorted2s are from previous runs or initialized

  void construct_stateb();	//fill in the a-states

  void construct_statef();	//diagonalize the final matrix, fill the f-state array

  void mat_ele_fill_derivs();	//fills the structure mat_ele
  void mat_ele_fill_orbital();  
  void mat_ele_convert(content *me1e, content *me2e, const char* mess); //converts single electron into 2e matrix elements using a-state->aij
  void mat_ele_fill_spinor();
  void mat_ele_fill_spin_imp();
  content mat_ele_from1e(int b1, int b2, content* m1e);//returns matrix element of SOC between states b1 and b2 using a-state->k
  content mat_ele_from2e(int b1, int b2, content* m2e);//returns matrix element of SOC between states b1 and b2 using a-state->aij
  content mat_ele_spin_imp(int b1, int b2);//matrix element of spin impurities between b-states <b1 b2>

};

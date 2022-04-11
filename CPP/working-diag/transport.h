
struct transport_class {
  int b;                //# of nearest neighbours
  int l;                //number of leads (only for two)
  prec hx,hy;           //grid step
  reg1* r_w;            //region with the leads included

                        //all xlines refer to the region without the leads - that is, the "inside"
                        //except of the following two numbers:
  int first, last;      //first left and last right xline of r_w which are inside r_wo
                        //that is, the rest are in the leads
     //xlines are labeled starting from the first active line to region.xlines (-1) being the last one

  reg1* r_wo;           //region without the leads
  int* dim;             //number of active points on xlines
  int xlines;           //number of xlines with active points
                      
  //wavefunction is indexed by (xline, sequential label, spin) 
  //operators (Hamiltonian,...) are then indexed by two such (multi)indexes
  // and stored as a matrix of matrices:

  //xxx(i,j)= matrix for points on xline(i) \times xline(j), that is
  //[i*xlines+j] is pointer to matrix[dim(i)*2*dim(j)*2]
  content **self_array;       //self energy
  content **g_array,**g_aux_array;  //green functions
  
  class lead_prototype lead[2];         //leads, only two (L=0 and R=1) implemented 
  struct s_matrix_prototype s_matrix;//the s-matrix

  //converts two multiindexes into a label into a matrix stored as a row of numbers (vector)
  int label(int i, int j, int yi, int yj, int si, int sj);

  //returns label in the matrixes space
  int labelM(int i, int j);

  //functions which return the pointer to the greens function / GF aux / self energy matrix corresponding to xlines i,j
  //and the values of this matrix specified by the y and spin coordinates
  content* g(int i, int j);
  content g(int i, int j, int yi, int yj, int si, int sj);
  content* g_aux(int i, int j);
  content* self(int i, int j);
  content self(int i, int j, int yi, int yj, int si, int sj); 

  //true if both the hamiltonian and self energy are zero for these two xline indexes
  bool iszeroh(int i, int j);
  //true if the Hamiltonian is zero
  bool iszeroH(reg1* r, int i, int j, int yi, int yj) ;
  //true if both the hamiltonian and self energy are zero for these two indexes
  bool iszerov(int i, int j, int yi, int yj);

  //matrix element of the Hamitlonian plus the self energy
  content HpS(int i, int j, int yi, int yj, int si, int sj);

  //fills the matrix representing the E-H-Sigma of the region without leads part (i,j)
  void fill_EmHmS(content* res, int i, int j, prec Ef);

  //multiplication of two matrixes with given multiindexes
  void MM(content* res,content* a, content* b,int i1,int j1, int i2, int j2);
  //add two matrixes res=res+in with multiindexes i,j
  void add(content* res, content* in, int i, int j, prec coef);
//add a constant to a diagonal matrix with multiindexes i,j res=res+x*id 
void addC(content* res, content x, int i, int j);
//compute norm of a matrix in
prec normM(content* in, int i, int j);

void check_invert(content* a, content *inv, int i);

//construct the region including leads
void region_w();
//construct the region excluding leads
void region_wo();
//to do once - the geometry dependent (energy independent) tasks
void prepare_geometry();
//the energy dependent tasks
void prepare_with_E();
//deallocate space from geometry dependent tasks
void clean_geometry();
//deallocate space from energy dependent tasks
void clean_with_E();
//returns the longitudinal part of the green function of a semiinfinite lead
content leadg_long(int i, int j, prec k, prec v);
//geometry dependent tasks for leads
void leads_geometry_init();
//energy dependent tasks for leads
void leads_with_E_init();
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
//checks that the greens function is the inverse of E-H-Sigma cut after line i. Computes the norm of
//[G (E-H-S)]_{jk} - delta(j,l) for each j,k from 0 to i (including)
//result is the list of residual norms into the log
void check_it(int i) ;
//print the upper left entry of the three matrixes
void print_matrixes();

void cycle(int i, int &j); 
//compute the green function between the first and the last xline using recursive method
//full -- compute the full inverse matrix (true) or just the four corners (false)
void recursive_green(bool full);
//returns y label of a point in the lead neighbouring with point y on xline inside
//if no, returns -1 - assumes no holes in the lead
int jcor(int xline, int y);
//tranmission amplitude from lead, mode, spin into a second combination of these three
content amplitude(int lead_from, int mode_from, int lead_to, int mode_to);

prec lead2lead(int from, int to);

prec trace(int from, int to);

void draw_potential(FILE* file, int where);

void check_GF(FILE* file, bool full);

//extract transmission and reflection amplitudes for a one dimensional case (all xlines have length (dim) 1)
//(note: valid only for spin along z)
void oneDamplitude(content &t, content &r);

//fills in the scattering matrix
content** scattering_matrix(int* M);


};

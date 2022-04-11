#include <math.h>
#include <string.h>
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <termios.h>
#include <unistd.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_bessel.h>

#define F77NAME(x) x ## _

using namespace std;

#define where_to_print 4
#define where_to_print_critical 1+4

#ifndef _AUXILS_H

#define _AUXILS_H

extern char buffer[1000];
extern FILE* logg;

 /*vypis retazec v poli co podla bitov v premennej kam

bit 1 na obrazovku
bit 2 do suboru f
bit 3 do suboru logg (definovany v vlastne.cpp
*/
void splachni(FILE* f,const char* co,int kam);

//simple string to the where_to_print_critical
void message(const char* message);

//simple string to the where_to (file is logg)
void message(const char* message, int where_to);

//simple string to where you want
void message(FILE* file,const char* message, int where_to);

//closes logg file and opens it under a given name
void other_name(const char * name);

//clears content of a file
void clear_file(const char * name);

class output_class {
  int NF;
  char** filenamebase;
  char** filename;
  FILE** files;
  
  char* logfilenamebase;
  char* logfilename;

public:  
  
  output_class(int NF_, const char* path, const char** name, const char* logname);
  ~output_class();
  void clear_files(int which=-1);
  void open_files();  
  void close_files(); //except the log
  void set_appendix(const char* tag, int p=-1);
  FILE* file(int i);
};

//#endif

#define max(x,y) ( (x>y) ? x : y )
#define min(x,y) ( (x>y) ? y : x )

/*gama funkcia s (polo)ciselnym argumentom
gama(i,j)=gama funkcia argumentu (i/j); j={1,2}*/
double gama(int i, int j);

//erfi
double erfi(double x);

//inverse of the error function: p=erf[ erfinv(p) ] 
double erfinv(double p);

//find first n solutions of equation J_l(x)=J_{l+1}(x)
void bessel_solver(int l, int n, double* res);

//factorial
double fac( int n);

//power with integer exponent
double pown(double x, int n);
 
//converts integer into ascii
void itoa(int i,char *num2);

//keeps x in the region <0, n-1>
inline int pos_mod(int x, int n)
{
  if (n==0) return(x);
  int aux=x % n;
  if (aux<0) aux+=n;
  return(aux);
}

//x1 and x2 (smaller/larger) are the solutions of the quadratic equation, returns 0 if root(s) exist, -1 if not
int quadratic(double p, double q, double&x1, double& x2);

//root of linear equation with complex coefficients c0 x + c1 x^* + c2 =0
//info = -1/0/1 for no/single/infinity solutions
int complex_bilinear(complex<double>* coef, complex<double>& x);

//finds roots of quadratic equation with complex coefficients c0 |x|^2 + c1 x + c2 x^* +c3 =0
//returns -1/0/1 for no/two/infinity solutions
int complex_quadratic(complex<double>* coefs, complex<double>& x1, complex<double>& x2);

//finds roots of cubic equation: y = sum_i={0..3} c[i] x^(3-i)
//returns the number of real roots (root[0..N-1] are real, the rest are complex)
int cubic(double* coef, complex<double>* root);

//mod with result both possitive and negative values: < or ( -n/2,n/2> 
//depending whether n is odd or even
int mod_pm(int i,int n);

//the nearest integer number
int my_round(double x);

//sort according to x in 'a'/'d'e scending order
void picsrt_inplace(int N, double* x, char order='d');
void picsrt_inplace(int N, int* x, char order='d');

//complex x: sort according to abs/real/imag('a'/'r'/'i'/'R'/'I') part of x [or absolute value of]
void picsrt_inplace(int N, complex<double>* x, char method='a', char order='d');

//sort according to x in 'a'/'d'e scending order
//create the index and rank arrays, the input array is not touched
void picsrt_index(int N, double* x, int * index, int *rank=0, char order='d');
void picsrt_index(int N, complex<double>* x, int * index, int *rank=0, char method='a', char order='d');

//sort the array x according to an index table index
void sort_from_index(int N, double * x, int * index);

//returns x coth(x) as (1+u)/(1+d)
void xcothx(double x, double &u, double &d);

//Brillouin function
double brillouin(int J, double x);



/*
//function fom complex numbers: norm(x)=x x*
double norm(const complex<double> & x);

//function fom complex numbers: abs(x)=sqrt(x x*)
double abs(const complex<double> & x);
*/

//returns a matrix representing a two dimensional differential operator
//dx^i dy^j using 2l+1 by 2l+1 points (including the central one)
//apart from the coefficient 1/ (hx^i hy^j)
void differential_operator(int i, int j, int l, double coef, double* out);

//returns coefficients defining the differential operator d^d/dx^d at the zero coordinate
//input: d (previous - the derivative), n = number of available points (=2b+1 previously)
//vector of coordinates of the available points (in internal length units)
void construct_differential_operator(int d, int n, double* coors, double* coefs);

class stopky {

 struct prvok {
   timeb cas;
   prvok *dalsi;
 };

  int bezi;
  prvok *prvy, *posledny;
public:
  stopky() {bezi=0; prvy=posledny=0;};
  ~stopky() {while (bezi) zastav();};
  void pusti();
  long int doba(int ktora=0);
  long int zastav (int ktora=0);
};


//precitaj znak z nebuferovanej klavesnice
int u_getch(void);

//stop executing - wait for pressing a key, ig it is 'esc' exit, otherwise continue
int stop();



/*generuj nahodne cislo z intervalu <0,1> */
double generuj();

//naplni generator aktualnou milisekundou, pripadne explicitnou hodnotou
int randomize(int seed=-1);

double ran1(int i=0);


class progress_bar_prototype {
  char* title;			//tag to begin the line - string, ends by \0

  double todo;			//goal in some units
  double done;			//units already done
  long int last_update;		//time of the last update (milli seconds)
  double eta;			//estimated time remaining (ms)

  char *mess;			//text (to be) printed on the screen
  int length;			//length of the message printed (if 0 - not out yet)

  char* aux, *aux2;		//string buffer
  stopky* timer;		//timer

  int file_desc_stdin;		//file descriptor of the standard output
  struct termios preserve, temp;//console atributes
  bool terminal;

public:
				//constructor, destructor
  progress_bar_prototype(const char* title_);
  ~progress_bar_prototype();

  void reset(const char* title_);//change the name, zero counters

  void start(bool on=true);			//start the timer

				//recomputes eta, deletes old message, if on, prints actual
  void print_message(bool on=true);	
				//update points todo, or done
  void add(double todo, double done, bool print=true);
				//change ETA to total time, print, put newline
  long int finished(bool on=true);

private:
				//converts integer to time string (h:)mm:ss, 
				//output into outside allocated res, terminated by \0
  void time2string(long int time, char* res); 
  

};


//check several nr routines
void nr_check1();

//polynomial interpolation of N points of data y at coordinates x
//value at x0 is computed by constructing a polynomial of order N-1 and the absolute error is estimated
complex<double> polint(double* x, complex<double>* y,  int N, double x0, double& err);
void polint_check(int N);
double trapzd( complex<double> f(double) , double a, double b, complex<double>& s, int n);
double midpnt( complex<double> f(double) , double a, double b, complex<double>& s, int n);
void trapzd_check(complex<double> f(double), complex<double> F(double));
complex<double> qromb(complex<double> f(double), double a, double b, int& maxsteps, int order, double& err);
complex<double> qromo(complex<double> f(double), double a, double b, int& maxsteps, int order, double& err);



/*class banded_matrix_prototype {


  int N;			//length (matrix linear dimension)
  int Nl,Nu;			//number of sub (lower) and super (upper) diagonals
  int* pntr;			//number of non-zero elements on the corresponding diagonal

  complex<double>** diagonal;   //diagonals stored one by one:
				//the matrix element (i,j) is on the (i-j)-th diagonal, with the sign encoded as the least significant bit, so that 
				//sub diagonals -1 to -(N-1) go to 1, 3, ..., 2(N-1)-1 and
				//super diagonals 1 to (N-1) go to 2, 4, ..., 2(N-1)
				//the main diagonal can be considered as super, 0 goes to 0
				//the sequential label of point (i,j)
				//if above the main diagonal (i<j), its seq label is max(i,j)=j
				//if below the main diagonal (i>j), its seq label is min(i,j)=j
  void cart2diag(int i, int j, int &d, int &seq);//conversion from cartesian into diagonal coordinates
  void diag2cart(int d, int seq, int &i, int &j);//conversion from cartesian into diagonal coordinates

  public:  
  complex<double>* AB;		//banded matrix stored in LAPACK form
  int KL, KU;			//number of diagonals of the AB matrix
  int LDAB;			//dimension of the banded matrix (=2*KL+KU+1)
  int *IPIV;			//auxiliary array
    
  void add(int i, int j, complex<double> x);//add the value, increase pointers if non-zero, possibly allocate the new diagonal
  complex<double> read(int i, int j);
  
  banded_matrix_prototype(int N);//constructor
  ~banded_matrix_prototype();//destructor
  
  void matrix_info(bool out=true);//into log
  void matrix2lapack();
  void LU_factorize();
  void LU_solve(complex<double>* x, complex<double> *y); //solve M x = y, do not overwrite y
  
};*/

class sparse_matrix_prototype {

public:
  int N;			//matrix linear dimension
	//the sparse storage: 
		//each line is stored as an array of labels (j coordinates) and values, so that the matrix elements are (the rest are zero)
		//value[i][j] = M(i,label[j]) for j=0,...,length[i]-1
  
  int max_length;		//should be set to expected number of elements in a line
  int* length;			//length[i] is the number of nonzero entries on line i
  int**label;			//label[i] is an array of labels of non zero entries
  complex<double>**value;	//value[i] is an array of values of non zero entries
private:	
	//the banded storage
			
		//the matrix element (i,j) is on the (i-j)-th diagonal, with the sign encoded as the least significant bit, so that 
		//sub diagonals -1 to -(N-1) go to 1, 3, ..., 2(N-1)-1 and
		//super diagonals 1 to (N-1) go to 2, 4, ..., 2(N-1)
		//the main diagonal can be considered as super, 0 goes to 0
		//the sequential label of point (i,j)
		//if above the main diagonal (i<j), its seq label is max(i,j)=j
		//if below the main diagonal (i>j), its seq label is min(i,j)=j
  inline void cart2diag(int i, int j, int &d, int &seq)//conversion from cartesian into diagonal coordinates
  { d=i-j;
    if (d<0) d=-1-2*d; else d*=2;
    seq=j; };
  inline void diag2cart(int d, int seq, int &i, int &j)//conversion from cartesian into diagonal coordinates
  { j=seq;
    if (d%2==1) d=(-1-d)/2; else d/=2;
    i=d+j;}
  
  int Nl,Nu;			//number of sub (lower) and super (upper) diagonals
  int* diagonal;		//number of non-zero elements on the corresponding diagonal
 
 public:
     
  complex<double>* AB;		//banded matrix stored in LAPACK form
  int KL, KU;			//number of diagonals of the AB matrix
  int LDAB;			//dimension of the banded matrix (=2*KL+KU+1)
  int *IPIV;			//auxiliary array of length N
			


  sparse_matrix_prototype(int N, int max_length);//constructor
  ~sparse_matrix_prototype();//destructor
  
  void add(int i, int j, complex<double> x);
  complex<double> read(int i, int j);
  void operate(complex<double>*x, complex<double>*y); //y = M x
  void matrix_structure(bool out=true);//into log
  void matrix2banded();
  void LU_banded_factorize();
  void LU_banded_solve(complex<double>* x, complex<double> *y); //solve M x = y, do not overwrite y
};

//void DeNoise(double delta, double v, double l, double t, int length, double* data, double &average, double& min, double& max, double &rms, double & maxd, FILE* log=0);
void RealFourierTransform(double *data, int N, int isign);
void FourierTransform(double *data, int N, int isign);

#endif

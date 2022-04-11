// replace all "z"  by "c" ("z" for double) and save under new name!!!
//and change typedef precisionon line 36 !!!

/*demand for dynamical memory allocation:

			dim*(ncv+nev+5)		of precision

			3*ncv*ncv+9*ncv+2*nev+2	of precision
			ncv+29			of int
*/

#ifndef _arpack_zn_H

#define _arpack_zn_H

//using namespace std;

#include "auxiliaries.h"

struct ret_from_arpack_zn;

class arpack_zn {

public:
typedef double precision;
typedef int integer;
typedef int logical;
typedef complex<precision> content;
typedef void (*pfAxy)(int ,void*, content*, content *);

//parameters for all the ARPACK top level routines
private:
integer ido;		//communication flag
char bmat;		//standard('I') or generalized ('G')
char *which;		//what kind of eigenvalues to return (f.e. "SM" = smallest)
integer nev;		//# of vectors asked
precision tol;		//relative accuracy
integer ncv;		//# of Lanczos vectors used
integer* iparam;	//lotsa
int max_iter; //max. numer of iterations

//parameters for snaupd

/*
iparam[0]=1;		//exact shifts
iparam[2]=999;		//max # of iterations
iparam[6]=1;		//mode 1 - A x = a x
*/

integer n;		//dimension of the problem = N
content *resid;		//vector with lenght N
content *V;		// array N times ncv
integer ldv;		//leading dimension of V
integer *ipntr;		//lotsa

content *workd;		//array 3*N
integer lworkl;		//=3ncv*ncv+5*ncv
content *workl;		//array of dim lworkl
precision *rwork;	//array ncv
integer info;		//communication flag


//parameters for sneupd

logical rvec;		//true = eigenvector also
char HowMny;		//'A' =all (=nev)     'S' -some, specified by select
logical* selekt;	// array of dim ncv : select[i] = true - eigenvector i wanted
content *d;		//array of dim nev+1 - this contains eigenvalues
content *z;		//array of dim N times nev+1 (can be set to V)
integer ldz;		//The leading dimension of the array Z. In any case,  LDZ .ge. 1.
content sigma;		//shift (not significant for mode 1)
content* workev;		//workspace of dim 2*ncv
// all others must be the same as in previous call to snaupd

class progress_bar_prototype* progress_bar;
bool with_progress_bar;

//checks if exit from ..upd was with zero error (i) after iter iterations, exit if input error, returns "ok"
bool make_sure(int i, int iter);

//free all allocated memory for arrays
void skonci();

  public:
//Prints to logfile (to where) the matrix <x_i A x_j> / d[j] for each of "for_hmn" vectors.
//returns the #of vectors where residual norm is more than "tolerance",
//d is array of eigenvalues
//z is 2d array of eigenvectors
//pf is user function for matrix times vector
int check_eigenset(int for_hmn,pfAxy pf,int n,void *A, content* d, content *z, int where);


public:

/*allocate memory for all variables needed for ARPACK routines, iniciate, parameters:
dim - dimensionality
howmany_eval - how many eigenvalues to find
tolerance - tolerance for convergence (f.e. 1e-10)
max_iter - maximum of iteration inside of ..aupd
eigen_vec - are eigenvec needed?
which_eig - "SM"-smallest,...
explic_ncv - set ncv by hand (default 0 -means let it on constructor)
*/
arpack_zn(int dim, int howmany_eval,precision tolerance,int max_iter,bool eigen_vec, const char* which_eig,int explic_ncv=0);
~arpack_zn();

/*recursive interface + checking for correctness of eigenvectors (if asked for):
INPUT - n- dimensionality
        pf - user function for A times x
	A - (some) information on matrix, it is user only through pf
	check - should the routine check validity of eigenvectors?
OUTPUT - eigenvalues are in 1d array d
       - eigenvectors in 2d array z
       - casy -array of 4 times in milisecs:
           casy[0] - time spend in ..aupd
	   casy[1] - times spend in ..eupd
	   casy[2] - time spend in user routine y=Ax
	   casy[3] - time spend in checking eigenvectors
*/
bool go_for_it(int n, pfAxy pf, void *A,ret_from_arpack_zn *output,bool check=false);
void show_stats(ret_from_arpack_zn* vysl,FILE* file, int where_to, bool newline,bool header);
void set_parameter(int, void*);
void read_parameter(int ,void*);
void set_progress_bar(progress_bar_prototype* pb);

};

struct ret_from_arpack_zn {
  struct  {
    int main,after,user,checking;
          } times;
  struct  {
    int iter,opx,reort,converged;
          } perform;
  arpack_zn::content *eigenvecs;
  arpack_zn::content *eigenvals;
  arpack_zn::content *eigenvecs_plaq; 
  arpack_zn::content *eigenvecs_orig;
  

};

//!lapack eigensys

class lapack_zgee_prototype {
                                        //I-input, O-output, W-workspace
  char jobvs, jobvl, jobvr;             //I:    V/N vectors (schur, left eig, right eig) are/are not computed
  char sort;                            //I:    S/N eigenvalues are/are not sorted
  int N, lda;                           //I:    order and leading dimension of the matrix A
  complex<double>* A;                   //I/O:  the matrix to be diagonalized / Schur matrix T
  int sdim;                             //O:    number of selected eigenvalues
  //int (select)(complex<double> x);    //I:    user supplied function used to select the eigenvalues to be sorted
  complex<double> *w;                   //O:    eigenvalues
  complex<double> *vs, *vl, *vr;        //O:    unitary matrix of vectors (Schur Z, left eig, right eig)
  int ldvs, ldvl, ldvr;                 //I:    leading order of matrix Z, VL, VR
  complex<double> *work;                        //W/O:  array od dimension lwork work[0] returns optimal lwork
  int lwork;                            //I:    workspace dimension (-1 for query)
  double* rwork;                        //W:    workspace double [2N]
  int* bwork;                           //W:
  int info;                             //O:    0-success, -N-wrong parameter#N, +N-only N eigenvalues converged
  
  public:
  lapack_zgee_prototype(int N);        //constructor
  ~lapack_zgee_prototype();            //destructor
  int diagonalize(complex<double> *H, complex<double> *eigval, complex<double> *eigvecsL=0, complex<double>* eigenvecsR=0, bool check=false); //diagonalize (and preserve?) the matrix, right/left eigenvectors optional, eigensystem check optional
};

class lapack_zhee_prototype {
                                        //I-input, O-output, W-workspace
  char jobz;                            //I:    V/N eigen vectors are/are not computed
  char uplo;                            //I:    U/L upper/lower triangul of A is provided on input
  int N, lda;                           //I:    order and leading dimension of the matrix A
  complex<double> *A;                   //I/O:  the matrix to be diagonalized
  double *w;                            //O:    eigenvalues
  complex<double> *work;                //W/O:  array od dimension lwork work[0] returns optimal lwork
  int lwork;                            //I:    workspace dimension (-1 for query)
  double* rwork;                        //W:    workspace double [3N-2]
  int info;                             //O:    0-success, -N-wrong parameter#N, +N-only N eigenvalues converged
  
  public:
  lapack_zhee_prototype(int N);        //constructor
  ~lapack_zhee_prototype();            //destructor
  int diagonalize(complex<double> *H, double *eigval, complex<double> *eigvecs=0, bool check=false); //diagonalize (and preserve?) the matrix, eigenvectors optional, eigensystem check optional
};

class lapack_dstev_prototype {
  //I:	input, O-output, W-workspace
  int N;				//I:	the matrix to decompose is N-rows by N-columns
  char jobz;			//I:	whether eigenvectors are to be computed
  double *D;		    //I/O:	diagonal entries of the matrix A, replaced by eigenvalues on output
  double *E;			//I/O:	subdiagonal entries of matrix A, on output overwritten
  double *Z;			//O:	eigenvectors of A
  int ldz;			    //I:	leading dimension of array Z(normally N)
  double *work;         //W:	workspace max(1,2*N-2)
  int info;				//O:	success
  
public:
  lapack_dstev_prototype(int N);        //constructor
  ~lapack_dstev_prototype();			//destructor
  int diagonalize(double *d, double *e, double *eigval, double *eigvecs=0, bool check=false);
};

class lapack_zgesvd_prototype {
  //I:	input, O-output, W-workspace
  int N,M;				//I:	the matrix to decompose is M-rows by N-columns
  char jobu, jobvt;			//I:	whether the matrices U and V^T are to be computed
  complex<double> *U, *VT;		//O:	U is MxM, VT is NxN matrix
  complex<double> *A;			//I/O:	input matrix MxN, may be overwritten
  int ldu, ldvt, lda;			//I:	matrices leading dimensions (normally M, N, min(M,N))
  complex<double> *work;		//W/O:	workspace
  int lwork;				//I:	work dimension, -1 for query
  double *S;				//O:	singular values, length=min(M,N)
  double *rwork;			//W/O:	length 5*min(M,N)
  int info;				//O:	success
  
public:
  lapack_zgesvd_prototype(int M, int N);	//constructor
  ~lapack_zgesvd_prototype();			//destructor
  int decompose(complex<double> *a, double *sv, complex<double> *vt=0, complex<double>* u=0, bool check=false, int M=-1, int N=-1);
  //int decompose2(complex<double> *a, double *sv, complex<double> *vt=0, complex<double>* u=0, bool check=false);
};




//interface for fortran library
extern "C"
{

void F77NAME(znaupd) ( arpack_zn::integer *ido, char *bmat, arpack_zn::integer *n, char *which, arpack_zn::integer *nev,arpack_zn::precision *tol, arpack_zn::content *resid, arpack_zn::integer *ncv, arpack_zn::content *V, arpack_zn::integer *ldv, arpack_zn::integer *iparam, arpack_zn::integer *ipntr, arpack_zn::content *workd, arpack_zn::content *workl, arpack_zn::integer *lworkl, arpack_zn::precision *rwork, arpack_zn::integer *info);

void F77NAME(zneupd) (arpack_zn::logical *rvec, char *HowMny, arpack_zn::logical *select, arpack_zn::content *d, arpack_zn::content *z, arpack_zn::integer *ldz, arpack_zn::content *sigma, arpack_zn::content* workev, char *bmat, arpack_zn::integer *n, char *which, arpack_zn::integer *nev, arpack_zn::precision *tol, arpack_zn::content *resid, arpack_zn::integer *ncv, arpack_zn::content *V, arpack_zn::integer *ldv, arpack_zn::integer *iparam, arpack_zn::integer *ipntr, arpack_zn::content *workd, arpack_zn::content *workl, arpack_zn::integer *lworkl,arpack_zn::precision *rwork, arpack_zn::integer *info);

void F77NAME(zgees) (char* JOBVS, char* SORT, int* (SELECT)(complex<double>*), int* N, complex<double>* A, int* LDA, int* SDIM, complex<double>* W, complex<double>* VS, int* LDVS, complex<double>* WORK, int* LWORK, double* RWORK, int* BWORK, int* INFO );

void F77NAME(zgeev) (char* JOBVL, char* JOBVR, int* N, complex<double>* A, int* LDA, complex<double>* W, complex<double>* VL,  int* LDVL, complex<double>* VR, int* LDVR, complex<double>* WORK, int* LWORK, double* RWORK, int* INFO );

void F77NAME(dstev) (char* JOBZ, int* N, double * D, double * E, double * Z, int* LDZ, double* WORK, int* INFO );
  
void F77NAME(dgesv)(int*, int*, double*, int*, int*, double*, int*, int*); 

void F77NAME(zheev)(char* JOBZ, char* UPLO, int* N, complex<double>* A, int* LDA, double* W, complex<double>* WORK, int *LWORK, double* RWORK, int* INFO );

void F77NAME(zgesvd)(char* JOBU, char* JOBVT, int* M, int* N, complex<double>* A, int* LDA, double* S, complex<double>* U, int* LDU, complex<double>* VT, int* LDVT, complex<double>* WORK, int* LWORK, double* RWORK, int* INFO );

}

#endif




#include <iostream>
#include <fstream>

#include "main.h"


int main()
{
  const char* name[]={"noname"};
  output_class output(0,"data/",name,"log");
  output.clear_files();
  output.open_files();
  
  int N=3;
  double k=N/2.0, theta=0, delta=1, t=1;
  
  content* auxval=new content[N];
  content* auxvec=new content[N*N];
  content* H=new content[N*N];
  content* eigenvecs=new content[N*N];
  
  lapack_zgee_prototype lapack(N);

  
  H[0+0] = content(1,0);H[1+0*N] = content(0,1);H[2+0*N] = content(0,1);
  H[0+1*N] = content(0,1);H[1+1*N] = content(1,0);H[2+1*N] = content(0,2);
  H[0+2*N] = content(1,1);H[1+2*N] = content(1,2);H[2+2*N] = content(3,0);
  lapack.diagonalize(H, auxval, auxvec, eigenvecs, true);
  
  cout << "Hamiltonian" << endl;
  for(int i = 0; i < N; i++)
  {
    for(int j= 0; j < N; j++)
    {
      cout << H[j+i*N] << '\t';
    }
    cout << endl;
  }

  cout << "eigenvalue1 : " << auxval[0] << '\t' << "eigenvalue2 : " << auxval[1] << '\t' << "eigenvalue3 : " << auxval[2] << endl;
  cout << "eigenvector1 :" << endl;
  cout << eigenvecs[0 + 0*N] << '\t' << eigenvecs[1 + 0*N] << '\t' << eigenvecs[2 + 0*N] << endl;
  cout << "eigenvector2 :" << endl;
  cout << eigenvecs[0 + 1*N] << '\t' << eigenvecs[1 + 1*N] << '\t' << eigenvecs[2 + 1*N] << endl;
  cout << "eigenvector3 :" << endl;
  cout << eigenvecs[0 + 2*N] << '\t' << eigenvecs[1 + 2*N] << '\t' << eigenvecs[2 + 2*N] << endl;

  cout << "conjugate(eigenvector1) :" << endl;
  cout << conj(eigenvecs[0 + 0*N]) << '\t' << conj(eigenvecs[1 + 0*N]) << '\t' << conj(eigenvecs[2 + 0*N]) << endl;


  delete[] H;
  delete[] auxval;
  delete[] auxvec;
  delete[] eigenvecs;
  return 0;

}











#include <iostream>
#include <fstream>

#include "main.h"


int main()
{
  /*coding("ponorit sa, alebo radsej plavat?",1);
  coding("plavat ako lod, ktora sa klze po tvojej vonavej a vlhkej lesklej hladine.",2);
  coding("ale musime pritom uzavriet elektricky okruh, medzi styrmi rukami a styrmi ocami.",3);
  coding("len vtedy sa mozeme do seba spravne zacvaknut v podpalubi.",1);
  coding("len vtedy sa lod spravne knise a vzdychaju jej pritom vsetky plachty.",2);
  coding("napinaju a uvolnuju sa lana a vydavaju pritom zvuky, ktore sa podobaju na starodavnu piesen o modrych ociach.",3);
  coding("pribeh o tom ako sa lod vynara a zase ponara, stale tam a spat, von a dnu, do tvojho zvlneneho oceanu.",1);*/
  /*
  const char* name[]={"noname"};
  output_class output(0,"data/",name,"log");
  output.clear_files();
  output.open_files();
  
  int N=5;
  double k=N/2.0, theta=0, delta=1, t=1;
  
  lapack_zhee_prototype lapack(N);
  complex<double> *H = new complex<double>[N*N];
  double *val = new double[N];
  complex<double> *vec = new complex<double>[N*N];
  complex<double> ii=complex<double>(0,1);
  
  
  val[0] =9.0;
  H[0+0] = complex<double>(1,0);
  H[0+1] = complex<double>(0,1);
  H[1+1] = complex<double>(0,1);
  H[2+1] = complex<double>(0,1);
  H[2+2] = complex<double>(0,1);
  H[2+3] = complex<double>(0,1);
  H[2+4] = complex<double>(0,1);
  H[2+5] = complex<double>(0,1);
  vec[0] = H[0];
  vec[1] = H[1];
  vec[2] = H[2];
  vec[3] = H[3];
  
  cout << H[0].imag() << endl;

  lapack.diagonalize(H, val, vec, false);

  cout << val[0] << '\t' << val[1] << endl;
  cout << vec[0] << '\t' << vec[1] << '\t'<< vec[2] << '\t' << vec[3] << endl;

  delete[] H;
  delete[] val;
  delete[] vec;
*/

  const char* name[]={"noname"};
  output_class output(0,"data/",name,"log");
  output.clear_files();
  output.open_files();
  
  int N=2;
  double k=N/2.0, theta=0, delta=1, t=1;
  
  lapack_dstev_prototype lapack(N);
  double *d=new double[N];
  double *e=new double[N];
  double *vecs=new double[N*N];
  double *vals=new double[N];
  content* auxval=new content[N];
  content* auxvec=new content[N*N];
  content* H=new content[N*N];
  content* eigenvecs=new content[N*N];
  
  progress_bar_prototype progress_bar("lapack dstev");
  progress_bar.start(true);
    for (int i=0;i<N;i++) {d[i]=delta*cos(M_PI/N*k*i+theta); e[i]=t;}
    lapack.diagonalize(d, e, vals, vecs, true);
  for (int i=0;i<N;i++) printf("eigval nr %i is %e \n",i,vals[i]);
  progress_bar.finished(true);

  cout << vals[0] << '\t' << vals[1] << '\t' << vals[2] << '\t' << vals[3] << '\t' << vals[3] << endl;
  
  lapack_zgee_prototype lapack2(N);

  H[0] = content(1,0);
  H[1] = content(0,1);
  H[2] = content(0,1);
  H[3] = content(0,1);
  lapack2.diagonalize(H, auxval, auxvec, eigenvecs, false);

  //cout << H[0] << endl;
  cout << auxval[0] << '\t' << auxval[1] << endl;
  cout << auxvec[0] << '\t' << auxvec[1] << endl;
  cout << eigenvecs[0] << '\t' << eigenvecs[1] << endl;
  cout << eigenvecs[2] << '\t' << eigenvecs[3] << endl;


  delete[] d;
  delete[] e;
  delete[] vecs;
  delete[] vals;
  exit(0);


  //check_de_noise();
  //exit(0);
  //printf("hello\n");
  //dices();
  //give_fermi_contact_coupling();
  //give_res_gap_renorm();
  //exit(0);
  //check_ipr();
  //average(10); exit(0);
  //for (int i=0;i<20;i++) convert_for_xmgrace(i); exit(0);
  //normalize();
  
  //printf("units.length=%e\n",units.length);
  //!aux: units calculation
  //give_res_helical2();
  //exit(0);
  //return(0);


}









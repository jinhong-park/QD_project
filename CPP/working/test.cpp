#include "main.h"

prec kx,ky;
int pntr;

content func(prec x, prec y)
{
  prec valr=sin(kx*x)*sin(ky*y);
  prec vali=0.8*kx*x*x*x+0.7*ky*y*y*y+0.6*kx*x*x*y+0.5*ky*y*y*x+ 0.4*x*x+0.3*y*y+0.2*x*y+x+y+0.1;
  return(content(valr,vali));
}

//d
content der(prec x, prec y)
{
  prec valr[10],vali[10];//dx dy dd dxdy dxdx dydy dxdxdy dxdydy
                         //dxdxdx dydydy
  valr[0]=kx*cos(kx*x)*sin(ky*y);
  valr[1]=ky*sin(kx*x)*cos(ky*y);
  valr[3]=kx*ky*cos(kx*x)*cos(ky*y);
  valr[4]=-kx*kx*sin(kx*x)*sin(ky*y);
  valr[5]=-ky*ky*sin(kx*x)*sin(ky*y);
  valr[2]=valr[4]+valr[5];
  valr[6]=-kx*kx*ky*sin(kx*x)*cos(ky*y);
  valr[7]=-kx*ky*ky*cos(kx*x)*sin(ky*y);
  valr[8]=-kx*kx*kx*cos(kx*x)*sin(ky*y);
  valr[9]=-ky*ky*ky*sin(kx*x)*cos(ky*y);
  vali[0]=3*0.8*kx*x*x+2*0.6*kx*x*y+0.5*ky*y*y+2*0.4*x+0.2*y+1;
  vali[1]=3*0.7*ky*y*y+0.6*kx*x*x+2*0.5*ky*y*x+2*0.3*y+0.2*x+1;
  vali[3]=2*0.6*kx*x+2*0.5*ky*y+0.2;
  vali[4]=6*0.8*kx*x+2*0.6*kx*y+2*0.4;
  vali[5]=6*0.7*ky*y+2*0.5*ky*x+2*0.3;
  vali[2]=vali[4]+vali[5];
  vali[6]=2*0.6*kx;
  vali[7]=2*0.5*ky;
  vali[8]=6*0.8*kx;
  vali[9]=6*0.7*ky;
  return(content(valr[pntr],vali[pntr]));
}

void test2(int dimx, int dimy, reg1::operators op,bool symmetrized, prec *outcome, FILE* file, int where_to)
{
  prec vysl[15];
  static prec lengthx,lengthy;
  if (! symmetrized) {//so calling of test with symetrized after non-symetrized will
  			//use the same parameters
    randomize();
    lengthx=1.23*ran1()+1;
    lengthy=1.23*ran1()+1;
    kx=ran1()*1.5+0.75;
    ky=ran1()*1.35+1.1;
  }

  pntr=op;
  reg1 r1(dimx,-lengthx,2*lengthx,false,dimy,-lengthy,2*lengthy,false,inside,2);

    for (int j=reg1::two;j<reg1::all_prec;j++) {

      r1.act_prec_set((reg1::operators_precision) j);
      r1.symmetrize_ops(symmetrized);
      //r1.ant(potential,units.length/units.nm);
      r1.test(op,func,der,vysl);
      outcome[j]=vysl[4*3+0]/vysl[1*3+0];	//norma residuum/norma derivacia

      sprintf(buffer,"operator nr. %i with symmetrization %i measured with control values:\n",op,symmetrized);
      splachni(file,buffer,where_to);
      sprintf(buffer,"lx=%e ly=%e kx=%e ky=%e\ndx=%i dy=%i\n",lengthx, lengthy, kx, ky, dimx, dimy);
      splachni(file,buffer,where_to);
      sprintf(buffer,"at precision %i\n with results:\n",j);
      splachni(file,buffer,where_to);
      message(file,"\t\tfunction\tderivation\tmeasured\tder-m\t\t(der-m)'/der'\n",where_to);

      message(file,"norm:\t\t",where_to);
      for (int i=0;i<5;i++) {
        sprintf(buffer,"%e\t",vysl[i*3+0]);
        splachni(file,buffer,where_to);
      }
      message(file,"\n",where_to);
      message(file,"max(Re):\t",where_to);
      for (int i=0;i<5;i++) {
        sprintf(buffer,"%e\t",vysl[i*3+1]);
        splachni(file,buffer,where_to);
      }
      message(file,"\n",where_to);
      message(file,"max(Im):\t",where_to);
      for (int i=0;i<5;i++) {
        sprintf(buffer,"%e\t",vysl[i*3+2]);
        splachni(file,buffer,where_to);
      }
      message(file,"\n",where_to);
    }

}

//#include "region2d_template.cpp"

//template void reg1::test(operators op, content func(prec,prec), content der(prec, prec), prec *output);
//void reg1::test(operators op, content func(prec,prec), content der(prec, prec), prec *output);

void test1(int dimx, int dimy, FILE* file, int where_to)
{
  prec vysl[2][reg1::all_op][reg1::all_prec];
  for (int i=0;i<reg1::all_op;i++)
    for (int j=0;j<2;j++)
      test2(dimx,dimy,(reg1::operators) i,(bool) j,vysl[j][i],logg,2);

  message(file,"operators control - - - results:\n",where_to);

  sprintf(buffer,"at dimensions x:%i\t y:%i\n",dimx, dimy);
  splachni(file, buffer, where_to);

  message(file,"precision level\t",where_to);
  for (int i=0;i<reg1::all_op;i++) {
    sprintf(buffer,"op # %i\t\t",i);
    splachni(file, buffer, where_to);
  }
  message(file, "\n", where_to);

  for (int j=0;j<reg1::all_prec;j++) {
    for (int k=0;k<2;k++) {
      sprintf(buffer,"prec %i sym %i\t",j,k);
      splachni(file, buffer,where_to);
        for (int i=0;i<reg1::all_op;i++) {
          sprintf(buffer,"%.2e\t",vysl[k][i][j]);
          splachni(file, buffer, where_to);
        }
      message(file, "\n", where_to);
    }
  }
}

content fnc1(prec x, prec y)
{
  return(content(3*sin(x),1/(1+y*y)));
}

content fnc2(prec x, prec y)
{
  return(content(4-sin(1/(1+y*y)),11.2*(x-y)));
}

content fnc3(prec x, prec y)
{
  return(content(5*(x+y),2*sin(x+y)));
}

content fnc4(prec x, prec y)
{
  return(content(y-sin(x),x*y*x));
}


void check_for_hermitivity(reg1& r1, FILE* file, int where_to)
{
  randomize();
  content *v, *w;
  for (int test=0;test<10;test++) {
    /*int dim_x=(int) (generuj()*100)+10;
    dim_x=50;
    int dim_y=(int) (generuj()*100)+10;
    dim_y=50;
    prec length_x=parameters.values[parameters_class::length_x];
    prec length_y=parameters.values[parameters_class::length_y];
 
    length_x*=(1+generuj()/10);
    length_y*=(1+generuj()/10);
    length_x=1;length_y=1;
    reg1 r1(dim_x,-length_x,2*length_x,false,dim_y,-length_y,2*length_y,false,inside);*/
    int sets=(int) parameters.read("sets","-");
    v=new content[2*sets*r1.Nw];
    w=new content[2*sets*r1.Nw];

    for (int i=0;i<sets*r1.Nw;i++) {
      v[i]=content(generuj(),generuj());
      w[i]=content(generuj(),generuj());
    }

    r1.operate_tot(v,v+sets*r1.Nw);
    r1.operate_tot(w,w+sets*r1.Nw);

    content val1=0,val2=0;
    for (int i=0;i<sets*r1.Nw;i++) {
      val1+=conj(w[i])*v[i+sets*r1.Nw];
      val2+=conj(conj(v[i])*w[i+sets*r1.Nw]);
    }

    //sprintf(buffer, "test # %i with dimension x: dim=%i && l=%e\t y:dim=%i,l=%e", test, dim_x, length_x,dim_y,length_y);
    //splachni(file, buffer, where_to);

    sprintf(buffer, "test #%i for Hamiltonian hermitivity: <v1Hv2>=%e%+e, <v2Hv1>^*=%e%+e, subtracted:%e%+e\n", test, val1.real(),val1.imag(),val2.real(),val2.imag(), val1.real()-val2.real(),val1.imag()-val2.imag());
    splachni(file, buffer, where_to);
    delete v, w;
  }
}



void ONO_table(ret_from_arpack_zn* pvysl,int size,int eig, FILE* file, int where_to)
{
  for (int i=0;i<eig;i++) {
    for (int j=0;j<eig;j++) {
      content* s1=pvysl->eigenvecs+2*i*size;
      content* s2=pvysl->eigenvecs+2*j*size;
      content val=0;
      for (int spin=0;spin<2;spin++) {
        s1+=spin*size;
	s2+=spin*size;
	for (int k=0;k<size;k++) val+=conj(*(s1+k))*(*(s2+k));
      }
      sprintf(buffer,"%e%+e\t",val.real(),val.imag());
      splachni(file,buffer,where_to);
    }
    message(file,"\n",where_to);
  }
}

void check_diag(FILE* file)
{
//!lapack vs arpack dense matrix diagonalization
  for (int dim=1;dim<100;dim++) {
    complex<double>* A=new complex<double>[dim*dim];
    complex<double>* eigvecs=new complex<double> [dim*dim];
    complex<double>* eigvals=new complex<double> [dim];
    complex<double>* eigvals2=new complex<double> [dim];
    for (int i=0;i<dim*dim;i++) A[i]=complex<double>(generuj(), generuj()); //generate the matrix
    for (int type=0;type<2;type++) {
      char msg[100];
      if (type==0) sprintf(msg,"Hermitian complex"); else sprintf(msg,"general complex  ");
      if (type==0) for (int i=0;i<dim;i++) {//Hermitian matrix
        for (int j=i+1;j<dim;j++) A[i*dim+j]=conj(A[j*dim+i]);
        A[i*dim+i]=A[i*dim+i].real();
      }
      stopky t=stopky();
      int pntr;
      double time[4]={0};
      lapack_zgee_prototype* lapack=0;
      
      //diagonalize with arpack and lapack
      for (int which=0;which<3;which++) {
        t.pusti();
        pntr=0;
        do {
          switch (which) {
            case 0 : {diagonalize_arp(dim,dim,A,eigvecs,eigvals,"SM"); break;}
            case 1 : {
              if (lapack!=0) delete lapack;
              lapack=new lapack_zgee_prototype(dim); 
              break;
            }
            case 2 : {
              lapack->diagonalize(A,eigvals2);
              for (int i=0;i<dim;i++) for (int j=0;j<dim-i-1;j++) {
                if (abs(eigvals2[j])>abs(eigvals2[j+1])) {
                  complex<double> aux=eigvals2[j]; eigvals2[j]=eigvals2[j+1];eigvals2[j+1]=aux;
                }
              } 
              break;
            }
          }
          pntr++;
        } while (t.doba()<10);
        time[which]=double(t.zastav())/pntr;
      }
      delete lapack;
      
      //compare the eigenvalues
      double diff=0;
      for (int i=0;i<dim;i++) diff+=abs(eigvals[i]-eigvals2[i]);
      if (diff>1e-5) {
        fprintf(logg,"at dim=%i eigenvalues differ by too much (%e):\neigval#\t\tarpack\t\tlapack\t\tdifference\n",dim,diff);
        for (int i=0;i<dim;i++) fprintf(logg,"(%.5e%+.5e)\t(%.5e%+.5e)\t(%.5e%+.5e)\n", eigvals[i].real(), eigvals[i].imag(), eigvals2[i].real(), eigvals2[i].imag(), (eigvals[i]-eigvals2[i]).real(), (eigvals[i]-eigvals2[i]).imag());
      }
      
      //put out the results
      fprintf(file, "for %s matrix of dimension %i, times: arpack:%.1e lapack:%.1e (once:%.1e, diag %.1e) [ms] (eigvals differ by %e)\n", msg, dim, time[0], time[1]+time[2], time[1], time[2], diff);
    }
    delete[] A;
    delete[] eigvecs;
    delete[] eigvals;
    delete[] eigvals2;
  }
}

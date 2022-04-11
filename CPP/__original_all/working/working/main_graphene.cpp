#include "main.h"

int main()
{
  parameters.set("system",'g',"-");	//graphene
  
  //!graphene parameters
  parameters.set("vf",1e+6,"m/s");	//fermi velocify
  parameters.set("gx",2,"-");		//g factors
  parameters.set("gy",2,"-");
  parameters.set("gz",2,"-");
  
  parameters.set("Dg",1e-3*0,"meV");	//Delta gap
  parameters.set("Di",12e-3*0,"meV");	//Delta intrinsic
  parameters.set("Dr",48e-3*0,"meV");	//Delta Bychkov-Rashba
  
  parameters.set("B",0.0,"Tesla");	//magnetic field (earth's magn field ~5E-5T)
  parameters.set("theta_B",0.0,"pi");	//0->perpendicular field; 0.5->in-plane field
  parameters.set("phi_B",0.0,"pi");	//in-plane orientation with respect to crystallographic x axis

  parameters.set("gr_U",0,"meV");
  //Martin R like
  parameters.set("d",0,"nm");
  parameters.set("phi_d",0,"pi");
  parameters.set("gr_dd",4,"nm");
  parameters.set("gr_w",-100,"nm");
  parameters.set("gr_r0",30,"nm");
  parameters.set("gr_V0",190,"meV");

  parameters.set("extent",5,"-");
  parameters.set("lv",1,"nm");
  parameters.set("grid_step",1,"-"); 
  
  //Christian E like
  /*parameters.set("gr_r0",70,"nm");
  parameters.set("gr_V0",2000,"meV");
  parameters.set("gr_dd",0.011*parameters.read("gr_r0","nm"),"nm");
  parameters.set("gr_w",parameters.read("gr_r0","nm")+10,"nm");

  parameters.set("extent",5,"-");
  parameters.set("lv",1,"nm");
  parameters.set("grid_step",1,"-"); */
  
  //!computational parameters
  parameters.set("eig",10,"-");		//number of single electron states to obtain in diagonalization
  parameters.set("eigmax",100,"-");	//max. number of states in sorting for the previous - must be large enough to accomodate all crossings.
                          
  //!sorting
  parameters.set("threshold",3E-1,"-");	//tolerated unprecision in sorting accoring to the symmetries.
  parameters.set("sorting_method",6,"-");
  //0->no sorting; 1->I; 2->Ix,Iy; 3->lz; 4->distance; 5->FD-states; 6->g-jz *-1->adds spin.

  parameters.set("sets",dof.length,"-");//number of components of the wavefunction - defined in struct dof in classes.h
  parameters.set("peierls",1,"-");

  const char* name[]={"measuredA","computedA","plottedA"};
  output_class output(3,"data/",name,"log");

  output.clear_files();                     //clear the output files
  output.open_files();                     //clear the output files
  
  fprintf(logg,"dof class was set with Nc=%i, length=%i\n",dof.Nc, dof.length);
  
  //TraceOn('z', 1, 0, 0, 0);
  
  //-------------------------------------------------
  //!the outer loop
  for (prec p1=90;p1<101;p1+=100) {
  //parameters.set("gr_w",p1,"nm");
  //fprintf(output.file(0),"%e ",p1);

  
  //!initialize the sorting class
  sorting_class sorting1e;
  sorting1e.init((int) parameters.read("eigmax","-"), (int) parameters.read("eig","-"));

    //-------------------------------------------------
    //!the inner loop
    //for (prec p2=2;p2>0.99;p2-=0.1) {
      //parameters.set("lv",p2,"nm");
      //parameters.set("grid_step",1,"-"); 


    for (prec p2=0;p2<8;p2+=0.1) {
      parameters.set("B",p2,"Tesla");

      fprintf(output.file(0),"%e ",p2);

      //!quantum dot geometry
      //double dot_radius=parameters.read("gr_r0","nm");  //dot radius in nm, at which the hard wall boundary is imposed 
      //double grid_step=1.8;    //distance between grid points in nm
      //conversion of the dot dimensions into GaAs parametrization
      //parameters.set("lv",grid_step,"nm");
      //parameters.set("extent",dot_radius/grid_step+1,"-");
      //parameters.set("grid_step",1,"-"); 
      
      
      fprintf(logg,"parameters set: p1=%f\n",p1);
      parameters.recompute();
      state_info *states=diagonalize();
      
      //#######################################
      //check for algorithm failure inspecting the imaginary part of the energy
      
      int eig=(int) parameters.read("eig","-");
      int eigmax=(int) parameters.read("eigmax","-");
      
      double imEmax=0, imEtol=1e-5;
      int failed=0;
      for (int i=0;i<eig;i++) {
	double MMmax=0, VVmax=0;
	int jMMmax=i, jVVmax=i;
 	for (int j=0;j<eig;j++) {
	  if (j==i) continue;
          //ONO failure measures:
          //M|v1> with M|v2>
          content MM=0;
	  int length=states->r1->g_length;
	  content *i0=states->vysl->eigenvecs_plaq+i*length, *j0=states->vysl->eigenvecs_plaq+j*length;
	  for (int k=0;k<length;k++) MM+=conj(i0[k])*j0[k];
	  if (abs(MM)>MMmax) {MMmax=abs(MM);jMMmax=j;}

	  //|v1> with |v2>
	  content VV=0;
	  length=states->r1->Nw*states->r1->sets;
	  i0=states->vysl->eigenvecs+i*length;
	  j0=states->vysl->eigenvecs+j*length;
	  for (int k=0;k<length;k++) VV+=conj(i0[k])*j0[k];
	  if (abs(VV)>VVmax) {VVmax=abs(VV);jVVmax=j;}
	}

	//imaginary part of the energy as a failure measure
      	double imE=states->vysl->eigenvals[i].imag()/states->vysl->eigenvals[i].real();
	if (abs(imE)>abs(imEtol)) {//failed
          failed++;
	  printf("imaginary part of state %i too high (rel %e), failure probable\n",i,imE);
	}
        if (abs(imE)>imEmax) imEmax=abs(imE);
	fprintf(logg, "state %2i failure stat: imE=%.3e, MMmax=%.3e (state %2i), VVmax=%.3e (state %2i)\n", i, imE, MMmax, jMMmax, VVmax, jVVmax);
      }
      printf("highest imaginary energy %e, number of probably wrong eigenvectors %i\n",imEmax,failed);
   
      //fprintf(output.file(0),"%e %e %e %e ",auxMM, auxMMlowest, auxvv, auxvvlowest);
      //#######################################
      
      
      //#######################################
      //sort the states
      if (failed==0) {
        message(logg,"sorting...",1);
        sorting1e.inspect_new_1e_states(states);
        sorting1e.sort(p2, states->sorted2un, states->unsorted2s);
        message(logg,"...done\n",1);
      }
      //#######################################
      
      
      //#######################################
      //some info about the single electron states into log and an output file
      states->show_set("clear");
      states->show_set("energy");    
      //states->show_set("energyI");
      states->show_set("g-jz");    
      //states->show_set("sz");    
      states->show_set("sigmaz");    
      //states->show_set("lz");    
      //states->show_set("s");
      //states->show_set("I");
      //states->show_set("Ix");
      //states->show_set("Iy");
      states->show_set("Isigmaz");
      statistics_multi(*states, 0, 1+2+4+8+16, logg,4);
      //#######################################

      
      //#######################################
      //numerical energies into file 1

      for (int i=0;i<eigmax;i++) {
	//fprintf(output.file(0),"%.14e%+ei ",states->read(states->sorted2un[i],"energy"), states->vysl->eigenvals[i].imag()*units.energy/units.meV);
	fprintf(output.file(0),"%.14e ",states->read(states->sorted2un[i],"energy"));
        if (states->sorted2un[i]==-1) continue;
	for (int j=0;j<states->r1->sets;j++) {
	  states->r1->draw(states->vysl->eigenvecs+states->sorted2un[i]*states->r1->Nw*states->r1->sets,j,-1,content(0,0),false,output.file(2),2);
	  fprintf(output.file(2),"\n\n");
	}
      }
      fprintf(output.file(0),"\n");
      //fprintf(output.file(0),"%e %e %e\n", imEsum/eig, sqrt(imEsum2/eig-pow(imEsum/eig,2)), imEmax);
      //#######################################

      /*content psipsi1=0, psipsi2=0, psiIpsi1=0, psiIpsi2=0;
      for (int n=0; n<states->r1->Nw;n++) {
	int ni=states->r1->invertconseq(n, region::inv);
	psipsi1+=conj(states->vysl->eigenvecs[n])*states->vysl->eigenvecs[n];
	psipsi2+=conj(states->vysl->eigenvecs[n+states->r1->Nw])*states->vysl->eigenvecs[n+states->r1->Nw];
	psiIpsi1+=conj(states->vysl->eigenvecs[n])*states->vysl->eigenvecs[ni];
	psiIpsi2+=conj(states->vysl->eigenvecs[n+states->r1->Nw])*states->vysl->eigenvecs[ni+states->r1->Nw];
      }
      printf("lowest eigenstate:\n |psi_1|^2 = %e%+e, |psi_1|^2 = %e%+e, |psi|^2 = %e%+e\n", psipsi1.real(), psipsi1.imag(), psipsi2.real(), psipsi2.imag(), (psipsi1+psipsi2).real(), (psipsi1+psipsi2).imag());
      printf("psi_1^* I psi_1 = %e%+e\n",psiIpsi1.real(), psiIpsi1.imag());
      printf("psi_2^* I psi_2 = %e%+e\n",psiIpsi2.real(), psiIpsi2.imag());
      */
      
      /* //#######################################
      //analytical energies into file 2
      prec Escale=units.hbar*parameters.read("vf","m/s")/(dot_radius*units.nm)/units.meV;
      prec* res=new prec[25];
      prec * ind_n=new prec[25*2];
      prec * ind_l=ind_n+25;
      int * index_rank=new int[25*2];
      for (int l=0;l<5;l++) {
	bessel_solver(l, 5, res+l*5);
	for (int n=0;n<5;n++) {
	  ind_n[l*5+n]=n;
	  ind_l[l*5+n]=l;
	}
      }
      picsrt_index_rank(25, res, index_rank, 'a');
      sort_from_index(25, res, index_rank);
      sort_from_index(25, ind_n, index_rank);
      sort_from_index(25, ind_l, index_rank);
	  
      fprintf(output.file(1),"%.14e ",p1);
      for (int i=0;i<10;i++) fprintf(output.file(1),"%e ",res[i]*Escale);
      fprintf(output.file(1),"\n");
      
      fprintf(logg,"\n\nanalytical energies:\n");
      for (int i=0;i<25;i++) fprintf(logg,"state %i: energy=%e meV, n=%i, l(up)=%i\n",i,res[i]*Escale, (int) ind_n[i], (int) ind_l[i]);
      fprintf(logg,"\n\n");
      delete res;
      delete ind_n;
      delete index_rank;
      //####################################### */
    
      fflush(0);
      clean_up(states);
    }  //!second parameter (p2) loop ends here 
    fprintf(output.file(0),"\n");
  } //!first parameter (p1) loop ends here 
  return(0);
}


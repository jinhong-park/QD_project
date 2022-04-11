#include <iostream>
#include <fstream>

#include "main.h"


int main()
{

  //!GaAs parameters
  parameters.set("eps0",12.9,"-");          //static dielectric constant
  parameters.set("density",5300,"kg/m3");   //density
  parameters.set("cl",5290,"m/s");          //longitudinal velocity
  parameters.set("ct",2480,"m/s");          //transversal velocity
  parameters.set("def",7,"eV");             //deformation potential
  parameters.set("piezo",1.4E9,"eV/m");     //piezoelectric constant
  parameters.set("m_eff",0.067,"m_e");      //effective mass (bulk: 0.067).
    
  parameters.set("gx",-0.39,"-");           //g factors
  parameters.set("gy",-0.39,"-");
  parameters.set("gz",-0.39,"-");
  
  parameters.set("br",1.58/sqrt(8440)*cos(55*M_PI/180),"meVA");            //bychkov-rashba
  parameters.set("dress",1.58/sqrt(8440)*sin(55*M_PI/180),"meVA");         //dresselhaus
  //parameters.set("br",3.3,"meVA");            //bychkov-rashba
  //parameters.set("dress",4.5,"meVA");         //dresselhaus
  parameters.set("dress3",27.5*0,"eVA^3");    //dresselhaus cubic
  parameters.set("T",0,"K");                //temperature
  
/*
  //!Si parameters
  //using (001)-grown ([001] plane confinement, 2DEG symmetry: C2v) SiGe/Si/SiGe (x=0.25) structure
  parameters.set("eps0",11.9,"-");          //static dielectric constant, see Landolt-Boernstein
  parameters.set("density",2330,"kg/m3");   //density  http://jas.eng.buffalo.edu/education/solid/unitCell/home.html
  parameters.set("cl",9150,"m/s");          //longitudinal velocity. SOURCE: Landolt-Boernstein, Grp IV-Elements, Si
  parameters.set("ct",5000,"m/s");          //transversal velocity. SOURCE: Landolt-Boernstein, Grp IV-Elements, Si
  double Xi_d=5.0;
  parameters.set("def",Xi_d,"eV");          //deformation potential Xi_d (5eV and -10eV possible). SOURCE: Cardona, arXiv:0801.4898
  parameters.set("def2",9.0,"eV");          //deformation potential Xi_u. SOURCE: Cardona, Landolt-Boernstein
  parameters.set("piezo",0,"eV/m");         //piezoelectric constant
  parameters.set("m_eff",0.198,"m_e");      //effective mass (bulk: 0.191). Slightly increased for strained Silicon:
                                            //sSi/Si(1-x)Ge(x): 0.198..0.199 for x=0.1..0.4. SOURCE: Phys. Rev. B 48 14276

  parameters.set("gx",2,"-");               //g factors
  parameters.set("gy",2,"-");
  parameters.set("gz",2,"-");

  parameters.set("br",0.05*1,"meVA");       //bychkov-rashba. SOURCE: Nestoklon, PRB 77 155328 (2008)
  parameters.set("dress",0.15*1,"meVA");    //dresselhaus for 5nm wide QW. SOURCE: Nestoklon, PRB 77 155328 (2008)
  parameters.set("dress3",0.,"eVA^3");      //dresselhaus cubic                                          
*/

  //!quantum dot setup
  //parameters.set("width_z",11.0,"nm");
  //parameters.set("lv",34,"nm");             //confinement length l_0
  
  parameters.set("B",3.0,"Tesla");         //magnetic field (earth's magn field ~5E-5T)
  parameters.set("theta_B",0.5*(90.0/90.0),"pi");        //0->perpendicular field; 0.5->in-plane field
  parameters.set("phi_B",0.0,"pi");          //in-plane orientation with respect to crystallographic x axis

  parameters.set("ALHO",1,"-");            //elliptic dot
  prec aux=2.86;
  parameters.set("lv",34*sqrt(aux),"nm");
  parameters.set("d",34*(aux-sqrt(aux)),"nm");
  //parameters.set("ALHO",0,"-");
  //parameters.set("d",0,"nm");  		        //half of the distance between quantum dot minima
  //parameters.set("lv",32,"nm");
  parameters.set("phi_d",0,"pi");         //in-plane orientation with respect to crystallographic x axis
  

  //!computational parameters
  int J = 20;                               //[./.]
  parameters.set("eig",J,"-");              //[J] number of single electron states (including spin) to obtain in diagonalization. Must be large enough because we will kick states when sorting.
  parameters.set("eigmax",2*J,"-");         //[4*eig] max. number of states in sorting for the previous - must be large enough to accomodate all crossings. Must be at least 3*eig.
  parameters.set("selected2e",4,"-");		//number of single electron states (no spin) to built the two electron basis; due to the spin it must be at least eig >= 2*selected, but if the states are sorted, it should be larger to accomodate all crossings
  parameters.set("eig2e",4,"-");    		//goal - number of the two electron states to get

  //!grids
  parameters.set("extent",5,"-");		        //[5] Don't change!
  parameters.set("grid_step",0.2,"-"); 	    //[0.2/0.15] determines grid dimension and single electron functions precision. (For comparison: dim = 2*extent/grid_step (approx) => dim = 50 @ grid_step=0.2).
  parameters.set("CE_exp",20,"-");          //[10/20] Sets precision of the relaxation integrals.
  parameters.set("sets",2,"-");             //=2S+1, where S is the spin of the single electron case (sets=2 => electrons with spin 1/2, which is needed for single-e relax)).
                          
  //!sorting
  parameters.set("threshold",8E-1,"-");		  //[1E-5] tolerated unprecision in sorting accoring to the symmetries.
  parameters.set("sorting_method",-4,"-"); 	//0->no sorting; 1->I; 2->Ix,Iy; 3->lz; 4->distance; 5->FD-states; *-1->adds spin. NOTE: 1e and b states MUST BE fine! f states need only to be fine up to the states needed.
                                            //example: use -1 for double dot (vs magnetic field) spectrum
  
  const char* name[]={"measured-relax.dat","measured-1e.dat","measured-relax-ave.dat"};
  output_class output(2,"data/",name,"log");
  
  
  //!the outer loop
  //int configs=10;
  //prec total_sum[100]={0};
  //for (prec p1=0;p1<configs;p1++)
  for (prec p1=0;p1<1;p1+=1)
  {
    //parameters.set("theta_B",0.5*(90-p1)/90,"pi");
    //parameters.set("phi_d",p1,"pi");
    
    //output.set_appendix("phi_d",(int) round(p1*100));      //appendix to the output files
    //output.set_appendix("iter",app);      //appendix to the output files
    output.clear_files();                          //clear the output files
    output.open_files();
    
    //!impurities
    /*parameters.set("imp_str",1,"meV");
    parameters.set("imp_ext",5,"nm");
    parameters.set("imp_num",10,"-");
    impurities.generate();*/
    
    //!hyperfine spins
    parameters.recompute();
    parameters.set("nuclear_state",0,"-");
    //single impurity volume
    double v0=pow(units.aGaAs/units.nm,3)/8;
    //coupling constant
    double betarho0=2e-3/v0;
    parameters.set("exch",-betarho0,"meV");
    fprintf(logg,"coupling constant: A=%e meV\n", -betarho0);
    //effective volume and number of impurities in a grid cell
    double hx=2*parameters.values[parameters_class::length_x]/parameters.values[parameters_class::dim_x];
    double hy=2*parameters.values[parameters_class::length_y]/parameters.values[parameters_class::dim_y];
    double width_z=parameters.read("width_z","nm");
    double V=(hx*units.length/units.nm)*(hy*units.length/units.nm)*(2*width_z/3); 	//grid cell volume
    double Ni=V/v0;
    fprintf(logg,"impurity volumes: V=%e v0=%e Ni=%e\n",V, v0, Ni);
    //scale of the effective grid cell impurity spin vector
    double Ieff=(2.0/5.0)*sqrt(15.0/4.0/Ni); //random spins
                                             //double Ieff=(2.0/5.0)*(3.0/2.0);// uniform polarization
    parameters.set("xMn",Ieff,"-");
    

    state_info *states=0, *states_previous=0;

    //!initialize the sorting classes for single
    sorting_class sorting1e;
    sorting1e.init((int) parameters.read("eigmax","-"), (int) parameters.read("eig","-"),3);
    
    //!and for double electron states
    sorting_class sorting2e;
    int Ns=(int) parameters.read("selected2e","-");
    int* s2un2e=new int [Ns*Ns];
    int* un2s2e=new int [Ns*Ns];
    for (int a=0;a<Ns*Ns;a++) un2s2e[a]=s2un2e[a]=a;
    sorting2e.init(Ns*Ns,Ns*Ns,1);
    
    //!the inner loop
    prec step_max=0.01, step_min=0.01, step=0.01;
    //class stack minimum(3,3,"min");
    
    for (prec p2=0;p2<1;p2+=step*2) {
      
      //parameters.set("lv",34*sqrt(p2),"nm");
      //parameters.set("d",34*(p2-sqrt(p2)),"nm");
      //parameters.set("B",p2,"Tesla");
      parameters.set("phi_B",p2,"pi");
      //parameters.set("br",p2*cos(39*M_PI/180),"meVA");            //bychkov-rashba
      //parameters.set("dress",p2*sin(39*M_PI/180),"meVA");         //dresselhaus
      //if (p2<134.01 || p2>140) step=2; else step=0.2;
      
      //prec Ebias=p2*1e-3*units.meV/(2*units.e*parameters.read("d","nm")*units.nm);
      //fprintf(logg,"bias energy of %f mueV converted to field %.2e V/m\n",p2,Ebias);
      //parameters.set("Ebias",Ebias,"V/m");

      sprintf(buffer,"parameters set: p1=%f, p2=%f\n",p1,p2);
      splachni(logg,buffer,1+4);
      
      //!the single electron diagonalization
      parameters.recompute();
      states=diagonalize();

      message(logg,"sorting...",1);
      sorting1e.inspect_new_1e_states(states);
      sorting1e.sort(p2, states->sorted2un, states->unsorted2s);
      message(logg,"...done\n",1);

      //some info about the single electron states into log and an output file
      fprintf(output.file(1),"%.14e ",p2);
      states->show_set("clear");
      states->show_set("energy");
      statistics_multi(*states, 0, 2+8, output.file(1),2);
  
      states->show_set("Ix");states->show_set("Iy");
      states->show_set("I");
      states->show_set("s");states->show_set("lz");
      //states->show_set("sx");states->show_set("sy");
      //states->show_set("sz");
      statistics_multi(*states, 0, 1+2+4+8+16, logg,4);
  
      //!relax rates
      //the looping parameter into output file
      sprintf(buffer,"%.14e ",p2);
      splachni(output.file(0),buffer,2);
      
      //for the lowest spin down state, first excited spin up to the ground state - phonon/interaction type
      content* temp=new content[states->r1->Nw*states->r1->sets];
      int from_last=4;
      for (int from=0; from<from_last; from++ ) {//initial state
        if (states->read(from,"s")*states->read(0,"s")>0) continue;
        from_last=from;
        for (int to=0; to<1; to++) {    //final state
          prec rate_sum=0;
          for (int which=0; which<3; which++) {//channel
            //#0: Def-LA (def=Xi_d,def2=Xi_u); #1,#2: Piezo(LA,TA); #3: Def-TA (def2=Xi_u)
            prec rate=0;
            if (from!=to) rate=transition_rate(*states,states->sorted2un[from],states->sorted2un[to],which);
            //fprintf(output.file(0),"%e ",rate);
            rate_sum+=rate;
          }
          fprintf(output.file(0),"%e ",rate_sum);
          //total_sum[(int) (p2*100)]+=rate_sum;
          //minimum.another_in(p2,rate_sum,"min");
          
          //dipole matrix elements
          /*
          prec phi_d=parameters.read("phi_d","pi")*M_PI;
          content* ket=states->vysl->eigenvecs+states->r1->Nw*states->r1->sets*states->sorted2un[to];
          content* bra=states->vysl->eigenvecs+states->r1->Nw*states->r1->sets*states->sorted2un[from];
          states->r1->op_ket(ket, temp, states->r1->r_x, states->r1->set);
          content x=0;
          for (int i=0;i<states->r1->sets;i++) x+=states->r1->braket(bra,i,temp,i); 
          states->r1->op_ket(ket, temp, states->r1->r_y, states->r1->set);
          content y=0;
          for (int i=0;i<states->r1->sets;i++) y+=states->r1->braket(bra,i,temp,i);
          content rd=x*cos(phi_d)+y*sin(phi_d), rcd=x*sin(phi_d)-y*cos(phi_d);
          prec c=units.length/units.nm;
          fprintf(output.file(0),"%e %e %e ",rate_sum, abs(rd)*c, abs(rcd)*c);
           */
        }
      }
      prec lbr=parameters.read("lbr","nm");
      prec ld=parameters.read("ldress","nm");
      prec lso=1/sqrt(1/ld/ld+1/lbr/lbr)/1000;
      fprintf(output.file(0),"%e \n", lso);
      //fprintf(output.file(0),"\n");
      
      /*if (minimum.filled>2) {//estimate where the minimum occurs
        prec estimate=minimum.coordinate(0);
        step=estimate-p2;
        if (step>step_max) step=step_max;
        if (step<0) {
          //fprintf(output.file(2),"%e %e %e\n", p1, estimate, minimum.value(estimate));
          //break;
          step=step_max;
        } else if (step<step_min) step=step_min;
      }*/
      
      
      //!two electron diagonalization
      /*two_electron_state_prototype state2e(states, &sorting2e);
      state2e.diagonalize(p2, s2un2e, un2s2e);
      
      prec exchange=state2e.read('b',2,"energy")-state2e.read('b',0,"energy");
      prec charging=state2e.read('b',0,"energy")-states->read(0,"energy");
      printf("2e ground state energy: %.2f meV, 1e ground state energy: %.2f meV\ncharging energy %.2f meV\n",state2e.read('b',0,"energy"),states->read(0,"energy"),charging);
      return(0);*/

      //step=sorting1e.adjust_step();
      if (states_previous!=0) clean_up(states_previous);
      states_previous=states;
      fflush(0);
    } //!second parameter (p2) loop ends here
    clean_up(states);
    sorting1e.deallocate();
  } //!first parameter (p1) loop ends here
    //for (int i=0;i<100;i++) fprintf(output.file(2),"%e %e\n",i/100.0,total_sum[i]/configs);

  return(0);
}


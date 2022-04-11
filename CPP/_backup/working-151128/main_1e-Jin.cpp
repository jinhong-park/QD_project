#include <iostream> 
#include <fstream>

#include "main.h"

//sets the parameters of the eliptic potential accoridng to the excitation energies
//input the excitation energies along the main axes in meV
void set_ALHO(double hwx, double hwy)
{
    double r=hwy/hwx, gmean=sqrt(hwx*hwy);
    if (logg!=0) {
        if (r<0) {fprintf(logg, "Warning: hwy/hwx negative, swapping sign in set_ALHO\n"); r=-r;}
        if (r<1) {fprintf(logg, "Warning: hwx > hwy, swapping them in set_ALHO\n"); r=1/r;}
    }
    double meff=parameters.read("m_eff","m_e");
    double lv=sqrt(units.hbar*units.hbar/(gmean*units.meV*meff*units.m_e))/units.nm;
    parameters.set("lv",lv,"nm");
    double d=lv*(sqrt(sqrt(r))-1.0);
    parameters.set("d",d,"nm");
    if (logg!=0) fprintf(logg,"from hwx=%.2e, hwy=%.2e, ratio=%.2e, geometric mean=%.2e the ALHO parameters lv=%e, d=%e\n",hwx,hwy,r,gmean,lv,d);
}

void interpolate_hw(double V, double& hwx, double& hwy)
{
    double V_data[]={-200,-150,-100,-50,0,50,100};
    content hwx_data[]={1.57, 1.54, 1.39, 1.40, 1.25, 1.11, 0.89};
    content hwy_data[]={2.65, 2.61, 2.42, 2.33, 2.21, 1.95, 1.68};
    int offset=0;
    while (V_data[offset+1]<V && offset<4) offset++;
    double err;
    hwx=polint(V_data+offset, hwx_data+offset, 3, V, err).real();
    hwy=polint(V_data+offset, hwy_data+offset, 3, V, err).real();
}

//varphi in deg, lso in mu m
void set_so(double varphi, double lso)
{
    varphi*=M_PI/180.0;
    lso*=units.nm*1000;
    double so=units.hbar*units.hbar/(units.m_e*parameters.read("m_eff","m_e")*lso)/(units.meV*units.nm/10);
    parameters.set("br",cos(varphi)*so,"meVA");            //bychkov-rashba
    parameters.set("dress",sin(varphi)*so,"meVA");         //dresselhaus
}

int main()
{
    int n;

    n=2;

    dof.reset_Ndof(n);   // 1
    
    //!GaAs parameters
    parameters.set("eps0",12.9,"-");          //static dielectric constant
    parameters.set("density",5300,"kg/m3");   //density
    parameters.set("cl",5290,"m/s");          //longitudinal velocity
    parameters.set("ct",2480,"m/s");          //transversal velocity
    parameters.set("def",7,"eV");             //deformation potential
    parameters.set("piezo",1.4E9,"eV/m");     //piezoelectric constant
    parameters.set("m_eff",0.067,"m_e");      //effective mass (bulk: 0.067).
    
    parameters.set("gx",-0.38,"-");           //g factors
    parameters.set("gy",-0.38,"-");
    parameters.set("gz",-0.38,"-");
    
    //set_so(60,20);
    parameters.set("br",0.0,"meVA");            //bychkov-rashba   // 3.3
    parameters.set("dress",0.0,"meVA");         //dresselhaus      // 4.5
    parameters.set("dress3",0.0,"eVA^3");      //dresselhaus cubic     // 11
    //parameters.set("T",0,"K");                //temperature
    
    //!quantum dot setup
    parameters.set("width_z",11.0,"nm");       //width of the 2DEG
    //parameters.set("lv",34,"nm");             //confinement length l_0
    
    parameters.set("B",1.0E-5,"Tesla");         //magnetic field (earth's magn field ~5E-5T)
    parameters.set("theta_B",90/180.0,"pi"); //out-of-plane orientation 0:perpendicular; 0.5:in-plane
    parameters.set("phi_B",0,"pi");          //in-plane orientation with respect to crystallographic x axis
    
    parameters.set("ALHO",1,"-");            //elliptic dot
    {
        double hwx,hwy, Vshape=0;               //values for given Vshape
        interpolate_hw(Vshape,hwx,hwy);
        set_ALHO(hwx,hwy);
    }
    //parameters.set("ALHO",0,"-");
    //parameters.set("d",0,"nm");  		        //half of the distance between quantum dot minima
    //parameters.set("lv",32,"nm");
    parameters.set("phi_d",60.0/180.0,"pi");    //dot potential in-plane orientation with respect to crystallographic x axis
    
    
    //!computational parameters
    int J = 20;                               //[./.]
    parameters.set("eig",J,"-");              //[J] number of single electron states (including spin) to obtain in diagonalization. Must be large enough because we will kick states when sorting.
    parameters.set("eigmax",2*J,"-");         //[4*eig] max. number of states in sorting for the previous - must be large enough to accomodate all crossings. Must be at least 3*eig.
    parameters.set("selected2e",10,"-");	  //number of single electron states (no spin) to built the two electron basis; due to the spin it must be at least eig >= 2*selected, but if the states are sorted, it should be larger to accomodate all crossings
    parameters.set("eig2e",4,"-");    		//goal - number of the two electron states to get
    
    //!grids
    parameters.set("extent",5,"-");		        //[5] Don't change!
    parameters.set("grid_step",0.2,"-"); 	    //[0.2/0.15] determines grid dimension and single electron functions precision. (For comparison: dim = 2*extent/grid_step (approx) => dim = 50 @ grid_step=0.2).
    parameters.set("CE_exp",20,"-");          //[10/20] Sets precision of the relaxation integrals.
    parameters.set("sets",pow(2,n),"-");             //=2S+1, where S is the spin of the single electron case (sets=2 => electrons with spin 1/2, which is needed for single-e relax)).
     dof.reset_Ndof(n);
    
    //!sorting
    parameters.set("threshold",8E-1,"-");		  //[1E-5] tolerated unprecision in sorting accoring to the symmetries.
    parameters.set("sorting_method",-4,"-"); 	//0->no sorting; 1->I; 2->Ix,Iy; 3->lz; 4->distance; 5->FD-states; *-1->adds spin. NOTE: 1e and b states MUST BE fine! f states need only to be fine up to the states needed.
    //example: use -1 for double dot (vs magnetic field) spectrum
    
    const char* name[]={"energies","g-factor","parameters"};
    output_class output(3,"data/Jin/2-3D-confinement/",name,"log");
    
    char app[1000];
    sprintf(app, "_phi_d_phi_B_%i",1);
    //output.set_appendix(app);                    //appendix to the output files
    output.clear_files();                          //delete the output files
    output.open_files();

    {
        parameters.recompute();
        double lbr=parameters.read("lbr","nm");
        double ldress=parameters.read("ldress","nm");
        double lz=parameters.read("width_z","nm");
        double m_eff=parameters.read("m_eff","-");
        double d=parameters.read("d","nm");
        double lv=parameters.read("lv","nm");
        double lx=lv*(1+d/lv);
        double ly=lv/(1+d/lv);
        double hwx=units.hbar*units.hbar/(units.m_e*m_eff*lx*lx*units.nm*units.nm)/units.meV;
        double hwy=units.hbar*units.hbar/(units.m_e*m_eff*ly*ly*units.nm*units.nm)/units.meV;

        fprintf(output.file(2),"lbr=%.0f nm, ldress=%.0f nm, lz=%.0f nm, lx=%.1f nm, ly=%.1f nm, hwx=%.1f meV, hwy=%.1f meV\n", lbr, ldress, lz, lx, ly, hwx, hwy);
    }
    
    //!the outer loop
   
    for (int wz=1; wz < 51; wz++){

    for (prec p1=89;p1<180;p1+=200)
    {
        parameters.set("phi_d",p1/180.0,"pi");
        
        state_info *states=0, *states_previous=0;
        
        //!initialize the sorting classes for single
        sorting_class sorting1e;
        sorting1e.init((int) parameters.read("eigmax","-"), (int) parameters.read("eig","-"),3);
        
        //!the inner loop
        for (prec p2=45;p2<180;p2+=200) {
            
            parameters.set("phi_B",p2/180.0,"pi");
            sprintf(buffer,"parameters set: p1=%f, p2=%f\n",p1,p2);
            splachni(logg,buffer,1+4);
            
            

            double dwz;
            dwz = 1.0+2.0*(((double) wz)-1.0);

            parameters.set("width_z",dwz,"nm"); 
                
            //!the single electron diagonalization
            parameters.recompute();
            states=diagonalize();
            
            message(logg,"sorting...",1);
            sorting1e.inspect_new_1e_states(states);
            sorting1e.sort(p2, states->sorted2un, states->unsorted2s);
            message(logg,"...done\n",1);
            
            //some info about the single electron states into log and an output file
            fprintf(output.file(0),"%.14e ",p2);
            states->show_set("clear");
            states->show_set("energy");
            statistics_multi(*states, 0, 2+8, output.file(0),2);
            
            states->show_set("Ix");states->show_set("Iy");
            states->show_set("I");
            states->show_set("s");states->show_set("lz");
            //states->show_set("sx");states->show_set("sy");
            //states->show_set("sz");
            statistics_multi(*states, 0, 1+2+4+8+16, logg,4);

            
            double e0=states->read(states->sorted2un[0],"energy");
            double e1=states->read(states->sorted2un[1],"energy");
            double B=parameters.read("B","Tesla");
            double g=(e1-e0)*units.meV/(B*units.muB);
            fprintf(output.file(1),"%e \t %e \t %e \t %e \t %e\n",dwz,e0,e1,e1-e0,g);
            
            
            if (states_previous!=0) clean_up(states_previous);
            states_previous=states;
            fflush(0);

       

        } //!second parameter (p2) loop ends here
        clean_up(states);
        sorting1e.deallocate();
    } //!first parameter (p1) loop ends here    


         } // width_z loop
    return(0);
}


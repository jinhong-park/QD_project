#include "main.h"

int main()
{

  parameters.set("precision",0,"-");	
  
  parameters.set("gx",-0.44,"-");               //g factors
  parameters.set("gy",-0.44,"-");
  parameters.set("gz",-0.44,"-");
 
  parameters.set("fermi_e",8,"meV");
  
  parameters.set("br",4.5*(1.26*0+0),"meVA");		//bychkov-rashba
  parameters.set("dress",3.3*1.26*0,"meVA");		//dresselhaus
  parameters.set("dress3",27.5*0,"eVA^3");	//dresselhaus cubic 
  
  parameters.set("grid_step",0.2,"-"); 
  parameters.set("lv",30.0,"nm");
  //parameters.set("d",10.0,"nm");		//interdot distance
  parameters.set("phi_d",0.0,"pi");  

  parameters.set("d",4.0,"nm");      //shielding bump range
  parameters.set("d_x",80.0,"nm");   //QPC length
  parameters.set("d_y",1000.0,"nm");   //QPC width

  parameters.set("B",5.25*2,"Tesla");		//magnetic field
  parameters.set("phi_B",0.0,"pi");
  parameters.set("theta_B",0.0,"pi");

  parameters.set("eig",10,"-");		//number of states
  parameters.set("eigmax",100,"-");      //max. number of states in sorting

  parameters.set("geometry",2.0,"-");
  parameters.set("flat_x",100.0,"nm");
  parameters.set("flat_y",100.0,"nm");

  //parameters.set("imp_str",50,"meV");
  //parameters.set("imp_ext",12,"nm");
  //parameters.set("imp_num",10,"nm");
  //impurities.generate();


  transport_class transport;

  const char* name[]={"measured", "measured_sz"};
  output_class output(2,"data/",name,"log");

  //output.set_appendix("varGBxdot100x1000Ef8all",0);  //appendix to the output files
  output.clear_files();                            //clear the output files
  output.open_files();                               //open all output files


  //!QPC probe = charge current derivative and spin current correlations
  /*parameters.set("leads",3.0,"-");
  prec QPC_length=1257+1124*0+2000*0;
  transport.set_lead(1,1,-parameters.read("flat_x","nm"),60,50,-6);
  transport.set_lead(3,2,parameters.read("flat_x","nm")+QPC_length,6,8,-31.6);
  transport.set_lead(2,1,-parameters.read("flat_x","nm"),-60,50,-6);

  extern prec GF_check_Eb;
  for (prec p1=80;p1>79;p1-=1) {
    //if (p1>2) p1=5;
    parameters.set("d_x",p1,"nm");

    output.set_appendix("-covardotNEWalldx",(int) (p1));                   //appendix to the output files
    //output.set_appendix("-try-both-notinQPC-",(int) p1);
    output.clear_files();                            //clear the output files
    output.open_files();                             //open all output files

    parameters.recompute();

    prec dx=parameters.read("d_x","nm");
    prec E0=pow(units.hbar,2)/(2*units.m_e*parameters.read("m_eff","m_e")*pow(dx*units.nm,2))/units.meV;
    prec Emax=E0*pow(QPC_length/2/dx,2);
    prec xmax=dx*sqrt(10/E0);
    fprintf(logg, "QPC length %e [nm], QPC curvature %e [nm]: at QPC end the energy offset reaches %e [meV], drop of 10 meV is reached at x=%e [nm]\n",QPC_length, dx, Emax, xmax);

    prec offset=-pow(parameters.read("dress","meVA"),2)-pow(parameters.read("br","meVA"),2);
    offset*=pow(units.meV*1e-10,2)*units.m_e*parameters.read("m_eff","m_e")/units.hbar/units.hbar/2;
    offset/=units.meV;
    fprintf(logg,"s-o offset is %e meV\n",offset);

    prec nso_phi=M_PI/2;
    if (parameters.read("dress","meVA")!=0) nso_phi=atan(parameters.read("br","meVA")/parameters.read("dress","meVA"));
    fprintf(logg, "spin-orbit field inside the QPC (along x) azimutal vector is %e [pi]\n", nso_phi/M_PI);

    prec ldress=parameters.read("ldress","nm");
    prec lbr=parameters.read("lbr","nm");
    prec leff=0;
    if (ldress!=0) leff=1/pow(ldress,2);
    if (lbr!=0) leff+=1/pow(lbr,2);
    if (leff!=0) leff=2*M_PI/sqrt(leff);
    fprintf(logg, "effective spin-orbit length (corresponding to 2 pi rotation) inside the QPC (along x) is %e [nm], angle of rotation inside QPC is %e [pi]\n", leff, QPC_length/leff*2);

    prec dI[2*2]={0}, dI2[2*2]={0}, Is[2*2]={0}, Is2[2*2]={0}, product[2*2]={0};
    int pntr[2]={0};
    for (prec p2=3;p2<10;p2+=0.01) {
      int which=0;
      if (p2<6.95) which=0; else
        if (p2>7.05) which=1; else continue;
      parameters.set("fermi_e",p2,"meV");
      extern prec GF_check_Eb;

      prec der, der0, t, t0;
  
      //Bstep in Tesla, corresponding to an (Zeeman) energy step of Estep [meV]
      prec Bstep=0.001;
      prec dEdB=parameters.read("gz","-")/4*units.hbar/units.m_e*units.e/units.meV;
      //prec Estep=abs(Bstep*dEdB);
  
      GF_check_Eb=parameters.read("fermi_e","meV")-offset;
      //der0=t0=0;
      //der0=transport.QPC_derivative(Estep,t0);
      der=-units.meV*M_PI/2.0*2.0;
      der/=pow(units.hbar,2)/(units.m_e*parameters.read("m_eff","m_e")*pow(parameters.read("d_x","nm")*units.nm,2));
      //printf("the ratio is %e\n",der0*dEdB);
      //exit(1);

      fprintf(logg,"Bstep for the derivatives is %e [meV], dE/dB is %e [meV/Tesla], and the conversion factor is der*dEdB=%e \n", Bstep, dEdB, der*dEdB);
  
      //der=der0;t=t0;
      //der=transport.set_QPC_sensitive(Estep,Vstep1,threshold,step1,Vstep2, threshold, step2, estimate, t);
      //sprintf(buffer,"fermi=%e V_QPC=%e der=%e at t=%e in %i+%i steps (vs V_QPC_est=%e der0=%e at t0=%e)\n", p1, GF_check_Eb, der, t, step1, step2, estimate, der0, t0);
      //splachni(logg, buffer, 1+4);
      
  
      prec V[3+1];
      V[1]=0.5;
      V[2]=-0.5;
      prec results[8];
  
      transport.derivative(V, Bstep, results);
      sprintf(buffer,"derivative: QPC voltage=%e, dI/dBx=%e (vs theory %e), dI/dBy=%e (vs theory %e), I_spinx=%e, I_spiny=%e, I_spinz=%e\n", V[3], results[0], der*results[2]*dEdB, results[1], der*results[3]*dEdB, results[2], results[3], results[4]);
      splachni(logg, buffer, 1+4);
  
      fprintf(output.file(0),"%e %e %e %e %e %e %e %e %e %e %e %e\n", p2, der*dEdB, results[0], results[2], der*results[2]*dEdB, results[1], results[3], der*results[3]*dEdB, results[4], results[5], results[6], results[7]);
      fflush(0);

      int IorV=0;//!!! 0 for current 6 for voltage
      dI[which]+=results[IorV];
      dI2[which]+=pow(results[IorV],2);
      Is[which]+=results[2];
      Is2[which]+=pow(results[2],2);
      product[which]+=results[IorV]*results[2];

      dI[which+2]+=results[IorV+1];
      dI2[which+2]+=pow(results[IorV+1],2);
      Is[which+2]+=results[3];
      Is2[which+2]+=pow(results[3],2);
      product[which+2]+=results[IorV+1]*results[3];

      pntr[which]++;
    }
    prec vardI[2*2], varIs[2*2], covar[2*2], vardI_tot[2], varIs_tot[2], covar_tot[2];
    int pntr_tot=pntr[0]+pntr[1];
    for (int what=0;what<2;what++) {
      vardI_tot[what]=(dI2[0+what*2]+dI2[1+what*2])/pntr_tot-pow((dI[0+what*2]+dI[1+what*2])/pntr_tot,2);
      varIs_tot[what]=(Is2[0+what*2]+Is2[1+what*2])/pntr_tot-pow((Is[0+what*2]+Is[1+what*2])/pntr_tot,2);
      covar_tot[what]=(product[0+what*2]+product[1+what*2])/pntr_tot-(dI[0+what*2]+dI[1+what*2])/pntr_tot*(Is[0+what*2]+Is[1+what*2])/pntr_tot;
      for (int which=0;which<2;which++) {
        dI[which+what*2]/=pntr[which];
        dI2[which+what*2]/=pntr[which];
        Is[which+what*2]/=pntr[which];
        Is2[which+what*2]/=pntr[which];
        product[which+what*2]/=pntr[which];
      }
    }
    for (int i=0;i<4;i++) {
      vardI[i]=dI2[i]-dI[i]*dI[i];
      varIs[i]=Is2[i]-Is[i]*Is[i];
      covar[i]=product[i]-dI[i]*Is[i];
    }

    fprintf(output.file(0),"XXX %e %i %i ",p1, pntr[0], pntr[1]);
    for (int i=0;i<4;i++) {
      fprintf(output.file(0),"%e %e %e %e %e ", dI[i], Is[i], vardI[i], varIs[i], covar[i]/sqrt(vardI[i]*varIs[i]));
    }
    fprintf(output.file(0),"%e %e %e %e\n", sqrt(vardI_tot[0]/varIs_tot[0]), sqrt(vardI_tot[1]/varIs_tot[1]), covar_tot[0]/sqrt(vardI_tot[0]*varIs_tot[0]), covar_tot[1]/sqrt(vardI_tot[1]*varIs_tot[1]));
    fflush(0);
  }
  return(0);
}*/



  //!QPC probe - as a function of the magnetic field
  /*parameters.set("leads",3.0,"-");
  transport.set_lead(1,1,-parameters.read("flat_x","nm"),60,50,-6);
  transport.set_lead(3,2,parameters.read("flat_x","nm")+1257,6,8,-31.6);
  transport.set_lead(2,1,-parameters.read("flat_x","nm"),-60,50,-6);

  extern prec GF_check_Eb;
  for (prec p1=8.00;p1<10;p1+=10.01) {
    parameters.set("fermi_e",p1,"meV");

    //bias the QPC such that the current derivative is maximal
    prec offset=-pow(parameters.read("dress","meVA"),2)-pow(parameters.read("br","meVA"),2);
    offset*=pow(units.meV*1e-10,2)*units.m_e*parameters.read("m_eff","m_e")/units.hbar/units.hbar/2;
    offset/=units.meV;
    fprintf(logg,"s-o offset is %e meV\n",offset);

    GF_check_Eb=parameters.read("fermi_e","meV")-offset;

    //prec der, der0, t, t0;
    //int step1=10, step2=20;
    //prec threshold=1e-2;
    //prec Vstep1=0.03, Vstep2=0.01, Estep=0.001;
    //der0=transport.QPC_derivative(Estep,t0);
    //der=transport.set_QPC_sensitive(Estep,Vstep1,threshold,step1,Vstep2, threshold, step2, estimate, t);
    //sprintf(buffer,"fermi=%e V_QPC=%e der=%e at t=%e in %i+%i steps (vs V_QPC_est=%e der0=%e at t0=%e)\n", p1, GF_check_Eb, der, t, step1, step2, estimate, der0, t0);
    //splachni(logg, buffer, 1+4);

    prec V[3+1];
    //V[2]=-V[1]=0.5;
    V[2]=V[1]=0;V[3]=1;
    prec T[4*4*16];
    prec I[4];
    int Nl=4;
    int M[4];

    int QPC=3;

    int pntr=0;
    prec step;
    for (prec p2=0;p2<10;p2+=step) {
      parameters.set("B",p2,"Tesla");
      transport.prepare_geometry(false);
      transport.prepare_with_E();
      transport.recursive_green(false);
      transport.scattering_matrix(M);
      for (int i=1;i<Nl;i++) for (int j=1;j<Nl;j++) for (int a=0;a<4;a++) for (int b=0;b<4;b++) {
        T[(i*Nl+j)*16+a*4+b]=transport.lead2lead_spinresolved(i,a,j,b);
        //fprintf(logg,"T_{%i%i}^{%i%i}=%e\n",i,j,a,b,T[(i*Nl+j)*16+a*4+b]);
      }

      //set QPC current to zero at the beggining
      if (pntr==0 && false) {
        V[QPC]=(T[(QPC*Nl+1)*16+0*4+0]*V[1]+T[(QPC*Nl+2)*16+0*4+0]*V[2])/(M[QPC]-T[(QPC*Nl+QPC)*16+0*4+0]);
        sprintf(buffer, "QPC potential is %e [mV]\n\tV[1]=%e, V[2]=%e, T_31^00=%e, T_32^00=%e (prop. channels:%i)\n",V[QPC], V[1], V[2], T[(QPC*Nl+1)*16+0*4+0], T[(QPC*Nl+2)*16+0*4+0], M[QPC]);
        splachni(logg, buffer, 1+4);
        //transport.check_unitarity();
        //transport.check_trs();
      }

      //compute the spin and charge current
      for (int a=0;a<4;a++) {
        I[a]=0;
        if (a==0) I[a]=M[QPC]*V[QPC];
        if (a==2 && M[QPC]==1) I[a]=M[QPC]*V[QPC]*((p2>0)*2-1);
        for (int i=1;i<=QPC;i++) I[a]-=T[(QPC*Nl+i)*16+a*4+0]*V[i];
      }
      prec current_unit=units.e*units.e/(2*M_PI*units.hbar)*1e-3*1e+9;
      for (int a=0;a<4;a++) I[a]*=current_unit;

      fprintf(logg, "currents through QPC: I[0]=%e, I[1]=%e, I[2]=%e, I[3]=%e\n", I[0], I[1], I[2], I[3]);

      fprintf(output.file(0),"%e ", p2);
      for (int a=0;a<4;a++) fprintf(output.file(0),"%e ", I[a]);
      fprintf(output.file(0),"\n");

      transport.clean_with_E();
      transport.clean_geometry();

      //prec results[6];
  
      //transport.derivative(V, 0.001, 0.001, results);
      //fprintf(logg, "results from derivative: %e %e %e %e %e %e\n", results[0], results[1], results[2], results[3], results[4], results[5]);

      
      if (pntr==0) {step=1e-3;} else {
        if (abs(p2)<0.1) step=0.01; else
          if (abs(p2)<1) step=0.1; else
            if (abs(p2)<10) step=0.1; else
             if (abs(p2)<100) step=0.1; else
               step=10;
      }

      if (p2>0) p2+=step;
      p2*=-1;
      pntr++;
      fflush(0);
    }
  }
return(0);
}*/

  /*
  //!QPC - dot with three leads, no QPC stick
  parameters.set("leads",3.0,"-");
  transport.set_lead(1,1,-parameters.read("flat_x","nm"),50,90,0);
  transport.set_lead(3,2,parameters.read("flat_x","nm")+0,9,20,0);
  transport.set_lead(2,1,-parameters.read("flat_x","nm"),-50,90,0);
  int M[4];

  extern prec GF_check_Eb;

  //output.set_appendix("Ef",(int) p1);            //appendix to the output files
  output.clear_files();                            //clear the output files
  output.open_files();                             //open all output files

  transport.prepare_geometry(false);

  for (prec p1=7;p1<15;p1+=0.1) {
    parameters.set("fermi_e",p1,"meV");

    int M[4];
    transport.prepare_with_E();
    transport.recursive_green(false);
    transport.scattering_matrix(M);
    prec R=transport.lead2lead_spinresolved(3,0,3,0);
    transport.clean_with_E();

    fprintf(output.file(0),"Ef=%e [meV] R=%e\n", p1, R);
    fflush(0);
  }
  transport.clean_geometry();

  return(0);
}

  //!var G check - quenching with the magnetic field
  /*parameters.set("leads",2.0,"-");

  prec step=1;
  for (prec p1=0;p1<50;p1+=step) {
    parameters.set("B",p1,"Tesla");

    prec G[4*2],G2[4*2];
    for (int i=0;i<4*2;i++) G[i]=G2[i]=0;
    int pntr[2]={0,0};

    int grid=(int) round(parameters.read("grid_step","-")*parameters.read("lv","nm"));
    int width=50/2;
    int ymax=((int) parameters.read("flat_y","nm")/grid) - width/grid-1;
    int lead1y=ymax;
    int lead2y=lead1y-2*width/grid-2;
    printf("grid=%i, ymax=%i ,lead1y=%i, lead2y=%i\n",grid, ymax, lead1y, lead2y);
    lead2y+=2;
    do {
      printf("y positions: lead1=%i, lead2=%i\n",lead1y, lead2y);
      lead2y-=2;
      if (lead2y<-ymax) {
        lead1y-=2;
        lead2y=lead1y-2*width/grid-2;
        if (lead1y<-ymax || lead2y<-ymax) break;
      }
      printf("reset to: lead1=%i, lead2=%i\n",lead1y, lead2y);
      transport.set_lead(1,1,-parameters.read("flat_x","nm"),lead1y*grid,50,-6);
      transport.set_lead(2,2,parameters.read("flat_x","nm"),lead2y*grid,50,-6);

      transport.prepare_geometry(false);
      for (prec p2=3;p2<11;p2+=2) {
        int which;
        parameters.set("fermi_e",p2,"meV");
        transport.prepare_with_E();
        if (lead[1].propagating_modes==4) which=0; else
          if (lead[1].propagating_modes==6) which=1; else
            continue;
  
        transport.recursive_green(false);
        int M[4];
        transport.scattering_matrix(M);
        prec G_act[4];
        for (int i=0;i<4;i++) G_act[i]=transport.lead2lead_spinresolved(2,i,1,0);
        transport.clean_with_E();
        fprintf(logg,"Ef=%e which=%i I0=%e Ix=%e Iy=%e Iz=%e\n", p2, which, G_act[0],G_act[1],G_act[2],G_act[3]);
        //fflush(0);

        for (int i=0;i<4;i++) {
          G[which*4+i]+=G_act[i];
          G2[which*4+i]+=G_act[i]*G_act[i];
        }
        pntr[which]++;
      }
      transport.clean_geometry();
    } while (true);

    fprintf(output.file(0),"%e %i %i ", p1, pntr[0], pntr[1]);
    for (int which=0; which<2; which++) {
      for (int i=0;i<4;i++) {
        G[which*4+i]/=pntr[which];
        G2[which*4+i]/=pntr[which];
        fprintf(output.file(0),"%e %e ",G[which*4+i], sqrt(G2[which*4+i]-G[which*4+i]*G[which*4+i]));
      }
    }
    fprintf(output.file(0),"\n");
    fflush(0);

    /*if (abs(p1)<0.1) step=0.01; else
      if (abs(p1)<1) step=0.1; else
        if (abs(p1)<10) step=1; else
         if (abs(p1)<100) step=10; else
           step=10;*/
  /*}
  return(0);
}*/

  //!QPC - two leads + channel of the same width
  parameters.set("leads",2.0,"-");
  parameters.set("flat_x",1000,"nm");
  parameters.set("flat_y",4.0,"nm");
  parameters.set("grid_step",0.2,"-");
  //parameters.set("flat_y",100.0,"nm");
  parameters.set("d",0.0,"nm");      //shielding bump range
  transport.set_lead(1,1,-parameters.read("flat_x","nm"),0,8,-31.6);
  transport.set_lead(2,2,parameters.read("flat_x","nm"),0,8,-31.6);
  //!
  //parameters.set("leads",3.0,"-");
  //transport.set_lead(1,1,-parameters.read("flat_x","nm"),60,50,-6);
  //transport.set_lead(2,1,-parameters.read("flat_x","nm"),-60,50,-6);
  //parameters.set("d_x",10000.0,"nm");
  //!

  extern prec GF_check_Eb;
  GF_check_Eb=9;

  for (prec p1=80;p1<81;p1+=1) {
    parameters.set("d_x",p1*1,"nm");
    prec QPC_length=0;

    prec dx=parameters.read("d_x","nm");
    prec E0=pow(units.hbar,2)/(2*units.m_e*parameters.read("m_eff","m_e")*pow(dx*units.nm,2))/units.meV;
    prec Emax=E0*pow(QPC_length/2/dx,2);
    prec xmax=dx*sqrt(10/E0);
    fprintf(logg,"QPC length %e [nm], QPC curvature %e [nm]: at QPC end the energy offset reaches %e [meV], drop of 10 meV is reached at x=%e [nm]\n",QPC_length, dx, Emax, xmax);
    //parameters.set("flat_x",xmax,"nm");
    //transport.set_lead(1,1,-parameters.read("flat_x","nm"),0,8,-31.6);
    //transport.set_lead(2,2,parameters.read("flat_x","nm"),0,8,-31.6);

    //!
    //transport.set_lead(3,2,parameters.read("flat_x","nm")+p1,6,8,-31.6);

    //fprintf(output.file(0),"\n\n");
    //int lead1width=70;
    //prec lead1offset=-9;
    //int lead1posmax=66;
    //int lead1pos=0;
    prec G[4],G2[4];
    for (int i=0;i<4;i++) G[i]=G2[i]=0;
    //int pntr=0;
    //do {
      //lead1pos-=(int) (parameters.read("grid_step","-")*parameters.read("lv","-")*2);
      //transport.set_lead(1,1,-parameters.read("flat_x","nm"),lead1pos,lead1width,lead1offset);
      //lead1pos*=-1;
      //if (lead1pos>=0) lead1pos+=(int) round(parameters.read("grid_step","-")*parameters.read("lv","nm"));
      //if (lead1pos>lead1posmax) break;
      //printf("injecting lead position updated to %i nm\n",lead1pos);
    
      transport.prepare_regions();
  
                //!hyperfine spins
    parameters.set("T", 0.0, "K");         //temperature
    parameters.set("width_z",10,"nm");     //width of the 2DEG	- used to be 3 nm

      //single impurity volume
      double v0=pow(units.aGaAs/units.nm,3)/8;
      //coupling constant
      double betarho0=2e-3/v0;
      parameters.set("exch",-betarho0,"meV");
      fprintf(logg,"hyperfine coupling constant: A=%e meV\n", -betarho0);
      //effective volume and number of impurities in a grid cell
      double hx,hy,width_z=parameters.read("width_z","nm");
      transport.region->give_par(hx,hy);
      double V=(hx*units.length/units.nm)*(hy*units.length/units.nm)*(2*width_z/3); 	//grid cell volume
      double Ni=V/v0;
      fprintf(logg,"impurity volumes: V=%e v0=%e Ni=%e\n",V, v0, Ni);
      //scale of the effective grid cell impurity spin vector
      //double Ieff=(2.0/5.0)*sqrt(15.0/4.0/Ni); //random spins
      double Ieff=(2.0/5.0)*(3.0/2.0);// uniform polarization
      parameters.set("xMn",Ieff,"-");
      spin_impurities.reinitialize(transport.region,0,0);
  
    bool full_GF=true, print_density=true;
    
      transport.prepare_geometry(full_GF);
    
      prec step=0.001, step_min=0.001, step_max=0.001, goal=0.2;
      prec G_old[4];
      for (int i=0;i<4;i++) G_old[i]=0;
      int refused=0;
      prec der=0, dermax=0, p2max=0, dermid=0;
      int pntr=0;
      for (prec p2=GF_check_Eb-0.5*0.2;p2<GF_check_Eb+0.5*1;p2+=step*1) {
        parameters.set("fermi_e",p2,"meV");

        //sprintf(buffer,"parameters set: p1=%.9e, p2=%.9e\n",p1,p2);
        //splachni(logg,buffer,1);
    
        transport.prepare_with_E();
        transport.recursive_green(full_GF);
        int M[4];
        transport.scattering_matrix(M);
        //transport.check_unitarity();
        //transport.check_trs();
        prec G_act[4]={0},G_act3[4]={0};
        for (int i=0;i<4;i++) {
          //G_act[i]=transport.lead2lead_spinresolved(2,0,i,0);
          G_act[i]=transport.lead2lead_spinresolved(2,i,1,0);
          //fprintf(logg,"lead2lead(2,%i,1,0)=%e\n",i,G_act[i]);
        }
        if (print_density) {
          print_density=false;
          transport.oneDspindensity(output.file(1));
          break;
        }
        G_act[2]=(G_act[0]+G_act[3])/2;
        G_act[3]=(G_act[0]-G_act[3])/2;
        der=(G_act[1]-G_old[1])/step;
        if (pntr>1) {if (der>dermax) {dermax=der;p2max=p2;}};
        pntr++;
        if (abs(p2-GF_check_Eb)<1e-6) {dermid=der;}
        //for (int i=0;i<4;i++) G_act3[i]=transport.lead2lead_spinresolved(2,i,3,0);
        transport.clean_with_E();
        prec dif=abs(G_act[0]-G_old[0]);
        if (dif!=0) dif/=max(abs(G_act[0]),abs(G_old[0]));
        fprintf(logg, "actual step: %e, transmissions: new=%e, old=%e, rel. difference: %e\n",step, G_act[0],G_old[0],dif);
        if (dif>goal && step>step_min && refused<0) {p2-=step;refused++;} //refuse
        else {//accept
          refused=0;
          fprintf(output.file(0),"%e %e %e %e %e\n", p2, G_act[0], G_act[1], G_act[2], G_act[3]);
          //fprintf(output.file(0),"%e %e %e %e %e\n", p2, G_act3[0], G_act3[1], G_act3[2], G_act3[3]);
          for (int i=0;i<4;i++) G_old[i]=G_act[i];
        }
        if (dif!=0) step/=dif/goal;
        if (step<step_min) step=step_min;
        if (step>step_max) step=step_max;
        //if (G21<1e-8 || G21>2-1e-8) step=step_max*10;
        //fprintf(logg,"new step set at %e\n",step);
        //fflush(0);
        //for (int i=0;i<4;i++) {
          //G[i]+=G_act[i];
          //G2[i]+=G_act[i]*G_act[i];
        //}
        //pntr++;
      }
      transport.clean_geometry();
    prec theory=M_PI/2.0*units.meV;
    theory/=pow(units.hbar,2)/(units.m_e*parameters.read("m_eff","m_e")*pow(parameters.read("d_x","nm")*units.nm,2));
    //fprintf(output.file(0),"%e %e %e %e %e\n", p1, theory, dermax, p2max, dermid);

    //}
    //while (true);
    //fprintf(output.file(0),"%e %e %i ", p1, parameters.read("ldress","nm"), pntr);
    //for (int i=0;i<4;i++) {
      //G[i]/=pntr;
      //G2[i]/=pntr;
      //fprintf(output.file(0),"%e %e ", G[i], sqrt(G2[i]-G[i]*G[i]));
    //}
    //fprintf(output.file(0),"\n");
  }
  return(0);
}

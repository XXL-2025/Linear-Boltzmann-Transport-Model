///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the key function of LBT
// Update elastic and inelastic evolution of the particle list for one time step
//
///////////////////////////////////////////////////////////////////////////////////////////////////
double Vx_cal(double milvx, double etas, double vzeta)
            {
                double cx;

                cx = milvx/(cosh(etas)+vzeta*sinh(etas));

                return cx;

            }

            double Vy_cal(double milvy, double etas, double vzeta)
            {
                double cy;

                cy = milvy/(cosh(etas)+vzeta*sinh(etas));

                return cy;

            }

            double Vz_cal(double etas, double vzeta) 
            {
            double cz;
            double uz,ut;

            uz = vzeta*cosh(etas)+sinh(etas);
            ut = cosh(etas)+vzeta*sinh(etas);
            cz = uz/ut;

            return cz;
            }

void bulklinear(double tau, double x,double y,double eta, double &temp,double &VX,double &VY,double &VZ,int &fraction)
{


	  
  double dxh=0.3;
  double dyh=0.3;
  double detah=0.3;  
  double dtauh=0.3;

  double x0h=-10.8;
  double y0h=-10.8;
  double eta0h=-4.8;  
  double tau0h=0.6;

/*
  double dxh=dxh_r;
  double dyh=dyh_r;
  double detah=detah_r;  
  double dtauh=dtauh_r;

  double x0h=x0h_r;
  double y0h=y0h_r;
  double eta0h=eta0h_r;  
  double tau0h=tau0h_r;
*/

  
 ///////////////////////////////////////////////////////////////////////////// 
 //cout<< tau <<" "<< x <<" "<< y <<" "<< eta <<" "<< ed <<" "<< temp <<" "<< fraction <<" "<< VX <<" "<< VY <<" "<< VZ << endl;
 ///////////////////////////////////////////////////////////////////////////// 


     //cout << tauhydro0 << "            " << nxhydro << "           " << nyhydro << "           " << netahydro << "             " << dtauhydro <<endl;

  
  
  int ix1=floor((x-x0h)/dxh)+1;
  int ix2=ix1+1;
  int iy1=floor((y-y0h)/dyh)+1;
  int iy2=iy1+1;
  int ieta1=floor((eta-eta0h)/detah)+1;
  int ieta2=ieta1+1;  
  int itau1=floor((tau-tau0h)/dtauh)+1;
  int itau2=itau1+1;  
  
  double x1=x0h+(ix1-1)*dxh;
  double x2=x1+dxh;
  double y1=y0h+(iy1-1)*dyh;
  double y2=y1+dyh;  
  double eta1=eta0h+(ieta1-1)*detah;
  double eta2=eta1+detah;  
  double tau1=tau0h+(itau1-1)*dtauh;
  double tau2=tau1+dtauh;  

  double xs1=(x2-x)/(x2-x1);
  double ys1=(y2-y)/(y2-y1);
  double etas1=(eta2-eta)/(eta2-eta1);
  double taus1=(tau2-tau)/(tau2-tau1);

  double xs2=(x-x1)/(x2-x1);
  double ys2=(y-y1)/(y2-y1);
  double etas2=(eta-eta1)/(eta2-eta1);
  double taus2=(tau-tau1)/(tau2-tau1);

 ///////////////////////////////////////////////////////////////////////////// 
 //cout<<"hahahahahahahaha "<<endl;
 /////////////////////////////////////////////////////////////////////////////  
  
  
  
  temp = xs1*ys1*etas1*taus1*temphydro[itau1][ix1][iy1][ieta1];
  temp+= xs1*ys1*etas1*taus2*temphydro[itau2][ix1][iy1][ieta1];
  temp+= xs1*ys1*etas2*taus1*temphydro[itau1][ix1][iy1][ieta2];
  temp+= xs1*ys1*etas2*taus2*temphydro[itau2][ix1][iy1][ieta2];
  temp+= xs1*ys2*etas1*taus1*temphydro[itau1][ix1][iy2][ieta1];
  temp+= xs1*ys2*etas1*taus2*temphydro[itau2][ix1][iy2][ieta1];
  temp+= xs1*ys2*etas2*taus1*temphydro[itau1][ix1][iy2][ieta2];
  temp+= xs1*ys2*etas2*taus2*temphydro[itau2][ix1][iy2][ieta2];
  temp+= xs2*ys1*etas1*taus1*temphydro[itau1][ix2][iy1][ieta1];
  temp+= xs2*ys1*etas1*taus2*temphydro[itau2][ix2][iy1][ieta1];
  temp+= xs2*ys1*etas2*taus1*temphydro[itau1][ix2][iy1][ieta2];
  temp+= xs2*ys1*etas2*taus2*temphydro[itau2][ix2][iy1][ieta2];
  temp+= xs2*ys2*etas1*taus1*temphydro[itau1][ix2][iy2][ieta1];
  temp+= xs2*ys2*etas1*taus2*temphydro[itau2][ix2][iy2][ieta1];
  temp+= xs2*ys2*etas2*taus1*temphydro[itau1][ix2][iy2][ieta2];
  temp+= xs2*ys2*etas2*taus2*temphydro[itau2][ix2][iy2][ieta2];

/*
cout<<"fraction--------------------------------------------------"<<endl;

  cout<<fractionhydro[itau1][ix1][iy1][ieta1]<<endl;
  cout<<fractionhydro[itau2][ix1][iy1][ieta1]<<endl;
  cout<<fractionhydro[itau1][ix1][iy1][ieta2]<<endl;
  cout<<fractionhydro[itau2][ix1][iy1][ieta2]<<endl;
  cout<<fractionhydro[itau1][ix1][iy2][ieta1]<<endl;
  cout<<fractionhydro[itau2][ix1][iy2][ieta1]<<endl;
  cout<<fractionhydro[itau1][ix1][iy2][ieta2]<<endl;
  cout<<fractionhydro[itau2][ix1][iy2][ieta2]<<endl;
  cout<<fractionhydro[itau1][ix2][iy1][ieta1]<<endl;
  cout<<fractionhydro[itau2][ix2][iy1][ieta1]<<endl;
  cout<<fractionhydro[itau1][ix2][iy1][ieta2]<<endl;
  cout<<fractionhydro[itau2][ix2][iy1][ieta2]<<endl;
  cout<<fractionhydro[itau1][ix2][iy2][ieta1]<<endl;
  cout<<fractionhydro[itau2][ix2][iy2][ieta1]<<endl;
  cout<<fractionhydro[itau1][ix2][iy2][ieta2]<<endl;
  cout<<fractionhydro[itau2][ix2][iy2][ieta2]<<endl;
*/


  fraction = xs1*ys1*etas1*taus1*fractionhydro[itau1][ix1][iy1][ieta1];
  fraction+= xs1*ys1*etas1*taus2*fractionhydro[itau2][ix1][iy1][ieta1];
  fraction+= xs1*ys1*etas2*taus1*fractionhydro[itau1][ix1][iy1][ieta2];
  fraction+= xs1*ys1*etas2*taus2*fractionhydro[itau2][ix1][iy1][ieta2];
  fraction+= xs1*ys2*etas1*taus1*fractionhydro[itau1][ix1][iy2][ieta1];
  fraction+= xs1*ys2*etas1*taus2*fractionhydro[itau2][ix1][iy2][ieta1];
  fraction+= xs1*ys2*etas2*taus1*fractionhydro[itau1][ix1][iy2][ieta2];
  fraction+= xs1*ys2*etas2*taus2*fractionhydro[itau2][ix1][iy2][ieta2];
  fraction+= xs2*ys1*etas1*taus1*fractionhydro[itau1][ix2][iy1][ieta1];
  fraction+= xs2*ys1*etas1*taus2*fractionhydro[itau2][ix2][iy1][ieta1];
  fraction+= xs2*ys1*etas2*taus1*fractionhydro[itau1][ix2][iy1][ieta2];
  fraction+= xs2*ys1*etas2*taus2*fractionhydro[itau2][ix2][iy1][ieta2];
  fraction+= xs2*ys2*etas1*taus1*fractionhydro[itau1][ix2][iy2][ieta1];
  fraction+= xs2*ys2*etas1*taus2*fractionhydro[itau2][ix2][iy2][ieta1];
  fraction+= xs2*ys2*etas2*taus1*fractionhydro[itau1][ix2][iy2][ieta2];
  fraction+= xs2*ys2*etas2*taus2*fractionhydro[itau2][ix2][iy2][ieta2];
  

/*
  
cout<<"                         "<<endl;
cout<<"VX--------------------------------------------------"<<endl;
cout<<"                         "<<endl;

  cout<<VXhydro[itau1][ix1][iy1][ieta1]<<endl;
  cout<<VXhydro[itau2][ix1][iy1][ieta1]<<endl;
  cout<<VXhydro[itau1][ix1][iy1][ieta2]<<endl;
  cout<<VXhydro[itau2][ix1][iy1][ieta2]<<endl;
  cout<<VXhydro[itau1][ix1][iy2][ieta1]<<endl;
  cout<<VXhydro[itau2][ix1][iy2][ieta1]<<endl;
  cout<<VXhydro[itau1][ix1][iy2][ieta2]<<endl;
  cout<<VXhydro[itau2][ix1][iy2][ieta2]<<endl;
  cout<<VXhydro[itau1][ix2][iy1][ieta1]<<endl;
  cout<<VXhydro[itau2][ix2][iy1][ieta1]<<endl;
  cout<<VXhydro[itau1][ix2][iy1][ieta2]<<endl;
  cout<<VXhydro[itau2][ix2][iy1][ieta2]<<endl;
  cout<<VXhydro[itau1][ix2][iy2][ieta1]<<endl;
  cout<<VXhydro[itau2][ix2][iy2][ieta1]<<endl;
  cout<<VXhydro[itau1][ix2][iy2][ieta2]<<endl;
  cout<<VXhydro[itau2][ix2][iy2][ieta2]<<endl;


cout<<"                         "<<endl;
cout<<"bulk_end--------------------------------------------------"<<endl;
cout<<"                         "<<endl;

*/










  
  VX = xs1*ys1*etas1*taus1*VXhydro[itau1][ix1][iy1][ieta1];
  VX+= xs1*ys1*etas1*taus2*VXhydro[itau2][ix1][iy1][ieta1];
  VX+= xs1*ys1*etas2*taus1*VXhydro[itau1][ix1][iy1][ieta2];
  VX+= xs1*ys1*etas2*taus2*VXhydro[itau2][ix1][iy1][ieta2];
  VX+= xs1*ys2*etas1*taus1*VXhydro[itau1][ix1][iy2][ieta1];
  VX+= xs1*ys2*etas1*taus2*VXhydro[itau2][ix1][iy2][ieta1];
  VX+= xs1*ys2*etas2*taus1*VXhydro[itau1][ix1][iy2][ieta2];
  VX+= xs1*ys2*etas2*taus2*VXhydro[itau2][ix1][iy2][ieta2];
  VX+= xs2*ys1*etas1*taus1*VXhydro[itau1][ix2][iy1][ieta1];
  VX+= xs2*ys1*etas1*taus2*VXhydro[itau2][ix2][iy1][ieta1];
  VX+= xs2*ys1*etas2*taus1*VXhydro[itau1][ix2][iy1][ieta2];
  VX+= xs2*ys1*etas2*taus2*VXhydro[itau2][ix2][iy1][ieta2];
  VX+= xs2*ys2*etas1*taus1*VXhydro[itau1][ix2][iy2][ieta1];
  VX+= xs2*ys2*etas1*taus2*VXhydro[itau2][ix2][iy2][ieta1];
  VX+= xs2*ys2*etas2*taus1*VXhydro[itau1][ix2][iy2][ieta2];
  VX+= xs2*ys2*etas2*taus2*VXhydro[itau2][ix2][iy2][ieta2];



  VY = xs1*ys1*etas1*taus1*VYhydro[itau1][ix1][iy1][ieta1];
  VY+= xs1*ys1*etas1*taus2*VYhydro[itau2][ix1][iy1][ieta1];
  VY+= xs1*ys1*etas2*taus1*VYhydro[itau1][ix1][iy1][ieta2];
  VY+= xs1*ys1*etas2*taus2*VYhydro[itau2][ix1][iy1][ieta2];
  VY+= xs1*ys2*etas1*taus1*VYhydro[itau1][ix1][iy2][ieta1];
  VY+= xs1*ys2*etas1*taus2*VYhydro[itau2][ix1][iy2][ieta1];
  VY+= xs1*ys2*etas2*taus1*VYhydro[itau1][ix1][iy2][ieta2];
  VY+= xs1*ys2*etas2*taus2*VYhydro[itau2][ix1][iy2][ieta2];
  VY+= xs2*ys1*etas1*taus1*VYhydro[itau1][ix2][iy1][ieta1];
  VY+= xs2*ys1*etas1*taus2*VYhydro[itau2][ix2][iy1][ieta1];
  VY+= xs2*ys1*etas2*taus1*VYhydro[itau1][ix2][iy1][ieta2];
  VY+= xs2*ys1*etas2*taus2*VYhydro[itau2][ix2][iy1][ieta2];
  VY+= xs2*ys2*etas1*taus1*VYhydro[itau1][ix2][iy2][ieta1];
  VY+= xs2*ys2*etas1*taus2*VYhydro[itau2][ix2][iy2][ieta1];
  VY+= xs2*ys2*etas2*taus1*VYhydro[itau1][ix2][iy2][ieta2];
  VY+= xs2*ys2*etas2*taus2*VYhydro[itau2][ix2][iy2][ieta2];
  
  VZ = xs1*ys1*etas1*taus1*VZhydro[itau1][ix1][iy1][ieta1];
  VZ+= xs1*ys1*etas1*taus2*VZhydro[itau2][ix1][iy1][ieta1];
  VZ+= xs1*ys1*etas2*taus1*VZhydro[itau1][ix1][iy1][ieta2];
  VZ+= xs1*ys1*etas2*taus2*VZhydro[itau2][ix1][iy1][ieta2];
  VZ+= xs1*ys2*etas1*taus1*VZhydro[itau1][ix1][iy2][ieta1];
  VZ+= xs1*ys2*etas1*taus2*VZhydro[itau2][ix1][iy2][ieta1];
  VZ+= xs1*ys2*etas2*taus1*VZhydro[itau1][ix1][iy2][ieta2];
  VZ+= xs1*ys2*etas2*taus2*VZhydro[itau2][ix1][iy2][ieta2];
  VZ+= xs2*ys1*etas1*taus1*VZhydro[itau1][ix2][iy1][ieta1];
  VZ+= xs2*ys1*etas1*taus2*VZhydro[itau2][ix2][iy1][ieta1];
  VZ+= xs2*ys1*etas2*taus1*VZhydro[itau1][ix2][iy1][ieta2];
  VZ+= xs2*ys1*etas2*taus2*VZhydro[itau2][ix2][iy1][ieta2];
  VZ+= xs2*ys2*etas1*taus1*VZhydro[itau1][ix2][iy2][ieta1];
  VZ+= xs2*ys2*etas1*taus2*VZhydro[itau2][ix2][iy2][ieta1];
  VZ+= xs2*ys2*etas2*taus1*VZhydro[itau1][ix2][iy2][ieta2];
  VZ+= xs2*ys2*etas2*taus2*VZhydro[itau2][ix2][iy2][ieta2]; 
  
  

}





void LBT(int &n, double &ti){
		
  int CT;                        //collision type 1 2 3 13 4 5 6 7 8
  int KATTC0;                    //flavor code 1/d 2/u 3/s -1/dbar -2/ubar -3/sbar 21/gluon 
  int KATT2;
  int KATT3;
  int KATT30;
		
  double RTE;                    //scattering rate (energy temperature)
  double striRTE;                    //scattering rate (energy temperature)

  double E;                      //parton energy
  double PLen;                   //parton momentum
  double T;                      //local temperature		

  double T1;                     //index for difference method
  double T2;
  double E1;
  double E2;
  int iT1;
  int iT2;
  int iE1;
  int iE2;	
		
  int nrad;
  int idlead;
  int idlead1;
  int idlead2;
  double Xtau;                   //....................main & titau	
  double Vx;
  double Vy;		
  double Veta;

  double tcar;                   //....................t x y z in Cartesian coordinate this is for the radiation process
  double xcar;		
  double ycar;
  double zcar;					

  double tcar0;                   //....................t x y z in Cartesian coordinate this is for the radiation process
  double xcar0;		
  double ycar0;
  double zcar0;

  
  double rans;
		
  int KATTx;
  int KATTy;

  double Ejp;
  double Elab;
		
  double qt;
		
  double px0;
  double py0;		
  double pz0;

  double Ncoll22;                    //...................average number of elastic scattering for particle i in the current step  
  int Nscatter;                   //...................number of elastic scattering for particle i in the current step

  int nnpp=np;                    //...................number of particles in the current np loop
  int np0=np;                     //...................number of particles in the current step (np0)
                                  //...................number of particles in the beginning of current step (np)

  int free=0;
  int free0=0;
  double fraction=0.0;
  double fraction0=0.0;
  double vc0b[4]={0.0};             //flow velocity     
  double pMag,vMag,flowFactor;

  double kt2;
  
  double vp0[4]={0.0};
  double p0temp[4]={0.0}; 

  double p0temp1=0.0; 
  double p0temp2=0.0; 
		
  double pcx[4]={0.0};
  double pcy[4]={0.0};		

  double pcx1[4]={0.0};
  double pcy1[4]={0.0};		

// for heavy quark
  double dt_lrf,maxFncHQ,Tdiff,lim_low,lim_high,lim_int;
  double probCol,probRad,probTot;

//...........................................................................................................NCUT
  ncut=0;
  ncut0=0;
  //...........................................................................................................NCUTEND		  
		  
  //...collision of the jet parton (i<=nj), recoiled partons, radiated gluons with the background thermal partons
  //...thermal partons
  //...np number of partons in current step
  //...nj number of jet partons 

  for(int i=1;i<=np;i++) {
//	cout<<"np:"<<"--------------------------------"<<np<<endl;
//		P[1][i]=0.0;
//		P[2][i]=0.0;
//		P[3][i]=sqrt(ener*ener-amss*amss);
//		P[0][i]=ener;
//		P[4][i]=amss;
//		P[5][i]=sqrt(P[1][i]*P[1][i]+P[2][i]*P[2][i]);
//		WT[i]=1.0;
//		KATT1[i]=Kjet;
      // if the production time of a parton is ahead of the current time, do nothing
    //  cout<<"i:"<<i<<"---------"<<"Vfrozen[0][i]:"<<Vfrozen[0][i]<<"---------------"<<"ti"<<ti<<endl;
      if(Vfrozen[0][i]>=ti) continue;  
//	cout<<"==================================="<<endl;
      icl22=0;
      qt=0;
      idlead=i;
		  
//......propagation of each positive particle and its negative counterpart in the current time step (for all particles)		
		
      if(P[0][i]<epsilon) {  //anti-overflow NOT NECESSARY
          V[1][i]=0.0;
          V[2][i]=0.0;
          V[3][i]=0.0;
          CAT[i]=1;		  
      }
      
      if(P0[0][i]<epsilon) {  //anti-overflow NOT NECESSARY
          V0[1][i]=0.0;
          V0[2][i]=0.0;
          V0[3][i]=0.0;
          CAT0[i]=1;
      }
      
 	CAT0[i]=1;
      if(CAT0[i] != 1) {	  
//...........................................t-z coordinates	
          if(tauswitch==0) {
		//cout<<"==================================================="<<endl;
              V0[1][i]=V0[1][i]+dt*P0[1][i]/P0[0][i];
              V0[2][i]=V0[2][i]+dt*P0[2][i]/P0[0][i];
	      V0[3][i]=V0[3][i]+dt*P0[3][i]/P0[0][i];

	      tcar0=ti;
	      zcar0=V0[3][i];
	      xcar0=V0[1][i];
	      ycar0=V0[2][i];

              double ed00,sd00,VX00,VY00,VZ00;
              int hydro_ctl0;



		  bulklinear(tcar0, xcar0, ycar0, zcar0, temp00, VX00, VY00, VZ00, hydro_ctl0);
	//	  cout<<"tcar0"<<" "<<tcar0<<" "<<"xcar0"<<" "<<xcar0<<" "<<"ycar0"<<" "<<ycar0<<" "<<"zcar0"<<" "<<zcar0<<" "<<endl;
		  
	//	  cout<<"temp00"<<" "<<temp00<<" "<<"hydro_ctl0"<<" "<<hydro_ctl0<<" "<<"VX00"<<" "<<VX00<<" "<<"VY00"<<" "<<VY00<<" "<<"VZ00"<<" "<<VZ00<<endl;

VX00=Vx_cal(VX00,zcar0,VZ00);
                 VY00=Vy_cal(VY00,zcar0,VZ00);
                 VZ00=Vz_cal(zcar0,VZ00);
                 
          //       cout<<"VX00"<<" "<<VX00<<" "<<"VY00"<<" "<<VY00<<" "<<"VZ00"<<" "<<VZ00<<endl;
		  
	      free0=0;

              if(bulkFlag==1) { // read OSU hydro
                  readhydroinfoshanshan_(&tcar0,&xcar0,&ycar0,&zcar0,&ed00,&sd00,&temp00,&VX00,&VY00,&VZ00,&hydro_ctl0);
              } else if(bulkFlag==2) { // read CCNU hydro
                 // hydroinfoccnu_(&tcar0, &xcar0, &ycar0, &zcar0, &temp00, &VX00, &VY00, &VZ00, &hydro_ctl0);
                 bulklinear(tcar0, xcar0, ycar0, zcar0, temp00, VX00, VY00, VZ00, hydro_ctl0);
                 VX00=Vx_cal(VX00,zcar0,VZ00);
                 VY00=Vy_cal(VY00,zcar0,VZ00);
                 VZ00=Vz_cal(zcar0,VZ00);
              } else if(bulkFlag==0) { // static medium
                  VX00=0.0;
                  VY00=0.0;
                  VZ00=0.0;
                  hydro_ctl0=0;
              }

              if(hydro_ctl0==0 && temp00>=hydro_Tc) {

	          qhat00=DebyeMass2(Kqhat0,alphas,temp00);
                  fraction0=1.0;

                  Vfrozen0[0][i]=ti;
                  Vfrozen0[1][i]=V0[1][i];
                  Vfrozen0[2][i]=V0[2][i];
                  Vfrozen0[3][i]=V0[3][i];
                  Tfrozen0[i]=temp00;
                  vcfrozen0[1][i]=VX00;
                  vcfrozen0[2][i]=VY00;
                  vcfrozen0[3][i]=VZ00;

              } else {
                  free0=1;
                  fraction0=0.0;
              }

	  }//if(tauswitch==0)
		  
      }//if(CAT0[i] != 1)

      if(CAT[i]!=1) {
  //    cout<<"-------------------1--------------------------"<<endl;
//...........................................t-z coordinates	
          if(tauswitch==0) {
//	cout<<"-------------------2--------------------------"<<endl;
	//cout<<"tcar"<<" "<<tcar<<" "<<"xcar"<<" "<<xcar<<" "<<"ycar"<<" "<<ycar<<" "<<"zcar"<<" "<<zcar<<" "<<endl;
//		cout<<"P[0][i]:"<<P[0][i]<<endl;
//		cout<<"V[1][i]:"<<"-----1-----"<<V[1][i]<<endl;
              V[1][i]=V[1][i]+dt*P[1][i]/P[0][i];
  //           cout<<"P[1][i]/P[0][i]:"<<P[1][i]/P[0][i]<<" "<<"dt:"<<dt<<endl;
    //         cout<<"V[1][i]:"<<"-----2-----"<<V[1][i]<<endl;
              
      //       cout<<"V[2][i]:"<<"-----1-----"<<V[2][i]<<endl;
              V[2][i]=V[2][i]+dt*P[2][i]/P[0][i];
        //     cout<<"P[2][i]/P[0][i]:"<<P[2][i]/P[0][i]<<" "<<"dt:"<<dt<<endl;
          //   cout<<"V[2][i]:"<<"-----2-----"<<V[2][i]<<endl;
              
            // cout<<"V[3][i]:"<<"-----1-----"<<V[3][i]<<endl;
              V[3][i]=V[3][i]+dt*P[3][i]/P[0][i];
            // cout<<"P[3][i]/P[0][i]:"<<P[3][i]/P[0][i]<<" "<<"dt:"<<dt<<endl;
             // cout<<"V[3][i]:"<<"-----2-----"<<V[3][i]<<endl;
cout<<"V:"<<V[1][i]<<V[1][2]<<V[3][i]<<endl;

	      tcar=ti;
	      zcar=V[3][i];
	      xcar=V[1][i];
	      ycar=V[2][i];
	     // cout<<"tcar"<<" "<<tcar<<" "<<"xcar"<<" "<<xcar<<" "<<"ycar"<<" "<<ycar<<" "<<"zcar"<<" "<<zcar<<" "<<endl;

          //...........................................Hydro part		

              double ed,sd,VX,VY,VZ;
              int hydro_ctl;

/*
		  bulklinear(tcar, xcar, ycar, zcar, temp0, VX, VY, VZ, hydro_ctl);
		  cout<<"tcar"<<" "<<tcar<<" "<<"xcar"<<" "<<xcar<<" "<<"ycar"<<" "<<ycar<<" "<<"zcar"<<" "<<zcar<<" "<<endl;
		  
		  cout<<"temp0"<<" "<<temp0<<" "<<"hydro_ctl"<<" "<<hydro_ctl<<" "<<"VX"<<" "<<VX<<" "<<"VY"<<" "<<VY<<" "<<"VZ"<<" "<<VZ<<endl;

		 VX=Vx_cal(VX,zcar,VZ);
                 VY=Vy_cal(VY,zcar,VZ);
                 VZ=Vz_cal(zcar,VZ);
                 
                 cout<<"VX"<<" "<<VX<<" "<<"VY"<<" "<<VY<<" "<<"VZ"<<" "<<VZ<<endl;
 */                
	      free=0;

              if(free==0) {

                  if(bulkFlag==1) { // read OSU hydro
                      readhydroinfoshanshan_(&tcar,&xcar,&ycar,&zcar,&ed,&sd,&temp0,&VX,&VY,&VZ,&hydro_ctl);
                  } else if(bulkFlag==2) { // read CCNU hydro
                    //  cout<<"tcar"<<" "<<tcar<<" "<<"xcar"<<" "<<xcar<<" "<<"ycar"<<" "<<ycar<<" "<<"zcar"<<" "<<zcar<<" "<<endl;
                      bulklinear(tcar, xcar, ycar, zcar, temp0, VX, VY, VZ, hydro_ctl);
                    //  cout<<"tcar"<<" "<<tcar<<" "<<"xcar"<<" "<<xcar<<" "<<"ycar"<<" "<<ycar<<" "<<"zcar"<<" "<<zcar<<" "<<endl;
                     // cout<<"temp0"<<" "<<temp0<<" "<<"hydro_ctl"<<" "<<hydro_ctl<<" "<<"VX"<<" "<<VX<<" "<<"VY"<<" "<<VY<<" "<<"VZ"<<" "<<VZ<<endl;
                  VX=Vx_cal(VX,zcar,VZ);
              //    cout<<"VX:"<<VX<<endl;
                 VY=Vy_cal(VY,zcar,VZ);
               //  cout<<"VY:"<<VY<<endl;
                 VZ=Vz_cal(zcar,VZ);
               //  cout<<"VZ:"<<VZ<<endl;
                // cout<<"VX"<<" "<<VX<<" "<<"VY"<<" "<<VY<<" "<<"VZ"<<" "<<VZ<<endl;
                  } else if(bulkFlag==0) { // static medium
                      VX=0.0;
                      VY=0.0;
                      VZ=0.0;
                      hydro_ctl=0;
                  }

                  if(hydro_ctl==0 && temp0>=hydro_Tc) {

                      vc0[1]=VX;
                      vc0[2]=VY;
                      vc0[3]=VZ;

//                      vc0[1]=0.0;
//                      vc0[2]=0.0;
//                      vc0[3]=0.0;

                      vc0b[1]=VX;
                      vc0b[2]=VY;
                      vc0b[3]=VZ;

        
                      //...alphas!	
                      alphas=alphas0(Kalphas,temp0);
                      //...Debye Mass square
                      qhat0=DebyeMass2(Kqhat0,alphas,temp0);	
	
                      fraction=1.0;
                      Vfrozen[0][i]=ti;
                      Vfrozen[1][i]=V[1][i];
                      Vfrozen[2][i]=V[2][i];
                      Vfrozen[3][i]=V[3][i];
                      Tfrozen[i]=temp0;
                      vcfrozen[1][i]=VX;
                      vcfrozen[2][i]=VY;
                      vcfrozen[3][i]=VZ;

                  } else {
                      free=1;
                      fraction=0.0;
                  }
  
       	          pc0[1]=P[1][i];
                  pc0[2]=P[2][i];
       	          pc0[3]=P[3][i];
       	          pc0[0]=P[0][i];

                  vMag=sqrt(vc0b[1]*vc0b[1]+vc0b[2]*vc0b[2]+vc0b[3]*vc0b[3]);
                  flowFactor=(1.0-(pc0[1]*vc0b[1]+pc0[2]*vc0b[2]+pc0[3]*vc0b[3])/pc0[0])/sqrt(1.0-vMag*vMag);
 
	          trans(vc0,pc0);
			  
	          if(pc0[0]<sqrt(qhat0)) {
	              free=1;
	          }
	  
	          if(i>nj && pc0[0]<Ecmcut) {
	              free=1;
                  }   
		  
	          transback(vc0,pc0);

              }//if(free==0)
	  } //if(tauswitch==0)
				

//...........................................tau-eta coordinates		
	  if(tauswitch==1) {
		
		//????????????? why I need to get the medium information before propagation
		
          //...........................................Hydro part		

		  /////////////////////////////////////////////////////////////
		  //           GetInfo( CD & tau, CD & x, CD &y, CD &etas, double & Ed, double & T, double & Vx, double & Vy, double & Veta, double & phase_ratio )        
		  /////////////////////////////////////////////////////////////

          //...........................................Hydro part end
		  
		  vp[0]=0.0;
		  vp[1]=V[1][i];
		  vp[2]=V[2][i];
		  vp[3]=V[3][i];

		  pc0[1]=P[1][i];
		  pc0[2]=P[2][i];
		  pc0[3]=P[3][i];
		  pc0[0]=P[0][i];
		  
		  titau(ti,vc0,vp,pc0,Vx,Vy,Veta,Xtau);

		  V[1][i]=V[1][i]+dt*Vx;
		  V[2][i]=V[2][i]+dt*Vy;
		  V[3][i]=V[3][i]+dt*Veta;

		  ///////////////////////////////////////time and position in Cartesian coordinate			 
		  tcar=ti*cosh(V[3][i]);
		  zcar=ti*sinh(V[3][i]);
		  xcar=V[1][i];
		  ycar=V[2][i];

                  alphas=alphas0(Kalphas,temp0);
		  qhat0=DebyeMass2(Kqhat0,alphas,temp0);


		  if(i>nj)
		  {
		  //......................................negative particle propagation
		  
		  /////////////////////////////////////////////////////////////
		  //           GetInfo( CD & tau, CD & x, CD &y, CD &etas, double & Ed, double & T, double & Vx, double & Vy, double & Veta, double & phase_ratio )        
		  /////////////////////////////////////////////////////////////
		  
		  vp0[0]=0.0;
		  vp0[1]=V0[1][i];
		  vp0[2]=V0[2][i];
		  vp0[3]=V0[3][i];

		  pc00[1]=P0[1][i];
		  pc00[2]=P0[2][i];
		  pc00[3]=P0[3][i];
		  pc00[0]=P0[0][i];
		  
		  titau(ti,vc0,vp,pc00,Vx,Vy,Veta,Xtau);

		  V0[1][i]=V0[1][i]+dt*Vx;
		  V0[2][i]=V0[2][i]+dt*Vy;
		  V0[3][i]=V0[3][i]+dt*Veta;

		  ///////////////////////////////////////time system in Cartesian coordinate			 
		  tcar=ti*cosh(V0[3][i]);
		  zcar=ti*sinh(V0[3][i]);
		  xcar=V0[1][i];
		  ycar=V0[2][i];
			
//                temp00= temp0;
			
//		  qhat00=DebyeMass2(Kqhat0,alphas,temp00);
		  
		  } // if(i>nj)
          } // if(tauswitch==1)			

//......propagation of each positive particle and its negative counterpart in the current time step	end	

//......CAT[i]=1 free streaming particle

//......free streaming when energy  < Ecut
//......free streaming when energy0 < Ecut0......in the next step

          if(free==0) {

              double qhatTP,RTE1,RTE2;
              double striqhatTP,striRTE1,striRTE2;

// modification: interplotate rate in l.r.f. 
       	      pc0[1]=P[1][i];
       	      pc0[2]=P[2][i];
       	      pc0[3]=P[3][i];
       	      pc0[0]=P[0][i];
              trans(vc0,pc0);
              E=pc0[0]; //  p4-the initial 4-momentum of the jet parton in the local rest frame
              PLen=sqrt(pc0[1]*pc0[1]+pc0[2]*pc0[2]+pc0[3]*pc0[3]);
              transback(vc0,pc0);

	      T=temp0;
	      KATTC0=KATT1[i];

//              if(E<T) preKT=4.0*pi/9.0/log(scaleAK*T*T/0.04)/alphas/fixKT;
//              else preKT=4.0*pi/9.0/log(scaleAK*E*T/0.04)/alphas/fixKT;

              runKT=4.0*pi/9.0/log(2.0*E*T/0.04)/0.3;
              preKT=fixAlphas/0.3;

	      lam(KATTC0,RTE,striRTE,PLen,T,T1,T2,E1,E2,iT1,iT2,iE1,iE2); //modified: use P instead

//              preKT=alphas/0.3;

// calculate p,T-dependence K factor
              KPfactor=1.0+KPamp*exp(-PLen*PLen/2.0/KPsig/KPsig);
              KTfactor=1.0+KTamp*exp(-pow((temp0-hydro_Tc),2)/2.0/KTsig/KTsig);
//              Kfactor=KPfactor*KTfactor*KTfactor*preKT*preKT;      
              Kfactor=KPfactor*KTfactor*KTfactor*runKT*preKT;      

// get qhat from table
              if(KATTC0==21) {
                  RTE1=(qhatG[iT2][iE1]-qhatG[iT1][iE1])*(T-T1)/(T2-T1)+qhatG[iT1][iE1];
                  RTE2=(qhatG[iT2][iE2]-qhatG[iT1][iE2])*(T-T1)/(T2-T1)+qhatG[iT1][iE2];
              } else if(KATTC0==4||KATTC0==-4) {
                  RTE1=(qhatHQ[iT2][iE1]-qhatHQ[iT1][iE1])*(T-T1)/(T2-T1)+qhatHQ[iT1][iE1];
                  RTE2=(qhatHQ[iT2][iE2]-qhatHQ[iT1][iE2])*(T-T1)/(T2-T1)+qhatHQ[iT1][iE2];

                  striRTE1=(striqhatHQ[iT2][iE1]-striqhatHQ[iT1][iE1])*(T-T1)/(T2-T1)+striqhatHQ[iT1][iE1];
                  striRTE2=(striqhatHQ[iT2][iE2]-striqhatHQ[iT1][iE2])*(T-T1)/(T2-T1)+striqhatHQ[iT1][iE2];
              } else {
                  RTE1=(qhatLQ[iT2][iE1]-qhatLQ[iT1][iE1])*(T-T1)/(T2-T1)+qhatLQ[iT1][iE1];
                  RTE2=(qhatLQ[iT2][iE2]-qhatLQ[iT1][iE2])*(T-T1)/(T2-T1)+qhatLQ[iT1][iE2];
              }

              qhatTP=(RTE2-RTE1)*(PLen-E1)/(E2-E1)+RTE1;
              striqhatTP=(striRTE2-striRTE1)*(PLen-E1)/(E2-E1)+striRTE1;

              qhatTP=qhatTP*Kfactor+striqhatTP;

              qhat_over_T3=qhatTP;  // what is read in is qhat/T^3 of quark/gluon

              // D2piT is diffusion coefficient "just for quark"
              if(KATTC0==21) D2piT=8.0*pi/(qhat_over_T3/CA*CF); else D2piT=8.0*pi/qhat_over_T3;

              qhat=qhatTP*pow(T,3); // for light quark and gluon, need additional table to make it correct

// update interaction time and average gluon number for heavy quark radiation 
// (as long as it's in medium, even though no collision)
//              dt_lrf=dt*sqrt(1.0-vc0[1]*vc0[1]-vc0[2]*vc0[2]-vc0[3]*vc0[3]);
              dt_lrf=dt*flowFactor;
              Tint_lrf[i]=Tint_lrf[i]+dt_lrf;
              Tdiff=Tint_lrf[i];

// get radiation probablity by reading tables -- same for heavy and light partons
              if(KATTC0==21) radng[i]+=nHQgluon(KATT1[i],dt_lrf,Tdiff,temp0,E,maxFncHQ)*KTfactor/2.0*runKT;
              else radng[i]+=nHQgluon(KATT1[i],dt_lrf,Tdiff,temp0,E,maxFncHQ)*KTfactor*runKT;
//              if(KATTC0==21) radng[i]+=nHQgluon(KATT1[i],dt_lrf,Tdiff,temp0,E,maxFncHQ)/2.0*KTfactor*preKT;
//              else radng[i]+=nHQgluon(KATT1[i],dt_lrf,Tdiff,temp0,E,maxFncHQ)*KTfactor*preKT;
              lim_low=sqrt(6.0*pi*alphas)*temp0/E;
//              lim_low=sqrt(6.0*pi*0.3)*temp0/E;
              if(abs(KATT1[i])==4) lim_high=1.0;
              else lim_high=1.0-lim_low;
              lim_int=lim_high-lim_low;
              if(lim_int>0.0) probRad=1.0-exp(-radng[i]);
              else probRad=0.0;

//...........................................t-z coordinates
              if(tauswitch==0) {

		  pc0[1]=P[1][i];
		  pc0[2]=P[2][i];
		  pc0[3]=P[3][i];
		  pc0[0]=P[0][i];
		  
//                V[0][i]=V[0][i]-dt*RTE/0.1970*sqrt(pow(P[1][i],2)+pow(P[2][i],2)+pow(P[3][i],2))/P[0][i];	
//		  Ncoll22=dt*RTE/0.197*((pc0[0]*1-pc0[1]*vc0[1]-pc0[2]*vc0[2]-pc0[3]*vc0[3])/E);		  

//??????@@@
//////////////////
                  V[0][i]=V[0][i]-fraction*dt*RTE/0.1970*sqrt(pow(P[1][i],2)+pow(P[2][i],2)+pow(P[3][i],2))/P[0][i];
                  probCol=fraction*dt_lrf*RTE/0.1970;
              }
//...........................................tau-eta coordinates
              if(tauswitch==1) {

                  V[0][i]=V[0][i]-dt*RTE*Xtau/0.1970*sqrt(pow(P[1][i],2)+pow(P[2][i],2)+pow(P[3][i],2))/P[0][i];
                  probCol=dt_lrf*RTE*Xtau/0.1970;  // ?? fraction
//		  Ncoll22=dt*RTE/0.197*Xtau*((pc0[0]*1-pc0[1]*vc0[1]-pc0[2]*vc0[2]-pc0[3]*vc0[3])/E);  		  
              }


//          cout<<"t, T, PLen, rate, lambda:  "<<ti<<"  "<<T<<"  "<<PLen<<"  "<<RTE<<"  "<<1/RTE*0.197<<endl;
////////////////////////////////////////...new elastic part	  
 
//          if(ran0(&NUM1)<1-exp(-Ncoll22)) // !Yes, collision!
////////////////////////////////////////		  

              if(KINT0==0) probRad=0.0;  // switch off radiation
//              probCol=probCol*KPfactor*KTfactor*preKT;
              probCol=probCol*preKT*runKT+fraction*dt_lrf*striRTE/0.1970;
              probCol=(1.0-exp(-probCol))*(1.0-probRad);  // probability of pure elastic scattering
              if(KINT0==2) probCol=0.0;
              probTot=probCol+probRad;

	              		  
              if(ran0(&NUM1)<probTot) { // !Yes, collision! Either elastic or inelastic.

//			  Nscatter=KPoisson(Ncoll22);

                  Nscatter=1; 
			  
                  for(int nsca=1;nsca<=Nscatter;nsca++) {			
			
//                          cout << "scatter!" << endl;
			  np0=np0+1;
			  pc0[1]=P[1][i];
			  pc0[2]=P[2][i];
			  pc0[3]=P[3][i];
			  pc0[0]=P[0][i];

// use PLen instead			  
                          RTE=RTE*preKT*runKT+striRTE;
			  flavor(CT,KATTC0,KATT2,KATT3,RTE,PLen,T,T1,T2,E1,E2,iT1,iT2,iE1,iE2);	
			 
			  //...........collision between a jet parton or a recoiled parton and a thermal parton from the medium 			 

                          if(CT==11) { // for Qq->Qq scattering
			      collHQ22Qq(CT,temp0,qhat0,vc0,pc0,pc2,pc3,pc4,qt);
                          } else if (CT==12) { // for Qg->Qg scattering
			      collHQ22Qg(CT,temp0,qhat0,vc0,pc0,pc2,pc3,pc4,qt);
                          } else { // for light parton scattering
			      colljet22(CT,temp0,qhat0,vc0,pc0,pc2,pc3,pc4,qt);
                          }

                          icl22 = 1;
                          n_coll22+=1;
			 
			  KATT1[i]=KATTC0;
			  KATT1[np0]=KATT2;
			  KATT10[np0]=KATT3;

			  if(pc0[0]<pc2[0] && abs(KATTC0)!=4 && KATTC0==KATT2) { //disable switch for heavy quark, only allow switch for identical particles
			      for(int k=0;k<=3;k++) {
				  p0temp[k]=pc2[k];
				  pc2[k]=pc0[k];
				  pc0[k]=p0temp[k];
			      }
			      KATT1[i]=KATT2;
     			      KATT1[np0]=KATTC0;
			  }

			  for(int j=0;j<=3;j++) {
//			    if(i<=nj) // always make sure the jet parton carries the largest momentum (maybe now all the parton includes recoiled for there could be multiple scattering in a single time step)
				
				P[j][i]=pc0[j];
				P[j][np0]=pc2[j];
				V[j][np0]=V[j][i];

				
//		        if(i<=nj || nsca!=1) // then for each recoiled parton there is a negative parton
//				{
				P0[j][np0]=pc3[j];
				V0[j][np0]=V[j][i];

                                Vfrozen[j][np0]=V[j][i];
                                Vfrozen0[j][np0]=V[j][i];
                                Tfrozen[np0]=temp0;
                                Tfrozen0[np0]=temp0;

                                if(j!=0) {                                 
                                    vcfrozen[j][np0]=vc0b[j];
                                    vcfrozen0[j][np0]=vc0b[j];
                                }
//				}
			  }
			  //..........................test

			//	cout<<"i:"<<i<<"-------------1----------- "<<"P[0][i]:"<<P[0][i]<<endl;
			//	cout<<"np0:"<<np0<<" "<<"P[0][np0]:"<<P[0][np0]<<endl;

                          // pass initial pT and weight information from jet parton to recoil and negative parton
                          P[5][np0]=P[5][i];
                          P0[5][np0]=P[5][i];
                          WT[np0]=WT[i];
                          WT0[np0]=WT[i];
                          P[4][np0]=0.0;  // recoiled parton is always massless
                          P0[4][np0]=0.0;


//CAT!!!			  
/*	       				
*/
                          CAT[np0]=2;

//.............................................................. collision between "negative" partons

//...background thermal parton scattering with the same thermal parton
//...the associated parton has scattered

//...negative parton can only get the rank i>nj and it can only act on the first collision in the current time step. 
//...This is not physical only because the recoiled parton and its negative partner appear in pairs 

//...if Ecut0
//

//.............................................................. collision between "negative" partons

//...background thermal parton scattering with the same thermal parton the associated parton has scattered
//...negative parton can only get the rank i>nj and it can only act on the first collision in the current time step. 
//...This is not physical only because the recoiled parton and its negative partner appear in pairs 

////...Scattering between P0[i] and the newest negative parton  // commented out by Shanshan, just let negative particles stream freely
//
//            if(CAT0[i] == 0)
//            {
////???
////                free0=1;
////???                
//                if(free0 == 0)
//                {
//		    if(i>nj && nsca==1)
//		    {
//
//                        temp00=temp0;
//
//		        qhat00=DebyeMass2(Kqhat0,alphas,temp00);
//
//                        pc00[1]=P0[1][i];
//                        pc00[2]=P0[2][i];
//                        pc00[3]=P0[3][i];
//                        pc00[0]=P0[0][i];
//				
//			pc30[1]=pc3[1];
//			pc30[2]=pc3[2];
//			pc30[3]=pc3[3];
//			pc30[0]=pc3[0];
//
//	                E=pc00[0];
//	                T=temp00;
//	                KATTC0=KATT10[i];
//			KATT30=KATT3;
//        
//			twflavor(CT,KATTC0,KATT30,E,T);		 			
//					   				   
//                        twcoll(CT,qhat00,vc0,pc00,pc30);		   
//				   				   				   
//	                KATT1[np0]=KATT2;
//	                KATT10[i]=KATTC0;
//	                KATT10[np0]=KATT30;
//
//                        for(int j=0;j<=3;j++)
//			{
//                            P0[j][i]=pc00[j];
//                            P0[j][np0]=pc30[j];
//                            V0[j][np0]=V[j][i];
//                        }
//               	        V0[0][np0]=-log(1.0-ran0(&NUM1));
//		    } //if(i>nj && nsca==1)
//                } // if(free0==0)
//            } //.............................................................. "negative" end
											
             
	          }//for(int nsca=1;nsca<=Nscatter;nsca++)			  			  
			  
//...............................................................................................the point where the elastic scattering end			  

                  if(qt<epsilon) { 
                      KINT=0;
                  } else KINT=1;
        
		  if(qt>epsilon)
		  {
			qperp[nqperp]=qt;
			nqperp=nqperp+1;
		  }
        
                  if(KINT0!=0 && KINT!=0) { // Switch on the radiation processes
                      for(int j=0;j<=3;j++) {
                          p0temp[j]=P[j][i];
                      }							
                	
                      px0=pc4[1];
                      py0=pc4[2];
                      pz0=pc4[3];
        
        	      Elab=pc4[0]; // p4-the initial 4-momentum of the jet parton in the lab frame
                      trans(vc0,pc4);
                      Ejp=pc4[0]; //  p4-the initial 4-momentum of the jet parton in the local rest frame
                      transback(vc0,pc4);
        
                      if(Ejp>2*sqrt(qhat0)) { // !radiation or not?
        
        		  Vtemp[1]=xcar;
        		  Vtemp[2]=ycar;
        		  Vtemp[3]=zcar;
        
                          if(KATT1[i] != 21) {
                              isp=1;
                              n_sp1 +=1;
                          } else {
                              isp=2;
                              n_sp2 +=1;
                          }
         
                          if(ran0(&NUM1)<probRad/probTot) { // radiation -- either heavy or light parton:
        
        	              for(int j=0;j<=3;j++) { // prepare for the collHQ23
        			  pc01[j]=pc4[j]; // 2->3 process
        			  pb[j]=pc4[j];
                      	      }
        
        //		      colljet23(temp0,qhat0,vc0,pc01,pc2,pc3,pc4,qt,icl23,tcar,tiscatter[i],tirad[i],Elab,Ejp);
        		      collHQ23(KATT1[i],temp0,qhat0,vc0,pc01,pc2,pc3,pc4,qt,icl23,Tdiff,Ejp,maxFncHQ,lim_low,lim_int);
                              n_coll23+=1;
        
        		      if(icl23!=1) {
        
                                  int ctGluon=1;
        
                                  ng_coll23+=1; 
        			  nrad=KPoisson(radng[i]);       
			       	cout<<"nrad"<<nrad<<endl;	  
                                  int nrad0=nrad;
        
        //                          nrad=1; // only allow single gluon emission right now
                                  
                                  ng_nrad+=nrad;
        			  np0=np0+1;
        			  KATT1[np0]=21;				  
        			  
                                  for(int j=0;j<=3;j++) {
            			      P[j][i]=pc01[j];            // heavy quark after colljet23
        			      P[j][np0-1]=pc2[j];         // replace the final state of 22
        			      V[j][np0-1]=V[j][i];
        			      P0[j][np0-1]=pc3[j];
        			      V0[j][np0-1]=V[j][i];
        			      P[j][np0]=pc4[j];           //radiated gluon from colljet23
        			      V[j][np0]=V[j][i];
        			      P0[j][np0]=0.0;
        			      V0[j][np0]=0.0;

				
                                      Vfrozen[j][np0]=V[j][i];
                                      Vfrozen0[j][np0]=0.0;
                                      if(j!=0) {                                 
                                          vcfrozen[j][np0]=vc0b[j];
                                          vcfrozen0[j][np0]=0.0;
                                      }
        			  }
        			  //..........................test
				
			//	cout<<"i:"<<i<<"--------------2------------ "<<"P[0][i]:"<<P[0][i]<<endl;
			//	cout<<"np0:"<<np0<<" "<<"P[0][np0]:"<<P[0][np0]<<endl;
			//	cout<<"np0-1:"<<np0-1<<" "<<"P[0][np0-1]:"<<P[0][np0-1]<<endl;
        
                                  Tfrozen[np0]=temp0;
                                  Tfrozen0[np0]=0.0;
        
                                  // pass initial pT and weight information
                                  P[5][np0]=P[5][i];
                                  P0[5][np0]=0.0;
                                  WT[np0]=WT[i];
                                  WT0[np0]=0.0;   // doesn't exist
                                  P[4][np0]=0.0;  // radiated gluons are massless
                                  P0[4][np0]=0.0;
        
                                  Tint_lrf[i]=0.0;  //reset radiation infomation for heavy quark
        
                                  eGluon=eGluon+pc4[0];
                                  nGluon=nGluon+1.0;
        
                       	          V[0][np0]=-log(1.0-ran0(&NUM1));
        
                                  // add multiple radiation for heavy quark
                                  while(nrad>1) {
        									
                                      for(int j=0;j<=3;j++) pc2[j]=pc4[j];
        
        //			      radiation(qhat0,vc0,pc01,pc2,pc4,pb,iclrad,tcar,tiscatter[i],tirad[i],Elab,Ejp);		   
        			      radiationHQ(KATT1[i],qhat0,vc0,pc2,pc01,pc4,pb,iclrad,Tdiff,Ejp,maxFncHQ,temp0,lim_low,lim_int);  
                                      n_radiation+=1;
        
        			      if(iclrad!=1) {
        			          np0=np0+1;
        			       	  ng_radiation+=1;
        			          KATT1[np0]=21;					  
        
        //                                  cout << "add one gluon" << endl;
                                          ctGluon++;
        
                                          for(int j=0;j<=3;j++) {
                    			      P[j][i]=pc01[j];            // heavy quark after one more gluon
                			      P[j][np0-1]=pc2[j];         // update the previous gluon
                			      P[j][np0]=pc4[j];           // record the new gluon
                			      V[j][np0]=V[j][i];
                			      P0[j][np0]=0.0;
                			      V0[j][np0]=0.0;
                			     
			
        
        
                                              Vfrozen[j][np0]=V[j][i];
                                              Vfrozen0[j][np0]=0.0;
                                              if(j!=0) {                                 
                                                  vcfrozen[j][np0]=vc0b[j];
                                                  vcfrozen0[j][np0]=0.0;
                                              }
                                          }
                                           //..........................test
                	//	cout<<"i:"<<i<<"-----------------------3------------------- "<<"P[0][i]:"<<P[0][i]<<endl;
			//	cout<<"np0:"<<np0<<" "<<"P[0][np0]:"<<P[0][np0]<<endl;
			//	cout<<"np0-1:"<<np0-1<<" "<<"P[0][np0-1]:"<<P[0][np0-1]<<endl;
        
                                          Tfrozen[np0]=temp0;
                                          Tfrozen0[np0]=0.0;
         
                                          // pass initial pT and weight information
                                          P[5][np0]=P[5][i];
                                          P0[5][np0]=0.0;
                                          WT[np0]=WT[i];
                                          WT0[np0]=0.0;   // doesn't exist
                                          P[4][np0]=0.0;  // radiated gluons are massless
                                          P0[4][np0]=0.0;
        
                                          eGluon=eGluon+pc4[0];
                                          nGluon=nGluon+1.0;
        
                              	          V[0][np0]=-log(1.0-ran0(&NUM1));
        
                                      } else { //end multiple radiation
                                          break;
                                      }  //if(icl!=1)radiation
        		   
                                      nrad=nrad-1;
        			  } //multiple radiation for heavy quark ends
        
        //                          cout<<"radiate! <Ng>: "<<nrad0<<"  real Ng H: "<< ctGluon << endl;
                              } // icl23 == 1, 1st gluon radiation from heavy quark
        
                          } // if(ran0(&NUM1)<probRad/probTot) 
        
        	      }  //if(Ejp>2*sqrt(qhat0))				  
        				  				  
        	  }  //if(KINT!=0)
        
        
        	  //...........collision end				  
        	  //...........determine the leading partons and set radng[i]
                      
        
                  //....tiscatter information
                  if(abs(KATT1[i])==4) radng[i]=0.0; // do it below
                  tiscatter[i]=tcar;
                  V[0][i]=-log(1.0-ran0(&NUM1));
        
        //          for(unsigned ip=np+1; ip<=np0; ++ip)
                  for(unsigned ip=nnpp+1; ip<=np0; ++ip)
                  {
                        tiscatter[ip]=tcar;
                    	V[0][ip]=-log(1.0-ran0(&NUM1));
                    	V0[0][ip]=-log(1.0-ran0(&NUM1));
                        tirad[ip]=tcar;
                        Tint_lrf[ip]=0.0;
                        radng[ip]=0.0;
                  }
        
        // SC: only do possible switch for g->ggg... here
                  if(KATT1[i]==21 && np0-nnpp>=2) {
                      idlead=nnpp+2;
                      p0temp1=P[0][idlead];
                      for(unsigned ip=nnpp+2; ip<=np0-1; ++ip) {
                          if(p0temp1 < P[0][ip+1]) {
                              p0temp1=P[0][ip+1];
                              idlead=ip+1;
                          }
                      }
        	      if(P[0][idlead]>P[0][i]) {
        	          KATTx=KATT1[i];
        	          KATTy=KATT1[idlead];
        	          KATT1[i]=KATTy;
        	          KATT1[idlead]=KATTx;					
                          for(int j=0;j<=3;j++) {
        	              pcx[j]=P[j][i];
        	              pcy[j]=P[j][idlead];
        	              P[j][i]=pcy[j];
        	              P[j][idlead]=pcx[j];
        	          }
                      }
                  }

	  //...........sorting block i np~np0(positive) np+1(negative)

	  //........................................................................................................

//CAT!!!			  
/*			  
*/			  
                  for(int m=nnpp+1;m<=np0;m++) {
                      if(CAT[i]==2) {
                          CAT[m]=2;             				 
                      }
                  }

                  if(P[0][nnpp+1]<Ecut) {
//                      cout << "cut recoil   " << nnpp << "    " << P[0][nnpp+1] << endl; 
                      CAT[nnpp+1]=1;
                      CAT0[nnpp+1]=1;
                      for(int j=0; j<=3; j++) {
                          P[j][nnpp+1]=0.0;
                          P0[j][nnpp+1]=0.0;
                      }
                  }

                  if(P[0][i]<Ecut) {
//                      cout << "cut leading" << endl;
                      CAT[i]=1;
                      for(int j=0; j<=3; j++) P[j][i]=0.0;
                  }

                  for(int m=nnpp+2; m<=np0; m++) {
                      if(P[0][m]<Ecut) {
//                          cout << "cut radiation   " << m << "    " << P[0][m] << endl;
                          CAT[m]=1;
                          for(int j=0; j<=3; j++) P[j][m]=0.0;
                      }
                  }

			  
//CAT!!!			  
/*			  
			  if((i>nj) && (P[0][i]<Ecut))
				{
				  CAT[i]=1;
				  ncut=ncut+1;
//				  for(int j=0;j<=3;j++)
//					{
//					  PP[j][ncut]=P[j][i];
//					  VV[j][ncut]=V[j][i];
//					  P[j][i]=0.0;
//					}
				}			 			 
*/
								
				
//		     cout << i << " " << np << endl;
//             for(unsigned ip=1; ip<=np0; ++ip)
//             {
//                 cout << ip << " " << P[0][ip] << endl;
//	     }
//		     cout << endl;
 
			  //........................................................................................................sorting block for p2 and radiated gluons
              
//CAT!!!			  
/*				  
			  for(int m=np+1;m<=np0;m++)
				{
				  if(P[0][m]<=Ecut)
					{			
					  CAT[m]=1;             
					  ncut=ncut+1;
//					  for(int j=0;j<=3;j++)
//						{
//						  PP[j][ncut]=P[j][m];
//						  VV[j][ncut]=V[j][m];
//						  P[j][m]=0.0;			 
//						}				 
					}
				}
*/

			  //........................................................................................................dump all the negative particles
			  //..................................................................CAT0=0/1  stay in LBT/dumped into hydro (negative partons)
//			  CAT0[np0]=1;             
//			  ncut0=ncut0+1;
//			  for(int j=0;j<=3;j++)
//                {
//				  PP0[j][ncut0]=P0[j][np+1];
//				  VV0[j][ncut0]=V0[j][np+1];
//				  P0[j][np0]=0.0;				 
//                }				
			  //...........................................................................................................++++++++++++++++++++++++END			 				
              } //...!Yes, collision!			

// even if there is no scattering, one should reset last scatter time now in the new framework
              tiscatter[i]=tcar;
              radng[i]=0.0;

          } //....for if(free==0)
      } //..........if CAT[i]=0 end	

//      } //..........if tcar<tiform[i]

      nnpp=np0;  

  } //......for np end

  //........time step end, np: the number of particles at this point  		
  np=np0;
  
  if(Kprimary==1) {
     np=nj;
  }

  if(np>=499999) {
      cout << "np exceeds the grid size ... terminate " << endl;
      exit (EXIT_FAILURE); 
  }

}
		
		
////////////////////////////////////////////////////////////////////////////////////////////////////



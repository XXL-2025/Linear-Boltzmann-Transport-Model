     char Hydroprofile_bulk3D_route[1024];
     sprintf(Hydroprofile_bulk3D_route,"hydroProfile/bulk.dat");     	
     ifstream f2(Hydroprofile_bulk3D_route);
	 cout<<Hydroprofile_bulk3D_route<<endl;

     //cout<<"ntauhydro"<<" "<<ntauhydro<<" "<<"nxhydro"<<" "<<nxhydro<<" "<<"nyhydro"<<" "<<nyhydro<<" "<<"netahydro"<<" "<<netahydro<<endl;
     int ntauhydro=31;
     int nxhydro=73;
     int nyhdro=73;
     int netahydro=33;
     double tauh,ed;
	 
     //ifstream f2("readindatafile/Hydroprofile/bulk3D.dat");
     if(!f2.is_open())
     {
      cout<<"Erro openning date file2!\n";
     }
     else
     {
      for(int itau=1;itau<=ntauhydro;itau++)
      {
      for(int ix=1;ix<=nxhydro;ix++)
      {
      for(int iy=1;iy<=nyhydro;iy++)
      {
      for(int ieta=1;ieta<=netahydro;ieta++)
      {
	  
	  f2>>tauh>>ed>>temphydro[itau][ix][iy][ieta]>>fractionhydro[itau][ix][iy][ieta]>>VXhydro[itau][ix][iy][ieta]>>VYhydro[itau][ix][iy][ieta]>>VZhydro[itau][ix][iy][ieta];
	 // cout<<tauh<<" "<<ed<<" "<<temphydro[itau][ix][iy][ieta]<<" "<<fractionhydro[itau][ix][iy][ieta]<<" "<<VXhydro[itau][ix][iy][ieta]<<" "<<VYhydro[itau][ix][iy][ieta]<<" "<<VZhydro[itau][ix][iy][ieta]<<endl;
	  
      }
      }
      }
      }

     }
     f2.close();
     
     
void bulklinear(double tau, double x,double y,double eta, double &temp,double &VX,double &VY,double &VZ,double &fraction)
{


  if(switchmedium==0)
  {
/*	  
  double dxh=0.3;
  double dyh=0.3;
  double detah=0.3;  
  double dtauh=0.3;

  double x0h=-10.8;
  double y0h=-10.8;
  double eta0h=-4.8;  
  double tau0h=0.6;
*/

  double dxh=dxh_r;
  double dyh=dyh_r;
  double detah=detah_r;  
  double dtauh=dtauh_r;

  double x0h=x0h_r;
  double y0h=y0h_r;
  double eta0h=eta0h_r;  
  double tau0h=tau0h_r;


  
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

}

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
if(corswitch==1)
          {

            VX=Vx_cal(VX,eta,VZ);
            VY=Vy_cal(VY,eta,VZ);
            VZ=Vz_cal(eta,VZ);

          }


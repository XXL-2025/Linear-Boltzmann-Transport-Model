#include"LBT.h"

//////////////////////////////////////////
float ran33(long *idum);
void rotate(double px,double py,double pz,double pr[4],int icc);
//////////////////////////////////////////

int main(int argc, char **argv) {	
  
    //.... LBTclass. Initialization
  LBTclass *LBT = new LBTclass();
  LBT->LBTinitialize();
  
  //.... The program starts
  struct tm *local_start;
  time_t time_start;
  time_start = time(NULL);
  local_start = localtime(&time_start);

  char buf1[80];
  strftime(buf1, 80, "Current Time: %Y-%m-%d %H:%M:%S", local_start);
  cout << "the program starts at:" << endl;
  cout << buf1 << endl;
  //.... Time counts


	 char readinpp_route[1024];	

  //.... Output files in AA

  ifstream infile("./file.in");
 //   infile >> LBT->NUM1;
 //   cout << "NUM1" << endl;
 //   cout << LBT->NUM1 << endl;

    infile >> readinpp_route;
    cout << "readinpp_route" << endl;
    cout << readinpp_route << endl;
    infile.close();
  ofstream positiveLBT("./outputdatafile/positive.dat");
  ofstream negativeLBT("./outputdatafile/negative.dat");

    ofstream numpositiveLBT("./outputdatafile/numptpositive.dat");
    ofstream numnegativeLBT("./outputdatafile/numptnegative.dat");

    ofstream f55("./outputdatafile/f55.dat");
     ofstream f_position("./outputdatafile/position.dat");
    ofstream f_position_initial("./outputdatafile/position_initial.dat");

    ofstream positiveInfo("positive_info.dat");//xx-25-12-17

  //...read in initial shower information from Pythia8 or Pythia6 or Hijing    (Dijet or Gammajet or Singlejet)

	 char numpt_route[1024];
     sprintf(numpt_route,"%s/numpt.dat",readinpp_route);   
     ifstream f3(numpt_route);
	 cout<<numpt_route<<endl;

	 char parton_route[1024];
     sprintf(parton_route,"%s/parton.dat",readinpp_route);   
     ifstream f4(parton_route);
	 cout<<parton_route<<endl;

	 char gamma_route[1024];
     sprintf(gamma_route,"%s/gamma.dat",readinpp_route);   
     ifstream f5(gamma_route);
	 cout<<gamma_route<<endl;
	

   if(!f3.is_open())
    {
    cout<<"Erro openning date file3!\n";
    }
    //ifstream f4("parton.dat");
    if(!f4.is_open())
    {
    cout<<"Erro openning date file4!\n";
    }
    //ifstream f5("gamma.dat");
    if(!f5.is_open())
    {
    cout<<"Erro openning date file5!\n";
    }




  int    ndrphoton;
    int    KATTgamma;
    int    numnucleon=201*201;
	double numjettotal;
	double randomxy=0.0;
	double R1=0.0;
	double XXX=0.0;
	double YYY=0.0;
  for(int i=1;i<=numnucleon;i++)
	    {
		f_position_initial<<LBT->numjet[i]<<" "<<LBT->Xnucleon[i]<<" "<<LBT->Ynucleon[i]<<endl;
        }
    


    double etagamma;
    double phi1gamma;
    double tanphigamma;
    double phiagamma;
    double energygamma;
    double ptgamma;
    double Egamma;

    double P0gamma;
    double P1gamma;
    double P2gamma;
    double P3gamma;
	
    int numGamma = 0;









  int numEvent=0;

  for (int n = 1; ; ++n)
  {

      //cout<<"n"<<"---------------------"<<n<<endl;


      if (numEvent >= LBT->ncall) break;


      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

      //.............gammajet information

      ndrphoton = 0;

      P0gamma = 0.0;
      P1gamma = 0.0;
      P2gamma = 0.0;
      P3gamma = 0.0;

      //......read in initial Gamma 4-momentum	
      int ne;
      int IDgamma;
      double asd1, asd2;
      double PGm[4] = { 0.0 };
      double weight1;


      //f5 >> ne >> PGm[1] >> PGm[2] >> PGm[3] >> PGm[0] >> asd1 >> asd2;	


//..............0122

      if (LBT->switchsingle == 0) {


          f5 >> ne >> weight1 >> IDgamma >> PGm[1] >> PGm[2] >> PGm[3] >> PGm[0];

          ndrphoton = ndrphoton + 1;

          P0gamma = PGm[0];
          P1gamma = PGm[1];
          P2gamma = PGm[2];
          P3gamma = PGm[3];

          etagamma = 1.0 / 2.0 * log((PGm[0] + PGm[3]) / (PGm[0] - PGm[3]));

          tanphigamma = abs(P2gamma / P1gamma);
          phiagamma = atan(tanphigamma);

          if (P1gamma > 0.0 && P2gamma > 0.0)
          {
              phi1gamma = phiagamma;
          }
          if (P1gamma < 0.0 && P2gamma>0.0)
          {
              phi1gamma = 3.1415926 - phiagamma;
          }
          if (P1gamma < 0.0 && P2gamma < 0.0)
          {
              phi1gamma = 3.1415926 + phiagamma;
          }
          if (P1gamma > 0.0 && P2gamma < 0.0)
          {
              phi1gamma = 2 * 3.1415926 - phiagamma;
          }

          energygamma = P0gamma;
          ptgamma = sqrt(pow(P1gamma, 2) + pow(P2gamma, 2));

          cout << "    ptgamma   " << ptgamma << "   etagamma  " << etagamma << endl;

      }//..............0122

//...................................................................................................................................isolation cut


//......CMS gamma-jet cut		
        //if(ndrphoton==0 || ptgamma<ptcutgammamin || ptgamma>ptcutgammamax || abs(etagamma)>etacutgamma) continue;		

//        numGamma=numGamma+1;




//......fastjet vector initialize
      int ui = -123456;
      int npart = 100000;

      vector<double> px(npart), py(npart), pz(npart), E(npart);
      vector<int> pid(npart);
      //......fastjet vector initialize end


      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

              //numEvent++; //put in wherever needed

      LBT->LBTclear();  //clear particle list in LBT

      //......Geometry sampling

      if (LBT->switchsingle == 0) {  //..............0122

          randomxy = LBT->ran0(&LBT->NUM1);

          R1 = 0.0;



          for (int i = 1; i <= numnucleon; i++)
          {
              R1 = R1 + LBT->numjet[i] / LBT->numjettotal;
              //cout<<"R1:"<<R1<<endl;
              if (randomxy < R1) {
                  XXX = LBT->Xnucleon[i];
                  YYY = LBT->Ynucleon[i];
                  break;
              }
          }

          //XXX=0.0;
          //YYY=0.0;

          cout << "XXX" << " " << XXX << " " << "YYY" << " " << YYY << endl;

          f_position << n << " " << XXX << " " << YYY << endl;

      }
      //.................................................................................803		

      double deltaR;

      double tanphiparton;
      double phiaparton;
      double phi1parton;

      double energycone = 0.0;

      //.................................................................................803		



      //......Parton information read in

              //int ne;		
      int njetshower;
      double energy_weight;
      int fjetshower;
      double weight;
      double leading_pT = 0.0;
      int ileading;

      int Aleading;

      double timeplus, timeplus_0, timeplus_0_1;

      double pythiap0;

      int nj = 0;

      double pxjetshower, pyjetshower, pzjetshower, p0jetshower, vxjetshower, vyjetshower, vzjetshower, vtjetshower, etajetshower, ptjetshower,jetcolor,jetanticolor; //??-1-11


      if (LBT->switchsingle == 0) {  //..............0122


          f3 >> ne >> njetshower >> energy_weight;

          //cout<<"ne"<<" "<<ne<<" "<<"njetshower"<<" "<<njetshower<<endl;

          //cout<<"LBT->nj"<<" "<<LBT->nj<<" "<<"nj"<<" "<<nj<<endl;		

          for (unsigned ip = 1; ip <= njetshower; ++ip)
          {

              //f4 >> ne >> fjetshower >> pxjetshower >> pyjetshower >> pzjetshower >> pythiap0 >> vxjetshower >> vyjetshower >> vzjetshower >> vtjetshower >> Aleading >> timeplus >> timeplus_0 >> timeplus_0_1;			

             // f4 >> ne >> fjetshower >> pxjetshower >> pyjetshower >> pzjetshower >> pythiap0 >> Aleading >> vxjetshower >> vyjetshower >> vzjetshower >> vtjetshower;	
              f4 >> ne >> fjetshower >> weight >> pxjetshower >> pyjetshower >> pzjetshower >> pythiap0 >> Aleading >> vxjetshower >> vyjetshower >> vzjetshower >> vtjetshower>> jetcolor>> jetanticolor;//??-1-11

              if (abs(fjetshower) == 4 || abs(fjetshower) == 5) { p0jetshower = pythiap0; } //??-1-11
              else { //xx-25-12-14

                  p0jetshower = sqrt(pow(pxjetshower, 2) + pow(pyjetshower, 2) + pow(pzjetshower, 2));
              }   //current calculation with light partons assume that they are massless 

              //p0jetshower=pythiap0;

              etajetshower = 1.0 / 2.0 * log((p0jetshower + pzjetshower) / (p0jetshower - pzjetshower));

              ptjetshower = sqrt(pow(pxjetshower, 2) + pow(pyjetshower, 2));

              //cout<<"ip"<<" "<<ip<<" "<<"fjetshower"<<" "<<fjetshower<<endl;

              //if(abs(Aleading) == 23) continue;

              if (abs(fjetshower) != 21 && abs(fjetshower) != 1 && abs(fjetshower) != 2 && abs(fjetshower) != 3 && abs(fjetshower) != 4 && abs(fjetshower) != 5) continue;   //??-1-11


              //.................................................................................803			

              tanphiparton = abs(pyjetshower / pxjetshower);
              phiaparton = atan(tanphiparton);

              if (pxjetshower > 0.0 && pyjetshower > 0.0)
              {
                  phi1parton = phiaparton;
              }
              if (pxjetshower < 0.0 && pyjetshower>0.0)
              {
                  phi1parton = 3.1415926 - phiaparton;
              }
              if (pxjetshower < 0.0 && pyjetshower < 0.0)
              {
                  phi1parton = 3.1415926 + phiaparton;
              }
              if (pxjetshower > 0.0 && pyjetshower < 0.0)
              {
                  phi1parton = 2 * 3.1415926 - phiaparton;
              }

              deltaR = sqrt(pow(abs(etajetshower - etagamma), 2) + pow(abs(phi1parton - phi1gamma), 2));

              if (abs(deltaR) < 0.4)
              {

                  energycone = energycone + p0jetshower;

              } //if(abs(deltaR)<0.4)


//.................................................................................803

            //cout<<ip<<" "<<nj<<endl;
            //cout<<ip<<" "<<etajetshower<<" "<<ptjetshower<<endl;

              if (abs(etajetshower) > 4.8) continue;     //track cut

              //if(ptjetshower<0.5) continue;     //track cut

              nj = nj + 1;

              LBT->KATT1[nj] = fjetshower;

              LBT->P[1][nj] = pxjetshower;
              LBT->P[2][nj] = pyjetshower;
              LBT->P[3][nj] = pzjetshower;
              LBT->P[0][nj] = p0jetshower;
              LBT->color[nj] = jetcolor; //??-1-11
              LBT->anticolor[nj] = jetanticolor;//??-1-11
              if (abs(LBT->color[nj]) > LBT->maxColor) LBT->maxColor = abs(LBT->color[nj]);//??-1-12
              if (abs(LBT->anticolor[nj] > LBT->maxColor)) LBT->maxColor = abs(LBT->anticolor[nj]);//??-1-12

              if (ptjetshower > leading_pT)
              {
                  ileading = nj;
              }


              //nj=1;

              //LBT->KATT1[nj]=1;

              //LBT->P[1][nj]=-P1gamma;
              //LBT->P[2][nj]=-P2gamma;
              //LBT->P[3][nj]=-P3gamma;
              //LBT->P[0][nj]=P0gamma;


  /*
              nj=1;

              LBT->KATT1[nj]=1;

              LBT->P[1][nj]=0.0;
              LBT->P[2][nj]=80.0;
              LBT->P[3][nj]=0.0;
              LBT->P[0][nj]=80.0;
  */


              LBT->tjp[nj] = vtjetshower;

              //cout << ne << " " << njetshower << " " << nj <<" " << LBT->KATT1[nj] << " " << LBT->P[1][nj] << " " << LBT->P[2][nj] << " " << LBT->P[3][nj] << " " << LBT->P[0][nj] << " " << etajetshower << endl;

  //..........formation time of the jet parton		  
              //LBT->tiform[nj]=2.0*p0jetshower/(pow(pxjetshower,2)+pow(pyjetshower,2));

              //LBT->tiform[nj]=2.0*LBT->P[0][nj]/(pow(LBT->P[1][nj],2)+pow(LBT->P[2][nj],2));
              //LBT->tiform[nj]=timeplus;
              LBT->tiform[nj] = 0.0;

              LBT->V[1][nj] = XXX;
              LBT->V[2][nj] = YYY;
              LBT->V[3][nj] = 0.0;

              LBT->tirad[nj] = LBT->tiform[nj];

              LBT->nj = nj;

          }//for(unsigned ip=1; ip<=njetshower; ++ip)



      }  //..............0122

//......CMS gamma-jet cut		
        //if(ndrphoton==0 || ptgamma<ptcutgammamin || ptgamma>ptcutgammamax || abs(etagamma)>etacutgamma) continue;	//..............0122	???

//		cout<<"energycone"<<" "<<energycone<<endl;

        //if(energycone>5.0) continue;		

      numGamma = numGamma + 1;

      //.................................................................................803

              //h100Ngamma->Fill(ptgamma,ww10);

      //.................................................................................803







      if (LBT->switchsingle == 1) {  //..............0122

          LBT->nj = 1;
          LBT->KATT1[1] = 5;
          LBT->P[1][1] = 0;
          LBT->P[3][1] = 0;
          LBT->P[0][1] = 100;
        //  LBT->P[2][1] = 100;
//LBT->P[2][1]=sqrt(LBT->P[0][1]*LBT->P[0][1] - 1.27 * 1.27);
LBT->P[2][1] = 99.91218;
LBT->V[0][1] = 0.6;
LBT->V[1][1] = 0;
LBT->V[2][1] = 0;
LBT->V[3][1] = 0;
LBT->color[1] = 0;//??-1-11
LBT->anticolor[1] = 110;//??-1-11

if (abs(LBT->color[1]) > LBT->maxColor) LBT->maxColor = abs(LBT->color[1]);//??-1-12
if (abs(LBT->anticolor[1] > LBT->maxColor)) LBT->maxColor = abs(LBT->anticolor[1]);//??-1-12


      }



      //cout<<LBT->nj<<"===================="<<nj<<endl;		

      if (LBT->switchsingle == 0) {
          LBT->nj = nj;
      }

      //cout<<LBT->nj<<"===================="<<nj<<endl;		
      //exit(1);
/*
        if(LBT->switchsingle==1)
        {
        LBT->LBTclear();  //clear particle list in LBT
        }
*/

      double epsilon = 0.000000001;

      //......time evolution in LBT

      double ti = LBT->time0;	//initial time inside or outside class?

      LBT->np = LBT->nj;


      int ntimestep = floor(LBT->timend / LBT->dt);

      if (LBT->LBTswitch == 0)
      {
          ntimestep = 1;
      }


      //......test		
      LBT->ntest22 = 0;
      LBT->ntestrad = 0;
      //......test
      cout << "ntimestep:" << ntimestep << endl;

      // xx-25-12-17
      for (int i = 1; i <= LBT->np; i++) {
          if (LBT->P[0][i] == 0) continue;
          // ��������
          double mass = 0.0;
          if (abs(LBT->KATT1[i]) == 4) {
              mass = 1.27; // c�������
          }
          else if (abs(LBT->KATT1[i]) == 5) { //??-1-12
              mass = 4.19; 
          }
          else {
              mass = sqrt(LBT->P[0][i] * LBT->P[0][i] - LBT->P[1][i] * LBT->P[1][i] -
                  LBT->P[2][i] * LBT->P[2][i] - LBT->P[3][i] * LBT->P[3][i]);
          }
          if (i == 1) {
              positiveInfo << "INIT " << n << " " << LBT->time0 << " " << LBT->KATT1[i] << "    "
                  << LBT->P[1][i] << "    " << LBT->P[2][i] << "    " << LBT->P[3][i] << "    "
                  << LBT->P[0][i] << "    " << mass << "    "
                  << LBT->V[1][i] << "    " << LBT->V[2][i] << "    " << LBT->V[3][i] << "    "
                  << LBT->V[0][i] << "    " << 0.0 << "    " << 0.0 << "    " << 0.0 << "    "
                  << 0.0 << "    " << sqrt(LBT->P[1][i] * LBT->P[1][i] + LBT->P[2][i] * LBT->P[2][i]) << "    "
                  << 0.0 << endl;
          }
      }
    for(int timestep=1;timestep<=ntimestep;++timestep) 
    {
      if(LBT->LBTswitch==1) 
      {		
        
        ti=ti+LBT->dt;
        cout << "ti:" << ti << endl;
        LBT->LinearBoltzmannTransport(n,ti);

        for (int i = 1; i <= LBT->np; i++) {  //xx-25-12-17
            if (LBT->P[0][i] == 0) continue;
            // ��������
            double mass = 0.0;
            if (abs(LBT->KATT1[i]) == 4) {
                mass = 1.27; // c�������
            }
            else if (abs(LBT->KATT1[i]) == 5) { //??-1-12
                mass = 4.19;
            }
            else {
                mass = sqrt(LBT->P[0][i] * LBT->P[0][i] - LBT->P[1][i] * LBT->P[1][i] -
                    LBT->P[2][i] * LBT->P[2][i] - LBT->P[3][i] * LBT->P[3][i]);
            }
            if (i == 1) {
                positiveInfo << "EVOL " << n << " " << ti << " " << LBT->KATT1[i] << "    "
                    << LBT->P[1][i] << "    " << LBT->P[2][i] << "    " << LBT->P[3][i] << "    "
                    << LBT->P[0][i] << "    " << mass << "    "
                    << LBT->Vfrozen[1][i] << "    " << LBT->Vfrozen[2][i] << "    " << LBT->Vfrozen[3][i] << "    "
                    << LBT->Vfrozen[0][i] << "    " << 0.0 << "    " << 0.0 << "    " << 0.0 << "    "
                    << 0.0 << "    " << sqrt(LBT->P[1][i] * LBT->P[1][i] + LBT->P[2][i] * LBT->P[2][i]) << "    "
                    << 0.0 << endl;
            }
        }


        if(LBT->Reachtauend==1 && LBT->switchmedium==0)  //!!! for uniform test  0122???
		{
		timestep=ntimestep;	
		}	

    }
     //.... time loop end
            int switchoutput=1;
            if(switchoutput>0)
            {
              if(timestep==ntimestep)
              {
                numEvent = numEvent + 1;

            //cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

                f55 << n <<" "<< PGm[1] <<" "<< PGm[2] <<" "<< PGm[3] <<" "<< PGm[0] <<" "<< endl;

    int ip = 0;
    for(int i=1;i<=LBT->np;i++)
    {
      if(LBT->P[0][i] == 0) ip+=1;
    }	
    int nnn=LBT->np-ip;
    numpositiveLBT<<n<<" "<<nnn<<" "<<LBT->nj<<endl;
    for(int i=1;i<=LBT->np;i++)
    {
      if(LBT->P[0][i] == 0) continue;
       positiveLBT<<n<<" "<<LBT->KATT1[i]<<" "<<LBT->P[1][i]<<" "<<LBT->P[2][i]<<" "<<LBT->P[3][i]<<" "<<LBT->P[0][i] <<" "<<LBT->Vfrozen[1][i]<<" "<<LBT->Vfrozen[2][i]<<" "<<LBT->Vfrozen[3][i]<<" "<<LBT->Vfrozen[0][i]<<" "<<LBT->CAT[i]<<endl;
      // positiveLBT<<LBT->KATT1[i]<<" "<<LBT->P[1][i]<<" "<<LBT->P[2][i]<<" "<<LBT->P[3][i]<<" "<<LBT->P[0][i]<<" "<<LBT->CAT[i]<<endl;
      //positiveLBT<<LBT->KATT1[i]<<" "<<LBT->P[1][i]<<" "<<LBT->P[2][i]<<" "<<LBT->P[3][i]<<" "<<LBT->P[0][i]<<" "<<LBT->Vfrozen[1][i]<<" "<<LBT->Vfrozen[2][i]<<" "<<LBT->Vfrozen[3][i]<<" "<<LBT->Vfrozen[0][i]<<" "<<LBT->CAT[i]<<endl;
    }	
   

    int ip0=0;
    for(int i=LBT->nj+1;i<=LBT->np;i++)
    {
      if(LBT->P0[0][i] == 0) ip0+=1;
    }	
    int nnn0=LBT->np-ip0-LBT->nj;
     numnegativeLBT<<n<<" "<<nnn0<<endl;	
    for(int i=LBT->nj+1;i<=LBT->np;i++)
    {
      if(LBT->P0[0][i] == 0) continue;
      negativeLBT<<n<<" "<<LBT->KATT10[i]<<" "<<-LBT->P0[1][i]<<" "<<-LBT->P0[2][i]<<" "<< -LBT->P0[3][i]<<" "<<-LBT->P0[0][i]<<" "<<LBT->Vfrozen0[1][i]<<" "<<LBT->Vfrozen0[2][i]<<" "<<LBT->Vfrozen0[3][i]<<" "<<LBT->Vfrozen0[0][i]<<" "<<LBT->CAT0[i]<<endl;
      //negativeLBT<<LBT->KATT10[i]<<" "<<-LBT->P0[1][i]<<" "<<-LBT->P0[2][i]<<" "<< -LBT->P0[3][i]<<" "<<-LBT->P0[0][i]<<" "<<endl;
      //negativeLBT<<LBT->KATT10[i]<<" "<<-LBT->P0[1][i]<<" "<<-LBT->P0[2][i]<<" "<< -LBT->P0[3][i]<<" "<<-LBT->P0[0][i]<<" "<<LBT->Vfrozen0[1][i]<<" "<<LBT->Vfrozen0[2][i]<<" "<<LBT->Vfrozen0[3][i]<<" "<<LBT->Vfrozen0[0][i]<<endl;
    }

     }//if(timestep==ntimestep)		
		
		    }//if(switchoutput>0)

    if(LBT->nj==0) continue;
  
    }//.... time evolution end


   // cout<<"numEvent"<<" "<<numEvent<<" "<<"numGamma"<<" "<<numGamma<<" "<<"n"<<" "<<n<<endl;


    if(numEvent>=LBT->ncall) break;



    //.... Output control
    int print=n%LBT->nprint;
	    //if(print==0)
    {
    //  cout << "n" << "    " << "ntest22" << "    " << "ntestrad" << endl;      
     // cout << n << "  " << LBT->ntest22 << " " << LBT->ntestrad << endl;	

     // cout << "n" << "    " << "np" << "    " << "LBT->nj" << endl;      
     // cout << n << "  " << LBT->np << " " << LBT->nj << endl;			

      //exit(1);
    }

  }
  //.... Event loop end



  

  positiveLBT.close();	
  negativeLBT.close();	

  numpositiveLBT.close();
  numnegativeLBT.close();

  f55.close();
	
    f3.close();
    f4.close();
    f5.close();
    positiveInfo.close(); //xx-25-12-17
  //.... Time counts
  struct tm *local_end;
  time_t time_end;
  time_end = time(NULL);
  local_end = localtime(&time_end);

  char buf2[80];
  strftime(buf2, 80, "Current Time: %Y-%m-%d %H:%M:%S", local_end);
  cout << "the program ends at:" << endl;
  cout << buf2 << endl;

  int cost, nh, nm, ns;
  cost = difftime(time_end, time_start);

  nh = cost / 3600;
  nm = (cost % 3600) / 60;
  ns = (cost % 3600) % 60;

  cout << "the program costs:" << endl;
  cout << cost << "s:" << " " << nh << "h" << " " << nm << "m" << " " << ns << "s" << endl;
  //.... The program ends
}




#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran33(long *idum)

{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0) { 
        if (-(*idum) < 1) *idum=1; 
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) { 
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1; 
    *idum=IA1*(*idum-k*IQ1)-k*IR1; 
    if (*idum < 0) *idum += IM1; 
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2; 
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV; 
    iy=iv[j]-idum2; 
    iv[j] = *idum; 
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX; 
    else return temp;
}


void rotate(double px,double py,double pz,double pr[4],int icc){
    //     input:  (px,py,pz), (wx,wy,wz), argument (i)
    //     output: new (wx,wy,wz)
    //     if i=1, turn (wx,wy,wz) in the direction (px,py,pz)=>(0,0,E)
    //     if i=-1, turn (wx,wy,wz) in the direction (0,0,E)=>(px,py,pz)


    double wx,wy,wz,E,pt,w,cosa,sina,cosb,sinb;
    double wx1,wy1,wz1;	   	   

    wx=pr[1];
    wy=pr[2];
    wz=pr[3];

    E=sqrt(px*px+py*py+pz*pz);
    pt=sqrt(px*px+py*py);

    w=sqrt(wx*wx+wy*wy+wz*wz);

    //  if(pt==0)
    if(pt<1e-6)
    {
        cosa=1;
        sina=0;
    } 
    else
    {
        cosa=px/pt;
        sina=py/pt;
    }

    if(E>1e-6) {

        cosb=pz/E;
        sinb=pt/E;

        if(icc==1) {
            wx1=wx*cosb*cosa+wy*cosb*sina-wz*sinb;
            wy1=-wx*sina+wy*cosa;
            wz1=wx*sinb*cosa+wy*sinb*sina+wz*cosb;
        } else {
            wx1=wx*cosa*cosb-wy*sina+wz*cosa*sinb;
            wy1=wx*sina*cosb+wy*cosa+wz*sina*sinb;
            wz1=-wx*sinb+wz*cosb;
        }
        wx=wx1;
        wy=wy1;
        wz=wz1;
    } else {
        cout << "warning: small E in rotation" << endl;
    }

    pr[1]=wx;
    pr[2]=wy;
    pr[3]=wz;      

    //  pr[0]=sqrt(pr[1]*pr[1]+pr[2]*pr[2]+pr[3]*pr[3]);

}


